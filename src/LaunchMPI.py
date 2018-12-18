import sys
import imp
import pyipopt
import logging
import time
import support
import os
import numpy as np
from mpi4py import MPI

import cantera as ct

import Global
import Parallelism as Par
import Slave
import Parameters as Param
import OptimFunctions as OptFn
import ReactorCstVolume as RCV
import OneDFlame as Flame


def eval_jac_g_wrapper(x, flag, user_data=None):
    noptim = Global.params.noptim
    if flag:
        logging.debug('Jacobian Initialization')
        tab = []
        for itcon in range(ncon):
            for k in range(noptim):
                tab.append(itcon)
        tab2 = []
        for itcon in range(ncon):
            for k in range(noptim):
                tab2.append(k)
        return (np.array(tab), np.array(tab2))
    else:
        return OptFn.gradient_constraints(x)


if __name__ == '__main__':
    # Main process, common to all tasks (master and slaves)
    # Reading input
    up = imp.load_source('user_param', sys.argv[1])

    # Cantera gas phase
    Global.gas = ct.Solution(up.mechanism)
    Global.valid_species = [line.rstrip(' \n') for line in open(up.spec_list)]

    # Define MPI message tags
    tags = Par.enum('READY', 'DONE', 'EXIT', 'START', 'SLEEP', 'WAKEUP')

    # Initializations and preliminaries
    Global.MPI = MPI
    # Maximum ignition delay time
    Global.tMax = 0.1
    comm = MPI.COMM_WORLD   # get MPI communicator object
    size = comm.size        # total number of processes
    rank = comm.rank        # rank of this process
    status = MPI.Status()   # get MPI status object

    Global.partitionid = -1

    if rank > 0:
        # Slave Process (waiting for task from Master)
        Slave.slave_execution()
    else:
        print "Rank is : ", rank
        # Main Process
        if not os.path.exists(up.directory):
            os.makedirs(up.directory)

        logfile = 'out.log'
        verbosity = 10
        if hasattr(up, 'verbosity'):
            if up.verbosity == 'DEBUG':
                verbosity = 10
            elif up.verbosity == 'INFO':
                verbosity = 20
        logging.basicConfig(filename=logfile, level=verbosity)
        
        logging.info('Start')
        # Initializing gas phase and damping
        gas = Global.gas

        # Load list of valid species
        valid_species = Global.valid_species
        print valid_species
        num_valid = len(valid_species)
        print num_valid

        logging.info('Species for optimization: %s' % str(num_valid))
        # Last species assumed to remove equality constraint
        noptim = num_valid - 1
        logging.info('Number of variables for optimization: %s' % noptim)

        # Species index build
        # species_index_exclude_init = []
        # for species in up.species_exclude_init:
        #     species_index_exclude_init.append(gas.species_index(species))

        # # species_index_exclude_zero = []
        # # for species in up.species_exclude_zero:
        # #     species_index_exclude_zero.append(gas.species_index(species))

        # # Create excluded species list
        # species_exclude_all = up.species_exclude_init 
        # species_index_exclude_all = []
        # for species in species_exclude_all:
        #     species_index_exclude_all.append(gas.species_index(species))

        # species_damp_init = ()
        # for species in gas.species_names:
        #     if not(species in up.species_exclude_init) and not(species in up.species_major):
        #         species_damp_init = species_damp_init + (species,)
        # species_index_damp_init = []
        # for species in species_damp_init:
        #     species_index_damp_init.append(gas.species_index(species))

        # species_damp = ()
        # for species in gas.species_names:
        #     if not(species in species_exclude_all) and not(species in up.species_major):
        #         species_damp = species_damp + (species,)
        # species_index_damp = []
        # for species in species_damp:
        #     species_index_damp.append(gas.species_index(species))

        cases = []
        quantityrefs = []
        tolerances = []

        # Group procs in for flames
        procs = range(0, comm.size-1)
        partitioning = []
        counting = -1
        for k in range(len(procs)):
          partitioning.append([-1])

        logging.debug('Partitioning: %s' % str(partitioning))

        Global.partitioning = partitioning

        # Read reference composition
        test_comp = up.test_comp
        palette = up.palette       
 
        # Lightweight constraints reference quantities
        num_lwcon = 3
        # MW
        mw_vec = []
        mw_valid_vec = []
        mw_full_vec = gas.molecular_weights
        # Reference MW
        for species in palette:
            sp_index = gas.species_index(species)
            mw_vec.append(mw_full_vec[sp_index])
        quantityrefs.append(np.dot(mw_vec,test_comp))
        # Valid species MW
        for species in valid_species:
            sp_index = gas.species_index(species)
            mw_valid_vec.append(mw_full_vec[sp_index])
        tolerances.append(up.tolerance_mw)
        print "Constructed MW" 

        # HC
        h_vec = []
        c_vec = []
        h_valid_vec = []
        c_valid_vec = []
        # Reference HC
        for species in palette:
            sp_index = gas.species_index(species)
            h_vec.append(gas.n_atoms(sp_index,'H'))
            c_vec.append(gas.n_atoms(sp_index,'C'))
        target_hc = np.dot(h_vec,test_comp)/np.dot(c_vec,test_comp)
        quantityrefs.append(target_hc)
        # Valid species HC
        for species in valid_species:
            sp_index = gas.species_index(species)
            h_valid_vec.append(gas.n_atoms(sp_index,'H'))
            c_valid_vec.append(gas.n_atoms(sp_index,'C'))
        tolerances.append(up.tolerance_hc)
        print "Constructed HC" 

        # Sum constraint
        quantityrefs.append(0.0)
        tolerances.append(0.0)

        # Constructing AI cases
        case_id = 0
        nai = len(up.P_ai)*len(up.T_ai)*len(up.phi_ai)
        for P in up.P_ai:
            for T in up.T_ai:
                for phi in up.phi_ai:
                    case_id += 1
                    cur = RCV.ReactorCstVolume(
                        Global.gas, T, P, phi, up.fuel, up.n2_o2_ratio)
                    cur.case_id = case_id
                    cur.isflame = False
                    cases.append(cur)
                    tolerances.append(up.tolerance_ai)
        print "Constructed AI cases" 

        # Computing reference cases
        tasks = []
        taskindex = -1
        print "Cases: ", cases

        # Construct reference x
        x0 = np.zeros(len(Global.valid_species))

        for case in cases:
            taskindex += 1
            # TODO : Modify this for new problem
            tasks.append((case, x0, taskindex))
        print "Beginning launched"
        quantityrefs_new = np.array(Par.tasklaunch(tasks))
        quantityrefs.append(quantityrefs_new[0])
        print "Computed reference cases"

        # RHS of constraint inequality
        tolerances = np.array(quantityrefs) * np.array(tolerances)
        tolerances[-1] = 1.0

        # TODO remove these global variables
        params = Param.Parameters()
        params.cases = cases
        # params.species_index_damp = species_index_damp
        # params.species_index_exclude_all = species_index_exclude_all
        params.num_valid = num_valid
        params.mw_valid_vec = mw_valid_vec
        params.h_valid_vec = h_valid_vec
        params.c_valid_vec = c_valid_vec
        params.noptim = noptim
        params.num_lwcon = num_lwcon
        params.quantityrefs = np.array(quantityrefs)
        params.tolerances = tolerances
        # params.species_damp = species_damp
        params.threshold = up.threshold
        Global.params = params

        print 'Number of cases: %i' % len(cases)
        print 'Reference quantities: %s' % quantityrefs
        logging.info('Number of cases: %i' % len(cases))
        logging.info('Reference quantities: %s' % quantityrefs)

        # Setting up IPOpt parameters
        nvar = noptim  # Number of variables
        x_L = np.ones(nvar)*0.0  # Vector of lower bounds
        x_U = np.ones(nvar)*1.0  # Vector of upper bounds
        ncon = len(cases) + num_lwcon  # Number of constraints
        logging.info('Number of constraints: %i' % ncon)

        # Building lower and upper bound vectors of tolerance
        # Dummy value since constraints are positive 
        g_L = -5.0 * np.ones(ncon)
        g_U = tolerances
        logging.info('Lower bounds: %s' % g_L)
        logging.info('Upper bounds: %s' % g_U)

        # Non-zero Jacobian elements
        nnzj = ncon * noptim
        # Non-zero Hessian (set to 0 since not evaluated)
        nnzh = 0

        # Initial coefficient vector
        if hasattr(up, 'x0'):
            x0 = np.array(up.x0)
            logging.info('Initializing x0 from input file')
        else:
            x0 = support.vec_from_palette(valid_species,palette,test_comp)
            print "x0:", x0
            logging.info('Initializing x0 from palette')

        if not(len(x0) == noptim):
            logging.error('Wrong length for x0')
            sys.exit(0)

        # TODO : Log initial x0
        # logging.info('Initial x0 vector')
        # for k in range(len(species_damp)):
        #     logging.info('%s %s' %
        #                  (gas.species_names[species_index_damp[k]], x0[k]))
        # logging.info('x0 passed initial: %s' % str(x0))

        type_calc = 'OPT'
        if hasattr(up, 'type_calc'):
            type_calc = up.type_calc

        if type_calc == 'OPT':
            # Optimization case
            logging.debug('IPOpt initialization')
            nlp = pyipopt.create(nvar, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh,
                                 OptFn.eval_f, OptFn.eval_grad_f, OptFn.eval_g, eval_jac_g_wrapper)

            #nlp.str_option('derivative_test','first-order')
            #nlp.num_option('derivative_test_tol',1e-2)

            logging.info('Starting optimization')
            x, zl, zu, constraint_multipliers, obj, status = nlp.solve(x0)
            nlp.close()
            logging.info('End optimization')

            # Threshold
            x_last = 1.0 - np.sum(x)
            x = np.append(x,x_last)
            x_unthreshold = x
            beta = np.array([0 if y < up.threshold else 1 for y in x])
            x = beta*x
  
            # MW
            mw = np.sum(mw_valid_vec*x)
            logging.info('Final MW: %s' % mw) 
            
            # HC
            n_h = np.sum(h_valid_vec*x)
            n_c = np.sum(c_valid_vec*x)
            hc = n_h/n_c
            logging.info('Final HC: %s' % hc) 

            # Sum = 1
            x_sum = np.sum(x)
            logging.info('Final sum: %s' % x_sum) 

            # Log final solution
            for k in range(num_valid):
                logging.info('%s %s %s' %
                            (valid_species[k], beta[k], x_unthreshold[k]))

        elif type_calc == 'VAL':
            # Validation case
            logging.info('Validation calculation')

            x1 = np.array(x0)
            xreduced = []
            zero_species = []
            remaining_species = []

            for k in range(len(species_damp)):
                if (x1[k] < params.threshold):
                    x1[k] = 0.0
                    zero_species.append(
                        gas.species_names[species_index_damp[k]])
                else:
                    x1[k] = 1.0
                    xreduced.append(x1[k])
                    remaining_species.append(
                        gas.species_names[species_index_damp[k]])
            logging.info('Eliminated species: %s' % str(zero_species))
            for k, species in enumerate(species_damp):
                if species in zero_species:
                    logging.info('%s %s' % (species, x0[k]))
            logging.info('Remaining species: %s' % str(remaining_species))
            for k, species in enumerate(species_damp):
                if species in remaining_species:
                    logging.info('%s %s' % (species, x0[k]))

            # Computing reference quantities
            tasks = []
            taskindex = -1
            for case in cases:
                taskindex += 1
                tasks.append((case, species_index_exclude_init, species_index_damp_init, np.ones(
                    len(species_index_damp_init)), taskindex, 'detailed'))

            for case in cases:
                taskindex += 1
                tasks.append((case, species_index_exclude_all,
                              species_index_damp, x1, taskindex, 'reduced_zero'))

            results = np.array(Par.tasklaunch(tasks))
            refs = results[0:len(cases)]
            optimized = results[len(cases):2*len(cases)]

            errortab = np.abs((optimized-refs)/refs)
            logging.info('Reference quantities: %s' % str(refs))
            logging.info('Errors: %s' % str(errortab))
            if Global.nflames > 0:
                logging.info('Maximum error flame speed: %s' %
                             np.max(errortab[0:len(flametab)]))
            if nai > 0:
                logging.info('Maximum error AI delay: %s' %
                             np.max(errortab[len(flametab):]))
            logging.info('Maximum error: %s' % np.max((optimized-refs)/refs))

        elif type_calc == 'SA':
            # Sensivity analysis case
            tasks = []
            taskindex = -1
            perturb = 0.01
            nvars = len(species_index_damp_init)
            for k in range(nvars):
                x0 = np.ones(nvars)
                x0[k] = 1.0-perturb
                for case in cases:
                    taskindex += 1
                    tasks.append((case, species_index_exclude_init,
                                  species_index_damp_init, x0, taskindex))

            results = np.array(Par.tasklaunch(tasks))
            results = results.reshape((nvars, len(cases)))
            results = np.abs(
                (results[:, :] - quantityrefs[:])/quantityrefs[:]) / perturb
            maxi = np.max(results, axis=1)

            sorting = np.argsort(maxi)
            for indsort in sorting:
                logging.info('%s %s' %
                             (gas.species_names[species_index_damp_init[indsort]], maxi[indsort]))

            for currentk in range(1, len(sorting)+1):
                toremove = sorting[0:currentk]
                x1 = np.ones(nvars)
                x1[toremove] = 0.0

                tasks = []
                taskindex = -1
                for case in cases:
                    taskindex += 1
                    tasks.append((case, species_index_exclude_init,
                                  species_index_damp_init, x1, taskindex))

                results = np.array(Par.tasklaunch(tasks))
                refs = quantityrefs
                logging.info('Last species removed: %s' %
                             gas.species_names[species_index_damp_init[toremove[-1]]])
                logging.info('Total error with %s removed species: %s' %
                             (currentk, np.max(np.abs((results-refs)/refs))))
                logging.info('Errors: %s' % str(np.abs((results-refs)/refs)))
