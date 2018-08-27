import numpy as np
import Global
import logging
import sys
import Parallelism as Par


def getQuantityCase(case, species_index_valid, x):
    case.storeMultiplier(species_index_valid, x)
    return case.getQuantity()

def eval_f(xf, user_data=None):
    logging.debug('F evaluation')
    noptim = Global.params.noptim
    x = np.array(xf)
    costfunc = 2.0 * np.sum(x*(1-x)) + np.sum(x*x)
    return costfunc

def eval_grad_f(xf, user_data=None):
    logging.debug('Grad F evaluation')
    noptim = Global.params.noptim
    x = np.array(xf)
    grad_f = 2.0 * (np.ones(noptim)-2*x) + 2*x
    return grad_f

def eval_g(xf, user_data=None):
    logging.debug('G evaluation')
    # Retrieve global parameters
    # species_index_damp = Global.params.species_index_damp
    # species_index_exclude = Global.params.species_index_exclude_all
    # species_index_valid = Global.params.species_index_valid    
    mw_valid_vec = Global.params.mw_valid_vec
    h_valid_vec = Global.params.h_valid_vec
    c_valid_vec = Global.params.c_valid_vec
    num_valid = Global.params.num_valid

    quantities = []
    quantityrefs = Global.params.quantityrefs
  
    # Lightweight constraints
    x = np.array(xf)
    x_last = 1.0 - np.sum(x)
    x = np.append(x,x_last)
    # MW
    mw = np.sum(mw_valid_vec*x)
    quantities.append(mw)

    # HC
    n_h = np.sum(h_valid_vec*x)
    n_c = np.sum(c_valid_vec*x)
    hc = n_h/n_c
    quantities.append(hc)

    # Sum constraint
    quantities.append(np.sum(x))

    # # Asynchronous constraints
    # cases = Global.params.cases
    # quantityrefs = Global.params.quantityrefs

    # # Building Tasks
    # tasks = []
    # taskindex = -1
    # for case in cases:
    #     taskindex += 1
    #     tasks.append((case, species_index_valid, x, taskindex))

    # # Retrieve array of results
    # quantities = np.array(Par.tasklaunch(tasks))
    return np.abs(quantities - quantityrefs)

def gradient_constraints(xf):
    logging.debug('Grad G evaluation')
    # Retrieve global variables
    noptim = Global.params.noptim
    cases = Global.params.cases
    # species_index_damp = Global.params.species_index_damp
    # species_index_exclude = Global.params.species_index_exclude_all
    # species_index_valid = Global.params.species_index_valid    
    quantityrefs = Global.params.quantityrefs
    tolerances = Global.params.tolerances
    # species_damp = Global.params.species_damp
    mw_valid_vec = Global.params.mw_valid_vec
    h_valid_vec = Global.params.h_valid_vec
    c_valid_vec = Global.params.c_valid_vec
    num_valid = Global.params.num_valid
    num_lwcon = Global.params.num_lwcon

    # Grad G for lightweight constraints
    # TODO : Do this analytically

    # Compute unperturbed
    cur_res = []

    # Lightweight constraints
    x = np.array(xf)
    x_last = 1.0 - np.sum(x)
    x = np.append(x,x_last)
    # MW
    mw = np.sum(mw_valid_vec*x)
    cur_res.append(mw)

    # HC
    n_h = np.sum(h_valid_vec*x)
    n_c = np.sum(c_valid_vec*x)
    hc = n_h/n_c
    cur_res.append(hc)

    # Sum constraint
    cur_res.append(np.sum(x))

    initquantities = cur_res

    # Compute perturbations
    perturbquantities = np.empty((num_lwcon,0),float)
    step = 1e-2
    results = []
    for kx in range(noptim):
        xperturb = np.array(xf)
        xperturb[kx] += step
        x_last = 1.0 - np.sum(xperturb)
        xperturb = np.append(xperturb,x_last)
       
        cur_res = []

        # Lightweight constraints
        x = xperturb
        # MW
        mw = np.sum(mw_valid_vec*x)
        cur_res.append(mw)

        # HC
        n_h = np.sum(h_valid_vec*x)
        n_c = np.sum(c_valid_vec*x)
        hc = n_h/n_c
        cur_res.append(hc)

        # Sum constraint
        cur_res.append(np.sum(x))

        cur_res = np.transpose(np.atleast_2d(cur_res))

        perturbquantities = np.hstack((perturbquantities,cur_res))

    # Collecting results
    stack = np.empty((0,noptim),float)
    for kcase in range(num_lwcon):
        grad = (np.abs(perturbquantities[kcase, :]-quantityrefs[kcase])
                - np.abs(initquantities[kcase]-quantityrefs[kcase]))/step
        stack = np.vstack((stack, grad))

    return stack
    # sys.exit(0)

    # # Building tasks
    # tasks = []
    # taskindex = -1
    # for case in cases:
    #     taskindex += 1
    #     tasks.append((case, species_index_valid, x, taskindex))

    # Perturbed quantities
    step = 1e-2
    for kcase, case in enumerate(cases):
        for kx in range(noptim):
            xperturb = np.array(x)
            xperturb[kx] += step
            taskindex += 1
            tasks.append((case, species_index_valid, xperturb, taskindex))

    # Collecting gradient results
    logging.debug('Collecting Grad G results')

    results = np.array(Par.tasklaunch(tasks))
    initquantities = results[0:len(cases)]
    perturbquantities = results[len(cases):].reshape((len(cases), noptim))

    for kcase, case in enumerate(cases):
        grad = (np.abs(perturbquantities[kcase, :]-quantityrefs[kcase])
                - np.abs(initquantities[kcase]-quantityrefs[kcase]))/step
        stack = np.hstack((stack, grad))

    # Logging info
    logging.info('Gradient G done')
    logging.info('Current x')
    logging.info(str(x).replace('  ', ','))
    logging.info('Current sum(x)')
    logging.info(str(np.sum(x)))
    logging.info('Current distance vector')
    logging.info(str(np.abs(initquantities-quantityrefs)/quantityrefs))

    candidatecheck = True
    errortab = np.abs(initquantities-quantityrefs)
    for kerror, error in enumerate(errortab):
        # logging.info( 'Error %s %s %s'%( kerror, error, tolerances[kerror] ) )
        if error > tolerances[kerror]*1.1:
            candidatecheck = False
    if candidatecheck:
        logging.info('Current solution satisfies constraint')
        logging.info('Current vector')
        logging.info(str(x).replace('  ', ','))
        # TODO : Modify debugging for eliminated species
        # logging.info('Number of eliminated species')
        # count = np.where(x < Global.params.threshold)
        # logging.info(str(len(count[0])))
        # logging.info('List of eliminated species')
        # for kdamp in range(len(species_damp)):
        #     if x[kdamp] < Global.params.threshold:
        #         logging.info(
        #             Global.gas.species_names[species_index_damp[kdamp]] + ' ' + str(x[kdamp]))
    else:
        logging.info('Current solution does not satisfy constraint')
        logging.info('Current vector')
        logging.info(str(x).replace('  ', ','))

    logging.info('End Gradient PP')

    return stack
