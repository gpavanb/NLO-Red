import numpy as np
import Global
import logging
import sys
import Parallelism as Par


def getQuantityCase(case, x):
    #case.storeMultiplier(species_index_valid, x)
    return case.getQuantity(x)

def eval_f(xf, user_data=None):
    logging.debug('F evaluation')
    noptim = Global.params.noptim
    x = np.array(xf)/sum(np.array(xf))
    costfunc = 2.0 - np.sum(x*x) # 1.0 - np.sum(x*x) - np.sum(x)*np.sum(x) + 2*np.sum(x)
    #costfunc = 2.0 * np.sum(x*(1-x)) + np.sum(x*x)
    return costfunc

def eval_grad_f(xf, user_data=None):
    logging.debug('Grad F evaluation')
    noptim = Global.params.noptim
    num_lwcon = Global.params.num_lwcon
 
    x = xf/np.sum(xf)
    cost = 2.0 - np.sum(x*x)

    # Evaluate gradient using perturbation
    perturbquantities = np.empty(noptim,float)
    step = Global.step
    results = []
    for kx in range(noptim):
        xp = np.copy(xf)
        xp[kx] += step
        xp = xp/np.sum(xp)
       
        perturbquantities[kx] = 2.0 - np.sum(xp*xp)

    grad_f = (perturbquantities-cost)/step 

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
 
    # Lightweight constraints
    x = np.copy(xf)
    x = x/np.sum(x)
    # MW
    mw = np.sum(mw_valid_vec*x)
    quantities.append(mw)

    # HC
    n_h = np.sum(h_valid_vec*x)
    n_c = np.sum(c_valid_vec*x)
    hc = n_h/n_c
    quantities.append(hc)

    # Asynchronous constraints
    cases = Global.params.cases
    quantityrefs = Global.params.quantityrefs

    # Building Tasks
    tasks = []
    taskindex = -1
    for case in cases:
        taskindex += 1
        tasks.append((case, x, taskindex))

    # Retrieve array of results
    quantities_new = np.array(Par.tasklaunch(tasks))
    quantities.extend(quantities_new)
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
    x_orig = np.copy(xf)
    x_orig = x_orig/np.sum(x_orig)
    x = np.copy(x_orig)
    # MW
    mw = np.sum(mw_valid_vec*x)
    cur_res.append(mw)

    # HC
    n_h = np.sum(h_valid_vec*x)
    n_c = np.sum(c_valid_vec*x)
    hc = n_h/n_c
    cur_res.append(hc)

    initquantities = cur_res
 
    # Compute perturbations
    perturbquantities = np.empty((num_lwcon,0),float)
    step = Global.step
    results = []
    for kx in range(noptim):
        xperturb = np.copy(xf)
        xperturb[kx] += step
        xperturb = xperturb/np.sum(xperturb)
       
        cur_res = []

        # Lightweight constraints
        x = np.copy(xperturb)
        # MW
        mw = np.sum(mw_valid_vec*x)
        cur_res.append(mw)

        # HC
        n_h = np.sum(h_valid_vec*x)
        n_c = np.sum(c_valid_vec*x)
        hc = n_h/n_c
        cur_res.append(hc)

        # Append data
        cur_res = np.transpose(np.atleast_2d(cur_res))
        perturbquantities = np.hstack((perturbquantities,cur_res))

    # Collecting results
    stack = np.empty((0,noptim),float)
    for kcase in range(num_lwcon):
        grad = (np.abs(perturbquantities[kcase, :]-quantityrefs[kcase])
                - np.abs(initquantities[kcase]-quantityrefs[kcase]))/step
        stack = np.vstack((stack, grad))

    # Building tasks
    tasks = []
    taskindex = -1
    for case in cases:
        taskindex += 1
        tasks.append((case, x_orig, taskindex))

    # Perturbed quantities
    step = Global.step
    for kcase, case in enumerate(cases):
        for kx in range(noptim):
            xperturb = np.copy(xf)
            xperturb[kx] += step
            xperturb = xperturb/np.sum(xperturb)
            taskindex += 1
            tasks.append((case, xperturb, taskindex))

    # Collecting gradient results
    logging.debug('Collecting Grad G results')

    results = np.array(Par.tasklaunch(tasks))
    initquantities_new = results[0:len(cases)]
    initquantities.extend(initquantities_new)
    perturb_new = results[len(cases):].reshape((len(cases), noptim))
    if (len(cases) != 0):
        perturbquantities = np.vstack((perturbquantities,perturb_new))

    for kcase in xrange(len(cases)):
        idx = kcase + num_lwcon # Only nonlinear constraints computed here
        grad = (np.abs(perturbquantities[idx, :]-quantityrefs[idx])
                - np.abs(initquantities[idx]-quantityrefs[idx]))/step
        stack = np.vstack((stack, grad))

    # Logging info
    logging.info('Gradient G done')
    logging.info('Current x')
    logging.info(str(x).replace('  ', ','))
    logging.info('Current sum(x)')
    logging.info(str(np.sum(x)))
    logging.info('Current distance vector')
    logging.info(str(np.abs(initquantities-quantityrefs)))

    candidatecheck = True
    errortab = np.abs(initquantities-quantityrefs)
    for kerror, error in enumerate(errortab):
        # logging.info( 'Error %s %s %s'%( kerror, error, tolerances[kerror] ) )
        # Skip check for sum constraint
        if (tolerances[kerror] == 0):
            continue
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

    logging.info('End Gradient constraints')

    return stack
