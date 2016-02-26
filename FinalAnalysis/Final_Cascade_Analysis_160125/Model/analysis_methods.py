from .model import *
import os

def evaluateOnlySS(a, k, b, u, n):
    '''
    A method that only evaluates the steady state points for a
    repression cascade system.

    :param a:
    :param k:
    :param b:
    :param u:
    :param n:
    :return:
    '''
    assert len(a) == len(k)
    ss = []
    ss_0 = []
    dr_array = []
    ss = cascade_ss(len(a), a, k, b, u, n)
    ss_0 = cascade_ss(len(a), a, k, b, 0.0, n)
    for ss_value, ss_0_value in zip(ss, ss_0):
        values = [ss_value, ss_0_value]
        dr = np.log10(max(values)/min(values))
        dr_array.append(dr)
    p = pd.DataFrame( {"ss": ss,
                      "ss_0": ss_0,
                      "dr": dr_array})
    p = p.unstack()
    p = pd.DataFrame(p).transpose()
    return p

def evaluateSystem(a, k, b, u, n, kmodel=False, dosemodel=False):
    '''
    Evaluates a repression cascade system by simulating a set of ODE
    as nested hill-functions.
    :param a: promoter strength parameters
    :param k: repression strength parameters
    :param b: degradation/dilution parameters
    :param u: input strength parameters
    :param n: hill-coefficient parameter
    :param kmodel: whether to evaluate the kinetics model trajectories. The kinetics model
    makes a couple of assumptions about how the parameters for each layer are related to
    each other. This is important in our kinetic induction experiments.
    :param dosemodel: whether to evaluate the steady state properties. Will evaluate the
    maximum slope and location of the inflection point in the dose response curve.
    :return:
    '''
    global t_final, dt
    assert len(a) == len(k)
    dr_array = []
    ttss_array = []
    ss = []
    ss_0 = []
    dose_max_rate = []
    dose_inflection = []
    max_rate = []
    traj = None
    if kmodel:
        traj = kineticsModel(a, k, b, u, n)
        traj.columns = [1, 2, 3, 4, 5]
        ss_0 = traj.iloc[0]
        ss = traj.iloc[-1]
    else:
        ss = cascade_ss(len(a), a, k, b, u, n)
        ss_0 = cascade_ss(len(a), a, k, b, 0.0, n)
        traj = integrateSystem(len(a), a, k, b, u, n)
        traj = pd.DataFrame(traj)
        traj = traj.set_index(0) #only for integrate system, comment if using kinetics model
        traj.reindex()
        traj = pd.DataFrame(traj)
        i_array = 0.006*2**np.linspace(0,20,100)
    for ss_value, ss_0_value in zip(ss, ss_0):
        values = [ss_value, ss_0_value]
        dr = np.log10(max(values)/min(values))
        dr_array.append(dr)
    ttss_array = getTimeToHalfMaximum(traj, ss_0, ss, 0.001)
    max_rate = traj.apply(np.gradient).apply(np.abs).apply(np.max)/dt
    p = pd.DataFrame( {"ss": ss,
            "ss_0": ss_0,
            "dr": dr_array,
            "ttss": ttss_array,
            "max_rate": list(max_rate)} )
    if dosemodel:
        dose = ssModel(i_array, a, k, b, u, n)
        dose.drop('rel exp', axis=1, inplace=True)
        dose_deriv = dose.apply(np.gradient).apply(lambda x: x/np.gradient(dose.index))
        dose_max_rate = dose_deriv.apply(np.abs).apply(np.max)
        dose_inflection = dose_deriv.apply(np.abs).apply(np.argmax)
        p["dose_max_rate"] = list(dose_max_rate)
        p["dose_inflection"] = list(dose_inflection)
    p = p.unstack()
    p = pd.DataFrame(p).transpose()
    return p

message = ""
describe = ""

def modelbootstrap(input_filename, rounds=10, outfile=None, overwrite=False):
    '''
    This method boot straps the parameter values by resampling each fit.
    Note: this is obsolete and no longer used in the paper.
    :param input_filename:
    :param rounds:
    :param outfile:
    :param overwrite:
    :return:
    '''
    paramset = pd.read_csv(input_filename, header=[0,1], index_col=0)
    ran_params = []
    global message, describe, t_final, dt
    t_final = 100.0
    dt = 0.1
    final_results = pd.DataFrame()
    for j in range(rounds):
        results = pd.DataFrame()
        for i in range(len(paramset)):
            p = paramset.iloc[np.random.randint(0,len(paramset))]
            a = np.array(p['a'])
            k = np.array(p['k'])
            b = [float(p['b'])]
            u = [float(p['u'][0])]
            n = [float(p['n'])]
            params = np.concatenate([a, k, b, u, n])
            index = []
            for i, key in zip([4, 4, 1, 1, 1], list("akbun")):
                index = index + zip([key]*i, range(i))
            index = pd.MultiIndex.from_tuples(index)
            params = pd.DataFrame(params, index=index).transpose() #, columns = pd.MultiIndex.from_tuples(index))
            b = b[0]
            u = u[0]
            n = n[0]
            metrics = evaluateSystem(a, k, b, u, n, kmodel=True)
            joined = params.join(metrics)
            results = pd.concat([results, joined], ignore_index=True)
        results = pd.DataFrame(results.apply(np.mean)).transpose()
        final_results = pd.concat([final_results, results])
        final_results = final_results.reindex()
        message = final_results
    if not outfile == None:
        if not os.path.isfile(outfile) or overwrite:
            final_results.to_csv(outfile)
        else:
            savedresults = pd.read_csv(outfile, header=[0,1], index_col=0)
            savedresults.columns = final_results.columns
            message = "%d params read from file" % len(savedresults)
            message = message + '\n' + "Appending to file %s" % os.path.basename(outfile)
            newresults = pd.concat([savedresults, final_results], ignore_index=True)
            newresults.to_csv(outfile)
            describe = newresults.describe()
    return final_results

def resample_parameters_and_evaluate(numlayers,
                                     input_filename = None,
                                     resample = True,
                                     parameter_bounds = None,
                                     num_rounds=500,
                                     outfile=None,
                                     overwrite=False,
                                     n_bounds=None,
                                     n_force=None,
                                     dosemodel=True,
                                     cleanup=None):
    '''
    This method projects the repression cascade model and randomly chooses a parameter from the parameter sets
    located in tne 'input_filename' file. It then simulates the dynamics of the system and saves the output
    to 'outfile.'

    :param numlayers: number of repression cascade layers to evaluate
    :param input_filename: input_filename (csv) from which to resample the parameter values
    :param resample: whether to resample from the parameter sets
    :param parameter bounds: bounds from which to resample parameters from a uniform distribution. Overrides
    resample
    :param num_rounds: number of rounds to resample parameters and evaluate system
    :param outfile: file location to save the output data
    :param overwrite: overwrites the current file locates at 'outfile'
    :param n_bounds: bounds which to examine the hill-coefficient. Overrides parameter resampling.
    :param n_force: forces the system to evaluate at hill-coefficient n_force. Overrides n_bounds.
    :param dosemodel:
    :return:
    '''
    global t_final, dt
    dt = 0.8
    t_final = 750.0
    ran_params = []
    results = pd.DataFrame()
    global message
    for i in range(num_rounds):
        print "Round {}".format(i)
        a, k, b, u, n = 0, 0, 0, 0, 0
        if resample:
            paramset = pd.read_csv(input_filename, header=[0,1], index_col=0)
            a = np.random.choice(paramset.a.unstack(), size=numlayers)
            k = np.random.choice(paramset.k.unstack(), size=numlayers)
            b = np.random.choice(paramset.b.unstack(), size=1)
            u = np.random.choice(paramset.u.unstack(), size=1)
            n = np.random.choice(paramset.n.unstack(), size=1)
        if parameter_bounds is not None:
            a = np.random.uniform(parameter_bounds['a'][0], parameter_bounds['a'][1], size=numlayers)
            k = np.random.uniform(parameter_bounds['k'][0], parameter_bounds['k'][1], size=numlayers)
            b = np.random.uniform(parameter_bounds['b'][0], parameter_bounds['b'][1], size=1)
            u = np.random.uniform(parameter_bounds['u'][0], parameter_bounds['u'][1], size=1)
            n = np.random.uniform(parameter_bounds['n'][0], parameter_bounds['n'][1], size=1)
        if not n_bounds == None:
            n = np.random.uniform(*n_bounds)
            n = np.array([n])
        if not n_force == None:
            n = [n_force]
        params = np.concatenate([a, k, b, u, n])
        index = []
        for i, key in zip([numlayers, numlayers, 1, 1, 1], list("akbun")):
            index = index + zip([key]*i, range(i))
        index = pd.MultiIndex.from_tuples(index)
        params = pd.DataFrame(params, index=index).transpose() #, columns = pd.MultiIndex.from_tuples(index))
        b = float(b)
        u = float(u)
        try:
            n = float(n)
        except:
            n = float(n[0])
        metrics = evaluateSystem(a, k, b, u, n, dosemodel=dosemodel)
        #metrics = evaluateOnlySS(a, k, b, u, n)
        joined = params.join(metrics)
        results = pd.concat([results, joined], ignore_index=True)
    if cleanup:
        results = cleanup_data(results, numlayers)
    if not outfile == None:
        if not os.path.isfile(outfile) or overwrite:
            print "Saving to {}".format(outfile)
            results.to_csv(outfile)
        else:
            savedresults = pd.read_csv(outfile, header=[0,1], index_col=0)
            savedresults.columns = results.columns
            message = "%d params read from file" % len(savedresults)
            message = message + '\n' + "Appending to file %s" % os.path.basename(outfile)
            newresults = pd.concat([savedresults, results], ignore_index=True)
            newresults.to_csv(outfile)
            print "Appending to {}".format(outfile)
            results = newresults
    return results


###################### Cleanup methods ######################
def add_level1_columns(df, label, df_level1):
    for col in df_level1.columns:
        df.loc[:,(label,col)] = df_level1.loc[:,col]
    return df

def drop_multiindex_column(df, level1=None, level2=None, inplace=True):
    c = df.loc[:,(level1)]
    if not inplace:
        df = df.copy()
    df.drop(level1, inplace=True, axis=1, level=0)
    c.drop(level2, axis=1, inplace=True)
    add_level1_columns(df, level1, c)
    return df

def reverse_multiindex_columns(df):
    columns = set(df.columns.get_level_values(0))
    for col in columns:
        c = df.loc[:,col]
        cols_level2 = c.columns
        c.columns = cols_level2[::-1]
        c.sort_index(axis=1, inplace=True)
        df[col] = c
    df.sort(inplace=True, axis=1)
    return df

def cleanup_data(df, numlayers):
    print "Dropping ", numlayers
    df.drop(str(numlayers+1), axis=1, level=1, inplace=True)
    df.drop(numlayers+1, axis=1, level=1, inplace=True)
    df = reverse_multiindex_columns(df)
    return df
