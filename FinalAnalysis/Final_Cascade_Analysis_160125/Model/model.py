import numpy as np
from scipy import integrate
import pandas as pd

global t_final
global dt
t_final = 500.
dt = 0.01

def integrateSystem(numlayers, a, k, b, u, n):
    '''
    Calculates trajectories for a repression cascade model.
    :param numlayers:
    :param a:
    :param k:
    :param b:
    :param u:
    :param n:
    :return:
    '''
    ss = cascade_ss(numlayers, a, k, b, 0.0, n)
    r1 = integrate.ode(cascade_kinetics)
    r1.set_initial_value(ss).set_f_params(a, k, b, u, n)
    data1 = []
    while r1.successful() and r1.t < t_final:
        r1.integrate(r1.t+dt)
        data1.append( [r1.t] + [float(i) for i in list(r1.y)] )
    return data1

def ssModel(i_array, a, k, b, u, n):
    '''
    Calculates the steady state dose response curve for a set of
    beta-estradiol concentrations.
    :param i_array: beta-estradiol concentrations in uM (np.array)
    :param a:
    :param k:
    :param b:
    :param u:
    :param n:
    :return:
    '''
    assert len(a) == len(k)
    num_layers = len(a)+1
    dose = pd.DataFrame(index = i_array, columns = list(range(1,num_layers+1)) + ['rel exp'])
    for i in i_array:
        rel_exp = tfxn(100.0, i)
        rel_u = rel_exp * u
        dose.loc[i] = [cascade_ss(l, a, k, b, rel_u, n)[0] for l in range(num_layers)] + [rel_exp]
    return dose

def kineticsModel(a, k, b, u, n):
    '''
    The kinetics model makes a couple of assumptions about how the parameters
    for each layer are related to each other. This is important in our kinetic
    induction experiments. This method is mainly used for modeling our kinetics
    experiments and for fitting.
    :param a:
    :param k:
    :param b:
    :param u:
    :param n:
    :return:
    '''
    assert len(a) == len(k)
    num_layers = len(a)+1
    layers = {}
    time = []
    for l in range(num_layers):
        system = integrateSystem(l, a, k, b, u, n)
        unzipped = zip(*system)
        layers[l] = unzipped[1]
        time = unzipped[0]
    return pd.DataFrame(layers, index=time)

def getTimeToHalfMaximum(trajectory, init, ss, ds):
    half = ss - init
    half = half/2.0
    half = ss - half
    ttss = []
    g = np.abs(trajectory - half)
    ttss = g.apply(np.argmin)
    return ttss

def getValue(t, df):
    i = np.abs(df.index - t).argmin()
    return df.iloc[i]

def cascadesystem_steadystate_simple(numlayers, a, k, b, u, n):

    Y = []
    Y.append(u/b)
    for i in range(numlayers):
        Y.append((1.*a/(1+k*Y[i]**n))/b)
    return np.array(Y)

def dynamic_range_plot(numlayers, a, k, b, u, n):
    params = [a,k,b,u,n]
    ss = cascadesystem_steadystate_simple(numlayers, *params)
    params = [a, k, b, 0.0, n]
    ss_0 = cascadesystem_steadystate_simple(numlayers, *params)
    #plt.yscale('log')
    #plt.plot(np.log10(np.abs(ss-ss_0)))
    dr = []
    for ss_value, ss_0_value in zip(ss, ss_0):
        values = [ss_value, ss_0_value]
        dr.append(max(values)/min(values))
    #return np.abs(ss-ss_0)
    metric = 20.0*np.log10(dr)
    #metric = np.abs(ss-ss_0)
    plt.plot(metric)
    return metric

def cascade_kinetics(t, r, a, k, b, u, n):
    '''
    Equation used to integrate the repression cascade model.

    :param t:
    :param r:
    :param a:
    :param k:
    :param b:
    :param u:
    :param n:
    :return:
    '''
    layer_depth = len(r) - 1
    dR = [0] * len(r)
    dR[layer_depth] = u - b * r[layer_depth]
    for l in range(layer_depth)[::-1]:
        dR[l] = a[l] / (1+k[l] * r[l+1] ** n) - b * r[l]
    return dR

def cascade_ss(layer_depth, a, k, b, u, n):
    '''
    Calculates the steady state point for a repression cascade model.

    :param layer_depth:
    :param a: list of promoter strength values of length layer_depth
    :param k: list of repression strength values of length layer_depth
    :param b: dilution degradation term (float)
    :param u: input strength term (float)
    :param n: hill-coefficient term (float)
    :return:
    '''
    assert len(a) == len(k)
    assert len(a) >= layer_depth
    R = [0]*(layer_depth+1)
    R[layer_depth] = u/b
    for l in range(layer_depth)[::-1]:
        R[l] = a[l] / ((1+k[l] * R[l+1] ** n)*b)
    return np.array(R)

def tfxn(inducer_at_max, inducer):
    '''
    Transfer function between amount of beta-estradiol and gRNA determined
    experimentally.
    :param inducer_at_max: amount of beta-estradiol in uM that produces maximal activation (float)
    :param inducer: amount of beta-estradiol in uM (flaot)
    :return: relative induction amount between 0.0 and 1.0 (float)
    '''
    AB = 542.8745566
    K = 0.269635359
    n = 1.238933919
    def fxn(u):
        return (AB*K*u**n)/(1+K*u**n)
    rel_induction = fxn(inducer)/fxn(inducer_at_max)
    return rel_induction