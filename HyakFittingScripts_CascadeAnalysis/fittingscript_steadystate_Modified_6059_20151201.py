# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:50:48 2015

@author: Justin Vrana
"""
###########################################################################
########################### Import Modules ##############################
###########################################################################

import numpy as np
from scipy import integrate
from scipy import optimize
import pandas as pd
#import pylab as plt #optional
import time
from datetime import datetime
import os
###########################################################################
########################### Input Parameters ##############################
###########################################################################

t_final = 100.0
dt = 1.0

root = "/Users/klavinslab"
save_location = os.path.join(root, "Documents/CRISPRCascadeAnalysis/")
output_location = os.path.join(save_location, "fittingoutput_6059.csv")
exp_data = os.path.join(save_location, "Data/ExperimentalData")
print(exp_data)
logfile_location = os.path.join(save_location, "test.csv")
dose_response_location = os.path.join(exp_data, "DoseResponse.csv")
kinetics_layer1_location = os.path.join(exp_data, "KineticsLayer1.csv")
kinetics_layer2_location = os.path.join(exp_data, "KineticsLayer2.csv")
kinetics_layer3_location = os.path.join(exp_data, "KineticsLayer3.csv")
kinetics_layer4_location = os.path.join(exp_data, "KineticsLayer4.csv")

fitKinetics = False # Use kinetics model to fit the data
fitDoseResponse = True # Use the dose response model to fit the data
best_params_df = None

au_scale = 1600. # scaling factor for the model/exp data for fitting
t_scale = 25.0 # scaling factor (time) for the model/exp data for fitting
lowest_score = float('inf') # lowest_score fit
best_params = [] # list of all lowest scores and parameters found during fit



###########################################################################
######################### Modeling Functions ##############################
###########################################################################

def importDoseResponseData(location):
    print("Importing Dose Response Data:")
    exp_dose = pd.read_csv(location)
    AU = exp_dose['FL1.Amean']
    treatment = exp_dose['treatment']
    exp_dose = pd.pivot_table(exp_dose, index='treatment', columns=["strain"], values="FL1.Amean")
    return exp_dose

def importKineticsData(layer1loc, layer2loc, layer3loc, layer4loc):
    print("Importing Kinetics Data")
    filelocations = [layer1loc, layer2loc, layer3loc, layer4loc]
    alllayers = [pd.read_csv(x) for x in filelocations]
    header = ['time', 'au']
    ttrim = 3000
    for i, layer in enumerate(alllayers):
        layer.columns = header
        layer['layer'] = i + 1
        layer = layer[layer.time < ttrim]
    layers = pd.concat(alllayers)
    exp_kinetics = pd.pivot_table(layers, values='au', columns=['layer', 'time'])
    return exp_kinetics

# Define kinetics for cascade model
def cascade_kinetics(t, r, a, k, b, u, n):
    layer_depth = len(r) - 1
    dR = [0] * len(r)
    dR[layer_depth] = u[layer_depth-1] - b * r[layer_depth]
    for l in range(layer_depth)[::-1]:
        dR[l] = a[l] / (1+k[l] * r[l+1] ** n) - b * r[l]
    return dR

# Define steady state for cascade model
def cascade_ss(layer_depth, a, k, b, u, n):
    assert len(a) == len(k) == len(u)
    assert len(a) >= layer_depth
    u = u[layer_depth-1]
    R = [0]*(layer_depth+1)
    R[layer_depth] = u/b
    for l in range(layer_depth)[::-1]:
        R[l] = a[l] / ((1+k[l] * R[l+1] ** n)*b)
    return R

# Beta-estradiol transfer function
def tfxn(inducer_at_max, inducer):
    AB = 542.8745566
    K = 0.269635359
    n = 1.238933919
    def fxn(u):
        return (AB*K*u**n)/(1+K*u**n)
    rel_induction = fxn(inducer)/fxn(inducer_at_max)
    return rel_induction

# Integrating cascade model
def integrateSystem(numlayers, a, k, b, u, n):
    ss = cascade_ss(numlayers, a, k, b, np.zeros(len(a)), n)
    r1 = integrate.ode(cascade_kinetics)
    r1.set_initial_value(ss).set_f_params(a, k, b, u, n)
    data1 = []
    while r1.successful() and r1.t < t_final:
        r1.integrate(r1.t+dt)
        data1.append( [r1.t] + list(r1.y) )
    return data1

# Define model for dose response experiment
# (including Beta-estradiol transfer function)
def ssModel(i_array, a, k, b, u, n):
    assert len(a) == len(k) == len(u)
    num_layers = len(a)+1
    dose = pd.DataFrame(index = i_array, columns = list(range(1,num_layers+1)) + ['rel exp'])
    for i in i_array:
        rel_exp = tfxn(100.0, i)
        rel_u = rel_exp * u
        dose.loc[i] = [cascade_ss(l, a, k, b, rel_u, n)[0] for l in range(num_layers)] + [rel_exp]
    return dose

# Define model for kinetics experiment
def kineticsModel(a, k, b, u, n):
    assert len(a) == len(k) == len(u)
    num_layers = len(a)+1
    layers = {}
    time = []
    for l in range(num_layers):
        system = integrateSystem(l, a, k, b, u, n)
        unzipped = zip(*system)
        layers[l] = unzipped[1]
        time = unzipped[0]
    return pd.DataFrame(layers, index=time)


###########################################################################
################### Scoring and Fitting Functions #########################
###########################################################################

# Scoring function for four layer cascade
def fourlayer_scoring(prm):
    a = prm[:4]
    k = prm[4:8]
    b = prm[8:10]
    u = prm[10]
    n = prm[11]
    return scoring(a, k, b, np.array([u]*len(a)), n)

# Scoring expeirmental vs. model for dose response data
def scoreDose(a, k, b, u, n, num_params):
    params = [a, k, b, u, n]
    N_ss = np.sum(exp_dose.count())
    v = N_ss - num_params - 1
    sumsquares_ss = 0
    ss_model = ssModel(exp_dose.index, *params) * au_scale
    ss_model = ss_model.iloc[:,1:-1]
    ss_diff = np.subtract(exp_dose, ss_model)
    ss_diff = np.divide(ss_diff, 1000.0)
    sumsquares_ss = np.sum(np.square(ss_diff[6059]))
    chisquared_ss = np.sum(sumsquares_ss/v)
    return chisquared_ss, ss_model #+ chisquared_k

# Scoring experimental vs. model for kinetics data
def scoreKinetics(a, k, b, u, n, num_params):
    params = [a, k, b, u, n]
    N_k = exp_kinetics.count()
    v = N_k - num_params - 1
    sumsquares_k = 0 #sum of squares

    kinetics_model = kineticsModel(*params) * au_scale
    kinetics_model = kinetics_model.set_index(kinetics_model.index * t_scale)
    layers = [exp_kinetics[i] for i in range(1,5)]
    layers_diff = []
    for i, layer in enumerate(layers):
        layer_model = []
        for t in layer.index:
            layer_model.append(getValue(t, kinetics_model)[1:])
        layer_model = pd.DataFrame(layer_model)[i+1]
        layer_diff = np.subtract(layer_model, layer)
        layer_diff = np.divide(layer_diff, 1000.0)
        sumsquares_k += np.sum(np.square(layer_diff))
    chisquared_k = sumsquares_k/v
    return chisquared_k, kinetics_model

# Helper function
def getValue(t, df):
    i = np.abs(df.index - t).argmin()
    return df.iloc[i]

# Scoring function
def scoring(a, k, b, u, n):
    assert len(a) == len(k) == len(u)
    global lowest_score, best_params, best_params_df, fitKinetics, fitDoseResponse

    num_params = len(a) + len(k) + 3 + len(u) #num parameters

    b1, b2 = b
    params = [a, k, b1, b2, u, n]
    chitotal = 0
    ss_model = None
    k_model = None
    if fitDoseResponse:
        chiD, ss_model = scoreDose(a, k, b2, u, n, num_params)
        chitotal += chiD
    if fitKinetics:
        chiK, k_model = scoreKinetics(a, k, b1, u, n, num_params)
        chitotal += chiK

    paramkey = ['score', 'a', 'k', 'b', 'u', 'n']
    param_iter = {'score': 1, 'a':4, 'k':4, 'b':2, 'u':4, 'n':1}
    index = []
    for key in paramkey:
        index = index + list(zip([key]*param_iter[key], range(param_iter[key])))
    index = pd.MultiIndex.from_tuples(index)


    if chitotal < lowest_score:
        print("Lowest Score Found:", chitotal)
        best_params.append([chitotal] + params)
        best_params.append([chitotal] + params)
        lowest_score = chitotal
        paramlist = []
        for thisparam in best_params:
            score, a, k, b1, b2, u, n = thisparam
            paramlist.append([score] + list(a) + list(k) + [b1] + [b2] + list(u) + [n])
        best_params_df = pd.DataFrame(paramlist, columns=index)
        best_params_df.to_csv(output_location)
    return chitotal

def timestamp():
    ts = int(time.time()*1000000)
    return datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H-%M-%S') + "_" + str(ts)

def makelog():
    global save_location, dose_response_location, logfile_location, kinetics_layer1_location
    global kinetics_layer2_location, kinetics_layer3_location, kinetics_layer4_location
    global a_bounds, k_bounds, u_bounds, n_bounds, b_bounds, au_scale, t_scale

    tstmp = timestamp()

    lg = "Timestamp:"+ str(tstmp)
    lg = lg + "\n" + "Differential Evoluation Fit Settings"
    lg = lg + "\n" + "\tImport Settings"
    lg = lg + "\n" + "\t\tUse Kinetics Data:"+ str(fitKinetics)
    lg = lg + "\n" + "\t\tUse Dose Response Data:"+ str(fitDoseResponse)
    lg = lg + "\n" + "\tImport Files:"
    if fitDoseResponse:
        lg = lg + "\n" + "\t\t\tDose Data:"+ str(dose_response_location)
    if fitKinetics:
        lg = lg + "\n" + "\t\t\tKinetics Data (1):"+ str(kinetics_layer1_location)
        lg = lg + "\n" + "\t\t\tKinetics Data (2):"+ str(kinetics_layer2_location)
        lg = lg + "\n" + "\t\t\tKinetics Data (3):"+ str(kinetics_layer3_location)
        lg = lg + "\n" + "\t\t\tKinetics Data (4):"+ str(kinetics_layer4_location)
    lg = lg + "\n" + "\t\tSave Setings"
    lg = lg + "\n" + "\t\t\tSave Location:"+ str(save_location)
    lg = lg + "\n" + "\tIntegration Settings"
    lg = lg + "\n" + "\t\ttime_final:"+ str(t_final)
    lg = lg + "\n" + "\t\tdt:"+ str(dt)
    lg = lg + "\n" + "\tScaling Factors"
    lg = lg + "\n" + "\t\tau scale:"+ str(au_scale)
    lg = lg + "\n" + "\t\ttime scale:"+ str(t_scale)
    lg = lg + "\n" + "\tBounds"
    lg = lg + "\n" + "\t\ta_bounds:"+ str(a_bounds)
    lg = lg + "\n" + "\t\tk_bounds:"+ str(k_bounds)
    lg = lg + "\n" + "\t\tb_bounds:"+ str(b_bounds)
    lg = lg + "\n" + "\t\tu_bounds:"+ str(u_bounds)
    lg = lg + "\n" + "\t\tn_bounds:"+ str(n_bounds)
    log = open(logfile_location, 'w')
    log.write(lg)
    log.close()
    return lg

###########################################################################
######################### Execute Fitting Procedure #########################
###########################################################################

# Define bounds
a_bounds = [(0.1, 1.5)] * 4
k_bounds = [(0.1, 1.5)] * 4
b_bounds = [(0.03, 0.13)] * 2
u_bounds = [(0.05, 1.0)] * 1
n_bounds = [(0.5, 3.5)] * 1

# Import expeirmental data
exp_kinetics = None
exp_dose = None
if fitKinetics:
    exp_kinetics = importKineticsData(kinetics_layer1_location, kinetics_layer2_location, kinetics_layer3_location, kinetics_layer4_location)
if fitDoseResponse:
    exp_dose = importDoseResponseData(dose_response_location)

# Execute diff. evolution procedure
bounds = a_bounds + k_bounds + b_bounds + u_bounds + n_bounds
makelog()
print("*" * 30)
print("Running Differential Evolution...")
print("*" * 30)
res = optimize.differential_evolution(fourlayer_scoring, bounds, disp=True, popsize=100, mutation=(0.7, 1.0), recombination=0.7)
