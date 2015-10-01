# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 12:05:43 2015

@author: Justin
"""
import pylab
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import numpy as np
from scipy import integrate
from scipy import optimize
from matplotlib.pylab import *
import os
import csv
import math
import re
import pandas as pd
class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

k_param = []
a_param = []
u_param = -1
b_param = -1
n_param = 1.1 #1.4459

dt = 0.1
t_final = 3000.

def cascadesystem(t, y, numlayers):
    # Grab global parameter values
    global k_param, a_param, u_param, b_param, n_param
    k, a, u, b, n = k_param, a_param, u_param, b_param, n_param
    
    #Check parameter lengths
    if not numlayers == len(k) or not numlayers == len(a):
        raise MyError("Number of parameters does not match number of layers")
    #Apply differentials with single parameter "layers"
    dY = [] #array of differentials
    dY.append(u - y[0]*b)
    for i in range(numlayers):
        dY.append(1.*a[i]/(1+k[i]*y[i]**n) - y[i+1]*b) #add y_i+1
    return dY

def cascadesystem_steadystate_simple(numlayers, k, a, u, b, n):
    
    Y = []
    Y.append(u/b)
    for i in range(numlayers):
        Y.append((1.*a/(1+k*Y[i]**n))/b)
    return np.array(Y)


def cascadesystem_steadystate(numlayers, k, a, u, b, n):
    if not numlayers == len(k) or not numlayers == len(a):
        raise MyError("Number of parameters does not match number of layers")
    
    Y = []
    Y.append(u/b)
    for i in range(numlayers):
        Y.append((1.*a[i]/(1+k[i]*Y[i]**n))/b)
    return np.array(Y)



def integrateCascade(numlayers, init, k, a, u, b, n):
    if not numlayers == len(k) or not numlayers == len(a):
        raise MyError("Number of parameters does not match number of layers")
    if not numlayers == len(init)-1:
        raise MyError("Number of initial values does not match number of layers")
        
    global k_param, a_param, u_param, b_param, n_param
    k_param, a_param, u_param, b_param, n_param = k, a, u, b, n
    r1 = integrate.ode(cascadesystem)
    r1.set_initial_value(init).set_f_params(numlayers)
    data1 = []
    while r1.successful() and r1.t < t_final:
        r1.integrate(r1.t+dt)
        data1.append( [r1.t] + list(r1.y) )
    return np.array(zip(*data1))



def dynamicrangeplot(numlayers, n, savefolder, newfoldername):
    ##################################
    ## Examine parameter sensitivites
    ##################################
    layers = numlayers
    data = []
    for a_mod in np.linspace(0.1, 5.0, 150):
        for k_mod in np.linspace(0.1, 5.0, 150):
            #for n_mod in np.linspace(0.5, 2.0, 10):
            a = np.ones(layers) * a_mod
            k = np.ones(layers) * k_mod
            u = 0.4
            b = 0.112
            ss = cascadesystem_steadystate(layers, k, a, u, b, n)
            ss_0 = cascadesystem_steadystate(layers, k, a, 0, b, n)
            dyn_range = abs(ss[-1] - ss_0[-1])
            degradation = ss[-1] - ss[-3]
            degradation = ss[-4] - ss[-2] #alternative plot
            data.append([a_mod, k_mod, n, dyn_range, degradation])
    data = np.transpose(np.array(data))
    x = data[0]
    y = data[1]
    z = data[3]
    #z[-1] = 70
    gridsize=100
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 12}
    pylab.rc('font', **font)
    #fig = plt.figure(figsize(16,12), dpi=300)
    plt.subplot(111)
    plt.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
    cb = plt.colorbar()
    cb.set_label('Dynamic range between steady states')
    plt.xlabel("a")
    plt.ylabel("k")
    title = "Dynamic range plot at layer %s for promoter strength 'a' \n and repression strength 'k' at n==%s" % (str(layers), str(n))
    plt.title(title)
    savelocation = savefolder + "/" + newfoldername
    if not os.path.exists(savelocation):
        os.makedirs(savelocation)
    plt.savefig(savelocation + "/" + title + ".png")
    plt.show()


def signaldegplot(numlayers, n, savefolder, newfoldername, nullclines=False):
    ##################################
    ## Examine parameter sensitivites
    ##################################
    layers = numlayers
    data = []
    for a_mod in np.linspace(0.1, 5.0, 150):
        for k_mod in np.linspace(0.1, 5.0, 150):
            #for n_mod in np.linspace(0.5, 2.0, 10):
            a = np.ones(layers) * a_mod
            k = np.ones(layers) * k_mod
            u = 0.4
            b = 0.112
            ss = cascadesystem_steadystate(layers, k, a, u, b, n)
            ss_0 = cascadesystem_steadystate(layers, k, a, 0, b, n)
            dyn_range = abs(ss[-1] - ss_0[-1])
            degradation = ss[-1] - ss[-3]
            degradation = ss[-4] - ss[-2] #alternative plot
            if nullclines:
                if degradation > 0:
                    degradation = 1
                elif degradation < 0:
                    degradation = -1
                else:
                    degradation = 0
            data.append([a_mod, k_mod, n, dyn_range, degradation])
    data = np.transpose(np.array(data))
    x = data[0]
    y = data[1]
    z = data[4]
    gridsize=100
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 12}
    pylab.rc('font', **font)
    #fig = plt.figure(figsize(16,12), dpi=300)
    plt.subplot(111)
    plt.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
    cb = plt.colorbar()
    cb.set_label('Positive increase in dynamic \nrange at layer %s' % str(layers))
    plt.xlabel("a")
    plt.ylabel("k")
    title = "Sensitifity plot for increase in dynamic range\n at layer %s for promoter strength 'a' \n and repression strength 'k' at n==%s" % (str(layers), str(n))
    if nullclines:
        title = "Nullcline " + title
    plt.title(title)
    savelocation = savefolder + "/" + newfoldername
    if not os.path.exists(savelocation):
        os.makedirs(savelocation)
    plt.savefig(savelocation + "/" + title + ".png")
    plt.show()

def getTimeToSteadyState(time, values, ds):
    ss = zip(*values)[-1]
    ttss = []  
    for i in range(len(values)):
        t = time[abs(np.abs(np.array(values[i])-ss[i]) - ds).argmin()]
        ttss.append(t)
    return np.array(ttss)

def getTimeToHalfMaximum(time, values, ds):
    ss = np.array(zip(*values)[-1])
    init = np.array(zip(*values)[1])
    half = (np.abs(ss - init))/2.0
    print "Steady state", ss
    print "Initialize", init
    print "Half", half
    ttss = []
    for i in range(len(values)):
        t = time[abs(np.abs(np.array(values[i])-half[i]) - ds).argmin()]
        ttss.append(t)
    return np.array(ttss)


def writeCSV(filelocation, x, y):
    data = zip(x, y)
    with open(filelocation, 'w') as f:
        a = csv.writer(f, delimiter=',')
        for row in data:
            s = str(row[0])
            s = s + ','.join([str(i) for i in row[1]])
            f.write(s)
            f.write("\n")



def generateSamplingData(params, params_sd, layers, numsamples, filename, timetosteadystate=True):
    #Add Header to sampling file
    sampling_file = open(filename, 'w')
    sampling_file.write("Sampling Parameters,a,k,u,b,n\n")
    sampling_file.write("mean,")
    sampling_file.write(','.join([str(i) for i in params]))
    sampling_file.write('\n')
    sampling_file.write("mean,")
    sampling_file.write(','.join([str(i) for i in params_sd]))
    sampling_file.write('\nSamples,a,k,u,b,n,steady state,steady state (u=0),time to steady state\n')
    #sampling_file.close()
    
    for i in range(numsamples):
        if i % 10 == 0:
            print "Iteration %s" % (str(i))
        # Sample parameter values
        random_parameters = []
        for j in range(5):
            param_value = params[j]
            param_value_sd = params_sd[j]
            num_sds = 1
            param_min = param_value - num_sds * param_value_sd
            param_max = param_value + num_sds * param_value_sd
            ran = np.random.normal(param_value, param_value_sd, size=layers).clip(param_min, param_max)
            random_parameters.append(ran)
        a = random_parameters[0]
        k = random_parameters[1]
        u = random_parameters[2][0]
        b = random_parameters[3][0]
        n = random_parameters[4][0]
#        a = np.random.normal(params[0], params_sd[0], size=layers).clip(params[0] - 3* params_sd[0], params[0] + 3* params_sd[0])
#        k = np.random.normal(params[1], params_sd[1], size=layers).clip(params[1] - 3* params_sd[1], params[1] + 3* params_sd[1])
#        u = np.random.normal(params[2], params_sd[2], size=1)[0].clip(params[2] - 3* params_sd[2], params[2] + 3* params_sd[2])
#        b = np.random.normal(params[3], params_sd[3], size=1).clip(0, params[3] - 3* params_sd[3], params[3] + 3* params_sd[3])[0]
#        n = np.random.normal(params[4], params_sd[4], size=1).clip(params[4] - 3* params_sd[4], params[4] + 3* params_sd[4])[0]

        ss = cascadesystem_steadystate(layers, k, a, u, b, n)
        ss_0 = cascadesystem_steadystate(layers, k, a, 0, b, n)
        for s in ss:
            if math.isnan(s):
                print a, k, u, b, n
                break
        for s in ss_0:
            if math.isnan(s):
                print a, k, u, b, n
                break
        # Compute time to steady state
        Y = 0
        if timetosteadystate:
            d = integrateCascade(layers, ss_0, k, a, u, b, n)
            Y = getTimeToSteadyState(d[0], d[1:], 0.001)
            Y = getTimeToHalfMaximum(d[0], d[1:], 0.001)
            if 1 == 1:
                for l in d[1::1]:
                    pylab.plot(d[0], l)
                pylab.xlim(0,500)
                pylab.show()
                print a
                print k
                print u
                print b
                print n
                print Y
        parameters = np.array([a, k, u, b, n])
        parameter_string = []
        for p in parameters:
            p_str = str(p)
            try:
                p_str = '|'.join(str(j) for j in p)
            except:
                pass
            parameter_string.append(p_str)
        
        data = np.array([ss * 1600., ss_0 * 1600., Y * 25.])
        
        data_string = []
        for d in data:
            d_str = str(d)
            try:
                d_str = '|'.join(str(j) for j in d)
            except:
                pass
            data_string.append(d_str)
        
        sampling_data = [str(i)] + parameter_string + data_string
        sampling_file.write(','.join(sampling_data))
        sampling_file.write('\n')
    sampling_file.close()

def readSamplingFile(filename):
        # Read sampling file
    sampling_file = open(filename, 'rb')
    ttss_sampling = []
    ss_sampling = []
    ss_0_sampling = []
    params = []
    with sampling_file as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        beginReading = False
        for row in reader:
            if beginReading:
                try:
                    ttss = [float(pt) for pt in row[8].split('|')]
                    ttss_sampling.append(ttss)
                except:
                    pass
                try:
                    ss = [float(pt) for pt in row[6].split('|')]
                    ss_0 = [float(pt) for pt in row[7].split('|')]
                    ss_sampling.append(ss)
                    ss_0_sampling.append(ss_0)
                    params.append(row[1:6])
                except:
                    pass
            if row[0] == "Samples":
                beginReading = True
    sampling_file.close()
    ttss_sampling = np.array(ttss_sampling)
    ss_sampling = np.array(ss_sampling)
    ss_0_sampling = np.array(ss_0_sampling)
    return ss_sampling, ss_0_sampling, ttss_sampling, params

if __name__ == "__main__":
    
    ###################################
    ## Plot Trajectories
    ###################################
    
    ###################################
    ## Plot Time to Steady-State and Dynamic Range
    ###################################    
    layers = 50
    numsamples = 1
    params = [0.068482583, 0.839400642, 0.60500040687, 0.12188042516, 1.48342953330]
    params_sd = [0.146595632, 0.294241234, 0.06754616, 0.009419987, 0.092091958]
    params_sd = [0.000001]*len(params)
    #params_sd[3] = 0.2
#    params_sd = np.ones(len(params_sd))*0.00001
#    params = params - params_sd
#    params_sd = np.ones(len(params_sd)) * 0.00001
    params[2] = params[0]
    params_sd[2] = params_sd[0]
    print params
    print params_sd
    root_folder = "/Users/Justin/Google Drive/KlavinsLab/Projects/FSM/Analysis/9-12-15/"
    filename = root_folder + "largesample.csv"
    generateSamplingData(params, params_sd, layers, numsamples, filename)
    
    n_files = []
    for n in linspace(0.1, 2.0, 20):
        filename = root_folder + "sampling_n%s%s.csv" % (str(n).split('.')[0], str(n).split('.')[1])
        pars = params[:]
        pars_sd = params_sd[:]
        pars[-1] = n
        pars_sd[-1] = 0.00001
        generateSamplingData(pars, pars_sd, layers, numsamples, filename, timetosteadystate=False) 
        n_files.append(filename)
    b_files = []
    for b in linspace(.05, 0.2, 3):
        filename =  root_folder + "sampling_b%s%s.csv" % (str(b).split('.')[0], str(b).split('.')[1])
        pars = params[:]
        pars_sd = params_sd[:]
        pars[-2] = b
        pars_sd[-2] = 0.00001
        print pars
        generateSamplingData(pars, pars_sd, layers, numsamples, filename, timetosteadystate=False) 
        b_files.append(filename)
#
#    filename = "/Users/Justin/Google Drive/KlavinsLab/Projects/FSM/Analysis/sampling_b05.csv"
#    params = [1.219238834, 0.868357922, 0.579300889, 0.123494219, 1.43980423]
#    params[-2] = 0.05
#    params_sd = [0.290498903, 0.290551181, 0.059192661, 0.00001, 0.066092914]
#    generateSamplingData(params, params_sd, 30, 1000, filename)
#    
#    filename = "/Users/Justin/Google Drive/KlavinsLab/Projects/FSM/Analysis/sampling_b20.csv"
#    params = [1.219238834, 0.868357922, 0.579300889, 0.123494219, 1.43980423]
#    params[-2] = 0.2
#    params_sd = [0.290498903, 0.290551181, 0.059192661, 0.00001, 0.066092914]
#    generateSamplingData(params, params_sd, 30, 1000, filename)

    def plotDelay(files, savefile=""):
        fileout = None
        if not savefile == "":
            print "Writing to file %s" % savefile
            fileout = open(savefile, 'w')
        for filename in files:
            ss_sampling, ss_0_sampling, ttss_sampling, params = readSamplingFile(filename)
            for i, ttss in enumerate(ttss_sampling):
                pylab.scatter(range(layers+1), ttss)
                print params[i]
                
            pylab.show
#            ttss_mean = np.mean(ttss_sampling, axis=0)
#            ttss_sd = np.std(ttss_sampling, axis=0)
#            print filename
#            pylab.plot(range(layers+1), ttss_mean)
#            pylab.fill_between(range(layers+1), ttss_mean+ttss_sd, ttss_mean-ttss_sd, facecolor='yellow', alpha=0.5)
#            if not fileout == None:
#                data = np.array([range(layers+1), ttss_mean, ttss_sd])
#                for row in data.T:
#                    row = [str(s) for s in row]
#                    fileout.write(','.join(row))
#                    fileout.write("\n")
        fileout.close()
        pylab.legend()
        pylab.show()

    def plotDynamicRange(files, savefile=""):
        fileout = None
        if not savefile == "":
            print "Writing to file %s" % savefile
            fileout = open(savefile, 'w')
        X = []
        Y = []
        for filename in files:
            ss_sampling, ss_0_sampling, ttss_sampling, params = readSamplingFile(filename)
            #ttss_mean = np.mean(ttss_sampling, axis=0)
            #ttss_sd = np.std(ttss_sampling, axis=0)
            ss_mean = np.mean(ss_sampling, axis=0)
            ss_sd = np.std(ss_sampling, axis=0)
            ss_0_mean = np.mean(ss_0_sampling, axis=0)
            ss_0_sd = np.std(ss_0_sampling, axis=0)
            dynamic_range = np.abs(ss_sampling - ss_0_sampling)
            dynamic_range_mean = np.mean(dynamic_range, axis=0)
            dynamic_range_sd = np.std(dynamic_range, axis=0)
            print filename
            lab = re.search("[a-z](\d+).csv", filename).group(1)
            n = lab[0] + "." + lab[1]
            n = float(n)
            X.append(n)
            logslope = (np.log10(dynamic_range_mean[0])-np.log10(dynamic_range_mean[-1]))/layers
            Y.append(logslope)
            pylab.plot(range(layers+1), dynamic_range_mean, label = lab)
            yer = dynamic_range_sd #np.sqrt(dynamic_range_mean)
            print dynamic_range_sd
            print np.sqrt(dynamic_range_mean)
            print yer
            pylab.errorbar(range(layers+1), dynamic_range_mean, xerr=0.0, yerr=yer, label = lab)
            
            pylab.plot(range(layers+1), dynamic_range_mean, label = lab)
            if not fileout == None:
                data = np.array([[lab]*(layers+1), range(layers+1), dynamic_range_mean, dynamic_range_sd])
                for row in data.T:
                    row = [str(s) for s in row]
                    fileout.write(','.join(row))
                    fileout.write("\n")
        fileout.close()
        pylab.legend()
        pylab.show()
        pylab.plot(X, Y)
        pylab.show()
    ss_sampling, ss_0_sampling, ttss_sampling, params = readSamplingFile(root_folder + "largesample.csv")
    #plotDynamicRange(n_files, savefile = root_folder + "dynamic_range_out_n.csv")
    #plotDynamicRange(b_files, savefile = root_folder + "dynamic_range_out_b.csv")
    #plotDelay([root_folder + "largesample.csv"])#, savefile = root_folder + "time_delay_plot.csv")
    X = []
    Y = []
    X2 = []
    Y2 = []
    params = [0.616706, 0.506873, 0.389255, 0.097741, 1.598991]
    
    def getlogslope(n, params):
        params_copy = params[:]
        params_copy[-1] = n
        ss = cascadesystem_steadystate_simple(30, *params_copy)
        params_copy[2] = 0
        ss_0 = cascadesystem_steadystate_simple(30, *params_copy)
        dynamic_range = np.abs(ss - ss_0)
        #pylab.plot(range(30+1), dynamic_range)
        logslope = (np.log10(dynamic_range[10]) - np.log10(dynamic_range[1]))/(range(30+1)[10] - range(30+1)[1])
        return logslope
        
    for n in linspace(0.6, 10.0, 1000):
        logslope = getlogslope(n, params)
        X.append(n)
        Y.append(logslope)
    X = np.array(X)
    Y = np.array(Y)
    #find points
    def findY(x_array, y_array, x):
        i = np.abs(x_array - x).argmin()
        min_index = np.abs(x)
        return x_array[i], y_array[i]
    
    def findX(x_array, y_array, y):
        i = np.abs(y_array - y).argmin()
        min_index = np.abs(y)
        return x_array[i], y_array[i]
    
    data = []
    data.append(findY(X, Y, 1.5))
    data.append(findY(X, Y, 1.0))
    data.append(findY(X, Y, 2.0))
    
    data2 = []
    data2.append(findX(X, Y, np.log10(0.9)/5))
    data2.append(findX(X, Y, np.log10(0.5)/5))
    Y = 10**Y
    Y = Y - 1
    pylab.xlim([0,3])
    pylab.yscale("Log")
    pylab.show()
    pylab.plot(X, Y)
    pylab.title("Metric")
    pylab.xlim([0,3])
    pylab.show()
    pylab.xscale("Log")
    pylab.plot(X, Y)
    pylab.scatter(zip(*data)[0], zip(*data)[1])
    pylab.scatter(zip(*data2)[0], zip(*data2)[1], marker="v")
    pylab.title("Metric")
    pylab.xlim([0.5,3])
    pylab.show()
    pd.DataFrame(zip(X, Y), columns=['n', 'metric']).to_csv("/Users/Justin/Desktop/metricplot.csv")
   # plotDynamicRange(n_files, savefile = root_folder + "dynamic_range_out_n.csv")
    #plotDynamicRange(b_files, savefile = root_folder + "dynamic_range_out_b.csv")
   # plotDelay([root_folder + "largesample.csv"])#, savefile = root_folder + "time_delay_plot.csv")