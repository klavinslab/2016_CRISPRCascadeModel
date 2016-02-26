import os
import glob
import pylab as plt
import seaborn as sns

from Model_160124 import *


paramset = pd.read_csv(os.path.join(
                        parameterfolder,
                                    "parameter_table_dose_response_fit_unitless.csv"),
                       header=0,
                       index_col=[0,1])
paramset = paramset['mean'].transpose()
print paramset.a
print np.random.choice(paramset.a, size=2)


def flipflop(r, a, k, b, u, n):
    dr = [0, 0, 0, 0]
    dr[0] = float(u[0] - b * r[0])
    dr[1] = float(u[1] - b * r[1])
    dr[2] = float(a[0]/(1+k[0]*r[0]**n + k[3]*r[3]**n) - b * r[2])
    dr[3] = float(a[1]/(1 + k[1]*r[1]**n + k[2]*r[2]**n) - b * r[3])
    return np.array(dr)

dt = 0.1
def inputArray(inputs):
    time_list = []
    input0_list = []
    input1_list = []
    for i in inputs:
        num, input0, input1 = i
        num = int(num/dt)
        time_list = time_list + [dt] * num
        input0_list = input0_list + [input0] * num
        input1_list = input1_list + [input1] * num
    time_list = np.array(time_list).cumsum()
    input0_list = np.array(input0_list)
    input1_list = np.array(input1_list)
    dic = {}
    for i, v in enumerate(zip(time_list, input0_list, input1_list)):
        t, i0, i1 = v
        dic[str(t)] = [i0, i1]
    return dic

def manualIntegrate(init, a, k, b, i, n):
    global t_final, dt
    t_final = np.array(i).sum(axis=0)[0]
    t_array = []
    y_array = []
    inputdic = inputArray(i)
    y_array.append(np.array(init))

    t_array.append(0)
    for t in np.arange(dt, t_final, dt):
        t_array.append(t)
        u = inputdic[str(t)]
        dy = flipflop(y_array[-1], a, k, b, u, n)
        new_y = y_array[-1] + dy * dt
        y_array.append(new_y)
    return pd.DataFrame(data=np.array(y_array), index=t_array)

def randomFlipFlop():
    init = [0, 0, 0, 0]
    a = np.random.choice(paramset.a, size=2)
    k = np.random.choice(paramset.k, size=4)
    b = np.random.choice(paramset.b, size=1)
    n = np.random.choice(paramset.n, size=1)
    ts = 150
    u = float(np.random.choice(paramset.u, size=1))
    leak = 0.0
    inputs = [
        [ts, u, leak],
        [ts*5, leak, leak],
        [ts, leak, u],
        [ts*5, leak, leak],
        [ts, u, leak],
    ]
    d = manualIntegrate(init, a, k, b, inputs, n)
    d.columns = ["R", "S", "Q", "Qbar"]
    d['aQ'] = a[0]
    d['aQbar'] = a[1]
    d['u'] = u
    d['k0'] = k[0]
    d['k1'] = k[1]
    d['k2'] = k[2]
    d['k3'] = k[3]
    d['b'] = float(b)
    d['n'] = float(n)
    # display.display(a, k, b, u, n)
    return d



df_set = []
for i in range(1000):
    d = randomFlipFlop()
    # display.clear_output(wait=True)
    # display.display(i)
    t800 = d.iloc[int(800/dt)]
    t1600 = d.iloc[int(1600/dt)]
    bistable = False
    if t800.Q < t800.Qbar and t1600.Q > t1600.Qbar:
        bistable = True
    d['set'] = i
    d['bistable'] = bistable
    df_set.append(d)
pd.concat(df_set).to_csv("/Users/klavinslab/Desktop/flipflop_set1000.csv")

df = pd.concat(df_set)
S = df.groupby('set')
sm = pd.DataFrame()
for s in S:
    sm = sm.append(s[1].iloc[0], ignore_index=True)
G = sm.groupby('bistable')
# for g in G:
#     display.display(g[1].apply(np.mean))