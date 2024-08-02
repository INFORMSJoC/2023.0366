import pandas as pd 
import numpy as np
from ARMUL import ARMUL
from sklearn import metrics
from ARMUL import Baselines

pathload = 'Demo_training_data.csv'
data = np.loadtxt(pathload, delimiter=',')
M = 10
d = 200
n_list = np.array([200]*10)
X = list()
z = list()

X_test = list()
z_test = list()
for m in range(M):
    start = sum(n_list[:m])
    end = sum(n_list[:(m+1)])
    idx = range(start,end)
    X.append(data[start:end,1:])
    z.append(data[start:end,0])
data = [X,z]
eta = 0.01 # step size for the optimization algorithm (proximal gradient descent)
T = 500 # number of iterations in optimization
seed = 10000 # random seed

C = [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.05, 0.1, 0.5]
lbd_list = np.array([i*np.sqrt(d / n_list) for i in C ] )


base = Baselines('logistic') # link == 'linear' or 'logistic'
test = ARMUL('logistic') 

base.clustered_train(data, eta_B=eta,T=T)
base.lowrank_train(data, eta_B=eta,eta_Z=eta,T=T)
test.cv(data, lbd_list, model = 'vanilla', n_fold = 5, seed = seed, eta_global = eta, eta_local = eta, eta_B = eta, eta_Z = eta, T_global = 500)


save_path_CLUSTER = 'CLUSTER.csv'
save_path_LOWRANK = 'LOWRANK.csv'
save_path_VANILLA = 'VANILLA.csv'
pd.DataFrame(base.models['clustered'][:,:,0]).to_csv(save_path_CLUSTER,sep=',',index=False,header=None)
pd.DataFrame(base.models['lowrank'][:,:,0]).to_csv(save_path_LOWRANK,sep=',',index=False,header=None)
pd.DataFrame(test.models['vanilla'][:,:,0]).to_csv(save_path_VANILLA,sep=',',index=False,header=None)
