"""
Created on Mon Mar 18 13:30:01 2019

@author: Matej Privoznik, Michal Kobak, Alejandro Angeli
@Date: March 18, 2019
"""


import pandas as pd
import numpy as np
import os as os
import math
from sklearn.linear_model import Lasso
from sklearn.model_selection import cross_val_score

#%% Utils functions
def save_solution(csv_file,pred_prob):
	with open(csv_file, 'w') as csv:
		df = pd.DataFrame.from_dict({'Id':range(10000,10000+len(pred_prob)),'y': pred_prob})
		df.to_csv(csv,index = False)

def read_data(csv_file):
    with open(csv_file,'r') as csv:
        df = pd.read_csv(csv)
    return df

def get_data(csv_file):
    with open(csv_file, 'r') as csv:
        df = pd.read_csv(csv)
        df = df.loc[:, df.columns != 'y']
        data = df.drop('Id', axis = 1)
    return data

def get_target(csv_file):
	with open(csv_file, 'r') as csv:
		df = pd.read_csv(csv)
		y = df['y']
	y = np.array(y)
	return y

df_train = pd.read_csv("train.csv", index_col="Id")

X_train=df_train.iloc[:,1:6]
y_train=df_train['y']

x21 = pd.DataFrame(np.ones((1000,1)))

X_features=pd.concat([X_train,pow(X_train,2),pow(math.e,X_train),np.cos(X_train),x21],axis=1)
X_features.columns=['x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13','x14','x15','x16','x17','x18','x19','x20','x21']

lassoreg = Lasso(alpha=0.17, normalize=False, fit_intercept=False, max_iter=10000)
model = lassoreg.fit(X_features, y_train)

np.savetxt('output.csv',model.coef_,fmt='%.12f')

