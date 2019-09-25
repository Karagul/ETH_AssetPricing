# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 12:59:44 2019

@author: aangeli
"""
import pandas as pd
import numpy as np
import os as os
import math
import cvxopt as opt
from cvxopt import solvers
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from matplotlib.backends.backend_pdf import PdfPages
from numpy.polynomial.polynomial import polyval
from scipy.stats import f
from datetime import datetime



solvers.options['show_progress'] = False

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

def PolyCoefficients(x, coeffs):
    """ Returns a polynomial for ``x`` values for the ``coeffs`` provided.

    The coefficients must be in ascending order (``x**0`` to ``x**o``).
    """
    o = len(coeffs)
    print(f'# This is a polynomial of order {ord}.')
    y = 0
    for i in range(o):
        y += coeffs[-1*i]*x**i
    return y


def frontier(mean_vec, cov, mu_vec = None, shortselling = True):
    n = len(cov)
    if mu_vec == None:
        N = 100
        mu_vec = [(t-10)/2500 for t in range(N)]
    
    # Convert to cvxopt matrices
    Sigma = opt.matrix(cov)
    mubar = opt.matrix(mean_vec)
    ones_vec = np.ones(n);


    # Create constraint matrices
    G = -opt.matrix(np.eye(n)) 
    h = opt.matrix(0.0, (n ,1))
    q = opt.matrix(0.0, (n, 1))
    

# negative n x n identity matrix
    if not(shortselling):
        G = None #opt.matrix(0.0 ())
        h = None 
        
    A = opt.matrix(np.concatenate((np.matrix(ones_vec), np.matrix(mean_vec))))    
    portfolios = [solvers.qp(Sigma, q, G = G, h = h, 
                             A =  A, b = opt.matrix([1, mu]))['x'] for mu in mu_vec]



    ## CALCULATE RISKS AND RETURNS FOR FRONTIER
    returns = [np.asscalar(np.dot(mubar.T, x)) for x in portfolios]
    risks = [np.asscalar(np.sqrt(np.dot(x.T, Sigma*x))) for x in portfolios]
    ## CALCULATE THE 2ND DEGREE POLYNOMIAL OF THE FRONTIER CURVE
    m1 = np.polyfit(returns, risks, 2)
    x1 = np.sqrt(m1[2] / m1[0])
    # CALCULATE THE OPTIMAL PORTFOLIO
    wt = solvers.qp(opt.matrix(Sigma), q, G = G, h = h, A = A, b = opt.matrix([1, x1]))['x']
    coef = m1
    inv_cov = np.linalg.inv(cov)
    gmv = np.matmul(inv_cov, np.ones((len(cov), 1)))
    gmv = gmv/np.sum(gmv)
    tanp = np.matmul(inv_cov, mean_vec)
    tanp = tanp/np.sum(tanp)
    
    
    
    x2 = m1[1]/(2*m1[0])
    wtan   = solvers.qp(opt.matrix(Sigma), q, G = G, h = h, A = A, b = opt.matrix([1, x2]))['x']
    
    wtan = np.linalg.lstsq(Sigma, mubar)[0]
    wtan = wtan/np.matmul(np.ones(len(wtan)), wtan)
    return returns, risks, coef, gmv, tanp

def mean_std_sr(vec):
    m = np.mean(vec)
    s = np.std(vec)
    sr = m/s
    return m, s, sr


def array2latex(arr):
    dim = arr.shape
    s = ""
    for l in range(dim[0]):
        for i in range(dim[1]-1):
            s = s + str(arr[l, i]) + " & "
        s = s + str(arr[l, -1])
        
        s = s + " \\\ \n"
    return s

    
    
def GSR(alpha, residuals, factors):
    s_vec  = residuals.shape
    T = s_vec[0]
    N = s_vec[1]
    L = factors.shape[1]
    alpha = np.asmatrix(a.reshape((len(a), 1)))
    mean_vec = np.mean(factors, axis = 0)
    mean_vec = np.asmatrix(np.reshape(mean_vec, (len(mean_vec), 1)))
    cov_res  = np.cov(res)
    cov_fac  = np.cov(factors)
    inv_cov_fac = np.asmatrix(np.linalg.inv(cov_fac))
    inv_cov_res = np.asmatrix(np.linalg.inv(cov_res))
    F = (T/N) * ((T-N-L)/(T-L- 1)) * (alpha.T * inv_cov_res * alpha)/(1 + (mean_vec.T *inv_cov_fac * mean_vec));
    pval = 1 - f.cdf(F, N, T - N - L)
    
    return F, pval
    
def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)  

def findSmallestDifference(A, B): 
    # Input: 
    #   A: (Call) d by n array like
    #   B: (Put) d by m array like 
    
    # Column is the Strike Prices 
    
    # NOTE: Vector must be sorted in ASCENDING order
    
    # A: Call
    # B: Putt
    
    # Sort both arrays  
    # using sort function 
       
    a = 0
    b = 0
    pos_a, pos_b
    
    n = A.shape[0]
    m = B.shape[0]
    
    
    # Check if sorted correctly 
    sortedB = all(B[i, 1] <= B[i+1, 1] for i in range(m - 1))
    sortedA = all(A[i, 1] <= A[i+1, 1] for i in range(n-1))
    
    if not(sortedA and sortedB):
        print("Call or Put premiumns are not sorted correctly")
        return
    
  
    # Initialize result as max value 
    result = sys.maxsize 
  
    # Scan Both Arrays upto 
    # sizeof of the Arrays 
    while (a < m and b < n): 
      
        if (abs(A[a] - B[b]) < result): 
            result = abs(A[a] - B[b]) 
            pos_a = a
            pos_b = b
  
        # Move Smaller Value 
        if (A[a] < B[b]): 
            a += 1
  
        else: 
            b += 1
    # return final sma result 
    return pos_a, pos_b

def findFuturePrice(Call, Put):
    Call = np.array(Call)
    Put = np.array(Put)
    
    bound = Call.shape[0] - 1
    shift = bound/2 - 1
    x = 0
    y = 0
    dx = 0
    dy = -1
    
    for i in range(bound**2):
        if (-bound/2 < x <= bound/2):
            # Check if line intersect 
            intersect, fut = LineIntersect(Call[x + shift, :], Call[x + shift + 1, :], 
                                           Put[x+ shift, :], Put[x + shift + 1, :])
            if intersect:
                return fut
            
        if x == y or (x < 0 and x == -y) or (x > 0 and x == 1-y):
            dx, dy = -dy, dx
        x, y = x+dx, y+dy
    

def LineIntersect(a1, a2, b1, b2):
    X1 = a1[0]
    X2 = a2[0]
    X3 = b1[0]
    X4 = b2[0]
    
    Y1 = a1[1] 
    Y2 = a2[1]
    Y3 = b1[1] 
    Y4 = b2[1]
    
    # Calc abcisses 
    I1 = [min(X1,X2), max(X1,X2)]
    I2 = [min(X3,X4), max(X3,X4)]
    
    if (max(X1,X2) < min(X3,X4)) or (max(Y1,Y2) < min(Y3,Y4)):
        return false #There is no mutual abcisses
    
    # Get equations of the lines 
    # Y = AX + B
    A1 = (Y1-Y2)/(X1-X2)
    A2 = (Y3-Y4)/(X3-X4) 
    B1 = Y1-A1*X1
    B2 = Y3-A2*X3
    
    # Parallel 
    if (A1 == A2):
        return False, -1
    
    # Find point of intersecion of lines defined by segment
    Xa = (B2 - B1) / (A1 - A2)
    
    # Check if intersection is out of bounds 
    if (Xa < max( min(X1,X2), min(X3,X4) )) or \
        (Xa > min( max(X1,X2), max(X3,X4) )):
        return False, -1, 
    else: 
        return True, Xa
    
    

    
def spiral(X, Y):
    x = y = 0
    dx = 0
    dy = -1
    for i in range(12**2):
        if (-X/2 < x <= X/2) and (-Y/2 < y <= Y/2):
            print (x, y)
            # DO STUFF...
        if x == y or (x < 0 and x == -y) or (x > 0 and x == 1-y):
            dx, dy = -dy, dx
        x, y = x+dx, y+dy

def integrat(call, put):
    fut
#%%

np.random.seed(2019)

#
path_dir = os.path.dirname(os.path.realpath('C:/Users/aangeli/Dropbox/MQF/Spring 2019/Asset Pricing/Project/SVIX.py'))
path_data= os.path.join(path_dir + '/SVIX2.xls')
SVIX2= pd.read_excel(path_data,  header = None)
#d = pd.to_datetime(SVIX2.iloc[:, 0])
#dd = pd.concat([d, SVIX2.iloc[:, 0]], axis = 1)

#path_data= os.path.join(path_dir + '/option.csv')
#df= pd.read_csv(path_data)

path_data= os.path.join(path_dir + '/df2.csv')
df3 = pd.read_csv(path_data)

d0 = pd.to_datetime(20000101, format = '%Y%m%d')
df3 = df3.iloc[:, [0, 1, 2,3, 5, 6, 7, 8, 9, 10, 11]]
for i in range(1, 4):    
    df3.iloc[:, i] = pd.to_datetime(df3.iloc[:, i], format = '%Y%m%d')
for i in range(1, 4):
    df3.iloc[:, i] = (df3.iloc[:, i] - d0).dt.days + 1
df3.iloc[:, 2] = df3.iloc[:, 2] - df3.iloc[:, 1]+ 1
# turn dates to number and find difference 


# turn 

path_data= os.path.join(path_dir + '/volsmile.csv')
vol_raw = pd.read_csv(path_data)

code_vec = np.array([[102456, ]])

for i in range(20):
    s = 26*i
    callp = vol_raw.iloc[13 + s:26 + s]
    putp = vol_raw.iloc[0+ s:13 + s]
    plt.plot(callp.iloc[:, 5], callp.iloc[:, 6], '-o')
    plt.plot(putp.iloc[:, 5], putp.iloc[:, 6], '-o')
    plt.title(str(vol_raw.iloc[s, 2]))
    plt.show()
    
    dif = callp.iloc[:, 6].values - putp.iloc[:, 6].values
    plt.plot(np.abs(dif), '-o')
    plt.show()