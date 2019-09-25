import pandas as pd
import numpy as np
import os as os

import math
from sklearn.linear_model import Lasso
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import StratifiedKFold
import cvxopt as opt
from cvxopt import solvers
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from matplotlib.backends.backend_pdf import PdfPages
from numpy.polynomial.polynomial import polyval
from scipy.stats import f

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

    
    
def GSR(alpha, residuals):
    s_vec  = residuals.shape
    T = s_vec[0]
    N = 4

    alpha = np.asmatrix(alpha.reshape((len(alpha), 1)))
    mean_vec = np.mean(residuals[:, -1], axis = 0)
    cov_fac  = np.cov(residuals)
    inv_cov_fac = np.asmatrix(np.linalg.inv(cov_fac[0:4, 0:4]))
    
    F = (T - N - 1)/N* (alpha.T * inv_cov_fac * alpha)/(1 + (mean_vec/np.sqrt(cov_fac[4, 4]))**2)/T;
    pval = 1-f.cdf(F, N, T - N - 1)
    
    return F, pval

         
    
    
    
    
    
    
#%%

np.random.seed(2019)


path_dir = os.path.dirname(os.path.realpath('C:/Users/aangeli/Dropbox/MQF/Spring 2019/Asset Pricing/A3/FamaFrench.py'))
path_data= os.path.join(path_dir + '/Problem3.5_data(csv).csv')

df=read_data(path_data)
#
##Find indices to split the data
df.rename( columns={'Unnamed: 0':'Date'}, inplace=True )
split_loc = df.loc[df['Date'].isin([192607, 196312, 196401, 201606])].index

excess_df = df - df.iloc[:, -1]

bound = 5 # Bound = includes the market return, bound = 4 does not

returns = df.iloc[:, 1:bound + 1] 
means_ff = returns.iloc[:, 0:4].mean(axis = 0)
rf = df.iloc[:,6]
mean_rf = rf.mean()
colnames=returns.columns
ex_ret_complete = pd.DataFrame((returns.values - rf.values[:, np.newaxis])/100, columns = colnames)

for h in range(2): 
    # Split into the two samples
    ex_ret = ex_ret_complete.iloc[split_loc[0 + 2*h]:split_loc[1+ 2*h]+1]
    returns_split = returns.iloc[split_loc[0 + 2*h]:split_loc[1+ 2*h]+1]/100
    df_mean_exces_ret = ex_ret.iloc[:, 0:bound].mean(axis = 0)
    rf_split = returns_split.iloc[:, -1]
    rf_split_mean = rf_split.mean()
    
  
    cov  = returns_split.cov()
    ret, risks, coeffs, mw, tw = frontier(df_mean_exces_ret.iloc[ 0:4], cov.iloc[0:4, 0:4].values, mu_vec = None, shortselling = False)
    print(array2latex(np.round(np.asmatrix(df_mean_exces_ret.values*100), 2)))
    print(array2latex(np.round(np.asmatrix(cov)*10000)))
    
    x = np.linspace(np.max(ret), 0, 1000)
    y = polyval(x, coeffs[::-1])
    
    market_point = np.array([np.sqrt(np.asscalar(np.matmul(np.transpose(np.matmul(cov.iloc[0:4, 0:4].values, mw)), mw))), \
                                       np.matmul(df_mean_exces_ret[ 0:4].values, mw)])
        
    tangent_point = np.array([np.sqrt(np.asscalar(np.matmul(np.transpose(np.matmul(cov.iloc[0:4, 0:4].values, tw)), tw))), \
                              np.matmul(df_mean_exces_ret[0:4].values, tw)])
        
    
    port_vec = np.zeros((5, 2))
    port_names = ['MV Frontier', 'CAL', 'GMV', 'TP' , 'SL', 'SH', 'BL', 'BH', 'Market']
    port_colors = 'bymcg'
    
    for i in range(5):
        port_vec[i,0] = cov.iloc[i, i]
        port_vec[i,1] = df_mean_exces_ret.iloc[i]
        
    
    s = 2.3
    si = 80
    plt.plot(risks, ret)
    plt.plot([0, tangent_point[0]*s], [0, tangent_point[1]*s], '--b') #CAL
    
    plt.plot(market_point[0], market_point[1],  'r*') #GMV
    plt.plot(tangent_point[0], tangent_point[1],'k*') #Tang


    for i in range(5):
        plt.plot(np.sqrt(np.asscalar(port_vec[i, 0])), port_vec[i, 1], port_colors[i] + 'o')
        
    invS = np.linalg.inv(cov.iloc[0:4, 0:4].values)
    
    plt.legend(port_names)
    plt.savefig('capm' + str(h) +'.pdf')
    plt.show()
    
    print('(%i) SR of Tangent Portfolio: %.3f' % (h+ 1, tangent_point[0]/tangent_point[1]))
    print('(%i) SR of Market Portfolio: %.3f' % (h+1, market_point[0]/market_point[1]))
    
    #(ii)
    twv = tw.reshape((len(tw), 1))
    mwv = mw.reshape((len(mw), 1))
    ff = np.eye(4) 
    
    weights_arr = ff
    
    cov_ii = returns_split.iloc[:, 0:4].cov().values

    beta_vec = np.zeros(weights_arr.shape[1])
    alpha_vec =  np.zeros(weights_arr.shape[1])
    market_sigma_vec = np.matmul(cov_ii, mwv)
    market_var = (np.asscalar(np.matmul(mw[:, 0], market_sigma_vec)))
    port_returns  =  np.zeros(weights_arr.shape[1])
    ff_split_mean = returns_split.iloc[:, 0:4].mean()
    
    ret_means = returns_split.mean(axis = 0)
    for i in range(weights_arr.shape[1]):
        port_returns[i] = np.asscalar(np.matmul(weights_arr[:, i], ff_split_mean.values.reshape(4, 1)))
        beta_vec[i] = cov.iloc[i, 4]/cov.iloc[4, 4] #np.asscalar(np.matmul(weights_arr[:, i], market_sigma_vec)))/market_var
        alpha_vec[i] = ret_means.iloc[i]- rf_split_mean - beta_vec[i]*(ret_means.iloc[-1] - rf_split_mean)
                        
                        
    betas = pd.DataFrame(np.round(beta_vec, 2), index = ['SL', 'SH', 'BL', 'BH'], columns = [''])
    alphas = pd.DataFrame(np.round(alpha_vec, 2), index = ['SL', 'SH', 'BL', 'BH'], columns = [''])
    
    beta_ret_line = np.polyfit(beta_vec,  port_returns[0:4], 1)
    
    x_beta = np.linspace(np.min(beta_vec)*0.9, np.max(beta_vec)*1.1, 1000)
    y_beta = polyval(x_beta, beta_ret_line[::-1])
    
    print(array2latex(np.round(np.asmatrix(alpha_vec*100), 2)))
    plt.plot(x_beta, y_beta, '--')
    plt.scatter(beta_vec, port_returns)
    plt.savefig('alpha' + str(h) +'.pdf')

    plt.show()


    print('\nBeta Vector:', betas)  
    print('\nAlpha Vector:', alphas) 
    
    F, p = GSR(alpha_vec, returns_split.values - rf_split_mean)
    
    print('p = ' , p, 'F = ', F)
    
    



#%%
    
# *** (i) ***
ret  = df.iloc[:, 1:bound + 1] 
rf = df.iloc[:, -1]

split_loc2 = df.loc[df['Date'].isin([196401, 201606, 196401, 199312,  199401, 201606])].index
HML = ret.iloc[:, 1] - ret.iloc[:, 0]
HML_cumsum = np.cumsum(HML)

HML_ma100 = [(HML_cumsum.iloc[i + 100] - HML_cumsum.iloc[i])/100 \
                      for i in range(split_loc2[0], len(HML_cumsum)-100)]

HML_ma100_1 = HML_ma100[split_loc2[0] - split_loc2[0] : split_loc2[1] - split_loc2[0]+1]
HML_ma100_2 = HML_ma100[split_loc2[2] - split_loc2[0] : split_loc2[3] - split_loc2[0]+1]


plt.plot(HML_ma100_1)
plt.plot(HML_ma100_2)
plt.savefig('bi.pdf')
plt.show()

# *** (ii) ***
HML1 = HML.iloc[split_loc2[0]: split_loc2[1]+1]
HML2 = HML.iloc[split_loc2[2]: split_loc2[3]+1]
rf1 = rf.iloc[split_loc2[0]: split_loc2[1]+1]
rf2 = rf.iloc[split_loc2[2]: split_loc2[3]+1]

rf1m = rf1.mean()
rf2m = rf2.mean()

ERM = df.iloc[:, 5] - df.iloc[:, 6]
ERM1 = ERM.iloc[split_loc2[0]: split_loc2[1]+1]
ERM2 = ERM.iloc[split_loc2[2]: split_loc2[3]+1]

HML_mean, HML_std, HML_sr = mean_std_sr(HML)
HML1_mean, HML1_std, HML1_sr = mean_std_sr((HML1-rf1m)*100)
HML2_mean, HML2_std, HML2_sr = mean_std_sr((HML2-rf2m)*100)

ERM_mean,ERM_std, ERM_sr = mean_std_sr(ERM)
ERM1_mean,ERM1_std, ERM1_sr = mean_std_sr(ERM1)
ERM2_mean,ERM2_std, ERM2_sr = mean_std_sr(ERM2)


HML1_mss = np.array([HML1_mean, HML1_std, HML1_sr])
HML2_mss = np.array([HML2_mean, HML2_std, HML2_sr])
ERM1_mss = np.array([ERM1_mean,ERM1_std, ERM1_sr])
ERM2_mss = np.array([ERM2_mean,ERM2_std, ERM2_sr])

bii = np.round(np.stack((HML1_mss, HML2_mss, ERM1_mss, ERM2_mss)), 3)

print(np.round(HML1_mean, 2), np.round(HML1_std, 2), np.round(HML1_sr, 2))
print(np.round(HML2_mean, 2), np.round(HML2_std, 2), np.round(HML2_sr, 2))


# *** (iii) ***

HMLQ = pd.Series([sum(HML.iloc[4*i:4*(i+1)]) for i in range(int(len(HML)/4))])
ERMQ = pd.Series([sum(ERM.iloc[4*i:4*(i+1)]) for i in range(int(len(ERM)/4))])

HMLQ1 = pd.Series([sum(HML1.iloc[4*i:4*(i+1)]) for i in range(int(len(HML1)/4))])
ERMQ1 = pd.Series([sum(ERM1.iloc[4*i:4*(i+1)]) for i in range(int(len(ERM1)/4))])

HMLQ2 = pd.Series([sum(HML2.iloc[4*i:4*(i+1)]) for i in range(int(len(HML2)/4))])
ERMQ2 = pd.Series([sum(ERM2.iloc[4*i:4*(i+1)]) for i in range(int(len(ERM2)/4))])

HMLQ_ac12 = [HMLQ.autocorr(i ) for i in range(13)]
ERMQ_ac12 = [ERMQ.autocorr(i) for i in range(13)]


plot_acf(HMLQ1, lags = range(13))
plt.savefig('HML0.pdf')


plot_acf(HMLQ2, lags = range(13))
plt.savefig('HML1.pdf')


#plot_acf(HMLQ1, lags = range(13))
#plot_acf(HMLQ2, lags = range(13))
#
#
#plot_acf(ERMQ1, lags = range(13))
#plot_acf(ERMQ2, lags = range(13))
plt.savefig('test.pdf')









