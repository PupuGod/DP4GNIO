# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14  2021

@purpose: Runfile for both l1 gnio and l2 gnio
@Structure: 
0 functions (line 23--43)
1 test data generation (line 44--47)
2 optimization (line 48--57)
3 plot the solutions (line 58--60)

@Author: Xuyu Chen, PhD candidate
@Insititution: School of Mathematical Sciences, Fudan University
"""

import l1gnio as l1
import l2gnio as l2
import numpy as np
from time import *
import random
import matplotlib.pyplot as plt

#==================================================
def random_data_gen( n ):
    data = list(np.random.uniform( -10, 10, n ) )
    w    = np.random.uniform(0.5,0.5,n)
    lbd  =  np.random.uniform(1,1,n-1) # Set regularzation parameter
    mu   =  np.random.uniform(1,1,n-1)
    
    return data, w, lbd, mu    

def drawresult(data,solution):
    plt.figure()
    plt.xlabel('')
    plt.ylabel('value of data')
    plt.title(('Nearly isotonic regression results') )
    plt.plot(data,color = 'r',linestyle = '', marker = 'o',label = 'data')
    plt.plot(solution,color = 'k',linestyle = '-', marker ='.' ,label = 'Regression')
    plt.legend(loc = 'upper left')    
    plt.show()
    
    
    
#========= S1. Generating the data (default: random data) ==============
n = 1000000
data,w,lbd,mu = random_data_gen(n)

#========= S2. optimization =========

begin_time = time()
print('The computation has started, please wait... ')
f = l2.func()
solution = f.optimize(data,w,lbd,mu)
end_time = time()
print('The solution has been successfully generated! ')
print('The computation time is %f s' % (end_time - begin_time) )

#========= S3. Plot the solution =============
drawresult(data,solution)






