# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14  2021

@purpose: Test files for both l1 gnio and l2 gnio
@Structure: 
1 test data generation
2 draw regression results
3 extrat solutions and data


@Author: Xuyu Chen, PhD candidate
@Insititution: School of Mathematical Sciences, Fudan University
@The implementation is not perfect, for any bug, please email chenxy18@fudan.edu.cn
"""

import l1gnio as l1
import l2gnio as l2
import numpy as np
from time import *
import random
import matplotlib.pyplot as plt

#1 random data generator
def random_data_gen( n ):
    data = list(np.random.uniform( -10, 10, n ) )
    w    = np.random.uniform(0.5,0.5,n)
    lbd  =  np.random.uniform(1,1,n-1) # Set regularzation parameter
    mu   =  np.random.uniform(1,1,n-1)
    
    return data, w, lbd, mu    


#2   plot solutions for GNIO
def drawresult(data,solution):
    plt.figure()
    plt.xlabel('')
    plt.ylabel('value of data')
    plt.title(('Nearly isotonic regression results') )
    plt.plot(data,color = 'r',linestyle = '', marker = 'o',label = 'data')
    plt.plot(solution,color = 'k',linestyle = '-', marker ='.' ,label = 'Regression')
    plt.legend(loc = 'upper left')    
    plt.show()
    
#========= Generating the data ==============
n = 1000000
data,w,lbd,mu = random_data_gen(n)

'''
# Reading the real-world data
file = open(  'D:/陈旭宇/DP_GIR/project/gold_index_open.txt' )
lines = file.readlines()
file.close()
n = len(lines)
data = np.zeros( n )
for i in range(n):
    data[i] = float( lines[i]  )
w = np.random.uniform(0.5,0.5,n)
lbd = np.random.uniform(1,100,n-1)
mu = np.random.uniform(1,100,n-1)
'''

#========= The test on the gnio subroutine =========

begin_time = time()
print('The computation has started, please wait... ')
f = l2.func()
solution = f.optimize(data,w,lbd,mu)
end_time = time()
print('The solution has been successfully generated! ')
print('The computation time is %f s' % (end_time - begin_time) )

# drawresult(data,solution)


#========= Output solution for further uses ===========

filename = 'C:\\Users\\AI\\Desktop\\zzy\\random_e4.txt'
output_file = np.zeros(5*n-1)
for x in lbd:
    if x>=np.inf:
        x = -1
for x in mu:
    if x>=np.inf:
        x = -1

output_file[0:n] = data
output_file[n:2*n] = w
output_file[2*n:3*n-1] = lbd
output_file[3*n-1:4*n-2] = mu
output_file[4*n-2:5*n-2] = solution
output_file[5*n-2] = end_time - begin_time
file = open(filename, 'w')
for num in output_file:
    file.write(str(num))
    file.write('\n')
file.close()

