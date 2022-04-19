# -*- coding: utf-8 -*-
"""
Created on Wed Feb 2 11:31:00 2022

Python Library for DP4GNIO -- A fast and robust solver for Generalized 
Nearly Isotonic Optimization problems

Find documentations on Github: https://github.com/PupuGod/DP4GNIO"

Warning: this library only works on Linux OS, for other versions, please visit
our Github homepage

@author: Xuyu Chen 
@Insititution: School of Mathematical Sciences, Fudan University
@Email: chenxy18@fudan.edu.cn
"""

import numpy as np
from ctypes import *
import time 


def _check(data, w, lbd, mu):
    flag1 = int(isinstance(data, np.ndarray))
    flag2 = int(isinstance(w, np.ndarray))
    flag3 = int(isinstance(lbd, np.ndarray))
    flag4 = int(isinstance(mu, np.ndarray))
    if (flag1*flag2*flag3*flag4):
        pass
    else:
        raise TypeError("Inputs must be 'ndarrays' ")
    
    n1 = len(data)
    n2 = len(w)
    n3 = len(lbd)
    n4 = len(mu)
    
    if n1 != n2:
        raise ValueError("Length of input arrays (data,w,lbd,mu) are incompatible")
    elif (n1 - n3) != 1:
        raise ValueError("Length of input arrays (data,w,lbd,mu) are incompatible")
    elif n3 != n4:
        raise ValueError("Length of input arrays (data,w,lbd,mu) are incompatible")
    else:
        pass
    
    for i in range(n3):
        flag_lbd = (lbd[i] >= 0)
        flag_mu = (mu[i] >= 0)
        flag_w = (w[i] >= 0) and (w[i] < np.inf)
        flag_data = (np.abs(data[i]) < np.inf)
        
        flag_all = flag_lbd * flag_data * flag_w * flag_mu
        
        if flag_all:
            pass
        else:
            raise ValueError("Invalid values, possibly negative value in (w,lbd,mu) or infinity in 'data'")
    
    flag_w = (w[n3] >= 0) and (w[n3] < np.inf)
    flag_data = (np.abs(data[n3]) < np.inf)
    if (flag_w * flag_data):
        pass
    else:
        raise ValueError("Invalid values, possibly negative value in (w,lbd,mu) or infinity in 'data'")
    
    return True


def Copyright():
    print("Thanks for using our software")
    print("Find documentations on Github: https://github.com/PupuGod/DP4GNIO")
    print("Copyright: Xuyu Chen and Xudong Li")
    print("Contact us: chenxy18@fudan.edu.cn")

def l1gnio(data, w, lbd, mu, check = True, display = True):
    """
    Parameters: data: ndarray
                An array of finite real numbers, as the raw data to be fitted
    
                w: ndarray
                An array of positive real numbers, as the weight in the L1 loss function
    
                lbd: ndarray
                An array of non-negative real number, which controls the increasing trend
    
                mu: ndarray
                An array of non-negative real number, which controls the decreasing trend
    
                check: bool, optional
                If this is set to be True, the software will automatically check whether 
                the inputs data,w,lbd,mu is valid (may increase computation expenses for
                large-scale cases ).  
    
                display: bool. optional
                If this is set to be True, the computation results and timing results 
                will be displayed.
    
    Returns:    results:  tuple
                A tuple wih two elements, tuple[0] is the computed solution, tuple[1] 
                is the corresponding objective value
    """
    if display:
        print("==========  The DP4GNIO software  ==========")
    
    if check:

        check_result = _check(data, w, lbd, mu)
        # print("Valid input")
        # print("============================================")
    
    pDll = CDLL("./lib_DP4GNIO.so")
    n = len(data)
    solution = np.zeros(n)
    
    pDll.l1gnio.argtypes = [POINTER(c_double*n), POINTER(c_double*n), \
                        POINTER(c_double*(n-1)), POINTER(c_double*(n-1)), \
                            c_int , POINTER(c_double * n)]
    
    data_c = (c_double * n)(*data)    
    w_c =  (c_double * n)(*w)    
    lbd_c = (c_double * (n-1))(*lbd)    
    mu_c = (c_double * (n-1))(*mu)    
    solution_c = (c_double * n)(*solution)      
    n_c = c_int(n)
    if display:
        print("Problem size: " + str(n))
        print("Start computing L1GNIO solution ....")
    time_start=time.time()
    pDll.l1gnio(data_c,w_c,lbd_c,mu_c,n_c, solution_c)
    time_end=time.time()
    
    solution = np.array(solution_c)
    
    obj_val = 0.0
    pDll.obj_l1.argtypes = [POINTER(c_double*n), POINTER(c_double*n), \
                        POINTER(c_double*(n-1)), POINTER(c_double*(n-1)), \
                            c_int , POINTER(c_double * n), POINTER(c_double * 1)]
    obj_pointer = (c_double * 1)(obj_val)
    
    pDll.obj_l1(data_c,w_c,lbd_c,mu_c,n_c, solution_c,obj_pointer)
    obj_val = np.array(obj_pointer)[0]
    if display:
        print("Optimal solution obtained !")
        print("Objective value: %.9e" %obj_val )
    
        print("Computation time: %.4f" %(time_end-time_start), "s")    
    
    return (solution,obj_val)


def l2gnio(data, w, lbd, mu, check = True, display = True):
    """
    Parameters: data: ndarray
                An array of finite real numbers, as the raw data to be fitted
    
                w: ndarray
                An array of positive real numbers, as the weight in the L1 loss function
    
                lbd: ndarray
                An array of non-negative real number, which controls the increasing trend
    
                mu: ndarray
                An array of non-negative real number, which controls the decreasing trend
    
                check: bool, optional
                If this is set to be True, the software will automatically check whether 
                the inputs data,w,lbd,mu is valid (may increase computation expenses for
                large-scale cases ).  
    
                display: bool. optional
                If this is set to be True, the computation results and timing results 
                will be displayed.
    
    Returns:    results:  tuple
                A tuple wih two elements, tuple[0] is the computed solution, tuple[1] 
                is the corresponding objective value
    """
    if display:
        print("==========  The DP4GNIO software  ==========")
    
    if check:

        checkresult = _check(data, w, lbd, mu)
        # print("Valid input")
        # print("============================================")
    
    pDll = CDLL("./lib_DP4GNIO.so")
    n = len(data)
    solution = np.zeros(n)
    
    pDll.l1gnio.argtypes = [POINTER(c_double*n), POINTER(c_double*n), \
                        POINTER(c_double*(n-1)), POINTER(c_double*(n-1)), \
                            c_int , POINTER(c_double * n)]
    
    data_c = (c_double * n)(*data)    
    w_c =  (c_double * n)(*w)    
    lbd_c = (c_double * (n-1))(*lbd)    
    mu_c = (c_double * (n-1))(*mu)    
    solution_c = (c_double * n)(*solution)      
    n_c = c_int(n)
    if display:
        print("Problem size: " + str(n))
        print("Start computing L2GNIO solution ....")
    time_start=time.time()
    pDll.l2gnio(data_c,w_c,lbd_c,mu_c,n_c, solution_c)
    time_end=time.time()
    
    solution = np.array(solution_c)
    
    obj_val = 0.0
    pDll.obj_l1.argtypes = [POINTER(c_double*n), POINTER(c_double*n), \
                        POINTER(c_double*(n-1)), POINTER(c_double*(n-1)), \
                            c_int , POINTER(c_double * n), POINTER(c_double * 1)]
    obj_pointer = (c_double * 1)(obj_val)
    
    pDll.obj_l2(data_c,w_c,lbd_c,mu_c,n_c, solution_c,obj_pointer)
    obj_val = np.array(obj_pointer)[0]
    if display:
        print("Optimal solution obtained !")
        print("Objective value: %.9e" %obj_val )    
        print("Computation time: %.4f" %(time_end-time_start), "s")    
    
    return (solution, obj_val)

if __name__ == '__main__':
    Copyright()








