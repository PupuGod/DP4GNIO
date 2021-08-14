# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 2021

@Purpose: Solving the l1 GNIO optimization (as a subroutine)
@Version: 1.0

@Data structure: SortedList from 
https://github.com/grantjenks/python-sortedcontainers

@author: Xuyu Chen, PhD candidate
@Insititution: School of Mathematical Sciences, Fudan University

"""

from sortedcontainers import SortedList
from collections import deque
import numpy as np

class func( object ):
    
    def __init__(self):
        
        self.bp = SortedList( [[0,0]] )
        self.lslope = float(0)
        self.rslope = float(0)
        
    def htog(self,w,y):
        self.bp.add([y, 2*w])        
        self.rslope += w
        self.lslope -= w
        
    def gtoh(self,lbd,mu):
                        
        if -lbd <= self.lslope and mu >= self.rslope:
            leftposi = -np.inf
            rightposi = np.inf
        
        elif -lbd > self.lslope and mu < self.rslope:
            
            left = self.lslope
            while left < -lbd:
                z = self.bp.pop(0)
                leftposi = z[0] # real position 
                left = left + z[1]
            self.bp.add( [leftposi, left + lbd ] )
            self.lslope = -lbd
            
            right = self.rslope
            while right > mu:
                z = self.bp.pop()
                rightposi = z[0] # real position 
                right = right - z[1]
            self.bp.add( [rightposi, mu - right ] )
            self.rslope = mu
                   
        elif -lbd <= self.lslope and mu < self.rslope:
            
            leftposi =  - np.inf 
            
            right = self.rslope
            while right > mu:
                z = self.bp.pop()
                rightposi = z[0] # real position 
                right = right - z[1]
            self.bp.add( [rightposi, mu - right ] )
            self.rslope = mu
            
        else:
            rightposi = np.inf 
            
            left = self.lslope
            while left < -lbd:
                z = self.bp.pop(0)
                leftposi = z[0] # real position 
                left = left + z[1]
            self.bp.add( [leftposi, left + lbd ] )
            self.lslope = -lbd        

        return tuple([leftposi,rightposi])

        
    def fmin(self):
        s = self.rslope
        if s <= 0:
            z = self.bp.pop()
            fmin = z[0]
        else:     
            while s > 0:
                z = self.bp.pop()
                fmin = z[0] # real position 
                s = s - z[1]
        
        return fmin
             
    
    def optimize(self, data, weight, lbd, mu):
        n = len(data)
        b_store = []    
        for i in range(n - 1):
            self.htog( weight[i], data[i])
            posi = self.gtoh(lbd[i],mu[i])
            b_store.append(posi)
        self.htog(weight[-1], data[-1])
        xmin = self.fmin()
        solution = recover(xmin,b_store)
        
        return solution
    
    
def recover(xn,b_store):
    val = deque([xn])
    n = len(b_store)
    for i in range(1,n+1):
        b = b_store[-i]
        if xn < b[0]:
            x = b[0]
        elif xn > b[1]:
            x = b[1]
        else:
            x = xn
        val.appendleft(x)
        xn = x
    val = list(val)
    return val 
