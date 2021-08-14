# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 2021

@Purpose: Solving the l1 GNIO optimization
@Version: 1.0

@Data structure: SortedList from 
https://github.com/grantjenks/python-sortedcontainers

@author: Xuyu Chen, PhD candidate
@Insititution: School of Mathematical Sciences, Fudan University
@The implementation is not perfect, for any bug, please email chenxy18@fudan.edu.cn
"""

from collections import deque
import numpy as np

class func(object):
    
    def __init__(self):
        self.bp = [0] # breakpoints
        self.ds = [[0,0]] # difference of coefficients
        self.sl = [0,0] # leftmost function 
        self.sr = [0,0] # rightmost function


    def htog(self,w,y):
        self.sl = [ self.sl[0] + w, self.sl[1] - 2*w*y]
        self.sr = [ self.sr[0] + w, self.sr[1] - 2*w*y]
        
    def gtoh(self,lbd,mu):
        k = len(self.bp)
        
        if lbd < np.inf and mu < np.inf:
            # finding b^-
            left = self.sl
            for p in range(k):
                s = 2*left[0]*self.bp[p] + left[1]
                if s > -lbd:
                    break
                left = [left[0] + (self.ds[p][0]) , left[1] + (self.ds[p][1]) ]
            leftposi = self.bp[p] - (s + lbd)/(2*left[0])
            
            # finding b^+
            right = self.sr
            for q in range(k-1,-1,-1):
                t = 2*right[0]*self.bp[q] + right[1]
                if t < mu:
                    break
                right = [ right[0] - (self.ds[q][0]) , right[1] - (self.ds[q][1]) ]
            rightposi = self.bp[q] + (mu - t)/(2*right[0])
            
            # Operating sequences and update new function 
            if leftposi > self.bp[k-1] or rightposi < self.bp[0]:
                self.bp = [leftposi , rightposi]
                self.sl = [0 , -lbd]
                self.sr = [0 , mu]
                self.ds = [ [left[0] - self.sl[0] , left[1] - self.sl[1] ] , [self.sr[0] - right[0], self.sr[1] - right[1]] ]
                #print('out of range')
            elif q == (p-1):
                self.bp = [leftposi , rightposi]
                self.sl = [0 , -lbd]
                self.sr = [0 , mu]
                self.ds = [[left[0] - self.sl[0],left[1] - self.sl[1]] , [self.sr[0] - right[0], self.sr[1] - right[1]] ]
                #print('q = p -1 ')
            else:
                self.bp = self.bp[p:(q+1)]
                self.ds = self.ds[p:(q+1)]
                self.sr = [0 , mu]
                self.bp.append(rightposi)
                self.ds.append([self.sr[0] - right[0], self.sr[1] - right[1]])
                self.sl = [0 , -lbd]
                self.bp.insert(0,leftposi)
                self.ds.insert(0,[left[0] - self.sl[0],left[1] - self.sl[1]])
                #print('normal')
        
        elif lbd == np.inf and mu == np.inf:
            leftposi = -np.inf
            rightposi = np.inf
        
        elif lbd == np.inf and mu < np.inf:
            leftposi = -np.inf
            right = self.sr
            for q in range(k-1,-1,-1):
                t = 2*right[0]*self.bp[q] + right[1]
                if t < mu:
                    break
                right = [ right[0] - (self.ds[q][0]) , right[1] - (self.ds[q][1]) ]
            rightposi = self.bp[q] + (mu - t)/(2*right[0])
            if  rightposi < self.bp[0]:
                self.bp =  [rightposi]
                self.sr = [0 , mu]
                self.ds = [ [self.sr[0] - right[0], self.sr[1] - right[1] ] ]
            else:
                self.bp = self.bp[0:(q+1)]
                self.ds = self.ds[0:(q+1)]
                self.bp.append( rightposi)
                self.sr = [0 , mu]
                self.ds.append( [self.sr[0] - right[0], self.sr[1] - right[1] ] )
        else:
            rightposi = np.inf 
            left = self.sl
            for p in range(k):
                s = 2*left[0]*self.bp[p] + left[1]
                if s > -lbd:
                    break
                left = [left[0] + (self.ds[p][0]) , left[1] + (self.ds[p][1]) ]
            leftposi = self.bp[p] - (s + lbd)/(2*left[0])
            if leftposi > self.bp[-1]:
                self.bp = [ leftposi ]
                self.sl = [ 0,-lbd]
                self.ds = [ [left[0] - self.sl[0],left[1] - self.sl[1]] ]
            else:
                self.bp = self.bp[p:]
                self.ds = self.ds[p:]
                self.sl = [0 , -lbd]
                self.bp.insert(0,leftposi)
                self.ds.insert(0,[left[0] - self.sl[0],left[1] - self.sl[1]])

        
        
        return tuple([leftposi,rightposi])
    
    def fmin(self):
        k = len(self.bp)
        cur = self.sl
        for i in range(k):
            s = 2*cur[0]*self.bp[i] + cur[1]
            if s > 0 :
                x = - cur[1] / (2 * cur[0])
                break
            cur = [cur[0] + self.ds[i][0],cur[1] + self.ds[i][1]]
        if s <= 0:
           x = - cur[1] / (2 * cur[0]) 
        
        return x
    
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
        