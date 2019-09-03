#!/usr/bin/env python
# coding: utf-8

#Author: Nitesh Sukhwani

import numpy as np

#  Smith-Waterman Algorithm
#  Finding local Alignment between two DNA sequence

'''
Smith Waterman's Algorithm for Local alignment 
seq1 = first nt sequences
seq2 = second nt sequences
match = match score
match = mismatch score
gap = gap penalty
'''

def SW_Mat(seq1,seq2,match,mismatch,gap):
    n1 = len(seq1)
    n2 = len(seq2)
    m = 0
    mat = np.zeros((n1+1,n2+1))
    for i in range(1,n1+1):
        for j in range(1,n2+1):
            if seq1[i-1]==seq2[j-1]:
                mat[i,j]=max(0,mat[i-1,j-1]+match,mat[i-1,j]+gap,mat[i,j-1]+gap)
                if mat[i,j]>m:
                    m = mat[i,j]
                    col_i = i
                    col_j = j
            else:
                mat[i,j]=max(0,mat[i-1,j-1]+mismatch,mat[i-1,j]+gap,mat[i,j-1]+gap)
                if mat[i,j]>m:
                    m = mat[i,j]
                    col_i = i
                    col_j = j
    return mat,(col_i,col_j)


def local_alignment(seq1, seq2,match,mismatch,gap):
    s =""
    mat,index = SW_Mat(seq1,seq2,match,mismatch,gap)
    lst=[]
    i = index[0]
    j = index[1]
    n = mat[i,j]
    while n!=0:
        i-=1
        j-=1
        lst.append(seq1[i])
        n = mat[i,j]
    lst.reverse()
    s = ''.join(lst)
    return s


b = 'aacctatagct'
a = 'gcgatatabubjnk'
print('global alignment of "',a,'" and "',b,'" is = ',local_alignment(a,b,1,-1,-1))

