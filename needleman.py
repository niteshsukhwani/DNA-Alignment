#!/usr/bin/env python
# coding: utf-8

#Author: Nitesh Sukhwani


import random
import numpy as np


# Needleman-Wunsch Algorithm 
# Finding Global Alignment between two DNA sequence

'''
Needlemen Wunsch Algorithm for Global alignment 
seq1 = first nt sequences
seq2 = second nt sequences
match = match score
match = mismatch score
gap = gap penalty
'''

def needleman(seq1,seq2,match,mismatch,gap):
    n1 = len(seq1)
    n2 = len(seq2)
    mat = np.zeros((n1+1,n2+1))
    for i in range(n1+1):
        for j in range(n2+1):
            if i==0 and j==0:
                continue
            if i==0:
                mat[i,j]=mat[i,j-1]+gap
                continue
            if j==0:
                mat[i,j]=mat[i-1,j]+gap
                continue
            if seq1[i-1]==seq2[j-1]:
                mat[i,j]=max(mat[i-1,j-1]+match,mat[i-1,j]+gap,mat[i,j-1]+gap)
            else:
                mat[i,j]=max(mat[i-1,j-1]+mismatch,mat[i-1,j]+gap,mat[i,j-1]+gap)
    return mat


def needleman_alignment(seq1, seq2,match,mismatch,gap):
    mat = needleman(seq1,seq2,match,mismatch,gap)
    i = len(seq1)
    j = len(seq2)
    n = mat[i,j]
    lst1 =[]
    lst2 =[]
    while i!=0 and j!=0:
        if mat[i,j]==mat[i-1,j]+gap:
            i-=1
            lst1.append('-')
            lst2.append(seq1[i])
        elif mat[i,j]==mat[i,j-1]+gap:
            j-=1
            lst2.append('-')
            lst1.append(seq2[j])
        else:
            i-=1
            j-=1
            if mat[i,j]+match==mat[i+1,j+1] or mat[i,j]+mismatch==mat[i+1,j+1] :
                lst1.append(seq1[i])
                lst2.append(seq2[j])
        n = mat[i,j]
    lst1.reverse()
    lst2.reverse()
    s1 = ''.join(lst1)
    s2 = ''.join(lst2)
    return s1,s2


b = 'actcg'
a = 'acagtag'
print('global alignment of "',a,'" and "',b,'" is = ',needleman_alignment(a,b,1,0,-1),sep='')


