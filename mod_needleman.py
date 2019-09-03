#!/usr/bin/env python
# coding: utf-8

#Author: Nitesh Sukhwani

import numpy as np

#  Modified Needleman-Wunsch Algorithm
#  For finding Semi-Global Alignment

'''
Modified Needlemen Wunsch Algorithm for Semi-Global alignment 
seq1 = first nt sequences
seq2 = second nt sequences
match = match score
match = mismatch score
gap = gap penalty
'''

def modified_needleman(seq1, seq2, match = 1, mismatch = 0 ,gap = -1):
    m = len(seq1)+1
    n = len(seq2)+1
    mat = np.zeros((m,n))
    for i in range(1,m):
        for j in range(1,n):
            if(i != m-1 and j != n-1):
                if(seq1[i-1] == seq2[j-1]):
                    mat[i,j] = max(mat[i-1,j] + gap, mat[i,j-1] + gap,   mat[i-1,j-1] + match)
                else:
                    mat[i,j] = max(mat[i-1,j] + gap, mat[i,j-1] + gap,   mat[i-1,j-1] + mismatch)
            elif(i == m-1 and j != n-1):
                if(seq1[i-1] == seq2[j-1]):
                    mat[i,j] = max(mat[i-1,j] + gap, mat[i,j-1],   mat[i-1,j-1] + match)
                else:
                    mat[i,j] = max(mat[i-1,j] + gap, mat[i,j-1],   mat[i-1,j-1] + mismatch)
            elif(i != m-1 and j == n-1):
                if(seq1[i-1] == seq2[j-1]):
                    mat[i,j] = max(mat[i-1,j], mat[i,j-1] + gap,   mat[i-1,j-1] + match)
                else:
                    mat[i,j] = max(mat[i-1,j], mat[i,j-1] + gap,   mat[i-1,j-1] + mismatch)
            else:
                if(seq1[i-1] == seq2[j-1]):
                    mat[i,j] = max(mat[i-1,j], mat[i,j-1],   mat[i-1,j-1] + match)
                else:
                    mat[i,j] = max(mat[i-1,j], mat[i,j-1],   mat[i-1,j-1] + mismatch)
                    
    return mat



def mod_needleman_alignment(seq1, seq2,match,mismatch,gp):
    mat = modified_needleman(seq1,seq2,match,mismatch,gp)
    i = len(seq1)
    j = len(seq2)
    l1 = len(seq1)
    l2 = len(seq2)
    n = mat[i,j]
    lst1 =[]
    lst2 =[]
    while i!=0 and j!=0:
        if i==l1 or j==l2:
            gap = 0
        else:
            gap = gp
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
            if mat[i,j]+match==mat[i+1,j+1] or mat[i,j]+mismatch==mat[i+1,j+1]:
                lst1.append(seq1[i])
                lst2.append(seq2[j])
        n = mat[i,j]
    lst1.reverse()
    lst2.reverse()
    s1 = ''.join(lst1)
    s2 = ''.join(lst2)
    return s1,s2


a = 'acactgatcg'
b = 'acactg'
print('global alignment of "',a,'" and "',b,'" is = ',mod_needleman_alignment(b,a,1,0,-1))


