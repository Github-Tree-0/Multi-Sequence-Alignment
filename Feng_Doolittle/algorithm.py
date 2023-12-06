# Implementation of Feng-Doolittle algorithm

import numpy as np
from collections import Counter

class MSA:
    def __init__(self, alphabet, score_matrix, f=1):
        self.alphabet = alphabet
        # the score matrix is a numpy matrix
        # the rows and columns are in the same order as the alphabet
        self.score_matrix = score_matrix
        assert '-' in self.alphabet
        self.alphabet2index = {c: i for i, c in enumerate(self.alphabet)}
        self.f = f
        
    def compute_profile(self, alignment):
        profile = np.zeros((len(self.alphabet), len(alignment[0])))
        for seq in alignment:
            for i in range(len(seq)):
                profile[self.alphabet2index[seq[i]]][i] += 1
        profile /= len(alignment)

        return profile
    
    def align_profile(self, algn1, algn2, cal_gaps=False, backtrace_flag=False):
        profile1 = self.compute_profile(algn1)
        profile2 = self.compute_profile(algn2)
        
        len1, len2 = len(algn1[0]), len(algn2[0])
        M = np.zeros((len1 + 1, len2 + 1))
        backtrace = np.zeros((len1 + 1, len2 + 1))

        v_gap = np.zeros(len(self.alphabet))
        v_gap[self.alphabet2index['-']] = 1
        for i in range(len1 + 1):
            for j in range(len2 + 1):
                if i != 0 or j != 0:
                    v1, v2 = profile1[:,i-1], profile2[:,j-1]
                    left = float('-inf') if i == 0 else M[i - 1, j] + \
                        np.sum(v1 * (self.score_matrix @ v_gap))
                    up = float('-inf') if j == 0 else M[i, j - 1] + \
                        np.sum(v_gap * (self.score_matrix @ v2))
                    diag = float('-inf') if i == 0 or j == 0 else M[i - 1, j - 1] + \
                        np.sum(v1 * (self.score_matrix @ v2))
                    # use argmin, record the path
                    index = np.argmax([left, up, diag])
                    M[i, j] = [left, up, diag][index]
                    backtrace[i, j] = index

        if backtrace_flag or cal_gaps:
            if cal_gaps:
                num_gaps = 0

            alignment = ['' for _ in range(len(algn1) + len(algn2))]
            # backtrace
            i, j = len1, len2
            while i > 0 or j > 0:
                if backtrace[i, j] == 0:
                    if cal_gaps:
                        num_gaps += 1
                    for indx in range(len(algn1)):
                        alignment[indx] = algn1[indx][i-1] + alignment[indx]
                    for indx in range(len(algn1), len(algn1)+len(algn2)):
                        alignment[indx] = '-' + alignment[indx]
                    i -= 1
                elif backtrace[i, j] == 1:
                    if cal_gaps:
                        num_gaps += 1
                    for indx in range(len(algn1)):
                        alignment[indx] = '-' + alignment[indx]
                    for indx in range(len(algn1), len(algn1)+len(algn2)):
                        alignment[indx] = algn2[indx-len(algn1)][j-1] + alignment[indx]
                    j -= 1
                else:
                    for indx in range(len(algn1)):
                        alignment[indx] = algn1[indx][i-1] + alignment[indx]
                    for indx in range(len(algn1), len(algn1)+len(algn2)):
                        alignment[indx] = algn2[indx-len(algn1)][j-1] + alignment[indx]
                    i -= 1
                    j -= 1
            if cal_gaps:
                return alignment, M[-1, -1], num_gaps
            
            return alignment, M[-1, -1]

        return None, M[-1, -1]

    def compute_D(self, sequences):
        num_seq = len(sequences)
        D = np.zeros((num_seq, num_seq))
        for i in range(num_seq):
            for j in range(i+1, num_seq):
                _, score = self.align_profile([sequences[i]], [sequences[j]])

                _, score_aa = self.align_profile([sequences[i]], [sequences[i]])
                _, score_bb = self.align_profile([sequences[j]], [sequences[j]])
                score_max = 0.5 * (score_aa + score_bb)

                alignment, _, num_gaps = self.align_profile([sequences[i]], [sequences[j]], cal_gaps=True)
                L = len(alignment[0])
                counter_a = Counter(sequences[i])
                counter_b = Counter(sequences[j])
                score_rand = np.sum([self.score_matrix[self.alphabet2index[x]][self.alphabet2index[y]] * counter_a[x] * counter_b[y] 
                                     for x in self.alphabet for y in self.alphabet]) / L - num_gaps
                
                D[i, j] = -np.log((score - score_rand) / (score_max - score_rand) * self.f)
                D[j, i] = D[i, j]

        return D
                

    def compute_guide_tree(self, sequences):
        pass

    def compute_msa(self, sqeuences):
        pass