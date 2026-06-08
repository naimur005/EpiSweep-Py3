#!/usr/bin/python3
from functools import reduce
import numpy as np
import sys, random, itertools, operator

def distance(x1, x2): #Hausdorff distance
    def cal_h(x1, x2):
        h_dist_list = []
        for x in x1.crd:
            dist_list = []
            for y in x2.crd:
                dist = np.linalg.norm(np.array(x)-np.array(y))
                dist_list.append(dist)
            h_dist_list.append(min(dist_list))
        return max(h_dist_list)
    return max(cal_h(x1, x2), cal_h(x2, x1))

def distance_coverage(x1, x2):
    union = set(x1.in_contact).union(set(x2.in_contact))
    intersect = set(x1.in_contact).intersection(set(x2.in_contact))
    return len(union - intersect)

def recluster(crd, docks, center=None):
    pos = {}
    for d in crd:
        idx = "-".join(d.pos)
        if idx not in pos: 
            pos[idx] = d
        else:
            if pos[idx].episcore < d.episcore:
                pos[idx] = d

    all_crd = pos.values()
    if center is None:
        all_crd.sort(key=lambda x: (x.episcore, x.coverage, x.freq), reverse=True)
        all_crd.sort(key=lambda x: (x.episcore, x.freq, x.coverage), reverse=True)
        all_crd.sort(key=lambda x: (x.coverage, x.episcore, x.freq), reverse=True)
        all_crd.sort(key=lambda x: (x.coverage, x.freq, x.episcore), reverse=True)
        all_crd.sort(key=lambda x: (x.freq, x.episcore, x.coverage), reverse=True)
        all_crd.sort(key=lambda x: (x.freq, x.coverage, x.episcore), reverse=True)
        comb = [[all_crd[0]]]
    else:
        center_pos = "-".join(crd[center[0]].pos)
        comb = [[pos[center_pos]]]

    for i in range(1, 20):
        print ("%d representatives among %d designs"%(i, len(all_crd)))
        episcore_sum = comb[0][0].episcore
        coverage = set(comb[0][0].in_contact)
        comb = comb + [all_crd]*(i-1)
        all_covered = False
        covering_triplet = []
        all_product = itertools.product(*comb)
        len_comb = len(comb)
        for c in all_product:
            if len(set(map(lambda x: "-".join(x.design_pos), c))) != len_comb:
                pass
            else:
                total_episcore = sum(map(lambda x: x.episcore, c))
                total_coverage = reduce(lambda l, m: l.union(m), map(lambda x: set(x.in_contact), c))
                if total_episcore > episcore_sum and len(total_coverage) > len(coverage):
                    if len(total_coverage) == len(docks):
                        all_covered = True
                        covering_triplet.append([c, total_episcore, len(total_coverage)])
        if all_covered:
            break
        else:
            print ("Failed.")
    covering_triplet.sort(key=lambda x: x[1], reverse=True)
    return covering_triplet[0]

def create_distance_matrix(mat, crd, medoid, oth_elm):
    for i, n in enumerate(oth_elm):
        for j, m in enumerate(medoid):
            mat[i, j] = distance(crd[n], crd[m])
    return mat
            
class Kmedoids:
    def __init__(self, k):
        np.random.seed(0)
        self.k = k
    def cluster(self, crd):
        self.crd = crd
        max_iter = 100
        print ("Design indices:")
        print ()
        for n, i in enumerate(self.crd):
            print (n, i.design)
        print ()
        idx = np.arange(self.crd.shape[0])
        np.random.shuffle(idx)
        self.medoid = idx[0:self.k]
        self.oth_elm = idx[self.k::]
        self.d_mat = np.empty((self.oth_elm.shape[0], self.medoid.shape[0]), dtype=float)
        self.d_mat = create_distance_matrix(self.d_mat, self.crd, self.medoid, self.oth_elm)
        
        self.clst = np.argmin(self.d_mat, axis=1)
        self.cost = np.sum(self.d_mat[np.arange(self.d_mat.shape[0]), self.clst])
    
        for i in range(max_iter):
            print
            print ("Attempt(s): %d"%(i+1))
            print
            medoid, oth_elm, clst, cost = self.k_clustering()
            if (cost < self.cost):
                self.medoid, self.oth_elm, self.clst, self.cost = medoid, oth_elm, clst, cost
            else:
                break

    def k_clustering(self):
        d = np.empty((self.oth_elm.shape[0], self.medoid.shape[0]), dtype=float)
        med = np.empty_like(self.medoid)
        oth = np.empty_like(self.oth_elm)
        idx = np.arange(self.oth_elm.shape[0])

        print ("Initial medoids:", self.medoid)
        print ("%d + (%d) members:\n"%(len(self.oth_elm), len(self.medoid)), self.oth_elm)
        print ("Grouping:\n", self.clst)
        print()
        
        med_cur = self.medoid.copy()
        oth_cur = self.oth_elm.copy()
        clst_cur = self.clst.copy()
        cost_cur = self.cost
        
        for i, m in enumerate(self.medoid):
            nn = len(self.oth_elm[self.clst == i])
            print ("Cluster number: %d"%(i+1), "(Number of cluster memebrs: %d)"%nn)
            for j, n in enumerate(self.oth_elm[self.clst == i]):
                med, oth = self.medoid.copy(), self.oth_elm.copy()
                med[i] = n
                tmp = oth[self.clst == i]
                tmp[j] = m
                oth[self.clst == i] = tmp
                d = create_distance_matrix(d, self.crd, med, oth)

                clst = np.argmin(d, axis=1)
                cost = np.sum(d[idx, clst])

                if cost <= cost_cur:
                    med_cur = med.copy()
                    oth_cur = oth.copy()
                    clst_cur = clst.copy()
                    cost_cur = cost
        print ()
        print ("Current medoids:", med_cur)
        print ("%d + (%d) members:\n"%(len(oth_cur), len(med_cur)), oth_cur)
        print ("Grouping:\n", clst_cur)
        print ()
        return med_cur, oth_cur, clst_cur, cost_cur