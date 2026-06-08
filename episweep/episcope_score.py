import os, sys
import csv, itertools
import numpy as np
from functools import reduce

class Filter:
    def __init__(self, struct, contact_f, score_mtx, avr_dist, epitope=None):
        score = Disruption_score(contact_f, score_mtx)  
        self.crd = []
        self.repr_crd = []
        self.docks = set([])
        
        pos_freq = {}
        for design_f in score.contact_list:
            with open(design_f, newline='') as csvfile:  
                mut_list = csv.DictReader(csvfile)
                for mut in mut_list:
                    if mut["index"] == "0":
                        wt_e = float(mut["potential"])
                    else:
                        m = DisruPPI_CSV_line(mut, struct)  
                        if epitope is not None:
                            score.score_cal(m, score.contact_list, epitope)
                        else:
                            score.score_cal(m, score.contact_list, [])

                        if m.inter_dist < avr_dist and m.e < wt_e:
                            self.crd.append(m)
                            self.docks = self.docks.union(set(m.in_contact))
                            for pp in m.pos:
                                if pp not in pos_freq:  
                                    pos_freq[pp] = 1
                                else:
                                    pos_freq[pp] += 1
        for i in self.crd:
            freq_sum = sum([pos_freq[x] for x in i.pos])  
            i.freq = freq_sum
        pos = {}
        for d in self.crd:
            idx = tuple(d.pos)  
            if idx not in pos:  
                pos[idx] = d
            else:
                if pos[idx].episcore < d.episcore:
                    pos[idx] = d
        self.repr_crd = self.crd
    
    def cluster_coverage(self, med, oth, cluster):
        print("\nMaximal docking model coverage calculation.")  
        self.all_cover = []
        self.selected = []
        cluster_set = [[self.repr_crd[x]] for x in med]  
        for n, c in enumerate(oth):
            cluster_set[cluster[n]].append(self.repr_crd[c])
        max_cover = 0
        all_combination = list(itertools.product(*cluster_set))
        print("%d combinations are being considered." % len(all_combination))  
        for comb in all_combination:
            coverage = []
            for cov in [x.in_contact for x in comb]:  
                coverage += cov
            coverage_set = set(coverage)
            if len(coverage_set) > max_cover:
                max_cover = len(coverage_set)
                self.all_cover = [comb]
            elif len(coverage_set) == max_cover:
                self.all_cover.append(comb)
        for i in self.all_cover:
            result = [i, sum(x.episcore for x in i), sum(x.freq for x in i), self.docks]  
            self.selected.append(result)
        self.selected.sort(key=lambda x: x[1], reverse=True)

class DisruPPI_CSV_line:
    def __init__(self, line, crd):
        self.e = float(line["potential"])
        self.design = line["mutation"]
        self.pos = list(map(lambda x: x[1:-1], line["mutation"].split("|")))
        self.crd = list(map(lambda x: crd[x], self.pos))
        self.triangle = reduce(lambda x, y: x+y, self.crd)
        if len(self.crd) == 1:
            self.inter_dist = 0.0001
        elif len(self.crd) == 2:
            self.inter_dist = np.linalg.norm(np.array(self.crd[0])-np.array(self.crd[1]))
        else:
            self.inter_dist = [np.linalg.norm(np.array(self.crd[p1])-np.array(self.crd[p2]))
                               for p1 in range(len(self.crd)-1) for p2 in range(p1+1, len(self.crd))]  
            self.inter_dist = sum(self.inter_dist) / len(self.inter_dist)

class DisruPPI_contact:
    def __init__(self, f):
        base_dir = '/app/targets/episcope'
        filename = '2vxt.00.3_contacts.csv'
        file_path = os.path.join(base_dir, filename)

        with open(file_path, newline='') as csvfile: 
            contact = list(csv.reader(csvfile))
        #with open(f, newline='') as csvfile:  # Use with statement
            #contact = list(csv.reader(csvfile))
        # Convert map objects to lists and dict directly
        self.contact = {x[0]: x[1].split(";") for x in contact}
        self.pos = list(self.contact.keys())
        self.contact_res = list(self.contact.values())

class Disruption_score:
    def __init__(self, contact_f, mtx):
        letter_code = { ... }  
        with open(contact_f, newline='') as csvfile:  
            design_contact_pair = list(csv.reader(csvfile))
        des_f_list = [int(os.path.isfile(x)) for x in (row[0] for row in design_contact_pair)]
        con_f_list = [int(os.path.isfile(x)) for x in (row[1] for row in design_contact_pair)]
        if not (len(des_f_list) == len(design_contact_pair) and len(con_f_list) == len(design_contact_pair)):
            print("Some files do not exist. Please check %s.\n" % contact_f, file=sys.stderr)  
        self.contact_list = {x[0]: DisruPPI_contact(x[1]) for x in design_contact_pair}

        self.mtx = {}
        with open(mtx, newline='') as mtxfile:  
            for l in csv.reader(mtxfile):
                ab = letter_code[l[1]]
                ag = letter_code[l[0]]
                score = float(l[2])
                if ab not in self.mtx:
                    self.mtx[ab] = {ag: score}
                else:
                    self.mtx[ab][ag] = score

    def score_cal(self, d, contact_list, epitope):
        design_pos = list(map(lambda x: x[1:-1], d.design.split("|")))
        d.design_pos = design_pos
        d.pos = "-".join(design_pos)

        d.average_score_list = []
        d.in_contact = []
        for c in contact_list:
            contact = contact_list[c]
            con_res = contact.contact_res  
            con_pos = contact.pos
            overlap = set(con_pos) & set(design_pos)
            if overlap:
                wt_score = 0.0
                mut_score = 0.0
                hit = []
                for res in d.design.split("|"):
                    wt = res[0]
                    mut = res[-1]
                    pos = res[1:-1]
                    if pos in epitope:
                        hit.append(pos)
                    if pos in overlap:
                        for aa in contact.contact[pos]:  
                            wt_score += self.mtx[wt][aa]  
                            mut_score += self.mtx[mut][aa]
                d.hit = len(hit)
                d.hit_pos = "|".join(hit)
                if mut_score - wt_score > 0.0:
                    d.average_score_list.append(mut_score - wt_score)
                    d.in_contact.extend(overlap)  
        if d.average_score_list:  
            d.episcore = np.mean(d.average_score_list)
        else:
            d.episcore = 0  
        d.coverage = len(d.in_contact)