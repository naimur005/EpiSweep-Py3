from .PDBParser import *
from episweep import *
import sys, re
from operator import methodcaller

class rotamer_pruning:
    def __init__(self, sweeper, unpruned_rotamer_file=None):
        if not sweeper.pruning:
            sys.stderr.write("Warning: Rotamer pruning is not going to be performed.\n")
        else:
            self.pruned_rotamers = []
            if not unpruned_rotamer_file:
                self.Goldstein_round = 1
                while sweeper.pruning:
                    self.GoldsteinDEE(sweeper)
                self.steric_clash_check(sweeper)
            else:
                rotamer_remainder = map(methodcaller("split", " "), re.findall(r"\((.*?)\)", open(unpruned_rotamer_file).read()))
                for p in range(len(sweeper.energies.nrots)):
                    for j in sweeper.energies.nrots[p]:
                        a = j[0]
                        AA_in_number = sweeper.energies.AA.index(a)
                        if a == sweeper.wt[p]: # We want to leave at least one wild type rotamer
                            wt_rotamer = 0
                            for r in range(j[1]):
                                current_rotamer = ["%d"%p, "%d"%AA_in_number, "%d"%r]
                                if rotamer_remainder.count(current_rotamer) > 0:
                                    wt_rotamer += 1
                            if wt_rotamer > 0:
                                for r in range(j[1]):
                                    current_rotamer = ["%d"%p, "%d"%AA_in_number, "%d"%r]
                                    if rotamer_remainder.count(current_rotamer) > 0: pass
                                    else: self.pruned_rotamers.append((p, a, r))
                        else:
                            for r in range(j[1]):
                                current_rotamer = ["%d"%p, "%d"%AA_in_number, "%d"%r]
                                if rotamer_remainder.count(current_rotamer) > 0: pass
                                else:
                                    self.pruned_rotamers.append((p, a, r))

            print ("Total %d rotamers were pruned."%len(self.pruned_rotamers))
            for rot in self.pruned_rotamers:
                print ("%d%s%d"%(rot[0], rot[1], rot[2]))
                del sweeper.energies.onebody[rot]
            print()

        accepted_rotamers = {}
        for key in sweeper.energies.onebody.keys():
            pos = key[0]; aa = key[1]; rot = key[2]
            if pos not in accepted_rotamers: 
                accepted_rotamers[pos] = {aa: [rot]}
            else:
                if aa not in accepted_rotamers[pos]: 
                    accepted_rotamers[pos][aa] = [rot]
                else:
                    accepted_rotamers[pos][aa].append(rot)

        sweeper.choices = []
        sweeper.seq_choices = []
        for pos in range(sweeper.n):
            sweeper.choices.append([]); sweeper.seq_choices.append([])
            if pos in accepted_rotamers:
                for aa in accepted_rotamers[pos]:
                    accepted_rotamers[pos][aa].sort()
                    pchoices = [RotChoice(aa,rot) for rot in accepted_rotamers[pos][aa]]
                    sweeper.choices[pos].extend(pchoices)
                    aa_choice = AAChoice(aa, aa==sweeper.wt[pos])
                    aa_choice.choices = pchoices
                    sweeper.seq_choices[pos].append(aa_choice)

        print ("Rotamer choices at each position.")
        for pos, rots in enumerate(sweeper.choices):
            print (pos+sweeper.PDBstart)
            for rot in rots: print(rot)
            print()

    def GoldsteinDEE(self, sweeper):
        def min_pair_E(p, a, r, q, b, s):
            min_diff_Epair = 0
            starting = True
            for c in sweeper.energies.nrots[q]:
                aa = c[0]
                for t in range(c[1]):
                    if self.pruned_rotamers.count((q, aa, t)) == 0:
                        if q > p: # twobody energy matrix has only upper triangle elements
                            #if sweeper.energies.twobody.has_key((p, a, r, q, aa, t)) and sweeper.energies.twobody.has_key((p, b, s, q, aa, t)):
                            Epair_diff = sweeper.energies.twobody[(p, a, r, q, aa, t)] - sweeper.energies.twobody[(p, b, s, q, aa, t)]
                        else:
                            #if sweeper.energies.twobody.has_key((q, aa, t, p, a, r)) and sweeper.energies.twobody.has_key((q, aa, t, p, b, s)):
                            Epair_diff = sweeper.energies.twobody[(q, aa, t, p, a, r)] - sweeper.energies.twobody[(q, aa, t, p, b, s)]
                        if starting:
                            min_diff_Epair = Epair_diff
                            starting = False
                        elif Epair_diff < min_diff_Epair:
                            min_diff_Epair = Epair_diff
            return min_diff_Epair

        ########### Goldstein DEE main function  ##############
        print ("Simple Goldstein DEE (Round %d)."%self.Goldstein_round)
        self.num_Goldstein_pruned = 0
        for p in range(len(sweeper.energies.nrots)):
            pruned_here = False
            for j in sweeper.energies.nrots[p]:
                a = j[0]
                if a != sweeper.wt[p]: # don't prune wild type rotamers
                    for r in range(j[1]):
                        if self.pruned_rotamers.count((p, a, r)) < 1:
                            for k in sweeper.energies.nrots[p]:
                                b = k[0]
                                for s in range(k[1]):
                                    if (p, a, r) in sweeper.energies.onebody and (p, b, s) in sweeper.energies.onebody:
                                        if self.pruned_rotamers.count((p, b, s)) < 1 and self.pruned_rotamers.count((p, a, r)) < 1:
                                            if (a != b) or (a == b and r != s):
                                                Energy_r = sweeper.energies.onebody[(p, a, r)]
                                                Energy_s = sweeper.energies.onebody[(p, b, s)]
                                                PairE = 0.0
                                                for q in range(len(sweeper.energies.nrots)):
                                                    if p!=q and sweeper.cm[p][q] == 1:
                                                        PairE += min_pair_E(p, a, r, q, b, s)
                                                TotalE = Energy_r - Energy_s + PairE
                                                if TotalE > sweeper.deltaE:
                                                    self.pruned_rotamers.append((p, a, r)) 
                                                    self.num_Goldstein_pruned += 1
                                                    if TotalE > 1.0e+6:
                                                        print ("%d%s%d pruned (%.3e)"%(p+1, a, r, TotalE))
                                                    else:
                                                        print ("%d%s%d pruned (%.3f)"%(p+1, a, r, TotalE))
                                                    pruned_here = True
        print ("%d rotamers were pruned at this stage."%self.num_Goldstein_pruned)
        if self.num_Goldstein_pruned > 0:
            self.Goldstein_round += 1
            sweeper.pruning = True
        else:
            print ("Simple Goldstein DEE is done.")
            sweeper.pruning = False

    def steric_clash_check(self, sweeper):
        print ("Checking steric clashes.")
        self.num_steric_clashes = 0
        for p in range(len(sweeper.energies.nrots)):
            for j in sweeper.energies.nrots[p]:
                a = j[0]
                if a == sweeper.wt[p]: pass # don't prune wild type rotamers
                else:
                    for r in range(j[1]):
                        if (p, a, r) in sweeper.energies.onebody: 
                            selfE = sweeper.energies.onebody[(p, a, r)]
                            if selfE >= sweeper.stericE and self.pruned_rotamers.count((p, a, r)) < 1:
                                self.num_steric_clashes += 1
                                self.pruned_rotamers.append((p, a, r))
                                if selfE > 1.0e+6:
                                    print ("%d%s%d pruned (%.3e)"%(p+1, a, r, selfE))
                                else:
                                    print ("%d%s%d pruned (%.3f)"%(p+1, a, r, selfE))
        print ("%d rotamers were pruned at this stage."%self.num_steric_clashes)