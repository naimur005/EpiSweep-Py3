#!/usr/bin/python3

import os, sys, csv, string
import argparse, itertools
global NMER
NMER = 9 # nonamer default

class Germline:
    def __init__(self, f):
        self.seq = list(csv.DictReader(open(f)))
        num = self.seq[0].keys()
        ref_num = []
        for n in num:
            try:
                int(n)
                ref_num.append(n+"-")
            except:
                ref_num.append(n)
        ref_num.sort(key=lambda x: int(x[:-1]))
        last_num = ref_num[-1]
        if last_num[-1] == "-":
            self.last_num = last_num[:-1]
        else:
            self.last_num = str(int(last_num[:-1])+1)

class Design:
    def __init__(self, line, L, H):
        self.L = L
        self.H = H
        self.index = line["index"]
        self.curve = line["curve"]
        self.HSC = line["episcore"]
        self.e = line["potential"]
        self.mutation = line["mutation"]
        self.seq = line["sequence"]
        self.seqL = self.seq[:len(L.header_number)]
        self.seqH = self.seq[len(L.header_number):]
        self.numbered_seqL = dict(map(lambda x: [L.header_number[x[0]], x[1]], enumerate(self.seqL)))
        self.numbered_seqH = dict(map(lambda x: [H.header_number[x[0]], x[1]], enumerate(self.seqH)))
        self.LN9, Lcount, self.raw_HSC_L = self.score_cal(self.seqL, self.L)
        self.HN9, Hcount, self.raw_HSC_H = self.score_cal(self.seqH, self.H)
        self.HSC_L = 100 * self.raw_HSC_L / (NMER * Lcount)
        self.HSC_H = 100 * self.raw_HSC_H / (NMER * Hcount)
        self.HSC = 100 * (self.raw_HSC_L + self.raw_HSC_H) / (NMER * (Lcount + Hcount))
        self.print_line = [self.index, self.curve, self.HSC_L, self.HSC_H, self.HSC, self.LN9, self.HN9, self.LN9+self.HN9, self.e, self.mutation, self.seq]

    def score_cal(self, seq, chain):
        N9 = 0
        count = 0
        raw = 0.0
        for i in range(len(seq)-NMER + 1):
            pep = seq[i:i+NMER]
            score = HSC_score(pep, chain.header_number[i:i+NMER], chain.germ)
            if score == NMER:
                N9 += 1
            count += 1
            raw += score
        return N9, count, raw

class AbInput:
    def __init__(self, f):
        header = list(csv.reader(open(f)))[0]
        self.header_number = self.number(header)
        self.numbered_seq = list(csv.DictReader(open(f)))[0]
        self.seq = "".join(map(lambda x: self.numbered_seq[x], self.header_number))
        self.nmer_index = []
        self.nmer_pep = []
        self.germ = None
        for n in range(len(self.header_number) - NMER + 1):
            indices = self.header_number[n:n+NMER]
            pep = map(lambda x: self.numbered_seq[x], indices)
            self.nmer_index.append(indices)
            self.nmer_pep.append(pep)

    def number(self, header):
        num_list = []
        for n in header:
            try:
                int(n[0]) # select only numbered items from the CSV input
                num_list.append(n)
            except: pass
        return num_list
    
    def mutation_choice(self, num_seq, cdr_def):
        iden_seq = []
        for s in self.germ.seq:
            identity = 0
            for p in self.header_number:
                if p in s: 
                    if self.numbered_seq[p] == s[p]:
                        identity += 1
            iden_seq.append([identity, s])
        iden_seq.sort(reverse=True)
        if num_seq:
            iden_seq = iden_seq[:num_seq]
        self.choice = []
        for p in self.header_number:
            if self.numbered_seq[p] == "-":
                pass
            else:
                pos_choice = [self.numbered_seq[p]]
                if p not in cdr_def:
                    for s in map(lambda x: x[1], iden_seq):
                        if s.has_key(p):
                            if pos_choice.count(s[p]) == 0 and s[p] != "-":
                                pos_choice.append(s[p])
                self.choice.append(pos_choice)

    def precomputation(self):
        self.precomputed = []
        for p in range(len(self.choice)-NMER + 1):
            wt_pep = "".join(map(lambda x: self.numbered_seq[x], self.header_number[p:p+NMER]))
            for pep in itertools.product(*self.choice[p:p+NMER]):
                pep_str = "".join(list(pep))
                if wt_pep == pep_str:
                    pep_status = "wt"
                else:
                    pep_status = "var"
                self.precomputed.append([pep_status, pep_str, HSC_score(pep, self.header_number[p:p+NMER], self.germ)])

def total_HSC(seq, chain):
    N9 = 0
    count = 0
    raw = 0.0
    for i in range(len(seq)-NMER + 1):
        pep = seq[i:i+NMER]
        score = HSC_score(pep, chain.header_number[i:i+NMER], chain.germ)
        if score == NMER:
            N9 += 1
        count += 1
        raw += score
    return raw, count, float(raw) / (NMER*count), N9

def HSC_score(pep, idx, germ):
    score_list = []
    for g in germ.seq:
        score = 0
        for n, i in enumerate(idx):
            if i not in g: 
                aa = "-"
            else: aa = g[i]
            if pep[n] == aa:
                score += 1
        score_list.append(score)
    return max(score_list)

def flush_HSC(precomp1, precomp2, output):
    pep1 = precomp1[-1][1]
    pep2 = precomp2[0][1]
    for i in range(NMER-1):
        connect_pep = pep1[i+1:] + pep2[:i+1]
        precomp1.append(["wt", connect_pep, 0])
    output_f = csv.writer(open(output, "wb"))
    counter = 0
    for p in precomp1+precomp2:
        if p[0] == "wt":
            counter += 1
        p.insert(1, counter)
        output_f.writerow(p)

def flush_DEE(choice, pdb):
    target = os.path.splitext(os.path.basename(pdb))[0]
    AA1to3 = {"A": "ALA",
            "E": "GLU",
            "Q": "GLN", 
            "D": "ASP", 
            "N": "ASN", 
            "L": "LEU", 
            "G": "GLY", 
            "K": "LYS", 
            "S": "SER", 
            "V": "VAL", 
            "R": "ARG", 
            "T": "THR", 
            "P": "PRO", 
            "I": "ILE", 
            "M": "MET", 
            "F": "PHE", 
            "Y": "TYR", 
            "C": "CYS", 
            "W": "TRP", 
            "H": "HIP"} # <-- It's temporarily written as HIP
    
    DEE_cfg = """runName %s
addWT true
addWTRots true
onlysinglestrand 0
numMaxMut 2
algOption 1
doMinimize false
doDACS false
minRatiodiff 0.15
initDepth 10
subDepth 1
diffFact 6
splitFlags false
distrDACS false
minimizeBB false
doBackrubs false
backrubFile none
minEnergyMatrixName %s
maxEnergyMatrixName none
useEref true
initEw 0
dielectConst 8.0
doSolvationE true
solvScale 0.5
pruningE 100.0
stericE 30.0
approxMinGMEC false
lambda 10000.0
preprocPairs true
pairSt 100.0
scaleInt false
maxIntScale 0.1
genInteractionGraph false
distCutoff 30.0
eInteractionCutoff 0.0
EREFMATRIXNAME %s.eref
outputConfInfo c_%s
outputPruneInfo p_%s
"""%tuple([target]*5)

    system_cfg = """pdbName %s
numOfStrands 1
strand0 1 %d
strandMutNums %d 0 
strandMut0 %s
strandAA0 true
strandAA1 false
strandRotTrans0 false
strandRotTrans1 true
numCofRes 0
cofMap0 0
"""%(pdb, len(choice), len(choice), " ".join(map(str, map(lambda x: x+1, range(len(choice))))))
        
    print ("Preparing OSPREY configuration files...")
    DEE_f = "%s_DEE.cfg"%target
    system_f = "%s_System.cfg"%target
    DEE = open(DEE_f, "w")
    system = open(system_f, "w")
    DEE.write(DEE_cfg)
    system.write(system_cfg)
    for n, c in enumerate(choice):
        DEE.write("resAllowed0_%d %s\n"%(n, " ".join(map(lambda x: AA1to3[x].lower(), c))))
    print ("OSPREY System file: %s")%DEE_f
    print ("OSPREY DEE file: %s")%system_f
    print

if __name__== "__main__":
    parser = argparse.ArgumentParser(description='HSC calculator')
    parser.add_argument('-H', '--germline_heavy', dest="H", action="store", default=None, help='Heavy chain sequence (ANACRI annotated)')
    parser.add_argument('-L', '--germline_light', dest="L", action="store", default=None, help='Light chain germline (ANACRI annotated)')
    parser.add_argument('-D', '--design_file', dest="design_file", action="store", default=None, help='(Optional) Design output file to recalculate the HSC score (Standard output)')
    parser.add_argument('-n', '--number_sequence', dest="number_sequence", action="store", type=int, default=False, help='Number of best matched germline sequences to include (Default: All sequences)')
    parser.add_argument('-o', '--output_HSC', dest="output_HSC", action="store", default=None, help='(Optional) Precomputed HSC score output (CSV format)')
    parser.add_argument('-p', '--pdb', dest="pdb", action="store", default=None, help='(Optional) PDB file (OSPREY configuration files will be created)')
    parser.add_argument('-s', '--numbering', dest="numbering", action="store", default="Kabat", help='(Optional) Numbering system (Default: Kabat - or Chothia)')
    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit()
    if args.numbering == "Chothia":
        numbering = "Chothia"
        CDR_L = map(lambda x: ("%d%s"%(x[0], x[1])).strip(), list(itertools.product(*[range(24, 34+1) + range(50, 56+1) + range(89, 97+1), list(" "+string.uppercase)])))
        CDR_H = map(lambda x: ("%d%s"%(x[0], x[1])).strip(), list(itertools.product(*[range(26, 32+1) + range(52, 56+1) + range(95, 102+1), list(" "+string.uppercase)])))
    else:
           numbering = "Kabat"
        #CDR_L = [f"{i}{j}".strip() for i, j in itertools.product(range(24, 35) + range(50, 57) + range(89, 98), " " + ascii_uppercase)]
        # Combining ranges by converting them to lists first and using itertools.chain for more efficiency
    combined_range = list(itertools.chain(range(24, 35), range(50, 57), range(89, 98)))

        # Using string.ascii_uppercase for uppercase letters and a space
    CDR_L = [f"{i}{j}".strip() for i, j in itertools.product(combined_range, " " + string.ascii_uppercase)]
        #CDR_H = [f"{i}{j}".strip() for i, j in itertools.product(range(31, 35) + range(50, 66) + range(95, 103), " " + ascii_uppercase)] + ["35", "35A", "35B"]
        # Combining ranges by converting them to lists first and using itertools.chain for efficiency
    combined_range = list(itertools.chain(range(31, 35), range(50, 66), range(95, 103)))

    CDR_H = [f"{i}{j}".strip() for i, j in itertools.product(combined_range, " " + string.ascii_uppercase)] + ["35", "35A", "35B"]

    HV_data = os.path.join(os.environ["EPISWEEP"], "data", "IGHV_human_germline_KABAT.csv")
    LV_data = os.path.join(os.environ["EPISWEEP"], "data", "IGKV_human_germline_KABAT.csv")
    HV_germline = Germline(HV_data)
    LV_germline = Germline(LV_data)
    H = AbInput(args.H)
    L = AbInput(args.L)
    H.germ = HV_germline
    L.germ = LV_germline
    raw_H, count_H, HSC_H, HN9 = total_HSC(H.seq, H)
    raw_L, count_L, HSC_L, LN9 = total_HSC(L.seq, L)
    H.mutation_choice(args.number_sequence, CDR_H)
    L.mutation_choice(args.number_sequence, CDR_L)

    if args.pdb is not None:
        flush_DEE(L.choice+H.choice, args.pdb)

    if args.output_HSC is not None:
        print ("Precomputing the HSC score...")
        if not args.number_sequence:
            print ("All germline sequences available are going to be used for mutational choices")
        else:
            print ("Only %d best matched germline sequences are going to be used for mutational choices")%args.number_sequence
        H.precomputation()
        L.precomputation()
        flush_HSC(L.precomputed, H.precomputed, args.output_HSC)
        print ("Precomputation of the HSC score is done: %s")%args.output_HSC
        print()

    if args.design_file is not None:
        designs = csv.DictReader(open(args.design_file))
        new_f = os.path.splitext(args.design_file)[0] + "_HSC.csv"
        new_design = csv.writer(open(new_f, "wb"))
        header = ["index", "curve", 
                "HSC_L", "HSC_H", "HSC_total", 
                "N9_L", "N9_H", "N9_total", 
                "potential", "mutation", "sequence"]
        new_design.writerow(header)
        for design in designs:
            d = Design(design, L, H)
            new_design.writerow(d.print_line)
        print ("HSC calculated design output file: %s")%new_f
    else:
        print ("The numbering scheme is %s")%numbering
        print ()
        print ("Heavy chain HSC score: %.3f (N9: %d)")%(100*HSC_H, HN9)
        print ("Light chain HSC score: %.3f (N9: %d)")%(100*HSC_L, LN9)
        print ("Total HSC score: %.3f (N9: %d)")%(100*float(raw_H+raw_L) / (NMER * (count_H + count_L)), LN9+HN9)