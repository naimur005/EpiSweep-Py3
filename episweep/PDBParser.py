import sys, math, csv, os

class PDBLine:
    def __init__(self, line):
        self.atom_name = line[12:16].strip()
        self.res_name = line[17:20].strip()
        self.resnum = int(line[22:26])
        self.chain = line[21]
        X = float(line[30:38])
        Y = float(line[38:46])
        Z = float(line[46:54])
        self.coord = [X, Y, Z]
        self.line = line

class PDB:
    def structure_pair_calculator(self, cutoff=8):
        def atom_distance(a, b):
            distance = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
            return distance
        self.possible_pairs = []
        if self.structure is None:
            sys.stderr.write("Structure file was not provided.\n")
            sys.exit()
        else:
            pdb_structure = open(self.structure).readlines()
            coord = filter(lambda x: x[:4] == "ATOM", pdb_structure)
            coord = map(lambda x: PDBLine(x), coord)
            CB_coord = []
            for p in coord:
                if p.res_name == "GLY":
                    if p.atom_name == "CA":
                        CB_coord.append(p.coord)
                else:
                    if p.atom_name == "CB":
                        CB_coord.append(p.coord)
            for i in range(len(CB_coord)-1):
                for j in range(i+1, len(CB_coord)):
                    d = atom_distance(CB_coord[i], CB_coord[j])
                    if d < cutoff:
                        self.possible_pairs.append((i, j))
        
        ########################
        # Save pair information to a file
        if self.pairs is not None:
            with open(self.pairs, "w", newline='') as pair_file:  
                pair_file_line = csv.writer(pair_file)
                pair_lines = map(lambda p: [p[0]+1, p[1]+1], self.possible_pairs)
                # Pair starts from 0. This make residues start at 1
                for pair in pair_lines:
                    pair_file_line.writerow(pair)

    def OSPREY_config_maker(self, target=None):
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
"""%(self.structure, len(self.wt_seq), len(self.wt_seq), " ".join(map(str, map(lambda x: x+1, range(len(self.wt_seq))))))
        
        respos = sorted(self.AA_choices.keys())
        resallowed = []
        for p, pos in enumerate(respos):
            allowed_AA = " ".join(map(lambda x: AA1to3[x].lower(), self.AA_choices[pos]))
            resallowed.append("resAllowed0_%d %s"%(p, allowed_AA))
        OSPREY_DEE_content = DEE_cfg + "\n".join(resallowed) + "\n"

        DEE_output_file = target + ".DEE.cfg"
        system_output_file = target + ".System.cfg"
        
        ########################
        # Save as a OSPREY DEE configuration file
        f = open(DEE_output_file, "w")
        f.write(OSPREY_DEE_content)
        f.close()
        f = open(system_output_file, "w")
        f.write(system_cfg)
        f.close()

        sys.stdout.write("OSPREY DEE and System (Protein) configuration files are saved as:\n")
        sys.stdout.write("\n%s\n%s\n\n"%(DEE_output_file, system_output_file))

    def CLEVER_design_file_maker(self, target=None, output=None):
        self.CLEVER_flag=True ### <-- This information is passed to TUBEscores for CLEVER special treatment
        respos = self.AA_choices.keys()
        respos.sort()
        clever_content = ["#design_start"]
        for pos in respos:
            clever_choices = []
            clever_choices = self.AA_choices[pos].remove(str(self.wt_seq[pos-1]))
            clever_choices = [self.wt_seq[pos-1]]+self.AA_choices[pos]
            #AAs = " ".join(self.AA_choices[pos])
            AAs = " ".join(clever_choices)
            clever_content.append("%d %s"%(pos, AAs))
        clever_content += ["#design_end", "", "#cluster_start"]
        clever_content += map(str, respos)
        pair = map(lambda x: "%d %d"%(x[0]+1, x[1]+1), self.possible_pairs)
        clever_content += pair
        clever_content += ["#cluster_end"]
        if output is not None:
            f = open(output, "w")
            f.write("\n".join(clever_content))
            f.close()

        sys.stdout.write("CLEVER design file generated successfully.\n")
