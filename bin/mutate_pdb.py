#!/usr/bin/python3

import re,sys,optparse,csv
import copy
from episweep import *

class mutation:
    def __init__(self, pdb):
        self.AA1to3 = {"A": "ALA", "E": "GLU", "Q": "GLN", "D": "ASP", "N": "ASN", "L": "LEU", "G": "GLY", "K": "LYS", "S": "SER", "V": "VAL", "R": "ARG", "T": "THR", "P": "PRO", "I": "ILE", "M": "MET", "F": "PHE", "Y": "TYR", "C": "CYS", "W": "TRP", "H": "HIS"}
        self.readPDB(pdb)

    def readPDB(self, structure):
        f = open(structure).readlines()
        atom_lines = filter(lambda x: x[:4] == "ATOM", f)
        atom_lines = map(lambda x: PDBLine(x), atom_lines)
        self.hetatm = filter(lambda x: x[:6] == "HETATM", f)
        self.coord = {}
        self.resnums = []
        for l in atom_lines:
            if not self.coord.has_key(l.resnum):
                self.coord[l.resnum] = [l.line]
                self.resnums.append(l.resnum)
            else:
                self.coord[l.resnum].append(l.line)

    def mutate(self, wt_seq, seq_list):
        bkbne = ["N", "CA", "C", "O"]
        for str_num, seq in enumerate(seq_list):
            f = "%s.pdb"%(str_num+1)
            mutated_file = open(f, "w")
            mut_list = []
            new_coord = copy.deepcopy(self.coord)
            for p, aa in enumerate(wt_seq):
                pos = p+1
                if aa != seq[p]:
                    coord_line = new_coord[pos][:]
                    new_coord[pos] = []
                    mut = self.AA1to3[seq[p]]
                    mut_list.append("%s%d%s"%(aa, pos, seq[p]))
                    for atom_line in coord_line:
                        if atom_line[13:15].strip() in bkbne:
                            new_line = atom_line[:17] + mut + atom_line[20:]
                            new_coord[pos].append(new_line)
            for res in self.resnums:
                mutated_file.write("".join(new_coord[res]))
            sys.stdout.write("Mutations list: %s\n"%(" ".join(mut_list)))
            sys.stdout.write("This variant is saved in %s\n\n"%f)

parser = optparse.OptionParser(description='Epitope predictor')
parser.add_option('-i', '--design_spec',
        dest="design_spec",
        action="store",
        help='Input: Design spec file')
parser.add_option('-s', '--sequence_file',
        dest="sequence_file",
        action="store",
        help='Input: Clever sequence file')

args, opts = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit()

design = Design()
design.load_spec(args.design_spec)

mutated_sequences = map(lambda x: "".join(x.split()[1:]), open(args.sequence_file).readlines())

wt_seq = "".join(map(lambda x: x.strip("\n"), open(design.target).readlines()[1:]))
a = mutation(design.structure)
a.mutate(wt_seq, mutated_sequences)