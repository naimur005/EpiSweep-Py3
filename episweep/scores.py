import struct
import sys, os, csv, gzip
from .choice import *
from .TUBEscores import *

class ScorePotential (TUBE):
    def initialize_score_potential(self):
        if self.score_potential is None:
            print ("Score potential file is not found. Exiting the programme.")
            sys.exit()

        self.nrots = []
        self.choices = []
        self.seq_choices = []
        self.onebody_score = {}
        self.twobody_score = {}
        self.filtered_seq_list = []
        self.all_score = []
        self.penalty = 0.0
        self.wt_onebody_score = {}
        self.wt_twobody_score = {}

        print ("Score potential %s is being loaded."%self.score_potential)
    
    def load_tube_score_potential(self):
        self.initialize_tube()
        with open(self.possible_choices, "rb") as choice_csv:
            choice = csv.reader(choice_csv)
            choices = []
            seq_choices = []
            for row in choice:
                p = int(row[0]) - 1
                tube = row[1]
                print (p, tube)
                while len(choices) <= p:
                    choices.append([]); seq_choices.append([])
                tube = row[1]
                pchoices = [tube]
                choices[p].extend(pchoices)
                if tube == self.wt_seq[p]:
                    tube_choice = TubeChoice(tube, list(self.wt_seq[p]), True)
                elif tube.isupper(): # For combinatorial method
                    tube_choice = TubeChoice(tube, list(tube), False)
                else: # For degenerate method
                    tube_choice = TubeChoice(tube, self.all_tubes_dict[tube], False)
                tube_choice.choices = pchoices
                seq_choices[p].append(tube_choice)
            self.choices = choices
            self.seq_choices = seq_choices
            self.choices = self.seq_choices

        with open(self.tube_score_potential, "rb") as score_csv:
            score = csv.reader(score_csv)
            for score_line in score:
                this_line_is_onebody = False
                p1, a1, p2, a2, e = score_line
                pos1 = int(p1) - 1
                try:
                    e = float(e)
                    pos2 = int(p2) - 1
                except:
                    e = float(p2)
                    this_line_is_onebody = True
                if this_line_is_onebody:
                    self.onebody_score[(pos1, a1, 0)] = e
                else:
                    self.twobody_score[(pos1, a1, 0, pos2, a2, 0)] = e
                self.all_score.append(e)
        self.penalty = max(self.all_score)

    def filtered_choice(self):
        res_pos = self.AA_choices.keys()
        res_pos.sort()
        for p in res_pos:
            pos = p-1
            self.choices.append([]); self.seq_choices.append([])
            for aa in self.AA_choices[p]:
                pchoices = [AAChoice_idx(aa, 0)]
                self.choices[pos].extend(pchoices)
                aa_choice = AAChoice(aa, aa==self.wt_seq[pos])
                aa_choice.choices = pchoices
                self.seq_choices[pos].append(aa_choice)
    
    def load_csv_score_potential(self):
        if not self.library_gen:
            self.AA_choices = {}
            with open(self.possible_choices, "r") as choice_csv:
                choice = csv.reader(choice_csv)
                for row in choice:
                    pos = int(row[0])
                    aa = row[1]
                    if pos not in self.AA_choices:
                        self.AA_choices[pos] = [aa]
                    else:
                        self.AA_choices[pos].append(aa)


        with open(self.score_potential, "r") as score_csv:  
            score = csv.reader(score_csv)
            for score_line in score:
                this_line_is_onebody = False
                p1, a1, p2, a2, e = score_line
                pos1 = int(p1) - 1
                try:
                    e = float(e)
                    pos2 = int(p2) - 1 
                except ValueError:
                    # Assuming e = float(p2) was meant to handle lines with only 4 elements where p2 is actually e
                    e = float(p2)  # p2 or a2, a2 may b correct # Corrected to 'a2' assuming 'p2, a2, e' pattern
                    this_line_is_onebody = True
                if this_line_is_onebody:
                    self.onebody_score[(pos1, a1, 0)] = e
                else:
                    self.twobody_score[(pos1, a1, 0, pos2, a2, 0)] = e
                self.all_score.append(e)
        self.filtered_choice()
        self.penalty = max(self.all_score)
    
    def load_binary_score_potential(self):
        OSPREY_AA = ["A","V","L","I","F","Y","W","C","M","S",
                "T","K","R","H","H","H","D","E","N","Q","G","P"]
        rot_nums = {"A": 1, "R": 34, "N": 7, "D": 5, "C": 3,
                "Q": 9, "E": 8, "G": 1, "H": 8, "I": 7,
                "L": 5, "K": 27, "M": 13, "F": 4, "P": 1,
                "S": 18, "T": 18, "W": 7, "Y": 8, "V": 3}
        possible_rot = []
        if self.possible_choices is not None:
            wt_seq_index = dict(enumerate(self.wt_seq))
            with open(self.possible_choices) as choice_csv_file:
                rot_reader = csv.reader(choice_csv_file)
                for row in rot_reader:
                    pos = int(row[0])-1
                    aa = row[1]
                    rot = int(row[2])
                    if pos in wt_seq_index:
                        if wt_seq_index[pos] == aa:
                            del wt_seq_index[pos]
                    row = [pos, aa, rot]
                    possible_rot.append(row)
            if len(wt_seq_index) != 0: # <- WT rotamers were pruned in this case
                sys.stdout.write("Some wild type rotamers were pruned.\n")
                sys.stdout.write("%s\n"%" ".join(map(lambda x: str(x+1)+self.wt_seq[x], wt_seq_index.keys())))
                sys.stdout.write("Rotamers for those positions are being added.\n")
                for pos in wt_seq_index:
                    for r in range(rot_nums[wt_seq_index[pos]]):
                        possible_rot.append([pos, wt_seq_index[pos], r])

        with open(self.score_potential, "rb") as f:
            # npos
            npos = struct.unpack(">I", f.read(4))[0]
            
            # nrots
            for p in range(npos):
                self.nrots.append([])
                naas = struct.unpack(">I", f.read(4))[0]
                for i in range(naas):
                    a = struct.unpack(">b", f.read(1))[0]
                    nrots = struct.unpack(">I", f.read(4))[0]
                    self.nrots[p].append((OSPREY_AA[a],nrots))
                    self.filtered_seq_list.append([p, OSPREY_AA[a]])

            # onebody
            # Construct choices from the score potential

            for p in range(0, npos):
                self.choices.append([]); self.seq_choices.append([])
                for (a,nr) in self.nrots[p]:
                    pchoices = []
                    for r in range(nr):
                        e = struct.unpack(">f", f.read(4))[0]
                        self.all_score.append(e)
                        if self.possible_choices is not None:
                            if [p, a, r] in possible_rot:
                                self.onebody_score[p,a,r] = e
                                pchoices.append(AAChoice_idx(a, r))
                        else:
                            self.onebody_score[p,a,r] = e
                            pchoices.append(AAChoice_idx(a, r))
                    self.choices[p].extend(pchoices)
                    aa_choice = AAChoice(a, a==self.wt_seq[p])
                    aa_choice.choices = pchoices
                    self.seq_choices[p].append(aa_choice)
            # twobody
            self.num_rotamers = {}
            for p in range(0, npos-1):
                for (a,nr) in self.nrots[p]:
                   #if not self.num_rotamers.has_key(a):
                   #    self.num_rotamers[a] = nr
                   #else:
                   #    if self.num_rotamers[a] != nr:
                   #        sys.stderr.write("The number of rotamers at %d must be incorrect. Check the energy matrix converter.\n"%p)
                   #        sys.exit()
                    for r in range(nr):
                        for q in range(p+1, npos):
                            for (b,ns) in self.nrots[q]:
                                for s in range(ns):
                                    e = struct.unpack(">f", f.read(4))[0]
                                    self.twobody_score[p,a,r,q,b,s] = e
                                    self.all_score.append(e)
        self.penalty = max(self.all_score)

    def load_epitope_score_list(self):
        if self.design_type == "disruption":
            self.max_episcore = 0.0
            self.interface_contact_list = dict(map(lambda x: [x[0], x[1].split(";")],csv.reader(open(self.interface_contacts))))
            for p in self.interface_contact_list:
                pos = int(p) - 1
                wt_aa = self.wt_seq[pos]
                in_contact = self.interface_contact_list[p]
                self.max_episcore += sum(map(lambda x: self.disruption_score_matrix[wt_aa][x], in_contact))
        else:
            if self.precomputed_epitope is not None:
                extension = os.path.splitext(self.precomputed_epitope)[1]
                if extension == ".gz":
                    episcore_f = gzip.open(self.precomputed_epitope, "rb")
                elif extension == ".csv":
                    episcore_f = open(self.precomputed_epitope, "rb")
                else:
                    sys.stdout.write("Unknown epitope score file extension: %s\n"%self.precomputed_epitope)
                    sys.stdout.write("Check the design spec file and the key 'precomputed_epitope'.\n")
                    sys.exit()
                epitope_score_list = csv.reader(episcore_f)
                self.max_episcore = 0.0
                self.precomputed_peptide_list = {}
                for i in epitope_score_list:
                    wt_mut = i[0]
                    pos = int(i[1]) - 1
                    peptide = i[2]
                    episcore = float(i[3])
                    if wt_mut == "wt":
                        self.max_episcore += episcore
                    if pos not in self.precomputed_peptide_list:
                        self.precomputed_peptide_list[pos] = [[peptide, episcore]]
                    else:
                        self.precomputed_peptide_list[pos].append([peptide, episcore])
                self.nmer = len(self.precomputed_peptide_list[0][0][0])
            else:
                self.max_episcore = self.profiles.num_hits(self.wt_seq, self.episcore_threshold)