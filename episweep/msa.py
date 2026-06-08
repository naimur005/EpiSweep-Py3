from collections import defaultdict
import math, sys, os, csv
from .epitopes import *

class MSA:
    def initialize_filter(self):
        if self.target is None:
            sys.stderr.write("Target sequence file is not given.\n")
            sys.exit()
        else:
            self.wt_seq = "".join(map(lambda x: x.strip("\n"), open(self.target).readlines()[1:]))
            self.n = len(self.wt_seq)
            if self.episcore_mtx is not None:
                self.profiles = PPProfiles(self.episcore_mtx)
            if self.msa is None:
                sys.stderr.write("MSA file is not given.\n")
                sys.exit()
            else:
                self.readMSA()

    def readMSA(self):
        print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print ('  EpiSweep Design Specification  ')
        print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

        self.AAs = "ACDEFGHIKLMNPQRSTVWY"
        self.penalty = 0.0
        self.seq_list=[]
        self.name_list=[]
        self.pos_aa_frequency = []

        self.mut_lock = []
        self.secondary_structure = None
        
        self.seq_weight_list = []
        self.onebody_score = {}
        self.twobody_score = {}

        self.coupling_count = defaultdict(int)
        
        self.AA_choices = {}  # contains all possible mutation choices at each position
        self.possible_pairs = [] # contains pair information

        self.choices = []
        self.seq_choices = [] # Those variables are fed into EpiSweep
        
        ################################
        # Reading MSA file 
        if not os.path.isfile(self.msa):
            sys.stderr.write("Error: %s does not exist,\n"%self.msa)
            sys.exit(0)

        in_file=open(self.msa,'r')
        data=''
        for line in in_file:
            line=line.strip()
            for c in line:
                if c=='>':
                    self.name_list.append(line)
                    if data:
                        self.seq_list.append(data)
                        data=''
                    break
                else:
                    line=line.upper()
            if all([k==k.upper() for k in line]):
                data=data+line
    
        ################################
        # Sometimes, sequence string contains header as well; so we remove it here
        for i in range(len(self.seq_list)):
            self.seq_list[i] = self.seq_list[i].replace(self.name_list[i], "")
            self.seq_weight_list.append(1) # This is a default sequence weight. If sequence_weight is called, it will be regenerated.
        self.filtered_seq_list = self.seq_list
        self.filtered_name_list = self.name_list
        self.penalty = (-1.0) * math.log(0.5/len(self.seq_list))
        self.remove_wt_gaps()

    ################################################################ SEQUENCE GAP FILTER ###
    #
    # This is MSA pre-filtering that removes columns from MSA if:
    #   1. Wildtype residue has a gap.
    #   2. Wildtype residue is a non aminoacid character.
    # Penalty is also calculated here based upon MSA size.
    def remove_wt_gaps(self):
        if self.wt_seq is None:
            self.wt_seq = self.seq_list[0]
            sys.stderr.write("Target sequence is not specified.\nThe first sequence in the MSA file will be used as a wild type sequence.\n")
            sys.stdout.write("Target wild type sequence:\n%s\n"%self.wt_seq)
            # if wild type sequence is not given, take the first one as the wild type sequence

        gap_pos = [p for p, aa in enumerate(self.wt_seq) if not aa in self.AAs]
        # The gapped positions are removed in all the sequences in the list
        if bool(gap_pos):
            for n, seq in enumerate(self.seq_list):
                seq_replace = list(seq)
                for p in gap_pos:
                    seq_replace[p] = ""
                self.seq_list[n] = "".join(seq_replace)

    ############################################### ONEBODY SCORE CALCULATION ###
    #
    # Weighted frequency calculation for each AA at each position
    # A sequence weight file is also required. If not provided, equal weights are
    #   assigned to all sequences (i.e. 1)
    # The frequencies are further converted to logs.
    def sequence_weight(self, seq_weights_file=None):
        print ('Sequence length is %d' % (len(self.wt_seq)))
        if seq_weights_file is not None:
            if os.path.isfile(seq_weights_file):
                self.seq_weight_list = []
                for weight in map(float, open(seq_weights_file,'r')):
                    self.seq_weight_list.append(float(weight))
            else:
                sys.stderr.write("Sequence weight file \"%s\" does not exist.\n"%seq_weights_file)
                sys.exit(0)
        else:
            sys.stdout.write('Sequence weight file not provided.\n')
            sys.stdout.write('Example file: Single column file with integer weights for each sequence in MSA).\n')
            sys.stdout.write('Assigning equal weights to all sequences...\n')

    ########################################################### IDENTITY THRESHOLD FILTER ###
    # Filter sequences here based upon
    #   1. Identity thresholds and
    #   2. Fraction of gaps allowed
    def msa_filter(self, low_thresh=0.35, high_thresh=0.90, fraction_gap_allowed=0.25):

        ################################
        # Wildtype sequence should be the first in MSA and is always included in the filtered list
        self.filtered_seq_list = []
        self.filtered_seq_list.append(self.wt_seq)
        wt_seq_length = len(self.wt_seq)
        print ("Remaining sequences:")
        for n, seq in enumerate(self.seq_list):
            gap_count = len(list(filter(lambda aa: aa not in self.AAs, list(seq))))
            if gap_count / wt_seq_length < fraction_gap_allowed:
                # not considered if there are too many gaps
                enu_seq = set(enumerate(seq))
                add_to_filtered = True
                for filtered_seq in self.filtered_seq_list:
                    enu_filtered_seq = set(enumerate(filtered_seq))
                    seq_id = len(enu_seq & enu_filtered_seq) / float(wt_seq_length)
                    if seq_id < low_thresh or seq_id > high_thresh:
                        add_to_filtered = False
                        break
                if add_to_filtered:
                    print (self.name_list[n].lstrip(">"))
                    self.filtered_seq_list.append(seq)

        print()
        print ('%d sequences selected from %d MSA input sequences' % (len(self.filtered_seq_list), len(self.seq_list)))
        print()

        self.onebody_score_calculator()

    def onebody_score_calculator(self):
        ################################
        # Assigning user provided weights at each position
        col_weight = []
        for i in range(len(self.wt_seq)):
            column_gaps = 0.0
            aa_dict = dict([(aa, 0.0) for aa in self.AAs])
            col_weight.append(0.0)
            for j in range(len(self.filtered_seq_list)):
                try:
                    aa_dict[self.filtered_seq_list[j][i]] = aa_dict[self.filtered_seq_list[j][i]] + self.seq_weight_list[j]
                    col_weight[-1] += self.seq_weight_list[j]
                except:
                    column_gaps += 1
                    col_weight[-1] += self.seq_weight_list[j]
            self.pos_aa_frequency.append(aa_dict)

        ################################
        # Calculate conservation log scores from AA frequecies
        max_conservation_score = 0.0
        for idx, val in enumerate(self.pos_aa_frequency):
            pos = idx + 1
            for aa in self.AAs:
                if val[aa] > 0.0 and aa != '-':
                    if pos not in self.AA_choices: 
                        self.AA_choices[pos] = [aa]
                    else: self.AA_choices[pos].append(aa)
                    conservation_log_score = (-1.0) * math.log(val[aa]/col_weight[idx])
                    self.onebody_score[idx,aa,0] = conservation_log_score
                    if conservation_log_score > max_conservation_score: max_conservation_score = conservation_log_score

        print ('AA conservation log scores calculated using weights ...')
        print ('max_conservation_score', max_conservation_score)
        print ('max penalty', self.penalty)
        print ()

    ################################################################ ONEBODY SCORE CALCULATION ###
    #
    # Weighted frequency calculation for AA pairs
    # Sequence weights are used from self.seq_weight_list variable.
    # Frequency for coupled 'a' and 'b' at positions 'i' and 'j' is stored as i, a, j, b
    def seq_potential_constructor(self):
        ################################
        # Calculating coupled frequencies
        pair_weight = {}
        pair_col_weight = {}

        for k in range(len(self.filtered_seq_list)):
            for i in range(len(self.filtered_seq_list[k])):
                a = self.filtered_seq_list[k][i]
                for j in range(i+1, len(self.filtered_seq_list[k])):
                    b = self.filtered_seq_list[k][j]
                    pair_col_weight_key = (i, j)
                    pair_weight_key = (i, a, j, b)
                    coupling_key = "%d%s%d%s"%(i, a, j, b)
                    if a != '-' and b != '-':
                        if pair_weight_key not in pair_weight: 
                            pair_weight[pair_weight_key] = self.seq_weight_list[k]
                            self.coupling_count[coupling_key] = 1
                        else:
                            pair_weight[pair_weight_key] += self.seq_weight_list[k]
                            self.coupling_count[coupling_key] += 1
                    if pair_col_weight_key not in pair_col_weight:
                        pair_col_weight[pair_col_weight_key] = self.seq_weight_list[k]
                    else:
                        pair_col_weight[pair_col_weight_key] += self.seq_weight_list[k]
        
        ################################
        # Calculating log scores from coupled frequencies
        max_coupled_score = 0.0
        for k in range(len(self.filtered_seq_list)):
            for i in range(len(self.filtered_seq_list[k])):
                a = self.filtered_seq_list[k][i]
                for j in range(i+1, len(self.filtered_seq_list[k])):
                    b = self.filtered_seq_list[k][j]
                    if a != '-' and b != '-':
                        couple_score = (-1.0) * math.log(float(pair_weight[(i,a,j,b)])/pair_col_weight[(i,j)])
                        if (i, a, 0) not in self.onebody_score: i_aa = self.penalty
                        else: i_aa = self.onebody_score[(i, a, 0)]
                        if (j, b, 0) not in self.onebody_score: j_bb = self.penalty
                        else: j_bb = self.onebody_score[(j, b, 0)]
                        
                        cov = couple_score - i_aa - j_bb
                        self.twobody_score[i,a,0,j,b,0] = cov
                        if self.twobody_score[i,a,0,j,b,0] > max_coupled_score: max_coupled_score = self.twobody_score[i,a,0,j,b,0]
    
        print ('max_coupled_score', max_coupled_score)
        print ('Coupled log scores step2 done ...')
        print ('Coupled log scores calculated for each pair ...')
        print ()

    def generate_seq_potential(self):
        self.CLEVER_flag = False
        res_pos = self.AA_choices.keys()
        res_pos.sort()
        if self.possible_choices is not None:
            with open(self.possible_choices, "wb") as choice_file:
                choice_file_writer = csv.writer(choice_file)
                for pos in res_pos:
                    for aa in self.AA_choices[pos]:
                        choice_file_writer.writerow([pos, aa])
            sys.stdout.write("\nAmino acid choices are is saved in \"%s\".\n\n"%self.possible_choices)
        if self.score_potential is not None:
            with open(self.score_potential, "wb") as score_potential_file:
                score_potential_file_writer = csv.writer(score_potential_file)
                # score potential starts from 0. This make residues start at 1
                for pos in res_pos:
                    p = pos - 1
                    for aa in self.AA_choices[pos]:
                        onebody_key = (p, aa, 0)
                        if onebody_key in self.onebody_score:
                            score_key = self.onebody_score[onebody_key]
                        else:
                            score_key = self.penalty
                        line = [pos, aa, score_key, '', '']
                        score_potential_file_writer.writerow(line)
                for p1 in range(len(res_pos)-1):
                    pos1 = res_pos[p1]
                    for aa1 in self.AA_choices[pos1]:
                        for p2 in range(p1+1, len(res_pos)):
                            pos2 = res_pos[p2]
                            for aa2 in self.AA_choices[pos2]:
                                twobody_key = (p1, aa1, 0 , p2, aa2, 0)
                                if twobody_key in self.twobody_score: 
                                    score_key = self.twobody_score[twobody_key]
                                else:
                                    score_key = self.penalty
                                line = [pos1, aa1, pos2, aa2, score_key]
                                score_potential_file_writer.writerow(line)
            sys.stdout.write("Score potential is saved in \"%s\".\n\n"%self.score_potential)

    ####################################################################### AA CHOICE FILTER ###
    #
    # Filter choices based upon expected AA frequencies:
    #   1. Pass a 2-column file with AA column and their frequency column.
    #     eg. McCaldon AA threshold
    #   2. A higher frequency value makes AA selection more stringent.
    #   3. If the file is missing, a low frequency cutoff is used, i.e., relaxed filtering. 
    def background_frequency_filter(self, input=None):
        sys.stdout.write("\nAmino acid background frequency filter removes mutations to each site.\n\n")
        background_frequency_file = os.path.join(self.episweep_data, input)
        ################################
        # Reading McCaldon expected AA frequencies (custom conservation threshold file)
        background_frequency_dict = {}
        if background_frequency_file is not None:
            if os.path.isfile(background_frequency_file):
                print ('Reading AA background frequency threshold data from AA background frequency file ...')
                back_freq_f = csv.reader(open(background_frequency_file, 'r'))
                #next(back_freq_f, None) # skip the header
                for line in back_freq_f:
                    background_frequency_dict[line[0]] = float(line[1])
            else:
                sys.stderr.write("Background amino acid frequency file \"%s\" does not exist.\n"%background_frequency_file)
                sys.exit(0)
        else:
            print ('Missing AA background frequency threshold file. (Example file: McCaldon AA frequencies)')
            print ('Considering equal background frequency threshold for all AA ...')
            for a in self.AAs:
                background_frequency_dict[a] = 10e-100 # nearly 0

        ################################
        # Filter using conservation threshold and store possible choices
        filtered = []
        for idx, val in enumerate(self.pos_aa_frequency):
            pos = idx + 1
            for aa in self.AAs:
                if val[aa] > 0.0 and aa != '-':
                    background_frequency_log_score = (-1.0) * math.log(background_frequency_dict[aa])
                    if not (self.onebody_score[idx,aa, 0] < background_frequency_log_score or self.wt_seq[idx] == aa):
                        filtered.append([pos, aa])
                        self.AA_choices[pos].remove(aa)
        print ("Filtered choices based on amino acid background frequency:")
        print (filtered)
        self.seq_choice_printer()

    ####################################################### ADD AA USING SECONDARY STRUCTURE INFO ###
    #
    # Add choices based upon secondary structure info:
    #   1. Provide secondary structure info for each wildtype residue.
    #       Using this information the module categorizes the positions in four groups: H, S, T and -
    #       H:Helix, S:Sheet, T:Turn, -:Coil
    #   2. Provide secondary structure amino acid propensity values. eg: Chou-Fasman
    #       Using the propensity cutoff value, more choices are added using seconadry structure info.
    #   Example: If position 20 is a helix, then more amino acids choices will be added
    #               using secondary structure propensities.
    #
    def secondary_structure_filter(self, ss_propensity_file=None, ss_propensity_cutoff=1.50):
        sys.stdout.write("\nSecondary structure filter adds mutations to each site.\n\n")
        if self.secondary_structure is None: 
            sys.stdout.write("Target secondary structure information is not given.\n")
            sys.stdout.write("This filter is going to be omitted.\n")
            sys.stdout.write("Please specify secondary structure information (e.g. msa.secondary_structure_information = \"HHHEEETTT\")\n")
            sys.stdout.write("H: Helix, E: Strand, T: Turn (Loop)\n")
        else:
            if ss_propensity_file is None:
                sys.stdout.write("Secondary structure propensity information is not given.\n")
                sys.stdout.write("This filter is going to be omitted.\n")
            elif not os.path.isfile(ss_propensity_file):
                sys.stderr.write("Secondary structure propensity information file \"%s\" does not exist. Please check the file.\n"%ss_propensity_file)
                sys.exit()
            else:
                ################################
                # Reading propensity file
                propensity_dict = {}
                try:
                    propensity_csv = csv.reader(open(ss_propensity_file, 'r'))
                    next(propensity_csv, None) # skip the header
                    for line in propensity_csv:
                        propensity_dict[line[0]] = map(float, line[1:])
                except:
                    print ('Missing aa propensity file.')
                    print ('Provide four column file, i.e., AA, Helix Propensity, Strand Propensity, Turn Propensity')
                print ('propensity_cutoff:', ss_propensity_cutoff)
                helixAAlist = [key for key,values in propensity_dict.iteritems() if values[0]>=ss_propensity_cutoff]
                strandAAlist = [key for key,values in propensity_dict.iteritems() if values[1]>=ss_propensity_cutoff]
                turnAAlist = [key for key,values in propensity_dict.iteritems() if values[2]>=ss_propensity_cutoff]
                print ('AA allowed for helix at the propensity cutoff:', ", ".join(helixAAlist))
                print ('AA allowed for strand at the propensity cutoff:', ", ".join(strandAAlist))
                print ('AA allowed for loop at the propensity cutoff:', ", ".join(turnAAlist))

                for p, ss in enumerate(self.secondary_structure):
                    pos = p+1
                    if ss == 'H':
                        for aa in helixAAlist:
                            if not self.AA_choices[pos].count(aa):
                                self.AA_choices[pos].append(aa)
                    elif ss == 'E':
                        for aa in strandAAlist:
                            if not self.AA_choices[pos].count(aa):
                                self.AA_choices[pos].append(aa)
                    elif ss == 'T':
                        for aa in turnAAlist:
                            if not self.AA_choices[pos].count(aa):
                                self.AA_choices[pos].append(aa)
                    else:
                        for aa in turnAAlist:
                            if not self.AA_choices[pos].count(aa):
                                self.AA_choices[pos].append(aa)

                self.seq_choice_printer()

    ########################################## Epitope score filter
    #
    # Filter possible AA choices based upon:
    #   Minimum number of epitopes that must be deleted my each mutation.
    def epitope_score_filter(self, min_epi_del=0):
        self.minEpiDel = min_epi_del ### <-- used in TUBEscores
        sys.stdout.write("\nEpitope score filter removes mutations to each site.\n\n")
        if self.episcore_mtx is None:
            sys.stdout.write("Epitope score matrix profile is not specified.\n")
            sys.stdout.write("This filter is going to be omitted.\n")
            sys.stdout.write("Please specify profiles if you want to activate this filter (e.g. msa.profiles=PPProfiles([episcore matrix file])\n")
        else:
            wt_episcore = self.profiles.num_hits(self.wt_seq, self.episcore_threshold)
            sys.stdout.write("Epitope score of this target: %d\n"%wt_episcore)
            for p in self.AA_choices:
                pos_muts = self.AA_choices[p][:]
                wt_aa = self.wt_seq[p-1]
                for aa in pos_muts:
                    mut = "%s%d%s"%(wt_aa, p, aa)
                    var_seq = list(self.wt_seq)
                    if (mut in self.mut_lock) or (aa == wt_aa):
                        pass
                    else:
                        var_seq[p-1] = aa
                        var_seq = "".join(var_seq)
                        var_episcore = self.profiles.num_hits(var_seq, self.episcore_threshold)
                        episcore_diff = wt_episcore - var_episcore
                        if episcore_diff >= min_epi_del:
                            print("Mutation %s (%d -> %d)" % (mut, wt_episcore, var_episcore)) 
                        else:
                            self.AA_choices[p].remove(aa)
            self.seq_choice_printer()

    ##################################### COUPLED PAIRS CALCULATION ###
    #
    # After paired log frequency calculation, coupled positions are identified
    # using Chi-squared test filtering at a significance level of 'N-choose-K'
    #
    # Acknowledgment:
    #
    # 'lzprob' and 'lchisqprob' functions were extracted
    # from Gary Strangman's 'stats.py' module.
    # "output" writes an output file in CSV format.
    def chisquare_pair_calculator(self, pair_f=None):
        def lzprob(z):
            """
            Returns the area under the normal curve 'to the left of' the given z value.
            Thus, 
                - for z<0, zprob(z) = 1-tail probability
                - for z>0, 1.0-zprob(z) = 1-tail probability
                - for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
            Adapted from z.c in Gary Perlman's |Stat.
        
            Usage:   lzprob(z)
            """
            Z_MAX = 6.0    # maximum meaningful z-value
            if z == 0.0:
                x = 0.0
            else:
                y = 0.5 * math.fabs(z)
                if y >= (Z_MAX*0.5):
                    x = 1.0
                elif (y < 1.0):
                    w = y*y
                    x = ((((((((0.000124818987 * w
                                -0.001075204047) * w +0.005198775019) * w
                            -0.019198292004) * w +0.059054035642) * w
                            -0.151968751364) * w +0.319152932694) * w
                        -0.531923007300) * w +0.797884560593) * y * 2.0
                else:
                    y = y - 2.0
                    x = (((((((((((((-0.000045255659 * y
                                    +0.000152529290) * y -0.000019538132) * y
                                -0.000676904986) * y +0.001390604284) * y
                                -0.000794620820) * y -0.002034254874) * y
                            +0.006549791214) * y -0.010557625006) * y
                            +0.011630447319) * y -0.009279453341) * y
                        +0.005353579108) * y -0.002141268741) * y
                        +0.000535310849) * y +0.999936657524
            if z > 0.0:
                prob = ((x+1.0)*0.5)
            else:
                prob = ((1.0-x)*0.5)
            return prob

        def lchisqprob(chisq,df):
            """
            Returns the (1-tailed) probability value associated with the provided
            chi-square value and df.  Adapted from chisq.c in Gary Perlman's |Stat.
        
            Usage:   lchisqprob(chisq,df)
            """
            BIG = 20.0
            def ex(x):
                BIG = 20.0
                if x < -BIG:
                    return 0.0
                else:
                    return math.exp(x)
            if chisq <=0 or df < 1:
                return 1.0
            a = 0.5 * chisq
            if df%2 == 0:
                even = 1
            else:
                even = 0
            if df > 1:
                y = ex(-a)
            if even:
                s = y
            else:
                s = 2.0 * lzprob(-math.sqrt(chisq))
            if (df > 2):
                chisq = 0.5 * (df - 1.0)
                if even:
                    z = 1.0
                else:
                    z = 0.5
                if a > BIG:
                    if even:
                        e = 0.0
                    else:
                        e = math.log(math.sqrt(math.pi))
                    c = math.log(a)
                    while (z <= chisq):
                        e = math.log(z) + e
                        s = s + ex(c*z-a-e)
                        z = z + 1.0
                    return s
                else:
                    if even:
                        e = 1.0
                    else:
                        e = 1.0 / math.sqrt(math.pi) / math.sqrt(a)
                    c = 0.0
                    while (z <= chisq):
                        e = e * (a/float(z))
                        c = c + e
                        z = z + 1.0
                    return (c*y+s)
            else:
                return s

        #################################
        # sequence potential must be constructed before calculating chi-square filter
        self.seq_potential_constructor()
        ################################
        # Bonferroni significance level calculation using "N choose K":
        n = int(len(self.wt_seq))+2
        k = int(2)+1
        table = [[0]*(n+2)]*(n+2)

        for i in range(1,n):
            table[i][i] = 1
        for i in range(1,n):
            for j in range(1,n-i):
                x = i+j
                if j == 1: table[x][j] = 1
                else: table[x][j] = table[x-1][j-1] + table[x-1][j]

        sigLevel = 0.01/table[n][k]
        print ('sigLevel for Chi-square comparison', sigLevel)
        
        total_pairs = 0.0
        dep_pairs = 0.0
        ################################
        # Performing chi-square test and filtering pairs using 'N choose K' criteria
        for i in range(len(self.wt_seq)):
            for j in range(i+1, len(self.wt_seq)):
                pos1 = i
                pos2 = j
                countPairs = {}
                # Count the weights for the a, b pair:
                for a in self.AAs:
                    for b in self.AAs:
                        coupling_key = "%d%s%d%s"%(i, a, j, b)
                        countPairs[(a, b)] = self.coupling_count[coupling_key]

                # Sum the rows and columns
                counta = {}; countb = {}
                for a in self.AAs:
                    for b in self.AAs:
                        if a in counta.keys(): counta[a] += countPairs[(a, b)]
                        else: counta[a] = countPairs[(a, b)]
                        if b in countb.keys(): countb[b] += countPairs[(a, b)]
                        else: countb[b] = countPairs[(a, b)]
                
                copya = {}; copyb = {}
                numa = 0; numb = 0
                
                # Count the number of rows greater than 0 and mark rows for copying
                for a in self.AAs:
                    if counta[a] > 0:
                        copya[a] = 1
                        numa += 1
                    else:
                        copya[a] = 0
                        
                # Count the number of cols greater than 0 and mark cols for copying
                for b in self.AAs:
                    if countb[b] > 0:
                        copyb[b] = 1
                        numb += 1
                    else:
                        copyb[b] = 0
                
                # Make the pairTable exactly the number of rows and cols with weight > 0
                pairTableA = {}; pairTable = {}

                # Copy the weights into the correct size table and transform into counts
                global_count = 0
                col_sum = {}
                nexta = -1
                for a in self.AAs:
                    if copya[a] == 1:
                        nexta += 1
                        count_col=0
                        for b in self.AAs:
                            pairTableA[(nexta, b)] = countPairs[(a, b)]
                            count_col = count_col + countPairs[(a, b)]
                        col_sum[nexta] = count_col
                        global_count = global_count + count_col

                row_sum = {}; 
                nextb = -1
                for b in self.AAs:
                    if copyb[b] == 1:
                        nextb += 1
                        count_row = 0
                        for a in range(numa):
                            pairTable[(a, nextb)] = pairTableA[(a, b)]
                            count_row = count_row + pairTableA[(a, b)]
                        row_sum[nextb] = count_row
                
                chi2 = False
                alpha = 0.01
                p = 1.0
                total_pairs += 1
                expectedtable = {}
                
                if numa<2 or numb<2:
                    chi2 = False
                else:
                    for key_row, value_row in row_sum.iteritems():
                        for key_col, value_col in col_sum.iteritems():
                            expectedtable[(key_col, key_row)] = (float(value_row)*float(value_col))/float(global_count)
                    chi_sum = 0.0
                    for k, v in expectedtable.iteritems():
                        chi_sum = chi_sum + (((float(pairTable[k]) - v)**2)/v)
                    dof = ((numa-1)*(numb-1))
                    p = lchisqprob(chi_sum, dof)
                    if p < sigLevel:
                        chi2 = True
                        self.possible_pairs.append((pos1, pos2))
                        dep_pairs +=1
        
        print ('Chi-squared based paired choices implemented ...')
        print ('Total pairs: ', total_pairs)
        print ('Dependent pairs: ', dep_pairs, (dep_pairs/total_pairs))
        print ()
        
        ########################
        # Save pair information to a file
        if pair_f is not None:
            f = pair_f
        else:
            if self.pairs is not None:
                if os.path.isfile(self.pairs):
                    f = self.pairs#+".new" <----- if necessary
                    sys.stdout.write("The pair file \"%s\" exists.\n"%f)
                    sys.stdout.write("New information will be saved in \"%s\".\n\n"%f)
                    sys.stdout.write("Check your pair file!!\n\n")
                else:
                    f = self.pairs
        with open(f, "wb") as pair_file:
            pair_file_line = csv.writer(pair_file)
            pair_lines = map(lambda x: [x[0]+1, x[1]+1], self.possible_pairs)
            # Pair starts from 0. This make residues start at 1
            for pair in pair_lines:
                pair_file_line.writerow(pair)

    #####################################################
    # If a mutational constraint file exists,
    # this filter adds/removes/forces mutations
    def apply_mutational_constraints(self):
        sys.stdout.write("Applying mutational constraints...\n")
        if self.mut_constraint is not None:
            mut_constraint = []
            with open(self.mut_constraint, "r") as f: 
                readAA = csv.reader(f)
                for row in readAA:
                    if len(row):
                        mut_constraint.append(row)
            
            mut_constraint = list(map(list, zip(*mut_constraint))) 
            positions = list(map(int, mut_constraint[0]))  
            exclusively_allowed = list(map(lambda x: list(x), mut_constraint[1])) 
            allowed = list(map(lambda x: list(x), mut_constraint[2]))  
            disallowed = list(map(lambda x: list(x), mut_constraint[3])) 

            warning_message = []
            for pos in self.AA_choices:
                if pos in positions:
                    p = positions.index(pos)
                    if len(exclusively_allowed[p]):
                        self.AA_choices[pos] = exclusively_allowed[p]
                        print ("Exclusively allowed amino acids at %d: %s"%(pos, ", ".join(exclusively_allowed[p])))
                        if self.wt_seq[p] not in exclusively_allowed[p]:
                            warning_message.append("WARNING!!! Wild type amino acid %s at %d was excluded."%(self.wt_seq[p], pos))
                    else:
                        if len(allowed[p]):
                            for aa in allowed[p]:
                                if aa not in self.AA_choices[pos]:
                                    self.AA_choices[pos].append(aa)
                            print ("Allowed amino acids at %d: %s"%(pos, ", ".join(allowed[p])))
                        if len(disallowed[p]):
                            for aa in disallowed[p]:
                                if aa in self.AA_choices[pos]:
                                    self.AA_choices[pos].remove(aa)
                                    if aa == self.wt_seq[p]:
                                        warning_message.append("WARNING!!! Wild type amino acid %s at %d was excluded."%(aa, pos))
                            print ("Disallowed amino acids at %d: %s"%(pos, ", ".join(disallowed[p])))
            print ("\n" + "\n".join(warning_message) + "\n")
            self.seq_choice_printer()
        else:
            sys.stdout.write("Mutational constraint file was not given.\n")
            sys.stdout.write("This filter is deactivated.\n\n")

    def seq_choice_printer(self):
        #pos = self.AA_choices.keys()
        #pos.sort()
        # Sort the keys of the dictionary
        pos = sorted(self.AA_choices.keys())
        num_mutable = 0
        sys.stdout.write("\nMutations at each position:\n\n")
        for p in pos:
            sys.stdout.write("%d: "%(p))
            mut = 0
            for aa in self.AA_choices[p]:
                sys.stdout.write("%s "%aa)
                mut += 1
            sys.stdout.write("\n")
            if mut > 1: num_mutable += 1
        sys.stdout.write("\nTotal mutable sites: %d\n\n"%num_mutable)
    
    ################ Filtered choice
    # This filter perhaps is called at the last stage of preprocess.
    # The variables 'choices' and 'seq_choices'
    # are going to be passed to EpiSweep.
    #
    # The two variables will be substituted if 
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