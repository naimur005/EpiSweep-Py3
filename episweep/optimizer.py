import cplex
import itertools
from .msa import *
from .scores import *
from .epitopes import *
from .TUBEscores import *
from .PDBParser import *
import os
import sys
import math

class EpiSweep (TUBE):
    # inputs
    # wt -- amino acid sequence (1-letter string)
    # n -- length of wt
    # nmuts -- number of positions at which to have choices other than just wt
    # max_episcore -- maximum epitope score
    # choices -- position-specific lists of allowed Choices
    # seq_choices -- position-specific lists of allowed SeqChoices
    # position_pairs -- list of (i,j) pairs for which to evaluate twobody

    # methods to be defined by subclasses for scoring
    # onebody -- position, choice -> score
    # twobody -- position i, choice i, position j, choice j -> score
    # episcore -- position, tuple of 9 choices -> score
    
    # internal state
    # opt -- the Cplex instance

    # results
    # chosen -- list of Choices (one per position)

    # create the optimization problem
    def initialize_optimizer(self):
        self.wt_seq = "".join(map(lambda x: x.strip("\n"), open(self.target).readlines()[1:]))
        if self.design_type == "disruption":
            pass
        else:
            self.profiles = PPProfiles(self.episcore_mtx)

        self.n = len(self.wt_seq)
        self.possible_pairs = map(lambda x: (int(x.split(",")[0])-1, int(x.split(",")[1])-1), open(self.pairs).readlines())
        self.initialize_score_potential()
        if self.library_gen:
            self.load_tube_score_potential()
        else:
            if os.path.splitext(self.score_potential)[1] in [".csv", ".CSV"]:
                self.load_csv_score_potential()
            else: self.load_binary_score_potential()
        self.load_epitope_score_list()

    def optimize(self):
        self.opt = cplex.Cplex()
        self.opt.objective.set_sense(self.opt.objective.sense.minimize)

        self.epsilon = 0.0001 # <- This tiny value fixes discrepancies of binary variables
        self.add_singles()
        self.add_pairs()
        if self.design_type == "disruption":
            self.add_disruption()
        elif self.design_type == "humanization": # This module will be implemented later
            self.add_humanness()
        else:
            self.add_epitopes()
        self.add_mutload(0) # <- calculate wild type
        self.sweep() # <- calculate wild type
        self.opt.linear_constraints.delete('mutload') # <- remove the constraint and reset
        self.add_mutload()
        if self.library_gen:
            self.add_library_size()
        self.sweep(self.num_curves)

        print()
        log_file = open(self.designs, "w")
        header = "index,curve,episcore,potential,mutation,sequence\n"
        sys.stdout.write(header)
        log_file.write(header)
        design_index = 0
        for curve_id, curve in enumerate(self.curve):
            for i in curve:
                l = [design_index, curve_id, i.episcore, i.score_potential, "|".join(i.mutation), i.seq]
                l = map(str, l)
                l = ",".join(l)+"\n"
                sys.stdout.write(l)
                log_file.write(l)
                design_index += 1
        log_file.close()
        sys.stdout.write("\nResults are saved in \"%s\".\n\n"%self.designs)

    # single-position variables, coefficients, and constraints
    def add_singles(self):
        for i in range(self.n):
            # handle all variables at position i together, so can enforce constraint of one choice per position
            # print self.choices[i]
            ivars = ['s%d_%s' % (i,c) for c in self.choices[i]]
            icoeffs = [self.onebody(i,c) for c in self.choices[i]]
            self.opt.variables.add(names=ivars, types=['I']*len(ivars), obj=icoeffs, lb=[0]*len(ivars), ub=[1]*len(ivars))
            self.opt.linear_constraints.add(lin_expr=[[ivars,[1]*len(ivars)]], senses=['E'], rhs=[1])

    # pairwise variables, coefficients, and constraints
    def add_pairs(self):
        pair_length = len(list(self.possible_pairs))
        sys.stdout.write("There are %d pairs.\n"%pair_length)
        count = 1
        for i, j in self.possible_pairs:
            # The first position in "possible_pairs" starts from 1.
            # Needs to be subtracted by 1.
            sys.stdout.write("%d out of %d pair variables have been added\n"%(count, pair_length))
            count += 1
            # handle all variables at positions i,j together, so can enforce consistency constraints
            ijvars = []; ijcoeffs = []
            cvars = {}; dvars = {}  # cvars: c->list of vars involving it, similarly for dvars
            #sys.stdout.write("%d %d %d %d\n"%(i+1, j+1, len(self.choices[i]), len(self.choices[j])))
            for c in self.choices[i]:
                for d in self.choices[j]:
                    var = 'p%d_%s_%d_%s' % (i,c,j,d)
                    ijvars.append(var)
                    ijcoeffs.append(self.twobody(i,c,j,d))
                    if c in cvars: cvars[c].append(var)
                    else: cvars[c] = [var]
                    if d in dvars: dvars[d].append(var)
                    else: dvars[d] = [var]
            # create binary variables
            self.opt.variables.add(names=ijvars, types=['I']*len(ijvars), obj=ijcoeffs, lb=[0]*len(ijvars), ub=[1]*len(ijvars))
            # consistency constraint: for each c, sum over its pair vars equals its singleton
            for c,vars in cvars.iteritems():
                svar = 's%d_%s'%(i,c)
                self.opt.linear_constraints.add(lin_expr=[[[svar]+vars,[-1]+[1]*len(vars)]], senses=['E'], rhs=[0])
            # consistency constraint: for each d, sum over its pair vars equals its singleton
            for d,vars in dvars.iteritems():
                svar = 's%d_%s'%(j,d)
                self.opt.linear_constraints.add(lin_expr=[[[svar]+vars,[-1]+[1]*len(vars)]], senses=['E'], rhs=[0])

    # constrain the total number of non-wt choices
    # constraint is named "mutload" so the value can be modified
    def add_mutload(self, nmut=None):
        if nmut is not None:
            mut_vars = ['s%d_%s' % (i,c) for i in range(self.n) for a in self.seq_choices[i] if not a.is_wt() for c in a.choices]
            self.opt.linear_constraints.add(lin_expr=[[mut_vars,[1]*len(mut_vars)]], senses=['E'], rhs=[nmut], names=['mutload'])
        else:
            if self.mut_load is None:
                sys.stderr.write("\nMutational loads were not specified.\n")
                sys.stderr.write("Free mutational loads will be assigned.\n")
            else:
                mut_vars = ['s%d_%s' % (i,c) for i in range(self.n) for a in self.seq_choices[i] if not a.is_wt() for c in a.choices]
                if type(self.mut_load) == int:
                    self.opt.linear_constraints.add(lin_expr=[[mut_vars,[1]*len(mut_vars)]], senses=['E'], rhs=[self.mut_load], names=['mutload'])
                elif type(self.mut_load) == list:
                    max_m = max(self.mut_load)
                    min_m = min(self.mut_load)
                    self.opt.linear_constraints.add(lin_expr=[[mut_vars,[1]*len(mut_vars)]], senses=['G'], rhs=[min_m], names=['mutload'])
                    self.opt.linear_constraints.add(lin_expr=[[mut_vars,[1]*len(mut_vars)]], senses=['L'], rhs=[max_m], names=['mutload'])
                else:
                    sys.stderr.write("\nMutation loads must be given as an integer number,")
                    sys.stderr.write("or list or tuple.\n")
                    sys.exit()

    # constrain the library size. Only non-wt choices are considered because size of wt choice is 1 => log(1) => 0
    def add_library_size(self):
        print (self.library_size)
        minlib = self.library_size[0]; maxlib = self.library_size[1]
        lib_vars = ['s%d_%s' % (i,c) for i in range(self.n) for a in self.seq_choices[i] if not a.is_wt() for c in a.choices]
        lib_coeff_max = [math.log(len(set(a.aas))) for i in range(self.n) for a in self.seq_choices[i] if not a.is_wt()]
        lib_coeff_min = [((-1.0)*math.log(len(set(a.aas)))) for i in range(self.n) for a in self.seq_choices[i] if not a.is_wt()]
        #lib_coeff_max = [math.log(len(set(a.aas))-1) for i in xrange(self.n) for a in self.seq_choices[i] if not a.is_wt()]
        #lib_coeff_min = [((-1.0)*math.log(len(set(a.aas))-1)) for i in xrange(self.n) for a in self.seq_choices[i] if not a.is_wt()]
        self.opt.linear_constraints.add(lin_expr=[[lib_vars, lib_coeff_max]], senses=['L'], rhs=[math.log(maxlib)], names=['maxlibsize'])
        self.opt.linear_constraints.add(lin_expr=[[lib_vars, lib_coeff_min]], senses=['L'], rhs=[(-1)*math.log(minlib)], names=['minlibsize'])

    # epitope variables and constraints
    # no contribution to objective function in episweep formulation; instead, have constraint named "episcore" that constrains the sum; set its rhs appropriately
    def tree_copy(self, i, depth=-1):
        print
        tuples = []
        def recurse(i, depth, aamut, t):
            if depth != 9:
                for aa in self.seq_choices[i+depth]:
                    if not aa.is_wt():
                        aamut += 1
                    print (aamut, depth, t, aa)
                    if aamut <= 2:
                        t.append(str(aa))
                        #t.append(str(aa))
                    else:
                        aamut -= 1

                    recurse(i, depth+1, aamut, t)
            else:
                if len(t) == 9:
                    tuples.append(t)
                    print ("success!!", aamut, "".join(t), depth)
                    t = []
                    depth=0
        recurse(i, 0, 0, []) # starts the recursion
        return tuples

    def tree_cut(self, i):
        tuples = []
        mut9merlimit=2
        possibleMuts = 0
        for aa in self.seq_choices[i+0]:
            aamut = 0
            if self.library_gen == True:
                if len(aa.aas)>1: aamut =1
            else:
                if not aa.is_wt(): aamut =1
            if aamut<=mut9merlimit:
                t = []
                t.append(aa)
                for bb in self.seq_choices[i+1]:
                    bbmut = 0
                    if self.library_gen == True:
                        if len(bb.aas)>1: bbmut =1
                    else:
                        if not bb.is_wt(): bbmut =1
                    if (aamut+bbmut)<=mut9merlimit:
                        t = t[0:1]
                        t.append(bb)
                        for cc in self.seq_choices[i+2]:
                            ccmut = 0
                            if self.library_gen == True:
                                if len(cc.aas)>1: ccmut =1
                            else:
                                if not cc.is_wt(): ccmut =1
                            if (aamut+bbmut+ccmut)<=mut9merlimit:
                                t = t[0:2]
                                t.append(cc)
                                for dd in self.seq_choices[i+3]:
                                    ddmut = 0
                                    if self.library_gen == True:
                                        if len(dd.aas)>1: ddmut =1
                                    else:
                                        if not dd.is_wt(): ddmut =1
                                    if (aamut+bbmut+ccmut+ddmut)<=mut9merlimit:
                                        t = t[0:3]
                                        t.append(dd)
                                        for ee in self.seq_choices[i+4]:
                                            eemut = 0
                                            if self.library_gen == True:
                                                if len(ee.aas)>1: eemut =1
                                            else:
                                                if not ee.is_wt(): eemut =1
                                            if (aamut+bbmut+ccmut+ddmut+eemut)<=mut9merlimit:
                                                t = t[0:4]
                                                t.append(ee)
                                                for ff in self.seq_choices[i+5]:
                                                    ffmut = 0
                                                    if self.library_gen == True:
                                                        if len(ff.aas)>1: ffmut =1
                                                    else:
                                                        if not ff.is_wt(): ffmut =1
                                                    if (aamut+bbmut+ccmut+ddmut+eemut+ffmut)<=mut9merlimit:
                                                        t = t[0:5]
                                                        t.append(ff)
                                                        for gg in self.seq_choices[i+6]:
                                                            ggmut = 0
                                                            if self.library_gen == True:
                                                                if len(gg.aas)>1: ggmut =1
                                                            else:
                                                                if not gg.is_wt(): ggmut =1
                                                            if (aamut+bbmut+ccmut+ddmut+eemut+ffmut+ggmut)<=mut9merlimit:
                                                                t = t[0:6]
                                                                t.append(gg)
                                                                for hh in self.seq_choices[i+7]:
                                                                    hhmut = 0
                                                                    if self.library_gen == True:
                                                                        if len(hh.aas)>1: hhmut =1
                                                                    else:
                                                                        if not hh.is_wt(): hhmut =1
                                                                    if (aamut+bbmut+ccmut+ddmut+eemut+ffmut+ggmut+hhmut)<=mut9merlimit:
                                                                        t = t[0:7]
                                                                        t.append(hh)
                                                                        for ii in self.seq_choices[i+8]:
                                                                            iimut = 0
                                                                            if self.library_gen == True:
                                                                                if len(ii.aas)>1: iimut =1
                                                                            else:
                                                                                if not ii.is_wt(): iimut =1
                                                                            if (aamut+bbmut+ccmut+ddmut+eemut+ffmut+ggmut+hhmut+iimut)<=mut9merlimit:
                                                                                t = t[0:8]
                                                                                t.append(ii)
                                                                                #print i, t, possibleMuts
                                                                                tuples.append(t)
        return tuples

    
    def tree_cut_2(self, i, mut_pepmer_limit=3):
        choices_in_nmer = []
        for pos in range(i,i+self.nmer):
            pos_choices = []
            for aa in self.seq_choices[pos]:
                pos_choices.append(aa)
            choices_in_nmer.append(pos_choices)
        tuples_temp = list(itertools.product(*choices_in_nmer))
        tuples = []
        for t in tuples_temp:
            codon_lengths = [1 if len(list(set(c.aas)))>1 else 0 for c in t]
            if sum(codon_lengths)<=mut_pepmer_limit: tuples.append(t)
        return tuples

    def disruption_score(self, i, t):
        score = 0
        i += 1
        contact_list = map(int, self.interface_contact_list.keys())
        contact_list.sort()
        for resn in range(len(t)):
            if resn+i in contact_list:
                for res in self.interface_contact_list[str(resn+i)]:
                    try:
                        score += self.disruption_score_matrix[str(t[resn])][res]
                    except KeyError as e:
                        raise Exception("Error on res %d with value %s contacting %s" % (resn+1, contact, str(t[resn])))
        print ("%s at %d has a score of %f" % ("".join([str(x) for x in t]), i, score))
        return score

    def add_disruption(self):
        epivars = []; episcores = []
        for i in range(self.n):               
            epitope_ranges = [0]
            tuples = [(opt,) for opt in self.seq_choices[i]]
            
            tvars = ['w%d_%d'%(i,t) for t in range(len(tuples))]
            tscores = [self.disruption_score(i, t) for t in tuples]
            
            if any(score != 0 for score in tscores):
 #FIXME: Should this be != 0 or just not here?
                print (len(tuples), "combinations")
                # keep track of those that are possible epitopes, in order to constrain total
                for t in range(len(tuples)):
                    if tscores[t] != 0:
 #FIXME: Should this be != 0? Or just not be here?
                        epivars.append(tvars[t])
                        episcores.append(tscores[t])
                # consistency constraint: for each position in tuple and possible seq choice there, sum over tuples with that choice equals sum over singletons for concrete choices for that seq choice

                # if no possible hits at this position, for any combination of mutations, then don't need to set up variables
                # (this generalizes Andrew's approach to ignoring positions without a possible hydrophobic at P1)
                self.opt.variables.add(names=tvars, types=['I']*len(tvars), lb=[-1]*len(tvars), ub=[1]*len(tvars))

                # note that this isn't what's written in the rotamer-based paper :( but it is consistent with the special case of sequence-based (IP^2)
                for h in epitope_ranges:
                    for a in self.seq_choices[i+h]:
                        wvars = ['w%d_%d'%(i,t) for t in range(len(tuples)) if tuples[t][h] == a]
                        svars = ['s%d_%s'%(i+h,c) for c in a.choices]
                        self.opt.linear_constraints.add(lin_expr=[[svars+wvars,[-1]*len(svars)+[1]*len(wvars)]], senses=['E'], rhs=[0])
            else:
                print (0, "combinations")
        # create episcore constraint: sum of epitope scores <= max_episcore
        self.opt.variables.add(names=['episcore'], types=[self.opt.variables.type.continuous], lb=[self.max_episcore])
        self.opt.linear_constraints.add(lin_expr=[[['episcore']+epivars,[-1]+episcores]], senses=['E'], rhs=[0])
 #episcore starts at 0?
        self.opt.linear_constraints.add(lin_expr=[[['episcore'],[1]]], senses=['G'], rhs=[self.max_episcore], names=['sweep_constraint'])
 #New episcore must be greater than wt. NOTE: changed L->G, so max_episcore is actually min_episcore

    def add_humanness(self):
        epivars = []; episcores = []
        for i in range(self.n-self.nmer + 1):
            tuples = map(lambda x: x[0], self.precomputed_peptide_list[i])
            epitope_ranges = range(self.nmer)
            tvars = ['w%d_%d'%(i,t) for t in range(len(tuples))]
            tscores = map(lambda x: x[1], self.precomputed_peptide_list[i])
            if sum(tscores)>0:
                print (len(tuples), "combinations")
                for t in range(len(tuples)):
                    if tscores[t]>0:
                        epivars.append(tvars[t])
                        episcores.append(tscores[t])
                self.opt.variables.add(names=tvars, types=['I']*len(tvars), lb=[0]*len(tvars), ub=[1]*len(tvars))
                for h in epitope_ranges:
                    for a in self.seq_choices[i+h]:
                        wvars = ['w%d_%d'%(i,t) for t in range(len(tuples)) if tuples[t][h] == str(a)]
                        svars = ['s%d_%s'%(i+h,c) for c in a.choices]
                        self.opt.linear_constraints.add(lin_expr=[[svars+wvars,[-1]*len(svars)+[1]*len(wvars)]], senses=['E'], rhs=[0])
            else:
                print (0, "combinations")
        self.opt.variables.add(names=['episcore'], types=[self.opt.variables.type.continuous], lb=[self.max_episcore])
        self.opt.linear_constraints.add(lin_expr=[[['episcore']+epivars,[-1]+episcores]], senses=['E'], rhs=[0])
 #episcore starts at 0?
        self.opt.linear_constraints.add(lin_expr=[[['episcore'],[1]]], senses=['G'], rhs=[self.max_episcore], names=['sweep_constraint'])
 #New episcore must be greater than wt. NOTE: changed L->G, so max_episcore is actually min_episcore

    def add_epitopes(self, propred_option=True):
        sys.stdout.write("\nAdding epitope variable.\n")

        multi_indexed = False
        if self.library_gen: multi_indexed = False
        for p in range(self.n):
            for aa in self.seq_choices[p]:
                if len(aa.choices) > 1:
                    multi_indexed = True
                    break

        epivars = []; episcores = []
        for i in range(self.n-self.nmer+1):
            #print i+1, "out of", self.n-8,
            if self.precomputed_epitope is not None:
                if self.library_gen:
                    tuples = self.tree_cut_2(i, mut_pepmer_limit=3)
                    epitope_ranges = range(self.nmer)
                else:
                    tuples = map(lambda x: x[0], self.precomputed_peptide_list[i])
                    epitope_ranges = range(self.nmer)
            else:
                if not multi_indexed:
                    tuples = self.tree_cut(i)
                    epitope_ranges = range(9)
                    #tuples = self.tree_copy(i, 9)
                    #print
                    #print map(lambda x: "".join(map(str, list(x))), tuples)
                else:
                    if propred_option:
                        seq_combination = []
                        tscores = []
                        tvars = []
                        epitope_ranges = []
                        for h in range(self.nmer):
                            if h == 4 or h == 7:
                                seq_combination.append([self.seq_choices[i+h][0]])
                            else:
                                seq_combination.append(self.seq_choices[i+h])
                                epitope_ranges.append(h)
                        tuples = list(itertools.product(*seq_combination))
                        
                    else:
                        tuples = list(itertools.product(*[self.seq_choices[i + h] for h in range(9)]))

                        epitope_ranges = range(self.nmer)

            tvars = ['w%d_%d'%(i,t) for t in range(len(tuples))]
            if self.precomputed_epitope is not None:
                if self.library_gen:
                    print (i, '...')
                    tscores = [self.episcore(i, t, 0) for t in tuples]
                else:
                    tscores = map(lambda x: x[1], self.precomputed_peptide_list[i])
            else:
                tscores = [self.episcore(i, t, self.episcore_threshold) for t in tuples]

            if sum(tscores)>0:
                #print len(tuples), "combinations"
                # keep track of those that are possible epitopes, in order to constrain total
                for t in range(len(tuples)):
                    #print i, t, tscores[t], len(tuples)
                    #print len(self.precomputed_peptide_list[i])
                    #print len(tscores), tscores[t]
                    if tscores[t]>0:
                        epivars.append(tvars[t])
                        episcores.append(tscores[t])
                # consistency constraint: for each position in tuple and possible seq choice there, sum over tuples with that choice equals sum over singletons for concrete choices for that seq choice

                # if no possible hits at this position, for any combination of mutations, then don't need to set up variables
                # (this generalizes Andrew's approach to ignoring positions without a possible hydrophobic at P1)
                self.opt.variables.add(names=tvars, types=['I']*len(tvars), lb=[0]*len(tvars), ub=[1]*len(tvars))

                # note that this isn't what's written in the rotamer-based paper :( but it is consistent with the special case of sequence-based (IP^2)
                #print len(tvars), tvars
                for h in epitope_ranges:
                    for a in self.seq_choices[i+h]:
                        if self.precomputed_epitope is not None:
                            if self.library_gen:
                                wvars = ['w%d_%d'%(i,t) for t in range(len(tuples)) if tuples[t][h] == a]
                            else:
                                wvars = ['w%d_%d'%(i,t) for t in range(len(tuples)) if tuples[t][h] == str(a)]
                        else:
                            wvars = ['w%d_%d'%(i,t) for t in range(len(tuples)) if tuples[t][h] == a]
                        svars = ['s%d_%s'%(i+h,c) for c in a.choices]
                        #print "wvars: %d, svars: %s"%(len(wvars), len(svars))
                        #print wvars, svars
                        self.opt.linear_constraints.add(lin_expr=[[svars+wvars,[-1]*len(svars)+[1]*len(wvars)]], senses=['E'], rhs=[0])
           #else:
           #    print 0, "combinations"
        # create episcore constraint: sum of epitope scores <= max_episcore
        self.opt.variables.add(names=['episcore'], types=[self.opt.variables.type.continuous])
        #self.opt.variables.add(names=['episcore'], types=[self.opt.variables.type.integer])
        self.opt.linear_constraints.add(lin_expr=[[['episcore']+epivars,[-1]+episcores]], senses=['E'], rhs=[0])
        self.opt.linear_constraints.add(lin_expr=[[['episcore'],[1]]], senses=['L'], rhs=[self.max_episcore], names=['sweep_constraint'])

    # solve the problem and extract the choices
    def solve(self):
        self.opt.solve()
        print ("Solution status =", self.opt.solution.get_status(), ":", self.opt.solution.status[self.opt.solution.get_status()])
        # visit the following address in order to check get_status() return values: http://pic.dhe.ibm.com/infocenter/cosinfoc/v12r2/index.jsp?topic=%2Filog.odms.cplex.help%2FContent%2FOptimization%2FDocumentation%2FCPLEX%2F_pubskel%2FCPLEX1210.html
        if self.opt.solution.get_status() == 103: return False # A curve finishes
        
        self.energy_value = self.opt.solution.get_objective_value()
        print ("Objective value =", self.energy_value)

        # extract choices and episcore
        vars = self.opt.variables.get_names()
        vals = self.opt.solution.get_values()
        self.chosen = [None]*self.n
        errors = []
        for v in range(len(vars)):
            # This checks any multiple choices, i.e. any singletons which have multiple 1s at a same position. Solutions may become non binary and so the epsilon value is needed.
            if (not (vals[v]>-self.epsilon and vals[v] <self.epsilon)) and vars[v][0]=='s':
                [pos,choice] = vars[v][1:].split('_')
                pos = int(pos)
                if self.chosen[pos] is not None:
                    errors.append('multiple choices at '+str(pos)+'%s. '%(self.chosen[pos]))
                self.chosen[pos] = choice

            elif (vals[v] > -self.epsilon or vals[v] < self.epsilon) and (vals[v] > 0 and vals[v] != 1.0) and vars[v][0] == 's':
                [pos,choice] = vars[v][1:].split('_')
                pos = int(pos)
                print ("Solution values are not binary at", pos, vars[v], vals[v])

            elif vars[v]=='episcore':
                self.episcore = vals[v]
        for pos in range(self.n):
            if self.chosen[pos] is None: errors.append('error: nothing chosen at '+str(pos))
        if len(errors)>0: raise Exception('.'.join(errors))

        print ('Episcore:', self.episcore)

        print ('solution:')
        for i,c in enumerate(self.chosen): print (str(i+1)+c)
        print()

        return True

    def sweep(self, ncurve=1, rotamer_degeneracy=True):
        if not hasattr(self, "curve"):
            self.curve = []
        for s in range(ncurve):
            sweep_line = []
            while self.solve():
                mvars = []
                rotvars = []
                mutated = 0
                var_seq = list(self.wt_seq)[:]
                pareto_design = plan()
                pareto_design.mutation = []
                pareto_design.score_potential = self.energy_value
                pareto_design.episcore = self.episcore
                for i in range(self.n):
                    if self.library_gen:
                        pass_condition = (self.chosen[i] != self.wt_seq[i])
                    else: pass_condition = (self.chosen[i][0] != self.wt_seq[i])
                    if pass_condition:
                        if self.library_gen: mutant = "%s%d%s"%(self.wt_seq[i], i+1, self.chosen[i])
                        else: mutant = "%s%d%s"%(self.wt_seq[i], i+1, self.chosen[i][0])
                        if self.library_gen: var_seq[i] = "*"
                        else: var_seq[i] = self.chosen[i][0]
                        pareto_design.mutation.append(mutant)
                        #if self.library_gen: print 'mutation @',i+1,': chosen',str(self.chosen[i]), ''.join(self.tube_dict[self.chosen[i]]), self.tube_ratio[i, self.chosen[i]],'; wt',self.wt_seq[i]
                        if self.library_gen:
                            if str(self.chosen[i]).isupper(): # For combinatorial method
                                print ('mutation @',i+1,': chosen', str(self.chosen[i]), str(self.chosen[i]), '; wt',self.wt_seq[i])
                            else: # For degenerate method
                                print ('mutation @',i+1,': chosen', str(self.chosen[i]), ''.join(self.all_tubes_dict[self.chosen[i]]), '; wt',self.wt_seq[i])

                        else: print ('mutation @',i+1,': chosen',str(self.chosen[i]),'; wt',self.wt_seq[i])
                        mutated += 1
                        mvars.append("s%d_%s"%(i, self.chosen[i]))
                        if (not self.library_gen) and rotamer_degeneracy:
                            for c in self.choices[i]:
                                if c.name[0] == self.chosen[i][0]:
                                    rotvars.append("s%d_%s"%(i, c.name))

                pareto_design.seq = "".join(var_seq)
                sweep_line.append(pareto_design)
                print ()
                print ("Mutations: %s"%(" ".join(pareto_design.mutation)))
                print ("Sequence: %s"%pareto_design.seq)
                print ()
                print ("Design output:", ",".join(map(str, [len(self.curve), pareto_design.episcore, pareto_design.score_potential, "|".join(pareto_design.mutation), pareto_design.seq])))
                print ()

                if len(mvars) > 0:
                    self.opt.linear_constraints.add(lin_expr=[[mvars,[1]*len(mvars)]], senses=['L'], rhs=[mutated-1]) # This constraint makes other solutions at least 1 mutated rotamer is different
                    if (not self.library_gen) and rotamer_degeneracy:
                        self.opt.linear_constraints.add(lin_expr=[[rotvars,[1]*len(rotvars)]], senses=['L'], rhs=[mutated-1]) # This constraint makes at least 1 amino acid choice different from other solutions
                
                if self.library_gen:
                    # Calculate library size of the solution variant and set constraint for the next round
                    self.final_libsize = 1
                    seq_final_list = []
                    for i,c in enumerate(self.chosen):
                        if c != self.wt_seq[i]:
                            if c.isupper(): # For combinatorial method
                                self.final_libsize = self.final_libsize * len(set(list(c)))
                                seq_final_list.append(list(list(c)))
                            else: # For degenerate method
                                self.final_libsize = self.final_libsize * len(set(self.all_tubes_dict[c]))
                                seq_final_list.append(list(self.all_tubes_dict[c]))
                        else:
                            seq_final_list.append(list(c))
                    #self.libraryNumHits(seq_final_list, self.episcore_threshold)
                    print ('Library Size:', self.final_libsize)
                    maxEpi = (((self.episcore*self.library_size[1])-(float(self.library_size[1])/self.library_size[0]))/float(self.library_size[1]))
                    print ('maxEpi for next run=', maxEpi)
                    self.opt.linear_constraints.set_rhs('sweep_constraint', maxEpi)
                else:
                    self.opt.linear_constraints.set_rhs('sweep_constraint', self.episcore+self.sweep_direction/self.score_scale)
            self.curve.append(sweep_line)

            if self.design_type == "disruption" or self.design_type == "humanization":
                self.opt.linear_constraints.set_rhs('sweep_constraint', self.max_episcore+self.sweep_direction/self.score_scale)
            else:
                self.opt.linear_constraints.set_rhs('sweep_constraint', self.max_episcore)
        sys.stderr.write("Sweeping done.\n")

    def libraryNumHits (self, final_seq, thresh):
        print ('Library NumHits:')
        best_min_score = []
        best_max_score = []
        print ('Position', 'AvgNumHit', 'NumHitStDev')
        for i in range(0, len(final_seq)-8):
            all_tubes = []
            norm = 1.0
            for j in range(i, i+9):
                all_tubes.append(list(set(final_seq[j])))
                norm = norm * len(final_seq[j])
            possible_peptides = list(itertools.product(*all_tubes))
            score = 0.0
            peptide_episcore_list = []
            for p in range(len(possible_peptides)):
                peptide = (''.join(c for c in possible_peptides[p]))
                peptide_episcore = self.profiles.peptide_hits(peptide, thresh)
                score = score + (peptide_episcore/norm)
                peptide_episcore_list.append(peptide_episcore)

            avg = self.average(peptide_episcore_list)
            variance = map(lambda x: (x - avg)**2, peptide_episcore_list)
            standard_deviation = math.sqrt(self.average(variance))
            best_min_score.append(min(peptide_episcore_list))
            best_max_score.append(max(peptide_episcore_list))
            print ((i+1), score, standard_deviation)
        print ('Best_min_episcore:', sum(best_min_score))
        print ('Best_max_episcore:', sum(best_max_score))

    def average(self, s): return sum(s) * 1.0 / len(s)

class plan:
    def __init__(self): pass
    