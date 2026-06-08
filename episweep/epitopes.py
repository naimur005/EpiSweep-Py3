# Fairly generic epitope stuff

# Note: epitope score of 0 means non-binder.  
# However, it's easier to do that than to have None floating around
# (esp. since it can be stuck in a database), and 0 isn't a worthwhile score anyway.
# Any letter other than the 20 AA types results in a 0 for the whole peptide.

import csv

# ----------------------------------------
# A peptide and a dictionary of scores for different alleles (by allele name)
class Epitope :
    def __init__(self, seq, scores):
        self.seq = seq
        self.scores = scores

    # Compare and hash just on sequence (e.g., for storing the epitope in a dictionary), as scores follow
    def __eq__(self, other):
        return self.seq == other.seq
    def __hash__(self):
        return hash(self.seq)

    def __getitem__(self, pos):
        return self.seq[pos]
    
class Profile :
    def __init__(self, allele, matrix, mean, std):
        self.allele = allele
        self.matrix = matrix
        self.mean = mean
        self.std = std

    # Returns the epitope score of the peptide starting at the position in the sequence
    # or 0 if doesn't satisfy the threshold (or has non-standard AA)
    def score(self, seq, pos=0, thresh=0):
        s = -self.mean
        try:
            for i in range(9):
                s += self.matrix[i][seq[pos+i]]
        except KeyError:
            return 0
        s /= self.std
        if s<thresh: return 0
        return s

# ProPred doesn't normalize, but instead compares sum to a percentile threshold.
class ProPred (Profile) :
    def __init__(self, allele, matrix, threshes):
        self.allele = allele
        self.matrix = matrix
        self.threshes = threshes

    # The score is the best threshold satisfied (if one is).
    # Note that lower is better for the threshold, and it can be at most 10.
    def score(self, seq, pos=0, thresh=10):
        s = 0
        try:
            for i in range(9):
                s += self.matrix[i][seq[pos+i]]
        except KeyError:
            return 0
        for t in range(1,thresh+1):
            if s >= self.threshes[t]: return t
        return 0

# ----------------------------------------
# A set of profiles of the same MHC class

class Profiles :
    # Reads a list of scoring matrices
    #   Header row: class:[1,2],by:[pos,aa]
    #   Each matrix has its own header row naming the allele: shortname,longname
    #   Then the scores, either 9x20 (by:pos) or 20x9 (by:aa)
    #   AAs are in alphabetical order by 1 letter codes
    #   Finally the mean and standard deviation
    # If alleles is provided, filters to those listed;
    #   sets all_alleles to indicate whether all from the file are kept
    def __init__(self, filename, alleles=None):
        self.profiles = []
        self.profiles_dict = {}
        self.all_alleles = True
        AAs = "ACDEFGHIKLMNPQRSTVWY"

        fp = csv.reader(open(filename,'rb'))

        # Overall header: MHC class and whether by pos or aa
        row = next(fp)
        self.set_mhc_class(int(row[0][-1]))
        byaa = row[1].endswith('aa')

        # Loop over alleles
        row = next(fp)
        while row:
            allele = row[0]
            matrix = []
            if byaa:
                for p in range(9): matrix.append({})
                for aa in AAs:
                    scores = next(fp)
                    for p in range(9): matrix[p][aa] = float(scores[p])
            else:
                for p in range(9):
                    matrix.append(dict(zip(AAs, map(float, next(fp)))))
            [mean, std] = map(float, fp.next())
            if alleles is None or allele in alleles: 
                profile = Profile(allele, matrix, mean, std)
                self.profiles.append(profile)
                self.profiles_dict[allele] = profile
            else:
                self.all_alleles = False
            row = next(fp,None)

    def set_mhc_class(self, c):
        self.mhc_class = c
        if c==1:
            self.up_positions = [3,4,5,6,7]
            self.down_positions = [0,1,2,8]
        elif c==2:
            self.up_positions = [1,2,4,6,7]
            self.down_positions = [0,3,5,8]
        else:
            raise Exception('unsupported mhc class')

    def __iter__(self):
        return self.profiles.__iter__()

    def __getitem__(self, allele):
        return self.profiles_dict[allele]

    def __len__(self):
        return len(self.profiles)

    # Convert percentile threshold to score threshold
    zthreshes = [0,2.326348,2.053749,1.880793,1.750686,1.644853,1.554773,1.475791,1.405071,1.340755,1.281551]
    def pct2thresh(self, pct):
        return Profiles.zthreshes[pct]
    
    # Returns a dictionary mapping a position (0-indexed) to a list of epitope scores
    # (one for each profile, in order by allele name) for those positions where
    # at least one epitope score meets the threshold.
    def episcan(self, seq, thresh):
        epitopes = {}
        for i in range(len(seq)-8):
            scores = dict((p.allele, p.score(seq, i, thresh)) for p in self.profiles)
            if sum(scores.itervalues())>0: epitopes[i]=Epitope(seq[i:i+9],scores)
        return epitopes

    # Returns a dictionary mapping a position (0-indexed) to the number of hits
    # (over all alleles) at that position, for those positions where at least
    # one epitope score meets the threshold.
    def num_hits_scan(self, seq, thresh):
        nhits = {}
        for i in range(len(seq)-8):
            nh = sum(p.score(seq, i, thresh) for p in self.profiles)
            if nh>0: nhits[i] = nh
        return nhits

    # Returns the total number of hits (over all positions and alleles) for the sequence.
    def num_hits(self, seq, thresh):
        total = 0
        for i in range(len(seq)-8):
            for p in self.profiles:
                total += p.score(seq, i, thresh)>0
        return total

    # Returns the total number of hits (over all alleles) for the single peptide.
    def peptide_hits(self, peptide, thresh):
        return sum(p.score(peptide, 0, thresh)>0 for p in self.profiles)

    # Returns a list (one per position in the sequence) of lists of scores (one per allele)
    def full_report(self, seq, thresh):
        return [[p.score(seq, i, thresh) for p in self.profiles] for i in range(len(seq)-8)]

    # Prints out a full report (see above) for testing
    def print_report(self, seq, thresh):
        for i in range(len(seq)-8):
            scores = [p.score(seq, i, thresh) for p in self.profiles]
            hits = sum(1 for s in scores if s>0)
            print (str(i+1)+','+seq[i:i+9]+','+str(i+9)+','+','.join(map(lambda s: '%.02f'%s,scores))+','+str(hits))

    # Returns whether or not the position (0-indexed) is TCR-facing ("up")
    def is_up(self, pos):
        return pos in self.up_positions

    # Returns a string of the TCR-facing residues ("ups") for the Epitope
    # TO-DO: assumes same ones across all alleles -- is that good enough?
    def up_residues(self, epitope):
        return ''.join(epitope[i] for i in self.up_positions)

# The ProPred version (with ProPred instances, and correspondingly different file format)
class PPProfiles (Profiles) :

    # Reads a list of scoring matrices, each in format
    #  Allele
    #  A1 A2 ... A9
    #  C1 C2 ... C9
    #  ...
    #  Y1 Y2 ... Y9
    #  t1 t2 ... t10
    # where each row has the contributions from one amino acid type for the 9 pockets,
    # and the final line has the thresholds for 1% up to 10%
    def __init__(self, filename, alleles=None):
        AAs = "ACDEFGHIKLMNPQRSTVWY"
        self.all_alleles = True
        self.profiles = []
        self.profiles_dict = {}
        self.allele_list = []
        self.set_mhc_class(2)
        self.filename = filename
        if filename is None:
            self.filename = None
        else:
            with open(filename,'r') as fp:
                epitope_mtx = csv.reader(fp)
                while True:
                    try:
                        allele = "-".join(filter(lambda x: x != "", next(epitope_mtx)))
                        matrix = [{} for p in range(9)]
                        for aa in AAs:
                            scores = list(filter(lambda x: x != "", next(epitope_mtx)))
                            for p in range(9):
                                matrix[p][aa] = float(scores[p])
                        threshes = list(map(float, filter(lambda x: x != "", next(epitope_mtx))))
                        threshes.insert(0, 0)
                        if alleles is None or allele in alleles: 
                            profile = ProPred(allele, matrix, threshes)
                            self.profiles.append(profile)
                            self.profiles_dict[allele] = profile
                            self.allele_list.append(allele)
                        else:
                            self.all_alleles = False
                    except StopIteration:
                        break

    def pct2thresh(self, pct): return pct