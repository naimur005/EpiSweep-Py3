from .scores import *
from .epitopes import *
from .PDBParser import *
from .DEE import *
from .msa import *
from .options import *
from .optimizer import *
from .choice import *
from .sausage_epitope import *
from .episcope_score import *
from .k_medoids import *

class Design (EpiSweep, ScorePotential, MSA, PDB, TUBE):
    def __init__(self):
        self.target = None
        self.episcore_mtx = None
        self.episcore_threshold = None
        self.precomputed_epitope = None
        self.nmer = 9

        self.pairs = None
        
        self.msa = None
        self.mut_constraint = None
        self.structure = None
        self.rotamer_choice = None

        self.mut_load = None
        self.num_curves = None

        self.library_gen = False
        self.library_size = None

        self.score_scale = 1.0

        self.designs = None # <- final result file

        self.episweep_path = os.environ["EPISWEEP"]
        self.episweep_data = os.path.join(self.episweep_path, "data")

        self.sweep_direction = -1 # Decrease episcore (For EpiSweep to remove T-cell epitope score)
        # For humanization and disruption, it will be set to be 1
		
        self.three_to_one_letter = {"GLY": "G","ALA": "A","VAL": "V","LEU": "L",
                "ILE": "I","MET": "M","PHE": "F","TRP": "W","PRO": "P","SER": "S",
                "THR": "T","CYS": "C","TYR": "Y","ASN": "N","GLN": "Q","ASP": "D",
                "GLU": "E","LYS": "K","ARG": "R","HIS": "H"}
    def onebody(self, i, c):
        key = (i, c.aa, c.idx)
        if key not in self.onebody_score:
            return 0
        else:
            return self.onebody_score[key]
        
    def twobody(self, i, c, j, d):
        key = (i, c.aa, c.idx, j, d.aa, d.idx)
        if key not in self.twobody_score:
            return 0
        else:
            return self.twobody_score[key]

    def episcore(self, i, t, thresh):
        all_tubes = []
        norm = 1.0
        for c in t:
            all_tubes.append(list(set(c.aas)))
            norm = norm * len(set(c.aas))
        possible_peptides = list(itertools.product(*all_tubes))
        score = 0.0
        if self.precomputed_epitope is not None:
            pep_epi_dict = {}
            for pep_epi in self.precomputed_peptide_list[i]: pep_epi_dict[pep_epi[0]]=pep_epi[1]
            for p in range(len(possible_peptides)):
                peptide = (''.join(c for c in possible_peptides[p]))
                peptide_episcore = pep_epi_dict[peptide]
                score += peptide_episcore/norm
        else:
            for p in range(len(possible_peptides)):
                peptide = (''.join(c for c in possible_peptides[p]))
                peptide_episcore = self.profiles.peptide_hits(peptide, thresh)
                score += peptide_episcore/norm
        return score
        
    def load_spec(self, f):
        spec = open(f, "r").readlines()
        spec = map(lambda x: map(lambda y: y.strip(), x.split(",")), spec)
        spec = dict(spec)
        spec_keywords = spec.keys()
        
        essential_keywords = ["design",
                "wt_seq",
                "epitope_score",
                "epitope_threshold",
                "pairs",
                "design_output"]
        design_type = ["sequence", "rotamer", "tube", "disruption", "humanization"]
        self.design_type = spec["design"]
        if spec["design"] == "humanization":
            essential_keywords.remove("epitope_score")
            essential_keywords.remove("epitope_threshold")
        
        ################
        # Error condition
        if set(spec_keywords) & set(essential_keywords) != set(essential_keywords):
            missing = list(set(essential_keywords) - (set(spec_keywords) & set(essential_keywords)))
            sys.stderr.write("Check the spec file. %s\n"%f)
            sys.stderr.write("Missing keywords: %s\n"%", ".join(missing))
            sys.exit()

        if spec["design"] not in design_type:
            sys.stderr.write("The design type you specified is %s\n"%spec["design"])
            sys.stderr.write("The design type must be one of these: %s\n"%", ".join(design_type))
            sys.exit()
        
        if not self.design_type == "humanization":
            if not os.path.isfile(os.path.join(self.episweep_data, spec["epitope_score"])):
                sys.stderr.write("The epitope score file %s does not exist in %s.\n"%(spec["epitope_score"], self.episweep_data))
                sys.stderr.write("Check your keyword: epitope_score.\n")
                sys.exit()

            self.episcore_mtx = os.path.join(self.episweep_path, "data", spec["epitope_score"])
            self.episcore_threshold = int(spec["epitope_threshold"])
        self.target = spec["wt_seq"]
        if "structure" in spec: self.structure = spec["structure"]
        if "msa" in spec: self.msa = spec["msa"]
        if "choices" in spec: self.possible_choices = spec["choices"]

        if "mut_constraint" in spec:
            if os.path.isfile(spec["mut_constraint"]):
                self.mut_constraint = spec["mut_constraint"]

        if "pairs" in spec: self.pairs = spec["pairs"]
        if "num_curves" in spec: self.num_curves = int(spec["num_curves"])

        # Handling 'mut_load' with potential range indicated by "-"
        if "mut_load" in spec:
            if spec["mut_load"].count("-"):
                self.mut_load = list(map(int, spec["mut_load"].split("-")))
            else:
                self.mut_load = int(spec["mut_load"])

        if "design_output" in spec: self.designs = spec["design_output"]
        if "precomputed_epitope" in spec: self.precomputed_epitope = spec["precomputed_epitope"]

        if spec["design"] == "sequence":
            if "seq_potential" not in spec:
                sys.stderr.write("You have not specified sequence potential file name.\n")
                sys.stderr.write("Check your keyword: seq_potential\n")
                sys.exit()
            else:
                self.score_potential = spec["seq_potential"]
        elif spec["design"] == "rotamer":
            if "rot_potential" not in spec:
                sys.stderr.write("You have not specified rotamer score potential file name.\n")
                sys.stderr.write("Check your keyword: rot_potential\n")
                sys.exit()
            else:
                self.score_potential = spec["rot_potential"]
        elif spec["design"] == "tube":
            if not ("tube_potential" in spec and "seq_potential" in spec):
                sys.stderr.write("If your choice is \"tube\", you have to specify both sequence and tube potential file names.\n")
                sys.stderr.write("Check your keyword: seq_potential and tube_potential\n")
                sys.exit()
            elif "library_size" not in spec:
                sys.stderr.write("If your choice is \"tube\", you have to specify library size.\n")
                sys.stderr.write("Check your keyword: library_size\n")
                sys.exit()
            else:
                self.library_gen = True
                self.score_potential = spec["seq_potential"]
                self.tube_score_potential = spec["tube_potential"]
                self.library_size = list(map(int, spec["library_size"].split("-")))
                self.CLEVER_flag = False

        elif spec["design"] == "humanization":
            if not ("rot_potential" in spec or "seq_potential" in spec):  
                sys.stderr.write("You have not specified sequence or rotamer score potential file name.\n")
                sys.stderr.write("Check your keyword: rot_potential or seq_potential\n")
                sys.exit()
            elif not ("precomputed_HSC" in spec or "precomputed_epitope" in spec):  
                sys.stderr.write("You must provide a precomputed human string content score file.\n")
                sys.stderr.write("Check your keyword: precomputed_HSC or precomputed_epitope\n")
                sys.exit()
            else:
                self.sweep_direction = 1
                if "precomputed_HSC" in spec:  
                    self.precomputed_epitope = spec["precomputed_HSC"]
                if "rot_potential" in spec:  
                    self.score_potential = spec["rot_potential"]
                else:
                    self.score_potential = spec["seq_potential"]
        elif spec["design"] == "disruption":
            if not ("contacts" in spec and "rot_potential" in spec): 
                sys.stderr.write("You have not specified rotamer score potential file name, or the contact file.\n")
                sys.stderr.write("Check your keyword: rot_potential, contacts\n")
                sys.exit()
            else:
                self.episcore_mtx = os.path.join(self.episweep_path, "data", spec["epitope_score"])
                self.score_potential = spec["rot_potential"]
                self.interface_contacts = spec["contacts"]
                self.score_scale = 10.0  # The default value is mentioned for clarity
                self.disruption_score_matrix = {}
                self.sweep_direction = 1
                for i in csv.reader(open(self.episcore_mtx)):
                    res1 = self.three_to_one_letter[i[0]]
                    res2 = self.three_to_one_letter[i[1]]
                    if res1 not in self.disruption_score_matrix: 
                        self.disruption_score_matrix[res1] = {res2: float(i[2])}
                    else:
                        self.disruption_score_matrix[res1][res2] = float(i[2]) 