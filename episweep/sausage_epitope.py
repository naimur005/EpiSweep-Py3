#from episweep import Design
from .epitopes import *
from .__init__ import *
import csv, sys, os
try:
    from pymol import *
except: pass

def sausage(spec=None, palette="white red", overlap=True, obj=None):
    __author__ = "Yoonjoo Choi, Miles Campbell, Chris Bailey-Kellogg"
    __copyright__ = "Copyleft 2015, Epi-Sausage project"

    __version__ = "0.2"
    __maintainer__ = "Yoonjoo Choi"
    __email__ = "yoonjoo@cs.dartmouth.edu"
    AA3to1 = {"ALA": "A", "GLU": "E", "GLN": "Q", "ASP": "D", 
            "ASN": "N", "LEU": "L", "GLY": "G", "LYS": "K", 
            "SER": "S", "VAL": "V", "ARG": "R", "THR": "T", 
            "PRO": "P", "ILE": "I", "MET": "M", "PHE": "F", 
            "TYR": "Y", "CYS": "C", "CYX": "C", "CYD": "C", 
            "TRP": "W", "HIS": "H", "HIP": "H", "HIE": "H", "HID": "H"}

    if spec is None:
        if obj is None:
            print ("Please specify a design spec file.")
            return False
    design = Design()
    design.load_spec(spec)
    if obj is None:
        obj = os.path.splitext(os.path.basename(design.structure))[0]

    if obj not in cmd.get_names():
        print ("Perhaps your structure is not in current PyMOL objects.")
        return False

    design.profiles = PPProfiles(design.episcore_mtx)
    protein = list(filter(lambda x: x[:4] == "ATOM", open(design.structure).readlines()))
    protein = list(filter(lambda x: x[13:15] == "CA", protein))

    episcore = {}
    resnums = {}
    seqs = {}
    max_episcore = 0.0
    
    #######################
    # Parsing PDB structure
    for c in protein:
        chain = c[21]
        res_num = c[22:26].strip()
        aa = AA3to1[c[17:20]]
        if not seqs.has_key(chain):
            resnums[chain] = [res_num]
            seqs[chain] = aa
            episcore[chain] = {res_num: 0}
        else:
            resnums[chain].append(res_num)
            seqs[chain] += aa
            episcore[chain][res_num] = 0

    for chain in seqs:
        for p in range(len(seqs[chain])-8):
            peptide = seqs[chain][p:p+9]
            res = resnums[chain][p]
            peptide_episcore = design.profiles.num_hits(peptide, design.episcore_threshold)
            if overlap:
                for q in resnums[chain][p:p+9]:
                    episcore[chain][q] += peptide_episcore
                    if episcore[chain][q] > max_episcore:
                        max_episcore = episcore[chain][q]
            else:
                episcore[chain][res] += peptide_episcore
                if episcore[chain][res] > max_episcore:
                    max_episcore = episcore[chain][res]

    #######################
    # Setting up colour palette
    color_vector = list(map(lambda x: cmd.get_color_tuple(x), palette.split()))
    scale = 60
    steps = scale / (len(color_vector) - 1)

    #######################
    # Epitope score setting
    cmd.alter(obj, 'b=0.0')
    for chain in episcore:
        for res in episcore[chain]:
            epi = str(episcore[chain][res])
            selected = obj+ " and n. CA and chain %s"%(chain)
            cmd.alter(selected + " and i. %s"%res, 'b='+epi)
            col_start = 0.0
            for p in range(len(color_vector)-1):
                for i in range(steps):
                    s = float(i) / steps
                    color_list = [color_vector[p + 1][j] * s + color_vector[p][j] * (1.0 - s) for j in range(3)]
                    color_name = '0x%02x%02x%02x' %(tuple(map(lambda x: x*255, color_list)))
                    cmd.color(color_name, "(%s) and b > %f"%(selected, col_start))
                    col_start = float(max_episcore) * (i+1+p*steps)/ scale

    
    ################
    # Other settings
    cmd.show('cartoon')
    cmd.cartoon('putty')
    cmd.set('cartoon_putty_radius', 0.2)
    cmd.set('ray_trace_mode', 0)

    cmd.color('white', obj+' and b=0')
    cmd.hide('lines')
    cmd.bg_color("white")

    