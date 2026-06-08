#!/usr/bin/python3

import sys
from episweep import *

design = Design() # Initialiser
design.load_spec("design_spec.csv")

design.initialize_filter()
design.msa_filter(low_thresh=0.30, high_thresh=0.95, fraction_gap_allowed=0.25)
design.background_frequency_filter(input='mccaldon.csv')
design.epitope_score_filter(min_epi_del=1)
design.apply_mutational_constraints()

design.structure_pair_calculator(cutoff=8)
design.OSPREY_config_maker(target="lst_cwbd")
