#!/usr/bin/python3

import sys
from episweep import *

if len(sys.argv) < 2:
    sys.stderr.write("EpiSweep - Deimmunization algorithm\n\n")
    sys.stderr.write("Design spec file is required.\n\n")
    sys.stderr.write("$ %s [Design spec csv file]\n"%sys.argv[0])
    sys.exit()

design = Design() # Initialiser: individual variant generator
design.load_spec(sys.argv[1])

# Optimise
design.initialize_optimizer()
design.optimize()

