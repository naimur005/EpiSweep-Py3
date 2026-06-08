#!/usr/bin/python3

import numpy as np
import os, sys, random, csv
from glob import glob
from episweep.k_medoids import Kmedoids
from episweep.episcope_score import Filter
import argparse

# Setup argument parser
parser = argparse.ArgumentParser(description='EpiScope')
parser.add_argument('-p', '--pdb', action='store', help='PDB structure')
parser.add_argument('-c', '--contact', action='store', help='Contact file [csv]')
parser.add_argument('-d', '--distance', action='store', type=int, default=12, help='Average distance cut-off (Default: 12)')
parser.add_argument('-n', '--cluster', action='store', type=int, help='Number of representatives')
parser.add_argument('-o', '--output', action='store', help='Output file')
args = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit()

# Ensure the PDB file exists
if not os.path.isfile(args.pdb):
    print(f"Error: PDB file '{args.pdb}' not found.")
    sys.exit(1)

# Read structure from PDB file
struct = [line for line in open(args.pdb).readlines() if line[:4] == "ATOM" and line[13:15] == "CA"]
struct = {line[22:26].strip(): list(map(float, [line[30:38], line[38:46], line[46:54]])) for line in struct}

# Load the scoring matrix from the expected environment variable path
score_mtx_path = os.path.join(os.getenv("EPISWEEP", "."), "data", "INT5.csv")
if not os.path.isfile(score_mtx_path):
    print(f"Error: Scoring matrix file '{score_mtx_path}' not found.")
    sys.exit(1)

# Initialize the Filter class with the provided arguments
a = Filter(struct, args.contact, score_mtx_path, args.distance)

# Perform clustering
print(f"Trying {args.cluster} cluster")
km = Kmedoids(k=args.cluster)
if not a.repr_crd:
    print("No data available for clustering. Exiting.")
    sys.exit(1)
km.cluster(np.array(a.repr_crd))

# Generate output
a.cluster_coverage(km.medoid, km.oth_elm, km.clst)

print("\nFinal selection:\n")
with open(args.output, "w", newline='') as output_file:
    output = csv.writer(output_file)
    header = ["design", "design_pos", "disruption_score", "number_of_dock", "dock_in_contact_list", "disruption_score_list"]
    output.writerow(header)
    for x in a.selected[0][0]:
        print(x.design)
        line = [
            x.design, 
            "|".join(x.design_pos),  # Assuming design_pos is a list of strings
            x.episcore, 
            len(x.in_contact), 
            "|".join(x.in_contact), 
            "|".join(map(str, x.average_score_list))
        ]
        output.writerow(line)

print("\nOutput file: %s" % args.output)
