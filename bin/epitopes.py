#!/usr/bin/python3

import re,sys,optparse,csv
from episweep import *

parser = optparse.OptionParser(description='Epitope predictor')
parser.add_option('-i', '--design_spec',
        dest="design_spec",
        action="store",
        help='Input: Design spec file')
parser.add_option('-n', '--index',
        dest="index", action="store",
        default="0",
        help='Design index (Default: 0; Wild type)')
parser.add_option('-o', '--output', 
        dest="output", 
        action="store", 
        help='Output: Epitope prediction results')

args, opts = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit()

design = Design()
design.load_spec(args.design_spec)
design.profiles = PPProfiles(design.episcore_mtx)

all_designs = map(lambda x: x.split(","), open(design.designs).readlines())
design_of_interest = list(filter(lambda x: x[0] == args.index, all_designs))[0] 

seq = design_of_interest[-1].strip("\n")

epitope_prediction = design.profiles.full_report(seq, design.episcore_threshold)


with open(args.output, "w", newline='') as file_out: 
    output = csv.writer(file_out)
    header = ["Start", "Peptide", "End"] + design.profiles.allele_list + ["Total"]
    output.writerow(header)  
    print(",".join(header))  
    
    for p, pred in enumerate(epitope_prediction):
        pred = list(map(lambda x: int(bool(x)), pred))  # Convert the map object to a list
        pred.append(sum(pred))
        line = [p+1, seq[p:p+9], p+9] + pred
        output.writerow(line)  # Writes each line to the CSV file
        print(",".join(map(str, line)))  