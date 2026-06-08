#!/usr/bin/python3

import re,sys,optparse,csv,random
from episweep import *

parser = optparse.OptionParser(description='Epitope predictor')
parser.add_option('-i', '--design_spec',
        dest="design_spec",
        action="store",
        help='Input: Design spec file')
parser.add_option('-n', '--index',
        dest="index", action="store",
        help='Design index (Default: 0; Wild type)')
parser.add_option('-e', '--number',
        dest="number", action="store",
        default=False,
        type=int,
        help='The number of random samples of library elements (Default: False; retrieve all)')
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
design.library_gen = False
design.initialize_optimizer()
design.all_tubes_dict = {}
design.tube_dict_generator()

all_designs = map(lambda x: x.split(","), open(design.designs).readlines())
wt_seq = map(lambda x: [x], design.wt_seq)
seq_len = len(design.wt_seq)
design_of_interest = next(filter(lambda x: x[0] == args.index, all_designs), None) 
design.load_csv_score_potential()

mutation = design_of_interest[-2].split("|")

for mut in mutation:
    mut_code = ""
    pos = ""
    for m in mut[1:]:
        try: int(m); pos += m
        except: mut_code += m
    p = int(pos) - 1
    if not mut_code.isupper():
        wt_seq[p] = design.all_tubes_dict[mut_code]
    else: wt_seq[p] = list(set(mut_code))

enumerated_seq = ["".join(x) for x in itertools.product(*wt_seq)] 
whole_num = len(enumerated_seq)
if args.number:
    random.shuffle(enumerated_seq)
    enumerated_seq = random.sample(enumerated_seq, args.number)
selected_num = len(enumerated_seq)

output = csv.writer(open(args.output, "wb"))
output.writerow(map(lambda x: x.strip("\n"), all_designs[0])) 
output.writerow([0,0]+map(lambda x: x.strip("\n"), design_of_interest)[2:]) 

for n, seq in enumerate(enumerated_seq):
    mut_list = []
    for p, a in enumerate(design.wt_seq):
        if a != seq[p]: mut_list.append("%s%d%s"%(a, p+1, seq[p]))
    episcore = design.profiles.num_hits(seq, design.episcore_threshold)
    score = 0.0
    for p1 in range(seq_len-1):
        a1 = seq[p1]
        key1 = (p1, a1, 0)
        if not design.onebody_score.has_key(key1): score += design.penalty
        else: score += design.onebody_score[key1]
        for p2 in range(p1+1, seq_len):
            a2 = seq[p2]
            key2 = (p2, a2, 0)
            pair_key = (p1, a1, 0, p2, a2, 0)
            if ((p1, p2) in design.possible_pairs) or ((p2, p1) in design.possible_pairs):
                if design.twobody_score.has_key(pair_key): score += design.twobody_score[pair_key]
                else:
                    if design.onebody_score.has_key(key1) and design.onebody_score.has_key(key2): score += 0
                    else: score += 2*design.penalty
    row = [n+1, 1, episcore, score, "|".join(mut_list), seq]
    output.writerow(row)
    print (" ").join(map(str, ["%d out of %d -"%(n+1, selected_num), "Episcore:", episcore, "Potential score:", score, mut_list]))

sys.stdout.write("\nSummary:\n\n")
sys.stdout.write("%d random samples (Total: %d) were retrieved.\n\n"%(selected_num, whole_num))
sys.stdout.write("The results are saved in \"%s\".\n"%args.output)