#!/usr/bin/python3

from optparse import OptionParser
import sys, re

parser = OptionParser()
parser.add_option('-i', '--input', dest='input', help="Clever log filename.", action='store_true')
parser.add_option('-o', '--output', dest='output', help="EpiSweep sequence potential filename.", action='store_true')
(options, args) = parser.parse_args()

mandatory = ['input','output']
for m in mandatory:
    if not options.__dict__[m]:
        print ("\nMissing options\n")
        parser.print_help()
        sys.exit()

if len(args) != 2:
    print ("\nMissing filenames\n")
    parser.print_help()
    sys.exit()

WT_reference = 0.0
onebody_aascores = {}
twobody_aascores = {}
wt = ''

file = args[0]
# Not a good way to read file, but works for reading Clever files.
# Read WT sequence first
with open(file) as input_data:
    # Skips text before the beginning of the interesting block:
    for line in input_data:
        if line.rstrip() == ' Num nDSeq   List of amino acids':  # Or whatever test is needed
            break
    # Reads text until the end of the block:
    for line in input_data:  # This keeps reading the file
        if line.rstrip() == 'Num: arbitrary index of clusters':
            break
        
        temp = re.split(' ',line.rstrip())  # Line is extracted (or block_of_lines.append(line), etc.)
        chars = filter(None, temp)
        if len(chars) >= 3: wt = wt + chars[2]


# Not a good way to read file, but works for reading Clever files.
with open(file) as input_data:
    # Skips text before the beginning of the interesting block:
    for line in input_data:
        if line.rstrip() == '  Num          ECI  NH-ratio     1/IPR CFInfo':  # Or whatever test is needed
            break
    # Reads text until the end of the block:
    for line in input_data:  # This keeps reading the file
        if line.rstrip() == '### Definition and other information ###':
            break
        
        temp = re.split(' ',line.rstrip())  # Line is extracted (or block_of_lines.append(line), etc.)
        scores = filter(None, temp)
        if len(scores) == 7 and scores[6] == 'CF': WT_reference = float(scores[1])
        elif len(scores) == 7: onebody_aascores[int(scores[6])-1,scores[5]] = float(scores[1])
        elif len(scores) == 9: twobody_aascores[int(scores[7])-1,scores[5],int(scores[8])-1,scores[6]] = float(scores[1])


for i in range(len(wt)):
    onebody_aascores[i,wt[i]] = 0.0

outputFH = open(args[1], 'w')

for key in sorted(onebody_aascores):
    outputFH.write(str(key[0]+1)+','+key[1]+','+str(onebody_aascores[key])+',,\n')

for key in sorted(twobody_aascores):
    outputFH.write(str(key[0]+1)+','+key[1]+','+str(key[2]+1)+','+key[3]+','+str(twobody_aascores[key])+'\n')

outputFH.close()




