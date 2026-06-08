#!/usr/bin/python3

import re,sys,argparse,csv

parser = argparse.ArgumentParser(description='Rotamer information extractor.')
parser.add_argument('-i', '--log', dest="log", action="store", help='Input: OSPREY energy calculation log file')
parser.add_argument('-o', '--rot', dest="rot", action="store", help='Output: Rotamer input file (rotamer-based EpiSweep)')

args = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit()

osprey_AA = ["A","V","L","I","F","Y","W","C","M","S","T","K","R","H","H","H","D","E","N","Q","G","P"]
accepted_rotamers = "("+re.findall(r"\((.*?)curIndexRed", open(args.log).read())[0]
rotamers = map(lambda x: x.split(), re.findall(r"\((.*?)\)", accepted_rotamers))
rotamers = map(lambda x: [int(x[0])+1, osprey_AA[int(x[1])], int(x[2])], rotamers) # position start from 0 in OSPREY2.1 and so move back to 1

with open(args.rot, "wb") as rot_file:
    rot_csv_file = csv.writer(rot_file)
    map(lambda row: rot_csv_file.writerow(row), rotamers)