Instruction to how to perform CoDAH humanization using EpiSweep

Document written by Yoonjoo Choi (3 Sep 2017)


=====
References 

Choi, Yoonjoo, Deeptak Verma, Karl E. Griswold, and Chris Bailey-Kellogg.
"EpiSweep: Computationally driven reengineering of therapeutic proteins to
reduce immunogenicity while maintaining function." Computational Protein
Design (2017): 375-398.

Choi, Yoonjoo, Christian Ndong, Karl E. Griswold, and Chris Bailey-Kellogg.
"Computationally driven antibody engineering enables simultaneous humanization
and thermostabilization." Protein Engineering, Design and Selection 29, no. 10
(2016): 419-426.

Choi, Yoonjoo, Casey Hua, Charles L. Sentman, Margaret E. Ackerman, and Chris
Bailey-Kellogg. "Antibody humanization by structure-based computational
protein design." MAbs, vol. 7, no. 6 (2015): 1045-1057.
=====

1. Annotate light and heavy chain sequences using ANACRI and save the results in CSV.

http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php

Currently CoDAH supports Chothia and Kabat numbering schemes.

See examples:

cetuximab_VH_kabat.csv
cetuximab_VK_kabat.csv


2. Calculate the wild type HSC score.

$ HSC_calculator.py -L cetuximab_VK_kabat.csv -H cetuximab_VH_kabat.csv 
The numbering scheme is Kabat

Heavy chain HSC score: 60.707 (N9: 1)
Light chain HSC score: 72.054 (N9: 7)
Total HSC score: 66.082 (N9: 8)

N9 indicates the number of exact nonamers found in human germline sequences.


3. The target antibody structure must match the annotated sequences.

In order to calculate an OSPREY-stryle energy matrix, the structure must be
single-chained. The chain order must be [Light] first and [Heavy] second.

Ex) cetuximab.pdb


4. Select mutational choices.

$ HSC_calculator.py -L cetuximab_VK_kabat.csv -H cetuximab_VH_kabat.csv -p
cetuximab.pdb -n 3 -s Kabat -o cetuximab_precomputed_HSC.csv

Options)
L: annotated light chain sequence (using ANACRI)
H: annotated heavy chain sequence
p: structure for the OSPREY energy calculation
n: number of best matched human germline sequences to include for mutational choices
o: precomputed HSC score output file name


5. Calculate the energy matrix. See an example script (osprey.bash)


6. Run EpiSweep

$ sweep.py design_spec.csv


7. Analyze and recalculate the design output

$ HSC_calculator.py -L cetuximab_VK_kabat.csv -H cetuximab_VH_kabat.csv -D
cetuximab_rot_designs.csv 

