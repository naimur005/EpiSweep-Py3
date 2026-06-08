Instruction to EpiScope

Document written by Yoonjoo Choi (7 Sep 2017)


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

1. Run a docking program and generate docking models

2. Run DisruPPI against all the docking models.

3. Place results files and contact files in the same folder. Create a results
- contact pair file. You need the structure file too.

Example ("design_contact_pair.csv")

2vxt.00.3.csv,2vxt.00.3_contacts.csv
2vxt.01.3.csv,2vxt.01.3_contacts.csv
2vxt.02.3.csv,2vxt.02.3_contacts.csv
...

4. Run EpiScope

$ episcope.py -p 2vxt_ag.pdb -d 12 -c design_contact_pair.csv -n 3 -o 2vxt.episcope.3.csv
