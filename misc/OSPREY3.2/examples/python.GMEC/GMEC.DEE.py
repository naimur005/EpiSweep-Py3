import osprey
import csv

osprey.start()

# Load the PDB file
strand = osprey.Strand('lst_cwbd_5.pdb')

# Configure flexibility of certain residues
#strand.flexibility['21'].setLibraryRotamers('ALA', 'GLY', 'SER')
#strand.flexibility['42'].setLibraryRotamers('VAL', 'LEU')

strand.flexibility['1'].setLibraryRotamers('LYS')
strand.flexibility['2'].setLibraryRotamers('THR')
# Create the configuration space
confSpace = osprey.ConfSpace(strand)

# Choose a forcefield
ffparams = osprey.ForcefieldParams()

# Setup the energy calculator
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# Define how energies of conformations are computed
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# Compute the energy matrix
emat = osprey.EnergyMatrix(confEcalc)

# Define the mapping from amino acid names to integer codes
aa_to_int = {
    "ALA": 0, "VAL": 1, "LEU": 2, "ILE": 3, "PHE": 4, "TYR": 5, "TRP": 6, "CYS": 7, "MET": 8,
    "SER": 9, "THR": 10, "LYS": 11, "ARG": 12, "HIP": 13, "HID": 14, "HIE": 15, "ASP": 16,
    "GLU": 17, "ASN": 18, "GLN": 19, "GLY": 20, "PRO": 21
}

with open('results.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    writer.writerow(['Number of Positions', len(confSpace.positions)])

    one_body_terms = []
    two_body_terms = []

    for position in confSpace.positions:
        posIndex = position.index
        numRotamers = len(position.resConfs)
        aas = set(position.resConfs[rc].template.name for rc in range(numRotamers))
        writer.writerow(['Position', posIndex + 1, 'Number of AAs', len(aas)])

        for aa in aas:
            aa_code = aa_to_int[aa]  # Get the integer code for the amino acid
            aa_rotamers = [rc for rc in range(numRotamers) if position.resConfs[rc].template.name == aa]
            writer.writerow(['AA', aa, 'Number of Rotamers', len(aa_rotamers)])

            for rc1 in aa_rotamers:
                single_energy = emat.getOneBody(posIndex, rc1)
                one_body_terms.append([posIndex, aa_code, rc1, single_energy])
                for pos2Index in range(posIndex + 1, len(confSpace.positions)):
                    numRotamers2 = len(confSpace.positions[pos2Index].resConfs)
                    for rc2 in range(numRotamers2):
                        aa2 = confSpace.positions[pos2Index].resConfs[rc2].template.name
                        aa2_code = aa_to_int[aa2]  # Get the integer code for the second amino acid
                        try:
                            pair_energy = emat.getPairwise(posIndex, rc1, pos2Index, rc2)
                            two_body_terms.append([posIndex, aa_code, rc1, pos2Index, aa2_code, rc2, pair_energy])
                        except Exception as e:
                            pass

    # one-body terms
    writer.writerow(['One-body terms'])
    for term in one_body_terms:
        posIndex, aa_code, rc1, single_energy = term
        writer.writerow([posIndex, aa_code, rc1, single_energy])

    # two-body terms
    writer.writerow(['Two-body terms'])
    for term in two_body_terms:
        posIndex, aa1_code, rc1, pos2Index, aa2_code, rc2, pair_energy = term
        writer.writerow([posIndex, aa1_code, rc1, pos2Index, aa2_code, rc2, pair_energy])

print("Results have been saved to 'results.csv'")

# run DEE with just steric pruning
pmat = osprey.DEE(confSpace, emat, showProgress=True)

# or run DEE with Goldstein pruning
#i0 = 10.0 # kcal/mol
#pmat = osprey.DEE(confSpace, emat, showProgress=True, singlesGoldsteinDiffThreshold=i0, pairsGoldsteinDiffThreshold=i0)

# Find the best sequence and rotamers using A* search
astar = osprey.AStarMPLP(emat, pmat)
gmec = osprey.GMECFinder(astar, confEcalc).find()

#######################################################################################################################
# Load the PDB file
strand = osprey.Strand('lst_cwbd_5.pdb')

# Define flexibility for heavy and light chains by specifying residues and chain identifiers
# Example: Heavy chain residues
strand.flexibility['H:21'].setLibraryRotamers('ALA', 'GLY', 'SER')  # Residue 21 on chain H
strand.flexibility['H:42'].setLibraryRotamers('VAL', 'LEU')

# Example: Light chain residues
strand.flexibility['L:1'].setLibraryRotamers('LYS')
strand.flexibility['L:2'].setLibraryRotamers('THR')
strand.flexibility['L:3'].setLibraryRotamers('ASN')

# Create the configuration space
confSpace = osprey.ConfSpace(strand)

# Choose a forcefield
ffparams = osprey.ForcefieldParams()

# Setup the energy calculator
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
####################################################################################################################
# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE)
