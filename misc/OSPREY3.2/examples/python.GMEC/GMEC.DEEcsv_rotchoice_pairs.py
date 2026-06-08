import osprey
import csv
import os

osprey.start()

# Input PDB file
pdb_file = 'lst_cwbd.pdb'

# Extract the base name (without extension) for dynamic file naming
base_name = os.path.splitext(os.path.basename(pdb_file))[0]

# Load the PDB file
strand = osprey.Strand(pdb_file)

# Configure flexibility of certain residues
strand.flexibility['21'].setLibraryRotamers('ALA', 'GLY', 'SER')
strand.flexibility['42'].setLibraryRotamers('VAL', 'LEU')

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

# Generate energy matrix CSV
energy_matrix_file = f"{base_name}_energy_matrix.csv"
with open(energy_matrix_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "AA", "Rotamer", "Single Energy", "Pair Energy"])
    for position in confSpace.positions:
        posIndex = position.index
        numRotamers = len(position.resConfs)
        for rc1 in range(numRotamers):
            aa1 = position.resConfs[rc1].template.name
            single_energy = emat.getOneBody(posIndex, rc1)
            writer.writerow([posIndex, aa1, rc1, single_energy, ""])
            for pos2Index in range(posIndex + 1, len(confSpace.positions)):
                numRotamers2 = len(confSpace.positions[pos2Index].resConfs)
                for rc2 in range(numRotamers2):
                    try:
                        aa2 = confSpace.positions[pos2Index].resConfs[rc2].template.name
                        pair_energy = emat.getPairwise(posIndex, rc1, pos2Index, rc2)
                        writer.writerow([posIndex, aa1, rc1, pos2Index, aa2, rc2, "", pair_energy])
                    except Exception:
                        writer.writerow([posIndex, aa1, rc1, pos2Index, aa2, rc2, "", "N/A"])

print(f"Energy matrix saved to '{energy_matrix_file}'")

# Generate rotamer choices file
rot_choices_file = f"{base_name}_rot_choices.csv"
with open(rot_choices_file, 'w', newline='') as choicesfile:
    writer = csv.writer(choicesfile)
    writer.writerow(["Position", "AA", "Rotamer"])
    for position in confSpace.positions:
        posIndex = position.index + 1  # Convert to 1-based index
        numRotamers = len(position.resConfs)
        for rcIndex in range(numRotamers):
            aa = position.resConfs[rcIndex].template.name
            writer.writerow([posIndex, aa, rcIndex])

print(f"Rotamer choices saved to '{rot_choices_file}'")

# Generate residue pairs file
pairs_file = f"{base_name}_pairs.csv"
with open(pairs_file, 'w', newline='') as pairsfile:
    writer = csv.writer(pairsfile)
    writer.writerow(["Position1", "Position2"])
    for position in confSpace.positions:
        posIndex = position.index + 1  # Convert to 1-based index
        for pos2Index in range(posIndex + 1, len(confSpace.positions) + 1):
            writer.writerow([posIndex, pos2Index])

print(f"Residue pairs saved to '{pairs_file}'")

# Run DEE with steric pruning
pmat = osprey.DEE(confSpace, emat, showProgress=True)

# Find the best sequence and rotamers using A* search
astar = osprey.AStarMPLP(emat, pmat)
gmec = osprey.GMECFinder(astar, confEcalc).find()

# Save GMEC results to CSV
gmec_results_file = f"{base_name}_gmec_results.csv"
with open(gmec_results_file, 'w', newline='') as gmecfile:
    writer = csv.writer(gmecfile)
    writer.writerow(["Rotamer numbers", "L", "L2", "L1", "Energy", "Score"])
    for conf in gmec:
        rotamer_numbers = conf[0]
        l = conf[1]
        l2 = conf[2]
        l1 = conf[3]
        energy = conf[4]
        score = conf[5]
        writer.writerow([rotamer_numbers, l, l2, l1, energy, score])

print(f"GMEC results saved to '{gmec_results_file}'")
