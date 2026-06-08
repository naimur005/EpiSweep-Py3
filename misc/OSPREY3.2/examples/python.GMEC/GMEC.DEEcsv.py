import osprey
import csv

osprey.start()

# Load the PDB file
strand = osprey.Strand('lst_cwbd.pdb')

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

# Open the CSV file for writing
with open('results.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "AA", "Rotamer", "Single Energy", "Pair Energy"])

    # Print and write to CSV
    print("Single and pair energy matrix's...")
    for position in confSpace.positions:
        posIndex = position.index
        numRotamers = len(position.resConfs)  # ResConfs holds rotamer or similar configurations
        for rc1 in range(numRotamers):
            aa1 = position.resConfs[rc1].template.name  # Get the amino acid type from the template
            single_energy = emat.getOneBody(posIndex, rc1)
            print(f"Single energy at position {posIndex}, AA {aa1}, rotamer {rc1}: {single_energy}")
            writer.writerow([posIndex, aa1, rc1, single_energy, ""])
            for pos2Index in range(posIndex + 1, len(confSpace.positions)):
                numRotamers2 = len(confSpace.positions[pos2Index].resConfs)
                for rc2 in range(numRotamers2):
                    try:
                        aa2 = confSpace.positions[pos2Index].resConfs[rc2].template.name  # Get the amino acid type for the rotamer
                        pair_energy = emat.getPairwise(posIndex, rc1, pos2Index, rc2)
                        print(f"Pair energy between pos {posIndex}, AA {aa1}, rot {rc1} and pos {pos2Index}, AA {aa2}, rot {rc2}: {pair_energy}")
                        writer.writerow([posIndex, aa1, rc1, pos2Index, aa2, rc2, "", pair_energy])
                    except Exception as e:
                        print(f"No pairwise energy available between pos {posIndex}, AA {aa1}, rot {rc1} and pos {pos2Index}, AA {aa2}, rot {rc2}")
                        writer.writerow([posIndex, aa1, rc1, pos2Index, aa2, rc2, "", "N/A"])

# Run DEE with just steric pruning, adjust thresholds here if set explicitly
pmat = osprey.DEE(confSpace, emat, showProgress=True)

# Find the best sequence and rotamers using A* search
astar = osprey.AStarMPLP(emat, pmat)
gmec = osprey.GMECFinder(astar, confEcalc).find()

# Save GMEC result to CSV
with open('gmec_results.csv', 'w', newline='') as gmecfile:
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

print("Results have been saved to 'results.csv' and 'gmec_results.csv'")
