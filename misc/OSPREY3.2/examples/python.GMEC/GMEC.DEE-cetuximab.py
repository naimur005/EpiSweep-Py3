import osprey
import csv
import os

osprey.start()

pdb_file = 'cetuximab.pdb'

# Extract the base name (without extension) for dynamic file naming
base_name = os.path.splitext(os.path.basename(pdb_file))[0]

# Load the PDB file
strand = osprey.Strand(pdb_file)

# Configure flexibility of certain residues
strand.flexibility['1'].setLibraryRotamers('ASP', 'GLU')
strand.flexibility['2'].setLibraryRotamers('ILE')
strand.flexibility['3'].setLibraryRotamers('LEU', 'VAL')
strand.flexibility['4'].setLibraryRotamers('LEU', 'MET')
strand.flexibility['5'].setLibraryRotamers('THR')
strand.flexibility['6'].setLibraryRotamers('GLN')
strand.flexibility['7'].setLibraryRotamers('SER')
strand.flexibility['8'].setLibraryRotamers('PRO')
strand.flexibility['9'].setLibraryRotamers('VAL', 'ASP', 'ALA')
strand.flexibility['10'].setLibraryRotamers('ILE', 'PHE', 'THR')
strand.flexibility['11'].setLibraryRotamers('LEU', 'GLN')
strand.flexibility['12'].setLibraryRotamers('SER')
strand.flexibility['13'].setLibraryRotamers('VAL')
strand.flexibility['14'].setLibraryRotamers('SER', 'THR')
strand.flexibility['15'].setLibraryRotamers('PRO')
strand.flexibility['16'].setLibraryRotamers('GLY', 'LYS')
strand.flexibility['17'].setLibraryRotamers('GLU')
strand.flexibility['18'].setLibraryRotamers('ARG', 'LYS')
strand.flexibility['19'].setLibraryRotamers('VAL', 'ALA')
strand.flexibility['20'].setLibraryRotamers('SER', 'THR')
strand.flexibility['21'].setLibraryRotamers('PHE', 'ILE', 'LEU')
strand.flexibility['22'].setLibraryRotamers('SER', 'THR')
strand.flexibility['23'].setLibraryRotamers('CYS')
strand.flexibility['24'].setLibraryRotamers('ARG')
strand.flexibility['25'].setLibraryRotamers('ALA')
strand.flexibility['26'].setLibraryRotamers('SER')
strand.flexibility['27'].setLibraryRotamers('GLN')
strand.flexibility['28'].setLibraryRotamers('SER')
strand.flexibility['29'].setLibraryRotamers('ILE')
strand.flexibility['30'].setLibraryRotamers('GLY')
strand.flexibility['31'].setLibraryRotamers('THR')
strand.flexibility['32'].setLibraryRotamers('ASN')
strand.flexibility['33'].setLibraryRotamers('ILE')
strand.flexibility['34'].setLibraryRotamers('HIP')
strand.flexibility['35'].setLibraryRotamers('TRP')
strand.flexibility['36'].setLibraryRotamers('TYR')
strand.flexibility['37'].setLibraryRotamers('GLN')
strand.flexibility['38'].setLibraryRotamers('GLN')
strand.flexibility['39'].setLibraryRotamers('ARG', 'LYS')
strand.flexibility['40'].setLibraryRotamers('THR', 'PRO')
strand.flexibility['41'].setLibraryRotamers('ASN', 'ASP', 'GLY')
strand.flexibility['42'].setLibraryRotamers('GLY', 'GLN')
strand.flexibility['43'].setLibraryRotamers('SER', 'ALA')
strand.flexibility['44'].setLibraryRotamers('PRO')
strand.flexibility['45'].setLibraryRotamers('ARG', 'LYS')
strand.flexibility['46'].setLibraryRotamers('LEU')
strand.flexibility['47'].setLibraryRotamers('LEU')
strand.flexibility['48'].setLibraryRotamers('ILE')
strand.flexibility['49'].setLibraryRotamers('LYS', 'TYR')
strand.flexibility['50'].setLibraryRotamers('TYR')
strand.flexibility['51'].setLibraryRotamers('ALA')
strand.flexibility['52'].setLibraryRotamers('SER')
strand.flexibility['53'].setLibraryRotamers('GLU')
strand.flexibility['54'].setLibraryRotamers('SER')
strand.flexibility['55'].setLibraryRotamers('ILE')
strand.flexibility['56'].setLibraryRotamers('SER')
strand.flexibility['57'].setLibraryRotamers('GLY')
strand.flexibility['58'].setLibraryRotamers('ILE', 'VAL')
strand.flexibility['59'].setLibraryRotamers('PRO')
strand.flexibility['60'].setLibraryRotamers('SER', 'ALA')
strand.flexibility['61'].setLibraryRotamers('ARG')
strand.flexibility['62'].setLibraryRotamers('PHE')
strand.flexibility['63'].setLibraryRotamers('SER')
strand.flexibility['64'].setLibraryRotamers('GLY')
strand.flexibility['65'].setLibraryRotamers('SER')
strand.flexibility['66'].setLibraryRotamers('GLY')
strand.flexibility['67'].setLibraryRotamers('SER')
strand.flexibility['68'].setLibraryRotamers('GLY')
strand.flexibility['69'].setLibraryRotamers('THR')
strand.flexibility['70'].setLibraryRotamers('ASP', 'GLU')
strand.flexibility['71'].setLibraryRotamers('PHE')
strand.flexibility['72'].setLibraryRotamers('THR')
strand.flexibility['73'].setLibraryRotamers('LEU')
strand.flexibility['74'].setLibraryRotamers('SER', 'THR')
strand.flexibility['75'].setLibraryRotamers('ILE')
strand.flexibility['76'].setLibraryRotamers('ASN', 'SER')
strand.flexibility['77'].setLibraryRotamers('SER')
strand.flexibility['78'].setLibraryRotamers('VAL', 'LEU')
strand.flexibility['79'].setLibraryRotamers('GLU', 'GLN')
strand.flexibility['80'].setLibraryRotamers('SER', 'ALA')
strand.flexibility['81'].setLibraryRotamers('GLU')
strand.flexibility['82'].setLibraryRotamers('ASP')
strand.flexibility['83'].setLibraryRotamers('ILE', 'ALA', 'PHE')
strand.flexibility['84'].setLibraryRotamers('ALA')
strand.flexibility['85'].setLibraryRotamers('ASP', 'THR', 'VAL')
strand.flexibility['86'].setLibraryRotamers('TYR')
strand.flexibility['87'].setLibraryRotamers('TYR')
strand.flexibility['88'].setLibraryRotamers('CYS')
strand.flexibility['89'].setLibraryRotamers('GLN')
strand.flexibility['90'].setLibraryRotamers('GLN')
strand.flexibility['91'].setLibraryRotamers('ASN')
strand.flexibility['92'].setLibraryRotamers('ASN')
strand.flexibility['93'].setLibraryRotamers('ASN')
strand.flexibility['94'].setLibraryRotamers('TRP')
strand.flexibility['95'].setLibraryRotamers('PRO')
strand.flexibility['96'].setLibraryRotamers('THR')
strand.flexibility['97'].setLibraryRotamers('THR')
strand.flexibility['98'].setLibraryRotamers('PHE')
strand.flexibility['99'].setLibraryRotamers('GLY')
strand.flexibility['100'].setLibraryRotamers('ALA')
strand.flexibility['101'].setLibraryRotamers('GLY')
strand.flexibility['102'].setLibraryRotamers('THR')
strand.flexibility['103'].setLibraryRotamers('LYS')
strand.flexibility['104'].setLibraryRotamers('LEU')
strand.flexibility['105'].setLibraryRotamers('GLU')
strand.flexibility['106'].setLibraryRotamers('LEU')
strand.flexibility['107'].setLibraryRotamers('LYS')
strand.flexibility['108'].setLibraryRotamers('GLN')
strand.flexibility['109'].setLibraryRotamers('VAL')
strand.flexibility['110'].setLibraryRotamers('GLN', 'THR')
strand.flexibility['111'].setLibraryRotamers('LEU')
strand.flexibility['112'].setLibraryRotamers('LYS', 'GLN')
strand.flexibility['113'].setLibraryRotamers('GLN', 'GLU')
strand.flexibility['114'].setLibraryRotamers('SER')
strand.flexibility['115'].setLibraryRotamers('GLY')
strand.flexibility['116'].setLibraryRotamers('PRO')
strand.flexibility['117'].setLibraryRotamers('GLY', 'VAL')
strand.flexibility['118'].setLibraryRotamers('LEU')
strand.flexibility['119'].setLibraryRotamers('VAL')
strand.flexibility['120'].setLibraryRotamers('GLN', 'LYS')
strand.flexibility['121'].setLibraryRotamers('PRO')
strand.flexibility['122'].setLibraryRotamers('SER', 'THR')
strand.flexibility['123'].setLibraryRotamers('GLN', 'GLU')
strand.flexibility['124'].setLibraryRotamers('SER', 'THR')
strand.flexibility['125'].setLibraryRotamers('LEU')
strand.flexibility['126'].setLibraryRotamers('SER', 'THR')
strand.flexibility['127'].setLibraryRotamers('ILE', 'LEU')
strand.flexibility['128'].setLibraryRotamers('THR')
strand.flexibility['129'].setLibraryRotamers('CYS')
strand.flexibility['130'].setLibraryRotamers('THR')
strand.flexibility['131'].setLibraryRotamers('VAL')
strand.flexibility['132'].setLibraryRotamers('SER')
strand.flexibility['133'].setLibraryRotamers('GLY')
strand.flexibility['134'].setLibraryRotamers('PHE', 'GLY')
strand.flexibility['135'].setLibraryRotamers('SER')
strand.flexibility['136'].setLibraryRotamers('LEU', 'ILE')
strand.flexibility['137'].setLibraryRotamers('THR', 'SER')
strand.flexibility['138'].setLibraryRotamers('ASN')
strand.flexibility['139'].setLibraryRotamers('TYR')
strand.flexibility['140'].setLibraryRotamers('GLY')
strand.flexibility['141'].setLibraryRotamers('VAL')
strand.flexibility['142'].setLibraryRotamers('HIP')
strand.flexibility['143'].setLibraryRotamers('TRP')
strand.flexibility['144'].setLibraryRotamers('VAL', 'ILE')
strand.flexibility['145'].setLibraryRotamers('ARG')
strand.flexibility['146'].setLibraryRotamers('GLN')
strand.flexibility['147'].setLibraryRotamers('SER', 'HIP', 'PRO')
strand.flexibility['148'].setLibraryRotamers('PRO')
strand.flexibility['149'].setLibraryRotamers('GLY')
strand.flexibility['150'].setLibraryRotamers('LYS')
strand.flexibility['151'].setLibraryRotamers('GLY', 'ALA')
strand.flexibility['152'].setLibraryRotamers('LEU')
strand.flexibility['153'].setLibraryRotamers('GLU')
strand.flexibility['154'].setLibraryRotamers('TRP')
strand.flexibility['155'].setLibraryRotamers('LEU', 'ILE')
strand.flexibility['156'].setLibraryRotamers('GLY', 'ALA')
strand.flexibility['157'].setLibraryRotamers('VAL')
strand.flexibility['158'].setLibraryRotamers('ILE')
strand.flexibility['159'].setLibraryRotamers('TRP')
strand.flexibility['160'].setLibraryRotamers('SER')
strand.flexibility['161'].setLibraryRotamers('GLY')
strand.flexibility['162'].setLibraryRotamers('GLY')
strand.flexibility['163'].setLibraryRotamers('ASN')
strand.flexibility['164'].setLibraryRotamers('THR')
strand.flexibility['165'].setLibraryRotamers('ASP')
strand.flexibility['166'].setLibraryRotamers('TYR')
strand.flexibility['167'].setLibraryRotamers('ASN')
strand.flexibility['168'].setLibraryRotamers('THR')
strand.flexibility['169'].setLibraryRotamers('PRO')
strand.flexibility['170'].setLibraryRotamers('PHE')
strand.flexibility['171'].setLibraryRotamers('THR')
strand.flexibility['172'].setLibraryRotamers('SER')
strand.flexibility['173'].setLibraryRotamers('ARG')
strand.flexibility['174'].setLibraryRotamers('LEU', 'VAL')
strand.flexibility['175'].setLibraryRotamers('SER', 'THR')
strand.flexibility['176'].setLibraryRotamers('ILE')
strand.flexibility['177'].setLibraryRotamers('ASN', 'SER')
strand.flexibility['178'].setLibraryRotamers('LYS', 'VAL')
strand.flexibility['179'].setLibraryRotamers('ASP')
strand.flexibility['180'].setLibraryRotamers('ASN', 'THR')
strand.flexibility['181'].setLibraryRotamers('SER')
strand.flexibility['182'].setLibraryRotamers('LYS')
strand.flexibility['183'].setLibraryRotamers('SER', 'ASN')
strand.flexibility['184'].setLibraryRotamers('GLN')
strand.flexibility['185'].setLibraryRotamers('VAL', 'PHE')
strand.flexibility['186'].setLibraryRotamers('PHE', 'SER', 'VAL')
strand.flexibility['187'].setLibraryRotamers('PHE', 'LEU')
strand.flexibility['188'].setLibraryRotamers('LYS', 'THR')
strand.flexibility['189'].setLibraryRotamers('MET', 'LEU')
strand.flexibility['190'].setLibraryRotamers('ASN', 'SER', 'THR')
strand.flexibility['191'].setLibraryRotamers('SER', 'ASN')
strand.flexibility['192'].setLibraryRotamers('LEU', 'VAL', 'MET')
strand.flexibility['193'].setLibraryRotamers('GLN', 'THR', 'ASP')
strand.flexibility['194'].setLibraryRotamers('SER', 'ALA', 'PRO')
strand.flexibility['195'].setLibraryRotamers('ASN', 'ALA', 'VAL')
strand.flexibility['196'].setLibraryRotamers('ASP')
strand.flexibility['197'].setLibraryRotamers('THR')
strand.flexibility['198'].setLibraryRotamers('ALA')
strand.flexibility['199'].setLibraryRotamers('ILE', 'VAL', 'THR')
strand.flexibility['200'].setLibraryRotamers('TYR')
strand.flexibility['201'].setLibraryRotamers('TYR')
strand.flexibility['202'].setLibraryRotamers('CYS')
strand.flexibility['203'].setLibraryRotamers('ALA')
strand.flexibility['204'].setLibraryRotamers('ARG')
strand.flexibility['205'].setLibraryRotamers('ALA')
strand.flexibility['206'].setLibraryRotamers('LEU')
strand.flexibility['207'].setLibraryRotamers('THR')
strand.flexibility['208'].setLibraryRotamers('TYR')
strand.flexibility['209'].setLibraryRotamers('TYR')
strand.flexibility['210'].setLibraryRotamers('ASP')
strand.flexibility['211'].setLibraryRotamers('TYR')
strand.flexibility['212'].setLibraryRotamers('GLU')
strand.flexibility['213'].setLibraryRotamers('PHE')
strand.flexibility['214'].setLibraryRotamers('ALA')
strand.flexibility['215'].setLibraryRotamers('TYR')
strand.flexibility['216'].setLibraryRotamers('TRP')
strand.flexibility['217'].setLibraryRotamers('GLY')
strand.flexibility['218'].setLibraryRotamers('GLN')
strand.flexibility['219'].setLibraryRotamers('GLY')
strand.flexibility['220'].setLibraryRotamers('THR')
strand.flexibility['221'].setLibraryRotamers('LEU')
strand.flexibility['222'].setLibraryRotamers('VAL')
strand.flexibility['223'].setLibraryRotamers('THR')
strand.flexibility['224'].setLibraryRotamers('VAL')
strand.flexibility['225'].setLibraryRotamers('SER')

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
