import sys
import os
import csv
import javaobj

def deserialize_java_object(file_path):
    with open(file_path, 'rb') as f:
        obj = javaobj.load(f)
    return obj

def main(args):
    valid_option = 0
    osprey_energy_matrix_file = None
    converted_energy_matrix_file = None
    osprey_ref_energy_matrix = None

    i = 0
    while i < len(args):
        if args[i][0] == '-':
            if len(args[i]) != 2:
                print(f"Not a valid argument: {args[i]}", file=sys.stderr)
                sys.exit(0)
            else:
                if args[i] == "-i":
                    valid_option += 1
                    osprey_energy_matrix_file = args[i + 1]
                elif args[i] == "-r":
                    valid_option += 1
                    osprey_ref_energy_matrix = args[i + 1]
                elif args[i] == "-o":
                    valid_option += 1
                    converted_energy_matrix_file = args[i + 1]
        i += 1

    if valid_option != 3:
        print("OSPREY Energy matrix converter.\n", file=sys.stderr)
        print("Three arguments are required:", file=sys.stderr)
        print("\t-i: OSPREY energy matrix file", file=sys.stderr)
        print("\t-r: OSPREY reference energy matrix file", file=sys.stderr)
        print("\t-o: Converted output CSV file", file=sys.stderr)
        sys.exit(0)

    if not os.path.isfile(osprey_energy_matrix_file):
        print(f'OSPREY energy matrix file "{osprey_energy_matrix_file}" does not exist.', file=sys.stderr)
        sys.exit(0)

    if not os.path.isfile(osprey_ref_energy_matrix):
        print(f'OSPREY reference energy file "{osprey_ref_energy_matrix}" does not exist.', file=sys.stderr)
        sys.exit(0)

    print("using OSPREY 2.0 AAs")
    int_to_aa = [
        "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET",
        "SER", "THR", "LYS", "ARG", "HIP", "HID", "HIE", "ASP", "GLU",
        "ASN", "GLN", "GLY", "PRO"
    ]

    # Deserialize the Java objects
    print("Deserializing the energy matrix...")
    energy_matrix = deserialize_java_object(osprey_energy_matrix_file)
    
    print("Deserializing the reference energy matrix...")
    ref_energy_matrix = deserialize_java_object(osprey_ref_energy_matrix)

    # Open converted energy matrix output file
    with open(converted_energy_matrix_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write the number of positions as the first row
        npos = len(energy_matrix) - 1
        print(f"There are {npos} positions in the energy matrix.")
        csv_writer.writerow(['Number of Positions', npos])

        # First part of output file gives for each position the number of AAs
        # and a list of the AA indices and their number of rotamers
        for p in range(npos):
            print(f"{p + 1} ", end="")
            naas = 0
            for a in range(len(energy_matrix[p])):
                if energy_matrix[p][a] is not None:
                    naas += 1
            csv_writer.writerow(['Position', p + 1, 'Number of AAs', naas])

            for a in range(len(energy_matrix[p])):
                if energy_matrix[p][a] is not None:
                    nrots = 0
                    for r in range(len(energy_matrix[p][a])):
                        if energy_matrix[p][a][r] is not None:
                            nrots += 1
                    print(f"{int_to_aa[a]}:{nrots} ", end="")
                    csv_writer.writerow(['AA', int_to_aa[a], 'Number of Rotamers', nrots])
            print()

        # Calculate Reference energy
        min_eaa = [float('inf')] * len(ref_energy_matrix[0])
        for p in range(npos):
            for a in range(len(ref_energy_matrix[p])):
                this_e = ref_energy_matrix[p][a]
                if this_e < min_eaa[a]:
                    min_eaa[a] = this_e

        for p in range(npos):
            for a in range(len(energy_matrix[p])):
                if energy_matrix[p][a] is not None:
                    for r in range(len(energy_matrix[p][a])):
                        if energy_matrix[p][a][r] is not None:
                            energy_matrix[p][a][r][p][0][1] -= ref_energy_matrix[p][a]

        # Second part of output file gives the one-body energy of each (position, AA, rotamer)
        # This includes the rotamer's internal energy [p][0][0] and that with
        # non-flexible template [p][0][1]
        print("One-body terms...")
        csv_writer.writerow(['One-body terms'])
        for p in range(npos):
            for a in range(len(energy_matrix[p])):
                if energy_matrix[p][a] is not None:
                    for r in range(len(energy_matrix[p][a])):
                        if energy_matrix[p][a][r] is not None:
                            csv_writer.writerow([p, a, r, energy_matrix[p][a][r][p][0][1]])

        # Third part of output file gives the two-body energy of each
        # (position, AA, rotamer, position', AA', rotamer')
        print("Two-body terms...")
        csv_writer.writerow(['Two-body terms'])
        for p in range(npos - 1):
            percent = (float(p + 1) / npos) * 100
            print(f"{percent:.1f}%", end="\r")
            for a in range(len(energy_matrix[p])):
                if energy_matrix[p][a] is not None:
                    for r in range(len(energy_matrix[p][a])):
                        if energy_matrix[p][a][r] is not None:
                            for q in range(p + 1, npos):
                                if energy_matrix[p][a][r][q] is not None:
                                    for b in range(len(energy_matrix[p][a][r][q])):
                                        if energy_matrix[p][a][r][q][b] is not None:
                                            for s in range(len(energy_matrix[p][a][r][q][b])):
                                                csv_writer.writerow([p, a, r, q, b, s, energy_matrix[p][a][r][q][b][s]])

        print("Done!")

if __name__ == "__main__":
    main(sys.argv[1:])
