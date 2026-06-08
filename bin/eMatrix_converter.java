import java.io.*;

public class eMatrix_converter {
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		int validOption = 0;
		String ospreyEnergyMatrixFile = null;
		String convertedEnergyMatrixFile = null;
		String ospreyRefEnergyMatrix = null;
		for (int i = 0; i < args.length; i++) {
			switch (args[i].charAt(0)) {
				case '-':
					if (args[i].length() != 2) {
						System.err.println("Not a valid argument: " + args[i]);
						System.exit(0);
					} else {
						if (args[i].equals("-i")) {
							validOption += 1;
							ospreyEnergyMatrixFile = args[i + 1];
						} else if (args[i].equals("-r")) {
							validOption += 1;
							ospreyRefEnergyMatrix = args[i + 1];
						} else if (args[i].equals("-o")) {
							validOption += 1;
							convertedEnergyMatrixFile = args[i + 1];
						}
					}
			}
		}
		if (validOption != 3) {
			System.err.println("OSPREY Energy matrix converter.\n");
			System.err.println("Three arguments are required:");
			System.err.println("\t-i: OSPREY energy matrix file");
			System.err.println("\t-r: OSPREY reference energy matrix file");
			System.err.println("\t-o: Converted output binary file");
			System.exit(0);
		}

		if (!(new File(ospreyEnergyMatrixFile).isFile())) {
			System.err.println("OSPREY energy matrix file \"" + ospreyEnergyMatrixFile + "\" does not exist.");
			System.exit(0);
		}

		if (!(new File(ospreyRefEnergyMatrix).isFile())) {
			System.err.println("OSPREY reference energy file \"" + ospreyRefEnergyMatrix + "\" does not exist.");
			System.exit(0);
		}

		String[] intToAA;
		System.out.println("using OSPREY 2.0 AAs");
		String[] o2 = { "ALA", "VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "CYS", "MET",
				"SER", "THR", "LYS", "ARG", "HIP", "HID", "HIE", "ASP", "GLU",
				"ASN", "GLN", "GLY", "PRO" };
		/*
		 * String[] o2 = { "A", "V", "L", "I", "F", "Y", "W", "C", "M",
		 * "S", "T", "K", "R", "H", "H", "H", "D", "E", "N", "Q", "G", "P" };
		 */
		intToAA = o2;

		// Load OSPREY energy matrix
		System.out.println("Loading the energy Matrix...");
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(ospreyEnergyMatrixFile));
		float energyMatrix[][][][][][] = (float[][][][][][]) in.readObject();
		in.close();

		// Open converted energy matrix output file
		DataOutputStream out = new DataOutputStream(new FileOutputStream(convertedEnergyMatrixFile));

		// CBK: why -1?
		int npos = energyMatrix.length - 1;
		System.out.println("There are " + npos + " positions in the energy matrix.");
		System.out.println();

		// Start output file with number of positions
		out.writeInt(npos);

		// First part of output file gives for each position the number of AAs
		// and a list of the AA indices and their number of rotamers
		for (int p = 0; p < npos; p++) {
			System.out.print(p + 1 + " ");
			byte naas = 0;
			for (byte a = 0; a < energyMatrix[p].length; a++) {
				if (energyMatrix[p][a] != null)
					naas++;
			}
			out.writeInt(naas);

			for (byte a = 0; a < energyMatrix[p].length; a++) {
				if (energyMatrix[p][a] != null) {
					byte nrots = 0;
					for (byte r = 0; r < energyMatrix[p][a].length; r++) {
						if (energyMatrix[p][a][r] != null)
							nrots++;
					}
					System.out.print(intToAA[a] + ":" + nrots + " ");
					out.writeByte(a);
					out.writeInt(nrots);
				}
			}
			System.out.println();
		}

		// Calculate Reference energy
		System.out.println("Loading the reference energy Matrix...");
		try (ObjectInputStream RefIn = new ObjectInputStream(new FileInputStream(ospreyRefEnergyMatrix))) {
			float RefEMatrix[][] = (float[][]) RefIn.readObject();
			in.close();
			float[] minEaa = new float[RefEMatrix[0].length];
			for (int p = 0; p < npos; p++) {
				for (int a = 0; a < RefEMatrix[p].length; a++) {
					float thisE = RefEMatrix[p][a];
					if (thisE < minEaa[a]) {
						minEaa[a] = thisE;
					}
				}
			}
			for (int p = 0; p < npos; p++) {
				for (int a = 0; a < energyMatrix[p].length; a++) {
					if (energyMatrix[p][a] != null) {
						for (int r = 0; r < energyMatrix[p][a].length; r++) {
							if (energyMatrix[p][a][r] != null) {
								// energyMatrix[p][a][r][p][0][1] -= minEaa[a];
								energyMatrix[p][a][r][p][0][1] -= RefEMatrix[p][a];
							}
						}
					}
				}
			}
		}
		// Second part of output file gives the one-body energy of each (position, AA,
		// rotamer)
		// This includes the rotamer's internal energy [p][0][0] and that with
		// non-flexible template [p][0][1]
		// CBK: check new OSPREY documentation -- this is based on Andrew's code from
		// earlier version
		System.out.println("One-body terms...");
		for (int p = 0; p < npos; p++) {
			// float percent = (float) (p + 1) / npos * 100;
			// System.out.printf("%2.1f%%", percent);

			for (byte a = 0; a < energyMatrix[p].length; a++) {
				if (energyMatrix[p][a] != null) {
					for (byte r = 0; r < energyMatrix[p][a].length; r++) {
						if (energyMatrix[p][a][r] != null) {
							// System.out.print("p= "+p+" a= "+intToAA[a]+" r= "+r+" ");
							// System.out.println(energyMatrix[p][a][r][p][0][1]);
							out.writeFloat(energyMatrix[p][a][r][p][0][1]);
						}
					}
				}
			}
			// System.out.print("\b\b\b\b\b\b");
		}

		// Third part of output file gives the two-body energy of each
		// (position, AA, rotamer, position', AA', rotamer')
		System.out.println("Two-body terms...");
		for (int p = 0; p < npos - 1; p++) {
			float percent = (float) (p + 1) / npos * 100;
			System.out.printf("%2.1f%%", percent);
			for (byte a = 0; a < energyMatrix[p].length; a++) {
				if (energyMatrix[p][a] != null) {
					for (byte r = 0; r < energyMatrix[p][a].length; r++) {
						if (energyMatrix[p][a][r] != null) {
							for (int q = p + 1; q < npos; q++) {
								if (energyMatrix[p][a][r][q] != null) {
									for (byte b = 0; b < energyMatrix[p][a][r][q].length; b++) {
										if (energyMatrix[p][a][r][q][b] != null) {
											for (byte s = 0; s < energyMatrix[p][a][r][q][b].length; s++) {
												out.writeFloat(energyMatrix[p][a][r][q][b][s]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
			System.out.print("\b\b\b\b\b\b");
		}

		System.out.println("Done!");
		out.close();
	}
}
