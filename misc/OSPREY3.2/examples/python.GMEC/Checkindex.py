import osprey

osprey.start()

# Load the PDB file
strand = osprey.Strand('lst_cwbd.pdb')

# Create the configuration space
conf_space = osprey.ConfSpace(strand)

# Check available attributes in conf_space
print("Attributes of conf_space:", dir(conf_space))

# Check if positions is an attribute of conf_space
if hasattr(conf_space, 'positions'):
    positions = conf_space.positions
    
    # Check if positions is an iterable and its length
    if isinstance(positions, (list, tuple)):
        print(f"Positions is a list or tuple with {len(positions)} elements.")
    else:
        print("Positions is not a list or tuple. Type:", type(positions))
    
    # If positions is iterable, iterate and print details
    try:
        print("Positions in conf_space:")
        for i, pos in enumerate(positions):
            print(f"Position {i}: {pos}")
            # Check if position has a residue number or other identifiable attribute
            try:
                res_number = pos.getResidueNumber()
                print(f"  Residue Number: {res_number}")
            except AttributeError:
                print("  No residue number attribute")
            
            # List all attributes of the position object
            print("  Attributes of position:", dir(pos))

    except Exception as e:
        print(f"Error iterating over positions: {e}")

else:
    print("No 'positions' attribute in conf_space.")

# We don't need to stop osprey as there's no stop method
# osprey.stop()
