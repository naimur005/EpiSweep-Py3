import optparse

def episweep_option(parser):
    parser.add_option('-t', '--target',
            action="store",
            dest="target_sequence",
            default=None,
            help="Target sequence (One letter code)")

    preprocessor_group = optparse.OptionGroup(parser,
            "Preprocessor",
            "These are options for preprocessor")
    
    preprocessor_group.add_option('-m', '--msa', 
            action="store", 
            dest="msa", 
            default=None, 
            help="Multiple sequence alignment file")
    preprocessor_group.add_option('-s', '--secondary_structure', 
            action="store", 
            dest="secondary_structure", 
            default=None, 
            help="[Optional] Secondary structure sequence (H: helix, E: strand, T: turn or loop) ex) TTTEETTHHHTTTT")
    preprocessor_group.add_option('-w', '--sequence_weight', 
            action="store", 
            dest="sequence_weight", 
            default=None, 
            help="[Optional] Sequence weight file. The number of entries must be the same as the number of sequence in the msa file input. A weight value ranges from 0 to 1.")
    preprocessor_group.add_option('-e', '--episcore_mtx', 
            action="store", 
            dest="episcore_mtx", 
            default=None, 
            help="[Optional] Epitope score matrix")
    preprocessor_group.add_option('-a', '--allowed_AA', 
            action="store", 
            dest="allowed_AA", 
            default=None, 
            help="[Optional] Output file for allowed amino acids. If omitted, results will be shown only on screen.")
    preprocessor_group.add_option('-p', '--pair_info', 
            action="store", 
            dest="pair_info_output", 
            default=None, 
            help="[Optional] Output file for pair information. If omitted, results will be shown only on screen.")

    parser.add_option_group(preprocessor_group)

# This function now returns the parsed arguments, which can be used directly in our script.
