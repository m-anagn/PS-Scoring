
"""

Description:
    LCD-Composer is a composition-based method for identifying low-complexity domains in protein sequences.

    Refer to Cascarina et. al. (2021) NAR Genomics & Bioinformatics (https://academic.oup.com/nargab/article/3/2/lqab048/6285187) 
    for a complete description of algorithm testing and application.

    Documentation for LCD-Composer is available at https://github.com/RossLabCSU/LCD-Composer.

============================================================================================

License_info:
    LCD-Composer is subject to the terms of the GPLv3.0 license. For a complete description of 
    license terms, please see the license at https://github.com/RossLabCSU/LCD-Composer.

"""

__author__ = 'Sean M Cascarina'
__copyright__ = 'Copyright 2020'
__credits__ = ['Sean M Cascarina']
__license__ = 'GPLv3.0'
__version__ = '1.0'
__maintainer__ = 'Sean M Cascarina'
__email__ = 'Sean.Cascarina@colostate.edu'



def main(args):

    all_aas = 'ACDEFGHIKLMNPQRSTVWY'
    win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, ign_disp = get_params(args)
    combined_aas = ''.join(amino_acids)
    start_time = str(datetime.datetime.now())
    print('\nStart time:', str(datetime.datetime.now()), amino_acids, '\n')
    
                
    output = open(args.results_file + '.tsv', 'w')
    output_params_header(output, args, ign_disp)
    
    frequency_output = open(args.results_file + '_BasicFrequencyData.tsv', 'w')
    frequency_output.write('\t'.join( ['File of Origin', 'Total # of Proteins in Proteome', 'Total # of LCDs Idenitified', '# of Proteins Containing at Least One LCD', 'Percentage of Proteins with LCD', 'ProtIDs with LCD (ampersand-delimited)'] ) + '\n')

    if args.verbose:
        output.write('\t'.join( ['File of Origin', 'Protein ID', '\t'.join([aa_group + '_Composition' for aa_group in amino_acids]), 'Combined Linear Dispersion'] ) + '\n')
    else:
        output.write('\t'.join( ['File of Origin', 'Protein ID','Domain Sequence','Domain Boundaries','Final Domain Composition','Final Domain Linear Dispersion','\t'.join([aa for aa in all_aas]) ]) + '\n')

    mins_df = {}
    maxs_df = {}
    
    completed = 0
    num_fastas = get_num_fastas(args)
    for (dirname, dir, files) in os.walk('.'):
        if not args.scan_subdirectories:
            dir.clear()
        for file in files:
            if args.file_keyword and args.file_keyword not in file:
                continue
            if file.endswith('.fasta') or file.endswith('.fsa') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fnn') or file.endswith('.faa') or file.endswith('.frn'):
                filename = os.path.join(dirname, file)
 
                h = open(filename)
                
                total_seqs = 0
                lcd_count = 0
                prots_with_lcd = set()

                for id, seq in fasta_parser(h):

                    #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
                    if seq.count('*') == 1 and seq[-1] == '*':
                        seq = seq[:-1]
                    elif seq.count('*') > 1:
                        if seq[-1] == '*':
                            seq = seq[:-1]

                    normed_stddevs = []
                    group_comps = {aa_group:[] for aa_group in amino_acids}
                    hit_positions = []
                    pos = -1
                    for i in range(len(seq) - win_size+1):
                        window = seq[i:i+win_size]
                            
                        pos += 1

                        #CALCULATE COMPOSITION OF EACH aa_group FOR THE WINDOW
                        comps, group_comps = calc_composition(window, amino_acids, group_comps)
                        
                        #MASK REGIONS THAT DON'T PASS COMPOSITION THRESHOLD
                        if 0 in [1 if comps[i] > comp_thresholds[i] else 0 for i in range(len(amino_acids))] and args.verbose==False:
                            normed_stddevs.append( -1 )
                            continue

                        #PERFORM FULL CALCULATIONS
                        else:
                            #QUICK CALCULATION FOR WINDOWS THAT ARE 100% AAs OF INTEREST
                            if sum([window.count(aa)/win_size*100 for aa in combined_aas]) > 100-0.000001:
                                normed_stddevs.append(1)
                                hit_positions.append(pos)
                                continue
                                
                            norm_stddev = calc_dispersion(window, combined_aas, other_aas, mins_df, maxs_df)
                            normed_stddevs.append( norm_stddev )

                            if sum(comps) >= ign_disp or norm_stddev > disp_threshold:
                                hit_positions.append( pos )

                    #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
                    domain_boundaries = merge_windows(hit_positions, win_size, id, seq)
                    seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

                    #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
                    trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)
                    
                    #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
                    trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
                    final_comps, final_stddevs = calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas)
                    
                    total_seqs += 1
                    
                    if len(trimmed_boundaries) > 0 and args.verbose==False:
                        lcd_count += len(trimmed_seqs)
                        prot_id, *desc = id.split(' ')
                        prots_with_lcd.add(prot_id)
                        for i in range(len(trimmed_seqs)):
                            output.write('\t'.join( [file, id, trimmed_seqs[i], trimmed_boundaries[i], str(round(final_comps[i], 2)), str(round(final_stddevs[i], 4)), '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
                    elif len(trimmed_boundaries) > 0 and args.verbose==True:
                        lcd_count += len(trimmed_seqs)
                        prot_id, *desc = id.split(' ')
                        prots_with_lcd.add(prot_id)
                        output.write('\t'.join( [file, id, '\t'.join( ['_'.join([str(round(x, 2)) for x in group_comps[aa_group]]) for aa_group in amino_acids] ), '_'.join([str(round(x, 4)) for x in normed_stddevs]) ] ) + '\n')
                    else:
                        continue
                        
                frequency_output.write('\t'.join( [file, str(total_seqs), str(lcd_count), str(len(prots_with_lcd)), str( len(prots_with_lcd) / total_seqs * 100 ), '&'.join(prots_with_lcd)] ) + '\n')
                completed += 1
                print(str(completed) + ' proteomes completed out of ' + str(num_fastas) + ' total proteomes\t', str(datetime.datetime.now()))
        
                h.close()
    output.close()
    
    print('\nStart time:', start_time)
    print('End time:', str(datetime.datetime.now()), '\n')
    
    
def get_num_fastas(args):
    """
    Determines the number of FASTA files that will be analyzed.
    
    Returns:
        Number of FASTA files (int)
    """
    num_fastas = 0
    for (dirname, dir, files) in os.walk('.'):
        if not args.scan_subdirectories:
            dir.clear()
        if args.file_keyword:
            num_fastas += sum( [1 if ( file.endswith('.fasta') or file.endswith('.fsa') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fnn') or file.endswith('.faa') or file.endswith('.frn') ) and args.file_keyword in file else 0 for file in files] )
        else:
            num_fastas += sum( [1 if ( file.endswith('.fasta') or file.endswith('.fsa') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fnn') or file.endswith('.faa') or file.endswith('.frn') ) else 0 for file in files] )
            
    return num_fastas

def calc_composition(window, amino_acids, group_comps):
    """
    Calculates the amino acid composition for each user-defined amino acid group in amino_acids.
    
    Returns:
        Compositions (list)
        Group compositions (dictionary)
    """
    comps = []
    for aa_group in amino_acids:
        comp = 0
        for aa in aa_group:
            comp += window.count(aa) / len(window) * 100
        comps.append( comp )
        group_comps[aa_group].append( comp )
        
    return comps, group_comps
    

def calc_dispersion(window, combined_aas, other_aas, mins_df, maxs_df):
    """
    Calculates the normalized spacing standard deviation for a given window sequence.
    
    Returns:
        Normalized spacing standard deviations (float)
    """
    
    res_spacing_stddev = calc_spacing(window, combined_aas, other_aas)
    num_muts = sum([window.count(aa) for aa in other_aas])

    min_stddev, max_stddev, mins_df, maxs_df = get_MinMax_stddevs(len(window), num_muts, mins_df, maxs_df)
    
    #HANDLES CASES WITH 100% COMPOSITION (min_stddev AND max_stddev ARE BOTH ZERO)
    if min_stddev == max_stddev:
        norm_stddev = 1.0
    else:
        norm_stddev = 1 - (res_spacing_stddev - min_stddev) / (max_stddev - min_stddev)
        
    return norm_stddev
                    

def calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas):
    """
    Calculates the final composition and normalized spacing standard deviation 
    for each merged & trimmed domain identified.
    
    Returns:
        Compositions (list)
        Normalized spacing standard deviations (list)
    """
    
    final_comps = []
    final_stddevs = []
    for seq in trimmed_seqs:
        comp = sum([seq.count(aa) for aa in combined_aas]) / len(seq) * 100
        norm_stddev = calc_dispersion(seq, combined_aas, other_aas, mins_df, maxs_df)
        num_muts = sum([seq.count(aa) for aa in other_aas])
        
        final_comps.append(comp)
        final_stddevs.append(norm_stddev)
        
    return final_comps, final_stddevs
    

def get_MinMax_stddevs(win_size, num_muts, mins_df, maxs_df):
    """
    Checks whether the minimum possible stddev and maximum possible stddev have been 
    calculated for a sequence of length win_size and a composition 
    corresponding to num_muts.
    
    If they have not already been calculated, generates sequences with minimum and 
    maximum standard deviations in spacings and adds them to lookup dictionaries.
    
    Returns:
        min_stddev (float)
        max_stddev (float)
        mins_df (dictionary)
        maxs_df (dictionary)
    """

    mins_df[win_size] = mins_df.get(win_size, {})
    maxs_df[win_size] = maxs_df.get(win_size, {})
    
    if num_muts not in mins_df[win_size]:
        min_seq = generate_min_seq(win_size, num_muts)
        max_seq = generate_max_seq(win_size, num_muts)

        min_stddev = calc_spacing(min_seq, 'Z', 'G')
        max_stddev = calc_spacing(max_seq, 'Z', 'G')

        mins_df[win_size][num_muts] = mins_df[win_size].get(num_muts, min_stddev)
        maxs_df[win_size][num_muts] = maxs_df[win_size].get(num_muts, max_stddev)

    return mins_df[win_size][num_muts], maxs_df[win_size][num_muts], mins_df, maxs_df
                    
                    
def trim_termini(seqs, amino_acids, domain_boundaries):
    """
    Trims domains that pass composition/separation thresholds until the termini match 
    one of the specified amino acids of interest.
    
    Returns:
        Trimmed sequences (list).
        Trimmed domains (list of tuples, where each tuple contains the start and 
                        end indices of the domain, respectively).
    """
    
    for i in range(len(seqs)):
        while seqs[i][0] not in amino_acids:
            seqs[i] = seqs[i][1:]
            domain_boundaries[i] = (domain_boundaries[i][0]+1, domain_boundaries[i][1])
            
        while seqs[i][-1] not in amino_acids:
            seqs[i] = seqs[i][:-1]
            domain_boundaries[i] = (domain_boundaries[i][0], domain_boundaries[i][1]-1)

    return seqs, domain_boundaries
    

def merge_windows(hit_positions, win_size, id, seq):
    """
    Merge neighboring windows that pass the composition/separation thresholds and
    lie within 1/2 the window size distance from each other.
    
    Returns:
        The boundaries (according to Python indexing, not protein positions) of all identified domains.
        This is a list of tuples, where each tuple contains the start and 
        end indices of the domain, respectively.
    """

    j = 0
    domain_boundaries = []
    
    while hit_positions:
        start = hit_positions[0]
        while j < (len(hit_positions) - 1) and ( (hit_positions[j] + win_size) >= ( hit_positions[j+1] ) ):
            j += 1
        
        #GETS THE ENDING POSITION FOR THE CURRENT WINDOW AND STORES THE 2-TUPLE
        end = hit_positions[j]
        domain_boundaries.append( (start , end+win_size) )

        #MODIFIES hit_positions TO DELETE ALL POSITIONS THAT WERE JUST MERGED, THEN RE-SETS j=0 TO START ITERATING AT THE FIRST POSITION IN THE NEW LIST
        hit_positions = hit_positions[j+1:]
        j = 0

    return domain_boundaries
    
        
def calc_spacing(seq, amino_acids, other_aas):
    """
    Calculates the spacings of the either the amino acid(s) of interest or all other amino acids (whichever is smallest).
    
    Returns:
        Standard deviation of the spacings (float)
    """

    #ADDS RESIDUE OF INTEREST TO THE ENDS SO THAT SPACINGS ARE ALSO CALCULATED BETWEEN THE N-TERM AND THE FIRST RESIDUE OF INTEREST, AND BETWEEN THE LAST RESIDUE OF INTEREST AND THE C-TERM
    amino_acids = ''.join(amino_acids)
    seq = amino_acids[0] + seq + amino_acids[0]
    res_pos = [x for x in range(len(seq)) if seq[x] in amino_acids]
    
    seq = other_aas[0] + seq[1:-1] + other_aas[0]
    nonres_pos = [x for x in range(len(seq)) if seq[x] not in amino_acids]
    res_spacings = sorted( [ res_pos[x] - res_pos[x-1] for x in range(1, len(res_pos))] +  [ nonres_pos[x] - nonres_pos[x-1] for x in range(1, len(nonres_pos))])
    
    return calc_stddev(res_spacings)
    
    
def calc_stddev(spacings):
    """
    Calculates the population standard deviation.
    
    Returns:
        Standard deviation (float)
    """
    
    mean = statistics.mean(spacings)
    var = sum( [(x - mean)**2 for x in spacings] ) / len(spacings)
    stddev = math.sqrt(var)
    
    return stddev
    
    
def fasta_parser(file):
    """Parses each instance in a FASTA formatted file into a gene id and a sequence.
    
    Yields:
        id, seq (both strings)
    """
    #INITIALIZES GENE AND SEQ FOR FIRST ITERATION
    gene, seq = '', []
    
    #LOOPING THROUGH FILE USING GENERATOR
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
            
        #YIELDS GENE AND SEQ IF THE NEXT ID IS REACHED
        if line.startswith('>'):
            if gene != '':
                yield (gene[1:], ''.join(seq))
            gene, seq = line, []
        else:
            seq.append(line)
            
    #YIELDS FINAL INSTANCE IN FASTA FILE
    if gene != '':
        yield (gene[1:], ''.join(seq))
                

def generate_min_seq(win_size, num_muts):
    """
    Generate a sequence with Z num_muts and maximum spacing between Z's.
    This will be used to calculate the minimum stddev value for normalization,
    since large spacing results in small stddevs.
    
    Returns:
        min_seq (string)
    """
    
    #CORNER CASE WITH SINGLE MUTATION
    if num_muts == 1:
        min_seq = 'G'*int(win_size/2) + 'Z'
        while len(min_seq) < win_size:
            min_seq += 'G'
            
    #CORNER CASES WHERE NUMBER OF MUTATIONS IS ~1/2 THE win_size
    elif num_muts == int(win_size / 2):
        num_muts = int(win_size / 2) + 1
        min_seq = 'ZG'*num_muts
            
    else:
        if num_muts >= int(win_size/2):
            num_muts = win_size-num_muts
            
        min_seq = ''
        block_size = int( (win_size - num_muts) / (num_muts+1) )
        remainder = (win_size-num_muts) % (num_muts+1)
        
        while remainder > 0:
            min_seq += 'Z'*(block_size+1) + 'G'
            remainder -= 1
            
        while len(min_seq) < win_size:
            min_seq += 'Z'*block_size + 'G'
            
    min_seq = min_seq[:win_size]
    
    return min_seq
    
    
def generate_max_seq(win_size, num_muts):
    """
    Generate a sequence with Z num_muts and all Z's clustered at one end (minimum spacing)
    This will be used to calculate the maximum stddev value for normalization.
    
    Returns:
        max_seq (string)
    """
    
    max_seq = num_muts*'Z' + 'G'*(win_size-num_muts)
    
    return max_seq
    
    
def get_params( args ):
    """Gather and define user parameters. Also generates an error message if an invalid reduced alphabet is passed in by the user.
    
    Returns:
        1) Window size (int)
        2) Amino acids of interest (str)
        3) Amino acids/characters that are not of interest (str)
        4) Composition threshold (float)
        5) Linear dispersion threshold (float)
        6) Ignore dispersion threshold (float)
    """
    
    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
    comp_thresholds = args.composition.split('_')
    comp_thresholds = [float(x)-0.000001 for x in comp_thresholds]
    
    if not args.ignore_dispersion_threshold:
        ign_disp = sum(comp_thresholds) + ( (100 - sum(comp_thresholds)) / 2 ) - 0.000001
    else:
        ign_disp = args.ignore_dispersion_threshold - 0.000001
    
    #RUN GATEKEEPER CHECKS=========================
    if sum(comp_thresholds) < 0-0.000001 or sum(comp_thresholds) > 100+0.000001:
        print('\n Invalid composition threshold. The composition threshold must be a number between 0-100 (inclusive)\n')
        exit()
    
    if args.dispersion < 0-0.000001 or args.dispersion > 1+0.000001:
        print('\n Invalid linear dispersion threshold. The linear dispersion threshold must be a number between 0-1 (inclusive)\n')
        exit()
        
    if ign_disp < 0-0.000001 or ign_disp > 100+0.000001:
        print('\n Invalid ignore dispersion threshold. The ignore dispersion threshold must be a number between 0-100 (inclusive)\n')
        exit()
    
    #GET USER-SPECIFIED PARAMETERS=================
    win_size = args.window_size
    amino_acids = args.amino_acids.split('_')
    amino_acids = [x.upper() for x in amino_acids]
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*_'
    other_aas = ''.join([x for x in chars if x not in ''.join(amino_acids)])
    disp_threshold = args.dispersion - 0.000001

    return win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, ign_disp
    

def output_params_header(output, args, ign_disp):
    
    output.write('*RUNTIME PARAMETERS*\n')
    output.write('Window Size: ' + str(args.window_size) + '\n')
    output.write('Amino Acid(s): ' + args.amino_acids + '\n')
    output.write('Composition Threshold(s): ' + args.composition + '\n')
    output.write('Linear Dispersion Threshold: ' + str(args.dispersion) + '\n')
    output.write('Composition to Ignore Dispersion: ' + str(round(ign_disp, 2)) + '\n\n')
    
    return output


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Identification of low-complexity domains on the basis of amino acid composition and linear dispersion', prog='LCD-Composer')

    parser.add_argument('results_file', help="""The name of your results file [tab-separated values (.tsv) format].""")
    
    parser.add_argument('-w', '--window_size', type=int, default=20, 
                        help="""Scanning window size used in the initial search. Default=20aa.
                        However, once a domain is initiated, the domain can be extended indefinitely.
                        Only integers (whole numbers) between 5-10000 are valid window sizes""")
                        
    parser.add_argument('-d', '--dispersion', type=float, default=0.5,
                        help="""Maximum amino acid dispersion threshold. Default=0.5.
                        
                        This parameter represents the degree of separation/de-mixing of amino acid(s) of interest.

                        All regions with normalized amino acid dispersion below this threshold are discarded. 
                        The value must be a decimal between 0.0 and 1.0 (inclusive).
                        
                        High dispersion values indicate well mixed sequences (e.g. GGGXGGGXGGG, GXGXGXGXGXG).
                        Low dispersion values indicate asymmetric/segregated sequences (e.g. XXGGGGGGGGG, XXXXXGGGGGG).
                        
                        Setting the dispersion threshold value = 0.0 will effectively turn off the dispersion threshold parameter.
                        Setting the dispersion threshold value = 1.0 will result in identification of ONLY perfectly mixed sequences.
                        """)
                        
    parser.add_argument('-c', '--composition', type=str, default='40',
                        help="""Composition threshold for defining X-rich regions, where X represents a particular amino acid or group of amino acids. Default = 40.0 (corresponding to 40% composition).
                        
                        Composition is calculated as the fraction of residues of interest in each window.
                        
                        Value must be between 0-100.
                        """)

    parser.add_argument('-a', '--amino_acids', type=str,
                        help="""Amino acid(s) of interest.
                        
                        Simple amino acid criteria should be an unbroken string consisting of the single-letter abbreviation for a single amino acid or a group of amino acids.
                        e.g.
                        Q       (Q composition and spacing vis-a-vis all other residues)
                        QN      (Composition and spacing of QN residues vis-a-vis all other residues)
                        QNST    (Composition and spacing of QNST residues vis-a-vis all other residues)
                        
                        Complex amino acid criteria allow for specification of different composition thresholds for distinct amino acids or distinct groups of amino acids.
                        Complex criteria use "AND" logic, meaning a sequence must contain each of the amino acid components at or above their respective composition thresholds to be considered.
                        Separate amino acids or groups of amino acids should be separated by an underscore '_'.
                        e.g.
                        Q_P     (must have specified Q composition AND must have specified P composition)
                        QN_ST   (must have specified QN composition AND must have specified ST composition)
                        G_R_Y   (must have specified G composition AND must have specified R composition AND must have specified Y composition)
                        
                        Composition thresholds for each of the separated amino acids or groups are specified by the -c parameter in a similar manner
                        e.g.
                        -a Q_P -c 40_20    (must have 40% Q composition AND must have 20% P composition)
                        -a QN_ST -c 60_25    (must have 60% QN composition AND must have 25% ST composition)
                        """)
                        
    parser.add_argument('-i', '--ignore_dispersion_threshold', type=float, default=None,
                        help="""Composition threshold at which the linear dispersion parameter is ignored when identifying and merging regions that pass the composition threshold.
                        
                        Defaults to half the distance between the user-specified composition threshold and 100%.
                        e.g.
                        A user specified composition threshold = 30% --> ignore_dispersion_threshold = 65%
                        A user specified composition threshold = 40% --> ignore_dispersion_threshold = 70%
                        A user specified composition threshold = 80% --> ignore_dispersion_threshold = 90%
                        
                        However, if the optional -i flag is used, it should be accompanied by a user-specified percentage value.
                        e.g. 
                        -i 60  --> ignore dispersion parameter if composition of the sequence is above 60%
                        -i 75  --> ignore dispersion parameter if composition of the sequence is above 75%
                        -i 100  --> ignore dispersion parameter for all sequences
                        """)

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="""Outputs per-position data (amino acid compositions and normalized spacing standard deviations) for all proteins analyzed.
                        
                        NOTE: using verbose output may increase the computation time, depending on the number of sequences, window size, composition threshold, and separation threshold.
                        """)
                        
    parser.add_argument('-k', '--file_keyword', type=str, default=None,
                        help="""Keyword required in file names in order to be analyzed by LCD-Composer.
                        Files not containing the specified keyword in the file name will be skipped.""")
                        
    parser.add_argument('-s', '--scan_subdirectories', action='store_true',
                        help="""Optional parameter to scan all subdirectories for FASTA files to analyze in
                        addition to those in the current directory.
                        
                        Default: False (subdirectories are NOT scanned)
                        
                        Use this flag in your command-line arguments in order to include 
                        files in all subdirectories for analysis.""")

                        
    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse, os, datetime, math, statistics, pickle
    args = get_args(sys.argv[1:])
    main(args)