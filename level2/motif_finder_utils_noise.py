def get_sequences_from_file(file_name):
    """Finding sequence length as well as returning all sequences from a given file"""
    sequences = []

    filename = 'data_set/file_' + file_name + "sequences.fa"
    file_obj = open(filename, 'r')
    for line in file_obj:
        sequences.append(line.rstrip('\n'))
        sl = len(line) - 1
    file_obj.close()

    return sl, sequences

def get_possible_motifs(sl,ml, sequence):
    """Derive all possible motifs from splitting the sequence"""
    possible_motifs = []

    for i in range (sl - ml):
        pos_motif = sequence[i:i + ml]
        possible_motifs.append(pos_motif)

    return possible_motifs


def build_motif_sites_dicts(sequences, sl, ml):
    motif_dicts = []
    sites_dicts = []
    for i in range (len(sequences)):
        motif_dict = dict()
        sites_dict = dict()

        possible_motifs = get_possible_motifs(sl,ml, sequences[i])

        for j in range(len(possible_motifs)):
            if possible_motifs[j] not in motif_dict:
                motif_dict[possible_motifs[j]] = 1
                sites_dict[possible_motifs[j]] = j
        motif_dicts.append(motif_dict)
        sites_dicts.append(sites_dict)

    return motif_dicts, sites_dicts


def get_best_motif(motif_dicts):
    motif_count_dict = dict()
    for i in range (len(motif_dicts)):
        for key in motif_dicts[i]:
            if key not in motif_count_dict:
                motif_count_dict[key] = 0
            else:
                motif_count_dict[key] +=1

    best_motif_count = 0
    best_motif = ""

    for key in motif_count_dict:
        if motif_count_dict[key] >= best_motif_count:
            best_motif_count = motif_count_dict[key]
            best_motif = key
    
    return best_motif


def get_sites(motif, sites_dicts):
    sites = []
    for i in range (len(sites_dicts)):
        if motif in sites_dicts[i]:
            sites.append(sites_dicts[i][motif])
        else:
            sites.append(-1)
    return sites


def find_motif(sequences, sl, ml):
    """
    Algorithm to find motif by searching through all sequences
    from a list of possible motifs
    """

    motif_dicts, sites_dicts = build_motif_sites_dicts(sequences, sl, ml)

    motif = get_best_motif(motif_dicts)

    sites = get_sites(motif, sites_dicts)

    return motif, sites



def create_output_files(motif, ml,sc,sl, index, positions_of_motif):
    """Returning predicted motif and predicted sites to output files """

    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)
    filename = start_name + "predictedmotif.txt"
    file_obj = open(filename, 'w')
    file_obj.write(motif)
    file_obj.close()

    filename = start_name + "predictedsites.txt"
    file_obj = open(filename, 'w')

    for i in range (sc):
        file_obj.write(str(positions_of_motif[i]) + "\n")
    file_obj.close()
                
                
