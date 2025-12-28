def get_ml_from_file(file_name):
    """Finding length of motif"""
    filename = 'data_set/file_' + file_name + "motif_length.txt"
    file_obj = open(filename, 'r')
    for line in file_obj:
	    ml = int(line)
    file_obj.close()
    return ml


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


def create_possible_motifs(sl,ml, sequences):
    """Derive all possible motifs from splitting the first sequence"""
    possible_motifs = []

    for i in range (sl - ml):
        pos_motif = sequences[0][i:i + ml]
        possible_motifs.append(pos_motif)

    return possible_motifs


def find_motif(possible_motifs, sequences, sl, ml):
    """
    Algorithm to find motif by searching through all sequences
    from a list of possible motifs
    """
    for i in range (len(sequences)):
        possible_motifs_copy = []
        for j in range(len(possible_motifs)):
            possible_motifs_copy.append(possible_motifs[j])
        for j in range (len(possible_motifs_copy)):
            check_motif = possible_motifs_copy[j]

            motif_there = False
            for k in range (sl - ml + 1):
                motif = sequences[i][k:k + ml]

                if motif == check_motif:
                    motif_there = True

            if motif_there == False:
                possible_motifs.remove(check_motif)

    matching_motif = possible_motifs[0]

    return matching_motif


def get_sites(sl,ml,motif, sequences):
    """Searching through each sequence to find the position of the motif in it"""
    positions_of_motif = []
    for k in range (len(sequences)):
        sequence = sequences[k]
        position = 0
        for i in range (sl - ml + 1):
            pos_motif = sequences[k][i:i + ml]
            if pos_motif == motif:
                positions_of_motif.append(i)
                break

    return positions_of_motif

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
                
                
