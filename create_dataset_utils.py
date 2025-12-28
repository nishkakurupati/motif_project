import random

nucleotide_list = ["A", "T", "C", "G"]


def create_motif(ml):
    """Create and return a random DNA motif of length ml."""
    motif = ""
    for i in range (ml):
        motif = motif + random.choice(nucleotide_list)
    return motif


def create_sequences(sc, sl):
    """Create sc sequences of length sl"""
    sequences = []
    for i in range (sc):
        sequence = ""
        for k in range (sl):
            sequence = sequence + random.choice(nucleotide_list)
        sequences.append(sequence)
    return sequences


def insert_motif(sequences, motif, ml, sc, sl):
    """Inside each sequence insert motif in random position"""
    positions_of_motif = []
    for i in range (sc):
        sequence = sequences[i]
        start_index = random.randint(0, sl-ml)
        positions_of_motif.append(start_index)
        end_index = start_index + ml

        first_half = sequence[:start_index] 
        second_half = sequence[end_index:]
        new_sequence = first_half + motif + second_half

        sequences[i] = new_sequence
    return sequences, positions_of_motif


def create_files(ml,sc,sl, index):
    """Create dataset files"""
    motif = create_motif(ml)
    first_s = create_sequences(sc, sl)
    sequences, positions_of_motif = insert_motif(first_s, motif, ml, sc, sl)

    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)

    filename = start_name + "sequences.fa"
    file_obj = open(filename, 'w')
    sum = 0
    for i in range (sc):
        file_obj.write((sequences[i]) + "\n")
    file_obj.close()

    filename = start_name  + "sites.txt"
    file_obj = open(filename, 'w')
    for i in range (sc):
        file_obj.write(str(positions_of_motif[i]) + "\n")
    file_obj.close()

    filename = start_name  + "motif.txt"
    file_obj = open(filename, 'w')
    file_obj.write(motif)
    file_obj.close()

    filename = start_name  + "motif_length.txt"
    file_obj = open(filename, 'w')
    file_obj.write(str(ml))
    file_obj.close()

