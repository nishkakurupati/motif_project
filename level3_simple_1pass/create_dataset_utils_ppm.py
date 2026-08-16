import random

nucleotide_list = ["A", "T", "C", "G"]

def create_sequences(sc, sl):
    """Create sc sequences of length sl"""
    sequences = []
    for i in range (sc):
        sequence = ""
        for k in range (sl):
            sequence = sequence + random.choice(nucleotide_list)
        sequences.append(sequence)
    return sequences

def create_motif(ml, p):
    """Create and return a random DNA motif of length ml."""
    motif = []
    for i in range(ml):
        preferred = random.choice(nucleotide_list)

        column = dict()
        for n in nucleotide_list:
            if n == preferred:
                column[n] = p
            else:
                column[n] = (1 - p) / 3
        
        motif.append(column)
    return motif


def create_random_motif(ppm):
    motif_instance = ""
    for i in ppm: #"i" is a column
        nucleotides = list(i.keys())
        probabilities = list(i.values())
        chosen = random.choices(nucleotides, probabilities)[0]
        motif_instance += chosen
    return motif_instance

def insert_motif(sequences, ppm, ml, sc, sl):
    """Insert a sampled motif (from PPM) into each sequence"""
    positions_of_motif = []
    planted_motifs = []

    for i in range(sc):
        sequence = sequences[i]
        motif_instance = create_random_motif(ppm)
        planted_motifs.append(motif_instance)
        start_index = random.randint(0, sl - ml)
        positions_of_motif.append(start_index)
        end_index = start_index + ml
        new_sequence = (sequence[:start_index] + motif_instance + sequence[end_index:])
        sequences[i] = new_sequence

    return sequences, positions_of_motif


def create_files(ml, sc, sl, index):
    """Create dataset files"""
    
    ppm = create_motif(ml, 0.8)
    sequences = create_sequences(sc,sl)
    sequences, positions_of_motif = insert_motif(sequences, ppm, ml, sc, sl)

    start_name = "data_set/file_" + str(ml) + "_" + str(sc) + "_" + str(sl) + "_" + str(index)

    filename = start_name + "sequences.fa"
    file_obj = open(filename, 'w')
    for i in range(sc):
        file_obj.write(sequences[i] + "\n")
    file_obj.close()

    filename = start_name + "sites.txt"
    file_obj = open(filename, 'w')
    for i in range(sc):
        file_obj.write(str(positions_of_motif[i]) + "\n")
    file_obj.close()
    
    filename = start_name + "motif.txt"
    file_obj = open(filename, 'w')
    for n in nucleotide_list:
        row = ""
        for col in range(len(ppm)):
            row += str(ppm[col][n]) + " "
        file_obj.write(n + " " + row.strip() + "\n")
    file_obj.close()

    filename = start_name + "motif_length.txt"
    file_obj = open(filename, 'w')
    file_obj.write(str(ml))
    file_obj.close()

