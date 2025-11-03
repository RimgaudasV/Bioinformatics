import os
import math

codontab = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',    # Serine
    'TTC': 'F', 'TTT': 'F',    # Phenylalanine
    'TTA': 'L', 'TTG': 'L',    # Leucine
    'TAC': 'Y', 'TAT': 'Y',    # Tirosine
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cisteine
    'TGA': '*',    # Stop
    'TGG': 'W',    # Tryptofan
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Proline
    'CAC': 'H', 'CAT': 'H',    # Histidine
    'CAA': 'Q', 'CAG': 'Q',    # Glutamine
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',    # Arginine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Threonine
    'AAC': 'N', 'AAT': 'N',    # Asparagine
    'AAA': 'K', 'AAG': 'K',    # Lysine
    'AGC': 'S', 'AGT': 'S',    # Serine
    'AGA': 'R', 'AGG': 'R',    # Arginine
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # Valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # Alanine
    'GAC': 'D', 'GAT': 'D',    # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # Glycine
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
START_CODONS = ['ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
AMINO_ACIDS = sorted(set(codontab.values()) - {'*'})

def read_file(filename):
    sequences = []
    with open(filename, 'r') as f:
        current_seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ''
            else:
                current_seq += line.upper()
        if current_seq:
            sequences.append(current_seq)
    return sequences


def reverse_complement(seq):
    complement = str.maketrans('ATGC', 'TACG')
    return seq.translate(complement)[::-1]



def find_orfs(seq):
    orfs = [] 
    seq_len = len(seq)

    starts = []
    for i in range(seq_len - 2):
        codon = seq[i:i+3]
        if codon in START_CODONS:
            starts.append(i)

    stops = []
    for i in range(seq_len - 2):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            stops.append(i)

    for stop_pos in stops:

        candidate_starts = []
        for s in starts:
            if s < stop_pos:
                candidate_starts.append(s)

        valid_starts = []
        for start_pos in candidate_starts:
            has_intermediate_stop = False
            for inter_stop in stops:
                if start_pos < inter_stop < stop_pos:
                    has_intermediate_stop = True
                    break
            if not has_intermediate_stop:
                valid_starts.append(start_pos)

        if valid_starts:
            farthest_start = max(valid_starts)
            orf_seq = seq[farthest_start:stop_pos + 3]

            if len(orf_seq) >= 100:
                orfs.append(orf_seq)

    return orfs


def translate_dna(seq):
    protein = ''
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        amino_acid = codontab.get(codon, 'X')
        if amino_acid == '*':
            break
        protein += amino_acid
    return protein


def codon_frequencies(protein):
    freqs = {aa: 0 for aa in AMINO_ACIDS}
    total = len(protein)
    for aa in protein:
        if aa in freqs:
            freqs[aa] += 1
    return {aa: freqs[aa] / total if total > 0 else 0 for aa in freqs}


def dicodon_frequencies(protein):
    freqs = {a1 + a2: 0 for a1 in AMINO_ACIDS for a2 in AMINO_ACIDS}
    total = len(protein) - 1 if len(protein) > 1 else 0
    for i in range(len(protein) - 1):
        dc = protein[i:i + 2]
        if dc in freqs:
            freqs[dc] += 1
    return {dc: freqs[dc] / total if total > 0 else 0 for dc in freqs}


def euclidean_distance(freqs1, freqs2):
    keys = sorted(freqs1.keys())
    s = 0
    for k in keys:
        s += (freqs1[k] - freqs2[k]) ** 2
    return math.sqrt(s)


def build_distance_matrix(objects, freq_list):
    n = len(objects)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            dist = euclidean_distance(freq_list[i], freq_list[j])
            matrix[i][j] = matrix[j][i] = dist
    return matrix


def save_to_file(filename, names, matrix):
    with open(filename, 'w') as f:
        f.write(f"{len(names)}\n")
        for i, name in enumerate(names):
            row = ' '.join(f"{matrix[i][j]:.3f}" for j in range(len(names)))
            f.write(f"{name[:10]:<10} {row}\n")

def main(file_names):
    codon_freqs = []
    dicodon_freqs = []
    data_dir = os.path.join(BASE_DIR, "data")
    for filename in file_names:
        file_path = os.path.join(data_dir, f"{filename}.fasta")
        sequences = read_file(file_path)
        all_orfs = []

        for seq in sequences:
            for s in [seq, reverse_complement(seq)]:
                orfs = find_orfs(s)
                all_orfs.extend(orfs)

        proteins = [translate_dna(orf) for orf in all_orfs if len(orf) >= 100]
        combined_protein = ''.join(proteins)

        cf = codon_frequencies(combined_protein)
        df = dicodon_frequencies(combined_protein)
        codon_freqs.append(cf)
        dicodon_freqs.append(df)

    codon_matrix = build_distance_matrix(file_names, codon_freqs)
    dicodon_matrix = build_distance_matrix(file_names, dicodon_freqs)

    results_dir = os.path.join(BASE_DIR, "results")
    os.makedirs(results_dir, exist_ok=True)

    save_to_file(os.path.join(results_dir, "codon_distance.phy"), file_names, codon_matrix)
    save_to_file(os.path.join(results_dir, "dicodon_distance.phy"), file_names, dicodon_matrix)


if __name__ == "__main__":
    files = [
        "mamalian1",
        "mamalian2",
        "mamalian3",
        "mamalian4",
        "bacterial1",
        "bacterial2",
        "bacterial3",
        "bacterial4"
    ]
    main(files)
    print("Finished")