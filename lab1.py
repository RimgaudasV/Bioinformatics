import os
import math

CODON_TABLE = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',   # Serine
    'TTC': 'F', 'TTT': 'F',    # Phenylalanine
    'TTA': 'L', 'TTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',   # Leucine
    'TAC': 'Y', 'TAT': 'Y',    # Tirosine
    'TAA': '*', 'TAG': '*', 'TGA': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cisteine
    'TGG': 'W',    # Tryptofan
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Proline
    'CAC': 'H', 'CAT': 'H',    # Histidine
    'CAA': 'Q', 'CAG': 'Q',    # Glutamine
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',   # Arginine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Threonine
    'AAC': 'N', 'AAT': 'N',    # Asparagine
    'AAA': 'K', 'AAG': 'K',    # Lysine
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # Valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # Alanine
    'GAC': 'D', 'GAT': 'D',    # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # Glycine
}

FILES = [
    "mamalian1", "mamalian2", "mamalian3", "mamalian4", 
    "bacterial1", "bacterial2", "bacterial3", "bacterial4"
]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
START_CODONS = ['ATG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']
AMINO_ACIDS = sorted(set(CODON_TABLE.values()) - {'*'})

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


# 1-3
def find_orfs(seq, min_len=100):
    orfs = []
    for frame in range(3):
        start_pos = None

        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]

            if codon in START_CODONS and start_pos is None:
                start_pos = i

            elif codon in STOP_CODONS and start_pos is not None:
                orf_seq = seq[start_pos:i+3]
                if len(orf_seq) >= min_len:
                    orfs.append(orf_seq)
                start_pos = None

    return orfs


# 4
def convert_to_protein_sequence(orf):
    protein = ''
    for i in range(0, len(orf) - 2, 3):
        codon = orf[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, 'X')
        if amino_acid == '*':
            break
        protein += amino_acid
    return protein

# 5
def codon_frequencies(protein):
    freqs = {aa: 0 for aa in AMINO_ACIDS}
    total = len(protein)
    for aa in protein:
        if aa in freqs:
            freqs[aa] += 1
    for aa in freqs:
        if total > 0:
            freqs[aa] = freqs[aa] / total
        else:
            freqs[aa] = 0
    return freqs

# 5
def dicodon_frequencies(protein):
    freqs = {a1 + a2: 0 for a1 in AMINO_ACIDS for a2 in AMINO_ACIDS}
    total = len(protein) - 1 if len(protein) > 1 else 0
    for i in range(len(protein) - 1):
        dc = protein[i:i + 2]
        if dc in freqs:
            freqs[dc] += 1
    for aa in freqs:
        if total > 0:
            freqs[aa] = freqs[aa] / total
        else:
            freqs[aa] = 0
    return freqs

# 6
def create_distance_matrix(n, freq_list):
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            s = 0
            for k in sorted(freq_list[i].keys()):
                s += (freq_list[i][k] - freq_list[j][k]) ** 2
            dist = math.sqrt(s)
            matrix[i][j] = matrix[j][i] = dist
    return matrix


def save_to_file(filename, names, matrix):
    with open(filename, 'w') as f:
        f.write(f"{len(names)}\n")
        for i, name in enumerate(names):
            row = ' '.join(f"{matrix[i][j]:.3f}" for j in range(len(names)))
            f.write(f"{name[:10]:<10} {row}\n")


if __name__ == "__main__":
    codon_freqs = []
    dicodon_freqs = []
    data_dir = os.path.join(BASE_DIR, "data")
    for filename in FILES:
        file_path = os.path.join(data_dir, f"{filename}.fasta")
        sequences = read_file(file_path)
        all_orfs = []

        for seq in sequences:
            for s in [seq, reverse_complement(seq)]:
                orfs = find_orfs(s)
                all_orfs.extend(orfs)

        proteins = [convert_to_protein_sequence(orf) for orf in all_orfs]
        combined_protein = ''.join(proteins)

        cf = codon_frequencies(combined_protein)
        df = dicodon_frequencies(combined_protein)
        codon_freqs.append(cf)
        dicodon_freqs.append(df)

    codon_matrix = create_distance_matrix(len(FILES), codon_freqs)
    dicodon_matrix = create_distance_matrix(len(FILES), dicodon_freqs)

    results_dir = os.path.join(BASE_DIR, "results")
    os.makedirs(results_dir, exist_ok=True)

    save_to_file(os.path.join(results_dir, "codon_distance.phy"), FILES, codon_matrix)
    save_to_file(os.path.join(results_dir, "dicodon_distance.phy"), FILES, dicodon_matrix)