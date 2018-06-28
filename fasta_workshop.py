import sys


def read_fasta(fasta_file):
    string_line = ''

    name = open(fasta_file).readline()
    with open(fasta_file) as f:
        next(f)
        for line in f:
            string_line += line.strip('\n')

    return string_line, name


def write_fasta(sequence, output_file, desc = ''):
    f = open(output_file, 'w')
    f.write(desc+'\n')
    line_len=0
    for x in range(len(sequence)):
        f.write(sequence[x])
        line_len = line_len + 1
        if line_len == 60:
            f.write('\n')
            line_len = 0
    f.close()
    return None

def read_codon_table(codon_table='codon_table.csv'):
    amino_acid_dict = {}
    with open(codon_table) as f:
        next(f)
        for line in f:
            striped_line =line.strip() 
            split_line = striped_line.split(',')
            if split_line[0] in amino_acid_dict:
                amino_acid_dict[split_line[0]].append(split_line[2])
            else:
                amino_acid_dict[split_line[0]] = [split_line[2]]

    return amino_acid_dict

def transcribe(dna_seq, direction='-'):
    rna_seq = ''
    length_seq = len(dna_seq) - 1
    rna_dict = {'T': 'A', 'A': 'U', 'G': 'C', 'C':'G'}

    if(direction == '-'):
        for i in range( len(dna_seq) ):
            rna_seq += rna_dict[dna_seq[length_seq]]
            length_seq -= 1
    else:
        for i in range( len(dna_seq) ):
            rna_seq += rna_dict[dna_seq[i]]
            
    return rna_seq


def translate(rna_seq, codon_to_amino):
    amino_acid_out = ''
    for x in range(0, len(rna_seq),3):
        if rna_seq[x:x+3] in codon_to_amino:
            amino_acid_out += ''.join(codon_to_amino[(rna_seq[x:x+3])] ) + ''
    return amino_acid_out



def main(input_fasta, output_fasta):
    dna_seq, descriptor = read_fasta(input_fasta)
    amino_dict = read_codon_table()
    rna_seq = transcribe(dna_seq)
    pro_seq = translate(rna_seq, amino_dict)
    write_fasta(pro_seq, output_fasta, descriptor)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])