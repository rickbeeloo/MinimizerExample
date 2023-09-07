import btllib  # https://www.theoj.org/joss-papers/joss.04720/10.21105.joss.04720.pdf
import collections 
from Bio import SeqIO

def get_minimizers(fasta_path, k, w):
    # we use btllib to get the minimzers.
    with btllib.Indexlr(fasta_path, k, w, btllib.IndexlrFlag.LONG_MODE) as indexer:
        minimizers = []
        for record in indexer: # Could also do this for each sequence in the FASTA seperate
            minimizers.extend([mx.out_hash for mx in record.minimizers])
        return dict(collections.Counter(minimizers))


def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = dna_sequence[::-1]  # Reverse the sequence
    complement_seq = ''.join([complement_dict[base] for base in reverse_seq])
    return complement_seq

# Adjusted from Heng li: https://github.com/lh3/kmer-cnt
def get_kmers(fasta_path, k):
    h = dict()
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = record.seq.upper()
        l = len(seq)
        if l < k:
            return
        for i in range(l - k + 1):
            kmer_for = seq[i:(i+k)]                                    
            if 'N' in kmer_for: continue # Skip N chars
            kmer_rev = reverse_complement(kmer_for)
            if kmer_for < kmer_rev: # get the smallest k-mer, we dont want to have both AAA and TTT as they are rev compl
                kmer = kmer_for
            else: 
                kmer = kmer_rev
            if kmer in h:
                h[kmer] += 1
            else: 
                h[kmer] = 1
    return h
                    

fasta = "GCA_001580015.1_ASM158001v1_genomic.fasta"

minimizers = get_minimizers(fasta, 31, 10).keys()
kmers = get_kmers(fasta, 31).keys() # This is slow, better use c++ to deal with this as well
print(len(minimizers))
print(len(kmers))
print("Saves: ", 100 - round(len(minimizers) / len(kmers) * 100, 0), "%") 

# 426154
# 2336989
# Saves:  82.0 %