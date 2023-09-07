import btllib  # https://www.theoj.org/joss-papers/joss.04720/10.21105.joss.04720.pdf
import collections 
from Bio import SeqIO


# Using minimizers is different from using smaller Ks. A smaller k-mer will decrease the number of k-mers but also the resolution you have.
# A minimizer will pick the smallest k-mer in a window but the k-mer size can still be large. 
# You might lose important k-mers this way but if it's more composition than exact k-mer that matter it wouldn't have much negative effect.
# Take for example:
# ACTACGAT
# The 3-mer would be:
# ACT
# CTA
# TAC
# ACG
# CGA
# GAT
# The minimizers with a window size of 4 would be:
# K-mer: ACT, Minimizer: ACT
# K-mer: CTA, Minimizer: ACT
# K-mer: TAC, Minimizer: ACT
# K-mer: ACG, Minimizer: ACG
# K-mer: CGA, Minimizer: ACG
# K-mer: GAT, Minimizer: ACG
# This would just summarize to ACT: 3, ACG: 3, which is just two kmers compared to 6. 
# If we then get another sequence and mutate it in two places:
# ACCACGGT
# We would get 
# K-mer: ACC, Minimizer: ACC
# K-mer: CCA, Minimizer: ACC
# K-mer: CAC, Minimizer: CAC
# K-mer: ACG, Minimizer: ACG
# K-mer: CGG, Minimizer: ACG
# K-mer: GGT, Minimizer: ACG
# So: ACC:2, CAC:1, ACG:3
# So you still notice the shared "ACG" minimizer among the sequences. 
# In total, we would now have ACT ACC CAC  ACG which is 4 kmers vs 11 ( 6 for ACTACGAT, 6 for ACCACGGT - 1 for ACG as it is shared)

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
# This is slow, better use c++ to deal with this as well
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

k_mer_size = 9
minimizers = get_minimizers(fasta, k_mer_size, 10).keys()
kmers = get_kmers(fasta, k_mer_size).keys() 
print(len(minimizers))
print(len(kmers))
print("Kmer size: {} saves: ".format(k_mer_size), 100 - round(len(minimizers) / len(kmers) * 100, 0), "%") 


k_mer_size = 31
minimizers = get_minimizers(fasta, k_mer_size, 10).keys()
kmers = get_kmers(fasta, k_mer_size).keys() 
print(len(minimizers))
print(len(kmers))
print("Kmer size: {} saves: ".format(k_mer_size), 100 - round(len(minimizers) / len(kmers) * 100, 0), "%") 


# 46009
# 127196
# Kmer size: 9 saves:  64.0 %

# 426154
# 2336989
# Kmer size: 31 saves:  82.0 %