
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return compute_min("".join(str_val))

def compute_reverse_comp(kmer):
    """
    Compute the reverse complement of a nucleotide sequence.
    :param str kmer: A k-mer of a DNA sequence
    :return str: returns the reverse complement of the DNA sequence
    """
    transcriber = {
                'A':'T',
                'C':'G',
                'T':'A',
                'G':'C'
                }
    new = ''
    for nuc in reversed(kmer):
        new+=transcriber[nuc]
    return new

def compute_min(kmer):
    """
    Compute the minimum binary value of a kmer, i.e. the canonical kmer
    :param str kmer: A k-mer of a DNA sequence
    :return str: the canonical version of the input k-mer
    """
    rev = compute_reverse_comp(kmer)

    if rev < kmer :
        return rev
    return kmer

def encode_nuc(letter):
    dico = {
            'A':0,
            'C':1,
            'T':2,
            'G':3
            }

    return dico[letter]

def encode_kmer(seq, k):
    kmer = 0
    for letter in seq[0:k]:
        kmer<<=2 #bit-shift
        kmer+=encode_nuc(letter)
    return kmer

def stream_kmers(seq, k):
    mask = (1<<(2*(k-1)))-1
    kmer = encode_kmer(seq, k)
    _km_list = []
    for i in range(k-2, len(seq)-k-1):
        _km_list.append(kmer)
        kmer&= mask
        kmer<<= 2
        kmer += encode_nuc(seq[i+k-1])

    return _km_list

"""
text, k = 'ATCGTGCTGATCGTATGC', 3
stream_kmers(text, int(k))
"""

"""
print(compute_min("GTGA"))
print(compute_min("ATCG"))
print(compute_min("GGGG"))
"""
