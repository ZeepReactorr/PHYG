import math, heapq
import random as rd 

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
    rev = compute_reverse_comp(kmer)

    if rev < kmer :
        return rev
    return kmer

def encode_nucl(nucl):
    """ 
    Encode a nucleotide into a 2-bit integer
    :param str nucl: The nucleotide to encode
    :return (int, int): The encoded nucleotide and its reverse complement
    """
    encoded = (ord(nucl) >> 1) & 0b11 # Extract the two bits of the ascii code that represent the nucleotide
    rencoded = (encoded + 2) & 0b11 # Complement encoding with bit tricks. Avoid slow if statement.

    return encoded, rencoded

def stream_kmers(seq, k):
    # Initialize the kmer and its reverse complement
    kmer = 0
    rkmer = 0

    # Add the first k-1 nucleotides to the first kmer and its reverse complement
    for i in range(k-1):
        nucl, rnucl = encode_nucl(seq[i])
        kmer |= nucl << (2*(k-2-i))
        rkmer |= rnucl << (2*(i+1))

    mask = (1 << (2*(k-1))) - 1

    # Yield the kmers
    for i in range(k-1, len(seq)):
        nucl, rnucl = encode_nucl(seq[i])
        # Remove the leftmost nucleotide from the kmer 
        kmer &= mask
        # Shift the kmer to make space for the new nucleotide
        kmer <<= 2
        # Add the new nucleotide to the kmer
        kmer |= nucl
        # Make space for the new nucleotide in the reverse kmer (remove the rightmost nucleotide by side effect)
        rkmer >>= 2
        # Add the new nucleotide to the reverse kmer
        rkmer |= rnucl << (2*(k-1))

        yield xorshift(min(kmer, rkmer))

def HMH(kmer : int):
    try :
        one = 64-f"{kmer:64b}".index('1',)-1 #find bit of highest index 
    except ValueError:
        one = 0
        
    n_zero = f'{64-one:06b}' #compute number of zeros before heaviest bit
    
    msk = (1 << (26))-1 #computation of the mask to obtain the 10 lower bits of the kmer
    kmer &= msk 
    kmer = n_zero + f"{kmer:026b}"
    return int(kmer)

def filter_smallest(seq : str, k : int, s: int):
    lst = [-math.inf for _ in range(s)]

    for kmer in stream_kmers(seq, k):
        if lst[0] < -kmer:
            heapq.heappushpop(lst, -kmer)

    return [HMH(-elt) for elt in lst if elt != -math.inf]

def xorshift(val):
    val ^= val << 13
    val &= 0xFFFFFFFFFFFFFFFFF
    val ^= val >> 7
    val ^= val << 17
    val &= 0xFFFFFFFFFFFFFFFFF
    return val

import time 

def bin_search_one(kmer : int):
    str_kmer = f"{kmer:64b}"

    size_mask = len(str_kmer)//2
    idx_one = 0

    while kmer != 1:      
        msk_right = (1 << (size_mask))-1
        msk_left = ~msk_right

        kmer_l = kmer&msk_left

        if kmer_l != 0:
            idx_one += size_mask
            kmer = kmer_l>>size_mask

        else:
            size_mask = size_mask//2

    return idx_one

def test_bin_search_one(n_test : int):
    """
    Testing function. Test ability to find the heaviest bit in 
    binarized kmer.

    And runtime comparison between a binary search of the heaviest bit
    and an O(n) built-in method of python : find()

    ==> result, find() is about twice as fast
    """
    for _ in range(n_test):
        val = rd.randint(1, 2**32)
        str_val = f"{val:32b}"
        h_bit = 32-str_val.index('1')-1

        assert(bin_search_one(val)==h_bit)
    
    st_bin = time.time()
    for _ in range(n_test):
        val = rd.randint(1, 2**64)
        bin_search_one(val)
    
    print(f"With bin search : {time.time() - st_bin}")

    st_fin = time.time()
    for _ in range(n_test):
        val = rd.randint(1, 2**64)
        h_bit = 64-f"{val:64b}".index('1')-1
    print(f"With bin search : {time.time() - st_fin}")
    
"""
test_bin_search_one(2000000)

"""