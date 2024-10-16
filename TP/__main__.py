from TP.loading import load_directory, clean_seq
from TP.kmers import filter_smallest
import time

def jaccard_from_sorted_lists(lstA, lstB):
    idxA = 0
    idxB = 0

    intersection = 0
    union = 0

    while idxA < len(lstA) and idxB < len(lstB):
        union += 1
        if lstA[idxA] == lstB[idxB]:
            intersection += 1
            idxA += 1
            idxB += 1
        elif lstA[idxA] < lstB[idxB]:
            idxA += 1
        else:
            idxB += 1

    union += len(lstA) - idxA
    union += len(lstB) - idxB

    return intersection / union

if __name__ == "__main__":
    st = time.time()
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("TP2_data")

    k = 21
    sample_size = 10240
    filenames = list(files.keys())

    # Create all the kmer lists (can be expensive in memory)
    print("Computing all kmer vectors")
    kmer_lists = {}
    for filename in filenames:
        kmer_lists[filename] = []

        for seq in clean_seq(files, filename):
            kmer_lists[filename].extend(filter_smallest(seq, k, sample_size))

        kmer_lists[filename].sort()

    print("Computing Jaccard similarity for all pairs of samples")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            jaccard = jaccard_from_sorted_lists(kmer_lists[filenames[i]], kmer_lists[filenames[j]])
            print(filenames[i], filenames[j], jaccard)

    et = time.time()
    print(f"Runtime : {(et-st)/60}")