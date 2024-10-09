from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str
from collections import Counter
import time

def jaccard(fileA, fileB, k):
    str_fileA = [kmer2str(km_bin, k) for km_bin in fileA]
    str_fileB = [kmer2str(km_bin, k) for km_bin in fileB]
    dicokmA, dicokmB = Counter(str_fileA), Counter(str_fileB)

    #get unique kmers
    outer = set(list(dicokmA.keys()))|set(list(dicokmB.keys()))
    #get common kmers
    inter = set(list(dicokmA.keys()))&set(list(dicokmB.keys()))

    jac = sum([min([dicokmA[km], dicokmB[km]]) for km in inter])/sum([max([dicokmA[km], dicokmB[km]]) for km in outer])
    return jac

if __name__ == "__main__":
    st = time.time()
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")

    k = 21
    filenames = list(files.keys())

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            fileA, fileB = [], []
            for index_seq in range(0, len(files[filenames[i]])):
                fileA += stream_kmers(files[filenames[i]][index_seq], k)

            for index_seq in range(0, len(files[filenames[j]])):
                fileB += stream_kmers(files[filenames[j]][index_seq], k)

            jaq = jaccard(fileA, fileB, k)
            print(filenames[i], filenames[j], jaq)

    et = time.time()
    print(f"Runtime : {(et-st)/60}")
