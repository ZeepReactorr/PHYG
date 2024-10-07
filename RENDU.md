# RENDU k-mer TP PHYG

Matteo BETTIATI

## Jaccard distance matrix

X | GCA_000069965.1 | GCA_000008865.2 | GCA_030271835.1 | GCA_000013265.1 | GCA_000005845.2
:---: | :---: | :---: | :---: | :---: | :---: | 
GCA_000069965.1 | X | 0.00231381248555095 | 0.031129655942524762 |0.0024370372626886356|0.0025674850826432915|
GCA_000008865.2 |0.00231381248555095 | X |0.002317306485120559 | 0.30704869735137047 | 0.4364853581146792 |
GCA_030271835.1 | 0.031129655942524762 |0.002317306485120559 |X|0.0024332528798458454|0.0025757914572849924|
GCA_000013265.1|0.0024370372626886356 |0.30704869735137047 |0.0024332528798458454|X|0.34100822682570797|
GCA_000005845.2|0.0025674850826432915 |0.4364853581146792|0.0025757914572849924|0.34100822682570797|X|

The Jaccard matrix indicates the proportion of shared kmer between pairs of bacterial genomes, and by extension an idea of how closely related the two bacteria strains are.

## Implementation

<p align="justify">
k-mer are produced using a sliding window large of k nucleotides, with a binary mask masking the first nucleotide and the next nucleotide being directly encoded. The entire k-mer is encoded in binary format, and then translated back to a string of nucleotides. At the same time, the canonical k-mer is determined and if lighter to encode, is kept in place of the one determined previously. Following this step, for both files for which we calculate the Jaccard distance, we make a dictionnary counting each kmer sequence using the `Counter` class from the standard library of Python. Then by leveraging the `set` type of python we compute the Jaccard.
</p>