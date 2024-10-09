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
First of all, a pair of sequence files are loaded. If multiple sequences are found in a file, their k-mer are constructed independently, and then the lists are concatenated. k-mer are generated using a sliding window large of k nucleotides, with a binary mask masking the first nucleotide and the next nucleotide being directly encoded. The entire k-mer is encoded in binary format. Following this, the canonical k-mer is determined as its reverse complement. Only the lighter one to encode is kept in the final k-mer list. Finally, the binarized k-mer list are computed as strings again, then using the properties of the <c>set</c> type in Python, the intersecting and union sets of k-mer are determined, and from there the Jaccard index is computed using the formula :
</p>

$$
J = \frac{A \cap B}{A \cup B}
$$
