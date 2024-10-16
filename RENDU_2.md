# RENDU 2 TP PHYG

Matteo BETTIATI

## Jaccard distance matrix

Jaccard Matrix between the three strains of bacteria using xorshift, sampling methods and the HyperMinHash. Sample size is set to 10240 for both datasets.

X| GCA_000069965.1 | GCA_000008865.2 | GCA_030271835.1 | GCA_000013265.1 | GCA_000005845.2
:---: | :---: | :---: | :---: | :---: | :---: | 
GCA_000069965.1 | 1 | 0.0003165 | 0.0172185 | 0.0003419 | 0.0004559 |
GCA_000008865.2 | 0.0003165 | 1 | 0.0003824 | 0.1391570 | 0.1848780 |
GCA_030271835.1 | 0.0172185 | 0.0003824 | 1 | 0.0003907 | 0.0006351 |
GCA_000013265.1 | 0.0003419 | 0.1391570 | 0.0003907 | 1 | 0.1833134 |
GCA_000005845.2 | 0.0004559 | 0.1848780 | 0.0006351 | 0.1833134 | 1 |

Runtime $\approx$ 15 seconds.

X| Human | Mice | Chimpanzee
:---: | :---: | :---: | :---: 
Human | 1 | 0.0339 | 0.0397 |
Mice| 0.0339 | 1 | 0.303 |
Chimpanzee | 0.04 | 0.303 | 1 |

Runtime $\approx$ 90 minutes.

The Jaccard matrix indicates the proportion of shared kmer between pairs of bacterial genomes, and by extension an idea of how closely related two bacteria strains are.

In the Jaccard matrix distance of Bacteria strains, while the values are pretty different from with the complete set, a similar distribution is found, with the strains GCA_000008865.2, GCA_000013265.1 and GCA_000005845.2 displaying the highest value of Jaccard distance, thus showing how close they are (which make sense since all three are *E. coli* strains). Similarily, the ranking of similarity between sequences is conserved overall.

In the Jaccard matrix containing the mammals Jaccard distances results, the results are more surprising : the highest similarity of kmer is found between the Mice and the Chimpanzee, while humans share about as many kmer with the Chimp as they do with the Mice. We were expecting Humans and Chimpanzee to share the most of their kmers, and it appears not to be the case. The reason for that is unclear at the moment, as the results from bacteria are in accordance with that of the previous program, which indicates that if there is an error it should not come from the program.

## Implementation

First thing was to download the complete reference genomes of human, mice and chimpanzee from NCBI. In the sequence can be found repeated and masked elements with "N" letters and lowercases nucleotides "atcg". The "N" were filtered out using a Regex and the masked sequence were kept.

The sample size was arbitrarily set to 10240. The idea is to obtain the 10240 smallest values among all computed kmers of each genome. To avoid overrepresentation of the smallest value nucleotide "A", the binarized k-mer are hashed with a xorshift, which randomizes and equalizes the distribution of kmer value, written on 64 bits. The, in order to save memory space, a HyperMinHash (HMH) is implemented. Typically a HMH is written on 16 bits, with the first 6 bits used for the number of 0 before the heaviest bit, and the last 10 bits are used to write the 10 weakest bits of the number. However, doing a straightforward implementation of the method leads to counting errors. In fact, most bits returned by the xorshift are encoded in 40-43 bits. This results in the first 6 bits of our HMH being able to take at most 3 or 4 values. Thus the total number of values possible to encode are $\text{Possible encoding} = 4*2^{10} = 4096$, which is much smaller than our sampling size, and thus bits counted as identical while not cannot be avoided, to the point where most of the kmer selected will be identical from one strain to another, which obviously does not reflect the biological reality. This problem also arises to a lesser extent with sample sizes inferior to 4096. With a sample size of 1024, similarity rankings is somewhat similar, but the similarity between strains that are far from one another are overestimated. <br>

To avoid this problem, the HMH is therefore encoding on 32 bits and not 16, with the first 6 dedicated to the heaviest bit position, and the last 26 for the 26 first bits of the xorshifted kmer. While the memory space taken is double that of the 16-bits encoding of typical HMH, original space taken by the kmers on 64 bits is thus halved, saving some space, and retaining most of the biological information needed for the Jaccard distance calculation. <br>

