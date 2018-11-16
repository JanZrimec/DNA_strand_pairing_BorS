Prediction of DNA strand pairing as used in [Zrimec et al. 2013: Band smearing of PCR amplified bacterial 16S rRNA genes: Dependence on initial PCR target diversity.](https://www.sciencedirect.com/science/article/pii/S0167701213002467?via%3Dihub)

The structures resulting from strand pairing of two nucleotide sequences were characterized with the procedure termed BorS. DNA strand pairing was based on rules of the NN model (Santalucia 2004) and a dynamic programming algorithm adapted from Smith and Waterman's local alignment algorithm (Smith & Waterman 1981). To accommodate analysis of large datasets, an implementation of DNA strand pairing was used, which was based on the NN model unified Watson-Crick (Santalucia 1998) and mismatch parameters (Santalucia 2004). The NN model in the BorS procedure was improved by addition of a gap opening parameter (gap).
In BorS, the individual optimal mating of two nucleotide sequences a and b was determined with the following algorithm that searched for the optimal path of pairing in a matrix, in which the interaction between the two sequences was evaluated:
i . A matrix H was built as follows:
![alt text](https://raw.githubusercontent.com/JanZrimec/Band_smear_algorithm_BorS/master/Figure1.png)
where:
m 	.. size of sequence a,
n 	.. size of sequence b,
b 	.. element of H, which represents the interaction between two nucleotides from the 	sequences a and b (ai-1ai/bj-1bj), 
NN 	.. free energy parameter ΔG°37, which evaluate the interaction (ai-1ai/bj-1bj), 
NNinit 	.. initialization parameter in the NN model (Santalucia 2004),
gap 	.. gap opening parameter, the value of which was optimized 

During building of the matrix H, pointers to previous elements, from which the consequent elements were calculated, were continuously stored.
ii. A traceback procedure was performed, where the greatest value in the last column or the last row in matrix H was found and the optimal path that led up to the maximum value was determined according to the array of pointers, stored while building H. The procedure was set to prefer a diagonal path without gap opening, since this corresponded to the NN strand pairing rules. The optimal path represented the optimal pairing of sequences a and b depending on the dissociation energy of the structures.

The BorS procedure tested all interactions between pairs of sequences that were determined as the most stable structural variants of the DNA with a lower dissociation energy. The following structural characteristics were determined: (i) the dissociation energy of DNA (dG), (ii) the entalpy (dH) and entropy (dS), (iii) the proportion of the number of base pairs with respect to length sequence (bp) (iv) of the average number of nucleotides in loops, which surround the base pairs, depending on the size of the area of overlap of the two sequences (loop), and (v) the average number of free nucleotides in both chains hanging outside the structure (dang).

An additional Needleman Wunsch (NW) global alingment implementation was added. The functions were tested in Matlab and Octave.

Use:
