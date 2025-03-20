
# load R-packages 
library(pwalign)

## Nucleotide global, local, and overlap alignments
s1 <- DNAString("ACTTCACCAGCTCCCTGGCGGTAAGTTGATCAAAGGAAACGCAAAGTTTTCAAG")
s2 <- DNAString("GTTTCACTACTTCCTTTCGGGTAAGTAAATATATAAATATATAAAAATATAATTTTCATC")

# First use a fixed substitution matrix
score_matrix <- nucleotideSubstitutionMatrix(match = 1, 
                                             mismatch = -3, 
                                             baseOnly = TRUE)

print(score_matrix)

# Global alignment
globalAlign <- pairwiseAlignment(s1, 
                                 s2, 
                                 substitutionMatrix = score_matrix,
                                 gapOpening = 0, 
                                 gapExtension = 3)
print(globalAlign)

# Ulempe med funksjonen: gir kun en alignment, selv om flere alignments gir samme score 