install.packages("ggseqlogo")
library(ggseqlogo)

data(ggseqlogo_sample)

p1 = ggseqlogo( seqs_dna$MA0001.1, method = 'bits' )
p2 = ggseqlogo( seqs_dna$MA0001.1, method = 'prob' )
gridExtra::grid.arrange(p1, p2)
