
###############################
# Bivalve Phylogenetic Tree
# 5/28/2020
# Alex Van Plantinga

library("rentrez")
library("Biostrings")
library("GenomicRanges")
library("clusterProfiler")
library("DOSE")
library("org.Hs.eg.db")
library("treeio")
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("msa")
library("seqinr")
library("phylobase")

species_list = c("Amblema plicata", "Cyrtonaias tampicoensis", 
                    "Lampsilis cardium", "Margaritifera margaritifera",
                    "Ellipsaria lineolata", "Elliptio complanata",
                    "Megalonaias nervosa", "Crassostrea virginica",
                    "Arctica islandica", # "Caenorhabditis elegans" , # nematode worm outgroup
                    "Conus ermineus", "Tridacna gigas")

ll = c()
for(i in 1:length(species_list)){
  #t = paste(species_list[i], "[Organism] AND ND1[Gene]", sep = "")
  t = paste(species_list[i], "[Organism] AND COI[Gene]", sep = "")
  m = entrez_search(db = "nuccore", term = t, retmax = 1)
  m = entrez_fetch(db="nuccore", id=m$ids, rettype="fasta")
  ll = append(ll, m)
}
setwd("C:\\Users\\Lenovo\\BioconductorAnimals")
#write(ll, "mussels_ND1.fasta", sep="\n")
write(ll, "mussels_COI.fasta", sep="\n")
mussels_COI_seqinr_format <- read.fasta("mussels_COI.fasta")

musselSeq <- readAAStringSet("mussels_COI.fasta")
musselAln <- msa(musselSeq)
## use default substitution matrix

musselAln

musselAln2 <- msaConvert(musselAln, type="seqinr::alignment")

d <- dist.alignment(musselAln2, "identity")

musselTree <- nj(d)

df = data.frame(color=sample(c('red', 'blue', 'green'), 
                             length(musselTree$tip.label), replace=T))

rownames(df) = musselTree$tip.label

#musselTree = phylo4d(as(musselTree, 'phylo4'), df)

for(i in 1:length(musselTree$tip.label)){
  musselTree$tip.label[i] = paste(strsplit(musselTree$tip.label, " ")[[i]][2:3], collapse = " ")
}

plot(musselTree, main="Phylogenetic Tree Based on Mollusk COI Sequences", cex = 0.8)

plot.phylo(musselTree, type="unrooted", cex = 0.5)


#musselTree = groupOTU(musselTree, musselTree$tip.label)

# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
#ggtree(musselTree,  layout='circular') + geom_tiplab(size=3, aes(angle=angle)) + theme(legend.position="none")+ aes(color=I(color))
ggtree(musselTree,  layout='slanted',  branch.length="none") + geom_tiplab(size=3, aes(angle=90)) + theme(legend.position="none") + 
  coord_flip() + ggplot2::xlim(0, 10)


tree <- groupClade(musselTree)

tree <- groupClade(musselTree, .node = c(12,13,14,18)+1)
tree$edge.length[2:10] = -tree$edge.length[2:10]

tree <- groupClade(musselTree, .node = c(12,13,14,18))

p = ggtree(tree, aes(color=group, linetype=group)) + 
  ggtitle("Unionid and Invertebrate Phylogenetic Tree with ND1 Sequences") +
  #geom_text(aes(label = node)) + # run this line to identify nodes to rotate
  geom_tiplab(aes(subset=(group==1))) +
  geom_tiplab(aes(subset=(group==2))) +
  geom_tiplab(aes(subset=(group==3))) +
  geom_tiplab(aes(subset=(group==4)), hjust = 1 ) + 
  ggplot2::xlim(-1.3, 1.5) + 
  scale_color_manual(labels = c("Uniodae1", "Unioda", "Uniodae2", "Other Invertebrates"), 
                     values = c("black", "purple", "navy", "magenta") ) +
  guides(color = guide_legend(override.aes = list(linetype = c('solid','solid', 'solid', 'solid'))),
         linetype = FALSE)

p
rotate(p, 19)
