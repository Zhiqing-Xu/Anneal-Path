if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("hyperdraw")
BiocManager::install("hypergraph")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

DHE0000 <- DirectedHyperedge(c("O=CC(O)C(O)C(O)C(O)CO"), c("D-Glucose"), "EC1.1.1_0")
DHE0001 <- DirectedHyperedge(c("Glycolaldehyde","D-Threose"), c("O=CC(O)C(O)C(O)C(O)CO"), "EC4.1.2_0")
DHE0002 <- DirectedHyperedge(c("Glycolaldehyde"), c("D-Threose"), "EC4.1.2_1")
DHE0003 <- DirectedHyperedge(c("O=C(CO)OP(=O)(O)O"), c("Glycolaldehyde"), "EC1.2.1_0")
DHE0004 <- DirectedHyperedge(c("Glycolate"), c("O=C(CO)OP(=O)(O)O"), "EC2.7.1_0")
DHE0005 <- DirectedHyperedge(c("Glycolate"), c("Glycolaldehyde"), "EC1.2.1_1")
Cnodes <- c("O=CC(O)C(O)C(O)C(O)CO","D-Glucose","Glycolaldehyde","D-Threose","O=C(CO)OP(=O)(O)O","Glycolate")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 1
edgeDataDefaults(testrabph, "color") <- "black"
nodeDataDefaults(testrabph, "margin") <- 'unit(1, "mm")'
nodeDataDefaults(testrabph, "shape") <- "box"
nodeDataDefaults(testrabph, "color") <- "black"
#plot(testrabph)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "")
save_path<-paste(other.name,"/pathways/pathway","0.png",sep="")
png(file = save_path, width = 3840, height = 2160, res = 300)
plot(testrabph)
dev.off()
