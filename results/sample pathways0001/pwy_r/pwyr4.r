if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("hyperdraw")
BiocManager::install("hypergraph")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

DHE0000 <- DirectedHyperedge(c("Serine"), c("Pyruvate"), "EC4.3.1_0")
DHE0001 <- DirectedHyperedge(c("Aldehyde","Glycine"), c("Serine"), "EC4.1.2_0")
DHE0002 <- DirectedHyperedge(c("D-Fructose 6-phosphate"), c("Aldehyde","D-Xylose-5-phosphate"), "EC4.1.2_1")
DHE0003 <- DirectedHyperedge(c("Glyoxylate"), c("Glycine"), "EC1.4.1_0")
DHE0004 <- DirectedHyperedge(c("beta-D-Fructose 6-phosphate"), c("D-Fructose 6-phosphate"), "EC1.1.1_0")
DHE0005 <- DirectedHyperedge(c("D-Allose 6-phosphate"), c("Glyoxylate","D-Xylose-5-phosphate"), "EC4.1.2_2")
DHE0006 <- DirectedHyperedge(c("D-Glucose 6-phosphate"), c("D-Allose 6-phosphate"), "EC1.1.1_1")
DHE0007 <- DirectedHyperedge(c("D-Glucose 6-phosphate"), c("beta-D-Fructose 6-phosphate"), "EC5.3.1_0")
Cnodes <- c("Serine","Pyruvate","Aldehyde","Glycine","D-Fructose 6-phosphate","D-Xylose-5-phosphate","Glyoxylate","beta-D-Fructose 6-phosphate","D-Allose 6-phosphate","D-Glucose 6-phosphate")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007)
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
save_path<-paste(other.name,"/pathways/pathway","4.png",sep="")
png(file = save_path, width = 1920, height = 1080, res = 300)
plot(testrabph)
dev.off()
