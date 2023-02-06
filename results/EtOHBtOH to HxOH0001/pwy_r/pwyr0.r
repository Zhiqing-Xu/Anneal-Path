#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

Cnodes <- c("C")
Rnodes <- list("E")
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
