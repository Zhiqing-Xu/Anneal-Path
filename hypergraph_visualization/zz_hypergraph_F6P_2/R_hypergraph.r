#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

DHEXXXX <- DirectedHyperedge(c("G6P"), c("F6P"), "")
#DHEYYYY <- DirectedHyperedge(c("DAHP"), c("YYY"), "RYYYY")
DHE0000 <- DirectedHyperedge(c("E4P","PEP"), c("DAHP"), "EC2.5.1")
DHE0001 <- DirectedHyperedge(c("2PG"), c("PEP"), "EC4.2.1")
DHE0002 <- DirectedHyperedge(c("3PG"), c("2PG"), "EC5.4.2")
DHE0003 <- DirectedHyperedge(c("GA3P","F6P"), c("E4P","X5P"), "EC2.2.1")
DHE0004 <- DirectedHyperedge(c("GA3P"), c("3PG"), "EC1.2.1")
DHE0005 <- DirectedHyperedge(c("F16bP"), c("GA3P","DHAP"), "EC4.1.2")
DHE0006 <- DirectedHyperedge(c("F6P"), c("F16bP"), "EC2.7.1")

Cnodes <- c("G6P","2PG","PEP","DAHP","3PG","E4P","X5P","F16bP","DHAP","GA3P","F6P")
Rnodes <- list(DHEXXXX, DHE0005, DHE0001, DHE0002, DHE0003, DHE0004, DHE0000, DHE0006)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)

testrabph <- graphLayout(hgbph)


edgeDataDefaults(testrabph, "lwd") <- 1
edgeDataDefaults(testrabph, "color") <- "black"

nodeDataDefaults(testrabph, "margin") <- 'unit(0.5, "mm")'
nodeDataDefaults(testrabph, "shape") <- "box"
nodeDataDefaults(testrabph, "color") <- "black"



plot(testrabph)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#



png(file = "test_1.png", width = 2000, height = 1000, res = 300)
plot(testrabph)
dev.off()







#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#



# nodeData(testrabph, c("AAA"), "fontsize") <- '50'
# testrabph@graph@AgNode[[8]]@color<-"pink"
# nodeData(testrabph, c("AAA"), "fillcolor") <- "blue"
# dev.new(width = 1920 * 2, height = 1920 * 2, unit = "px")


# pdf("mygraph.pdf", width = 11, height = 8)