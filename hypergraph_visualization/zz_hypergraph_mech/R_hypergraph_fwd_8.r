#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)


DHE0000 <- DirectedHyperedge(c("C0004","C0005"), c("C0002","C0003"), "R0000")
DHE0001 <- DirectedHyperedge(c("C0000","C0001"), c("C0004","C0005"), "R0001")
DHE0002 <- DirectedHyperedge(c("C0004","C0005"), c("C0006","C0007"), "R0002")
DHE0003 <- DirectedHyperedge(c("C0004","C0010"), c("C0008","C0009"), "R0003")
DHE0004 <- DirectedHyperedge(c("C0000","C0001"), c("C0011","C0010"), "R0004")
DHE0005 <- DirectedHyperedge(c("C0004","C0010"), c("C0013","C0012"), "R0005")
DHE0006 <- DirectedHyperedge(c("C0004","C0011"), c("C0015","C0014"), "R0006")
DHE0007 <- DirectedHyperedge(c("C0004","C0011"), c("C0017","C0016"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0005","C0010"), c("C0019","C0018"), "R0008")
DHE0009 <- DirectedHyperedge(c("C0005","C0010"), c("C0020","C0021"), "R0009")
DHE0010 <- DirectedHyperedge(c("C0005","C0011"), c("C0022","C0023"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0005","C0011"), c("C0024","C0025"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0011","C0010"), c("C0026","C0027"), "R0012")
DHE0013 <- DirectedHyperedge(c("C0011","C0010"), c("C0028","C0029"), "R0013")
Cnodes <- c("C0004","C0005","C0002","C0003","C0000","C0001","C0006","C0007","C0010","C0008","C0009","C0011","C0013","C0012","C0015","C0014","C0017","C0016","C0019","C0018","C0020","C0021","C0022","C0023","C0024","C0025","C0026","C0027","C0028","C0029")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007,DHE0008,DHE0009,DHE0010,DHE0011,DHE0012,DHE0013)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)




testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 2
edgeDataDefaults(testrabph, "color") <- "black"


nodeDataDefaults(testrabph, "margin") <- 'unit(-4, "mm")'
nodeDataDefaults(testrabph, "shape") <- "circle"
nodeDataDefaults(testrabph, "color") <- "transparent"

plot(testrabph)

png(file = "test.png", width = 1920 * 2, height = 1080 * 2, res = 300)
plot(testrabph)
dev.off()