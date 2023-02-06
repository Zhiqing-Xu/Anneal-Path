#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)


DHE0000 <- DirectedHyperedge(c("C0002"), c("C0001"), "R0000")
DHE0001 <- DirectedHyperedge(c("C0003"), c("C0002"), "R0001")
DHE0002 <- DirectedHyperedge(c("C0000"), c("C0003"), "R0002")
DHE0003 <- DirectedHyperedge(c("C0002"), c("C0004"), "R0003")
DHE0004 <- DirectedHyperedge(c("C0006"), c("C0005"), "R0004")
DHE0005 <- DirectedHyperedge(c("C0003"), c("C0006"), "R0005")
DHE0006 <- DirectedHyperedge(c("C0006"), c("C0007"), "R0006")
DHE0007 <- DirectedHyperedge(c("C0009"), c("C0008"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0010"), c("C0009"), "R0008")
DHE0009 <- DirectedHyperedge(c("C0000"), c("C0010"), "R0009")
DHE0010 <- DirectedHyperedge(c("C0009"), c("C0011"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0013"), c("C0012"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0010"), c("C0013"), "R0012")
DHE0013 <- DirectedHyperedge(c("C0013"), c("C0014"), "R0013")
DHE0014 <- DirectedHyperedge(c("C0017"), c("C0016"), "R0014")
DHE0015 <- DirectedHyperedge(c("C0018"), c("C0017"), "R0015")
DHE0016 <- DirectedHyperedge(c("C0015"), c("C0018"), "R0016")
DHE0017 <- DirectedHyperedge(c("C0017"), c("C0019"), "R0017")
DHE0018 <- DirectedHyperedge(c("C0021"), c("C0020"), "R0018")
DHE0019 <- DirectedHyperedge(c("C0018"), c("C0021"), "R0019")
DHE0020 <- DirectedHyperedge(c("C0021"), c("C0022"), "R0020")
DHE0021 <- DirectedHyperedge(c("C0024"), c("C0023"), "R0021")
DHE0022 <- DirectedHyperedge(c("C0025"), c("C0024"), "R0022")
DHE0023 <- DirectedHyperedge(c("C0015"), c("C0025"), "R0023")
DHE0024 <- DirectedHyperedge(c("C0024"), c("C0026"), "R0024")
DHE0025 <- DirectedHyperedge(c("C0028"), c("C0027"), "R0025")
DHE0026 <- DirectedHyperedge(c("C0025"), c("C0028"), "R0026")
DHE0027 <- DirectedHyperedge(c("C0028"), c("C0029"), "R0027")
Cnodes <- c("C0002","C0001","C0003","C0000","C0004","C0006","C0005","C0007","C0009","C0008","C0010","C0011","C0013","C0012","C0014","C0017","C0016","C0018","C0015","C0019","C0021","C0020","C0022","C0024","C0023","C0025","C0026","C0028","C0027","C0029")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007,DHE0008,DHE0009,DHE0010,DHE0011,DHE0012,DHE0013,DHE0014,DHE0015,DHE0016,DHE0017,DHE0018,DHE0019,DHE0020,DHE0021,DHE0022,DHE0023,DHE0024,DHE0025,DHE0026,DHE0027)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)


testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 2
edgeDataDefaults(testrabph, "color") <- rgb(192 / 256, 0, 0)


nodeDataDefaults(testrabph, "margin") <- 'unit(-5, "mm")'
nodeDataDefaults(testrabph, "shape") <- "circle"
nodeDataDefaults(testrabph, "color") <- "transparent"
plot(testrabph)


png(file = "mygraphic.png", width = 1920 * 3, height = 1080 * 3, res = 600)
plot(testrabph)
dev.off()