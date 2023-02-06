#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)


DHE0000 <- DirectedHyperedge(c("C0002"), c("C0001"), "R0000")
DHE0001 <- DirectedHyperedge(c("C0003"), c("C0002"), "R0001")
DHE0002 <- DirectedHyperedge(c("C0004"), c("C0003"), "R0002")
DHE0003 <- DirectedHyperedge(c("C0005"), c("C0004"), "R0003")
DHE0004 <- DirectedHyperedge(c("C0006"), c("C0005"), "R0004")
DHE0005 <- DirectedHyperedge(c("C0007"), c("C0006"), "R0005")
DHE0006 <- DirectedHyperedge(c("C0008"), c("C0007"), "R0006")
DHE0007 <- DirectedHyperedge(c("C0000"), c("C0008"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0009"), c("C0007"), "R0008")
DHE0009 <- DirectedHyperedge(c("C0000"), c("C0009"), "R0009")
DHE0010 <- DirectedHyperedge(c("C0010"), c("C0002"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0011"), c("C0010"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0005"), c("C0011"), "R0012")
DHE0013 <- DirectedHyperedge(c("C0012"), c("C0001"), "R0013")
DHE0014 <- DirectedHyperedge(c("C0003"), c("C0012"), "R0014")
DHE0015 <- DirectedHyperedge(c("C0013"), c("C0001"), "R0015")
DHE0016 <- DirectedHyperedge(c("C0010"), c("C0013"), "R0016")
DHE0017 <- DirectedHyperedge(c("C0014"), c("C0013"), "R0017")
DHE0018 <- DirectedHyperedge(c("C0011"), c("C0014"), "R0018")
Cnodes <- c("C0002","C0001","C0003","C0004","C0005","C0006","C0007","C0008","C0000","C0009","C0010","C0011","C0012","C0013","C0014")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007,DHE0008,DHE0009,DHE0010,DHE0011,DHE0012,DHE0013,DHE0014,DHE0015,DHE0016,DHE0017,DHE0018)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)



