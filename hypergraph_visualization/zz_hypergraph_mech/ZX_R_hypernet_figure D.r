#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

DHEXXXX <- DirectedHyperedge(c("G6P"), c("F6P"), "RXXXX")
#DHEYYYY <- DirectedHyperedge(c("DAHP"), c("YYY"), "RYYYY")
DHE0000 <- DirectedHyperedge(c("E4P","PEP"), c("DAHP"), "EC2.5.1")
DHE0001 <- DirectedHyperedge(c("2PG"), c("PEP"), "EC4.2.1")
DHE0002 <- DirectedHyperedge(c("3PG"), c("2PG"), "EC5.4.2")
DHE0003 <- DirectedHyperedge(c("GA3P","F6P"), c("E4P","X5P"), "EC2.2.1")
DHE0004 <- DirectedHyperedge(c("GA3P"), c("3PG"), "EC1.2.1")
DHE0005 <- DirectedHyperedge(c("F16bP"), c("GA3P","DHAP"), "EC4.1.2")
DHE0006 <- DirectedHyperedge(c("F6P"), c("F16bP"), "EC2.7.1")

DHE0007 <- DirectedHyperedge(c("E4P","C0010"), c("DAHP"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0011"), c("C0010"), "R0008")
DHE0009 <- DirectedHyperedge(c("C0012"), c("E4P"), "R0009")
DHE0010 <- DirectedHyperedge(c("C0014"), c("C0011"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0013"), c("C0012"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0016"), c("C0013"), "R0012")
DHE0013 <- DirectedHyperedge(c("C0017"), c("C0014","E4P"), "R0013")
DHE0014 <- DirectedHyperedge(c("F16bP"), c("C0015","C0016"), "R0014")
DHE0015 <- DirectedHyperedge(c("F6P"), c("C0017"), "R0015")
DHE0016 <- DirectedHyperedge(c("C0019"), c("C0016","C0018"), "R0016")
DHE0017 <- DirectedHyperedge(c("F6P"), c("C0019"), "R0017")
DHE0018 <- DirectedHyperedge(c("C0021"), c("C0016","C0020"), "R0018")
DHE0019 <- DirectedHyperedge(c("F6P"), c("C0021"), "R0019")
DHE0020 <- DirectedHyperedge(c("C0023"), c("C0016","C0022"), "R0020")
DHE0021 <- DirectedHyperedge(c("F6P"), c("C0023"), "R0021")
DHE0022 <- DirectedHyperedge(c("C0025"), c("C0024","C0014"), "R0022")
DHE0023 <- DirectedHyperedge(c("C0017"), c("C0025"), "R0023")
DHE0024 <- DirectedHyperedge(c("F6P"), c("C0026","C0016"), "R0024")
DHE0025 <- DirectedHyperedge(c("C0027"), c("C0025"), "R0025")
DHE0026 <- DirectedHyperedge(c("F6P"), c("C0027"), "R0026")
DHE0027 <- DirectedHyperedge(c("C0029"), c("C0014","C0028"), "R0027")
DHE0028 <- DirectedHyperedge(c("C0017"), c("C0029"), "R0028")
DHE0029 <- DirectedHyperedge(c("C0030"), c("C0029"), "R0029")
DHE0030 <- DirectedHyperedge(c("F6P"), c("C0030"), "R0030")
DHE0031 <- DirectedHyperedge(c("C0032"), c("C0014","C0031"), "R0031")
DHE0032 <- DirectedHyperedge(c("C0017"), c("C0032"), "R0032")
DHE0033 <- DirectedHyperedge(c("C0033"), c("C0032"), "R0033")
DHE0034 <- DirectedHyperedge(c("F6P"), c("C0033"), "R0034")
DHE0035 <- DirectedHyperedge(c("C0035"), c("C0014","C0034"), "R0035")
DHE0036 <- DirectedHyperedge(c("C0036"), c("C0035"), "R0036")
DHE0037 <- DirectedHyperedge(c("F6P"), c("C0036"), "R0037")
DHE0038 <- DirectedHyperedge(c("C0017"), c("C0035"), "R0038")
DHE0039 <- DirectedHyperedge(c("C0037"), c("C0010"), "R0039")
DHE0040 <- DirectedHyperedge(c("3PG"), c("C0037"), "R0040")
DHE0041 <- DirectedHyperedge(c("C0029"), c("C0038","3PG"), "R0041")
DHE0042 <- DirectedHyperedge(c("C0037"), c("PEP"), "R0042")
DHE0043 <- DirectedHyperedge(c("C0039"), c("PEP"), "R0043")
DHE0044 <- DirectedHyperedge(c("3PG"), c("C0039"), "R0044")
Cnodes <- c("G6P","2PG","C0036","C0034","PEP","C0015","C0014","C0017","C0016","C0011","C0010","C0013","C0012","DAHP","C0039","C0038","C0019","C0018","3PG","E4P","X5P","C0035","F16bP","C0030","C0028","C0029","C0024","C0025","C0026","C0027","C0020","C0021","C0022","C0023","C0033","C0032","C0037","C0031","DHAP","GA3P","F6P")
Rnodes <- list(DHEXXXX, DHE0007, DHE0008, DHE0009, DHE0010, DHE0011, DHE0012, DHE0013, DHE0014, DHE0015, DHE0016, DHE0017, DHE0018, DHE0019, DHE0020, DHE0021, DHE0022, DHE0023, DHE0024, DHE0025, DHE0026, DHE0027, DHE0028, DHE0029, DHE0030, DHE0031, DHE0032, DHE0033, DHE0005, DHE0001, DHE0002, DHE0003, DHE0004, DHE0000, DHE0006,DHE0034, DHE0035, DHE0036, DHE0037, DHE0038, DHE0039, DHE0040, DHE0041, DHE0042, DHE0043, DHE0044)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)

testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 1.5
edgeDataDefaults(testrabph, "color") <- "gray60"

edgeData(testrabph, c("F6P", "EC2.7.1"), c("EC2.7.1", "F16bP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6P", "EC2.7.1"), c("EC2.7.1", "F16bP"), "color") <- "firebrick"


edgeData(testrabph, c("F16bP", "EC4.1.2"), c("EC4.1.2", "GA3P"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F16bP", "EC4.1.2"), c("EC4.1.2", "DHAP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F16bP", "EC4.1.2"), c("EC4.1.2", "GA3P"), "color") <- rgb(85 / 256, 142 / 256, 213 / 256)
edgeData(testrabph, c("F16bP", "EC4.1.2"), c("EC4.1.2", "DHAP"), "color") <- rgb(85 / 256, 142 / 256, 213 / 256)

edgeData(testrabph, c("F6P", "EC2.2.1"), c("EC2.2.1", "X5P"), "lwd") <- c("3", "3")
edgeData(testrabph, c("GA3P", "EC2.2.1"), c("EC2.2.1", "E4P"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6P", "EC2.2.1"), c("EC2.2.1", "X5P"), "color") <- rgb(102 / 256, 51 / 256, 0 / 256)
edgeData(testrabph, c("GA3P", "EC2.2.1"), c("EC2.2.1", "E4P"), "color") <- rgb(102 / 256, 51 / 256, 0 / 256)

edgeData(testrabph, c("GA3P", "EC1.2.1"), c("EC1.2.1", "3PG"), "lwd") <- c("3", "3")
edgeData(testrabph, c("GA3P", "EC1.2.1"), c("EC1.2.1", "3PG"), "color") <- "firebrick"

edgeData(testrabph, c("3PG", "EC5.4.2"), c("EC5.4.2", "2PG"), "lwd") <- c("3", "3")
edgeData(testrabph, c("3PG", "EC5.4.2"), c("EC5.4.2", "2PG"), "color") <- "firebrick"

edgeData(testrabph, c("2PG", "EC4.2.1"), c("EC4.2.1", "PEP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("2PG", "EC4.2.1"), c("EC4.2.1", "PEP"), "color") <- "firebrick"

edgeData(testrabph, c("PEP", "EC2.5.1"), c("EC2.5.1", "DAHP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("PEP", "EC2.5.1"), c("EC2.5.1", "DAHP"), "color") <- rgb(112 / 256, 48 / 256, 160 / 256)
edgeData(testrabph, c("E4P", "EC2.5.1"), c("EC2.5.1", "DAHP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("E4P", "EC2.5.1"), c("EC2.5.1", "DAHP"), "color") <- rgb(112 / 256, 48 / 256, 160 / 256)

nodeDataDefaults(testrabph, "margin") <- 'unit(-3, "mm")'
nodeDataDefaults(testrabph, "shape") <- "circle"
nodeDataDefaults(testrabph, "color") <- "transparent"


nodeData(testrabph, c("F6P"), "color") <- "white"
nodeData(testrabph, c("F6P"), "margin") <- 'unit(2, "mm")'

nodeData(testrabph, c("F16bP"), "color") <- "black"
nodeData(testrabph, c("F16bP"), "margin") <- 'unit(0.5, "mm")'

nodeData(testrabph, c("DHAP"), "color") <- "black"
nodeData(testrabph, c("DHAP"), "margin") <- 'unit(0.5, "mm")'
nodeData(testrabph, c("GA3P"), "color") <- "black"
nodeData(testrabph, c("GA3P"), "margin") <- 'unit(1, "mm")'

nodeData(testrabph, c("X5P"), "color") <- "black"
nodeData(testrabph, c("X5P"), "margin") <- 'unit(2, "mm")'
nodeData(testrabph, c("E4P"), "color") <- "black"
nodeData(testrabph, c("E4P"), "margin") <- 'unit(2, "mm")'
nodeData(testrabph, c("3PG"), "color") <- "black"
nodeData(testrabph, c("3PG"), "margin") <- 'unit(2, "mm")'

nodeData(testrabph, c("2PG"), "color") <- "black"
nodeData(testrabph, c("2PG"), "margin") <- 'unit(2, "mm")'

nodeData(testrabph, c("PEP"), "color") <- "black"
nodeData(testrabph, c("PEP"), "margin") <- 'unit(2, "mm")'

nodeData(testrabph, c("DAHP"), "color") <- "black"
nodeData(testrabph, c("DAHP"), "margin") <- 'unit(0.5, "mm")'

nodeData(testrabph, c("EC2.7.1"), "color") <- "black"
nodeData(testrabph, c("EC4.1.2"), "color") <- "black"
nodeData(testrabph, c("EC2.2.1"), "color") <- "black"
nodeData(testrabph, c("EC1.2.1"), "color") <- "black"
nodeData(testrabph, c("EC5.4.2"), "color") <- "black"
nodeData(testrabph, c("EC4.2.1"), "color") <- "black"
nodeData(testrabph, c("EC2.5.1"), "color") <- "black"
#nodeData(testrabph, c("DAHP"), "margin") <- 'unit(1, "mm")'



plot(testrabph)


# nodeData(testrabph, c("AAA"), "fontsize") <- '50'
# testrabph@graph@AgNode[[8]]@color<-"pink"
# nodeData(testrabph, c("AAA"), "fillcolor") <- "blue"
# dev.new(width = 1920 * 2, height = 1920 * 2, unit = "px")
# pdf("mygraph.pdf", width = 11, height = 8)



png(file = "FigureD.png", width = 1080*4, height = 1920*4, res = 600)
plot(testrabph)
dev.off()