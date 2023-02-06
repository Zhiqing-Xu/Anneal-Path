#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)


DHE0000 <- DirectedHyperedge(c("C0006","C0004"), c("C0002","C0003"), "R0000")
DHE0001 <- DirectedHyperedge(c("C0008","C0009"), c("C0004","C0005"), "R0001")
DHE0002 <- DirectedHyperedge(c("C0011","C0009"), c("C0006","C0007"), "R0002")
DHE0003 <- DirectedHyperedge(c("C0000","C0001"), c("C0008","C0009"), "R0003")
DHE0004 <- DirectedHyperedge(c("C0000","C0001"), c("C0011","C0010"), "R0004")
DHE0005 <- DirectedHyperedge(c("C0014","C0004"), c("C0013","C0012"), "R0005")
DHE0006 <- DirectedHyperedge(c("C0011","C0010"), c("C0015","C0014"), "R0006")
DHE0007 <- DirectedHyperedge(c("C0005","C0019"), c("C0017","C0016"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0010","C0008"), c("C0019","C0018"), "R0008")
DHE0009 <- DirectedHyperedge(c("C0024","C0022"), c("C0020","C0021"), "R0009")
DHE0010 <- DirectedHyperedge(c("C0008","C0009"), c("C0022","C0023"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0010","C0009"), c("C0024","C0025"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0029","C0022"), c("C0026","C0027"), "R0012")
DHE0013 <- DirectedHyperedge(c("C0011","C0010"), c("C0028","C0029"), "R0013")
DHE0014 <- DirectedHyperedge(c("C0019","C0023"), c("C0031","C0030"), "R0014")
DHE0015 <- DirectedHyperedge(c("C0028","C0023"), c("C0033","C0032"), "R0015")
DHE0016 <- DirectedHyperedge(c("C0036","C0018"), c("C0035","C0034"), "R0016")
DHE0017 <- DirectedHyperedge(c("C0010","C0008"), c("C0037","C0036"), "R0017")
DHE0018 <- DirectedHyperedge(c("C0040","C0036"), c("C0039","C0038"), "R0018")
DHE0019 <- DirectedHyperedge(c("C0011","C0009"), c("C0040","C0041"), "R0019")
DHE0020 <- DirectedHyperedge(c("C0029","C0036"), c("C0042","C0043"), "R0020")
DHE0021 <- DirectedHyperedge(c("C0037","C0046"), c("C0044","C0045"), "R0021")
DHE0022 <- DirectedHyperedge(c("C0011","C0008"), c("C0046","C0047"), "R0022")
DHE0023 <- DirectedHyperedge(c("C0051","C0018"), c("C0048","C0049"), "R0023")
DHE0024 <- DirectedHyperedge(c("C0009"), c("C0051","C0050"), "R0024") ## F-arc
DHE0025 <- DirectedHyperedge(c("C0040","C0018"), c("C0053","C0052"), "R0025")
DHE0026 <- DirectedHyperedge(c("C0046","C0019"), c("C0055","C0054"), "R0026")
DHE0027 <- DirectedHyperedge(c("C0019","C0047"), c("C0057","C0056"), "R0027")
DHE0028 <- DirectedHyperedge(c("C0041","C0019"), c("C0059","C0058"), "R0028")
DHE0029 <- DirectedHyperedge(c("C0062","C0046"), c("C0060","C0061"), "R0029")
DHE0030 <- DirectedHyperedge(c("C0011","C0008"), c("C0062","C0063"), "R0030")
DHE0031 <- DirectedHyperedge(c("C0046","C0028"), c("C0064","C0065"), "R0031")
DHE0032 <- DirectedHyperedge(c("C0014","C0047"), c("C0066","C0067"), "R0032")
DHE0033 <- DirectedHyperedge(c("C0050","C0062"), c("C0068","C0069"), "R0033")
DHE0034 <- DirectedHyperedge(c("C0024","C0062"), c("C0071","C0070"), "R0034")
DHE0035 <- DirectedHyperedge(c("C0025","C0062"), c("C0073","C0072"), "R0035")
DHE0036 <- DirectedHyperedge(c("C0051","C0063"), c("C0075","C0074"), "R0036")
DHE0037 <- DirectedHyperedge(c("C0063"), c("C0076"), "R0037")# 
DHE0038 <- DirectedHyperedge(c("C0015","C0063"), c("C0079","C0078"), "R0038")
DHE0039 <- DirectedHyperedge(c("C0051","C0040"), c("C0080","C0081"), "R0039")
DHE0040 <- DirectedHyperedge(c("C0024","C0040"), c("C0082","C0083"), "R0040")
DHE0041 <- DirectedHyperedge(c("C0007","C0040"), c("C0084","C0085"), "R0041")
DHE0042 <- DirectedHyperedge(c("C0007","C0041"), c("C0086","C0087"), "R0042")
DHE0043 <- DirectedHyperedge(c("C0015","C0028"), c("C0088"), "R0043") # B-arc
DHE0044 <- DirectedHyperedge(c("C0014","C0029"), c("C0091","C0090"), "R0044")
DHE0045 <- DirectedHyperedge(c("C0015","C0014"), c("C0093","C0092"), "R0045")
Cnodes <- c("C0006","C0004","C0002","C0003","C0008","C0009","C0005","C0011","C0007","C0000","C0001","C0010","C0014","C0013","C0012","C0015","C0019","C0017","C0016","C0018","C0024","C0022","C0020","C0021","C0023","C0025","C0029","C0026","C0027","C0028","C0031","C0030","C0033","C0032","C0036","C0035","C0034","C0037","C0040","C0039","C0038","C0041","C0042","C0043","C0046","C0044","C0045","C0047","C0051","C0048","C0049","C0050","C0053","C0052","C0055","C0054","C0057","C0056","C0059","C0058","C0062","C0060","C0061","C0063","C0064","C0065","C0066","C0067","C0068","C0069","C0071","C0070","C0073","C0072","C0075","C0074","C0076","C0079","C0078","C0080","C0081","C0082","C0083","C0084","C0085","C0086","C0087","C0088","C0091","C0090","C0093","C0092")
Rnodes <- list(DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007,DHE0008,DHE0009,DHE0010,DHE0011,DHE0012,DHE0013,DHE0014,DHE0015,DHE0016,DHE0017,DHE0018,DHE0019,DHE0020,DHE0021,DHE0022,DHE0023,DHE0024,DHE0025,DHE0026,DHE0027,DHE0028,DHE0029,DHE0030,DHE0031,DHE0032,DHE0033,DHE0034,DHE0035,DHE0036,DHE0037,DHE0038,DHE0039,DHE0040,DHE0041,DHE0042,DHE0043,DHE0044,DHE0045)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)

#plot(hgbph)


testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 2
edgeDataDefaults(testrabph, "color") <- rgb(80 / 256, 51 / 256, 0 / 256)


nodeDataDefaults(testrabph, "margin") <- 'unit(-5, "mm")'
nodeDataDefaults(testrabph, "shape") <- "circle"
nodeDataDefaults(testrabph, "color") <- "transparent"


edgeData(testrabph, c("C0015", "R0043"), c("R0043", "C0088"), "lwd") <- c("4", "4")
edgeData(testrabph, c("C0015", "R0043"), c("R0043", "C0088"), "color") <- rgb(112 / 256, 48 / 256, 160 / 256)
edgeData(testrabph, c("C0028", "R0043"), c("R0043", "C0088"), "lwd") <- c("4", "4")
edgeData(testrabph, c("C0028", "R0043"), c("R0043", "C0088"), "color") <- rgb(112 / 256, 48 / 256, 160 / 256)

edgeData(testrabph, c("C0009", "R0024"), c("R0024", "C0051"), "lwd") <- c("4", "4")
edgeData(testrabph, c("C0009", "R0024"), c("R0024", "C0051"), "color") <- rgb(85 / 256, 142 / 256, 213 / 256)
edgeData(testrabph, c("C0009", "R0024"), c("R0024", "C0050"), "lwd") <- c("4", "4")
edgeData(testrabph, c("C0009", "R0024"), c("R0024", "C0050"), "color") <- rgb(85 / 256, 142 / 256, 213 / 256)


edgeData(testrabph, c("C0063", "R0037"), c("R0037", "C0076"), "lwd") <- c("4", "4")
edgeData(testrabph, c("C0063", "R0037"), c("R0037", "C0076"), "color") <- rgb(192 / 256,0,0)




plot(testrabph)


png(file = "mygraphic.png", width = 1920 * 3, height = 1080 * 3, res = 600)
plot(testrabph)
dev.off()

