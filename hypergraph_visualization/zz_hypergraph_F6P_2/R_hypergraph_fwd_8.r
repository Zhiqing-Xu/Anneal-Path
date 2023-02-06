#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(hypergraph)
library(hyperdraw)

DHEXXXX <- DirectedHyperedge(c("G6PXX"), c("F6PXX"), "RXXXX")
DHE0000 <- DirectedHyperedge(c("C0004","C0003"), c("C0001"), "R0000")
DHE0001 <- DirectedHyperedge(c("C0005"), c("C0002","C0003"), "R0001")
DHE0002 <- DirectedHyperedge(c("C0006"), c("C0004"), "R0002")
DHE0003 <- DirectedHyperedge(c("F6PXX"), c("C0005"), "R0003")
DHE0004 <- DirectedHyperedge(c("F6PXX"), c("C0006","C0007"), "R0004")
DHE0005 <- DirectedHyperedge(c("C0009"), c("C0004","C0008"), "R0005")
DHE0006 <- DirectedHyperedge(c("F6PXX"), c("C0009"), "R0006")
DHE0007 <- DirectedHyperedge(c("C0012","C0013"), c("C0010"), "R0007")
DHE0008 <- DirectedHyperedge(c("C0014"), c("C0012","C0011"), "R0008")
DHE0009 <- DirectedHyperedge(c("F6PXX"), c("C0013"), "R0009")
DHE0010 <- DirectedHyperedge(c("F6PXX"), c("C0014"), "R0010")
DHE0011 <- DirectedHyperedge(c("C0004","C0016"), c("C0015"), "R0011")
DHE0012 <- DirectedHyperedge(c("C0017"), c("C0016"), "R0012")
DHE0013 <- DirectedHyperedge(c("F6PXX"), c("C0017"), "R0013")
DHE0014 <- DirectedHyperedge(c("C0018"), c("C0016"), "R0014")
DHE0015 <- DirectedHyperedge(c("F6PXX"), c("C0018"), "R0015")
DHE0016 <- DirectedHyperedge(c("C0004","C0006"), c("C0019"), "R0016")
DHE0017 <- DirectedHyperedge(c("C0021"), c("C0006","C0020"), "R0017")
DHE0018 <- DirectedHyperedge(c("F6PXX"), c("C0021"), "R0018")
DHE0019 <- DirectedHyperedge(c("C0021"), c("C0006","C0022"), "R0019")
DHE0020 <- DirectedHyperedge(c("C0024"), c("C0023","C0006"), "R0020")
DHE0021 <- DirectedHyperedge(c("F6PXX"), c("C0024"), "R0021")
DHE0022 <- DirectedHyperedge(c("C0026"), c("C0006","C0025"), "R0022")
DHE0023 <- DirectedHyperedge(c("F6PXX"), c("C0026"), "R0023")
DHE0024 <- DirectedHyperedge(c("C0018"), c("C0006","C0027"), "R0024")
DHE0025 <- DirectedHyperedge(c("C0029"), c("C0028"), "R0025")
DHE0026 <- DirectedHyperedge(c("C0017"), c("C0029"), "R0026")
DHE0027 <- DirectedHyperedge(c("C0030"), c("C0029"), "R0027")
DHE0028 <- DirectedHyperedge(c("F6PXX"), c("C0030"), "R0028")
DHE0029 <- DirectedHyperedge(c("C0032"), c("C0031"), "R0029")
DHE0030 <- DirectedHyperedge(c("C0006","C0007"), c("C0032"), "R0030")
DHE0031 <- DirectedHyperedge(c("C0012","C0034"), c("C0033"), "R0031")
DHE0032 <- DirectedHyperedge(c("F16bP"), c("C0034","C0007"), "R0032")
DHE0033 <- DirectedHyperedge(c("F6PXX"), c("F16bP"), "R0033")
DHE0034 <- DirectedHyperedge(c("C0006"), c("C0034"), "R0034")
DHE0035 <- DirectedHyperedge(c("C0037"), c("C0036"), "R0035")
DHE0036 <- DirectedHyperedge(c("C0017"), c("C0037"), "R0036")
DHE0037 <- DirectedHyperedge(c("C0009"), c("C0037"), "R0037")
DHE0038 <- DirectedHyperedge(c("C0012","C0018"), c("C0038"), "R0038")
DHE0039 <- DirectedHyperedge(c("C0026"), c("C0018","C0039"), "R0039")
DHE0040 <- DirectedHyperedge(c("C0042","C0012"), c("C0040"), "R0040")
DHE0041 <- DirectedHyperedge(c("C0043"), c("C0041","C0042"), "R0041")
DHE0042 <- DirectedHyperedge(c("F6PXX"), c("C0043"), "R0042")
DHE0043 <- DirectedHyperedge(c("C0045"), c("C0044"), "R0043")
DHE0044 <- DirectedHyperedge(c("C0026"), c("C0045"), "R0044")
DHE0045 <- DirectedHyperedge(c("C0017"), c("C0045"), "R0045")
DHE0046 <- DirectedHyperedge(c("C0046"), c("C0044"), "R0046")
DHE0047 <- DirectedHyperedge(c("C0017"), c("C0046"), "R0047")
DHE0048 <- DirectedHyperedge(c("C0047"), c("C0046"), "R0048")
DHE0049 <- DirectedHyperedge(c("F6PXX"), c("C0047"), "R0049")
DHE0050 <- DirectedHyperedge(c("C0045"), c("C0048"), "R0050")
DHE0051 <- DirectedHyperedge(c("C0050"), c("C0049"), "R0051")
DHE0052 <- DirectedHyperedge(c("C0017"), c("C0050"), "R0052")
DHE0053 <- DirectedHyperedge(c("C0051"), c("C0050"), "R0053")
DHE0054 <- DirectedHyperedge(c("F6PXX"), c("C0051"), "R0054")
DHE0055 <- DirectedHyperedge(c("C0027","C0047"), c("C0052"), "R0055")
DHE0056 <- DirectedHyperedge(c("C0053"), c("C0047"), "R0056")
DHE0057 <- DirectedHyperedge(c("C0007"), c("C0027"), "R0057")
DHE0058 <- DirectedHyperedge(c("F6PXX"), c("C0053"), "R0058")
DHE0059 <- DirectedHyperedge(c("C0034","C0055"), c("C0054"), "R0059")
DHE0060 <- DirectedHyperedge(c("F6PXX"), c("C0055"), "R0060")
DHE0061 <- DirectedHyperedge(c("C0042","C0013"), c("C0056"), "R0061")
DHE0062 <- DirectedHyperedge(c("C0045"), c("C0057"), "R0062")
DHE0063 <- DirectedHyperedge(c("C0012","C0006"), c("C0058"), "R0063")
DHE0064 <- DirectedHyperedge(c("C0060","C0027"), c("C0059"), "R0064")
DHE0065 <- DirectedHyperedge(c("F6PXX"), c("C0060"), "R0065")
DHE0066 <- DirectedHyperedge(c("C0062","C0006"), c("C0061"), "R0066")
DHE0067 <- DirectedHyperedge(c("C0017"), c("C0062"), "R0067")
DHE0068 <- DirectedHyperedge(c("C0014"), c("C0062"), "R0068")
DHE0069 <- DirectedHyperedge(c("C0012"), c("C0063"), "R0069")
DHE0070 <- DirectedHyperedge(c("C0012"), c("C0064"), "R0070")
DHE0071 <- DirectedHyperedge(c("C0067"), c("C0066","E4P"), "R0071")
DHE0072 <- DirectedHyperedge(c("C0017"), c("C0067"), "R0072")
DHE0073 <- DirectedHyperedge(c("F16bP"), c("C0067"), "R0073")
DHE0074 <- DirectedHyperedge(c("C0069","C0006"), c("C0068"), "R0074")
DHE0075 <- DirectedHyperedge(c("C0018"), c("C0069"), "R0075")
DHE0076 <- DirectedHyperedge(c("C0051"), c("C0069"), "R0076")
DHE0077 <- DirectedHyperedge(c("C0062"), c("C0070"), "R0077")
DHE0078 <- DirectedHyperedge(c("C0006","C0034"), c("C0071"), "R0078")
DHE0079 <- DirectedHyperedge(c("C0073"), c("C0072"), "R0079")
DHE0080 <- DirectedHyperedge(c("C0006","C0007"), c("C0073"), "R0080")
DHE0081 <- DirectedHyperedge(c("C0027","C0013"), c("C0074"), "R0081")
DHE0082 <- DirectedHyperedge(c("C0006","C0076"), c("C0075"), "R0082")
DHE0083 <- DirectedHyperedge(c("C0009"), c("C0076"), "R0083")
DHE0084 <- DirectedHyperedge(c("C0018"), c("C0076"), "R0084")
DHE0085 <- DirectedHyperedge(c("C0067"), c("C0077"), "R0085")
DHE0086 <- DirectedHyperedge(c("C0007"), c("C0078"), "R0086")
DHE0087 <- DirectedHyperedge(c("C0027","C0079"), c("C0078"), "R0087")
DHE0088 <- DirectedHyperedge(c("C0007"), c("C0079"), "R0088")
DHE0089 <- DirectedHyperedge(c("C0051"), c("C0080","C0079"), "R0089")
DHE0090 <- DirectedHyperedge(c("C0030"), c("C0007","C0081"), "R0090")
DHE0091 <- DirectedHyperedge(c("C0060"), c("C0082","C0007"), "R0091")
DHE0092 <- DirectedHyperedge(c("C0055"), c("C0083","C0007"), "R0092")
DHE0093 <- DirectedHyperedge(c("C0032"), c("C0084"), "R0093")
DHE0094 <- DirectedHyperedge(c("C0086"), c("C0085"), "R0094")
DHE0095 <- DirectedHyperedge(c("F16bP"), c("C0086"), "R0095")
DHE0096 <- DirectedHyperedge(c("C0018"), c("C0086"), "R0096")
Cnodes <- c("G6PXX", "C0004","C0003","C0001","C0005","C0002","C0006","F6PXX","C0007","C0009","C0008","C0012","C0013","C0010","C0014","C0011","C0016","C0015","C0017","C0018","C0019","C0021","C0020","C0022","C0024","C0023","C0026","C0025","C0027","C0029","C0028","C0030","C0032","C0031","C0034","C0033","F16bP","C0037","C0036","C0038","C0039","C0042","C0040","C0043","C0041","C0045","C0044","C0046","C0047","C0048","C0050","C0049","C0051","C0052","C0053","C0055","C0054","C0056","C0057","C0058","C0060","C0059","C0062","C0061","C0063","C0064","C0067","C0066","E4P","C0069","C0068","C0070","C0071","C0073","C0072","C0074","C0076","C0075","C0077","C0078","C0079","C0080","C0081","C0082","C0083","C0084","C0086","C0085")
Rnodes <- list(DHEXXXX, DHE0000,DHE0001,DHE0002,DHE0003,DHE0004,DHE0005,DHE0006,DHE0007,DHE0008,DHE0009,DHE0010,DHE0011,DHE0012,DHE0013,DHE0014,DHE0015,DHE0016,DHE0017,DHE0018,DHE0019,DHE0020,DHE0021,DHE0022,DHE0023,DHE0024,DHE0025,DHE0026,DHE0027,DHE0028,DHE0029,DHE0030,DHE0031,DHE0032,DHE0033,DHE0034,DHE0035,DHE0036,DHE0037,DHE0038,DHE0039,DHE0040,DHE0041,DHE0042,DHE0043,DHE0044,DHE0045,DHE0046,DHE0047,DHE0048,DHE0049,DHE0050,DHE0051,DHE0052,DHE0053,DHE0054,DHE0055,DHE0056,DHE0057,DHE0058,DHE0059,DHE0060,DHE0061,DHE0062,DHE0063,DHE0064,DHE0065,DHE0066,DHE0067,DHE0068,DHE0069,DHE0070,DHE0071,DHE0072,DHE0073,DHE0074,DHE0075,DHE0076,DHE0077,DHE0078,DHE0079,DHE0080,DHE0081,DHE0082,DHE0083,DHE0084,DHE0085,DHE0086,DHE0087,DHE0088,DHE0089,DHE0090,DHE0091,DHE0092,DHE0093,DHE0094,DHE0095,DHE0096)
hg <- Hypergraph(Cnodes,Rnodes)
hgbph <- graphBPH(hg)
#plot(hgbph)


testrabph <- graphLayout(hgbph)
edgeDataDefaults(testrabph, "lwd") <- 1
edgeDataDefaults(testrabph, "color") <- "black"


nodeDataDefaults(testrabph, "margin") <- 'unit(-4, "mm")'
nodeDataDefaults(testrabph, "shape") <- "circle"
nodeDataDefaults(testrabph, "color") <- "transparent"







edgeData(testrabph, c("C0004", "R0000"), c("R0000", "C0001"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0004", "R0000"), c("R0000", "C0001"), "color") <- "blue"

edgeData(testrabph, c("C0003", "R0000"), c("R0000", "C0001"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0003", "R0000"), c("R0000", "C0001"), "color") <- "blue"

edgeData(testrabph, c("C0005", "R0001"), c("R0001", "C0002"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0005", "R0001"), c("R0001", "C0002"), "color") <- "blue"

edgeData(testrabph, c("C0005", "R0001"), c("R0001", "C0003"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0005", "R0001"), c("R0001", "C0003"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0002"), c("R0002", "C0004"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0002"), c("R0002", "C0004"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0003"), c("R0003", "C0005"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0003"), c("R0003", "C0005"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0004"), c("R0004", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0004"), c("R0004", "C0006"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0004"), c("R0004", "C0007"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0004"), c("R0004", "C0007"), "color") <- "blue"

edgeData(testrabph, c("C0009", "R0005"), c("R0005", "C0004"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0009", "R0005"), c("R0005", "C0004"), "color") <- "blue"

edgeData(testrabph, c("C0009", "R0005"), c("R0005", "C0008"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0009", "R0005"), c("R0005", "C0008"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0006"), c("R0006", "C0009"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0006"), c("R0006", "C0009"), "color") <- "firebrick"

edgeData(testrabph, c("C0012", "R0007"), c("R0007", "C0010"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0007"), c("R0007", "C0010"), "color") <- "blue"

edgeData(testrabph, c("C0013", "R0007"), c("R0007", "C0010"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0013", "R0007"), c("R0007", "C0010"), "color") <- "blue"

edgeData(testrabph, c("C0014", "R0008"), c("R0008", "C0012"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0014", "R0008"), c("R0008", "C0012"), "color") <- "blue"

edgeData(testrabph, c("C0014", "R0008"), c("R0008", "C0011"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0014", "R0008"), c("R0008", "C0011"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0009"), c("R0009", "C0013"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0009"), c("R0009", "C0013"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0010"), c("R0010", "C0014"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0010"), c("R0010", "C0014"), "color") <- "firebrick"

edgeData(testrabph, c("C0004", "R0011"), c("R0011", "C0015"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0004", "R0011"), c("R0011", "C0015"), "color") <- "blue"

edgeData(testrabph, c("C0016", "R0011"), c("R0011", "C0015"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0016", "R0011"), c("R0011", "C0015"), "color") <- "blue"

edgeData(testrabph, c("C0017", "R0012"), c("R0012", "C0016"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0017", "R0012"), c("R0012", "C0016"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0013"), c("R0013", "C0017"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0013"), c("R0013", "C0017"), "color") <- "firebrick"

edgeData(testrabph, c("C0018", "R0014"), c("R0014", "C0016"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0014"), c("R0014", "C0016"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0015"), c("R0015", "C0018"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0015"), c("R0015", "C0018"), "color") <- "firebrick"

edgeData(testrabph, c("C0004", "R0016"), c("R0016", "C0019"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0004", "R0016"), c("R0016", "C0019"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0016"), c("R0016", "C0019"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0016"), c("R0016", "C0019"), "color") <- "blue"

edgeData(testrabph, c("C0021", "R0017"), c("R0017", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0021", "R0017"), c("R0017", "C0006"), "color") <- "blue"

edgeData(testrabph, c("C0021", "R0017"), c("R0017", "C0020"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0021", "R0017"), c("R0017", "C0020"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0018"), c("R0018", "C0021"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0018"), c("R0018", "C0021"), "color") <- "blue"

edgeData(testrabph, c("C0021", "R0019"), c("R0019", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0021", "R0019"), c("R0019", "C0006"), "color") <- "blue"

edgeData(testrabph, c("C0021", "R0019"), c("R0019", "C0022"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0021", "R0019"), c("R0019", "C0022"), "color") <- "blue"

edgeData(testrabph, c("C0024", "R0020"), c("R0020", "C0023"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0024", "R0020"), c("R0020", "C0023"), "color") <- "blue"

edgeData(testrabph, c("C0024", "R0020"), c("R0020", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0024", "R0020"), c("R0020", "C0006"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0021"), c("R0021", "C0024"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0021"), c("R0021", "C0024"), "color") <- "blue"

edgeData(testrabph, c("C0026", "R0022"), c("R0022", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0026", "R0022"), c("R0022", "C0006"), "color") <- "blue"

edgeData(testrabph, c("C0026", "R0022"), c("R0022", "C0025"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0026", "R0022"), c("R0022", "C0025"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0023"), c("R0023", "C0026"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0023"), c("R0023", "C0026"), "color") <- "firebrick"

edgeData(testrabph, c("C0018", "R0024"), c("R0024", "C0006"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0024"), c("R0024", "C0006"), "color") <- "blue"

edgeData(testrabph, c("C0018", "R0024"), c("R0024", "C0027"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0024"), c("R0024", "C0027"), "color") <- "blue"

edgeData(testrabph, c("C0029", "R0025"), c("R0025", "C0028"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0029", "R0025"), c("R0025", "C0028"), "color") <- "firebrick"

edgeData(testrabph, c("C0017", "R0026"), c("R0026", "C0029"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0026"), c("R0026", "C0029"), "color") <- "firebrick"

edgeData(testrabph, c("C0030", "R0027"), c("R0027", "C0029"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0030", "R0027"), c("R0027", "C0029"), "color") <- "firebrick"

edgeData(testrabph, c("F6PXX", "R0028"), c("R0028", "C0030"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0028"), c("R0028", "C0030"), "color") <- "firebrick"

edgeData(testrabph, c("C0032", "R0029"), c("R0029", "C0031"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0032", "R0029"), c("R0029", "C0031"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0030"), c("R0030", "C0032"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0030"), c("R0030", "C0032"), "color") <- "blue"

edgeData(testrabph, c("C0007", "R0030"), c("R0030", "C0032"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0007", "R0030"), c("R0030", "C0032"), "color") <- "blue"

edgeData(testrabph, c("C0012", "R0031"), c("R0031", "C0033"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0031"), c("R0031", "C0033"), "color") <- "blue"

edgeData(testrabph, c("C0034", "R0031"), c("R0031", "C0033"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0034", "R0031"), c("R0031", "C0033"), "color") <- "blue"

edgeData(testrabph, c("F16bP", "R0032"), c("R0032", "C0034"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F16bP", "R0032"), c("R0032", "C0034"), "color") <- "blue"

edgeData(testrabph, c("F16bP", "R0032"), c("R0032", "C0007"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F16bP", "R0032"), c("R0032", "C0007"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0033"), c("R0033", "F16bP"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0033"), c("R0033", "F16bP"), "color") <- "firebrick"

edgeData(testrabph, c("C0006", "R0034"), c("R0034", "C0034"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0034"), c("R0034", "C0034"), "color") <- "blue"

edgeData(testrabph, c("C0037", "R0035"), c("R0035", "C0036"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0037", "R0035"), c("R0035", "C0036"), "color") <- "firebrick"

edgeData(testrabph, c("C0017", "R0036"), c("R0036", "C0037"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0036"), c("R0036", "C0037"), "color") <- "firebrick"

edgeData(testrabph, c("C0009", "R0037"), c("R0037", "C0037"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0009", "R0037"), c("R0037", "C0037"), "color") <- "firebrick"

edgeData(testrabph, c("C0012", "R0038"), c("R0038", "C0038"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0038"), c("R0038", "C0038"), "color") <- "blue"

edgeData(testrabph, c("C0018", "R0038"), c("R0038", "C0038"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0038"), c("R0038", "C0038"), "color") <- "blue"

edgeData(testrabph, c("C0026", "R0039"), c("R0039", "C0018"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0026", "R0039"), c("R0039", "C0018"), "color") <- "blue"

edgeData(testrabph, c("C0026", "R0039"), c("R0039", "C0039"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0026", "R0039"), c("R0039", "C0039"), "color") <- "blue"

edgeData(testrabph, c("C0042", "R0040"), c("R0040", "C0040"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0042", "R0040"), c("R0040", "C0040"), "color") <- "blue"

edgeData(testrabph, c("C0012", "R0040"), c("R0040", "C0040"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0040"), c("R0040", "C0040"), "color") <- "blue"

edgeData(testrabph, c("C0043", "R0041"), c("R0041", "C0041"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0043", "R0041"), c("R0041", "C0041"), "color") <- "blue"

edgeData(testrabph, c("C0043", "R0041"), c("R0041", "C0042"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0043", "R0041"), c("R0041", "C0042"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0042"), c("R0042", "C0043"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0042"), c("R0042", "C0043"), "color") <- "blue"

edgeData(testrabph, c("C0045", "R0043"), c("R0043", "C0044"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0045", "R0043"), c("R0043", "C0044"), "color") <- "firebrick"

edgeData(testrabph, c("C0026", "R0044"), c("R0044", "C0045"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0026", "R0044"), c("R0044", "C0045"), "color") <- "firebrick"

edgeData(testrabph, c("C0017", "R0045"), c("R0045", "C0045"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0045"), c("R0045", "C0045"), "color") <- "firebrick"

edgeData(testrabph, c("C0046", "R0046"), c("R0046", "C0044"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0046", "R0046"), c("R0046", "C0044"), "color") <- "firebrick"

edgeData(testrabph, c("C0017", "R0047"), c("R0047", "C0046"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0047"), c("R0047", "C0046"), "color") <- "firebrick"

edgeData(testrabph, c("C0047", "R0048"), c("R0048", "C0046"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0047", "R0048"), c("R0048", "C0046"), "color") <- "firebrick"

edgeData(testrabph, c("F6PXX", "R0049"), c("R0049", "C0047"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0049"), c("R0049", "C0047"), "color") <- "firebrick"

edgeData(testrabph, c("C0045", "R0050"), c("R0050", "C0048"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0045", "R0050"), c("R0050", "C0048"), "color") <- "firebrick"

edgeData(testrabph, c("C0050", "R0051"), c("R0051", "C0049"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0050", "R0051"), c("R0051", "C0049"), "color") <- "firebrick"

edgeData(testrabph, c("C0017", "R0052"), c("R0052", "C0050"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0052"), c("R0052", "C0050"), "color") <- "firebrick"

edgeData(testrabph, c("C0051", "R0053"), c("R0053", "C0050"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0051", "R0053"), c("R0053", "C0050"), "color") <- "firebrick"

edgeData(testrabph, c("F6PXX", "R0054"), c("R0054", "C0051"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F6PXX", "R0054"), c("R0054", "C0051"), "color") <- "firebrick"

edgeData(testrabph, c("C0027", "R0055"), c("R0055", "C0052"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0027", "R0055"), c("R0055", "C0052"), "color") <- "blue"

edgeData(testrabph, c("C0047", "R0055"), c("R0055", "C0052"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0047", "R0055"), c("R0055", "C0052"), "color") <- "blue"

edgeData(testrabph, c("C0053", "R0056"), c("R0056", "C0047"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0053", "R0056"), c("R0056", "C0047"), "color") <- "blue"

edgeData(testrabph, c("C0007", "R0057"), c("R0057", "C0027"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0007", "R0057"), c("R0057", "C0027"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0058"), c("R0058", "C0053"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0058"), c("R0058", "C0053"), "color") <- "blue"

edgeData(testrabph, c("C0034", "R0059"), c("R0059", "C0054"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0034", "R0059"), c("R0059", "C0054"), "color") <- "blue"

edgeData(testrabph, c("C0055", "R0059"), c("R0059", "C0054"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0055", "R0059"), c("R0059", "C0054"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0060"), c("R0060", "C0055"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0060"), c("R0060", "C0055"), "color") <- "blue"

edgeData(testrabph, c("C0042", "R0061"), c("R0061", "C0056"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0042", "R0061"), c("R0061", "C0056"), "color") <- "blue"

edgeData(testrabph, c("C0013", "R0061"), c("R0061", "C0056"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0013", "R0061"), c("R0061", "C0056"), "color") <- "blue"

edgeData(testrabph, c("C0045", "R0062"), c("R0062", "C0057"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0045", "R0062"), c("R0062", "C0057"), "color") <- "firebrick"

edgeData(testrabph, c("C0012", "R0063"), c("R0063", "C0058"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0063"), c("R0063", "C0058"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0063"), c("R0063", "C0058"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0063"), c("R0063", "C0058"), "color") <- "blue"

edgeData(testrabph, c("C0060", "R0064"), c("R0064", "C0059"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0060", "R0064"), c("R0064", "C0059"), "color") <- "blue"

edgeData(testrabph, c("C0027", "R0064"), c("R0064", "C0059"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0027", "R0064"), c("R0064", "C0059"), "color") <- "blue"

edgeData(testrabph, c("F6PXX", "R0065"), c("R0065", "C0060"), "lwd") <- c("5", "5")
edgeData(testrabph, c("F6PXX", "R0065"), c("R0065", "C0060"), "color") <- "blue"

edgeData(testrabph, c("C0062", "R0066"), c("R0066", "C0061"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0062", "R0066"), c("R0066", "C0061"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0066"), c("R0066", "C0061"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0066"), c("R0066", "C0061"), "color") <- "blue"

edgeData(testrabph, c("C0017", "R0067"), c("R0067", "C0062"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0067"), c("R0067", "C0062"), "color") <- "firebrick"

edgeData(testrabph, c("C0014", "R0068"), c("R0068", "C0062"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0014", "R0068"), c("R0068", "C0062"), "color") <- "firebrick"

edgeData(testrabph, c("C0012", "R0069"), c("R0069", "C0063"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0069"), c("R0069", "C0063"), "color") <- "blue"

edgeData(testrabph, c("C0012", "R0070"), c("R0070", "C0064"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0012", "R0070"), c("R0070", "C0064"), "color") <- "blue"

edgeData(testrabph, c("C0067", "R0071"), c("R0071", "C0066"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0067", "R0071"), c("R0071", "C0066"), "color") <- "blue"

edgeData(testrabph, c("C0067", "R0071"), c("R0071", "E4P"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0067", "R0071"), c("R0071", "E4P"), "color") <- "blue"

edgeData(testrabph, c("C0017", "R0072"), c("R0072", "C0067"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0017", "R0072"), c("R0072", "C0067"), "color") <- "firebrick"

edgeData(testrabph, c("F16bP", "R0073"), c("R0073", "C0067"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F16bP", "R0073"), c("R0073", "C0067"), "color") <- "firebrick"

edgeData(testrabph, c("C0069", "R0074"), c("R0074", "C0068"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0069", "R0074"), c("R0074", "C0068"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0074"), c("R0074", "C0068"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0074"), c("R0074", "C0068"), "color") <- "blue"

edgeData(testrabph, c("C0018", "R0075"), c("R0075", "C0069"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0075"), c("R0075", "C0069"), "color") <- "blue"

edgeData(testrabph, c("C0051", "R0076"), c("R0076", "C0069"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0051", "R0076"), c("R0076", "C0069"), "color") <- "blue"

edgeData(testrabph, c("C0062", "R0077"), c("R0077", "C0070"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0062", "R0077"), c("R0077", "C0070"), "color") <- "firebrick"

edgeData(testrabph, c("C0006", "R0078"), c("R0078", "C0071"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0078"), c("R0078", "C0071"), "color") <- "blue"

edgeData(testrabph, c("C0034", "R0078"), c("R0078", "C0071"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0034", "R0078"), c("R0078", "C0071"), "color") <- "blue"

edgeData(testrabph, c("C0073", "R0079"), c("R0079", "C0072"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0073", "R0079"), c("R0079", "C0072"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0080"), c("R0080", "C0073"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0080"), c("R0080", "C0073"), "color") <- "blue"

edgeData(testrabph, c("C0007", "R0080"), c("R0080", "C0073"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0007", "R0080"), c("R0080", "C0073"), "color") <- "blue"

edgeData(testrabph, c("C0027", "R0081"), c("R0081", "C0074"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0027", "R0081"), c("R0081", "C0074"), "color") <- "blue"

edgeData(testrabph, c("C0013", "R0081"), c("R0081", "C0074"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0013", "R0081"), c("R0081", "C0074"), "color") <- "blue"

edgeData(testrabph, c("C0006", "R0082"), c("R0082", "C0075"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0006", "R0082"), c("R0082", "C0075"), "color") <- "blue"

edgeData(testrabph, c("C0076", "R0082"), c("R0082", "C0075"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0076", "R0082"), c("R0082", "C0075"), "color") <- "blue"

edgeData(testrabph, c("C0009", "R0083"), c("R0083", "C0076"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0009", "R0083"), c("R0083", "C0076"), "color") <- "blue"

edgeData(testrabph, c("C0018", "R0084"), c("R0084", "C0076"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0018", "R0084"), c("R0084", "C0076"), "color") <- "blue"

edgeData(testrabph, c("C0067", "R0085"), c("R0085", "C0077"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0067", "R0085"), c("R0085", "C0077"), "color") <- "firebrick"

edgeData(testrabph, c("C0007", "R0086"), c("R0086", "C0078"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0007", "R0086"), c("R0086", "C0078"), "color") <- "blue"

edgeData(testrabph, c("C0027", "R0087"), c("R0087", "C0078"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0027", "R0087"), c("R0087", "C0078"), "color") <- "blue"

edgeData(testrabph, c("C0079", "R0087"), c("R0087", "C0078"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0079", "R0087"), c("R0087", "C0078"), "color") <- "blue"

edgeData(testrabph, c("C0007", "R0088"), c("R0088", "C0079"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0007", "R0088"), c("R0088", "C0079"), "color") <- "blue"

edgeData(testrabph, c("C0051", "R0089"), c("R0089", "C0080"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0051", "R0089"), c("R0089", "C0080"), "color") <- "blue"

edgeData(testrabph, c("C0051", "R0089"), c("R0089", "C0079"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0051", "R0089"), c("R0089", "C0079"), "color") <- "blue"

edgeData(testrabph, c("C0030", "R0090"), c("R0090", "C0007"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0030", "R0090"), c("R0090", "C0007"), "color") <- "blue"

edgeData(testrabph, c("C0030", "R0090"), c("R0090", "C0081"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0030", "R0090"), c("R0090", "C0081"), "color") <- "blue"

edgeData(testrabph, c("C0060", "R0091"), c("R0091", "C0082"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0060", "R0091"), c("R0091", "C0082"), "color") <- "blue"

edgeData(testrabph, c("C0060", "R0091"), c("R0091", "C0007"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0060", "R0091"), c("R0091", "C0007"), "color") <- "blue"

edgeData(testrabph, c("C0055", "R0092"), c("R0092", "C0083"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0055", "R0092"), c("R0092", "C0083"), "color") <- "blue"

edgeData(testrabph, c("C0055", "R0092"), c("R0092", "C0007"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0055", "R0092"), c("R0092", "C0007"), "color") <- "blue"

edgeData(testrabph, c("C0032", "R0093"), c("R0093", "C0084"), "lwd") <- c("5", "5")
edgeData(testrabph, c("C0032", "R0093"), c("R0093", "C0084"), "color") <- "blue"

edgeData(testrabph, c("C0086", "R0094"), c("R0094", "C0085"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0086", "R0094"), c("R0094", "C0085"), "color") <- "firebrick"

edgeData(testrabph, c("F16bP", "R0095"), c("R0095", "C0086"), "lwd") <- c("3", "3")
edgeData(testrabph, c("F16bP", "R0095"), c("R0095", "C0086"), "color") <- "firebrick"

edgeData(testrabph, c("C0018", "R0096"), c("R0096", "C0086"), "lwd") <- c("3", "3")
edgeData(testrabph, c("C0018", "R0096"), c("R0096", "C0086"), "color") <- "firebrick"

plot(testrabph)



png(file = "mygraphic_2.png", width = 1920 * 4, height = 1080 * 4, res = 600)
plot(testrabph)
dev.off()