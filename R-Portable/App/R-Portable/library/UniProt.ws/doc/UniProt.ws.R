### R code from vignette source 'UniProt.ws.Rnw'

###################################################
### code chunk number 1: loadPkg
###################################################
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)


###################################################
### code chunk number 2: help (eval = FALSE)
###################################################
## help("UniProt.ws")


###################################################
### code chunk number 3: show
###################################################
up


###################################################
### code chunk number 4: availSpecies
###################################################
availableUniprotSpecies(pattern="musculus")


###################################################
### code chunk number 5: setTaxID
###################################################
mouseUp <- UniProt.ws(10090)
mouseUp


###################################################
### code chunk number 6: columns
###################################################
head(keytypes(up))


###################################################
### code chunk number 7: keytypes
###################################################
head(columns(up))


###################################################
### code chunk number 8: keys (eval = FALSE)
###################################################
## egs = keys(up, "ENTREZ_GENE")


###################################################
### code chunk number 9: select
###################################################
keys <- c("1","2")
columns <- c("PDB","UNIGENE","SEQUENCE")
kt <- "ENTREZ_GENE"
res <- select(up, keys, columns, kt)
head(res)


