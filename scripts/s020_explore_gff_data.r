##s020_explore_gff_data.r
##2015-02-28 dmontaner@cipf.es
##Study of miRBase sequences

##gff-version 3
##date 2014-6-22
#
# Chromosomal coordinates of Homo sapiens microRNAs
# microRNAs:               miRBase v21
# genome-build-id:         GRCh38
# genome-build-accession:  NCBI_Assembly:GCA_000001405.15
#
# Hairpin precursor sequences have type "miRNA_primary_transcript". 
# Note, these sequences do not represent the full primary transcript, 
# rather a predicted stem-loop portion that includes the precursor 
# miRNA. Mature sequences have type "miRNA".
#

## GTF format: 
## http://www.ensembl.org/info/website/upload/gff.html

## 1. seqname - name of the chromosome or scaffold; chromosome names can be given
##    with or without the 'chr' prefix.
##    Important note: the seqname must be one used within Ensembl, i.e. a standard
##    chromosome name or an Ensembl identifier such as a scaffold ID, without any
##    additional content such as species or assembly. See the example GFF output below.
## 2. source - name of the program that generated this feature, or the data source
##    (database or project name)
## 3. feature - feature type name, e.g. Gene, Variation, Similarity
## 4. start - Start position of the feature, with sequence numbering starting at 1.
## 5. end - End position of the feature, with sequence numbering starting at 1.
## 6. score - A floating point value.
## 7. strand - defined as + (forward) or - (reverse).
## 8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the
##    feature is the first base of a codon, '1' that the second base is the first
##    base of a codon, and so on.
## 9. attribute - A semicolon-separated list of tag-value pairs, providing
##    additional information about each feature.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.0 (2015-04-16)"

try (source (".job.r")); try (.job)

options (width = 170)


###DATA
datos <- read.table (file.path (.job$dir$raw, "hsa.gff3"), as.is = TRUE)
colnames (datos) <- c ("chr", "source", "feature", "sta", "end", "score", "strand", "frame", "attribute")
sapply (datos, class)
dim (datos)
datos[1:3,]

summary (datos)

################################################################################

## Some descriptive
sort (table (datos[,"chr"]))
table (datos[,"strand"])
table (datos[,"chr"], datos[,"strand"])

table  (datos[,"feature"])
unique (datos[,"feature"])

mat <- datos[datos$feature== "miRNA",                    c ("chr", "sta", "end", "strand", "attribute")]
hpi <- datos[datos$feature== "miRNA_primary_transcript", c ("chr", "sta", "end", "strand", "attribute")]

dim (mat)
dim (hpi)

mat[1:3,]
hpi[1:3,]

################################################################################

## Explore duplicates

table (duplicated (mat)) ## OK no dups
table (duplicated (hpi)) ## OK no dups

table (duplicated (mat[,c("chr", "sta", "end")]))           ## Duplicated position
table (duplicated (mat[,c("chr", "sta", "end", "strand")])) ## But not the strand

table (duplicated (hpi[,c("chr", "sta", "end")]))    ## Duplicated position
table (duplicated (hpi[,c("chr", "sta", "strand")])) ## But not the strand

## duplicated in mature
dup.mat <- mat[duplicated (mat[,c("chr", "sta", "end")]),]
touse <- (mat$chr %in% dup.mat$chr) & (mat$sta %in% dup.mat$sta) & (mat$end %in% dup.mat$end)

dup.mat <- mat[touse,]
orden <- order (dup.mat$chr, dup.mat$sta, dup.mat$strand)
dup.mat <- dup.mat[orden,]
dup.mat


## duplicated in hairpin
dup.hpi <- hpi[duplicated (hpi[,c("chr", "sta", "end")]),]
touse <- (hpi$chr %in% dup.hpi$chr) & (hpi$sta %in% dup.hpi$sta) & (hpi$end %in% dup.hpi$end)

dup.hpi <- hpi[touse,]
orden <- order (dup.hpi$chr, dup.hpi$sta, dup.hpi$strand)
dup.hpi <- dup.hpi[orden,]
dup.hpi

################################################################################

## FORMAT IDS
spl <- strsplit (mat$attribute, split = ";")
lon <- sapply (spl, length)
table (lon)
##
spl[1:3]
##
mat[,"ID"]           <- sapply (spl, function (x) sub ("ID=",           "", x[1]))
mat[,"Alias"]        <- sapply (spl, function (x) sub ("Alias=",        "", x[2]))
mat[,"Name"]         <- sapply (spl, function (x) sub ("Name=",         "", x[3]))
mat[,"Derives_from"] <- sapply (spl, function (x) sub ("Derives_from=", "", x[4]))

spl <- strsplit (hpi$attribute, split = ";")
lon <- sapply (spl, length)
table (lon)
##
spl[1:3]
##
hpi[,"ID"]           <- sapply (spl, function (x) sub ("ID=",           "", x[1]))
hpi[,"Alias"]        <- sapply (spl, function (x) sub ("Alias=",        "", x[2]))
hpi[,"Name"]         <- sapply (spl, function (x) sub ("Name=",         "", x[3]))



## REVISE DUPLICATES

mat[1:3,]
hpi[1:3,]

table (duplicated (mat[,"ID"]))    ## OK UNIQUE IDS
table (duplicated (mat[,"Alias"]))
table (duplicated (mat[,"Name"]))

table (duplicated (hpi[,"ID"]))    ## OK UNIQUE IDS
table (duplicated (hpi[,"Name"]))  ## OK UNIQUE IDS

rownames (hpi) <- hpi[,"ID"]
hpi[1:3,]


## REORDER
orden <- order (mat[,"Derives_from"], mat[,"ID"])
mat <- mat[orden,]

orden <- order (hpi[,"ID"])
hpi <- hpi[orden,]

mat[1:3,]
hpi[1:3,]

################################################################################

### Mature from hairpin
table (duplicated (mat[,"Name"]))
table (duplicated (mat[,"Derives_from"]))
table (duplicated (mat[,c("Name", "Derives_from")]))  ## OK unique

table (table (mat[,"Derives_from"]))  ## OK a hairpin produces at most 2 matures

na <- names (sort (table (mat[,"Derives_from"])))[1:10] ## some of the hpi with a single mat
na
hpi[na,]


touse <- mat[,"Name"] %in% mat[duplicated (mat[,"Name"]), "Name"]
table (touse)

dup.mat <- mat[touse,]
orden <- order (dup.mat$Name, dup.mat$Derives_from)
dup.mat <- dup.mat[orden,]
dup.mat[1:10,]

table (table (dup.mat[,"Name"]))
sort (table (dup.mat[,"Name"]))

################################################################################

### LENGTHS
mat[,"lon"] <- mat[,"end"] - mat[,"sta"] + 1
hpi[,"lon"] <- hpi[,"end"] - hpi[,"sta"] + 1

summary (mat[,"lon"])  ## OK ALL POSITIVE
summary (hpi[,"lon"])  ## OK ALL POSITIVE




################################################################################
### Find RELATIVE position
################################################################################

## Find relative positions of the mature sequences over the hairpins

## case EXAMPLE

mat[1:10,]

h.p <- "MI0000060"  ## example of hairpin from POSITIVE strand
h.n <- "MI0000061"  ## example of hairpin from POSITIVE strand

hpi[h.p,]
mat[mat$Derives_from == h.p,]
mat[mat$Derives_from == h.p, c ("sta", "end")] - hpi[h.p, "sta"] + 1

hpi[h.n,]
mat[mat$Derives_from == h.n,]
hpi[h.n, "end"] - mat[mat$Derives_from == h.n, c ("end", "sta")] + 1


### Compute all relative positions
mat[1:3,]
table (mat$Derives_from %in% hpi$ID)  ## OK all in

mat[mat$strand == "+", "ref"] <- hpi[mat$Derives_from, "sta"][mat$strand == "+"]
mat[mat$strand == "-", "ref"] <- hpi[mat$Derives_from, "end"][mat$strand == "-"]

mat[,c ("Rsta", "Rend")] <- NA
mat[mat$strand == "+", c ("Rsta", "Rend")] <-   mat[mat$strand == "+", c ("sta", "end")] - mat[mat$strand == "+", "ref"] + 1
mat[mat$strand == "-", c ("Rsta", "Rend")] <- - mat[mat$strand == "-", c ("end", "sta")] + mat[mat$strand == "-", "ref"] + 1

mat[1:10, c ("strand", "ID", "Name", "Rsta", "Rend", "Derives_from")]

tail (mat[,c ("strand", "ID", "Name", "Rsta", "Rend", "Derives_from")])

hpi[mat[1:10, "Derives_from"],]



################################################################################
### Find the tail
################################################################################

### Fin the tail (or side: 5'; 3') of the hairpin
##  in which the mature miRNA is allocated

mat[,"hpi.center"] <- hpi[mat$Derives_from, "lon"] / 2
mat[1:3,]

mat[,"p"] <- NA
##
touse <- (mat[,"Rsta"] < mat[,"hpi.center"]) & (mat[,"Rend"] <= mat[,"hpi.center"])
table (touse)
mat[touse, "p"] <- 5
##
touse <- (mat[,"hpi.center"] <= mat[,"Rsta"]) & (mat[,"hpi.center"] < mat[,"Rend"])
table (touse)
mat[touse, "p"] <- 3
##
table (mat[,"p"], exclude = NULL)   ## almost all of them get a classification
##
mat[is.na (mat$p),]


## Use the middle position of the mature miRNA for the tail allocation
mat[,"mat.center"] <- (mat[,"Rend"] + mat[,"Rsta"]) / 2
##
table (mat$mat.center == mat$hpi.center) ## OK not equal centers
##
mat[,"tail"] <- NA
mat[mat$mat.center < mat$hpi.center, "tail"] <- 5
mat[mat$mat.center > mat$hpi.center, "tail"] <- 3
##
table (mat[,"tail"], exclude = NULL)

table (tail = mat[,"tail"], p = mat[,"p"], exclude = NULL)
mat[is.na (mat$p),]


## review by name
table (mat$tail, exclude = NULL)
mat.5p <- mat[mat$tail == 5,]
mat.3p <- mat[mat$tail == 3,]
dim (mat.5p)
dim (mat.3p)

grep ("3p", mat.5p[,"Name"], value = TRUE)
grep ("5p", mat.3p[,"Name"], value = TRUE)

malos <- mat[,"Name"] %in% grep ("5p", mat.3p[,"Name"], value = TRUE)
table (malos)
mat[malos,]    ## WRONG TAIL

mat[mat$Derives_from == "MI0006327",]
mat[mat$Derives_from == "MI0022631",]


###SAVE
setwd (.job$dir$proces)
save (mat, hpi, file = "info_from_gff.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
