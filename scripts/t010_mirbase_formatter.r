##t010_mirbase_formatter.r
##2015-06-12 dmontaner@cipf.es
##Download miRBase hairpin information and format to use it with aligners

## To Do: make a function with this pipeline

## igual seria mejor intentar construir un GFF para las secuencias...

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.0 (2015-04-16)"
library (Biostrings); packageDescription ("Biostrings", fields = "Version") #"2.36.1"
#help (package = Biostrings)

try (source (".job.r")); try (.job)

options (width = 170)


### GET DATA FROM MIRBASE
seq <- readRNAStringSet      ("ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz")
str <- readLines (gzcon (url ("ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz")))


################################################################################


### EXTRACT INFORMATION from sequences
datos <- data.frame (seq.names = names (seq), stringsAsFactors = FALSE)

nombres <- names (seq)
nombres <- sub (" +stem-loop", "", nombres)
nombres <- sub (" +stem +loop", "", nombres)
##
nombres <- strsplit (nombres, split = " +")
##
datos[,"id"]  <- sapply (nombres, "[", 1)
datos[,"acc"] <- sapply (nombres, "[", 2)
datos[,"lon"] <- width (seq)
datos[,"spe"] <- substr (datos[, "id"], 1, 3)

##sort (table (datos[,"spe"]))

if (any (duplicated (datos[,"id"]))) stop ("duplicated IDs in sequences")


################################################################################


### EXTRACT INFORMATION FORM LINES
str <- grep ("^>", str, value = TRUE)
str <- sub ("^>", "", str)

datos1 <- data.frame (str.names = str, stringsAsFactors = FALSE)

str <- strsplit (datos1[,"str.names"], split = " +")
 datos1[,"lon"] <- sapply (str, length)
datos1[,"id"]  <- sapply (str, "[", 1)
 datos1[,"parentesis"] <- sapply (str, "[", 2)
datos1[,"p5"] <- sapply (str, "[", 3)
datos1[,"p3"] <- sapply (str, "[", 4)

datos1[,"p5"] <- sub ("\\[", "", datos1[,"p5"])
datos1[,"p5"] <- sub ("]",   "", datos1[,"p5"])
##
datos1[,"p3"] <- sub ("\\[", "", datos1[,"p3"])
datos1[,"p3"] <- sub ("]",   "", datos1[,"p3"])

p5 <- strsplit (datos1[,"p5"], split = ":")
datos1[,"p5.id"] <- sapply (p5, "[", 1)
##
p5 <- strsplit (sapply (p5, "[", 2), split = "-")
datos1[,"p5.ini"] <- as.integer (sapply (p5, "[", 1))
datos1[,"p5.end"] <- as.integer (sapply (p5, "[", 2))
####
p3 <- strsplit (datos1[,"p3"], split = ":")
datos1[,"p3.id"] <- sapply (p3, "[", 1)
##
p3 <- strsplit (sapply (p3, "[", 2), split = "-")
datos1[,"p3.ini"] <- as.integer (sapply (p3, "[", 1))
datos1[,"p3.end"] <- as.integer (sapply (p3, "[", 2))

table (datos1$lon)
datos1[datos1$lon > 4,] ## MORE THAN 2 MATURE SEQUENCES !!!!!

################################################################################

### MERGE DATASETS
if (any (datos$id != datos1$id)) stop = "WRONG ORDER"

datos <- cbind (datos, datos1[,c ("p5.id", "p3.id", "p5.ini", "p5.end",  "p3.ini", "p3.end")])

##############################

head (datos)
sapply (datos, class)

table (is.na (datos$p5.end))
table (is.na (datos$p5.ini))

table (is.na (datos$p3.end))
table (is.na (datos$p3.ini))
table (is.na (datos$p3.id))

table (duplicated (datos[,"p5.id"]))
table (duplicated (datos[,"p3.id"]))

################################################################################

### OVERLAPPING p5 p3
## It seems that some p3 or p3 mature sequences may have different ISOFORMS
## Those are the overlapping ones
malos <- datos[,"p5.end"] > datos[,"p3.ini"]
table (malos, exclude = NULL)
malos <- which (malos)
malos

## There are few of this weird miRNAs
## and non of them in human
## We KEEP JUST THE FIRST of the mature sequences 
datos[malos, c ("p3.id",  "p3.ini", "p3.end")] <- NA

################################################################################

### DEFINE CENTRAL POSITION
datos[,"center"] <- round ((datos[,"p5.end"] + datos[,"p3.ini"]) / 2)
datos[is.na (datos$center), "center"] <- datos[is.na (datos$center), "lon"]  ## so that we will always assign the p5 mature seq

################################################################################

### EXPLORE ALPHABET
malos <- NULL
for (l in setdiff (alphabet (seq), c ("A", "C", "G", "U"))) {
    pos <- grep (l, seq, fixed = TRUE)  ##fixed = TRUE so that + and . are not taken as part of a regular expression
    malos <- c (malos, pos)
}
malos <- unique (malos)
length (malos)
datos[,"ambiguous"] <- FALSE
datos[malos, "ambiguous"] <- TRUE

### RNA to DNA
seq <- DNAStringSet (seq)


################################################################################

### RENAME SEQUENCES
table (names (seq) == datos$seq.names)
table (duplicated (datos$id))

names (seq) <- rownames (datos) <- datos$id



################################################################################
### FINAL FILTERING
################################################################################

### FILTERING out hairpins with ambiguous nucleotides
seq <- seq[!datos$ambiguous,]
datos <- datos[!datos$ambiguous,]

### FILTERING SPECIES
##slist <- c ("hsa", "mmu")
slist <- c ("hsa")
touse <- datos$spe %in% slist
table (touse)

datos <- datos[touse,]
seq <- seq[touse]

table (names (seq) == rownames (datos))


### SAVE
setwd (.job$dir$proces)
writeXStringSet (seq, filepath = "mirbase_hairpin_clean.fa", compress = FALSE, format = "fasta", width = 10000)

save (datos, file = "mirbase_hairpin_clean_info.RData")

###EXIT
warnings ()
sessionInfo ()
q ("no")
