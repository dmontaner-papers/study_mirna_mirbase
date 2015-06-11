##s050_hairping_vs_genome.r
##2015-03-01 dmontaner@cipf.es
##We compare miRBase sequences with those taken from the reference genome according to the human GFF file

##gff-version 3
##date 2014-6-22
#
# Chromosomal coordinates of Homo sapiens microRNAs
# microRNAs:               miRBase v21
# genome-build-id:         GRCh38
# genome-build-accession:  NCBI_Assembly:GCA_000001405.15

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.1 (2014-07-10)"
library (Biostrings); packageDescription ("Biostrings", fields = "Version") #
#help (package = Biostrings)

try (source (".job.r")); try (.job)

options (width = 170)

###DATA
setwd (.job$dir$proces)
load ("info_from_gff.RData")
load ("haiprin_dna_human.RData")
ls ()

mat[1:3,]
hpi[1:3,]
dnah

dim (mat)
dim (hpi)
length (dnah)  ## OK as hpi

##rename hairping sequences
names (dnah)[1:3]
names (dnah) <- sapply (strsplit (names (dnah), " "), "[", 2)
dnah

table (names (dnah) == rownames (hpi))  ## OK ALL IN ORDER

################################################################################


### READ GENOMIC DNA
system.time (ref  <- readDNAStringSet (filepath = file.path (.job$dir$annot, "GCA_000001405.15_GRCh38_genomic.fna.gz"), format="fasta"))
## system.time (ref0 <- readDNAStringSet (filepath = file.path (.job$dir$annot, "GCA_000001405.15_GRCh38_genomic.fna"),    format="fasta"))
## identical (ref, ref0)

ref
length (ref)
names (ref)

lon <- width (ref)
names (lon) <- names (ref)
sort (lon)
sum (as.numeric (lon)) #3,209,286,105


## rename (just) chromosomes 
names (ref)[1:3]

cromosomas <- c (grep ("\\d, GRCh38 reference primary assembly", names (ref), value = TRUE),
                 grep (  "X, GRCh38 reference primary assembly", names (ref), value = TRUE),
                 grep (  "Y, GRCh38 reference primary assembly", names (ref), value = TRUE))
cr <- sub (", GRCh38 reference primary assembly", "", cromosomas)
li <- strsplit (cr, split = " Homo sapiens chromosome ")
cr <- sapply (li, "[", 2)
cr
newnames <- paste0 ("chr", cr)
cbind (cromosomas, newnames)

##rename 
for (i in 1:length (cromosomas)) {
    print (newnames[i])
    touse <- names (ref) == cromosomas[i]
    print (sum (touse))
    names (ref)[touse] <- newnames[i]
}

names (ref)[1:3]

################################################################################

### Extract hairpin sequences from the genome
system.time ({
    seqs <- DNAStringSet (NULL)
    for (i in 1:nrow (hpi)) {
        sq <- subseq (ref[ hpi[i, "chr"]] , start = hpi[i, "sta"], end = hpi[i, "end"])
        ##sq <- subseq (ref[[hpi[i, "chr"]]], start = hpi[i, "sta"], end = hpi[i, "end"])  ## not good to concatenate with c
        names (sq) <- hpi[i, "ID"]
        seqs <- c (seqs, sq)
    }
})
seqs
table (rownames (hpi) == names (seqs))  ## OK ALL IN ORDER

## REVERSE COMPLEMENTARY
torev <- hpi[,"strand"] == "-"
table (torev)
seqs[torev] <- reverseComplement (seqs[torev])

################################################################################

### COMPARE with miRBase reads
table (names (dnah) == names (seqs)) ## OK same order

table (dnah == seqs) ## OK ALL THE SAME

dnah[1:2]
seqs[1:2]


###EXIT
warnings ()
sessionInfo ()
q ("no")
