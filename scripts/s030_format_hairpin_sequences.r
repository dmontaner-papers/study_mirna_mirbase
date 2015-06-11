##s030_format_hairpin_sequences.r
##2015-04-19 dmontaner@cipf.es
##Format hairpin sequences from miRBase

### NOTE:
## Hairpins are in RNA sequence (U)
## Alphabet is extended to Y K, N. This may be a problem for some aligners.
## Mean length for a hairpin 103

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
library (Biostrings); packageDescription ("Biostrings", fields = "Version") #"2.32.0"
#help (package = Biostrings)

try (source (".job.r")); try (.job)

options (width = 170)


###READ FASTA
setwd (file.path (.job$dir$raw))
seqs <- readRNAStringSet (filepath = "hairpin.fa.gz")
class (seqs)
seqs
length (seqs)

## sequence lengths
lon <- width (seqs)
names (lon) <- names (seqs)
lon[1:3]
summary (lon)


################################################################################

## ALPHABET
alphabet (seqs)

###grep ("Y", seqs, value = TRUE)
names (grep ("Y", seqs, value = TRUE))

## non standard letters
malos <- NULL
##
system.time ({
for (l in setdiff (alphabet (seqs), c ("A", "C", "G", "U"))) {
    cat  ("\n========== ", l, " =========\n")
    print (l)
    pos <- grep (l, seqs, fixed = TRUE)  ##fixed = TRUE so that + and . are not taken as part of a regular expression
    print (length (pos))
    malos <- c (malos, pos)
}
})
table (duplicated (malos)) ## some duplicates
malos <- unique (malos)
length (malos)
sort (names (seqs[malos]))

grep ("hsa", names (seqs[malos]))  ## OK NO HUMAN
grep ("mus", names (seqs[malos]))  ## OK NO mouse


## ELIMINATE sequences with non standard letters
table (malos %in% 1:length (seqs))
seqs <- seqs[-malos]
length (seqs)

################################################################################

### CONVERT TO DNA
dna <- DNAStringSet (seqs)
class (dna)
dna

table (as.character (dna) == gsub ("U", "T", as.character (seqs)))  ## OK

################################################################################

### KEEP JUST HUMANS
humanos <- grep ("hsa", names (dna))
length (humanos)

dnah <- dna[humanos]
length (dnah)


## SAVE 
setwd (.job$dir$proces)

save (dna,  file = "haiprin_dna.RData")
save (dnah, file = "haiprin_dna_human.RData")

writeXStringSet (dna,  filepath = "haiprin_dna.fa")
writeXStringSet (dnah, filepath = "haiprin_dna_human.fa")

###EXIT
warnings ()
sessionInfo ()
q ("no")
