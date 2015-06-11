##s010_download_mirbase_data.r
##2015-06-08 dmontaner@cipf.es
##Study of miRBase sequences

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.0 (2015-04-16)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

options (width = 170)


###DATA
setwd (file.path (.job$dir$rawdat))
dir ()

system.time (download.file (url = "ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz",    destfile = "hairpin.fa.gz"))
system.time (download.file (url = "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",     destfile = "mature.fa.gz"))
system.time (download.file (url = "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3", destfile = "hsa.gff3"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
