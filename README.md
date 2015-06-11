Exploring microRNA sequences from miRBase
================================================================================

In this scripts I use R/Bioconductor to explore some files downloaded from the
[miRBase](http://www.mirbase.org/)
[FPT site](ftp://mirbase.org/pub/mirbase/CURRENT/).

I mainly focus on __human__ miRNAs.

I used miRBase version __21__.


Objectives
----------------------------------------

- Understanding miRBase _fasta_ files.
- Understanding forward and reverse strand assignment of miRNAs.
- Understanding _mature_ and _hairpin_ relationship,
  including __relative position__ of mature miRNAs respect to their _pre-miRNAs_.


Analysis
----------------------------------------

1. Download miRBase files
1. Explore human GFF file:
    - IDs, Names and Aliases
	- mature to hairpin relationship
	- length distribution
	- strand orientation
	- derive relative position (of the mature over the hairpin)
	- revise whether the __middle base__ in the hairpin separates the _5p_ _3p_ mature miRNAs
1. Format hairpin sequences:
    - explore the _alphabet_
    - exclude sequences with __ambiguity codes__. See [Ensembl Glossary](http://www.ensembl.org/Help/Glossary).
	- transform RNA into DNA
1. Compare miRBase sequences to those extracted form the human genome according to the GFF file information.


Findings
--------------------------------------------------------------------------------


1. Fasta files in miRBase use the RNA _alphabet_ (U instead of T).
1. Not many but some miRNA sequences may have __ambiguity codes__ (letters other than A, C, U, G).

1. Each hairpin may produce just 1 or 2 mature molecules, not more.
1. The same mature molecules may produced from different hairpins. 

1. Some hairpins can be __expressed form both strands__ (same chromosomic position & different strand in the GFF file).  
   They are considered two different microRNAs and have different assigned IDs in miRBase.  
   This wont be found in NGS data if they are not _directional_, and may cause double heats when aligning reads against miRBase. 

1. miRNA sequences with associated __reverse__ strand match with the reverse complementary of their genomic DNA region (according to the provided GFF file).

1. The __middle base__ in the hairpin separates the _5p_ _3p_ mature miRNAs for most of human cases.
   Just 2 human microRNAs where wrongly classified using the _middle base_.


Miscellany
--------------------------------------------------------------------------------

About [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/):

miRBase hairpin sequences need to be transformed to DNA before they can be indexed for mots aligners.

Ideally FASTX-Toolkit may help you by: 

1. Shaping _multi-line_ FASTA files to _single-line_: 

    fasta_formatter -i hairpin.fa -o formatted.fa

1. Converting RNA to DNA sequences:

    fasta_nucleotide_changer -i formatted.fa -o hairpinDNA.fa -d

But `fasta_nucleotide_changer` does not work with __ambiguity codes__.


