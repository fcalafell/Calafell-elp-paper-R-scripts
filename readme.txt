This folder contains the R scripts used in the paper "An assessment of a massively parallel sequencing approach for the 
identification of individuals from mass graves of the Spanish Civil War (1936-1939)”, 
by Calafell et al., submitted to Electrophoresis, as well as example input files.

Basic R skills are required to run the scripts and they must be edited to use different parameters. The parameters that 
can be modified are commented on within the scripts.

The scripts are:
LR_true_father_offspring.R. Simulates the LR distribution between a parent and offspring for autosomal loci. It reads an 
allele frequency matrix in the csv format (Excel files can easily be exported as csv, but they use semicolons rather than 
comas. A text editor such as textpad should be used then to replace all semicolons with comas). Alleles are in rows and 
STRs or SNPs in columns, and each cell contains a relative allele frequency. Rows and columns should not be labeled. An 
example input file, alleles.csv, is included. For reference, an Excel file (autosomal loci reference.xlsx) containing the 
same information as alleles.csv is included, but with row and column labels; this cannot be used to run the script. Within 
the script, the number of loci and the maximum number of alleles per locus must be specified. Parameters defining the number
of iterations, whether heterozygotes only are considered, and the fraction of allele drop-out can be edited. LR summary 
statistics are printed to the screen, and the whole LR distribution is written to an output file.

LR_false_father_offspring.R. Simulates unrelated individuals and compares for autosomal loci them as if they were candidate 
parent and offspring. In each simulated case, paternity may be excluded or not. The frequency of exclusions is counted and 
printed to the screen; for non-exclusions, an LR is computed; LR summary statistics are printed to the screen, and the whole 
LR distribution is written to an output file. The input files and modifiable parameters are as in LR_true_father_offspring.R 
(see above).

LR_uncle_nephewniece.R. Simulates the LR distribution between second degree relations (uncle vs. nephew/niece; grandfather vs.
grandchild) for autosomal loci. It can produce either the LR distribution between second degree relations, or between unrelated
individuals compared as if they were second degree relations. Note that second degree relations cannot be excluded solely from
autosomal loci. To switch between these two modes, set the variable truepaternity to 1 or 0, respectively. Otherwise, the input
files and modifiable parameters are as in LR_true_father_offspring.R (see above). 

yhapstrue.R. Simulates the LR distribution of pairs of patrilineal relatives by using Y-chromosome haplotypes. Mutation is not 
considered. The example input file is yhaps.csv, in which haplotypes are in rows and STR alleles are in columns, with the last
column designating the absolute allele frequency of the haplotype. Rows and columns must not be labeled. For easier reference,
a version of yhaps.csv with labeled columns, yhapref.xlsx, is included.  Within the script, the numbers of loci, haplotypes and
individuals present in the input file must be specified. Parameters defining the number of iterations, and the fraction of allele
drop-out can be edited. LR summary statistics are printed to the screen, and the whole LR distribution is written to an output file.

yhapsfalse.R. Simulates unrelated individuals and compares them for Y-chromosome haplotypes as if they were patrilinearly related.
Mutation is not considered. In each simulated case, relatedness may be excluded or not. The frequency of exclusions is counted and
printed to the screen; for non-exclusions, an LR is computed; LR summary statistics are printed to the screen, and the whole LR
distribution is written to an output file. The input files and modifiable parameters are as in yhapstrue.R (see above)

LR_fatherdaughterX.R Simulates the LR distribution between father and daughter for X-linked loci (STRs or STR haplotypes). It can
produce either the LR distribution between father and daughter, or between unrelated individuals compared as if they were father and
daughter; in the latter case, it also provides the proportion of cases showing exclusion. To switch between these two modes, set the
variable truepaternity to 1 or 0, respectively. The example input file is xhaps.csv, in which haplotypes or alleles are in rows loci
are in columns. Rows and columns must not be labeled. For easier reference, a version of xhaps.csv with labeled columns,
xhapsreference.xlsx, is included.  Within the script, the numbers of loci and maximum number of alleles or haplotypes per locus must
be specified. Parameters defining the number of iterations, and the fraction of allele drop-out can be edited. LR summary statistics,
chance of exclusion, and mean number of exclusionary loci are printed to the screen, and the whole LR distribution is written to an 
output file.

LR_gf_maternalgrandsonX.R Simulates the LR distribution between a grandfather and a maternal grandson (a man’s daughter’s son) for 
X-linked loci (STRs or STR haplotypes). It can produce either the LR distribution between a grandfather and a maternal grandson, or
between unrelated individuals compared as if they were a grandfather and a maternal grandson. To switch between these two modes, set
the variable truepaternity to 1 or 0, respectively. The input files and modifiable parameters are as in LR_fatherdaughterX.R (see above)

LR_uncle_maternalnephewX.R Simulates the LR distribution between an uncle and a maternal nephew (a man’s sister’s son) for X-linked loci
(STRs or STR haplotypes). It can produce either the LR distribution between an uncle and a maternal nephew, or between unrelated individuals
compared as if they were an uncle and a maternal nephew. To switch between these two modes, set the variable truepaternity to 1 or 0,
respectively. The input files and modifiable parameters are as in LR_fatherdaughterX.R (see above)

LR_uncle_maternalnieceX.R Simulates the LR distribution between an uncle and a maternal niece (a man’s sister’s daughter) for X-linked loci
(STRs or STR haplotypes). It can produce either the LR distribution between an uncle and a maternal niece, or between unrelated individuals
compared as if they were an uncle and a maternal niece. To switch between these two modes, set the variable truepaternity to 1 or 0,
respectively. The input files and modifiable parameters are as in LR_fatherdaughterX.R (see above)

LR_uncle_paternalnieceX.R Simulates the LR distribution between an uncle and a paternal niece (a man’s brother’s daughter) for X-linked loci
(STRs or STR haplotypes). It can produce either the LR distribution between an uncle and a paternal niece, or between unrelated individuals
compared as if they were an uncle and a paternal niece. To switch between these two modes, set the variable truepaternity to 1 or 0,
respectively. The input files and modifiable parameters are as in LR_fatherdaughterX.R (see above)

homdistr.R generates the distribution of the number of expected homozygote genotypes given a set of loci (and their allele frequencies), a
fraction of drop-out alleles, and Hardy-Weinberg’s law. It reads an allele frequency matrix in the csv format (Excel files can easily be
exported as csv, but they use semicolons rather than comas. A text editor such as textpad should be used then to replace all semicolons with
comas). Alleles are in rows and STRs or SNPs in columns, and each cell contains a relative allele frequency. Rows and columns should not be
labeled. An example input file, alleles.csv, is included. For reference, an Excel file (autosomal loci reference.xlsx) containing the same
information as alleles.csv is included, but with row and column labels; this cannot be used to run the script. Within the script, the number
of loci and the maximum number of alleles per locus must be specified. Parameters defining the number of iterations, whether heterozygotes
only are considered, and the fraction of allele drop-out can be edited. It prints to the screen the mean, median, and 95% CI of expected
number of homozygotes, and it writes to a file the whole distribution of the number of expected homozygotes. The output file name can be
edited from within the script.   

  
