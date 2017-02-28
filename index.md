### Overview
Mendel Gamete Competition is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. The gamete competition model is an application of the Bradley-Terry model and can be considered a parametric form of the TDT for use with pedigree data because, besides getting p-values, we also get a measure of the strength of the allelic associations. The Bradley-Terry model was originally applied to problems such as ranking teams in a sports league based on the intra-league win/loss records. In genetics, alleles assume the role of teams, and transmission parameters (the *τs*) assume the role of the winning propensities [2](#2), [3](#3). As implemented in this version of Open Mendel, the gamete competition is an affected only association analysis and we assume the allele frequencies for the marker are known without error.  Let allele *i* be assigned a segregation parameter *τ<sub>i</sub>* , then the probability that a heterozygous parent with genotype *i*=j transmits allele *i* is the ratio *τ<sub>i</sub>/(τ<sub>i</sub>+τ<sub>j</sub>)*.  Because this ratio is invariant when *i* and *j* are multiplied by the same constant *c*, we need to impose the constraint that the most frequent allele *k* has segregation parameter *τ<sub>k</sub>* = 1. These propensities replace the normally used Mendelian segregation parameters for heterozygous parents' transmissions in the Elston-Stewart-Ott representation of the likelihood of the pedigrees.  In fact, under the null of no association between the marker and the trait in the gamete competition, Mendelian segregation ratios hold for heterozygous parents' transmissions so that *τ<sub>i</sub>*=1 is true for all alleles *i* and the likelihood reverts to the standard one. Note that the transmissions for homozygous parents always conform to standard Mendelian segregation ratios both under the null or alternative thus, like the TDT, only heterozygous parents are informative.  To test whether Mendelian segregation can be rejected, we estimate these *τs* by maximum likelihood and conduct a likelihood ratio test. P-values are calculated assuming that the likelihood ratio test statistic is asymptotically chi square distributed. The degrees of freedom are equal to the number of alleles minus 1.

### Appropriate Problems and Data Sets
The Gamete Competition model applies to pedigrees, including those with missing marker data. With too many marker alleles computational efficiency suffers and large sample statistical assumptions become suspect. We recommend consolidating alleles until at most eight alleles remain and each has a frequency of 0.05 or greater. If the fraction of missing data is large, ethnic stratification may come into play. One remedy is to limit analysis to a single ethnic group; another is to use ethnic-specific allele frequencies. If you opt for the latter strategy, then you cannot simultaneously estimate allele frequencies and transmission parameters.

### Installation
*Note: Three OpenMendel packages - [SnpArrays](https://github.com/OpenMendel/SnpArrays.jl), [Search](https://github.com/OpenMendel/Search.jl), and [MendelBase](https://github.com/OpenMendel/MendelBase.jl) must be installed before any Mendel analysis packages will run.*

Within Julia, use the package manager to install MendelGameteCompetition:

    Pkg.clone("https://github.com/OpenMendel/MendelGameteCompetition.jl.git")

This package supports Julia v0.4.

### Input Files
The MendelGameteCompetition analysis package uses the following input files. Example input files can be found in the [docs]( https://github.com/OpenMendel/MendelGameteCompetition.jl/tree/master/docs) subfolder of the MendelGameteCompetition project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File]( https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File]( https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File]( https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

### Control file<a id="control-file"></a>
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Gamete Competition:


	#
	# Input and Output files.
	#
	locus_file = gamete competition LocusFrame.txt
	pedigree_file = gamete competition PedigreeFrame.txt
	output_file = gamete competition Output.txt
	#
	# Analysis parameters for Gamete Competition option.
	#
	disease_status = ACE
	affected_designator = 1
	standard_errors = true

In the example above, there are six keywords. The three keywords specify the input and output files: *gamete competition LocusFrame.txt*, *gamete competition PedigreeFrame.txt*, and *gamete competition Output.txt*. The last three keywords specify the analysis parameters: *disease_status*, *affected_designator*, and *standard_errors*. The text after the '=' are the keyword values.

### Keywords<a id="keywords-table"></a>
This is a list of OpenMendel keywords specific to Gamete Competition. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)


 Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------      |  ----------------       |  ----------------      |  ----------------
GameteCompetition_output_file  |GameteCompetition_Output_File.txt | User defined output file name |   Creates a lod score table output file 
repetitions          |                   |                         |
xlinked_analysis  |  FALSE  |  TRUE, FALSE  |  Whether or not markers are on the X chromosome


### Data Files
Gamete Competition requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file]( https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File]( https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Gamete Competition [docs]( https://github.com/OpenMendel/MendelGameteCompetition.jl/tree/master/docs) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelGameteCompetition

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> GameteCompetition("Control_file.txt")

*Note: The package is called* MendelGameteCompetition *but the analysis function is called simply* GameteCompetition.

### Interpreting the results
There are two forms of output.  A table is output to the screen that corresponds to a data frame that can be used in other analyses as desired.  For the SNP data provided Open Mendel, the results are:

Row | Marker  | LowAllele | Low τ   | HighAllele | High τ | Pvalue  |
----|  -----  |-------    |  --------|--------   |  --------|----     |
 1  | SNP1  | 1      | 1.0       | 2       | 5.41918  | 2.80892e-5 |
 2  | SNP2  | 1       | 1.0      |  2      | 5.0539   | 2.57148e-5 |
 3         | SNP3|  1          | 1.0            |      2        |     5.23133 | 3.22822e-6 | 
 4         | SNP4|  1          | 1.0            |      2        |     4.52948 | 3.96425e-5 |
 5         | SNP5|  1          | 1.0            |      2        |     4.59438 | 3.25329e-5 |
 6         | SNP6|  1          | 0.654208 |      2        |     1.0          | 0.00527033|
 7         | SNP7|  2          | 1.0            |      1        |     8.04062 | 2.28436e-6 |
 8         | ID      | 1           | 1.0            |     2         |     6.67263 | 4.60103e-6 |
 9         | SNP9| 1           | 1.0            |     2         |     7.45069 | 7.69303e-6 |
 10       | CT     | 1           | 1.0            |    2          |     8.51002 | 7.94566e-7 | 

For each marker, the allele with the smallest transmission, its corresponding *τ*, the allele with the largest transmission, and its corresponding *τ* are provided along with the p-value for the test of association of the marker with the trait. In the example provided, the CT Marker is the most associated with the trait because it has the smallest p-value.  The most frequent allele is the 1 allele so it is assigned *τ1* = 1.  The 2 allele is ~8.51 times more likely to be transmitted from a 1/2 parent than the 1 allele.  Details of the analysis are provided in the output text file. In this text file, the iterations of the numeric loglikelihood maximization, the maximum likelihood estimates at the maximum log likelihood, their standard errors and their correlations are provided for each marker (see the example output, [gamete competition Output.txt](https://github.com/OpenMendel/MendelGameteCompetition.jl/blob/master/docs/gamete%20competition%20Output.txt))

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

<a id="1"></a> *Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<a id="2"></a> *Sinsheimer JS, Blangero J, Lange K (2000). Gamete competition models. American Journal of Human Genetics 66:1168-1172.*

<a id="3"></a> *Sinsheimer JS, McKenzie CA, Keavney B, Lange K (2001). SNPs and snails and puppy dogs' tails: Analysis of SNP data using the gamete competition model. Annals of Human Genetics 65:483-490.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.