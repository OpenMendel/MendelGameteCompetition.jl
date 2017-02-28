# MendelGameteCompetition

This [Julia](http://julialang.org/) package  implements a gamete competition analysis, which is a generalization of the TDT analysis.  MendelGameteCompetition is one component of the umbrella [OpenMendel](https://openmendel.github.io) project.

[![](https://img.shields.io/badge/docs-current-blue.svg)](https://OpenMendel.github.io/MendelGameteCompetition.jl)

## Installation

*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelGameteCompetition:

    Pkg.clone("https://github.com/OpenMendel/MendelGameteCompetition.jl.git")

This package supports Julia v0.4 and v0.5.

## Data Files

To run this analysis package you will need to prepare a Control file and have your data files available. The Control file holds the names of your data files and any optional parameters for the analysis. Details on the general format and contents of the Control and data files can be found on the MendelBase [documentation page](https://openmendel.github.io/MendelBase.jl). Descriptions of the specific options available within the MendelGameteCompetition analysis package are in its [documentation page](https://openmendel.github.io/MendelGameteCompetition.jl).

There are example data files in the "docs" subfolder of each Mendel package, for example, ~/.julia/v0.5/MendelGameteCompetition/docs.

## Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelGameteCompetition

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> GameteCompetition("Control_file.txt")

*Note: The package is called* MendelGameteCompetition *but the analysis function is called simply* GameteCompetition.

## Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

*Sinsheimer JS, Blangero J, Lange K (2000). Gamete competition models. American Journal of Human Genetics 66:1168-1172.*

*Sinsheimer JS, McKenzie CA, Keavney B, Lange K (2001). SNPs and snails and puppy dogs' tails: Analysis of SNP data using the gamete competition model. Annals of Human Genetics 65:483-490.*

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
