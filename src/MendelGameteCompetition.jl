"""
This module orchestrates a gamete competition analysis.
"""
module MendelGameteCompetition
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures           # Now in MendelBase.
# using ModelConstruction        # Now in MendelBase.
# using ElstonStewartPreparation # Now in MendelBase.
# using ElstonStewartEvaluation  # Now in MendelBase.
using Search
using SearchSetup
#
# Required external modules.
#
using DataFrames    # From package DataFrames.
using Distributions # From package Distributions.

export GameteCompetition

"""
This is the wrapper function for the Gamete Competition analysis option.
"""
function GameteCompetition(control_file = ""; args...)

  const GAMETE_COMPETITION_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Gamete Competition analysis option")
  println("        version ", GAMETE_COMPETITION_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["gamete_competition_table"] = "gamete competition Table Output.txt"
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "gametecompetition")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "GameteCompetition"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = false
  skipped_loci = gamete_competition_option(pedigree, person, nuclear_family,
    locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function GameteCompetition

"""
This function tests the gamete competition model at each locus.
It returns the number of skipped loci. Relevant keywords are
disease_status and affected_designator.
"""
function gamete_competition_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame, 
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{ASCIIString, Any})

  io = keyword["output_unit"]
  skipped_loci = 0
  #
  # Eliminate genotypes but do not lump alleles at each locus.
  #
  keyword["eliminate_genotypes"] = true
  #
  # Define a gamete competition data frame.
  #
  gamete_competition_frame = DataFrame(Marker = ASCIIString[],
    LowAllele = ASCIIString[], LowTau = Float64[],
    HighAllele = ASCIIString[], HighTau = Float64[], Pvalue = Float64[])
  #
  # Subject each locus to gamete competition analysis.
  #
  (model_loci, locus.model_loci) = (locus.model_loci, 1)
  model_locus = copy(locus.model_locus)
  locus.model_locus = zeros(Int, 1)
  for loc = 1:locus.loci
    locus.model_locus[1] = loc
    #
    # Define the parameter data structure.
    #
    keyword["constraints"] = 1
    keyword["goal"] = "maximize"
    keyword["parameters"] = locus.alleles[loc]
    keyword["title"] = "Gamete competition analysis for " * locus.name[loc]
    parameter = set_parameter_defaults(keyword)
    parameter = initialize_optimization(locus, parameter, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree, 
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println("One or more pedigrees exceeds the complexity threshold.")
      println("$locus.name[loc] is being skipped.")
      println(io, "One or more pedigrees exceeded the complexity threshold.")
      println(io, "$locus.name[loc] was skipped.")
      skipped_loci = skipped_loci + 1
      continue
    end
    #
    # Pass the variables to optimize for maximum likelihood estimation.
    #
    function fun(par)
      copy!(parameter.par, par)
      f = elston_stewart_loglikelihood(pedigree, person, locus, parameter, 
        instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = optimize(fun, parameter)
    (low_tau, low_allele) = findmin(best_par)
    (high_tau, high_allele) = findmax(best_par)
    #
    # Report the pvalue of likelihood ratio test statistic.
    #
    lrt = 2 * (parameter.function_value[2] - parameter.function_value[1])
    degrees_of_freedom = locus.alleles[loc] - 1
    pvalue = ccdf(Chisq(degrees_of_freedom), lrt)
    push!(gamete_competition_frame, [locus.name[loc], 
      locus.allele_name[loc][low_allele], low_tau, 
      locus.allele_name[loc][high_allele], high_tau, pvalue])
  end
  writetable(keyword["gamete_competition_table"], gamete_competition_frame)
  show(gamete_competition_frame)
  close(io)
  return skipped_loci
end # function gamete_competition_option

end # module MendelGameteCompetition

