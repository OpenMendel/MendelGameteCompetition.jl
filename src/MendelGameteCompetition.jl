__precompile__()

"""
This module orchestrates a gamete competition analysis.
"""
module MendelGameteCompetition
#
# Required OpenMendel packages and modules.
#
using MendelBase
# namely: DataStructures, ModelConstructions,
# ElstonStewartPreparations, ElstonStewartEvaluations
using MendelSearch
#
# Required external modules.
#
using CSV
using DataFrames
using Distributions
using LinearAlgebra

export GameteCompetition
"""
This is the wrapper function for the Gamete Competition analysis option.
"""
function GameteCompetition(control_file = ""; args...)
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Gamete Competition analysis option")
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["gamete_competition_output_table"] =
    "gamete competition Output Table.txt"
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
    locus_frame, phenotype_frame, person_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Check if SNP data were read.
  #
  if snpdata.snps != 0
    println(" \n\nERROR: This analysis does not use data from SNP files!\n")
  else
  #
  # Execute the specified analysis.
  #
    println(" \nAnalyzing the data.\n")
    execution_error = false
    skipped_loci = gamete_competition_option(pedigree, person, nuclear_family,
      locus, locus_frame, phenotype_frame, person_frame, keyword)
    if execution_error
      println(" \n \nERROR: Mendel terminated prematurely!\n")
    else
      println(" \n \nMendel's analysis is finished.\n")
    end
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
  phenotype_frame::DataFrame, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any})

  io = keyword["output_unit"]
  skipped_loci = 0
  #
  # Eliminate genotypes but do not lump alleles at each locus.
  #
  keyword["eliminate_genotypes"] = true
  #
  # Define a gamete competition data frame.
  #
  gamete_competition_frame = DataFrame(Marker = AbstractString[],
    LowAllele = AbstractString[], LowTau = Float64[],
    HighAllele = AbstractString[], HighTau = Float64[], Pvalue = Float64[])
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
    parameter =
      initialize_optimization_gamete_competition!(locus, parameter, keyword)
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
    # Pass the variables to search for maximum likelihood estimation.
    #
    function fun(par)
      copyto!(parameter.par, par)
      f = elston_stewart_loglikelihood(penetrance_gamete_competition,
        prior_gamete_competition, transmission_gamete_competition,
        pedigree, person, locus, parameter, instruction, person_frame, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = mendel_search(fun, parameter)
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
  gc_table_file = string(keyword["gamete_competition_output_table"])
  CSV.write(gc_table_file, gamete_competition_frame;
    writeheader = true, delim = keyword["output_field_separator"],
    missingstring = keyword["output_missing_value"])
  show(gamete_competition_frame)
  close(io)
  return skipped_loci
end # function gamete_competition_option

"""
Supply a penetrance for individual i.
"""
function penetrance_gamete_competition(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, endlocus::Int, i::Int)

  pen = 1.0
  return pen
end # function penetrance_gamete_competition

"""
Supply a prior probability for founder i.
"""
function prior_gamete_competition(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, endlocus::Int, i::Int)

  prior_prob = 1.0
  for l = startlocus:endlocus
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_gamete_competition

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
"""
function transmission_gamete_competition(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  person_frame::DataFrame, keyword::Dict{AbstractString, Any},
  startlocus::Int, endlocus::Int, i::Int, j::Int)
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  loc = locus.model_locus[startlocus]
  xlinked = locus.xlinked[loc]
  if xlinked && person.male[i] && person.male[j]
    return 1.0
  end
  #
  # Implement the gamete competition model.
  #
  k = multi_genotype[1, startlocus]
  l = multi_genotype[2, startlocus]
  m = gamete[startlocus]
  trans = 0.0
  if person.disease_status[j] == keyword["affected_designator"]
    if k == m; trans = trans + par[k] / (par[k] + par[l]); end
    if l == m; trans = trans + par[l] / (par[k] + par[l]); end
  else
    if k == m; trans = trans + 0.5; end
    if l == m; trans = trans + 0.5; end
  end
  return trans
end # function transmission_gamete_competition

"""
Initialize the optimization problem.
"""
function initialize_optimization_gamete_competition!(locus::Locus,
  parameter::Parameter, keyword::Dict{AbstractString, Any})
  #
  # Initialize, bound, and name the parameters.
  #
  for i = 1:parameter.parameters
    parameter.par[i] = 1.0
    parameter.min[i] = 1e-5
    parameter.name[i] = "tau" * " $i"
    parameter.name[i] = rpad(parameter.name[i], 8, ' ')
  end
  #
  # Constrain the parameter for the most frequent allele.
  #
  loc = locus.model_locus[1]
  n = 0
  sum_freq = 0.0
  for allele = 1:locus.alleles[loc]
    if sum(locus.frequency[loc][:, allele]) > sum_freq
       n = allele
       sum_freq = sum(locus.frequency[loc][:, allele])
    end
  end
  if n > 0; parameter.constraint[1, n] = 1.0; end
  parameter.constraint_level[1] = 1.0
  return parameter
end # function initialize_optimization_gamete_competition!
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelGameteCompetition
