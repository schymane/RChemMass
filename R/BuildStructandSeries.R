# Functions for building structures and series
# Compiled from various scripts.
# E. Schymanski, 15/2/2017

# This should go in DESCRIPTION
# library(enviPat)
# library(rcdk)
# library(RMassBank)
# library(mzR)
# library(curl)
# library(rsvg)
# source("C:/DATA/R/cachedEic.R")
# source("C:/DATA/R/cachedMsms.R")
# source("C:/DATA/R/plotMSMS.R")
data(isotopes, envir=asNamespace("RChemMass"))
data(adducts, envir=asNamespace("RChemMass"))
m_e <- 5.45799e-4

#' Build SMILES from R groups
#'
#' @description This function performs basic string manipulation on SMILES codes to generate
#' systematic homologues, facilitaing the generation of identifiers for structures
#' that may not exist in databases for use in MS workflows.
#' This requires a decent understanding of SMILES as the text manipulation can have
#' unexpected consequences. This does not perform any advanced chemical processing.
#' Using large ranges over many R groups will result in combinatorial explosion; for
#' these cases structure generation is recommended. SMILES are built with the
#' \code{RtoSMILES} function.
#'
#' @usage buildSmiles(genSmiles, R1toN, nR1toN, ExtraAtoms_R1toN=NULL, RDB_R1toN=NULL)
#'
#' @param genSmiles A "Generic" \code{SMILES} code to expand, with R groups encoded in square brackets.
#' The first group can be \code{"[R]"} or \code{"[R1]"}, subsequent groups should be numbered
#' systematically (e.g. \code{"[R2]", "[R3]", ...}). Examples \code{"O=S(O)(=O)OCCO[R1]C[R2]", "CCO[R]CCO"}
#' @param R1toN A vector containing 1 to n \code{SMILES} entries representing the repeating unit. One entry
#' per R group is required, comma separated and enclosed in brackets. Examples \code{"(C,C,C)", "(C,CCO)"}.
#' @param nR1toN The range for each R group (start-end,start-end, ...,start-end).
#' If a single range is given for n R groups, this should be split in advance with \code{splitRrange}
#' @param ExtraAtoms_R1toN (optional, to be removed from this function and handled in advance).
#' Extra atoms to be added to beginning or end of R group. Must be valid \code{SMILES} to generate
#' valid structures. Use with caution.
#' @param RDB_R1toN (optional, to be removed from this function and handled in advance).
#' Unsaturation to be added to beginning or end of R group. Must be valid \code{SMILES} to generate
#' valid structures. Use with caution - double-bond stereochemistry is not created.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{RtoSMILES}}, \code{\link{splitRrange}}, \code{\link{adjustRgroup}},
#' To view created SMILES: \code{\link{renderSMILES.CDKdepict}}, \code{\link{renderSMILES.rcdk}}.
#'
#' @return Returns a vector containing the resulting SMILES codes
#' @export
#'
#' @examples
#' test_smiles <- "P(=O)(OC[C@H](COP(=O)(O)OC[C@H](OC(=O)[R2])COC(=O)[R1])O)(O)OC[C@H](OC(=O)[R4])COC(=O)[R3]"
#' buildSmiles(test_smiles,"(C,C,C,C)","(1-3,1-3,1-3,1-3)")
#'
buildSmiles <- function(genSmiles, R1toN, nR1toN, ExtraAtoms_R1toN=NULL, RDB_R1toN=NULL,
                        atStart=TRUE, adjustRange=FALSE) {
  # process the inputs
  smiles <- as.character(genSmiles)
  #  r_groups <- strsplit(gsub("[/(/)]","",R1toN),",")[[1]]
  r_groups <- strsplit(substring(R1toN,2,(nchar(R1toN)-1)),",")[[1]]
  # n_r_groups <- strsplit(gsub("[/(/)]","",nR1toN),",")[[1]]
  n_r_groups <- strsplit(substring(nR1toN,2,(nchar(nR1toN)-1)),",")[[1]]
  n_r_groups <- matrix(as.numeric(unlist(strsplit(n_r_groups,"-"))),ncol=2,byrow=T)
  colnames(n_r_groups) <- c("Rmin","Rmax")
  if (!is.null(ExtraAtoms_R1toN)) {
    extraAtoms <- strsplit(substring(ExtraAtoms_R1toN,2,(nchar(ExtraAtoms_R1toN)-1)),",")[[1]]
    if (length(r_groups)!=length(extraAtoms)) {
      warning(paste0("Mismatch in number of R groups (", length(r_groups),") and extra atoms (",
              length(extraAtoms), "). Check results for ",smiles))
    }
  } else {
    extraAtoms <- "" # otherwise leave empty
  }
  if (!is.null(RDB_R1toN)) {
    RDB_entries <- as.numeric(strsplit(substring(RDB_R1toN,2,(nchar(RDB_R1toN)-1)),",")[[1]])
    if (length(r_groups)!=length(RDB_entries)) {
      warning(paste0("Mismatch in number of R groups (", length(r_groups),") and RDB entries (",
                     length(RDB_entries), "). Check results for ",smiles))
    }
  } else {
    RDB_entries <- "" # otherwise leave empty
  }  #
  # calculate total number of molecules covered
  n_r_groups_range <- n_r_groups[,2] - n_r_groups[,1] + 1
  n_smiles <- prod(n_r_groups_range)
  smiles_list <- vector("character",n_smiles)
  smiles_index <- 0
  smiles_index2 <- 0
  temp_smiles_list <- ""
  next_temp_list <- ""
  #   r2 <- grepl("[R2]",smiles,fixed=TRUE)
  r1 <- grepl("[R1]",smiles,fixed=TRUE)
  r0 <- grepl("[R]",smiles,fixed=TRUE)
  r <- grepl("R",smiles,fixed=TRUE)
  if (r1) {
    first_r_group <- "[R1]"
  } else {
    first_r_group <- "[R]"
  }
  r_group <- ""
  if (!r) {
    stop("SMILES does not contain an R group")
  }
  if (!r1 && !r0) {
    stop("SMILES does not contain an R1 or undefined R group")
  }
  # first add extra atoms, RDB etc to generic SMILES
  if (!is.null(ExtraAtoms_R1toN) || !is.null(RDB_R1toN)) {
    for (n in 1:length(n_r_groups_range)) {
      if (n == 1) {
        r_group <- first_r_group
      } else {
        r_group <- paste0("[R",n,"]")
      }
      if (adjustRange) {
        if (!is.null(ExtraAtoms_R1toN) && !is.null(RDB_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, extraAtoms[n],RDB_entries[n],
                                      atStart=atStart,R_range = c(n_r_groups[n,1],n_r_groups[n,2]))
        } else if (!is.null(ExtraAtoms_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, extraAtoms[n],atStart=atStart)
          temp_smiles <- c(temp_smiles[1], n_r_groups[n,1],n_r_groups[n,2])
        } else if (!is.null(RDB_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, RDB = RDB_entries[n],atStart=atStart,
                                      R_range = c(n_r_groups[n,1],n_r_groups[n,2]))
        } else {
          temp_smiles <- c(smiles, n_r_groups[n,1],n_r_groups[n,2])
        }
        smiles <- temp_smiles[1]
        n_r_groups[n,1] <- as.numeric(temp_smiles[2])
        n_r_groups[n,2] <- as.numeric(temp_smiles[3])
      } else {
        if (!is.null(ExtraAtoms_R1toN) && !is.null(RDB_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, extraAtoms[n],RDB_entries[n],atStart=atStart)
        } else if (!is.null(ExtraAtoms_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, extraAtoms[n],atStart=atStart)
        } else if (!is.null(RDB_R1toN)) {
          temp_smiles <- adjustRgroup(smiles, r_group, RDB=RDB_entries[n],atStart=atStart)
        } else {
          temp_smiles <- smiles
        }
        smiles <- temp_smiles[1]
      }
    }
  }

  # then expand the SMILES by R groups
  for (n in 1:length(n_r_groups_range)) {
    # first case: R=1. Either R or R1.
    if (n == 1) {
      if (r1) {
        r_group <- "[R1]"
      } else {
        r_group <- "[R]"
      }
      # adjust the SMILES by the R group
      temp_smiles_list <- RtoSMILES(smiles,r_group,r_groups[n],c(n_r_groups[n,1],n_r_groups[n,2]))
    } else {
      r_group <- paste0("[R",n,"]")
      next_temp_list <- vector("character",prod(n_r_groups_range[1:n]))
      for (j in 1:length(temp_smiles_list)) {
        smiles_index <- smiles_index+1
        smiles_index2 <- smiles_index+n_r_groups_range[n]-1
        next_temp_list[smiles_index:smiles_index2] <- RtoSMILES(temp_smiles_list[j],r_group,r_groups[n],
                                                                c(n_r_groups[n,1],n_r_groups[n,2]))
        smiles_index <- smiles_index2
      }
      temp_smiles_list <- next_temp_list
      next_temp_list <- ""
      smiles_index <- 0
      smiles_index2 <- 0
    }
  }
  return(temp_smiles_list)
}



#' Split a single R group range into n ranges
#'
#' @description This splits a range into n for the case where e.g. the total number of carbons is
#' known but multiple Rs exist. This assumes symmetrical filling, with remainders
#' added to the earlier groups. Used to feed into the \code{buildSmiles} function
#' This function may perform strangely and needs testing.
#'
#' @usage splitRrange(genSmiles, nR1toN)
#'
#' @param genSmiles A "Generic" \code{SMILES} code to expand, with R groups encoded in square brackets.
#' The first group can be \code{"[R]"} or \code{"[R1]"}, subsequent groups should be numbered
#' systematically (e.g. \code{"[R2]", "[R3]", ...}). Examples \code{"O=S(O)(=O)OCCO[R1]C[R2]", "CCO[R]CCO"}
#' @param nR1toN The range for each R group (start-end,start-end, ...,start-end) or a single range.
#' If a single range is given for n R groups, this is split to form a range for each R group.
#' Warnings are produced if mismatches occur, in which case the input is returned.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{RtoSMILES}}, \code{\link{buildSmiles}}.
#'
#' @return Returns the split range or original input for use in \code{\link{buildSmiles}}
#' @export
#'
#' @examples
#' test_smiles <- "P(=O)(OC[C@H](COP(=O)(O)OC[C@H](OC(=O)[R2])COC(=O)[R1])O)(O)OC[C@H](OC(=O)[R4])COC(=O)[R3]"
#' splitRrange(test_smiles,"(1-37)")
#' splitRrange(test_smiles,"(1-38)")
#' splitRrange(test_smiles,"(1-39)")
#' splitRrange(test_smiles,"(1-40)")
#' splitRrange(test_smiles,"(3-40)")
#' splitRrange(test_smiles,"(1-3)")
#' splitRrange(test_smiles,"(16-22)")
#'
splitRrange <- function(genSmiles, nR1toN) {
  # function tests whether Rrange needs to be split
  nR <- length(gregexpr("R",genSmiles,fixed=TRUE)[[1]])
  n_r_groups <- strsplit(substring(nR1toN,2,(nchar(nR1toN)-1)),",")[[1]]
  #n_r_groups <- strsplit(gsub("[/(/)]","",nR1toN),",")[[1]]
  n_r_groups <- matrix(as.numeric(unlist(strsplit(n_r_groups,"-"))),ncol=2,byrow=T)
  colnames(n_r_groups) <- c("Rmin","Rmax")
  n_ranges <- length(n_r_groups[,1])
  if (nR == n_ranges) {
    # no need to convert, just output
    nR1toN_out <- nR1toN
  } else if (n_ranges == 1 && nR > 1) {
    n_r_groups_range <- n_r_groups[,2] - n_r_groups[,1] + 1
    start_num <- vector("numeric",nR)
    end_num <- vector("numeric",nR)
    nR1toN_out <- "("
    start_rem <- n_r_groups[1,1]%%nR
    end_rem <- n_r_groups[1,2]%%nR
    for (i in 1:nR) {
      if (i == 1) {
        start_num[i] <- ceiling(n_r_groups[1,1]/nR)
        end_num[i] <- ceiling(n_r_groups[1,2]/nR)
        nR1toN_out <- paste0(nR1toN_out,start_num[i],"-",end_num[i],",")
      } else if (i == nR) {
        start_num[i] <- floor(n_r_groups[1,1]/nR)
        end_num[i] <- floor(n_r_groups[1,2]/nR)
        nR1toN_out <- paste0(nR1toN_out,start_num[i],"-",end_num[i],")")
      } else {
        if (i <= start_rem) {
          start_num[i] <- floor(n_r_groups[1,1]/nR) + 1
        } else {
          start_num[i] <- floor(n_r_groups[1,1]/nR)
        }
        if (i <= end_rem) {
          end_num[i] <- floor(n_r_groups[1,2]/nR) + 1
        } else {
          end_num[i] <- floor(n_r_groups[1,2]/nR)
        }
        nR1toN_out <- paste0(nR1toN_out,start_num[i],"-",end_num[i],",")
      }
    }
    #nR1toN_out <- paste0(nR1toN_out,")")
    range_out <- sum(end_num)-sum(start_num)+1
    if (!(n_r_groups_range == range_out)) {
      warning(paste0("Mismatch in range: input=",n_r_groups_range,", output=",range_out))
    }
    if (range_out < nR) {
      warning(paste0("Range (", range_out, ") is less than number of R groups (",nR,")"))
    }

  } else {
    warning(paste0("Mismatch in conditions, check nR1toN entry. nR=",nR," , nR1toN=",nR1toN))
    nR1toN_out <- nR1toN
  }
  return(nR1toN_out)
}



#' Build SMILES from single R groups
#'
#' @description This works behind the scenes of \code{buildSmiles} on individual R groups, but can be
#' used independently. This performs basic string manipulation on SMILES codes to generate
#' systematic homologues, facilitaing the generation of identifiers for structures
#' that may not exist in databases for use in MS workflows.
#' This requires a decent understanding of SMILES as the text manipulation can have
#' unexpected consequences.
#'
#' @usage RtoSMILES(genSmiles, R_format, R_smiles, R_range)
#'
#' @param genSmiles A "Generic" \code{SMILES} code to expand, with R groups encoded in square brackets.
#' The R group to replace is defined in R_format.
#' @param R_format The format of the R group to replace, including square brackets, e.g.
#' \code{[R], [R1], [R2]}.
#' @param R_smiles The valid \code{SMILES} code representing the repeating unit. Examples
#' \code{"C", "CCO"}. Must be valid \code{SMILES} to generate valid structures.
#' Use with caution - it is highly recommended to view the structures with
#' \code{\link{renderSMILES.CDKdepict}} or \code{\link{renderSMILES.rcdk}} to detect unwanted
#' behaviour.
#' @param R_range The range for each R group (start-end).
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{BuildSmiles}}, \code{\link{splitRrange}}, \code{\link{adjustRgroup}},
#' To view created SMILES: \code{\link{renderSMILES.CDKdepict}}, \code{\link{renderSMILES.rcdk}}.
#'
#' @return Returns a vector containing the resulting SMILES codes
#' @export
#'
#' @examples
#' RtoSMILES("OC(=O)[R]","[R]","C",c(1,5))
#' test_smiles <- "P(=O)(OC[C@H](COP(=O)(O)OC[C@H](OC(=O)[R2])COC(=O)[R1])O)(O)OC[C@H](OC(=O)[R4])COC(=O)[R3]"
#' RtoSMILES(test_smiles,"[R1]","C",c(1,5))
#'
RtoSMILES <- function(genSmiles, R_format, R_smiles, R_range) {
  # return either one smiles, or a vector of length R_smiles
  R_group <- R_format
  Rmin <- min(R_range)
  Rmax <- max(R_range)
  nRs <- Rmax-Rmin+1
  # check if square brackets are there, if not, add
  if (!grepl("[",R_group,fixed=TRUE) && !grepl("]",R_group,fixed=TRUE)) {
    R_group <- paste0("[",R_group,"]")
  }
  r_test <- grepl(R_group,genSmiles,fixed=TRUE)
  if (!r_test) {
    stop(paste0("Given R format not found in SMILES: ",R_group))
  }
  smiles_list <- vector("character",nRs)
  start_smiles <- vector("character",Rmin)
  start_smiles[1:Rmin] <- R_smiles
  start_smiles <- paste0(start_smiles[1:Rmin],collapse="")
  r_list <- vector("character",nRs)
  # now fill the smiles_list..
  for (i in 1:nRs) {
    if (i == 1) {
      r_list[i] <- start_smiles
      smiles_list[i] <- sub(R_group,r_list[i],genSmiles,fixed=TRUE)
    } else {
      r_list[i] <-  paste0(r_list[(i-1)],R_smiles)
      smiles_list[i] <- sub(R_group,r_list[i],genSmiles,fixed=TRUE)
    }
  }
  return(smiles_list)
}


#' Add extra atoms and unsaturations before or after an R group
#'
#' @description This works behind the scenes of \code{buildSmiles} on individual R groups, but can be
#' used independently. This performs basic string manipulation on SMILES codes to add
#' additional atoms (as valid SMILES codes) and unsaturations before or after R groups
#' prior to generation of systematic homologues with RtoSMILES. This facilitates the generation
#' of identifiers for structures that may not exist in databases for use in MS workflows.
#' This requires a decent understanding of SMILES as the text manipulation can have
#' unexpected consequences.
#' It is highly recommeded to view these structures to ensure the groups are added as intended.
#'
#' @usage adjustRgroup(genSmiles,R_format,ExtraAtoms=NULL,RDB=NULL,R_range=NULL,atStart=TRUE)
#'
#' @param genSmiles A "Generic" \code{SMILES} code to expand, with R groups encoded in square brackets.
#' The R group to replace is defined in R_format.
#' @param R_format The format of the R group to replace, including square brackets, e.g.
#' \code{[R], [R1], [R2]}.
#' @param ExtraAtoms A valid \code{SMILES} code representing the single unit to add.
#' Single-valent atoms should be included with their attached atom. Examples
#' \code{"C(=O)", "C(Br)", "C(O)", "O"}. Must be valid \code{SMILES} to generate valid structures.
#' Use with caution - it is highly recommended to view the structures with
#' \code{\link{renderSMILES.CDKdepict}} or \code{\link{renderSMILES.rcdk}} to detect unwanted
#' behaviour.
#' @param RDB A number (Ring and Double Bond count) indicating the degree of unsaturations to add.
#' This number determines how many \code{"C=C"} units are added. Note that no double-bond stereoisomerism
#' is added. If \code{R_range} is given, this will be adjusted according to the number of \code{"C=C"} groups
#' added, unless a negative range would be generated, in this case the original range is returned.
#' @param R_range The range for the R group (start-end). Will be corrected by \code{2*RDB} if given.
#' @param atStart Default \code{TRUE} adds extra atoms and unsaturated units before the R group,
#' \code{FALSE} adds afterwards. Triply substituted Cs should only be added after terminal Rs, for instance.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{BuildSmiles}}, \code{\link{splitRrange}}, \code{\link{RtoSMILES}},
#' To view created SMILES: \code{\link{renderSMILES.CDKdepict}}, \code{\link{renderSMILES.rcdk}}.
#'
#' @return Returns a vector containing the resulting SMILES codes and,
#' if range is given, the new min and max values for range.
#' @export
#'
#' @examples
#'
#' test_genSmiles <- adjustRgroup("OC(=O)[R]","[R]",ExtraAtoms="C(=O)")
#' test_genSmiles <- adjustRgroup("OC(=O)[R]","[R]",ExtraAtoms="C(Br)")
#' test_genSmiles <- adjustRgroup("OC(=O)[R]","[R]",RDB=2,atStart=FALSE)
#' buildCDKdepictURL(test_genSmiles)
#' # if the range is returned, need to only depict first entry
#' test_genSmiles <- adjustRgroup("OC(=O)[R]","[R]",ExtraAtoms="C(=O)",RDB=2,atStart=TRUE, R_range = c(4,6))
#' test_genSmiles <- adjustRgroup("OC(=O)[R]","[R]",RDB=2,atStart=TRUE, R_range = c(2,20))
#' buildCDKdepictURL(test_genSmiles[1])
#'
adjustRgroup <- function(genSmiles,R_format,ExtraAtoms=NULL,RDB=NULL,R_range=NULL,atStart=TRUE) {
  if (is.null(ExtraAtoms) && is.null(RDB)) {
    warning("No adjustments to be made, returning input SMILES")
    new_genSmiles <- genSmiles
  }
  if (length(grep(R_format,genSmiles,fixed=TRUE))==0) {
    error(paste0("R format ", R_format," not found in SMILES, please check and try again"))
  }
  RDB_adjust <- ""
  # adjust according to RDB
  if (!is.null(RDB) && is.numeric(RDB) && RDB>=1) {
    for (i in 1:RDB) {
      RDB_adjust <- paste0(RDB_adjust,"C=C")
    }
    if (!is.null(R_range)) {
      new_R_range <- R_range-2*RDB
      if (min(new_R_range) < 0) {
        warning("New range is negative, keeping original values")
      } else {
        R_range <- new_R_range
      }
    }
  }
  if (!is.null(RDB) && is.null(R_range)) {
    # cannot correct the R_range
    warning(paste0("No R_range given, cannot adjust range. Chain will grow by ",2*RDB,"Cs beyond range"))
  }
  # so, new RDB string is in RDB_adjust, now need ExtraAtoms
  ExtraAtoms_toAdd <- ""
  if (!is.null(ExtraAtoms)) {
    ExtraAtoms_toAdd <- paste0(ExtraAtoms)
  }

  if (atStart) {
    string_to_add <- paste0(ExtraAtoms_toAdd,RDB_adjust,R_format)
    new_genSmiles <- sub(R_format,string_to_add,genSmiles,fixed=TRUE)
  } else {
    string_to_add <- paste0(R_format,ExtraAtoms_toAdd,RDB_adjust)
    new_genSmiles <- sub(R_format,string_to_add,genSmiles,fixed=TRUE)
  }
  # if R_range is not defined, only one entry is returned
  return(c(new_genSmiles,R_range))
}


#' Trim SMILES list from buildSmiles to one SMILES per formula
#'
#' @description This function takes the output of \code{\link{buildSmiles}} and
#' returns the first SMILES entry per molecular formula.
#'
#' @usage trimBuiltSmiles(smiles_list)
#'
#' @param smiles_list A list of \code{SMILES} from \code{\link{buildSmiles}}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{BuildSmiles}}, \code{\link{MolFormFromSmiles.rcdk}}.
#' To view created SMILES: \code{\link{renderSMILES.CDKdepict}}, \code{\link{renderSMILES.rcdk}}.
#'
#' @return Returns a vector containing the resulting SMILES codes
#' @export
#'
#' @examples
#' LAS_gen_smiles <- "OS(=O)(=O)c1ccc(cc1)C([R1])[R2]"
#' LAS_smiles <- buildSmiles(LAS_gen_smiles, "{C,C}", splitRrange(LAS_gen_smiles, "{9-13}"))
#' trimBuiltSmiles(LAS_smiles)
#'
trimBuiltSmiles <- function(smiles_list) {
  formulas <- as.vector(sapply(smiles_list, MolFormFromSmiles.rcdk))
  smiles_trim <- LAS_smiles[!duplicated(formulas)]
  return(smiles_trim)
}


#' Build Labels for Homologue Series
#'
#' @description This function generates labels to plot homologue series with the inputs used
#' in the \code{buildSmiles} function.
#'
#' @usage buildSeriesLabel(genSmiles, SeriesName, R1toN, nR1toN, ExtraAtoms_R1toN=NULL,
#' RDB_R1toN=NULL, split_nR1toN=FALSE)
#'
#' @param genSmiles A "Generic" \code{SMILES} code to expand, with R groups encoded in square brackets.
#' @param SeriesName A general name for the series, preferably short (abbreviation).
#' @param R1toN A vector containing 1 to n \code{SMILES} entries representing the repeating unit. One entry
#' per R group is required, comma separated and enclosed in brackets. Examples \code{"(C,C,C)", "(C,CCO)"}.
#' @param nR1toN The range for each R group (start-end,start-end, ...,start-end).
#' If a single range is given for n R groups, this can optionally be split by setting \code{split_nR1toN}.
#' Different labels are produced depending on this setting.
#' @param ExtraAtoms_R1toN (optional) If present, is added to the end of the label.
#' @param RDB_R1toN (optional) If present, is added to the end of the label (after ExtraAtoms)
#' @param split_nR1toN If \code{TRUE}, splits \code{nR1toN} with \code{\link{splitRrange}}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{buildSmiles}}, \code{\link{splitRrange}}.
#'
#' @return Returns a label that can be used to annotate generic structures in CDK Depict
#' @export
#'
#' @examples
#' genSmiles <- "OS(=O)(=O)c1ccc(cc1)C([R1])[R2]"
#' buildSeriesLabel(genSmiles,"LAS","{C,C}","{9-13}")
#' buildSeriesLabel(genSmiles,"LAS","{C,C}","{9-13}", split_nR1toN=TRUE)
#'
buildSeriesLabel <- function(genSmiles, SeriesName, R1toN, nR1toN, ExtraAtoms_R1toN=NULL,
                             RDB_R1toN=NULL, split_nR1toN=FALSE) {
  # count the R groups
  nR <- length(gregexpr("R",genSmiles)[[1]])
  # extract the SMILES codes for each R
  Rsmiles <- strsplit(substring(R1toN,2,(nchar(R1toN)-1)),",")[[1]]
  # extract the range(s)
  if (split_nR1toN) {
    nR1toN_in <- splitRrange(genSmiles,nR1toN)
    Rrange <- strsplit(substring(nR1toN_in,2,(nchar(nR1toN_in)-1)),",")[[1]]
  } else {
    nR1toN_in <- nR1toN
    Rrange <- strsplit(substring(nR1toN_in,2,(nchar(nR1toN_in)-1)),",")[[1]]
  }
  # create the R part of the label
  if (length(Rrange)<nR) {
    R_label <- paste(paste0("R",as.character(seq(1:nR)),"=",Rsmiles),collapse=",")
    full_label <- paste0(SeriesName,"; ", R_label,"; TotalR=", Rrange)
  } else {
    R_label <- paste(paste0("R",as.character(seq(1:nR)),"=",Rsmiles,"(",Rrange,")"),collapse=",")
    full_label <- paste0(SeriesName,"; ", R_label)
  }
  # add extra atoms or RDB, if needed
  if (!is.null(ExtraAtoms_R1toN)) {
    ExtraAtom_text <- substring(ExtraAtoms_R1toN,2,(nchar(ExtraAtoms_R1toN)-1))
    full_label <- paste0(full_label, "; ExtraAtoms=",ExtraAtom_text)
  }
  if (!is.null(RDB_R1toN)) {
    RDB_text <- substring(RDB_R1toN,2,(nchar(RDB_R1toN)-1))
    full_label <- paste0(full_label, "; RDB=",RDB_text)
  }
  return(full_label)
}


#' Build Labels for Members of a Homologue Series
#'
#' @description This function generates labels for members of a homologue series generated with
#' the \code{buildSmiles} function (or also in general).
#'
#' @usage buildMemberLabels(SeriesName, n_start, n_total, smiles_list=NULL, addFormula=FALSE)
#'
#' @param SeriesName A general name for the series, preferably short (abbreviation).
#' @param n_start The starting number for the members in the series.
#' @param nR1toN The total number of members that need names.
#' @param smiles_list (optional) If present, used to calculate formulas for \code{addFormula=TRUE}.
#' @param RDB_R1toN (optional) If \code{TRUE}, adds molecular formula to the labels.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{buildSmiles}}, \code{\link{buildSeriesLabel}}.
#'
#' @return Returns a label that can be used to annotate individual structures in CDK Depict
#' @export
#'
#' @examples
#' genSmiles <- "OS(=O)(=O)c1ccc(cc1)C([R1])[R2]"
#' nR1toN_in <- splitRrange(genSmiles,"{9-13}")
#' LAS_smiles <- buildSmiles(genSmiles,"{C,C}",nR1toN_in)
#' LAS_smiles_trim <- trimBuiltSmiles(LAS_smiles)
#'
#' # now create labels
#' buildMemberLabels("LAS",9,5)
#' LAS_trim_labels <- buildMemberLabels("LAS",9,5, LAS_smiles_trim, addFormula = TRUE)
#' # to write a file that can be copied directly to CDK depict:
#' # write.table(cbind(LAS_smiles_trim, LAS_trim_labels),"LAS_trim_smiles_labels.smi",
#' # quote=F, row.names=F, col.names=F)
#'
buildMemberLabels <- function(SeriesName, n_start, n_total, smiles_list=NULL,
                              addFormula=FALSE) {
  # calculate formulas, if needed
  if (addFormula && !is.null(smiles_list)) {
    MemberFormulas <- as.vector(sapply(smiles_list, MolFormFromSmiles.rcdk))
  } else {
    MemberFormulas <- NULL
  }
  # create labels
  n_formulas <- length(MemberFormulas)
  if (addFormula && (n_formulas < n_total)) {
    warning("Number of formulas less than series length, setting addFormula=FALSE")
    addFormula=FALSE
  }
  if (addFormula) {
    memberLabels <- paste0(SeriesName, "-", seq(from=n_start, length.out=n_total), ": ", MemberFormulas)
  } else {
    memberLabels <- paste0(SeriesName, "-", seq(from=n_start, length.out=n_total))
  }
  return(memberLabels)
}


#' Build molecular formulae for a homologue series
#'
#' @description Given a starting formula and a "building block", this builds a series for formulas
#' from \code{n_start} to \code{n_max}, adding a retention time for direct use in
#' suspect screening. Naming is built on \code{series_name} and \code{n}. Formula
#' manipulations are performed with \code{\link{enviPat}} and \code{\link{RMassBank}}.
#' Can be used if generic structures are not available for \code{BuildSmiles}.
#'
#' @usage build.homol(start_formula, building_block, series_name, n_start, n_max, rt=15)
#'
#' @param start_formula Molecular formula to start building the series. Standardised with
#' \code{checkform}, along with \code{building_block}.
#' @param building_block Formula of the "building block" to add n times to \code{start_formula}.
#' @param series_name Codename for the series. Short names are recommended for downstream use.
#' @param n_start Starting number for member naming only. Formula building goes from n=1 (
#' \code{start_formula}) to n=(\code{n_max}-\code{n_start}).
#' @param n_max Maximum number to calculate.
#' @param rt Default \code{15} (minutes); can be adjusted for different systems and anticipated
#' series behaviour as desired. For first round screening, RT in the middle of the elution is
#' recommended, such that a symmetrical window over the whole program can be used.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{BuildSmiles}}, \code{\link{check_chemform}}, \code{\link{mergeform}},
#' \code{\link{enviPat}}, \code{\link{findMz.formula}}.
#'
#' @return Returns a \code{n_max-n_start} by \code{4} array containing member names, formulas, RTs
#' and neutral monoisotopic mass.
#' @export
#'
#' @examples
#' las <- build.homol("C16H26O3S","CH2","LAS",10,14)
#' spacs <- build.homol("C9H10O5S","CH2","SPAC",3,15)
#' spadcs <- build.homol("C9H8O7S","CH2","SPADC",1,15)
#' c12aes <- build.homol("C14H30O5S","C2H4O","C12AES",1,10)
#' c13aes <- build.homol("C15H32O5S","C2H4O","C13AES",1,6)
#' npeos <- build.homol("C17H28O5S","C2H4O","NPEO_SO4_",1,8)
#'
build.homol <- function(start_formula, building_block, series_name, n_start, n_max, rt=15) {
  formulas <- c(start_formula, building_block)
  checked <- check_chemform(isotopes,formulas)
  member_formulas <- checked[1,2]
  n <- 1
  new_formula <- ""
  neutral_monoiso_mass <- as.numeric(findMz.formula(member_formulas[1],"")[3])
  #series_MmHm <- as.numeric(findMz.formula(series_formulas[1],"mH")[3])
  member_names <- paste(series_name,as.character(n_start),sep="")
  rt <- rt # in minutes
  for(i in 2:(n_max-n_start)) {
    #n <- n+1
    formulas[1] <- member_formulas[n]
    checked <- check_chemform(isotopes,formulas)
    new_formula <- mergeform(checked[1,2],checked[2,2])
    member_formulas[i] <- new_formula
    neutral_monoiso_mass[i] <- as.numeric(findMz.formula(new_formula,"")[3])
    #series_MmHm[i] <- as.numeric(findMz.formula(new_formula,"mH")[3])
    member_names[i] <- paste(series_name,as.character(n_start+i-1),sep="")
    rt[i] <- rt[1]
    n <- n+1
  }
  homol_series <- cbind(member_names, member_formulas, rt, neutral_monoiso_mass)
  return(homol_series)
}


#' Build Exact Masses for Homologue Series Screening
#'
#' @description This builds a series of exact masses given a starting mass and a
#' mass difference and the range \code{n_start} to \code{n_end}, adding a retention
#' time for direct use in suspect screening. Naming is built on \code{series_name} and
#' \code{n_start+n}.
#'
#' @usage build.homol.byMass(start_mass, mass_diff, series_name, n_start, n_end, n_adj=0, rt=15)
#'
#' @param start_mass Starting mass for the series. \code{n_start} can be negative if this
#' is a "central" mass.
#' @param mass_diff Mass of the "building block" for the series.
#' @param series_name Codename for the series. Short names are recommended for downstream use.
#' @param n_start Starting number for series, if negative, \code{mass_diff} is subtracted.
#' @param n_end End number for series.
#' @param n_adj This adjusts the number for naming.
#' @param rt Default \code{15} (minutes); can be adjusted for different systems and anticipated
#' series behaviour as desired. For first round screening, RT in the middle of the elution is
#' recommended, such that a symmetrical window over the whole program can be used. This
#' can also be achieved with \code{rt_window=NULL} in \code{\link{plotEICs}}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{build.homol}}, \code{\link{screen.homol}}.
#'
#' @return Returns a \code{n_end-n_start} by \code{3} array containing member names, masses and RTs.
#'
#' @export
#'
#' @examples
#' findMz.formula("C16H26O3S1", mode="mH")
#' las_MmHm <- build.homol.byMass(297.156,14.01565,"LAS",0,4,10)
#' findMz.formula("C9H10O5S", mode="mH")
#' spacs_MmHm <- build.homol.byMass(229.0199,14.01565,"SPAC",0,12,3)
#' c12aes <- build.homol("C14H30O5S","C2H4O","C12AES",1,10)
#' findMz.formula("C14H30O5S", mode="mH")
#' c12aes_MmHm <- build.homol.byMass(309.1772,44.026215,"C12AES",0,9,1)
#'
build.homol.byMass <- function(start_mass, mass_diff, series_name, n_start,
                               n_end, n_adj=0, rt=15) {
  n <- 1
  exact_mass <- start_mass+n_start*mass_diff
  member_names <- paste(series_name,as.character(n_start+n_adj),sep="")
  rt <- rt # in minutes
  for(i in 2:(n_end-n_start)) {
    #n <- n+1
    exact_mass[i] <- start_mass+(n_start+n)*mass_diff
    member_names[i] <- paste(series_name,as.character(n_start+n+n_adj),sep="")
    rt[i] <- rt[1]
    n <- n+1
  }
  homol_series <- cbind(member_names, exact_mass, rt)
  return(homol_series)
}

