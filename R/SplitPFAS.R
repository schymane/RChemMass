# R Wrapper for the MetFrag Tools SplitPFAS jar
# SplitPFAS Author: Christoph Ruttkies
# R Wrapper: Emma Schymanski
# July 2018

# NOTE: This is under development for a PFAS Classification Project
# Initiated by Zhanyun Wang, ETHZ.
# 
# Jar available for download from:
# https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-Tools.jar

# SplitPFAS Background
# we have four cases
# 1: print help, no conversion - separate function
# 2: just convert smiles with default smarts
# 3: add smarts via a file
# 4: add image to show where bond is split. 

#' Display Help for SplitPFAS Jar
#' 
#' @description A small wrapper function to display the "help" text for SplitPFAS Jar.
#' The jar file is available from 
#' \url{https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-Tools.jar}
#' 
#' @usage splitPFAS.help(splitPFAStools_name, splitPFAStools_dir="")
#'
#' @param splitPFAStools_name The name of the SplitPFAS Jar file (e.g. 
#' \code{"MetFrag2.4.5-Tools.jar"}).
#' @param splitPFAStools_dir The directory containing the jar file; if empty assumes the
#' currect directory.
#'
#' @return Returns the help text in the console
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragTools)
#' and Zhanyun Wang (ETHZ, PFAS examples).
#' 
#' @seealso \code{\link{splitPFAS}}
#' 
#' @export
#'
#' @examples
#' splitPFAStools_name <- "MetFrag2.4.5-Tools.jar"
#' splitPFAS.help(splitPFAStools_name)
#' 
splitPFAS.help <- function(splitPFAStools_name, splitPFAStools_dir="") {
  # send command to the splitPFAS tool and output the help
  curr_dir <- getwd()
  if (nchar(splitPFAStools_dir)>0) {
    setwd(splitPFAStools_dir)
  }
  CommandPart1 <- paste0("java -cp ",splitPFAStools_name, " de.ipbhalle.metfrag.split.SplitPFAS")
  SplitPFASCommand <- CommandPart1
  CommandOut <- suppressWarnings(system(command=SplitPFASCommand,intern=FALSE,show.output.on.console=TRUE))
  setwd(curr_dir)
  #return(CommandOut)
}




#' Split Chemical (as SMILES) into Perfluorinated Part and Functional Group using SplitPFAS
#' 
#' @description A small wrapper function to run the SplitPFAS Jar to split perfluorinated compounds
#' (PFAS) into PFAS chain and functional group, by given SMARTS group(s). NOTE this is under 
#' development: Please contact the authors for more information. The jar is available from 
#' \url{https://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-Tools.jar}
#' 
#' @usage splitPFAS(smiles, splitPFAStools_name, incl_image=FALSE, smarts_file="", 
#' splitPFAStools_dir="")
#'
#' @param smiles The SMILES code of the PFAS to split. Preferably an "MS-ready" form (no salts etc)
#' @param splitPFAStools_name The name of the SplitPFAS Jar file (e.g. 
#' \code{"MetFrag2.4.5-Tools.jar"}).
#' @param incl_image If \code{TRUE}, a PNG image showing the bond broken is saved in a subdirectory 
#' \code{"SplitPFAS_PNGs"} in the current directory. The jar creates a file in the temp directory 
#' which is copied over and deleted from the temp dir. Default \code{FALSE} switches off image creation.
#' @param smarts_file If empty, PFAS are split at the end of the perfluorinated chain. If a path
#' is given, the SMILES is split according to the matching SMARTS in the file, where order of the 
#' SMARTS in the file indicates priority matching order (first match is processed)
#' @param splitPFAStools_dir The directory containing the jar file; if empty assumes the
#' currect directory.
#'
#' @return Returns a list containing the processed SMILES (\code{SMILES_in}), the (non-PFAS) 
#' functional group (\code{Func_Group}), the PFAS Group(s) split off (\code{PFAS_Groups}), 
#' both of which are pipe-separated (\code{"|"}) if different groups are present; the number of PFAS 
#' groups split off (\code{n_PFAS_Groups}), the SMARTS string used to do the splitting 
#' \code{SplitSMARTS}, the PNG location (\code{PNG_location}), which is empty if PNG creation is off
#' and finally \code{Error_msg}, which contains an error message if issues were encountered. 
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragTools)
#' and Zhanyun Wang (ETHZ, PFAS examples).
#' 
#' @seealso \code{\link{splitPFAS.help}}
#' 
#' @export
#'
#' @examples
#' splitPFAStools_name <- "MetFrag2.4.5-Tools.jar"
#' smarts_file <- system.file("extdata","SplitPFAS_smarts.txt",package="RChemMass")
#' pfas_smiles <- "CN(CCOC(=O)C(C)=C)S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
#' pfas_smiles <- "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)CS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
#' pfas_smiles <- "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOP(O)(OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
#' pfas_smiles <- "FCC(F)(F)C(F)C(F)(F)F"
#' pfas_smiles <- "CCC(CN)C(F)(F)C(F)(F)C(F)(F)F"
#' 
#' splitPFAS(pfas_smiles,splitPFAStools_name)
#' splitPFAS(pfas_smiles,splitPFAStools_name,smarts_file = smarts_file)
#' splitPFAS(pfas_smiles,splitPFAStools_name,smarts_file = smarts_file,incl_image = TRUE)
#' 
splitPFAS <- function(smiles, splitPFAStools_name, incl_image=FALSE, incl_smarts=TRUE,
                      smarts_file="", splitPFAStools_dir="") {
  # first, get set up
  curr_dir <- getwd()
  CommandPart1 <- paste0("java -cp ",splitPFAStools_name, " de.ipbhalle.metfrag.split.SplitPFAS")
  CommandPart1 <- paste0(CommandPart1," smiles=")
  CommandPart2 <- ""
  if (incl_smarts && nchar(smarts_file)>0) {
    CommandPart2 <- paste0(CommandPart2," smartspath=",smarts_file)
  } else if (incl_smarts) {
    smarts_file <- system.file("extdata","SplitPFAS_smarts.txt",package="RChemMass")
    CommandPart2 <- paste0(CommandPart2," smartspath=",smarts_file)
  }
  if (incl_image) {
    CommandPart2 <- paste0(CommandPart2," image=yes")
  }
  # build command
  SplitPFASCommand <- paste0(CommandPart1,smiles," ",CommandPart2)
  
  # run the command
  if (nchar(splitPFAStools_dir)>0) {
    setwd(splitPFAStools_dir)
  }
  #  CommandOut <- system(command=SplitPFASCommand,intern=TRUE,show.output.on.console=FALSE)
  CommandOut <- suppressWarnings(system(command=SplitPFASCommand,intern=TRUE,
                                        show.output.on.console=FALSE))
  setwd(curr_dir)
  
  # process the output
  #find index of SMILES
  smiles_index <- grep("Reading",CommandOut)
  #see if PNG exists:
  png_start_string <- "stored bond(s) to break highlighted in "
  png_index <- grep(png_start_string,CommandOut,fixed=TRUE)
  png_dir <- paste0(curr_dir,"/SplitPFAS_PNGs")
  if (!dir.exists(png_dir)) {
    dir.create(png_dir)
  }
  #length(png_index)=0 if it is not found
  #then find the results index
  results_index <- smiles_index+length(png_index)+1 #length png_index is either 0 or 1
  
  if (length(png_index)>0) {
    png_location <- sub(png_start_string,"",CommandOut[png_index],fixed=TRUE)
    #find file, copy it to current dir and delete from temp.
    file.exists(png_location)
    file.copy(png_location,png_dir)
    file.remove(png_location)
    png_location <- paste0(png_dir,"/",basename(png_location))
  } else {
    png_location <- ""
  }
  
  # Error messages to detect
  # "No splittable bond found for the input molecule FC(F)C(F)(F)C(F)(F)C(F)(F)F"
  # even if image=yes, no image is output. 
  # "Error: Problem occured. Check input molecule. More than one functional group left after split."
  # "Error: Problem occured. Check input."
  # "Error: The initial match for the PFAS chain ending seemed to cause an error. Check your input SMILES."
  error_msg_start_strings <- c("No splittable bond found for the input molecule ",
                               "Error: Problem occured. ",
                               "Error: The initial match for the ")
  error_index1 <- grep(error_msg_start_strings[1], CommandOut, fixed=TRUE)
  error_index2 <- grep(error_msg_start_strings[2], CommandOut, fixed=TRUE)
  error_index3 <- grep(error_msg_start_strings[3], CommandOut, fixed=TRUE)
  # if we have an error and then a successful result, need to check whether we 
  # have valid results ... 
  if (length(error_index1)>0) {
    error_msg <- CommandOut[error_index1]
    results_index <- results_index+1
  } else if (length(error_index2)>0) {
    error_msg <- CommandOut[error_index2]
    results_index <- results_index+1
  } else if (length(error_index3)>0) {
    error_msg <- CommandOut[error_index3]
    results_index <- results_index+1
  } else {
    error_msg <- ""
  }
  # if we have an error and then a successful result, need to check 
  
  smiles_in <- sub("Reading ","",CommandOut[smiles_index],fixed=TRUE)
  results_smiles <- CommandOut[results_index]
  # the results are four parts, func_smiles, PFAS_part, n_PFAS_parts and SMARTS
  if (length(results_smiles)>0) {
    #if (!is.na(results_smiles)) {
    results <- strsplit(results_smiles," ",fixed=TRUE)[[1]]
  } else {
    results <- c("","","","")
  }
  #results <- strsplit(results_smiles," ",fixed=TRUE)[[1]]
  
  # now, make output
  # test if there's an error, in which case output the error message and
  # empty results strings. Otherwise output results and empty error_msg string
  
  splitPFASterms <- list()
  splitPFASterms[['SMILES_in']] <- smiles_in
  
  # if (nchar(error_msg)>2) {
  #   splitPFASterms[['Func_Group']] <- ""
  #   splitPFASterms[['PFAS_Groups']] <- ""
  #   splitPFASterms[['n_PFAS_Groups']] <- ""
  #   splitPFASterms[['SplitSMARTS']] <- ""
  # } else {
  splitPFASterms[['Func_Group']] <- results[1]
  splitPFASterms[['PFAS_Groups']] <- results[2]
  splitPFASterms[['n_PFAS_Groups']] <- results[3]
  splitPFASterms[['SplitSMARTS']] <- gsub("'","",sub("smarts=","",results[4]),fixed=TRUE)
  # }
  splitPFASterms[['PNG_location']] <- png_location
  splitPFASterms[['Error_msg']] <- error_msg
  
  return(splitPFASterms)
}

#' Run SplitPFAS on a CSV file containing SMILES 
#' 
#' @description A little wrapper function to run SplitPFAS on a csv file
#' containing SMILES and other fields, optionally running InChIKey checks
#' NOTE: some options only work for very specific files. Use with caution. 
#'
#' @param pfas_file Name of the CSV file containing at least one SMILES column
#' @param splitPFAStools_name Name of the SplitPFAS jar. See \code{\link{SplitPFAS}}
#' @param smarts_file Name of the SMARTS file. See \code{\link{SplitPFAS}}
#' @param babel_dir Location of Open Babel dir. Needed for InChIKey checks. 
#' @param incl_image Whether to include image generation. See \code{\link{SplitPFAS}}
#' @param smiles_col_index Index of the column containing SMILES. Default 3. 
#' @param SMILES_KeyCheck Boolean to indicate if SMILES InChIKey check should be performed.
#' @param func_group_col_index Index of column containing functional group. For validation only. 
#' @param Func_Group_KeyCheck Default \code{FALSE}. Checks functional group InChIKeys if \code{TRUE} 
#' - only for test files that contain this information (advanced users only)  
#' @param pfas_part_col_index Index of column containing functional group. For validation only. 
#' @param PFAS_Part_KeyCheck  Default \code{FALSE}. Checks functional group InChIKeys if \code{TRUE} 
#' - only for test files that contain this information (advanced users only).
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @details Returns a csv file containing the SplitPFAS results. 
#' 
#' @seealso \code{\link{splitPFAS.help}}, \code{\link{splitPFAS}}
#' 
#' @export
#'
#' @examples
#' pfas_file <- system.file("extdata","PFAS_SplitSmiles_ExamplesExt.csv",package="RChemMass")
#' splitPFAStools_name <- "MetFrag2.4.5-Tools.jar"
#' smarts_file <- system.file("extdata","SplitPFAS_smarts.txt",package="RChemMass")
#' babel_dir <- "C:/Program Files/OpenBabel-2.4.1"
#' splitPFAS.csv(pfas_file, splitPFAStools_name, smarts_file, babel_dir,
#'               incl_image = TRUE, smiles_col_index = 3, SMILES_KeyCheck = TRUE,
#'               func_group_col_index = 5, Func_Group_KeyCheck = FALSE,
#'               pfas_part_col_index = 6, PFAS_Part_KeyCheck = FALSE)
#'               
#' splitPFAS.csv(pfas_file, splitPFAStools_name, smarts_file, babel_dir,
#'               incl_image = TRUE, smiles_col_index = 3, SMILES_KeyCheck = TRUE,
#'               func_group_col_index = 5, Func_Group_KeyCheck = TRUE, 
#'               pfas_part_col_index = 6, PFAS_Part_KeyCheck = TRUE)
#'               
#' pfas_file <- system.file("extdata","Fluroseedset_1to30.csv",package="RChemMass")
#' splitPFAS.csv(pfas_file, splitPFAStools_name, smarts_file, babel_dir,
#'               incl_image = TRUE, smiles_col_index = 5, SMILES_KeyCheck = TRUE,
#'               Func_Group_KeyCheck = FALSE, PFAS_Part_KeyCheck = FALSE)
#'               
splitPFAS.csv <- function(pfas_file, splitPFAStools_name, smarts_file, babel_dir,
                          incl_image = TRUE, smiles_col_index = 3, SMILES_KeyCheck = TRUE,
                          func_group_col_index = 5, Func_Group_KeyCheck = FALSE, 
                          pfas_part_col_index = 6, PFAS_Part_KeyCheck = FALSE) {
  #with InChIKey checks
  #pfas_file <- "PFAS_SplitSmiles_ExamplesExt.csv"
  pfas_info <- read.csv(pfas_file,stringsAsFactors = F)
  smiles_col_index <- smiles_col_index
  func_group_col_index <- func_group_col_index
  pfas_part_col_index <- pfas_part_col_index
  
  pfas_info$SMILES_in <- ""
  pfas_info$Func_Group_Out <- ""
  pfas_info$PFAS_Groups_Out <- ""
  pfas_info$N_PFAS_Groups_Out <- ""
  pfas_info$SplitSMARTS <- ""
  pfas_info$PNG_location <- ""
  pfas_info$SplitPFAS_errorMsg <- ""
  if (SMILES_KeyCheck) {
    pfas_info$SMILES_InChIKey <- ""
    pfas_info$SMILES_in_InChIKey <- ""
    pfas_info$SMILES_InChIKeyCheck <- ""
  }
  if (Func_Group_KeyCheck) {
    pfas_info$FUNC_GROUP_Key <- ""
    pfas_info$Func_Group_Out_Key <- ""
    pfas_info$FuncGroupKeyCheck <- ""
  }
  if (PFAS_Part_KeyCheck) {
    pfas_info$PFAS_GROUP_Key <- ""
    pfas_info$PFAS_Group_Out_Key <- ""
    pfas_info$PFASGroupKeyCheck <- ""
  }
  
  
  for (i in 1:length(pfas_info[,smiles_col_index])) {
    # get SMILES
    pfas_smiles <- pfas_info[i,smiles_col_index]
    SplitPFAS_output <- splitPFAS(pfas_smiles, splitPFAStools_name, incl_image = incl_image, 
                                  smarts_file = smarts_file)
    if (SMILES_KeyCheck) {
      SMILES_InChIKey <- getInChIKey.obabel(pfas_smiles, babel_dir)
      SMILES_in_InChIKey <- getInChIKey.obabel(SplitPFAS_output$SMILES_in, babel_dir)
    }
    if (Func_Group_KeyCheck) {
      FUNC_GROUP_Key <- getInChIKey.obabel(pfas_info$FUNC_GROUP[i], babel_dir)
      Func_Group_Out_Key <- getInChIKey.obabel(SplitPFAS_output$Func_Group, babel_dir)
    } 
    if (PFAS_Part_KeyCheck) {
      PFAS_GROUP_Key <- getInChIKey.obabel(pfas_info$PFAS_PART[i], babel_dir)
      PFAS_Group_Out_Key <- getInChIKey.obabel(SplitPFAS_output$PFAS_Groups, babel_dir)
    }
    
    
    #save the values in the output file
    pfas_info$SMILES_in[i] <- SplitPFAS_output$SMILES_in
    pfas_info$Func_Group_Out[i] <- SplitPFAS_output$Func_Group
    pfas_info$PFAS_Groups_Out[i] <- SplitPFAS_output$PFAS_Groups
    pfas_info$N_PFAS_Groups_Out[i] <- SplitPFAS_output$n_PFAS_Groups
    pfas_info$SplitSMARTS[i] <- SplitPFAS_output$SplitSMARTS
    pfas_info$PNG_location[i] <- SplitPFAS_output$PNG_location
    pfas_info$SplitPFAS_errorMsg[i] <- SplitPFAS_output$Error_msg
    if (SMILES_KeyCheck) {
      pfas_info$SMILES_InChIKey[i] <- SMILES_InChIKey
      pfas_info$SMILES_in_InChIKey[i] <- SMILES_in_InChIKey
      pfas_info$SMILES_InChIKeyCheck[i] <- (SMILES_InChIKey == SMILES_in_InChIKey)
    }
    if (Func_Group_KeyCheck) {
      pfas_info$FUNC_GROUP_Key[i] <- FUNC_GROUP_Key
      pfas_info$Func_Group_Out_Key[i] <- Func_Group_Out_Key
      pfas_info$FuncGroupKeyCheck[i] <- (FUNC_GROUP_Key==Func_Group_Out_Key)
    }
    if (PFAS_Part_KeyCheck) {
      pfas_info$PFAS_GROUP_Key[i] <- PFAS_GROUP_Key
      pfas_info$PFAS_Group_Out_Key[i] <- PFAS_Group_Out_Key
      pfas_info$PFASGroupKeyCheck[i] <- (PFAS_GROUP_Key == PFAS_Group_Out_Key) 
    }
  }
  curr_dir <- getwd()
  new_pfas_csv <- paste0(curr_dir,"/", sub(".csv","_wResultsChecks.csv",basename(pfas_file),fixed=TRUE))
  write.csv(pfas_info,new_pfas_csv,row.names=F)
  return_msg <- paste0("Output written in ", new_pfas_csv)
  return(return_msg)  
}
