# Chemical curation in R
# Collection of various scripts written over the years
# Taking advantage of several different packages
# E. Schymanski, 15/2/2017

# # package installation and dependencies (note rinchi is only on github)
# library(devtools)
# # install_github("cdkr", "CDK-R", subdir="rinchi")
# # install_github("CDK-R/cdkr", subdir="rinchi")
# library(rinchi)
# library(RMassBank)
# library(rcdk)
# library(enviPat)
m_electron <- 0.0005485799090


#' Extract Molecular Formula from InChI
#'
#' @description This basic function using string manipulation to extract a formula from
#' the front of an InChI string, returning the second field (the formula).
#' Note: this does NOT account for labeling, i.e. it is the formula associated
#' with the structural skeleton (the InChIKey first block).
#'
#' @usage MolFormFromInChI(InChI)
#'
#' @param InChI InChI string to extract the formula.
#'
#' @return Returns the molecular formula as is in the InChI string.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{MolFormFromSmiles.rcdk}}.
#'
#' @export
#'
#' @examples
#' inchi <- "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
#' MolFormFromInChI(inchi)
#'
MolFormFromInChI <- function(InChI) {
  molform <- strsplit(InChI,"/")[[1]][2]
  return(molform)
}


#' Calculate Molecular Formula from SMILES with rcdk
#'
#' @description A small wrapper function to calculate the molecular formula
#' from SMILES with the rcdk, as recommended. Note: current version does not handle labelling.
#'
#' @usage MolFormFromSmiles.rcdk(smiles)
#'
#' @param smiles Valid SMILES code used to calculate the formula.
#'
#' @return Returns the molecular formula only (the \code{rcdk} function contains additional text).
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{parse.smiles}}, \code{\link{get.mol2formula}}, \code{\link{MolFormFromInChI}}.
#'
#' @export
#'
#' @examples
#' MolFormFromSmiles.rcdk("OC(=O)C")
#'
MolFormFromSmiles.rcdk <- function(smiles) {
  mol <- parse.smiles(smiles)[[1]]
  convert.implicit.to.explicit(mol)
  formula <- get.mol2formula(mol, charge=0)
  return(formula@string)
}

#' Calculate Molecular Formula and Exact Masses from SMILES
#'
#' @description A small wrapper function to calculate the molecular formula
#' from SMILES with the rcdk, as well as selected MS-relevant exact masses for fixed adducts
#' using enviPat. This uses \code{getMolecule} instead of \code{parse.smiles}.
#' Note: current version does not handle labelling as rcdk is used to obtain the formula.
#' See \code{\link{getSuspectMasses}} for an extension handling all enviPat adducts.
#'
#' @usage getSuspectFormulaMass(smiles)
#'
#' @param smiles Valid SMILES code used to calculate the formula and from then on the masses.
#'
#' @return Returns the molecular formula, monoisotopic mass, [M+H]+, [M+NH4]+, [M+Na]+ and [M-H]-
#' masses in a named list. If no H is present for [M-H]-, NA is returned.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getMolecule}}, \code{\link{get.mol2formula}}, \code{\link{enviPat}},
#' \code{\link{getSuspectMasses}}.
#'
#' @export
#'
#' @examples
#' getSuspectFormulaMass("c1ccccc1")
#'
#'
getSuspectFormulaMass <- function(smiles) {
  mol <- try(getMolecule(smiles), silent=T)
  if (is.atomic(mol)) {
    mol_form <- ""
    Monoiso_mass <- ""
    MpHp_mass <- ""
    MpNH4_mass <- ""
    MpNa_mass <- ""
    MmHm_mass <- ""
    print(paste0("Invalid smiles: ", smiles))
  } else {
    mol_form <- strsplit(capture.output(get.mol2formula(mol, charge=0))," ")[[1]][3]
    mol_form_checked <- check_chemform(isotopes,mol_form)
    mol_form <- mol_form_checked$new_formula
    # calculate masses
    Monoiso_mass <- isopattern(isotopes, mol_form, verbose=FALSE)[[1]][1]
    MpHp_form <- mergeform(mol_form,as.character(adducts$Formula_add[1]))
    MpHp_mass <- isopattern(isotopes, MpHp_form, charge=adducts$Charge[1],verbose=FALSE)[[1]][1]
    MpNH4_form <- mergeform(mol_form,as.character(adducts$Formula_add[2]))
    MpNH4_mass <- isopattern(isotopes, MpNH4_form, charge=adducts$Charge[2],verbose=FALSE)[[1]][1]
    MpNa_form <- mergeform(mol_form,as.character(adducts$Formula_add[3]))
    MpNa_mass <- isopattern(isotopes, MpNa_form, charge=adducts$Charge[3],verbose=FALSE)[[1]][1]
    if (length(grep("H",mol_form))>0) {
      MmHm_form <- subform(mol_form,as.character(adducts$Formula_ded[6]))
      MmHm_mass <- isopattern(isotopes, MmHm_form, charge=adducts$Charge[6],verbose=FALSE)[[1]][1]
    } else {
      MmHm_mass <- NA
    }
    #}
    #  MmHm_form <- subform(mol_form,as.character(adducts$Formula_ded[6]))
    #  MmHm_mass <- isopattern(isotopes, MmHm_form, charge=adducts$Charge[6],verbose=FALSE)[[1]][1]
    #  neutral_monoiso_mass <- isopattern(isotopes, mol_form, charge=0,verbose=FALSE)[[1]][1]
  }
  SusFormMass <- list()
  SusFormMass[['MolForm']] <- mol_form
  SusFormMass[['Monoiso_mass']] <- Monoiso_mass
  SusFormMass[['MpHp_mass']] <- MpHp_mass
  SusFormMass[['MpNH4_mass']] <- MpNH4_mass
  SusFormMass[['MpNa_mass']] <- MpNa_mass
  SusFormMass[['MmHm_mass']] <- MmHm_mass
  #  }
  return(SusFormMass)
}


#' Calculate Exact Masses of Fixed Adducts from Molecular formula
#'
#' @description A small wrapper function to calculate selected MS-relevant exact masses
#' for fixed adducts from a molecular formula
#' using enviPat.
#'
#' @usage getAdductMassesFromFormula(MolForm)
#'
#' @param MolForm Molecular formula used to calculate the adduct masses.
#'
#' @return Returns the monoisotopic mass, [M+H]+, [M+NH4]+, [M+Na]+, [M-H]- and [M+FA-H]-
#' masses in a named list. If no H is present for [M-H]-, NA is returned.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectFormulaMass}}, \code{\link{enviPat}},
#' \code{\link{getSuspectMasses}}.
#'
#' @export
#'
#' @examples
#' getAdductMassesFromFormula("C12H12")
#' getAdductMassesFromFormula("C8H3D3Cl2O3")
#' getAdductMassesFromFormula("C9H6D3ClO3")
#' getAdductMassesFromFormula("C6H6")
#'
getAdductMassesFromFormula <- function(MolForm) {
  mol_form <- MolForm
  mol_form_checked <- check_chemform(isotopes,mol_form)
  mol_form <- mol_form_checked$new_formula
  if (mol_form_checked$warning) {
    mol_form <- MolForm
    Monoiso_mass <- ""
    MpHp_mass <- ""
    MpNH4_mass <- ""
    MpNa_mass <- ""
    MmHm_mass <- ""
    print(paste0("Invalid formula: ", MolForm))
  } else {
    # calculate masses
    Monoiso_mass <- isopattern(isotopes, mol_form, verbose=FALSE)[[1]][1]
    MpHp_form <- mergeform(mol_form,as.character(adducts$Formula_add[1]))
    MpHp_mass <- isopattern(isotopes, MpHp_form, charge=adducts$Charge[1],verbose=FALSE)[[1]][1]
    MpNH4_form <- mergeform(mol_form,as.character(adducts$Formula_add[2]))
    MpNH4_mass <- isopattern(isotopes, MpNH4_form, charge=adducts$Charge[2],verbose=FALSE)[[1]][1]
    MpNa_form <- mergeform(mol_form,as.character(adducts$Formula_add[3]))
    MpNa_mass <- isopattern(isotopes, MpNa_form, charge=adducts$Charge[3],verbose=FALSE)[[1]][1]
    MpFAm_form <- mergeform(mol_form,as.character(adducts$Formula_add[9]))
    MpFAm_form <- subform(MpFAm_form,as.character(adducts$Formula_ded[9]))
    MpFAm_mass <- isopattern(isotopes, MpFAm_form, charge=adducts$Charge[9],verbose=FALSE)[[1]][1]
    if (length(grep("H",mol_form))>0) {
      MmHm_form <- subform(mol_form,as.character(adducts$Formula_ded[6]))
      MmHm_mass <- isopattern(isotopes, MmHm_form, charge=adducts$Charge[6],verbose=FALSE)[[1]][1]
    } else {
      MmHm_mass <- NA
    }
  }

  SusFormMass <- list()
  SusFormMass[['MolForm']] <- mol_form
  SusFormMass[['Monoiso_mass']] <- Monoiso_mass
  SusFormMass[['MpHp_mass']] <- MpHp_mass
  SusFormMass[['MpNH4_mass']] <- MpNH4_mass
  SusFormMass[['MpNa_mass']] <- MpNa_mass
  SusFormMass[['MmHm_mass']] <- MmHm_mass
  SusFormMass[["MpFAm_mass"]] <- MpFAm_mass
  #  }
  return(SusFormMass)
}



#' Calculate Masses for various adducts from SMILES
#'
#' @description A small wrapper function to calculate adduct masses for any entries in
#' the \code{enviPat} \code{adducts} table from the SMILES code, using
#' the \code{getMolecule} function in \code{rcdk} to calculate the formula from SMILES first.
#' Note: current version does not handle labelling as rcdk is used to obtain the formula.
#'
#' @usage getSuspectMasses(smiles, adduct_list)
#'
#' @param smiles Valid SMILES code used to calculate the formula and from then on the masses.
#' @param adduct_list A list of adducts matching those in the \code{enviPat} \code{adducts} table.
#'
#' @return Returns the masses of the given adducts in a matrix with column names matching the given adduct.
#' If the adduct/deduct cannot be calculated for a given entry, NA is returned.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getMolecule}}, \code{\link{get.mol2formula}}, \code{\link{enviPat}},
#' \code{\link{check_ded2}}, \code{\link{getSuspectFormulaMass}}.
#'
#' @export
#'
#' @examples
#' getSuspectMasses("c1ccccc1", c("M+H","M+NH4"))
#' getSuspectMasses("c1ccccc1", c("M-H","M-H2O-H"))
#' getSuspectMasses("c1ccccc1", c("M+H","M+2H","M+3Na"))
#'
#'
getSuspectMasses <- function(smiles, adduct_list) {
  mol <- try(getMolecule(smiles), silent=T)
  adduct_list <- adduct_list
  n_adducts <- length(adduct_list)
  SusMasses <- matrix(ncol=n_adducts)
  colnames(SusMasses) <- adduct_list
  adduct <- ""
  adduct_formula <- ""
  to_add <- ""
  to_ded <- ""
  if (is.atomic(mol)) {
    stop(paste0("Invalid smiles: ", smiles))
  } else {
    mol_form <- strsplit(capture.output(get.mol2formula(mol, charge=0))," ")[[1]][3]
    mol_form_checked <- check_chemform(isotopes,mol_form)
    mol_form <- mol_form_checked$new_formula
    # calculate masses
    for (i in 1:length(adduct_list)) {
      #Monoiso_mass <- isopattern(isotopes, mol_form, verbose=FALSE)[[1]][1]
      adduct_entry <- which(adducts$Name == adduct_list[i])
      # restart from mol_form for new adduct
      adduct_formula <- mol_form
      #adduct_entry <- grep(adduct_list[i],adducts$Name, fixed=TRUE)[1]
      if (length(adduct_entry)== 0) {
        to_add <- FALSE
        to_ded <- FALSE
        adduct_formula <- NA
      } else {
        # check if an adduct, hence add:
        to_add <- as.character(adducts$Formula_add[adduct_entry])
        if (!(to_add == FALSE)) {
        adduct_formula <- mergeform(adduct_formula,as.character(adducts$Formula_add[adduct_entry]))
        }
        # check if a deduct entry, hence subtract
        to_ded <- as.character(adducts$Formula_ded[adduct_entry])
        if (!(to_ded == FALSE)) {
          # check that we can subtract the formula
          ded_OK <- check_ded2(adduct_formula,as.character(adducts$Formula_ded[adduct_entry]))
          if (ded_OK) {
            adduct_formula <- subform(adduct_formula, as.character(adducts$Formula_ded[adduct_entry]))
          } else {
            warning(paste0("Cannot subtract ", as.character(adducts$Formula_ded[adduct_entry]), " from formula ",
                           adduct_formula, " smiles ", smiles))
            adduct_formula <- NA
          }
        }
      }
      # now we need to calculate the mass for the formula
      if (!is.na(adduct_formula)) {
        SusMasses[1,i] <- isopattern(isotopes, adduct_formula, charge=adducts$Charge[adduct_entry],
                                     verbose=FALSE)[[1]][1]
      }

    }
  }
  return(SusMasses)
}

#' Check if "deduct" can be subtracted from given molecular formula
#'
#' @description This function uses the formula functionality in \code{RMassBank} to check
#' if a formula can be subtracted from another (i.e. is a subset of another), especially for
#' calculating adduct masses. This serves in the place of \code{check_ded} from \code{enviPat},
#' which returns the opposite result.
#'
#' @usage check_ded2(formula, deduct)
#'
#' @param formula Molecular formula to check.
#' @param deduct Deduct (molecular formula to subtract from \code{formula}).
#'
#' @return Returns \code{TRUE} if the \code{deduct} can be subtracted from \code{formula}, otherwise
#' \code{FALSE}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{check_ded}}, \code{\link{enviPat}}, \code{\link{getSuspectFormulaMass}},
#' \code{\link{formulastring.to.list}}.
#'
#' @export
#'
#' @examples
#' check_ded2("C6H6", "H2")
#' check_ded2("C6H6", "H3O")
#'
check_ded2 <- function(formula, deduct) {
  formula <- formula
  deduct <- deduct
  formula_list <- formulastring.to.list(formula)
  deduct_list <- formulastring.to.list(deduct)
  ded_elements <- names(deduct_list)
  formula_elements <- names(formula_list)
  el_in_formula <- vector(length=length(ded_elements))
  n_element_in_formula <- 0
  n_element_in_ded <- 0
  for (i in 1:length(ded_elements)) {
    element <- ded_elements[i]
    n_element_in_ded <- deduct_list[i]
    ind_el <- grep(element,formula_elements)
    if (length(ind_el)==0) {
      ind_el <- 0
      n_element_in_formula <- 0
    } else {
      n_element_in_formula <- as.numeric(formula_list[ind_el])
      if (n_element_in_formula >= n_element_in_ded) {
        el_in_formula[i] <- TRUE
      }
    }
  }
  # now test if deduct is completely within formula
  ded_OK <- FALSE
  if (length(grep("FALSE",el_in_formula,fixed=TRUE))==0) {
    ded_OK <- TRUE
  }
  return(ded_OK)
}

#' Retrieve an InChIKey by SMILES or Name with Cactus/CTS
#'
#' @description A small wrapper function to retrieve an InChIKey with webservices. Cactus
#' (\url{https://cactus.nci.nih.gov/chemical/structure}) is used to obtain a Key
#' from SMILES; if this fails the Chemical Translation Service (CTS,
#' \url{http://cts.fiehnlab.ucdavis.edu/} is queried by name.
#' Future extensions possible. If SMILES is available and Open Babel is installed,
#' use \code{\link{getSuspectInChIKey.babel}} instead.
#'
#' @usage getSuspectInChIKey(smiles,name)
#'
#' @param smiles Valid SMILES code used to retrieve the InChIKey.
#' @param name Chemical name used to query CTS. Strange names will likely fail.
#'
#' @return Returns the InChIKey retrieved, or \code{NA}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectInChIKey.babel}}.
#'
#' @export
#'
#' @examples
#' getSuspectInChIKey("c1ccccc1", "benzene")
#'
getSuspectInChIKey <- function(smiles,name) {
  InChIKey <- sub("InChIKey=","",getCactus(smiles,"stdinchikey")[1],fixed=TRUE)
  if (nchar(InChIKey)!=27) {
    InChIKey2 <- getCtsKey(name)
    if (nchar(InChIKey2)!=27) {
      InChIKey <- NA
    } else {
      InChIKey <- InChIKey2
    }
  }
  # if this doesn't work, use babel - in separate function
  return(InChIKey)
}


#' Convert SMILES to an InChIKey with OpenBabel (obabel)
#'
#' @description A small wrapper function to convert SMILES to an InChIKey with OpenBabel
#' (\url{http://openbabel.org/wiki/Main_Page}). Requires pre-installation of
#' OpenBabel. If this is not the case, use webservices instead via
#' \code{\link{getSuspectInChIKey}}. Note this supercedes \code{\link{getSuspectInChIKey.babel}}
#' This function uses Babel default InChI options; standard InChIKeys will be generated.
#'
#' @usage getInChIKey.obabel(smiles,babel_dir)
#'
#' @param smiles SMILES code to convert to the InChIKey.
#' @param babel_dir Location of folder containing \code{"obabel.exe"}.
#'
#' @return Returns the InChIKey retrieved, or alternative output from Babel. Run
#' \code{\link{InChIKey_test}} to determine if valid.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectInChIKey}}, \code{\link{InChIKey_test}}, \code{\link{getSuspectInChIKey.babel}}.
#'
#' @export
#'
#' @examples
#' babel_dir <- "C:/Program Files (x86)/OpenBabel-2.3.2"
#' getInChIKey.obabel("c1ccccc1", babel_dir)
#' Various failed conversions:
#' getInChIKey.obabel("blah", babel_dir)
#' getInChIKey.obabel("", babel_dir)
#' InChIKey_test(getInChIKey.obabel("blah",babel_dir))
#'
getInChIKey.obabel <- function(smiles,babel_dir) {
  dir <- getwd()
  if (file.exists(babel_dir)) {
    setwd(babel_dir)
    Babel_out <- system("obabel -ismi -oinchikey",input=smiles,show.output.on.console = F,intern=T)
  } else {
    Babel_out = ""
    warning("Incorrect Babel directory or directory does not exist")
  }
  setwd(dir)
  return(Babel_out[1])
}



#' Convert SMILES to an InChIKey with OpenBabel (superceded)
#'
#' @description A small wrapper function to convert SMILES to an InChIKey with OpenBabel
#' (\url{http://openbabel.org/wiki/Main_Page}). Requires pre-installation of
#' OpenBabel. If this is not the case, use webservices instead via
#' \code{\link{getSuspectInChIKey}}. Note this creates \code{"temp_smiles.smi"}
#' in the \code{temp_dir} directory and will overwrite any previous file of the same name.
#' This function uses Babel default InChI options; standard InChIKeys will be generated.
#' This is superceded by \code{\link{getInChIKey.obabel}}.
#'
#' @usage getSuspectInChIKey.babel(smiles,babel_dir,temp_dir)
#'
#' @param smiles Valid SMILES code used to convert to the InChIKey.
#' @param babel_dir Location of folder containing \code{"babel.exe"}.
#' @param temp_dir Location of folder to save \code{"temp_smiles.smi"}.
#'
#' @return Returns the InChIKey retrieved, or alternative output from Babel. Run
#' \code{\link{InChIKey_test}} to determine if valid.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectInChIKey}}, \code{\link{getInChIKey.obabel}}, \code{\link{InChIKey_test}}.
#'
#' @export
#'
#' @examples
#' getSuspectInChIKey.babel("c1ccccc1", "C:/Program Files (x86)/OpenBabel-2.3.2", "C:/DATA")
#' Various failed conversions:
#' getSuspectInChIKey.babel("blah", "C:/Program Files (x86)/OpenBabel-2.3.2", "C:/DATA")
#' getSuspectInChIKey.babel("", "C:/Program Files (x86)/OpenBabel-2.3.2", "C:/DATA")
#'
getSuspectInChIKey.babel <- function(smiles,babel_dir,temp_dir) {
  dir <- getwd()
  if (file.exists(babel_dir)) {
    setwd(babel_dir)
    temp_file_loc <- paste0(temp_dir,"/temp_smiles.smi")
    file.conn <- file(temp_file_loc)
    open(file.conn,open="wt")
    writeLines(smiles, con=file.conn)
    close(file.conn)
    BabelCommand <- paste("babel -ismi ", temp_file_loc, " -oinchikey",sep="")
    Babel_out <- system(command=BabelCommand,intern=TRUE,show.output.on.console=FALSE)
  } else {
    Babel_out = ""
    warning("Incorrect Babel directory or directory does not exist")
  }
  setwd(dir)
  return(Babel_out[1])
}


#' Convert SMILES to an InChI with OpenBabel (obabel)
#'
#' @description A small wrapper function to convert SMILES to an InChI with OpenBabel
#' (\url{http://openbabel.org/wiki/Main_Page}). Requires pre-installation of
#' OpenBabel. 
#' This function uses Babel default InChI options; standard InChIs will be generated.
#'
#' @usage getInChI.obabel(smiles,babel_dir)
#'
#' @param smiles SMILES code to convert to the InChI.
#' @param babel_dir Location of folder containing \code{"obabel.exe"}.
#'
#' @return Returns the InChI retrieved, or alternative output from Babel. 
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getInChI.obabel}}.
#'
#' @export
#'
#' @examples
#' babel_dir <- "C:/Program Files (x86)/OpenBabel-2.3.2"
#' getInChI.obabel("c1ccccc1", babel_dir)
#' Various failed conversions:
#' getInChI.obabel("blah", babel_dir)
#' getInChI.obabel("", babel_dir)
#'
getInChI.obabel <- function(smiles,babel_dir) {
  dir <- getwd()
  if (file.exists(babel_dir)) {
    setwd(babel_dir)
    Babel_out <- system("obabel -ismi -oinchi",input=smiles,show.output.on.console = F,intern=T)
  } else {
    Babel_out = ""
    warning("Incorrect Babel directory or directory does not exist")
  }
  setwd(dir)
  return(Babel_out[1])
}




#' Convert an InChI to SMILES with OpenBabel
#'
#' @description A small wrapper function to convert an InChI to SMILES with OpenBabel
#' (\url{http://openbabel.org/wiki/Main_Page}). Requires pre-installation of
#' OpenBabel. If this is not the case, use webservices instead via
#' \code{\link{getCactus}}.  Note this creates \code{"temp_inchi.inchi"}
#' in the \code{temp_dir} directory and will overwrite any previous file of the same name.
#'
#' @usage getSmilesFromInChI.babel(inchi,babel_dir,temp_dir)
#'
#' @param smiles Valid InChI code used to convert to SMILES. Default babel settings are used.
#' @param babel_dir Location of folder containing \code{"babel.exe"}.
#' @param temp_dir Location of folder to save \code{"temp_inchi.inchi"}.
#'
#' @return Returns the SMILES retrieved, or alternative output from Babel.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getCactus}}, \code{\link{getSuspectInChIKey.babel}}.
#'
#' @export
#'
#' @examples
#' inchi <- "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
#' getSmilesFromInChI.babel(inchi, "C:/Program Files (x86)/OpenBabel-2.3.2", "C:/DATA")
#'
getSmilesFromInChI.babel <- function(inchi,babel_dir, temp_dir) {
  dir <- getwd()
  if (file.exists(babel_dir)) {
    setwd(babel_dir)
    temp_file_loc <- paste0(temp_dir,"/temp_smiles.smi")
    file.conn <- file(temp_file_loc)
    open(file.conn,open="wt")
    writeLines(inchi, con=file.conn)
    close(file.conn)
    BabelCommand <- paste("babel -iinchi C:/DATA/temp_inchi.inchi -osmi",sep="")
    Babel_out <- system(command=BabelCommand,intern=TRUE,show.output.on.console=FALSE)
    Babel_out <- substring(Babel_out,1,(nchar(Babel_out)-1))
  } else {
    Babel_out = ""
    warning("Incorrect Babel directory or directory does not exist")
  }
  setwd(dir)
  return(Babel_out[1])
}

#' Get PubChem CID, ChemSpider ID and CAS from InChIKey via webservices
#'
#' @description A wrapper function to obtain identifiers from an InChIKey
#' using various webservices. PubChem CIDs are obtained with \code{\link{getPcId}} from
#' PubChem direct (\url{https://pubchem.ncbi.nlm.nih.gov/}), ChemSpider IDs
#' directly from ChemSpider (\url{http://www.chemspider.com/}) with
#' \code{\link{getCSID}}. CAS numbers are retrieved using \code{\link{getCactus}} from
#' Cactus (\url{https://cactus.nci.nih.gov/chemical/structure});
#' if this fails the Chemical Translation Service (CTS,
#' \url{http://cts.fiehnlab.ucdavis.edu/} is queried. Options are given to
#' return all CAS numbers, the "smallest" CAS number or the first entry.
#'
#' @usage getSuspectIdentifiers(InChIKey, allCAS=FALSE, sep=",", minCAS=TRUE)
#'
#' @param InChIKey Valid InChIKey used to retrieve information.
#' @param allCAS Default \code{FALSE} will result in only one CAS number.
#' If \code{TRUE}, all CAS are returned in one string, separated by \code{sep}.
#' This overrides \code{minCAS}.
#' @param sep Determines the separator used with \code{allCAS}, default \code{","}.
#' @param minCAS Default \code{TRUE} returns the smallest CAS number (shortest string) if
#' \code{allCAS=FALSE}. If \code{FALSE}, returns the first CAS number in the list.
#'
#' @return Returns a list with the resulting identifiers, with \code{NA} if the query failed.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getCactus}}, \code{\link{getPcId}}, \code{\link{getCSID}}
#' , \code{\link{getCtsRecord}}, \code{\link{getSuspectInChIKey}}.
#'
#' @export
#'
#' @examples
#' benz_key <- "UHOVQNZJYSORNB-UHFFFAOYSA-N"
#' getSuspectIdentifiers(benz_key)
#' getSuspectIdentifiers(benz_key,allCAS=TRUE)
#' getSuspectIdentifiers(benz_key,allCAS=FALSE,minCAS=TRUE)
#' getSuspectIdentifiers(benz_key,allCAS=FALSE,minCAS=FALSE)
#'
getSuspectIdentifiers <- function(InChIKey, allCAS=FALSE, sep=",", minCAS=TRUE) {
  PcId <- getPcId(InChIKey)
  CsId <- as.numeric(getCSID(InChIKey))
  CTS_info <- getCtsRecord(InChIKey)
  cas <- getCactus(InChIKey,"cas")
  if (!is.na(cas) && allCAS) {
    cas <- paste0(cas,collapse=sep)
  } else if (!is.na(cas) && !allCAS && minCAS) {
    cas <- cas[[which.min(nchar(cas))]] #take the smallest entry
  } else if (!is.na(cas) && !allCAS && !minCAS) {
    cas <- cas[1] #take the first entry
  }
  # if CAS was not retrieved from CACTUS, try CTS
  if ((CTS_info[1] == "Sorry, we couldn't find any matching results") || is.null(CTS_info[1])) {
    CTS_info <- NA
    #cas <- "NA"
  }
  # if we have something from CTS, get the CAS number. Otherwise stays as NA
  if(!is.na(CTS_info[1]) && !is.na(nchar(cas[1])<=2)){
    if("CAS" %in% CTS.externalIdTypes(CTS_info)) {
      # If multiple CAS, take the shortest one.
      cas <- CTS.externalIdSubset(CTS_info,"CAS")
      if (allCAS) {
        cas <- paste0(cas,collapse=sep)
      } else if (!allCAS && minCAS) {
        cas <- cas[[which.min(nchar(cas))]] #take the smallest entry
      } else if (!allCAS && !minCAS) {
        cas <- cas[1] #take the first entry
      }
    }
  }
  #no cross-linking to StoffIdent yet.
  SusIDs <- list()
  SusIDs[['PCID']] <- PcId
  SusIDs[['CSID']] <- CsId
  SusIDs[['CAS']] <- cas
  return(SusIDs)
}


#' Use rcdk to calculate AlogP and XlogP from SMILES
#'
#' @description A small wrapper function to calculate AlogP and XlogP
#' with rcdk from a SMILES code.
#'
#' @usage getCDKlogPs(smiles, kekulise=TRUE)
#'
#' @param smiles Valid SMILES code used to perform calculation.
#' @param kekulise Default \code{TRUE} performs aromaticity detection prior to calculation.
#' \code{FALSE} interprets the SMILES "as is". Can result in dramatically different logP values.
#'
#' @return Returns the resulting AlogP and XlogP values in a list.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{parse.smiles}}, \code{\link{get.alogp}}, \code{\link{get.xlogp}},
#' \code{\link{getCDKlogPsAndRT}}.
#'
#' @export
#'
#' @examples
#' getCDKlogPs("c1ccccc1")
#' getCDKlogPs("c1ccccc1",kekulise=FALSE)
#' getCDKlogPs("c1ccc2c(c1)[nH]nn2")
#' getCDKlogPs("c1ccc2c(c1)[nH]nn2",kekulise=FALSE)
#'
getCDKlogPs <- function(smiles, kekulise=TRUE) {
  smiles <- as.character(smiles)
  e <- simpleError("Invalid SMILES: try other kekulise setting or creating a unified canonical SMILES and retry")
  mol <- tryCatch(parse.smiles(smiles, kekulise), error= function(e) e)
  if (length(mol)>1) {
    stop(e)
  } else {
    mol <- mol[[1]]
  }
  # calculate values
  if (!is.null(mol)) {
    alogp_suspect <- get.alogp(mol)
    xlogp_suspect <- get.xlogp(mol)
  } else {
    alogp_suspect <- NA
    xlogp_suspect <- NA
  }
  logPs <- list()
  logPs[['CDKXlogP']] <- xlogp_suspect
  logPs[['CDKAlogP']] <- alogp_suspect
  return(logPs)
}


#' Use rcdk to calculate AlogP and XlogP from SMILES and estimate retention times
#'
#' @description A small wrapper function to calculate AlogP and XlogP
#' with rcdk from a SMILES code (as \code{\link{getCDKlogPs}}) plus estimate retention
#' times using a simple linear model. Default model was calculated on 810 Eawag reference
#' standards and is only valid for this Xbridge chromatographic program.
#'
#' @usage getCDKlogPsAndRT <- function(smiles,kekulise=TRUE,model=c(3.654,1.6225),quiet=TRUE)
#'
#' @param smiles Valid SMILES code used to perform calculation.
#' @param kekulise Default \code{TRUE} performs aromaticity detection prior to calculation.
#' \code{FALSE} interprets the SMILES "as is". Can result in dramatically different logP values.
#' @param model A vector of length 2 containing the intercept and slope for the model. Default values
#' \code{c(3.654,1.6225)} arise from the standard Eawag non-target chromatographic program on 810 standards.
#'
#' @return Returns the resulting AlogP, XlogP and RT estimated from XlogP values in a list.
#' Default model returns RT in minutes.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{parse.smiles}}, \code{\link{get.alogp}}, \code{\link{get.xlogp}},
#' \code{\link{getCDKlogPs}}.
#'
#' @export
#'
#' @examples
#' getCDKlogPsAndRT("c1ccccc1")
#' getCDKlogPsAndRT("c1ccccc1",kekulise=FALSE)
#' getCDKlogPsAndRT("c1ccc2c(c1)[nH]nn2")
#' getCDKlogPsAndRT("c1ccc2c(c1)[nH]nn2",kekulise=FALSE)
#'
getCDKlogPsAndRT <- function(smiles,kekulise=TRUE,model=c(3.654,1.6225),quiet=TRUE) {
  smiles <- as.character(smiles)
  e <- simpleError("Invalid SMILES: try other kekulise setting or creating a unified canonical SMILES and retry")
  mol <- tryCatch(parse.smiles(smiles, kekulise), error= function(e) e)
  if (length(mol)>1) {
    stop(e)
  } else {
    mol <- mol[[1]]
  }
  # calculate values
  alogp_suspect <- get.alogp(mol)
  xlogp_suspect <- get.xlogp(mol)
  RT_fromXlogP <- round(model[1]+model[2]*xlogp_suspect,digits = 4)

  logPs <- list()
  logPs[['CDKXlogP']] <- round(xlogp_suspect,digits=4)
  logPs[['CDKAlogP']] <- round(alogp_suspect,digits=4)
  logPs[['RT_fromXlogP']] <- RT_fromXlogP
  if (!quiet) {
    print("Default model: RTs (minutes) calulated on dataset of 810 Eawag standards using CDKXlogP; RSE=2.72")
  }
  return(logPs)
}

#' Test whether CAS number follows valid format
#'
#' @description This basic function using string manipulation to test whether a CAS number
#' is in the expected format. Useful to detect errors that occur in excel etc.
#' Potential future extension: to return corrected CAS numbers (not yet implemented).
#' This function checks if the string is longer than 7 and contains 3 dashes.
#'
#' @usage CAS_test(CAS_RN)
#'
#' @param CAS_RN Text string containing the CAS number to test.
#'
#' @return Returns \code{TRUE} or \code{FALSE} (fails test).
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectIdentifiers}}.
#'
#' @export
#'
#' @examples
#' CAS_test("95-14-7")
#' CAS_test("01.03.7440")
#'
CAS_test <- function(CAS_RN) {
  cas_test <- FALSE
  if (nchar(CAS_RN)>=7 && length(strsplit(CAS_RN,"-")[[1]]) == 3) {
    cas_test <- TRUE
  }
  return(cas_test)
}


#' Fix CAS numbers if they fail CAS_test
#'
#' @description This function uses basic string manipulation to split
#' CAS numbers that have been converted to dates in excel. The CAS number
#' should have failed \code{\link{CAS_test}} prior to fixing.
#'
#' @usage CAS_fix(CAS_RN)
#'
#' @param CAS_RN Text string containing the CAS number to fix.
#'
#' @return Returns the fixed CAS number.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{CAS_test}}.
#'
#' @export
#'
#' @examples
#' CAS_fix("95-14-7")
#' CAS_fix("01.03.7440")
#' CAS_fix("1/12/5392")
#' CAS_fix("NotCASNumber")
#'
CAS_fix <- function(CAS_RN) {
  if (CAS_test(CAS_RN)) {
    fixed_cas <- CAS_RN
  } else {
    split_cas <- strsplit(CAS_RN,"[[:punct:]]")[[1]]
    #fixed_cas <- ""
    n <- length(split_cas)
    if (n == 3) {
      fixed_cas <- ""
      for (i in seq_along(split_cas)) {
        if (i == 1) {
          fixed_cas <- split_cas[n-(i-1)]
        } else {
          fixed_cas <- paste(fixed_cas, split_cas[n-(i-1)],sep="-")
        }
      }
    } else {
      fixed_cas <- CAS_RN
      warning(paste0("CAS number does not fit pattern; returning input: ", CAS_RN))
    }
  }
  return(fixed_cas)
}


#' Trim text from start of Cactus InChIKey
#'
#' @description This basic function trims \code{"InChIKey="} from the start of
#' InChIKeys obtained via Cactus with \code{\link{getCactus}}. Returns the original
#' Key if the preceding text is not present.
#'
#' @usage trimKey(InChIKey)
#'
#' @param InChIKey InChIKey to trim (or not).
#'
#' @return Returns the trimmed key or the original input.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getCactus}}.
#'
#' @export
#'
#' @examples
#' CactusKey <- getCactus("benzene","stdinchikey")
#' trimKey(CactusKey)
#' trimKey("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#'
trimKey <- function(InChIKey) {
  key <- sub("InChIKey=","",InChIKey)
  return(key)
}


#' Test whether InChIKey follows valid format
#'
#' @description This basic function using string manipulation to test whether an InChIKey
#' is in the expected format. Useful to detect conversions and retrieval results.
#' This function checks if the string contains 27 characters and 2 dashes/3 blocks.
#'
#' @usage InChIKey_test(InChIKey)
#'
#' @param InChIKey InChIKey to test.
#'
#' @return Returns \code{TRUE} or \code{FALSE} (fails test).
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getSuspectInChIKey}}.
#'
#' @export
#'
#' @examples
#' CactusKey <- getCactus("benzene","stdinchikey")
#' InChIKey_test(CactusKey)
#' CactusKey <- trimKey(CactusKey)
#' InChIKey_test(CactusKey)
#' InChIKey_test("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#'
InChIKey_test <- function(InChIKey) {
  Key_test <- FALSE
  if (nchar(InChIKey)==27 && length(strsplit(InChIKey,"-")[[1]]) == 3) {
    Key_test <- TRUE
  }
  return(Key_test)
}


#' Retrieve InChIKey from PubChem CID
#' 
#' Retrieves InChIKey from PubChem CID using PUG REST.
#' 
#' This searches PubChem for the CID and returns the InChIKey. The function 
#' may work for other output, but this has not been tested.
#' 
#' @usage getPCInChIKey(query, from = "cid", to="InChIKey")
#' 
#' @param query ID to be converted
#' @param from Type of input ID (default \code{cid}, i.e. PubChem CID)
#' @return The InChIKey
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' 
#' @examples
#' getPCInChIKey("2519")
#' 
#' @export
getPCInChIKey <- function(query, from = "cid", to="InChIKey")
{
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
  url <- paste(baseURL, from, query, "property", to, "json",sep="/")
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula,InChIKey/
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=5),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault))
    return(NA)
  
  #titleEntry <- which(unlist(lapply(r$InformationList$Information, function(i) !is.null(i$Title))))
  
  #titleEntry <- titleEntry[which.min(sapply(titleEntry, function(x)r$InformationList$Information[[x]]$inchikey))]
  
  InChIKey <- as.character(unlist(r)[2])
  
  if(is.null(InChIKey)){
    return(NA)
  } else{
    return(InChIKey)
  }
}


#' Retrieve PubChem CIDs from SMILES or InChI
#' 
#' Retrieves PubChem CIDs from SMILES or InChI using PUG REST.
#' By "default", this also returns the CID of the prefered tautomer
#' according to current PubChem standardization procedures
#' 
#' This searches PubChem by SMILES and/or InChI and returns CIDs. The function 
#' may work for other output, but this has not been tested.
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCID.smiles(query, from = "smiles", to="cids", timeout=30)
#' 
#' @param query SMILES or InChI (or other) to be converted
#' @param from Type of input ID (default \code{"smiles"}, i.e. SMILES, alternative
#' \code{"inchi"} to retrieve via InChI, or others (untried))
#' @param timeout Timeout in seconds for the query
#' @return The matching CID
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' 
#' @examples
#' getPCID.smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
#' 
#' @export
getPCID.smiles <- function(query, from = "smiles", to="cids", timeout=30)
{
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", to, "/JSON?", from, "=", 
                URLencode(query,reserved=TRUE))
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON?smiles=CC(=O)OC1=CC=CC=C1C(=O)O
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(url,timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault))
    return(NA)
  
  PCID <- as.character(unlist(r)[1])
  
  if(PCID==0){
    return(NA)
  } else{
    return(PCID)
  }
}


# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/22206/xrefs/PatentID,RN,PubMedID/JSON



#' Retrieve PubMedID and PatentID Counts from PubChem
#' 
#' Retrieves PubMedIDs and PatentIDs from PubChem and 
#' returns a total count of each, using PUG REST.
#' Default behaviour is to retrieve by CID. 
#' 
#' For this function to work as expected, only a unique search query should be
#' used (e.g. CID). At this stage, only counts of PubMedIDs and Patent Counts 
#' will be returned. 
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCxrefs.count(query, from = "cid", xrefs="PatentID,PubMedID",
#' timeout=30)
#' 
#' @param query string of the identifier (CID, InChIKey, name) to be converted
#' @param from Type of input ID (default \code{"cid"}, i.e. PubChem Compound ID, 
#' alternative \code{"inchikey"} to retrieve via InChIKey, \code{"name"} 
#' to retrieve via name (caution), \code{"inchi"} to retrieve by InChI, 
#' or others (untried)).
#' @param xrefs The type of reference information to retrieve. Default 
#' \code{"PatentID,PubMedID"}. While any other string terms recognised by 
#' PUG REST can be used, these are currently ignored. 
#' @param timeout The timeout, in seconds. For records with many entries, 
#' this may need to be increased - e.g. if only one \code{NA} is returned. 
#' @return A list containing the PCID and total counts of patent IDs and PubMed IDs
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' PubChem PUG REST:
#' \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
#' 
#' @examples
#' getPCxrefs.count("22206")
#' getPCxrefs.count("FZXISNSWEXTPMF-UHFFFAOYSA-N",from="inchikey")
#' #a non-live record (returns NAs)
#' getPCxrefs.count("4644")
#' #a CID with many entries (tests timeout)
#' getPCxrefs.count("2244")
#' 
#' @export
getPCxrefs.count <- function(query, from = "cid", xrefs="PatentID,PubMedID",timeout=30)
{
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", query, "/xrefs/", xrefs, "/JSON")
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/22206/xrefs/PatentID,RN,PubMedID/JSON
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault)) {
    Counts <- list()
    Counts[['PCID']] <- as.numeric(query)
    Counts[['PatentCount']] <- 0
    Counts[['PubMedRefCount']] <- 0
    return(Counts)
  } else {
    PCID <- r$InformationList$Information[[1]]$CID
    PatentCount <- length(r$InformationList$Information[[1]]$PatentID)
    PubMedRefCount <- length(r$InformationList$Information[[1]]$PubMedID)
    
    if (is.null(PCID)) {
      PCID <- NA
    }
    if (is.null(PatentCount)) {
      PatentCount <- NA
    }
    if (is.null(PubMedRefCount)) {
      PubMedRefCount <- NA
    }
    
    Counts <- list()
    Counts[['PCID']] <- PCID
    Counts[['PatentCount']] <- PatentCount
    Counts[['PubMedRefCount']] <- PubMedRefCount
    return(Counts)
    
  }
}



#' Retrieve Property Information from PubChem
#' 
#' Retrieves property data for MetFrag from PubChem using PUG REST.
#' Default behaviour is to retrieve by CID or InChIKey. Properties
#' returned are MolecularFormula, IsomericSMILES, InChI, InChIKey,
#' IUPACName and MonoisotopicMass
#' 
#' For this function to work as expected, only one search entry should be
#' used (e.g. one CID or InChIKey). The URL can accept comma separated CIDs, 
#' but this is currently ignored downstream. 
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCproperty.MF(query, from = "cid", timeout=10)
#' 
#' @param query string of the identifier (CID, InChIKey) to be converted
#' @param from Type of input ID (default \code{"cid"}, i.e. PubChem Compound ID, 
#' alternative \code{"inchikey"} to retrieve via InChIKey, or others (untried)).
#' @param timeout The timeout, in seconds.  
#' @return A list containing the related property information
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' PubChem PUG REST:
#' \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
#' 
#' @examples
#' getPCproperty.MF("22206")
#' getPCproperty.MF("FZXISNSWEXTPMF-UHFFFAOYSA-N",from="inchikey")
#' getPCproperty.MF("2244")
#' #a non-live record 
#' getPCproperty.MF("4644")
#' 
#' @export
getPCproperty.MF <- function(query, from = "cid", timeout=10)
{
  property="MolecularFormula,IsomericSMILES,InChI,InChIKey,IUPACName,MonoisotopicMass"
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", query, "/property/", property, "/JSON")
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/22206/xrefs/PatentID,RN,PubMedID/JSON
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault)) {
    Properties <- list()
    Properties[['PCID']] <- NA
    Properties[['MolecularFormula']] <- NA
    Properties[['IsomericSMILES']] <- NA
    Properties[['InChI']] <- NA
    Properties[['InChIKey']] <- NA
    Properties[['IUPACName']] <- NA
    Properties[['MonoisotopicMass']] <- NA
    return(Properties)
  } else {
    PCID <- r$PropertyTable$Properties[[1]]$CID
    MolecularFormula <- r$PropertyTable$Properties[[1]]$MolecularFormula
    IsomericSMILES <- r$PropertyTable$Properties[[1]]$IsomericSMILES
    InChI <- r$PropertyTable$Properties[[1]]$InChI
    InChIKey <- r$PropertyTable$Properties[[1]]$InChIKey
    IUPACName <- r$PropertyTable$Properties[[1]]$IUPACName
    MonoisotopicMass <- r$PropertyTable$Properties[[1]]$MonoisotopicMass
    
    if (is.null(PCID)) {
      PCID <- NA
    }
    if (is.null(MolecularFormula)) {
      MolecularFormula <- NA
    }
    if (is.null(IsomericSMILES)) {
      IsomericSMILES <- NA
    }
    if (is.null(InChI)) {
      InChI <- NA
    }
    if (is.null(InChIKey)) {
      InChIKey <- NA
    }
    if (is.null(IUPACName)) {
      IUPACName <- NA
    }
    if (is.null(MonoisotopicMass)) {
      MonoisotopicMass <- NA
    }
    
    Properties <- list()
    Properties[['PCID']] <- PCID
    Properties[['MolecularFormula']] <- MolecularFormula
    Properties[['IsomericSMILES']] <- IsomericSMILES
    Properties[['InChI']] <- InChI
    Properties[['InChIKey']] <- InChIKey
    Properties[['IUPACName']] <- IUPACName
    Properties[['MonoisotopicMass']] <- MonoisotopicMass
    return(Properties)
    
  }
}


#' Retrieve Record Title (Compound Name) from PubChem
#' 
#' Retrieves record title (title name) from PubChem using PUG REST.
#' Default behaviour is to retrieve by CID or InChIKey. 
#' 
#' For this function to work as expected, only one search entry should be
#' used (e.g. one CID or InChIKey). The URL can accept comma separated CIDs, 
#' but this is currently ignored downstream. 
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCdesc.title(query, from = "cid", timeout=10)
#' 
#' @param query string of the identifier (CID, InChIKey) to be converted
#' @param from Type of input ID (default \code{"cid"}, i.e. PubChem Compound ID, 
#' alternative \code{"inchikey"} to retrieve via InChIKey, or others (untried)).
#' @param timeout The timeout, in seconds.  
#' @return A list containing the CID and related record title (Name)
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' PubChem PUG REST:
#' \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
#' 
#' @examples
#' getPCdesc.title("22206")
#' getPCdesc.title("FZXISNSWEXTPMF-UHFFFAOYSA-N",from="inchikey")
#' getPCdesc.title("2244")
#' getPCdesc.title("1983")
#' #a non-live record 
#' getPCdesc.title("4644")
#' 
#' @export
getPCdesc.title <- function(query, from = "cid", timeout=10)
{
  #description="Title"
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", query, "/description/JSON")
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1983/description/JSON
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault)) {
    Desc <- list()
    Desc[['PCID']] <- NA
    Desc[['Title']] <- NA
    return(Desc)
  } else {
    PCID <- r$InformationList$Information[[1]]$CID
    Title <- r$InformationList$Information[[1]]$Title

    if (is.null(PCID)) {
      PCID <- NA
    }
    if (is.null(Title)) {
      Title <- NA
    }

    Desc <- list()
    Desc[['PCID']] <- PCID
    Desc[['Title']] <- Title
    return(Desc)
    
  }
}

#' Retrieve SMILES from Compound Name from PubChem
#' 
#' Retrieves the Isomeric SMILES via name (synonyms) from PubChem using PUG REST.
#' Note that more than one SMILES can be returned and the "best match" is
#' usually at the top (i.e. n=1). 
#' 
#' For this function to work as expected, only one search entry should be
#' used with the default "name" and "isomericsmiles". n should only be adjusted
#' if you know what you are doing (i.e. trust n=1 first).
#' The URL can accept comma separated CIDs, but this is currently 
#' ignored downstream. 
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCproperty.IsoSMILES(query, from = "name", to = "isomericsmiles", n=1, timeout=30)
#' 
#' @param query string of the compound name for which you need SMILES
#' @param from Type of input ID (default \code{"name"} should be kept for this 
#' function to work as expected).
#' @param to Type of output desired (default \code{"isomericsmiles"} should be 
#' kept for this function to work as expected).
#' @param n Default \code{n=1}. For some names, multiple matches will be 
#' returned and n>1 can be returned by defining n. \code{n=1} is best match.  
#' @param timeout The timeout, in seconds.  
#' @return A list containing the CID and the Isomeric SMILES
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' PubChem PUG REST:
#' \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
#' 
#' @examples
#' getPCproperty.IsoSMILES("aspirin")
#' #an example with multiple results available (up to 3)
#' getPCproperty.IsoSMILES("carnitine")
#' getPCproperty.IsoSMILES("carnitine",n=1)
#' getPCproperty.IsoSMILES("carnitine",n=3)
#' # a nonsense example
#' getPCproperty.IsoSMILES("blah")
#' 
#' @export
getPCproperty.IsoSMILES <- function(query, from = "name", to = "isomericsmiles", n=1, timeout=30)
{
  #build URL
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", query, "/property/", to, "/JSON")
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/isomericsmiles/JSON
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault)) {
    Properties <- list()
    Properties[['PCID']] <- NA
    Properties[['IsomericSMILES']] <- NA
    return(Properties)
  } else {
    PCID <- r$PropertyTable$Properties[[n]]$CID
    IsomericSMILES <- r$PropertyTable$Properties[[n]]$IsomericSMILES
    
    if (is.null(PCID)) {
      PCID <- NA
    }
    if (is.null(IsomericSMILES)) {
      IsomericSMILES <- NA
    }
    Properties <- list()
    Properties[['PCID']] <- PCID
    Properties[['IsomericSMILES']] <- IsomericSMILES
    return(Properties)
    
  }
}


#' Retrieve various CID types from CID via PubChem
#' 
#' Retrieves the various types of related CIDs from a query CID from 
#' PubChem using PUG REST. See details. 
#' 
#' For this function to work as expected, only one search entry should be
#' used. The URL can accept comma separated CIDs, but this is currently 
#' ignored downstream. 
#' Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance. 
#' 
#' @usage getPCIDs.CIDtype(query, type="parent", from = "cid", to = "cids", timeout=30)
#' 
#' @param query Input CID (as string) to search
#' @param from Type of input ID (default \code{"cid"} should be kept for this 
#' function to work as expected).
#' @param to Type of output desired (default \code{"cids"} should be 
#' kept for this function to work as expected).
#' @param timeout The timeout, in seconds.  
#' @return A list containing the related CIDs of the desired type
#' 
#' @details 
#' PubChem have a lot of related CIDs, for instance for this record: 
#' \url{https://pubchem.ncbi.nlm.nih.gov/compound/1234#section=Related-Compounds}.
#' This function enables you to retrieve CIDs by these different types:
#' original, parent, component, similar_2d, similar_3d, same_stereo, 
#' same_isotopes, same_connectivity, same_tautomer, same_parent, 
#' same_parent _stereo, same_parent _isotopes, same_parent _connectivity, 
#' same_parent _tautomer (if more options are available but not implemented 
#' here, this will fail the input tests, pls post an issue).
#' If you have a mixture, "component" gives you the components of the mixture. 
#' If you have an individual component, "component" gives you all the mixtures 
#' containing this component. 
#' If the CID is a "parent", it is the neutralized form. 
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#' 
#' @references 
#' PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' PubChem PUG REST:
#' \url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
#' 
#' @examples
#' # The original returns the input
#' getPCIDs.CIDtype("3053",type="original")
#' # Find the parent CID (neutral version) for a salt
#' getPCIDs.CIDtype("167781",type="parent")
#' # If the parent is not available, go to "component" and get the bits
#' getPCIDs.CIDtype("104265",type="parent")
#' # [1] NA
#' getPCIDs.CIDtype("104265",type="component")
#' # [1] 13360  1004
#' # For deprecated records, only "original" works
#' getPCIDs.CIDtype("4644",type="parent")
#' # [1] NA
#' getPCIDs.CIDtype("4644",type="original")
#' # [1] 4644
#' getPCIDs.CIDtype("4644",type="component")
#' # [1] NA
#' 
#' @export
getPCIDs.CIDtype <- function(query, type="parent", from = "cid", to = "cids", timeout=30)
{
  # test type parameters
  if(!(type %in% c("original","parent","component", "preferred",
                   "similar_2d", "similar_3d", 
                   "same_stereo", "same_isotopes", "same_connectivity", "same_tautomer",
                   "same_parent", "same_parent_stereo", "same_parent_isotopes", 
                   "same_parent_connectivity", "same_parent_tautomer"))) {
    stop("Incorrect type: select one of original, parent, component, preferred, similar_2d, 
          similar_3d, same_stereo, same_isotopes, same_connectivity, same_tautomer,
          same_parent, same_parent_stereo, same_parent_isotopes, 
          same_parent_connectivity, same_parent_tautomer")
  }
  #build URL
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
  url <- paste0(baseURL, from, "/", query, "/", to, "/JSON?cids_type=", type)
  #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2735102/cids/JSON?cids_type=component
  
  
  errorvar <- 0
  currEnvir <- environment()
  
  tryCatch(
    url_data <- getURL(URLencode(url),timeout=timeout),
    error=function(e){
      currEnvir$errorvar <- 1
    })
  
  if(errorvar){
    return(NA)
  }
  
  # This happens if the PCID is not found:
  r <- fromJSON(url_data)
  
  if(!is.null(r$Fault)) {
    CIDs <- NA
    return(CIDs)
  } else {
    CIDs <- r$IdentifierList$CID
    
    if (is.null(CIDs)) {
      CIDs <- NA
    }
    return(CIDs)
  }
}





# # Loops that need to be turned into functions:
# for (i in 1:length(cmpds$ID)) {
#   cas_rn <-  as.character(cmpds$Dat_PSMV_Anhang_1.CASNr[i])
#   cas_rn2 <- as.character(cmpds$Dat_Bioz_Infos.CASNr[i])
#   cmpd_name <- as.character(cmpds$Name[i])
#   cmpd_name2 <- as.character(cmpds$IUPAC.Bezeichnung[i])
#   #retreive keys for all entries if key is not yet filled
#   if (!is.na(cmpds$InChIKey[i]) && InChIKey_test(cmpds$InChIKey[i])) {
#     cmpds$InChIKey[i] <- cmpds$InChIKey[i]
#   } else  if (CAS_test(cas_rn)) {
#     inchikey1 <- trimKey(getCactus(cas_rn,"stdinchikey"))
#   } else {
#     inchikey1 <- "NA"
#   }
#   if (CAS_test(cas_rn2)) {
#     inchikey2 <- trimKey(getCactus(cas_rn2,"stdinchikey"))
#   } else {
#     inchikey2 <- "NA"
#   }
#   inchikey3 <- trimKey(getCactus(cmpd_name,"stdinchikey"))
#   inchikey4 <- trimKey(getCactus(cmpd_name2,"stdinchikey"))
#   # export only one (you could also export all, just create 4 InChIKey fields...)
#   if (InChIKey_test(inchikey1)) {
#     cmpds$InChIKey[i] <- inchikey1
#   } else if (InChIKey_test(inchikey2)) {
#     cmpds$InChIKey[i] <- inchikey2
#   } else if (InChIKey_test(inchikey3)) {
#     cmpds$InChIKey[i] <- inchikey3
#   } else if (InChIKey_test(inchikey4)) {
#     cmpds$InChIKey[i] <- inchikey4
#   } else {
#     cmpds$InChIKey[i] <- "NA"
#   }
# }
#
# inchi <- getCactus(trimKey(inchikey),"stdinchi")
# inchi
# smiles <- getCactus(trimKey(inchikey),"smiles")
# smiles
#
# for (i in 1:length(cmpds$ID_Substanz)) {
#   inchikeyi <- as.character(cmpds$InChIKey[i])
#   if (InChIKey_test(inchikeyi)) {
#     inchi <- getCactus(inchikeyi,"stdinchi")
#     cmpds$InChI[i] <- inchi
#     smiles <- getCactus(inchikeyi,"smiles")[1]
#     cmpds$SMILES[i] <- smiles
#   } else {
#     inchi <- NA
#     smiles <- NA
#   }
#   # get the formula from the InChI
#   if ((nchar(inchi) > 3)&& !is.na(inchi)) {
#     molform <- MolFormFromInChI(inchi)
#     cmpds$MolFormula[i] <- molform
#   } else {
#     molform <- NA
#   }
#   # get exact masses etc.
#   if (!is.na(smiles) && nchar(smiles)>=1) {
#     mass_data <- getSuspectFormulaMass(smiles)
#     cmpds$MolFormula2[i] <- mass_data$MolForm
#     cmpds$NeutralExactMass[i] <- mass_data$Monoiso_mass
#     cmpds$MpHp_mass[i] <- mass_data$MpHp_mass
#     cmpds$MpNH4_mass[i] <- mass_data$MpNH4_mass
#     cmpds$MpNa_mass[i] <- mass_data$MpNa_mass
#     cmpds$MmHm_mass[i] <- mass_data$MmHm_mass
#     logPs <- getCDKlogPs(smiles)
#     cmpds$CDKXlogP[i] <- logPs$CDKXlogP
#     cmpds$CDKAlogP[i] <- logPs$CDKAlogP
#   } else {
#     mass_data <- NA
#     logPs <- NA
#   }
#
# }

# from Barbara toxins:
# info_to_check <- read.csv("PhytotoxinsTest_ES_OnlyWInfo.csv")
# info_to_check$Formula_InChI_std <- ""
# info_to_check$Formula_Match <- ""
#
# formula_to_std <- as.character(info_to_check$Formula_InChI[17])
# checked_formula <- check_chemform(isotopes, formula_to_std)$new_formula
# check_text <- checked_formula == as.character(info_to_check$MolForm[17])
# checked_formula == formula_to_std
#
# for (i in 1:length(info_to_check$Formula_InChI)) {
#   formula_to_std <- as.character(info_to_check$Formula_InChI[i])
#   checked_formula <- check_chemform(isotopes, formula_to_std)$new_formula
#   check_text <- checked_formula == as.character(info_to_check$MolForm[i])
#   #checked_formula == formula_to_std
#   info_to_check$Formula_InChI_std[i] <- checked_formula
#   info_to_check$Formula_Match[i] <- check_text
#
# }




# log Kow from CDK
# CDK_XlogP <- get.xlogp(mol)

# inchi <- get.inchi('CCC')
# inchik <- get.inchi.key('CCC')



# Comments:
# standardization:
# http://ambit.sourceforge.net/download_ambitcli.html

# remove stereochemistry in SMILES:
# babel -i smi -o can -xi
#
# remove charges etc in InChI:
#C:\Program Files (x86)\OpenBabel-2.3.2>obabel -:"[Na+].CC(C)C(=O)C([O-])=O" -oinchi -r -xT /nochg


