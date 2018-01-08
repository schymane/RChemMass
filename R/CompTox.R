# Functions for the Dashboard
# E. Schymanski, 10/5/2017

#library(RCurl)

#' Build URL for CompTox Chemistry Dashboard
#'
#' @description This builds a URL for searching the US EPA CompTox Chemistry Dashboard, optionally
#' sending this direct to the default browser.
#'
#' @usage buildCompToxURL(search_string, type="simple", add_utf8=FALSE, browse=FALSE,
#' mass_entry="", single_component=TRUE, no_isotopes=TRUE)
#'
#' @param search_string A string to search CompTox. Must be compatible with one of the \code{type}
#' entries to perform a valid search.
#' @param type Type of search to perform, one of \code{("simple", "name", "cas","DTXSID","InChIKey",
#' "skeleton", "similarity","mass_ppm", "mass_Da", "mass_range", "formula")}. See Details.
#' @param add_utf8 If \code{TRUE}, adds utf8 character into the URL for compatability with IE8.
#' If this is added, it is preferable to copy/paste the link into the browser rather than use
#' \code{browse=TRUE} to avoid encoding issues.
#' @param browse Default \code{FALSE} returns just a URL for this function, if \code{TRUE}, the
#' URL is opened in the default browser as well.
#' @param mass_entry This must be non-empty for all mass searches; for \code{mass_ppm} it is the ppm
#' value to use, for \code{mass_Da} the absolute mass error and for \code{mass_range} it is the
#' higher mass entry in the range (and \code{search_string} must be the lower entry).
#' @param single_component Default \code{TRUE} means only single component records are searched;
#'  \code{FALSE} includes mixtures. Only used in combination with mass and formula searches.
#' @param no_isotopes Default \code{TRUE} eliminates non-standard isotope species. \code{FALSE}
#' includes these. Used only in mass searches.
#'
#' @details
#' Currently the CompTox Dashboard supports name, CAS, DTXSID, InChIKey and InChIKey skeleton
#' searches as a simple text-based search requiring only \code{search_string}.
#'
#' The \code{type="simple"}
#' covers all these categories, alternatively the respective \code{type} can be used as follows:
#'
#' name = \code{"name"}, CAS reference number = \code{"cas"},
#' DSSToxID (Dashboard identifier) = \code{"DTXSID"} (including the letters DTXSID!), InChIKey =
#' \code{InChIKey} and InChIKey skeleton (first block) search = \code{"skeleton"}.
#'
#' The similarity search (\code{type="similarity"}) requires a DTXSID as \code{search_string}.
#' Mass-based searches are a combination of \code{search_string} and \code{mass_entry} (see above).
#'
#' Formula searches require a molecular formula as \code{search_string} and type=\code{"formula"}.
#' Note that a SMILES search is not yet supported in the Dashboard.
#'
#' @author Emma Schymanski (R wrapper, <emma.schymanski@@uni.lu>), Antony Williams (CompTox Dashboard)
#'
#' @return Returns a URL for direct use in a browser, optionally opening it in the default browser
#' @export
#'
#' @examples
#' buildCompToxURL("benzene")
#' buildCompToxURL("nicotine", browse=TRUE)
#'
#' # if you add utf8, better to copy and paste link into browswer rather than use browse=TRUE.
#' buildCompToxURL("1912-24-9", type = "cas", add_utf8 = TRUE)
#'
#' buildCompToxURL("MXWJVTOOROXGIU-UHFFFAOYSA-N", browse=TRUE)
#' buildCompToxURL("SNICXCGAKADSCV", type="skeleton", browse=TRUE)
#'
#' # this searches just for the exact ID
#' buildCompToxURL("DTXSID9020112", browse=TRUE)
#' # this searches for similar substances, appearing lower down in the record
#' buildCompToxURL("DTXSID9020112", type="similarity", browse=TRUE)
#'
#' # Molecular formula search: (default with single component)
#' buildCompToxURL("C8H14ClN5", type="formula", browse=TRUE)
#' buildCompToxURL("C8H14ClN5", type="formula", browse=TRUE, single_component = FALSE)
#' # These two links return the same result, this shouldn't be the case.
#' # Single component part does not seem to be working (also in web interface). Comment submitted.
#' buildCompToxURL("C8H15Cl2N5", type="formula", browse=TRUE, single_component = TRUE)
#' buildCompToxURL("C8H15Cl2N5", type="formula", browse=TRUE, single_component = FALSE)
#'
#' # Mass search with ppm and Da
#' buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="5")
#' buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="5", single_component = FALSE,
#'                 no_isotopes = FALSE)
#' buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="100")
#' buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.01")
#' buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.1")
#' buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.1",
#'                 single_component = FALSE, no_isotopes = FALSE)
#' # Search by mass range
#' buildCompToxURL("162.115", type="mass_range", browse=TRUE, mass_entry="162.116",
#'                 single_component = FALSE, no_isotopes = FALSE)
#'
buildCompToxURL <- function(search_string, type="simple", add_utf8=FALSE, browse=FALSE,
                            mass_entry="", single_component=TRUE, no_isotopes=TRUE) {
  # this builds the URL with varying complexity depending on input
  url_base1 <- "https://comptox.epa.gov/dashboard/dsstoxdb/results?"
  if (add_utf8) {
    url_base <- paste0(url_base1,"utf8=%E2%9C%93&search=")
    #url_base <- "https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=<U+2713>&search="
    #url_base <- paste0(url_base1, "utf8=✓&search=")
    #only use this option if internet explorer is in use...
    # browseURL(paste0(url_base, "benzene"),browser="C:/Program Files/Internet Explorer/iexplore.exe")
  } else {
    url_base <- paste0(url_base1,"search=")
  }

  # search_type tells what we need to add to the URL
  # values: "simple" means nothing to add to search_string
  # otherwise; "name", "cas", "DTXSID", "InChIKey", "skeleton" are also simple searches
  # "similarity" needs DTXSID and additional text
  # "mass" needs additional conditions
  # "formula" needs additional conditions
  #   c("simple","cas","DTXSID","InChIKey","skeleton", "similarity", "mass", "formula")
  type_entries <- c("simple","name", "cas","DTXSID","InChIKey","skeleton", "similarity",
                    "mass_ppm", "mass_Da", "mass_range", "formula")
  if (length(grep(type, type_entries))<1) {
    warning("Type entry mismatch, defaulting to \"simple\", search may not work correctly")
    type <- "simple"
  }
  simple_types <- c("simple","name", "cas","DTXSID","InChIKey","skeleton")
  # if simple search, do this
  if (length(grep(type, simple_types))>=1) {
    url_base <- paste0(url_base, search_string)
  } else if (grepl(type, "similarity")) {
    # similarity search, check if we have a DTXSID in search_string
    if (grepl("DTXSID",search_string)) {
      url_base <- paste0(url_base, search_string,"#similarity")
    } else {
      warning("Search string does not contain DTXSID for similarity search. Defaulting to simple search")
      url_base <- paste0(url_base, search_string)
    }
  } else if (grepl("mass",type)) {
    # mass search - note no "search="
    # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1&single_component=1
    if (add_utf8) {
      url_base <- paste0(url_base1,"utf8=%E2%9C%93&search=")
      #url_base <- paste0(url_base1,"utf8=✓&")
    } else {
      url_base <- url_base1
    }
    # now add the mass conditions
    # exact mass & ppm: &mass=2&mass1=162.115698&mass2=5&ppm=1
    if (grepl("mass_ppm", type)) {
      mass_entry <- as.character(mass_entry)
      url_base <- paste0(url_base, "&mass=2&mass1=",search_string,"&mass2=", mass_entry,"&ppm=1")
      # mass entry has to be the ppm number, as a string
    } else if (grepl("mass_Da", type)) {
      # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=0.01&ppm=0
      mass_entry <- as.character(mass_entry)
      url_base <- paste0(url_base, "&mass=2&mass1=",search_string,"&mass2=", mass_entry,"&ppm=0")
    } else if (grepl("mass_range", type)) {
      # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=1&mass1=162.115&mass2=162.116
      mass_entry <- as.character(mass_entry)
      url_base <- paste0(url_base, "&mass=1&mass1=",search_string,"&mass2=", mass_entry)
      # so here the mass_entry has to be the highest mass
    } else {
      # shouldn't need this case anyway ...
      url_base <- paste0(url_base1,"search=",search_string)
    }

    # add component, isotope options
    if (single_component) {
      url_base <- paste0(url_base, "&single_component=1")
    }
    if (no_isotopes) {
      url_base <- paste0(url_base, "&isotopes=1")
    }

  } else if (grepl(type, "formula")) {
    # formula search - note this needs a different order!
    # https://comptox.epa.gov/dashboard/dsstoxdb/results?formula=1&search=C8H14ClN5&single_component=1
    url_base <- paste0(url_base1, "formula=1&search=", search_string)
    # add component options
    if (single_component) {
      url_base <- paste0(url_base, "&single_component=1")
    }
    # questions: do I have to deal with utf8 (doesn't seem so?), other values for formula?
  } else {
    # shouldn't need this case, it's the same as "simple"
    url_base <- paste0(url_base, search_string)
  }

  if (browse) {
    browseURL(url_base)
  }

  return(url_base)
}

# buildCompToxURL("benzene")
# buildCompToxURL("nicotine", browse=TRUE)
# buildCompToxURL("1912-24-9", type = "cas", add_utf8 = TRUE)
# buildCompToxURL("MXWJVTOOROXGIU-UHFFFAOYSA-N", browse=TRUE)
# buildCompToxURL("SNICXCGAKADSCV", type="skeleton", browse=TRUE)
# # this searches just for the exact ID
# buildCompToxURL("DTXSID9020112", browse=TRUE)
# # this searches for similar substances, appearing lower down in the record
# buildCompToxURL("DTXSID9020112", type="similarity", browse=TRUE)
# # Molecular formula search: (default with single component)
# buildCompToxURL("C8H14ClN5", type="formula", browse=TRUE)
# buildCompToxURL("C8H14ClN5", type="formula", browse=TRUE, single_component = FALSE)
# # These two links return the same result, this shouldn't be the case.
# # Single component part does not seem to be working (also in web interface)
# buildCompToxURL("C8H15Cl2N5", type="formula", browse=TRUE, single_component = TRUE)
# buildCompToxURL("C8H15Cl2N5", type="formula", browse=TRUE, single_component = FALSE)
# # Mass search with ppm
# buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="5")
# # Searched by Mass and single component chemicals, ignoring isotopes: Found 135 results for '162.115698 ± 5 ppm'.
# buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="5", single_component = FALSE,
#                 no_isotopes = FALSE)
# # Searched by Mass: Found 137 results for '162.115698 ± 5 ppm'.
# buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="10")
# # Searched by Mass and single component chemicals, ignoring isotopes: Found 135 results for '162.115698 ± 10 ppm'.
# # coincidence that there's nothnig additional at 10 ppm; it's all either same mass or >11 ppm
# buildCompToxURL("162.115698", type="mass_ppm", browse=TRUE, mass_entry="100")
# # Searched by Mass and single component chemicals, ignoring isotopes: Found 481 results for '162.115698 ± 100 ppm'.
# buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.01")
# # Searched by Mass and single component chemicals, ignoring isotopes: Found 227 results for '162.115698 ± 0.01 amu'
# buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.1")
# # Searched by Mass and single component chemicals, ignoring isotopes: Found 2030 results for '162.115698 ± 0.1 amu'
# buildCompToxURL("162.115698", type="mass_Da", browse=TRUE, mass_entry="0.1",
#                 single_component = FALSE, no_isotopes = FALSE)
# # Searched by Mass: Found 2116 results for '162.115698 ± 0.1 amu'.
# buildCompToxURL("162.115", type="mass_range", browse=TRUE, mass_entry="162.116",
#                 single_component = FALSE, no_isotopes = FALSE)
# # Searched by Mass: Found 137 results for '162.115 to 162.116 amu'


# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=atrazine
# # 1912-24-9 - no changes needed
# # DTXSID9020112 - no changes needed
# # MXWJVTOOROXGIU-UHFFFAOYSA-N - no changes needed
# # MXWJVTOOROXGIU - no changes needed
# #
#
# # C8H14ClN5 - needs additional options!
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?formula=1&search=C8H14ClN5&single_component=1
#
# # similarity
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?search=DTXSID90374722#similarity
#
# # CCNC1=NC(NC(C)C)=NC(Cl)=N1 - SMILES doesn't work
#
#
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=nicotine
#
# # InChI and SMILES don't work
# # CN1CCC[C@H]1C1=CN=CC=C1
# #  InChI=1S/C10H14N2/c1-12-7-3-5-10(12)9-4-2-6-11-8-9/h2,4,6,8,10H,3,5,7H2,1H3/t10-/m0/s1
# # SNICXCGAKADSCV-JTQLQIEISA-N
# # 162.115698
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1&single_component=1
#
# # browseURL("https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene")
# # returns https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%C3%A2%C5%93%E2%80%9C&search=benzene
# URLencode("https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene")
# URLdecode("https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene")
#
# # browseURL("https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene", encodeIfNeeded=TRUE)
# # browseURL("https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&search=benzene", encodeIfNeeded=FALSE)
#
#
# # mass searches
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1&single_component=1
# # no single component option
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1
# # no isotopes
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1&isotopes=1
# # both options
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=5&ppm=1&single_component=1&isotopes=1
#
# # Da ..
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=2&mass1=162.115698&mass2=0.01&ppm=0
#
# # Mass Range:
# # https://comptox.epa.gov/dashboard/dsstoxdb/results?utf8=%E2%9C%93&mass=1&mass1=162.115&mass2=162.116&single_component=1
