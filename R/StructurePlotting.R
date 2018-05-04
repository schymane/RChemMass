# Plotting functions for chemical structures
# Compiled from various scripts
# E. Schymanski, 15/2/2017

# Package dependencies (should go in description)
# library(rcdk)
# library(rcdklibs)
# library(curl)
# library(rsvg)


#' Render SMILES into 2D image for plotting via rcdk
#'
#' @description This function uses the rcdk to parse the smiles into a mol, with options
#' to switch kekulise (aromaticity detection) on or off and to define desired
#' coordinates. Output requires that plot.new has been called, i.e. this is
#' designed to be used directly during plotting.
#' This function uses default depiction options. For more options with
#' rcdk>3.4.1, use \code{\link{renderSMILES.rcdk}}. Best results with rcdk>3.4.1 and
#' rcdklibs>2.0.
#'
#' @usage renderSMILES.rcdk.default(smiles, kekulise=TRUE, coords=c(0,0,100,100))
#'
#' @param smiles A valid SMILES code for rendering (e.g. \code{"c1ccccc1"}).
#' @param kekulise If \code{TRUE}, performs CDK aromaticiy detection, which is recommended.
#' Setting \code{FALSE} can force rendering of invalid SMILES with undefined bond orders.
#' Older rcdk versions may behave differently.
#' @param coords This is used to control the size of the image within the plot. Values
#' \code{c(xmin,ymin,xmax,ymax)}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @details More information about aromaticity: \url{https://github.com/rajarshi/cdkr/issues/49}
#'
#' @seealso \code{\link{renderSMILES.CDKdepict}}, \code{\link{renderSMILES.rcdk}}
#'
#' @return Returns an image for use during plotting
#' @export
#'
#' @examples
#' smiles <- "OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O"
#'   plot.new()
#'   plot.window(xlim=c(0,200), ylim=c(0,100))
#'   renderSMILES.rcdk(smiles,kekulise=FALSE)
#'   renderSMILES.rcdk(smiles,kekulise=TRUE)
#'
renderSMILES.rcdk.default <- function(smiles, kekulise=TRUE, coords=c(0,0,100,100)) {
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise=kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not rendered: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      print("Replotting with kekulise=FALSE")
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    # img <- view.image.2d(mol, width = 150, height = 150)
    # this is just to make sure I see something
    if (length(img)>2) {
      rasterImage(img, coords[1],coords[2], coords[3],coords[4])
    }
  }
}

#' Render SMILES into 2D image for plotting via rcdk
#'
#' @description This function uses the rcdk to parse the smiles into a mol, with options
#' to switch aromaticity detection on or off and to define desired
#' coordinates, as well as depiction control. Compatible only with rdck above 3.4.1. Default
#' parameters produce the same as \code{\link{renderSMILES.rcdk.default}}.
#' Output requires that plot.new has been called, i.e. this is
#' designed to be used directly during plotting.
#' Users with rcdk<3.4.1 should use \code{\link{renderSMILES.rcdk.default}} or
#' \code{renderSMILES.CDKdepict} to view latest CDK updates.
#' As opposed to older versions, kekulise=FALSE now produces strange behaviour for some
#' SMILES.
#'
#' @usage renderSMILES.rcdk(smiles, kekulise=TRUE, coords=c(0,0,100,100),
#' width=200, height=200, zoom=1.3,style="cow", annotate="off", abbr="on",suppressh=TRUE,
#' showTitle=FALSE, smaLimit=100, sma=NULL)
#'
#' @param smiles A valid SMILES code for rendering (e.g. \code{"c1ccccc1"}).
#' @param kekulise If \code{TRUE}, performs CDK aromaticiy detection, which can lead to
#' undesirable plotting results with older versions. Setting \code{FALSE} renders
#' the SMILES code as is but will plot aromatic structures with undefined bond orders.
#' @param coords This is used to control the size of the image within the plot. Values
#' \code{c(xmin,ymin,xmax,ymax)}.
#' @param width Defines width of depiction object
#' @param height Defines height of depiction object
#' @param zoom Controls the size of the image, default 2, increasing to 5 results in a large
#' structure in the web browser. Warnings at zoom>=10.
#' @param style The plotting style, one of 5 options \code{("cow","cob","bow","wob","nob")},
#' i.e. color on white, color on black, black on white, white on black, neon on black. The
#' appearance of \code{SMARTS} and code{annotate} options change with style.
#' @param annotate (optional) Features to interpret the structure. Default \code{"none"}, other
#' options Atom Numbers, Atom Mapping, Color Map, Atom Value
#' (\code{"number", "mapidx", "colmap", "atomvalue"}) respectively.
#' @param abbr Default \code{"off"}, this controls whether the structure is displayed "as is"
#' \code{"off"} or whether groups (\code{"groups"}), reagents (\code{"reagents"}) or both reagents
#' and groups (\code{"on"}) are abbreviated.
#' @param suppressh Default \code{TRUE} suppresses the display of Hs. \code{FALSE} does the oppposite.
#' @param showTitle Default \code{FALSE} turns title display off (no title is included)
#' @param smaLimit sets the limit for \code{SMARTS} patterns to highlight
#' @param sma (optional) \code{SMARTS} codes can be entered to highlight parts of the molecule.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @details More information about aromaticity: \url{https://github.com/rajarshi/cdkr/issues/49}
#'
#' @seealso \code{\link{view.molecule.2d}} \code{\link{renderSMILES.CDKdepict}}
#'
#' @return Returns an image for use during plotting
#' @export
#'
#' @examples
#' smiles <- "OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O"
#' plot.new()
#' plot.window(xlim=c(0,200), ylim=c(0,100))
#' renderSMILES.rcdk(smiles,kekulise=FALSE)
#' renderSMILES.rcdk(smiles,kekulise=TRUE)
#' renderSMILES.rcdk(smiles,kekulise=TRUE, abbr="off")
#'
renderSMILES.rcdk <- function(smiles, kekulise=TRUE, coords=c(0,0,100,100), width=200, height=200,
                              zoom=1.3,style="cow", annotate="off", abbr="on",suppressh=TRUE,
                              showTitle=FALSE, smaLimit=100, sma=NULL) {
  dep <- get.depictor(width = width, height = height, zoom = zoom, style = style, annotate = annotate,
                      abbr = abbr, suppressh = suppressh, showTitle = showTitle, smaLimit = smaLimit,
                      sma = NULL)
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise=kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol, depictor=dep))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not rendered: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol, depictor=dep))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    if (length(img)>2) {
      rasterImage(img, coords[1],coords[2], coords[3],coords[4])
    }
  }
}



#' Build URL for CDK Depict
#'
#' @description This builds a URL for the CDK Depict website for direct viewing or for
#' retrieving svgs for \code{renderSMILES.CDKdepict}. Contains all the options
#' in Version 0.2, reverting to default options if invalid entries are specified.
#'
#' @usage buildCDKdepictURL(smiles, style="bow",title=NULL, abbr="off", smarts=NULL,
#' suppressh=TRUE,showtitle=TRUE,zoom=2,annotate="none")
#'
#' @param smiles A valid \code{SMILES} code or a \code{CxSMILES} for advanced plotting.
#' The latter is not yet accessible in R and must be compiled manually or externally.
#' @param style The plotting style, one of 5 options \code{("cow","cob","bow","wob","nob")},
#' i.e. color on white, color on black, black on white, white on black, neon on black. The
#' appearance of \code{SMARTS} and code{annotate} options change with style.
#' @param title (optional) title to display under structure.
#' @param abbr Default \code{"off"}, this controls whether the structure is displayed "as is"
#' \code{"off"} or whether groups (\code{"groups"}), reagents (\code{"reagents"}) or both reagents
#' and groups (\code{"on"}) are abbreviated.
#' @param smarts (optional) \code{SMARTS} codes can be entered to highlight parts of the molecule.
#' @param suppressh Default \code{TRUE} suppresses the display of Hs. \code{FALSE} does the oppposite.
#' @param showtitle Default \code{TRUE} displays the title below the molecule. \code{FALSE} turns
#' this off.
#' @param zoom Controls the size of the image, default 2, increasing to 5 results in a large
#' structure in the web browser. Warnings at zoom>=10.
#' @param annotate (optional) Features to interpret the structure. Default \code{"none"}, other
#' options Atom Numbers, Atom Mapping, Color Map, Atom Value
#' (\code{"number", "mapidx", "colmap", "atomvalue"}) respectively.
#'
#' @author Emma Schymanski (R wrapper, <emma.schymanski@@uni.lu>), John Mayfield (CDK Depict)
#'
#' @seealso \code{\link{renderSMILES.CDKdepict}} to use the resulting structure in plotting,
#' \code{\link{renderSMILES.rcdk}} for alternative rendering.
#'
#' @return Returns a URL for direct use in a browser, or for rendering in R.
#' @export
#'
#' @examples
#' buildCDKdepictURL("c1ccccc1")
#' buildCDKdepictURL("CN1C=NC2=C1C(=O)N(C(=O)N2C)C",style="cow",title="caffeine",
#'                   abbr="off",smarts="c1cncnc1",suppressh=FALSE,showtitle=TRUE,zoom=5,annotate="number")
#' buildCDKdepictURL("CCOCCOCCO",style="cow")
#' buildCDKdepictURL("CCOCCOCCO |Sg:n:3,4,5::ht|",style="cow", smarts="CCO")
#' buildCDKdepictURL("CCOCCOCCO |Sg:n:3,4,5::ht|",style="cow", smarts="CCO",title="PEGn")
#' buildCDKdepictURL("OS(=O)(=O)c1ccc(cc1)C([R1]C(=O)O)[R2]C(=O)O", style="cow",title="SPADCs",
#' abbr="off",suppressh=TRUE,showtitle=TRUE,zoom=5,annotate="number",smarts="C(=O)O")
#' buildCDKdepictURL("OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O |Sg:n:15:m:ht,Sg:n:11:n:ht|", style="cow",
#' title="SPADCs, n+m=0-5",abbr="off",suppressh=TRUE,showtitle=TRUE,zoom=5,annotate="number",smarts="C(=O)O")
#'
buildCDKdepictURL <- function(smiles, style="bow",title=NULL, abbr="off", smarts=NULL,
                              suppressh=TRUE,showtitle=TRUE,zoom=2,annotate="none") {
  # this builds the URL with varying complexity depending on input
  url_base <- "https://cdkdepict-openchem.rhcloud.com/depict/"
  # style is front of the URL, before SMILES
  # values: color on white, color on black, black on white, white on black, neon on black
  #   c("cow","cob","bow","wob","nob")
  style_entries <- c("cow", "cob", "bow", "wob", "nob")
  if (length(grepl(style, style_entries))<1) {
    warning("Style setting mismatch, defaulting to \"black on white\"")
    style <- "bow"
  }
  # bow/svg?smi=
  url_base <- paste0(url_base,style,"/svg?smi=")
  # now get to SMILES
  smiles_string <- as.character(smiles)
  # if "CxSmiles", this will have the coords as part of the smiles...

  #   # see if coords exist - this comes after smiles. If null, do nothing.
  #   if (!is.null(coords)) {
  #     # this assumes "coords" has the contents needed without seperators
  #     smiles_string <- paste0(smiles_string, " |",as.character(coords),"|")
  #   }
  # now include the title
  # NOTE: may have to check that spacing is correct with CxSMILES
  if (!is.null(title)) {
    # if title has a plus, substitute this first
    title <- sub("+","%2B",title,fixed=TRUE)
    smiles_string <- paste0(smiles_string, " ", as.character(title))
  }
  url_base <- paste0(url_base,URLencode(smiles_string))
  # now define the options at the end
  url_end <- "&abbr="
  #abbr: 4 settings: None, Groups, Reagents, Reagents and Groups
  # URL code: ( "off", "groups", "reagents", "on")
  abbr_entries <- c("off", "groups", "reagents" , "on")
  if (length(grepl(abbr, abbr_entries))<1) {
    warning("abbr setting mismatch, defaulting to \"off\"")
    abbr <- "off"
  }
  # Suppress Hs or not...
  if (!suppressh) {
    url_end <- paste0(url_end,abbr,"&suppressh=false","&showtitle=")
  } else {
    url_end <- paste0(url_end,abbr,"&suppressh=true","&showtitle=")
  }
  # Show the title
  if (!showtitle) {
    url_end <- paste0(url_end,"false")
  } else {
    url_end <- paste0(url_end,"true")
  }
  # if smarts are defined, add this field
  if (!is.null(smarts)) {
    url_end <- paste0(url_end,"&sma=",URLencode(as.character(smarts)))
  }
  # calculate the zoom - warn if too large
  if (!is.numeric(zoom) || zoom<0) {
    zoom <- 2
    warning("Non-numeric or negative value for zoom, defaulting to 2")
  } else if (zoom>10) {
    warning("Large value for zoom, normal range is 2-5, warning >10")
  }
  url_end <- paste0(url_end,"&zoom=", as.character(zoom),"&annotate=")
  # now finish with the annotation
  # Options: None, Atom Numbers, Atom Mapping, Color Map, Atom Value
  # c("none", "number", "mapidx", "colmap", "atomvalue")
  annotate_entries <- c("none", "number", "mapidx", "colmap", "atomvalue")
  if (length(grepl(annotate, annotate_entries))<1) {
    warning("Annotate setting mismatch, defaulting to \"none\"")
    annotate <- "none"
  }
  url_end <- paste0(url_end,annotate)

  url_base <- paste0(url_base, url_end)
  # replace spaces with %20
  url_base <- gsub(" ","%20",url_base)
  return(url_base)
}


#' Render SMILES into 2D image for plotting via CDK Depict
#'
#' @description This function uses the URL built with \code{buildCDKdepictURL} to create an
#' image for direct use in plotting, via a temporary svg file and given
#' coordinates. Output requires that plot.new has been called.
#'
#' @usage renderSMILES.CDKdepict(depictURL,coords=c(0,0,100,100), filename="tmp.svg")
#'
#' @param depictURL A URL created with \code{buildCDKdepictURL}
#' @param coords This is used to control the size of the image within the plot. Values
#' \code{c(xmin,ymin,xmax,ymax)}.
#' @param filename Name of the temporary svg file for saving the image. Can be used
#' to save all image files if renamed for every instance.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{buildCDKdepictURL}} to create the URL, \code{\link{renderSMILES.rcdk}}
#' for alternative plotting direct from SMILES.
#'
#' @return Returns an image for use during plotting
#' @export
#'
#' @examples
#' plot.new()
#' plot.window(xlim=c(0,200), ylim=c(0,100))
#' test_url <- buildCDKdepictURL("OS(=O)(=O)c1ccc(cc1)C(CC(=O)O)CC(=O)O |Sg:n:15:m:ht,Sg:n:11:n:ht|",
#' style="cow",title="SPADCs, n+m=0-5", abbr="off",suppressh=TRUE,showtitle=TRUE,zoom=5,annotate="off",
#' smarts="C(=O)O")
#' renderSMILES.CDKdepict(test_url)
#'
renderSMILES.CDKdepict <- function(depictURL,coords=c(0,0,100,100), filename="tmp.svg") {
  curl_download(depictURL,filename)
  img <- rsvg(filename)
  if (length(img)>2) {
    rasterImage(img,xleft=coords[1],ybottom=coords[2],xright=coords[3],ytop=coords[4])
  }
}

