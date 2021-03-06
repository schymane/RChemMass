% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChemicalCuration.R
\name{getPCxrefs.count}
\alias{getPCxrefs.count}
\title{Retrieve PubMedID and PatentID Counts from PubChem}
\usage{
getPCxrefs.count(query, from = "cid", xrefs="PatentID,PubMedID",
timeout=30)
}
\arguments{
\item{query}{string of the identifier (CID, InChIKey, name) to be converted}

\item{from}{Type of input ID (default \code{"cid"}, i.e. PubChem Compound ID, 
alternative \code{"inchikey"} to retrieve via InChIKey, \code{"name"} 
to retrieve via name (caution), \code{"inchi"} to retrieve by InChI, 
or others (untried)).}

\item{xrefs}{The type of reference information to retrieve. Default 
\code{"PatentID,PubMedID"}. While any other string terms recognised by 
PUG REST can be used, these are currently ignored.}

\item{timeout}{The timeout, in seconds. For records with many entries, 
this may need to be increased - e.g. if only one \code{NA} is returned.}
}
\value{
A list containing the PCID and total counts of patent IDs and PubMed IDs
}
\description{
Retrieves PubMedIDs and PatentIDs from PubChem and 
returns a total count of each, using PUG REST.
Default behaviour is to retrieve by CID.
}
\details{
For this function to work as expected, only a unique search query should be
used (e.g. CID). At this stage, only counts of PubMedIDs and Patent Counts 
will be returned. 
Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance.
}
\examples{
getPCxrefs.count("22206")
getPCxrefs.count("FZXISNSWEXTPMF-UHFFFAOYSA-N",from="inchikey")
#a non-live record (returns NAs)
getPCxrefs.count("4644")
#a CID with many entries (tests timeout)
getPCxrefs.count("2244")

}
\references{
PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 

PubChem PUG REST:
\url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
}
\author{
Emma Schymanski <emma.schymanski@uni.lu>
}
