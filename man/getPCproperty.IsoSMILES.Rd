% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChemicalCuration.R
\name{getPCproperty.IsoSMILES}
\alias{getPCproperty.IsoSMILES}
\title{Retrieve SMILES from Compound Name from PubChem}
\usage{
getPCproperty.IsoSMILES(query, from = "name", to = "isomericsmiles", n=1, timeout=30)
}
\arguments{
\item{query}{string of the compound name for which you need SMILES}

\item{from}{Type of input ID (default \code{"name"} should be kept for this 
function to work as expected).}

\item{to}{Type of output desired (default \code{"isomericsmiles"} should be 
kept for this function to work as expected).}

\item{n}{Default \code{n=1}. For some names, multiple matches will be 
returned and n>1 can be returned by defining n. \code{n=1} is best match.}

\item{timeout}{The timeout, in seconds.}
}
\value{
A list containing the CID and the Isomeric SMILES
}
\description{
Retrieves the Isomeric SMILES via name (synonyms) from PubChem using PUG REST.
Note that more than one SMILES can be returned and the "best match" is
usually at the top (i.e. n=1).
}
\details{
For this function to work as expected, only one search entry should be
used with the default "name" and "isomericsmiles". n should only be adjusted
if you know what you are doing (i.e. trust n=1 first).
The URL can accept comma separated CIDs, but this is currently 
ignored downstream. 
Thanks to Paul Thiessen and Evan Bolton from PubChem team for assistance.
}
\examples{
getPCproperty.IsoSMILES("aspirin")
#an example with multiple results available (up to 3)
getPCproperty.IsoSMILES("carnitine")
getPCproperty.IsoSMILES("carnitine",n=1)
getPCproperty.IsoSMILES("carnitine",n=3)
# a nonsense example
getPCproperty.IsoSMILES("blah")

}
\references{
PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 

PubChem PUG REST:
\url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
}
\author{
Emma Schymanski <emma.schymanski@uni.lu>
}
