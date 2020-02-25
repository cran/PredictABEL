\name{ORunivariate}
\alias{ORunivariate}
\title{Function to compute univariate ORs for genetic predictors.}
\usage{ORunivariate(data, cOutcome, cGenPreds, filenameGeno, filenameAllele)}
\description{The function computes the univariate ORs with 95\% CIs for genetic predictors.}
\details{The function computes the univariate ORs with 95\% CIs for the specified 
genetic variants both per allele and per genotype. The ORs are saved with the data from which they are  
calculated. Genotype frequencies are provided for      
persons with and without the outcome  
of interest. The genotype or allele that is coded as \code{'0'} is considered 
as the reference to computes the ORs.}
\value{The function returns two different tables. One table contains genotype frequencies 
and univariate ORs with 95\% CIs and the other contains allele frequencies and 
univariate ORs with 95\% CIs.}
\keyword{manip}
\seealso{\code{\link{ORmultivariate}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{cOutcome}{Column number of the outcome variable. \code{cOutcome=2}  
means that the second column of the dataset is the outcome variable.}
\item{cGenPreds}{Column numbers of genetic variables for which the ORs 
are calculated.}
\item{filenameGeno}{Name of the output file in which the univariate ORs  
and frequencies per genotype will be saved. The file is saved in the working directory as 
a txt file. When no \code{filenameGeno} is specified, the output is not saved.}
\item{filenameAllele}{Name of the output file in which the univariate ORs and 
frequencies per allele will be saved. The file is saved in the working 
directory as a txt file. When no \code{filenameAllele} is specified, the output is not saved.}}
\examples{# specify dataset with outcome and predictor variables 
data(ExampleData)
# specify column number of the outcome variable
cOutcome <- 2 
# specify column numbers of genetic predictors
cGenPreds <- c(11:13,16)

# compute univariate ORs
ORunivariate(data=ExampleData, cOutcome=cOutcome, cGenPreds=cGenPreds)}
