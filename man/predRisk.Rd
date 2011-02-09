\name{predRisk}
\alias{predRisk}
\title{Function to compute predicted risks for all individuals in the dataset.}
\usage{predRisk(riskModel, data, cID, filename)}
\description{Function to compute predicted risks for all individuals in the dataset.}
\details{The function computes predicted risks from a specified logistic regression model.   
The function \code{\link{fitLogRegModel}} can be used to construct such a model.}
\value{The function returns a vector of predicted risks.}
\keyword{htest}
\seealso{\code{\link{fitLogRegModel}}, \code{\link{plotCalibration}}, 
\code{\link{plotROC}}, \code{\link{plotPriorPosteriorRisk}}}
\arguments{\item{riskModel}{Name of logistic regression model that can be fitted using  
the function \code{\link{fitLogRegModel}}.}
\item{data}{Data frame or matrix that includes the outcome, ID number and 
predictor variables.}
\item{cID}{Column number of ID variable. The ID number and predicted risks 
will be saved under \code{filename}. When \code{cID} is not specified, the output is not saved.}
\item{filename}{Name of the output file in which the ID number and 
estimated predicted risks will be saved. The file is saved in the working  
directory as a txt file. When no \code{filename} is specified, the output is not saved.}}
\examples{# specify dataset with outcome and predictor variables
data(ExampleData)
# specify column number of the outcome variable
cOutcome <- 2
# specify column number of ID variable
cID <- 1
# specify column numbers of non-genetic predictors
cNonGenPred <- c(3:10)
# specify column numbers of non-genetic predictors that are categorical
cNonGenPredCat <- c(6:8)
# specify column numbers of genetic predictors
cGenPred <- c(11,13:16)
# specify column numbers of genetic predictors that are categorical
cGenPredCat <- c(0)

# fit logistic regression model
riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome,
cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)

# obtain predicted risks 
predRisk <- predRisk(riskModel=riskmodel, filename="name.txt")}
