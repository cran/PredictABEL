\name{fitLogRegModel}
\alias{fitLogRegModel}
\title{Function to fit a logistic regression model.}
\usage{fitLogRegModel(data, cOutcome, cNonGenPreds, cNonGenPredsCat, 
cGenPreds, cGenPredsCat)}
\description{The function fits a standard GLM function for the logistic regression model.}
\details{The function fits a standard GLM function for the logistic regression model.
This function can be used to construct a logistic regression model based on genetic and non-genetic
predictors. The function also allows to enter the genetic predictors 
as a single risk score. For that purpose, the function requires that 
the dataset additionally includes the risk score.
A new dataset can be constructed using
"\code{NewExampleData <- cbind(ExampleData,riskScore)}". 
The genetic risk scores can be obtained
using the function \code{\link{riskScore}} in this package or be
imported from other methods.}
\value{No value returned.}
\keyword{models}
\seealso{\code{\link{predRisk}}, \code{\link{ORmultivariate}}, \code{\link{riskScore}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictor variables.}
\item{cOutcome}{Column number of the outcome variable. \code{cOutcome=2}
means that the second column of the dataset is the outcome variable.
To fit the logistic regression model, the outcome variable needs to be 
(re)coded as \code{1} for the presence and \code{0} for the absence of the
outcome of interest.}
\item{cNonGenPreds}{Column numbers of the non-genetic predictors that are  
included in the model. An example to denote column numbers is  
\code{c(3,6:8,10)}. Choose \code{c(0)} when no non-genetic predictors 
are considered.}
\item{cNonGenPredsCat}{Column numbers of the non-genetic predictors that  
are entered as categorical variables in the model. When non-genetic 
predictors are not specified as being categorical they are treated as     
continuous variables in the model. If no non-genetic predictors are 
categorical, denote \code{c(0)}.}
\item{cGenPreds}{Column numbers of the genetic predictors or genetic risk score.
Denote \code{c(0)}
when the prediction model does not consider 
genetic predictors or genetic risk score.}
\item{cGenPredsCat}{Column numbers of the genetic predictors that are 
entered as categorical variables in the model. When SNPs are considered as 
categorical, the model  
will estimate effects per genotype. Otherwise, SNPs are considered as  
continuous variables for which  the model will estimate an allelic effect.  
Choose c(0) when no genetic predictors are considered as categorical 
or when genetic predictors are entered as a risk score into the model.}}
\examples{# specify dataset with outcome and predictor variables 
data(ExampleData)
# specify column number of outcome variable
cOutcome <- 2 
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

# show summary details for the fitted risk model
summary(riskmodel)}
