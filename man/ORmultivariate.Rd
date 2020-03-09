\name{ORmultivariate}
\alias{ORmultivariate}
\title{Function to obtain multivariate odds ratios from a logistic regression model.}
\usage{ORmultivariate(riskModel, filename)}
\description{The function estimates multivariate (adjusted) odds ratios (ORs) with    
95\% confidence intervals (CIs) for all the genetic and non-genetic variables 
in the risk model.}
\details{The function requires that first a logistic regression  
model is fitted either by using \code{GLM} function or the function 
\code{\link{fitLogRegModel}}. In addition to the multivariate ORs, 
the function returns summary statistics of model performance, namely the Brier 
score and the Nagelkerke's \eqn{R^2} value.  
The Brier score quantifies the accuracy of risk predictions by comparing 
predicted risks with observed outcomes at individual level (where outcome 
values are either 0 or 1). The Nagelkerke's \eqn{R^2} value indicates the percentage of variation   
of the outcome explained by the predictors in the model.}
\value{The function returns: 
\item{Predictors Summary}{OR with 95\% CI and corresponding  p-values for 
each predictor in the model}
\item{Brier Score}{Brier score}
\item{Nagelkerke Index}{Nagelkerke's \eqn{R^2} value}}
\keyword{htest}
\references{Brier GW. Verification of forecasts expressed in terms of probability. 
Monthly weather review 1950;78:1-3.   


Nagelkerke NJ. A note on a general definition of the coefficient
of determination. Biometrika 1991;78:691-692.}
\seealso{\code{\link{fitLogRegModel}}}
\arguments{\item{riskModel}{Name of logistic regression model that can be fitted using  
the function \code{\link{fitLogRegModel}}.}
\item{filename}{Name of the output file in which the multivariate    
ORs will be saved. If no directory is specified, the file is 
saved in the working directory as a txt file.   
When \code{filename} is not specified, the output is not saved.}}
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

# obtain multivariate OR(95\% CI) for all predictors of the fitted model
ORmultivariate(riskModel=riskmodel)}
