\name{riskScore}
\alias{riskScore}
\title{Function to compute genetic risk scores.}
\usage{riskScore(weights, data, cGenPreds, Type)}
\description{The function computes unweighted 
or weighted genetic risk scores. The relative effects (or weights) of 
genetic variants can either come from beta coefficients of a risk model 
or from a vector of beta coefficients imported into R, e.g., when beta cofficients are obtained from meta-analysis.}
\details{The function calculates unweighted 
or weighted genetic risk scores. The unweighted genetic risk score is a simple
risk allele count assuming that all alleles have the same effect. For this
calculation, it is required that the genetic variables are coded as the number of risk
alleles. Beta coefficients are used to determine which allele is the risk
allele. When the sign of the beta coefficient is negative, the allele coding
is reversed. The weighted risk score is a sum of the number of risk alleles
multiplied by their beta coefficients.

The beta coefficients can come from two different sources, either beta coefficients of a risk model 
or a vector of beta coefficients imported into R, e.g., when beta cofficients are obtained from meta-analysis.
This vector of beta coefficients 
should be a named vector containing the same names as mentioned in genetic variants.
A logistic regression model can be constructed using \code{\link{fitLogRegModel}}
from this package.}
\note{When a vector of beta coefficients is imported, it should be checked 
whether the DNA strands and the coding of the risk alleles are the same 
as in the study data. The functions are available in the package \code{GenABEL}
to accurately compute risk scores when the DNA strands are different or the risk 
alleles are coded differently in the study data and the data used in meta-analysis.}
\value{The function returns a vector of risk scores.}
\keyword{htest}
\seealso{\code{\link{plotRiskDistribution}}, \code{\link{plotRiskscorePredrisk}}}
\arguments{\item{weights}{The vector that includes the weights given to the genetic 
variants. See details for more informations.}
\item{data}{Data frame or matrix that includes the outcome  
and predictors variables.}
\item{cGenPreds}{Column numbers of the genetic variables on the basis of   
which the risk score is computed.}
\item{Type}{Specification of the type of risk scores that will be computed.
Type can be weighted (\code{Type="weighted"}) or 
unweighted (\code{Type="unweighted"}).}}
\examples{# specify dataset with outcome and predictor variables
data(ExampleData) 
# specify column numbers of genetic predictors
cGenPred <- c(11:16)

# fit a logistic regression model
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel <- ExampleModels()$riskModel2

# compute unweighted risk scores 
riskScore <- riskScore(weights=riskmodel, data=ExampleData, 
cGenPreds=cGenPred, Type="unweighted")}
