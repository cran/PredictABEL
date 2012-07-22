\name{reclassification}
\alias{reclassification}
\title{Function for reclassification table and statistics.}
\usage{reclassification(data, cOutcome, predrisk1, predrisk2, cutoff)}
\description{The function creates a reclassification table and provides statistics.}
\details{The function creates a reclassification table and computes the 
categorical and continuous net reclassification improvement (\code{NRI}) and 
integrated discrimination improvement (\code{IDI}). A reclassification table 
indicates the number of individuals who move to another risk category or remain 
in the same risk category as a result of updating the risk model. Categorical \code{NRI} equal to 
\code{x\%} means that compared with individuals without outcome, 
individuals with outcome were almost \code{x\%} more likely to move up a category than down.
The function also computes continuous \code{NRI}, which does not require any discrete 
risk categories and relies on the proportions of individuals with outcome 
correctly assigned a higher probability and individuals without outcome 
correctly assigned a lower probability by an updated model compared with the 
initial model.
\code{IDI} equal to \code{x\%} means that the difference in average 
predicted risks between the individuals with and without the outcome  
increased by \code{x\%} in the updated model.
The function requires predicted risks estimated by using two separate risk 
models. Predicted risks can be obtained using the functions
\code{\link{fitLogRegModel}} and \code{\link{predRisk}}
or be imported from other methods or packages.}
\value{The function returns the reclassification table, separately  
for individuals with and without the outcome of interest and the following measures: 
\item{NRI (Categorical)}{Categorical Net Reclassification Improvement with 95\% CI and \code{p-value} of the test}
\item{NRI (Continuous)}{Continuous Net Reclassification Improvement with 95\% CI and \code{p-value} of the test}
\item{IDI}{Integrated Discrimination Improvement with 95\% CI and \code{p-value} 
of the test}}
\keyword{htest}
\references{Cook NR. Use and misuse of the receiver operating characteristic
curve in risk prediction. Circulation 2007;115(7):928-935.      

Pencina MJ, D'Agostino RB Sr, D'Agostino RB Jr, Vasan RS.
Evaluating the added predictive ability of a new marker: from
area under the ROC curve to reclassification and beyond. Stat
Med 2008;27(2):157-172; discussion 207-212.}
\seealso{\code{\link{plotDiscriminationBox}}, \code{\link{predRisk}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{cOutcome}{Column number of the outcome variable.}
\item{predrisk1}{Vector of predicted risks of all individuals using initial
model.}
\item{predrisk2}{Vector of predicted risks of all individuals using updated
model.}
\item{cutoff}{Cutoff values for risk categories. 
Define the cut-off values as \code{c(0,...,1)}. 
Multiple values can be defined and always specify \code{0} and \code{1}. 
Example: \code{c(0,.20,.30,1)}}}
\examples{# specify dataset with outcome and predictor variables 
data(ExampleData)
# specify column number of the outcome variable
cOutcome <- 2

# fit logistic regression models
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel1 <- ExampleModels()$riskModel1
riskmodel2 <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)
# specify cutoff values for risk categories
cutoff <- c(0,.10,.30,1)    

# compute reclassification measures
reclassification(data=ExampleData, cOutcome=cOutcome, 
predrisk1=predRisk1, predrisk2=predRisk2, cutoff)}
