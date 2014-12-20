\name{plotROC}
\alias{plotROC}
\title{Function for a receiver operating characteristic curve (ROC) plot and area under the ROC curve (AUC) value.}
\usage{plotROC(data, cOutcome, predrisk, labels, plottitle, xlabel, ylabel, 
fileplot, plottype)}
\description{The function produces ROC curve and corresponding AUC value with 95\% CI. 
The function can plot one or multiple ROC curves in a single plot.}
\details{The function requirs predicted risks or risk scores and the outcome of 
interest for all individuals. 
Predicted risks can be obtained using the functions
\code{\link{fitLogRegModel}} and \code{\link{predRisk}}
or be imported from other methods or packages.}
\value{The function creates ROC plot and returns AUC value with 95\% CI.}
\keyword{htest}
\keyword{hplot}
\references{Hanley JA, McNeil BJ. The meaning and use of the area under a  
receiver operating characteristic (ROC) curve. Radiology 1982;143:29-36.


Tobias Sing, Oliver Sander, Niko Beerenwinkel, Thomas Lengauer.
ROCR: visualizing classifier performance in R.
Bioinformatics 2005;21(20):3940-3941.}
\seealso{\code{\link{predRisk}}, \code{\link{plotRiskDistribution}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{cOutcome}{Column number of the outcome variable.}
\item{predrisk}{Vector of predicted risk. When multiple curves need to
be presented in one plot, specify multiple vectors of predicted
risks as \code{predrisk=cbind(predrisk1, predrisk2,...,predriskn)}.}
\item{labels}{Label(s) given to the ROC curve(s). Specification of \code{labels} is optional.
When specified, the \code{labels} should be in the same order as specified in \code{predrisk}.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "ROC plot".}
\item{xlabel}{Label of x-axis. Specification of \code{xlabel} is optional. Default is "1- Specificity".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Sensitivity".}
\item{fileplot}{Name of the output file that contains the plot. The file is
saved in the working directory in the format specified under \code{plottype}. Example:
\code{fileplot="plotname"}. Note that the extension is not specified here.
When \code{fileplot} is not specified, the plot is not saved.}
\item{plottype}{The format in which the plot is saved. Available formats are
wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
When \code{plottype} is not specified, the plot will be saved in jpg format.}}
\examples{# specify the arguments in the function to produce ROC plot
# specify dataset with outcome and predictor variables
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

# specify label of the ROC curve
labels <- c("without genetic factors", "with genetic factors")

# produce ROC curve
plotROC(data=ExampleData, cOutcome=cOutcome, 
predrisk=cbind(predRisk1,predRisk2), labels=labels)}
