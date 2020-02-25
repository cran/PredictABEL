\name{plotRiskscorePredrisk}
\alias{plotRiskscorePredrisk}
\title{Function to plot predicted risks against risk scores.}
\usage{plotRiskscorePredrisk(data, riskScore, predRisk, plottitle, xlabel, 
ylabel, rangexaxis, rangeyaxis, filename, fileplot, plottype)}
\description{This function is used to make a plot of predicted risks against risk scores.}
\details{The function creates a plot of predicted risks against risk scores. 
Predicted risks can be obtained using the functions
\code{\link{fitLogRegModel}} and \code{\link{predRisk}}
or be imported from other methods or packages. 
The function \code{\link{riskScore}} can be
used to compute unweighted or weighted risk scores.}
\value{The function creates a plot of predicted risks against risk scores.}
\keyword{hplot}
\seealso{\code{\link{riskScore}}, \code{\link{predRisk}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{riskScore}{Vector of (weighted or unweighted) genetic risk scores.}
\item{predRisk}{Vector of predicted risks.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "Risk score predicted risk plot".}
\item{xlabel}{Label of x-axis. Specification of \code{xlabel} is optional. Default is "Risk score".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Predicted risk".}
\item{rangexaxis}{Range of the x axis. Specification of \code{rangexaxis} is optional.}
\item{rangeyaxis}{Range of the y axis. Specification of \code{rangeyaxis} is optional. Default is \code{c(0,1)}.}
\item{filename}{Name of the output file in which risk scores and
predicted risks for each individual will be saved. If no directory is 
specified, the file is saved in the working directory as a txt file. 
When no \code{filename} is specified, the output is not saved.}
\item{fileplot}{Name of the output file that contains the plot. The file is
saved in the working directory in the format specified under \code{plottype}. Example:
\code{fileplot="plotname"}. Note that the extension is not specified here. 
When \code{fileplot} is not specified, the plot is not saved.}
\item{plottype}{The format in which the plot is saved. Available formats are
wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
eps or pdf. For example, \code{plottype="eps"} will save the plot in eps format.
When \code{plottype} is not specified, the plot will be saved in jpg format.}}
\examples{# specify dataset with outcome and predictor variables
data(ExampleData)

# fit a logistic regression  model
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk <- predRisk(riskmodel)

# specify column numbers of genetic predictors
cGenPred <- c(11:16)

# function to compute unweighted genetic risk scores
riskScore <- riskScore(weights=riskmodel, data=ExampleData, 
cGenPreds=cGenPred, Type="unweighted")

# specify range of x-axis
rangexaxis <- c(0,12)   
# specify range of y-axis
rangeyaxis <- c(0,1)     
# specify label of x-axis
xlabel <- "Risk score"     
# specify label of y-axis
ylabel <- "Predicted risk" 
# specify title for the plot
plottitle <- "Risk score versus predicted risk"

# produce risk score-predicted risk plot
plotRiskscorePredrisk(data=ExampleData, riskScore=riskScore, predRisk=predRisk, 
plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis, 
rangeyaxis=rangeyaxis)}
