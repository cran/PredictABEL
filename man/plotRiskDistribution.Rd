\name{plotRiskDistribution}
\alias{plotRiskDistribution}
\title{Function to plot histogram of risks separated for individuals with and 
without the outcome of interest.}
\usage{plotRiskDistribution(data, cOutcome, risks, interval, rangexaxis, 
rangeyaxis, plottitle, xlabel, ylabel, labels, fileplot, plottype)}
\description{Function to plot histogram of risks separated for individuals with 
and without the outcome of interest.}
\value{The function creates the histogram of risks separated for individuals 
with and without the outcome of interest.}
\keyword{hplot}
\seealso{\code{\link{plotROC}}, \code{\link{riskScore}}}
\arguments{\item{data}{Data frame or numeric matrix that includes the outcome and 
predictor variables.}
\item{cOutcome}{Column number of the outcome variable.}
\item{risks}{Risk of each individual. It
is specified by either a vector of risk scores or a vector of predicted risks.}
\item{interval}{Size of the risk intervals. For example, \code{interval=.1} 
will construct the following intervals for predicted risks: 
\code{0-0.1, 0.1-0.2,..., 0.9-1}.}
\item{rangexaxis}{Range of the x-axis. Specification of \code{rangexaxis} is optional.}
\item{rangeyaxis}{Range of the y-axis.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "Histogram of risks".}
\item{xlabel}{Label of x-axis. Specification of \code{xlabel} is optional. Default is "Risk score".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Percentage".}
\item{labels}{Labels given to the groups of individuals without and 
with the outcome of interest. Specification of \code{labels} is optional. Default is
\code{c("Without outcome", "With outcome")}.}
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
# specify column number of the outcome variable
cOutcome <- 2

# fit a logistic regression model
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk <- predRisk(riskmodel)

# specify the size of each interval
interval <- .05
# specify label of x-axis
xlabel <- "Predicted risk"
# specify label of y-axis
ylabel <- "Percentage"
# specify range of x-axis
xrange <- c(0,1)
# specify range of y-axis
yrange <- c(0,40)
# specify title for the plot
maintitle <- "Distribution of predicted risks"
# specify labels
labels <- c("Without outcome", "With outcome")

# produce risk distribution plot
plotRiskDistribution(data=ExampleData, cOutcome=cOutcome,
risks=predRisk, interval=interval, plottitle=maintitle, rangexaxis=xrange,
rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)}
