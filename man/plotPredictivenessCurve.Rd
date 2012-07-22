\name{plotPredictivenessCurve}
\alias{plotPredictivenessCurve}
\title{Function for predictiveness curve.}
\usage{plotPredictivenessCurve(predrisk, rangeyaxis, labels, plottitle, 
xlabel, ylabel, fileplot, plottype)}
\description{The function creates a plot of cumulative percentage 
of individuals to the predicted risks.}
\details{The Predictiveness curve is a plot of cumulative percentage 
of individuals to the predicted risks. Cumulative percentage indicates 
the percentage of individual that has a predicted risk equal or lower 
than the risk value. 
Predicted risks can be obtained using the functions 
\code{\link{fitLogRegModel}} and \code{\link{predRisk}} 
or be imported from other methods or packages.}
\value{The function creates a predictiveness curve.}
\keyword{hplot}
\seealso{\code{\link{predRisk}}}
\references{Pepe MS, Feng Z, Huang Y, et al. Integrating the predictiveness   
of a marker with its performance as a classifier. 
Am J Epidemiol 2008;167:362-368.}
\arguments{\item{predrisk}{Vector of predicted risk. When multiple curves need to 
be presented in one plot, specify multiple vectors of predicted 
risks as \code{predrisk=cbind(predrisk1, predrisk2,...,predriskn)}.}
\item{rangeyaxis}{Range of the y axis. Default \code{rangeyaxis} is \code{c(0,1)}.}
\item{labels}{Label(s) given to the predictiveness curve(s). Specification of \code{labels} is optional.
When specified, the \code{labels} should be in the same order as specified in \code{predrisk}.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "Predictiveness curve".}
\item{xlabel}{Label of x-axis. Specification of \code{xlabel} is optional. Default is "Cumulative percentage".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Predicted risks".}
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

# fit logistic regression models
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel1 <- ExampleModels()$riskModel1
riskmodel2 <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)

# specify range of y-axis
rangeyaxis <- c(0,1) 
# specify labels of the predictiveness curves
labels <- c("without genetic factors", "with genetic factors")

# produce predictiveness curves
plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2),
rangeyaxis=rangeyaxis, labels=labels)}
