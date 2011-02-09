\name{plotPriorPosteriorRisk}
\alias{plotPriorPosteriorRisk}
\title{Function to plot posterior risks against prior risks.}
\usage{plotPriorPosteriorRisk(data, priorrisk, posteriorrisk, cOutcome, 
plottitle, xlabel, ylabel, rangeaxis, plotAll=TRUE, labels, 
filename, fileplot, plottype)}
\description{Function to plot posterior risks against prior risks.}
\details{The function creates a plot of posterior risks (predicted risks using 
the updated model) against prior risks (predicted risks using the initial 
model). Predicted risks can be obtained using the functions 
\code{\link{fitLogRegModel}} and \code{\link{predRisk}} or be             
imported from other packages or methods.}
\value{The function creates a plot of posterior risks against prior risks.}
\keyword{hplot}
\seealso{\code{\link{predRisk}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{priorrisk}{Vector of predicted risks based on initial model.}
\item{posteriorrisk}{Vector of predicted risks based on updated model.}
\item{cOutcome}{Column number of the outcome variable.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "PriorPosteriorRisk plot".}
\item{xlabel}{Label of x-axis. Specification of \code{xlabel} is optional. Default is "Prior risk".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Posterior risk".}
\item{rangeaxis}{Range of x-axis and y-axis. Specification of \code{rangeaxis} is optional. Default is \code{c(0,1)}.}
\item{plotAll}{\code{plotAll=TRUE} will create one plot for the 
total population. When  \code{plotAll=FALSE} separate plots will be created 
for individuals with and without the outcome of interest.   
means two separate plots for with and without outcome of interest.}
\item{labels}{Labels given to the groups of individuals without and with 
the outcome of interest. Default \code{labels} is
\code{c("without outcome", "with outcome")}. Note that when 
\code{plotAll=TRUE}, specification of \code{labels} is not necessary.}
\item{filename}{Name of the output file in which prior and posterior 
risks for each individual with the outcome will be saved. If no directory is 
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
# specify column number of outcome variable
cOutcome <- 2

# fit logistic regression models
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel1 <- ExampleModels()$riskModel1
riskmodel2 <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)

# specify label of x-axis 
xlabel <- "Prior risk"                   
# specify label of y-axis
ylabel <- "Posterior risk"               
# specify title for the plot
titleplot <- "Prior versus posterior risk" 
# specify range of the x-axis and y-axis 
rangeaxis <- c(0,1)                       
# labels given to the groups without and with the outcome of interest
labels<- c("without outcome", "with outcome")           

# produce prior risks and posterior risks plot
plotPriorPosteriorRisk(data=ExampleData, priorrisk=predRisk1,
posteriorrisk=predRisk2, cOutcome=cOutcome, xlabel=xlabel, ylabel=ylabel,   
rangeaxis=rangeaxis, plotAll=TRUE, plottitle=titleplot, labels=labels)}
