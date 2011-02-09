\name{plotCalibration}
\alias{plotCalibration}
\title{Function for calibration plot and Hosmer-Lemeshow goodness of fit test.}
\usage{plotCalibration(data, cOutcome, predRisk, groups, 
rangeaxis, plottitle, xlabel, ylabel, filename, fileplot, plottype)}
\description{The function produces a calibration plot and provides Hosmer-Lemeshow 
goodness of fit test statistics.}
\details{Hosmer-Lemeshow test statistic is a measure of the fit 
of the model, comparing observed and predicted risks across subgroups of 
the population. The default number of groups is 10. 

The function requires the outcome of interest and predicted risks of 
all individuals. Predicted risks can be obtained from the 
functions \code{\link{fitLogRegModel}} and \code{\link{predRisk}} or 
be imported from other packages or methods.}
\value{The function creates a calibration plot and returns the following measures: 
\item{Chi_square}{Chi square value of  Hosmer-Lemeshow test}
\item{df}{Degrees of freedom, which is \code{(groups-2)} where \code{groups:} number 
of groups}
\item{p_value}{p-value of Hosmer-Lemeshow test for goodness of fit}}
\keyword{htest}
\keyword{hplot}
\references{Hosmer DW, Hosmer T, Le Cessie S, Lemeshow S. A comparison of
goodness-of-fit tests for the logistic regression model. Stat Med 1997;
16:965-980.}
\seealso{\code{\link{predRisk}}}
\arguments{\item{data}{Data frame or numeric matrix that includes the outcome and 
predictor variables.}
\item{cOutcome}{Column number of the outcome variable.}
\item{predRisk}{Vector of predicted risks of all individuals in the dataset.}
\item{groups}{Number of groups considered in 
Hosmer-Lemeshow test. Specification of \code{groups} is optional (default \code{groups} is 10).}
\item{rangeaxis}{Range of x-axis and y-axis. Specification of \code{rangeaxis} is optional. Default is \code{c(0,1)}.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} is optional. Default is "Calibration plot".}
\item{xlabel}{Label of x-axis Default. Specification of \code{xlabel} is optional. Default is "Predicted risk".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} is optional. Default is "Observed risk".}
\item{filename}{Name of the output file in which the calibration table is saved. 
The file is saved as a txt file in the working directory. When  no 
\code{filename} is specified, the output is not saved.}
\item{fileplot}{Name of the file that contains the calibation plot.
The file is saved in the working directory in the format specified under \code{plottype}. Example:
\code{fileplot="plotname"}. Note that the extension is not specified here.
When \code{fileplot} is not specified, the plot is not saved.}
\item{plottype}{The format in which the plot is saved. Available formats are
wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
eps or pdf. Foe example, \code{plottype="eps"} will save the plot in eps format.
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

# specify range of x-axis and y-axis  
rangeaxis <- c(0,1) 
# specify number of groups for Hosmer-Lemeshow test            
groups <- 10 

# compute calibration measures and produce calibration plot                          
plotCalibration(data=ExampleData, cOutcome=cOutcome, predRisk=predRisk, 
groups=groups, rangeaxis=rangeaxis, filename="calibration.txt")}
