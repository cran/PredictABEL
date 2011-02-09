\name{plotDiscriminationBox}
\alias{plotDiscriminationBox}
\title{Function for box plots of predicted risks separately for individuals with
 and without the outcome of interest.}
\usage{plotDiscriminationBox(data, cOutcome, predrisk, labels, plottitle, 
ylabel, fileplot, plottype)}
\description{The function produces box plots of predicted risks for individuals with 
and without the outcome of interest and calculates the discrimination slope.}
\details{The discrimination slope is the difference between the mean predicted risks 
of individuals with and without the outcome of interest. Predicted risks 
can be obtained using the 
\code{\link{fitLogRegModel}} and \code{\link{predRisk}} or be 
imported from other programs. The difference between discrimination 
slopes of two separate risk models is equivalent 
to (\code{IDI}) which is discussed 
in details in the \code{\link{reclassification}} function.}
\value{The function creates a box plots of predicted risks for  
individuals with and without the outcome of interest and returns the discrimination slope.}
\keyword{hplot}
\references{Yates JF. External correspondence: decomposition of the mean 
probability score. Organizational Behavior and Human Performance 1982;30:132-156.}
\seealso{\code{\link{reclassification}}, \code{\link{predRisk}}}
\arguments{\item{data}{Data frame or matrix that includes the outcome and 
predictors variables.}
\item{cOutcome}{Column number of the outcome variable.}
\item{predrisk}{Vector of predicted risks.}
\item{labels}{Labels given to the groups of individuals without and with 
the outcome of interest. Specification of \code{label} is optional. 
Default is \code{c("Without disease", "With disease")}.}
\item{plottitle}{Title of the plot. Specification of \code{plottitle} 
is optional. Default is "Box plot".}
\item{ylabel}{Label of y-axis. Specification of \code{ylabel} 
is optional. Default is "Predicted risks".}
\item{fileplot}{Name of the file that contains the plot. The file is
saved in the working directory in the format specified under \code{plottype}. 
Example: \code{fileplot="name"}. Note that the extension is not specified here. 
When \code{fileplot} is not specified, the plot is not saved.}
\item{plottype}{The format in which the plot is saved. Available formats are
wmf, emf, png, jpg, jpeg, bmp, tif, tiff, ps,
eps or pdf. Foe example, \code{plottype="eps"} will save the plot in eps format.
When \code{plottype} is not specified, the plot will be saved in jpg format.}}
\examples{# specify dataset with outcome and predictor variables
data(ExampleData)
# specify column number of outcome variable
cOutcome <- 2

# fit a logistic regression model
# all steps needed to construct a logistic regression model are written in a function
# called 'ExampleModels', which is described on page 4-5
riskmodel <- ExampleModels()$riskModel2

# obtain predicted risks
predRisk <- predRisk(riskmodel)
# specify labels for the groups without and with the outcome of interest 
labels <- c("Without disease", "With disease")    

# produce discrimination box plot     
plotDiscriminationBox(data=ExampleData, cOutcome=cOutcome, predrisk=predRisk, 
labels=labels)}
