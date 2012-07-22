\name{ExampleData}
\alias{ExampleData}
\docType{data}
\title{A hypothetical dataset that is used to demonstrate all functions.}
\description{
ExampleData is a hypothetical dataset constructed to demonstrate all functions 
in the package. ExampleData is a data frame containing a binary outcome variable 
(e.g., disease present/absent) and genetic and non-genetic predictor variables 
for 10,000 persons. In this dataset, column 1 is the ID number, column 2 is the
outcome variable, columns 3-10 are non-genetic variables and columns 11-16 are  
genetic variables.
  
 }
\usage{data(ExampleData)}
%\details{
%}
%\source{
%}
%\references{
%}
\examples{

data(ExampleData)
# show first 5 records (rows) of the dataset
head(ExampleData,5)
}
\keyword{datasets}
