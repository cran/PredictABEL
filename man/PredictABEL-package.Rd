\name{PredictABEL-package}
\alias{PredictABEL-package}
\title{An R package for the analysis of (genetic) risk prediction studies.}
\description{An R package for the analysis of (genetic) risk prediction studies.}
\details{Fueled by the substantial gene discoveries from genome-wide association 
studies, there is increasing interest in investigating the predictive 
ability of genetic risk models. To assess the performance of genetic risk 
models, PredictABEL includes functions for the various measures and plots 
that have been used in empirical studies, including univariate and 
multivariate odds ratios (ORs) of the predictors, the c-statistic (or AUC), 
Hosmer-Lemeshow goodness of fit test, reclassification table, net 
reclassification improvement (NRI) and integrated discrimination 
improvement (IDI). The plots included are the ROC plot, calibration plot, 
discrimination box plot, predictiveness curve, and several risk distributions. 


These functions can be applied to predicted risks that are obtained using 
logistic regression analysis, to weighted or unweighted risk scores, for 
which the functions are included in this package. The functions can also be 
used to assess risks or risk scores that are constructed using other methods, e.g., Cox Proportional 
Hazards regression analysis, which are not included in the current version. 
Risks obtained from other methods can be imported into R for assessment 
of the predictive performance.


The functions to construct the risk models using logistic regression analyses 
are specifically written for models that include genetic variables, 
eventually in addition to non-genetic factors, but they can also be applied 
to construct models that are based on non-genetic risk factors only. \cr 


Before using the functions \code{\link{fitLogRegModel}} for constructing 
a risk model or \code{\link{riskScore}} for computing risk 
scores, the following checks on the dataset are advisable to be done: 

(1) Missing values: The logistic regression analyses and computation of 
the risk score are done only for subjects that have no missing data. In case 
of missing values, individuals with missing data can be removed from the 
dataset or imputation strategies can be used to fill in missing data. 
Subjects with missing data can be removed with the R function \code{na.omit} 
(available in \code{stats} package). 
Example: \code{DataFileNew <- na.omit(DataFile)} 
will make a new dataset (\code{DataFileNew}) with no missing values;


(2) Multicollinearity: When there is strong correlation between the 
predictor variables, regression coefficients may be estimated imprecisely 
and risks scores may be biased because the assumption of independent effects 
is violated. In genetic risk prediction studies, problems with 
multicollinearity should be expected when single nucleotide polymorphisms 
(SNPs) located in the same gene are                                                    
in strong linkage disequilibrium (LD). For SNPs in LD it is common to select 
the variant with the lowest p-value in the model;  


(3) Outliers: When the data contain significant outliers, either clinical 
variables with extreme values of the outcomes or extreme values resulting 
from errors in the data entry, these may impact the construction of the risk models and 
computation of the risks scores. Data should be carefully checked and outliers 
need to be removed or replaced, if justified;

(4) Recoding of data: In the computation of unweighted risk scores, it is assumed 
that the genetic variants are coded \code{0,1,2}  
representing the number of alleles carried. When variants 
are coded \code{0,1} representing a dominant or recessive effect of the alleles, 
the variables need to be recoded before unweighted risk scores can be computed. \cr


To import data into R several alternative strategies can be used. Use the 
\code{Hmisc} package for importing SPSS and SAS data into R. 
Use "\code{ExampleData <- read.table("DataName.txt", header=T, sep="\t")}" for text 
files where variable names are included as column headers and data are 
separated by tabs.  
Use "\code{ExampleData <- read.table("Name.csv", sep=",", header=T)}" 
for comma-separated files with variable names as column headers. 
Use \code{"setwd(dir)"} to set the working directory to "dir". The datafile 
needs to be present in the working directory. \cr 


To export datafiles from R tables to a tab-delimited textfile with the first row as
the name of the variables, 
use "\code{write.table(R_Table, file="Name.txt", row.names=FALSE, sep="\t")}"  and 
when a comma-separated textfile is requested and variable names are provided in the first row,
use "\code{write.table(R_Table, file="Name.csv", row.names=FALSE, sep=",")}".  
When the directory is not specified, the file will be
saved in the working directory. For exporting R data into SPSS, SAS and 
Stata data, use functions in the the \code{foreign} package. \cr

Several functions in this package depend on other R packages:

(1) \code{Hmisc}, is used to compute NRI and IDI;

(2) \code{ROCR}, is used to produce ROC plots;

(3) \code{epitools}, is used to compute  univariate odds ratios;

(4) \code{PBSmodelling}, is used to produce predictiveness curve.}
\note{The current version of the package includes the basic measures 
and plots that are used in the assessment of (genetic) risk prediction models and the
function to construct a simulated dataset that contains individual genotype 
data, estimated genetic risk and disease status, used for the evaluation of 
genetic risk models (see Janssens et al, Genet Med 2006). 
Planned extensions of the package include functions to construct risk 
models using Cox Proportional Hazards analysis for prospective data and 
assess the performance of risk models for time-to-event data.}
\section{Acknowledgements}{The authors would like to acknowledge Lennart Karssen, 
Maksim Struchalin and Linda Broer from the Department of Epidemiology, 
Erasmus Medical Center, Rotterdam for their valuable comments and 
suggestions to make this package.
}
\author{Suman Kundu  

Yurii S. Aulchenko

A. Cecile J.W. Janssens}
\keyword{package}
\references{S Kundu, YS Aulchenko, CM van Duijn, ACJW Janssens. PredictABEL:
an R package for the assessment of risk prediction models.
Eur J Epidemiol. 2011;26:261-4. \cr

ACJW Janssens, JPA Ioannidis, CM van Duijn, J Little, MJ Khoury.
Strengthening the Reporting of Genetic Risk Prediction Studies: The GRIPS
Statement Proposal. Eur J Epidemiol. 2011;26:255-9. \cr

ACJW Janssens, JPA Ioannidis, S Bedrosian, P Boffetta, SM Dolan, N Dowling,
I Fortier, AN. Freedman, JM Grimshaw, J Gulcher, M Gwinn, MA Hlatky, H Janes,
P Kraft, S Melillo, CJ O'Donnell, MJ Pencina, D Ransohoff, SD Schully,
D Seminara, DM Winn, CF Wright, CM van Duijn, J Little, MJ Khoury.
Strengthening the reporting of genetic risk prediction studies
(GRIPS)-Elaboration and explanation. Eur J Epidemiol. 2011;26:313-37. \cr   

Aulchenko YS, Ripke S, Isaacs A, van Duijn CM. GenABEL: an R package for genome-wide
association analysis. Bioinformatics 2007;23(10):1294-6.}
