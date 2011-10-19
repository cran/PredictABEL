\name{simulatedDataset}
\alias{simulatedDataset}
\title{Function to construct a simulated dataset containing individual  
genotype data, estimated genetic risk and disease status.}
\usage{simulatedDataset(ORfreq, poprisk, popsize, filename)}
\description{Construct a dataset that contains individual genotype data, genetic risk,   
and disease status for a hypothetical population. 
The dataset is constructed using simulation in such a way that the frequencies 
and odds ratios (ORs) of the genetic variants and the population disease risk 
computed from this dataset are the same as specified by the input parameters.}
\details{The function will execute when the matrix with odds ratios and frequencies, 
population disease risk and the number of individuals are specified. \cr 

The simulation method is described in detail in the references. \cr 


The method assumes that (i) the combined effect of the genetic variants 
on disease risk follows a multiplicative (log additive) risk model; 
(ii) genetic variants inherit independently, that is no linkage disequilibrium 
between the variants; (iii) genetic variants have independent effects on the 
disease risk, which indicates no interaction among variants; and (iv) all 
genotypes and allele proportions are in Hardy-Weinberg equilibrium. 
Assumption (ii) and (iv) are used to generate the genotype data, and assumption 
(ii) and (iii) are used to calculate disease risk.


Simulating the dataset involves three steps: (1) modelling genotype data, 
(2) modelling disease risks, and (3) modelling disease status. Brief 
descriptions of these steps are as follows:


(1) Modelling genotype data: For each variant the genotype 
frequencies are either specified or calculated from the allele frequencies 
using Hardy-Weinberg equilibrium. Then, the genotypes for each genetic 
variant are randomly distributed without replacement over all individuals.


(2) Modelling disease risks: For the calculation of the individual disease 
risk, Bayes' theorem is used, which states that the posterior odds of disease 
are obtained by multiplying the prior odds by the likelihood ratio (LR) of 
the individual genotype data. The prior odds are calculated from the 
population disease risk or disease prevalence 
(prior odds= prior risk/ (1- prior risk)) and the posterior odds are converted 
back into disease risk (disease risk= posterior odds/ (1+ posterior odds)). 
Under the no linkage disequilibrium (LD) assumption, the LR of a genetic profile  
is obtained by multiplying the LRs of the single genotypes that are included in 
the risk model. The LR of a single genotype is calculated using frequencies 
and ORs of genetic variants and population disease risk. See references 
for more details.


(3) Modelling disease status: To model disease status, we used a procedure 
that compares the estimated disease risk of each subject to a randomly drawn value 
between 0 and 1 from a uniform distribution. A subject is assigned to the 
group who will develop the disease when the disease risk is higher than the 
random value and to the group who will not develop the disease when the risk 
is lower than the random value. 

This procedure ensures that for each genomic profile, the percentage of 
people who will develop the disease equals the population disease risk 
associated with that profile, when the subgroup of individuals with that 
profile is sufficiently large.}
\value{The function returns:
\item{Dataset}{A data frame or matrix that includes genotype data, 
genetic risk and disease status for a hypothetical population. 
The dataset contains (4 + number of genetic variants included) columns,  
in which the first column is the un-weighted risk score, which is the sum 
of the number of risk alleles for each individual, the third column is the 
estimated genetic risk, the forth column is the individual disease status expressed 
as '0' or '1', indicating without or with the outcome of interest, and the fifth until 
the end column are genotype data for the variants expressed as '0', '1' or '2', 
which indicate the number of risk alleles present in each individual for the genetic variants.}}
\keyword{models}
\references{Janssens AC, Aulchenko YS, Elefante S, Borsboom GJ, Steyerberg EW, 
van Duijn CM. Predictive testing for complex diseases using multiple genes: 
fact or fiction? Genet Med. 2006;8:395-400.


Janssens AC, Moonesinghe R, Yang Q, Steyerberg EW, van Duijn CM, Khoury MJ. 
The impact of genotype frequencies on the clinical validity of genomic 
profiling for predicting common chronic diseases. Genet Med. 2007;9:528-35.


van der Net JB, Janssens AC, Sijbrands EJ, Steyerberg EW. Value of genetic 
profiling for the prediction of coronary heart disease. 
Am Heart J. 2009;158:105-10.


van Zitteren M, van der Net JB, Kundu S, Freedman AN, van Duijn CM, 
Janssens AC. Genome-based prediction of breast cancer risk in the general 
population: a modeling study based on meta-analyses of genetic associations. 
Cancer Epidemiol Biomarkers Prev. 2011;20:9-22.}
\arguments{\item{ORfreq}{Matrix with ORs and frequencies of the genetic variants. 
The matrix contains four columns in which the first two describe ORs and the 
last two describe the corresponding frequencies. The number of rows in this 
matrix is same as the number of genetic variants included. Genetic variants 
can be specified as per genotype, per allele, or as dominant/ recessive 
effect of the risk allele. When per genotype data are used, OR of the 
heterozygous and homozygous risk genotypes are mentioned in the first two 
columns and the corresponding genotype frequencies are mentioned in the last 
two columns. When per allele data are used, the OR and frequency of the risk 
allele are specified in the first and third column and the remaining two cells 
are coded as '1'.  Similarly, when dominant/ recessive effects of the risk 
alleles are used, the 
OR and frequency of the dominant/ recessive variant are specified in the first 
and third column, and the remaining two cells are coded as '0'.}
\item{poprisk}{Population disease risk (expressed in proportion).}
\item{popsize}{Total number of individuals included in the dataset.}
\item{filename}{Name of the file in which the dataset will be saved. 
The file is saved in the working directory as a txt file. When no filename 
is specified, the output is not saved.}}
\examples{# specify the matrix containing the ORs and frequencies of genetic variants 
# In this example we used per allele effects of the risk variants
ORfreq<-cbind(c(1.35,1.20,1.24,1.16), rep(1,4), c(.41,.29,.28,.51),rep(1,4))

# specify the population disease risk
popRisk <- 0.3
# specify size of hypothetical population
popSize <- 10000

# Obtain the simulated dataset
Data <- simulatedDataset(ORfreq=ORfreq, poprisk=popRisk, popsize=popSize)

# Obtain the AUC and produce ROC curve
plotROC(data=Data, cOutcome=4, predrisk=Data[,3])}
