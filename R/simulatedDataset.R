#' Function to construct a simulated dataset containing individual genotype data, genetic risks and disease status.
#' Construct a dataset that contains genotype data, estimated risk based on
#' genetic variants, and disease status for a hypothetical population.
#' The dataset is constructed using simulation in such a way that the frequencies
#' and odds ratios (OR) of the genetic variants and the population disease risk
#' computed from this dataset are the same as specified by the input parameters.
#'
#' The function will execute when the matrix with odds ratios and frequencies,
#' population disease risk and the number of individuals are specified. \cr
#'
#' The simulation method is described in detail in the references. \cr
#'
#'
#' The method assumes that (i) the combined effect of the genetic variants
#' on disease risk follows a multiplicative (log additive) risk model;
#' (ii) genetic variants inherit independently, that is no linkage disequilibrium
#' between the variants; (iii) genetic variants have independent effects on the
#' disease risk, which indicates no interaction among variants; and (iv) all
#' genotypes and allele proportions are in Hardy-Weinberg equilibrium.
#' Assumption (ii) and (iv) are used to generate the genotype data, and assumption
#'(i), (ii) and (iii) are used to calculate disease risk.
#'
#'
#' Simulating the dataset involves three steps: (1) modelling genotype data,
#' (2) modelling disease risks, and (3) modelling disease status. Brief
#' descriptions of these steps are as follows:
#'
#'
#' (1) Modelling genotype data: For each genetic variant the genotype
#' frequencies are either specified or calculated from the allele frequencies
#' using Hardy-Weinberg equilibrium. Then, the genotypes for each genetic
#' variant are randomly distributed without replacement over all individuals.
#'
#'
#' (2) Modelling disease risks: For the calculation of the individual disease
#' risk, Bayes' theorem is used, which states that the posterior odds of disease
#' are obtained by multiplying the prior odds by the likelihood ratio (LR) of
#' the individual genotype data. The prior odds are calculated from the
#' population disease risk or disease prevalence
#' (prior odds= prior risk/ (1- prior risk)) and the posterior odds are converted
#' back into disease risk (disease risk= posterior odds/ (1+ posterior odds)).
#' Under the no linkage disequilibrium (LD) assumption, the LR is obtained
#' by multiplying the LRs of all individual genotypes that are included in
#' the risk model. The LR of each genotype is calculated using frequencies
#' and ORs of genetic variants and population disease risk. See references
#' for more details.
#'
#'
#' (3) Modelling disease status: To model disease status, we used a procedure
#' that compares the disease risk of each subject to a randomly drawn value
#' between 0 and 1 from a uniform distribution. A subject was assigned to the
#' group who will develop the disease when the disease risk was higher than the
#' random value and to the group who will not develop the disease when the risk
#' was lower than the random value.
#'
#'
#' This procedure ensures that for each genomic profile, the percentage of
#' people who will develop the disease equals the population disease risk
#' associated with that profile, when the subgroup of individuals with that
#' profile is sufficiently large.
#'
#'
#' @param ORfreq Matrix with ORs and frequencies of the genetic variants.
#' The matrix contains four columns in which the first two describe ORs and the
#' last two describe the corresponding frequencies. The number of rows in this
#' matrix is same as the number of genetic variants included. Genetic variants
#' can be specified as per genotype, per allele, or per dominant/ recessive
#' effect of the risk allele. When per genotype data are used, OR of the
#' heterozygous and homozygous risk genotypes are mentioned in the first two
#' columns and the corresponding genotype frequencies are mentioned in the last
#' two columns. When per allele data are used, the OR and frequency of the risk
#' allele are specified in the first and third column and the remaining two cells
#' are coded as '1'.  Similarly, when per dominant/ recessive data are used, the
#' OR and frequency of the dominant/ recessive variant are specified in the first
#' and third column, and the remaining two cells are coded as '0'. \cr
#' Note that, when OR of a genetic variant is less than 1, modify the reference
#' group such that the OR for the new reference group is 1 and above 1 for
#' other groups. Also, change the corresponding frequencies accordingly.
#' @param poprisk Population disease risk (expressed in proportion).
#' @param popsize  Total number of individuals included in the dataset.
#' @param filename Name of the file in which the dataset will be saved.
#' The file is saved in the working directory as a txt file. When no filename
#' is specified, the output is not saved.
#'
#'  @return
#'   The function returns:
#'   \item{Dataset}{A data frame or matrix that includes genotype data,
#' estimated genetic risk and disease status for a hypothetical population.
#' The dataset contains (4+number of genetic variants included) columns. The
#' first column of this dataset is the unweighted risk score, which is the sum
#' of the number of risk alleles for each individual, the third column is the
#' estimated genetic risk, the forth column is the individual disease status with '1'
#' indicates with and '0' as without the outcome of interest, and the fifth until
#' the end column are genotype data for the variants expressed as '0', '1' or '2',
#' which indicate the number of risk alleles for that genetic variant.}
#'
#'
#' @keywords models
#'
#'
#' @references Hanley JA, McNeil BJ. The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. Radiology 1982;143:29-36.
#'
#'
#'  Janssens AC, Aulchenko YS, Elefante S, Borsboom GJ, Steyerberg EW,
#' van Duijn CM. Predictive testing for complex diseases using multiple genes:
#' fact or fiction? Genet Med. 2006;8:395-400.
#'
#'
#' Janssens AC, Moonesinghe R, Yang Q, Steyerberg EW, van Duijn CM, Khoury MJ.
#' The impact of genotype frequencies on the clinical validity of genomic
#' profiling for predicting common chronic diseases. Genet Med. 2007;9:528-35.
#'
#'
#' van der Net JB, Janssens AC, Sijbrands EJ, Steyerberg EW. Value of genetic
#' profiling for the prediction of coronary heart disease.
#' Am Heart J. 2009;158:105-10.
#'
#'
#' van Zitteren M, van der Net JB, Kundu S, Freedman AN, van Duijn CM,
#' Janssens AC. Genome-based prediction of breast cancer risk in the general
#' population: a modeling study based on meta-analyses of genetic associations.
#' Cancer Epidemiol Biomarkers Prev. 2011;20:9-22.
#'
#'
#' @examples
#' # specify the matrix containing the ORs and frequencies of genetic variants.
#' # In this example we used per allele genetic variants
#' ORfreq<-cbind(c( 1.35,1.20,1.24,1.16), rep(1,4), c(.41,.29,.28,.51),rep(1,4))
#'
#' # Obtain the dataset
#' Data <- simulatedDataset(ORfreq=ORfreq, poprisk=.3, popsize=1000)
#'
#' # Obtain the AUC and produce ROC curve
#' plotROC(data=Data, cOutcome=4, predrisk=Data[,3])
#'
"simulatedDataset" <- function(ORfreq, poprisk, popsize, filename)
{
if (missing(poprisk)) {stop("Population disease risk is not specified")}
if (missing(popsize)) {stop("Total number of individuals is not mentioned")}

g <- nrow(ORfreq)
reconstruct.2x2table <- function(p,d,OR,s)
{
a <- 0
b <- 0
c <- (OR*p*s*(1-d)*d*s)/((1-p)*s*(1-d)+OR*p*s*(1-d))
dd <- p*s-c
e <- d*s-c
f <- (1-p)*s-e
tabel <- cbind(a,b,c,dd,e,f,g,OR)
tabel
}
###################################################################
# Reconstruct 2*3 table from OR - no rare disease assumption
###################################################################
reconstruct.2x3table <- function(OR1,OR2,p1,p2,d,s){
        a	<- 1
        eOR	<- 0
        while (eOR<=OR2){
                b	<- p2*s*(1-d)
                snew <- s-a-b
                p1new <-p1/(1-p2)
                dnew <- (d-(a/s))/((d-(a/s))+ ((1-d)-b/s))
                c	<-	(OR1*p1new*snew*(1-dnew)*dnew*snew)/((1-p1new)*snew*(1-dnew)+OR1*p1new*snew*(1-dnew))
                dd	<-	p1new*((1-d)-b/s)*s
                e	<-	(d-(a/s))*s-c
                f	<- ((1-d)-b/s)*s-dd
                eOR	<- (a*f)/(b*e)
                tabel <- cbind(a,b,c,dd,e,f,g,OR1,OR2)
                a	<- a+1
                tabel
                }
        tabel
}
###################################################################
# Reconstruct 2*3 table from OR - based on HWE - no rare disease assumption
###################################################################
reconstruct.2x3tableHWE <- function(OR,p,d,s){
  OR1 <- OR
  OR2 <- OR^2
        p1 <-  2*p*(1-p)
        p2 <-  p*p

        a	<- 1
        eOR	<- 0
        while (eOR<=OR2){
                b	<- p2*s*(1-d)
                snew <- s-a-b
                p1new <-p1/(1-p2)
                dnew <- (d-(a/s))/((d-(a/s))+ ((1-d)-b/s))
                c	<-	(OR1*p1new*snew*(1-dnew)*dnew*snew)/((1-p1new)*snew*(1-dnew)+OR1*p1new*snew*(1-dnew))
                dd	<-	p1new*((1-d)-b/s)*s
                e	<-	(d-(a/s))*s-c
                f	<- ((1-d)-b/s)*s-dd
                eOR	<- (a*f)/(b*e)
                tabel <- cbind(a,b,c,dd,e,f,g,OR1,OR2)
                a	<- a+1
                tabel
                }
        tabel
}
###################################################################
# Adjust such that mean (postp) = pd
###################################################################
# this correction is needed when the number of genes or ORs get large to ensure that the mean (postp) equals the prior risk

adjust.postp <- function (pd, LR){
        odds.diff <- 0
        prior.odds <- pd/(1-pd)
        for (i in (1:100000)) {
        Postp <- (prior.odds*LR)/(1+(prior.odds*LR))
        odds.diff <- (pd-mean(Postp))/ (1-(pd-mean(Postp)))
        prior.odds	<- prior.odds+odds.diff
        if (odds.diff < .0001) break
        }
        Postp
}

################################################################################

func.data <- function(p,d,OR,s,g){
  Data <- matrix (NA,s,4+g)
  Data[,1] <- rep(0,s)
        Data[,2] <- rep(1,s)
        Data[,3] <- rep(0,s)
        i <- 0
        while (i < g){
    i <- i+1
    cells2x3 <- rep(NA,9)
    cells2x3 <- if(p[i,2]==0) {reconstruct.2x2table(p=p[i,1],d,OR=OR[i,1],s)} else {if(p[i,2]==1) {reconstruct.2x3tableHWE(OR=OR[i,1],p=p[i,1],d,s)}
  else {reconstruct.2x3table(OR1=OR[i,1],OR2=OR[i,2],p1=p[i,1],p2=p[i,2],d,s)}}   # reconstruct table for calculation of likelihood ratios for genotypes
      LREE        <- ((cells2x3[1]/d*s)/(cells2x3[2]/(1-d)*s))			# calculate likelihood ratios
      LREe        <- ((cells2x3[3]/d*s)/(cells2x3[4]/(1-d)*s))
      LRee        <- ((cells2x3[5]/d*s)/(cells2x3[6]/(1-d)*s))

 Gene <- if(p[i,2]==0){c(rep(0,((1-p[i,1]-p[i,2])*s)),rep(1,p[i,1]*s),rep(2,p[i,2]*s))}
     else {c(rep(0,((1-p[i,1]*p[i,1]-2*p[i,1]*(1-p[i,1]))*s)),rep(1,2*p[i,1]*(1-p[i,1])*s),rep(2,p[i,1]*p[i,1]*s))}		# create vector of genotypes for all subjects based on hardy-weinberg distribution of alleles
                Filler <- s-length(Gene)                               #soms is Gene 1 te subject te kort en dan werkt het niet
                Gene <- sample(c(Gene,rep(0,Filler)),s,rep=FALSE)
    Data[,4+i] <- Gene
    GeneLR <- ifelse(Gene==0,LRee,ifelse(Gene==1,LREe,LREE))

    Data[,1] <- Data[,1]+Gene
    Data[,2] <- Data[,2]*GeneLR

#        cat(i,"")										# report proces on screen
                }

                Data[,3] <- adjust.postp(pd=d, LR=Data[,2])
                Data[,4]  <- ifelse(runif(s)<=(Data[,3]), 1, 0)
    Data <- as.data.frame(Data)
    Data
    }

 simulatedData <- func.data (p=ORfreq[,c(3,4)],d=poprisk,OR=ORfreq[,c(1,2)],s=popsize,g=nrow(ORfreq))


if (!missing(filename))
        {write.table( simulatedData,file=filename, row.names=TRUE,sep = "\t")  }

 return(simulatedData)
}
