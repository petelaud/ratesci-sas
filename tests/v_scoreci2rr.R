#   Validation of SAS SCORECI macro against R functions for RR
#   PropCIs::riskscoreci  (for Koopman interval)
#   DescTools (for Koopman intervals)
#   sasLM (for Koopman intervals)
#   contingencytables (for MN intervals)
#   gsDesign (for Koopman & MN intervals)
#   pairwiseCI (demonstrating mismatch for both)
#   and ratesci::scoreci  (for SCAS interval)
#

install.packages('PropCIs')
library(PropCIs)
install.packages('gsDesign')
library(gsDesign)
install.packages('DescTools')
library(DescTools)
install.packages('contingencytables')
library(contingencytables)
install.packages('sasLM')
library(sasLM)
install.packages('ratesci')
library(ratesci)
install.packages('pairwiseCI')
library(pairwiseCI)

SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2rr.csv")
SASCIs <- read.csv("C:/Mac/Home/Documents/Main/GitHub/ratesci-sas/tests/sasval2rr.csv")
nsamp <- dim(SASCIs)[[1]]
head(SASCIs)

SASCIs[5710:5720, ]

SASCIs$y1 <- SASCIs$e1
SASCIs$y2 <- SASCIs$e0
SASCIs$n2 <- SASCIs$n0

RCIs <- R2CIs <- R3CIs <- R4CIs <- array(NA,dim=c(nsamp,2))
dimnames(RCIs)[[2]]<-c("lclR","uclR")
dimnames(R2CIs)[[2]]<-c("lclR2","uclR2")
dimnames(R3CIs)[[2]]<-c("lclR3","uclR3")
dimnames(R4CIs)[[2]]<-c("lclR4","uclR4")

j <- 5719

# One or more of these is very slow!
# 
for (j in 1:nsamp) {
	try(RCIs[j,] <- unlist(PropCIs::riskscoreci(x1 = SASCIs[j,"e1"],
	                                            n1 = SASCIs[j,"n1"],
	                                            x2 = SASCIs[j,"e0"],
	                                            n2 = SASCIs[j,"n0"],
	                          conf.level = SASCIs[j,"CONFLEV"])))
  try(R2CIs[j,] <- unlist(sasLM::RRmn1(y1 = SASCIs[j,"e1"],
                                       n1 = SASCIs[j,"n1"],
                                       y2 = SASCIs[j,"e0"],
                                       n2 = SASCIs[j,"n0"], 
                                       conf.level = SASCIs[j, "CONFLEV"])[4:5]))
  try(R3CIs[j,] <- unlist(DescTools::BinomRatioCI(x1 = SASCIs[j,"e1"],
                                       n1 = SASCIs[j,"n1"],
                                       x2 = SASCIs[j,"e0"],
                                       n2 = SASCIs[j,"n0"], 
                                       conf.level = SASCIs[j, "CONFLEV"])[2:3]))
  try(R4CIs[j,] <- unlist(gsDesign::ciBinomial(x1 = SASCIs[j,"e1"],
                                       n1 = SASCIs[j,"n1"],
                                       x2 = SASCIs[j,"e0"],
                                       n2 = SASCIs[j,"n0"], 
                                       scale = "RR",
                                       alpha = 1 - SASCIs[j, "CONFLEV"])))
 # try(R2CIs[i,] <- unlist(RRmn(SASCIs[i, c("y1", "n1", "y2", "n2")], conf.level = SASCIs[i, "CONFLEV"])[4:5]))
  #('try' function allows the loop to continue in the event of an error)
}
#warnings()

# PropCIs::riskscoreci does not apply the 'N-1' variance bias correction. 
#      See code at https://github.com/shearer/PropCIs/blob/master/R/riskscoreci.R
# sasLM::RRmn1 likewise, but also appears to be calculated with less precision. See code at https://github.com/cran/sasLM/blob/master/R/RRmn1.R
# Furthermore RRmn1 produces (NA,NA) for 214/291 vs 0/1
# DescTools::BinomRatioCI large discrepancies in LCL e.g. 81/101 vs 0/278
# Collect Sys.time() data to assess which packages are slowest

allCIs <- cbind(SASCIs, RCIs, R2CIs, R3CIs, R4CIs)
head(allCIs)
attach(allCIs)

#identify any instances of missing intervals from riskscoreci
allCIs[is.na(lclR),]
allCIs[is.na(uclR),]
lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR
ucld[U_BOUND == Inf & uclR == Inf] <- 0
lcld2 <- L_BOUND - lclR2
ucld2 <- U_BOUND - uclR2
ucld2[U_BOUND == Inf & uclR2 == Inf] <- 0
lcld3 <- L_BOUND - lclR3
ucld3 <- U_BOUND - uclR3
ucld3[U_BOUND == Inf & uclR3 == Inf] <- 0
lcld4 <- L_BOUND - lclR4
ucld4 <- U_BOUND - uclR4
ucld4[U_BOUND == Inf & uclR4 == Inf] <- 0

summary(lcld)
summary(ucld)
summary(lcld2)
summary(ucld2)
summary(lcld3)
summary(ucld3)
summary(lcld4)
summary(ucld4)


allCIs[is.na(lcld2),]
allCIs[(ucld2 > 1),]
allCIs[(lcld3 > 10),]



# hist(lcld)
# hist(ucld)

m1 <- max(abs(lcld[!is.na(lclR)])) #maximum discrepancy between SAS macro & riskscoreci version
allCIs[abs(lcld) == m1,] # 99/99 vs 50/52
m2 <- max(abs(ucld[!is.na(uclR)])) #maximum discrepancy between SAS macro & riskscoreci version
allCIs[abs(ucld) == m2,] # 230/233 vs 237/237

m3 <- max(abs(lcld2[!is.na(lclR2)])) #maximum discrepancy between SAS macro & sasLM::RRmn1 version
allCIs[abs(lcld2[!is.na(uclR2)]) == m3,] # 256/256 vs 1/259
m4 <- max(abs(ucld2[!is.na(uclR2)])) #maximum discrepancy between SAS macro & sasLM::RRmn1 version
allCIs[abs(ucld2[!is.na(uclR2)]) == m4,] # 256/256 vs 1/259
allCIs[abs(ucld2) > 1,] # 256/256 vs 1/259
#NB precision in the diffscoreci code is +/-(1e-07) so discrepancies of that order are to be expected

m5 <- max(abs(lcld3[!is.na(lclR3)])) #maximum discrepancy between SAS macro & DescTools version
allCIs[abs(lcld3) == m5,] # 
m6 <- max(abs(ucld3[!is.na(uclR3)])) #maximum discrepancy between SAS macro & DescTools version
allCIs[abs(ucld3) == m6,] # 

m7 <- max(abs(lcld4[!is.na(lclR4)])) #maximum discrepancy between SAS macro & gsDesign version
allCIs[abs(lcld4) == m7,] # 
m8 <- max(abs(ucld4[!is.na(uclR4)])) #maximum discrepancy between SAS macro & gsDesign version
allCIs[abs(ucld4) == m8,] # 


# PropCIs::riskscore bug: LCL doesn't change when conf.level changes, if x1=n1 and x2 is nearly n2
PropCIs::riskscoreci(x1=99,
                     n1=99,
                     x2=45,
                     n2=52,
                     conf.level=0.9)
# And UCL doesn't change when conf.level changes, if x2=n2
# Warning about acos function
PropCIs::riskscoreci(x1=230,
                     n1=233,
                     x2=237,
                     n2=237,
                     conf.level=0.9)


sasLM::RRmn1(y1=256,
             n1=256,
             y2=1,
             n2=259, 
             conf.level = 0.99)
PropCIs::riskscoreci(x1=256,
             n1=256,
             x2=1,
             n2=259, 
             conf.level = 0.99)
ratesci::scoreci(256,256,1,259, contrast="RR", bcf=FALSE, skew = FALSE, level=0.99)


detach(allCIs)


#Now do the same with skewness corrected intervals (also including 'N-1' bcf)
# To check against ratesci output
SASCIs2 <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2rrskew.csv")
SASCIs2 <- read.csv("C:/Mac/Home/Documents/Main/GitHub/ratesci-sas/tests/sasval2rrskew.csv")
nsamp <- dim(SASCIs2)[[1]]
head(SASCIs2)

RCIs2 <- ratesci::scoreci(x1=SASCIs2[,"e1"], 
                          n1=SASCIs2[,"n1"],
                          x2=SASCIs2[,"e0"], 
                          n2=SASCIs2[,"n0"],
                          level=SASCIs2[,"CONFLEV"], 
                          contrast = "RR",
                          skew=TRUE, 
                          precis=10)$estimates[,c(1,3)]

allCIs2 <- cbind(SASCIs2, RCIs2)
head(allCIs2)

attach(allCIs2)

rm(n1)

#identify any instances of missing intervals from diffscoreci
allCIs2[is.na(lower),]
allCIs2[is.na(upper),]
lcld <- L_BOUND - lower
ucld <- U_BOUND - upper
ucld[U_BOUND == Inf & upper == Inf] <- 0

summary(lcld)
summary(ucld)

allCIs2[is.na(lcld),]
allCIs2[ucld>0.01,]

# hist(lcld)
# hist(ucld)

m1 <- max(abs(lcld[!is.na(lower)])) #maximum discrepancy between SAS macro & scoreci version
allCIs2[abs(lcld) == m1,] # 
m2 <- max(abs(ucld[!is.na(upper)])) #maximum discrepancy between SAS macro & scoreci version
allCIs2[abs(ucld) == m2,] # mystery
# Check uncorrected interval
allCIs[abs(ucld) == m2,] # 
ratesci::scoreci(71, 94, 1, 280, contrast="RR", bcf=FALSE, skew=FALSE, level=0.99)
ratesci::scoreci(71, 94, 1, 280, contrast="RR", bcf=TRUE, skew=TRUE, level=0.99, precis=-1)

detach(allCIs2)


