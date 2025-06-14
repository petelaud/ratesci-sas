#   Validation of SAS SCORECI macro against R functions
#   PropCIs::orscoreci  (for MN interval)
#   and ratesci::scoreci  (for SCAS interval)
#

install.packages('PropCIs')
library(PropCIs)
install.packages('gsDesign')
library(gsDesign)
install.packages('sasLM')
library(sasLM)
install.packages('ratesci')
library(ratesci)

SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2or.csv")
SASCIs <- read.csv("C:/Mac/Home/Documents/Main/GitHub/ratesci-sas/tests/sasval2or.csv")
nsamp <- dim(SASCIs)[[1]]
head(SASCIs)

SASCIs[5710:5720, ]

# SASCIs$y1 <- SASCIs$e1
# SASCIs$y2 <- SASCIs$e0
# SASCIs$n2 <- SASCIs$n0

RCIs <- R2CIs <- R3CIs <- array(NA,dim=c(nsamp,2))
dimnames(RCIs)[[2]]<-c("lclR","uclR")
dimnames(R2CIs)[[2]]<-c("lclR2","uclR2")
dimnames(R3CIs)[[2]]<-c("lclR3","uclR3")

j <- 5719

for (j in 1:nsamp) {
	try(RCIs[j, ] <- unlist(PropCIs::orscoreci(x1=SASCIs[j,"e1"],
	                                            n1=SASCIs[j,"n1"],
	                                            x2=SASCIs[j,"e0"],
	                                            n2=SASCIs[j,"n0"],
	                          conf.level=SASCIs[j,"CONFLEV"])))
  try(R2CIs[j, ] <- unlist(sasLM::ORmn1(y1=SASCIs[j,"e1"],
                                       n1=SASCIs[j,"n1"],
                                       y2=SASCIs[j,"e0"],
                                       n2=SASCIs[j,"n0"], 
                                       conf.level = SASCIs[j, "CONFLEV"])[4:5]))
 # try(R2CIs[i,] <- unlist(RRmn(SASCIs[i, c("y1", "n1", "y2", "n2")], conf.level = SASCIs[i, "CONFLEV"])[4:5]))
  #('try' function allows the loop to continue in the event of an error)
}
R3CIs <- ratesci::scoreci(x1=SASCIs[,"e1"], n1=SASCIs[,"n1"],
                          x2=SASCIs[,"e0"], n2=SASCIs[,"n0"],
                          level=SASCIs[,"CONFLEV"], contrast = "OR",
                          skew=FALSE, ORbias=FALSE, precis=12)$estimates[,c(1,3)]
#warnings()

# PropCIs::orscoreci has precision to ~3dps. sasLM::ORmn1 ~4dps 
#      See code at https://github.com/shearer/PropCIs/blob/master/R/orscoreci.R
#                  https://github.com/cran/sasLM/blob/master/R/ORmn1.R
# Furthermore ORmn1 produces (NA,NA) for 214/291 vs 0/1

allCIs <- cbind(SASCIs, RCIs, R2CIs, R3CIs)
head(allCIs)
attach(allCIs)

#identify any instances of missing intervals from orscoreci
allCIs[is.na(lclR),]
allCIs[is.na(uclR),]
lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR
ucld[U_BOUND == Inf & uclR == Inf] <- 0
lcld2 <- L_BOUND - lclR2
ucld2 <- U_BOUND - uclR2
ucld2[U_BOUND == Inf & uclR2 == Inf] <- 0
lcld3 <- L_BOUND - Lower
ucld3 <- U_BOUND - Upper
ucld3[U_BOUND == Inf & Upper == Inf] <- 0

summary(lcld)
summary(ucld)
summary(lcld2)
summary(ucld2)
summary(lcld3)
summary(ucld3)


allCIs[is.na(lcld2),]



# hist(lcld)
# hist(ucld)

m1 <- max(abs(lcld[!is.na(lclR)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs[abs(lcld) == m1,] # 75/75 vs 0/109, 
ratesci::scoreci(75, 75, 0, 109, contrast = "OR", skew = FALSE, or_bias = FALSE, level=0.9, precis=12)
ratesci::scoreci(75, 75, 0, 109, contrast = "OR", skew = TRUE, or_bias = TRUE, level=0.9, precis=12)
m2 <- max(abs(ucld[!is.na(uclR)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs[abs(ucld) == m2,] # 201/202 vs 7/157
ratesci::scoreci(201, 202, 7, 157, contrast = "OR", skew = FALSE, or_bias = FALSE, level=0.95, precis=12)
m3 <- max(abs(lcld2[!is.na(lclR2)])) #maximum discrepancy between SAS macro & sasLM::ORmn1 version
allCIs[abs(lcld2[!is.na(uclR2)]) == m3,] # 201/202 vs 7/157
m4 <- max(abs(ucld2[!is.na(uclR2)])) #maximum discrepancy between SAS macro & sasLM::ORmn1 version
allCIs[abs(ucld2[!is.na(uclR2)]) == m4,] # 256/256 vs 1/259
allCIs[abs(ucld2) > 1,] 
#NB precision in the diffscoreci code is +/-(1e-03) so discrepancies of that order are to be expected

m5 <- max(abs(lcld3[!is.na(Lower)])) #maximum discrepancy between SAS macro & sasLM::ORmn1 version
allCIs[abs(lcld3[!is.na(Lower)]) == m5,] # 201/202 vs 7/157
m6 <- max(abs(ucld3[!is.na(Upper)])) #maximum discrepancy between SAS macro & sasLM::ORmn1 version
allCIs[abs(ucld3[!is.na(Lower)]) == m6,] # 201/202 vs 7/157

detach(allCIs)


# To check against ratesci output
SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2or.csv")
nsamp <- dim(SASCIs)[[1]]
head(SASCIs)

RCIs <- ratesci::scoreci(x1 = SASCIs2[,"e1"], 
                          n1 = SASCIs2[,"n1"],
                          x2 = SASCIs2[,"e0"], 
                          n2 = SASCIs2[,"n0"],
                          level = SASCIs2[,"CONFLEV"], 
                          contrast = "OR",
                          skew = FALSE,
                          or_bias = FALSE,
                          precis = 10)$estimates[,c(1,3)]

allCIs <- cbind(SASCIs, RCIs)
head(allCIs)

attach(allCIs)

rm(n1)

#identify any instances of missing intervals from diffscoreci
allCIs[is.na(lower),]
allCIs[is.na(upper),]
lcld <- L_BOUND - lower
ucld <- U_BOUND - upper
ucld[U_BOUND == Inf & upper == Inf] <- 0

summary(lcld)
summary(ucld)

m1 <- max(abs(lcld[!is.na(lower)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs[abs(lcld) == m1,] # 
m2 <- max(abs(ucld[!is.na(upper)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs[abs(ucld) == m2,] # 


detach(allCIs)


#Now do the same with skewness corrected intervals (also including OR bias correction)
# To check against ratesci output
SASCIs2 <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2orskew.csv")
nsamp <- dim(SASCIs2)[[1]]
head(SASCIs2)

RCIs2 <- ratesci::scoreci(x1 = SASCIs2[,"e1"], 
                          n1 = SASCIs2[,"n1"],
                          x2 = SASCIs2[,"e0"], 
                          n2 = SASCIs2[,"n0"],
                          level = SASCIs2[,"CONFLEV"], 
                          contrast = "OR",
                          skew = TRUE, 
                          or_bias = TRUE,
                          precis = 10)$estimates[,c(1,3)]

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


m1 <- max(abs(lcld[!is.na(lower)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs2[abs(lcld) == m1,] # 
m2 <- max(abs(ucld[!is.na(upper)])) #maximum discrepancy between SAS macro & orscoreci version
allCIs2[abs(ucld) == m2,] # 216/218 vs 1/217, can't win 'em all - small discrepancy relative to the size of UCL


allCIs2[is.na(lcld),]
allCIs2[abs(ucld) > 500,]
# 19/19 vs 0/207 - SAS macro doesnt give UCL=Inf

# hist(lcld)
# hist(ucld)

detach(allCIs2)


