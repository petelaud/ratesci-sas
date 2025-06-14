#   Validation of SAS SCORECI macro against R function
#    ratesci::scoreci  (for stratified MN and SCAS intervals)
#

install.packages('ratesci')
library(ratesci)
install.packages("pairwiseCI")
library(pairwiseCI)

SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3rr.csv")
nsamp <- max(SASCIs$sample)
head(SASCIs)

RCIs <- array(NA, dim = c(nsamp, 2))
dimnames(RCIs)[[2]] <- c("lclR", "uclR")
i <- 1
for (i in 1:nsamp) {
  onesample <- SASCIs[SASCIs$sample == i, ]
  onesample$y1 <- onesample$e1
  onesample$y2 <- onesample$e0
  onesample$n2 <- onesample$n0
  
  try(RCIs[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"], contrast = "RR",
                                    level = onesample[,"CONFLEV"][1], stratified = TRUE,
                                    skew = FALSE, precis=10, weighting = "MH")$estimates[,c(1,3)])
#  try(R2CIs[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
#                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"], contrast = "RR",
#                                    level = onesample[,"CONFLEV"][1], stratified = TRUE, bcf = FALSE,
#                                    skew = FALSE, precis=10, weighting = "MN")$estimates[,c(1,3)])
  
  # RRmn function cannot cope with zero event counts!
#  sasLM::RRmn(onesample[, c("y1", "n1", "y2", "n2")], conf.level = onesample[1, "CONFLEV"]) #$Common[4:5]
  
  
  #('try' function allows the loop to continue in the event of an error)
}

allCIs <- cbind(SASCIs[SASCIs$stratum == 1,], RCIs)
head(allCIs)
attach(allCIs)

# Summarise differences between SAS and R
lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

summary(lcld)
summary(ucld)

# hist(lcld)
# hist(ucld)

max(abs(lcld)) #maximum discrepancy between SAS macro & ratesci version
max(abs(ucld)) #maximum discrepancy between SAS macro & ratesci version

detach(allCIs)


# Check example dataset from Lee & Bae 2022, Table 3
# Some issues with this paper/function:
# Claimed "The MN score method has not been previously implemented in R software for data with stratification"
# when the ratesci package had been available since 2016 (MN weighting refined in Dec2021)
# Output estimates for RR are inconsistent with estimates of p1 and p2
# Output estimates for RR are p2/p1 if p2>p1
# zero cell counts produce an error
# 'N-1' correction has been omitted from code since the paper was published
# - this seems to be the reason for the discrepancy vs PropCIs functions noted in the paper
#   (the authors of PropCIs have included the BCF for RD only)
#   but even comparing like-for-like excluding BCF, RRmn1 still manages to mismatch slightly
# "Two different sets of functions were necessary" for with and without stratification
# Example from Bernal et al has UCL closer to zero when including BCF and SKEW
# Point out stratified calculation not available in SAS PROC FREQ, but macro on GitHub

d1 = data.frame(matrix(c(25, 339, 28, 335, 23, 370, 40, 364), nrow=2, byrow=TRUE))
colnames(d1) =  c("y1", "n1", "y2", "n2")
d1out <- RDmn(d1)
d1out$Common[1] - d1out$Common[2]
d1out2 <- scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RD", skew=FALSE, bcf=FALSE, weighting="MN", stratified=TRUE)
d1out2a <- scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RD", skew=TRUE, bcf=TRUE, weighting="MH", stratified=TRUE)
d1out <- RRmn(d1)
d1out$Common[1]/d1out$Common[2]
d1out2 <- scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RR", skew=FALSE, bcf=FALSE, weighting="MN", stratified=TRUE)
d1out2$estimates[,7]/d1out2$estimates[,8]


RDmn(d1)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RD", skew=FALSE, bcf=FALSE, weighting="MN", stratified=TRUE)

ORmn(d1)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="OR", skew=FALSE, bcf=FALSE, weighting="MN", stratified=TRUE)

# Figure 2
library(sasLM)
RDmn1(y1=28, n1=385, y2=53, n2=377)
scoreci(28, 385, 53, 377, contrast="RD", skew=FALSE, bcf=FALSE, precis=10)
PropCIs::diffscoreci(28, 385, 53, 377, 0.95)
RRmn1(y1=28, n1=385, y2=53, n2=377)
scoreci(28, 385, 53, 377, contrast="RR", skew=FALSE, bcf=FALSE, precis=10)
pairwiseCI::Prop.ratio(c(28,385-28), c(53,377-53), CImethod="Score")
scoreci(28, 385, 53, 377, contrast="RR", skew=FALSE, bcf=TRUE, precis=10)
pairwiseCI::Prop.ratio(c(28,385-28), c(53,377-53), CImethod="MNScore")
PropCIs::riskscoreci(28, 385, 53, 377, 0.95)
ORmn1(y1=28, n1=385, y2=53, n2=377)
#Try pairwiseci functions for good measure?


# Check example from Gart 1988. scoreci doesn't quite match, but RRmn is way off.
d1 = data.frame(matrix(c(4,16,5,79, 2,16,3,87, 4,18,10,90, 1,15,3,82), nrow=4, byrow=TRUE))
colnames(d1) =  c("y1", "n1", "y2", "n2")
RRmn(d1)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RR", skew=FALSE, bcf=FALSE, weighting="MN", stratified=TRUE)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RR", skew=TRUE, bcf=TRUE, weighting="MN", stratified=TRUE)


# Check example from Gart 1990. scoreci doesn't quite match, but RRmn is further off.
d1 = data.frame(matrix(c(34,46,1,9, 30,38,2,10), nrow=2, byrow=TRUE))
colnames(d1) =  c("y1", "n1", "y2", "n2")
RDmn(d1)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RD", skew=FALSE, bcf=FALSE, weighting="MH", stratified=TRUE)
scoreci(d1$y1, d1$n1, d1$y2, d1$n2, contrast="RD", skew=TRUE, bcf=FALSE, weighting="MH", stratified=TRUE)



# Now do the same with skewness corrected intervals
SASCIs2<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3rrskew.csv")
nsamp <- max(SASCIs2$sample)
head(SASCIs2)

RCIs2 <- array(NA, dim = c(nsamp, 2))
dimnames(RCIs2)[[2]] <- c("lclR", "uclR")
for (i in 1:nsamp) {
  onesample <- SASCIs2[SASCIs2$sample == i, ]
  try(RCIs2[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"], contrast="RR",
                                    level = onesample[,"CONFLEV"][1], stratified = TRUE,
                                    skew = TRUE, precis=10)$estimates[,c(1,3)])
  #('try' function allows the loop to continue in the event of an error)
}

allCIs2 <- cbind(SASCIs2[SASCIs2$stratum == 1,], RCIs2)
head(allCIs2)

attach(allCIs2)

lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

summary(lcld)
summary(ucld)

# hist(lcld)
# hist(ucld)

max(abs(lcld)) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld)) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs2)


# Now do the same with Poisson intervals
SASCIsp<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3prr.csv")
nsamp <- max(SASCIsp$sample)
head(SASCIsp)

RCIsp <- array(NA, dim = c(nsamp, 2))
dimnames(RCIsp)[[2]] <- c("lclR", "uclR")
for (i in 1:nsamp) {
  onesample <- SASCIsp[SASCIsp$sample == i, ]
  try(RCIsp[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"], contrast="RR",
                                    level = onesample[,"CONFLEV"][1], stratified = TRUE,
                                    skew = FALSE, distrib = 'poi', precis=10)$estimates[,c(1,3)])
  #('try' function allows the loop to continue in the event of an error)
}

allCIsp <- cbind(SASCIsp[SASCIsp$stratum == 1,], RCIsp)
head(allCIsp)

attach(allCIsp)

lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

max(abs(lcld)) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld)) #maximum discrepancy between SAS macro & scoreci version

detach(allCIsp)


# Now do the same with skewness corrected intervals
SASCIsp2<-read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval3prrskew.csv")
nsamp <- max(SASCIsp2$sample)
head(SASCIsp2)

RCIsp2 <- array(NA, dim = c(nsamp, 2))
dimnames(RCIsp2)[[2]] <- c("lclR", "uclR")
for (i in 1:nsamp) {
  onesample <- SASCIsp2[SASCIsp2$sample == i, ]
  try(RCIsp2[i, ] <- ratesci::scoreci(x1 = onesample[,"e1"], n1 = onesample[,"n1"],
                                    x2 = onesample[,"e0"], n2 = onesample[,"n0"], contrast="RR",
                                    level = onesample[,"CONFLEV"][1], stratified = TRUE,
                                    skew = TRUE, distrib = 'poi', precis=10)$estimates[,c(1,3)])
  #('try' function allows the loop to continue in the event of an error)
}

allCIsp2 <- cbind(SASCIsp2[SASCIsp2$stratum == 1,], RCIsp2)
head(allCIsp2)

attach(allCIsp2)

lcld <- L_BOUND - lclR
ucld <- U_BOUND - uclR

max(abs(lcld)) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld)) #maximum discrepancy between SAS macro & scoreci version

detach(allCIsp2)



