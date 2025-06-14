#   Validation of SAS SCORECI macro against R functions
#   ratesci::scoreci  (for MN & SCAS interval for Poisson RR)
#
#   SAS sometimes fails to converge for large UCL beyond ~2dps for 99% CI - why?

install.packages('ratesci')
library(ratesci)

install.packages("pak")
pak::pak("petelaud/ratesci")

SASCIs <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2prr.csv")
#SASCIs <- read.csv("C:/Mac/Home/Documents/Main/GitHub/ratesci-sas/tests/sasval2prr.csv")
nsamp <- dim(SASCIs)[[1]]
head(SASCIs)

RCIs <- ratesci::scoreci(x1=SASCIs[,"e1"], n1=SASCIs[,"n1"],
                          x2=SASCIs[,"e0"], n2=SASCIs[,"n0"],
                          level=SASCIs[,"CONFLEV"], contrast = 'RR',
                          skew=FALSE, precis=10, distrib='poi')$estimates[,c(1,3)]

allCIs <- cbind(SASCIs,RCIs)
head(allCIs)
attach(allCIs)

#identify any instances of missing intervals from ratesci
allCIs[is.na(lower),]
allCIs[is.na(upper),]
lcld <- L_BOUND - lower
ucld <- U_BOUND - upper
ucld[U_BOUND == Inf & upper == Inf] <- 0

summary(lcld)
summary(ucld)


allCIs[is.na(lcld),]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(lower)])) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld[!is.na(upper)])) #maximum discrepancy between SAS macro & scoreci version
allCIs[(abs(ucld) > 10),]

detach(allCIs)


#Now do the same with skewness corrected intervals
SASCIs2 <- read.csv("/Users/ssu/Documents/Main/GitHub/ratesci-sas/tests/sasval2prrskew.csv")
SASCIs2 <- read.csv("C:/Mac/Home/Documents/Main/GitHub/ratesci-sas/tests/sasval2prrskew.csv")
nsamp <- dim(SASCIs2)[[1]]
head(SASCIs2)

RCIs2 <- ratesci::scoreci(x1=SASCIs2[,"e1"], n1=SASCIs2[,"n1"],
                          x2=SASCIs2[,"e0"], n2=SASCIs2[,"n0"],
                          level=SASCIs2[,"CONFLEV"], contrast = "RR",
                          skew=TRUE, precis=12, distrib='poi')$estimates[,c(1,3)]

allCIs2 <- cbind(SASCIs2, RCIs2)
head(allCIs2)

attach(allCIs2)

#identify any instances of missing intervals from ratesci
allCIs2[is.na(lower),]
allCIs2[is.na(upper),]
lcld <- L_BOUND - lower
ucld <- U_BOUND - upper
ucld[U_BOUND == Inf & upper == Inf] <- 0

summary(lcld)
summary(ucld)
summary(ucld[CONFLEV < 0.99])

allCIs2[is.na(lcld),]
allCIs2[abs(ucld)>0.01,]

# hist(lcld)
# hist(ucld)

max(abs(lcld[!is.na(Lower)])) #maximum discrepancy between SAS macro & scoreci version
max(abs(ucld[!is.na(Upper)])) #maximum discrepancy between SAS macro & scoreci version

detach(allCIs2)


