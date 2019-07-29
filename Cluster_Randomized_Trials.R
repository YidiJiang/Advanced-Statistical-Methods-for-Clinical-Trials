# Cluster level analysis of binary data 

# Use standard t-test to compare the average proportions between groups
# Aggregated the binary outcome data into cluster level event proportions
# Compute the mean square error within clusters (MSW), the mean square error among clusters (MSC)
# and the intra-class correlation coefficient (ICC)

clussum <- crct_bin %>% # take a data frame
  group_by(rx,site) %>%  #drop into groupby #groupby: to create praticular structure to a group
  summarize(m=length(y),Yij=sum(y),Pij=mean(y)) %>%  #mij: # of clusters per group Yij: Pij:
  group_by(rx) %>% 
  mutate(Pi=sum(Yij)/sum(m)) 

K <- nrow(clussum)
M <- sum(clussum$m)

SS <- clussum %>%
  summarize(SSC=sum(m*(Pij-Pi)^2),SSW=sum(m*Pij*(1-Pij)),mm=sum(m^2/sum(m)),
            Pbar=mean(Pij),S2=sum((Pij-Pbar)^2),k=length(m))

(MSC <- sum(SS$SSC)/(K-2))
(MSW <- sum(SS$SSW)/(M-K))
(m0 <- (M - sum(SS$mm))/(K-2))

(rhohat <- (MSC - MSW)/(MSC + (m0-1)*MSW))

# Perform the cluster level T-test
S2 <- sum(SS$S2)/(K-2)
tu <- (SS$Pbar[1] - SS$Pbar[2])/sqrt(S2*sum(1/SS$k)) # computing a two-sided T-test
2*pt(tu,K-2)

# Use GEE (Generalized Estimating Equation) model to perform an individual level analysis
library(geepack)
# The data should be sorted by cluster before if GEE is used
# "id" is needed to specify which "i" is your data
fit.gee <- geese(y~rx,id=site,data=crct_bin,family=binomial,corstr="exch")
summary(fit.gee)




# Cluster level analysis of continuous data 

# Use standard t-test to compare the average proportions between groups
# Aggregated the binary outcome data into cluster level event proportions
# Compute the mean square error within clusters (MSW), the mean square error among clusters (MSC)
# and the intra-class correlation coefficient (ICC)

# Need the cluster means with the individual data for MSW calc
crct_cont_unbal1 <- crct_cont_unbal %>%
  group_by(trtGrp,site) %>%
  mutate(Yij=mean(y1))

clussum2 <- crct_cont_unbal %>%
  group_by(trtGrp,site) %>%
  summarize(m=length(y1),Yij=mean(y1)) %>%
  group_by(trtGrp) %>%
  mutate(Yi = sum(m*Yij)/sum(m))

K <- nrow(clussum2)
M <- sum(clussum2$m)

SSW <- crct_cont_unbal1 %>%
  group_by(trtGrp) %>%
  summarize(SSW = sum((y1-Yij)^2))

SSC <- clussum2 %>%
  summarize(SSC = sum(m*(Yij-Yi)^2),mm=sum(m^2/sum(m)),SS=sum((Yij-Yi)^2),
            k=length(m),Ybar=mean(Yi))

(MSC <- sum(SSC$SSC)/(K-2))
(MSW <- sum(SSW$SSW)/(M-K))
(m0 <- (M - sum(SSC$mm))/(K-2))

(rhohat <- (MSC - MSW)/(MSC + (m0-1)*MSW))

# Cluster level
S2 <- sum(SSC$SS)/(K-2)
delta <- with(SSC, Ybar[2]-Ybar[1])
(tu <- delta/sqrt(S2*sum(1/SSC$k)))
2*pt(tu,K-2,lower.tail=FALSE)

# Use linear mixed effect model to conduct individual level analysis of cluster data with continuous outcome
library(lme4)
fit.lme <- lmer(y1~trtGrp + (1|site), data=crct_cont_unbal)
summary(fit.lme)
anova(fit.lme)
# Extracting the beta coefficients and CI is easy enough
fixef(fit.lme)
confint(fit.lme)
# Extract the variance components to obtain the ICC
(vcomp <- VarCorr(fit.lme))
print(vcomp,comp=c("Variance"))
vc <- as.data.frame(vcomp)[,"vcov"]
(icc.lme <- vc[1]/sum(vc))

# Returning to the matter of p-values, the lmerTest package provides various options.
# In this case there are no differences because the data are balanced.
library(lmerTest)
fit.lme <- update(fit.lme)
anova(fit.lme)
anova(fit.lme,ddf="Kenward-Roger")
anova(fit.lme,ddf="lme4")
summary(fit.lme)
summary(fit.lme,ddf="Kenward-Roger")
summary(fit.lme,ddf="lme4")
fixef(fit.lme)
confint(fit.lme)



# Rearrange the sample size formula for a continuous outcome
# For a hypothetical cluster trial scenario, the power at various cluster size and cluster number values can be calculated and plotted
mcid <- 0
simga_parameter <- 0
type1_zscore <- 1.96

type2_zscore <- function(k, m, rho, mcid, sigma, type1_zscore) {
  z <- sqrt((k*m*mcid^2) / (2*sigma^2*(1+rho*(m-1)))) - type1_zscore
  return(z)
}

k <- seq(2,100,1)
m <- seq(1,100,1)
rho <- 0

data <- matrix(NA, nrow=length(k)*length(m)*length(rho), ncol=5)
colnames(data) <- c("k","m","rho","type2_zscore","power")

row_num <- 1
for(i in k) {
  for(j in m) {
    for(l in rho) {
      
      data[row_num, "k"] <- i
      data[row_num, "m"] <- j
      data[row_num, "rho"] <- l
      data[row_num, "type2_zscore"] <- type2_zscore(k=i, m=j, rho=l,
                                                    mcid=mcid, sigma=simga_parameter, type1_zscore=type1_zscore)
      data[row_num, "power"] <- 1 - pnorm(data[row_num, "type2_zscore"])
      
    }
  }
}