# This project is to plan the sample size of the trial by using simulation.

# Before simulating the power, verify the analysis method is valid
# Compare three analysis methods: Naive model, GEE model and Random Intercept model (for bias, standard 
# deviation and 95% confidence interval coverage probability)
library(geepack)
library(nlme)
n <- 100
m <- 5000
theta <- -5
pointest <- matrix(NA, m, 3)
varest <- matrix(NA, m, 3)
coverage <- matrix(NA, m, 3)
power <- matrix(NA, m, 3)
for (j in 1:m) {
  if(j %% 100 == 0) {print(j)}
  # Generate dataset:
  z <- rbinom(n, 1, 0.5)
  a <- rnorm(n, mean=50, sd=sqrt(50))
  e1 <- rnorm(n, sd=sqrt(50))
  e2 <- rnorm(n, sd=sqrt(50))
  e3 <- rnorm(n, sd=sqrt(50))
  e4 <- rnorm(n, sd=sqrt(50))
  y1 <- a + e1
  y2 <- a + theta*z + e2
  y3 <- a + theta*z + e3
  y4 <- a + theta*z + e4
  id <- 1:n
  simulation_data <- data.frame("id"=rep(id, 3),
                                "z"=rep(z, 3),
                                "y1"=rep(y1, 3),
                                "y"=c(y2, y3, y4),
                                "time"=rep(c(2,4,6), each=n))
  simulation_data <- simulation_data[order(simulation_data$id),]
  # Naive model:
  naive_model <- lm(y~z, data=simulation_data)
  pointest[j,1] <- coef(naive_model)[2]
  varest[j,1] <- vcov(naive_model)[2,2]
  cil <- pointest[j,1] + qnorm(0.025) * sqrt(varest[j,1])
  ciu <- pointest[j,1] + qnorm(0.975) * sqrt(varest[j,1])
  coverage[j,1] <- (theta >= cil) & (theta <= ciu)
  power[j,1] <- (cil > 0.0) | (ciu < 0.0)
  # GEE:
  gee_model <- geeglm(y~z, corstr="independence", std.err="san.se", id=id, data=simulation_data)
  pointest[j,2] <- coef(gee_model)[2]
  varest[j,2] <- gee_model$geese$vbeta[2,2]
  cil <- pointest[j,2] + qnorm(0.025) * sqrt(varest[j,2])
  ciu <- pointest[j,2] + qnorm(0.975) * sqrt(varest[j,2])
  coverage[j,2] <- (theta >= cil) & (theta <= ciu)
  power[j,2] <- (cil > 0.0) | (ciu < 0.0)
  # Random intercept:
  random_inter_model <- lme(y~z, random=~1|id, data=simulation_data)
  pointest[j,3] <- random_inter_model$coefficients$fixed[2]
  varest[j,3] <- vcov(random_inter_model)[2,2]
  cil <- pointest[j,3] + qnorm(0.025) * sqrt(varest[j,3])
  ciu <- pointest[j,3] + qnorm(0.975) * sqrt(varest[j,3])
  coverage[j,3] <- (theta >= cil) & (theta <= ciu)
  power[j,3] <- (cil > 0.0) | (ciu < 0.0)
}
results <- cbind(colMeans(pointest), colMeans(pointest) - theta, apply(pointest, 2, sd),
                 sqrt(colMeans(varest)),
                 100*(apply(pointest, 2, var) + (colMeans(pointest) - theta)^2), sqrt(colMeans(varest)/m),
                 colMeans(coverage), colMeans(power))
rownames(results) <- c('Naive', 'GEE', 'RandomIntercept')
colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', 'MCE', 'Coverage', 'Power')
round(results, 3)

# add the baseline volume measurement influence as a covariate to three models to check how baseline volume
# measurement influences the results
n <- 100
m <- 5000
theta <- 5
pointest <- matrix(NA, m, 3)
varest <- matrix(NA, m, 3)
coverage <- matrix(NA, m, 3)
power <- matrix(NA, m, 3)
for (j in 1:m) {
  if(j %% 100 == 0) {print(j)}
  # Generate dataset:
  z <- rbinom(n, 1, 0.5)
  a <- rnorm(n, mean=50, sd=sqrt(50))
  e1 <- rnorm(n, sd=sqrt(50))
  e2 <- rnorm(n, sd=sqrt(50))
  e3 <- rnorm(n, sd=sqrt(50))
  e4 <- rnorm(n, sd=sqrt(50))
  y1 <- a + e1
  y2 <- a + theta*z + e2
  y3 <- a + theta*z + e3
  y4 <- a + theta*z + e4
  id <- 1:n
  simulation_data <- data.frame("id"=rep(id, 3),
                                "z"=rep(z, 3),
                                "y1"=rep(y1, 3),
                                "y"=c(y2, y3, y4),
                                "time"=rep(c(2,4,6), each=n))
  simulation_data <- simulation_data[order(simulation_data$id),]
  # Naive model:
  naive_model <- lm(y~z+y1, data=simulation_data)
  pointest[j,1] <- coef(naive_model)[2]
  varest[j,1] <- vcov(naive_model)[2,2]
  cil <- pointest[j,1] + qnorm(0.025) * sqrt(varest[j,1])
  ciu <- pointest[j,1] + qnorm(0.975) * sqrt(varest[j,1])
  coverage[j,1] <- (theta >= cil) & (theta <= ciu)
  power[j,1] <- (cil > 0.0) | (ciu < 0.0)
  # GEE:
  gee_model <- geeglm(y~z+y1, corstr="independence", std.err="san.se", id=id,
                      data=simulation_data)
  pointest[j,2] <- coef(gee_model)[2]
  varest[j,2] <- gee_model$geese$vbeta[2,2]
  cil <- pointest[j,2] + qnorm(0.025) * sqrt(varest[j,2])
  ciu <- pointest[j,2] + qnorm(0.975) * sqrt(varest[j,2])
  coverage[j,2] <- (theta >= cil) & (theta <= ciu)
  power[j,2] <- (cil > 0.0) | (ciu < 0.0)
  # Random intercept:
  random_inter_model <- lme(y~z+y1, random=~1|id, data=simulation_data)
  pointest[j,3] <- random_inter_model$coefficients$fixed[2]
  varest[j,3] <- vcov(random_inter_model)[2,2]
  cil <- pointest[j,3] + qnorm(0.025) * sqrt(varest[j,3])
  ciu <- pointest[j,3] + qnorm(0.975) * sqrt(varest[j,3])
  coverage[j,3] <- (theta >= cil) & (theta <= ciu)
  power[j,3] <- (cil > 0.0) | (ciu < 0.0)
}
results <- cbind(colMeans(pointest), colMeans(pointest) - theta, apply(pointest, 2, sd),
                 sqrt(colMeans(varest)),
                 100*(apply(pointest, 2, var) + (colMeans(pointest) - theta)^2), sqrt(colMeans(varest)/m),
                 colMeans(coverage), colMeans(power))
rownames(results) <- c('Naive', 'GEE', 'RandomIntercept')
colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', 'MCE', 'Coverage', 'Power')
round(results, 3)

# For random intercept model without adjustment for baseline volume, the power curve to reject the null is
sim <- function(n, m, theta) {
  pointest <- matrix(NA, m, 3)
  varest <- matrix(NA, m, 3)
  coverage <- matrix(NA, m, 3)
  power <- matrix(NA, m, 3)
  for (j in 1:m) {
    if(j %% 100 == 0) {print(j)}
    # Generate dataset:
    z <- rbinom(n, 1, 0.5)
    a <- rnorm(n, mean=50, sd=sqrt(50))
    e1 <- rnorm(n, sd=sqrt(50))
    e2 <- rnorm(n, sd=sqrt(50))
    e3 <- rnorm(n, sd=sqrt(50))
    e4 <- rnorm(n, sd=sqrt(50))
    y1 <- a + e1
    y2 <- a + theta*z + e2
    y3 <- a + theta*z + e3
    y4 <- a + theta*z + e4
    id <- 1:n
    simulation_data <- data.frame("id"=rep(id, 3),
                                  "z"=rep(z, 3),
                                  "y1"=rep(y1, 3),
                                  "y"=c(y2, y3, y4),
                                  "time"=rep(c(2,4,6), each=n))
    simulation_data <- simulation_data[order(simulation_data$id),]
    # Random intercept (without adjustment for baseline):
    random_inter_model <- lme(y~z, random=~1|id, data=simulation_data)
    pointest[j,3] <- random_inter_model$coefficients$fixed[2]
    varest[j,3] <- vcov(random_inter_model)[2,2]
    cil <- pointest[j,3] + qnorm(0.025) * sqrt(varest[j,3])
    ciu <- pointest[j,3] + qnorm(0.975) * sqrt(varest[j,3])
    coverage[j,3] <- (theta >= cil) & (theta <= ciu)
    power[j,3] <- (cil > 0.0) | (ciu < 0.0)
  }
  results <- cbind(colMeans(pointest), colMeans(pointest) - theta, apply(pointest, 2, sd),
                   sqrt(colMeans(varest)),
                   100*(apply(pointest, 2, var) + (colMeans(pointest) - theta)^2), sqrt(colMeans(varest)/
                                                                                          m), colMeans(coverage), colMeans(power))
  rownames(results) <- c('Naive', 'GEE', 'RandomIntercept')
  colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', 'MCE', 'Coverage', 'Power')
  return(round(results, 3))
}
# Estimate the power to reject the null at several different sample sizes
# Use sample sizes from 100 to 200 and the jump size equals 5
sample_size <- c(seq(100,200,5))
sample_power <- matrix(NA, 21, 2)
sample_power[,1] <- sample_size
colnames(sample_power) <- c("n", "power")
for(i in 1:nrow(sample_power)) {
  sample_power[i,"power"] <- sim(n=sample_power[i,"n"], m=1000, theta=-5)
  ["RandomIntercept","Power"]
}
plot(x=sample_power[,"n"], y=sample_power[,"power"], xlab="Sample Size", ylab="Power",
     type="b", lwd=1.5, lty=3)
# From the power curve, the required minimum sample size needed for 90% power is around 108


# Simulate random dropout with probability of 10% before each follow-up measurement
sim <- function(n, m, theta) {
  pointest <- matrix(NA, m, 3)
  varest <- matrix(NA, m, 3)
  coverage <- matrix(NA, m, 3)
  power <- matrix(NA, m, 3)
  for (j in 1:m) {
    if(j %% 100 == 0) {print(j)}
    # Generate dataset:
    z <- rbinom(n, 1, 0.5)
    a <- rnorm(n, mean=50, sd=sqrt(50))
    e1 <- rnorm(n, sd=sqrt(50))
    e2 <- rnorm(n, sd=sqrt(50))
    e3 <- rnorm(n, sd=sqrt(50))
    e4 <- rnorm(n, sd=sqrt(50))
    y1 <- a + e1
    y2 <- a + theta*z + e2
    y3 <- a + theta*z + e3
    y4 <- a + theta*z + e4
    id <- 1:n
    simulation_data <- data.frame("id"=rep(id, 3),
                                  "z"=rep(z, 3),
                                  "y1"=rep(y1, 3),
                                  "y"=c(y2, y3, y4),
                                  "dropout"=rep(0, n*3),
                                  "time"=rep(c(2,4,6), each=n))
    simulation_data <- simulation_data[order(simulation_data$id),]
    # Simulate dropout:
    id1 <- id
    dropout1 <- id1*rbinom(length(id1), 1, 0.1)
    dropout1 <- dropout1[dropout1!=0]
    id2 <- setdiff(id,dropout1)
    dropout2 <- id2*rbinom(length(id2), 1, 0.1)
    dropout2 <- dropout2[dropout2!=0]
    id3 <- setdiff(id,c(dropout1,dropout2))
    dropout3 <- id3*rbinom(length(id3), 1, 0.1)
    dropout3 <- dropout3[dropout3!=0]
    for(k in dropout1) {
      simulation_data$dropout[simulation_data$id==k] <- c(1,1,1)
    }
    for(k in dropout2) {
      simulation_data$dropout[simulation_data$id==k] <- c(0,1,1)
    }
    for(k in dropout3) {
      simulation_data$dropout[simulation_data$id==k] <- c(0,0,1)
    }
    simulation_data <- simulation_data[simulation_data$dropout==0,]
    # Random intercept:
    random_inter_model <- lme(y~z, random=~1|id, data=simulation_data)
    pointest[j,3] <- random_inter_model$coefficients$fixed[2]
    varest[j,3] <- vcov(random_inter_model)[2,2]
    cil <- pointest[j,3] + qnorm(0.025) * sqrt(varest[j,3])
    ciu <- pointest[j,3] + qnorm(0.975) * sqrt(varest[j,3])
    coverage[j,3] <- (theta >= cil) & (theta <= ciu)
    power[j,3] <- (cil > 0.0) | (ciu < 0.0)
  }
  results <- cbind(colMeans(pointest), colMeans(pointest) - theta, apply(pointest, 2, sd),
                   sqrt(colMeans(varest)),
                   100*(apply(pointest, 2, var) + (colMeans(pointest) - theta)^2), sqrt(colMeans(varest)/
                                                                                          m), colMeans(coverage), colMeans(power))
  rownames(results) <- c('Naive', 'GEE', 'RandomIntercept')
  colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', 'MCE', 'Coverage', 'Power')
  return(round(results, 3))
}
# Re-estimate the power to reject the null at several different sample sizes 
# Use sample sizes from 100 to 200 and the jump size equals 5;
sample_size <- c(seq(100,200,5))
sample_power <- matrix(NA, 21, 2)
sample_power[,1] <- sample_size
colnames(sample_power) <- c("n", "power")
for(i in 1:nrow(sample_power)) {
  sample_power[i,"power"] <- sim(n=sample_power[i,"n"], m=1000, theta=-5)
  ["RandomIntercept","Power"]
}
plot(x=sample_power[,"n"], y=sample_power[,"power"], xlab="Sample Size", ylab="Power",
     type="b", lwd=1.5, lty=3)
# From the power curve, the required minimum sample size needed for 90% power is around 130


# Simulate random treatment discontinuation in intervention arm with probability of 10% before each follow-up measurement
sim <- function(n, m, theta) {
  pointest <- matrix(NA, m, 3)
  varest <- matrix(NA, m, 3)
  coverage <- matrix(NA, m, 3)
  power <- matrix(NA, m, 3)
  for (j in 1:m) {
    if(j %% 100 == 0) {print(j)}
    # Generate dataset:
    z <- rbinom(n, 1, 0.5)
    a <- rnorm(n, mean=50, sd=sqrt(50))
    e1 <- rnorm(n, sd=sqrt(50))
    e2 <- rnorm(n, sd=sqrt(50))
    e3 <- rnorm(n, sd=sqrt(50))
    e4 <- rnorm(n, sd=sqrt(50))
    y1 <- a + e1
    y2 <- a + theta*z + e2
    y3 <- a + theta*z + e3
    y4 <- a + theta*z + e4
    alt_y2 <- a + theta*(1-z) + e2
    alt_y3 <- a + theta*(1-z) + e3
    alt_y4 <- a + theta*(1-z) + e4
    id <- 1:n
    simulation_data <- data.frame("id"=rep(id, 3),
                                  "y1"=rep(y1, 3),
                                  "z"=rep(z, 3),
                                  "y"=c(y2, y3, y4),
                                  "alt_z"=rep(1-z, 3),
                                  "alt_y"=c(alt_y2, alt_y3, alt_y4),
                                  "discontinue"=rep(0, n*3),
                                  "time"=rep(c(2,4,6), each=n),
                                  "final_z"=rep(NA, n*3),
                                  "final_y"=rep(NA, n*3))
    simulation_data <- simulation_data[order(simulation_data$id),]
    # Simulate treatment discontinuation:
    treatment_id <- unique(simulation_data$id[simulation_data$z==1])
    id1 <- treatment_id
    discontinue1 <- id1*rbinom(length(id1), 1, 0.1)
    discontinue1 <- discontinue1[discontinue1!=0]
    id2 <- setdiff(treatment_id,discontinue1)
    discontinue2 <- id2*rbinom(length(id2), 1, 0.1)
    discontinue2 <- discontinue2[discontinue2!=0]
    id3 <- setdiff(treatment_id,c(discontinue1,discontinue2))
    discontinue3 <- id3*rbinom(length(id3), 1, 0.1)
    discontinue3 <- discontinue3[discontinue3!=0]
    for(k in discontinue1) {
      simulation_data$discontinue[simulation_data$id==k] <- c(1,1,1)
    }
    for(k in discontinue2) {
      simulation_data$discontinue[simulation_data$id==k] <- c(0,1,1)
    }
    for(k in discontinue3) {
      simulation_data$discontinue[simulation_data$id==k] <- c(0,0,1)
    }
    for(l in 1:nrow(simulation_data)) {
      if(simulation_data$discontinue[l] == 1){
        simulation_data$final_z[l] <- simulation_data$alt_z[l]
        simulation_data$final_y[l] <- simulation_data$alt_y[l]
      } else {
        simulation_data$final_z[l] <- simulation_data$z[l]
        simulation_data$final_y[l] <- simulation_data$y[l]
      }
    }
    # Random intercept:
    random_inter_model <- lme(final_y~final_z, random=~1|id, data=simulation_data)
    pointest[j,3] <- random_inter_model$coefficients$fixed[2]
    varest[j,3] <- vcov(random_inter_model)[2,2]
    cil <- pointest[j,3] + qnorm(0.025) * sqrt(varest[j,3])
    ciu <- pointest[j,3] + qnorm(0.975) * sqrt(varest[j,3])
    coverage[j,3] <- (theta >= cil) & (theta <= ciu)
    power[j,3] <- (cil > 0.0) | (ciu < 0.0)
  }
  results <- cbind(colMeans(pointest), colMeans(pointest) - theta, apply(pointest, 2, sd),
                   sqrt(colMeans(varest)),
                   100*(apply(pointest, 2, var) + (colMeans(pointest) - theta)^2), sqrt(colMeans(varest)/
                                                                                          m), colMeans(coverage), colMeans(power))
  rownames(results) <- c('Naive', 'GEE', 'RandomIntercept')
  colnames(results) <- c('Mean', 'Bias', 'SD', 'Mean SE', '100xMSE', 'MCE', 'Coverage', 'Power')
  return(round(results, 3))
}
# Re-estimate the power to reject the null at several different sample sizes
# Use sample sizes from 100 to 200 and the jump size equals 5
sample_size <- c(seq(100,200,5))
sample_power <- matrix(NA, 21, 2)
sample_power[,1] <- sample_size
colnames(sample_power) <- c("n", "power")
for(i in 1:nrow(sample_power)) {
  sample_power[i,"power"] <- sim(n=sample_power[i,"n"], m=1000, theta=-5)
  ["RandomIntercept","Power"]
}
plot(x=sample_power[,"n"], y=sample_power[,"power"], xlab="Sample Size", ylab="Power",
     type="b", lwd=1.5, lty=3)
# From the power curve, the required minimum sample size needed for 90% power is around 132