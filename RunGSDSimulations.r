# This script runs the various simulations used in the ESD paper by Eaton, Moore and MacKenzie
# It requires that the file GSD_Data.RData be loaded, and it creates (and WILL OVERWRITE)
# the file GSD_simulations.RData. The original version of the simulation data on GitHub produces
# the exact results in the paper, whereas re-running the simulation will produce similar but
# not identical results.

### run the monte carlo simulation #################################################
k = 100  #set the sample size
m = 5000  #set the number of Monte Carlo simulations to run
qntls = c(0.1, 0.16, 0.25, 0.5, 0.75, 0.84, 0.90)

ci = WolmanCI(cfd1all, k, 
             probs = qntls,
             alpha = 0.05,
             Np = length(all1.data$D))

ci50 = WolmanCI(cfd1all, k, 
              probs = qntls,
              alpha = 0.5,
              Np = length(all1.data$D))

results = matrix(data = NA, ncol = length(ci$probs), nrow = m)

for(i in seq_along(results[,1])){
  sample.data = sample(all1.data$D, replace = T, size = k)
  cfdtmp = MakeCFD(sample.data)
  results[i,] = WolmanCI(cfdtmp, k, probs = qntls)[[2]]
}

mclimits = matrix(data = NA, ncol = 2, nrow = length(ci$probs))
for(i in seq_along(mclimits[,1])){
  mclimits[i,] = as.numeric(quantile(results[,i], probs = c(0.005, 0.995)))
}

### perform the simulation of log normal error #####################################


NSim = 3000
SErr84 = array(data = NA, dim =  c(NSim,4))
SErr50 = array(data = NA, dim =  c(NSim,4))

for(i in seq(1,NSim,1)){
  log50 = runif(1, min = 4.5, max = 6.5)  #vary d50 from 22.6 mm to 90.5 mm
  logSD = runif(1, min = 0.5, 2.0)  #vary D16 to D84 range from factor of 2 to factor of 16
  #logSD = 0.5
  n = round(runif(1, min = 50, max = 1000))  #generate sample sizes from 100 to 500
  obs = 2^rnorm(n,log50,logSD)  #take n samples from a log normal distribution
  gsd = MakeCFD(obs, plot = F)  #extract the cumulative distribution from the sample
  
  #estimate the error for the D84 from the distribution
  p = c(0.50,0.84)
  tmp = WolmanCI(gsd,n,p)  
  err = 0.5*(tmp$upper - tmp$lower)/tmp$estimate
  SErr50[i,] = c(n,logSD,log50,err[1])
  SErr84[i,] = c(n,logSD, log50, err[2])
}

SErr84 = as.data.frame(SErr84)
colnames(SErr84) = c("N", "logSD", "log50", "err")
SErr50 = as.data.frame(SErr50)
colnames(SErr50) = c("N", "logSD", "log50", "err")

powfit50 = lm(log(SErr50$err) ~ log(SErr50$N) + SErr50$logSD)
rsq50 = as.numeric(summary(powfit50)[9])
tmp = as.matrix(anova(powfit50)[2])
var50 = as.numeric(tmp[1:3]/sum(tmp[1:3]))
coeff50 = as.numeric(powfit50$coefficients)

powfit84 = lm(log(SErr84$err) ~ log(SErr84$N) + SErr84$logSD)
rsq84 = as.numeric(summary(powfit84)[9])
tmp = as.matrix(anova(powfit84)[2])
var84 = as.numeric(tmp[1:3]/sum(tmp[1:3]))
coeff84 = as.numeric(powfit84$coefficients)


### clean up the workspace #########################################################

save(m, mclimits, results, ci, ci50,
     NSim, SErr84, SErr50, 
     powfit50, rsq50, var50, coeff50,
     powfit84, rsq84, var84, coeff84,
     file = paste(mydir, "GSD_simulations.RData", sep = ""))

