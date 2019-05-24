### run the monte carlo simulation #################################################
# run the analysis for k = 100, 200, 400, keep alpha at 0.05

m = 10000 #set the number of Monte Carlo simulations to run

# present the analysis graphically using confidence bounds
qntls = seq(0.05, 0.95,0.05)
cfdlab = MakeCFD(all1.data$D)
ci100 = WolmanCI(cfdlab, 100, qntls, alpha = 0.05)
#ci200 = WolmanCI(cfd1all, 200, qntls, alpha = 0.05)
ci400 = WolmanCI(cfdlab, 400, qntls, alpha = 0.05)

#sample some data from the distribution
k = 100
sample.100 = matrix(data = NA, nrow = k, ncol = 25)
for(i in seq(1,25,1)){
  sample.100[,i] = sort(sample(all1.data$D, replace = T, size = k))
}

#sample some data from the distribution
k = 400
sample.400 = matrix(data = NA, nrow = k, ncol = 25)
for(i in seq(1,25,1)){
  sample.400[,i] = sort(sample(all1.data$D, replace = T, size = k))
}

# PERFORM THE BOOTSTRAP ANALYSIS USING LAB DATA
k = 100  #set the sample size
tmp = matrix(data=NA, nrow = m, ncol = length(qntls))
for(i in seq(1,m,1)){
  tmp2 = sort(sample(all1.data$D, replace = T, size = k))
  tmp[i,] = tmp2[k*qntls]
}
tmp = apply(tmp,2,sort) #sort each quantile estimate
upper = tmp[m*0.975,]  #estimated upper bounds by bootstrap
estimate = tmp[m*0.5,]
lower = tmp[m*0.025,]  #estimated upper bounds by bootstrap
bootstrap100 = data.frame(qntls,estimate, lower,upper)

# k = 200  #set the sample size
# for(i in seq(1,m,1)){
#   tmp2 = sort(sample(all1.data$D, replace = T, size = k))
#   tmp[i,] = tmp2[k*qntls]
# }
# tmp = apply(tmp,2,sort) #sort each quantile estimate
# upper = tmp[m*0.95,]  #estimated upper bounds by bootstrap
# estimate = tmp[m*0.5,]
# lower = tmp[m*0.05,]  #estimated upper bounds by bootstrap
# bootstrap200 = data.frame(qntls,estimate, lower,upper)

k = 400  #set the sample size
for(i in seq(1,m,1)){
  tmp2 = sort(sample(all1.data$D, replace = T, size = k))
  tmp[i,] = tmp2[k*qntls]
}
tmp = apply(tmp,2,sort) #sort each quantile estimate
upper = tmp[m*0.975,]  #estimated upper bounds by bootstrap
estimate = tmp[m*0.5,]
lower = tmp[m*0.025,]  #estimated upper bounds by bootstrap
bootstrap400 = data.frame(qntls,estimate, lower,upper)

# results = matrix(data = NA, ncol = length(ci$probs), nrow = m)
# 
# for(i in seq_along(qntls)){
#   sample.data = sample(all1.data$D, replace = T, size = k)
#   cfdtmp = MakeCFD(sample.data)
#   results[i,] = WolmanCI(cfdtmp, k, probs = qntls)[[2]]
# }
# 
# mclimits = matrix(data = NA, ncol = 2, nrow = length(ci$probs))
# for(i in seq_along(mclimits[,1])){
#   mclimits[i,] = as.numeric(quantile(results[,i], probs = c(0.005, 0.995)))
# }

# repeat above for simulated data from a known distribution
sim.gsd = 2^rnorm(1000000, mean = 5.5, sd = 1)
cfdsim = MakeCFD(sim.gsd)
sim100 = WolmanCI(cfdsim, 100, qntls, alpha = 0.05)
sim400 = WolmanCI(cfdsim, 400, qntls, alpha = 0.05)

k = 100  #set the sample size
tmp = matrix(data=NA, nrow = m, ncol = length(qntls))
for(i in seq(1,m,1)){
  tmp2 = sort(sample(sim.gsd, replace = T, size = k))
  tmp[i,] = tmp2[k*qntls]
}
tmp = apply(tmp,2,sort) #sort each quantile estimate
upper = tmp[m*0.975,]  #estimated upper bounds by bootstrap
estimate = tmp[m*0.5,]
lower = tmp[m*0.025,]  #estimated upper bounds by bootstrap
simstrap100 = data.frame(qntls,estimate, lower,upper)

k = 400  #set the sample size
tmp = matrix(data=NA, nrow = m, ncol = length(qntls))
for(i in seq(1,m,1)){
  tmp2 = sort(sample(sim.gsd, replace = T, size = k))
  tmp[i,] = tmp2[k*qntls]
}
tmp = apply(tmp,2,sort) #sort each quantile estimate
upper = tmp[m*0.975,]  #estimated upper bounds by bootstrap
estimate = tmp[m*0.5,]
lower = tmp[m*0.025,]  #estimated upper bounds by bootstrap
simstrap400 = data.frame(qntls,estimate, lower,upper)

### perform the simulation of log normal error #####################################


NSim = 3000
SErr84 = array(data = NA, dim =  c(NSim,4))
SErr50 = array(data = NA, dim =  c(NSim,4))

for(i in seq(1,NSim,1)){
  log50 = runif(1, min = 4.5, max = 6.5)  #vary d50 from 22.6 mm to 90.5 mm
  logSD = runif(1, min = 0.5, 2.5)  #vary D16 to D84 range from factor of 2 to factor of 16
  #logSD = 0.5
  n = round(runif(1, min = 50, max = 1000))  #generate sample sizes from 100 to 500
  obs = 2^rnorm(n,log50,logSD)  #take n samples from a log normal distribution
  gsd = MakeCFD(obs, plot = F)  #extract the cumulative distribution from the sample
  
  #estimate the error for the D84 from the distribution
  p = c(0.50,0.84)
  tmp = WolmanCI(gsd,n,p)  
  err = 0.5*(tmp$upper - tmp$lower)/tmp$estimate
  SErr50[i,] = c(n,2*logSD,log50,err[1])
  SErr84[i,] = c(n,2*logSD, log50, err[2])
}

SErr84 = as.data.frame(SErr84)
colnames(SErr84) = c("N", "si", "log50", "err")
SErr50 = as.data.frame(SErr50)
colnames(SErr50) = c("N", "si", "log50", "err")

powfit50 = lm(log(SErr50$err) ~ log(SErr50$N) + SErr50$si)
rsq50 = as.numeric(summary(powfit50)[9])
tmp = as.matrix(anova(powfit50)[2])
var50 = as.numeric(tmp[1:3]/sum(tmp[1:3]))
coeff50 = as.numeric(powfit50$coefficients)

powfit84 = lm(log(SErr84$err) ~ log(SErr84$N) + SErr84$si)
rsq84 = as.numeric(summary(powfit84)[9])
tmp = as.matrix(anova(powfit84)[2])
var84 = as.numeric(tmp[1:3]/sum(tmp[1:3]))
coeff84 = as.numeric(powfit84$coefficients)


### clean up the workspace #########################################################

save(m, ci100, ci400, sample.100, sample.400, bootstrap100, bootstrap400,
     sim100, sim400, simstrap100, simstrap400,
     NSim, SErr84, SErr50, 
     powfit50, rsq50, var50, coeff50,
     powfit84, rsq84, var84, coeff84,
     file = paste(mydir, "GSD_simulations.RData", sep = ""))

