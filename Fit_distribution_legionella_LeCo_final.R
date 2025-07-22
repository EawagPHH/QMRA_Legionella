graphics.off() 
rm(list=ls(all=TRUE)) 
source("DBDA2E-utilities.R") #this provides essential functions for JAGS
fileNameRoot="A"
library(hdrcde)

### Import the dataset named Legionella_data
###
mydataframe <-Legionella_data
data = dplyr::filter(mydataframe, Temperature=="Hot",Sample_type == "First draw") #Here you can chose first draw sample, flushed sample, cold water and hot water
y = unlist(data[, "Count"])
x = unlist(data[, "Volume"]*1000)
conc=log((y + 1)/(0.001 * x)) # concentrations using a LOD of 1 CFU/L

# Model Poisson-Gamma  ---------------------------

#Packaging data for JAGS:
dataList = list(
  y = y,
  x = x,
  N = length(y)
)

# Building the model:
modelstring = "
model {
for ( i in 1:N ) {
  #Likelihood
    y[i] ~ dpois(lamu[i]*x[i]/1000)
    lamu[i] ~ dgamma(shape,rate)
}
# Priors
shape ~ dunif(0, 1)
rate ~ dunif(1.0E-12, 1.0E-1)
}
"
# Write model to a file, and send to JAGS:
writeLines(modelstring,con="model.txt")

# MCMC simulation
parameters = c("shape", "rate")  
adaptSteps = 1.0E+3  # Number of steps to "tune" the samplers.           
burnInSteps = 1.0E+3  # Burn-in period.            
nChains = 3    # Number of chain.       
numSavedSteps = 1.0E+5 # Number of iterations after burn-in.         
thinSteps = 1   # Thinning factor.
nPerChain = ceiling(( numSavedSteps * thinSteps ) / nChains) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel9 = jags.model( "model.txt" , data=dataList, 
                         n.chains=nChains, n.adapt=adaptSteps)
# Burn-in:
cat("Burning in the MCMC chain...\n")
update(jagsModel9 , n.iter=burnInSteps)

# The saved MCMC chain:
cat("Sampling final MCMC chain...\n")
mcmcCoda9 = coda.samples(jagsModel9 , variable.names=parameters,
                          n.iter=nPerChain, thin=thinSteps)
mcmcChain9 = as.matrix(mcmcCoda9)
chainLength9 = NROW(mcmcChain9)

# Display diagnostics of chain.
parameterNames = varnames(mcmcCoda9)
for ( parName in parameterNames ) {
  diagMCMC(codaObject=mcmcCoda9, parName=parName) 
}


# Contour plot of parameter values
df = as.data.frame(mcmcChain9)
o <- hdr.2d(x=df[,2], y=log10(1/df[,1]), prob = c(0.01, 0.05, 0.1))
openGraph(height=5,width=7)
plot(o, xlab = "shape",ylab = "log10(scale)")
legend("topleft", legend = o$alpha, title = "Credible region",
       fill = gray((length(o$alpha):1)/(length(o$alpha)+1)))

# Model Poisson-lognormal ---------------------------

#Packaging data for JAGS:
dataList = list(
  y = y ,
  x = x,
  N = length(y)
)

# Building the model:
modelstring = "
model {
for (i in 1:N) {
y[i] ~ dpois(lambda[i] * x[i] / 1000)
lambda[i] ~ dlnorm(meanlog, 1/sdlog^2)
}
# Priors
meanlog ~ dunif(-10,10)
sdlog ~ dexp(1)
}
" 
# Write model to a file, and send to JAGS:
writeLines(modelstring,con="model.txt")

# MCMC simulation
parameters = c( "meanlog", "sdlog")  
adaptSteps = 1.0E+3 # Number of steps to "tune" the samplers.              
burnInSteps = 1.0E+3 # Burn-in period.          
nChains = 3 # Number of chain.                   
numSavedSteps = 1.0E+5 # Number of iterations after burn-in.            
thinSteps = 1  # Thinning factor.                  
nPerChain = ceiling( (numSavedSteps * thinSteps) / nChains) 
# Create, initialize, and adapt the model:
jagsModel2 = jags.model("model.txt" , data=dataList, 
                         n.chains=nChains , n.adapt=adaptSteps)
# Burn-in:
cat("Burning in the MCMC chain...\n")
update(jagsModel2 , n.iter=burnInSteps)

# The saved MCMC chain:
cat("Sampling final MCMC chain...\n")
mcmcCoda2 = coda.samples(jagsModel2, variable.names=parameters ,
                          n.iter=nPerChain, thin=thinSteps)
mcmcChain2 = as.matrix(mcmcCoda2)
chainLength2 = NROW(mcmcChain2)

# Display diagnostics of chain.
parameterNames = varnames(mcmcCoda2)
for (parName in parameterNames) {
  diagMCMC(codaObject = mcmcCoda2 , parName = parName)
}


# Contour plot for parameter values
df = as.data.frame(mcmcChain2)
o <- hdr.2d(x=df[,1], y=df[,2], prob = c(0.01, 0.05, 0.1,0.2))
openGraph(height=5,width=7)
plot(o, xlab = "meanlog",ylab = "sdlog")
legend("topleft", legend = o$alpha, title = "Credible region",
       fill = gray((length(o$alpha):1)/(length(o$alpha)+1)))

#########################
# Plot CCDF ------------------------------------------------------------------------------

tiff(file = paste("Legio_HOT_CCDF.tiff", sep = ""), width = 3543, height = 3543, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
par(mar=c(4.1, 4.1, 1.1, 1.1))

plot(point_x, point_y, col="white", pch = 19, xlab = expression(paste(bolditalic("L. pneumophila"),
    bold(" (CFU/L)"))), ylab = "Exceedance probability", xaxt="n", xlim = range(myTicks),ylim=c(0.002, 1),
    font.lab = 2, log='xy')
x_ticks <- 10^seq(0, ceiling(log10(max(lnorm_xMLE))), by = 1)  
axis(1, at = x_ticks, labels = parse(text = paste0("10^", log10(x_ticks))))

#Plot log-normal distribution ------------------------------------------------------------------------------

#Simulating uncertainty interval:
HDI_output = matrix(1:80, nrow = 40, ncol = 2)
pb = txtProgressBar(min = 0, max = 40, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 40, 1)) {
  output <- numeric(length = nrow(mcmcChain2))
  for (i in 1:nrow(mcmcChain2)) {            
    output[[i]] = 1 - plnorm(10 ^ (j / 4), meanlog = mcmcChain2[i,"meanlog"], sdlog = mcmcChain2[i,"sdlog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 10, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(1, 0, 0, 0.2),border = NA) 

# Simulation and representation of the best fit distribution:
Mean_meanlog = mean(mcmcChain2[,"meanlog"])
Mean_sdlog = mean(mcmcChain2[, "sdlog"])
lnorm_dataMLE = rlnorm(n = 1.0E+6, meanlog = Mean_meanlog, sdlog = Mean_sdlog)
lnorm_xMLE = sort(lnorm_dataMLE)
lnorm_yMLE = 1 - ecdf(lnorm_dataMLE)(sort(lnorm_dataMLE))
lines(lnorm_xMLE, lnorm_yMLE, col="red", lwd = 2,lty = 1)

# Plot gamma distribution ------------------------------------------------------------------------------

# Simulating uncertainty interval:
HDI_output_B = matrix(1:80, nrow = 40, ncol = 2)
pb = txtProgressBar(min = 0, max = 40, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1,40,1)) {
  output <- numeric(length = nrow(mcmcChain9))
  for (i in 1:nrow(mcmcChain9)) {            
    output[[i]] = 1-pgamma(10^(j/4), shape = mcmcChain9[i,"shape"], rate = mcmcChain9[i,"rate"])
  }
  HDI_output_B[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixD = cbind(c(10 ^ (seq(0.25,10,0.25))), HDI_output_B)
MatrixE = MatrixD[(MatrixD[,2] > 0),]
MatrixF = MatrixD[(MatrixD[,3] > 0),]
polygon(c(MatrixF[,1], rev(MatrixE[,1])), c(MatrixF[,3], rev(MatrixE[,2])), col = rgb(0, 0, 1, 0.3), border = NA)

# Simulation and representation of the best fit distribution:
Mean_shapeNB = mean(mcmcChain9[,"shape"])
Mean_rateNB = mean(mcmcChain9[,"rate"])
gamma_dataMLE = rgamma(n = 1.0E+6, shape = Mean_shapeNB, rate = Mean_rateNB)
gamma_xMLE = sort(gamma_dataMLE)
gamma_yMLE = 1 - ecdf(gamma_dataMLE)(sort(gamma_dataMLE))
lines(gamma_xMLE,gamma_yMLE, col = "blue", lwd = 2,lty = 1)

# Data point ------------------------------------------------------------------------------
N = length(y)
point_x = sort(1000 * (y / x))
point_y = (1 + 1 / N) - ecdf(1000 * y / x)(sort(1000 * y / x) )
points(point_x, point_y ,col = "black", bg = "grey", pch = 19, lwd = 1, cex = 1)

legend("topright", 
       legend = c("PLA", "PGA", "Observed data"),
       col = c("red", "blue", "black"),
       lty = c(1, 1, NA),  
       pch = c(NA, NA, 19), 
       lwd = c(2, 2, NA),   
       pt.bg = c(NA, NA, "grey"), 
       pt.cex = c(NA, NA, 1),     
       inset = 0.02,              
       bg = "white",             
       bty = "n")               

dev.off()

