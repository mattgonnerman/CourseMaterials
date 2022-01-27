lapply(c("dplyr","ggplot2","R2jags","RMark", "gtools"),
       require, #Use the require function to load packages
       character.only = T, #Necessary for lapply
       quietly = T) #Whether to report successful loading of packages 


### Describing the model system ###
## Simulated turkey survival
# Base daily survival rate
mean.DSR <- 0.999  

## Regression Coefficients
# Intercept
alpha.s <- logit(mean.DSR) #Use logit link to create an intercept

# For every 1 sd increase in snow accumulation, the logit-linked probability of surviving a given day changes by ...
beta.snow <- -.002 

# For every 1 sd increase in max daily temperature, the logit-linked probability of surviving a given day changes by ...
beta.temp <- -.25 

# For every 1 sd increase in body condition, the logit-linked probability of surviving a given day changes by ...
beta.body <- -.003 
gtools::inv.logit(alpha.s) - gtools::inv.logit(alpha.s+beta.body)

# In comparison to suburban landscapes, turkeys in a agricultural and forested landscape differ in logit-linked probability of surviving a given day by ...
beta.ag <- .1
beta.for <- -.5

#########################################################################################################
### Simulating our Model System ###

#### Simulate individual turkey characteristics and associated weather information

w.days <- 30  #the length of winter
n.turkey <- 200  #number of simulated turkeys

#Simulate individual turkey characteristics
turkey.sim <- matrix(NA, nrow = n.turkey, ncol = 3)
for(i in 1:n.turkey){
  turkey.sim[i,1] <- rnorm(1, mean = 0, sd = 1.5) #Body Condition Metric, referenced to mean body condition (pre-z-standardized)
  #Categorical covariate for Landscape, need dummy code it
  h <- sample(1:3, 1) #1 = Suburban, 2 = Agriculture, 3 = Forested
  turkey.sim[i,2] <-  ifelse(h == 2, 1, 0)
  turkey.sim[i,3] <-  ifelse(h == 3, 1, 0)
}

#Simulate daily weather characteristics
weather.sim <- matrix(NA, nrow = 2, ncol = w.days) #Create an empty matrix to place daily weather information into.
#For loop to populate the matrix
for(j in 1:w.days){
  #Snow Accumulation
  weather.sim[1,j] <- rlnorm(1, meanlog = log(1), sdlog = log(2)) #continuous lognormal distribution
  #Temperature
  weather.sim[2,j] <- rnorm(1, mean = 20, sd = 10) #sample discrete unifrom. option to provide a vector of probabilities to weight samples.
}

#Since we created our weather data to match how it is normally measured, lets next z-standardize the simulated values so they make sense with our beta coefficients.
for(k in 1:2){
  weather.sim[k,] <- scale(weather.sim[k,], center = T, scale = T)
}

# Each turkey will also have unique pressures and differences in behavior which we did not measure but will still cause changes in survival rates, so we need to simulate these influences for each bird at each time step.
random.sim <- matrix(NA, nrow = n.turkey, ncol = w.days)
for(i in 1:n.turkey){
  for(j in 1:w.days){
    random.sim[i,j] <- rnorm(1, mean = 0, sd = 1)
  }
}

#We can now simulate each turkeys individual survival outcomes
#First lets find combine the weather information with the betas we defines to create and unlinked probability
ind.daily.prob <- matrix(NA, nrow = n.turkey, ncol = w.days)
for(i in 1:n.turkey){
  for(t in 1:w.days){
    ind.daily.prob[i,t] <- alpha.s + beta.snow*weather.sim[1,t] + beta.temp*weather.sim[2,t] + beta.body*turkey.sim[i,1] + beta.ag*turkey.sim[i,2] + beta.for*turkey.sim[i,3] + random.sim[i,t]
    
    #We can now use the logit link to tranform into probabilities of surviving a given day
    ind.daily.prob[i,t] <- gtools::inv.logit(ind.daily.prob[i,t])
  }
}

#We can use a bernoulli trial (weighted coin flip) to simulate survival, assuming a bird was alive, need to make sure that they all start alive at capture, all birds assumed caught Dec 31. 

ind.surv.hist <- matrix(NA, nrow = n.turkey, ncol = w.days + 1)
ind.surv.hist[,1] <- 1

#Combine if/else statement with forloop to simulate

for(i in 1:n.turkey){
  for(t in 2:(w.days+1)){
    if(ind.surv.hist[i, t-1] == 1){
      ind.surv.hist[i, t] <- rbinom(1,1, ind.daily.prob[i,(t-1)])
    }else{
      break #ends the loop for this bird if it has died
    }
  }
}

#Before we move on, lets change the NAs in the matrix to 0, as once a bird dies it is dead for good.
ind.surv.hist[is.na(ind.surv.hist)] <- 0

#Encounter history for first bird
ind.surv.hist[1,]

#Lets see how many survived
sum(ind.surv.hist[,w.days + 1], na.rm = T)


####Simulate Resulting Survival history
obs.hist.mat <- matrix(NA, nrow = n.turkey, ncol = w.days + 1)
obs.hist.mat[,1] <- 1 #We know the bird was alive at capture
#Now we can simulate our technicians going out and doing the work to collect our data.
for(i in 1:n.turkey){
  for(t in 2:(w.days+1)){
    succ.obs <- rbinom(1,1, obs.rate)
    
    if(succ.obs == 1){
      obs.hist.mat[i, t] <- ind.surv.hist[i, t]
      
      if(ind.surv.hist[i, t] == 0){ #if the bird is dead we would stop checking on it
        break #so we exit loop for this bird when we identify a mortality
      }
    }
  }
}

JAGS.turkey.S.model <- function(){
  ###Priors
  base.S ~ dbeta(1,1) 
  a.S <- logit(base.S) 
  b.snow ~ dnorm(0, 1/1000)
  b.temp ~ dnorm(-1, 1/10)
  b.bc ~ dnorm(-1, 1/2) 
  b.ag ~ dnorm(0, 1/1000)
  b.for ~ dnorm(0, 1/1000)
  
  
  ###LIKELIHOOD
  for(i in 1:n.ind){
    for(t in 1:n.days){
      logit(DSR[i,t]) <- a.S + b.snow*x.snow[t] + b.temp*x.temp[t] + b.bc*x.bc[i] + b.ag*x.ag[i] + b.for*x.for[i]
    }
  }
  
  for(j in 1:n.obs){
    obs.S[j] ~ dbern(period.S[j]) #Probability of survival since previous check
    period.S[j] <- prod(DSR[turkey.id[j], (day.found[j] - (1:days.since[j]))]) #Aggregrate survival across missed days
  }
}