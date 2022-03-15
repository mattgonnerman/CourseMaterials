lapply(c("dplyr",
         # "ggplot2", 
         "R2jags",
         "RMark", 
         "gtools"),
       require,
       character.only = T,
       quietly = T) #Whether to report successful loading of packages 


##Simulated turkey survival
w.days <- 90  #number of days turkeys are monitored for
mean.DSR <- 0.99  #base daily survival rate
alpha.s <- logit(mean.DSR) #Use logit link to create an intercept

# Effects on daily survival
beta.body <- .003 #Body Condition
beta.ag <- .5 #Agriculture
beta.for <- -1 #Forested


#Simulate individual turkey characteristics
n.turkey <- 250  #number of simulated turkeys

turkey.sim <- matrix(NA, nrow = n.turkey, ncol = 5)
for(i in 1:n.turkey){
  #Body Condition Metric, referenced to mean body condition (pre-z-standardized)
  turkey.sim[i,1] <- rnorm(1, mean = 0, sd = 1) 
  
  #Categorical covariate for Landscape, need dummy code it
  h <- sample(1:3, 1) #1 = Suburban, 2 = Agriculture, 3 = Forested
  turkey.sim[i,2] <-  ifelse(h == 2, 1, 0)
  turkey.sim[i,3] <-  ifelse(h == 3, 1, 0)
  
  #Random variation associated with individual turkeys
  turkey.sim[i,4] <- rnorm(1, mean = 0, sd = .1) 
  
  #We can now combine individual betas to define individual survival rates
  turkey.sim[i,5] <- gtools::inv.logit(alpha.s + beta.body*turkey.sim[i,1] + beta.ag*turkey.sim[i,2] +
                                             beta.for*turkey.sim[i,3] + turkey.sim[i,4])
    
}



#We can use a bernoulli trial (weighted coin flip) to simulate survival, assuming a bird was alive, need to make sure that they all start alive at capture, all birds assumed caught Dec 31. 
ind.surv.hist <- matrix(NA, nrow = n.turkey, ncol = w.days+1)
ind.surv.hist[,1] <- 1 #First day is day of capture, set to 1 since known alive

#Combine if/else statement with forloop to simulate
for(i in 1:n.turkey){
  for(t in 2:(w.days+1)){
    if(ind.surv.hist[i, t-1] == 1){ #If the bird was alive at previous visit
      ind.surv.hist[i, t] <- rbinom(1,1, turkey.sim[i,5]) 
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


obs.rate <- .95 #We can't expect to observe each turkey every day, so this just randomizes whether we successfully collect data on a bird on a given day. Let's assume that on average we find a bird every 4 days.

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
  b.bc ~ dnorm(0, 1/100) 
  b.ag ~ dnorm(0, 1/100)
  b.for ~ dnorm(0, 1/100)
  
  
  ###LIKELIHOOD
  for(i in 1:n.ind){
      logit(DSR[i]) <- a.S + b.bc*x.bc[i] + b.ag*x.ag[i] + b.for*x.for[i]
    }
  
  for(j in 1:n.obs){
    period.S[j] <- pow(DSR[turkey.id[j]], interval[j]) # Total probability of survival since previous visit
    obs.S[j] ~ dbern(period.S[j]) # Observation process follows bernoulli distribution
  }
}


### Format data to make compatible with model specification
# Remove first column (capture day) from encounter history, will not be used
obs.hist.mat.noDay1 <- obs.hist.mat[,-1]

# JAGS doesn't like NAs, so we need to simplify observation history. 
# An easy way to do this is to convert it to a vector
obs.hist.vec <- as.vector(t(obs.hist.mat.noDay1))
obs.hist <- obs.hist.vec[!is.na(obs.hist.vec)]

# Vector format requires a bit more information to function correctly
# For example, we need a vector to tell the model which bird we are referencing for each visit
ID <- matrix(1:n.turkey, nrow = n.turkey, ncol = w.days)
ID <- as.vector(t(ID))
ID <- ID[!is.na(obs.hist.vec)]      # ID marker

# We also need to know how long since the last visit
visit_ind <- matrix(NA, ncol = ncol(obs.hist.mat), nrow = nrow(obs.hist.mat))
get.last <- function(x) max(which(!is.na(x)))

for(i in 1:n.turkey){
  for(j in 2:(w.days+1)){
    if(is.na(obs.hist.mat[i,j]) == FALSE){
      visit_ind[i,j] <- j - get.last(obs.hist.mat[i,1:(j-1)])
    }
  }
}

interval <- as.vector(t(visit_ind))
interval <- interval[!is.na(interval)]





#Data to supply to the model
data.list <- list(
  #Observation Information
  obs.S = obs.hist, #Vector containing observed survival history
  turkey.id = ID, #Vector describing which turkey an observation was for
  interval = interval,#Time since last visit
  
  n.ind = nrow(turkey.sim), #Number of individual turkeys monitored
  n.obs = length(obs.hist), #Number of total observations
  
  #Covariate Information
  x.bc = turkey.sim[,1], #Vector describing body condition for a given turkey
  x.ag = turkey.sim[,2], #Vector describing whether a turkey primarily uses agriculture landscape
  x.for = turkey.sim[,3] #Vector describing whether a turkey primarily uses forested landscape
)

#Parameter Monitors
parameters.null <-c(
  #Regression Coefficients
  "a.S",
  "b.bc",
  "b.ag",
  "b.for"
)

#Initial Values
inits.null <- function(){
  list(
    base.S = .998 #Set initial value for base survival rate
  )
}


#Make final MCMC related decisions
ni <- 2000 #Number of simulations to run
nb <- 1000 #Length of burnin in period
nt <- 1 #How much thinning should there be
nc <- 3 #Number of chains to run simultaneously

#Run the model
Surv_JAGS_output <- R2jags::jags(data = data.list,
                         parameters.to.save = parameters.null,
                         inits = inits.null,
                         model.file = JAGS.turkey.S.model,
                         n.iter = ni,
                         n.burnin = nb,
                         n.thin = nt,
                         n.chains = nc) 


R2jags::traceplot(Surv_JAGS_output, #bugs output object
          mfrow = c(3,2), #Rows and Columns in figure
          ask = F, #Ask before plotting?
          varname = c("b.ag", "b.for", "b.bc")) #Which parameters to plot

#Basic Parameter Estimates
print(Surv_JAGS_output)



######################################
### RMARK
#######################################
require(RMark)
EH.rmark <- as.data.frame(obs.hist.mat)[,] %>%
  mutate(TurkID = row_number()) %>%
  rowwise() %>%
  mutate(FirstFound = 1, #Conveniently, all of our simulated birds were caught on the same day
         LastPresent = max(which(c_across(1:(w.days)) == 1)), #When was the turkey last observed alive
         LastChecked = max(which(!is.na(c_across(1:(w.days))))), #When was the turkey first observed dead
         Fate = min(c_across(everything()), na.rm = T), #What was the final known fate of a bird through the period of itnerest?
         Freq = 1) %>% #Number of repetitions of a given encounter history, I don't usually aggregate
  mutate(LastChecked = ifelse(Fate == 1, w.days + 1, LastChecked)) %>%
  # mutate(LastChecked = ifelse(FirstDead > LastAlive)) %>%
  dplyr::select(TurkID,
                FirstFound, 
                LastPresent,
                LastChecked,
                Fate,
                Freq)

indcov.rmark <- as.data.frame(turkey.sim[,1:3]) %>%
  rename(BC = V1, Ag = V2, For = V3) %>%
  mutate(TurkID = row_number(), 
         Ag = as.factor(Ag), #RMark will convert automatically, but can preemptively do this to avoid messages
         For = as.factor(For))

data.rmark <- merge(EH.rmark, indcov.rmark, by = "TurkID")

WintSurv.process = process.data(data.rmark,
                                model="Nest", #We will run this in a nest survival framework
                                nocc=w.days + 1, #Number of occassions in the EH
                                groups=c("Ag", "For")) #Identify categorical covariates here
WintSurv.ddl = make.design.data(WintSurv.process)

RMark.model <- mark(WintSurv.process, 
                    WintSurv.ddl, 
                    model.parameters=list(S=list(formula =~ BC + Ag + For, 
                                                 link="logit")))

summary(RMark.model)
