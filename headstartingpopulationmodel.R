rm(list=ls())

## Required Packages ##
library(progress)       # progress bar
library(tidyverse)      # ggplot
library(pracma)         # maths for stochasticity
library(purrr)          # programming tools
library(hrbrthemes)     # plot themes
library(reshape)        # data wrangling
library(RColorBrewer)   # colour palettes


## Set the current population in Chitwan ##
# wild females 1-9, adult females, wild males 1-18, adult males, head-started females 1r 1n 2-9, head-started males 1r 1n 2-18 #
currentpop <- c(20,2, rep(0, times=7), 30,    7, 1, rep(0, times=16), 3,    65,5, (21:18), 8  ,(16:14),     (9:1), rep(0, times=10))
length(currentpop)
sum(currentpop)
currentpop

### Parameters ###
# Set the baseline parameter values #

B <- 0.8*39            # birthrate 0.8*39 reproductive frequency * clutch size
seB <- 0.02            # birthrate standard error 0.02
EM <- 0.5              # early maturity rate 0.5
Sa <- 0.95             # adult survival 0.95
seSa <- 0.003          # adult survival standard error 0.003
XR <- 0.75             # sex ratio 1:3
Sh <- 0.1              # hatchling survival 0.1
Sim2 <- 0.2            # immature survival 0.2, increasing by 0.1 each year
Sim3 <- Sim2+0.1       #0.3
Sim4 <- Sim2+0.2       #0.4
Sim5 <- Sim2+0.3       #0.5
Sim6 <- Sim2+0.4       #0.6
Sim7 <- Sim2+0.5       #0.7
Sim8 <- Sim2+0.6       #0.8
Sim9 <- Sim2+0.7       #0.9
seSim <- 0.002         # immature survival standard error 0.002

catprob <- 0.2         # catastrophe probability
cateffect <- 0.1       # catastrophe effect 0.1*hatchling survival = 0.01

# head-starting parameters #

currentreleases <- 447 # eggs collected for first five years 447
secr <- 57             # eggs collected standard error 57
removal <- 0.9         # proportion of eggs removed for headstarting 0.9
DN <- (1-0.45)         # dispersal survival Narayani 0.55
seDN <- 0.031          # standard error^ 0.031
DR <- (1-0.15)         # dispersal survival Rapti 0.85
seDR <- 0.006          # standard error^ 0.006
Sc <- 0.2              # captivity survival 0.2
riversplit <- 0.8      # headstart location release proportion
releaseXR <- 0.8       # head-start releases sex ratio

# head-start survival (years post release)
releaseSY1r <- 0.6     # Rapti 0.6
sereleaseSY1r <- 0.02  # standard error 0.02
releaseSY1n <- 0.3     # Narayani 0.3
sereleaseSY1n <- 0.07  # standard error 0.07
releaseSY2 <- 0.63     # Year 2 post release, all individuals 0.63
sereleaseSY2 <- 0.012  # standard error 0.012
releaseSY3 <- 0.68     # Year 3 post release 0.68
sereleaseSY3 <- 0.073  # standard error 0.073
releaseSY4 <- releaseSY3+0.055    # Year 4 post release 0.735
releaseSY5 <- releaseSY3+0.11     # Year 5 post release 0.79     
releaseSY6 <- releaseSY3+0.165    # Year 6 post release 0.845
releaseSY7 <- releaseSY3+0.22     # Year 7 post release 0.9


## Demographic Stochasticity functions ##

# norminv, returns the inverse of the normal cumulative distribution 
norminv <- function(p,mu=0,sigma=1) { 
  x0 <- -sqrt(2)*pracma::erfcinv(2*p)
  x0*sigma+mu
}

# Survival/dispersal function
Survival <- function(Estimate = 0.8, SE = 0.01){
  
  LogEst <- log(Estimate/(1-Estimate))
  SEt <- SE/(Estimate*(1-Estimate))
  Random <- norminv(runif(1), LogEst, SEt)
  RealRandom <- exp(Random)/(1 + exp(Random))
  
  data.frame(Estimate, SE, LogEst, SEt, Random, RealRandom)
  
}

# Reproduction function
Reproduction <- function(Estimate = 0.8, SE = 0.01){
  
  LogEst <- log(Estimate)
  SEt <- SE/Estimate
  Random <- norminv(runif(1), LogEst, SEt)
  RealRandom <- exp(Random)
  
  data.frame(Estimate, SE, LogEst, SEt, Random, RealRandom)
  
  RealRandom
  
}

################################################################################
############# DETERMENISTIC VERSION ############################################
################################################################################

## one time step function ##
step <- function(pop){ #input initial population
  
  matrix <- matrix(0, ncol = 58, nrow = 58)          # create the transition matrix

## Wild females ##  
  # female offspring 
  matrix[1,10] <- B * (1-removal) * XR * Sh
  # male offspring
  matrix[11,10] <- B * (1-removal) * (1-XR) * Sh
  # female early maturity
  matrix[10,5] <- Sim5 * EM
  matrix[6,5] <- Sim5 * (1-EM)
  # immature survival each year
  matrix[2,1] <- Sh
  matrix[3,2] <- Sim2
  matrix[4,3] <- Sim3
  matrix[5,4] <- Sim4
  matrix[7,6] <- Sim6
  matrix[8,7] <- Sim7
  matrix[9,8] <- Sim8
  matrix[10,9] <- Sim9
  # adult survival
  matrix[10,10] <- Sa
  
## Wild males ##
  # early maturity
  matrix[29,20] <- Sim9 * EM
  matrix[21,20] <- Sim9 * EM
  # immature survival each year
  matrix[12,11] <- Sh
  matrix[13,12] <- Sim2
  matrix[14,13] <- Sim3
  matrix[15,14] <- Sim4
  matrix[16,15] <- Sim5
  matrix[17,16] <- Sim6
  matrix[18,17] <- Sim7
  matrix[19,18] <- Sim8
  matrix[20,19] <- Sim9
  matrix[22,21] <- Sim9
  matrix[23,22] <- Sim9
  matrix[24,23] <- Sim9
  matrix[25,24] <- Sim9
  matrix[26,25] <- Sim9
  matrix[27,26] <- Sim9
  matrix[28,27] <- Sim9
  matrix[29,28] <- Sim9
  # adult survival
  matrix[29,29] <- Sa
  
## FEMALE HEADSTARTS ##
  matrix[32,30] <- releaseSY1r
  matrix[32,31] <- releaseSY1n
  matrix[33,32] <- releaseSY2
  matrix[34,33] <- releaseSY3
  matrix[35,34] <- releaseSY4
  matrix[10,35] <- releaseSY5 * EM
  matrix[36,35] <- releaseSY5 * (1-EM)
  matrix[37,36] <- releaseSY6
  matrix[38,37] <- releaseSY7
  matrix[39,38] <- releaseSY7
  matrix[10,39] <- releaseSY7
  
## MALE HEADSTARTS ##
  matrix[42,40] <- releaseSY1r
  matrix[42,41] <- releaseSY1n
  matrix[43,42] <- releaseSY2
  matrix[44,43] <- releaseSY3
  matrix[45,44] <- releaseSY4
  matrix[46,45] <- releaseSY5
  matrix[47,46] <- releaseSY6
  matrix[48,47] <- releaseSY7
  matrix[49,48] <- releaseSY7
  matrix[50,49] <- releaseSY7
  matrix[51,50] <- releaseSY7 * EM
  matrix[29,50] <- releaseSY7 * (1-EM)
  matrix[52,51] <- releaseSY7
  matrix[53,52] <- releaseSY7
  matrix[54,53] <- releaseSY7
  matrix[55,54] <- releaseSY7
  matrix[56,55] <- releaseSY7
  matrix[57,56] <- releaseSY7
  matrix[58,57] <- releaseSY7
  matrix[29,58] <- releaseSY7
  
  names <- c("1F","2F","3F","4F","5F","6F","7F","8F","9F","AF","1M","2M","3M","4M","5M","6M","7M","8M","9M","10M","11M","12M","13M","14M","15M","16M","17M","18M","AM","hs1rF","hs1nF","hs2F","hs3F","hs4F","hs5F","hs6F","hs7F","hs8F","hs9F","hs1rM","hs1nM","hs2M","hs3M","hs4M","hs5M","hs6M","hs7M","hs8M","hs9M","hs10M","hs11M","hs12M","hs13M","hs14M","hs15M","hs16M","hs17M","hs18M")
  rownames(matrix) <- names
  colnames(matrix) <- names

# density dependence of the number of hatchlings 
  if (pop[1]+pop[11] > 1000){
    pop[1] <- 750
    pop[11] <- 250
  }

# no adult males = no reproduction
  if (pop[29] < 1){
    matrix[1,10] <- 0
    matrix[11,10] <- 0
  }
  

  # multiply the initial population by the transition matrix to get the population after a year  
  newpop <- matrix %*% pop
  
  return(newpop[,1])
}



## Function for the deterministic model with releases ##
run <- function(initpop, runs=50){ #INPUT- initial population vector, number of runs
  pop <- initpop
  
  # create a dataframe for each years population output
  outputdf <- data.frame(matrix(0, nrow = runs + 1, ncol = length(pop)))
  colnames(outputdf) <- c("1F","2F","3F","4F","5F","6F","7F","8F","9F","AF","1M","2M","3M","4M","5M","6M","7M","8M","9M","10M","11M","12M","13M","14M","15M","16M","17M","18M","AM","hs1rF","hs1nF","hs2F","hs3F","hs4F","hs5F","hs6F","hs7F","hs8F","hs9F","hs1rM","hs1nM","hs2M","hs3M","hs4M","hs5M","hs6M","hs7M","hs8M","hs9M","hs10M","hs11M","hs12M","hs13M","hs14M","hs15M","hs16M","hs17M","hs18M")       #### change seq to initpop length -1 
  outputdf[1,] <- pop
  
  # create a matrix for for the number of collected eggs 
  collectedmatrix <- matrix(1, nrow=4)
  
   
  for (i in 1:runs){
    
    # run one step of the model
    pop <- step(pop)
    
    #Number of collected for the current step (will be added to population in 5 years)
    collected <- (pop[1]+pop[11]) / (1-removal) / Sh * removal
    # collected density dependence
    collected = pmin(collected, 800)

    collectedmatrix <- rbind(collectedmatrix, collected)

    # Add releases to the population
    if (i <= 4){
      pop[30] <- currentreleases * Sc * releaseXR * riversplit * DR
      pop[31] <- currentreleases * Sc * releaseXR * (1-riversplit) * DN
      pop[40] <- currentreleases * Sc * (1-releaseXR) * riversplit * DR
      pop[41] <- currentreleases * Sc* (1-releaseXR) * (1-riversplit) * DN
    }

    if (i > 4){
      pop[30] <- collectedmatrix[i] * releaseXR * Sc * riversplit * DR
      pop[31] <- collectedmatrix[i] * releaseXR * Sc * (1-riversplit) * DN
      pop[40] <- collectedmatrix[i] * (1-releaseXR) * Sc * riversplit * DR
      pop[41] <- collectedmatrix[i] * (1-releaseXR) * Sc * (1-riversplit) * DN
    }
    
    pop
    outputdf[i + 1, ] <- pop
    #print(round(pop, 2))  
    
  }
  return(outputdf)
}

## deterministic run of the model + total population at each step + plotted ##
dm <- run(currentpop, 50)
dmdf <- as.data.frame(dm)
dmdf$totalpop <- rowSums(dmdf)
dmdf$year <- 1:nrow(dmdf)
ggplot(dmdf, aes(y=totalpop, x=year)) + geom_line()
tail(dmdf$totalpop)
tail(dmdf$AF)
tail(dmdf$AM)
(tail(dmdf$AM,1)+tail(dmdf$AF,1))/tail(dmdf$totalpop,1)*100


################################################################################
##################STOCHASTIC VERSION############################################
################################################################################

# single time step of the stochastic model #
stochstep <- function(pop){ #input initial population
  
  # create the transition matrix
  matrix <- matrix(0, ncol = 58, nrow = 58)
  
  # add parameter values

## Wild females ##
  # female offspring 
  matrix[1,10] <- Reproduction(B, seB) * (1-removal) * Survival(XR)$RealRandom * Survival(Sh)$RealRandom
  # male offspring
  matrix[11,10] <- Reproduction(B, seB) * (1-removal) * (1-Survival(XR)$RealRandom) * Survival(Sh)$RealRandom
  # early maturity
  matrix[10,5] <- Survival(Sim5, seSim)$RealRandom * Survival(EM)$RealRandom
  matrix[6,5] <- Survival(Sim5, seSim)$RealRandom * (1-Survival(EM)$RealRandom)
  # immature survival each year
  matrix[2,1] <- Survival(Sh)$RealRandom
  matrix[3,2] <- Survival(Sim2, seSim)$RealRandom
  matrix[4,3] <- Survival(Sim3, seSim)$RealRandom
  matrix[5,4] <- Survival(Sim4, seSim)$RealRandom
  matrix[7,6] <- Survival(Sim6, seSim)$RealRandom
  matrix[8,7] <- Survival(Sim7, seSim)$RealRandom
  matrix[9,8] <- Survival(Sim8, seSim)$RealRandom
  matrix[10,9] <- Survival(Sim9, seSim)$RealRandom
  # adult survival
  matrix[10,10] <- Survival(Sa, seSa)$RealRandom
  
## Wild males ##
  # early maturity
  matrix[29,20] <- Survival(Sim9, seSim)$RealRandom * Survival(EM)$RealRandom
  matrix[21,20] <- Survival(Sim9, seSim)$RealRandom * (1-Survival(EM)$RealRandom)
  # immature survival each year
  matrix[12,11] <- Survival(Sh)$RealRandom
  matrix[13,12] <- Survival(Sim2, seSim)$RealRandom
  matrix[14,13] <- Survival(Sim3, seSim)$RealRandom
  matrix[15,14] <- Survival(Sim4, seSim)$RealRandom
  matrix[16,15] <- Survival(Sim5, seSim)$RealRandom
  matrix[17,16] <- Survival(Sim6, seSim)$RealRandom
  matrix[18,17] <- Survival(Sim7, seSim)$RealRandom
  matrix[19,18] <- Survival(Sim8, seSim)$RealRandom
  matrix[20,19] <- Survival(Sim9, seSim)$RealRandom
  matrix[22,21] <- Survival(Sim9, seSim)$RealRandom
  matrix[23,22] <- Survival(Sim9, seSim)$RealRandom
  matrix[24,23] <- Survival(Sim9, seSim)$RealRandom
  matrix[25,24] <- Survival(Sim9, seSim)$RealRandom
  matrix[26,25] <- Survival(Sim9, seSim)$RealRandom
  matrix[27,26] <- Survival(Sim9, seSim)$RealRandom
  matrix[28,27] <- Survival(Sim9, seSim)$RealRandom
  matrix[29,28] <- Survival(Sim9, seSim)$RealRandom
  # adult survival
  matrix[29,29] <- Survival(Sa, seSa)$RealRandom
  
## FEMALE HEADSTARTS ##
  matrix[32,30] <- Survival(releaseSY1r, sereleaseSY1r)$RealRandom
  matrix[32,31] <- Survival(releaseSY1n, sereleaseSY1n)$RealRandom
  matrix[33,32] <- Survival(releaseSY2, sereleaseSY2)$RealRandom
  matrix[34,33] <- Survival(releaseSY3, sereleaseSY3)$RealRandom
  matrix[35,34] <- Survival(releaseSY4)$RealRandom
  matrix[10,35] <- Survival(releaseSY5)$RealRandom*Survival(EM)$RealRandom
  matrix[36,35] <- Survival(releaseSY5)$RealRandom*(1-Survival(EM)$RealRandom)
  matrix[37,36] <- Survival(releaseSY6)$RealRandom
  matrix[38,37] <- Survival(releaseSY7)$RealRandom
  matrix[39,38] <- Survival(releaseSY7)$RealRandom
  matrix[10,39] <- Survival(releaseSY7)$RealRandom
  
## MALE HEADSTARTS ##
  matrix[42,40] <- Survival(releaseSY1r, sereleaseSY1r)$RealRandom
  matrix[42,41] <- Survival(releaseSY1n, sereleaseSY1n)$RealRandom
  matrix[43,42] <- Survival(releaseSY2, sereleaseSY2)$RealRandom
  matrix[44,43] <- Survival(releaseSY3, sereleaseSY2)$RealRandom
  matrix[45,44] <- Survival(releaseSY4)$RealRandom
  matrix[46,45] <- Survival(releaseSY5)$RealRandom
  matrix[47,46] <- Survival(releaseSY6)$RealRandom
  matrix[48,47] <- Survival(releaseSY7)$RealRandom
  matrix[49,48] <- Survival(releaseSY7)$RealRandom
  matrix[50,49] <- Survival(releaseSY7)$RealRandom
  matrix[51,50] <- Survival(releaseSY7)$RealRandom*Survival(EM)$RealRandom
  matrix[29,50] <- Survival(releaseSY7)$RealRandom*(1-Survival(EM)$RealRandom)
  matrix[52,51] <- Survival(releaseSY7)$RealRandom
  matrix[53,52] <- Survival(releaseSY7)$RealRandom
  matrix[54,53] <- Survival(releaseSY7)$RealRandom
  matrix[55,54] <- Survival(releaseSY7)$RealRandom
  matrix[56,55] <- Survival(releaseSY7)$RealRandom
  matrix[57,56] <- Survival(releaseSY7)$RealRandom
  matrix[58,57] <- Survival(releaseSY7)$RealRandom
  matrix[29,58] <- Survival(releaseSY7)$RealRandom
  
  names <- c("1F","2F","3F","4F","5F","6F","7F","8F","9F","AF","1M","2M","3M","4M","5M","6M","7M","8M","9M","10M","11M","12M","13M","14M","15M","16M","17M","18M","AM","hs1rF","hs1nF","hs2F","hs3F","hs4F","hs5F","hs6F","hs7F","hs8F","hs9F","hs1rM","hs1nM","hs2M","hs3M","hs4M","hs5M","hs6M","hs7M","hs8M","hs9M","hs10M","hs11M","hs12M","hs13M","hs14M","hs15M","hs16M","hs17M","hs18M")
  rownames(matrix) <- names
  colnames(matrix) <- names
  
  # no males = no reproduction
    if (pop[29] < 1){
      matrix[1,10] <- 0
      matrix[11,10] <- 0
    }


  # catastrophe
   if (runif(1) < catprob){
     matrix[1,10] <- matrix[1,10] * cateffect
     matrix[11,10] <- matrix[11,10] * cateffect
   }
  
  pop[is.na(pop)] <- 0
  

  # multiply the initial population by the transition matrix to get the population after a year  
  newpop <- matrix %*% matrix(pop)

  return(newpop[,1])
}



## Function for the stochastic model with releases ##
stochrun <- function(initpop, runs){ #INPUT- initial population vector, number of runs
  pop <- initpop
  
  #create dataframe for each years population output
  outputdf <- data.frame(matrix(0, nrow = runs + 1, ncol = length(pop)))
  colnames(outputdf) <- c("1F","2F","3F","4F","5F","6F","7F","8F","9F","AF","1M","2M","3M","4M","5M","6M","7M","8M","9M","10M","11M","12M","13M","14M","15M","16M","17M","18M","AM","hs1rF","hs1nF","hs2F","hs3F","hs4F","hs5F","hs6F","hs7F","hs8F","hs9F","hs1rM","hs1nM","hs2M","hs3M","hs4M","hs5M","hs6M","hs7M","hs8M","hs9M","hs10M","hs11M","hs12M","hs13M","hs14M","hs15M","hs16M","hs17M","hs18M")       #### change seq to initpop length -1 
  outputdf[1,] <- pop
  
  # create matrix for number of eggs collected
  collectedmatrix <- matrix(1, nrow=4)
  
   for (i in 1:runs){

     # run one step of the model
     pop <- stochstep(pop)


    collected <- (pop[1]+pop[11]) / (1-removal) / Survival(Sh)$RealRandom * removal
    #collected density dependence
    collected <- Reproduction(collected, secr)
    collected = pmin(collected, 800)

    # add to collected matrix
    collectedmatrix <- rbind(collectedmatrix, collected)
    collectedmatrix[is.na(collectedmatrix)] <- 0


    # Add releases to the population
    if (i <= 5){
      pop[30] <- Reproduction(currentreleases, secr) * Sc * releaseXR * riversplit * Survival(DR, seDR)$RealRandom
      pop[31] <- Reproduction(currentreleases, secr) * Sc * releaseXR * (1-riversplit) * Survival(DN, seDN)$RealRandom
      pop[40] <- Reproduction(currentreleases, secr) * Sc * (1-releaseXR) * riversplit * Survival(DR, seDR)$RealRandom
      pop[41] <- Reproduction(currentreleases, secr) * Sc * (1-releaseXR) * (1-riversplit) * Survival(DN, seDN)$RealRandom
    }

    else{
     pop[30] <- collectedmatrix[i] * releaseXR * Survival(Sc)$RealRandom * riversplit * Survival(DR, seDR)$RealRandom
     pop[31] <- collectedmatrix[i] * releaseXR * Survival(Sc)$RealRandom * (1-riversplit) * Survival(DN, seDN)$RealRandom
     pop[40] <- collectedmatrix[i] * (1-releaseXR) * Survival(Sc)$RealRandom * riversplit * Survival(DR, seDR)$RealRandom
     pop[41] <- collectedmatrix[i] * (1-releaseXR) * Survival(Sc)$RealRandom * (1-riversplit) * Survival(DN, seDN)$RealRandom
     }

    pop[is.na(pop)] <- 0
    
    pop
    outputdf[i + 1, ] <- pop

  }
  return(outputdf)
  }

# Run stochastic model + total population at each step + plotted
gg <- stochrun(currentpop, 50)
df <- as.data.frame(gg)
df$totalpop <- rowSums(df)
df$year <- 1:nrow(df)
ggplot(df, aes(y=totalpop, x=year)) + geom_line()


################################################################################
############MULTIPLE REPLICATES OF THE STOCHASTIC MODEL######################### 
################################################################################

# create blank list for the results
results <- list()
# set the number of replicates
reps <- 1000
# create a blank list for the total population at each step
population_values <- vector("list", reps)
# counter for the number of replicates with extinction (no adult males or no adult females)
count_adult <- 0

# set progress bar 
pb <- progress_bar$new(total=reps)
pb$tick(0)

# Run the stochastic model for the number of replicates
for (i in 1:reps) {
  
  # set the input population and number of runs
  initpop <- currentpop
  runs <- 50
  
  # progress bar
  pb$tick(1)
  
  # return results and population values
  output <- stochrun(initpop, runs)
  results[[i]] <- output
  population_values[[i]] <- rowSums(output)
  
  # Check for extinction at any step
  if (any(output[, "AM"] < 1) || any(output[, "AF"] < 1)) {
    count_adult <- count_adult + 1
  }
}


################################################################################
## Results ##

# Calculate the mean population at each step across all replicates
mean_population_at_steps <- sapply(1:runs, function(step) {
  mean_population <- mean(sapply(population_values, function(pop) pop[step]))
  return(mean_population)
})

# Calculate the standard deviation for each step
sd_population_at_steps <- sapply(1:runs, function(step) {
  sd <- sd(sapply(population_values, function(pop) pop[step]))
  return(sd)
})


# Calculate the standard error for each step
se_population_at_steps <- sapply(1:runs, function(step) {
  se <- sd(sapply(population_values, function(pop) pop[step])) / sqrt(reps)
  return(se)
})


# Calculate the 95% confidence intervals for each step
confidence_intervals <- sapply(1:runs, function(step) {
  mean_pop <- mean_population_at_steps[step]
  se <- se_population_at_steps[step]
  #margin_of_error <- qt(0.975, df = reps - 1) * se
  margin_of_error <- 1.96 * se
  
  lower_ci <- mean_pop - margin_of_error
  upper_ci <- mean_pop + margin_of_error
  return(c(lower_ci, upper_ci))
})

# Create a data frame with time step, mean population, standard deviation, standard error, and confidence intervals
population_stats_df <- data.frame(
  step = 1:runs,
  mean_population = mean_population_at_steps,
  sd_population = sd_population_at_steps,
  se_population = se_population_at_steps,
  lower_ci = confidence_intervals[1, ],
  upper_ci = confidence_intervals[2, ]
)

# 0-50 instead of 1-51
population_stats_df$step <- population_stats_df$step-1

############

## calculate extinction probability
# Calculate the percentage of repetitions where condition is met
percentage_adult <- (count_adult / reps) * 100
# Print the results
cat("Percentage of repetitions with less than 1 adult males or females:", percentage_adult, "%\n")

############

## plot the mean total population over time 
ggplot(population_stats_df, aes(x = step, y = mean_population)) +
  geom_line() +
  labs(x = "Step", y = "Mean Population", title = "Mean Population Over Steps")


## plot the mean total population over time with confidence intervals
plot(population_stats_df$step, population_stats_df$mean_population, type = "l", ylim = c(0, (max(population_stats_df$upper_ci)+550)),
     xlab = "Year", ylab = "Population Size")
lines(population_stats_df$step, population_stats_df$upper_ci, col = "red", lty = 2)
lines(population_stats_df$step, population_stats_df$lower_ci, col = "blue", lty = 2)
## add legend 
legend("topleft", legend = c("Mean Population", "Upper 95% CI", "Lower 95% CI", "Mean Adult Population"), 
       col = c("black", "red", "blue", "springgreen3"), lty = c(1, 2, 2,1))
## add total adult population at each step
lines(0:(runs), (output[["AF"]]+output[["AM"]]), col = "springgreen3")
## add sample of 25 replicates
for (i in 1:25) {
  total_pop <- rowSums(results[[i]])
  lines(0:(runs), total_pop, lwd = 0.5, col="lightgrey")
}
# make lines clearer
lines(population_stats_df$step, population_stats_df$mean_population, col = "black")
lines(population_stats_df$step, population_stats_df$upper_ci, col = "red", lty = 2)
lines(population_stats_df$step, population_stats_df$lower_ci, col = "blue", lty = 2)


# Calculate mean final population size
final_pop_sizes <- sapply(results, function(output) tail(rowSums(output), 1))
mean_final_pop <- mean(final_pop_sizes)
median_final_pop <- median(final_pop_sizes)

# calculate mean final adult population size 
final_Fpops <- sapply(results, function(output) tail(output[["AF"]], 1))
mean_Fpop <- mean(final_Fpops)
final_Mpops <- sapply(results, function(output) tail(output[["AM"]], 1))
mean_Mpop <- mean(final_Mpops)


# Print the results
cat("Mean final population size:", mean_final_pop, "\n", 
    "Mean final adult population size:", mean_Fpop, mean_Mpop,"\n",
    "Adult percentage:", (mean_Fpop+mean_Mpop)/mean_final_pop*100, "\n", 
    "Extinction risk:", percentage_adult, "\n")

################################################################################
############## SENSITIVITY ANALYSIS ############################################
################################################################################

## List of parameter values to test sensitivity # changed by 10%
full_parameter_values <- list(
  B = c((B-B*0.1), B, (B+B*0.1)),                                                               # birthrate
  Sa = c((Sa-Sa*0.1), Sa, (Sa+Sa*0.1)),                                                         # adult survival
  removal = c((removal-removal*0.1), removal, (removal+removal*0.1)),                           # removal for head-starting
  Sh = c((Sh-Sh*0.1), Sh, (Sh+Sh*0.1)),                                                         # hatchling survival
  Sc = c((Sc-Sc*0.1), Sc, (Sc+Sc*0.1)),                                                         # captivity survival
  DN = c((DN-DN*0.1), DN, (DN+DN*0.1)),                                                         # dispersal (narayani)
  releaseSY1n = c((releaseSY1n-releaseSY1n*0.1), releaseSY1n, (releaseSY1n+releaseSY1n*0.1)),   # postrelease survival (year 1 narayani)
  releaseSY1r = c((releaseSY1r-releaseSY1r*0.1), releaseSY1r, (releaseSY1r+releaseSY1r*0.1)),   # postrelease survival (year 1 rapti)
  releaseSY2 = c((releaseSY2-releaseSY2*0.1), releaseSY2, (releaseSY2+releaseSY2*0.1)),         # postrelease survival (year 2)
  EM = c((EM-EM*0.1), EM, (EM+EM*0.1)),                                                         # early maturity rate
  XR = c((XR-XR*0.1), XR, (XR+XR*0.1)),                                                         # wild sex ratio
  Sim2 = c((Sim2-Sim2*0.1), Sim2, (Sim2+Sim2*0.1)),                                             # immature survival (year 1)
  Sim3 = c((Sim3-Sim3*0.1), Sim3, (Sim3+Sim3*0.1)),                                             # immature survival (year 2)
  releaseXR = c((releaseXR-releaseXR*0.1), releaseXR, (releaseXR+releaseXR*0.1)),               # release sex ratio
  riversplit = c((riversplit-riversplit*0.1), riversplit, (riversplit+riversplit*0.1)),         # release site proportion
  DR = c((DR-DR*0.1), DR, (DR+DR*0.1))                                                          # dispersal (rapti)
  )

# Create dataframe for results of sensitivity analysis (upper, lower baseline values)
ubl_df <- data.frame(full_parameter_values)
row.names(ubl_df) <- c("lower", "baseline", "upper")

##############################################################################
### Adult determintistic  ###
Sasensitivity_values <- c((Sa-Sa*0.1), Sa, (0.99))       # set values to test
results_sensitivity <- list()                            # blank lists for results
population_values <- list()

# loop over the sensitivity values
for (i in 1:length(Sasensitivity_values)) {
  # set the sensitivity value for adult survival
  Sa <- Sasensitivity_values[i]

  # run the simulation for the current sensitivity value
  initpop <- currentpop
  runs <- 50

  output <- run(initpop, runs)
  # final total population
  population_values[[i]] <- rowSums(output)
  
}

SADF <- as.data.frame(do.call(cbind, population_values))
colnames(SADF) <- Sasensitivity_values
SADF$run <- 1:(runs+1)

# add results to dataframe
ubl_df$Sa <- c(SADF[51,1], SADF[51,2], SADF[51,3])
# reset parameter
Sa<-0.95

################################################################################
################## stochastic version Sa sensitivity ###########################
### Sa 0.1-0.9 (0.95) ###
Sasensitivity_values <- c(seq(0.1, 0.9, 0.1), 0.71, 0.95, 0.99)  # set values to test
results_sensitivity <- list()                                    # blank lists for results
population_values <- list()
reps_per_sensitivity <- 10
results <- list()
runs <-100

# Initialize a progress bar
pb <- progress_bar$new(total = length(Sasensitivity_values))

# Loop over the sensitivity values
for (i in 1:length(Sasensitivity_values)) {
  pb$tick(1)
  
  Sa <- Sasensitivity_values[i]
  
  for (j in 1:reps_per_sensitivity) {
    initpop <- currentpop
    runs <- 50
    
    Sa <- Sasensitivity_values[i]
    
    output <- stochrun(initpop, runs)
    
    population_values[[j]] <- rowSums(output)
    mean_population_df <- data.frame(Sasensitivity_values = rowMeans(sapply(population_values, "[", 1:runs)))
    
  }
  results[[i]] <- mean_population_df
  
}

sSaSADF <- as.data.frame(do.call(cbind, results))
colnames(sSaSADF) <- Sasensitivity_values
sSaSADF$run <- 1:(runs)

#ubl
ubl_df$Sa <- c(sSaSADF[50,10], sSaSADF[50,11], sSaSADF[50,12])
Sa<-0.95

##############################################################################
### Hatchling survival ###
Shsensitivity_values <- c((Sh-Sh*0.1), Sh, (Sh+Sh*0.1))      # set values to test
 
# empty list to store the results for each sensitivity value
results_sensitivity <- list()
Shpopulation_values <- list()

# Set the number of runs for each sensitivity value
# runs_per_sensitivity <- 10

# Initialize a progress bar
pb <- progress_bar$new(total = length(Shsensitivity_values))

# Loop over the sensitivity values
for (i in 1:length(Shsensitivity_values)) {
  # Set the sensitivity value for adult survival
  Sh <- Shsensitivity_values[i]
  
  # Run the simulation for the current sensitivity value
  initpop <- currentpop
  runs <- 50
  Sh <- Shsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  Shpopulation_values[[i]] <- rowSums(output)
  
}

ShSADF <- as.data.frame(do.call(cbind, Shpopulation_values))
colnames(ShSADF) <- Shsensitivity_values
ShSADF$run <- 1:(runs+1)

#ubl
ubl_df$Sh <- c(ShSADF[51,1], ShSADF[51,2], ShSADF[51,3])
Sh <- 0.1

#### repeat for all parameters ####
#################################################################################
### Dispersal Narayani ###
DNsensitivity_values <- c((DN-DN*0.1), DN, (DN+DN*0.1))
DNresults_sensitivity <- list()
DNpopulation_values <- list()
pb <- progress_bar$new(total = length(DNsensitivity_values))

for (i in 1:length(DNsensitivity_values)) {
  DN <- DNsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  DN <- DNsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  DNpopulation_values[[i]] <- rowSums(output)
  
}

DNSADF <- as.data.frame(do.call(cbind, DNpopulation_values))
colnames(DNSADF) <- DNsensitivity_values
DNSADF$run <- 1:(runs+1)
#ubl
ubl_df$DN <- c(DNSADF[51,1], DNSADF[51,2], DNSADF[51,3])
DN <- (1-0.45)

#################################################################################
### y1 headstart survival narayani  ###
releaseSY1nsensitivity_values <- c((releaseSY1n-releaseSY1n*0.1), releaseSY1n, (releaseSY1n+releaseSY1n*0.1))
releaseSY1nresults_sensitivity <- list()
releaseSY1npopulation_values <- list()
pb <- progress_bar$new(total = length(releaseSY1nsensitivity_values))

for (i in 1:length(releaseSY1nsensitivity_values)) {
  releaseSY1n <- releaseSY1nsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  releaseSY1n <- releaseSY1nsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  releaseSY1npopulation_values[[i]] <- rowSums(output)
  
}

releaseSY1nSADF <- as.data.frame(do.call(cbind, releaseSY1npopulation_values))
colnames(releaseSY1nSADF) <- releaseSY1nsensitivity_values
releaseSY1nSADF$run <- 1:(runs+1)

ubl_df$releaseSY1n <- c(releaseSY1nSADF[51,1], releaseSY1nSADF[51,2], releaseSY1nSADF[51,3])
releaseSY1n <- 0.3

#################################################################################
### y1 headstart survival rapti ###
releaseSY1rsensitivity_values <- c((releaseSY1r-releaseSY1r*0.1), releaseSY1r, (releaseSY1r+releaseSY1r*0.1))
releaseSY1rresults_sensitivity <- list()
releaseSY1rpopulation_values <- list()
pb <- progress_bar$new(total = length(releaseSY1rsensitivity_values))

for (i in 1:length(releaseSY1rsensitivity_values)) {
  releaseSY1r <- releaseSY1rsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  releaseSY1r <- releaseSY1rsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  releaseSY1rpopulation_values[[i]] <- rowSums(output)
  
}

releaseSY1rSADF <- as.data.frame(do.call(cbind, releaseSY1rpopulation_values))
colnames(releaseSY1rSADF) <- releaseSY1rsensitivity_values
releaseSY1rSADF$run <- 1:(runs+1)
#ubl
ubl_df$releaseSY1r <- c(releaseSY1rSADF[51,1], releaseSY1rSADF[51,2], releaseSY1rSADF[51,3])
releaseSY1r <- 0.6

#########################################
### y2 headstart survival ###
releaseSY2sensitivity_values <- c((releaseSY2-releaseSY2*0.1), releaseSY2, (releaseSY2+releaseSY2*0.1))
releaseSY2results_sensitivity <- list()
releaseSY2population_values <- list()
pb <- progress_bar$new(total = length(releaseSY2sensitivity_values))

for (i in 1:length(releaseSY2sensitivity_values)) {
  releaseSY2 <- releaseSY2sensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  releaseSY2 <- releaseSY2sensitivity_values[i]
  
  
  output <- run(initpop, runs)
  releaseSY2population_values[[i]] <- rowSums(output)
  
}

releaseSY2SADF <- as.data.frame(do.call(cbind, releaseSY2population_values))
colnames(releaseSY2SADF) <- releaseSY2sensitivity_values
releaseSY2SADF$run <- 1:(runs+1)
#ubl
ubl_df$releaseSY2 <- c(releaseSY2SADF[51,1], releaseSY2SADF[51,2], releaseSY2SADF[51,3])
releaseSY2 <- 0.63

############################################
### Early Maturity ###
EMsensitivity_values <- c((EM-EM*0.1), EM, (EM+EM*0.1))
EMresults_sensitivity <- list()
EMpopulation_values <- list()
pb <- progress_bar$new(total = length(EMsensitivity_values))

for (i in 1:length(EMsensitivity_values)) {
  EM <- EMsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  EM <- EMsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  EMpopulation_values[[i]] <- rowSums(output)
  
}

EMSADF <- as.data.frame(do.call(cbind, EMpopulation_values))
colnames(EMSADF) <- EMsensitivity_values
EMSADF$run <- 1:(runs+1)

ubl_df$EM <- c(EMSADF[51,1], EMSADF[51,2], EMSADF[51,3])
EM <- 0.5

####################################################################
### wild sex ratio ###
XRsensitivity_values <- c((XR-XR*0.1), XR, (XR+XR*0.1))
XRresults_sensitivity <- list()
XRpopulation_values <- list()
pb <- progress_bar$new(total = length(XRsensitivity_values))

for (i in 1:length(XRsensitivity_values)) {
  XR <- XRsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  XR <- XRsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  XRpopulation_values[[i]] <- rowSums(output)
  
}

XRSADF <- as.data.frame(do.call(cbind, XRpopulation_values))
colnames(XRSADF) <- XRsensitivity_values
XRSADF$run <- 1:(runs+1)

ubl_df$XR <- c(XRSADF[51,1], XRSADF[51,2], XRSADF[51,3])
XR <- 0.75

###################################################################
### immature survival  age 2 ###
Sim2sensitivity_values <- c((Sim2-Sim2*0.1), Sim2, (Sim2+Sim2*0.1))
Sim2results_sensitivity <- list()
Sim2population_values <- list()

for (i in 1:length(Sim2sensitivity_values)) {
  Sim2 <- Sim2sensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  Sim2 <- Sim2sensitivity_values[i]
  
  
  output <- run(initpop, runs)
  Sim2population_values[[i]] <- rowSums(output)
  
}

Sim2SADF <- as.data.frame(do.call(cbind, Sim2population_values))
colnames(Sim2SADF) <- Sim2sensitivity_values
Sim2SADF$run <- 1:(runs+1)

ubl_df$Sim2 <- c(Sim2SADF[51,1], Sim2SADF[51,2], Sim2SADF[51,3])
Sim2 <- 0.2
###################################################################
### immature survival age 3 ###
Sim3sensitivity_values <- c((Sim3-Sim3*0.1), Sim3, (Sim3+Sim3*0.1))
Sim3results_sensitivity <- list()
Sim3population_values <- list()

for (i in 1:length(Sim3sensitivity_values)) {
  Sim3 <- Sim3sensitivity_values[i]
  initpop <- currentpop
  runs <- 50

  output <- run(initpop, runs)
  Sim3population_values[[i]] <- rowSums(output)
  
}

Sim3SADF <- as.data.frame(do.call(cbind, Sim3population_values))
colnames(Sim3SADF) <- Sim3sensitivity_values
Sim3SADF$run <- 1:(runs+1)

ubl_df$Sim3 <- c(Sim3SADF[51,1], Sim3SADF[51,3], Sim3SADF[51,3])
Sim3 <- 0.3

#####################################################################
### headstart sex ratio ###
releaseXRsensitivity_values <- c((releaseXR-releaseXR*0.1), releaseXR, (releaseXR+releaseXR*0.1))

# empty list to store the results for each sensitivity value
releaseXRresults_sensitivity <- list()
releaseXRpopulation_values <- list()

# Loop over the sensitivity values
for (i in 1:length(releaseXRsensitivity_values)) {
  # Set the sensitivity value for adult survival
  releaseXR <- releaseXRsensitivity_values[i]
  
  # Run the simulation for the current sensitivity value
  initpop <- currentpop
  runs <- 50
  releaseXR <- releaseXRsensitivity_values[i]

  output <- run(initpop, runs)
  releaseXRpopulation_values[[i]] <- rowSums(output)
}

releaseXRSADF <- as.data.frame(do.call(cbind, releaseXRpopulation_values))
colnames(releaseXRSADF) <- releaseXRsensitivity_values
releaseXRSADF$run <- 1:(runs+1)

ubl_df$releaseXR <- c(releaseXRSADF[51,1], releaseXRSADF[51,2], releaseXRSADF[51,3])
releaseXR <- 0.8
################################################################################
### release site ###
riversplitsensitivity_values <- c((riversplit-riversplit*0.1), riversplit, (riversplit+riversplit*0.1))

# empty list to store the results for each sensitivity value
riversplitpopulation_values <- list()

# Loop over the sensitivity values
for (i in 1:length(riversplitsensitivity_values)) {
  # Set the sensitivity value for adult survival
  riversplit <- riversplitsensitivity_values[i]
  
  # Run the simulation for the current sensitivity value
  initpop <- currentpop
  runs <- 50
  riversplit <- riversplitsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  riversplitpopulation_values[[i]] <- rowSums(output)
}

riversplitSADF <- as.data.frame(do.call(cbind, riversplitpopulation_values))
colnames(riversplitSADF) <- riversplitsensitivity_values
riversplitSADF$run <- 1:(runs+1)

ubl_df$riversplit <- c(riversplitSADF[51,1], riversplitSADF[51,2], riversplitSADF[51,3])
riversplit <- 0.8

################################################################################
### dispersal survival rapti  ###
DRsensitivity_values <- c((DR-DR*0.1), DR, (DR+DR*0.1))

# empty list to store the results for each sensitivity value
DRpopulation_values <- list()

# Loop over the sensitivity values
for (i in 1:length(DRsensitivity_values)) {
  # Set the sensitivity value for adult survival
  DR <- DRsensitivity_values[i]
  
  # Run the simulation for the current sensitivity value
  initpop <- currentpop
  runs <- 50
  DR <- DRsensitivity_values[i]
  
  output <- run(initpop, runs)
  DRpopulation_values[[i]] <- rowSums(output)
  
}

DRSADF <- as.data.frame(do.call(cbind, DRpopulation_values))
colnames(DRSADF) <- DRsensitivity_values
DRSADF$run <- 1:(runs+1)

ubl_df$DR <- c(DRSADF[51,1], DRSADF[51,2], DRSADF[51,3])
DR <- (1-0.15)


################################################################################
### Removal ###
removalsensitivity_values <- c((removal-removal*0.1), removal, (removal+removal*0.1))
# 
# Initialize an empty list to store the results for each sensitivity value
population_values <- list()

# Loop over the sensitivity values
for (i in 1:length(removalsensitivity_values)) {
  # Set the sensitivity value for removal
  removal <- removalsensitivity_values[i]
  initpop <- currentpop
  runs <- 50

  output <- run(initpop, runs)
  population_values[[i]] <- rowSums(output)
}

removalSADF <- as.data.frame(do.call(cbind, population_values))
colnames(removalSADF) <- removalsensitivity_values
removalSADF$run <- 1:(runs+1)

ubl_df$removal <- c(removalSADF[51,1], removalSADF[51,2], removalSADF[51,3])
removal <- 0.9
################################################################################
### Captivity survival ###
Scsensitivity_values <- c((Sc-Sc*0.1), Sc, (Sc+Sc*0.1))

# Initialize an empty list to store the results for each sensitivity value
Scpopulation_values <- list()

# Loop over the sensitivity values
for (i in 1:length(Scsensitivity_values)) {
  # Set the sensitivity value for adult survival
  Sc <- Scsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  
  output <- run(initpop, runs)
  Scpopulation_values[[i]] <- rowSums(output)
}

ScSADF <- as.data.frame(do.call(cbind, Scpopulation_values))
colnames(ScSADF) <- Scsensitivity_values
ScSADF$run <- 1:(runs+1)

ubl_df$Sc <- c(ScSADF[51,1], ScSADF[51,2], ScSADF[51,3])
Sc <- 0.2
################################################################################
### Birth rate ###
Bsensitivity_values <- c((B-B*0.1), B, (B+B*0.1))
# empty list to store the results for each sensitivity value
population_values <- list()

# Loop over the sensitivity values
for (i in 1:length(Bsensitivity_values)) {
  # Set the sensitivity value for adult survival
  B <- Bsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  output <- run(initpop, runs)
  population_values[[i]] <- rowSums(output)
}

BSADF <- as.data.frame(do.call(cbind, population_values))
colnames(BSADF) <- Bsensitivity_values
BSADF$run <- 1:(runs+1)

ubl_df$B <- c(BSADF[51,1], BSADF[51,2], BSADF[51,3])
B <- 0.8*39

################################################################################

# Format the dataframe for plotting
ubl_df <- t(ubl_df)
ubl_df <- as.data.frame(ubl_df)
parameters <- c("Birth Rate", "Adult Survival", "Removal", "Hatchling Survival", "Captivity Survival",
                "Dispersal (Narayani)", "Year 1 HS Survival (Narayani)", "Year 1 HS Survival (Rapti)", 
                "Year 2 HS Survival",
                "Early Maturity", "Wild Sex Ratio",
                "Immature Survival (Year 1)", "Immature Survival (Year 2)", "Release Sex Ratio",
                "Release Site", "Dispersal (Rapti)")
ubl_df$parameter <- parameters
ubl_df$range <- ubl_df$upper-ubl_df$lower
# convert absolute population value to percentage change from baseline
ubl_df$lowchange <- -(ubl_df$baseline-ubl_df$lower)/ubl_df$baseline*100
ubl_df$highchange <- (ubl_df$upper-ubl_df$baseline)/ubl_df$baseline*100
  
  
# Plot the sensitivity analysis
ubl_df %>%
  arrange(range) %>%    
  mutate(parameter=factor(parameter, levels=parameter)) %>% 
  ggplot( aes(x=parameter, y=(-100:100))) +
  geom_segment( aes(x=parameter, xend=parameter, y=0, yend=highchange), color="black") +
  geom_segment( aes(x=parameter, xend=parameter, y=lowchange, yend=0), color="black") +
  geom_point( aes(x=parameter, y=lowchange), color=rgb(0.7,0.2,0.1,0.5), size=4 ) +
  geom_point( aes(x=parameter, y=highchange), color=rgb(0.2,0.7,0.1,0.5), size=4) +
  coord_flip() +
  theme_ipsum() +
  theme(legend.position = "none", )+
  geom_abline(aes(intercept=0, slope=0))+
  xlab("Parameter") +
  ylab("Change in Population Size (%)") 




##################### Further sensitivity analysis #############################
### Adult survival ###

Sasensitivity_values <- c(seq(0.1,0.9,0.1),0.95)            # set values to test
Sapopulation_values <- list()                               # blank list for results

for (i in 1:length(Sasensitivity_values)) {                 # run model
  Sa <- Sasensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  Sa <- Sasensitivity_values[i]*
  
  
  output <- run(initpop, runs)
  Sapopulation_values[[i]] <- rowSums(output)
  
}

# format data
SaSADF <- as.data.frame(do.call(cbind, Sapopulation_values))
colnames(SaSADF) <- Sasensitivity_values
SaSADF$run <- 1:(runs+1)

# colour scheme
color_palette <- viridis(9, "plasma")

# plot !
ggplot(SaSADF, aes(x = run)) +
  geom_line(aes(y = `0.95`), color = "black", lty=2, size = 1.2) +
  geom_line(aes(y = `0.9`), color = color_palette[1], size = 1.2) +
  geom_line(aes(y = `0.8`), color = color_palette[2], size = 1.2) +
  geom_line(aes(y = `0.7`), color = color_palette[3], size = 1.2) +
  geom_line(aes(y = `0.6`), color = color_palette[4], size = 1.2) +
  geom_line(aes(y = `0.5`), color = color_palette[5], size = 1.2) +
  geom_line(aes(y = `0.4`), color = color_palette[6], size = 1.2) +
  geom_line(aes(y = `0.3`), color = color_palette[7], size = 1.2) +
  geom_line(aes(y = `0.2`), color = color_palette[8], size = 1.2) +
  geom_line(aes(y = `0.1`), color = color_palette[9], size = 1.2) +
  ylim(0, 650) +
  labs(x = "Year", y = "Population Size", title = "Adult Survival") +
  theme_classic() +
  scale_color_viridis_d("plasma")

# reset value
Sa <- 0.95

# Repeat for other parameters \/ #
################################################################################
### Captivity survival ###

Scsensitivity_values <- seq(0.1,0.9,0.1)
Scresults_sensitivity <- list()
Scpopulation_values <- list()

for (i in 1:length(Scsensitivity_values)) {
  Sc <- Scsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  Sc <- Scsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  Scpopulation_values[[i]] <- rowSums(output)
  
}

ScSADF <- as.data.frame(do.call(cbind, Scpopulation_values))
colnames(ScSADF) <- Scsensitivity_values
ScSADF$run <- 1:(runs+1)
ggplot(ScSADF, aes(x = run)) +
  geom_line(aes(y = `0.9`), color = color_palette[1], size = 1.2) +
  geom_line(aes(y = `0.8`), color = color_palette[2], size = 1.2) +
  geom_line(aes(y = `0.7`), color = color_palette[3], size = 1.2) +
  geom_line(aes(y = `0.6`), color = color_palette[4], size = 1.2) +
  geom_line(aes(y = `0.5`), color = color_palette[5], size = 1.2) +
  geom_line(aes(y = `0.4`), color = color_palette[6], size = 1.2) +
  geom_line(aes(y = `0.3`), color = color_palette[7], size = 1.2) +
  geom_line(aes(y = `0.2`), color = "Black", lty=2, size = 1.2) +
  geom_line(aes(y = `0.1`), color = color_palette[9], size = 1.2) +
  ylim(0, 3000) +
  labs(x = "Year", y = "Population Size", title = "Survival in Captivity") +
  theme_classic() 

Sc <- 0.2
################################################################################
### removal ###

removalsensitivity_values <- seq(0.1,0.9,0.1)
removalresults_sensitivity <- list()
removalpopulation_values <- list()

for (i in 1:length(removalsensitivity_values)) {
  removal <- removalsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  removal <- removalsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  removalpopulation_values[[i]] <- rowSums(output)
  
}

removalSADF <- as.data.frame(do.call(cbind, removalpopulation_values))
colnames(removalSADF) <- removalsensitivity_values
removalSADF$run <- 1:(runs+1)
ggplot(removalSADF, aes(x = run)) +
  geom_line(aes(y = `0.9`), color = "Black", lty=2, size = 1.2) +
  geom_line(aes(y = `0.8`), color = color_palette[2], size = 1.2) +
  geom_line(aes(y = `0.7`), color = color_palette[3], size = 1.2) +
  geom_line(aes(y = `0.6`), color = color_palette[4], size = 1.2) +
  geom_line(aes(y = `0.5`), color = color_palette[5], size = 1.2) +
  geom_line(aes(y = `0.4`), color = color_palette[6], size = 1.2) +
  geom_line(aes(y = `0.3`), color = color_palette[7], size = 1.2) +
  geom_line(aes(y = `0.2`), color = color_palette[8], size = 1.2) +
  geom_line(aes(y = `0.1`), color = color_palette[9], size = 1.2) +
  ylim(0, 1150) +
  labs(x = "Year", y = "Population Size", title = "Removal") +
  theme_classic() +
  theme(legend.position = "top", legend.direction = "horizontal")

removal <- 0.9
################################################################################
### dispersal survival ###

DRsensitivity_values <- c(seq(0.1,0.9,0.1), 0.85)
DRresults_sensitivity <- list()
DRpopulation_values <- list()

for (i in 1:length(DRsensitivity_values)) {
  DR <- DRsensitivity_values[i]
  initpop <- currentpop
  runs <- 50
  DR <- DRsensitivity_values[i]
  
  
  output <- run(initpop, runs)
  DRpopulation_values[[i]] <- rowSums(output)
  
}

DRSADF <- as.data.frame(do.call(cbind, DRpopulation_values))
colnames(DRSADF) <- DRsensitivity_values
DRSADF$run <- 1:(runs+1)
ggplot(DRSADF, aes(x = run)) +
  geom_line(aes(y = `0.9`), color = color_palette[1], size = 1.2) +
  geom_line(aes(y = `0.85`), color = "Black", lty=2, size = 1.2) +
  geom_line(aes(y = `0.8`), color = color_palette[2], size = 1.2) +
  geom_line(aes(y = `0.7`), color = color_palette[3], size = 1.2) +
  geom_line(aes(y = `0.6`), color = color_palette[4], size = 1.2) +
  geom_line(aes(y = `0.5`), color = color_palette[5], size = 1.2) +
  geom_line(aes(y = `0.4`), color = color_palette[6], size = 1.2) +
  geom_line(aes(y = `0.3`), color = color_palette[7], size = 1.2) +
  geom_line(aes(y = `0.2`), color = color_palette[8], size = 1.2) +
  geom_line(aes(y = `0.1`), color = color_palette[9], size = 1.2) +
  ylim(0, 700) +
  labs(x = "Year", y = "Population Size", title = "Dispersal Survival (Rapti)") +
  theme_classic() 

DR <- 0.85
#################
# legend for sensitivity analysis plots
legend("topleft", legend = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","Baseline"), 
       col = c(color_palette[9], color_palette[8], color_palette[7], color_palette[6], color_palette[5], 
               color_palette[4], color_palette[3], color_palette[2],color_palette[1],
               "black"), lty = c(rep(1, 9), 2), horiz = TRUE, lwd = 5, title = "Parameter Value")
################################################################################
####################### HEATMAPs ################################################
################################################################################
par(mfrow=c(2,2))

################################################################################
############# hatchling survival heat map ######################################

# create blank lists for results
heat_population_values <- list()
heat_resultSh <- list()
# set removal and immature survival valued to be tested
Shsensitivity_values <- seq(0.1, 0.9, 0.1)
removalsensitivity_values <- seq(0.1, 0.9, 0.1)

#run with each removal value against each Sh value
for (i in 1:length(removalsensitivity_values)) {
  for (j in 1:length(Shsensitivity_values)) {
    removal <- removalsensitivity_values[i]
    Sh <- Shsensitivity_values[j]
    
    initpop <- currentpop
    runs <- 50
    
    output <- run(initpop, runs)
    heat_population_values[[j]] <- rowSums(output)
  }
  heat_resultSh[[i]] <- sapply(heat_population_values, tail,1)
}

# plot with base r
heat_matrixSh <- do.call(rbind, heat_resultSh)
rownames(heat_matrixSh) <- removalsensitivity_values
colnames(heat_matrixSh) <- Shsensitivity_values
heatmap(heat_matrixSh, Colv = NA, Rowv = NA, xlab = "Hatchling Survival", ylab = "Removal") 
heatmap(heat_matrixSh, Colv = NA, Rowv = NA, xlab = "Hatchling Survival", ylab = "Removal", revC = TRUE) 

                                             
# plot with ggplot
  # reorder 
Sh_melt <- melt(heat_matrixSh)    
  # plot
ggplot(Sh_melt, aes(X1, X2)) +                           
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "pink", high = "orangered4",
                       name = "Population size") +
  labs(x="Removal", y="Hatchling Survival") +
  theme_minimal() +
  scale_x_discrete(name ="Removal", 
                     limits=seq(0.1,0.9,0.1)) +
  scale_y_discrete(name ="Hatchling Survival", 
                   limits=seq(0.1,0.9,0.1)) +
  theme(axis.text.x = element_text(vjust = 0.5),  # Adjust vertical text placement
        axis.text.y = element_text(hjust = 0.5))
  
 Sh <- 0.1
 removal <-0.9
 
################################################################################
############# Immature survival (year1) heat map ###############################

# create blank lists for results
 heat_population_values <- list()
 heat_resultSim2 <- list()
# set removal and immature survival valued to be tested
 removalsensitivity_values <- seq(0.1, 0.9, 0.1)
 Sim2sensitivity_values <- seq(0.1, 0.9, 0.1)
 
#run with each removal value against each Sim2 value
 for (i in 1:length(removalsensitivity_values)) {
   for (j in 1:length(Sim2sensitivity_values)) {
     removal <- removalsensitivity_values[i]
     Sim2 <- Sim2sensitivity_values[j]
     
     initpop <- currentpop
     runs <- 50
     
     output <- run(initpop, runs)
     heat_population_values[[j]] <- rowSums(output)
   }
   heat_resultSim2[[i]] <- sapply(heat_population_values, tail,1)
 }
 
# plot with base r
 heat_matrixSim2 <- do.call(rbind, heat_resultSim2)
 rownames(heat_matrixSim2) <- removalsensitivity_values
 colnames(heat_matrixSim2) <- Sim2sensitivity_values
 heatmap(heat_matrixSim2, Colv = NA, Rowv = NA, xlab = "Immature Survival (Age 2)", ylab = "Removal")
 

# plot with ggplot
 # reorder data
Sim2_melt <- melt(heat_matrixSim2)                                           
 # plot
ggplot(Sim2_melt, aes(X1, X2)) +                     
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "pink", high = "orangered4",
                      name = "Population size") +
  labs(x="Removal", y="Immature Survival (Age 2)") +
  theme_minimal() +
  scale_x_discrete(name ="Removal", 
                   limits=seq(0.1,0.9,0.1)) +
  scale_y_discrete(name ="Immature Survival (Age 2)", 
                   limits=seq(0.1,0.9,0.1)) +
  theme(axis.text.x = element_text(vjust = 0.5),  
        axis.text.y = element_text(hjust = 0.5))

 Sim2 <- 0.2
 removal <- 0.9
 
################################################################################
######### catastrophe sensitivity using heatmap ################################
 
# create blank lists for the results 
 heat_population_values <- list()
 heat_resultcat <- list()
# set catastrophe sensitivity values for probability and effect
 catprobsensitivity_values <- seq(0.1, 0.9, 0.05)
 cateffectsensitivity_values <- seq(0.05, 0.13, 0.005)
 
# run the model for each catastrophe probability value against each catastrophe effect value
 for (i in 1:length(catprobsensitivity_values)) {
   for (j in 1:length(cateffectsensitivity_values)) {
     catprob <- catprobsensitivity_values[i]
     cateffect <- cateffectsensitivity_values[j]
     
     initpop <- currentpop
     runs <- 50
     
     output <- run(initpop, runs)
     heat_population_values[[j]] <- rowSums(output)
   }
   heat_resultcat[[i]] <- sapply(heat_population_values, tail,1)
 }
 
# plot using base R
 heat_matrixcat <- do.call(rbind, heat_resultcat)
 rownames(heat_matrixcat) <- catprobsensitivity_values
 colnames(heat_matrixcat) <- cateffectsensitivity_values
 heatmap(heat_matrixcat, Colv = NA, Rowv = NA, xlab = "Catastrophe Probability", ylab = "Catastrophe Effect (Sh)") #, col=brewer.pal(n = 9, name = "OrRd")
 
# plot using ggplot
 # reorder data
 cat_melt <- melt(heat_matrixcat)  
 # plot
 ggplot(cat_melt, aes(X1, X2)) +        
   geom_tile(aes(fill = value)) + 
   scale_fill_gradient2(low = "pink", high = "orangered4",
                        name = "Population size") +
   labs(x="Catastrophe Probability", y="Catastrophe Effect (Hatchlings)") +
   theme_minimal() 
   scale_x_discrete(name ="Catastrophe Probability",
                    limits=seq(0.1,0.9,0.1)) +
   scale_y_discrete(name ="Catastrophe Effect (Sh)",
                    limits=seq(0.05, 0.15, 0.01)) +
   theme(axis.text.x = element_text(vjust = 0.5),
         axis.text.y = element_text(hjust = 0.5))

catprob <- 0.2
cateffect <- 0.1
###########################################################################





