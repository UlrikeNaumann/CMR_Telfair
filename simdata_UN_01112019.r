# HW, - adapted by UN - Estimating Abundance from multiple sampling capture recatpure data via multi-state
# simulating capture recapture data

require(pracma)
# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data for 2 genders
# ------------------------------------------------------------------------------
# Name: sim.data.gender
# Objective: To generate multiple years of multi-state CMR data using the given parameter values
# Inputs: Nmale - total number of male  
#         Nfemale - total number of female    
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified, assume constant for genders 
#           - can be either one value or a vector with the first value for male, second for female 
#         smale - male survival, constant -  next: dependent on t
#         sfemale - female survival, constant - next: dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         pmale - capture probability for each state, constant, length R, for male . next: dependent on t 
#         pfemale - capture probability for each state, constant, length R, for female . next: dependent on t 
#         # psi - transition probabilities matrix - seems we don't need it, maybe use again for model with temp emigration
#         pclassmal - classification probabilities for true male. vector of 3 values, classified as male/female/unknown given male
#         pclassfem - classification probabilities for true female. vector of 3 values, classified as male/female/unknown given female
#         pclassmal needs to sum to 1, pclassfem needs to sum to 1, last class unknown can be missing, then = 0
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history matrix for each year, truth of which state individual was in whilst present at site
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = male, present, 2 = female, present, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders

sim.data.gender <- function(Nmale, Nfemale ,T, K=1, r, smale,sfemale , 
                            pmale, pfemale, pclassmal, pclassfem, random = FALSE, psi = matrix(c(1,0,0,1), ncol=2,nrow=2))  {

    if(sum(pclassmal) != 1 ) stop("classification probabilities must sum up to 1!")
    if(sum(pclassfem) != 1 ) stop("classification probabilities must sum up to 1!")
  
    if(length(pclassmal) == 2) pclassmal <- c(pclassmal, 0)
    if(length(pclassfem) == 2) pclassfem <- c(pclassfem, 0)
  
    if(smale < 0 | sfemale <0 | pmale < 0 | pfemale < 0 |smale >1 | sfemale >1 | pmale >1 | pfemale >1  ) stop(
      "probabilities must have values between 0 and 1"  )
    
    # now run the simulation for male
    # possibly change how N is coded in sim.data! to avoid having to provide 2 values
    testsimmale = sim.data(N=c(Nmale,0), T=T, K=K, r=rep(0.2, times=5), s=smale, p = pmale)
    
    # adjust the state for the animals where state was falsely coded
    # we need testsimmale$state.laake, testsimmale$encounter,  pclassmal
    # state.obsm: 0 = not observed. 1 = observed as male, 2 = observed as female, 3 = observed as unknown gender
    
    # get multinomial draws
    tmp <- rmultinom(n=sum(testsimmale$attendance), size=1, prob=pclassmal)
    # make these into a vector of 1,2, or 3
    obsgender <- apply(tmp , MARGIN=2, FUN = function(x) which(x == 1))
    state.obsm <- matrix(0, nrow = Nmale, ncol = T)
    state.encm <- matrix(0, nrow = Nmale, ncol = 3)
    noencounters <- rowSums(testsimmale$encounter)
    running <- 1
    for(i in 1:length(noencounters)){
      
      if(noencounters[i]!= 0){ # if 0, go to next animal
        # get no. of encounters per animal
        wh <- which(testsimmale$encounter == 1)
        
        state.obsm[i,wh] <- obsgender[running:(running+noencounters[i] - 1)]
        # how many of state.obs[i,wh] are respectively 1, 2 or 3
        state.encm[i,] <- c(sum(state.obsm[i,wh]==1) , sum(state.obsm[i,wh]==2), sum(state.obsm[i,wh]==3)  ) 
        # code state as 3 numbers for each indiv.: of all obs. times obs as male/female/unknown
        running <-   running+noencounters[i] 
      }
    }
    
     
    testsimfemale = sim.data(N=c(0, Nfemale), T=T, K=K, r=rep(0.2, times=5), s=sfemale, p = pfemale)
    
    # adjust the state for the animals where state was falsely coded
    # we need testsimmale$state.laake, testsimmale$encounter,  pclassmal
    # state.obsf: 0 = not observed. 1 = observed as male, 2 = observed as female, 3 = observed as unknown gender
    
    # get multinomial draws
    tmp <- rmultinom(n=sum(testsimfemale$attendance), size=1, prob=pclassmal)
    # make these into a vector of 1,2, or 3
    obsgender <- apply(tmp , MARGIN=2, FUN = function(x) which(x == 1))
    state.obsf <- matrix(0, nrow = Nmale, ncol = T)
    state.encf <- matrix(0, nrow = Nmale, ncol = 3)
    noencounters <- rowSums(testsimfemale$encounter)
    running <- 1
    for(i in 1:length(noencounters)){
      
      if(noencounters[i]!= 0){ # if 0, go to next animal
        # get no. of encounters per animal
        wh <- which(testsimfemale$encounter == 1)
        
        state.obsf[i,wh] <- obsgender[running:(running+noencounters[i] - 1)]
        # how many of state.obs[i,wh] are respectively 1, 2 or 3
        state.encf[i,] <- c(sum(state.obsf[i,wh]==1) , sum(state.obsf[i,wh]==2), sum(state.obsf[i,wh]==3)  ) 
        # code state as 3 numbers for each indiv.: of all obs. times obs as male/female/unknown
        running <-   running+noencounters[i] 
      }
    }
    
    # TO DO
    # now combine the results from testsimfemale and testsimmale
    
    
    # observed gender proportion at each time point. HERE, we ignore whether the same animal is assessed with different gender
    # at different time point. we have male, female, uncertain.
    p.state.obs <- matrix(0, ncol = T, row = 3)
    p.state.obsN <- matrix(0, ncol = T, row = 3)
    
    state.enc <- rbind( state.encm , state.encf ) 
    state.obs <- rbind(state.obsm, state.obsf)
    for(t in 1:T){
      p.state.obsN[,i] <- c( sum( state.obs[,i] == 1), sum( state.obs[,i] == 2),sum( state.obs[,i] == 3)  )
      p.state.obs[,i]  <- p.state.obsN[,i] / sum( p.state.obsN[,i])
    }
    
      # what is the estimated gender proportion after readjusting the gender for animals with uncertain gender assessment
    #  p.state.obs.corrected # this term will be needed for decoding algorithm! 
    #TO DO
    # eg for each animal where gender is unclear, we could use an extra category: unclear
      
    if(random){
      
      # randomly mix together male and female, ie. mix up all observations
      neworder <- order( rnorm(sum(Nmale, Nfemale)) )
      
      # Nmale, Nfemale
       # frand = rbinom(sum(Nmale, Nfemale),1,.5) # works only if the same number of male and female
       # does not work, the number needs to be exact
       
       # neworder[frand == 0] <- 1:Nmale
       # neworder[frand == 1] <-(Nmale+1):(Nmale+Nfemale) 
       # use neworder to reorder all outputs 
      
      # now reorder all outputs
      attendance <- attendance[neworder ,]
      encounter  <- encounter[neworder ,]
      state      <- state[neworder]
      state.laake <- state.laake[neworder ,];
      # p.state # true gender proportion
      for(j in 1:T){
      capture.full[[j]] <- capture.full[[j]][neworder]
      capture.0 <- which(capture.full[[j]] == 0)
        capture[[j]] <- capture.full[[j]][- capture.0] 
      }    
    }
    # N.obs = number of individuals seen each year (observed gender) 
    N.t <- matrix(0, ncol = T, nrow = 3)
    N.t[2,] <- testsimmale$N.attend
    N.t[3,] <- testsimfemale$N.attend
    N.t[1,] <- N.t[3,] + N.t[2,];
      #  N.c = N.encounter - number of individuals seen each year (true gender) 
      
    N.c <- matrix(0, ncol = T, nrow = 3)
    N.c[2,] <- testsimmale$N.encounter
    N.c[3,] <- testsimfemale$N.encounter
    N.c[1,] <- N.c[3,] + N.c[2,];
    
    N.obs <- 
    nsee <- c(testsimmale$N.seen , testsimfemale$N.seen) 
    # state.laake: 0 not yet recruited. 1 recruited male. 2 recruited female. 3 departed
    
    # to do: still need to refine this
    return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, 'p.state' = p.state,
                'p.state.obs' = p.state.obs, 'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 
                'N.encounter'=N.c, 'N.obs' = N.obs , 'N.seen'= c(sum(nsee),nsee ) ))
}
#head(testsim1$state.laake)
 #head(testsim1$encounter)
# dim(testsim1$capture); head(testsim1$capture)
# testsim1$N.encounter

# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data
# ------------------------------------------------------------------------------
# Name: sim.data
# Objective: To generate multiple years of multi-state stopover data using the given parameter values
# Inputs: N - total number of individuals # to do: 2 numbers - male and female
#             if only one number specified, assume 1/2 male, half female
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified, assume constant for genders - to do:  make dependent on gender
#         s - survival, constant - to do: make dependent on gender. next: dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         p - capture probability for each state, constant, length R -to do: make dependent on gender. next: dependent on t 
#         # psi - transition probabilities matrix - seems we don't need it, maybe use again for model with temp emigration
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history vector for each year, truth of which state individual was in whilst present at site (male/female)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = male, present, 2 = female, present, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders

sim.data <- function(N, T, K=1, r, s, p, psi = matrix(c(1,0,0,1), ncol=2,nrow=2))  {

  if (length(N)==1) Ngender <- c(floor(N/2) , ceiling(N/2)) else if(length(N)==2) Ngender <- N else {
    stop("Only 2 different states (gender). Please specifify N correctly.") }
  # if N is one number but not even, the first state will have 1 less individuum then the second state
  
  # reassign N
  N <- Ngender[1] + Ngender[2] ;
  
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # next: array with gender as dimension
  encounter <- matrix(1, nrow=N, ncol=T)  # next: array with gender as dimension
  p.state <- rep(0, times = T)
  capture.full <- NULL
  capture <- NULL
  
  if(length(r)==1) r<- rep(r, times = T) else if(length(r) != T) stop("r needs to have length T");
  
  # attendance history
  # generate recruitment distribution

  
  # recruitment <- rmultinom(1, N, r) # to do. with N and r  as 2dim vector (gender) 
  recruitment1 <- rmultinom(1, Ngender[1], r)
  recruitment2 <- rmultinom(1, Ngender[2], r)
  # make vector of arrival times
  # recruitment.i <- rep(1:T, recruitment) # adapt for 2-gender version
  recruitment1.i <- rep(1:T, recruitment1)
  recruitment2.i <- rep(1:T, recruitment2) 
  recruitment.i <- c(recruitment1.i, recruitment2.i);
  # state <-   # true state (gender), vector of length N
  state <- c(rep(1, times = Ngender[1]), rep(2, times = Ngender[2]) ); 
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]] <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- state[i]
    # state.laake: 0 not yet recruited. 1 recruited male. 2 recruited female. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year
          survival <- rbinom(1,1,s) # adapt with survival prob s dependent on gender
          # store 0 or 1 following year
          attendance[i,t+1] <- survival
          if (survival == 0) state.laake[i,t+1] <- 3 else state.laake[i,t+1] <- state[i];
        } else state.laake[i,t+1] <- 3 # ie departed
      }
    }  
  }
  
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  N1.t <- colSums(attendance[state == 1 , ]);
  N2.t <- colSums(attendance[state == 2 , ]);
  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration

  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    
    capture.t <- matrix(0, nrow=N, ncol=K)

    # I will need this when the states are also temporary emigration states - then adapt this
    # loop over individuals and assign states
#    for (i in 1:N) {
      # find the occasions when present
#      present <- which(presence.t[i,] == 1)
      # assign state at initial capture
#      state.init <- 1+rbinom(1,1,alpha[2])  # if necessary, add alpha again as argument, if needed for temp. emigration modelling
#      state.t[i,present[1]] <- state.init
      # loop over remaining occasions when present and transition between states
#      if (length(present) > 1)  {
#        for (k in present[2]:present[length(present)])  {
#          if (state.t[i,k-1] == 1)  {
#            state.k <- 1+rbinom(1,1,psi[1,2])  # if necessary, add psi again as argument, if needed for temp. emigration modelling
#          } else if (state.t[i,k-1] == 2)  {
#            state.k <- 1+rbinom(1,1,psi[2,2])
#          }
#          state.t[i,k] <- state.k
#        }
#      }
#    }
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
         # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
      # loop over occasions
      for (k in 1:K)  {
        # determine whether they were captured or not
        # capture.t[i,k] <- state.t[i,1]*rbinom(1,1,p[state.t[i,1]])
        capture.t[i,k] <- state[i]*rbinom(1,1,p[state[i] ])  # needs to be adapted if temp emigration possible
      } }
    }
    
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0)
    
    # update encounter histories
    for (c in capture.0)  {
      encounter[c,t] <- 0
    }
    
    # store state history for year t
   #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <- c( capture.t)  #  capture.t
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}

# example
# testsim1 = sim.data(N=1000, T=5, K=1, r=rep(0.2, times=5), s=0.3, p = c(0.3,0.5))
# testsim = sim.data(N=c(1000,0), T=5, K=1, r=rep(0.2, times=5), s=0.3, p = 0.4) # only male
# testsim1$p.state; head(testsim1$encounter)
# head(testsim1$state); tail(testsim1$state)
# head(testsim1$state.laake); tail(testsim1$state.laake)
# head(testsim1$attendance)
# dim(testsim$attendance); 
# head(testsim$state.laake); tail(testsim$state.laake)


# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data WITHOUT gender
# ------------------------------------------------------------------------------
# Name: sim.data.0
# Objective: To generate multiple years of capture recapture data using the given parameter values, constant survival and capture prob
# Inputs: N - total number of individuals 
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified
#         s - survival, constant - . next: dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         p - capture probability for each state, constant, length R -to do: make dependent on gender. next: dependent on t 
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history vector for each year, truth of which state individual was in whilst present at site (male/female)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = present,  3 = dead  # not yet: 1 = male, present, 2 = female, present,
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders

sim.data.0 <- function(N, T, K=1, r, s, p, psi = matrix(c(1,0,0,1), ncol=2,nrow=2))  {
  
  #  if (length(N)==1) Ngender <- c(floor(N/2) , ceiling(N/2)) else if(length(N)==2) Ngender <- N else {
  #    stop("Only 2 different states (gender). Please specifify N correctly.") }
  # if N is one number but not even, the first state will have 1 less individuum then the second state
  
  # reassign N
  # N <- Ngender[1] + Ngender[2] ;
  
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # next: array with gender as dimension
  encounter <- matrix(1, nrow=N, ncol=T)  # next: array with gender as dimension
  # p.state <- rep(0, times = T)
  capture.full <- NULL
  capture <- NULL
  
  if(length(r)==1) r<- rep(r, times = T) else if(length(r) != T) stop("r needs to have length T");
  
  # attendance history
  # generate recruitment distribution
  
  
  recruitment <- rmultinom(1, N, r) # produces a matrix wiht rows = T = lenght r , with the number of elements first recruited in each row
  # to do. with N and r  as 2dim vector (gender) 
  # recruitment1 <- rmultinom(1, Ngender[1], r)
  # recruitment2 <- rmultinom(1, Ngender[2], r)
  # make vector of arrival times
  recruitment.i <- rep(1:T, recruitment) # a number for each individual # adapt for 2-gender version
  # recruitment1.i <- rep(1:T, recruitment1)
  # recruitment2.i <- rep(1:T, recruitment2) 
  # recruitment.i <- c(recruitment1.i, recruitment2.i);
  # state <-   # true state (gender), vector of length N
  # state <- c(rep(1, times = Ngender[1]), rep(2, times = Ngender[2]) ); 
  #   state <- rep(1, times = N) 
  # now use
  state <- matrix(0, nrow=N, ncol=T) # state = 0, and 1 if alive
  
  # instead of state, I now use state.laake. state was in HW a list with attendance vector as elements for each secondary period
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]] <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- 1 # state[i]
    state[i,recruitment.i[i]] <- 1 
    # state.laake: 0 not yet recruited. 1 recruited male. 2 recruited female. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year
          survival <- rbinom(1,1,s) # adapt with survival prob s dependent on gender
          # store 0 or 1 following year
          attendance[i,t+1] <- survival
          state[i,t+1] <- survival;
          if (survival == 0) { state.laake[i,t+1] <- 3 
          # state[i,t+1] <- 0 # it is defaulted to 0 anyway
          } else { state.laake[i,t+1] <- 1; 
          #  state[i,t+1] <- 1;
          }#  state[i];
          
        } else { # ie, between first recruitment and end of study are some time points, but has died now
          state.laake[i,t+1] <- 3 # ie departed - 
          state[i,t+1] <- 0;
          # since this is developed from gender version, we have states: 0 (unborn)), 1 (alive), 3 (departed)
        }
      }
    }  
  }
  
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  # N1.t <- colSums(attendance[state == 1 , ]);
  # N2.t <- colSums(attendance[state == 2 , ]);
  #  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration
  
  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    
    capture.t <- matrix(0, nrow=N, ncol=K)
    
    # I will need this when the states are also temporary emigration states - then adapt this
    # loop over individuals and assign states
    #    for (i in 1:N) {
    # find the occasions when present
    #      present <- which(presence.t[i,] == 1)
    # assign state at initial capture
    #      state.init <- 1+rbinom(1,1,alpha[2])  # if necessary, add alpha again as argument, if needed for temp. emigration modelling
    #      state.t[i,present[1]] <- state.init
    # loop over remaining occasions when present and transition between states
    #      if (length(present) > 1)  {
    #        for (k in present[2]:present[length(present)])  {
    #          if (state.t[i,k-1] == 1)  {
    #            state.k <- 1+rbinom(1,1,psi[1,2])  # if necessary, add psi again as argument, if needed for temp. emigration modelling
    #          } else if (state.t[i,k-1] == 2)  {
    #            state.k <- 1+rbinom(1,1,psi[2,2])
    #          }
    #          state.t[i,k] <- state.k
    #        }
    #      }
    #    }
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
      # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
        # loop over occasions
        capture.t[i,] <- rbinom(K,1,p)  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob! in primary time  - possible would also be time dept within secondary
        #    for (k in 1:K)  {
        # determine whether they were captured or not
        ## capture.t[i,k] <- state.t[i,1]*rbinom(1,1,p[state.t[i,1]])
        ## capture.t[i,k] <- state[i]*rbinom(1,1,p[state[i] ])  # needs to be adapted if temp emigration possible
        #  capture.t[i,k] <- rbinom(1,1,p)  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob!
        # } 
      }
    }
    
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0) # ie for which year are all secondary capt histories = 0
    
    # update encounter histories
    for (c in capture.0)  {
      encounter[c,t] <- 0
    }
    
    # store state history for year t
    #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <- c( capture.t)  #  capture.t
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, # 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}

# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data WITH temp. emigration
# ------------------------------------------------------------------------------
# this is now out of date, and replaced by sim.data.tau.time
# Name: sim.data.tau
# Objective: To generate multiple years of capture recapture data using the given parameter values, constant survival and capture prob
# Inputs: N - total number of individuals 
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified
#         s - survival, constant - . next: dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         p - capture probability for each state, constant, length R -to do: make dependent on gender. next: dependent on t 
#         tau - probability of temporary emigration, given 1) not temp. emigrated in previous year, 
#                                                          2) temp. emigrated in previous year, 
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history vector for each year, truth of which state individual was in whilst present at site (male/female)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = alive, 2 = temp emig, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders
# temp = sim.data.tau(N=100, T=10, K=1, r= 0.6, s=0.5, p=0.5, tau = c(0.2, 0.3)) ; temp$capture.full; temp$N.attend ; temp$N.encounter; temp$N.seen
# temp$state; temp$state.laake
sim.data.tau <- function(N, T, K=1, r, s, p, tau)  {
  
#  if (length(N)==1) Ngender <- c(floor(N/2) , ceiling(N/2)) else if(length(N)==2) Ngender <- N else {
#    stop("Only 2 different states (gender). Please specifify N correctly.") }
  # if N is one number but not even, the first state will have 1 less individuum then the second state
  
  # reassign N
  # N <- Ngender[1] + Ngender[2] ;
  if(s + tau[1] > 1) {
    warning("s + tau[1] must be <= 1")
    sp = s + tau[1];
    s <- s/sp;
    tau[1]/sp 
  }  
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # alive an not temp emigrated
  encounter <- matrix(1, nrow=N, ncol=T)  # next: array with gender as dimension
 # p.state <- rep(0, times = T)
  capture.full <- NULL
  capture <- NULL
  
  if(length(r)==1) r<- rep(r, times = T) else if(length(r) != T) stop("r needs to have length T");
  
  # attendance history
  # generate recruitment distribution
  
  
  recruitment    <- rmultinom(1, N, r) # produces a matrix with rows = T = lenght r , with the number of elements first recruited in each row

  # make vector of arrival times
   recruitment.i <- rep(1:T, recruitment) # a number for each individual # adapt for 2-gender version
   # now use
   state <- matrix(0, nrow=N, ncol=T) # state = 0, and 1 if alive
   
   # instead of state, I now use state.laake. state was in HW a list with attendance vector as elements for each secondary period
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
  # state = 1 means temp. emigrated or alive
  # attendance is 1 if alive and not temp emigrated!
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]]  <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- 1 
    state[i,recruitment.i[i]] <- 1 # state should be used later to for male/female, here: just one gender
    # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year, or temp emigrate or die
          survival <- rmultinom(1,1,c(s,tau[1], 1-s-tau[1] ) ) # adapt with survival prob s dependent on gender
          surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
          # store 0 or 1 following year
          
          # state[i,t+1] <- survival;
          state.laake[i,t+1] <- surv
          
          if (surv == 3) {  # state.laake[i,t+1] <- 3 
                       state[i,t+1] <- 0 # it is defaulted to 0 anyway
          } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
                  state[i,t+1] <- 1;
                  attendance[i,t+1] <- 1
          } else if(surv == 2){ # temp. emigrated
            #  attendance[i,t+1] <- 0 # is anyway the default
            state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
            # here: state = 1 means temp. emigrated or alive
          }
        # attendance is 1 if alive and not temp emigrated!
          } else { # is either dead or temp emigrated
            if (state.laake[i,t] == 2){ # temp. emigrated
              # now probability of emigrationg again is tau[2]
              emigagain <- rbinom(1,1, tau[2]) # 1 = emig 0 = alive
              state[i,t+1] <- 1 # either way from temp. emig. no death
              state.laake[i,t+1] <- emigagain + 1  #  1 = alive, 2 = temp. emig, 3 = dead
              attendance[i,t+1] == 1 - emigagain
            } else {  # ie dead -
            
            state.laake[i,t+1] <- 3  
            state[i,t+1] <- 0;
        # since this is developed from gender version, we have states: 0 (unborn)), 1 (alive), 3 (departed)
          }
      }
    }  
  }
  }
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  # N1.t <- colSums(attendance[state == 1 , ]);
  # N2.t <- colSums(attendance[state == 2 , ]);
 #  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration
  
  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    
    capture.t <- matrix(0, nrow=N, ncol=K)
    
    # I will need this when the states are also temporary emigration states - then adapt this
    # loop over individuals and assign states
    #    for (i in 1:N) {
    # find the occasions when present
    #      present <- which(presence.t[i,] == 1)
    # assign state at initial capture
    #      state.init <- 1+rbinom(1,1,alpha[2])  # if necessary, add alpha again as argument, if needed for temp. emigration modelling
    #      state.t[i,present[1]] <- state.init
    # loop over remaining occasions when present and transition between states
    #      if (length(present) > 1)  {
    #        for (k in present[2]:present[length(present)])  {
    #          if (state.t[i,k-1] == 1)  {
    #            state.k <- 1+rbinom(1,1,psi[1,2])  # if necessary, add psi again as argument, if needed for temp. emigration modelling
    #          } else if (state.t[i,k-1] == 2)  {
    #            state.k <- 1+rbinom(1,1,psi[2,2])
    #          }
    #          state.t[i,k] <- state.k
    #        }
    #      }
    #    }
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
      # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
        # loop over occasions
        capture.t[i,] <- rbinom(K,1,p)  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob! in primary time  - possible would also be time dept within secondary
    #    for (k in 1:K)  {
          # determine whether they were captured or not
          ## capture.t[i,k] <- state.t[i,1]*rbinom(1,1,p[state.t[i,1]])
         ## capture.t[i,k] <- state[i]*rbinom(1,1,p[state[i] ])  # needs to be adapted if temp emigration possible
        #  capture.t[i,k] <- rbinom(1,1,p)  # needs to be adapted if temp emigration possible
          # next: use time dependent capture prob!
       # } 
        }
    }
    
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0) # ie for which year are all secondary capt histories = 0
    
    # update encounter histories
    for (c in capture.0)  {
      encounter[c,t] <- 0
    }
    
    # store state history for year t
    #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <- c( capture.t)  #  capture.t
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, # 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}


# to do: should we allow the first time point of attendance, that the animal is already temp. emigrated?
# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data WITH temp. emigration
# ------------------------------------------------------------------------------
# Name: sim.data.tau.t
# Objective: To generate multiple years of capture recapture data using the given parameter values, constant survival and capture prob
# Inputs: N - total number of individuals 
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified
#         s - survival, dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         p - capture probability for each state, dependent on t 
#         tau - probability of temporary emigration, given 1) not temp. emigrated in previous year, 
#                                                          2) temp. emigrated in previous year, 
#         this is now changed to be a factor which is multiplicated (1-tau)*s to the survival
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history matrix (vector for each year), truth of which state individual was in whilst present at site (0 = not present, 1 = present)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = alive, 2 = temp emig, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders
# temp = sim.data.tau.time(N=100, T=10, K=1, r= 0.6, s=rep(0.5,10)-seq(0.01, 0.1, length=10), p=rep(0.5,10)+seq(0.01, 0.1, length=10), tau = c(0.2, 0.3)) ; temp$capture.full; temp$N.attend ; temp$N.encounter; temp$N.seen
# temp$state; temp$state.laake
# 1/11.2019: this is now recoded. survprob = s*(1-tau), emig prob = s*tau, dead prob = 1-s

sim.data.tau.time <- function(N, T, K=1, r, s, p, tau)  {

  if(length(r)==1) r<- rep(r, times = T) else if(length(r) != T) stop("r needs to have length T");
  if(length(s)==1) s<- rep(s, times = T) else if(length(s) != T) stop("s needs to have length T");
  if(length(p)==1) p<- rep(p, times = T) else if(length(p) != T) stop("p needs to have length T");
    
  if(any(s  > 1)) {
    stop("s must be <= 1")  }
  if(any(p  > 1)) {
    stop("p must be <= 1")  }
  if(any(r  > 1)) {
    stop("r must be <= 1")  }
  if(any(tau  > 1)) {
    stop("tau must be <= 1")  }
    # sp = s + tau[1];
    # s <- s/sp;
    # tau[1]<- tau[1]/sp 
   
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # alive an not temp emigrated
  encounter <- matrix(1, nrow=N, ncol=T)  # next: array with gender as dimension
  # p.state <- rep(0, times = T)
  capture.full <- NULL
  capture <- NULL
  
  
  # attendance history
  # generate recruitment distribution
  
  
  recruitment    <- rmultinom(1, N, r) # produces a matrix with rows = T = lenght r , with the number of elements first recruited in each row
  
  # make vector of arrival times
  recruitment.i <- rep(1:T, recruitment) # a number for each individual # adapt for 2-gender version
  # now use
  state <- matrix(0, nrow=N, ncol=T) # state = 0, and 1 if alive
  
  # instead of state, I now use state.laake. state was in HW a list with attendance vector as elements for each secondary period
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
  # state = 1 means temp. emigrated or alive
  # attendance is 1 if alive and not temp emigrated!
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]]  <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- 1 
    state[i,recruitment.i[i]] <- 1 # state should be used later to for male/female, here: just one gender
    # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year, or temp emigrate or die
          # survival <- rmultinom(1,1,c(s[t],tau[1], 1-s[t]-tau[1] ) ) 
          # recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
          survival <- rmultinom(1,1,c(s[t]*(1-tau[1]),s[t]*tau[1], 1-s[t] ) ) 
          surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
          # store 0 or 1 following year
          
          # state[i,t+1] <- survival;
          state.laake[i,t+1] <- surv
          
          if (surv == 3) {  # state.laake[i,t+1] <- 3 
            state[i,t+1] <- 0 # it is defaulted to 0 anyway
          } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
            state[i,t+1] <- 1;
            attendance[i,t+1] <- 1
          } else if(surv == 2){ # temp. emigrated
            #  attendance[i,t+1] <- 0 # is anyway the default
            state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
            # here: state = 1 means temp. emigrated or alive
          }
          # attendance is 1 if alive and not temp emigrated!
        } else { # is either dead or temp emigrated
          if (state.laake[i,t] == 2){ # temp. emigrated
            # lets build in a death prob.: recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
            survival <- rmultinom(1,1,c(s[t]*(1-tau[2]),s[t]*tau[2], 1-s[t] ) ) 
            surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
            # store 0 or 1 following year
            
            # state[i,t+1] <- survival;
            state.laake[i,t+1] <- surv
            
            if (surv == 3) {  # state.laake[i,t+1] <- 3 
              state[i,t+1] <- 0 # it is defaulted to 0 anyway
            } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
              state[i,t+1] <- 1;
              attendance[i,t+1] <- 1
            } else if(surv == 2){ # temp. emigrated
              #  attendance[i,t+1] <- 0 # is anyway the default
              state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
              # here: state = 1 means temp. emigrated or alive
            }
            
            # now probability of emigrationg again is tau[2]
            # emigagain <- rbinom(1,1, tau[2]) # 1 = emig 0 = alive
            # state[i,t+1] <- 1 # either way from temp. emig. no death
            # assumption here is: no death after temp emig? 
            # state.laake[i,t+1] <- emigagain + 1  #  1 = alive, 2 = temp. emig, 3 = dead
            # attendance[i,t+1] == 1 - emigagain
          } else {  # ie dead -
            
            state.laake[i,t+1] <- 3  
            state[i,t+1] <- 0;
            # since this is developed from gender version, we have states: 0 (unborn)), 1 (alive), 3 (departed)
          }
        }
      }  
    }
  }
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  # N1.t <- colSums(attendance[state == 1 , ]);
  # N2.t <- colSums(attendance[state == 2 , ]);
  #  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration
  
  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    
    capture.t <- matrix(0, nrow=N, ncol=K)
    
    # I will need this when the states are also temporary emigration states - then adapt this
    # loop over individuals and assign states
    #    for (i in 1:N) {
    # find the occasions when present
    #      present <- which(presence.t[i,] == 1)
    # assign state at initial capture
    #      state.init <- 1+rbinom(1,1,alpha[2])  # if necessary, add alpha again as argument, if needed for temp. emigration modelling
    #      state.t[i,present[1]] <- state.init
    # loop over remaining occasions when present and transition between states
    #      if (length(present) > 1)  {
    #        for (k in present[2]:present[length(present)])  {
    #          if (state.t[i,k-1] == 1)  {
    #            state.k <- 1+rbinom(1,1,psi[1,2])  # if necessary, add psi again as argument, if needed for temp. emigration modelling
    #          } else if (state.t[i,k-1] == 2)  {
    #            state.k <- 1+rbinom(1,1,psi[2,2])
    #          }
    #          state.t[i,k] <- state.k
    #        }
    #      }
    #    }
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
      # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
        # loop over occasions
        capture.t[i,] <- rbinom(K,1,p[t])  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob! in primary time  - possible would also be time dept within secondary
        #    for (k in 1:K)  {
        # determine whether they were captured or not
        ## capture.t[i,k] <- state.t[i,1]*rbinom(1,1,p[state.t[i,1]])
        ## capture.t[i,k] <- state[i]*rbinom(1,1,p[state[i] ])  # needs to be adapted if temp emigration possible
        #  capture.t[i,k] <- rbinom(1,1,p)  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob!
        # } 
      }
    }
    
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0) # ie for which year are all secondary capt histories = 0
    
    # update encounter histories
    for (c in capture.0)  {
      encounter[c,t] <- 0
    }
    
    # store state history for year t
    #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <- c( capture.t)  #  capture.t
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, # 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}


# to do: should we allow the first time point of attendance, that the animal is already temp. emigrated?

# ------------------------------------------------------------------------------
# Code to generate robust multi-state (multi-period) CRM data WITH temp. emigration
# ------------------------------------------------------------------------------
# Name: sim.data.robust.tau.time
# Objective: To generate multiple years of capture recapture data using the given parameter values, constant survival and capture prob
# Inputs: N - total number of individuals 
#         T - number of years
#         K - number of occasions each year, allow to vary over time, can be length = 1 or T
#         r - recruitment probabilities, length T, if not state dept specified
#         s - survival, dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         p.sec - capture probability, dependent on t, assumed to be constant within primary
#             here: p.sec = capt prob for each sec. point k[j] in t, prob of being obs at all in the primary time t
#                       is then p.prim = 1 - dbinom(x=0, size=K, prob=p.sec)
#         p.prim - probability of being captured at least once in the primary occassion t, time dept
#         from p.prim, can derive p.sec = 1 - nthroot(x=1 - p.prim,n=K)
#           - next: could allow for additional option. eg trap shyness, or linear increase/decrease within primary
#           - however, these would be only to test the general program for robustness
#         for specification of p.sec or p.prim: p.sec overwrites p.prim, if both are specified!
#               if K differs for each t, the same p.sec for each prim time point, does not correspond to the same p.prim
#         tau - probability of temporary emigration, given 1) not temp. emigrated in previous year, 
#                                                          2) temp. emigrated in previous year, 
#         this is now changed to be a factor which is multiplicated (1-tau)*s to the survival
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, number of times an individuals was seen in each year
#          state - state history vector for each year, truth of which state individual was in whilst present at site (male/female)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = alive, 2 = temp emig, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a matrix for each year for all individuals, no of times in t that individuals were seen
#          ? capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders
# temp = sim.data.robust.tau.time(N=100, T=10, K=5, r= 0.6, s=rep(0.5,10)-seq(0.01, 0.1, length=10), p.prim=rep(0.5,10)+seq(0.01, 0.1, length=10), tau = c(0.2, 0.3)) ; 
# temp = sim.data.robust.tau.time(N=100, T=10, K=c(3,3,3,3,2,2,2,2,2,2), r= 0.6, s=rep(0.5,10)-seq(0.01, 0.1, length=10), p.prim=rep(0.5,10)+seq(0.01, 0.1, length=10), tau = c(0.2, 0.3)) ; 
# temp$capture.full; temp$N.attend ; temp$N.encounter; temp$N.seen
# temp$state; temp$state.laake; temp$encounter; temp$capture
# 1/11.2019: this is now recoded. survprob = s*(1-tau), emig prob = s*tau, dead prob = 1-s

sim.data.robust.tau.time <- function(N, T, K=1, r, s, p.prim=NULL, p.sec=NULL, tau)  {
  

  if(length(r)==1) r<- rep(r, times = T) else if(length(r) != T) stop("r needs to have length T");
  if(length(s)==1) s<- rep(s, times = T) else if(length(s) != T) stop("s needs to have length T");
#   if(length(p)==1) p<- rep(p, times = T) else if(length(p) != T) stop("p needs to have length T");
  if(length(K)==1) K<- rep(K, times = T) else if(length(K) != T) stop("K needs to have length T");
  
  if(is.null(p.prim)) { 
    if(is.null(p.sec)) stop("either p.prim or p.sec need to be specified") else {
      p.prim = 1 - dbinom(x=0, size=K, prob=p.sec) # returns vector if K is a vector
    }
  }
  if(is.null(p.sec)) p.sec = 1 -  nthroot(1-p.prim,K);
  # from here on, we use only p = p.sec for calculations  
  p <- p.sec
  if(length(p)==1) p<- rep(p, times = T) else if(length(p) != T) stop("p.prim / p.sec needs to have length T");
  
  if(any(s  > 1)) {
    stop("s must be <= 1")  }
  if(any(p  > 1)) {
    stop("p must be <= 1")  }
  if(any(r  > 1)) {
    stop("r must be <= 1")  }
  if(any(tau  > 1)) {
    stop("tau must be <= 1")  }
  # sp = s + tau[1];
  # s <- s/sp;
  # tau[1]<- tau[1]/sp 
  
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # alive an not temp emigrated
  encounter <- matrix(0, nrow=N, ncol=T)  # encounter is a matrix with T columns
  # p.state <- rep(0, times = T)
  capture.full <- NULL                    # capture.full is a list of T matrizes with each K columns
  capture <- NULL
  
  
  # attendance history
  # generate recruitment distribution
  
  
  recruitment    <- rmultinom(1, N, r) # produces a matrix with rows = T = lenght r , with the number of elements first recruited in each row
  
  # make vector of arrival times
  recruitment.i <- rep(1:T, recruitment) # a number for each individual # adapt for 2-gender version
  # now use
  state <- matrix(0, nrow=N, ncol=T) # state = 0, and 1 if alive
  
  # instead of state, I now use state.laake. state was in HW a list with attendance vector as elements for each secondary period
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
  # state = 1 means temp. emigrated or alive
  # attendance is 1 if alive and not temp emigrated!
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]]  <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- 1 
    state[i,recruitment.i[i]] <- 1 # state should be used later to for male/female, here: just one gender
    # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year, or temp emigrate or die
          # survival <- rmultinom(1,1,c(s[t],tau[1], 1-s[t]-tau[1] ) ) 
          # recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
          survival <- rmultinom(1,1,c(s[t]*(1-tau[1]),s[t]*tau[1], 1-s[t] ) ) 
          surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
          # store 0 or 1 following year
          
          # state[i,t+1] <- survival;
          state.laake[i,t+1] <- surv
          
          if (surv == 3) {  # state.laake[i,t+1] <- 3 
            state[i,t+1] <- 0 # it is defaulted to 0 anyway
          } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
            state[i,t+1] <- 1;
            attendance[i,t+1] <- 1
          } else if(surv == 2){ # temp. emigrated
            #  attendance[i,t+1] <- 0 # is anyway the default
            state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
            # here: state = 1 means temp. emigrated or alive
          }
          # attendance is 1 if alive and not temp emigrated!
        } else { # is either dead or temp emigrated
          if (state.laake[i,t] == 2){ # temp. emigrated
            # lets build in a death prob.: recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
            survival <- rmultinom(1,1,c(s[t]*(1-tau[2]),s[t]*tau[2], 1-s[t] ) ) 
            surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
            # store 0 or 1 following year
            
            # state[i,t+1] <- survival;
            state.laake[i,t+1] <- surv
            
            if (surv == 3) {  # state.laake[i,t+1] <- 3 
              state[i,t+1] <- 0 # it is defaulted to 0 anyway
            } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
              state[i,t+1] <- 1;
              attendance[i,t+1] <- 1
            } else if(surv == 2){ # temp. emigrated
              #  attendance[i,t+1] <- 0 # is anyway the default
              state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
              # here: state = 1 means temp. emigrated or alive
            }
            
            # now probability of emigrationg again is tau[2]
            # emigagain <- rbinom(1,1, tau[2]) # 1 = emig 0 = alive
            # state[i,t+1] <- 1 # either way from temp. emig. no death
            # assumption here is: no death after temp emig? 
            # state.laake[i,t+1] <- emigagain + 1  #  1 = alive, 2 = temp. emig, 3 = dead
            # attendance[i,t+1] == 1 - emigagain
          } else {  # ie dead -
            
            state.laake[i,t+1] <- 3  
            state[i,t+1] <- 0;
            # since this is developed from gender version, we have states: 0 (unborn)), 1 (alive), 3 (departed)
          }
        }
      }  
    }
  }
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  # N1.t <- colSums(attendance[state == 1 , ]);
  # N2.t <- colSums(attendance[state == 2 , ]);
  #  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration
  
  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    capture.t <- matrix(0, nrow=N, ncol=K[t])  
    #capture.tt <- matrix(0, nrow=N, ncol=1) # i am only interested in no of captures
    # I will need this when the states are also temporary emigration states - then adapt this
    # loop over individuals and assign states
    #    for (i in 1:N) {
    # find the occasions when present
    #      present <- which(presence.t[i,] == 1)
    # assign state at initial capture
    #      state.init <- 1+rbinom(1,1,alpha[2])  # if necessary, add alpha again as argument, if needed for temp. emigration modelling
    #      state.t[i,present[1]] <- state.init
    # loop over remaining occasions when present and transition between states
    #      if (length(present) > 1)  {
    #        for (k in present[2]:present[length(present)])  {
    #          if (state.t[i,k-1] == 1)  {
    #            state.k <- 1+rbinom(1,1,psi[1,2])  # if necessary, add psi again as argument, if needed for temp. emigration modelling
    #          } else if (state.t[i,k-1] == 2)  {
    #            state.k <- 1+rbinom(1,1,psi[2,2])
    #          }
    #          state.t[i,k] <- state.k
    #        }
    #      }
    #    }
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
      # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
        # loop over occasions
        capture.t[i,] <- rbinom(K[t],1,p[t])  # outputs K values of 0 and 1
        
        # capture.tt[i,] <- sum(capture.ti[i,]) # rbinom(1,K,p[t])  # outputs the number between 0 and K
        # next: use time dependent capture prob! in primary time  - possible would also be time dept within secondary
        #    for (k in 1:K)  {
        # determine whether they were captured or not
        ## capture.t[i,k] <- state.t[i,1]*rbinom(1,1,p[state.t[i,1]])
        ## capture.t[i,k] <- state[i]*rbinom(1,1,p[state[i] ])  # needs to be adapted if temp emigration possible
        #  capture.t[i,k] <- rbinom(1,1,p)  # needs to be adapted if temp emigration possible
        # next: use time dependent capture prob!
        # } 
        encounter[i,t] <- sum(capture.t[i,]) # i am only interested in no of captures # outputs the number between 0 and K
          
      }
    }
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0) # ie for which year are all secondary capt histories = 0
    
    # update encounter histories
    for (c in capture.0)  {       encounter[c,t] <- 0    } # should work anyway

    # store state history for year t
    #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <-  capture.t # a matrix with N rows and K[t] columns
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter >= 1)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter  >= 1) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, # 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}



# ------------------------------------------------------------------------------
# Code to generate multi-state (multi-period) CRM data for 2 genders
# ------------------------------------------------------------------------------
# Name: sim.data.gender
# Objective: To generate multiple years of multi-state CMR data using the given parameter values
# Inputs: Nmale - total number of male  
#         Nfemale - total number of female    
#         T - number of years
#         K - number of occasions each year, constant
#         r - recruitment probabilities, length T, if not state dept specified, assume constant for genders 
#           - can be either one value or a vector with the first value for male, second for female 
#         smale - male survival, constant -  next: dependent on t
#         sfemale - female survival, constant - next: dependent on t
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         pmale - capture probability for each state, constant, length R, for male . next: dependent on t 
#         pfemale - capture probability for each state, constant, length R, for female . next: dependent on t 
#         # psi - transition probabilities matrix - seems we don't need it, maybe use again for model with temp emigration
#         pclassmal - classification probabilities for true male. vector of 3 values, classified as male/female/unknown given male
#         pclassfem - classification probabilities for true female. vector of 3 values, classified as male/female/unknown given female
#         pclassmal needs to sum to 1, pclassfem needs to sum to 1, last class unknown can be missing, then = 0
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, which years individuals were or were not seen
#          state - state history matrix for each year, truth of which state individual was in whilst present at site
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = male, present, 2 = female, present, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a vector for each year for all individuals, when individuals were or were not seen
#          capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders

# to do: should we allow the first time point of attendance, that the animal is already temp. emigrated?

# in preparation
# ------------------------------------------------------------------------------------------------
# Code to generate robust multi-state (multi-period) CRM data WITH temp. emigration for 2 genders
# -------------------------------------------------------------------------------------------------
# Name: sim.data.robust.tau.time.gender
# Objective: To generate multiple years of capture recapture data using the given parameter values, constant survival and capture prob
# Inputs: Nmale - total number of male  
#         Nfemale - total number of female    
#         T - number of years
#         K - number of occasions each year, constant, next step: allow not constant
#         rmale - recruitment probabilities for male, length T (or length 1 for constant)
#         rfemale - recruitment probabilities for male, length T (or length 1 for constant)
#         smale - male survival, length T
#         sfemale - male survival, length T
#         # alpha - initial discrete state probability (length R) - seems we don't need it, maybe use again for model with temp emigration
#         pmale.sec - capture probability for male, dependent on t, assumed to be constant within primary
#             here: p.sec = capt prob for each sec. point k[j] in t, prob of being obs at all in the primary time t
#                       is then p.prim = 1 - dbinom(x=0, size=K, prob=p.sec)
#         pfemale.sec - capture probability for female
#         pmale.prim - probability for male of being captured at least once in the primary occassion t, time dept
#         pfemale.prim - probability for female 
#         from p.prim, can derive p.sec = 1 - nthroot(x=1 - p.prim,n=K)
#           - next: could allow for additional option. eg trap shyness, or linear increase/decrease within primary
#           - however, these would be only to test the general program for robustness
#         taumale - probability of temporary emigration for male (length = 2), 
#                                   given 1) not temp. emigrated in previous year, 
#                                   2) temp. emigrated in previous year, 
#                   this is now changed to be a factor which is multiplicated (1-tau)*s to the survival
#         taufemale - probability of temporary emigration for female
#         pclassmal - classification probabilities for true male. vector of 3 values, classified as male/female/unknown given male
#         pclassfem - classification probabilities for true female. vector of 3 values, classified as male/female/unknown given female
#         pclassmal needs to sum to 1, pclassfem needs to sum to 1, last class unknown can be missing, then = 0
#         random - whether to randomly mix together male and female, ie. mix up all observations
# Outputs: attendance - attendance matrix, truth of which years individuals are present at the site
#          encounter - encounter matrix, number of times an individuals was seen in each year, dim = N x T
#          gencounter - encounter array, number of times an individuals was captured and classified as male [,,1], female [,,2], unknown [,,3] in each year, dim = N x T x 3
#          state - state history vector for each year, truth of which state individual was in whilst present at site (0 = not present 1 = male, 2 = female)
#          state.laake - state history matrix for each year, truth of which state individual was in whilst present at site,
#                        style as in Laake paper. 0 = not yet born, 1 = alive, 2 = temp emig, 3 dead
#          p.state - yearly true state (gender) proportion | attendance
#          capture.full - capture history list with a matrix for each year for all individuals, no of times in t that individuals were seen
#          ? capture - capture history matrix for each year, when individuals were or were not seen
#          N.attend - actual number of individuals at the site each year - to do: for all + both genders
#          N.encounter - number of individuals seen each year (+ gender specific)  - to do: for all + both genders
#          N.seen - number of individuals seen over entire study (+ gender specific)  - to do: for all + both genders
# temp = sim.data.robust.tau.time(N=100, T=10, K=5, r= 0.6, s=rep(0.5,10)-seq(0.01, 0.1, length=10), p.prim=rep(0.5,10)+seq(0.01, 0.1, length=10), tau = c(0.2, 0.3)) ; 
# temp$capture.full; temp$N.attend ; temp$N.encounter; temp$N.seen
# temp$state; temp$state.laake; temp$encounter; temp$capture
# 1/11.2019: this is now recoded. survprob = s*(1-tau), emig prob = s*tau, dead prob = 1-s

sim.data.gender <- function(Nmale, Nfemale ,T, K=1, rmale, rfemale, smale,sfemale ,pmale.prim=NULL, pfemale.prim=NULL, pmale.sec=NULL,pfemale.sec=NULL, 
                            taumale , taufemale, pmale, pfemale, pclassmal, pclassfem, random = FALSE)  {
  
sim.data.robust.tau.time.gender <- function(N, T, K=1, r, s, p.prim=NULL, p.sec=NULL, tau)  {
  
  if(sum(pclassmal) != 1 ) stop("classification probabilities must sum up to 1!")
  if(sum(pclassfem) != 1 ) stop("classification probabilities must sum up to 1!")
  
  if(length(pclassmal) == 2) pclassmal <- c(pclassmal, 0)
  if(length(pclassfem) == 2) pclassfem <- c(pclassfem, 0)
  
  if(any(smale < 0) | any(sfemale <0 ) | any(pmale < 0) | any(pfemale < 0) | any(smale >1) | any(sfemale >1) | any(pmale >1) | any(pfemale >1)  ) stop(
    "probabilities must have values between 0 and 1"  )
  
  if(is.null(pfemale.prim)) { 
    if(is.null(pfemale.sec)) stop("either pfemale.prim or p.sec need to be specified") else {
      pfemale.prim = 1 - dbinom(x=0, size=K, prob=pfemale.sec)
    }
  }
  if(is.null(pfemale.sec)) pfemale.sec = 1 -  nthroot(1-pfemale.prim,K);
  # from here on, we use only p = p.sec for calculations  
  #p <- pfemale.sec
    
  if(is.null(pmale.prim)) { 
    if(is.null(pmale.sec)) stop("either pmale.prim or pmale.sec need to be specified") else {
      pmale.prim = 1 - dbinom(x=0, size=K, prob=pmale.sec)
    }
  }
  if(is.null(pmale.sec)) pmale.sec = 1 -  nthroot(1-pmale.prim,K);
  # from here on, we use only p = p.sec for calculations  
  #p <- p.sec
  if(length(rmale)==1) rmale<- rep(rmale, times = T) else if(length(rmale) != T) stop("rmale needs to have length T");
  if(length(rfemale)==1) rfemale<- rep(rfemale, times = T) else if(length(rfemale) != T) stop("rfemale needs to have length T");
  
  if(length(smale)==1) smale<- rep(smale, times = T) else if(length(smale) != T) stop("smale needs to have length T");
  if(length(sfemale)==1) sfemale<- rep(sfemale, times = T) else if(length(sfemale) != T) stop("sfemale needs to have length T");
  if(length(pfemale.sec)==1) pfemale.sec<- rep(pfemale.sec, times = T) else if(length(pfemale.sec) != T) stop("pfemale.sec needs to have length T");
  if(length(pmale.sec)==1) pmale.sec<- rep(pmale.sec, times = T) else if(length(pmale.sec) != T) stop("pmale.sec needs to have length T");
  

  if(any(rmale  > 1)) {
    stop("r must be <= 1")  }
  if(any(rfemale  > 1)) {
    stop("r must be <= 1")  }
  if(any(taumale  > 1)) {
    stop("taumale must be <= 1")  }
  if(any(taufemale  > 1)) {
    stop("taufemale must be <= 1")  }
  
  # sp = s + tau[1];
  # s <- s/sp;
  # tau[1]<- tau[1]/sp 
  
  
  testsimmale = sim.data.robust.tau.time(N=Nmale, T=T, K=K, r=rmale, s=smale, p.sec=pmale.sec, tau = taumale)
  testsimfemale = sim.data.robust.tau.time(N=Nfemale, T=T, K=K, r=rfemale, s=sfemale, p.sec=pfemale.sec, tau = taufemale)
 
  # go on here!
   
   # sim.data.robust.tau.time(N, T, K=1, r, s, p.prim=NULL, p.sec=NULL, tau)
  
  # create storage variables
  attendance <- matrix(0, nrow=N, ncol=T) # alive an not temp emigrated
  encounter <- matrix(0, nrow=N, ncol=T)  # encounter is a matrix with T columns
  # p.state <- rep(0, times = T)
  capture.full <- NULL                    # capture.full is a list of T matrizes with each K columns
  capture <- NULL
  
  
  # attendance history
  # generate recruitment distribution
  
  
  recruitment    <- rmultinom(1, N, r) # produces a matrix with rows = T = lenght r , with the number of elements first recruited in each row
  
  # make vector of arrival times
  recruitment.i <- rep(1:T, recruitment) # a number for each individual # adapt for 2-gender version
  # now use
  state <- matrix(0, nrow=N, ncol=T) # state = 0, and 1 if alive
  
  # instead of state, I now use state.laake. state was in HW a list with attendance vector as elements for each secondary period
  state.laake <- matrix(0, nrow=N, ncol=T) # is C_t in Laake for each individual 
  # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
  # state = 1 means temp. emigrated or alive
  # attendance is 1 if alive and not temp emigrated!
  
  # loop over individuals
  for (i in 1:N)  {
    # record first attendance
    attendance[i,recruitment.i[i]]  <- 1  # recruitment.i[i] is the first atttendance time point
    state.laake[i,recruitment.i[i]] <- 1 
    state[i,recruitment.i[i]] <- 1 # state should be used later to for male/female, here: just one gender
    # state.laake: 0 not yet recruited. 1 recruited . 2 temp. emigrated. 3 departed
    # if not recruited in final year
    if (recruitment.i[i] < T)  {
      # loop over remaining occasions
      for (t in recruitment.i[i]:(T-1))  {
        # if present at site
        if (attendance[i,t] == 1)  {
          # work out whether they survive until the next year, or temp emigrate or die
          # recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
          survival <- rmultinom(1,1,c(s[t]*(1-tau[1]),s[t]*tau[1], 1-s[t] ) ) 
          surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
          # store 0 or 1 following year
          
          # state[i,t+1] <- survival;
          state.laake[i,t+1] <- surv
          
          if (surv == 3) {  # state.laake[i,t+1] <- 3 
            state[i,t+1] <- 0 # it is defaulted to 0 anyway
          } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
            state[i,t+1] <- 1;
            attendance[i,t+1] <- 1
          } else if(surv == 2){ # temp. emigrated
            #  attendance[i,t+1] <- 0 # is anyway the default
            state[i,t+1] <- 1; # state: 0 - dead. 1 male or no gender measured
            # here: state = 1 means temp. emigrated or alive
          }
          # attendance is 1 if alive and not temp emigrated!
        } else { # is either dead or temp emigrated
          if (state.laake[i,t] == 2){ # temp. emigrated
            # lets build in a death prob.: recoded 1/11/2019: death prob = 1-s[t],emig prob = s[t]*tau[2] 
            survival <- rmultinom(1,1,c(s[t]*(1-tau[2]),s[t]*tau[2], 1-s[t] ) ) 
            surv <- which(survival == 1)  #  1 = alive, 2 = temp. emig, 3 = dead
            # store 0 or 1 following year
            
            # state[i,t+1] <- survival;
            state.laake[i,t+1] <- surv
            
            if (surv == 3) {  # state.laake[i,t+1] <- 3 
              state[i,t+1] <- 0 # it is defaulted to 0 anyway
            } else if(surv == 1) {#  state.laake[i,t+1] <- 1; 
              state[i,t+1] <- 1;
              attendance[i,t+1] <- 1
            } else if(surv == 2){ # temp. emigrated
              #  attendance[i,t+1] <- 0 # is anyway the default
              state[i,t+1] <- 1;
              # here: state = 1 means temp. emigrated or alive
            }

          } else {  # ie dead -
            
            state.laake[i,t+1] <- 3  
            state[i,t+1] <- 0;
            # since this is developed from gender version, we have states: 0 (unborn)), 1 (alive), 3 (departed)
          }
        }
      }  
    }
  }
  # find number that are present each year
  N.t <- colSums(attendance)  # + add for each gender
  # N1.t <- colSums(attendance[state == 1 , ]);
  # N2.t <- colSums(attendance[state == 2 , ]);
  #  p.state <- N1.t / N.t  # true state proportion for each primary capture occassion
  # state.t <- matrix(0, nrow=N, ncol=1) # needed when we include the temp. emigration
  
  # capture histories for each year
  # loop over years
  for (t in 1:T)  {
    
    # storage variables for this year
    capture.t <- matrix(0, nrow=N, ncol=K)  
    
    # loop over individuals and check whether they were captured or not
    for (i in 1:N)  {
      # capture.t is defaulted to 0
      # find the primary occassions where they were present
      if(attendance[i,t]==1) {
        # loop over occasions
        capture.t[i,] <- rbinom(K,1,p[t])  # outputs K values of 0 and 1
        encounter[i,t] <- sum(capture.t[i,]) # i am only interested in no of captures # outputs the number between 0 and K
        
      }
    }
    # find which rows are all zeros in capture histories
    capture.0 <- which(rowSums(capture.t) == 0) # ie for which year are all secondary capt histories = 0
    
    # update encounter histories
    for (c in capture.0)  {       encounter[c,t] <- 0    } # should work anyway
    
    # store state history for year t
    #  state[[t]] <- state.t - will be needed if temp. emigration is also a state
    # store full capture history for year t
    capture.full[[t]] <-  capture.t # a matrix with N rows and K columns
    # store capture history for year t
    capture[[t]] <- capture.t[-capture.0,]
  }
  
  # count those encountered each year
  N.c <- colSums(encounter >= 1)  # to do: count for both male and female and all
  # count those who are seen at least once
  N.seen <- N-length(which(rowSums(encounter  >= 1) == 0))  # to do: count for both male and female and all
  
  # return
  return(list('attendance'=attendance, 'encounter'=encounter, 'state'=state, 'state.laake' = state.laake, # 'p.state' = p.state,
              'capture.full'=capture.full,  'capture'=capture, 'N.attend'=N.t, 'N.encounter'=N.c, 'N.seen'=N.seen))
}


