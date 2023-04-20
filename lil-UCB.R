library(rlist) #Library for working on lists, which will be how bandit problems will be stored
library(tictoc) #Library for timer

g <<- 0.1^2 #Declaring gamma as a global variable
#Alpha vectors for experiments
#alphas <<- c(1.5,1.5,1.5)
alphas <<- c(1,1.5,2)
#alphas <<- c(1.5,1.5,1.5)
#alphas <<- c(2,1.5,2,2.5,1)
#alphas <<- rep(1.5,5)

K <<- length(alphas) #Number of arms
#Generating a bandit arm in the rare events setting. Inputs are rarity parameters gamma and alpha, upper bound B, and number of nonzero support points n.
generate_rare_bandit <- function(gamma,alpha,B,n){
  s_v <- seq(1,(0.5*B),by = 0.1^4)
  v <- sort(sample(s_v,n,replace = FALSE))*(gamma^(-alpha))
  val <- c(0,v)
  
  s_p <- seq(0.1^3,0.3+0.1^5,by=0.1^4)
  p <- sample(s_p,n)*gamma^(alpha)
  prob <- c(1-sum(p),p)
  
  bandit <- rbind(val,prob)
  return(bandit)
}
#Computing the mean reward of an arm
mean_reward <- function(arm){
  return(sum(arm[1,]*arm[2,]))
}
#Generating a bandit instance, taking rarity parameters as input. Also returns the vector of upper bounds
generate_bandit_problem <- function(gamma,alpha_vector){
  bandits <- list()
  K <- length(alpha_vector)
  UBs <- c(K+7,K+4,seq(K+1,4,by=-1))
  for(i in 1:K){
    new <- generate_rare_bandit(gamma,alpha_vector[i],UBs[i],floor(1.1*(K+1-i)+20))
    name <- paste(i)
    bandits <- list.append(bandits,new)
  }
  bandits <- list.append(bandits,"UpperBounds"=UBs)
  return(bandits)
}
#Sampling an arm, taking as input the bandit model, the arm index and the number of samples
sample_arm <- function(bandit,arm_idx,num){
  val <- bandit[[arm_idx]][1,]
  probs <- bandit[[arm_idx]][2,]
  
  ans <- sample(val,num,replace=TRUE,prob=probs)
  
  return(ans)
}
#Given a vector of samples, the following function returns the empirical probability distribution
emp_dist <- function(x){
  v <- sort(unique(x))
  probs <- c()
  for (value in v){
    probs <- c(probs,mean(value==x))
  }
  return(rbind(v,probs))
}
#Combining a list of bandit arms and upper bounds to generate a bandit instance
combine_arms <- function(arm_list,UBs){
  return(list.append(arm_list,"UpperBounds"=UBs))
}
#Updating the empirical bandit model based on the newest samples seen. This is done ore arm at a time.
update_emp_bandit <- function(emp_bandit,num_samples,samples,arm_idx){
  arm <- emp_bandit[[arm_idx]]
  freq_list <- (arm[2,])*num_samples[arm_idx]
  v <- arm[1,]
  
  for(old in v){
    updated_old_freq <- sum(samples==old)
    idx <- which(v==old)[1]
    freq_list[idx] <- freq_list[idx] + updated_old_freq
  }
  fresh_obs <- unique(samples)
  new_obs <- setdiff(fresh_obs,v)
  for(new in new_obs){
    v <- c(v,new)
    new_freq <- sum(samples==new)
    freq_list <- c(freq_list,new_freq)
  }
  
  updated_num_samples <- num_samples
  updated_num_samples[arm_idx] <- num_samples[arm_idx]+length(samples)
  probs <- (1/updated_num_samples[arm_idx])*freq_list
  
  updated_arm <- rbind(v,probs)    
  #updated_arm <- updated_arm[,order(updated_arm[1,],decreasing=FALSE)]
  
  updated_emp_bandit <- emp_bandit
  updated_emp_bandit[[arm_idx]] <- updated_arm
  return(list("Bandit"=updated_emp_bandit,"NumSamples"=updated_num_samples))
}

#Setting a seed number for replicability
Seed <- 304
set.seed(Seed)
cat("\n Seed for Bandit Model = ",Seed,"\n")

bandit_prob <<- generate_bandit_problem(g,alphas) #Generating the bandit problem
B <<- bandit_prob$UpperBounds #The upper bounds

cat("\nMean of arm 1 = ",mean_reward(bandit_prob[[1]]),"\n")
for(i in 2:K){
  #print(bandit_prob[[i]])
  cat("\nMean of arm", i, "= ",mean_reward(bandit_prob[[i]]),"\n")
}

delta <<- 0.01 #Confidence level
#epsilon, beta and a parameters for lil-UCB
epsilon <<- 18.19
Beta <<- 1
a <<- ((Beta+2)/Beta)^2

#We'll run lil-UCB on one sample path of our bandit instance
set.seed(NULL)
seed_set <- 200:3000
Seed <- sample(seed_set,1)
set.seed(Seed)
cat("\n Seed for lil-UCB = ",Seed,"\n")

init_arms <- list() #List to store initial empirical arm distributions
num_samples <- c()
#Initial exploration
for(i in 1:K){
  alpha <- alphas[i]
  n <- 2
  num_samples <- c(num_samples,length(x))
  arm <- emp_dist(x)
  init_arms <- list.append(init_arms,arm)
}
t <- sum(num_samples)

emp_bandit <- combine_arms(init_arms,B)
#Setting the stopping rule
stop_thresholds <- 1 + a*(sum(num_samples)-num_samples)
stop_rule <- (sum(num_samples>=stop_thresholds)==0)
tic() #Timer starts
#lil-UCB iterations
while(stop_rule){
  x1 <- log(log((1+epsilon)*num_samples)/delta)
  x2 <- 2*(max(B*g^(-alphas)))*(1+epsilon)*x1/4
  x3 <- x2/num_samples
  ucb <- (1+Beta)*(1+sqrt(epsilon))*sqrt(x3)
  emp_means <- c()
  for(i in 1:K){
    emp_means <- c(emp_means,mean_reward(emp_bandit[[i]]))
  }
  indices <- emp_means+ucb
  best_ucb_arm <- sample(which(indices==max(indices)),1)
  #num_samp <- g^(-alphas[best_ucb_arm])
  num_samp <- 1
  x <- sample_arm(bandit_prob,best_ucb_arm,num_samp)
  u <- update_emp_bandit(emp_bandit,num_samples,x,best_ucb_arm)
  emp_bandit <- u$Bandit
  num_samples <- u$NumSamples
  t <- sum(num_samples)
  stop_thresholds <- 1 + a*(sum(num_samples)-num_samples)
  stop_rule <- (sum(num_samples>=stop_thresholds)==0)
}
toc() #Timer ends
