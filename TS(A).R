library(rlist) #Library for working on lists, which will be how bandit problems will be stored
library(rootSolve) #Library to solve equations using bisection method
library(tictoc) #Library for timer

#Computation of KL divergence exactly
KL_exact <- function(p,q){
  if(sum(q<=0)==0){
    return(sum(p*log(p/q))) #num of operations = O(length(p)=length(q))
  }else{cat("ERROR: Invalid arguments in KL divergence","\n")}   
} 
#Computation of K_inf(L). Inputs are the bandit model, the vector of arms' powers, the arm index, and the number with respect to which computation is to be done.
#Dual form is used for computation
K_inf_L <- function(bandit,powers,i,x){
  arm <- bandit[[i]]
  ai <- arm[1,]
  pi <- arm[2,]
  #p0 <- as.numeric(arm[2,1])
  dual_mean_eq <- function(lambda){
    t1 <- (ai-x)*pi #O(length(ai)) many computations
    t2 <- 1-lambda*(x-ai) #O(length(ai)) many computations
    return(-sum(t1/t2))
  }
  #print(dual_mean_eq(1/x))
  if(x>=mean_reward(arm)){
    return(0)
  }else{
    u_smart <- smart_ub(dual_mean_eq,0,1/x)
    #The uniroot function implements the bisection method
    opt.lambda <- uniroot(dual_mean_eq,c(0,u_smart),tol=tolerance)$root
    #print(opt.lambda)
    ans <- sum(pi*log(1-opt.lambda*(x-ai))) #O(length(ai))*(-log(tolerance)) many computations
    if(ans<0 && ans>-(tolerance)){
      return(0)
    }else{
      return(ans)
    }
  }
} 
#Computation of K_inf(U). Inputs are the bandit model, the vector of arms' powers, the arm index, and the number with respect to which computation is to be done.
#Dual form is used for computation
K_inf_U <- function(bandit,powers,i,x){
  arm <- bandit[[i]]
  ai <- arm[1,]
  pi <- arm[2,]
  dual_mean_eq <- function(lambda){
    t1 <- pi*(x-ai) #O(length(ai)) many computations
    t2 <- 1+lambda*(x-ai) #O(length(ai)) many computations
    return(sum(t1/t2))
  }
  upper.bd <- 1/(B[i]*(g^(-powers[i]))-x)
  if(x<=mean_reward(arm)){
    opt.lambda <- 0
  }else if(dual_mean_eq(0)*dual_mean_eq(upper.bd)>=0){
    opt.lambda <- upper.bd
  }else{
    opt.lambda <- uniroot(dual_mean_eq,c(0,upper.bd),tol=tolerance)$root #O(length(ai))*(-log(tolerance)) many computations
  }
  #print(opt.lambda)
  ans <- sum(pi*log(1+opt.lambda*(x-ai))) #O(length(ai)) many computations
  
  if(ans<0 && ans>-(tolerance)){
    return(0)
  }else{
    return(ans)
  }
  
}
#Computation of KL divergence approximately
KL_approx <- function(p,q){
  px <- p[-1]
  qx <- q[-1]
  if(sum(q<=0)==0){
    return(abs(sum(px*log(px/qx)) + sum(qx-px))) #O(length(ai)) many computations
  }else{cat("ERROR: Invalid arguments in KL divergence","\n")}   
} 
#Computation of K_inf(L) approximately using C constants
KL_approx_L <- function(gamma,alpha,arm,C_1){
  ai <- arm[1,-1]*gamma^(alpha)
  pi <- arm[2,-1]*gamma^(-alpha)
  tilde_pi <- (pi/(1+(C_1*ai)))*gamma^(alpha)
  tilde_pi_dist <- c((1-sum(tilde_pi)),tilde_pi)
  return(KL_approx(arm[2,],tilde_pi_dist)) #O(length(ai)) many computations
} 
#Computation of K_inf(U) approximately using C constants
KL_approx_U <- function(gamma,alpha1,arm1,i,alphai,armi,C_1,C,Bounds){#to be called only after find_Ci is used
  a1 <- arm1[1,-1]*gamma^(alpha1)
  p1 <- arm1[2,-1]*gamma^(-alpha1)
  ai <- armi[1,-1]*gamma^(alphai)
  pi <- armi[2,-1]*gamma^(-alphai)
  tilde_p1i <- p1/(1+(C_1*a1))
  tilde_mu <- sum(a1*tilde_p1i)
  if(C==(1/Bounds[i])){
    tilde_pi <- (pi/(1-(C*ai)))
    ans <- sum(pi*(log(pi/tilde_pi)))+C*(tilde_mu) #O(length(ai)) many computations
    return(abs(ans)*gamma^(alphai))
  }else{
    tilde_pi <- abs((pi/(1-(C*ai)))*gamma^(alphai))
    if(sum(tilde_pi)<=1){
      tilde_pi_dist <- c(1-sum(tilde_pi),tilde_pi)
    }else{
      cat("ERROR! CANNOT FINF ALTERNATIVE DISTRIBUTION FOR K_inf_U approx!! \n")
    }
    
    return(KL_approx(armi[2,],tilde_pi_dist)) #O(length(ai)) many computations
  } 
}
#smart selection of upper bounds for bisection method when function value at theoretical upper bound tends to infinity
smart_ub <- function(f_name,lb,ub){
  #flag variable to see extra output text while debugging
  repair <- FALSE
  if(repair){
    cat("smart_ub START \n")
  }
  #print(f_name)
  lb_sign <- sign(do.call(f_name,list(lb)))
  ub_value <- do.call(f_name,list(ub))
  ub_sign <- sign(ub_value)
  if(lb_sign*ub_sign<=0){
    if(ub_value<Inf){
      if(repair){
        cat("smart_ub END \n")
      }
      return(ub)
    }else{
      pow <- 1
      mid <- (0.5^pow)*(lb+ub)
      mid_sign <- sign(do.call(f_name,list(mid)))
      while(lb_sign*mid_sign>0){
        pow <- 1
        mid <- (1-0.5^pow)*mid+(0.5^pow)*ub
        mid_sign <- sign(do.call(f_name,list(mid)))
      }
      if(repair){
        cat("smart_ub END \n")
      }
      return(mid)
    }
  }else{
    cat("ROOT DOES NOT EXIST IN INTERVAL PROVIDED, CANNOT FIND POINT OF SIGN CHANGE", "\n")
    if(repair){
      cat("smart_ub END \n")
    }
    print(f_name)
  }            
}
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
    idx <- which(v==old)
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
  updated_arm <- updated_arm[,order(updated_arm[1,],decreasing=FALSE)]
  
  updated_emp_bandit <- emp_bandit
  updated_emp_bandit[[arm_idx]] <- updated_arm
  return(list("Bandit"=updated_emp_bandit,"NumSamples"=updated_num_samples))
}

#Alpha vectors for experiments
#alphas <<- c(1,1.5,2)
#alphas <<- c(1.5,1.5,1.5)
#alphas <<- c(2,1.5,2,2.5,1)
#alphas <<- rep(1.5,5)
#alphas <<- c(2,1.5,1)

#Setting a seed number for replicability
Seed <- 304
set.seed(Seed)
cat("\n Seed for Bandit Model = ",Seed,"\n")
g <<- 0.1^2 #Declaring gamma as a global variable

alphas <<- c(1,1.5,2)

tolerance <<- 0.1^12 #Tolerance level for solving equations by bisection method

K <<- length(alphas) #Number of arms
bandit_prob <<- generate_bandit_problem(g,alphas) #Generating the bandit problem
B <<- bandit_prob$UpperBounds #The upper bounds
F_0 <<- c() #Thresholds to determine whether or not the extra support point would be necessary for K_inf(U)
for (i in 2:K){
  arm <- bandit_prob[[i]]
  ai <- arm[1,-1]*(g^(alphas[i]))
  pi <- arm[2,-1]*(g^(-alphas[i]))
  xi <- B[i]/((sum((ai*pi)/(B[i]-ai)))^(-1))
  F_0 <<- c(F_0,xi)
}

cat("\nMean of arm 1 = ",mean_reward(bandit_prob[[1]]),"\n")
for(i in 2:K){
  #print(bandit_prob[[i]])
  cat("\nMean of arm", i, "= ",mean_reward(bandit_prob[[i]]), " Threshold = ",F_0[i-1],"\n")
}

#The three functions that follow will be useful in case the first two empirical arms are not empirically the best and the second best
#As per our code, solving the lower bound problem requires the first two empirical arms to be empirically the best and the second best
#We will be required to solve the lower bound problem for the empirical bandits we see after each round of sampling

#Switching the indices of i1-th and i2-th arms in the bandit model that has been input
swap_arm <- function(bandit,i1,i2){        
  switched_bandit <- bandit
  switched_UBs <- bandit$UpperBounds
  
  switched_bandit[[i1]] <- bandit[[i2]]
  switched_bandit[[i2]] <- bandit[[i1]]
  
  switched_UBs[i1] <- bandit$UpperBounds[i2]
  switched_UBs[i2] <- bandit$UpperBounds[i1]
  
  switched_bandit$UpperBounds <- switched_UBs
  return(switched_bandit)
}
#Given a vector v, the next function makes switches between the pairs of indices (idx11,idx12) and (idx21,idx22)
switch_vector <- function(v,idx11,idx12,idx21,idx22){    
  new_v <- v
  
  if((idx11 == idx21)||(idx12 == idx22)){
    cat("INCORRECT INDICES ENTERED FOR switch_vector. \n")
    return()
  }
  if(idx11!=idx12){
    new_v[idx11] <- v[idx12]
    new_v[idx12] <- v[idx11]
  }
  two_cycle <- (idx21==idx12)&&(idx22==idx11)
  if(two_cycle){
    return(new_v)
  }else{
    if(idx21!=idx22){
      if(idx21==idx12){
        dummy <- new_v[idx12]
        new_v[idx12] <- new_v[idx22]
        new_v[idx22] <- dummy
      }else if(idx22==idx11){
        dummy <- new_v[idx12]
        new_v[idx12] <- new_v[idx21]
        new_v[idx21] <- dummy   
      }else{
        new_v[idx21] <- v[idx22]
        new_v[idx22] <- v[idx21]
      }
      return(new_v)
    }
  } 
  return(v)
}
#Putting the two empirically best arms at indices 1 and 2.
#After solving the lower bound problem, the indices are switched back to what they were in the given bandit model.
put_best_first <- function(emp_bandit){#ensures that emp_bandit has first and second arms as best and second best          
  emp_means <- c()
  for(i in 1:K){
    emp_means <- c(emp_means,mean_reward(emp_bandit[[i]]))
  }
  
  best_emp_mean <- max(emp_means)
  best_emp_mean_idx <- which(emp_means==best_emp_mean)
  swap1 <- swap_arm(emp_bandit,1,best_emp_mean_idx)
  #emp_means[1] <- mean_reward(swap1[[1]])
  #emp_means[best_emp_mean_idx] <- mean_reward(swap1[[best_emp_mean_idx]])
  
  scd_best_emp_mean <- max(emp_means[-best_emp_mean_idx])
  scd_best_emp_mean_idx <- which(emp_means==scd_best_emp_mean)
  
  dummy_emp_means <- c()
  for(i in 1:K){
    dummy_emp_means <- c(dummy_emp_means,mean_reward(swap1[[i]]))
  }
  dummy_idx <- which(dummy_emp_means==scd_best_emp_mean)
  swap2 <- swap_arm(swap1,2,dummy_idx)
  
  sorted_alphas <- switch_vector(alphas,1,best_emp_mean_idx,2,scd_best_emp_mean_idx)
  
  ans <- list("SortedBandit"=swap2,"BestEmpArm"=best_emp_mean_idx,"SecondBestEmpArm"=scd_best_emp_mean_idx,"Alphas"=sorted_alphas)
  return(ans)
}

#Solving the inner minimization problem
#Inputs are the vector of weights, the bandit instance, the vector of rarities and the index of the arm with which the empirical best arm will contest
solve_inner_opt <- function(weights,bandit,alphas,idx){
  
  arm1 <- bandit[[1]]
  alpha1 <- alphas[1]
  arm2 <- bandit[[idx]]
  alpha2 <- alphas[idx]
  
  bounds <- bandit$UpperBounds
  
  a1 <- arm1[1,-1]*g^(alpha1)
  p1 <- arm1[2,-1]*g^(-alpha1)
  
  a2 <- arm2[1,-1]*g^(alpha2)
  p2 <- arm2[2,-1]*g^(-alpha2)
  
  eval <- function(C_1i){
    t1 <- sum((a1*p1)/(1+C_1i*a1))
    Ci <- (C_1i*weights[1]*g^(alpha1-alpha2))/weights[2]
    t2 <- sum((a2*p2)/(1-Ci*a2))
    
    return(t1/t2 - 1)        
  }
  
  l <- 0
  u <- min(sum(p1)/sum(a2*p2),(weights[2]*g^(alpha2-alpha1))/(weights[1]*max(a2)))
  C_1i <- uniroot(eval,c(l,u),tol=tolerance)$root
  Ci <- (C_1i*weights[1]*g^(alpha1-alpha2))/weights[2]
  
  if(Ci > 1/bounds[idx]){
    Ci <- 1/bounds[idx]
    C_1i <- (Ci*weights[2]*g^(alpha2-alpha1))/weights[1]
  }
  return(c(C_1i,Ci))
}
#Computing the weighted sum of K_inf quantities
wt_KL_exact <- function(w,bandit_problem,powers,i,x){
  ans <- w[1]*K_inf_L(bandit_problem,powers,1,x)+w[i]*K_inf_U(bandit_problem,powers,i,x) #O(length(ai)*(-log(tolerance)) many computations
  return(ans)
}
#Solving the inner minimization problem exactly
#Inputs are the vector of weights, the bandit instance, the vector of rarities and the index of the arm with which the empirical best arm will contest
#This will be required in the computation of theoretical lower bound on sample complexity
solve_inner_opt_exact <- function(w,bandit_problem,powers,i){
  eval <- function(x){
    return(wt_KL_exact(w,bandit_problem,powers,i,x)) #O(length(ai)*(-log(tolerance)) many computations
  }
  mu1 <- mean_reward(bandit_problem[[1]])
  mui <- mean_reward(bandit_problem[[i]])
  init <- (mu1+mui)/2
  ans <- optim(init,eval,method="Brent",lower=mui,upper=mu1)
  #O(length(ai)*(-log(tolerance))*(log(mu1-mui)-log(tolerance))^2) many computations
  if(ans$value <= min(w[1]*K_inf_L(bandit_problem,powers,1,mui),w[i]*K_inf_U(bandit_problem,powers,i,mu1),wt_KL_exact(w,bandit_problem,powers,i,(mu1+mui)/2))){
    return(list("value"=ans$value,"optx"=ans$par))
  }else{
    cat("Inner optimization not solved correctly!!")
    return()
  }
}

#Computing likelihood ratio exactly
Z_exact <- function(emp_bandit,num_samples){
  sorted <- put_best_first(emp_bandit)
  sorted_emp_bandit <- sorted$SortedBandit
  sorted_alphas <- sorted$Alphas
  sorted_UBs <- sorted_emp_bandit$UpperBounds
  
  switched_num_samples <- switch_vector(num_samples,1,sorted$BestEmpArm,2,sorted$SecondBestEmpArm)
  
  wt <- switched_num_samples/sum(switched_num_samples)
  
  ans <- c()
  for(i in 2:K){
    x <- solve_inner_opt_exact(wt,sorted_emp_bandit,sorted_alphas,i)$optx
    zee <- wt[1]*K_inf_L(sorted_emp_bandit,sorted_alphas,1,x)+wt[i]*K_inf_U(sorted_emp_bandit,sorted_alphas,i,x)
    zee <- zee*sum(switched_num_samples)
    ans <- c(ans,zee)
  }
  
  return(min(ans))
}
#Stopping rule threshold
beta_alt <- function(N,delta){
  ans <- log((K-1)/delta)+4*log(log(N+1))+2
  return(ans)
}
#Solving the outer maximization problem
solve_LB <- function(emp_bandit){
  repair <- FALSE
  if(repair){cat("\n Lower bound solving START \n")}
  sorted <- put_best_first(emp_bandit)
  sorted_emp_bandit <- sorted$SortedBandit
  sorted_alphas <- sorted$Alphas
  sorted_UBs <- sorted_emp_bandit$UpperBounds
  
  sorted_temp_F_0 <- c()
  for (i in 2:K){
    arm <- sorted_emp_bandit[[i]]
    ai <- arm[1,-1]*(g^(sorted_alphas[i]))
    pi <- arm[2,-1]*(g^(-sorted_alphas[i]))
    xi <- sorted_UBs[i]/((sum((ai*pi)/(sorted_UBs[i]-ai)))^(-1))
    sorted_temp_F_0 <- c(sorted_temp_F_0,xi)
  }
  
  h <- function(w1,w2,p1,p,q1,q){
    return(w1*KL_approx(p1,p)+w2*KL_approx(q1,q))
  }
  
  wt_KL_C <- function(gamma,arm1,alpha1,armi,alphai,C_1,C){
    a_1 <- arm1[1,-1]*gamma^(alpha1)
    p_1 <- arm1[2,-1]*gamma^(-alpha1)
    a_i <- armi[1,-1]*gamma^(alphai)
    p_i <- armi[2,-1]*gamma^(-alphai)
    if(C>0){
      wt_sum <- sum(p_1*log(1+C_1*a_1))+(C_1/C)*sum(p_i*log(1-C*a_i)) #O(length(ai)) many computations
    }else{
      wt_sum <- sum(p_1*log(1+C_1*a_1))-C_1*sum(p_i*a_i)
    }
    
    if(abs(wt_sum)<2*tolerance){
      wt_sum <- abs(wt_sum)
    }
    
    return(wt_sum)
  }
  
  C_1i_ub <- function(gamma,arm1,alpha1,armi,alphai){
    eval<- function(x){
      mu_i <- mean_reward(armi)
      a1 <- arm1[1,-1]*(gamma^(alpha1))
      p1 <- arm1[2,-1]*(gamma^(-alpha1))
      LHS <- sum((a1*p1)/(1+x*a1))
      RHS <- mu_i
      return(LHS-RHS)
    }
    ub <- 1/tolerance
    #pow <- 1
    #while(eval(ub)>=0){
    #    pow <- pow+1
    #    ub <- ub^pow
    #}
    ans <- uniroot(eval,c(0,ub),tol=tolerance)$root #O(length(ai)*(-log(tolerance))) many computations
    return(ans)
  }
  
  find_Ci <- function(gamma,arm1,alpha1,i,armi,alphai,C1){
    if(repair){
      cat("find_Ci START \n")
      cat("\t C_1i = ",C1,"\n")
    }
    a1 <- arm1[1,-1]*(gamma^(alpha1))
    ai <- armi[1,-1]*(gamma^(alphai))
    p1 <- arm1[2,-1]*(gamma^(-alpha1))
    pi <- armi[2,-1]*(gamma^(-alphai))
    
    tilde_p1i <- (p1/(1+(C1*a1)))
    tilde_mu <- (sum(a1*tilde_p1i))
    #print(c(F_0[i-1],tilde_mu))
    
    if(sorted_temp_F_0[i-1]<min(sum(a1*p1),tilde_mu)){
      if(repair){cat("find_Ci END \n")}
      return(1/sorted_UBs[i])
    }else{
      chk <- C_1i_ub(gamma,arm1,alpha1,armi,alphai)
      if(repair){
        cat("\t C_1i upper bound = ",chk,"\n")
      }
      if(chk<C1){ #O(length(ai)*(-log(tolerance))) many computations
        if(repair){
          cat("\t","Invalid C_1i entered! The upper bound was ",chk,"\n")
          cat("find_Ci END \n")
        }
        return()
      }
      
      C.min <- 0
      C.max <- min(1/ai)
      Ci_solver <- function(C){
        
        if(C1==chk && C == 0){
          return(0)
        }else{
          LHS <- sum((a1*p1)/(1+C1*a1))
          RHS <- sum((ai*pi)/(1-C*ai))
          return(RHS-LHS) #O(length(ai)) many computations
        }
        
      }
      
      if(repair){
        cat("\t Value at lower bound, i.e. zero =",Ci_solver(C.min),"\n")            
        cat("\t Upper bound on C_i =",C.max,"\n")
        cat("\t Value at upper bound =",Ci_solver(C.max),"\n")
      }
      
      if(Ci_solver(C.max)<=Inf){
        C.ub <- smart_ub(Ci_solver,C.min,C.max)
        C <- uniroot(Ci_solver,c(C.min,C.ub),tol=tolerance)$root 
        #O(length(ai)*(-log(tolerance))) many computations
        if(repair){cat("find_Ci END \n")}
        return(C)
      }else{
        cat("\t Possibla NaNs!! Check Ci_solver \n")
        if(repair){
          cat("find_Ci END \n")
        }
        return()
      }
      
    }
  }
  
  wt_KL_C1 <- function(gamma,arm1,alpha1,i,armi,alphai,C_1){
    if(repair){cat("wt_KL_C1 START \n")}
    C <- find_Ci(gamma,arm1,alpha1,i,armi,alphai,C_1) #O(length(ai)*(-log(tolerance))) many computations
    if(repair){cat("\t C_i = ",C,"\n")}
    ans <- wt_KL_C(gamma,arm1,alpha1,armi,alphai,C_1,C) #O(length(ai)) many computations
    if(repair){cat("wt_KL_C1 END \n")}
    return(ans)
  }
  
  find_C1i <- function(gamma,arm1,alpha1,arm2,alpha2,i,armi,alphai,C_12){
    if(repair){cat("find_C1i START \n")}
    LHS <- wt_KL_C1(gamma,arm1,alpha1,2,arm2,alpha2,C_12) #O(length(ai)*(-log(tolerance))) many computations
    if(repair){cat("\t LHS of C_1i = ",LHS,"\n")}
    C_1i_solver <- function(x){            
      RHS <- wt_KL_C1(gamma,arm1,alpha1,i,armi,alphai,x) #O(length(ai)*(-log(tolerance))) many computations
      if(repair){cat("\t RHS of C_1i = ",RHS,"\n")}
      return(RHS-LHS)
    }
    C_1i.min <- 0
    C_1i.max <- C_1i_ub(gamma,arm1,alpha1,armi,alphai)
    chk_C_1i_solver <- C_1i_solver(C_1i.max)
    if(repair){
      cat("\t Theoretical upper bound on C_1i is",C_1i.max,"\n")
      cat("\t Value at theoretical upper bound is: ",chk_C_1i_solver,"\n")
    }
    
    if(chk_C_1i_solver<=Inf){
      C_1i.ub <- smart_ub(C_1i_solver,C_1i.min,C_1i.max)
      if(repair){
        cat("\t Selected upper bound on C_1i is",C_1i.ub,"\n")
        cat("\t Values at 0 and selected bound are: ",c(C_1i_solver(C_1i.min),C_1i_solver(C_1i.ub)),"\n")
      }
      C_1i <- uniroot(C_1i_solver,c(C_1i.min,C_1i.ub),tol=tolerance)$root 
      #O(length(ai)*(-log(tolerance))^2) many computations
      if(repair){cat("find_C1i END \n")}
      return(C_1i)
    }else{
      cat("\t Possibla NaNs!! Check C_1i_solver \n")
      if(repair){cat("find_C1i END \n")}
    }        
  }
  
  KL_ratio <- function(gamma,arm1,alpha1,i,armi,alphai,C_1i){
    if(repair){cat("KL_ratio START \n")}
    C_i <- find_Ci(gamma,arm1,alpha1,i,armi,alphai,C_1i) #O(length(ai)*(-log(tolerance))) many computations
    up <- KL_approx_L(gamma,alpha1,arm1,C_1i) #O(length(ai)) many computations
    down <- KL_approx_U(gamma,alpha1,arm1,i,alphai,armi,C_1i,C_i,sorted_UBs) #O(length(ai)) many computations
    if(repair){cat("KL_ratio END \n")}
    return(up/down)
  }
  
  KL_ratio_sum <- function(C_12){
    if(repair){cat("KL_ratio_sum START \n")}
    
    arm1 <- sorted_emp_bandit[[1]]
    alpha1 <- sorted_alphas[1]
    arm2 <- sorted_emp_bandit[[2]]
    alpha2 <- sorted_alphas[2]
    s <- KL_ratio(g,arm1,alpha1,2,arm2,alpha2,C_12) #O(length(ai)*(-log(tolerance))) many computations
    if(repair){cat("Current C_12 = ",C_12,"\n")}
    for(i in 3:K){
      if(repair){cat("\t Arm index = ",i,"\n")}
      armi <- sorted_emp_bandit[[i]]
      alphai <- sorted_alphas[i]
      C_1i <- find_C1i(g,arm1,alpha1,arm2,alpha2,i,armi,alphai,C_12)
      if(repair){cat("\t C_1", i, "=",C_1i,"\n")}
      #O(length(ai)*(-log(tolerance))^2) many computations
      summand <- KL_ratio(g,arm1,alpha1,i,armi,alphai,C_1i)
      #O(length(ai)*(-log(tolerance))) many computations
      s <- s+summand
    } #O(K*length(ai)*(-log(tolerance))^2) many computations
    if(repair){cat("KL_ratio_sum END \n")}
    return(s-1)
  }
  
  
  C_12.min <- 0
  arm1 <- sorted_emp_bandit[[1]]
  alpha1 <- sorted_alphas[1]
  p_1 <- arm1[2,-1]*g^(-alpha1)
  mu_2 <- mean_reward(sorted_emp_bandit[[2]])
  C_12.max <- C_1i_ub(g,arm1,alpha1,sorted_emp_bandit[[2]],sorted_alphas[2])
  if(repair){cat("Upper bound on C_12 is",C_12.max,"\n")}
  if(KL_ratio_sum(C_12.max)<=Inf){
    C_12.ub <- smart_ub(KL_ratio_sum,C_12.min,C_12.max)
  }else{
    cat("Possible NaNs!! Check KL_ratio_sum \n")
  } 
  #C_12.max <- 0.63*C_1i_ub(g,arm1,alpha1,sorted_emp_bandit[[2]],sorted_alphas[2])
  #print(KL_ratio_sum((1-g)*C_1i_ub(g,arm1,alpha1,sorted_emp_bandit[[2]],sorted_alphas[2])))
  
  #tic()
  opt_C_12 <- uniroot(KL_ratio_sum,c(C_12.min,C_12.ub),tol=tolerance)$root
  #O(K*length(ai)*(-log(tolerance))^3) many computations
  opt_C1 <- c(opt_C_12)
  opt_C <- c(find_Ci(g,arm1,alpha1,2,sorted_emp_bandit[[2]],sorted_alphas[2],opt_C_12)) #O(length(ai)*(-log(tolerance))) many computations
  
  for(i in 3:K){
    y1 <- find_C1i(g,arm1,alpha1,sorted_emp_bandit[[2]],sorted_alphas[2],i,sorted_emp_bandit[[i]],sorted_alphas[i],opt_C_12)
    #O(length(ai)*(-log(tolerance))^2) many computations
    opt_C1 <- c(opt_C1,y1) 
    y2 <- find_Ci(g,arm1,alpha1,i,sorted_emp_bandit[[i]],sorted_alphas[i],y1)
    #O(length(ai)*(-log(tolerance))) many computations
    opt_C <- c(opt_C,y2)
  }#O(K*length(ai)*(-log(tolerance))^2) many computations
  
  wi_by_w1 <- (opt_C1/opt_C)*g^(sorted_alphas[1]-sorted_alphas[-1])
  opt_w1 <- 1/(1+sum(wi_by_w1))
  opt_w <- c(opt_w1,opt_w1*wi_by_w1)
  #toc()
  
  final_opt_w <- switch_vector(opt_w,sorted$BestEmpArm,1,sorted$SecondBestEmpArm,2)
  
  
  if(repair){cat("\n Lower bound solving END \n")}
  return(final_opt_w)
}
#Computing the theoretical lower bound on sample complexity at a given confidence level delta
compute_cpx <- function(emp_bandit,delta){
  repair <- FALSE
  if(repair){cat("\n Complexity computation START \n")}
  sorted <- put_best_first(emp_bandit)
  sorted_emp_bandit <- sorted$SortedBandit
  sorted_alphas <- sorted$Alphas
  sorted_UBs <- sorted_emp_bandit$UpperBounds
  
  
  Z <- function(i,y){
    if(repair){cat("\n Z function START\n")}
    mu1 <- mean_reward(sorted_emp_bandit[[1]])
    mui <- mean_reward(sorted_emp_bandit[[i]]) 
    
    if(y==0){
      return(list("value"=0,"xstar"=mu1))
    }
    if(y==Inf){
      return(list("value"=K_inf_L(sorted_emp_bandit,sorted_alphas,1,mui),"xstar"=mui)) #O(length(ai)*(-log(tolerance)) many computations
    }
    
    eval_x <- function(x){
      return(K_inf_L(sorted_emp_bandit,sorted_alphas,1,x)+y*K_inf_U(sorted_emp_bandit,sorted_alphas,i,x))
    }
    init <- (mu1+mui)/2 
    ans <- optim(init,eval_x,method="Brent",lower=mui,upper=mu1,control=list(pgtol=g^(max(alphas)+1)))
    #O(length(ai)*(-log(tolerance))*(log(mu1-mui)-log(tolerance))^2) many computations
    #print(ans)
    if(repair){cat("\n Z function END\n")}
    return(list("value"=ans$value,"xstar"=ans$par))
  }
  
  Z_inv <- function(i,t){
    if(repair){cat("\n Z_inv function START\n")}
    eval_y <- function(y){
      return(g^(-max(alphas))*(Z(i,y)$value-t))
      #O(length(ai)*(-log(tolerance))*(log(mu1-mui)-log(tolerance))^2) many computations
    }
    eval_y_tfm <- function(y2){
      y <- g^(-max(alphas))*y2
      return(eval_y(y))
    }
    mui <- mean_reward(sorted_emp_bandit[[i]])
    if(t>K_inf_L(sorted_emp_bandit,sorted_alphas,1,mui)){
      cat("Invalid t entered for Z^(-1)!! \n")
    }else if(t==K_inf_L(sorted_emp_bandit,sorted_alphas,1,mui)){
      return(Inf)
    }else{
      #print(c(eval_y(0),eval_y(0.9*K_inf_L(sorted_emp_bandit,sorted_alphas,1,mui))))
      
      pow <- -15
      while(eval_y(10^(pow))<=0){
        pow <- pow + 1
      }
      u <- 10^(pow)
      u_tfm <- u*g^max(alphas)
      #print(pow)
      
      ans <- uniroot(eval_y,c(0,u),tol=tolerance)$root
      #O(length(ai)*(-log(tolerance))^2*(log(mu1-mui)-log(tolerance))^2) many computations
      #ans <- uniroot(eval_y,c(0,1),tol=tolerance)$root
      if(repair){cat("\n Z_inv function END\n")}
      return(ans)
    }
  }
  
  KL_ratio_sum <- function(wt_KL_value){
    if(repair){cat("KL_ratio_sum START \n")}
    sum <- 0
    for(i in 2:K){
      #print(i)
      yi <- Z_inv(i,wt_KL_value)
      xi <- Z(i,yi)$xstar
      sum <- sum + K_inf_L(sorted_emp_bandit,sorted_alphas,1,xi)/K_inf_U(sorted_emp_bandit,sorted_alphas,i,xi)
      #cat(xi,",",yi,",",sum,"\n")
    }#O(K*length(ai)*(-log(tolerance))^2*(log(mu1-mui)-log(tolerance))^2) many computations  
    if(repair){cat("KL_ratio_sum END \n")}
    return(sum-1)
  }
  
  mu2 <- mean_reward(sorted_emp_bandit[[2]])
  
  tic()
  
  theoretical_ub <- K_inf_L(sorted_emp_bandit,sorted_alphas,1,mu2)
  if(KL_ratio_sum(theoretical_ub)<= Inf){
    u <- smart_ub(KL_ratio_sum,0,theoretical_ub)
  }else{
    cat("Possible NaNs!! Check KL_ratio_sum \n")
  } 
  
  opt_wt_KL_value <- uniroot(KL_ratio_sum,c(0,u),tol=tolerance)$root #might need smart_UB
  #O(K*length(ai)*(-log(tolerance))^3*(log(mu1-mui)-log(tolerance))^2) many computations
  
  wi_by_w1 <- c()
  for(i in (2:K)){
    opt_Z_inv <- Z_inv(i,opt_wt_KL_value)
    #print(g^(-alphas[i])*(Z(i,opt_Z_inv)$value-opt_wt_KL_value))
    wi_by_w1 <- c(wi_by_w1,opt_Z_inv)
  }#O(K*length(ai)*(-log(tolerance))^2*(log(mu1-mui)-log(tolerance))^2) many computations 
  opt_w1 <- 1/(1+sum(wi_by_w1))
  opt_w <- c(opt_w1,opt_w1*wi_by_w1)
  final_opt_w <- switch_vector(opt_w,1,sorted$BestEmpArm,2,sorted$SecondBestEmpArm)
  toc()
  
  
  x <- solve_inner_opt_exact(opt_w,sorted_emp_bandit,sorted_alphas,2)$optx
  T_star <- (wt_KL_exact(opt_w,sorted_emp_bandit,sorted_alphas,2,x))^(-1)
  
  if(repair){cat("\n Complexity computation END \n")}
  return(T_star*log(1/delta))
}

delta <<- 0.01 #Confidence level
sol <- solve_LB(bandit_prob)
T <- compute_cpx(bandit_prob,delta)

cat("Optimal weights for current current bandit problem: \n",sol,"\n")

#The TS(A) algorithm
#Initial exploration is done until at least one nonzero reward is seen in each arm
#Inner optimization is solved approximately (in batch sizes of $\gamma^(-\alpha_{max})$), but likelihood ratio is computed exactly
ApproxAlgo <- function(Seed){
  set.seed(NULL)
  set.seed(Seed)
  cat("\n Seed for Approx Algo = ",Seed,"\n")
  init_arms <- list()
  num_samples <- c()
  successes <- rep(0,K)
  #Initial Exploration
  for(i in 1:K){
    alpha <- alphas[i]
    n <- as.integer((g^(-alpha)+1))
    x <- sample_arm(bandit_prob,i,n)
    successes[i] <- successes[i]+sum(x!=0)
    while(sum(x)==0){
      x <- c(x,sample_arm(bandit_prob,i,n)) 
      successes[i] <- successes[i]+sum(x!=0)
    }
    num_samples <- c(num_samples,length(x))
    arm <- emp_dist(x)
    init_arms <- list.append(init_arms,arm)
  }
  emp_bandit <- combine_arms(init_arms,B)
  num_samples
  emp_bandit
  m <- as.integer(g^(-max(alphas)))+1
  l <- 1
  stop <- FALSE
  #ALGORITHM
  while(stop==FALSE){
    w_hat <- solve_LB(emp_bandit)  
    stopping_index <- Z_exact(emp_bandit,num_samples)
    threshold <- beta_alt((l*m),delta)
    
    if(stopping_index>threshold){
      stop <- TRUE
    }else{    
      stv <- c()
      for(i in 1:K){
        nn <- as.integer(sqrt((l+1)*m)-num_samples[i])
        stv <- c(stv,max(0,nn))
      }
      if(sum(stv)<=m){
        for(i in 1:K){
          if(stv[i]>0){
            x <- sample_arm(bandit_prob,i,stv[i])
            successes[i] <- successes[i]+sum(x!=0)
            u <- update_emp_bandit(emp_bandit,num_samples,x,i)
            emp_bandit <- u$Bandit
            #print(emp_bandit)
            num_samples <- u$NumSamples
            #print(num_samples)
          }
        }
        if(sum(stv)<m){
          diff <- m - sum(stv)
          s <- sample(1:K,diff,prob=w_hat,replace=TRUE)
          for(j in 1:K){
            num <- sum(s==j)
            if(num>0){
              x <- sample_arm(bandit_prob,j,num)
              successes[j] <- successes[j]+sum(x!=0)
              u <- update_emp_bandit(emp_bandit,num_samples,x,j)
              emp_bandit <- u$Bandit
              num_samples <- u$NumSamples
            }
          }
        }
      }else{
        stv_wt <- stv/sum(stv)
        s <- m*stv_wt
        for(i in 1:K){
          num <- s[i]
          if(num>0){
            x <- sample_arm(bandit_prob,i,num)
            successes[i] <- successes[i]+sum(x!=0)
            u <- update_emp_bandit(emp_bandit,num_samples,x,i)
            emp_bandit <- u$Bandit
            num_samples <- u$NumSamples
          }
        }
      }
    }
    l <- l+1
  }
  final_emp_means <- c()
  for(i in 1:K){
    final_emp_means <- c(final_emp_means,mean_reward(emp_bandit[[i]]))
  }
  correctness <- TRUE
  if(final_emp_means[1]<max(final_emp_means)){
    correctness <- FALSE
  }
  cat("Optimal weights for current empirical distribution at l = ",l,": \n",w_hat,"\n",
      "Stopping index = ",stopping_index,"Threshold = ",threshold,"\n")
  cat("\n Means of arms: \n",final_emp_means,"\n")
  cat("\n Number of samples drawn for each arm: \n", num_samples)
  cat("\n Overall sampling complexity: ", sum(num_samples))
  return(list("Complexity" = sum(num_samples),
              "Means" = final_emp_means,
              "Correctness"= correctness,
              "Successes"=successes))
}

#The following code will run 50 sample paths of TS(A) on our chosen bandit instance
set.seed(NULL)
seed_set <- sample(200:3000,10,replace=FALSE)

complexities <- c() #Vector to store sample complexities
correctnesses <- c() #Vector to store Boolean variables that state whether or not the correct best arm is recommended
successful_draws <- c() #Vector to store how many draws of each arm return nonzero rewards
tic() #Timer starts
for(pip in seed_set){
  solution <- ApproxAlgo(pip)
  complexities <- c(complexities,solution$Complexity)
  correctnesses <- c(correctnesses,solution$Correctness)
  successful_draws <- rbind(successful_draws,solution$Successes)
}
toc() #Timer ends
cat("Lower bound on expected sampling complexity at ",(delta*100),"% level of confidence = ",T,"\n")
cat("Average complexity = ",mean(complexities),"\n","Average correctness = ",mean(correctnesses),"\n")
#The following code can be un-commented to get output in the form of .csv spreadsheets
#data <- as.data.frame(cbind(complexities,correctnesses,successful_draws))
#write.csv(data,"Homogeneous output approx 5 arms.csv")