setwd("/home/user5/R/WYY/research2")
require(nimble)
require(coda)
library(doParallel)
library(nimbleHMC)
library(ggplot2)
library(MCMCvis)
library(zoo)
library(reshape2)
library(dplyr)
library(rlang)
library(gridExtra)
load("dengue.RData") 

n_v_ext <- 11.6667*N_h_ext
i_v_ext <- 40*y_ext
s_v_ext <- n_v_ext-i_v_ext
s_h_ext <- N_h_ext-y_ext
i_h_ext <- y_ext
r_h_ext <- data.frame(matrix(0, ncol = nR, nrow = nT ))
gamma_h <- 0.7903
nu_i <- rep(0.1, nR)
n_v_int <- 11.6667*N_h_int
y_int[1, 8] <- 1
y_int[1, 9] <- 1
i_v_int <- 40*y_int
s_v_int <- n_v_int-i_v_int
s_h_int <- N_h_int-y_int
i_h_int <- y_int
r_h_int <- data.frame(matrix(0, ncol = nR, nrow = nT ))

# constants
nimbleConsts <- list(nR, nT, 
                     as.matrix(N_h_ext), as.matrix(n_v_ext),
                     as.matrix(s_h_ext), as.matrix(i_h_ext), as.matrix(r_h_ext),
                     as.matrix(s_v_ext), as.matrix(i_v_ext), 
                     as.matrix(N_h_int), as.matrix(n_v_int),
                     as.matrix(s_h_int), as.matrix(i_h_int), as.matrix(r_h_int),
                     as.matrix(s_v_int), as.matrix(i_v_int), 
                     as.matrix(mu_v), #as.numeric(rho), as.numeric(mu_h), 
                     as.matrix(a), as.matrix(c), gamma_h, as.matrix(T_ij_ex), as.matrix(T_ij_in), as.numeric(nu_i))
names(nimbleConsts) <- c("nR", "nT", 
                         "N_h_ext", "N_v_ext", "s_h_ext", "i_h_ext", "r_h_ext", 
                         "s_v_ext", "i_v_ext", 
                         "N_h_int", "N_v_int", "s_h_int", "i_h_int", "r_h_int", 
                         "s_v_int", "i_v_int", 
                         "mu_v", #"rho", "mu_h",
                         "a", "c", "gamma_h", "T_ij_ex", "T_ij_in", "nu_i") 

# nimble data
nimbleData <- list(as.matrix(y_ext), as.matrix(y_int))
names(nimbleData) <- c("y_ext", "y_int")


nbModel <- nimbleCode({
  # 综合更新
  # 外部更新
  for (i in 1:nR) {
    S_v_ext[1,i] <- s_v_ext[1,i]
    I_v_ext[1,i] <- i_v_ext[1,i]
    S_h_ext[1,i] <- s_h_ext[1,i]
    I_h_ext[1,i] <- i_h_ext[1,i]
    R_h_ext[1,i] <- r_h_ext[1,i]
    lambda_h_ext[1,i] <- i_h_ext[1,i]
    
    # p_ext[1,i] <- theta_ext[i]/(theta_ext[i]+lambda_h_ext[1,i])
    # y_ext[1,i] ~ dnegbin(p_ext[1,i], theta_ext[i])
    # predy_ext[1,i] ~ dnegbin(p_ext[1,i], theta_ext[i])
    y_ext[1,i] ~ dpois(lambda_h_ext[1,i])
    predy_ext[1,i] ~ dpois(lambda_h_ext[1,i])
    
    S_v_int[1,i] <- s_v_int[1,i]
    I_v_int[1,i] <- i_v_int[1,i]
    S_h_int[1,i] <- s_h_int[1,i]
    I_h_int[1,i] <- i_h_int[1,i]
    R_h_int[1,i] <- r_h_int[1,i]
    lambda_h_int[1,i] <- i_h_int[1,i]
    
    # p_int[1,i] <- theta_int[i]/(theta_int[i]+lambda_h_int[1,i])
    # y_int[1,i] ~ dnegbin(p_int[1,i], theta_int[i])
    # predy_int[1,i] ~ dnegbin(p_int[1,i], theta_int[i])
    y_int[1,i] ~ dpois(lambda_h_int[1,i])
    predy_int[1,i] ~ dpois(lambda_h_int[1,i])
    
    # 考虑外部流入
    # ext_S_flow[1] <- sum(S_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6])
    # ext_I_flow[1] <- sum(I_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6])
    # ext_R_flow[1] <- sum(R_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6])
    ext_S[1, i] <- (sum(S_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6]))/N_h_int[1, i]
    ext_I[1, i] <- (sum(I_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6]))/N_h_int[1, i]
    ext_R[1, i] <- (sum(R_h_ext[1,1:10] * T_ij_ex[1:10,6] / N_h_int[1, 6]))/N_h_int[1, i]
    
    
    S_v_ext[2,i] <- S_v_ext[1,i] + mu_v[1,i] * N_v_ext[1,i] - a[1,i] * c[1,i] * S_v_ext[1,i] * I_h_ext[1,i] / N_h_ext[1,i] - mu_v[1,i] * S_v_ext[1,i]
    I_v_ext[2,i] <- I_v_ext[1,i] + a[1,i] * c[1,i] * S_v_ext[1,i] * I_h_ext[1,i] / N_h_ext[1,i] - mu_v[1,i] * I_v_ext[1,i]
    
    lambda_h_ext[2,i] <- beta_h_ext[1,i] * S_h_ext[1,i] * I_v_ext[1,i] / N_h_ext[1,i]
    # p_ext[2,i] <- theta_ext[i] / (theta_ext[i] + lambda_h_ext[2,i])
    # y_ext[2,i] ~ dnegbin(p_ext[2,i], theta_ext[i])# rnbinom(1, size = theta_ext[i], prob = p[2,i])
    # predy_ext[2,i] ~ dnegbin(p_ext[2,i], theta_ext[i])
    y_ext[2,i] ~ dpois(lambda_h_ext[2,i])
    predy_ext[2,i] ~ dpois(lambda_h_ext[2,i])
    
    S_h_ext[2,i] <- S_h_ext[1,i] - y_ext[2,i] - nu_i[i] * S_h_ext[1,i] + sum(S_h_ext[1,i] / N_h_ext[1,i] * T_ij_ex[,i])
    I_h_ext[2,i] <- I_h_ext[1,i] + y_ext[2,i] - gamma_h * I_h_ext[1,i] - nu_i[i] * I_h_ext[1,i] + sum(I_h_ext[1,i] / N_h_ext[1,i] * T_ij_ex[,i]) 
    R_h_ext[2,i] <- R_h_ext[1,i] + gamma_h * I_h_ext[1,i] - nu_i[i] * R_h_ext[1,i] + sum(R_h_ext[1,i] / N_h_ext[1,i] * T_ij_ex[,i])
    
    # 内部更新 (第六个地区)
    S_v_int[2,i] <- S_v_int[1,i] + mu_v[1,6] * N_v_int[1,i] - a[1,6] * c[1,6] * S_v_int[1,i] * I_h_int[1,i] / N_h_int[1,i] - mu_v[1,6] * S_v_int[1,i]
    I_v_int[2,i] <- I_v_int[1,i] + a[1,6] * c[1,6] * S_v_int[1,i] * I_h_int[1,i] / N_h_int[1,i] - mu_v[1,6] * I_v_int[1,i]
    
    lambda_h_int[2,i] <- beta_h_int[1,i] * S_h_int[1,i] * I_v_int[1,i] / N_h_int[1,i]
    # p_int[2,i] <- theta_int[i] / (theta_int[i] + lambda_h_int[2,i])
    # y_int[2,i] ~ dnegbin(p_int[2,i], theta_int[i])
    # predy_int[2,i] ~ dnegbin(p_int[2,i], theta_int[i])
    y_int[2,i] ~ dpois(lambda_h_int[2,i])
    predy_int[2,i] ~ dpois(lambda_h_int[2,i])
    
    # 考虑外部流入
    # ext_S_flow[2] <- sum(S_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6])
    # ext_I_flow[2] <- sum(I_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6])
    # ext_R_flow[2] <- sum(R_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6])
    
    ext_S[2, i] <- (sum(S_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6]))/N_h_int[2, i]
    ext_I[2, i] <- (sum(I_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6]))/N_h_int[2, i]
    ext_R[2, i] <- (sum(R_h_ext[2,1:10] * T_ij_ex[1:10,6] / N_h_int[2, 6]))/N_h_int[2, i]
    
    S_h_int[2,i] <- S_h_int[1,i] - y_int[2,i] - nu_i[i] * S_h_int[1,i] + sum(S_h_int[1,i] / N_h_int[1,i] * T_ij_in[,i]) + ext_S[1, i]
    I_h_int[2,i] <- I_h_int[1,i] + y_int[2,i] - gamma_h * I_h_int[1,i] - nu_i[i] * I_h_int[1,i] + sum(I_h_int[1,i] / N_h_int[1,i] * T_ij_in[,i]) + ext_I[1, i]
    R_h_int[2,i] <- R_h_int[1,i] + gamma_h * I_h_int[1,i] - nu_i[i] * R_h_int[1,i] + sum(R_h_int[1,i] / N_h_int[1,i] * T_ij_in[,i]) + ext_R[1, i]
  }
  
  for (t in 3:nT) {
    # 外部更新
    for (i in 1:nR) {
      S_v_ext[t,i] <- S_v_ext[t-1,i] + mu_v[t-1,i] * N_v_ext[t-1,i] - a[t-1,i] * c[t-1,i] * S_v_ext[t-1,i] * I_h_ext[t-1,i] / N_h_ext[t-1,i] - mu_v[t-1,i] * S_v_ext[t-1,i]
      I_v_ext[t,i] <- I_v_ext[t-1,i] + a[t-1,i] * c[t-1,i] * S_v_ext[t-1,i] * I_h_ext[t-1,i] / N_h_ext[t-1,i] - mu_v[t-1,i] * I_v_ext[t-1,i]
      
      lambda_h_ext[t,i] <- beta_h_ext[t-1,i] * S_h_ext[t-1,i] * I_v_ext[t-1,i] / N_h_ext[t-1,i]
      # p_ext[t,i] <- theta_ext[i] / (theta_ext[i] + lambda_h_ext[t,i])
      # y_ext[t,i] ~ dnegbin(p_ext[t,i], theta_ext[i])
      # predy_ext[t,i] ~ dnegbin(p_ext[t,i], theta_ext[i])
      y_ext[t,i] ~ dpois(lambda_h_ext[t,i])
      predy_ext[t,i] ~ dpois(lambda_h_ext[t,i])
      
      S_h_ext[t,i] <- S_h_ext[t-1,i] - y_ext[t,i] - nu_i[i] * S_h_ext[t-1,i] + sum(S_h_ext[t-1,i] / N_h_ext[t-1,i] * T_ij_ex[,i]) - sum(S_h_ext[t-2,i] / N_h_ext[t-2,i] * T_ij_ex[,i]) + nu_i[i] * S_h_ext[t-2,i]
      I_h_ext[t,i] <- I_h_ext[t-1,i] + y_ext[t,i] - gamma_h * I_h_ext[t-1,i] - nu_i[i] * I_h_ext[t-1,i] + sum(I_h_ext[t-1,i] / N_h_ext[t-1,i] * T_ij_ex[,i]) - sum(I_h_ext[t-2,i] / N_h_ext[t-2,i] * T_ij_ex[,i]) + nu_i[i] * I_h_ext[t-2,i]
      R_h_ext[t,i] <- R_h_ext[t-1,i] + gamma_h * I_h_ext[t-1,i] - nu_i[i] * R_h_ext[t-1,i] + sum(R_h_ext[t-1,i] / N_h_ext[t-1,i] * T_ij_ex[,i]) - sum(R_h_ext[t-2,i] / N_h_ext[t-2,i] * T_ij_ex[,i]) + nu_i[i] * R_h_ext[t-2,i]
    }
    
    # 内部更新 (第六个地区)
    for (i in 1:nR) {
      S_v_int[t,i] <- S_v_int[t-1,i] + mu_v[t-1,6] * N_v_int[t-1,i] - a[t-1,6] * c[t-1,6] * S_v_int[t-1,i] * I_h_int[t-1,i] / N_h_int[t-1,i] - mu_v[t-1,6] * S_v_int[t-1,i]
      I_v_int[t,i] <- I_v_int[t-1,i] + a[t-1,6] * c[t-1,6] * S_v_int[t-1,i] * I_h_int[t-1,i] / N_h_int[t-1,i] - mu_v[t-1,6] * I_v_int[t-1,i]
      
      lambda_h_int[t,i] <- beta_h_int[t-1,i] * S_h_int[t-1,i] * I_v_int[t-1,i] / N_h_int[t-1,i]
      # p_int[t,i] <- theta_int[i] / (theta_int[i] + lambda_h_int[t,i])
      # y_int[t,i] ~ dnegbin(p_int[t,i], theta_int[i])
      # predy_int[t,i] ~ dnegbin(p_int[t,i], theta_int[i])
      y_int[t,i] ~ dpois(lambda_h_int[t,i])
      predy_int[t,i] ~ dpois(lambda_h_int[t,i])
      
      # 考虑外部流入
      # ext_S_flow[t] <- sum(S_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6])
      # ext_I_flow[t] <- sum(I_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6])
      # ext_R_flow[t] <- sum(R_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6])
      ext_S[t, i] <- (sum(S_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6]))/N_h_int[t, i]
      ext_I[t, i] <- (sum(I_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6]))/N_h_int[t, i]
      ext_R[t, i] <- (sum(R_h_ext[t,1:10] * T_ij_ex[1:10,6] / N_h_int[t, 6]))/N_h_int[t, i]
      
      
      S_h_int[t,i] <- S_h_int[t-1,i] - y_int[t,i] - nu_i[i] * S_h_int[t-1,i] + sum(S_h_int[t-1,i] / N_h_int[t-1,i] * T_ij_in[,i]) - sum(S_h_int[t-2,i] / N_h_int[t-2,i] * T_ij_in[,i]) + nu_i[i] * S_h_int[t-2,i] + ext_S[t-1, i] - ext_S[t-2, i]
      I_h_int[t,i] <- I_h_int[t-1,i] + y_int[t,i] - gamma_h * I_h_int[t-1,i] - nu_i[i] * I_h_int[t-1,i] + sum(I_h_int[t-1,i] / N_h_int[t-1,i] * T_ij_in[,i]) - sum(I_h_int[t-2,i] / N_h_int[t-2,i] * T_ij_in[,i]) + nu_i[i] * I_h_int[t-2,i] + ext_I[t-1, i] - ext_I[t-2, i]
      R_h_int[t,i] <- R_h_int[t-1,i] + gamma_h * I_h_int[t-1,i] - nu_i[i] * R_h_int[t-1,i] + sum(R_h_int[t-1,i] / N_h_int[t-1,i] * T_ij_in[,i]) - sum(R_h_int[t-2,i] / N_h_int[t-2,i] * T_ij_in[,i]) + nu_i[i] * R_h_int[t-2,i] + ext_R[t-1, i] - ext_R[t-2, i]
    }
  }
  
  
  for (t in 1:nT-1){
    for (i in 1:nR){
      # beta_h_ext[t, i] <- alpha_ext*exp(-((t-u_ext[i])^2)/(2*sig_ext[i]^2))
      # beta_h_int[t, i] <- alpha_ext*exp(-((t-u_int[i])^2)/(2*sig_int[i]^2))
      # beta_h_ext[t, i] ~ dnorm(0, sd=10)
      # beta_h_int[t, i] ~ dnorm(0, sd=10)
      beta_h_ext[t, i] <- exp(b_ext1+b_ext2[i]*t+b_ext3[i]*t^2+b_ext4[i]*t^3)#+b_ext4[i]*t^4
      beta_h_int[t, i] <- exp(b_int1+b_int2[i]*t+b_int3[i]*t^2+b_int4[i]*t^3)#+b_int4[i]*t^4
      # beta_h_ext[t, i] <- exp(alpha_ext+b_ext[i]*cos(2*3.14*t/52) + c_ext[i]*sin(2*3.14*t/52) + phi_ext[t, i])
      # beta_h_int[t, i] <- exp(alpha_int+b_int[i]*cos(2*3.14*t/52) + c_int[i]*sin(2*3.14*t/52) + phi_int[t, i])
      # beta_h_ext[t, i] <- exp(alpha_ext+b_ext[i]*cos(2*3.14*t/52) + c_ext[i]*sin(2*3.14*t/52)) #+ phi_ext[t, i] #alpha_ext+b_ext[i]*cos(2*3.14*(t-phi_ext_26[i])/52) + c_ext[i]*cos(2*3.14*(t-phi_ext_52[i])/52)
      # beta_h_int[t, i] <- exp(alpha_int+b_int[i]*cos(2*3.14*t/52) + c_int[i]*sin(2*3.14*t/52)) #+ phi_int[t, i] #alpha_int+b_int[i]*cos(2*3.14*(t-phi_int_26[i])/52) + c_int[i]*cos(2*3.14*(t-phi_int_52[i])/52)
      # phi_ext[t, i] ~ dnorm(-3, sd=0.001)#dnorm(0.0001, sd=0.001)
      # phi_int[t, i] ~ dnorm(-3, sd=0.001)#dnorm(0.0001, sd=0.001)
    }
  }
  # theta_ext ~ dinvgamma(0.4, 0.3)
  # theta_int ~ dinvgamma(0.4, 0.3)
  b_ext1 ~ dnorm(0.03, sd=10) #dinvgamma(2,0.5)#dnorm(0, 0.001) # 均值为0，方差较大（非信息性先验）
  b_int1 ~ dnorm(0.03, sd=10) #dinvgamma(2,0.5)#dnorm(0, 0.001) # 均值为0，方差较大（非信息性先验）
  for (i in 1:nR) {
    # theta_ext[i] ~ dinvgamma(0.4, 0.3)
    # theta_int[i] ~ dinvgamma(0.4, 0.3)
    # u_ext[i] ~ dnorm(0, sd=10) # 均值为0，方差较大（非信息性先验）
    # sig_ext[i] ~ dnorm(0, sd=10) # 均值为0，方差较大（非信息性先验）
    # u_int[i] ~ dnorm(0, sd=10) # 均值为0，方差较大（非信息性先验）
    # sig_int[i] ~ dnorm(0, sd=10) # 均值为0，方差较大（非信息性先验）
    #b_ext1[i] ~ dnorm(0, sd=10)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_ext2[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    #b_int1[i] ~ dnorm(0, sd=10)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_int2[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_ext3[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_ext4[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_int3[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    b_int4[i] ~ dnorm(0, sd=1)#dinvgamma(2,0.5) # 均值为0，方差较大（非信息性先验）
    # phi_ext_26[i] ~ dunif(0, 26)
    # phi_int_26[i] ~ dunif(0, 26)
    # phi_ext_52[i] ~ dunif(0, 52)
    # phi_int_52[i] ~ dunif(0, 52)
  }
  # maxphi_ext <- max(phi_ext[1:nT,1:nR])
  # constrain_max_phi_ext ~ dconstraint(maxphi_ext <= 0.0003)
  # # minphi_ext <- min(phi_ext[1:nT,1:nR])
  # # constrain_min_phi_ext ~ dconstraint(minphi_ext <= -8)
  # maxphi_int <- max(phi_int[1:nT,1:nR])
  # constrain_max_phi_int ~ dconstraint(maxphi_int <= 0.0003)
  # # minphi_int <- min(phi_int[1:nT,1:nR])
  # # constrain_min_phi_int ~ dconstraint(minphi_int <= -8)
})

inits <- list(#alpha_ext=0.05,alpha_int=0.05,
  # beta_h_int=structure(.Data=rep(0.1, nR*(nT-1)),.Dim=c(nT-1, nR)),
  # beta_h_ext=structure(.Data=rep(0.005, nR*(nT-1)),.Dim=c(nT-1, nR)),
  #theta_ext=2,theta_int=2,#theta_ext=rep(2,10),theta_int=rep(2,10),
  # constrain_max_phi_ext=1,constrain_max_phi_int=1,
  # constrain_min_phi_ext=1,constrain_min_phi_int=1,
  # phi_ext=structure(.Data=rep(0, nR*(nT-1)),.Dim=c(nT-1, nR)),
  # phi_int=structure(.Data=rep(0, nR*(nT-1)),.Dim=c(nT-1, nR)),
  # u_ext=rep(10,10),u_int=rep(13,10),sig_ext=rep(15,10),sig_int=c(rep(20,8),80,20),
  b_ext1=-2,b_int1=-2,
  b_ext2=rep(0,10),b_int2=rep(0,10),
  b_ext3=rep(0,10),b_int3=rep(0,10),
  b_ext4=rep(0,10),b_int4=rep(0,10),
  #phi_ext_26=rep(10,10),phi_ext_52=rep(10,10),phi_int_26=rep(10,10),phi_int_52=rep(10,10),
  #beta_h_int=structure(.Data=rep(0.01, nR*nT),.Dim=c(nT, nR)),theta_int[i]=2,
  lambda_h_int=as.matrix(i_h_int),
  S_h_int=as.matrix(s_h_int),I_h_int=as.matrix(i_h_int),R_h_int=as.matrix(r_h_int),
  I_v_int=as.matrix(i_v_int),S_v_int=as.matrix(s_v_int),
  #beta_h_ext=structure(.Data=rep(0.005, nR*nT),.Dim=c(nT, nR)),theta_ext[i]=2,
  predy_ext=structure(.Data=rep(0, nR*nT),.Dim=c(nT, nR)),
  predy_int=structure(.Data=rep(0, nR*nT),.Dim=c(nT, nR)),
  lambda_h_ext=as.matrix(i_h_ext),
  S_h_ext=as.matrix(s_h_ext),I_h_ext=as.matrix(i_h_ext),R_h_ext=as.matrix(r_h_ext),
  I_v_ext=as.matrix(i_v_ext),S_v_ext=as.matrix(s_v_ext))
flow <- nimbleModel(code=nbModel, constants=nimbleConsts, data=nimbleData, 
                    inits = inits, buildDerivs = TRUE)
flow$initializeInfo()
flow$b_ext2
flow$beta_h_int
flow$beta_h_ext
flow$lambda_h_ext
flow$lambda_h_int
flow$y_ext
flow$beta_h_int

nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)#"phi_ext","phi_int",
nimbleOptions(getDependenciesIncludesPredictiveNodes=TRUE)# "theta","deta","beta_h","pred_y","alpha","u","sig""pred_y","S_h_int","I_h","R_h""u_ext", "sig_ext", "u_int", "sig_int",
MCMCconfig <- configureMCMC(flow, monitors=c("beta_h_ext","beta_h_int","b_ext1", "b_ext2", "b_ext3", "b_int1", "b_int2", "b_int3",  #"phi_ext_26", "phi_ext_52", "phi_int_26", "phi_int_52", #"theta_ext[i]","theta_int[i]","S_h_ext","I_h_ext","R_h_ext","S_v_ext","I_v_ext"
                                             "predy_ext","predy_int"),useConjugacy = FALSE, enableWAIC = TRUE)#"S_h_int","I_h_int","R_h_int","S_v_int","I_v_int",


# Compile

# Compile
nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(flow)#,resetFunctions = TRUE)
CnbMCMC <- compileNimble(nbMCMC, project = flow,resetFunctions = TRUE)


seed <- c(1111)
start_time <- Sys.time()
ful_4 <- runMCMC(CnbMCMC,inits=inits,setSeed=c(seed,seed*2),nchains = 2, nburnin=100000, niter = 130000, samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC = TRUE, thin=20)
end_time <- Sys.time()
end_time-start_time
save(ful_4, file = "/home/user5/R/WYY/research2/model_4_polynomial.RData")

# summ <- full_samples$summary$all.chains
# summ1 <- data.frame(summ)
# 
results  <- as.mcmc.list(ful_4$samples)
# 
# mcsum <- MCMCsummary(object = results, round = 6)
# View(mcsum)
# pdf("tracec.pdf", width = 8, height = 6)
# MCMCtrace(object = results,
#           pdf = TRUE, # no export to PDF
#           ind = TRUE,
#           filename = "tracec.pdf") # separate density lines per chain
# # 
plot(results, ask=T)

# plot(results[,"predy_ext[1, 1]"], ask=T)
# plot(results[,"beta_h_int[10, 7]"], ask=T)
# plot(results[,"I_h_int[10, 7]"], ask=T)
