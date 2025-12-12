
rm(list=ls()); gc()



simData = function(N = 300, hypothesis = 1){
  
  require(psyphy)
  require(MASS)
  require(dplyr)
  require(ordinal)
  
  X = data.frame(b_buoy=NA, se_buoy=NA, p_buoy=NA, 
                 b_avoid=NA, se_avoid=NA, p_avoid=NA,
                 b_worry=NA, se_worry=NA, p_worry=NA,
                 b_sc=NA, se_sc=NA, p_sc=NA)
  
  tryCatch({
    
    n_item <- 60            
    n_item_session <- 30   
    N_per_Class <- 20       
    N_Classes <- N / N_per_Class
    ID <- 1:N
    
    Class <- rep(1:ceiling(N_Classes), each = N_per_Class)[1:N]
    Grade <- ifelse(Class < (round(N_Classes) / 2), 4, 5)
    Age <- round(runif(N, 9, 10), 1)
    Age <- ifelse(Grade == 5, Age + 1, Age)
    
    Class <- rep(Class, each=n_item)
    Grade <- rep(Grade, each=n_item)
    Age <- rep(Age, each=n_item) 
    Gender <- sample(c("M", "F"), size = N, replace = T) %>% rep(each=n_item)
    
    ability <- rnorm(N, mean = 0, sd = 1)
    difficulty <- rnorm(n_item, mean = 0, sd = 1)
    
    
    ########################################
    # TRAIT MEASURERES SIMULATION
    ########################################
    sigma = 0.5
    var_names <- c("MA", "buoy", "avoid", "sc_trait", "worry", "sc_state")
    n_vars <- length(var_names)
    Sigma <- matrix(0.5,nrow=n_vars,ncol=n_vars) + diag(n_vars)*(0.5)
    rownames(Sigma) <- colnames(Sigma) <- var_names
    
    Sigma["MA", "buoy"] <- -0.5
    Sigma["MA", "sc_trait"] <- -0.4
    Sigma["MA", "sc_state"] <- -0.4
    
    Sigma["buoy", "avoid"] <- -0.5
    Sigma["buoy", "worry"] <- -0.5
    
    Sigma["avoid", "sc_trait"] <- -0.5
    Sigma["avoid", "sc_state"] <- -0.5
    
    Sigma["sc_trait", "worry"] <- -0.4
    
    Sigma["worry", "sc_state"] <- -0.4
    
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
    rIntTraits <- mvrnorm(n=N, mu = rep(0, n_vars), Sigma = Sigma)
    
    
    
    
    # MA
    math_anxiety <- rIntTraits[, "MA"] + rnorm(length(rIntTraits[, "MA"]), 0, sigma)
    math_anxiety <- round(plogis(math_anxiety) * (45-9)+9) %>%
      rep(each=n_item)
    
    # buoyancy
    buoyancy <- rIntTraits[, "buoy"] + rnorm(length(rIntTraits[, "buoy"]), 0, sigma)
    buoyancy <- round(plogis(buoyancy) * (28-4)+4) %>%
      rep(each=n_item)
    
    # avoidance
    avoidance <- rIntTraits[, "avoid"] + rnorm(length(rIntTraits[, "avoid"]), 0, sigma)
    avoidance <- round(plogis(avoidance) * (40-10)+10) %>%
      rep(each=n_item)
    
    # perceived competence (trait)
    sc <- rIntTraits[, "sc_trait"] + rnorm(length(rIntTraits[, "sc_trait"]), 0, sigma)
    sc_trait <- round(plogis(sc) * (40-10)+10) %>%
      rep(each=n_item)
    
    ########################################
    # STATE MEASURES SIMULATION
    ########################################
    
    # worry
    worry <- rIntTraits[, "worry"] + rnorm(length(rIntTraits[, "worry"]), 0, sigma)
    worry <- round(plogis(worry) * (12-3)+3) %>%
      rep(each=n_item)
    
    
    # perceived competence (state)
    sc_state <- rIntTraits[, "sc_state"] + rnorm(length(rIntTraits[, "sc_state"]), 0, sigma)
    sc_state <- round(plogis(sc_state) * (12-3)+3) %>%
      rep(each=n_item)
    
    
    if (hypothesis == 1) {
      # Buoyancy --> + self
      # Avoidance/Anxiety --> + skip
      beta_ord <- c(wor = -0.5, avoid = -0.6, buoy = 0.5)
    } else {
      beta_ord <- c(wor = 0, avoid = 0, buoy = 0) #H0
    }
    
    tot_rows <- N * n_item
    id <- rep(ID, each = n_item)
    item <- replicate(N, sample(1:n_item, size = n_item, replace = F)) %>% as.vector()
    session <- rep(c(rep(1, n_item_session), rep(2, n_item_session)), N)
      
    wor_z <- scale(worry)
    avoid_z <- scale(avoidance)
    buoy_z <- scale(buoyancy)
      
    df <- data.frame(id, Gender, Age, Grade, Class, item, session, 
                       math_anxiety, buoyancy, avoidance, sc_trait,
                       worry, sc_state,
                       wor_z, avoid_z, buoy_z)
      
    prob_corr <- function(theta, b) { 1 / (1 + exp(-(theta - b))) }
    df$b <- difficulty[df$item]
    df$theta <- ability[df$id]
    df$p_correct <- prob_corr(df$theta, df$b)
    df$accuracy <- rbinom(tot_rows, 1, df$p_correct)
      
    i_s2 <- which(df$session == 2)
    
    sdRand <- .5
    subjInt <- rnorm(ID, 0, sdRand)
    randEff <- subjInt[df$id[i_s2]]
    
    # generate linear component
    yLinear <- beta_ord["wor"] * df$wor_z[i_s2] + 
      beta_ord["avoid"] * df$avoid_z[i_s2] + 
      beta_ord["buoy"] * df$buoy_z[i_s2] + 
      randEff +
      rlogis(length(i_s2))
    
    # cut points to separate skip/help/self
    tau1 <- -1
    tau2 <- 1
    
    df$choice <- NA
    df$choice[i_s2] <- case_when(
      yLinear < tau1 ~ "skip",
      yLinear >= tau1 & yLinear < tau2 ~ "help",
      TRUE ~ "self"
    )
    
    df$accuracy[i_s2] <- ifelse(df$choice[i_s2] == "skip", 0, df$accuracy[i_s2])
    
    
    ma_df <- df %>% filter(session == 1) %>% group_by(id) %>% summarise(math_ability_T0 = sum(accuracy), .groups='drop')
    
    # filtering only the part of the dataset we need
    d <- df %>%
      left_join(ma_df, by = "id") %>%
      filter(session == 2)
    
    ########################################
    # ANALYSIS
    ########################################
    
    d$choice_c  <- factor(d$choice, levels = c("skip", "help", "self"))
    
    # standardizing predictors
    d$buoyancy <- scale(d$buoyancy)
    d$avoidance <- scale(d$avoidance)
    d$worry <- scale(d$worry)
    d$math_anxiety <- scale(d$math_anxiety)
    
    
    ########################################
    # RQ2 (Predictors on categorical model)
    ########################################
    d$choice_ord <- ordered(d$choice, levels = c("skip", "help", "self"))
    my_formula_ord <- choice_ord ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id)
    fit_ord = clmm(formula = my_formula_ord, 
                   data = d, 
                   link = "logit", 
                   control = clmm.control(maxIter = 200))
    
    summary <- summary(fit_ord)
    coefs <- summary$coefficients
    
    # store results
    X$b_buoy = coefs["buoyancy", "Estimate"]
    X$se_buoy = coefs["buoyancy", "Std. Error"]
    X$p_buoy = coefs["buoyancy", "Pr(>|z|)"]
    X$b_avoid <- coefs["avoidance", "Estimate"]
    X$se_avoid <- coefs["avoidance", "Std. Error"]
    X$p_avoid <- coefs["avoidance", "Pr(>|z|)"]
    X$b_worry <- coefs["worry", "Estimate"]
    X$se_worry <- coefs["worry", "Std. Error"]
    X$p_worry <- coefs["worry", "Pr(>|z|)"]
    X$b_sc <- coefs["sc_state", "Estimate"]
    X$se_sc <- coefs["sc_state", "Std. Error"]
    X$p_sc <- coefs["sc_state", "Pr(>|z|)"]
  }, error=function(e){})
  
  return(X)
}

simData(N=200, hypothesis = 0)

##############################################################
# H1

library(parallel)

num_cores = detectCores() - 1
cl = makeCluster(num_cores)

hypothesis = 1
N = 300
clusterExport(cl, varlist = c("simData", "N", "hypothesis"))
start_time <- Sys.time()
results = parLapply(cl, 1:1000, function(x) simData(N, hypothesis))
end_time <- Sys.time()
results
test_results <- (end_time - start_time)
test_results

stopCluster(cl)

save(results,file="results1000_H1.RData")

##############################################################
# H0

library(parallel)

num_cores = detectCores() - 1
cl = makeCluster(num_cores)

hypothesis = 0
N = 300
clusterExport(cl, varlist = c("simData", "N", "hypothesis"))
results = parLapply(cl, 1:1000, function(x) simData(N, hypothesis))
results

stopCluster(cl)

save(results,file="results1000_H0.RData")

##############################################################
# RUN POWER ANALYSIS AND CHECK PARAMETERS

library(dplyr)

load("results1000_H1.RData")
x = do.call(rbind,results)
lapply(x,median, na.rm=T)
lapply(x[, grep("b_|se_", names(x))], sd, na.rm=T)

p_unadj = x[,grep("p_",names(x))]
apply(p_unadj,2,function(x) mean(x<.05, na.rm=T))
mean(apply(p_unadj[,c("p_buoy","p_avoid","p_worry")],1,function(x) mean(x<.05)==1), na.rm=T)

p_adj = t(apply(p_unadj, 1, function(x) p.adjust(x, method = "fdr")))
apply(p_adj,2,function(x) mean(x<.05, na.rm=T))
mean(apply(p_adj[,c("p_buoy","p_avoid","p_worry")],1,function(x) mean(x<.05)==1), na.rm=T)







