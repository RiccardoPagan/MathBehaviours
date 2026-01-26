rm(list = ls()); gc()

#########################################################################
simData = function(N=300, target_coef = "muself_buoyancy", hypothesis = 1){
  
  require(psyphy)
  require(MASS)
  require(dplyr)
  require(brms) 
  
  X = data.frame(estimated_B = NA, SE = NA)
  
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
    
    beta_ord <- c(wor = -0.25, avoid = -0.30, buoy = 0.30, sc = 0.25)
    
    if (hypothesis == 1) {
      beta_skip <- c(b0 = 0.13, avoid =  0.30, wor =  0.25, buoy = -0.3, sc = -0.25)
      beta_hint <- c(b0 = -0.34, avoid = 0.15, wor =  0.10, buoy = -0.15, sc = -0.10)
    } else {
      beta_skip <- c(b0 = 0.13, avoid =  0, anx =  0, buoy = 0)
      beta_hint <- c(b0 = -0.34, avoid =  0, anx =  0, buoy = 0)
    }
    
    tot_rows <- N * n_item
    id<- rep(ID, each = n_item)
    item <- replicate(N, sample(1:n_item, size = n_item, replace = F)) %>% as.vector()
    session <- rep(c(rep(1, n_item_session), rep(2, n_item_session)), N)
    
    wor_z <- scale(worry)
    avoid_z <- scale(avoidance)
    buoy_z <- scale(buoyancy)
    sc_z <- scale(sc_state)
    
    df <- data.frame(id, Gender, Age, Grade, Class, item, session, 
                     math_anxiety, buoyancy, avoidance, sc_trait,
                     worry, sc_state,
                     wor_z, avoid_z, buoy_z, sc_z)
    
    prob_corr <- function(theta, b) { 1 / (1 + exp(-(theta - b))) }
    df$b <- difficulty[df$item]
    df$theta <- ability[df$id]
    df$p_correct <- prob_corr(df$theta, df$b)
    df$accuracy <- rbinom(tot_rows, 1, df$p_correct)
    
    i_s2 <- which(df$session == 2)
    
    TauSubj <- 1
    subjInt <- rnorm(N, 0, TauSubj)
    rSubj   <- subjInt[df$id[i_s2]]
    
    tauItem <- .5
    ItemInt <- rnorm(n_item, 0, tauItem)
    rItem <- ItemInt[df$item[i_s2]]
    
    lp_skip <- beta_skip["b0"] + beta_skip["wor"] * df$wor_z[i_s2] + beta_skip["avoid"] * df$avoid_z[i_s2] + beta_skip["buoy"] * df$buoy_z[i_s2] + beta_skip["sc"] * df$sc_z[i_s2] + rSubj + rItem
    lp_hint <- beta_hint["b0"] + beta_hint["wor"] * df$wor_z[i_s2] + beta_hint["avoid"] * df$avoid_z[i_s2] + beta_hint["buoy"] * df$buoy_z[i_s2] + beta_hint["sc"] * df$sc_z[i_s2] + rSubj + rItem
    lp_self <- rep(0, length(i_s2))
    
    logits <- cbind(lp_skip, lp_hint, lp_self)
    logits <- logits - apply(logits, 1, max)
    probs <- exp(logits)
    probs <- probs / rowSums(probs)
    
    # sampling behaviour
    pick <- runif(length(i_s2))
    p_skip <- probs[,1]
    p_hint <- probs[,2]
    df$skip[i_s2] <- as.numeric(pick < p_skip)
    df$help[i_s2] <- as.numeric(pick >= p_skip & pick < (p_skip + p_hint))
    df$self[i_s2] <- as.numeric(pick >= (p_skip + p_hint))
    
    df$accuracy[i_s2] <- ifelse(df$skip[i_s2] == 1, 0, df$accuracy[i_s2])
    df[df$session == 1, c("skip", "help", "self")] <- NA
    
    ma_df <- df %>% filter(session == 1) %>% group_by(id) %>% summarise(math_ability_T0 = sum(accuracy), .groups='drop')
    
    # filtering only the part of the dataset we need
    d <- df %>%
      left_join(ma_df, by = "id") %>%
      filter(session == 2)
    
    ########################################
    # ANALYSIS
    ########################################
    d$choice = paste0(d$skip,d$help,d$self)
    d$choice = case_when(
      d$choice=="100" ~ "skip",
      d$choice=="010" ~ "help",
      d$choice=="001" ~ "self",
      TRUE ~ NA
    )
    
    d$choice  <- factor(d$choice, levels = c("skip", "help", "self"))
    
    # standardizing predictors
    d$buoyancy <- scale(d$buoyancy)
    d$avoidance <- scale(d$avoidance)
    d$worry <- scale(d$worry)
    d$math_anxiety <- scale(d$math_anxiety)
    
    # categorical model
    my_formula <- choice ~ buoyancy + avoidance + worry + sc_state + math_anxiety + sc_trait + math_ability_T0 + Gender + (1|id)
    fit_cat <- brm(
      formula = my_formula,
      data = d,
      family = categorical(link = "logit"), 
      cores = 7
    )
    
    # check if the power of the priors (by default) is weak
    # get_prior(
    #   choice ~ buoyancy + avoidance + worry + sc_state + math_anxiety + 
    #     sc_trait + math_ability_T0 + Gender + (1|id),
    #   data = d,
    #   family = categorical(link = "logit"),
    #   cores = 7
    # )
    # 
    # seePrior <- prior_summary(fit_cat)
    
   
    
    list(fit = fit_cat)  
    
  },error=function(e){
    message("Error in simData (N=", N, "): ", e$message)
    list(fit = NULL)
  })
}
#########################################################################

#########################################################################
N_grid <- seq(100, 500, by = 50)
n_reps <- 5


# main parameters
pars_interest <- c( "muself_buoyancy", "muself_sc_state", "muself_avoidance", 
                    "muhelp_buoyancy")

vcov_results <- list()

for (i in seq_along(N_grid)) {
  N <- N_grid[i]
  message("Simulating N = ", N, " (", i, "/", length(N_grid), ")")
  
  # array to save the matrix for this N (4×4×n_reps)
  V_array <- array(NA, dim = c(length(pars_interest), length(pars_interest), n_reps))
  dimnames(V_array) <- list(pars_interest, pars_interest, paste0("rep", 1:n_reps))
  
  for (rep in 1:n_reps) {
    message("  Rep ", rep, "/", n_reps)
    out <- simData(N = N, hypothesis = 1)
    
    if (!is.null(out$fit)) {
      V <- vcov(out$fit, pars = pars_interest)
      V_array[, , rep] <- V
    }
  }
  
  vcov_results[[i]] <- list(
    N = N,
    vcov_array = V_array,
    vcov_mean = apply(V_array, c(1, 2), mean, na.rm = TRUE),
    SE_mean = sqrt(diag(apply(V_array, c(1, 2), mean, na.rm = TRUE)))
  )
}


# saving the output
saveRDS(vcov_results, "vcov_matrices.rds")
message("Matrices saved in the file named vcov_matrices.rds")
#########################################################################


#########################################################################
# COMPARISON MATRICES FOR VARIANCES AND COVARIANCES
#########################################################################

# variance matrix for each N
variance_matrix <- matrix(NA, nrow = length(N_grid), ncol = length(pars_interest))
rownames(variance_matrix) <- paste0("N=", N_grid)
colnames(variance_matrix) <- pars_interest

for (i in seq_along(vcov_results)) {
  variance_matrix[i, ] <- diag(vcov_results[[i]]$vcov_mean)
}

print(round(variance_matrix, 6))

# covariance matrix for main pairs across N
cov_pairs <- list(
  "buoy_sc" = c("muself_buoyancy", "muself_sc_state"),
  "buoy_avoid" = c("muself_buoyancy", "muself_avoidance"),
  "sc_avoid" = c("muself_sc_state", "muself_avoidance")
)

covariance_matrix <- matrix(NA, nrow = length(N_grid), ncol = length(cov_pairs))
rownames(covariance_matrix) <- paste0("N=", N_grid)
colnames(covariance_matrix) <- names(cov_pairs)

for (i in seq_along(vcov_results)) {
  V <- vcov_results[[i]]$vcov_mean
  for (j in seq_along(cov_pairs)) {
    pair <- cov_pairs[[j]]
    covariance_matrix[i, j] <- V[pair[1], pair[2]]
  }
}

print(round(covariance_matrix, 6))

# correlation matrix for the same pairs
correlation_matrix <- matrix(NA, nrow = length(N_grid), ncol = length(cov_pairs))
rownames(correlation_matrix) <- paste0("N=", N_grid)
colnames(correlation_matrix) <- names(cov_pairs)

for (i in seq_along(vcov_results)) {
  V <- vcov_results[[i]]$vcov_mean
  cor_mat <- cov2cor(V)
  for (j in seq_along(cov_pairs)) {
    pair <- cov_pairs[[j]]
    correlation_matrix[i, j] <- cor_mat[pair[1], pair[2]]
  }
}

print(round(correlation_matrix, 4))


# final output 
comparison_results <- list(
  N_grid = N_grid,
  variance_matrix = variance_matrix,
  covariance_matrix = covariance_matrix,
  correlation_matrix = correlation_matrix,
  vcov_full = vcov_results  
)

saveRDS(comparison_results, "vcov_comparison.rds")
message("Comparison matrices saved in vcov_comparison.rds")









