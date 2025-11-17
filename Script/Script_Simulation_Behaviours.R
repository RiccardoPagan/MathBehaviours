# empty workspace
rm(list=ls())
gc()

#load required packages
library(psyphy)
library(MASS)
library(ggplot2)

# for reproducibility
set.seed(42)

########################################
# PARAMETERS
########################################
N <- 300                
n_item <- 60            
n_item_session <- 30   
N_per_Class <- 20       
N_Classes <- N / N_per_Class
Class <- rep(1:ceiling(N_Classes), each = N_per_Class)[1:N]
Grade <- ifelse(Class < (round(N_Classes) / 2), 4, 5)
Age <- round(runif(N, 9, 10), 1)
Age <- ifelse(Grade == 5, Age + 1, Age)
Gender <- sample(c("M", "F"), size = N, replace = T)

# Child ability (theta) ~ N(0,1)
ability <- rnorm(N, mean = 0, sd = 1)
# Item difficulty (b) ~ N(0,1)
difficulty <- rnorm(n_item, mean = 0, sd = 1)

# empty dataframe
df <- data.frame(
  ID = numeric(0),
  gender = character(0),
  age = numeric(0),
  grade = numeric(0),
  Class = numeric(0),
  item = numeric(0),
  session = numeric(0),
  accuracy = numeric(0),
  skip = numeric(0),
  help = numeric(0),
  self = numeric(0),
  math_anxiety = numeric(0),
  buoyancy = numeric(0),
  avoidance = numeric(0),
  pc = numeric(0),
  worry = numeric(0),
  metacognition= numeric(0)
)

########################################
# TRAIT MEASURERES SIMULATION
########################################

# MA
ma <- rnorm(N, 0, 1)
sigma_ma <- 0.5
math_anxiety <- ma + rnorm(N, 0, sigma_ma)
math_anxiety <- round(plogis(math_anxiety) * (45-9)+9)

# buoyancy
buo <- rnorm(N, 0, 1)
sigma_buo <- 0.5
buoyancy <- buo + rnorm(N, 0, sigma_buo)
buoyancy <- round(plogis(buoyancy) * (28-4)+4)

# avoidance
av <- rnorm(N, 0, 1)
sigma_av <- 0.5
avoidance <- av + rnorm(N, 0, sigma_av)
avoidance <- round(plogis(avoidance) * (40-10)+10)

# perceived competence
hart <- rnorm(N, 0, 1)
sigma_hart <- 0.5
pc <- hart + rnorm(N, 0, sigma_hart)
pc <- round(plogis(pc) * (40-10)+10)


########################################
# STATE MEASURES SIMULATION
########################################

# worry
wor <- rnorm(N, 0, 1)
sigma_wor <- 0.5
worry <- wor + rnorm(N, 0, sigma_wor)
worry <- round(plogis(worry) * (40-10)+10)

# metacognition
meta <- rnorm(N, 0, 1)
sigma_meta <- 0.5
metacognition <- meta + rnorm(N, 0, sigma_meta)
metacognition <- round(plogis(metacognition) * (4-1)+1)


########################################
# function based on Rasch's model
########################################
prob_corr <- function(theta, b) {
  1 / (1 + exp(-(theta - b)))
}

########################################
# SIMULATE BEHAVIOURS FOR 2ND SESSION
########################################

anx_z <- scale(math_anxiety)
avoid_z <- scale(avoidance)
buoy_z <- scale(buoyancy)

softmax <- function(mx) {
  mx <- sweep(mx, 1, apply(mx,1,max), "-")
  ex <- exp(mx)
  ex / rowSums(ex)
}

simulate_from_betas <- function(beta_skip, beta_hint, anx_z, avoid_z, buoy_z) {
  lp_skip <- beta_skip["b0"] + beta_skip["avoid"]*avoid_z + beta_skip["anx"]*anx_z + beta_skip["buoy"]*buoy_z
  lp_hint  <- beta_hint["b0"]  + beta_hint["avoid"]*avoid_z  + beta_hint["anx"]*anx_z  + beta_hint["buoy"]*buoy_z
  lp_self <- rep(0, length(lp_skip))
  M <- cbind(lp_skip, lp_hint, lp_self)
  probs <- softmax(M)
  P1 <- probs[,1]; P2 <- probs[,2]; P3 <- probs[,3]
  weighted_avoid <- P1*beta_skip["avoid"] + P2*beta_hint["avoid"]
  deriv_P3_avoid <- - P3 * weighted_avoid
  weighted_buoy  <- P1*beta_skip["buoy"]  + P2*beta_hint["buoy"]
  deriv_P3_buoy  <- - P3 * weighted_buoy
  list(probs = probs,
       deriv_P3_avoid = deriv_P3_avoid,
       deriv_P3_buoy = deriv_P3_buoy)
}

beta_skip <- c(b0 = -0.5, avoid =  1.0, anx =  0.9, buoy = -0.8)
beta_hint <- c(b0 = -0.2, avoid = -0.4, anx =  0.9, buoy = -0.3)
res <- simulate_from_betas(beta_skip, beta_hint, anx_z, avoid_z, buoy_z)

# check if P3 is working as we want
cat("Proporzione osservazioni dove P3 (self) diminuisce con avoidance (deriv < 0):\n")
print(mean(res$deriv_P3_avoid < 0))

cat("Proporzione osservazioni dove P3 (self) aumenta con buoyancy (deriv > 0):\n")
print(mean(res$deriv_P3_buoy > 0))

summary(res$deriv_P3_avoid)
summary(res$deriv_P3_buoy)

########################################
# MATH SIMULATION AND DATASET POPULATION
########################################
for (j in 1:N) {
  
  ID <- j
  gender <- Gender[j]
  age <- Age[j]
  grade <- Grade[j]
  Class_j <- Class[j]
  theta_j <- ability[j]
  
  anx_j <- math_anxiety[j]
  buoy_j <- buoyancy[j]
  av_j <- avoidance[j]
  pc_j <- pc[j]
  
  worry_j  <- worry[j]
  metacognition_T1_j <- metacognition[j]
  
  anx_j_z <- anx_z[j]
  buoy_j_z <- buoy_z[j]
  av_j_z <- avoid_z[j]
  
  # Session 1
  item_sess1 <- sample(1:n_item, size = n_item_session, replace = F)
  item_sess2 <- setdiff(1:n_item, item_sess1)
  
  for (i in item_sess1) {
    b_i <- difficulty[i]
    p <- prob_corr(theta_j, b_i)
    accuracy <- rbinom(1, 1, p)
    
    df <- rbind(df, data.frame(
      ID = ID, 
      gender = gender, 
      age = age, 
      grade = grade, 
      Class = Class_j,
      item = i, 
      session = 1, 
      accuracy = accuracy,
      skip = NA, 
      help = NA, 
      self = NA,
      math_anxiety = anx_j,
      buoyancy = buoy_j, 
      avoidance = av_j, 
      pc = pc_j,
      worry = worry_j, 
      metacognition = metacognition_T1_j
    ))
  }
  
  # Session 2
  for (i in item_sess2) {
    b_i <- difficulty[i]
    lp_skip <- beta_skip["b0"] + beta_skip["anx"]*anx_j_z + beta_skip["avoid"]*av_j_z + beta_skip["buoy"]*buoy_j_z
    lp_hint <- beta_hint["b0"] + beta_hint["anx"]*anx_j_z + beta_hint["avoid"]*av_j_z + beta_hint["buoy"]*buoy_j_z
    lp_self <- 0
    
    logits <- c(lp_skip, lp_hint, lp_self)
    probs <- exp(logits - max(logits))
    probs <- probs / sum(probs)
    
    behaviour <- sample(c("skip", "hint", "self"), 1, prob = probs)
    
    choice_skip <- as.numeric(behaviour == "skip")
    choice_hint <- as.numeric(behaviour == "hint")
    choice_self <- as.numeric(behaviour == "self")
    
    p <- prob_corr(theta_j, b_i)
    accuracy <- ifelse(behaviour == "skip", 0, rbinom(1, 1, p))
    
    df <- rbind(df, data.frame(
      ID = ID,
      gender = gender,
      age = age,
      grade = grade,
      Class = Class_j,
      item = i,
      session = 2,
      accuracy = accuracy,
      skip = choice_skip,
      help = choice_hint,
      self = choice_self,
      math_anxiety = anx_j,
      buoyancy = buoy_j, 
      avoidance = av_j, 
      pc = pc_j,
      worry = worry_j,
      metacognition = metacognition_T1_j
    ))
  }
}

########################################
# OUTPUT
########################################
df <- df[order(df$ID, df$session, df$item), ]
write.csv(df, "Data/DatasetSimulato.csv", row.names = F)





#cat("Simulazione terminata ed esportato in file 'DatasetSimulato.csv'")
table(df$gender)/n_item


# pr1 <- subset(df, session == 1)
# pr2 <- subset(df, session == 2)
# 
# cor(pr1$worry_pre, pr2$worry_pre)
# cor(pr1$worry_post, pr2$worry_post)
# summary(pr1$worry_post)
