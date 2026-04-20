# empty workspace
rm(list=ls()); gc()

#load required packages
library(psyphy)
library(MASS)
library(dplyr)

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
ID <- 1:N
Class <- rep(1:ceiling(N_Classes), each = N_per_Class)[1:N]
Grade <- ifelse(Class < (round(N_Classes) / 2), 4, 5)
Age <- round(runif(N, 9, 10), 1)
Age <- ifelse(Grade == 5, Age + 1, Age)
Class <- rep(Class, each=n_item)
Grade <- rep(Grade, each=n_item)
Age <- rep(Age, each=n_item) 
Gender <- sample(c("M", "F"), size = N, replace = T) %>%
  rep(each=n_item)


# Child ability (theta) ~ N(0,1)
ability <- rnorm(N, mean = 0, sd = 1)
# Item difficulty (b) ~ N(0,1)
difficulty <- rnorm(n_item, mean = 0, sd = 1)


########################################
# TRAIT MEASURERES SIMULATION
########################################
sigma <- 0.5
var_names <- c("MA", "buoy", "avoid", "sc_trait", "worry", "sc_state")
n_vars <- length(var_names)
Sigma <- matrix(0.5,nrow=n_vars,ncol=n_vars) + diag(n_vars)*(0.5)
rownames(Sigma) <- colnames(Sigma) <- var_names

Sigma["MA", "buoy"] <- -0.5
Sigma["MA", "sc_trait"] <- -0.4
Sigma["MA", "sc_state"] <- -0.4
Sigma["MA", "worry"] <- 0.2

Sigma["buoy", "avoid"] <- -0.5
Sigma["buoy", "worry"] <- -0.5

Sigma["avoid", "sc_trait"] <- -0.5
Sigma["avoid", "sc_state"] <- -0.5

Sigma["sc_trait", "worry"] <- -0.4
Sigma["sc_trait", "sc_state"] <- 0.2

Sigma["worry", "sc_state"] <- -0.4

Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
rIntTraits <- MASS::mvrnorm(n=N, mu = rep(0, n_vars), Sigma = Sigma)

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


########################################
# function based on Rasch's model
########################################
prob_corr <- function(theta, b) {
  1 / (1 + exp(-(theta - b)))
}


# beta standardized
anx_z <- scale(math_anxiety)
avoid_z <- scale(avoidance)
buoy_z <- scale(buoyancy)

# beta parameters for 2nd session
beta_skip <- c(b0 = -0.5, avoid =  1.0, anx =  0.9, buoy = -0.8)
beta_hint <- c(b0 = -0.2, avoid = -0.4, anx =  0.9, buoy = -0.3)


########################################
# PREPARE DATAFRAME
########################################

tot_rows <- N * n_item
id <- rep(ID, each = n_item)

# matrix with N columns with the item in a random order --> vectorized
item <- replicate(N, sample(1:n_item, size = n_item, replace = F)) %>%
  as.vector()

# splitting the vector in 2 sessions
session <- rep(c(rep(1, n_item_session), rep(2, n_item_session)), N)

df <- data.frame(id, Gender, Age, Grade, Class, item, session, math_anxiety, 
                 buoyancy, avoidance, sc_trait, worry, sc_state, anx_z, avoid_z, buoy_z)

df$b <- difficulty[df$item]
df$theta <- ability[df$id]
df$p_correct <- prob_corr(df$theta, df$b)
df$accuracy <- rbinom(tot_rows, 1, df$p_correct)

# linear logit

i_s2 <- which(df$session == 2)
lp_skip <- beta_skip["b0"] + 
  beta_skip["anx"] * df$anx_z[i_s2] + 
  beta_skip["avoid"] * df$avoid_z[i_s2] + 
  beta_skip["buoy"] * df$buoy_z[i_s2]

lp_hint <- beta_hint["b0"] + 
  beta_hint["anx"] * df$anx_z[i_s2] + 
  beta_hint["avoid"] * df$avoid_z[i_s2] + 
  beta_hint["buoy"] * df$buoy_z[i_s2]

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

# setting accuracy = 0 if skip
df$accuracy[i_s2] <- ifelse(df$skip[i_s2] == 1, 0, df$accuracy[i_s2])
df[df$session == 1, c("skip", "help", "self")] <- NA

########################################
# OUTPUT
########################################

df <- df %>%
  select(id, Gender, Age, Grade, Class, item, session, accuracy,
         skip, help, self, math_anxiety, buoyancy, avoidance, sc_trait, 
         worry, sc_state) %>%
  arrange(id, session, item, Age)

ma_df <- df %>%
  filter(session == 1) %>%
  group_by(id) %>%
  summarise(math_ability_T0 = sum(accuracy))

df <- df %>%
  left_join(ma_df, by = "id")

write.csv(df, "SimulatedDatasetBehaviours.csv", row.names = F)


