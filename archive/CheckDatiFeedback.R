library(readxl)
DATA_MERGED_FINAL_FEEDBACK <- read_excel("~/Downloads/DATA_MERGED_FINAL_FEEDBACK.xlsx")
View(DATA_MERGED_FINAL_FEEDBACK)
d = data.frame(DATA_MERGED_FINAL_FEEDBACK )
rm (DATA_MERGED_FINAL_FEEDBACK)


names(d)

d = subset(d, Gruppo == "CG")

d$trait_MA = rowSums(d[, grep("AMAS.block_1.", names(d), value = T)], na.rm = T)
d$state_sc = rowSums(d[, grep("Pre_perception.block_1.", names(d), value = T)], na.rm = T)

d$trait_sc = rowSums(d[, grep("Self_Comp.", names(d), value = T)], na.rm = T)
d$trait_sc = rowSums(d[, c("Self_Comp.question1_reversed", "Self_Comp.question2", "Self_Comp.question3_reversed",
                       "Self_Comp.question4", "Self_Comp.question5_reversed", "Self_Comp.question6", 
                       "Self_Comp.question7_reversed", "Self_Comp.question8", "Self_Comp.question9_reversed",
                       "Self_Comp.question10")], na.rm = T)

mtx = cor(d[, c("trait_MA","worry_pre", "worry_post", "state_sc", "trait_sc")], use.comple)
