library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(lmerTest)
library(stringi)
library(stringr)
library(dplyr)
library(purrr)
library(tidyverse)
library(scales)
library(optimx)
library(lme4)

options(scipen = 999)


#path<-'E:/MSUPE/piansrann/df_28/df_400_800'
path<-'C:/Users/User/Downloads/df_0_400'
path<-'C:/Users/User/Downloads/df_bl'


read_plus <- function(flnm) {
  read.csv(flnm) %>% 
    mutate(filename = flnm)
}

######### Load dataframes #########
path<-'C:/Users/User/Downloads/df_bl'

tbl_with_sources <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~fread(.))

#bl experiments
tbl_with_sources <- tbl_with_sources[tbl_with_sources$round!= "1", ]
tbl_with_sources <- tbl_with_sources [tbl_with_sources$round!= "2", ] # table with common bl for stim_type, label, round, subj like in original python version
tbl_with_sources_bl = tbl_with_sources %>%
  group_by(subject, label, stim_type, round) %>%
  summarise(theta_bl = mean(theta_power)) %>%
  ungroup() 
df_bl <-tbl_with_sources_bl
#tbl_with_sources$mean <- rowMeans(tbl_with_sources[, 2:40])
#tbl_with_sources$mean
colnames(tbl_with_sources)
tbl_with_sources_bl = tbl_with_sources %>%
  group_by(subject, label, stim_type, round) %>%
  summarise(theta_bl = mean(theta_power)) %>%
  ungroup() 

tbl_with_sources_bl$US_type <- ifelse(grepl("_", tbl_with_sources_bl$stim_type), sapply(strsplit(tbl_with_sources_bl$stim_type, "_"), `[`, 2), "CS_minus")
tbl_with_sources_bl$CS_type <- sub("_.*", "", tbl_with_sources_bl$stim_type)

df_bl <- tbl_with_sources_bl[tbl_with_sources_bl$round!= "1", ]
df_bl <- df_bl [df_bl$round!= "2", ] # table with common bl for stim_type, label, round, subj like in original python version

#tbl_with_sources$US_type <- ifelse(grepl("_", tbl_with_sources$stim_type), sapply(strsplit(tbl_with_sources$stim_type, "_"), `[`, 2), "CS_minus")
#tbl_with_sources$CS_type <- sub("_.*", "", tbl_with_sources$stim_type)


#df_bl <- tbl_with_sources[tbl_with_sources$round!= "1", ]
#df_bl <- df_bl [df_bl$round!= "2", ]

path<-'C:/Users/User/Downloads/df_0_400'

tbl_with_sources <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~fread(.))
tbl_with_sources$US_type <- ifelse(grepl("_", tbl_with_sources$stim_type), sapply(strsplit(tbl_with_sources$stim_type, "_"), `[`, 2), "CS_minus")
tbl_with_sources$CS_type <- sub("_.*", "", tbl_with_sources$stim_type)

df1 <- tbl_with_sources[tbl_with_sources$round!= "1", ]
df1 <- df1 [df1 $round!= "2", ]

merged_df <- merge(df_bl, df1, by = c('subject', 'label', 'stim_type', 'round')) ########exp bl

#test <- df1$theta_power - df_bl$theta_power
df1$test <- merged_df$theta_power - merged_df$theta_bl #bl correction 
#df1 <- df1 [df1 $round!= "run2", ]
######### if you need, remove outliers #######here i use column with bl correction 'test'
data_beta_sum = df1 %>%
  group_by(subject, label) %>%
  summarise(theta_mean = mean(test),
            theta_sd = sd(test)) %>%
  ungroup() %>%
  mutate(theta_high = theta_mean + (2 * theta_sd)) %>%
  mutate(theta_low = theta_mean - (2 * theta_sd))

data_accuracy_clean = df1 %>%
  inner_join(data_beta_sum) %>%
  filter(test < theta_high) %>%
  filter(test > theta_low)

######### if you need, remove outliers ####### i dont use this part
#data_beta_sum = df1 %>%
# group_by(subject, label) %>%
 # summarise(theta_mean = mean(theta_power),
#            theta_sd = sd(theta_power)) %>%
  #ungroup() %>%
 # mutate(theta_high = theta_mean + (2 * theta_sd)) %>%
#  mutate(theta_low = theta_mean - (2 * theta_sd))

#data_accuracy_clean = df1 %>%
#  inner_join(data_beta_sum) %>%
#  filter(theta_power < theta_high) %>%
#  filter(theta_power > theta_low)
###################################################################

elem<- filter(data_accuracy_clean, CS_type=='pict'|CS_type=='sound')
elem$stim_type<- 'element'
comp<- filter(data_accuracy_clean, CS_type=='comb')
comp$stim_type<- 'complex'
df2<- rbind(comp,elem)
colnames(df2)
colnames(df2)[11]<-'theta power'

labels<- unique(tbl_with_sources$label)
labels<- na.omit(labels)
p_vals<- data.table()
for (l in 1:length(labels)){
  temp<- subset(df2,label == labels[l])
  print(temp)
  m <-  lmer(theta_power ~ stim_type*US_type +(1|subject), data = temp,REML = FALSE)
  summary(m)
  an <- anova(m)
  #print(an)m
  an <- data.table(an,keep.rownames = TRUE)
  an_cols <- c('rn','Pr(>F)') 
  an <- an[, ..an_cols]
  an$`Pr(>F)` <- format(an$`Pr(>F)`, digits = 3)
  an$interval <- "theta_power"
  an$interval <- gsub('theta power','',an$interval)
  an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
  an$label <- unique(temp$label)
  p_vals <- rbind(p_vals,an)
}
p_vals <- subset(p_vals, !grepl("\\?\\?\\?-lh|\\?\\?\\?-rh", label)) # to remove ??? roi
setwd("E:/MSUPE/piansrann/theta_main")
write.csv(p_vals, "main_effects_0_400_ms_d2.csv")

p <- p_vals$`stim_type:US_type`
p_vals$fdr <- p.adjust(p, method = 'fdr', n = length(p))

library(permutes)
library(buildmer)
library(lmPerm)

labels<- unique(tbl_with_sources$label)
labels<- na.omit(labels)


labels<- unique(tbl_with_sources$label)
labels<- na.omit(labels)
p_vals<- data.table()
for (l in 1:length(labels)){
  temp<- subset(df2,label == labels[l])
  tm<- t(data.table(m))
  m <-  perm.lmer(theta_power ~ stim_type*US_type +(1|subject), data = temp, type='anova')
  print(m)
  m$label <- labels[l]
  p_vals <- rbind(p_vals,m)
}

temp<- subset(df2,label == labels[1])
tm<- t(data.table(m))
m <-  perm.lmer(theta_power ~ stim_type*US_type +(1|subject), data = temp, type='anova')
print(m)
m$label <- labels[1]

colnames(tm) <- tm[1, ]
tm <- tm[-1, ]
tm <- data.table(m)
p_vals <- subset(p_vals, !grepl("\\?\\?\\?-lh|\\?\\?\\?-rh", label)) # to remove ??? roi
setwd("E:/MSUPE/piansrann/theta_main")
write.csv(p_vals, "main_effects_only_comb_0_400_ms_d2.csv")





p_type <- filter(p_vals, Factor=='stim_type:US_type')
p_type$fdr <- p.adjust(p, method = 'fdr',length(p))



