
library(ggplot2)
library(ggsignif)
library(stringr)
library(tidyverse)
library(ggpubr)

# Exp9 <- read.csv("Z:\\Taylor\\Experiments\\009_DrugPathways\\CellProfiler_Exp9_Plate01\\cellData.csv")
Exp14 <- read_csv("cellData14.csv")
Exp15 <- read_csv("cellData15.csv")

se <- function(x) sqrt(var(x)/length(x))

theme_Bryan <- function () {
  theme_bw(base_size=14*96/72) %+replace%
    theme(
      #text = element_text(family = "sans"),
      plot.title = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text.x = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5,
      strip.text.y = element_text(size = 7*96/72, hjust = 0.5), #hjust = 0.5, #10 == 2.57mm
      strip.background = element_blank(),
      axis.title =element_text(size=6*96/72), #, face = "bold"
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(size = 5*96/72),
      legend.title=element_text(size=6*96/72),
      legend.text=element_text(size=5*96/72),
      axis.ticks.length =  unit(2.75, "pt"))
}

df <- data.frame(isolation=c(rep("1",30),rep("2",30)),
                 conditions=c(rep("Media",3),rep("Sarpogrelate",3),rep("TGFb",3),rep("PE",3),
                              rep("TGFb + Escitalopram",3),rep("TGFb + Sarpogrelate 2uM",3),rep("TGFb + Sarpogrelate 5uM",3),
                              rep("PE + Escitalopram",3),rep("PE + Sarpogrelate 2uM",3),rep("PE + Sarpogrelate 5uM",3),
                              rep("Media",3),rep("Sarpogrelate",3),rep("TGFb",3),rep("PE",3),
                              rep("TGFb + Escitalopram",3),rep("TGFb + Sarpogrelate 2uM",3),rep("TGFb + Sarpogrelate 5uM",3),
                              rep("PE + Escitalopram",3),rep("PE + Sarpogrelate 2uM",3),rep("PE + Sarpogrelate 5uM",3)),
                 medians=c(Exp14$Media,Exp14$Sar4,Exp14$TGFb,Exp14$PE,Exp14$TGFb_Esc,Exp14$TGFb_Sar3,Exp14$TGFb_Sar4,Exp14$PE_Esc,Exp14$PE_Sar3,Exp14$PE_Sar4,
                           Exp15$Media,Exp15$Sar4,Exp15$TGFb,Exp15$PE,Exp15$TGFb_Esc,Exp15$TGFb_Sar3,Exp15$TGFb_Sar4,Exp15$PE_Esc,Exp15$PE_Sar3,Exp15$PE_Sar4))

# Take the means of the group median values
df_means_1 <- df %>% 
  filter(isolation == "1") %>% 
  group_by(conditions) %>% 
  summarise_at(vars(medians), list(means=mean, sem=se)) %>% 
  mutate(isolation=c(rep("1",10)))
df_means_2 <- df %>% 
  filter(isolation == "2") %>% 
  group_by(conditions) %>% 
  summarise_at(vars(medians), list(means=mean, sem=se)) %>% 
  mutate(isolation=c(rep("2",10)))
df_means <- rbind(df_means_1,df_means_2)

# Separate by condition
df_tgfb <- df_means %>% 
  filter(!grepl('PE',conditions))
df_pe <- df_means %>% 
  filter(!grepl('TGFb',conditions))

# PE data
df_pe %>% 
  ggplot(aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", color='black') +
  geom_errorbar(aes(ymin=means-sem, ymax=means+sem), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Media","Sarpogrelate","PE","PE + Escitalopram","PE + Sarpogrelate 2uM","PE + Sarpogrelate 5uM"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() +
  theme_Bryan()

# TGFb data
df_tgfb %>% 
  ggplot(aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", color='black') +
  geom_errorbar(aes(ymin=means-sem, ymax=means+sem), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Media","Sarpogrelate","TGFb","TGFb + Escitalopram","TGFb + Sarpogrelate 2uM","TGFb + Sarpogrelate 5uM"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() +
  theme_Bryan()


## Statistics test

library(emmeans)

# Separate by condition
Traw <- df %>% 
  filter(!grepl('PE',conditions))
Praw <- df %>% 
  filter(!grepl('TGFb',conditions))

# releveling conditions to reference T_1
Traw <- Traw %>% 
  mutate(conditions=factor(conditions,levels=c("TGFb",
                                               "TGFb + Escitalopram",
                                               "TGFb + Sarpogrelate 2uM",
                                               "TGFb + Sarpogrelate 5uM",
                                               "Media",
                                               "Sarpogrelate")))

Praw <- Praw %>% 
  mutate(conditions=factor(conditions,levels=c("PE",
                                               "PE + Escitalopram",
                                               "PE + Sarpogrelate 2uM",
                                               "PE + Sarpogrelate 5uM",
                                               "Media",
                                               "Sarpogrelate")))

modelT <- lm(medians~isolation+conditions,data=Traw)
summary(modelT)  #these are RAW t-tests with RAW p-values

modelP <- lm(medians~isolation+conditions,data=Praw)
summary(modelP)  #these are RAW t-tests with RAW p-values

# post-hoc
Tmeans <- emmeans(modelT,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Tmeans$contrasts

Pmeans <- emmeans(modelP,specs= trt.vs.ctrl ~ conditions) # dunnet's, compare treatments to control
Pmeans$contrasts
