rm(list=ls())

library(ggplot2)
library(ggsignif)
library(stringr)

Exp7PE <- read.csv("cellData7pe.csv")
Exp7TGFb <- read.csv("cellData7tgf.csv")
Exp8PE <- read.csv("cellData8pe.csv")
Exp8TGFb <- read.csv("cellData8tgf.csv")

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

# TGFb data
df <- data.frame(isolation=c(rep("1",7),rep("2",7)),
                     conditions=c("Control","TGFb","TGFb + Escitalopram","TGFb + Tirbanibulin","TGFb + OSI-930","TGFb + Rifabutin","TGFb + Mifepristone",
                                  "Control","TGFb","TGFb + Escitalopram","TGFb + Tirbanibulin","TGFb + OSI-930","TGFb + Rifabutin","TGFb + Mifepristone"),
                     means=c(mean(Exp7TGFb$Media), mean(Exp7TGFb$TGFb), mean(Exp7TGFb$Esc_24), mean(Exp7TGFb$Tir_40), mean(Exp7TGFb$OSI_160), mean(Exp7TGFb$Rif_80), mean(Exp7TGFb$Mif_6), 
                             mean(Exp8TGFb$Media), mean(Exp8TGFb$TGFb), mean(Exp8TGFb$Esc_24), mean(Exp8TGFb$Tir_40), mean(Exp8TGFb$OSI_160), mean(Exp8TGFb$Rif_80), mean(Exp8TGFb$Mif_6)),
                     sm=c(se(Exp7TGFb$Media), se(Exp7TGFb$TGFb), se(Exp7TGFb$Esc_24), se(Exp7TGFb$Tir_40), se(Exp7TGFb$OSI_160), se(Exp7TGFb$Rif_80), se(Exp7TGFb$Mif_6), 
                          se(Exp8TGFb$Media), se(Exp8TGFb$TGFb), se(Exp8TGFb$Esc_24), se(Exp8TGFb$Tir_40), se(Exp8TGFb$OSI_160), se(Exp8TGFb$Rif_80), se(Exp8TGFb$Mif_6)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Control","TGFb","TGFb + Escitalopram","TGFb + Tirbanibulin","TGFb + OSI-930","TGFb + Rifabutin","TGFb + Mifepristone"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

# PE data
df <- data.frame(isolation=c(rep("1",7),rep("2",7)),
                 conditions=c("Control","PE","PE + Escitalopram","PE + Tirbanibulin","PE + OSI-930","PE + Rifabutin","PE + Mifepristone",
                              "Control","PE","PE + Escitalopram","PE + Tirbanibulin","PE + OSI-930","PE + Rifabutin","PE + Mifepristone"),
                 means=c(mean(Exp7PE$Media), mean(Exp7PE$PE), mean(Exp7PE$Esc_24), mean(Exp7PE$Tir_40), mean(Exp7PE$OSI_160), mean(Exp7PE$Rif_80), mean(Exp7PE$Mif_6), 
                         mean(Exp8PE$Media), mean(Exp8PE$PE), mean(Exp8PE$Esc_24), mean(Exp8PE$Tir_40), mean(Exp8PE$OSI_160), mean(Exp8PE$Rif_80), mean(Exp8PE$Mif_6)),
                 sm=c(se(Exp7PE$Media), se(Exp7PE$PE), se(Exp7PE$Esc_24), se(Exp7PE$Tir_40), se(Exp7PE$OSI_160), se(Exp7PE$Rif_80), se(Exp7PE$Mif_6), 
                      se(Exp8PE$Media), se(Exp8PE$PE), se(Exp8PE$Esc_24), se(Exp8PE$Tir_40), se(Exp8PE$OSI_160), se(Exp8PE$Rif_80), se(Exp8PE$Mif_6)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Control","PE","PE + Escitalopram","PE + Tirbanibulin","PE + OSI-930","PE + Rifabutin","PE + Mifepristone"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()



