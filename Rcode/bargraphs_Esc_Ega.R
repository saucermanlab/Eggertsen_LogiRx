
library(ggplot2)
library(ggsignif)
library(stringr)

Exp12 <- read.csv("cellData12.csv")
Exp13 <- read.csv("cellData13.csv")

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

## Ega comparison

# TGFb data
df <- data.frame(isolation=c(rep("1",4),rep("2",4)),
                 conditions=c("Escitalopram","TGFb","TGFb + Escitalopram","TGFb + Eganelisib",
                              "Escitalopram","TGFb","TGFb + Escitalopram","TGFb + Eganelisib"),
                 means=c(mean(Exp12$Esc), mean(Exp12$TGFb), mean(Exp12$TGFb_Esc), mean(Exp12$TGFb_Ega),  
                         mean(Exp13$Esc),mean(Exp13$TGFb), mean(Exp13$TGFb_Esc), mean(Exp13$TGFb_Ega)), 
                 sm=c(se(Exp12$Esc), se(Exp12$TGFb), se(Exp12$TGFb_Esc), se(Exp12$TGFb_Ega),
                      se(Exp13$Esc), se(Exp13$TGFb), se(Exp13$TGFb_Esc), se(Exp13$TGFb_Ega)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Escitalopram","TGFb","TGFb + Escitalopram","TGFb + Eganelisib"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()


# PE data
df <- data.frame(isolation=c(rep("1",4),rep("2",4)),
                 conditions=c("Escitalopram","PE","PE + Escitalopram","PE + Eganelisib",
                              "Escitalopram","PE","PE + Escitalopram","PE + Eganelisib"),
                 means=c(mean(Exp12$Esc), mean(Exp12$PE), mean(Exp12$PE_Esc), mean(Exp12$PE_Ega), 
                         mean(Exp13$Esc), mean(Exp13$PE), mean(Exp13$PE_Esc), mean(Exp13$PE_Ega)),
                 sm=c(se(Exp12$Esc), se(Exp12$PE), se(Exp12$PE_Esc), se(Exp12$PE_Ega), 
                      se(Exp13$Esc), se(Exp13$PE), se(Exp13$PE_Esc), se(Exp13$PE_Ega)))

ggplot(df, aes(x=conditions, y=means, fill=isolation)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=means-sm, ymax=means+sm), width=0.2, position=position_dodge(0.9)) +
  scale_x_discrete(limits = c("Escitalopram","PE","PE + Escitalopram","PE + Eganelisib"),
                   labels = function(x) str_wrap(x,width = 10)) +
  scale_fill_grey() + 
  theme_Bryan()

