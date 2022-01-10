library(tidyverse)
library(readxl)
library(reshape2)
library(gridExtra)
library(ggbeeswarm)
library(ggpubr)
library(EnvStats)  
library(multcompView)
library(data.table)
library(multcomp) # package to add compact letter display for significance

setwd("~/Documents/University/Manuscripts/Brochet_SCFAs/Datasets/Gut_Metabolomics")

# CFUs data-analysis ----

mydata_CFUs <- read.csv("CFUs_exp1.csv")

mydata_CFUs$group <- as.factor(mydata_CFUs$group)
mydata_CFUs$group <- factor(mydata_CFUs$group, levels = (unique(mydata_CFUs$group)))

colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000",
  "MD" = "black")

#CFUs_res <- compare_means(value ~ species,  data = mydata_CFUs, method="wilcox.test")
pdf("CFUs.pdf")

ggplot(data=mydata_CFUs, aes(x=group, y=log10(value), color=group))+
  scale_color_manual(values = colors) +
  geom_boxplot()+
  #\stat_compare_means(aes(label = ..p.signif..),comparisons=my_comparisons, label.y = c(9,9.5,9)) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
  theme(axis.title.x = element_blank()) +
  #scale_y_log10()+
  ylab("CFUs/rectum") +
  ylim(c(0,10)) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

dev.off()

# Stats on CFUs

# Test for normality data-distribution 

mydata_CFUs <- read.csv("CFUs_exp1.csv")
mydata_CFUs <- subset(mydata_CFUs, species!="MD")
mydata_CFUs$value <- log10(mydata_CFUs$value)
qqnorm(mydata_CFUs$value)

# Perform ANOVA

mydata_CFUs$group <- as.factor(mydata_CFUs$group)
anova_CFUs <- aov(value~group, mydata_CFUs)
library(multcomp) # package to add compact letter display for significance
TukeyHSD_CFUs <- glht(anova_CFUs, linfct=mcp(group="Tukey"))
cld(TukeyHSD_CFUs)

# Weight data analysis ----

mydata_weight <- read.csv("rectum_weight_exp1.csv")

mydata_weight$Sample <- as.factor(mydata_weight$Sample)
mydata_weight$Sample <- factor(mydata_weight$Sample, levels = (unique(mydata_weight$Sample)))

colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000",
  "MD" = "black")

weight_res <- compare_means(Rectum_weight ~ species,  data = mydata_weight, method="wilcox.test")
pdf("weight.pdf")

ggplot(data=mydata_weight, aes(x=Sample, y=log10(Rectum_weight), color=Sample))+
  scale_color_manual(values = colors) +
  geom_boxplot()+
  #stat_compare_means(aes(label = ..p.signif..),comparisons=my_comparisons, label.y = c(9,9.5,9)) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
  theme(axis.title.x = element_blank()) +
  #scale_y_log10()+
  ylab("mg/rectum") +
  #ylim(c(0,10)) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

dev.off()

mydata_weight$Sample <- as.factor(mydata_weight$Sample)
anova_weight <- aov(Rectum_weight~Sample, mydata_weight)
library(multcomp) # package to add compact letter display for significance
TukeyHSD_weights <- glht(anova_weight, linfct=mcp(Sample="Tukey"))
cld(TukeyHSD_weights)


# Stats on weight

# Test for normality data-distribution 

mydata_weight <- subset(mydata_weight, species!="MD")
qqnorm(mydata_weight$Rectum_weight)

# Perform ANOVA

mydata_weight$Sample <- as.factor(mydata_weight$Sample)
anova_weight <- aov(Rectum_weight~Sample, mydata_weight)
TukeyHSD_weight <- glht(anova_weight, linfct=mcp(Sample="Tukey"))
cld(TukeyHSD_weight)


# Read GC-MS raw-data ----

Responses = read.csv(file = "gut_met_new.csv")
names(Responses)[names(Responses) == "Sample"] <- "Acq_time"
names(Responses)[names(Responses) == "X"] <- "Sample"
Responses= Responses[-1, ]

# Visualise ISTD responses across dataset

ggplot(Responses, aes(Acq_time, as.numeric(isovalerate.Results )))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Responses= filter(Responses, as.numeric(isovalerate.Results) >= 100) ### Based on the plot remove failed samples with no response for ISTD

ggplot(Responses, aes(Acq_time, as.numeric(isovalerate.Results )))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Reformat table (transpose)

data_long = Responses %>% pivot_longer(!c("Sample", "Acq_time"), names_to = "Metabolite", values_to = "Response") # turn data to
data_long$Sample <- str_sub(data_long$Sample, end = -3) # Remove .D from samples
data_long$Metabolite=gsub(".Results", "", data_long$Metabolite) # Remove ".Result" suffix from data_file column
data_long$Response= as.numeric(data_long$Response)

# Normalisation with ISTD 

Isovalerate_istd = subset(data_long, Metabolite =="isovalerate") # create df with the values for ISTD
Isovalerate_med = median(Isovalerate_istd$Response)

setDT(data_long) ## Turn df into dt
setkey(data_long, Sample)
setDT(Isovalerate_istd)
setkey(Isovalerate_istd, Sample)

data_long= 
  data_long[Isovalerate_istd, on = c("Sample", "Acq_time"),  Norm_resp := Response/i.Response*Isovalerate_med][]

`%notlike%` =negate (`%like%`)

# Create dataframe without Standard curves (data_long_samples)

data_long_samples=data_long[data_long$Sample %notlike% "STD", ]

write.csv(data_long_samples, "ISTD_normalised_responses.csv")

# Need to normalise them by weight  as well
# Sample 4_10 no weight measurement = exclude

data_long_samples <- subset(data_long_samples, Sample!="4_10")
rectum_weight <- read.csv("rectum_weight.csv")

write.csv(data_long_samples, "ISTD_normalised_responses.csv")

data_long_samples <- data_long_samples %>%
  add_column(Rectum_weight = c(1:(length(data_long_samples$Norm_resp))))

data_long_samples$Rectum_weight <- as.numeric(data_long_samples$Rectum_weight)

for (i in 1:nrow(data_long_samples)){    
  for (j in 1:nrow(rectum_weight)){ 
    if(data_long_samples[i,2]==(rectum_weight[j,1])){  # if variable match
      data_long_samples[i,6]<-rectum_weight[j,2]
      print(i)
      break
    }
    
  }
}

data_long_samples$weight_norm <- data_long_samples$Norm_resp * data_long_samples$Rectum_weight

write.csv(data_long_samples, "ISTD_weight_normalised_responses.csv")
data_long_samples <- read.csv("ISTD_weight_normalised_responses.csv")

# Plot normalised responses vs. non normalised responses

ggplot(data_long_samples, aes(Response, weight_norm))+
  geom_point()+
  stat_cor() +
  stat_regline_equation(label.x = 3, label.y = 32)+
  geom_smooth(method="lm")

data_long_samples <- data_long_samples[data_long_samples$Metabolite != "isovalerate", ] 

# Z-score conversion ----

# need to substract from all samples the mean of all samples and divide by the standard deviation

data_long_samples$z_score <- (data_long_samples$weight_norm - mean(data_long_samples$weight_norm)) / (sd(data_long_samples$weight_norm))


# PCA ----

group_file <- read.table("group_file.txt", header = TRUE, check.names = FALSE)
data_long_samples1 <- data_long_samples %>%
  add_column(group = NA)
data_long_samples1$group <- as.factor(data_long_samples1$group)
data_long_samples1 <- data_long_samples1 %>%
  add_column(SDP = NA)
data_long_samples1$SDP <- as.factor(data_long_samples1$SDP)
#group_file$group <- as.factor(group_file$group)


for (i in 1:nrow(data_long_samples1)){    
  for (j in 1:nrow(group_file)){ 
    if(data_long_samples1[i,2]==(group_file[j,1])){  # if variable match
      data_long_samples1[i,9]<-group_file[j,2]
      data_long_samples1[i,10]<-group_file[j,3]
      print(i)
      break
    }
    
  }
}



library(dplyr)
PCA_data <- subset(data_long_samples1, select=c("group", "Metabolite","z_score", "SDP"))

library(tidyverse)

PCA_data1 <- PCA_data %>% pivot_wider(names_from = Metabolite, 
            values_from = z_score, 
            values_fn = list(z_score=list))


PCA_data1 <- unnest(PCA_data1, cols = c(formate, acetate, propanoate, butyrate, lactic.acid, 
                                         succinate))

write.table(PCA_data1,"PCA_data1.txt")
PCA_data1 <- read.table("PCA_data1.txt", check.names = FALSE, header = TRUE)

pca <- prcomp(PCA_data1[c(-1,-2)],
              center = TRUE,
             scale. = TRUE)

colour <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000",
  "MD" = "black")

colors1 <- c(
  "SDP1" = "#FF87AB",
  "SDP4" = "#12BBEF",
  "SDP2" = "#FFD072",
  "SDP3" = "#3CB694",
  "MD" = "black")


size1 <- c(
  "SDP1" = 5,
  "SDP4" = 5,
  "SDP2" = 5,
  "SDP3" = 5,
  "MD" = 5)



library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
PCA_data1$SDP <- as.factor (PCA_data1$SDP)

pdf("rectum_PCA2.pdf")

ggbiplot(pca,
              obs.scale = 1,
              var.scale = 1,
              groups = PCA_data1$SDP,
              ellipse = TRUE,
              circle = FALSE,
         var.axes=TRUE)+
#geom_line(aes(size = 2)) +
scale_color_manual(values = colors1) +
theme(legend.position = 'none') +
geom_point(aes(fill=PCA_data1$group), shape=21, size=5)+
scale_fill_manual(values = colour) 
  
dev.off()

# Boxplots ----

mydatarel <-  data_long_samples1 %>% 
  extract(Sample, into = c("variable", "replicate"), "(.*)_([^_]+)$")
mydatarel$group <- as.factor(mydatarel$group)
mydatarel$Metabolite <- as.factor(mydatarel$Metabolite)


colors <- c(
  "ESL0185" = "#FF87AB",
  "ESL0186" = "#12BBEF",
  "ESL0183" = "#FFD072",
  "ESL0184" = "#3CB694",
  "ESL0353" = "#C9184A",
  "wkb8" = "#E6830B",
  "ESL0354" = "#E6A100",
  "ESL0263" = "#F14670",
  "ESL0351" = "#064789",
  "ESL0350" = "#007200",
  "ESL0261" = "#0089AF",
  "ESL0260" = "#38B000",
  "MD" = "black")


pdf("boxplots_rel3.pdf")

# Iterate over samples (=variable)

for (metaboliteName in levels(mydatarel$Metabolite)) {
  
  # Get sample specific data
  
  metaboliteData <- mydatarel %>% 
    filter(Metabolite == metaboliteName) %>% 
    mutate(plot_strain= as.factor(1:n()))
  
  p <- metaboliteData %>%
    mutate(group = fct_relevel(group, 
                              "ESL0185", "ESL0263", "ESL0353", 
                              "ESL0183", "wkb8", "ESL0354", 
                              "ESL0184", "ESL0260", "ESL0350", "ESL0186", "ESL0261", "ESL0351")) %>%
    
    ggplot(metaboliteData, mapping = aes(x=group, y=z_score, fill=group))+
    scale_fill_manual(values = colors) +
    geom_boxplot()+
    ggtitle(metaboliteName) +
    stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "MD") +
    theme_bw()+
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
    theme(axis.title.x = element_blank()) +
    ylab("z-score")#+
    #scale_y_continuous(limits = c(0, 2.5e07))
  

  print(p)
  
}

dev.off()



