##### Code used to generate plots and statistics for phenotypic analyses
#Plotting and analyzing changes in pigmentation scores across latitudinal, seasonal, and experimental populations

##### Skyler Berardi; contact: berardis@sas.upenn.edu

#Packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(nlme)
library(readxl)
library(emmeans)
library(plotrix)


#####
##### East Coast Latitudinal Cline: Plot and Statistics

## Read dataset
pigmentation_latitude_table <- read_excel("Pigmentation_latitudinal_scoring.xlsx")

latitude_levels <- c("1","2","3","4","5","6")
pigmentation_latitude_table <- pigmentation_latitude_table %>% mutate(latitude_pop = factor(latitude_pop, levels = latitude_levels, ordered = TRUE))
pigmentation_latitude_table <- pigmentation_latitude_table %>% mutate_at('isofemale_line', as.factor)

view(pigmentation_latitude_table)


## Statistics: East Coast Latitudinal Cline
# Linear mixed effects model: pigmentation scores ~ latitude, random = ~1 | isofemale line
LM_latitude_1 <- lme(fixed = pigmentation_score ~ latitude_pop, random = ~1 | isofemale_line, data = pigmentation_latitude_table)
anova(LM_latitude_1)
summary(LM_latitude_1)
LM_latitude_1$coefficients
LM_latitude_1$residuals


## Extracting means and standard errors from model for plotting
LM_latitude_means_1 <- emmeans(LM_latitude_1, ~ latitude_pop)
latitude_plot_1 <- as.data.frame(LM_latitude_means_1)
latitude_plot_1$dataset <- c(rep("latitudinal",6))
view(latitude_plot_1)


## Plotting latitudinal pigmentation data
#jitter
jitter_cline <- position_jitter(width = 0.2, height = 0)
#raw data
ggplot(pigmentation_latitude_table, aes(x=latitude_pop, y=pigmentation_score)) +
  geom_point(size=0.05, position=jitter_cline, alpha=1/3) +
  #Adding LM slope line
  geom_smooth(data=pigmentation_latitude_table, aes(x=latitude_pop, y=pigmentation_score, group=dataset, color=dataset), method='lm', se=FALSE,  color="#823830") +
  #line plot
  geom_line(data=latitude_plot_1, aes(x=latitude_pop, y=emmean, group=dataset, color=dataset), linetype = 3, linewidth=.75, inherit.aes = FALSE) +
  xlab("Latitude (°N)") +
  ylab("Mean Pigmentation Score") +
  ggtitle("East Coast Latitudinal Cline")+
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("25.28","30.2","33.95","38.02","39.88","42.45"))+
  scale_color_manual(values=c("#823830")) +
  #Adding mean points
  geom_point(data=latitude_plot_1, aes(x=latitude_pop, y=emmean, group=dataset, color=dataset), color="#823830",shape=15,size=2, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=latitude_plot_1, aes(x=latitude_pop, ymin=emmean-SE, ymax=emmean+SE, group=dataset, color=dataset), linewidth=.45,color="#823830",width=0.05, inherit.aes = FALSE)


## Export plot
tiff("latitude_plot_pig", units="in", width=6, height=4, res=300)

# insert ggplot code
ggplot(pigmentation_latitude_table, aes(x=latitude_pop, y=pigmentation_score)) +
  geom_point(size=0.05, position=jitter_cline, alpha=1/3) +
  #Adding LM slope line
  geom_smooth(data=pigmentation_latitude_table, aes(x=latitude_pop, y=pigmentation_score, group=dataset, color=dataset), method='lm', se=FALSE,  color="#823830") +
  #line plot
  geom_line(data=latitude_plot_1, aes(x=latitude_pop, y=emmean, group=dataset, color=dataset), linetype = 3, linewidth=.75, inherit.aes = FALSE) +
  xlab("Latitude (°N)") +
  ylab("Mean Pigmentation Score") +
  ggtitle("East Coast Latitudinal Cline")+
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("25.28","30.2","33.95","38.02","39.88","42.45"))+
  scale_color_manual(values=c("#823830")) +
  #Adding mean points
  geom_point(data=latitude_plot_1, aes(x=latitude_pop, y=emmean, group=dataset, color=dataset), color="#823830",shape=15,size=2, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=latitude_plot_1, aes(x=latitude_pop, ymin=emmean-SE, ymax=emmean+SE, group=dataset, color=dataset), linewidth=.45,color="#823830",width=0.05, inherit.aes = FALSE)

dev.off()



#####
##### Seasonal pigmentation patterns in wild populations: Plot and Statistics
#Populations collected from 2010 - 2015 in Linvilla Orchards, Media, Pennsylvania

## Read dataset
pigmentation_seasonal_table <- read_excel("Pigmentation_seasonal_scoring.xlsx")

seasonal_year_levels <- c("2010","2011","2012","2013","2014","2015")
seasonal_season_levels <- c("1_spring","2_fall")

pigmentation_seasonal_table <- pigmentation_seasonal_table %>% mutate(year = factor(year, levels = seasonal_year_levels, ordered = TRUE))
pigmentation_seasonal_table <- pigmentation_seasonal_table %>% mutate(season = factor(season, levels = seasonal_season_levels, ordered = TRUE))
pigmentation_seasonal_table <- pigmentation_seasonal_table %>% mutate_at('timepoint_ID', as.factor)
pigmentation_seasonal_table <- pigmentation_seasonal_table %>% mutate_at('isofemale_line', as.factor)
pigmentation_seasonal_table$year_season <- paste(pigmentation_seasonal_table$year, "-", pigmentation_seasonal_table$season)

view(pigmentation_seasonal_table)

#Calculating mean pigmentation score for each isofemale line
seasonal_isofemale_means <- pigmentation_seasonal_table %>%
  group_by(isofemale_line) %>%
  mutate(isofemale_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(year, season, year_season, isofemale_line, isofemale_mean_pigscore)

view(seasonal_isofemale_means)


## Statistics: Wild Seasonal Populations (Media, PA)
# Linear mixed effects model: pigmentation scores ~ year*season, random = ~1 | isofemale line
LM_seasonal_1 <- lme(fixed = isofemale_mean_pigscore ~ year*season, random = ~1 | isofemale_line, data=seasonal_isofemale_means)

anova(LM_seasonal_1)
summary(LM_seasonal_1)
LM_seasonal_1$coefficients
LM_seasonal_1$residuals


## Extracting means and standard errors from model for plotting
LM_seasonal_means_1 <- emmeans(LM_seasonal_1, ~ year*season)
print(LM_seasonal_means_1)
seasonal_plot_1 <- as.data.frame(LM_seasonal_means_1)
seasonal_plot_1$year_season <- paste(seasonal_plot_1$year, "-", seasonal_plot_1$season)
view(seasonal_plot_1)


## Plotting seasonal pigmentation patterns in wild populations
ggplot(seasonal_plot_1, aes(x= year_season, y=emmean, group=1)) +
  geom_line(color="#823830")+
  geom_point(color="#823830")+
  theme_classic() +
  ggtitle("Seasonal Patterns (Media, PA)") +
  xlab("Timepoint")+
  ylab("Mean Pigmentation Score") +
  scale_x_discrete(labels=c("June\n2010","Nov.\n2010","June\n2011","Nov.\n2011","June\n2012","Nov.\n2012","June\n2013","Nov.\n2013","June\n2014","Nov.\n2014","June\n2015","Nov.\n2015")) +
  #Adding standard error
  geom_errorbar(aes(x=year_season, ymin=emmean-SE, ymax=emmean+SE), linewidth=.45,color="#823830",width=0.10) +
  #data points for individual females
  geom_point(data=seasonal_isofemale_means, aes(x=year_season, y=isofemale_mean_pigscore), alpha=1/4, size=0.05)


# Export plot
tiff("seasonal_plot_pig", units="in", width=6, height=4, res=300)

# insert ggplot code
ggplot(seasonal_plot_1, aes(x= year_season, y=emmean, group=1)) +
  geom_line(color="#823830")+
  geom_point(color="#823830")+
  theme_classic() +
  ggtitle("Seasonal Patterns (Media, PA)") +
  xlab("Timepoint")+
  ylab("Mean Pigmentation Score") +
  scale_x_discrete(labels=c("June\n2010","Nov.\n2010","June\n2011","Nov.\n2011","June\n2012","Nov.\n2012","June\n2013","Nov.\n2013","June\n2014","Nov.\n2014","June\n2015","Nov.\n2015")) +
  #Adding standard error
  geom_errorbar(aes(x=year_season, ymin=emmean-SE, ymax=emmean+SE), linewidth=.45,color="#823830",width=0.10) +
  #data points for individual females
  geom_point(data=seasonal_isofemale_means, aes(x=year_season, y=isofemale_mean_pigscore), alpha=1/4, size=0.05)

dev.off()



#####
##### Seasonal pigmentation patterns in experimental populations: Plot and Statistics
#Experiment run 2016 in Philadelphia, Pennsylvania (Schmidt Lab experimental orchard)
#Mesocosms were labled "S1 - S10", but flies from S2 were not available to score (N=9 mesocosms).

## Read dataset
pigmentation_experimental_table <- read_excel("Pigmentation_experimental_scoring.xlsx")

experimental_levels = c("0","1","2")
pigmentation_experimental_table <- pigmentation_experimental_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = experimental_levels, ordered = TRUE))
pigmentation_experimental_table <- pigmentation_experimental_table %>% mutate_at('mesocosm', as.factor)
pigmentation_experimental_table$timepointID_mesocosm <- paste(pigmentation_experimental_table$timepoint_ID, "-", pigmentation_experimental_table$mesocosm)

view(pigmentation_experimental_table)

#Make data frame with S mesocosms only (no Founder; for linear mixed effects model)

pigmentation_experimental_table_S <- pigmentation_experimental_table %>% filter(population != 'Founder')
view(pigmentation_experimental_table_S)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

pigmentation_experimental_table_founder <- pigmentation_experimental_table %>% filter(population == 'Founder')
view(pigmentation_experimental_table_founder)

founder_mean <- mean(pigmentation_experimental_table_founder$pigmentation_score)
print(founder_mean)

founder_SE <- std.error(pigmentation_experimental_table_founder$pigmentation_score)
print(founder_SE)


#Calculating mean pigmentation score for each mesocosm
experimental_mesocosm_means <- pigmentation_experimental_table_S %>%
  group_by(timepointID_mesocosm) %>%
  mutate(mesocosm_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_mesocosm, timepoint_ID, mesocosm, mesocosm_mean_pigscore)

view(experimental_mesocosm_means)


## Statistics: Experimental Orchard 2016
# Linear mixed effects model: pigmentation scores ~ timepoint, random = ~1 | mesocosm/timepoint
LM_experimental_1 <- lme(mesocosm_mean_pigscore ~ timepoint_ID, random=~1|mesocosm/timepoint_ID, data = experimental_mesocosm_means)
anova(LM_experimental_1)
summary(LM_experimental_1)
LM_experimental_1$coefficients
LM_experimental_1$residuals


## Extracting means and standard errors from model for plotting
LM_experimental_means_1 <- emmeans(LM_experimental_1, ~ timepoint_ID)
print(LM_experimental_means_1)
experimental_plot_1 <- as.data.frame(LM_experimental_means_1)

view(experimental_plot_1)


## Experimental Orchard 2016: plotted mean pigmentation scores across seasons
ggplot(experimental_mesocosm_means, aes(x=timepoint_ID, y=mesocosm_mean_pigscore, group=mesocosm, color=mesocosm)) +
  geom_line() +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Experimental Orchard: 2016") +
  theme_classic() +
  ylim(7,14.7) +
  scale_x_discrete(labels=c("August 5","October 5"))+
  scale_color_manual(values=c("#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F"))+
  theme(legend.position = "none") +
  #Adding average points
  geom_line(data=experimental_plot_1, aes(x=timepoint_ID, y=emmean, group=1), linewidth=2, color="#823830", inherit.aes = FALSE) +
  geom_point(data=experimental_plot_1, aes(x=timepoint_ID, y=emmean, group=1), size=2.5, color="#823830", inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=experimental_plot_1, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE), linewidth=.45, color="#823830", width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=0.6, y=founder_mean), color="#823830", size=2.5) +
  geom_errorbar(aes(x=0.6, ymin=founder_mean-founder_SE, ymax=founder_mean+founder_SE), linewidth=.35, color="#823830", width=0.05)


# Export plot
tiff("experimental_plot_pig", units="in", width=6, height=4, res=300)

# insert ggplot code
ggplot(experimental_mesocosm_means, aes(x=timepoint_ID, y=mesocosm_mean_pigscore, group=mesocosm, color=mesocosm)) +
  geom_line() +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Experimental Orchard: 2016") +
  theme_classic() +
  ylim(7,14.7) +
  scale_x_discrete(labels=c("August 5","October 5"))+
  scale_color_manual(values=c("#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F"))+
  theme(legend.position = "none") +
  #Adding average points
  geom_line(data=experimental_plot_1, aes(x=timepoint_ID, y=emmean, group=1), linewidth=2, color="#823830", inherit.aes = FALSE) +
  geom_point(data=experimental_plot_1, aes(x=timepoint_ID, y=emmean, group=1), size=2.5, color="#823830", inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=experimental_plot_1, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE), linewidth=.45, color="#823830", width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=0.6, y=founder_mean), color="#823830", size=2.5) +
  geom_errorbar(aes(x=0.6, ymin=founder_mean-founder_SE, ymax=founder_mean+founder_SE), linewidth=.35, color="#823830", width=0.05)

dev.off()


