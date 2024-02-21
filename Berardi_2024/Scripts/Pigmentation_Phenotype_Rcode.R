##### Code used to generate plots and statistics for phenotypic analyses
#Plotting and analyzing changes in pigmentation scores across latitudinal, seasonal, and experimental populations

##### Skyler Berardi; contact: skylerberardi@gmail.com

#Packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(nlme)
library(lme4)


#####
##### East Coast Latitudinal Cline: Plot and Statistics
#Entered pigmentation scores from Excel scoring sheet and assigned to latitudinal populations. 
#The mean pigmentation scores of 17-30 isofemale lines per latitudinal population were entered here.
#Each isofemale line score was the mean of 10 individual females that were scored per line.

#Enter pigmentation scores: 
HomeFL <- c(18,13.6,12.8,9.6,11.7,17.5,11,15.1,22.9,15.9,15,17,10.6,18,19.8,15,16.6,10.1,14.6,16.8,14.9,16.8,9.4,11.1,14.5,14.5,16.4,12.5,12.2)
JackFL <- c(15.5,18,20.9,14.1,24,14.1,20.56,20.4,10.8,11.4,21.5,22.3,20,15.7,11.3,18,14.8,12.1,7.2,10.1,23,19.5,23.1,16.9,17.6,14.5,11.9,23.8,11.6,17.3)
AtheGA <- c(14.25,16.75,19.65,12.85,22.75,12.85,19.31,19.15,9.55,10.15,20.25,21.05,18.75,14.45,10.05,16.75,13.55)
CharVA <- c(15.2,11.6,16.7,14.3,16.2,19.7,26.5,23.6,15.5,16,19.8,21,21.2,19.7,11.7,22.6,18.1,20.8,20.7,29.1,14.3,18.1,15.2,22.3,19.5,15,22.4,19.2,17.7)
LinvPA <- c(15.4,30.3,19.2,25.3,22.1,28.7,24.5,15,23.6,23.9,27.9,19.6,21.2,18.6,18.8,21.8,6.8,12.8,26,15.6,22.4,22.7,18.4,18.1,14.1,13.6,18.2,15.2,20.5)
LancMA <- c(11.8,22.3,19.6,30.1,21.6,22.6,24.9,18.6,29.4,16.8,23.6,26.1,21.2,19.6,21.4,21.9,27.5,20,21.6,26.1,22.1,21.5,16.5,25.9,24.2,21.2,22,21.8,14.8)

#Made a list of each latitudinal population and computed the average score per population from raw data for each population
East_Coast_Cline_line <- list(HomeFL,JackFL,AtheGA,CharVA,LinvPA,LancMA)
East_Coast_Cline_avg_line <- sapply(East_Coast_Cline_line, mean)

#Created a dataset for average pigmentation scores for each population: common garden temperature (25C for all populations), latitudinal population (numbered 1-6), average pigmentation scores
CommonGardenTemp_Cline_line <- c(rep("25°C" , 6) )
latitudes_Cline_line <- rep(c("1", "2", "3", "4", "5", "6") , 1)
scores_Cline_line <- East_Coast_Cline_avg_line
data_Cline_line <- data.frame(CommonGardenTemp_Cline_line,latitudes_Cline_line,scores_Cline_line)

#Created a dataset showing mean pigmentation scores for each isofemale line at each latitude: common garden temperature, latitudinal population (numbered 1-6), individual data points (pts) for the isofemale line pigmentation scores. 
#Will use this to plot raw data points for isofemale lines on the latitudinal graph and run a linear model. 
Temp_points_Cline <- c(rep("25°C" , 163) )
Lat_points_Cline <- c(rep("1",29),rep("2",30),rep("3",17),rep("4",29),rep("5",29),rep("6",29))
Pts_points_Cline <- c(HomeFL,JackFL,AtheGA,CharVA,LinvPA,LancMA)
Data_points_Cline <- data.frame(Temp_points_Cline,Lat_points_Cline,Pts_points_Cline)


#Plotting East Coast Latitudinal Cline: plotted average pigmentation scores, raw data points, and a regression line (pigmentation scores ~ latitude)
jitter_Cline <- position_jitter(width = 0.2, height = 0)
#raw data
ggplot(Data_points_Cline, aes(x=Lat_points_Cline, y=Pts_points_Cline)) +
  geom_point(size=0.05, position=jitter_Cline, alpha=1/3) + 
  #line plot
  geom_line(data=data_Cline_line, aes(x = latitudes_Cline_line, y = scores_Cline_line, group=CommonGardenTemp_Cline_line, color=CommonGardenTemp_Cline_line), linetype = 3, size=.75) +
  xlab("Latitude (°N)") +
  ylab("Average Pigmentation Score") +
  ggtitle("East Coast Latitudinal Cline")+
  theme_classic() +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("25.28","30.2","33.95","38.02","39.88","42.45"))+
  scale_color_manual(values=c("#823830")) +
  #Adding points
  geom_point(aes(x=1,y=14.61),colour="#823830",shape=15,size=2)+
  geom_point(aes(x=2,y=16.73),colour="#823830",shape=15,size=2)+
  geom_point(aes(x=3,y=16),colour="#823830",shape=15,size=2)+
  geom_point(aes(x=4,y=18.74),colour="#823830",shape=15,size=2)+
  geom_point(aes(x=5,y=20.01),colour="#823830",shape=15,size=2)+
  geom_point(aes(x=6,y=21.95),colour="#823830",shape=15,size=2)+
  #Adding error bars with standard error calculated in Excel
  geom_errorbar(aes(x=1, ymin = 14.0208413, ymax = 15.2136387), size=.15,color="#823830",width=0.05) +
  geom_errorbar(aes(x=2, ymin = 15.8818757, ymax = 17.5821243), size=.15,color="#823830",width=0.05) +
  geom_errorbar(aes(x=3, ymin = 15.0017836, ymax = 17.0111564), size=.15,color="#823830",width=0.05) +
  geom_errorbar(aes(x=4, ymin = 17.9918035, ymax = 19.5047565), size=.15,color="#823830",width=0.05) +
  geom_errorbar(aes(x=5, ymin = 19.0282925, ymax = 20.9923875), size=.15,color="#823830",width=0.05) +
  geom_errorbar(aes(x=6, ymin = 21.2040036, ymax = 22.7063364), size=.15,color="#823830",width=0.05) +
  #Adding LM slope line
  geom_smooth(data=data_Cline_line, aes(x = latitudes_Cline_line, y = scores_Cline_line, group=CommonGardenTemp_Cline_line, color=CommonGardenTemp_Cline_line), method='lm', se=FALSE,  color="#823830")


#Statistics: East Coast Latitudinal Cline
#Linear model: pigmentation scores ~ latitude
LM_Cline <- lm(formula = Pts_points_Cline ~ Lat_points_Cline, data = Data_points_Cline)
summary(LM_Cline)



#####
#####Seasonal pigmentation patterns in wild populations: Plot and Statistics
#Populations collected from 2010 - 2015 in Linvilla Orchards ("Lin"), Media, Pennsylvania
#Entered pigmentation scores from Excel scoring sheet and assigned to seasonal populations (Sp=spring/early season, Fa=fall/late season) in each year. 
#The pigmentation scores of 10-20 isofemale lines per timepoint/year were entered here.
#Each isofemale line score was the mean of 5-10 individual females that were scored per line.

#Enter pigmentation scores
Lin_Sp_2010 <- c(36.9,38.6,29.2,27.7,25.2,34.85714286,31.66666667,38.5,35.3,27.6)
Lin_Fa_2010 <- c(21.6,16.7,26,23.22222222,24.66666667,18.7,22,22.4,23.5,31.5)
Lin_Sp_2011 <- c(29.625,32.2,21.4,25.8,32.4,24.1,30.1,29.6,32.6,22.2)
Lin_Fa_2011 <- c(20.6,11.8,13.3,18.7,17.9,21.3,9.7,18.2,18.9,12.1)
Lin_Sp_2012 <- c(20.9,25.2,21.4,26.33333333,27.5,22.9,29.55555556,31.2,29.875,28.71428571)
Lin_Fa_2012 <- c(12.1,12.6,11.8,13.8,15.375,13,15,20.4,16.1,14.6)
Lin_Sp_2013 <- c(20.9,24.1,23.5,17.3,23.6,30.9,20.6,19.7,24.6,23.7)
Lin_Fa_2013 <- c(13.66666667,19.1,16.1,15.6,14.3,14,13.7,18.8,16.3,16.7)
Lin_Sp_2014 <- c(23.6,28.8,25.6,19.6,21.6,19.4,14.8,17.4,26.6,27.8,17,22.8,24.2,22,25,20,23.6,20.2,24.6,24.6)
Lin_Fa_2014 <- c(24.8,20.4,24.8,23.2,21.2,19.4,24.4,22.6,20.2,21.6,24.2,19.6,15.2,17.8,18.2,23,21.2,20,25,32.2)
Lin_Sp_2015 <- c(31.6,24.4,29.6,30.4,30.2,34.6,29.4,35,27,34.6,33.8,35.2,31.4,20.6,31,27.6,32.4,32.6,34.2,32)
Lin_Fa_2015 <- c(20.8,13.8,21.6,18.2,19.6,12.4,12.4,15,24,20.4,11.4,14.2,16.6,21.4,25.8,25.8,18.2,28.4,20.2,20.6)

#Made a list of each seasonal population and computed average score for each population
Lin_line <- list(Lin_Sp_2010,Lin_Fa_2010,Lin_Sp_2011,Lin_Fa_2011,Lin_Sp_2012,Lin_Fa_2012,Lin_Sp_2013,Lin_Fa_2013,Lin_Sp_2014,Lin_Fa_2014,Lin_Sp_2015,Lin_Fa_2015)
Lin_avg_line <- sapply(Lin_line, mean)

#Created a dataset for average pigmentation scores for each population: year collected, spring/fall timepoint (numbered 1-12), average pigmentation scores
Lin_years <- c("2010Sp", "2010Fa", "2011Sp", "2011Fa", "2012Sp", "2012Fa", "2013Sp", "2013Fa", "2014Sp", "2014Fa", "2015Sp", "2015Fa")
Lin_timepoints <- c(1,2,3,4,5,6,7,8,9,10,11,12)
Lin_data <- data.frame(Lin_years,Lin_timepoints,Lin_avg_line)

#Created a dataset showing pigmentation scores for each isofemale line: year collected, spring/fall timepoint (numbered 1-12), individual data points (pts) for isofemale line pigmentation scores. 
#Will use this to plot raw data points for isofemale lines on the seasonal graph.
Year_points_Lin <- c(rep("2010Sp" , 10) , rep("2010Fa" , 10) , rep("2011Sp" , 10) , rep("2011Fa" , 10) , rep("2012Sp" , 10) , rep("2012Fa" , 10), rep("2013Sp" , 10) , rep("2013Fa" , 10), rep("2014Sp" , 20) , rep("2014Fa" , 20), rep("2015Sp" , 20) , rep("2015Fa" , 20))
Tmpt_points_Lin <- c(rep("1",10), rep("2",10), rep("3",10), rep("4",10), rep("5",10), rep("6",10), rep("7",10), rep("8",10), rep("9",20), rep("10",20), rep("11",20), rep("12",20))
Pts_points_Lin <- c(Lin_Sp_2010,Lin_Fa_2010,Lin_Sp_2011,Lin_Fa_2011,Lin_Sp_2012,Lin_Fa_2012,Lin_Sp_2013,Lin_Fa_2013,Lin_Sp_2014,Lin_Fa_2014,Lin_Sp_2015,Lin_Fa_2015)
Data_points_Lin <- data.frame(Year_points_Lin,Tmpt_points_Lin,Pts_points_Lin)


#Plotting seasonal pigmentation patterns in wild populations: plotted average pigmentation scores and isofemale line data points
ggplot(Lin_data, aes(x = factor(Lin_years, level = Lin_years), y = Lin_avg_line, group=1)) +
  geom_line(color="#823830")+
  geom_point(color="#823830")+
  theme_classic()+
  ggtitle("Seasonal Patterns (Media, PA)") +
  xlab("Timepoint")+
  ylab("Average Pigmentation Score") +
  #Adding error bars: standard error calculated in Excel
  geom_errorbar(aes(x=1, ymin = 31.8506469, ymax = 33.2541131), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=2, ymin = 22.3066231, ymax = 23.7511569), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=3, ymin = 27.360592, ymax = 28.644408), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=4, ymin = 15.532629, ymax = 16.967371), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=5, ymin = 25.7426731, ymax = 26.9729669), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=6, ymin = 13.9457301, ymax = 15.0092699), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=7, ymin = 22.3138941, ymax = 23.4661059), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=8, ymin = 15.3868551, ymax = 16.2664849), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=9, ymin = 21.8922791, ymax = 23.0277209), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=10, ymin = 21.3854988, ymax = 22.5145012), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=11, ymin = 30.2925419, ymax = 31.4674581), size=.45,color="#823830",width=0.05) +
  geom_errorbar(aes(x=12, ymin = 18.4038944, ymax = 19.6761056), size=.45,color="#823830",width=0.05) +
  #data points for individual isofemale lines
  geom_point(data=Data_points_Lin, aes(x=Year_points_Lin, y=Pts_points_Lin), alpha=1/4, size=0.05)


#Statistics: 2010-2015 Linvilla Seasonal Linear Model
#Made a dataset of mean pigmentation scores for isofemale lines at each year/timepoint, and then ran a linear model
#pigmentation scores ~ year * timepoint
Year_points_LinLM <- c(rep("2010" , 20) , rep("2011" , 20) , rep("2012" , 20) , rep("2013" , 20), rep("2014" , 40) , rep("2015" , 40))
Tmpt_points_LinLM <- c(rep("1",10), rep("2",10), rep("1",10), rep("2",10), rep("1",10), rep("2",10), rep("1",10), rep("2",10), rep("1",20), rep("2",20), rep("1",20), rep("2",20))
Pts_points_LinLM <- c(Lin_Sp_2010,Lin_Fa_2010,Lin_Sp_2011,Lin_Fa_2011,Lin_Sp_2012,Lin_Fa_2012,Lin_Sp_2013,Lin_Fa_2013,Lin_Sp_2014,Lin_Fa_2014,Lin_Sp_2015,Lin_Fa_2015)
Data_points_LinLM <- data.frame(Year_points_LinLM,Tmpt_points_LinLM,Pts_points_LinLM)
#Seasonal Linvilla Linear Model
LM_LinLM <- lm(formula = Pts_points_LinLM ~ Year_points_LinLM * Tmpt_points_LinLM, data = Data_points_LinLM)
summary(LM_LinLM)



#####
#####Seasonal pigmentation patterns in experimental populations: Plot and Statistics
#Experiment run 2016 in Philadelphia, Pennsylvania (Schmidt lab experimental orchard)
#Entered pigmentation scores from Excel scoring sheet and assigned to seasonal populations in each year. 
#20 females per mesocosm per timepoint were scored and entered here. 
#Mesocosms were labled "S1 - S10", and flies from S2 were not available to score (N=9 mesocosms).

#Enter pigmentation scores
Founder_2016 <- c(4,12,13,5,9,7,9,8,13,13,14,15,12,10,10,12,13,14,10,10,14,12,13,12,10,12,13,8,12,8,12,10,12,13,9,15,11,12,12,13,13,15,7,12,13,11,8,15,12,14,15,14,10,15,13,17,11,19,11,14,12,13,12,12,11,14,14,12,13,11,13,17,12,15,12,13,14,12,14,11)
##Note: calculated founder average pigmentation score = 12.025
S1_Summer_16 <- c(15,14,8,9,11,11,5,7,9,12,10,15,12,11,16,11,7,10,9,14)
S1_Fall_16 <- c(8,7,6,7,8,11,7,3,8,13,4,11,10,11,11,11,12,7,8,12)
S3_Summer_16 <- c(12,7,12,13,12,10,12,14,17,13,12,7,19,10,10,12,15,14,14,9)
S3_Fall_16 <- c(12,7,11,9,13,12,7,9,12,9,9,10,11,6,12,5,13,8,6,7)
S4_Summer_16 <- c(12,16,14,12,17,13,9,13,10,11,9,15,12,17,16,14,14,16,16,13)
S4_Fall_16 <- c(8,7,9,10,11,7,8,9,7,12,7,7,11,11,8,7,15,10,9,9)
S5_Summer_16 <- c(12,12,11,13,13,12,13,21,11,9,17,17,8,12,12,12,16,15,13,11)
S5_Fall_16 <- c(15,8,8,6,14,12,9,12,10,10,5,10,15,11,9,12,11,9,6,9)
S6_Summer_16 <- c(13,13,15,6,18,16,10,9,19,11,15,12,12,12,8,8,11,10,10,15)
S6_Fall_16 <- c(10,12,6,9,11,11,8,8,10,8,16,4,10,8,10,5,6,5,6,8)
S7_Summer_16 <- c(12,19,16,10,13,11,14,15,14,11,10,15,16,13,18,21,9,14,14,13)
S7_Fall_16 <- c(11,8,13,8,11,9,9,9,6,9,12,9,4,12,7,5,7,10,6,8)
S8_Summer_16 <- c(8,12,10,15,7,14,6,18,14,9,9,11,16,7,12,17,10,11,14,17)
S8_Fall_16 <- c(5,8,8,6,6,10,8,8,8,15,12,5,9,9,9,7,14,13,5,10)
S9_Summer_16 <- c(12,8,11,20,14,20,11,10,15,19,7,8,20,11,3,12,16,10,12,12)
S9_Fall_16 <- c(13,8,6,13,11,4,3,7,8,4,4,9,8,11,9,5,7,10,6,9)
S10_Summer_16 <- c(11,16,12,10,13,20,14,14,15,6,10,7,8,20,10,18,7,12,20,19)
S10_Fall_16 <- c(8,9,7,9,9,6,13,4,14,8,15,9,12,13,8,11,6,13,11,5)

#Entered average summer and fall pigmentation scores calculated in Excel
Avg_2016S_summer_point = 12.56
Avg_2016S_fall_point = 8.94

#Made a list of each seasonal population and computed average score from raw data for each population
Orchard_2016_line <- list(S1_Summer_16, S1_Fall_16, S3_Summer_16, S3_Fall_16, S4_Summer_16, S4_Fall_16, S5_Summer_16, S5_Fall_16, S6_Summer_16, S6_Fall_16, S7_Summer_16, S7_Fall_16, S8_Summer_16, S8_Fall_16, S9_Summer_16, S9_Fall_16, S10_Summer_16, S10_Fall_16)
Orchard_2016_avg_line <- sapply(Orchard_2016_line, mean)

#Created a dataset for average pigmentation scores for each mesocosm population: mesocosm, timepoint ("2" = end of summer, "3" = end of fall), average pigmentation scores
Orchard_2016_mesocosms <- c(rep("S1" , 2) , rep("S3" , 2) , rep("S4" , 2) , rep("S5" , 2) , rep("S6" , 2) , rep("S7" , 2) , rep("S8" , 2) , rep("S9" , 2) , rep("S10" , 2) )
Orchard_2016_timepoints <- rep(c("2" , "3") , 9)
Orchard_2016_scores <- Orchard_2016_avg_line
Orchard_2016_data <- data.frame(Orchard_2016_mesocosms,Orchard_2016_timepoints,Orchard_2016_scores)

#Created a dataset for raw data points showing individual pigmentation scores in each population: mesocosm collected, timepoint (end of summer = "2", end of fall = "3"), individual data points (pts) for pigmentation scores. 
#Did not plot these points; used this dataset for the linear model.
Mesocosm_points_Orch2016 <- c(rep("S1" , 40) , rep("S3" , 40) , rep("S4" , 40) , rep("S5" , 40) , rep("S6" , 40) , rep("S7" , 40) , rep("S8" , 40) , rep("S9" , 40) , rep("S10" , 40) )
Tmpt_points_Orch2016 <- rep(c("2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3"), 9)
Pts_points_Orch2016 <- c(S1_Summer_16,S1_Fall_16,S3_Summer_16,S3_Fall_16,S4_Summer_16,S4_Fall_16,S5_Summer_16,S5_Fall_16,S6_Summer_16,S6_Fall_16,S7_Summer_16,S7_Fall_16,S8_Summer_16,S8_Fall_16,S9_Summer_16,S9_Fall_16,S10_Summer_16,S10_Fall_16)
Data_points_Orch2016 <- data.frame(Mesocosm_points_Orch2016,Tmpt_points_Orch2016,Pts_points_Orch2016)


#Experimental Orchard 2016: plotted average pigmentation scores across seasons
#Founder average pigmentation score: 12.025
#Founder standard error: 0.282940594
#September average pigmentation score: 12.55555556
#September standard error: 0.2631073238
#November average pigmentation score: 8.944444444
#November standard error: 0.2046514079
ggplot(Orchard_2016_data, aes(x = Orchard_2016_timepoints, y = Orchard_2016_scores, group=Orchard_2016_mesocosms, color=Orchard_2016_mesocosms)) +
  geom_line() +
  xlab("Timepoint") +
  ylab("Average Pigmentation Score") +
  ggtitle("Experimental Orchard: 2016")+
  theme_classic() +
  ylim(7,14.7) +
  scale_x_discrete(labels=c("August","October"))+
  scale_color_manual(values=c("#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F","#34434F"))+
  theme(legend.position = "none")+
  #Adding average points
  #geom_point(aes(x=1,y=12.025),colour="#823830",size=4,pch=17) +
  geom_point(aes(x=1,y=Avg_2016S_summer_point),colour="#823830",size=2.5)+
  geom_point(aes(x=2,y=Avg_2016S_fall_point),colour="#823830",size=2.5)+
  #Adding line segment between summer and fall points
  geom_segment(mapping=aes(x=1, y=Avg_2016S_summer_point, xend=2, yend=Avg_2016S_fall_point), arrow=NULL, size=2, color="#823830") +
  #Adding founder point
  geom_point(aes(x=0.6, y=12.025),colour="#823830",size=2.5) +
  #Founder standard error bars
  geom_errorbar(aes(x=0.6, ymin = 11.74205941, ymax = 12.30794059), size=.45,color="#823830",width=0.03) +
  #September standard error bars
  geom_errorbar(aes(x=1, ymin = 12.29244823, ymax = 12.81866288), size=.45,color="#823830",width=0.03) +
  #November standard error bars
  geom_errorbar(aes(x=2, ymin = 8.739793037, ymax = 9.149095852), size=.45,color="#823830",width=0.03)


#Statistics: Experimental Orchard 2016
#Linear mixed effects model: lme(pigmentation scores ~ timepoint, random=~1|mesocosm/timepoint)
#Tested for significance using ANOVA
LM_Orch2016 <- lme(Pts_points_Orch2016 ~ Tmpt_points_Orch2016, random=~1|Mesocosm_points_Orch2016/Tmpt_points_Orch2016, data = Data_points_Orch2016)
anova(LM_Orch2016)


