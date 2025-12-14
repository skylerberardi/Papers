##### Skyler Berardi; contact: berardis@sas.upenn.edu

##### pigmentation measured from flies directly sampled from mesocosms (no common garden treatment)
##### 2019 interspecific competition experiment


#Packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(lme4)
library(lmerTest)
library(readxl)
library(emmeans)
library(plotrix)


#####
##### Interspecific Competition (2019)

## Read competition 2019 dataset
field19_pig_table <- read_excel("Competition2019_FieldPigmentation_R.xlsx")

pigmentation_levels <- c("0","1","2")
field19_pig_table <- field19_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
field19_pig_table <- field19_pig_table %>% mutate_at('cage', as.factor)
field19_pig_table <- field19_pig_table %>% mutate_at('treatment', as.factor)
field19_pig_table$timepointID_cage <- paste(field19_pig_table$timepoint_ID, "-", field19_pig_table$cage)

view(field19_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

field19_pig_table_orchard <- field19_pig_table %>% filter(treatment != 'Founder')
view(field19_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

field19_pig_table_founder <- field19_pig_table %>% filter(treatment == 'Founder')
view(field19_pig_table_founder)

field19_pig_founder_mean <- mean(field19_pig_table_founder$pigmentation_score)
print(field19_pig_founder_mean)

field19_pig_founder_SE <- std.error(field19_pig_table_founder$pigmentation_score)
print(field19_pig_founder_SE)


#Calculating mean pigmentation score for each cage
field19_pig_cage_means <- field19_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(field19_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_field19_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = field19_pig_table_orchard)
anova(LM_field19_pig)
summary(LM_field19_pig)
confint(LM_field19_pig)
plot(LM_field19_pig)
qqnorm(resid(LM_field19_pig))
qqline(resid(LM_field19_pig))

# Estimated marginal means
emm_field19 <- emmeans(LM_field19_pig, ~ timepoint_ID*treatment)
print(emm_field19)
plot(emm_field19)

# Pairwise comparisons

# All comparisons
pairs(emm_field19)

# Make list of planned comparisons
field19_contrasts <- list(
  #control common garden over time
  E_T1.T2 = c(1, -1, 0, 0),
  #control field over time
  F_T1.T2 = c(0, 0, 1, -1),
  #TP1 control common garden vs. control field
  T1_E.F = c(1, 0, -1, 0),
  #TP2 control common garden vs. control field
  T2_E.F = c(0, 1, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_field19 = contrast(emm_field19, method = field19_contrasts, adjust = "none")
print(contrasts_field19)
confint(contrasts_field19)

# Test contrasts with a Holm correction
contrasts_holm_field19 = contrast(emm_field19, method = field19_contrasts, adjust = "holm")
print(contrasts_holm_field19)


## Extracting means and standard errors from field19 model for plotting
LM_field19_pig_means <- emmeans(LM_field19_pig, ~ timepoint_ID*treatment)
print(LM_field19_pig_means)
LM_field19_pig_means <- as.data.frame(LM_field19_pig_means)

view(LM_field19_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_field19_pig_plot <- subset(LM_field19_pig_means, select = -c(df, lower.CL, upper.CL))
LM_field19_pig_plot <- LM_field19_pig_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", emmean = field19_pig_founder_mean, SE = field19_pig_founder_SE, .before = 1)
LM_field19_pig_plot <- LM_field19_pig_plot %>% add_row(timepoint_ID = "0", treatment = "E", emmean = field19_pig_founder_mean, SE = field19_pig_founder_SE, .before = 1)

view(LM_field19_pig_plot)

# Update individual cage dataframe to include founder
field19_pig_cage_plot <- subset(field19_pig_cage_means, select = -c(timepointID_cage))
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E23.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E22.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E20.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E16.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E11.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E7.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E.field", cage = "E3.field", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E23", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E22", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E20", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E16", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E11", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E7", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)
field19_pig_cage_plot <- field19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E3", cage_mean_pigscore = field19_pig_founder_mean, .before = 1)

view(field19_pig_cage_plot)

# Plot
ggplot(field19_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ylim(8.5,24.5) +
  ggtitle("Control Populations (2019)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 9","September 11","November 8")) +
  scale_color_manual(name="Treatment", values=c("green3","darkgreen"), labels=c("Common Garden \nTreated","Field Collected")) +
  #Adding average points
  geom_line(data=LM_field19_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_field19_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_field19_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=field19_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=field19_pig_founder_mean-field19_pig_founder_SE, ymax=field19_pig_founder_mean+field19_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)








