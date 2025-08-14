##### Code used to generate plots and statistics for phenotypic analyses
##### Skyler Berardi; contact: berardis@sas.upenn.edu

##### Abiotic and biotic drivers of pigmentation evolution
# Intraspecific Competition (2022)
# Temperature (2022)
# Interspecific Competition (2019)
# Intraspecific Competition (2017)
# Diet (2020)
# Microbial Additions (2020)

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
##### Intraspecific competition (reduced population density; 2022)

## Read density 2022 dataset
den22_pig_table <- read_excel("Density2022_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","1","2")
den22_pig_table <- den22_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
den22_pig_table <- den22_pig_table %>% mutate_at('cage', as.factor)
den22_pig_table <- den22_pig_table %>% mutate_at('treatment', as.factor)
den22_pig_table$timepointID_cage <- paste(den22_pig_table$timepoint_ID, "-", den22_pig_table$cage)

view(den22_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

den22_pig_table_orchard <- den22_pig_table %>% filter(treatment != 'Founder')
view(den22_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

den22_pig_table_founder <- den22_pig_table %>% filter(treatment == 'Founder')
view(den22_pig_table_founder)

den22_pig_founder_mean <- mean(den22_pig_table_founder$pigmentation_score)
print(den22_pig_founder_mean)

den22_pig_founder_SE <- std.error(den22_pig_table_founder$pigmentation_score)
print(den22_pig_founder_SE)


#Calculating mean pigmentation score for each cage
den22_pig_cage_means <- den22_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(den22_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_den22_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = den22_pig_table_orchard)
anova(LM_den22_pig)
summary(LM_den22_pig)
confint(LM_den22_pig)
plot(LM_den22_pig)
qqnorm(resid(LM_den22_pig))
qqline(resid(LM_den22_pig))

# Estimated marginal means
emm_den22 <- emmeans(LM_den22_pig, ~ timepoint_ID*treatment)
print(emm_den22)
plot(emm_den22)

# Pairwise comparisons

# All comparisons
pairs(emm_den22)

# Make list of planned comparisons
den22_contrasts <- list(
  #control over time
  E_T1.T2 = c(1, -1, 0, 0),
  #density over time
  N_T1.T2 = c(0, 0, 1, -1),
  #TP1 control vs. density
  T1_E.N = c(1, 0, -1, 0),
  #TP2 control vs. density
  T2_E.N = c(0, 1, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_den22 = contrast(emm_den22, method = den22_contrasts, adjust = "none")
print(contrasts_den22)
confint(contrasts_den22)


# Test contrasts with a Holm correction
contrasts_holm_den22 = contrast(emm_den22, method = den22_contrasts, adjust = "holm")
print(contrasts_holm_den22)


## Extracting means and standard errors from den22 model for plotting
LM_den22_pig_means <- emmeans(LM_den22_pig, ~ timepoint_ID*treatment)
print(LM_den22_pig_means)
LM_den22_pig_means <- as.data.frame(LM_den22_pig_means)

view(LM_den22_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_den22_pig_plot <- subset(LM_den22_pig_means, select = -c(df, lower.CL, upper.CL))
LM_den22_pig_plot <- LM_den22_pig_plot %>% add_row(timepoint_ID = "0", treatment = "N", emmean = den22_pig_founder_mean, SE = den22_pig_founder_SE, .before = 1)
LM_den22_pig_plot <- LM_den22_pig_plot %>% add_row(timepoint_ID = "0", treatment = "E", emmean = den22_pig_founder_mean, SE = den22_pig_founder_SE, .before = 1)

view(LM_den22_pig_plot)

# Update individual cage dataframe to include founder
den22_pig_cage_plot <- subset(den22_pig_cage_means, select = -c(timepointID_cage))
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N9", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N8", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N7", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N6", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N5", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N4", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N3", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N2", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "N", cage = "N1", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E9", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E8", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E7", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E6", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E5", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E4", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E3", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E2", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)
den22_pig_cage_plot <- den22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E1", cage_mean_pigscore = den22_pig_founder_mean, .before = 1)

view(den22_pig_cage_plot)

# Plot
ggplot(den22_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Intraspecific Competition (2022)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  ylim(8,16) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 6","September 7","November 8")) +
  scale_color_manual(name="Treatment", breaks = , values=c("darkorange4","orange"), labels = c("Control", "Reduced Density")) +
  #Adding average points
  geom_line(data=LM_den22_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_den22_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_den22_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=den22_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=den22_pig_founder_mean-den22_pig_founder_SE, ymax=den22_pig_founder_mean+den22_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)





#####
##### Temperature (2022)

## Read temperature 2022 dataset
tmp22_pig_table <- read_excel("Warming2022_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","1","2")
tmp22_pig_table <- tmp22_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
tmp22_pig_table <- tmp22_pig_table %>% mutate_at('cage', as.factor)
tmp22_pig_table <- tmp22_pig_table %>% mutate_at('treatment', as.factor)
tmp22_pig_table$timepointID_cage <- paste(tmp22_pig_table$timepoint_ID, "-", tmp22_pig_table$cage)

view(tmp22_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

tmp22_pig_table_orchard <- tmp22_pig_table %>% filter(treatment != 'Founder')
view(tmp22_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

tmp22_pig_table_founder <- tmp22_pig_table %>% filter(treatment == 'Founder')
view(tmp22_pig_table_founder)

tmp22_pig_founder_mean <- mean(tmp22_pig_table_founder$pigmentation_score)
print(tmp22_pig_founder_mean)

tmp22_pig_founder_SE <- std.error(tmp22_pig_table_founder$pigmentation_score)
print(tmp22_pig_founder_SE)


#Calculating mean pigmentation score for each cage
tmp22_pig_cage_means <- tmp22_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(tmp22_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_tmp22_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = tmp22_pig_table_orchard)
anova(LM_tmp22_pig)
summary(LM_tmp22_pig)
confint(LM_tmp22_pig)
plot(LM_tmp22_pig)
qqnorm(resid(LM_tmp22_pig))
qqline(resid(LM_tmp22_pig))

# Estimated marginal means
emm_tmp22 <- emmeans(LM_tmp22_pig, ~ timepoint_ID*treatment)
print(emm_tmp22)
plot(emm_tmp22)

# Pairwise comparisons

# All comparisons
pairs(emm_tmp22)

# Make list of planned comparisons
tmp22_contrasts <- list(
  #control over time
  E_T1.T2 = c(1, -1, 0, 0),
  #warming over time
  W_T1.T2 = c(0, 0, 1, -1),
  #TP1 control vs. warming
  T1_E.W = c(1, 0, -1, 0),
  #TP2 control vs. warming
  T2_E.W = c(0, 1, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_tmp22 = contrast(emm_tmp22, method = tmp22_contrasts, adjust = "none")
print(contrasts_tmp22)
confint(contrasts_tmp22)

# Test contrasts with a Holm correction
contrasts_holm_tmp22 = contrast(emm_tmp22, method = tmp22_contrasts, adjust = "holm")
print(contrasts_holm_tmp22)


## Extracting means and standard errors from tmp22 model for plotting
LM_tmp22_pig_means <- emmeans(LM_tmp22_pig, ~ timepoint_ID*treatment)
print(LM_tmp22_pig_means)
LM_tmp22_pig_means <- as.data.frame(LM_tmp22_pig_means)

view(LM_tmp22_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_tmp22_pig_plot <- subset(LM_tmp22_pig_means, select = -c(df, lower.CL, upper.CL))
LM_tmp22_pig_plot <- LM_tmp22_pig_plot %>% add_row(timepoint_ID = "0", treatment = "W", emmean = tmp22_pig_founder_mean, SE = tmp22_pig_founder_SE, .before = 1)
LM_tmp22_pig_plot <- LM_tmp22_pig_plot %>% add_row(timepoint_ID = "0", treatment = "E", emmean = tmp22_pig_founder_mean, SE = tmp22_pig_founder_SE, .before = 1)

view(LM_tmp22_pig_plot)

# Update individual cage dataframe to include founder
tmp22_pig_cage_plot <- subset(tmp22_pig_cage_means, select = -c(timepointID_cage))
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W9", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W8", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W7", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W6", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W5", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W4", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W3", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W2", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "W", cage = "W1", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E9", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E8", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E7", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E6", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E5", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E4", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E3", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E2", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)
tmp22_pig_cage_plot <- tmp22_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E1", cage_mean_pigscore = tmp22_pig_founder_mean, .before = 1)

view(tmp22_pig_cage_plot)

# Plot
ggplot(tmp22_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ylim(8,16) +
  ggtitle("Increased Temperature (2022)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 6","September 7","November 8")) +
  scale_color_manual(name="Treatment", values=c("blue2","red2"), labels = c("Control", "Warming")) +
  #Adding average points
  geom_line(data=LM_tmp22_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_tmp22_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_tmp22_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=tmp22_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=tmp22_pig_founder_mean-tmp22_pig_founder_SE, ymax=tmp22_pig_founder_mean+tmp22_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)





#####
##### Interspecific Competition (2019)

## Read competition 2019 dataset
cmp19_pig_table <- read_excel("Competition2019_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","1","2")
cmp19_pig_table <- cmp19_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
cmp19_pig_table <- cmp19_pig_table %>% mutate_at('cage', as.factor)
cmp19_pig_table <- cmp19_pig_table %>% mutate_at('treatment', as.factor)
cmp19_pig_table$timepointID_cage <- paste(cmp19_pig_table$timepoint_ID, "-", cmp19_pig_table$cage)

view(cmp19_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

cmp19_pig_table_orchard <- cmp19_pig_table %>% filter(treatment != 'Founder')
view(cmp19_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

cmp19_pig_table_founder <- cmp19_pig_table %>% filter(treatment == 'Founder')
view(cmp19_pig_table_founder)

cmp19_pig_founder_mean <- mean(cmp19_pig_table_founder$pigmentation_score)
print(cmp19_pig_founder_mean)

cmp19_pig_founder_SE <- std.error(cmp19_pig_table_founder$pigmentation_score)
print(cmp19_pig_founder_SE)


#Calculating mean pigmentation score for each cage
cmp19_pig_cage_means <- cmp19_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(cmp19_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_cmp19_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = cmp19_pig_table_orchard)
anova(LM_cmp19_pig)
summary(LM_cmp19_pig)
confint(LM_cmp19_pig)
plot(LM_cmp19_pig)
qqnorm(resid(LM_cmp19_pig))
qqline(resid(LM_cmp19_pig))

# Estimated marginal means
emm_cmp19 <- emmeans(LM_cmp19_pig, ~ timepoint_ID*treatment)
print(emm_cmp19)
plot(emm_cmp19)

# Pairwise comparisons

# All comparisons
pairs(emm_cmp19)

# Make list of planned comparisons
cmp19_contrasts <- list(
  #control over time
  E_T1.T2 = c(0, 0, 1, -1),
  #competition over time
  C_T1.T2 = c(1, -1, 0, 0),
  #TP1 control vs. competition
  T1_E.C = c(-1, 0, 1, 0),
  #TP2 control vs. competition
  T2_E.C = c(0, -1, 0, 1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_cmp19 = contrast(emm_cmp19, method = cmp19_contrasts, adjust = "none")
print(contrasts_cmp19)
confint(contrasts_cmp19)

# Test contrasts with a Holm correction
contrasts_holm_cmp19 = contrast(emm_cmp19, method = cmp19_contrasts, adjust = "holm")
print(contrasts_holm_cmp19)


## Extracting means and standard errors from cmp19 model for plotting
LM_cmp19_pig_means <- emmeans(LM_cmp19_pig, ~ timepoint_ID*treatment)
print(LM_cmp19_pig_means)
LM_cmp19_pig_means <- as.data.frame(LM_cmp19_pig_means)

view(LM_cmp19_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_cmp19_pig_plot <- subset(LM_cmp19_pig_means, select = -c(df, lower.CL, upper.CL))
LM_cmp19_pig_plot <- LM_cmp19_pig_plot %>% add_row(timepoint_ID = "0", treatment = "C", emmean = cmp19_pig_founder_mean, SE = cmp19_pig_founder_SE, .before = 1)
LM_cmp19_pig_plot <- LM_cmp19_pig_plot %>% add_row(timepoint_ID = "0", treatment = "E", emmean = cmp19_pig_founder_mean, SE = cmp19_pig_founder_SE, .before = 1)

view(LM_cmp19_pig_plot)

# Update individual cage dataframe to include founder
cmp19_pig_cage_plot <- subset(cmp19_pig_cage_means, select = -c(timepointID_cage))
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C28", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C24", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C19", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C15", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C10", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "C", cage = "C4", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E23", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E22", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E20", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E16", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E11", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E7", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)
cmp19_pig_cage_plot <- cmp19_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "E", cage = "E3", cage_mean_pigscore = cmp19_pig_founder_mean, .before = 1)

view(cmp19_pig_cage_plot)

# Plot
ggplot(cmp19_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ylim(8,16) +
  ggtitle("Interspecific Competition (2019)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 9","September 11","November 8")) +
  scale_color_manual(name="Treatment", values=c("darkgreen","green3"), labels=c("Competition","Control")) +
  #Adding average points
  geom_line(data=LM_cmp19_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_cmp19_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_cmp19_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=cmp19_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=cmp19_pig_founder_mean-cmp19_pig_founder_SE, ymax=cmp19_pig_founder_mean+cmp19_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)





#####
##### Intraspecific Competition (reduced population density; 2017)

## Read density 2017 dataset
den17_pig_table <- read_excel("Density2017_Pigmentation_R.xlsx")

pigmentation_levels <- c("1","2","3","4")
den17_pig_table <- den17_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
den17_pig_table <- den17_pig_table %>% mutate_at('cage', as.factor)
den17_pig_table <- den17_pig_table %>% mutate_at('treatment', as.factor)
den17_pig_table$timepointID_cage <- paste(den17_pig_table$timepoint_ID, "-", den17_pig_table$cage)

view(den17_pig_table)

#Note: no founder flies available to score.

#Calculating mean pigmentation score for each cage
den17_pig_cage_means <- den17_pig_table %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(den17_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_den17_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = den17_pig_table)
anova(LM_den17_pig)
summary(LM_den17_pig)
confint(LM_den17_pig)
plot(LM_den17_pig)
qqnorm(resid(LM_den17_pig))
qqline(resid(LM_den17_pig))

# Estimated marginal means
emm_den17 <- emmeans(LM_den17_pig, ~ timepoint_ID*treatment)
print(emm_den17)
plot(emm_den17)

# Pairwise comparisons

# All comparisons
pairs(emm_den17)

# Make list of planned comparisons
den17_contrasts <- list(
  #control over time (end summer to end fall)
  E_T2.T4 = c(0, 1, 0, -1, 0, 0, 0, 0),
  #density over time (end summer to end fall)
  D_T2.T4 = c(0, 0, 0, 0, 0, 1, 0, -1),
  #TP2 control vs. density (end summer)
  T2_E.D = c(0, 1, 0, 0, 0, -1, 0, 0),
  #TP4 control vs. density (end fall)
  T4_E.D = c(0, 0, 0, 1, 0, 0, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_den17 = contrast(emm_den17, method = den17_contrasts, adjust = "none")
print(contrasts_den17)
confint(contrasts_den17)

# Test contrasts with a Holm correction
contrasts_holm_den17 = contrast(emm_den17, method = den17_contrasts, adjust = "holm")
print(contrasts_holm_den17)


## Extracting means and standard errors from den17 model for plotting
LM_den17_pig_means <- emmeans(LM_den17_pig, ~ timepoint_ID*treatment)
print(LM_den17_pig_means)
LM_den17_pig_means <- as.data.frame(LM_den17_pig_means)

view(LM_den17_pig_means)


## Pigmentation plot
ggplot(den17_pig_cage_means, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Intraspecific Competition (2017)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("August 8","September 22", "October 18", "November 12")) +
  scale_color_manual(values=c("darkorange4","orange"), labels=c("Control","Reduced Density")) +
  #Adding average points
  geom_line(data=LM_den17_pig_means, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_den17_pig_means, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_den17_pig_means, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) 





#####
##### Diet (2020)

## Read nutritional quality 2020 dataset
nq20_pig_table <- read_excel("Nutrition2020_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","1","2","3","4","5")
nq20_pig_table <- nq20_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
nq20_pig_table <- nq20_pig_table %>% mutate_at('cage', as.factor)
nq20_pig_table <- nq20_pig_table %>% mutate_at('treatment', as.factor)
nq20_pig_table$timepointID_cage <- paste(nq20_pig_table$timepoint_ID, "-", nq20_pig_table$cage)

view(nq20_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

nq20_pig_table_orchard <- nq20_pig_table %>% filter(treatment != 'Founder_control')
nq20_pig_table_orchard <- nq20_pig_table_orchard %>% filter(treatment != 'Founder_apple')

view(nq20_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

#control
nq20_pig_table_founder_control <- nq20_pig_table %>% filter(treatment == 'Founder_control')
view(nq20_pig_table_founder_control)

nq20_pig_founder_control_mean <- mean(nq20_pig_table_founder_control$pigmentation_score)
print(nq20_pig_founder_control_mean)

nq20_pig_founder_control_SE <- std.error(nq20_pig_table_founder_control$pigmentation_score)
print(nq20_pig_founder_control_SE)

#apple
nq20_pig_table_founder_apple <- nq20_pig_table %>% filter(treatment == 'Founder_apple')
view(nq20_pig_table_founder_apple)

nq20_pig_founder_apple_mean <- mean(nq20_pig_table_founder_apple$pigmentation_score)
print(nq20_pig_founder_apple_mean)

nq20_pig_founder_apple_SE <- std.error(nq20_pig_table_founder_apple$pigmentation_score)
print(nq20_pig_founder_apple_SE)


#Calculating mean pigmentation score for each cage
nq20_pig_cage_means <- nq20_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(nq20_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_nq20_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = nq20_pig_table_orchard)
anova(LM_nq20_pig)
summary(LM_nq20_pig)
confint(LM_nq20_pig)
plot(LM_nq20_pig)
qqnorm(resid(LM_nq20_pig))
qqline(resid(LM_nq20_pig))

# Estimated marginal means
emm_nq20 <- emmeans(LM_nq20_pig, ~ timepoint_ID*treatment)
print(emm_nq20)
plot(emm_nq20)

# Pairwise comparisons

# All comparisons
pairs(emm_nq20)

# Make list of planned comparisons
nq20_contrasts <- list(
  #control over time (end of summer to end of fall)
  E_T2.T5 = c(0, 0, 0, 0, 0, 0, 1, 0, 0, -1),
  #low quality (apple) over time (end of summer to end of fall)
  A_T2.T5 = c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0),
  #TP2 control vs. low quality (apple) (end of summer)
  T2_E.A = c(0, -1, 0, 0, 0, 0, 1, 0, 0, 0),
  #TP5 control vs. low quality (apple) (end of fall)
  T5_E.A = c(0, 0, 0, 0, -1, 0, 0, 0, 0, 1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_nq20 = contrast(emm_nq20, method = nq20_contrasts, adjust = "none")
print(contrasts_nq20)
confint(contrasts_nq20)

# Test contrasts with a Holm correction
contrasts_holm_nq20 = contrast(emm_nq20, method = nq20_contrasts, adjust = "holm")
print(contrasts_holm_nq20)


## Extracting means and standard errors from nq20 model for plotting
LM_nq20_pig_means <- emmeans(LM_nq20_pig, ~ timepoint_ID*treatment)
print(LM_nq20_pig_means)
LM_nq20_pig_means <- as.data.frame(LM_nq20_pig_means)

view(LM_nq20_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_nq20_pig_plot <- subset(LM_nq20_pig_means, select = -c(df, lower.CL, upper.CL))
LM_nq20_pig_plot <- LM_nq20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "Control", emmean = nq20_pig_founder_control_mean, SE = nq20_pig_founder_control_SE, .before = 1)
LM_nq20_pig_plot <- LM_nq20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", emmean = nq20_pig_founder_apple_mean, SE = nq20_pig_founder_apple_SE, .before = 1)

view(LM_nq20_pig_plot)

# Update individual cage dataframe to include founder
nq20_pig_cage_plot <- subset(nq20_pig_cage_means, select = -c(timepointID_cage))
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B6b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B5b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B4b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B3b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B2b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "B1b", cage_mean_pigscore = nq20_pig_founder_control_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A6a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A5a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A4a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A3a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A2a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)
nq20_pig_cage_plot <- nq20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Apple", cage = "A1a", cage_mean_pigscore = nq20_pig_founder_apple_mean, .before = 1)

view(nq20_pig_cage_plot)

# Plot
ggplot(nq20_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ylim(8,16) +
  ggtitle("Diet (2020)") +
  theme_classic() +
  theme(text = element_text(size=16), axis.text.x = element_text(size=10)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 15","August 13","September 7","September 30","November 9","November 24")) +
  scale_color_manual(name="Treatment", values=c("#FFC107","#CC3A7C"), labels = c("Apple", "Control")) +
  #Adding average points
  geom_line(data=LM_nq20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_nq20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_nq20_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) 





#####
##### Resident Microbe Additions (2020)

## Read microbial additions 2020 dataset
mic20_pig_table <- read_excel("Microbes2020_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","2","4","5")
mic20_pig_table <- mic20_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
mic20_pig_table <- mic20_pig_table %>% mutate_at('cage', as.factor)
mic20_pig_table <- mic20_pig_table %>% mutate_at('treatment', as.factor)
mic20_pig_table$timepointID_cage <- paste(mic20_pig_table$timepoint_ID, "-", mic20_pig_table$cage)

view(mic20_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

mic20_pig_table_orchard <- mic20_pig_table %>% filter(treatment != 'Founder')
view(mic20_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

mic20_pig_table_founder <- mic20_pig_table %>% filter(treatment == 'Founder')
view(mic20_pig_table_founder)

mic20_pig_founder_mean <- mean(mic20_pig_table_founder$pigmentation_score)
print(mic20_pig_founder_mean)

mic20_pig_founder_SE <- std.error(mic20_pig_table_founder$pigmentation_score)
print(mic20_pig_founder_SE)


#Calculating mean pigmentation score for each cage
mic20_pig_cage_means <- mic20_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(mic20_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_mic20_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = mic20_pig_table_orchard)
anova(LM_mic20_pig)
summary(LM_mic20_pig)
confint(LM_mic20_pig)
plot(LM_mic20_pig)
qqnorm(resid(LM_mic20_pig))
qqline(resid(LM_mic20_pig))

# Estimated marginal means
emm_mic20 <- emmeans(LM_mic20_pig, ~ timepoint_ID*treatment)
print(emm_mic20)
plot(emm_mic20)

# Pairwise comparisons

# All comparisons
pairs(emm_mic20)

# Make list of planned comparisons
mic20_contrasts <- list(
  #control over time
  E_T2.T5 = c(0, 0, 0, 1, 0, -1, 0, 0, 0),
  #At over time
  At_T2.T5 = c(1, 0, -1, 0, 0, 0, 0, 0, 0),
  #Lb over time
  Lb_T2.T5 = c(0, 0, 0, 0, 0, 0, 1, 0, -1),
  #TP2 control vs. At
  T2_E.At = c(-1, 0, 0, 1, 0, 0, 0, 0, 0),
  #TP2 control vs. Lb
  T2_E.Lb = c(0, 0, 0, 1, 0, 0, -1, 0, 0),
  #TP2 At vs. Lb
  T2_At.Lb = c(1, 0, 0, 0, 0, 0, -1, 0, 0),
  #TP5 control vs. At
  T5_E.At = c(0, 0, -1, 0, 0, 1, 0, 0, 0),
  #TP5 control vs. Lb
  T5_E.Lb = c(0, 0, 0, 0, 0, 1, 0, 0, -1),
  #TP5 At vs. Lb
  T5_At.Lb = c(0, 0, 1, 0, 0, 0, 0, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_mic20 = contrast(emm_mic20, method = mic20_contrasts, adjust = "none")
print(contrasts_mic20)
confint(contrasts_mic20)

# Test contrasts with a Holm correction
contrasts_holm_mic20 = contrast(emm_mic20, method = mic20_contrasts, adjust = "holm")
print(contrasts_holm_mic20)


## Extracting means and standard errors from mic20 model for plotting
LM_mic20_pig_means <- emmeans(LM_mic20_pig, ~ timepoint_ID*treatment)
print(LM_mic20_pig_means)
LM_mic20_pig_means <- as.data.frame(LM_mic20_pig_means)

view(LM_mic20_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_mic20_pig_plot <- subset(LM_mic20_pig_means, select = -c(df, lower.CL, upper.CL))
LM_mic20_pig_plot <- LM_mic20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "AT", emmean = mic20_pig_founder_mean, SE = mic20_pig_founder_SE, .before = 1)
LM_mic20_pig_plot <- LM_mic20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "LB", emmean = mic20_pig_founder_mean, SE = mic20_pig_founder_SE, .before = 1)
LM_mic20_pig_plot <- LM_mic20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "Control", emmean = mic20_pig_founder_mean, SE = mic20_pig_founder_SE, .before = 1)

view(LM_mic20_pig_plot)

# Update individual cage dataframe to include founder
mic20_pig_cage_plot <- subset(mic20_pig_cage_means, select = -c(timepointID_cage))
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN12", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN10", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN8", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN6", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN4", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN2", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN11", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN9", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN7", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN5", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN3", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN1", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A6a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A5a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A4a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A3a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A2a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)
mic20_pig_cage_plot <- mic20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A1a", cage_mean_pigscore = mic20_pig_founder_mean, .before = 1)

view(mic20_pig_cage_plot)

# Plot
ggplot(mic20_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ylim(8,16) +
  ggtitle("Resident Microbe Additions (2020)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  scale_x_discrete(labels=c("July 15","September 7","November 9","November 24")) +
  scale_color_manual(name="Treatment", values=c("purple2","#FFC107","#1E88E5"), labels = c("Apple + At", "Apple", "Apple + Lb")) +
  #Adding average points
  geom_line(data=LM_mic20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_mic20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_mic20_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=mic20_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=mic20_pig_founder_mean-mic20_pig_founder_SE, ymax=mic20_pig_founder_mean+mic20_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)





#####
##### Resident Microbe Additions: just AT (2020)

## Read microbial additions 2020 dataset
mic20_pig_table <- read_excel("Microbes2020_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","2","4","5")
mic20_pig_table <- mic20_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
mic20_pig_table <- mic20_pig_table %>% mutate_at('cage', as.factor)
mic20_pig_table <- mic20_pig_table %>% mutate_at('treatment', as.factor)
mic20_pig_table$timepointID_cage <- paste(mic20_pig_table$timepoint_ID, "-", mic20_pig_table$cage)

view(mic20_pig_table)

#Filter out LB -- just comparing AT and control
micAT20_pig_table <- mic20_pig_table %>% filter(treatment!='LB')
view(micAT20_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

micAT20_pig_table_orchard <- micAT20_pig_table %>% filter(treatment != 'Founder')
view(micAT20_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

micAT20_pig_table_founder <- micAT20_pig_table %>% filter(treatment == 'Founder')
view(micAT20_pig_table_founder)

micAT20_pig_founder_mean <- mean(micAT20_pig_table_founder$pigmentation_score)
print(micAT20_pig_founder_mean)

micAT20_pig_founder_SE <- std.error(micAT20_pig_table_founder$pigmentation_score)
print(micAT20_pig_founder_SE)


#Calculating mean pigmentation score for each cage
micAT20_pig_cage_means <- micAT20_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(micAT20_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_micAT20_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = micAT20_pig_table_orchard)
anova(LM_micAT20_pig)
summary(LM_micAT20_pig)
confint(LM_micAT20_pig)
plot(LM_micAT20_pig)
qqnorm(resid(LM_micAT20_pig))
qqline(resid(LM_micAT20_pig))

# Estimated marginal means
emm_micAT20 <- emmeans(LM_micAT20_pig, ~ timepoint_ID*treatment)
print(emm_micAT20)
plot(emm_micAT20)

# Pairwise comparisons

# All comparisons
pairs(emm_micAT20)

# Make list of planned comparisons
micAT20_contrasts <- list(
  #control over time (end of summer to end of fall)
  E_T2.T5 = c(0, 0, 0, 1, 0, -1),
  #At over time (end of summer to end of fall)
  At_T2.T5 = c(1, 0, -1, 0, 0, 0),
  #TP2 control vs. At (end of summer)
  T2_E.At = c(-1, 0, 0, 1, 0, 0),
  #TP5 control vs. At (end of fall)
  T5_E.At = c(0, 0, -1, 0, 0, 1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_micAT20 = contrast(emm_micAT20, method = micAT20_contrasts, adjust = "none")
print(contrasts_micAT20)
confint(contrasts_micAT20)

# Test contrasts with a Holm correction
contrasts_holm_micAT20 = contrast(emm_micAT20, method = micAT20_contrasts, adjust = "holm")
print(contrasts_holm_micAT20)


## Extracting means and standard errors from micAT20 model for plotting
LM_micAT20_pig_means <- emmeans(LM_micAT20_pig, ~ timepoint_ID*treatment)
print(LM_micAT20_pig_means)
LM_micAT20_pig_means <- as.data.frame(LM_micAT20_pig_means)

view(LM_micAT20_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_micAT20_pig_plot <- subset(LM_micAT20_pig_means, select = -c(df, lower.CL, upper.CL))
LM_micAT20_pig_plot <- LM_micAT20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "AT", emmean = micAT20_pig_founder_mean, SE = micAT20_pig_founder_SE, .before = 1)
LM_micAT20_pig_plot <- LM_micAT20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "Control", emmean = micAT20_pig_founder_mean, SE = micAT20_pig_founder_SE, .before = 1)

view(LM_micAT20_pig_plot)

# Update individual cage dataframe to include founder
micAT20_pig_cage_plot <- subset(micAT20_pig_cage_means, select = -c(timepointID_cage))
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN12", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN10", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN8", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN6", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN4", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "AT", cage = "ATN2", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A6a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A5a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A4a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A3a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A2a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)
micAT20_pig_cage_plot <- micAT20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A1a", cage_mean_pigscore = micAT20_pig_founder_mean, .before = 1)

view(micAT20_pig_cage_plot)

# Plot
ggplot(micAT20_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Microbial Additions: AT (2020)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  ylim(8,16) +
  scale_x_discrete(labels=c("July 15","September 7","November 9","November 24")) +
  scale_color_manual(name="Treatment", values=c("purple2","#FFC107"), labels=c("Apple + At", "Apple")) +
  #Adding average points
  geom_line(data=LM_micAT20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_micAT20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_micAT20_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=micAT20_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=micAT20_pig_founder_mean-micAT20_pig_founder_SE, ymax=micAT20_pig_founder_mean+micAT20_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)





#####
##### Resident Microbe Additions: just LB (2020)

## Read microbial additions 2020 dataset
mic20_pig_table <- read_excel("Microbes2020_Pigmentation_R.xlsx")

pigmentation_levels <- c("0","2","4","5")
mic20_pig_table <- mic20_pig_table %>% mutate(timepoint_ID = factor(timepoint_ID, levels = pigmentation_levels, ordered = TRUE))
mic20_pig_table <- mic20_pig_table %>% mutate_at('cage', as.factor)
mic20_pig_table <- mic20_pig_table %>% mutate_at('treatment', as.factor)
mic20_pig_table$timepointID_cage <- paste(mic20_pig_table$timepoint_ID, "-", mic20_pig_table$cage)

view(mic20_pig_table)

#Filter out AT -- just comparing LB and control
micLB20_pig_table <- mic20_pig_table %>% filter(treatment!='AT')
view(micLB20_pig_table)

#Make data frame with TP1 and TP2 scores only (no Founder; for linear mixed effects model)

micLB20_pig_table_orchard <- micLB20_pig_table %>% filter(treatment != 'Founder')
view(micLB20_pig_table_orchard)

#Make data frame with Founder only (for plotting)
#Calculate raw mean and SE from founder data for adding founder point to plot

micLB20_pig_table_founder <- micLB20_pig_table %>% filter(treatment == 'Founder')
view(micLB20_pig_table_founder)

micLB20_pig_founder_mean <- mean(micLB20_pig_table_founder$pigmentation_score)
print(micLB20_pig_founder_mean)

micLB20_pig_founder_SE <- std.error(micLB20_pig_table_founder$pigmentation_score)
print(micLB20_pig_founder_SE)


#Calculating mean pigmentation score for each cage
micLB20_pig_cage_means <- micLB20_pig_table_orchard %>%
  group_by(timepointID_cage) %>%
  mutate(cage_mean_pigscore = mean(pigmentation_score, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(timepointID_cage, timepoint_ID, treatment, cage, cage_mean_pigscore)

view(micLB20_pig_cage_means)


## Pigmentation Statistics

# Linear mixed effects model (LMM)
LM_micLB20_pig <- lmer(pigmentation_score ~ timepoint_ID*treatment + (1|cage), data = micLB20_pig_table_orchard)
anova(LM_micLB20_pig)
summary(LM_micLB20_pig)
confint(LM_micLB20_pig)
plot(LM_micLB20_pig)
qqnorm(resid(LM_micLB20_pig))
qqline(resid(LM_micLB20_pig))

# Estimated marginal means
emm_micLB20 <- emmeans(LM_micLB20_pig, ~ timepoint_ID*treatment)
print(emm_micLB20)
plot(emm_micLB20)

# Pairwise comparisons

# All comparisons
pairs(emm_micLB20)

# Make list of planned comparisons
micLB20_contrasts <- list(
  #control over time (end of summer to end of fall)
  E_T2.T5 = c(1, 0, -1, 0, 0, 0),
  #Lb over time (end of summer to end of fall)
  Lb_T2.T5 = c(0, 0, 0, 1, 0, -1),
  #TP2 control vs. Lb (end of summer)
  T2_E.Lb = c(1, 0, 0, -1, 0, 0),
  #TP5 control vs. Lb (end of fall)
  T5_E.Lb = c(0, 0, 1, 0, 0, -1)
)

# Test contrasts with no p value correction
# Effect size (estimate) and 95% confidence interval 
# 95% CI = [lower.CL, upper.CL]

contrasts_micLB20 = contrast(emm_micLB20, method = micLB20_contrasts, adjust = "none")
print(contrasts_micLB20)
confint(contrasts_micLB20)

# Test contrasts with a Holm correction
contrasts_holm_micLB20 = contrast(emm_micLB20, method = micLB20_contrasts, adjust = "holm")
print(contrasts_holm_micLB20)


## Extracting means and standard errors from micLB20 model for plotting
LM_micLB20_pig_means <- emmeans(LM_micLB20_pig, ~ timepoint_ID*treatment)
print(LM_micLB20_pig_means)
LM_micLB20_pig_means <- as.data.frame(LM_micLB20_pig_means)

view(LM_micLB20_pig_means)


## Plot with founder connected
# Update linear model means and SE dataframe to include founder
LM_micLB20_pig_plot <- subset(LM_micLB20_pig_means, select = -c(df, lower.CL, upper.CL))
LM_micLB20_pig_plot <- LM_micLB20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "LB", emmean = micLB20_pig_founder_mean, SE = micLB20_pig_founder_SE, .before = 1)
LM_micLB20_pig_plot <- LM_micLB20_pig_plot %>% add_row(timepoint_ID = "0", treatment = "Control", emmean = micLB20_pig_founder_mean, SE = micLB20_pig_founder_SE, .before = 1)

view(LM_micLB20_pig_plot)

# Update individual cage dataframe to include founder
micLB20_pig_cage_plot <- subset(micLB20_pig_cage_means, select = -c(timepointID_cage))
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN11", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN9", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN7", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN5", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN3", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "LB", cage = "LBN1", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A6a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A5a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A4a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A3a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A2a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)
micLB20_pig_cage_plot <- micLB20_pig_cage_plot %>% add_row(timepoint_ID = "0", treatment = "Control", cage = "A1a", cage_mean_pigscore = micLB20_pig_founder_mean, .before = 1)

view(micLB20_pig_cage_plot)

# Plot
ggplot(micLB20_pig_cage_plot, aes(x=timepoint_ID, y=cage_mean_pigscore, group=cage, color=treatment)) +
  geom_line(alpha=0.25) +
  xlab("Timepoint") +
  ylab("Mean Pigmentation Score") +
  ggtitle("Microbial Additions: LB (2020)") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  labs(color='Treatment') +
  ylim(8,16) +
  scale_x_discrete(labels=c("July 15","September 7","November 9","November 24")) +
  scale_color_manual(name="Treatment", values=c("#FFC107","#1E88E5"), labels=c("Apple", "Apple + Lb")) +
  #Adding average points
  geom_line(data=LM_micLB20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), linewidth=2, inherit.aes = FALSE) +
  geom_point(data=LM_micLB20_pig_plot, aes(x=timepoint_ID, y=emmean, group=treatment, color=treatment), size=2.5, inherit.aes = FALSE) +
  #Adding standard error
  geom_errorbar(data=LM_micLB20_pig_plot, aes(x=timepoint_ID, ymin=emmean-SE, ymax=emmean+SE, group=treatment, color=treatment), linewidth=.45, width=0.05, inherit.aes = FALSE) +
  #Adding founder point and standard error based on raw means
  geom_point(aes(x=1, y=micLB20_pig_founder_mean), color="gray25", size=2.5) +
  geom_errorbar(aes(x=1, ymin=micLB20_pig_founder_mean-micLB20_pig_founder_SE, ymax=micLB20_pig_founder_mean+micLB20_pig_founder_SE), linewidth=.35, color="gray25", width=0.05)












