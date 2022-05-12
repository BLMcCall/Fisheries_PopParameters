#*****************************************
#
#     CPUE comparison among drainage
#             06 April 2022
#
#
#*****************************************


# CPUE averaged for each season per drainage ----
C_CPUE <- c(0.43, 0.48, 1.13, 1.17, 0.20, 0.62, 1.08, 0.67)
O_CPUE <- c(0.40, 0.30, 0.90, 0.46, 0.14, 0.51, 0.88, 0.41)

# Are datasets normally distributed?
shapiro.test(C_CPUE)   # W = 0.91 | p = 0.34
shapiro.test(O_CPUE)   # W = 0.90 | p = 0.27

# Performing a paired t-test
t.test(C_CPUE, O_CPUE, paired = TRUE)   
# Results: t = 2.95 | df = 7 | p = 0.02 *

# Visualize difference, if informative consider incorporating plots
if(!require(devtools)) install.packages("devtools")
install.packages("ggpubr")

CPUE <- data.frame(Rivers = rep(c("Caddo", "Ouachita"), each = 8),
                   weight = c(C_CPUE, O_CPUE))

library(ggpubr)
ggboxplot(CPUE, x = "Rivers", y = "weight",
          color = "Rivers", palette = c("#B10026", "#084594"),
          order = c("Caddo", "Ouachita"), add = "jitter",
          ylab = "CPUE", xlab = "River Drainage",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.15,
          legend = "none")

# CPUE compared across seasons ----
# building dataframe
Winter <- c(0.43, 0.40, 0.20, 0.14)
Spring <- c(0.48, 0.30, 0.62, 0.51)
Summer <- c(1.13, 0.90, 1.08, 0.88)
Fall <- c(1.17, 0.14, 0.67, 0.41)

seasonal <- data.frame(id = 1:4)
seasonal$winter <- Winter
seasonal$spring <- Spring
seasonal$summer <- Summer
seasonal$fall <- Fall

library(reshape2)
seasonal_long <- melt(seasonal,
                      id.vars = "id",
                      measure.vars = c("winter", "spring", "summer", "fall"),
                      variable.name = "seasons",
                      value.name = "cpue")

## Normally distributed?
shapiro.test(seasonal$CPUE) # W = 0.92 | p =0.19

# ANOVA
anova <- aov(cpue~seasons, data = seasonal_long)
summary(anova)
# Results: F = 5.78 | df = 3,12 | p = 0.01 *

## Post-hoc tests
TukeyHSD(aov(cpue~seasons, data = seasonal_long))
# Results: 

# CPUE for the sexes averaged for each season per drainage ----
## FEMALES
FC_CPUE <- c(0.27, 0.20, 0.87, 0.61, 0.14, 0.26, 0.68, 0.23)
FO_CPUE <- c(0.15, 0.17, 0.55, 0.26, 0.10, 0.22, 0.50, 0.16)

# Normally distributed?
shapiro.test(FC_CPUE) # W = 0.85 | p = 0.09
shapiro.test(FO_CPUE) # W = 0.82 | p = 0.04 *

# Nonparametric paired t-test
wilcox.test(FC_CPUE, FO_CPUE, paired = TRUE)
# Results: V = 36 | p = 0.01 *

# Boxplots
Female_CPUE <- data.frame(Rivers = rep(c("Caddo", "Ouachita"), each = 8),
                   weight = c(FC_CPUE, FO_CPUE))

ggboxplot(Female_CPUE, x = "Rivers", y = "weight",
          color = "Rivers", palette = c("#B10026", "#084594"),
          order = c("Caddo", "Ouachita"), add = "jitter",
          ylab = "CPUE", xlab = "River Drainage",
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.15,
          legend = "none")

## MALES
MC_CPUE <- c(0.10, 0.07, 0.24, 0.31, 0.06, 0.13, 0.39, 0.09)
MO_CPUE <- c(0.50, 0.12, 0.38, 0.04, 0.07, 0.14, 0.38, 0.12)

# Normally distributed?
shapiro.test(MC_CPUE) # W = 0.85 | p = 0.12
shapiro.test(MO_CPUE) # W = 0.85 | p = 0.09

# Paired t-test
t.test(MC_CPUE, MO_CPUE, paired = TRUE)
# Results: t = -0.69 | df = 7 | p = 0.51

## JUVENILES
JC_CPUE <- c(0.05, 0.43, 0.03, 0.05, 0.00, 0.23, 0.00, 0.36)
JO_CPUE <- c(0.00, 0.02, 0.00, 0.30, 0.00, 0.16, 0.00, 0.10)

# Normally distributed?
shapiro.test(JC_CPUE) # W = 0.80 | p = 0.03 *
shapiro.test(JO_CPUE) # W = 0.75 | p = 0.01 *

# Nonparametric paired t-test
wilcox.test(JC_CPUE, JO_CPUE, paired = TRUE)
# Results: V = 17 | p = 0.21

# SL comparison across drainages between seasons ----
SLC <- c(46, 40, 25, 37, 44, 44, 49, 42)
OLC <- c(35, 40, 28, 22, 41, 38, 25, 41)

# Normally distributed?
shapiro.test(SLC) # W = 0.87 | p = 0.14
shapiro.test(OLC) # W = 0.86 | p = 0.12

# Paired t-test
t.test(SLC, OLC, paired = TRUE)
# Results: t = 2.23 | df = 7 | p = 0.06