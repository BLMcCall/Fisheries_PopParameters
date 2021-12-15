#***********************************************************************
#
#                        PBD manuscript analyses
#                             10 Feb 2021
#
#
#***********************************************************************

### Relative abundance and sex percentages datasets ----

## Importing dataset from excel
library(readxl)

# To plot data as a barplot with multiple variables, has to be in the format 
# "X-axis | Condition | Y-axis"
#   If more than three variables, then will need to melt data with [melt]
ra <- read_excel("F:/Publications/SFCProceedings_Epallididorsum_Survey/
                 RA_SA_TropFishR.xlsx", sheet = "RA_Drainage")

sex <- read_excel("F:/Publications/SFCProceedings_Epallididorsum_Survey/
                  RA_SA_TropFishR.xlsx", sheet = "Sex")

## Plotting relative abundance and sex ratio
library(ggplot2)
library(tidyverse)
library(ggsci) #color palletes

# Need to reorder the x-axis so seasons are in order; convert x-axis to factors and arrange the
# order how I like using 'levels'
ra$year <- factor(ra$year, 
                  levels = c("2016W", "2016S", "2016Su", "2016F", "2017W", "2017S", "2017Su", "2017F"))

sex$year <- factor(sex$year, 
                  levels = c("2016W", "2016S", "2016Su", "2016F", "2017W", "2017S", "2017Su", "2017F"))

ggplot(data = ra, aes(x = year, y = value, fill = drainage)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_grey()  + # color theme
  theme_bw() + # setting the background and axis lines
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "Seasonal Year", y = "Relative Abundance (CPUE)", fill = "Drainages")

ggplot(data = sex, aes(x = year, y = value, fill = sex)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_grey()  + # color theme
  theme_bw() + # setting the background and axis lines
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(x = "Seasonal Year", y = "Sex Ratio (%)", fill = "Sex")

### Fish growth data arrangement ----
## Following https://cran.r-project.org/web/packages/TropFishR/vignettes/tutorial.html
## AND https://www.fishbase.in/manual/key%20facts.htm 
library(TropFishR)
library(devtools)

## Import LFQ dataset and arrange data for analysis format
lfq <- read.csv(choose.files())

caddo_lfq <- read.csv(choose.files())

ouachita_lfq <- read.csv(choose.files())

# Arranging the data for analyses following
#  https://cran.r-project.org/web/packages/TropFishR/vignettes/lfqData.html

##********************
# One population \ species-level
##********************
dates <- colnames(lfq)[-1]
dates <- strsplit(dates, "X")
dates <- unlist(lapply(dates, function(x) x[2]))
dates <- as.Date(dates, "%Y.%d.%m")

lfq.2 <- list(dates = dates,
              midLengths = lfq$lengthClass,
              catch = as.matrix(lfq[,-1]))

class(lfq.2) <- "lfq"  # Matrix with a dataframe as 'lfq' required for 
#                       [TropFishR]
lfq.2$catch[is.na(lfq.2$catch)] <- 0
plot(lfq.2, Fname="catch", ylab = "Length classes (mm)", date.axis = "modern")

## Raw data v reconstructed data
# Set seed for reproducible results
set.seed(123)

# Adjust bin size (need to do research to understand what bin size is)
## Wang et al. 2020 Selecting optimal bin size to account for growth variability in ELEFAN
## DOI: 10.1016/j.fishres.2019.105474
## Wang et al. 2020 developed a rule of thumb to calculate bin size, OBS = 1.86 x(Lmax/A)^0.45
lfq.3 <- lfqModify(lfq.2, bin_size = 8)

# Comparing raw and reconstructed lfq data
lfqbin <- lfqRestructure(lfq.3, MA = 5, addl.sqrt = FALSE)

plot(lfqbin, Fname = "rcounts", date.axis = "modern")


##*********************
# Caddo population
##*********************
dates <- colnames(caddo_lfq)[-1]
dates <- strsplit(dates, "X")
dates <- unlist(lapply(dates, function(x) x[2]))
dates <- as.Date(dates, "%Y.%d.%m")

lfq.c2 <- list(dates = dates,
              midLengths = caddo_lfq$lengthClass,
              catch = as.matrix(caddo_lfq[,-1]))

class(lfq.c2) <- "lfq"  # Matrix with a dataframe as 'lfq' required for 
#                       [TropFishR]
lfq.c2$catch[is.na(lfq.c2$catch)] <- 0
plot(lfq.c2, Fname="catch", ylab = "Length classes (mm)", date.axis = "modern")

## Raw data v reconstructed data
# Set seed for reproducible results
set.seed(123)

# Adjust bin size
lfq.c3 <- lfqModify(lfq.c2, bin_size = 8) 


##*******************
# Ouachita population
##*******************
dates <- colnames(ouachita_lfq)[-1]
dates <- strsplit(dates, "X")
dates <- unlist(lapply(dates, function(x) x[2]))
dates <- as.Date(dates, "%Y.%d.%m")

lfq.o2 <- list(dates = dates,
               midLengths = ouachita_lfq$lengthClass,
               catch = as.matrix(ouachita_lfq[,-1]))

class(lfq.o2) <- "lfq"  # Matrix with a dataframe as 'lfq' required for 
#                       [TropFishR]
lfq.o2$catch[is.na(lfq.o2$catch)] <- 0
plot(lfq.o2, Fname="catch", ylab = "Length classes (mm)", date.axis = "modern")


## Raw data v reconstructed data
# Set seed for reproducible results
set.seed(123)

# Adjust bin size
lfq.o3 <- lfqModify(lfq.o2, bin_size = 8) 


### Powell-Wetherall method ----

pw <- powell_wetherall(param = lfq.3,
                       catch_columns = 1:ncol(lfq.3$catch),
                       reg_int = NULL)
paste("Linf =", round(pw$Linf_est), "±", round(pw$se_Linf)) # 46 ± 4


##******************
# Caddo population
##******************
caddo_pw <- powell_wetherall(param = lfq.c3,
                       catch_columns = 1:ncol(lfq.c3$catch),
                       reg_int = NULL)
paste("Linf =", round(caddo_pw$Linf_est), "±", round(caddo_pw$se_Linf)) # 49 ± 9

##*****************
# Ouachita population
##*****************
ouachita_pw <- powell_wetherall(param = lfq.o3,
                       catch_columns = 1:ncol(lfq.o3$catch),
                       reg_int = NULL)
paste("Linf =", round(ouachita_pw$Linf_est), "±", round(ouachita_pw$se_Linf)) # 46 ± 5

# Getting a quick K using the pre-determimned Linf value
# ELEFAN = Electronical Length Frequency Analysis
kscan <- ELEFAN(lfq.3, Linf_fix = pw$Linf_est, MA = 5, 
                addl.sqrt = TRUE, hide.progressbar = FALSE)

kscan$par; kscan$Rn_max # Linf 45.64 (CI 33.31 - 57.98)| K 0.74 (CI 0.66 - 0.84)| 
#                         t_anchor 0.83 | C 0 | ts 0 | agemax 3 | ncohort 5 | Rn 0.48

### Simulated annealing and genetic algorithm ----

# Previous analysis cannot test for different combinations of Linf and K
# but there other methods that can like the genetic algorithm and 
# a simulated annealing. The optimized analyses disrupt the assumption of
# constant growth and by allowing flexibility to the parameters of
# "ts" (temperature | describes probability of accepting worse conditions) 
# and "C" (w/in year growth rate change).

# Linf values will be the largest size class (51-mm) +/- 10-mm for upper
# and lower limit. Initial K-values range from 0.5-1.0 with lower and upper
# limits (0.01-2.0). "t anchor" anchors growth curves corresponding to peak
# spawning month (October in our data, but may spawn late summer). "t anchor" 
# is a fraction of the year, so to get Oct 1st, 273 days into the year 
# (273/365 = 0.75) so "t anchor" is 0.75.

##***************************
# Simulated annealing | one population
##***************************

#  Five minute
sa_5m <- ELEFAN_SA(lfq.3, SA_time = 60*5, SA_temp = 6e5, MA = 5, 
                seasonalised = TRUE, addl.sqrt = FALSE,
                init_par = list(Linf = 51, K = 0.5, t_anchor = 0.75, C = 0.5, ts = 0.5),
                low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)

#ncohort 5 | agemax 2 | Linf 41.61 | K 0.90 | t_anchor 0.46 | C 0.85 | ts 0.67 | phil 3.19 | Rn 0.83

##***********************
# Simulated annealing | two populations
##***********************
# Caddo
# 5 minutes
caddo_sa <- ELEFAN_SA(lfq.c3, SA_time = 60*5, SA_temp = 6e5, MA = 5, 
                 seasonalised = TRUE, addl.sqrt = FALSE,
                 init_par = list(Linf = 51, K = 0.5, t_anchor = 0.75, C = 0.5, ts = 0.5),
                 low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                 up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)
#ncohort 5 | agemax 2 | Linf 53.02 | K 0.87 | t_anchor 0.02 | C 0.97 | ts 0.31 | phil 3.39 | Rn 0.73




# Ouachita
# 5 minute
ouachita_sa <- ELEFAN_SA(lfq.o3, SA_time = 60*5, SA_temp = 6e5, MA = 5, 
                          seasonalised = TRUE, addl.sqrt = FALSE,
                          init_par = list(Linf = 51, K = 0.5, t_anchor = 0.75, C = 0.5, ts = 0.5),
                          low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                          up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)

#ncohort 5 | agemax 2 | Linf 55.51 | K 0.85 | t_anchor 0.77 | C 0.92 | ts 0.51 | phil 3.42 | Rn 0.74


# Can use jack knife mehods to construct CI with a for loop of the [ELEFAN]
# Change the dataset name to repeat for each population
jk <- vector("list", length(lfq.3$dates))
for(i in 1:length(lfq.3$dates)){
  loop_data <- list(dates = lfq.3$dates[-i],
                    midLengths = lfq.3$midLengths,
                    catch = lfq.3$catch[,-i])
  tmp <- ELEFAN_SA(loop_data, SA_time = 60*5, SA_temp = 6e5, MA = 5, 
                   seasonalised = TRUE, addl.sqrt = FALSE, agemax = 2,
                   init_par = list(Linf = 51, K = 0.5, t_anchor = 0.75, C = 0.5, ts = 0.5),
                   low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                   up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1),
                   plot = FALSE)
  jk[[i]] <- unlist(c(tmp$par, list(Rn_max = tmp$Rn_max)))
}
jkRes <- do.call(cbind, jk)

jk_mean <- apply(as.matrix(jkRes), MARGIN = 1, FUN = mean)         # mean
jk_CI <- apply(as.matrix(jkRes), MARGIN = 1, FUN = function(x) 
  quantile(x, probs = c(0.025, 0.975)))                            # CI
jk_CI <- t(jk_CI)                                                  # transposing matrix
colnames(jk_CI) <- c("lower", "upper")


##***************************
# Genetic algoritihm | one population
##***************************
# 450 iterations
ga <- ELEFAN_GA(lfq.3, MA = 5, seasonalised = TRUE, maxiter = 450, addl.sqrt = FALSE,
                 low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                 up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)
#ncohort 5 | agemax 2 | Linf 54.35 | K 0.69 | t_anchor 0.56 | C 0.46 | ts 0.31 | phil 3.31 | Rn 0.67


##***********************
# Genetic algorithm | Two populations
##***********************
# Caddo
# 450 iterations
caddo_ga <- ELEFAN_GA(lfq.c3, MA = 5, seasonalised = TRUE, maxiter = 450, addl.sqrt = FALSE,
                low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)
#ncohort 4 | agemax 2 | Linf 55.99 | K 0.79 | t_anchor 0.83 | C 0.47 | ts 0.32 | phil 3.39 | Rn 0.73


# Ouachita
# 850 iterations
ouachita_ga <- ELEFAN_GA(lfq.o3, MA = 5, seasonalised = TRUE, maxiter = 850, addl.sqrt = FALSE,
                      low_par = list(Linf = 41, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = list(Linf = 61, K = 1, t_anchor = 1, C = 1, ts = 1), agemax = 2)
#ncohort 5 | agemax 2 | Linf 55.61 | K 0.87 | t_anchor 0.76 | C 0.85 | ts 0.57 | phil 3.43 | Rn 0.74


### Plotting growth curves ----

# Seasonality was included in growth curve plots because seasonality did seem to play a role with 
# a high 'C' value for all.

##***************************
# One population
##***************************
# Simmulated annealing results
plot(lfqbin, Fname = "rcounts", date.axis = "modern", ylim = c(0,60))
lfqFitCurves(lfq.3, par = sa_5m$par, draw = TRUE, col = "black", lty = 1, lwd = 1, agemax = 2)

##***************************
# Two populations
##***************************
# Caddo
plot(caddo_lfqbin, Fname = "rcounts", date.axis = "modern", ylim = c(0,60))
lfqFitCurves(lfq.c3, par = caddo_sa$par, draw = TRUE, col = "black", lty = 1, lwd = 1, agemax = 2)

# Ouachita
plot(ouachita_lfqbin, Fname = "rcounts", date.axis = "modern", ylim = c(0,60))
lfqFitCurves(lfq.o3, par = ouachita_sa$par, draw = TRUE, col = "black", lty = 1, lwd = 1, agemax = 2)


### Natural mortality estimates ----

## Going to move forward with the simulated annealing algorithm results since they are best fit (high Rn values)
# Assign sa algorithm estimates to the data list
lfq_sa <- c(lfq.3, sa_5m$par)
class(lfq_sa) <- "lfq"

## Caddo River
lfq_sa_caddo <- c(lfq.c3, caddo_sa$par)
class(lfq_sa_caddo) <- "lfq"

## Ouachita River
lfq_sa_ouachita <- c(lfq.o3, ouachita_sa$par)
class(lfq_sa_ouachita) <- "lfq"

## Estimates of natural mortality, life span, and maturity estimates
# When no controlled experiments or tagged data are available, the main approach for
# this estimation is to use empirical formulas that can be found with ?M_empirical 
# and in the detail section describing 'method'

# I only have data representing length (mm) and Linf and K estimates of length from previous
# analyses; therefore, "Then_growth" is most appropriate
nm <- M_empirical(Linf = sa_5m$par$Linf, K_l = sa_5m$par$K, method = "Then_growth")
lfq_sa$M <- as.numeric(nm)  # natural mortality of 1.11 year^-1

## Caddo River
nm_caddo <- M_empirical(Linf = caddo_sa$par$Linf, K_l = caddo_sa$par$K,
                        method = "Then_growth")
lfq_sa_caddo$M <- as.numeric(nm_caddo) # natural mortality of 1.01 year^-1


# Ouachita River
nm_ouachita <- M_empirical(Linf = ouachita_sa$par$Linf, K_l = ouachita_sa$par$K,
                           method = "Then_growth")
lfq_sa_ouachita$M <- as.numeric(nm_ouachita) # natural mortality of 0.98 year^-1

# Life span
# Following a commonly cited equation from Taylor (1958) tmax = t0 + 2.996/k
lspan <- 0 + (2.996/0.60) # This equals 4.9 (~5), which is the same estimate provided
#                           in the ga elefan for agemax. T0 is set to zero b/c cannot
#                           derive t0 from lfq data alone and is a constant in the 
#                           statistic with no biological significance. Statistic also
#                           referred to as A_0.95

