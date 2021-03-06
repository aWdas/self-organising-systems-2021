# Experiment 1: Inertia effect on dispersion
install.packages("tidyverse")
library(dplyr)
library(ggplot2)

install.packages("caret")
library(caret)


rms <- function(X, Y) {
  x_mean <- mean(X)
  y_mean <- mean(Y)
  n <- length(X)
  
  return(sqrt(1/n * sum((X - x_mean)^2 + (Y - y_mean)^2)))
}

## Schaffer
schaffer <- read.csv("netlogo_inertia_penalty_constraints_schaffer_coordinates.csv", 
                     header = FALSE, col.names = c("fitness_function", "run_nr", "inertia", "particle_nr", "x", "y"))

dispersion_schaffer <- schaffer %>% group_by(fitness_function, inertia, run_nr) %>% summarize(rms = rms(x,y)) 

lm_schaffer <- lm(rms ~ inertia, dispersion_schaffer)
summary(lm_schaffer)

## Shubert
shubert <- read.csv("netlogo_inertia_penalty_constraints_shubert_coordinates.csv", 
                     header = FALSE, col.names = c("fitness_function", "run_nr", "inertia", "particle_nr", "x", "y"))

dispersion_shubert <- shubert %>% group_by(fitness_function, inertia, run_nr) %>% summarize(rms = rms(x,y)) 

lm_shubert <- lm(rms ~ inertia, dispersion_shubert)
summary(lm_shubert)

## Eggholder
eggholder <- read.csv("netlogo_inertia_penalty_constraints_eggholder_coordinates.csv", 
                     header = FALSE, col.names = c("fitness_function", "run_nr", "inertia", "particle_nr", "x", "y"))

dispersion_eggholder <- eggholder %>% group_by(fitness_function, inertia, run_nr) %>% summarize(rms = rms(x,y)) 

lm_eggholder <- lm(rms ~ inertia, dispersion_eggholder)
summary(lm_eggholder)

## Together
experiment1_all <- bind_rows(dispersion_schaffer, dispersion_shubert, dispersion_eggholder)
experiment1_all$fitness_function <- factor(experiment1_all$fitness_function, 
                                          levels=c("Fitness function Schaffer", "Fitness function Shubert", "Fitness function Eggholder"),
                                          labels=c("Schaffer", "Shubert", "Eggholder"))

lm_all <- lm(rms ~ inertia + fitness_function + inertia:fitness_function, experiment1_all)
summary(lm_all)

ggplot(experiment1_all, aes(x=inertia, y=rms)) + 
  geom_point(aes(colour = fitness_function))


# Experiment 2: Effect of population-size on global-best-val

experiment2 <- read.csv("netlogo_population_size_best_vals.csv", 
                        header = FALSE, col.names = c("run_nr", "fitness_function", "population_size", "global_best_val"))

experiment2$fitness_function <- factor(experiment2$fitness_function, 
                                          levels=c("Fitness function Schaffer", "Fitness function Shubert", "Fitness function Eggholder"),
                                          labels=c("Schaffer", "Shubert", "Eggholder"))

lm_experiment2 <- lm(global_best_val ~ population_size, experiment2)
summary(lm_experiment2)

ggplot(experiment2, aes(x=population_size, y=global_best_val)) + 
  geom_point(aes(colour = fitness_function), position = position_jitter(width = 1, height = 0.02))
