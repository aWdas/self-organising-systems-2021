install.packages("tidyverse")
library(dplyr)
library(ggplot2)

####### Experiment 3 ##############
experiment3 = read.csv("netlogo_inertia_convergence_best_vals.csv",
    header = FALSE, col.names = c("run_nr", "fitness_function", "inertia","iterations", "global_best_val"))
experiment3$fitness_function <- factor(experiment3$fitness_function, 
                                       levels=c("Fitness function Schaffer", "Fitness function Shubert", "Fitness function Eggholder"),
                                       labels=c("Schaffer", "Shubert", "Eggholder"))
#Plotting iterations before optimum
ggplot(experiment3, aes(x=inertia, y=iterations)) + 
  geom_point(aes(colour = fitness_function), position = position_jitter(width = 0.005, height = 2)) +
  stat_summary(size=1, aes(y = iterations,color=fitness_function), fun=mean, geom="line") +
  labs (x = 'Particle Inertia', y = 'Iterations before optimum')

#Inspecting cases that did not converge
nc = filter(experiment3,experiment3$iterations==100)
ggplot(nc, aes(x=inertia, y=global_best_val)) +
  geom_point(aes(colour = fitness_function), position = position_jitter(width = 0.005, height = 0.005)) +
  labs( x = "Partical Inertia", y = "Best value found before timeout")

####### Experiment 4 ##############
experiment4 = read.csv("netlogo_speed_limit_best_vals.csv",
                  header = FALSE, col.names = c("run_nr", "fitness_function", "maxspeed","iterations", "global_best_val"))
experiment4$fitness_function <- factor(experiment4$fitness_function, 
                                  levels=c("Fitness function Langermann", "Fitness function Schwefel", "Fitness function Shubert", "Fitness function Schaffer", "Fitness function Eggholder", "Fitness function Easom", "Fitness function Booth"),
                                  labels=c("Langermann", "Schwefel","Shubert", "Schaffer", "Eggholder", "Easom", "Booth"))
#Plotting iterations before optimum
ggplot(experiment4, aes(x=maxspeed, y=iterations)) + 
#  geom_point(aes(colour = fitness_function), position = position_jitter(width = 0.1, height = 1)) +
  stat_summary(size=1, aes(y = iterations,color=fitness_function), fun=median, geom="line") +
  stat_summary(size=1, aes(y = iterations), fun=median, geom="line") +
    labs (x = 'Maximum Speed', y = 'Iterations before optimum', colour = "Fitness Function")

#Plotting as deviation from total mean
ff = aggregate(experiment4$iterations,list(experiment4$maxspeed,experiment4$fitness_function),mean)
total = aggregate(experiment4$iterations,list(experiment4$maxspeed),mean)
ff$x = ff$x-total$x
ggplot(ff, aes(x=Group.1, y=x)) +
  geom_hline(yintercept=0)+
  stat_summary(size=1, aes(y = x,color=Group.2), fun=median, geom="line") +
  labs(x = "Maximum Speed", y = "Deviation from total mean", colour="Fitness Function")
  