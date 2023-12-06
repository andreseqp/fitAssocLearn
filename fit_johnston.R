## Fitting a RW model to the data of Johnston 2023 -------------------------
# Johnston, M., K. F. Brecht, and A. Nieder. 2023. Crows flexibly apply 
# statistical inferences based on previous experience. 
# Current Biology 33:3238-3243.e3.


## Import important packages
library(data.table)
library(readxl)
library(here)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(shinystan)
library(stringi)
library(tidyverse)
theme_set(new = theme_classic())

## Read the data 
full_train <- read_excel(here("data","Johnston_2023","Data.xlsx"),
                             sheet = 'Initial training')
full_test <- read_excel(here('data','Johnston_2023','Data.xlsx'),
                        sheet = 'Testing')

full_train <- full_train %>% select(c(crow,session,trial,condition,reward))
full_test  <- full_test %>% select(c(crow,session,trial,condition,left_p,right_p,
                                     reward,choice))

## Need to write the stan model that includes the training session