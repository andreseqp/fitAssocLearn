## Let´s try with cmdstanr -----------------------------------------------------
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

# Install cmdstan
set_cmdstan_path("M:\\Projects\\cmdstan")
#install_cmdstan()
cmdstan_path()
cmdstan_version()

cmdstanr::install_cmdstan(cores =4, overwrite = TRUE)

rebuild_cmdstan()


# Run the example

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)

data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)
