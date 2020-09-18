
library(rstan)
library(ggplot2)
library(plotly)
library(dplyr) ## for mutate/select
library(broom)
library(quantmod)
#install.packages("matrixStats")

library(matrixStats)
library(matrixStats)
library(bayesrules)
rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())

#getSymbols("^GSPC",from = "2007-01-01",to = "2020-09-03")
#saveRDS(GSPC,"gspc.rds")
GSPC <- readRDS("GSPC.rds")
SP500              <- GSPC$GSPC.Adjusted
SP500Ret           <- diff(log(SP500), lag = 1)
SP500Ret[is.na(SP500Ret)] <- 0
colnames(SP500Ret) <- c('LogReturn')


## Viz of log return
tidySP500Ret <- tidy(SP500Ret)
plot<- ggplot(tidySP500Ret, aes(x = index, y = value, color = series)) + 
  geom_line() + 
  theme_bw() +
  labs(title = "SP500 Returns Returns from 200 to 2020", x = "")
plot
ggplotly(plot, tooltip = c("city"))


unwantedObs <- SP500Ret["2007-01-01/2019-12-31", which.i=TRUE]

# Remove all data up to 2019 December
SP500Ret2020 <- SP500Ret[-unwantedObs,]


ggplot(tidy(SP500Ret2020), aes(x = index, y = value, color = series)) + 
  geom_line() + 
  theme_bw() +
  labs(title = "SP500 Log Returns from January to Sept 02,2020", x = "")



garchmodel <- "data {
  int<lower=0> T;  //total number of days
  int<lower=0> H; //number of days ahead forecast
  real y[T];  // log return at time T
  real<lower=0> sigma1; //3 sqrt volatiltiy
}
parameters {
  real mu; 
  real<lower=0> alpha0;          
  real<lower=0,upper=1> alpha1;  
  real<lower=0, upper=(1-alpha1)> beta1; 
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0 +
                      + alpha1 * square(y[t - 1] - mu)
                    + beta1  * square(sigma[t - 1]));
}
model {
  // prior
  alpha0 ~ normal(2.47e-6,0.5);
  alpha1 ~ normal(0.1256,0.5);
  beta1 ~ normal(0.8543,0.5);
  // likelihood
  y ~ normal(mu,sigma);
}
generated quantities {
  real<lower=0> sigma_pred[T];
  real<lower=0> sigma_bfore[H];
  real y_bfore[H];
  real y_pred[T];

  for (h in 1:H) {
    sigma_bfore[h] = sqrt(
      alpha0
      + alpha1 * pow(y[T + h - 1], 2)
      + beta1 * pow(sigma[T + h - 1], 2)
      );
    y_bfore[h] = normal_rng(0, sigma_bfore[h]);
  }
  
  sigma_pred[1] = sigma1;
  y_pred[1] = normal_rng(0, sigma1);

  for (t in 2:T) {
    sigma_pred[t] = sqrt(
    alpha0
    + alpha1 * pow(y_pred[t-1], 2)
    + beta1 * pow(sigma_pred[t-1], 2)
    );
    y_pred[t] = normal_rng(0, sigma_pred[t]);
}
}"

## priors : https://www.mdpi.com/1911-8074/13/3/51
SP500Ret[is.na(SP500Ret)] <- 0


y <-  list(y = as.numeric(SP500Ret2020), T = length(SP500Ret2020),H=1, 
           sigma1 = sd(SP500Ret2020))
SP500_model <- stan(model_code = garchmodel, data = y,
                    chains = 4,
                    iter = 1250,
                    warmup = 250,
)

## Posterior traceplot and distribution


params <- rstan::extract(SP500_model)
bayesplot::mcmc_trace(SP500_model,pars = c("alpha0","alpha1","beta1"))
par(mfrow=c(1,3))
hist(params$alpha0,xlab = "Posterior alpha0",main = paste("Histogram of alpha0"))
hist(params$alpha1,xlab = "Posterior alpha1",main = paste("Histogram of alpha1"))
hist(params$beta1,xlab = "Posterior beta1",main = paste("Histogram of beta1"))


## Posterior summary statistitcs

mymat<- summary(SP500_model,pars = c("alpha0","alpha1","beta1"))$ summary

posterior_summary <- data.frame(c(as.vector(summary(SP500_model,pars = c("alpha0"))$ summary)),
                                  c(as.vector(summary(SP500_model,pars = c("alpha1"))$ summary)),
                                  c(as.vector(summary(SP500_model,pars = c("beta1"))$ summary)),
                                row.names = colnames(mymat))

colnames(posterior_summary) <- c("alpha0","alpha1","beta1")
posterior_summary






# Why are prior and posteriors different?
params <- rstan::extract(SP500_model)
mu = median(params$mu)
alpha0 = median(params$alpha0)
alpha1 = median(params$alpha1)
beta1 <- median(params$beta1)

Priors <- c(2.47e-6,0.1256,0.8543)
Posteriors <-c(alpha0,alpha1,beta1)


priorVPosterior <- data.frame(Priors, Posteriors,row.names = c("alpha0","alpha1","beta1"))
priorVPosterior


## credible interval and predictions
spred <- extract(SP500_model, pars = "sigma_pred")[[1]]
ypred <- extract(SP500_model, pars = "y_pred")[[1]]
sfore <- extract(SP500_model, pars = "sigma_bfore")[[1]]
yfore <- extract(SP500_model, pars = "y_bfore")[[1]]



# median volatility and logreturn values from model
medVolatility <- median(colMedians(spred))
medLogRet <- median(colMedians(ypred))


# Median volatility, logreturn forecast  
nextDayvolatility <- median(sfore)
nextDayLogRet <- median(yfore)



# Error analysis
mypred <- colMedians(ypred)
mytrue <- as.vector(SP500Ret2020)

mae <- median(abs(mypred-mytrue))
mae


