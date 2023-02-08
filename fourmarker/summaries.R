###############################
###############################
## Program to read in AFT model,
## extract measures of interest,
## and compute summary
## statistics from model fit4.
###############################
###############################

###############################
## Load libraries
###############################
library(tidyverse)
library(rstan)

###############################
## Load the simulated data
###############################
simdata <- readRDS("simdata.rds")

# Recreate the list object passed to stan
# Recall
# id: 1, 2, ... with no holes
# outcome: 1, 2, 3, or 4
# center education at 16
uid <- sort(unique(simdata$id))
if (any(simdata$id != match(simdata$id, uid))) stop("id is not sequential")

stan.list <- list(id = simdata$id, y= simdata$y, 
               outcome= as.numeric(simdata$outcome),
               age = simdata$age,
               x = cbind(simdata$male, simdata$apoepos, simdata$apoemiss,
                         simdata$educ-16),
               adrc = simdata$adrc
)

## df for the t distribution defined in stanfit.R
dft <- 10

###############################
## Load the AFT model
###############################

## This line can be slow.
## use this trick so model is not re-loaded
if (!exists("fit4")) load("fit4save.rda", verbose=T)

## 's4' = summary of model on 4 outcomes
if (!exists("s4")) s4 <- summary(fit4) ## can be slow to execute

## used to extract certain parameters by name
sname <- rownames(s4[[1]])

###############################
## Review log-liklihoods for each chain
###############################
## More thorough checks that the model has converged are recommended and
## are documented at
## https://mc-stan.org/docs/2_25/reference-manual/convergence.html

## Get the log-liklihood
round(s4[[2]]["lp__", c(3,1,7),])

## visualize the log-likelihood to assess how much they overlapped or not
s4[[2]]["lp__", c(3,1,7),] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var='chains') %>%
ggplot +
    geom_pointrange(aes(y=chains, x=`mean`, xmin=`2.5%`, xmax=`97.5%`),
                    size=1, color="blue", fill="white", shape=22)
## Does it appears though the chains have converged to the same place?




#################################
## Review per-chain random intercepts, residual error
#################################
atau <- s4[[2]][grepl("atau", sname),,]
sigma <- s4[[2]][grepl("sigma", sname),,]

round(atau[,"mean",],1)

round(sigma[,"mean",],3)

## We want to see that the random intercepts (alpha) for all
## chains are essentially identical.
## The std of the random intercept (atau) for all chains are
## essentially identical.
## The residual standard error (sigma) for all chains are
## essentially identical (sigma can be robust to differences
## in atau, if they had existed)




###############################
###############################
## Process stan output to get summary statistics
###############################
###############################

###############################
## Some helpers
###############################
pctls <- c("mean", "2.5%", "97.5%")

## Another package seems to be overwriting the extract function. fix that.
extract <- function(...) rstan::extract(...)

myci <- function(x) {
    setNames(c(mean(x), quantile(x, c(0.025, 0.975))), c("mean", "lcl", "ucl"))
}
formatCI <- function (x, digits = 2, format = "f", ci.sep = ", ", ...) {
    Est <- formatC(x[1], digits = digits, format = format)
    LCL <- formatC(x[2], digits = digits, format = format)
    UCL <- formatC(x[3], digits = digits, format = format)
    paste(Est, " (", LCL, ci.sep, UCL, ")", sep = "")
}



###############################
###############################
## Covariate effects
###############################
###############################

## Q: How do covariates affect the onset of accumulation?

## A: The covariate effects report the estimated acceleration.
## For example, the estimate of 6.1 for APOE 𝜀4+
## implies an amyloid accumulation that begins, on average,
## 6.1 years earlier in APOE 𝜀4 carriers than non-carriers.


## -------------------------
## Male sex effect
## =========================
formatCI(s4[[1]]["beta[1,1]", pctls], digits = 1) ## Male sex effect on pib
formatCI(s4[[1]]["beta[1,2]", pctls], digits = 1) ## Male sex effect on tau
formatCI(s4[[1]]["beta[1,3]", pctls], digits = 1) ## Male sex effect on wmh
formatCI(s4[[1]]["beta[1,4]", pctls], digits = 1) ## Male sex effect on fa

## Posterior samples
male.mcmc <- extract(fit4, c("beta[1,1]", "beta[1,2]", "beta[1,3]", "beta[1,4]"))

## The coefficients are in terms age adjustments. 
(formatCI(myci(male.mcmc[[1]]), digits = 1))
(formatCI(myci(male.mcmc[[2]]), digits = 1))
(formatCI(myci(male.mcmc[[3]]), digits = 1))
(formatCI(myci(male.mcmc[[4]]), digits = 1))

## sign switch to get female sex effect
(formatCI(myci(-male.mcmc[[1]]), digits = 1))
(formatCI(myci(-male.mcmc[[2]]), digits = 1))
(formatCI(myci(-male.mcmc[[3]]), digits = 1))
(formatCI(myci(-male.mcmc[[4]]), digits = 1))


## -------------------------
## APOE effects
## =========================

## e4+
formatCI(s4[[1]]["beta[2,1]", pctls], digits = 1)
formatCI(s4[[1]]["beta[2,2]", pctls], digits = 1)
formatCI(s4[[1]]["beta[2,3]", pctls], digits = 1)
formatCI(s4[[1]]["beta[2,4]", pctls], digits = 1)


## -------------------------
## Education 1-y effect
## =========================

formatCI(s4[[1]]["beta[4,1]", pctls], digits = 1)
formatCI(s4[[1]]["beta[4,2]", pctls], digits = 1)
formatCI(s4[[1]]["beta[4,3]", pctls], digits = 1)
formatCI(s4[[1]]["beta[4,4]", pctls], digits = 1)





###############################
###############################
## Correlation coefficient between individual-level adjustments
###############################
###############################

## individual adjustment correlations between outcomes
omega <- matrix(s4[[1]][grepl("Omega", sname), "mean"], ncol=4, byrow=TRUE)
dimnames(omega) <- list(c("PIB", "Tau", "WMH", "FA"),
                       c("PIB", "Tau", "WMH", "FA"))

round(omega,2)

## The strongest correlations are PiB to tau, and WMH to FA



###############################
###############################
## Assemble a dataset for plotting
###############################
###############################

## The estimated fixed effects per iteration
betas <- extract(fit4, "beta")[[1]]

## The covariate values & indicator of referral to ADRC
x.subj <- cbind(stan.list$x[!duplicated(stan.list$id),],
                'adrc'=stan.list$adrc[!duplicated(stan.list$id)])

## The random effects aka individual adjustments aka alphas
alphas <- extract(fit4, "alpha")$alpha

## The estimated referral effect per iteration
adrc.shift.mcmc <- extract(fit4, "adrc_shift")$adrc_shift

## Compute each subject's xTbeta + ADRC effect
xTbeta.pib <- rowMeans(x.subj %*% t(cbind(betas[, , 1], adrc.shift.mcmc)))
xTbeta.tau <- rowMeans(x.subj %*% t(cbind(betas[, , 2], adrc.shift.mcmc)))
xTbeta.wmh <- rowMeans(x.subj %*% t(cbind(betas[, , 3], adrc.shift.mcmc)))
xTbeta.fa  <- rowMeans(x.subj %*% t(cbind(betas[, , 4], adrc.shift.mcmc)))

## The estimated random effects
## aka individual-level adjustments aka alphas
alpha.pib <- colMeans(alphas$alpha[, , 1])
alpha.tau <- colMeans(alphas$alpha[, , 2])
alpha.wmh <- colMeans(alphas$alpha[, , 3])
alpha.fa <- colMeans(alphas$alpha[, , 4])

## Put data together as dataframe for subsetting and plotting
d <- data.frame(id = as.vector(stan.list$id),
                age = stan.list$age,
                stan.list$x,
                outcome = as.vector(stan.list$outcome),
                y = stan.list$y,
	        alpha.pib=alpha.pib[stan.list$id],
	        alpha.tau=alpha.tau[stan.list$id],
	        alpha.wmh=alpha.wmh[stan.list$id],
	        alpha.fa=alpha.fa[stan.list$id],
	        xTbeta.pib=xTbeta.pib[stan.list$id],
	        xTbeta.tau=xTbeta.tau[stan.list$id],
	        xTbeta.wmh=xTbeta.wmh[stan.list$id],
	        xTbeta.fa=xTbeta.fa[stan.list$id]
)

## Indicate whether a person ever had a given outcome measured
## This is useful for plotting going forward because the model will
## provide an estimated alpha for each person, each outcome
## regardless of them having ever been measured on that outcome.
## These estimates do no affect the fit, but we would like exclude
## those from the plots below.
d$has.pib <- ifelse(d$id %in% unique(d$id[d$outcome == 1]), 1, 0)
d$has.tau <- ifelse(d$id %in% unique(d$id[d$outcome == 2]), 1, 0)
d$has.wmh <- ifelse(d$id %in% unique(d$id[d$outcome == 3]), 1, 0)
d$has.fa  <- ifelse(d$id %in% unique(d$id[d$outcome == 4]), 1, 0)

d <- mutate(d,
            ## Compute adjusted ages
            aaage = age + alpha.pib + xTbeta.pib, ## 'aaage' is adjusted amyloid age
            atage = age + alpha.tau + xTbeta.tau, ## adjusted tau age
            awage = age + alpha.wmh + xTbeta.wmh, ## adjusted WMH age
            afage = age + alpha.fa + xTbeta.fa)   ## adjusted FA age

###############################
## Scatterplots biomarker by age
###############################

par(mfrow=c(2,2))
dat <- filter(d, outcome==1) %>% filter(!duplicated(id))
plot(exp(y) ~ age, data = dat, ylab='Amyloid, SUVR', xlab='Age, years')

dat <- filter(d, outcome==2) %>% filter(!duplicated(id))
plot(exp(y) ~ age, data = dat, ylab='Tau, SUVR', xlab='Age, years')

dat <- filter(d, outcome==3) %>% filter(!duplicated(id))
plot(y ~ age, data = dat, ylab='WMH, %', xlab='Age, years')

dat <- filter(d, outcome==4) %>% filter(!duplicated(id))
plot(I(1-y) ~ age, data = dat, ylab='FA', xlab='Age, years')




###############################
## Scatterplot biomarker by adjusted age
###############################

par(mfrow=c(2,2))
dat <- filter(d, outcome==1) %>% filter(!duplicated(id))
plot(exp(y) ~ aaage, data = dat, ylab='Amyloid, SUVR', xlab='Adjusted age, years', log='y')

dat <- filter(d, outcome==2) %>% filter(!duplicated(id))
plot(exp(y) ~ atage, data = dat, ylab='Tau, SUVR', xlab='Adjusted age, years', log='y')

dat <- filter(d, outcome==3) %>% filter(!duplicated(id))
plot(y ~ awage, data = dat, ylab='WMH, %', xlab='Adjusted age, years')

dat <- filter(d, outcome==4) %>% filter(!duplicated(id))
plot(I(1-y) ~ afage, data = dat, ylab='FA', xlab='Adjusted age, years')





###############################
## Pairwise scatterplots of biomarkers
###############################

d.wide <- select(d, id, age, outcome, y) %>%
    arrange(outcome) %>%
    pivot_wider(names_from=outcome, values_from=y, names_prefix='y') %>%
    mutate(pib=exp(y1),
           tau=exp(y2),
           wmh=y3,
           fa=1-y4)

## All observations. Ns will differ in each panel.
pairs(select(d.wide, pib,tau,wmh,fa), log=c(1,2))

## Most recent visit per person.
filter(d.wide, !duplicated(id)) %>%
    select(pib,tau,wmh,fa) %>%
    pairs(., log=c(1,2))




###############################
## Pairwise scatterplots of individual-level adjustments
###############################

# The model will estimate four alphas for each subject - one per outcome -
# whether or not a given subject was measured on all 4 outcomes.
# This is expected but a nuissance for plotting.
# We choose not to plot estimated alphas for those not measured for a given outcome.

# Reverse limits so bottom left corner is earlier onset, top right later onset

add.lines <- function(...){
    abline(h = 0, lty = 3) #horizontal line
    abline(v = 0, lty = 3) #vertical line
    }
    
add.ellipse <- function(x,y,...){
    ## draw an 80% confidence ellipse to show density.
    .t80 <- qt(0.9, df=dft) ## get dft from the model
    .mu <- c(mean(x, na.rm=T), mean(y, na.rm=T))
    .var <- var( cbind(x,y), use='pairwise')
    exy <- ellipse::ellipse(.var, center = .mu, t=.t80)
    lines(exy[,1], exy[,2], col='black', lwd=2, lty=2)
    }

add.hints <- function(){ 
    ## Add hints
    text(xpd=T, x=grconvertX(0.05,'npc'),y=grconvertY(0.05,'npc'), labels='earlier onset', adj=0)
    text(xpd=T, x=grconvertX(0.95,'npc'),y=grconvertY(0.95,'npc'), labels='later onset', adj=1)
}

## plotting function 
pfun1 <- function(x,y,...){
    ## set up plot with reversed limits so earlier onset is bottom left
    ## increase range of limits by 10% making space for text hints
    plot(y ~ x, xlim = rev(range(x))*1.1, ylim = rev(range(y))*1.1, ..., type = 'n')
    grid()
    add.lines()
    ## Add some transparency to the point color since there is overlap
    points(x=x, y=y, col=alpha('gray30', 0.4))
    add.ellipse(x=x, y=y)
    add.hints()
    }
    
par(mfcol=c(3,2))
dat <- filter(d, !duplicated(id) & has.pib == 1 & has.tau == 1)
pfun1(x=dat$alpha.pib, y=dat$alpha.tau, xlab = "Amyloid adjustment, y", ylab = "Tau adjustment, y")

dat <- filter(d, !duplicated(id) & has.pib == 1 & has.wmh == 1)
pfun1(x=dat$alpha.pib, y=dat$alpha.wmh, xlab = "Amyloid adjustment, y", ylab = "WMH adjustment, y")

dat <- filter(d, !duplicated(id) & has.pib == 1 & has.fa == 1)
pfun1(x=dat$alpha.pib, y=dat$alpha.fa,  xlab = "Amyloid adjustment, y", ylab = "FA adjustment, y")

dat <- filter(d, !duplicated(id) & has.tau == 1 & has.wmh == 1)
pfun1(x=dat$alpha.tau, y=dat$alpha.wmh,  xlab = "Tau adjustment, y", ylab = "WMH adjustment, y")

dat <- filter(d, !duplicated(id) & has.tau == 1 & has.fa == 1)
pfun1(x=dat$alpha.tau, y=dat$alpha.fa,  xlab = "Tau adjustment, y", ylab = "FA adjustment, y")

dat <- filter(d, !duplicated(id) & has.wmh == 1 & has.fa == 1)
pfun1(x=dat$alpha.wmh, y=dat$alpha.fa,  xlab = "WMH adjustment, y", ylab = "FA adjustment, y")



