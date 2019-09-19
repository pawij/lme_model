# example power analysis of treatment effect using a mixed effects model
set.seed(42)
library(lme4)
options(warn=-1)

### USER INPUT START
###
# number of sites, subjects, treatment effect sizes, and times to simulate
site.values <- c(6, 8, 10, 15, 20, 40, 60, 100)
subject_site.values <- c(20)
effect.values <- c(0.4,0.6,0.8,1.0)
time.values <- c(3)
# desired biomarker parameters
beta.0 <- 6.
beta.time <- -0.2
###
### USER INPUT END

# simulate data
create_fake <- function(J, K, L, t){
    # J : number of sites
    # K : number of subjects / site
    # L : number of years
    time <- rep(seq(0,2,length=L), J*K)
    subject <- rep(1:(J*K), each=L)
    site0 <- sample(rep (1:J, K))
    site <- factor(site0[subject])
    treatment0 <- sample(rep (0:1, J*K/2))
    treatment <- treatment0[subject]
    # time coefficient
    g.0.true <- beta.time#-0.21
    # treatment coefficient
    g.1.true <- -1.0*as.numeric(t)*g.0.true
    # intercept
    mu.a.true <- beta.0#6.89
    # fixed effects
    b.true <- (g.0.true + g.1.true*treatment0) # length J*K; multiply across by time later
    # random effects
    sigma.y.true <- 0.2    
    sigma.site <- 0.35
    sigma.subject <- 0.59
    sigma.site.time <- 0.13
    sigma.subject.time <- 0.01
    re.site <- rnorm(J, 0, sigma.site)
    re.subject <- rnorm(J*K, 0, sigma.subject)
    re.site_time <- rnorm(J, 0, sigma.site.time)
    re.subject_time <- rnorm(J*K, 0, sigma.subject.time)
    # simulate
    a <- mu.a.true + re.site[site] + re.subject[subject]
    b <- b.true[subject] + re.site_time[site] + re.subject_time[subject] + re.site_time[site]*treatment + re.subject_time[subject]*treatment
    y <- rnorm(J*K*L, a + b*time, sigma.y.true)
    return(data.frame( y, time, subject, treatment, site ))
}
# simulate hypothetical treatment
trial.power <- function(n_site, n_subject_site, n_time, effect, n.sims=1000){
  signif <- rep (NA, n.sims)
  for (s in 1:n.sims){
    fake <- create_fake(n_site, n_subject_site, n_time, effect)    
    # if we try to fit to y, but do not sufficiently capture the variance of either site or subject, then the std dev of time:treatment will be large
    lme.power <- lmer(
                      y ~ time + time:treatment + (1 + time | site/subject),
    		      data=fake
		     )
    theta.hat <- coef(summary(lme.power))['time:treatment', 'Estimate']
    theta.se <- coef(summary(lme.power))['time:treatment', 'Std. Error']
    signif[s] <- ifelse (theta.hat - 2*theta.se > 0, 1, 0)
  }
  power <- mean(signif)
  return(power)
}
# print stuff
for(i1 in 1:length(effect.values)){
 for(i2 in 1:length(site.values)){
  for(i3 in 1:length(subject_site.values)){
   for(i4 in 1:length(time.values)){
    power <- trial.power(site.values[i2], subject_site.values[i3], time.values[i4], effect.values[i1])
    cat("Power =", power ,", for effect =", effect.values[i1] ,", n_sites =", site.values[i2], ", n_subject_site =", subject_site.values[i3], ", n_time =", time.values[i4], "\n")
   }
  }
 }
}
