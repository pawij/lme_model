# Power analysis used in Predict, Track, Image-HD meta-analysis
# Author: Peter Wijeratne (p.wijeratne@ucl.ac.uk)
set.seed(42)
library(nlme)
library(lme4)
options(warn=-1)

###################################################################################################

isLME4=TRUE
isImaging=TRUE
data.path = "/home/paw/code/Python-EBM/pti_data_combined_test.csv"

###################################################################################################

if(isImaging){
    vnames = c('Caudate','Pallidum','Putamen','Insula.White.Matter','Non.ventricular.CSF')
} else{
    vnames = c('sdmt','swrt','tms')
}

for(vIdx in 1:length(vnames)){

    volname=vnames[vIdx]

    for(simtype in 1:2){
        print(volname)
        if(simtype==1){
            site.values <- c(4) # N sites
            subject_site.values <- c(30, 40, 50, 75, 100, 200, 300, 500) # N subs / site
        } else{
            site.values <- c(6, 8, 10, 15, 20, 40, 60, 100) # N sites
            subject_site.values <- c(20) # N subs / site
        }
        effect.values <- c(0.2,0.4)
        time.values <- c(3)

        lme.data <- read.table(data.path,header=TRUE,sep=",",na.strings="NA",dec=".")

        ############################################################################## cut controls - here untreated HD is the control sample
        lme.data <- lme.data[lme.data$group!=0,]
        ##############################################################################

        lme.data$time <- as.numeric(lme.data$time)
        lme.data$age <- as.numeric(lme.data$age)
        lme.data$tiv <- as.numeric(lme.data$tiv)
        lme.data$cag <- as.numeric(lme.data$cag)
        lme.data$sex <- factor(lme.data$sex)
        lme.data$subject <- factor(lme.data$subject)
        lme.data$site <- factor(lme.data$site)
        lme.data$study <- factor(lme.data$study)
        lme.data$scanner <- factor(lme.data$scanner)
        lme.data$group <- factor(lme.data$group)
        lme.data$field <- factor(lme.data$field)

        lme.data$tms <- as.numeric(lme.data$tms)
        lme.data$swrt <- as.numeric(lme.data$swrt)
        lme.data$sdmt <- as.numeric(lme.data$sdmt)
        # convert volumes to z-scores
        vars <- colnames(lme.data)
        startIdx = match('Non.ventricular.CSF',vars)
        for(i in startIdx:length(vars)){
            lme.data[,vars[i]] <- lme.data[,vars[i]] / sd(lme.data[,vars[i]])
        }
        # extract volume
        y = lme.data[,volname]

        if(isLME4){
            if(isImaging){
                mod <- lmer(
                    y ~ time + age + sex + tiv + group + cag + voxel + scanner + time:age + (1 + time | site/subject), # for imaging biomarkers
                    data = lme.data
                )
            } else{
                mod <- lmer(
                    y ~ time + age + sex + group + (1 + time | site/subject), # for clinical biomarkers
                    data = lme.data
                )
            }
        } else{
            if(isImaging){
                mod <- lme(
                    y ~ time + age + sex + tiv + group + cag + voxel + scanner + time:age, # for imaging biomarkers
                    random = ~1 + time | site/subject,
                    data = lme.data,
                    control = lmeControl(opt = "optim",
                                         maxIter = 800, msMaxIter = 800)
                )
            } else{
                mod <- lme(
                    y ~ time + age + sex + group, # for clinical biomarkers
                    random = ~1 + time | site/subject,
                    data = lme.data,
                    control = lmeControl(opt = "optim",
                                         maxIter = 800, msMaxIter = 800)
                )
            }
        }
        #        print(summary(mod))

        getHyperparam<-function(x){
            # fixed effects
            sum_obj <- summary(x)
            mu.a.true <- as.numeric(fixef(x)[1]) # intercept
            g.0.true <- as.numeric(fixef(x)[2]) # time
            # random effects
            if(isLME4){
                vc <- as.data.frame(VarCorr(x))
                sigma.a0.true <- as.numeric(vc[1,'sdcor']) # subject:site	
                sigma.b0.true <- as.numeric(vc[2,'sdcor']) # subject:site:time
                sigma.a.true <- as.numeric(vc[4,'sdcor']) # site
                sigma.b.true <- as.numeric(vc[5,'sdcor']) # site:time
                sigma.y.true <- as.numeric(vc[7,'sdcor']) # Residual
            } else{
                vc <- VarCorr(x)
                sigma.a0.true <- as.numeric(vc[5,2]) # subject:site	
                sigma.b0.true <- as.numeric(vc[6,2]) # subject:site:time
                sigma.a.true <- as.numeric(vc[2,2]) # site
                sigma.b.true <- as.numeric(vc[3,2]) # site:time
                sigma.y.true <- as.numeric(vc[7,2]) # Residual
            }
    
            hp<-c(g.0.true, mu.a.true, sigma.y.true, sigma.a.true, sigma.a0.true, sigma.b.true, sigma.b0.true)
            names(hp)<-c("g.0.true", "mu.a.true", "sigma.y.true", "sigma.a.true", "sigma.a0.true", "sigma.b.true", "sigma.b0.true")

            return(hp)
        }

        HP <- getHyperparam(mod)

        create_fake <- function(J, K, L, t, HP){
            # J : number of sites
            # K : number of subjects / site
            # L : number of years
            # t : treatment effect
            # HP: hyper-parameters
            time <- rep(seq(0,2,length=L), J*K)
            subject <- rep(1:(J*K), each=L)
            site0 <- sample(rep (1:J, K))
            site <- factor(site0[subject])
            treatment0 <- sample(rep (0:1, J*K/2))
            treatment <- treatment0[subject]
            # parameters from model fit
            slope <- as.numeric( HP['g.0.true'] ) # time coefficient            
            effect <- -1.0*as.numeric(t)*slope # treatment coefficient
            intercept <- as.numeric( HP['mu.a.true'] ) # intercept
            # fixed effects
            b.true <- (slope + effect*treatment0) # length J*K; multiply across by time later
            # random effects
            sigma.y.true <- as.numeric( HP['sigma.y.true'] )
            sigma.site <- as.numeric( HP['sigma.a.true'] )
            sigma.subject <- as.numeric( HP['sigma.a0.true'] )
            sigma.site.time <- as.numeric( HP['sigma.b.true'] )
            sigma.subject.time <- as.numeric( HP['sigma.b0.true'] )
            re.site <- rnorm(J, 0, sigma.site)
            re.subject <- rnorm(J*K, 0, sigma.subject)
            re.site_time <- rnorm(J, 0, sigma.site.time)
            re.subject_time <- rnorm(J*K, 0, sigma.subject.time)
            a <- intercept + re.site[site] + re.subject[subject]
            b <- b.true[subject] + re.site_time[site] + re.subject_time[subject] + re.site_time[site]*treatment + re.subject_time[subject]*treatment
            # simulate
            y <- rnorm(J*K*L, a + b*time, sigma.y.true)

            return(data.frame( y, time, subject, treatment, site ))
        }
        
        trial.power <- function(n_site, n_subject_site, n_time, effect, n.sims=1000){
            signif <- rep (NA, n.sims)
            for (s in 1:n.sims){
                fake <- create_fake(n_site, n_subject_site, n_time, effect, HP)
                
                # if we try to fit to y, but do not sufficiently capture the variance of either site or subject, then the std dev of time:treatment will be large
                if(isLME4){
                    lme.power <- lmer(
                        y ~ time + time:treatment + (1 + time | site/subject),
                        data=fake
                    )
                    theta.hat <- coef(summary(lme.power))['time:treatment', 'Estimate']
                    theta.se <- coef(summary(lme.power))['time:treatment', 'Std. Error']
                } else{
                    lme.power <- lme(
                        y ~ time + time:treatment,
                        random = ~1 + time | site/subject,
                        data=fake,
                        control = lmeControl(opt = "optim",
                                             maxIter = 800, msMaxIter = 800)
                    )
                    theta.hat <- coef(summary(lme.power))['time:treatment', 'Value']
                    theta.se <- coef(summary(lme.power))['time:treatment', 'Std.Error']
                }
                # want CSF volume not to increase - hence we want significantly negative treatment
                if(volname=='Non.ventricular.CSF' || volname=='tms' || volname=='sdmt'){
                    signif[s] <- ifelse (theta.hat + 2*theta.se < 0, 1, 0)
                }
                # want solid regions' volume not to decrease - hence we want significantly positive treatment
                else{
                    signif[s] <- ifelse (theta.hat - 2*theta.se > 0, 1, 0)
                }
            }
            power <- mean(signif)
            return(power)
        }

        power.out <- data.frame(matrix(ncol = 5, nrow = 0))
        names(power.out) <- c('power','effect','sites','subs_site','time')

        for(i1 in 1:length(effect.values)){
            for(i2 in 1:length(site.values)){
                for(i3 in 1:length(subject_site.values)){
                    for(i4 in 1:length(time.values)){
                        power <- trial.power(site.values[i2], subject_site.values[i3], time.values[i4], effect.values[i1])
                        cat("Power =", power ,", for effect =", effect.values[i1] ,", n_sites =", site.values[i2], ", n_subject_site =", subject_site.values[i3], ", n_time =", time.values[i4], "\n")
                        newRow <- data.frame(power, effect.values[i1], site.values[i2], subject_site.values[i3], time.values[i4])
                        names(newRow) <- c('power','effect','sites','subs_site','time')
                        power.out <- rbind(power.out,newRow)
                    }
                }
            }
        }
        simname='_site'
        if(simtype==1){
            simname='_subs'
        }
        write.csv(power.out, file=paste(volname,simname,'_power.csv',sep=''), row.names=FALSE)
    }
}
