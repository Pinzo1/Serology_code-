#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White                                                   ##
## Institut Pasteur                                                ##
## michael.white@pasteur.fr                                        ##
## m.white08@imperial.ac.uk                                        ##
##                                                                 ##
#####################################################################
#####################################################################


rm(list=ls())

#setwd("")

## call the libraries that are needed for analysis

library(MASS)
library(compiler)

###############################################
###############################################
##          ##                               ##
##   ####   ##  ####    ####  ######  ####   ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ##  ##  ##  ## ######   ##   ######  ##
##  ##  ##  ##  ## ##  ##  ##   ##   ##  ##  ##
##   ####   ##  ####   ##  ##   ##   ##  ##  ##
##          ##                               ##
###############################################
###############################################
 
###############################################
## 0.1 Read in data



data <- read.csv(".csv", header = TRUE)					## call in the csv datafile that needs to be analysed 

head(data)												## check that the data has been read in correctly 

data <- data[,]											## You may need to subset data from the dataframe

data_2 <- as.data.frame(cbind(data$age, data$conc))

data_3 <- na.omit(data_2) 								## make sure there are no NAs

AB_data <- data_3										## assign the object to a new name that will be evalauated by the function

###############################################
## 0.2 bin the data by age group to plot

age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)							## age bins will be set from the youngest age group
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )		## in your data to the maximum age by = user defined

###############################################
## 0.2 Prepare data for plotting

N_bins <- length(age_bins) - 1 


GMT_bins      <- rep(NA, N_bins)

AB_range_bins <- matrix(NA, nrow=N_bins, ncol=3)									## generate an empty matrix to be filled with antibody levels for each age bin and the binomial CIs 
colnames(AB_range_bins) <- c("med", "low", "high")

for(i in 1:N_bins)
{
	index <- intersect( which(AB_data[,1]>age_bins[i]), which(AB_data[,1]<=age_bins[i+1]) ) 
	temp  <- AB_data[index,2]

	GMT_bins[i] <- exp(mean(log(temp)))

	AB_range_bins[i,] <- quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


###############################################
## 0.3 Plot data
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,60), ylim=c(0.1,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Cross-sectional antibody data"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}



###################################################
###################################################
##        ##                                     ##
##   ##   ##  #     #  ####  ####   ##### ##     ##
##  ###   ##  ##   ## ##  ## ## ##  ##    ##     ##
##   ##   ##  ####### ##  ## ##  ## ####  ##     ##
##   ##   ##  ## # ## ##  ## ## ##  ##    ##     ##
##  ####  ##  ##   ##  ####  ####   ##### #####  ##
##        ##                                     ##
###################################################
###################################################


################################################### 
## 1.1 MODEL  specifiy the AA model that may have
## generated the data  

par_MC <- c(5, 0.1, 0.1, 10, 0.6)  ## (alpha_0, gamma, rr, time_c, sigma)


model_M2 <- function(a, par)
{
	alpha_0 <- par[1]
	gamma   <- par[2]
      rr      <- par[3]
	time_c  <- par[4]
	sigma   <- par[5]

	alpha_c <- gamma*alpha_0

	age_xx  <- a - time_c 
 
	if( age_xx<=0 ){
		AB_titre <- ( alpha_c/rr )*( 1 - exp(-rr*a) )
	}

	if( (age_xx>0) && (age_xx<a) ){
		AB_titre <- ( alpha_0/rr )*exp(-rr*a)*( exp(age_xx*rr) - 1 ) +
			     ( alpha_c/rr )*( 1 - exp( -rr*(a-age_xx) ) )
	}

	if( a <= age_xx ){
		AB_titre <- ( alpha_0/rr )*( 1 - exp(-rr*a) )
	} 

	AB_titre
}

model_M2 <- cmpfun(model_M2, options=list(optimize=3)) 



###################################################
## 1.2 define and compute the LIKELIHOOD 
 
loglike_M2 <- function( par )
{
	alpha_0 <- par[1]
	gamma   <- par[2]
      rr      <- par[3]
	time_c  <- par[4]
	sigma   <- par[5]

	AB_model <- sapply(AB_data[,1], model_M2, par=par)

	mu <- log(AB_model) 

	loglike <- -log(AB_data[,2]) - log(2.506628*sigma) - 0.5*( (log(AB_data[,2])-mu)/sigma )^2

      sum( loglike )
}

loglike_M2 <- cmpfun(loglike_M2, options=list(optimize=3))

###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M2 <- function( par )
{
	alpha_0 <- par[1]
	gamma   <- par[2]
      rr      <- par[3]
	time_c  <- par[4]
	sigma   <- par[5]

	######################################
	## Uniform prior on alpha_0 ~ U(0,1000)

	if( alpha_0>0 && alpha_0<1000 )
	{
		prior_alpha_0 <- log(1/1000)
	}else{
		prior_alpha_0 <- -LARGE
	}

	######################################
	## Uniform prior on gamma ~ U(0,1)

	if( gamma>0 && gamma<1 )
	{
		prior_gamma <- log(1/1)
	}else{
		prior_gamma <- -LARGE
	}

	######################################
	## Uniform prior on rr ~ U(0,10)

	if( rr>0 && rr<10 )
	{
		prior_rr <- log(1/10)
	}else{
		prior_rr <- -LARGE
	}

	######################################
	## Uniform prior on time_c ~ U(0,60)

	if( time_c>0 && time_c<60 )
	{
		prior_time_c <- log(1/60)
	}else{
		prior_time_c <- -LARGE
	}

	######################################
	## Uniform prior on sigma ~ U(0,10)

	if( sigma>0 && sigma<10 )
	{
		prior_sigma <- log(1/10)
	}else{
		prior_sigma <- -LARGE
	}

	prior <- prior_alpha_0 + prior_gamma + prior_rr + prior_time_c + prior_sigma

	prior
}

prior_M2 <- cmpfun(prior_M2, options=list(optimize=3))


#################################################
#################################################
##          ##                                 ##
##   ####   ##  #     #  ####  #     #  ####   ##
##  ##  ##  ##  ##   ## ##  ## ##   ## ##  ##  ##
##     ##   ##  ####### ##     ####### ##      ##
##    ##    ##  ## # ## ##  ## ## # ## ##  ##  ##
##   #####  ##  ##   ##  ####  ##   ##  ####   ##
##          ##                                 ##
#################################################
#################################################


N_mcmc       <- 50000      ## Number of MCMC iterations
N_tune_start <- 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   <- 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      <- 6000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  <- 1           ## Scaler for step size
MCMC_accept <- 0           ## Track the MCMC acceptance rate

max_corr    <- 0.75        ## Maximum degree of correlation


#################################################
## 2.1 Robbins-munro step scaler


rm_scale <- function(step_scale, mc, log_prob){

	dd <- exp(log_prob)
	if( dd < -30 ){ dd <- 0 }
	dd <- min( dd, 1 )

	rm_temp <- ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
		
	out <- step_scale*exp(rm_temp)
	
	out <- max( out, 0.05 )
	out <- min( out, 5)
	out
}


#################################################
## 2.2 Prepare object for MCMC fitting
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par) <- c("alpha_0", "gamma", "rr", "time_c", "sigma", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations


par_MC <- c(50, 0.5, 0.1, 10, 0.6)       ## (alpha_0, gamma, rr, time_c, sigma)


Sigma_MC <- diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

                  

loglike_MC <- loglike_M2( par_MC ) + prior_M2( par_MC )




for(mc in 1:N_mcmc)
{
	par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)
 
	if( prior_M2(par_MCp1) > -0.5*LARGE ){
 
		loglike_MCp1 <- loglike_M2( par_MCp1 ) + prior_M2( par_MCp1 )


		log_prob = min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC <- par_MCp1
			
			loglike_MC  <- loglike_MCp1
			MCMC_accept <- MCMC_accept + 1                       
		}

		#######################################
		## RM scaling of proposal step size

		if( mc < N_adapt ){
			step_scale <- rm_scale( step_scale, mc, log_prob)
		}

		#######################################
		## Adaptive tuning of covariance matrix

		if( (mc > N_tune_start) && (mc < N_tune_end) )
		{
			cov_MC <- cov( MCMC_par[1:(mc-1),1:5] )

			###########################
			## Checks for tuning

			if( min(diag(cov_MC)) > 1e-6 )
			{
				###########################
				## Check for high degree of correlation

				sd_MC_inv <- 1/sqrt(diag(cov_MC))

				corr_MC <- t(t(cov_MC*sd_MC_inv)*sd_MC_inv)

				corr_MC[intersect( which( corr_MC > max_corr ), which(corr_MC<0.99999) )] <- max_corr
				corr_MC[  which( corr_MC < -max_corr )] <- -max_corr

				cov_MC <- t(t(corr_MC*(1/sd_MC_inv))*(1/sd_MC_inv))
 
				Sigma_MC <- cov_MC
			}
		}
	}


	MCMC_par[mc,1:5] <- par_MC
	MCMC_par[mc,6]   <- loglike_MC
	MCMC_par[mc,7]   <- prior_M2( par_MC )
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,3))



#####################################
## PANEL 1: alpha_0 MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,1], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="alpha_0", 
main="alpha_0")



#####################################
## PANEL 2: gamma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,2], 
ylim=c(0,1),
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="gamma", 
main="gamma")



#####################################
## PANEL 3 rr MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,3], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="rr", 
main="rr")



#####################################
## PANEL 4: time_c MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,4], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="time_c", 
main="time_c" )



#####################################
## PANEL 5: sigma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,5], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="sigma", 
main="sigma" )



#####################################
## PANEL 6: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,6], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
ylim=quantile( MCMC_par[,6], prob=c(0.01,1)),
main="likelihood" )


	




#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,3))


#####################################
## PANEL 1: alpha_0 MCMC posterior


DEN <- density( MCMC_burn[,1] )
	
QUANT <- quantile( MCMC_burn[,1], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="alpha_0", ylab="", 
main="posterior: alpha_0" )

	
low_index  <- which(DEN$x<QUANT[1])
mid_index  <- intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index <- which(DEN$x>QUANT[3])

polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
	   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
             col="pink")

polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
	   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
	   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
	   col="pink")

points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)





#####################################
## PANEL 2: gamma MCMC posterior


DEN <- density( MCMC_burn[,2] )
	
QUANT <- quantile( MCMC_burn[,2], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="gamma", ylab="", 
main="posterior: gamma" )

	
low_index  <- which(DEN$x<QUANT[1])
mid_index  <- intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index <- which(DEN$x>QUANT[3])

polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
	   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
             col="pink")

polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
	   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
	   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
	   col="pink")

points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)




#####################################
## PANEL 3: rr MCMC posterior


DEN <- density( MCMC_burn[,3] )
	
QUANT <- quantile( MCMC_burn[,3], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="rr", ylab="", 
main="posterior: rr" )

	
low_index  <- which(DEN$x<QUANT[1])
mid_index  <- intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index <- which(DEN$x>QUANT[3])

polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
	   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
             col="pink")

polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
	   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
	   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
	   col="pink")

points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)



#####################################
## PANEL 4: time_c MCMC posterior


DEN <- density( MCMC_burn[,4] )
	
QUANT <- quantile( MCMC_burn[,4], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, 60),
xlab="sigma", ylab="", 
main="posterior: sigma" )

	
low_index  <- which(DEN$x<QUANT[1])
mid_index  <- intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index <- which(DEN$x>QUANT[3])

polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
	   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
             col="pink")

polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
	   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
	   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
	   col="pink")

points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)



#####################################
## PANEL 5: sigma MCMC posterior


DEN <- density( MCMC_burn[,5] )
	
QUANT <- quantile( MCMC_burn[,5], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="sigma", ylab="", 
main="posterior: sigma" )

	
low_index  <- which(DEN$x<QUANT[1])
mid_index  <- intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index <- which(DEN$x>QUANT[3])

polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
	   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
             col="pink")

polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
	   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
	   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
	   col="pink")

points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)




#############################################
#############################################
##          ##                             ##
##   ####   ##  ###### #####  ###  ######  ##
##  ##  ##  ##    ##   ##    ##      ##    ##
##     ##   ##    ##   ####   ###    ##    ##
##  ##  ##  ##    ##   ##       ##   ##    ##
##   ####   ##    ##   #####  ###    ##    ##
##          ##                             ##
#############################################
#############################################

#############################################
## 3.1 Extract posterior medians and 
##     calculate model prediction

#load("Model2_Rombo_CT694.RData")

par_median <- apply(X=MCMC_burn[,1:5], MARGIN=2, FUN=median)


age_seq <- seq(from=0, to=90, by=1)

M2_predict <- sapply(age_seq, model_M2, par=par_median)




###############################################
## 3.1 Plot data and model prediction
jpeg("AA_Model2_Rombo_CT694.jpg")

par(mfrow=c(1,1))


plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,80), ylim=c(0.1,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Antibody acquisition Model 2 - Rombo CT694"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}


points(x=age_seq, y=M2_predict, 
type='l', lwd=3, col="green")

######################################
######################################
##          ##                      ##
##  ##      ##  ####   ####  ####   ##
##  ## ##   ##  ## ##   ##  ##  ##  ##
##  ######  ##  ##  ##  ##  ##      ##
##     ##   ##  ## ##   ##  ##  ##  ##
##     ##   ##  ####   ####  ####   ## 
##          ##                      ##
######################################
######################################
##
## Estimate the Deviance Information Criterion.
## Note that the deviance is -2*log(likelihood)

######################################
## Mean of the posterior

theta_bar = apply( X=MCMC_burn[,1:5], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*mean( MCMC_burn[,6] - MCMC_burn[,7] )


######################################
## Effective number of model parameters
	
pD = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC = pD + D_bar

######################################
## Median of the posterior

theta_bar = apply( X=MCMC_burn[,1:5], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*median( MCMC_burn[,6] - MCMC_burn[,7] )


######################################
## Effective number of model parameters
	
pD_median = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC_median = pD_median + D_bar



