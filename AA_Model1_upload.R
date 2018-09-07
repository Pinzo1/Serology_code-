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


#setwd("")

## call the libraries that are needed for analysis

library(MASS)
library(compiler)
library(binom)

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

data <- read.csv("", header = TRUE)  					 ## call in the csv datafile that needs to be analysed 

head(data)						 						 ## check that the data has been read in correctly 

#data <- data[",]										 ## You may need to subset data from the dataframe

data_2 <- as.data.frame(cbind(data$age, data$conc))

data_3 <- na.omit(data_2) 								## make sure there are no NAs

AB_data <- data_3										## assign the object to a new name that will be evalauated by the function

###############################################
## 0.2 bin the data by age group to plot

age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)							 ## age bins will be set from the youngest age group
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )		 ## in your data to the maximum age by = user defined

###############################################
## 0.2 Prepare data for plotting


N_bins <- length(age_bins) - 1 


GMT_bins      <- rep(NA, N_bins)

AB_range_bins <- matrix(NA, nrow=N_bins, ncol=3)						## generate an empty matrix to be filled with antibody levels for each age bin and the binomial CIs 
colnames(AB_range_bins) <- c("med", "low", "high")

for(i in 1:N_bins)
{
	index <- intersect( which(AB_data[,1]>age_bins[i]), which(AB_data[,1]<=age_bins[i+1]) ) 
	temp  <- AB_data[index,2]

	GMT_bins[i] = exp(mean(log(temp)))

	AB_range_bins[i,] <- quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


###############################################
## 0.3 Plot data
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,90), ylim=c(0.1,100),
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
## 1.1 MODEL   specifiy the AA model that may have
## generated the data 

par_MC <- c(2, 0.1, 0.55)  ## (alpha, rr, sigma)
 
model_M1 <- function(a, par)
{
	alpha <- par[1]
     rr    <- par[2]
	sigma <- par[3]
 
	AB_titre <- ( alpha/rr )*( 1 - exp(-rr*a) )
 
	AB_titre
}


model_M1 <- cmpfun(model_M1, options=list(optimize=3)) 


###################################################
## 1.2 LIKELIHOOD evaluate the binomial likelihood
 
loglike_M1 <- function( par )
{
	alpha <- par[1]
    rr    <- par[2]
	sigma <- par[3]

	AB_model <- model_M1( AB_data[,1], par )

	mu <- log(AB_model) 

	loglike <- - log(AB_data[,2]) - log(2.506628*sigma) - 0.5*( (log(AB_data[,2])-mu)/sigma )^2

      sum( loglike )
}

loglike_M1 <- cmpfun(loglike_M1, options=list(optimize=3))


###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M1 <- function( par )
{
	alpha <- par[1]
      rr    <- par[2]
	sigma <- par[3]

	######################################
	## Uniform prior on alpha ~ U(0,100)

	if( alpha>0 && alpha<100 )
	{
		prior_alpha <- log(1/100)
	}else{
		prior_alpha <- -LARGE
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
	## Uniform prior on sigma ~ U(0,10)

	if( sigma>0 && sigma<10 )
	{
		prior_sigma <- log(1/10)
	}else{
		prior_sigma <- -LARGE
	}

	prior <- prior_alpha + prior_rr + prior_sigma

	prior
}

prior_M1 <- cmpfun(prior_M1, options=list(optimize=3))


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


N_mcmc       <- 10000      ## Number of MCMC iterations
N_tune_start <- 300        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   <- 3000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      <- 4000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  <- 1           ## Scaler for step size
MCMC_accept <- 0           ## Track the MCMC acceptance rate

max_corr    <- 0.75        ## Maximum degree of correlation

#################################################
## 2.1 Robbins-munro step scaler


rm_scale <- function(step_scale, mc, log_prob)
{
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
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par) <- c("alpha", "rr", "sigma", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations

par_MC <- c(50, 0.1, 0.5)                 ## Initial guess: (alpha, rr, sigma)

Sigma_MC <- diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

loglike_MC <- loglike_M1( par_MC ) + prior_M1( par_MC )



for(mc in 1:N_mcmc)
{
	par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)



	if( prior_M1(par_MCp1) > -0.5*LARGE  ){
 
		loglike_MCp1 <- loglike_M1( par_MCp1 ) + prior_M1( par_MCp1 )


		log_prob <- min( loglike_MCp1-loglike_MC, 0 )           
                   
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
			cov_MC <- cov( MCMC_par[1:(mc-1),1:3] )

			if( min(diag(cov_MC))>1e-6 )
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

	MCMC_par[mc,1:3] <- par_MC
	MCMC_par[mc,4]   <- loglike_MC
	MCMC_par[mc,5]   <- prior_M1( par_MC )
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,2))



#####################################
## PANEL 1: alpha MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,1], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="alpha", 
main="alpha")




#####################################
## PANEL 2: rr MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,2], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="rr", 
main="rr")




#####################################
## PANEL 3: sigma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,3], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="sigma", 
main="sigma" )




#####################################
## PANEL 4: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,4], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
main="likelihood" )







#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(1,3))


#####################################
## PANEL 1: alpha MCMC posterior


DEN <- density( MCMC_burn[,1] )
	
QUANT <- quantile( MCMC_burn[,1], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="alpha", ylab="", 
main="posterior: alpha" )

	
low_index  = which(DEN$x<QUANT[1])
mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index = which(DEN$x>QUANT[3])

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
## PANEL 2: rr MCMC posterior


DEN <- density( MCMC_burn[,2] )
	
QUANT <- quantile( MCMC_burn[,2], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="rr", ylab="", 
main="posterior: rr" )

	
low_index  = which(DEN$x<QUANT[1])
mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index = which(DEN$x>QUANT[3])

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
## PANEL 3: sigma MCMC posterior


DEN <- density( MCMC_burn[,3] )
	
QUANT <- quantile( MCMC_burn[,3], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="sigma", ylab="", 
main="posterior: sigma" )

	
low_index  = which(DEN$x<QUANT[1])
mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
high_index = which(DEN$x>QUANT[3])

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


par_median <- apply(X=MCMC_burn[,1:3], MARGIN=2, FUN=median)



age_seq <- seq(from=0, to=90, by=1)

M1_predict <- model_M1(age_seq, par_median )





###############################################
## 3.1 Plot data and model prediction
 
par(mfrow=c(1,1))


plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,60), ylim=c(0.1,500),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Antibody acquisition Model 1 fit"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}


points(x=age_seq, y=M1_predict, 
type='l', lwd=3, col="blue")


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

theta_bar = apply( X=MCMC_burn[,1:3], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M1( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*mean( MCMC_burn[,4] - MCMC_burn[,5] )


######################################
## Effective number of model parameters
	
pD = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC = pD + D_bar


######################################
## Median of the posterior

theta_bar = apply( X=MCMC_burn[,1:3], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M1( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*median( MCMC_burn[,4] - MCMC_burn[,5] )


######################################
## Effective number of model parameters
	
pD_median = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC_median = pD_median + D_bar

