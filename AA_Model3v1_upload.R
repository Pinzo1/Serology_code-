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
##


age_SD = 18									## define the age of sexual debut


data <- read.csv("", header = TRUE)			## call in the csv datafile that needs to be analysed

head(data)									## check that the data has been read in correctly 


data <- data[,]								## you may need to subset the data

data_2 <- as.data.frame(cbind(data$age, data$conc))

data_3 <- na.omit(data_2) 					## make sure there are no NAs

AB_data <- data_3

AB_cs1 <- data_3							## assign the object to a new name that will be evalauated by the function


###############################################
## 0.1 Prepare data for plotting

age_bins     = seq(from=0, to=90, by=5)													## age bins will be set from the youngest age group
age_bins_mid = 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )	## in your data to the maximum age by = user defined

N_bins = length(age_bins) - 1 


GMT_bins_cs1 = rep(NA, N_bins)

AB_range_bins_cs1 = matrix(NA, nrow=N_bins, ncol=3)
colnames(AB_range_bins_cs1) = c("med", "low", "high")

for(i in 1:N_bins)
{
	index = intersect( which(AB_cs1[,1]>age_bins[i]), which(AB_cs1[,1]<=age_bins[i+1]) ) 
	temp  = AB_cs1[index,2]

	GMT_bins_cs1[i] = exp(mean(log(temp)))

	AB_range_bins_cs1[i,] = quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


###############################################
## 0.5 Plot data
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins_cs1, 
pch=15, cex=2,
log="y", xlim=c(0,70), ylim=c(0.1,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Cross-sectional antibody data"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins_cs1[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins_cs1[i,3], 
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

time_0  = 80  ## time when lambda = lambda_0 (i.e. 20 years ago)
age_max = 80   ## maximum age in population - dataset dependent

################################################### 
## 1.1 MODEL specifiy the AA model that may have
## generated the data 

model_M3 <- function(a, t_survey, par_M3)
{
	alpha_0 = par_M3[1]
	gamma   = par_M3[2]
      rr      = par_M3[3]
	sigma   = par_M3[4]

	alpha_c = gamma*alpha_0

	AB_titre = (1/rr)*( alpha_c + (alpha_0-alpha_c)*((a + t_survey)/time_0) + (alpha_0 - alpha_c)/(rr*time_0) )*(1-exp(-rr*a)) - ((alpha_0-alpha_c)/(rr*time_0))*a
}

model_M3 = cmpfun(model_M3, options=list(optimize=3)) 

###################################################
## 1.2 define and compute the LIKELIHOOD 

loglike_M3_cs1 <- function( par_M3 )
{
	alpha_0 = par_M3[1]
	gamma   = par_M3[2]
      rr      = par_M3[3]
	sigma   = par_M3[4]

	AB_model = sapply(AB_cs1[,1], model_M3, t_survey=15, par_M3=par_M3)

	mu = log(AB_model) 

	loglike = -log(AB_cs1[,2]) - log(2.506628*sigma) - 0.5*( (log(AB_cs1[,2])-mu)/sigma )^2

      sum( loglike )
}

loglike_M3_total = cmpfun(loglike_M3_cs1, options=list(optimize=3)) 



###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M3 = function( par_M3 )
{
	alpha_0 = par_M3[1]
	gamma   = par_M3[2]
      rr      = par_M3[3]
	sigma   = par_M3[4]

	######################################
	## Uniform prior on alpha_0 ~ U(0,1000)

	if( alpha_0>0 && alpha_0<1000 )
	{
		prior_alpha_0 = log(1/1000)
	}else{
		prior_alpha_0 = -LARGE
	}

	######################################
	## Uniform prior on gamma ~ U(0,1)

	if( gamma>0 && gamma<1 )
	{
		prior_gamma = log(1/1)
	}else{
		prior_gamma = -LARGE
	}

	######################################
	## Uniform prior on rr ~ U(0,10)

	if( rr>0 && rr<10 )
	{
		prior_rr = log(1/10)
	}else{
		prior_rr = -LARGE
	}


	######################################
	## Uniform prior on sigma ~ U(0,10)

	if( sigma>0 && sigma<10 )
	{
		prior_sigma = log(1/10)
	}else{
		prior_sigma = -LARGE
	}

	prior = prior_alpha_0 + prior_gamma + prior_rr + prior_sigma

	prior
}

prior_M3 = cmpfun(prior_M3, options=list(optimize=3))


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


N_mcmc       <- 30000      ## Number of MCMC iterations
N_tune_start <- 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   <- 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      <- 6000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  <- 1           ## Scaler for step size
MCMC_accept <- 0           ## Track the MCMC acceptance rate

max_corr    <- 0.75        ## Maximum degree of correlation


#################################################
## 2.1 Robbins-munro step scaler


rm_scale = function(step_scale, mc, log_prob)
{
	dd = exp(log_prob)
	if( dd < -30 ){ dd <- 0 }
	dd = min( dd, 1 )

	rm_temp = ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
		
	out = step_scale*exp(rm_temp)
	
	out = max( out, 0.02 )
	out = min( out, 2)
	out
}


#################################################
## 2.2 Prepare object for MCMC fitting
 
MCMC_par           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par) = c("alpha_0", "gamma", "rr", "sigma", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations


par_MC = c(40, 0.5, 0.1, 0.6)       ## (alpha_0, gamma, rr, sigma)


Sigma_MC = diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

                  
prior_MC = prior_M3( par_MC )

loglike_MC = loglike_M3_total( par_MC ) + prior_MC




for(mc in 1:N_mcmc)
{
	par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)

	prior_MCp1 = prior_M3( par_MCp1 )
 
	if( prior_MCp1 > -0.5*LARGE )
	{
 		loglike_MCp1 = loglike_M3_total( par_MCp1 ) + prior_MCp1


		log_prob = min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC <- par_MCp1
			
			loglike_MC  = loglike_MCp1
			prior_MC    = prior_MCp1

			MCMC_accept = MCMC_accept + 1                       
		}

		#######################################
		## RM scaling of proposal step size

		if( mc < N_adapt )
		{
			step_scale = rm_scale( step_scale, mc, log_prob)
		}

		#######################################
		## Adaptive tuning of covariance matrix

		if( (mc > N_tune_start) && (mc < N_tune_end) )
		{
			cov_MC = cov( MCMC_par[1:(mc-1),1:4] )
		}
	}


	MCMC_par[mc,1:4] = par_MC
	MCMC_par[mc,5]   = loglike_MC
	MCMC_par[mc,6]   = prior_MC
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,3))

for(k in 1:4)
{
	#####################################
	## PANEL k: 

	plot(x=1:N_mcmc, y=MCMC_par[,k], 
	pch=19, col="grey", cex=0.25,
	xlab="MCMC iteration", ylab="", 
	main=colnames(MCMC_par)[k])
}


#####################################
## PANEL 6: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,5], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
ylim=quantile( MCMC_par[,5], prob=c(0.01,1)),
main="likelihood" )


	




#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn = MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,2))


#####################################
## PANEL 1: MCMC posterior

for(k in 1:4)
{
	DEN = density( MCMC_burn[,k] )
	
	QUANT = quantile( MCMC_burn[,k], prob=c(0.025, 0.5, 0.975) )

	plot(x=DEN$x, y=DEN$y, type='l',
	xlim=c(0, max(DEN$x)),
	xlab="alpha_0", ylab="", 
	main=colnames(MCMC_par)[k] )


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
}





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

par_median = apply(X=MCMC_burn[,1:4], MARGIN=2, FUN=median)


age_seq = seq(from=0, to=90, by=1)

M3_predict_cs1 = sapply(age_seq, model_M3, t_survey=15, par=par_median)


###############################################
## 0.2 bin the data by age group to plot

age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )

###############################################
## 0.2 Prepare data for plotting

 
N_bins <- length(age_bins) - 1 


GMT_bins      <- rep(NA, N_bins)

AB_range_bins <- matrix(NA, nrow=N_bins, ncol=3)
colnames(AB_range_bins) <- c("med", "low", "high")

for(i in 1:N_bins)
{
	index <- intersect( which(AB_data[,1]>age_bins[i]), which(AB_data[,1]<=age_bins[i+1]) ) 
	temp  <- AB_data[index,2]

	GMT_bins[i] <- exp(mean(log(temp)))

	AB_range_bins[i,] <- quantile( temp, prob=c(0.5, 0.025, 0.975) )
}


#############################################
## 3.2 Posterior prediction intervals 

N_sam = 100
sam_seq = round(seq(from=1, to=nrow(MCMC_burn), length=N_sam))


M3_sam_cs1 = matrix(NA, nrow=N_sam, ncol=length(age_seq))
for(k in 1:N_sam)
{
	M3_sam_cs1[k,] = sapply(age_seq, model_M3, t_survey=15, par=MCMC_burn[sam_seq[k],1:4])
}

M3_quant_cs1 = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
	M3_quant_cs1[,j] = quantile( M3_sam_cs1[,j], prob=c(0.025, 0.5, 0.975) )
}




###############################################
## 3.3 Plot data and model prediction
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=GMT_bins, 
pch=15, cex=2,
log="y", xlim=c(0,90), ylim=c(0.1,200),
xlab="age (years)", ylab="Geometric mean antibody titre", 
main="Antibody acquisition Model 3 - Rombo CT694"  )

	
for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=AB_range_bins[i,2], 
             x1=age_bins_mid[i], y1=AB_range_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}



points(x=age_seq, y=M3_predict_cs1, 
type='l', lwd=3, col="green")


polygon(x=c(age_seq, rev(age_seq)), 
y=c( M3_quant_cs1[1,], rev(M3_quant_cs1[3,]) ),
col=rgb(144/256,144/256,144/256,0.4), border=NA)



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

theta_bar = apply( X=MCMC_burn[,1:4], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M3_total( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*mean( MCMC_burn[,5] - MCMC_burn[,6] )


######################################
## Effective number of model parameters
	
pD = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC = pD + D_bar


######################################
## Mean of the posterior

theta_bar = apply( X=MCMC_burn[,1:4], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M3_total( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar_median = -2*median( MCMC_burn[,5] - MCMC_burn[,6] )


######################################
## Effective number of model parameters
	
pD_median = D_bar_median - D_theta_bar


######################################
## Estimate of DIC

DIC = pD_median + D_bar_median



