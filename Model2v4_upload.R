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

setwd("")

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

data <- data[data$province == ,]			## if there are multiple study sites in the data set subdivide them


PGP3_pre_bin <- rep(NA, (dim(data)[1]-1))	## define an empty vector to store whether someone is sero-pos or not


PGP3_cut <- 0.80							## specify the cut-off vlaue which defines people as +ve or -ve

for(i in 1:nrow(data))
{
	if(data$OD[i] > PGP3_cut)				## loop through the data to fill in the vector to define people as +ve or -ve
	{
		PGP3_pre_bin[i] <- 1
	}else{
		PGP3_pre_bin[i] <- 0 
	}

}




Kin_pre_binned <- cbind(data$age, PGP3_pre_bin)	## make a dataframe of all individuals age and sero-status for analysis 



colnames(Kin_pre_binned)[1] <- "age"
colnames(Kin_pre_binned)[2] <- "PGP3_bin"


AB_data <- Kin_pre_binned						## assign the object to a new name that will be evalauated by the function 
SP_cs <- Kin_pre_binned

###############################################
## 0.2 bin the data by age group to plot


age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)						## age bins will be set from the youngest age group 
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] ) ## in your data to the maximum age by = user defined
 
 
N_bins <- length(age_bins) - 1 


SP_bins      <- rep(NA, N_bins)

SP_range_bins_cs <- matrix(NA, nrow=N_bins, ncol=3)
colnames(SP_range_bins_cs) <- c("med", "low", "high")		## generate an empty matrix to be filled with sero-prev for each age bin and the binomial CIs

for(i in 1:N_bins)
{
	index <- intersect( which(Kin_pre_binned[,1]>age_bins[i]), which(Kin_pre_binned[,1]<=age_bins[i+1]) ) 
	temp  <- Kin_pre_binned[index,2]

	SP_bins[i] <- sum(Kin_pre_binned[index,2])/length(index)

	SP_range_bins_cs[i,] <- qbinom( c(0.5,0.025,0.975), size=length(index), prob=sum(Kin_pre_binned[index,2])/length(index) )/length(index)
}


###############################################
## 0.3 Plot data
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=SP_bins, 
pch=15, cex=2,
xlim=c(0,90), ylim=c(0,1),
xlab="age (years)", ylab="Proportion seropositive", 
main="Cross-sectional serological data"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=SP_range_bins_cs[i,2], 
             x1=age_bins_mid[i], y1=SP_range_bins_cs[i,3], 
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
## 1.1 specificy the serological model which may
## have generated our sero-prev data 

rho = 0.02611204 ## defined the fix value of the SRR

model_M2_STI = function(a, t_survey, par_M2)
{
	lambda_0   = par_M2[1]
	gamma      = par_M2[2]
	time_c     = par_M2[3]
	lambda_STI = par_M2[4]

	lambda_c = gamma*lambda_0

	age_xx = a - time_c + t_survey
	

	if( (age_SD <= age_xx) && (age_xx <= a) )
	{
		SP_prop = ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*age_SD) )
		
		SP_prop = SP_prop*exp(-(lambda_0+lambda_STI+rho)*(age_xx-age_SD)) + ( (lambda_0+lambda_STI)/(lambda_0+lambda_STI+rho) )*( 1 - exp(-(lambda_0+lambda_STI+rho)*(age_xx-age_SD)) )

		SP_prop = SP_prop*exp(-(lambda_c+lambda_STI+rho)*(a-age_xx)) + ( (lambda_c+lambda_STI)/(lambda_c+lambda_STI+rho) )*( 1 - exp(-(lambda_c+lambda_STI+rho)*(a-age_xx)) )
	}

	if( (age_SD <= a) && (a <= age_xx) )
	{
		SP_prop = ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*age_SD) )
		
		SP_prop = SP_prop*exp(-(lambda_0+lambda_STI+rho)*(a-age_SD)) + ( (lambda_0+lambda_STI)/(lambda_0+lambda_STI+rho) )*( 1 - exp(-(lambda_0+lambda_STI+rho)*(a-age_SD)) )
	}

	if( (age_xx > 0) && (age_xx <= age_SD) && (age_SD <= a) )
	{
		SP_prop = ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*age_xx) )
		
		SP_prop = SP_prop*exp(-(lambda_c+rho)*(age_SD-age_xx)) + ( (lambda_c)/(lambda_c+rho) )*( 1 - exp(-(lambda_c+rho)*(age_SD-age_xx)) )

		SP_prop = SP_prop*exp(-(lambda_c+lambda_STI+rho)*(a-age_SD)) + ( (lambda_c+lambda_STI)/(lambda_c+lambda_STI+rho) )*( 1 - exp(-(lambda_c+lambda_STI+rho)*(a-age_SD)) )
	}

	if( (age_xx > 0) && (age_xx <= a) && (a <= age_SD) )
	{
		SP_prop = ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*age_xx) )
		
		SP_prop = SP_prop*exp(-(lambda_c+rho)*(a-age_xx)) + ( (lambda_c)/(lambda_c+rho) )*( 1 - exp(-(lambda_c+rho)*(a-age_xx)) )
	}	

	if( (age_xx > 0) && (a <= age_xx) && (a <= age_SD) )
	{
		SP_prop = ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*a) )
	}	

	if( (age_xx <= 0) && (a <= age_SD) )
	{
		SP_prop = ( lambda_c/(lambda_c+rho) )*( 1 - exp(-(lambda_c+rho)*a) )
	}	

	if( (age_xx <= 0) && (age_SD <= a) )
	{
		SP_prop = ( lambda_c/(lambda_c+rho) )*( 1 - exp(-(lambda_c+rho)*age_SD) )

		SP_prop = SP_prop*exp(-(lambda_c+lambda_STI+rho)*(a-age_SD)) + ( (lambda_c+lambda_STI)/(lambda_c+lambda_STI+rho) )*( 1 - exp(-(lambda_c+lambda_STI+rho)*(a-age_SD)) )
	}	

	SP_prop
}

model_M2_STI = cmpfun(model_M2_STI, options=list(optimize=3)) ## wrap the transmission model into a function


###################################################
## 1.2 calculate the likelihood of the paramater set
## using a binomial likelihhod


loglike_M2_STI_cs = function( par_M2 )
{
	SP_model = sapply(SP_cs[,1], model_M2_STI, t_survey=0, par_M2=par_M2)

	loglike = SP_cs[,2]*log(SP_model) + (1-SP_cs[,2])*log(1-SP_model)

      sum( loglike )
}

loglike_M2_STI_cs = cmpfun(loglike_M2_STI_cs, options=list(optimize=3)) 



###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
LARGE = 1e10     ## Large value for rejecting parameters with prior

prior_M2 = function( par_M2 )
{
	lambda_0   = par_M2[1]
	gamma      = par_M2[2]
	time_c     = par_M2[3]
	lambda_STI = par_M2[4]


	######################################
	## Uniform prior on lambda_0 ~ U(0,10)

	if( lambda_0>0 && lambda_0<10 )
	{
		prior_lambda_0 = log(1/10)
	}else{
		prior_lambda_0 = -LARGE
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
	## Uniform prior on time_c ~ U(0,60)

	if( time_c>0 && time_c<60 )
	{
		prior_time_c = log(1/60)
	}else{
		prior_time_c = -LARGE
	}


	######################################
	## Uniform prior on lambda_STI ~ U(0,10)

	if( lambda_STI>0 && lambda_STI<10 )
	{
		prior_lambda_STI = log(1/10)
	}else{
		prior_lambda_STI = -LARGE
	}


	prior = prior_lambda_0 + prior_gamma + prior_time_c + prior_lambda_STI

	prior
}

prior_M2 = cmpfun(prior_M2, options=list(optimize=3))


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


N_mcmc       = 30000      ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate



#################################################
## 2.1 Robbins-munro step scaler


rm_scale = function(step_scale, mc, log_prob)
{
	dd = exp(log_prob)
	if( dd < -30 ){ dd = 0 }
	dd = min( dd, 1 )

	rm_temp = ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
		
	out = step_scale*exp(rm_temp)
	
	out = max( out, 0.02 )
	out = min( out, 2)
	out
}


#################################################
## 2.2 Prepare object to store MCMC fitting output
 
MCMC_par           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par) = c("lambda_0", "gamma", "time_c", "lambda_STI", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations


par_MC = c(0.1, 0.5, 20, 0.2)       ## (lambda_0, gamma, rho, time_c, lambda_STI)


Sigma_MC = diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist

prior_MC = prior_M2( par_MC )             

loglike_MC = loglike_M2_STI_cs( par_MC ) + prior_MC



for(mc in 1:N_mcmc)
{
	par_MCp1 = mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)

	prior_MCp1 = prior_M2( par_MCp1 )
 
	if( prior_MCp1 > -0.5*LARGE )
	{
 		loglike_MCp1 = loglike_M2_STI_cs( par_MCp1 ) + prior_MCp1


		log_prob = min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC = par_MCp1
			
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

	
	MCMC_par[mc,1:4] = par_MC			## store the output parameter values from each MCMC realisation 
	MCMC_par[mc,5]   = loglike_MC		## return the associated loglikelihood
	MCMC_par[mc,6]   = prior_MC			## return the prior
}



#########################################################
## 2.4 Examine MCMC chains
 

par(mfrow=c(2,4))

for(k in 1:4)
{
	#####################################
	## PANEL k

	plot(x=1:N_mcmc, y=MCMC_par[,k], 
	pch=19, col="grey", cex=0.25,
	xlab="MCMC iteration", ylab=colnames(MCMC_par)[k], 
	main=colnames(MCMC_par)[k])
}


#####################################
## PANEL 6: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,6], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
ylim=quantile( MCMC_par[,6], prob=c(0.01,1)),
main="likelihood" )


#####################################
## PANEL 7: prior

plot(x=1:N_mcmc, y=MCMC_par[,7], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="prior", 
ylim=quantile( MCMC_par[,7], prob=c(0.01,1)),
main="prior" )



#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,3))

for(k in 1:4)
{
	#####################################
	## PANEL k: MCMC posterior

	DEN = density( MCMC_burn[,k] )
	
	QUANT = quantile( MCMC_burn[,k], prob=c(0.025, 0.5, 0.975) )

	plot(x=DEN$x, y=DEN$y, type='l',
	xlim=c(0, max(DEN$x)),
	xlab="", ylab="", 
	main=colnames(MCMC_burn)[k] )

	
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

M2_STI_predict_cs = sapply(age_seq, model_M2_STI, t_survey=0, par=par_median)


#############################################
## 3.2 Posterior prediction intervals 


N_sam = 500
sam_seq = round(seq(from=1, to=nrow(MCMC_burn), length=N_sam))



M2_STI_sam_cs = matrix(NA, nrow=N_sam, ncol=length(age_seq))
for(k in 1:N_sam)
{
	M2_STI_sam_cs[k,] = sapply(age_seq, model_M2_STI, t_survey=0, par=MCMC_burn[sam_seq[k],1:5])
}

M2_STI_quant_cs = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
	M2_STI_quant_cs[,j] = quantile( M2_STI_sam_cs[,j], prob=c(0.025, 0.5, 0.975) )
}



###############################################
## 3.3 Plot data and model prediction
 
par(mfrow=c(1,1))

plot(x=age_bins_mid, y=SP_range_bins_cs[,1], 
pch=15, cex=2,
xlim=c(0,90), ylim=c(0,1),
xlab="age (years)", ylab="Proportion seropositive", 
main="Cross-sectional serological data"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=SP_range_bins_cs[i,2], 
             x1=age_bins_mid[i], y1=SP_range_bins_cs[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}

points(x=age_seq, y=M2_STI_predict_cs, 
type='l', lwd=3, col="black")

polygon(x=c(age_seq, rev(age_seq)), 
y=c( M2_STI_quant_cs[1,], rev(M2_STI_quant_cs[3,]) ),
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
## We also calacuate the DIC using the median
## instead of the mean as orginally done

######################################
## Mean of the posterior

theta_bar = apply( X=MCMC_burn[,1:4], FUN=mean, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2_STI_cs( theta_bar )


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
## Median of the posterior

theta_bar = apply( X=MCMC_burn[,1:4], FUN=median, MARGIN=2)


######################################
## Deviance at mean of the posterior

D_theta_bar = -2*loglike_M2_STI_cs( theta_bar )


######################################
## Mean deviance (averaged over posterior)

D_bar = -2*median( MCMC_burn[,5] - MCMC_burn[,6] )


######################################
## Effective number of model parameters
	
pD_median = D_bar - D_theta_bar


######################################
## Estimate of DIC

DIC_median = pD_median + D_bar

