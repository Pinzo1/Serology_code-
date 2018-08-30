#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White                                                   ##
## Imperial College Lonodon                                        ##
## m.white08@imperial.ac.uk                                        ##
##                                                                 ##
#####################################################################
#####################################################################

#####################################################################
rm(list=ls())

setwd("")

## call the libraries that are needed for analysis

library(MASS)
library(compiler)
library(binom)
 

data <- read.csv(".csv", header = TRUE)					 ## call in the csv datafile that needs to be analysed 

head(data)						 						 ## check that the data has been read in correctly 


data <- data[data$province == "",]						 ## if there are multiple study sites in the data set subdivide them




PGP3_pre_bin <- rep(NA, (dim(data)[1]-1))				 ## define an empty vector to store whether someone is sero-pos or not


PGP3_cut <- 0.80										 ## specify the cut-off vlaue which defines people as +ve or -ve

for(i in 1:nrow(data))									 ## loop through the data to fill in the vector to define people as +ve or -ve
{
	if(data$OD[i] > PGP3_cut)
	{
		PGP3_pre_bin[i] <- 1
	}else{
		PGP3_pre_bin[i] <- 0 
	}

}



Kin_pre_binned <- cbind(data$age, PGP3_pre_bin)			 ## make a dataframe of all individuals age and sero-status for analysis 



colnames(Kin_pre_binned)[1] <- "age"
colnames(Kin_pre_binned)[2] <- "PGP3_bin"

colnames(Kin_pre_binned) 

AB_data <- Kin_pre_binned								## assign the object to a new name that will be evalauated by the function  



###############################################
## 0.2 bin the data by age group to plot
				

age_bins <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90)							## age bins will be set from the youngest age group 
age_bins_mid <- 0.5*( age_bins[1:(length(age_bins)-1)] + age_bins[2:length(age_bins)] )		## in your data to the maximum age by = user defined
 
 
N_bins <- length(age_bins) - 1 


SP_bins      <- rep(NA, N_bins)

SP_range_bins <- matrix(NA, nrow=N_bins, ncol=3)
colnames(SP_range_bins) <- c("med", "low", "high")						## generate an empty matrix to be filled with sero-prev for each age bin and the binomial CIs

for(i in 1:N_bins)
{
	index <- intersect( which(Kin_pre_binned[,1]>age_bins[i]), which(Kin_pre_binned[,1]<=age_bins[i+1]) ) 
	temp  <- Kin_pre_binned[index,2]

	SP_bins[i] <- sum(Kin_pre_binned[index,2])/length(index)

	SP_range_bins[i,] <- qbinom( c(0.5,0.025,0.975), size=length(index), prob=sum(Kin_pre_binned[index,2])/length(index) )/length(index)
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
	arrows(x0=age_bins_mid[i], y0=SP_range_bins[i,2], 
             x1=age_bins_mid[i], y1=SP_range_bins[i,3], 
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

par_MC <- c(5, 0.1, 10)  ## (lambda_0, gamma, time_c)


rho = 0.02611204 			## specify the fixed value of the SRR


model_M2 <- function(a, par)
{
	lambda_0 <- par[1]
	gamma    <- par[2]
	time_c   <- par[3]

	lambda_c <- gamma*lambda_0

	age_xx  <- a - time_c 
 
	if( age_xx<=0 ){
		SP_prop <- ( lambda_c/(lambda_c+rho) )*( 1 - exp(-(lambda_c+rho)*a) )
	}

	if( (age_xx>0) && (age_xx<a) ){
		SP_prop <- lambda_c/(lambda_c+rho) + 
			( (lambda_0-lambda_c)*rho/( (lambda_0+rho)*(lambda_c+rho) ) )*exp( -(lambda_c+rho)*(a-age_xx) ) -
			( lambda_0/(lambda_0+rho) )*exp( -(lambda_c+rho)*a )*exp( -(lambda_0-lambda_c)*age_xx )
	}

	if( a <= age_xx ){
		SP_prop <- ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*a) )
	} 
	
	
	for(k in 1:length(SP_prop))
 	{
 	
 	if(SP_prop[k] > 0.9999999){
 		SP_prop[k] = 0.9999990
 	}
 	
 }


	SP_prop
}


###################################################
## 1.2 calculate the likelihood of the paramater set
## using a binomial likelihhod
## 
 
loglike_M1 <- function( par )
{
	lambda_0 <- par[1]
	gamma    <- par[2]
	time_c   <- par[3]

	SP_model <- sapply(Kin_pre_binned[,1], model_M2, par=par)

	mu <- log(SP_model) 

	loglike <- Kin_pre_binned[,2]*log(SP_model) + (1-Kin_pre_binned[,2])*log(1-SP_model)

      sum( loglike )
}



###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters
 
prior_M1 <- function( par ){
 
	lambda_0 <- par[1]
	gamma    <- par[2]
	time_c   <- par[3]


	######################################
	## Uniform prior on lambda_0 ~ U(0,10)

	if( lambda_0>0 && lambda_0<10 )
	{
		prior_lambda_0 <- 1/10
	}else{
		prior_lambda_0 <- -1e6
	}

	######################################
	## Uniform prior on gamma ~ U(0,1)

	if( gamma>0 && gamma<1 )
	{
		prior_gamma <- 1/1
	}else{
		prior_gamma <- -1e6
	}



	######################################
	## Uniform prior on time_c ~ U(0,60)

	if( time_c>0 && time_c<50 )
	{
		prior_time_c <- 1/50
	}else{
		prior_time_c <- -1e6
	}



	prior <- prior_lambda_0 + prior_gamma +  prior_time_c

	prior
}


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

N_adapt      <- 5000

#################################################
## 2.1 Robbins-munro step scaler


step_scale  <- 1
MCMC_accept <- 0


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
## 2.2 Prepare object to store MCMC fitting output
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par) <- c("lambda_0", "gamma","time_c", "loglike", "prior")

 
#########################################################
## 2.3 Implement MCMC iterations


par_MC <- c(0.2, 0.1, 2)  ## (lambda_0, gamma, rho, time_c)

Sigma_MC <- diag( c(0.01, 0.01, 2) )


## For Sigma_MC Ideally fill in your existing best guess
## based on the estimated posterior from a previous burn-in.
## The most important thing is to make sure that things are
## approximately the correct order of magnitude.



loglike_MC <- loglike_M1( par_MC ) + prior_M1( par_MC )



for(mc in 1:N_mcmc)
{
	par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)

 
	if( par_MCp1[1] > 0 &&
	    par_MCp1[2] > 0 &&
	    par_MCp1[2] < 1 &&
	    par_MCp1[3] > 0 &&
	    par_MCp1[3] < 50  ){
 
		loglike_MCp1_data <- loglike_M1( par_MCp1 ) 
		
		loglike_MCp1 <- loglike_MCp1_data + prior_M1( par_MCp1 )


		log_prob = min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC <- par_MCp1
		
			loglike_MC  <- loglike_MCp1
			MCMC_accept <- MCMC_accept + 1                       
		}

		if( mc < N_adapt ){
			step_scale <- rm_scale( step_scale, mc, log_prob)
		}

		
	}

	MCMC_par[mc,1:3] <- par_MC					## store the output parameter values from each MCMC realisation 
	MCMC_par[mc,4]   <- loglike_MC				## return the associated loglikelihood
	MCMC_par[mc,5]   <- prior_M1( par_MC )		## return the prior

}



#########################################################
## 2.4 Examine MCMC chains
 


par(mfrow=c(2,3))



#####################################
## PANEL 1: alpha_0 MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,1], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="lambda_0", 
main="lambda_0")



#####################################
## PANEL 2: gamma MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,2], 
ylim=c(0,1),
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="gamma", 
main="gamma")


#####################################
## PANEL 4: time_c MCMC chain

plot(x=1:N_mcmc, y=MCMC_par[,3], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="time_c", 
main="time_c" )



#####################################
## PANEL 5: likelihood

plot(x=1:N_mcmc, y=MCMC_par[,4], 
pch=19, col="grey", cex=0.25,
xlab="MCMC iteration", ylab="likelihood", 
main="likelihood" )







#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]


par(mfrow=c(2,2))


#####################################
## PANEL 1: lambda_0 MCMC posterior


DEN <- density( MCMC_burn[,1] )
	
QUANT <- quantile( MCMC_burn[,1], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, max(DEN$x)),
xlab="lambda_0", ylab="", 
main="posterior: lambda_0" )

	
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
## PANEL 3: time_c MCMC posterior


DEN <- density( MCMC_burn[,3] )
	
QUANT <- quantile( MCMC_burn[,3], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN$x, y=DEN$y, type='l',
xlim=c(0, 60),
xlab="time_c", ylab="", 
main="posterior: time_c" )

	
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


par_median <- apply(X=MCMC_burn[,1:3], MARGIN=2, FUN=median)


age_seq <- seq(from=0, to=90, by=5)

M2_predict <- sapply(age_seq, model_M2, par=par_median)




###############################################
## 3.1 Plot data and model prediction
 



par(mfrow=c(1,1))


plot(x=age_bins_mid, y=SP_bins, 
pch=15, cex=2,
xlim=c(0,90), ylim=c(0,1),
xlab="age (years)", ylab="Proportion seropositive", 
main="Sero-catalytic Model 2 fit RenBel"  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=SP_range_bins[i,2], 
             x1=age_bins_mid[i], y1=SP_range_bins[i,3], 
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
## We also calacuate the DIC using the median
## instead of the mean as orginally done

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
######################################
######################################
## Mean of the posterior

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

DIC_median = pD + D_bar















