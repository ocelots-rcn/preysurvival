# R script for Modeling Exercise: Prey Clumping & Density Dependent Survival
# by Sabrina Russo (srusso2@unl.edu) & Glenn Ledder, University of Nebraska - Lincoln

# This modeling exercise is meant for students in an upper-level undergraduate course who do not have very much math background
# For a published study corresponding to these modeling exercises, 
#	see Russo & Augspurger 2004, Ecology Letters &  Russo et al. 2006, Ecology

# computes seed survival from granivores with seeds normally-dispersed in the x direction

# sigma is the standard deviation for the normal pdf, which is used as the seed dispersal curve

# h is the handling time of the granivore consuming seeds
# a parameter for positive density-dependence caused by limitations to how fast the granivores can eat seeds
# h>0 means seed survival is enhanced at high population densities because the granivores need time to digest
# This effect is increased for larger h.

# tf is the total dimensionless time for the simulation
# As a rough guide, note that there is about 60% survival with sigma=1, h=0, tf=2

# xmax is the maximum distance for calculations
# It should be no smaller than 2.8*sigma, but should not be excessive

## The core of this script is the function survivors, which outputs a list of quantities of interest.
## Later functions produce graphs depicting relationships to illustrate the combined effects of 
#	positive density dependence caused by handling time and 
#	negative density dependence caused by seed crowding (due to limited seed dispersal)

# run the function using the command 
#   ans <- survivors(sigma,h,tf,xmax)
# and then unpack the desired results

# S <- ans[1] is the total survival at time tf
# x <- unlist(ans[2]) is a list of x coordinates needed for spatial plots
# s <- unlist(ans[3]) is a list of survival probabilities keyed to x
# p <- unlist(ans[4]) is a list of population density function values, used to plot the final seed distribution

## SET SCENARIO PARAMETERS

#sigma = 1		# sigma must be < 10 and > 0.1 or else the curves below get funky!
#h = 0
#tf = 0:4
#xmax = 2.8*sigma

require(pracma)

####################################################################################################################

### DO NOT CHANGE THE FOLLOWING FUNCTION CODE!!
## DEFINE FUNCTION TO CALCULATE SURVIVAL PROBABILITIES s(x) and S


survivors.S = function(sigma,h,tf)
{
  xmax = 2.8*sigma
  ss = sqrt(2)*sigma
  sp = ss*sqrt(pi)
  H = h/sp
  C = (tf-h)/sp
  umax = xmax/ss
  du = .02*umax
  u = seq(0,umax,by=du)
  u2 = u^2
  G = u2+C*exp(-u2)
  y = exp(-u2-C*exp(-u2))
  if (H>0)
    for (k in 1:51)
    {
      range <- c(0,y[k])
      ans <- uniroot(function(Y) log(Y)+H*Y+G[k], range)
      y[k] <- ans$root
    }
  n <- 1:51
  alt <- cos(n*pi)
  odd <- 1-alt
  even <- 1+alt
  trap <- (odd %*% y-y[1]-y[51])*du
  mid <- (even %*% y)*du
  int <- (2*mid+trap)/3
  S <- 2*int/sqrt(pi)
  x <- ss*u
  s <- y*exp(u2)
  p <- y/sp
  return(as.numeric(S))
}



survivors = function(sigma,h,tf,xmax)
{
  ss = sqrt(2)*sigma
  sp = ss*sqrt(pi)
  H = h/sp
  C = (tf-h)/sp
  umax = xmax/ss
  du = .02*umax
  u = seq(0,umax,by=du)
  u2 = u^2
  G = u2+C*exp(-u2)
  y = exp(-u2-C*exp(-u2))
  if (H>0)
    for (k in 1:51)
    {
      range <- c(0,y[k])
      ans <- uniroot(function(Y) log(Y)+H*Y+G[k], range)
      y[k] <- ans$root
    }
  n <- 1:51
  alt <- cos(n*pi)
  odd <- 1-alt
  even <- 1+alt
  trap <- (odd %*% y-y[1]-y[51])*du
  mid <- (even %*% y)*du
  int <- (2*mid+trap)/3
  S <- 2*int/sqrt(pi)
  x <- ss*u
  s <- y*exp(u2)
  p <- y/sp
  return(list(S=S,x=x,s=s,p=p))
}

### DO NOT CHANGE THE PRECEDING FUNCTION CODE!!

# run the function by assigning it to a variable name:
# S is ans[1]
# x is unlist(ans[2])
# s(x) is unlist(ans[3])
# p(x) is unlist(ans[4])

# ans = survivors(sigma,h,tf[2])


####################################################################################################################


# Student worksheet Question #1
# Seed dispersal curve based on normal distribution (used above)
#	and sigma set as above
#	mean = 0, the location of the mother tree
# Note that sigma < 4 gives strange plots for the seed density distribution in the next graph (Question 2)


plot.curve = function(sigmas=c(4, 6, 8, 10), mean=0, max.x=NULL, max.y=NULL)  	{

	colors.use = colorRampPalette(c("blue", "red"))(length(sigmas))

	if (is.null(max.x)) {
	  max.x <- 2.8*max(sigmas)+mean
	}
	
	if (is.null(max.y)) {
	  if (min(sigmas) < 0.1) max.y = 4 
	  if (min(sigmas) >= 4) max.y = 0.12 else max.y = 1
	}

	curve(dnorm(x, mean=0, sd=sigmas[1]), from = -1*max.x, to = max.x,
		col=colors.use[1], lwd = 2, ylim=c(0,max.y), type="l", 
		main="Seed dispersal curve with varying sigma", xlab="Seed Dispersal Distance from Mother Tree", 
		ylab="Probability of Seed at Distance")

	if (length(sigmas) > 1 ) {
		for (i in 2:length(sigmas)) {

		curve(dnorm(x, mean=0, sd=sigmas[i]), from = -2.8*max(sigmas)+mean, to = 2.8*max(sigmas)+mean, 
			col=colors.use[i], lwd = 2, type="l", add=T) 
		}
	}
	
	text(x=15, y = seq(0.75*max.y, max.y, length.out=length(sigmas)), labels = paste("sigma =", sigmas), 
		col=colors.use, pos=4)

 }


# Example Usage: plot.curve(sigmas= c(4, 6, 8, 10))


####################################################################################################################


# Student worksheet Question #2
# plot and return distribution of seed densities for a thousand seeds dispersed 
#	according to seed dispersal curves with varying sigmas
# note that there is no randomness at all in this function
# xx are evenly-distributed positions along the maximum possible x-axis range (given the largest sigma)
# xx are fed to dnorm, which simply uses the formula for the probability density function of a normal distribution
#	to find the y-value (probability)
# then trapz takes the area under the curve and gives that for each bin (AUC), since that corresponds to the number of seeds (integration of the pdf)
#	and find the number of seeds per bin (AUC*no.seeds)
# plot Empirical PDFs http://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
# note that the bins are not quite correct because the last bin is smaller than 2 units.  


dist.seed.dens = function(sigmas, mean=0, adj=1)  	{

	require(pracma)
	is.odd <- function(x) x %% 2 != 0 
	blue.col <- rgb(0,0,1,0.3)

	no.seeds = 1e3
	dens.by.sigma =list()
	Pdens.by.sigma=list()
	max.dens = rep(NA, length(sigmas))
	mean.by.sigma = SD.by.sigma = rep(NA, length(sigmas))
	xlim.max.by.sigma = ylim.max.by.sigma = rep(NA, length(sigmas))

	xx = seq(-2.8*max(sigmas)+mean, 2.8*max(sigmas)+mean, length.out=1000)
	bins = seq(-2.8*max(sigmas)+mean, 2.8*max(sigmas)+mean, by=2)
	bins=c(bins, 1000)
	xx.binned = cut(x=xx, breaks=bins, include.lowest=T, right=F)

	AUC = rep(NA, (length(bins)-1))
	no.seeds.bin  = rep(NA, (length(bins)-1))

	for (i in 1:length(sigmas))  {
		hx = dnorm(xx, mean=0, sd=sigmas[i])

		for (j in 1:(length(bins)-1))  {

			xx.bin = xx[xx.binned==levels(xx.binned)[j]]
			hx.bin = hx[xx.binned==levels(xx.binned)[j]]

			AUC[j] = trapz(xx.bin, hx.bin)
			no.seeds.bin[j] = AUC[j]*no.seeds

			}

		dens.by.sigma[[i]] = round(no.seeds.bin,2)
		SD.by.sigma[i] = sd(no.seeds.bin)
		mean.by.sigma[i] = mean(no.seeds.bin)


# for histogram version
#		max.dens[i] = max(no.seeds.bin, na.rm=T)
#		if (i==1) max.to.use = ceiling(max.dens[i])
#		if (i > 1) max.to.use = ceiling(max(max.dens, na.rm=T))
#		if ( is.odd(max.to.use) )  max.to.use = max.to.use+1
#		hist(no.seeds.bin, breaks=seq(0, max.to.use, by=2), xlab = "Seed Density", ylab = "No. Patches", 
#			main=paste("sigma (dispersal curve) =", sigmas[i], "\n", "Std. deviation of seed density =", ceiling(SD.by.sigma[i])), 
#			ylim=c(0,26), xlim=c(0,max.to.use), cex.main=1)

		# calculate probability density

		dens = density(no.seeds.bin, from=0, adjust=adj)

		Pdens.by.sigma[[i]] <- dens
		xlim.max.by.sigma[i] <- max(dens$x)
		ylim.max.by.sigma[i] <- max(dens$y)

		}

	x.lim = max(xlim.max.by.sigma)
	y.lim = max(ylim.max.by.sigma)

	par(mfcol = c(2, length(sigmas)/2) )

	for (i in 1:length(sigmas))  {

		plot(Pdens.by.sigma[[i]], xlim = c(0,x.lim), ylim=c(0,y.lim), xlab = "Seed Density in Patch", ylab = "Probability",
		    main=paste("sigma (dispersal curve) =", sigmas[i], "\n", "Std. deviation of seed density =", ceiling(SD.by.sigma[i])),
		    cex.lab=1.1, cex.axis=1.0, cex.main=1.2, type="l", lwd=2, col="blue")

		# put our density plots in

		polygon(x=c(0,Pdens.by.sigma[[i]]$x), y=c(0,Pdens.by.sigma[[i]]$y), density = -1, col = blue.col)


		}

	return(list(dens.by.sigma = dens.by.sigma, SD.by.sigma = SD.by.sigma, mean.by.sigma=mean.by.sigma))

  }



# Example Usage: dens.by.sigma = dist.seed.dens(sigmas=c(4, 6, 8, 10))
# to get histogram (replace 4 with whichever sigma you want to look at):
	# sab2=hist(sab$dens.by.sigma[[4]], breaks=brks)


####################################################################################################################


# Student worksheet Question #3
# Holling Type II functional response curve with different values of handling time h


plot.Holl.T2 = function(hs=c(0, 0.5, 1, 2, 3, 4))  	{

	colors.use = colorRampPalette(c("blue", "red"))(length(hs))

	curve(x/(1 + hs[1]*x), from = 0, to = 2, col=colors.use[1], lwd = 2, type="l", 
		main="Holling Type II curve with varying h", xlab="Seed Density in the Patch", 
		ylab="Consumption Rate of Seeds in the Patch")

	if (length(hs) > 1 ) {
		for (i in 2:length(hs)) {

		curve(x/(1 + hs[i]*x), from = 0, to = 2, col=colors.use[i], lwd = 2, type="l", add=T) 
		}

	}
	
	space=0.1

	text(x=0, y = seq(2-(space*length(hs)), 2, length.out=length(hs)), labels = paste("h =", hs), 
		col=colors.use, pos=4)

 }

# Example Usage: plot.Holl.T2(hs=c(0, 0.5, 1, 2, 3, 4))


####################################################################################################################


# Student worksheet Question #5
# plot S as a function of sigma for various values of h
# S is the total proportion of population surviving
# sigmas must be vector; hs can be a vector; tf is a scalar


S.sigma.curve = function(hs = 0:4, sigmas = seq(0.05,5,by=0.05), tf=2)  {

	J <- length(hs)
	K <- length(sigmas)
	S <- rep(0,K)

	xmaxs = 2.8*sigmas

	colors.use = colorRampPalette(c("blue", "red"))(J)

	for (k in 1:K)
		{
		ans <- survivors(sigmas[k], hs[1], tf, xmaxs[k])
		S[k] <- ans[1]
		}
	plot(sigmas, S, type='l', ylim=c(0,1), col=colors.use[1], lwd=2, xlab = "Standard deviation of seed dispersal curve (sigma)",
	  ylab = "Total proportion of seeds in the population surviving (S)",
	  main="Proportion of all seeds surviving with respect to sigma with varying h", cex.main=0.9)

	for (j in 2:J)
		{
		for (k in 1:K)
			{
			ans <- survivors(sigmas[k], hs[j], tf, xmaxs[k])
			S[k] <- ans[1]
			}
		lines(sigmas, S, col=colors.use[j], lwd=2)

		text(x=4, y = seq(0.3-(0.05*J), 0.3, length.out=J), labels = paste("h =", hs), col=colors.use, pos=4)
	}

  }

# Example Usage: S.sigma.curve(hs = 0:4, sigmas = seq(0.05,5,by=0.05))


# Student worksheet Question #5
# plot S as a function of h for various values of sigmas
# S is the total proportion of population surviving
# hs must be vector; sigmas can be a vector; tf is a scalar


S.h.curve = function(hs = seq(0, 4, by=0.05), sigmas =c(0.1, 1:5), tf=2)  {

	J <- length(sigmas)
	K <- length(hs)
	S <- rep(0,K)

	xmaxs = 2.8*sigmas

	colors.use = colorRampPalette(c("blue", "red"))(J)

	for (k in 1:K)	# loop for x-axis
		{
		ans <- survivors(sigmas[1], hs[k], tf, xmaxs[1])
		S[k] <- ans[1]
		}

	plot(hs, S, type='l', ylim=c(0,1), col=colors.use[1], lwd=2, xlab = "Handling time (h)",
	  ylab = "Total proportion of seeds in the population surviving (S)",
	  main="Proportion of all seeds surviving with respect to h with varying sigma", cex.main=0.9)

	for (j in 2:J)
		{
		for (k in 1:K)
			{
			ans <- survivors(sigmas[j], hs[k], tf, xmaxs[j])
			S[k] <- ans[1]
			}
		lines(hs, S, col=colors.use[j], lwd=2)

		text(x=3, y = seq(0.3-(0.05*J), 0.3, length.out=J), labels = paste("sigma =", sigmas), col=colors.use, pos=4)
	}

  }

# Example Usage: S.h.curve(hs = seq(0, 4, by=0.05), sigmas =c(0.1, 1:5))


####################################################################################################################


