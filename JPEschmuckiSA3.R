#################################################################################
#
#	Appendix S2
#	
#	- PLEASE READ CAREFULLY -
#
# This script contains two worked examples:
#
# 1. showing how to generate simulated butterfly count data, impose a
# specific degradation level, and compute the regional GAM abundance indices using the 
# regionalGAM method presented in:
# 
# Schmucki, R., Pe’er, G., Roy, D. b., Stefanescu, C., Van Swaay, C. a. m., Oliver, T. h., Kuussaari, M.,
# Van Strien, A. j., Ries, L., Settele, J., Musche, M., Carnicer, J., Schweiger, O., Brereton, T. m., 
# Harpke, A., Heliölä, J., Kühn, E. & Julliard, R. (2015) A regionally informed abundance index for supporting 
# integrative analyses across butterfly monitoring schemes. Journal of Applied Ecology. DOI: 10.1111/1365-2664.12561
#
# 2. Installation and usage of the R package regionalGAM to compute abundance indices with the 
# regional GAM method, using the dataset collected for Gatekeeper Pyronia tithonus (Linnaeus 1767). 
# These data are available in Appendix S4 in Supporting Information and from the R package RegionalGAM 
# that can be installed from the GitHub repository https://github.com/RetoSchmucki/regionalGAM.
#
#
# ==================================================================================
#
# Version 1.0
# R version 3.2.1
# 
# Requiere packages
# 	- mgcv
#	- sampling
#
# This script was developed by Reto Schmucki - reto.schmucki[at]mail.mcgill.ca
#
# The two-stage model was first presented in: 
# Dennis, E.B., Freeman, S.N., Brereton, T. & Roy, D.B. (2013) Indexing butterfly abundance 
# whilst accounting for missing counts and variability in seasonal pattern. 
# Methods in Ecology and Evolution, 4, 637–645.
#
#
# NOTE: The latest version of this script is available for download from my GitHub repository
#
# => https://github.com/RetoSchmucki/regionalGAMpaper.git
#
# The functions necessary to compute the regional GAM index are part of an R package that 
# is available from installation my GitHub repository https://github.com/RetoSchmucki/regionalGAM
#
# For installation and usage, please refer to the instruction available in the README.md document
# that is available here: https://github.com/RetoSchmucki/regionalGAM/tree/master
#
######################################################################################


################################################################
#
# Example 1. Data Simulation and Degradation,
#            using the functions provided in Appendix S2 (JPEschmuckiSA2.R) 
#			 available in online supporting information.
#
################################################################

    setwd("PATH TO /JPEschmuckiSA2.R")
    
    # load the necessary function found in "JPEschmuckiSA2.R"
    source("JPEschmuckiSA2.R")

	################################################################################################################
	#
	# Generate data for a univoltine species over 100 sites (nsite x miter) and two years with a 10% decline
	#
	# - long_trend : the trend over the full period defined by nyear.
	# - nyear : the number of year to simulate.
	# - sd_peak_shift : standard deviation allow around the peak week.
	# - bivoltine : if TRUE, generate a flight curve with two peaks, if FALSE a univoltine flight cure is generated
	#
	################################################################################################################   
    
    new.data <- data_gen(miter=2,nsite=50,long_trend=-0.1,nyear=2,sd_peak_shift=2,bivoltine='FALSE')

    # Compute flight curve the each year
    pheno_full <- flight_curve(new.data)

    # plot the flight curve
    plot(pheno_full$DAYNO[pheno_full$year==1],pheno_full$nm[pheno_full$year==1],type="l", col='red', xlab="day",ylab="relative abundance")
    points(pheno_full$DAYNO[pheno_full$year==2],pheno_full$nm[pheno_full$year==2],type="l", col='blue')
    legend('topright',legend=c("year 1","year 2"),fill=c("red","blue"),bty='n')



	########################################################################################################
	#
	# Generate degradation patterns to create a dataset with specific proportion of missing observations per monitoring season
	# with or without applying a constraint on peak week.
	#
	# - perc : define the proportion of week sampled (to keep)
	# - fully_random_subsample : if TRUE, no constraint are imposed on the peak week, only stratification.
	#
	#	IF fully_random_subsample is FALSE, then:
	#
	# - keep_peak_in_subsample : if TRUE, peak week is kept in the series of sampled week.
	# - keep_peak_in_subsample : if FALSE, peak week is systematically excluded from the sampled week.
	#
	########################################################################################################
     
    degra <- degradation_prop(new.data, pheno = pheno_full, perct = 1, fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE)

    # Apply the degradation pattern over the newly created data and compute abundance indices
    index_full <- data_degradation(new.data, degra = degra, perct = 1)

    # Example with 40% missing [perct = 0.6] with the constraint to keep peak weeks in the degradation process
    degra <- degradation_prop(new.data, pheno = pheno_full, perct = 0.6,fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE)
	index_degrad60 <- data_degradation(new.data, degra = degra, perct = 0.6)

    # Compute percent error for the degraded dataset, 40% missing [perct = 0.6]
    perct_error06rg <- data.frame(degrad=0.6,perct_pheno= index_degrad60$prop_pheno_sampled,
	pct_error=(index_degrad60$regional_gam-index_full$regional_gam)/index_full$regional_gam)

    # produce an histogram
    hist(perct_error06rg$pct_error)
        
      
################################################################
#
# Example 2. Install the R package "regionalGAM" and compute
#            abundance indices for the Gatekeeper count data
#
################################################################
	
	# clean working environment
	rm(list=ls())
	
	# for installing the packages from github, you will need to install the devtools package
	install.packages("devtools")

	# install and load the regionalGAM package
	library(devtools)
	install_github("RetoSchmucki/regionalGAM")
	library(RegionalGAM)


	# load count data for the gatekeeper in the Cold temperate and Moist region
	data("gatekeeper_CM")
	head(gatekeeper_CM)
	
	# format data to compute the flight curve
	dataset1 <- gatekeeper_CM[,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]

	# compute the annual flight curve for the regional dataset
	# WARNING, this may take some time.
	pheno <- flight_curve(dataset1)
	
	# plot pheno for year 2005
	plot(pheno$DAYNO[pheno$year==2005],pheno$nm[pheno$year==2005],pch=19,cex=0.7,type='o',col='red',xlab="day",ylab="relative abundance")

	# format data to compute abundance indices
		# Note: here we restrict index computation to the sites used in the trend analysis in the paper entitled,
		# A Regionally informed abundance index for supporting integrative analyses across butterfly monitoring schemes
		# But this can be extended to all sites
	dataset2 <- gatekeeper_CM[gatekeeper_CM$TREND==1,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]
	
	# compute the annual abundance indices
	data.index <- abundance_index(dataset2, pheno)

	# save your results (pheno and data.index
	write.csv(pheno, file="flight_pheno.csv",row.names=FALSE)
	write.csv(data.index, file="abundance_indices.csv",row.names=FALSE)

###############################################################################################################
# Temporal Trend Estimate
# Trend can also be estimated in TRIM
# for TRIM, refer to http://www.cbs.nl/en-GB/menu/themas/natuur-milieu/methoden/trim/default.htm?Languageswitch=on
###############################################################################################################

	# load packages used for trend estimate
	library(nlme)
	library(MASS)

	# compute collated annual indices
	glmm.mod_fullyear <- glmmPQL(regional_gam~ as.factor(YEAR)-1,data=data.index,family=quasipoisson,random=~1|SITE, correlation = corAR1(form = ~ YEAR | SITE),verbose = FALSE)
	summary(glmm.mod_fullyear)

	# extract values and plot
	year <- unique(data.index$YEAR)
	col.index <- as.numeric(glmm.mod_fullyear$coefficients$fixed)
	plot(year,col.index,type='o', xlab="year",ylab="collated index")

	# model temporal trend with a simple linear regression
	mod1 <- gls(col.index ~ year)
	summary(mod1)

	# check for temporal autocorrelation in the residuals
	acf(residuals(mod1,type="normalized"))

	# adjust the model to account for autocorrelation in the residuals
	mod2 <- gls(col.index ~ year, correlation = corARMA(p=2))
	summary(mod2)
	
	# check for remaining autocorrelation in the residuals
	acf(residuals(mod2,type="normalized"))

	# plot abundance with trend line
	plot(year,col.index, type='o',xlab="year",ylab="collated index")
	abline(mod2,lty=2,col="red")

	# NOTE: bootstrap procedure can be used to estimate confidence intervals around collated indices.
