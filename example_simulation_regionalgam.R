#################################################################################
#
# This is a working example to generate simulated butterfly count data, impose
# specific degradation level and compute the regional GAM abundance indices using the 
# regionalGAM method described in our paper subimitted for publication in Journal of Applied Ecology
# 
# Title:
# Regionally informed abundance index for supporting integrative analyses across butterfly monitoring schemes
# 
# Authors:
# Schmucki, Reto; Pe'er, Guy; Roy, David; Stefanescu, Constanti; Van Swaay, Chris; Oliver, Tom; 
# Kuussaari, Mikko; Van Strien, Arco; Ries, Leslie; Settele, Josef; Musche, Martin; Carnicer, 
# Jofre; Schweiger, Oliver; Brereton, Tom; Heliölä, Janne; Harpke, Alexander; Kühn, Elisabeth; 
# Julliard, Romain
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
# and use function available in "functions_regionalgam_simulation.R"
#
#
# NOTE: latest version of this script available for download from my github repository
#
# => https://github.com/RetoSchmucki/regionalGAMpaper.git
#
######################################################################################

    setwd("PATH TO functions_regionalgam_simulation.R) 

    # load the necessary function available in "functions_regionalgam_simulation.R"
    source("functions_regionalgam_simulation.R")

    # Generate data for a univoltine species over 200 sites (nsite x miter) and two years
    new.data <- data_gen(miter=2,nsite=100,long_trend=-0.1,nyear=2,sd_peak_shift=2,bivoltine='FALSE')

       # Compute flight curve the each year
       pheno_full <- flight_curve(new.data)

       # plot the flight curve
       plot(pheno_full$DAYNO[pheno_full$year==1],pheno_full$nm[pheno_full$year==1],type="l", col='red', xlab="day",ylab="relative abundance")
       points(pheno_full$DAYNO[pheno_full$year==2],pheno_full$nm[pheno_full$year==2],type="l", col='blue')
       legend('topright',legend=c("year 1","year 2"),fill=c("red","blue"),bty='n')


    # Generate degradation pattern to create a dataset with no missing observation per monitoring season
    # with no contraint on to keep the week of peak abundance  
    degra <- degradation_prop(new.data, pheno = pheno_full, perct = 1, fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE)

    # Apply the degradation pattern over the newly created data and compute abundance indices
       index_full <- data_degradation(new.data, degra = degra, perct = 1)

        # Example with 40% missing with a constraint to keep peak weeks in the degradation process
        degra <- degradation_prop(new.data, pheno = pheno_full, perct = 0.3,fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE)
		
        index_degrad60 <- data_degradation(new.data, degra = degra, perct = 0.3)

        # Compute percent error for the degraded dataset (60%)
        perct_error06rg <- data.frame(degrad=0.6,perct_pheno= index_degrad60$prop_pheno_sampled,pct_error=(index_degrad60$regional_gam-index_full$regional_gam)/index_full$regional_gam)

        # produce an histogram
        hist(perct_error06rg$pct_error)
