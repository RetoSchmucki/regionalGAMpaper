#################################################################################
# 
#	Appendix S1
#	
#	- PLEASE READ CAREFULLY -
#
# R-script corresponding to the method presented in: 
#
# Schmucki, R., Pe’er, G., Roy, D. b., Stefanescu, C., Van Swaay, C. a. m., Oliver, T. h., Kuussaari, M.,
# Van Strien, A. j., Ries, L., Settele, J., Musche, M., Carnicer, J., Schweiger, O., Brereton, T. m., 
# Harpke, A., Heliölä, J., Kühn, E. & Julliard, R. (2015) A regionally informed abundance index for supporting 
# integrative analyses across butterfly monitoring schemes. Journal of Applied Ecology. DOI: 10.1111/1365-2664.12561
#
# NOTE: The latest version of this script is available for download from my GitHub repository
#
# => https://github.com/RetoSchmucki/regionalGAMpaper.git
#
# The functions required to compute the regional GAM index are part of an R package that 
# is available for installation from this GitHub repository https://github.com/RetoSchmucki/regionalGAM
#
# For installation and usage, please refer to the second example in Appendix S3 or in follow
# the instructions in the README.md document available on the GitHub repository:
# https://github.com/RetoSchmucki/regionalGAM/tree/master
#
# ==================================================================================
#
# Version: 1.0
# R version 3.2.1
#
# Requiere packages
# 	- mgcv
#	- sampling
#
# Function available
#	- year_day_func()
#	- trap_area()
#	- trap_index()
#	- flight_count_sim()
#	- data_gen()
#	- flight_curve()
#	- flight_curve1y()
#	- degradation_prop()
#	- data_degradation()
#
#
# This script was developed by Reto Schmucki - reto.schmucki[at]mail.mcgill.ca
# Functions trap_area(), trap_index(), and flight_curve() where adapted from the script initially developed by
# Colin A. Harrower at the NERC Centre for Ecology and Hydrology (CEH), Crowmarsh Gifford, Oxfordshire, UK
#
# The two-stage model was first presented in: 
# Dennis, E.B., Freeman, S.N., Brereton, T. & Roy, D.B. (2013) Indexing butterfly abundance 
# whilst accounting for missing counts and variability in seasonal pattern. 
# Methods in Ecology and Evolution, 4, 637–645.
#
# A worked example for data simulation and degradation is provided in Appendix S3 available in 
# online supporting information.
# 
######################################################################################


#################################################
# 1. compute year, month and day full dataset
#################################################

year_day_func = function(sp_data) {
    
    year <- unique(sp_data$YEAR)
    origin.d <- paste(year, "01-01", sep = "-")
    
    # adjust for leep year with real data (i.e. years > 1000)
    if (year > 1000) {
    
    	if ((year%%4 == 0) & ((year%%100 != 0) | (year%%400 == 0))) {
        	nday <- 366
    	} else {
       		nday <- 365
    	}
    
    } else {
    nday <- 365
    }
    
    date.serie <- as.POSIXlt(seq(as.Date(origin.d), length = nday, by = "day"), format = "%Y-%m-%d")
    
    dayno <- as.numeric(julian(date.serie, origin = as.Date(origin.d)) + 1)
    month <- as.numeric(strftime(date.serie, format = "%m"))
    week <- as.numeric(strftime(date.serie, format = "%W"))
    week_day <- as.numeric(strftime(date.serie, format = "%u"))
    day <- as.numeric(strftime(date.serie, format = "%d"))
    
    site_list <- sp_data[!duplicated(sp_data$SITE), c("SITE")]
    
    all_day_site <- data.frame(SPECIES = sp_data$SPECIES[1], SITE = rep(site_list, 
        rep(nday, length(site_list))), YEAR = sp_data$YEAR[1], MONTH = month, WEEK = week, 
        DAY = day, DAY_WEEK = week_day, DAYNO = dayno, COUNT = NA)
    
    count_index <- match(paste(sp_data$SITE, sp_data$DAYNO, sep = "_"), paste(all_day_site$SITE, 
        all_day_site$DAYNO, sep = "_"))
    all_day_site$COUNT[count_index] <- sp_data$COUNT

    # add count number to ease extraction of single count
    site_count_length <- aggregate(sp_data$COUNT, by = list(sp_data$SITE), function(x) list(1:length(x)))
    names(site_count_length$x) <- as.character(site_count_length$Group.1)
    site_countno <- stack(site_count_length$x)
    all_day_site$COUNTNO <- NA
    all_day_site$COUNTNO[count_index] <- site_countno$values
    
    # Add zeros (ANCHORS) to close the season before and after the first and the last observation
    first_obs <- min(all_day_site$DAYNO[!is.na(all_day_site$COUNT)])
    last_obs <- max(all_day_site$DAYNO[!is.na(all_day_site$COUNT)])
    
    closing_season <- c((first_obs - 11):(first_obs - 7), (last_obs + 7):(last_obs + 
        11))
    
    # If closing season is before day 1 or day 365(366), use the five first and last days as ANCHORS
    if (min(closing_season) < 1) 
        closing_season[1:5] <- c(1:5)
    if (max(closing_season) > nday) 
        closing_season[6:10] <- c((nday - 4):nday)
    
    all_day_site$COUNT[all_day_site$DAYNO %in% closing_season] <- 0
    all_day_site$ANCHOR <- 0
    all_day_site$ANCHOR[all_day_site$DAYNO %in% closing_season] <- 1
    
    all_day_site <- all_day_site[order(all_day_site$SITE, all_day_site$DAYNO), ]
    
    return(all_day_site)
}

#################################################
# 2. compute area under the curve trapezoide method
#################################################

trap_area = function(x, y = NULL) {
    # If y is null and x has multiple columns then set y to x[,2] and x to x[,1]
    if (is.null(y)) {
        if (length(dim(x)) == 2) {
            y = x[, 2]
            x = x[, 1]
        } else {
            stop("ERROR: need to either specifiy both x and y or supply a two column data.frame/matrix to x")
        }
    }


    
    # Check x and y are same length
    if (length(x) != length(y)) {
        stop("ERROR: x and y need to be the same length")
    }
    
    # Need to exclude any pairs that are NA for either x or y
    rm_inds = which(is.na(x) | is.na(y))
    if (length(rm_inds) > 0) {
        x = x[-rm_inds]
        y = y[-rm_inds]
    }
    
    # Determine values of trapezoids under curve Get inds
    inds = 1:(length(x) - 1)
    # Determine area using trapezoidal rule Area = ( (b1 + b2)/2 ) * h where b1 and
    # b2 are lengths of bases (the parallel sides) and h is the height (the
    # perpendicular distance between two bases)
    areas = ((y[inds] + y[inds + 1])/2) * diff(x)
    
    # total area is sum of all trapezoid areas
    tot_area = sum(areas)
    
    # Return total area
    return(tot_area)
}

######################################
# 3. compute Abundance index
######################################

trap_index = function(sp_data, data_col = "IMP", time_col = "DAYNO", by_col = c("SPECIES", 
    "SITE", "YEAR")) {
    
    # Build output data.frame
    out_obj = unique(sp_data[, by_col])
    # Set row.names to be equal to collapsing of output rows (will be unique, you
    # need them to make uploading values back to data.frame will be easier)
    row.names(out_obj) = apply(out_obj, 1, paste, collapse = "_")
    
    # Using this row.names from out_obj above as index in by function to loop through
    # values all unique combs of by_cols and fit trap_area to data
    ind_dat = by(sp_data[, c(time_col, data_col)], apply(sp_data[, by_col], 1, paste, 
        collapse = "_"), trap_area)
    
    # Add this data to output object
    out_obj[names(ind_dat), "SINDEX"] = round(ind_dat/7, 1)
    
    # Set row.names to defaults
    row.names(out_obj) = NULL
    
    # Return output object
    return(out_obj)
}

######################################
# 4. simulate flight curve and counts
######################################

flight_count_sim <- function(fully_random_subsample = FALSE, nyear = 10, nsite = 100, 
    long_trend = -0.5, sd_peak_shift = 2, bivoltine = FALSE, keep_peak_in_subsample = TRUE, 
    iter = 1) {
    
    if (bivoltine == FALSE) {
        speciesGAM <- "Univoltine species"
    } else {
        speciesGAM <- "Bivoltine species"
    }
    
    pow_result_comb <- data.frame()
    
    annual_decline <- -1 * log(long_trend + 1)/nyear
    declinetrend <- round(annual_decline, 3)
    
    # generate a year day data.frame
    
    origin.d <- paste(1, "01-01", sep = "-")
    nday <- 365
    date.serie <- as.POSIXlt(seq(as.Date(origin.d), length = nday, by = "day"), format = "%Y-%m-%d")
    
	dayno <- as.numeric(julian(date.serie, origin = as.Date(origin.d)) + 1)
    month <- as.numeric(strftime(date.serie, format = "%m"))
    week <- as.numeric(strftime(date.serie, format = "%W"))
    week_day <- as.numeric(strftime(date.serie, format = "%u"))
    day <- as.numeric(strftime(date.serie, format = "%d"))
    
    year_day <- data.frame(JULIAN_D = dayno, MONTH = month, WEEK = week, WEEK_DAY = week_day, DAY=day)
    year_day_sample_cummul <- data.frame()
    
    t <- year_day$JULIAN_D
    
    my.list_yz <- list()
    
    for (z in (((iter - 1) * nsite) + 1):(nsite * iter)) {
        
        # generate lambda for population abundance with limits between 0.1 and 20
        lamb <- exp(abs(rnorm(1, 5, 1.5)))/length(2:51)
        if (lamb >= 1) {
            lamb <- lamb
        } else {
            lamb <- 1
        }
        if (lamb >= 20) {
            lamb <- 20
        } else {
            lamb <- lamb
        }
        
        for (y in 1:nyear) {
            
            #################################################### 
            
            # generate temporal trend in the abundance indices at site 'z'
            lamb <- lamb - ((lamb * declinetrend) * (min(1, (y - 1))))
            
            if (bivoltine == TRUE) {
                
                speciesGAM <- "Polyommatus bellargus"
                
                # EMERGENCE 1
                
                peak <- 150
                peak_shift <- rnorm(1, 0, sd_peak_shift)
                u1 <- peak + peak_shift
                sig <- 0.15
                beta <- 3
                
                fE1 <- (sig * (exp((t - u1)/beta)))/(beta * (1 + exp((t - u1)/beta)))^(sig + 
                  1)  # Logistic with skew
                sdfe <- fE1/sum(fE1)  # standardized to 1
                prop <- sdfe/mean(sdfe)
                
                # Generate count across the sampling season using a poisson process
                N <- c()
                for (i in t) {
                  N <- c(N, rpois(n = 1, lambda = (lamb * 0.6) * prop[i]))
                }
                
                #################### EMERGENCE 2
                
                peak2 <- u1 + round(rnorm(1, 60, 5))
                peak_shift2 <- rnorm(1, 0, 1)
                u <- peak2 + peak_shift2
                sig2 <- 0.15
                beta2 <- 3
                
                fE12 <- (sig2 * (exp((t - u)/beta2)))/(beta2 * (1 + exp((t - u)/beta2)))^(sig2 + 
                  1)  # Logistic with skew
                sdfe2 <- fE12/sum(fE12)  # standardized to 1
                prop2 <- sdfe2/mean(sdfe2)
                
                
                
                # Generate count across the sampling season using a poisson process for gen 2
                N2 <- c()
                for (i in t) {
                  N2 <- c(N2, rpois(n = 1, lambda = lamb * prop2[i]))
                }
                
                N <- N + N2
                
            } else {
                
                # EMERGENCE 1
                
                peak <- 190
                peak_shift <- rnorm(1, 0, sd_peak_shift)
                u <- peak + peak_shift
                sig <- 0.15
                beta <- 3
                
                fE1 <- (sig * (exp((t - u)/beta)))/(beta * (1 + exp((t - u)/beta)))^(sig + 
                  1)  # Logistic with skew
                sdfe <- fE1/sum(fE1)  # standardized to 1
                prop <- sdfe/mean(sdfe)
                
                ##################################################### 
                
                # Generate count across the sampling season using a poisson process
                N <- c()
                for (i in t) {
                  N <- c(N, rpois(n = 1, lambda = lamb * prop[i]))
                }
            }
            
            year_day$COUNT <- N
            

            ######################################################## 
            
            # Simulate sampling day pattern within weeks
            sampling_day <- data.frame(WEEK = 2:51, DAY = sample(2:6, length(2:51), 
                replace = T))
            year_day_sample <- year_day[paste(year_day$WEEK, year_day$WEEK_DAY) %in% 
                paste(sampling_day$WEEK, sampling_day$DAY), ]
            year_day_sample$SPECIES <- speciesGAM
            year_day_sample$SITE <- z
            year_day_sample$YEAR <- y
            year_day_sample$peak <- u
            year_day_sample$sig <- sig
            year_day_sample$beta <- beta
            year_day_sample$lamb <- lamb
            
            list.name_yz <- as.character(paste(y, z, sep = "_"))
            my.list_yz[[list.name_yz]] <- year_day_sample
            
            # year_day_sample_cummul <- rbind(year_day_sample_cummul,year_day_sample)
            
        }  # end of year loop
    }  # end of site loop
    
    year_day_sample_cummul <- do.call(rbind, my.list_yz)
    
    # Keep 26 weeks from April to September
    year_day_sample_cummul_26w <- year_day_sample_cummul[year_day_sample_cummul$WEEK >= 
        14 & year_day_sample_cummul$WEEK <= 39, ]
    
    return(year_day_sample_cummul_26w)
}



##################################################################
# 5. function to generate fake count datasets for multiple iteration
##################################################################

data_gen <- function(miter = 2, nsite = 1, long_trend = -0.05, nyear = 10, sd_peak_shift = 3, 
    bivoltine = FALSE) {
    my.list <- list()
    for (iter in 1:miter) {
        new.data <- flight_count_sim(nsite = nsite, iter = iter, long_trend = long_trend, 
            nyear = nyear, sd_peak_shift = sd_peak_shift, bivoltine = bivoltine)
        list.name <- as.character(iter)
        my.list[[list.name]] <- new.data
    }
    new.data <- do.call(rbind, my.list)
    return(new.data)
}

######################################################################
# 6. function to generate flight curve with 200 sites for multiple years
######################################################################

flight_curve <- function(year_day_sample_cummul_26w, y = 1) {
    
    require(mgcv) 

    dataset <- year_day_sample_cummul_26w[, c("SPECIES", "SITE", "YEAR", "MONTH", 
        "DAY","JULIAN_D", "COUNT")]
    dataset$X_COORD <- 1
    dataset$Y_COORD <- 1
    dataset <- data.frame(dataset, stringsAsFactors = FALSE)
    
    nsite <- length(unique(dataset$SITE))
    
    # Rename columns to standard names
    names(dataset) = c("SPECIES", "SITE", "YEAR", "MONTH", "DAY", "DAYNO", "COUNT", 
        "X_COORD", "Y_COORD")
    
    for (y in unique(dataset$YEAR)){
    
    dataset_y <- dataset[dataset$YEAR == y, ]
    if(nsite > 200){
	dataset_y <- dataset_y[as.character(dataset_y$SITE)%in%as.character(unique(dataset_y$SITE)[sample(1:nsite,200,replace=F)]),]
    } else { dataset_y <- dataset_y 
    }


    # Determine missing days and add to dataset
    sp_data_all <- year_day_func(dataset_y)
    sp_data_all <- sp_data_all
    
    # Expand ANCHORS before and after the season
    sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
    sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0

    sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
    gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), data = sp_data_all, 
        family = poisson(link = "log")), silent = TRUE)
    
    
    # second try if the GAM did not converge
    if (class(gam_obj_site)[1] == "try-error") {
        
        dataset_y <- dataset[dataset$YEAR == y, ]
        if(nsite > 200){
	    dataset_y <- dataset_y[as.character(dataset_y$SITE)%in%as.character(unique(dataset_y$SITE)[sample(1:nsite,200,replace=F)]),]
        } else { dataset_y <- dataset_y 
        }

        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)
        sp_data_all <- sp_data_all
        sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
        sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), 
            data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)

        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"

        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)

        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        
    } else {
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)
        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)
        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"
        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)
        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
    }
    
    # retry if the GAM produce peak in the first of last 5 weeks
    if (min(sp_data_all$WEEK[sp_data_all$NM == max(sp_data_all$NM)]) %in% c(1:5) | 
        max(sp_data_all$WEEK[sp_data_all$NM == max(sp_data_all$NM)]) %in% c(47:52)) {
        
        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)
        sp_data_all <- sp_data_all
        sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
        sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), 
            data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)

        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"

        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)

        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        
        sp_data_filled <- sp_data_all
        
        flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR, 
            week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO, 
            nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR, sp_data_filled$DAYNO, 
            sep = "_")), ]
        
        flight_curve <- flight_curve[order(flight_curve$DAYNO), ]
        
    } else {
        
        sp_data_filled <- sp_data_all
        
        flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR, 
            week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO, 
            nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR, sp_data_filled$DAYNO, 
            sep = "_")), ]
        
        flight_curve <- flight_curve[order(flight_curve$DAYNO), ]
    }
    
    
    # bind if exist else create
    if ("flight_pheno" %in% ls()) {
        flight_pheno <- rbind(flight_pheno, flight_curve)
    } else {
        flight_pheno <- flight_curve
    }
    
    } # end of year loop
    
    return(flight_pheno)
}

################################################################
# 7. function to generate flight curve with 200 sites for one year
################################################################

flight_curve1y <- function(dataset_y) {
    
    require(mgcv) 

    nsite <- length(unique(dataset_y$SITE))
    
    if(nsite > 200){
	dataset_y <- dataset_y[as.character(dataset_y$SITE)%in%as.character(unique(dataset_y$SITE)[sample(1:nsite,200,replace=F)]),]
    } else { dataset_y <- dataset_y 
    }


    # Determine missing days and add to dataset
    sp_data_all <- dataset_y
    sp_data_all <- sp_data_all

    # Expand ANCHORS before and after the season
    sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
    sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0

    sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
    gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), data = sp_data_all, 
        family = poisson(link = "log")), silent = TRUE)
    
    
    # second try if the GAM did not converge
    if (class(gam_obj_site)[1] == "try-error") {
        
        dataset_y <- dataset[dataset$YEAR == y, ]
        if(nsite > 200){
	    dataset_y <- dataset_y[as.character(dataset_y$SITE)%in%as.character(unique(dataset_y$SITE)[sample(1:nsite,200,replace=F)]),]
        } else { dataset_y <- dataset_y 
        }

        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)
        sp_data_all <- sp_data_all
        sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
        sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), 
            data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)

        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"

        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)

        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        
    } else {
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)
        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)
        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"
        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)
        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
    }
    
    # retry if the GAM produce peak in the first of last 5 weeks
    if (min(sp_data_all$WEEK[sp_data_all$NM == max(sp_data_all$NM)]) %in% c(1:5) | 
        max(sp_data_all$WEEK[sp_data_all$NM == max(sp_data_all$NM)]) %in% c(47:52)) {
        
        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)
        sp_data_all <- sp_data_all
        sp_data_all[sp_data_all$trimDAYNO < 70, "COUNT"] <- 0
        sp_data_all[sp_data_all$trimDAYNO > 305, "COUNT"] <- 0
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        gam_obj_site <- try(gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE), 
            data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        
        # Generate a list of values for all days from the aditive model and use these
        # value to fill the missing obserations
        sp_data_all[, "FITTED"] <- predict.gam(gam_obj_site, newdata = sp_data_all[, 
            c("trimDAYNO", "SITE")], type = "response")
        sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
        sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
        
        # Define the flight curve from the fitted values and append them over years (this
        # is one flight curve per year for all site)        
        site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE), 
            FUN = sum)

        # Rename sum column
        names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"

        # Add data to sp_data data.frame (ensure merge does not sort the data!)
        sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE, sort = FALSE)

        # Calculate normalised values
        sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        
        sp_data_filled <- sp_data_all
        
        flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR, 
            week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO, 
            nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR, sp_data_filled$DAYNO, 
            sep = "_")), ]
        
        flight_curve <- flight_curve[order(flight_curve$DAYNO), ]
        
    } else {
        
        sp_data_filled <- sp_data_all
        
        flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR, 
            week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO, 
            nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR, sp_data_filled$DAYNO, 
            sep = "_")), ]
        
        flight_pheno <- flight_curve[order(flight_curve$DAYNO), ]
    }
    
    
    return(flight_pheno)
}


################################################
# 8. compute degradation pattern for count dataset
################################################

degradation_prop <- function(year_day_sample_cummul_26w, pheno = flight_pheno, perct = 0.6, 
    fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE) {
    
    require(sampling)
    
    if ("degradation_pattern_cummul" %in% ls()) 
        rm(degradation_pattern_cummul)
    
    # start year loop
    for (y in unique(year_day_sample_cummul_26w$YEAR)) {
        
        year_pheno <- pheno[pheno$year == y, ]
        regional_peak_week <- year_pheno$week[year_pheno$nm == max(year_pheno$nm)]
        
        sp_data_year <- year_day_sample_cummul_26w[year_day_sample_cummul_26w$YEAR == 
            y, ]
        
        for (iteration in unique(as.character(sp_data_year$SITE))) {
            
            sp_data_site <- sp_data_year[as.character(sp_data_year$SITE) == iteration, 
                ]
            
            peak_index1 <- which(abs(sp_data_site$WEEK - regional_peak_week) == min(abs(sp_data_site$WEEK - 
                regional_peak_week)))[1]
            
            if (fully_random_subsample == TRUE) {
                
                tot_length <- length(sp_data_site$WEEK)
                week_vect <- sp_data_site$WEEK
                sample_size <- round(tot_length * perct)
                strat_size <- (tot_length)%/%sample_size
                strat_1 <- rep(c(1:sample_size), strat_size)
                strat <- c(strat_1, sample(unique(strat_1), (tot_length) - length(strat_1), 
                  replace = F))
                strat <- strat[order(strat)]
                data_strat <- data.frame(ind = c(week_vect), strat = strat, prob = rep(c(0.05, 
                  0.25, 0.7, 0.05, 0.25, 0.7)[sample(1:3, 1):5][1:3], 15)[1:(tot_length)])
                subset_weekno <- data_strat$ind[strata(data_strat, stratanames = "strat", 
                  size = rep(1, tot_length), method = c("systematic"), pik = data_strat$prob)$ID_unit]
                subset_weekno <- unique(subset_weekno)
                subset_weekno <- subset_weekno[order(subset_weekno)]
                
            } else {
                
                
                if (keep_peak_in_subsample == TRUE) {
                  
                  # force the peak
                  if (perct >= 0.5) {
                    tot_length <- length(sp_data_site$WEEK)
                    week_vect <- sp_data_site$WEEK
                    sample_size <- round(tot_length * perct) - 1
                    strat_size <- (tot_length - 1)%/%sample_size
                    strat_1 <- rep(c(1:sample_size), strat_size)
                    strat <- c(strat_1, sample(unique(strat_1), (tot_length) - length(strat_1) - 
                      1, replace = F))
                    strat <- strat[order(strat)]
                    data_strat <- data.frame(ind = c(week_vect)[-peak_index1], strat = strat, 
                      prob = rep(c(0.05, 0.25, 0.7, 0.05, 0.25, 0.7)[sample(1:3, 
                        1):5][1:3], 15)[1:(tot_length)][-peak_index1])
                    subset_weekno <- c(data_strat$ind[strata(data_strat, stratanames = "strat", 
                      size = rep(1, tot_length), method = c("systematic"), pik = data_strat$prob)$ID_unit], 
                      week_vect[peak_index1])
                    subset_weekno <- unique(subset_weekno)
                    subset_weekno <- subset_weekno[order(subset_weekno)]
                  } else {
                    tot_length <- length(sp_data_site$WEEK)
                    week_vect <- sp_data_site$WEEK
                    sample_size <- round(tot_length * perct) - 1
                    strat_size <- (tot_length - 3)%/%sample_size
                    strat_1 <- rep(c(1:sample_size), strat_size)
                    strat <- c(strat_1, sample(unique(strat_1), (tot_length) - length(strat_1) - 
                      3, replace = F))
                    strat <- strat[order(strat)]
                    data_strat <- data.frame(ind = c(week_vect)[-c((peak_index1 - 
                      1):(peak_index1 + 1))], strat = strat, prob = rep(c(0.05, 0.25, 
                      0.7, 0.05, 0.25, 0.7)[sample(1:3, 1):5][1:3], 15)[1:(tot_length)][-c((peak_index1 - 
                      1):(peak_index1[1] + 1))])
                    subset_weekno <- data_strat$ind[strata(data_strat, stratanames = "strat", 
                      size = rep(1, tot_length), method = c("systematic"), pik = data_strat$prob)$ID_unit]
                    subset_weekno <- c(week_vect[peak_index1], subset_weekno)
                    subset_weekno <- unique(subset_weekno)
                    subset_weekno <- subset_weekno[order(subset_weekno)]
                  }
                  # force one around the peak -> this should be around the peak value and not
                  # include the peak
                } else {
                  if (perct >= 0.5) {
                    tot_length <- length(sp_data_site$WEEK)
                    week_vect <- sp_data_site$WEEK
                    sample_size <- round(tot_length * perct)
                    strat_size <- (tot_length - 1)%/%sample_size
                    strat_1 <- rep(c(1:sample_size), strat_size)
                    strat <- c(strat_1, sample(unique(strat_1), (tot_length) - length(strat_1) - 
                      1, replace = F))
                    strat <- strat[order(strat)]
                    data_strat <- data.frame(ind = week_vect[-peak_index1], strat = strat, 
                      prob = rep(c(0.05, 0.25, 0.7, 0.05, 0.25, 0.7)[sample(1:3, 
                        1):5][1:3], 15)[1:(tot_length)][-peak_index1])
                    subset_weekno <- data_strat$ind[strata(data_strat, stratanames = "strat", 
                      size = rep(1, tot_length), method = c("systematic"), pik = data_strat$prob)$ID_unit]
                    subset_weekno <- unique(subset_weekno)
                    subset_weekno <- subset_weekno[order(subset_weekno)]
                  } else {
                    tot_length <- length(sp_data_site$WEEK)
                    week_vect <- sp_data_site$WEEK
                    sample_size <- round(tot_length * perct) - 1
                    strat_size <- (tot_length - 3)%/%sample_size
                    strat_1 <- rep(c(1:sample_size), strat_size)
                    strat <- c(strat_1, sample(unique(strat_1), (tot_length) - length(strat_1) - 
                      3, replace = F))
                    strat <- strat[order(strat)]
                    data_strat <- data.frame(ind = c(week_vect)[-c((peak_index1 - 
                      1):(peak_index1 + 1))], strat = strat, prob = rep(c(0.05, 0.25, 
                      0.7, 0.05, 0.25, 0.7)[sample(1:3, 1):5][1:3], 15)[1:(tot_length)][-c((peak_index1 - 
                      1):(peak_index1 + 1))])
                    subset_weekno <- data_strat$ind[strata(data_strat, stratanames = "strat", 
                      size = rep(1, tot_length), method = c("systematic"), pik = data_strat$prob)$ID_unit]
                    subset_weekno <- c(sample(week_vect[c(peak_index1 - 1, peak_index1 + 
                      1)], 1), subset_weekno)
                    subset_weekno <- unique(subset_weekno)
                    subset_weekno <- subset_weekno[order(subset_weekno)]
                  }
                }
            }
            
            weeksnumbermat <- as.data.frame(matrix(NA, ncol = 26, nrow = 1))
            weeksnumbermat[1:length(subset_weekno)] <- subset_weekno
            
            degradation_pattern <- cbind(data.frame(year = y, peak_week = regional_peak_week, 
                prop_week_sampled = perct, iteration = iteration), weeksnumbermat)
            
            # bind if exist else create
            if ("degradation_pattern_cummul" %in% ls()) {
                degradation_pattern_cummul <- rbind(degradation_pattern_cummul, degradation_pattern)
            } else {
                degradation_pattern_cummul <- degradation_pattern
            }
            
        }  # close iteration loop (nsite)
    }  # close year loop
    
    return(list(degradation_pattern_cummul = degradation_pattern_cummul, pheno = pheno))
}

#####################################################
# 9. compute Abundance indices for degraded datasets
#####################################################

data_degradation <- function(year_day_sample_cummul_26w =new.data, degra = degra, perct = 1,
    fully_random_subsample = FALSE, keep_peak_in_subsample = TRUE) {

     if ("cumullated_indices" %in% ls()) 
        rm(cumullated_indices)
    cumullated_indices <- data.frame()
    
    for (y in unique(year_day_sample_cummul_26w$YEAR)) {
        
        # year_pheno <- degra$pheno[degra$pheno$year == y, ]
        sp_data_site <- year_day_sample_cummul_26w[year_day_sample_cummul_26w$YEAR == 
            y, ]
        
        perct_deg <- degra$degradation_pattern_cummul[degra$degradation_pattern_cummul$year == 
            y, ]
        
        site_deg_mat <- as.matrix(perct_deg[, 5:30])
        
        site_deg_long <- data.frame(site = rep(unique(sp_data_site$SITE), rep(ncol(site_deg_mat), 
            nrow(site_deg_mat))), week = c(t(site_deg_mat)))
        
        site_deg_clean <- site_deg_long[!is.na(site_deg_long$week), ]
        
    	dataset <- sp_data_site[, c("SPECIES", "SITE", "YEAR", "MONTH", 
            "DAY", "JULIAN_D", "COUNT")]
        dataset$X_COORD <- 1
        dataset$Y_COORD <- 1
        dataset <- data.frame(dataset, stringsAsFactors = FALSE)
    
        nsite <- length(unique(dataset$SITE))
    
        # Rename columns to standard names
        names(dataset) = c("SPECIES", "SITE", "YEAR", "MONTH", "DAY", "DAYNO", "COUNT", 
            "X_COORD", "Y_COORD")

        sp_data_site <- year_day_func(dataset)
        
        sp_data_site$trimDAYNO <- sp_data_site$DAYNO - min(sp_data_site$DAYNO) + 1
        
        sp_data_site$COUNT[!paste(sp_data_site$SITE, sp_data_site$WEEK, sep = "_") %in% 
            paste(site_deg_clean$site, site_deg_clean$week, sep = "_")] <- NA
        sp_data_site$COUNT[sp_data_site$ANCHOR == 1] <- 0

        year_pheno <- flight_curve1y(sp_data_site)
     
        sp_data_site <- merge(sp_data_site, year_pheno[, c("DAYNO", "nm")], by = c("DAYNO"), 
            all.x = TRUE, sort = FALSE)
        
        # compute proportion of the flight curve sampled due to missing visits
        pro_missing_count <- data.frame(SITE = sp_data_site$SITE, WEEK = sp_data_site$WEEK, 
            NM = sp_data_site$nm, COUNT = sp_data_site$COUNT, ANCHOR = sp_data_site$ANCHOR)
        pro_missing_count$site_week <- paste(pro_missing_count$SITE, pro_missing_count$WEEK, 
            sep = "_")
        siteweeknocount <- aggregate(pro_missing_count$COUNT, by = list(pro_missing_count$site_week), 
            function(x) sum(!is.na(x)) == 0)
        pro_missing_count <- pro_missing_count[pro_missing_count$site_week %in% siteweeknocount$Group.1[siteweeknocount$x == 
            TRUE], ]
        pro_count_agg <- aggregate(pro_missing_count$NM, by = list(pro_missing_count$SITE), 
            function(x) 1 - sum(x, na.rm = T))
        names(pro_count_agg) <- c("SITE", "PROP_PHENO_SAMPLED")
        
        # Compute the regional gam index
        glm_obj_site <- glm(COUNT ~ factor(SITE) + offset(log(nm)) - 1, data = sp_data_site, 
            family = quasipoisson(link = "log"), control = list(maxit = 100))
        
        sp_data_site[, "FITTED"] <- predict.glm(glm_obj_site, newdata = sp_data_site, 
            type = "response")
        sp_data_site[, "COUNT_IMPUTED"] <- sp_data_site$COUNT
        sp_data_site[is.na(sp_data_site$COUNT), "COUNT_IMPUTED"] <- sp_data_site$FITTED[is.na(sp_data_site$COUNT)]
        
        # add fitted value for missing mid-week data
        sp_data_site <- sp_data_site[!paste(sp_data_site$DAY_WEEK, sp_data_site$COUNT) %in% 
            c("1 NA", "2 NA", "3 NA", "5 NA", "6 NA", "7 NA"), ]
        
        # remove all added mid-week values for weeks with real count (observation)
        sp_data_site$site_week <- paste(sp_data_site$SITE, sp_data_site$WEEK, sep = "_")
        siteweekcount <- aggregate(sp_data_site$COUNT, by = list(sp_data_site$site_week), 
            function(x) sum(!is.na(x)) > 0)
        sp_data_site <- sp_data_site[!(is.na(sp_data_site$COUNT) & (sp_data_site$site_week %in% 
            siteweekcount$Group.1[siteweekcount$x == TRUE])), names(sp_data_site) != 
            "site_week"]
        
        # Compute the regional gam index
        regional_gam_index <- trap_index(sp_data_site, data_col = "COUNT_IMPUTED", 
            time_col = "DAYNO", by_col = c("SPECIES", "SITE", "YEAR"))
        year_count_index <- trap_index(sp_data_site, data_col = "COUNT", time_col = "DAYNO", 
            by_col = c("SPECIES", "SITE", "YEAR"))
        
        cumu_index <- merge(regional_gam_index, pro_count_agg, by = c("SITE"), all.x = TRUE, 
            sort = FALSE)
        cumu_index <- merge(cumu_index, year_count_index[, c("SITE", "SINDEX")], 
            by = c("SITE"), all.x = TRUE, sort = FALSE)
        names(cumu_index) <- c("SITE", "SPECIES", "YEAR", "regional_gam", "prop_pheno_sampled", 
            "linear_interp_index")
        cumu_index$regional_peak_week <- degra$degradation_pattern_cummul[degra$degradation_pattern_cummul$year == 
            y, ][1, "peak_week"]
        cumu_index <- cumu_index[,c("SITE", "SPECIES", "YEAR", "regional_gam","linear_interp_index","regional_peak_week","prop_pheno_sampled")] 
        cumu_index$prop_week_sampled <- perct
        cumu_index <- cumu_index[order(cumu_index$SITE), ]
        
        
        # bind if exist else create
        if ("cumullated_indices" %in% ls()) {
            cumullated_indices <- rbind(cumullated_indices, cumu_index)
        } else {
            cumullated_indices <- cumu_index
        }
        
    }  # end of year loop
    
    return(cumullated_indices)
}



