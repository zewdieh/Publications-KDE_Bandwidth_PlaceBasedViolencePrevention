#-------------------------------------------------------------------------------  
# Title: Code for manuscript: Bandwidth selection in neighborhood violence 
#        prevention research 
# R Script Author: Hiwot Y. Zewdie
# Affiliation: University of Washington
# Date: Sept 2025
#
# Description: Code to generate results for two manuscript objectives: 
#               1) Quantify variability in data-driven bandwidth selectors
#               2) Simulate the changes in intervention effects across 
#                  bandwidth sizes
#-------------------------------------------------------------------------------

#load libraries
library(spatstat)
library(raster)
library(terra)
library(sf)
library(tidyverse)
library(tidycensus)
library(mapview)
library(ggpubr)
library(here)

#-------------------------------------------------------------------------------  
#                            Objective 0: Set up 
#-------------------------------------------------------------------------------

  #set working directory
  setwd(here())

  #load Philadelphia county geographic extent  
        philly <- get_acs(
          state = "PA",
          county = "Philadelphia",
          geography = "county",
          variables = "B19013_001",
          geometry = T,
          year = 2015)%>%
          dplyr::select(-c(variable,estimate,moe))%>%
          st_transform(phily,crs = 2272) #philly coordinate system, units = survey ft
        
        philly_owin <- as.owin(philly) #observation window (owin) object for KDE
        
        #load crime 
        vcrime <- readRDS("data/kde_manuscript_crime_data.RDS")

#-------------------------------------------------------------------------------  
#     Objective 1: Quantify variability in data-driven bandwidth selectors
#-------------------------------------------------------------------------------
        #set up parameters for loop        
        years <- c(2013:2023) 
        crime_type <-c("all", "all_firearm","assault", "assault_firearm", 
                       "homicide", "homicide_firearm", "robbery",
                       "robbery_firearm")
        bw <- c("bw_default_cvl", "bw_default_dig") 
        
        
        #create empty dataframe to store results 
        kde_results <- data.frame()
        
        #for loop
        for (b in bw) {
          for (i in years) {
            for (c in crime_type) {
              
              #filter all violent crime data by year
              vcrime_filtered <- vcrime%>%
                filter(year == i)
              
              #filter data by crime types
              if (c != "all" & c != "all_firearm") { 
                
                crime_base <- gsub("_.*", "", c)
                firearm <- grepl("_f", c)
                
                vcrime_filtered <- vcrime_filtered%>%
                  filter(vcrime_cat == crime_base)
                
                if (firearm) {
                  vcrime_filtered <- vcrime_filtered%>%
                    filter(grepl("_f", crime_type))
                }
                
              }else if (c == "all_firearm") {
                vcrime_filtered <- vcrime_filtered%>%
                  filter(grepl("_f", crime_type))
              }
              
              #transform resulting crime data to spatial points
              vcrime_sf <- vcrime_filtered%>%
                st_as_sf(coords = c("Y_2272", "X_2272"), crs = 4269)%>% 
                st_transform(crs = 2272)%>% 
                dplyr::select(geometry)
              
              #create ppp object from crime data
              vcrime_ppp <- ppp(
                x = st_coordinates(vcrime_sf)[, 1],
                y = st_coordinates(vcrime_sf)[, 2],
                window = philly_owin)
              
              #check for rejected points - outside philly extent 
              x <- attr(vcrime_ppp, "rejects")
              if (!is.null(x) && x$n > 0) {
                print(paste("Rejected points for bandwidth:", b, 
                            "| Crime type:", c, "| Year:", i, ":", x$n, "points"))
              } 
              
              #set bandwidth (data driven) 
              if (b == "bw_default_cvl") {
                sig <- bw.CvL(vcrime_ppp)
              } else if (b == "bw_default_dig") {
                sig <- bw.diggle(vcrime_ppp)
              } 
              
              temp <- data.frame(
                year = i,
                crime_type = c,
                n_crime = nrow(vcrime_filtered),
                bandwidth_name = b,
                bandwidth_size = sig[1]*sqrt(7)) #to get sigma in terms of bw (in ft)
              
              #bind to og kde_results dataframe
              kde_results <- rbind(kde_results, temp)
            }
          }
        }

        #-----------------------------------------------------------------------
        #                   Objective 1: Visualizations 
        #-----------------------------------------------------------------------
        
        #Table 1
        kde_results%>%
          filter(bandwidth_name == "bw_default_cvl")%>% #swap out as needed 
          group_by(crime_type)%>%
          summarise(mean_bw = mean(bandwidth_size), 
                    sd_bw = sd(bandwidth_size),
                    min_bw = min(bandwidth_size),
                    max_bw = max(bandwidth_size),
                    mean_n = mean(n_crime))
        
        #Figure 2. Default bandwidth sizes across years and crime types 
        kde_results%>%
          filter(bandwidth_name == "bw_default_cvl")%>% #swap out as needed 
          mutate(crime_type_pretty = case_when(
            crime_type == "all" ~ "Overall Violent Crime", 
            crime_type == "all_firearm" ~ "Overall Firearm Violent Crime", 
            crime_type == "assault" ~ "Aggravated Assault", 
            crime_type == "assault_firearm" ~ "Firearm Aggravated Assault", 
            crime_type == "homicide" ~ "Homicide", 
            crime_type == "homicide_firearm" ~ "Firearm Homicide", 
            crime_type == "robbery" ~ "Robbery", 
            crime_type == "robbery_firearm" ~ "Firearm Robbery"), 
            crime_type_pretty = factor(crime_type_pretty, 
                                       levels = c("Overall Violent Crime", 
                                                  "Overall Firearm Violent Crime", 
                                                  "Aggravated Assault", 
                                                  "Firearm Aggravated Assault", 
                                                  "Homicide",
                                                  "Firearm Homicide",
                                                  "Robbery", 
                                                  "Firearm Robbery")))%>%
          ggplot(aes(x = year, 
                     y = bandwidth_size))+
          geom_line()+
          labs(x = "Year", y = "Bandwidth size (ft)", 
               title = "a) Cronie and van Lieshout's criterion")+ #swap out title as needed 
          scale_x_continuous(breaks = as.numeric(seq(2000, 2025, by = 2))) +
          scale_y_continuous(limits=c(0,60000))+ #change scale for diggle 
          facet_wrap(~crime_type_pretty, 4)+
          theme_minimal()

        
        
#-----------------------------------------------------------------------
#                   Objective 2: Simulations
#-----------------------------------------------------------------------
        
        #load vacant lot data 
        lots <- readRDS("data/kde_manuscript_lot_data.RDS")%>%
          st_as_sf(coords=c("long", "lat"), crs = 4269)%>%
          st_transform(lots, crs = 2272) #philly crs
        
        #filter crime to just assaults for simulations
        vcrime_sf <- vcrime%>%
          filter(vcrime_cat == "assault")%>%
          st_as_sf(coords = c("Y_2272", "X_2272"), crs = 4269)%>% 
          st_transform(crs = 2272)
        
        #set up simualtion parameters
        n_iterations <- 30  #number of randomization iterations
        set_of_seeds <- 1:n_iterations #random seeds
        years <- 2013:2023  #year range
        bw <- seq(100, 10000, by = 100) #bandwidth values(ft)

        #set up output directories to save sim results
        dir.create("analysis/sim_rds", showWarnings = F, recursive = T)
        
        for (i in set_of_seeds) {
          seed <- set_of_seeds[i]
          message("Starting iteration ", i, " (seed = ", seed, ")") #tracking
          set.seed(seed)
          
          #-------------------------------------------------------------
          #    Objective 2: create adjusted crime data 
          #-------------------------------------------------------------
          
          #randomize treated lots
          lots <- lots%>%
            mutate(pb_intervention = sample(c(rep(1, floor(nrow(.) * 0.5)), #random assignment of place-based intervention
                                              rep(0, ceiling(nrow(.) * 0.5)))))%>%
            mutate(cluster = as.character(cluster))
          
          #buffers + subset crimes
          treated_lots <- lots%>%
            filter(pb_intervention == 1)%>%
            st_buffer(dist = 500) #spatial extent of place-based intervention
          intersects_itsx <- st_intersects(vcrime_sf, treated_lots) 
          near_flag <- lengths(intersects_itsx) > 0 #flag assaults within 500 ft of treated lot
          
          vcrime_near_treated_post2016 <- vcrime_sf%>%
            mutate(near_treated = near_flag)%>%
            filter(year > 2016 & near_treated) #keep assaults after 2016 and within 500ft of treated lot
          
          set.seed(seed)
          
          #adjust crime data 
          keep_near <- sample(seq_len(nrow(vcrime_near_treated_post2016)),
                              size = floor(nrow(vcrime_near_treated_post2016) * 0.7)) #remove 30% of assaults 
          vcrime_near_treated_post2016_filtered <- vcrime_near_treated_post2016[keep_near, ]
         
          vcrime_pre2016 <- vcrime_sf%>%filter(year <= 2016)
          
          vcrime_outside_buffer_post2016 <- vcrime_sf%>%
            filter(year > 2016 & !near_flag)
          
          #bind together adjusted crime for simulation 
          vcrime_final_iter <- bind_rows(vcrime_near_treated_post2016_filtered,
                                         vcrime_pre2016,
                                         vcrime_outside_buffer_post2016)%>%
            mutate(year = as.integer(year))
          
          #---------------------------------------------------------------------
          #    Objective 2: simulate over bandwidths - kernel density estimation 
          #---------------------------------------------------------------------
          #set up storage for this iteration
          kde_results_iter <- data.frame()
          model_results_iter <- data.frame()
          
          #loop over bandwidths and years
          for (b in bw) {
            sig <- b / sqrt(7) #to set bw in terms of quartic kernel 
            
            for (yr in years) {
              vcrime_filtered <- vcrime_final_iter%>%
                filter(year == yr)
              
              print(paste("Iteration:", seed,", Bandwidth size:",b," Year:",yr)) #tracking 
              
              if (nrow(vcrime_filtered) == 0) {
                temp <- data.frame(cluster = lots$cluster,
                                   year = yr,
                                   treat = lots$pb_intervention,
                                   n_crime = 0,
                                   bw_value = b,
                                   kde_value = NA_real_,
                                   iter = i)
                kde_results_iter <- rbind(kde_results_iter, temp)
                next
              }
              
              coords <- st_coordinates(vcrime_filtered)
              vcrime_ppp <- ppp(x = coords[,1], y = coords[,2], window = philly_owin)
              
              #fit kde
              kde <- density.ppp(x = vcrime_ppp, 
                                 sigma = sig, 
                                 kernel = "quartic", 
                                 diggle = T, 
                                 eps = 100) 
              
              kde_raster <- raster(kde)
              crs(kde_raster) <- "EPSG:2272"
              kde_raster[] <- kde_raster[] * 27878400 #transform density values from sq ft to sq mile 
              
              kde_values <- raster::extract(kde_raster, lots)
              
              temp <- data.frame(cluster = lots$cluster,
                                 year = yr,
                                 treat = lots$pb_intervention,
                                 n_crime = nrow(vcrime_filtered),
                                 bw_value = b,
                                 kde_value = kde_values,
                                 iter = i)
              kde_results_iter <- rbind(kde_results_iter, temp)
            }
            
            #---------------------------------------------------------------------
            #    Objective 2: simulate over bandwidths - regression model
            #---------------------------------------------------------------------
            kde_bw <- kde_results_iter%>%
              filter(bw_value == b)%>%
              mutate(post = ifelse(year > 2016, 1, 0))%>%
              filter(!is.na(kde_value))
            
              did_model <- try(glm(kde_value ~ treat * post,
                                   data = kde_bw%>%
                                     mutate(kde_value = round(kde_value, 0)), 
                                   family = "poisson"), silent = TRUE)
              if (!inherits(did_model, "try-error")) {
                sm <- summary(did_model)
                  coef_val <- sm$coefficients["treat:post", "Estimate"]
                  se_val   <- sm$coefficients["treat:post", "Std. Error"]
                
                model_results_iter <- rbind(model_results_iter,
                                            data.frame(iter = i,
                                                       bw_value = b,
                                                       coef = coef_val,
                                                       se = se_val))
              }
            }
          
          #save outputs from this set of iterations
          saveRDS(kde_results_iter,
                  file = sprintf("/analysis/sim_rds/kde_results_iter%02d.rds", i))
          saveRDS(model_results_iter,
                  file = sprintf("/analysis/sim_rds/model_results_iter%02d.rds", i))
        }
        
        #-----------------------------------------------------------------------
        #    Objective 2: simulate over bandwidths - average results for viz
        #-----------------------------------------------------------------------
        
        model_files <- list.files("/analysis/sim_rds", pattern = "model_results_iter.*rds", full.names = T)
        model_results_all <- map_dfr(model_files, readRDS)
        
        #summarize across iterations
        bandwidth_effects <- model_results_all%>%
          group_by(bw_value)%>%
          summarise(
            mean_coef = mean(coef, na.rm = T),
            sd_coef = sd(coef, na.rm = T),
            n_iter = n(),
            .groups = "drop")%>%
          mutate(
            exp_mean = exp(mean_coef),
            exp_lower = exp(mean_coef - sd_coef),
            exp_upper = exp(mean_coef + sd_coef))
        
        #-----------------------------------------------------------------------
        #    Objective 2: simulate over bandwidths - plots 
        #-----------------------------------------------------------------------
        
        #prep iteration-level data for spaghetti lines
        bandwidth_effects_iter <- model_results_all%>%
          mutate(exp_est = exp(coef))
        
        #Fig 3
        ggplot() +
          # spaghetti lines (each iteration faint)
          geom_line(data = bandwidth_effects_iter,
                    aes(x = bw_value, y = exp_est, group = iter),
                    color = "grey60", alpha = 0.2)+
          geom_vline(aes(xintercept = 500), linetype = 3)+ #spatial extent of place-based intervention 
          geom_line(aes(x = bandwidth_effects$bw_value, y = 0.7), linetype = 3)+ #true effect 
          geom_text(aes(x=5000, label="'True effect'", y=0.71), 
                    colour="blue", angle=0, size=3) +
          geom_line(aes(x = bandwidth_effects$bw_value, y = 1), linetype = 1)+ #null 
          geom_line(data = bandwidth_effects, aes(x = bw_value, y = exp_mean), size = 1) + #average coef (line)
          geom_point(data = bandwidth_effects, aes(x = bw_value, y = exp_mean), size = 2) + #average coef (point)
          geom_ribbon(data = bandwidth_effects, aes(x = bw_value, ymin = exp_lower, #sd
                                                    ymax = exp_upper), alpha = 0.5)+ 
          labs(x = "Bandwidth (ft)", y = "Intervention effect") +
          theme_minimal()
        
        
        
       