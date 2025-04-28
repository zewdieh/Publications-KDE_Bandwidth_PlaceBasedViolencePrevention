#-------------------------------------------------------------------------------  
# Title: Code for manuscript: Bandwidth selection in neighborhood violence 
#        prevention research 
# R Script Author: Hiwot Y. Zewdie
# Affiliation: University of Washington
# Date: May 2025
#
# Description: Provides code for the two manuscript objectives: 
#               1) Quantify variability in data-driven bandwidth selectors
#               2) Simulate the changes in interventional effects across 
#                  bandwidth sizes
#-------------------------------------------------------------------------------

#load libraries
library(spatstat)
library(raster)
library(tidyverse)
library(tidycensus)
library(terra)
library(mapview)
library(sf)
library(ggpubr)

#-------------------------------------------------------------------------------  
#     Objective 1: Quantify variability in data-driven bandwidth selectors
#-------------------------------------------------------------------------------

        #load Philadelphia county geographic extent  
        philly <- get_acs(
          state = "PA",
          county = "Philadelphia",
          geography = "county",
          variables = "B19013_001",
          geometry = T,
          year = 2015)%>%
          select(-c(variable,estimate,moe))%>%
          st_transform(phily,crs = 2272)
        
        #load crime 
        vcrime <- readRDS("kde_manuscript_crime_data.RDS")
        
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
                st_transform(crs = 2272)%>% #philly coordinate system, units = survey ft
                dplyr::select(geometry)
              
              #set philly extent (Philly county shapefile from census)
              philly_owin <- as.owin(philly)
              
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
              
              #set bandwidth - 500, 5000, data driven 
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
          filter(bandwidth_name == "bw_default_dig")%>% #swap out as needed 
          group_by(crime_type)%>%
          summarise(mean_bw = mean(bandwidth_size), 
                    sd_bw = sd(bandwidth_size),
                    min_bw = min(bandwidth_size),
                    max_bw = max(bandwidth_size),
                    mean_n = mean(n_crime))
        
        #Figure 3. Default bandwidth sizes across years and crime types 
        kde_results%>%
          filter(bandwidth_name == "bw_default_dig")%>% #swap out as needed 
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
          scale_y_continuous(limits=c(0,60000))+
          facet_wrap(~crime_type_pretty, 4)+
          theme_minimal()

        
        
#-----------------------------------------------------------------------
#                   Objective 2: Simulations
#-----------------------------------------------------------------------
        
        #load vacant lot data 
        lots <- readRDS("kde_manuscript_lot_data.RDS")%>%
          st_as_sf(coords=c("long", "lat"), crs = 4269)%>%
          st_transform(lots, crs = 2272) #philly crs
        
        set.seed(4)
        lots <- lots%>% 
          #randomly assign lots to place based intervention (i.e., community grill)
          mutate(community_grill = sample(c(rep(1, floor(nrow(.) * 0.5)), 
                                            rep(0, ceiling(nrow(.) * 0.5)))))  
        
        #filter crime to just assaults for simulations
        vcrime_sf <- vcrime%>%
          filter(vcrime_cat == "assault")%>%
          st_as_sf(coords = c("Y_2272", "X_2272"), crs = 4269)%>% 
          st_transform(crs = 2272)
        
                  #-------------------------------------------------------------
                  #    Objective 2: create adjusted crime data for simulations
                  #-------------------------------------------------------------
                  
                  #1. define 'nearby'
                  treated_lots <- lots%>%
                    filter(community_grill == 1)%>%
                    st_buffer(dist = 500)
                  
                  #2. identify crime incidents occurring after 2016 that are within 500ft of treated lots
                  vcrime_near_treated_post2016 <- vcrime_sf%>%
                    filter(year > 2016 & lengths(st_intersects(geometry, treated_lots)) > 0)
                  
                  #3. randomly remove 30% of these crime incidents nearby - 'true effect'
                  set.seed(123)
                  vcrime_near_treated_post2016_filtered <- vcrime_near_treated_post2016%>%
                    sample_frac(0.7)#keep 70% (remove 30%)
                  
                  #4. keep all crimes from 2016 and earlier (unchanged)
                  vcrime_pre2016 <- vcrime_sf%>%
                    filter(year <= 2016)
                  
                  #5. keep crimes after 2016 that were not in the buffer (unchanged)
                  vcrime_outside_buffer_post2016 <- vcrime_sf%>%
                    filter(year > 2016 & lengths(st_intersects(geometry, treated_lots)) == 0)
                  
                  #6. combine all datasets
                  vcrime_final <- bind_rows(
                    vcrime_near_treated_post2016_filtered,  #post-2016 crimes near treated lots (30% removed)
                    vcrime_pre2016,  #all crimes from 2016 and earlier (unchanged)
                    vcrime_outside_buffer_post2016)  #post-2016 crimes outside buffer (unchanged)

                  #-------------------------------------------------------------
                  #    Objective 2: estimate crime density surfaces (annual)
                  #-------------------------------------------------------------
                  
                  years <- c(2013:2023)
                  bw <- seq(100,5000, by = 100) #5000ft for fig
                  
                  #set up empty dataframe to store results
                  kde_results <- data.frame()
                  
                  #for loop
                  for (b in bw) {
                    for (i in years) {
                      
                      #filter adjusted crime data by year
                      vcrime_filtered <- vcrime_final%>%
                        filter(year == i)
                      
                      #set philly extent (Philly county shapefile from census)
                      philly_owin <- as.owin(philly)
                      
                      #create ppp object from adjusted crime data
                      vcrime_ppp <- ppp(
                        x = st_coordinates(vcrime_filtered)[, 1],
                        y = st_coordinates(vcrime_filtered)[, 2],
                        window = philly_owin)
                      
                      #check for rejected points - outside philly extent 
                      x <- attr(vcrime_ppp, "rejects")
                      if (!is.null(x) && x$n > 0) {
                        print(paste("Rejected points for bandwidth:", 
                                    b, "| Year:", i, ":", x$n, "points"))
                      }
                      
                      #set bandwidth in terms of SD of quartic kernel
                      sig <- b/sqrt(7) 
                      
                      #perform kde
                      kde <- density.ppp(
                        x = vcrime_ppp,
                        sigma = sig,
                        kernel = "quartic",
                        diggle = T, #edge correction
                        eps = 100)
                      
                      kde_raster <- raster(kde)
                      crs(kde_raster) <- "EPSG:2272" #philly crs
                      
                      #transform density values from sq ft to sq mile 
                      kde_raster@data@values <- kde_raster@data@values*27878400 
                      
                      #extract estimated crime density values for each lot (centroid)
                      kde_values <- extract(kde_raster, lots)
                      
                      temp <- data.frame(
                        cluster = lots$cluster,
                        year = i,
                        treat = lots$community_grill,
                        n_crime = nrow(vcrime_filtered),
                        bw_value = b, #bandwidth size (in ft)
                        kde_value = kde_values) #crime density
                      
                      #bind to og kde_results dataframe
                      kde_results <- rbind(kde_results, temp)
                    }
                  }
                  
                  #-------------------------------------------------------------
                  # Objective 2: intervention effect estimation across bandwidths 
                  #-------------------------------------------------------------
                  
                  kde_results <- kde_results2%>%
                    mutate(post = ifelse(year > 2016, 1, 0))  #intervention indicator

                  #empty list to store model results across bandwidth values
                  bandwidth_results <- list()  
                  
                  #loop through estimation models for each bandwidth value 
                  for (b in bw) {
                    kde_results_filtered <- kde_results%>%
                      filter(bw_value == b)
                    
                    did_model<- glm(kde_value ~ treat*post, 
                                    data = kde_results_filtered%>%
                                       mutate(kde_value = round(kde_value,0)), #make integer 
                                    family = "poisson")
                    
                    bandwidth_results[[as.character(b)]] <- summary(did_model)
                  }
                  

                  #extract estimated intervention effects across bandwidths
                  bandwidth_effects <- data.frame(
                    bw_value = numeric(),
                    treat_post_coef = numeric(),
                    treat_post_se = numeric()
                  )
                  
                  for (b in bw) {
                    model_summary <- bandwidth_results[[as.character(b)]]
                    treat_post_row <- model_summary$coefficients["treat:post", ]
                    
                    bandwidth_effects <- rbind(bandwidth_effects, 
                                               data.frame(bw_value = b,
                                                          treat_post_coef = treat_post_row[1],
                                                          treat_post_se = treat_post_row[2]))
                  }
                  
                  #exponentiate to get relative risk 
                  bandwidth_effects <- bandwidth_effects%>%
                    mutate(exp_est = exp(treat_post_coef),
                           exp_ci_lower = exp(treat_post_coef - 1.96 * treat_post_se), 
                           exp_ci_upper = exp(treat_post_coef + 1.96 * treat_post_se))
                  
                  
                  #-------------------------------------------------------------
                  #               Objective 2: visualizations 
                  #-------------------------------------------------------------
                  
                  #plot change in estimated intervention effect by bandwidth size 
                  ggplot(bandwidth_effects, aes(x = bw_value, y = exp_est)) +
                    geom_point()+
                    geom_line() +
                    geom_line(aes(y=1), linetype = 2) +
                    geom_ribbon(aes(ymin = exp_ci_lower,
                                    ymax = exp_ci_upper), alpha = 0.2) +
                    labs(x = "Estimated effect (RR)",
                         y = "Bandwidth size (ft)")+
                    theme_minimal()
                  
                  