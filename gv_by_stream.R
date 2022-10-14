
library(tidyverse)
library(rgdal)
library(raster)
library(ncdf4)

# Extract some basic stats from GV phenology image
getAnnualMeanGV <- function(gv_pheno)
{
  mean <- calc(gv_pheno, fun=mean)
  #max <- calc(gv_pheno, fun=max)
  #min <- calc(gv_pheno, fun=min)
  #max_month <- calc(gv_pheno, fun=which.max)
  #min_month <- calc(gv_pheno, fun=which.min)
  
  #return( stack(mean, max, min, max_month, min_month) )
  return(mean)
}

# Masks INPUT_IMAGE to riparian vegetation around one target stream
#   STREAM_NAME - string name for target stream, assumes in NHD+ Flowlines variable "GNIS_Name"
#   RIPARIAN_MASK - 1 at riparian cells, NA at non-riparian cells
#   FLOWLINES - SpatialLineDataFrame object containing vector info on streams, assuming NHD+ format
maskToStream <- function(stream_name, riparian_mask, flowlines, input_image, buffer_width)
{
  target_stream <- flowlines[flowlines$GNIS_Name == stream_name,]
  target_stream_buffer <- buffer(target_stream, width=buffer_width)
  output_image <- mask(input_image, target_stream_buffer) * riparian_mask
  output_image <- crop(output_image, target_stream_buffer)
  return(output_image)
}

# Process One Year to phenoseries and mask to riparian area
#    save last year to use in differencing
last_phenology_image <- raster()
maskAnnualRaster <- function(year, stream_name, riparian_mask, flowlines, get_change=FALSE, last_stream_raster=raster(), buffer_width=1000)
{
  # Load image from target year
  phenology_image <- stack(paste("D:/landsat_vafb/pheno_ndvi_vafb/pheno_NDVI_",
                                 year,
                                 "_annual.tif", 
                                 sep=""))
  # Mask to Riparian Areas
  stream_data <- maskToStream(stream_name, riparian_mask, flowlines, phenology_image, buffer_width)
  # Generate Monthly Mean Values across Riparian
  monthly_data <- data.frame(year = numeric(0),
                             month = numeric(0),
                             mean_gv = numeric(0),
                             median_gv = numeric(0),
                             max_gv = numeric(0),
                             min_gv = numeric(0),
                             sd_gv = numeric(0),
                             mean_gv_change = numeric(0),
                             median_gv_change = numeric(0),
                             max_gv_change = numeric(0),
                             min_gv_change = numeric(0),
                             sd_gv_change = numeric(0),
                             stream_name = character(0))
  for(ind in 1:nlayers(stream_data))
  {
    mean_gv <- cellStats(stream_data[[ind]], stat=mean)
    median_gv <- cellStats(stream_data[[ind]], stat=median)
    max_gv <- cellStats(stream_data[[ind]], stat=max)
    min_gv <- cellStats(stream_data[[ind]], stat=min)
    sd_gv <- cellStats(stream_data[[ind]], stat=sd)
    print(paste(month.name[[ind]], " ", year,
                " mean GV across ",
                stream_name, 
                " is ",
                mean_gv),
          sep="")
    if(get_change)
    {
      mean_gv_change <- cellStats(stream_data[[ind]] - last_stream_raster[[ind]], stat=mean)
      median_gv_change <- cellStats(stream_data[[ind]] - last_stream_raster[[ind]], stat=median)
      max_gv_change <- cellStats(stream_data[[ind]] - last_stream_raster[[ind]], stat=max)
      min_gv_change <- cellStats(stream_data[[ind]] - last_stream_raster[[ind]], stat=min)
      sd_gv_change <- cellStats(stream_data[[ind]] - last_stream_raster[[ind]], stat=sd)
    }
    else 
    {
      mean_gv_change <- NA
      median_gv_change <- NA
      max_gv_change <- NA
      min_gv_change <- NA
      sd_gv_change <- NA
    }
    print(paste(month.name[[ind]], " ", year,
                " mean GV change is ",
                mean_gv_change),
          sep="")
    monthly_data <- rbind(monthly_data,
                    data.frame(year = year,
                               month = ind,
                               mean_gv = mean_gv,
                               median_gv = median_gv,
                               max_gv = max_gv,
                               min_gv = min_gv,
                               sd_gv = sd_gv,
                               mean_gv_change = mean_gv_change,
                               median_gv_change = median_gv_change,
                               max_gv_change = max_gv_change,
                               min_gv_change = min_gv_change,
                               sd_gv_change = sd_gv_change,
                               stream_name = stream_name))
  }
  # Generate annual mean GV values for each pixel
  stream_mean_gv <- getAnnualMeanGV(stream_data)
  return(list(stream_mean_gv,
              monthly_data,
              stream_data))
}

# Load GV Image
gv_pheno_scene <- stack("D:/landsat_vafb/data/stacked_scenes/pheno_NDVI_1986_annual.tif")
# Flowlines Vector File 
flowlines <- readOGR("D:/SERDP/Vandenberg/hydrology/derived/flowlines_clipped_narrow.shp")
flowlines <- spTransform(flowlines, crs(gv_pheno_scene))
# Subset to a few major streams we are interested in
target_streams <- c("Santa Ynez River", "Santa Maria River","San Antonio Creek")
flowlines <- flowlines[flowlines$GNIS_Name %in% target_streams,]

# Choose a Riparian Mask 
#   Classifier was run in timeseries
#   Keep something as 'riparian' if it was classified that way in at least 90% of models
classified_years <- c(1984,1986,1988,1995,1996,1998,2000,2002,2013,2015,2017,2019)
riparian_mask <- (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
                         classified_years[1],
                         ".tif",
                         sep=""))==1)
for(ind in classified_years)
{
  riparian_mask <- riparian_mask +
                  (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
                                ind,
                                ".tif",
                                sep=""))==1)
}
riparian_mask <- riparian_mask >= 10
riparian_mask[riparian_mask == 0] <- NA

# Process All Scenes
years <- c(seq(1984,2011),seq(2013,2021))
monthly_data <- data.frame(year = numeric(0),
                           month = numeric(0),
                           mean_gv = numeric(0),
                           median_gv = numeric(0),
                           stream_name = character(0))
annual_gv_raster_list <- list()
previous_stream_image <- raster()
for(stream in target_streams)
{
  annual_gv_raster_new <- raster()
  for(year in years)
  {
    # Pick Which Years for Riparian Class Mask
    #   Originally tried this to flexibly change the riparian extent between years
    #   This approach is fraught... for now, pick only the consistently-riparian areas 
    #     to use for the whole timeseries (above, outside loop)
    # class_index <- max(which(year >= classified_years)) # last classified year BEFORE target
    # riparian_mask <- (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
    #                                classified_years[class_index],
    #                                ".tif",
    #                                sep=""))==1)
    # if(class_index < length(classified_years))
    # {
    #   riparian_mask <- riparian_mask + 
    #                    (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
    #                                 classified_years[class_index+1],
    #                                 ".tif",
    #                                 sep=""))==1)
    # }
    # riparian_mask[riparian_mask==0] <- NA
    
    # For first year, no difference image is made
    if(which(year == years) < 2)
    {
      new_data <- maskAnnualRaster(year, stream, riparian_mask, flowlines, get_change=FALSE, last_stream_raster=raster(), buffer_width=1000)
    }
    else
    {
      new_data <- maskAnnualRaster(year, stream, riparian_mask, flowlines, get_change=TRUE, last_stream_raster=previous_stream_image, buffer_width=1000)
    }
    previous_stream_image <- new_data[[3]]
    annual_gv_raster_new <- addLayer(annual_gv_raster_new, new_data[[1]])
    monthly_data <- rbind(monthly_data, new_data[[2]])
  }
  annual_gv_raster_list <- c(annual_gv_raster_list, annual_gv_raster_new)
}



# Load SPEI information
spei <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_12.csv")
names(spei) <- c("date", "spei_12")
# 1-month
spei_1 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_1.csv")
spei$spei_1 <- spei_1$spei
# 2-month
spei_2 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_2.csv")
spei$spei_2 <- spei_2$spei
# 3-month
spei_3 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_3.csv")
spei$spei_3 <- spei_3$spei
# 4-month
spei_4 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_4.csv")
spei$spei_4 <- spei_4$spei
# 5-month
spei_5 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_5.csv")
spei$spei_5 <- spei_5$spei
# 6-month
spei_6 <- read_csv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/spei_6.csv")
spei$spei_6 <- spei_6$spei
# Reformat date to Day, Month, Year, Day of Year
spei <- spei %>%
  mutate(year = as.numeric(substr(date,1,4)),
         month = as.numeric(substr(date,6,7)),
         day = as.numeric(substr(date,9,10)))


# Annual Summaries
#   GV
gv_annual <- monthly_data %>% 
  filter(month %in% seq(3,9),       # Keep only growing season
         ) %>% 
  group_by(year, stream_name) %>% 
  summarize(gv_mean = mean(mean_gv),
            gv_max = max(mean_gv),
            gv_mean_change = mean(mean_gv_change),
            gv_max_change = max(mean_gv_change))
#   SPEI
spei_annual <- spei %>%
  filter(year %in% gv_annual$year,   # Filter to target years  
         month %in% c(3,9)) %>%           # Filter to growing season (Mar. - Sep.)
  group_by(year) %>%
  summarize(spei_12_end = sum(spei_12*(month %in% c(8,9))/2),
            spei_1 = mean(spei_1),
            spei_2 = mean(spei_2),
            spei_3 = mean(spei_3),
            spei_4 = mean(spei_4),
            spei_5 = mean(spei_5),
            spei_6 = mean(spei_6),
            spei_12 = mean(spei_12))
# Combine SPEI and GV
annual_summary <- merge(gv_annual, spei_annual, by="year")

# Load Streamflow
san_antonio_flow <- read_tsv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/san_antonio_annual_flow", skip=34)
santa_maria_flow <- read_tsv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/san_antonio_annual_flow", skip=34)
santa_ynez_flow <- read_tsv("D:/SERDP/scripts/riparian_groundwater_dependence_vafb/hydrology/san_antonio_annual_flow", skip=34)
san_antonio_flow <- san_antonio_flow[2:nrow(san_antonio_flow),] %>%
  filter(year_nu > 1983) %>%
  mutate(mean_va = as.numeric(mean_va))
santa_maria_flow <- santa_maria_flow[2:nrow(santa_maria_flow),] %>%
  filter(year_nu > 1983) %>%
  mutate(mean_va = as.numeric(mean_va))
santa_ynez_flow <- santa_ynez_flow[2:nrow(santa_ynez_flow),] %>%
  filter(year_nu > 1983) %>%
  mutate(mean_va = as.numeric(mean_va))
names(san_antonio_flow) <- c("agency_cd", "site_no", "parameter_cd", "ts_id", "year", "flow")
names(santa_maria_flow) <- c("agency_cd", "site_no", "parameter_cd", "ts_id", "year", "flow")
names(santa_ynez_flow) <- c("agency_cd", "site_no", "parameter_cd", "ts_id", "year", "flow")



# ******************* Linear Models and Statistics *******************
#  SYR
santa_ynez_annual <- merge(annual_summary %>% filter(stream_name == "Santa Ynez River"), 
                           santa_ynez_flow, by="year")
santa_ynez_model <- lm(data=santa_ynez_annual %>% drop_na(), gv_mean_change ~ spei_12_end+flow)
summary(santa_ynez_model)
santa_ynez_model_results <- data.frame(gv = santa_ynez_annual$gv_mean_change,
                                       spei = santa_ynez_annual$spei_12_end,
                                       flow = santa_ynez_annual$flow,
                                       gv_pred = predict(santa_ynez_model, newdata=santa_ynez_model))
ggplot(santa_ynez_model_results %>% drop_na()) + 
  geom_point(aes(x=gv,y=gv_pred)) + 
  theme_bw() + 
  xlab("Greenness Observation") + 
  ylab("Greenness Prediction") + 
  geom_abline(intercept=0, slope=1, col="gray") + 
  ggtitle("Santa Ynez - Linear Model Prediction Results")
#  SMR
santa_maria_annual <- merge(annual_summary %>% filter(stream_name == "Santa Maria River"), 
                            santa_maria_flow, by="year")
santa_maria_model <- lm(data=santa_maria_annual %>% drop_na(), gv_mean_change ~ spei_12_end+flow)
summary(santa_maria_model)
santa_maria_model_results <- data.frame(gv = santa_maria_annual$gv_mean_change,
                                       spei = santa_maria_annual$spei_12_end,
                                       flow = santa_maria_annual$flow,
                                       gv_pred = predict(santa_maria_model, newdata=santa_maria_annual))
ggplot(santa_maria_model_results %>% drop_na()) + 
  geom_point(aes(x=gv,y=gv_pred)) + 
  theme_bw() + 
  xlab("Greenness Observation") + 
  ylab("Greenness Prediction") + 
  geom_abline(intercept=0, slope=1, col="gray") + 
  ggtitle("Santa Maria - Linear Model Prediction Results")
#  SAC
#   Note - for this one SPEI_3 seems to be helping meaningfully (more so than for others, too)
san_antonio_annual <- merge(annual_summary %>% filter(stream_name == "San Antonio Creek"), 
                            san_antonio_flow, by="year")
summary(lm(data=san_antonio_annual %>% drop_na(), gv_mean_change ~ spei_12_end+flow))
#  SYR + SMR
#   And actually SPEI_3 helps this one too
syr_smr_annual <- rbind(santa_ynez_annual, santa_maria_annual)
syr_smr_model <- lm(data=syr_smr_annual, gv_mean_change ~ spei_12_end+flow+stream_name)
summary(syr_smr_model)
santa_maria_model_results <- data.frame(gv = syr_smr_annual$gv_mean_change,
                                        spei = syr_smr_annual$spei_12_end,
                                        flow = syr_smr_annual$flow,
                                        stream = syr_smr_annual$stream_name,
                                        gv_pred = predict(syr_smr_model, newdata=syr_smr_annual))
ggplot(santa_maria_model_results %>% drop_na()) + 
  geom_point(aes(x=gv,y=gv_pred)) + 
  theme_bw() + 
  xlab("Greenness Observation") + 
  ylab("Greenness Prediction") + 
  geom_abline(intercept=0, slope=1, col="gray") + 
  ggtitle("Santa Maria - Linear Model Prediction Results")

# ******************* Linear Models and Statistics *******************





# ******************* Visualization *******************
#   Changes in GV over Time by River
ggplot(data=monthly_data %>% filter(year %in% seq(2014,2018))) + 
  geom_line(aes(x=month, y=mean_gv, group=year, col=year)) + 
  facet_wrap(~stream_name) + 
  xlab("Month") + 
  scale_x_continuous(breaks=seq(1,6)*2, expand=c(0,0)) + 
  ylab("Mean Greenness") + 
  theme_bw() + 
  ggtitle("Seasonal Phenology of Riparian Plants")
# Annual change
#   For some reason, all three of these show really strong positive change in GV with time
#   Not sure if this is due to recovery from some major event?
#   Or maybe there's something wrong with my data? 
#   Difference between Landsat 5 and 8, maybe? idk
ggplot(data=annual_summary) + 
  geom_line(aes(x=year, y=gv_mean_change),col="green3") + 
  geom_line(aes(x=year, y=spei_growing/100),col="red") +
  facet_wrap(~stream_name)




# ******************* Monthly Modeling *******************
monthly_subset <- merge(monthly_data, spei, by=c("year","month")) %>%
  filter(stream_name %in% c("Santa Ynez River","Santa Maria River"))
# Add in a single annual SPEI value for each year
annual_spei_simple <- spei_annual %>% 
  dplyr::select(year,spei_12_end)
monthly_subset <- merge(monthly_subset, annual_spei_simple, by="year")
monthly_subset$month2 <- monthly_subset$month^2
monthly_subset$month3 <- monthly_subset$month^3
monthly_model <- lm(data=monthly_subset, mean_gv ~ spei_1+stream_name+month+month2+month3)
summary(monthly_model)
monthly_subset$pred <- predict(monthly_model, newdata=monthly_subset)
santa_maria_model_results <- data.frame(gv = syr_smr_annual$gv_mean_change,
                                        spei = syr_smr_annual$spei_12_end,
                                        flow = syr_smr_annual$flow,
                                        stream = syr_smr_annual$stream_name,
                                        gv_pred = predict(monthly_model, newdata=syr_smr_annual))
ggplot(santa_maria_model_results %>% drop_na()) + 
  geom_point(aes(x=gv,y=gv_pred)) + 
  theme_bw() + 
  xlab("Greenness Observation") + 
  ylab("Greenness Prediction") + 
  geom_abline(intercept=0, slope=1, col="gray") + 
  ggtitle("Santa Maria - Linear Model Prediction Results")



# SPEI Visualizationggplot(spei %>% filter(year > 2010)) + 
geom_line(aes(x=date, y=spei_12), col="red") + 
  geom_line(aes(x=date, y=spei_3), col="cyan3") + 
  ggtitle("SPEI 3- and 12-month Trends on Vandenberg") + 
  xlab("Year") + 
  ylab("SPEI") + 
  theme_bw()




