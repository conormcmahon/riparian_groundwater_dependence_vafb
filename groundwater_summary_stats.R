
library(tidyverse)
library(raster)
library(janitor)
library(sp)

# Water Elevation Data
well_data <- read_csv("D:/SERDP/Vandenberg/hydrology/well_data/GroundwaterElevation.csv") %>%
  janitor::clean_names() %>%
  mutate(date = as.Date(msmt_date)) %>%
  filter(date > as.Date("1983-12-31"),
         !is.infinite(wse))

# Load Station Data 
site_information <- read_csv("D:/SERDP/Vandenberg/hydrology/well_data/Station.csv") %>%
  janitor::clean_names()

# Combine Latitude/Longitude with water elevation data
well_data <- merge(well_data, site_information[,c(1,6,7)], by=c("site_code"))

well_summary <- well_data %>%
  group_by(site_code) %>%
  summarize(start_date = min(date, na.rm=TRUE),
            end_date = max(date, na.rm=TRUE),
            min_elevation = min(wse, na.rm=TRUE),
            max_elevation = max(wse, na.rm=TRUE), 
            mean_elevation = mean(wse, na.rm=TRUE),
            stdev_elevation = sd(wse, na.rm=TRUE),
            ground_elevation = mean(gse_wse+wse, na.rm=TRUE),
            latitude = mean(latitude, na.rm=TRUE),
            longitude = mean(longitude, na.rm=TRUE),
            count = n()) %>%
  mutate(variation = max_elevation - min_elevation) %>%
  drop_na()

# Add well-based mean, min, max, stdev, etc. to well elevation dataset
well_data <- merge(well_data, well_summary[,c(1,4,5,6,7,11,12)], by="site_code")

well_summary_spatial <- well_summary
well_summary_spatial$start_year <- as.numeric(substr(well_summary_spatial$start_date,1,4))
well_summary_spatial$end_year <- as.numeric(substr(well_summary_spatial$end_date,1,4))
coordinates(well_summary_spatial) <- ~longitude+latitude
plot(well_summary_spatial)
projection(well_summary_spatial) <- CRS("+init=epsg:4267")
well_summary_spatial <- spTransform(x=well_summary_spatial, CRS("+init=epsg:32610"))

# Load Raster Background
class_raster <- raster("D:/SERDP/Vandenberg/Landsat/classified/classes_2019.tif")
riparian_raster_sub <- crop(class_raster==1, extent(well_summary_spatial))
# Plot raster and well points together
plot(riparian_raster_sub, col=(colorRampPalette(c("white","red")))(2))
plot(well_summary_spatial, add=TRUE)

ggplot(well_summary %>% filter(count > 10) %>% filter(wse != 0)) +
  geom_histogram(aes(x=max_elevation-min_elevation)) + 
#  scale_x_continuous(trans="log1p", expand=c(0,0), breaks=c(0,10,100,1000)) + 
  xlab("Water Table Elevation Variation (m)") +
  ylab("Frequency") + 
  ggtitle("Variation in Water Table Height at Individual Wells")

ggplot(well_summary) +
  geom_histogram(aes(x=count)) + 
  scale_x_continuous(trans="log1p", expand=c(0,0)) + 
  xlab("Number of Survey Datapoints at Well") +
  ylab("Frequency") + 
  ggtitle("Range in Number of Samples per Well")


# Which wells have the most records?
most_counts <- well_summary %>%
  arrange(-count)

# Plot of water elevation timeseries for four example wells (with lots of data)
ggplot(well_data %>% filter(site_code %in% most_counts[c(1,4,5,6),]$site_code)) + 
  geom_line(aes(x=date, y=(wse-mean_elevation)/stdev_elevation, group=site_code, col=site_code)) +
  xlab("Year") + 
  ylab("Water Elevation Z-score") + 
  ggtitle("Water Elevation Variation over Time")

# Plot of annual fluctuation in number of well datapoints available
ggplot(well_data) +
  geom_histogram(aes(x=as.numeric(substr(as.character(date),1,4)))) + 
  #  scale_x_continuous(trans="log1p", expand=c(0,0), breaks=c(0,10,100,1000)) + 
  xlab("Water Table Elevation Variation (m)") +
  ylab("Frequency") + 
  ggtitle("Variation in Water Table Height at Individual Wells")

