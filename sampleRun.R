# Read data
bucket.data <- read.csv("experimental-data.csv", stringsAsFactors = FALSE)
bucket.data$date <- as.POSIXct(bucket.data$Measurement.Time, format = "%m/%d/%Y %I:%M %p")
colnames(bucket.data)[which(names(bucket.data) == "PrecipitationIn")] <- "rainfall"

# Apply HyperSTL and plot result
source('HyperSTL.R')
y <- HyperSTL(bucket.data, 'D..10cm')
plot(y)
