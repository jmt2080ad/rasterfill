library(RCurl)
library(XML)
library(rlist)
library(parallel)

## get states
shpDir <- tempdir()
download.file("http://www2.census.gov/geo/tiger/GENZ2017/shp/cb_2017_us_state_20m.zip", file.path(shpDir, "US.zip"))
unzip(file.path(shpDir, "US.zip"), exdir = shpDir)
states <- read_sf(file.path(shpDir, "cb_2017_us_state_20m.shp"))

## get attributes
attr <- readHTMLTable(getURI("https://en.wikipedia.org/wiki/United_States_presidential_election,_2016"))
attr <- list.clean(attr, fun = is.null, recursive = FALSE)
attrMain  <- attr[[31]][,1:14]
attrPres  <- unlist(lapply(c("HC", "DT", "GJ", "JS"), rep, times = 3))
attrNames <- tolower(gsub(" ", "_", gsub("\n", " ", unlist(attrMain[1,]))))
attrNames <- c(attrNames[1:2], paste(attrPres, attrNames[3:length(attrNames)], sep = "_"))
attrMain  <- attrMain[2:nrow(attrMain),]
names(attrMain) <- gsub("%", "percent", gsub("#", "count", attrNames))
attrMain  <- cbind(attrMain[1:2],
                   data.frame(lapply(attrMain[3:ncol(attrMain)],
                                     function(x) as.numeric(gsub(",|%", "", x)))))
attrMain  <- attrMain[!grepl(", 1st|, 2nd|, 3rd", attrMain$state_or_district),]
attrMain$state_or_district <- gsub(" \\(at-lg\\)", "", attrMain$state_or_district) 

## merge attributes with shapes
states <- merge(states, attrMain, by.x="NAME", by.y="state_or_district")

## get only lower 48
states <- states[!states$NAME %in% c("Puerto Rico", "Alaska", "Hawaii"),]

## transfrom to US friendly projection 
states <- st_transform(states, crs = 2163)

## write out data
st_write(states, "./data/UsElectionStates.shp")
