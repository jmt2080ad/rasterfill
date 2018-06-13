library(raster)
library(sf)
library(RCurl)
library(XML)
library(rlist)
library(parallel)
library(data.table)

#### get data
## states
shpDir <- tempdir()
download.file("http://www2.census.gov/geo/tiger/GENZ2017/shp/cb_2017_us_state_20m.zip", file.path(shpDir, "US.zip"))
unzip(file.path(shpDir, "US.zip"), exdir = shpDir)
states <- read_sf(file.path(shpDir, "cb_2017_us_state_20m.shp"))

## attributes, takes forever on my machine...
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

## merge meta data with shapes
states <- merge(states, attrMain, by.x="NAME", by.y="state_or_district")

## get only lower 48
states <- states[!states$NAME %in% c("Puerto Rico", "Alaska", "Hawaii"),]

## transfrom to utm 
states <- st_transform(states, crs = 2163)

## now on to the symbology
buildRasters <- function(recordNum, outcellsize = 5000){
    record <- states[recordNum,]
    print(recordNum)
    buff <- st_buffer(st_cast(record$geometry, "MULTILINESTRING"), 3000)
    stat <- st_difference(record$geometry, buff)
    rbox <- st_bbox(record$geometry)

    rast <- raster(vals = 0,
                   xmn = rbox["xmin"],
                   xmx = rbox["xmax"],
                   ymn = rbox["ymin"],
                   ymx = rbox["ymax"],
                   resolution = outcellsize)

    rama <- mask(rast,
                 as(st_difference(record$geometry, buff), "Spatial"))

    raex <- extract(rast,
                    as(st_difference(record$geometry, buff), "Spatial"),
                    cellnumber = T)

    hc <- record$HC_percent/100
    dt <- record$DT_percent/100
    gj <- record$GJ_percent/100
    js <- record$JS_percent/100
    other <- 1 - sum(c(hc, dt, gj, js), na.rm = T)
    modList <- list(hc, dt, gj, js, other)

    cellVec <- unlist(mapply(function(x, y){x <- x[!is.na(x)]; y <- y[!is.na(x)]; rep(y, times = ceiling(x * length(rama[!is.na(rama)])))},
                             modList, seq(modList), SIMPLIFY = F))[1:length(rama[!is.na(rama)])]
   
    invisible(mapply(function(x, y) rama[x] <<- y, raex[[1]][,1], cellVec))
    return(rama)
}

cl <- makeCluster(detectCores())
clusterExport(cl, list("states"))
clusterEvalQ(cl, {library(raster); library(sf); library(sp)})
stateRast <- parLapply(cl, seq(nrow(states)), buildRasters)
stopCluster(cl)

mergeVec <- 2:length(stateRast)

outRast <- extend(stateRast[[1]], extent(states))

lapply(mergeVec,
       function(x){
           outRast <<- merge(outRast, resample(stateRast[[x]], outRast))
       })

jpeg("./outImg.jpg", width = 1000, height = 600)
plot(outRast, col = c("lightblue", "salmon", "gold", "lightgreen", "purple"))
plot(states, add = T, col = NA)
dev.off()
