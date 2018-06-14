library(raster)
library(sf)
library(parallel)

states <- st_read("./data/UsElectionStates.shp")

## define function for filling rasters proportionally
buildRasters <- function(recordNum, outcellsize = 5000){

    ## make blank raster as defined by bounding box
    rast <- raster(vals = 0,
                   xmn = rbox["xmin"],
                   xmx = rbox["xmax"],
                   ymn = rbox["ymin"],
                   ymx = rbox["ymax"],
                   resolution = outcellsize)

    ## make all cells outside of polygon NA
    rama <- mask(rast, as(st_difference(record$geometry, buff), "Spatial"))

    ## extract cell ids from raster inside polygon
    raex <- extract(rast, as(st_difference(record$geometry, buff), "Spatial"), cellnumber = T)

    ## process proportions
    hc <- record$HC_percent/100
    dt <- record$DT_percent/100
    gj <- record$GJ_percent/100
    js <- record$JS_percent/100
    other <- 1 - sum(c(hc, dt, gj, js), na.rm = T)
    modList <- list(hc, dt, gj, js, other)

    ## assign proportions based on total valid cells
    cellVec <- unlist(mapply(function(x, y){
        x <- x[!is.na(x)]
        y <- y[!is.na(x)]
        rep(y, times = ceiling(x * length(rama[!is.na(rama)])))
    }, modList, seq(modList), SIMPLIFY = F))[1:length(rama[!is.na(rama)])]

    ## assign cell values to raster based on cell ID
    invisible(mapply(function(x, y) rama[x] <<- y, raex[[1]][,1], cellVec))
    return(rama)
}

## run build raster in parallel 
cl <- makeCluster(detectCores())
clusterExport(cl, list("states"))
clusterEvalQ(cl, {library(raster); library(sf); library(sp)})
stateRast <- parLapply(cl, seq(nrow(states)), buildRasters)
stopCluster(cl)

## merge each state into one raster
mergeVec <- 2:length(stateRast)
outRast <- extend(stateRast[[1]], extent(states))
lapply(mergeVec,
       function(x){
           outRast <<- merge(outRast, resample(stateRast[[x]], outRast))
       })

## export map
jpeg("./outImg.jpg", width = 1000, height = 600)
plot(outRast, col = c("lightblue", "salmon", "gold", "lightgreen", "purple"))
plot(states, add = T, col = NA)
dev.off()
