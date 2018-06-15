library(raster)
library(sf)
library(parallel)

states <- st_read("./data/UsElectionStates.shp")

## define function for building buffers based on turn out
buildBuffers <- function(recordNum, dist){
    message(recordNum)
    record <- states[recordNum,]
    buff <- st_cast(record$geometry, "MULTILINESTRING") %>%
        st_buffer(dist = dist) %>%
        st_union() %>%
        st_intersection(y = record$geometry)

## define function for filling rasters proportionally
buildRasters <- function(recordNum, outcellsize = 5000){
    message(recordNum)
    ## make blank raster as defined by bounding box
    record <- states[recordNum,]
    rbox <- st_bbox(record$geometry)
    rast <- raster(vals = 0,
                   xmn = rbox["xmin"],
                   xmx = rbox["xmax"],
                   ymn = rbox["ymin"],
                   ymx = rbox["ymax"],
                   resolution = outcellsize)

    ## make all cells outside of polygon NA
    rama <- mask(rast, as(record, "Spatial"))

    ## extract cell ids from raster inside polygon
    raex <- extract(rast, as(record, "Spatial"), cellnumber = T)[[1]]

    ## process proportions
    hc <- record$HC_prcn/100
    dt <- record$DT_prcn/100
    gj <- record$GJ_prcn/100
    js <- record$JS_prcn/100
    other <- 1 - sum(c(hc, dt, gj, js), na.rm = T)
    modList <- list(hc, dt, gj, js, other)

    ## assign proportions based on total valid cells
    cellVec <- unlist(mapply(function(x, y, rastLen){
        x <- x[!is.na(x)]
        y <- y[!is.na(x)]
        rep(y, times = round(x * rastLen))
    },
    modList, seq(modList),
    MoreArgs = list(rastLen = nrow(raex)), SIMPLIFY = F))

    ## replace extract table with values
    raex[,2] <- cellVec[1:nrow(raex)]

    ## replace cells
    rama[raex[,1]] <- raex[,2]
    
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
plot(outRast, col = c("dodgerblue", "firebrick2", "gold", "lightgreen", "darkorchid"))
plot(states["geometry"], add = T, col = NA)
dev.off()
