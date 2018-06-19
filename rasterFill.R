library(raster)
library(sf)
library(parallel)
library(data.table)

states <- st_read("./data/UsElectionStates.shp")

## define function for building buffers based on turn out
## root finding equation: (1 - buffArea/stateArea) - turnout_ppo
rootFun <- function(x, record, stateArea, turnout_ppo){
    buffArea <- st_cast(record$geometry, "MULTILINESTRING") %>%
        st_union() %>%
        st_buffer(dist = x) %>%
        st_intersection(y = record$geometry) %>%
        st_area() %>%
        as.numeric()
    (1 - buffArea/stateArea) - turnout_ppo
}

buildBuff <- function(record){
    stateArea   <- as.numeric(st_area(record$geometry))
    turnout_ppo <- record$turn_pp
    buffDist    <- uniroot(rootFun,
                           record = record,
                           stateArea = stateArea,
                           turnout_ppo = turnout_ppo,
                           lower = 1,
                           upper = stateArea)
    return(buffArea <- st_cast(record$geometry, "MULTILINESTRING") %>%
               st_union() %>%
               st_buffer(dist = buffDist$root))
}

## define function for filling rasters proportionally
buildRasters <- function(recordNum, outcellsize = 1000){
    message(recordNum)
    ## make blank raster as defined by bounding box
    record <- states[recordNum,]
    buff <- buildBuff(record)
    stat <- st_difference(record$geometry, buff)
    rbox <- st_bbox(record$geometry)
    rast <- raster(vals = NA,
                   xmn = rbox["xmin"],
                   xmx = rbox["xmax"],
                   ymn = rbox["ymin"],
                   ymx = rbox["ymax"],
                   resolution = outcellsize)

    ## extract cell ids from raster inside polygon
    raex <- extract(rast, as(stat, "Spatial"), cellnumber = T)[[1]]

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
    rast[raex[,1]] <- raex[,2]
    
    return(rast)
}

## run build raster in parallel 
cl <- makeCluster(detectCores())
clusterExport(cl, list("states", "buildBuff", "rootFun"))
clusterEvalQ(cl, {library(raster); library(sf); library(sp)})
stateRast <- parLapply(cl, seq(nrow(states)), buildRasters)
stopCluster(cl)

## merge each state into one raster
mergeVec <- 2:length(stateRast)
outRast <- extend(stateRast[[1]], extent(states))
lapply(mergeVec,
       function(x){
           outRast <<- merge(outRast, resample(stateRast[[x]], outRast))
           gc()
       })


## export map
jpeg("./outImg.jpg", width=900, height=600, quality=100)
cols <- c("dodgerblue", "firebrick2", "gold", "green", "darkorchid")
plot(outRast, col = cols, legend = F, axes = F)
plot(states["geometry"], add = T, col = NA, lwd = 0.5)
legend(x="bottomleft", legend = c("Dems", "GOP", "Libs", "Greens", "Others"), fill= cols)
dev.off()

writeRaster(outRast, "outRast.tif", overwrite = T)
