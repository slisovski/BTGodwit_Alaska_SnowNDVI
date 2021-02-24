
<style type="text/css">

body, td {
   font-size: 17px;
}
code.r{
  font-size: 12px;
}
pre {
  font-size: 12px
}
</style>

## Conklin, Lisovski & Battley 2021 - Environmental data and analysis

In Arctic-breeding shorebirds, timing of breeding closely follows the
retreat of snow-cover from tundra nest sites in the spring. We compared
two indices of breeding phenology (timing of snowmelt and spring
‘green-up’), summarized separately for two regions (North, South) of the
Alaska breeding range (Figure 1) of bar-tailed godwits. We did this for
two temporal periods: long-term trends encompassing the entire study
period (2008–2020), and a shorter term related directly to the period in
which we tracked individuals with geolocators (2008–2014).

<center>
<img src="images/Fig01_breedingRange.png"></img>
<figcaption>
Figure 1: Breeding range of Bar-tailed godwits in Alaska. Map data from
Natural Earth <https://www.naturalearthdata.com/>. Breeding range
supplied by BirdLife International and Handbook of the Birds of the
World (2017) Bird species distribution maps of the world. Version 7.0.
Available at <http://datazone.birdlife.org/species/requestdis>
</figcaption>
</center>

## The Datasets

### IMS Daily Northern Hemisphere Snow and Ice Analysis

Remotely sensed IMS Daily Northern Hemisphere Snow and Ice Analysis data
for the period 2008–2020 on a scale of 4 x 4 km were downloaded from the
National Snow & Ice Data Center \[1\].

1.  List all available downloaded IMS scenes (paths to files) and
    extract dates.

``` r
## list all downloaded IMS ASCI files
fls.gz <- list.files("~/Desktop/4km", pattern = ".gz", recursive = T,  full.names = T)

## get dates for ASCI files
dates  <- as.Date(as.POSIXct(unlist(lapply(strsplit(fls.gz, "ims"), function(x) {
  strsplit(x[[2]], "_4km")}))[c(TRUE, FALSE)], format = "%Y%j"))
```

1.  Initialize the file structure (projection, raster indices)

``` r
prj <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"

## read random file to initialize raster
rastID  <- sample(1:length(fls.gz), 1)
asciDat <- readLines(fls.gz[rastID])
## delete non-data
tab     <- asciDat[-which(unlist(suppressWarnings(lapply(asciDat, function(x) is.na(as.numeric(gsub(" ", "", x)))))))]
  
  z <- unlist(lapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]))) ## snow data to vector
  m <- matrix(z, ncol = 6144, nrow = 6144, byrow = T)[6144:1,]                   ## to matrix (and flip)                         

  snowR <- raster(m, crs = CRS(prj))                                             ## to raster
  extent(snowR)      <- c(-12288000, 12288000, -12288000, 12288000)              ## define extent

## get raster indices for cells inside the breeding range  
rInd <- extract(snowR, as(btg %>% st_transform(prj), "Spatial"), cellnumbers=TRUE)[[1]][,1]
  
# png("images/Fig02_init.png", width = 2000, height = 1000)
#   opar <- par(mfrow = c(1,2), mar = c(0,0,0,0), bty = "n")
#   plot(snowR, legend = FALSE, breaks = seq(-1, 4),
#      col = c("transparent", "lightsteelblue3", "palegreen4", "grey80", "grey99"),
#      xaxt = "n", yaxt = "n")
#   plot(btg %>% st_transform(prj) %>% st_geometry(), add = T, 
#      col = adjustcolor("purple4", alpha.f = 0.5), border = NA)
# 
#   plot(btg %>% st_transform(prj) %>% st_geometry())
#   points(coordinates(snowR)[rInd,], pch = 16, cex = 0.1)
#   par(opar)
# dev.off()
```

<center>
<img src="images/Fig02_init.png"></img>
<figcaption>
Figure 2: Random IMS scene (green = open land, white = snow covered
land, grey = ice covered ocean) in the left and breeding range polygon
with center coordinates of cells that intersect with the breeding range
on the right (n = 13,442 cells).
</figcaption>
</center>

1.  Create matrix with each pixel that intersects the breeding range in
    rows and for each scene in columns. Change `mclapply` to `lapply`
    (and delete the `mc.cores` option) on systems running windows.

``` r
snowM <- do.call("cbind", lapply(1:length(dates), function(x) {  
    
    tmp <- readLines(fls.gz[x])
    tab <- tmp[-which(unlist(suppressWarnings(parallel::mclapply(tmp, function(x) is.na(as.numeric(gsub(" ", "", x))), mc.cores = parallel::detectCores()-1))))]
    
    z <- unlist(parallel::mclapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]), mc.cores = parallel::detectCores()-1))
    m <- matrix(z, ncol = 6144, nrow = 6144, byrow = T)[6144:1,] 
    
    raster(m)[][rInd]
    
}))
  
snowRaw <- list(crds = coordinates(snowR)[rInd,], dates = dates, snow = snowM)
save(snowRaw, file = "results/snowRaw_4km_2004_2020.RData")
```

1.  Define (a) asymmetric gaussian model, (b) log-likelihood function
    and (c) function that calculates the dates for specified thresholds
    in the model prediction.

``` r
library(bbmle)
library(zoo)

### a. assymetric gaussian model
gauss.curve <- function(parms, tab) {
  t <- 1:nrow(tab)
  parms <- as.list(parms)
  fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
  fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
  c(fit1, fit2[-1])
}

### b. log-likelihood function with binomial error distribution.
fitGauss <- function(tab) {
  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), tab)  
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    -sum(dbinom(x = tab[,2]*100, size = rep(100, length(fit)), prob = fit, log = TRUE), na.rm=T)
  } 
  
  mle <- suppressWarnings(bbmle::mle2(gauss.loglik, method="L-BFGS-B",
                                      start=list(a1 = 200,     a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
                                      lower=list(a1 = 120,  a2 = 5,   a3 = 0.5,a4 = 5,  a5 = 0.5),
                                      upper=list(a1 = 240,  a2 = Inf,  a3 = Inf, a4 =  Inf, a5 =  Inf)
  ))
  
  coef(mle)
} 

### c. Curve intersection
curveIntersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {

  curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
  curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)

  point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                    c(min(curve1$x), max(curve1$x)))$root

  point_y <- curve2_f(point_x)
  
  list(x = point_x, y = point_y)
} 
```

1.  Fit asymmetric gaussian model and extract date of 1/3 snow (model
    prediction intersects falls below 66.66) free per pixel and year. In
    addition dates will be extracted when the pixel was snow free for
    the first time in each year and the first date when the pixel
    remained snow free for at least 60 days.

``` r
smM <- do.call("rbind", parallel::mclapply(1:nrow(snowRaw$crds), function(j) {

  ## only interested in snow free land vs. snow covered land (different years have different values e.g. 4 or 165 for snow free)
  y <- ifelse(snowRaw$snow[j,]%in%c(4,165), 4, ifelse(snowRaw$snow[j,]%in%c(3,164), 3, snowRaw$snow[j,]))

  if(sum(!y%in%c(2,4)) < length(y)/5) {

    tab0 <- data.frame(year = as.numeric(format(snowRaw$dates, "%Y")),
                       s    = ifelse(y%in%c(0,1), 0, ifelse(y==2, 0, 1)),
                       doy =  as.numeric(format(snowRaw$dates, "%j")))

    ## split by year
    spl <- split(tab0, f = as.character(tab0$year))

    sm1 <- do.call("rbind", lapply(spl, function(x) {
      
      tab <- merge(data.frame(day = 1:365), data.frame(day = as.numeric(x[,3]), p = as.numeric(x[,2])), all.x = T)

      mle <- fitGauss(tab)
      fit <- gauss.curve(mle, tab)

      sm <- tryCatch(curveIntersect(data.frame(x = tab[,1], y = fit)[1:mle[1],], data.frame(x = tab[,1], y = 0.666))$x,
                     error = function(x) NA)
      
      open0 <- sapply(which(tab$p<0.5), function(s) {

        if(any(!is.na(tab$p[s:nrow(tab)]) & tab$p[s:nrow(tab)]>0.5)) {
          nextSnow <- min(which(tab$p[s:nrow(tab)]>0.5))+s
        } else nextSnow <- nrow(tab)
        
        diff(c(s, nextSnow))

      })
      
      cbind(year = median(as.numeric(as.numeric(as.character(x[,1])))), sm = sm, 
            first = min(tab$day[which(tab$p<0.5)]), 
            open  = tab$day[which(tab$p<0.5)[min(which(open0>=60))]])

    }))

    matrix(unlist(c(merge(data.frame(year = 2004:2020), sm1, all.x = T)[,-1])), ncol = 3)

  } else {
    matrix(NA, nrow = length(2004:2020), ncol = 3)
  }

}, mc.cores = parallel::detectCores()-1))

snow <- list(crds = snowRaw$crds, smM = smM)
save(snow, file = "results/snowMelt_4km.RData")
```

### Noise-removed NDVI from NOAA STAR

### References

\[1\] IMS Daily Northern Hemisphere Snow and Ice Analysis at 1 km, 4 km,
and 24 km Resolutions, Version 1. <https://doi.org/10.7265/N52R3PMC>
