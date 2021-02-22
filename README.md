
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
‘green-up’), summarized separately for two regions (North, South) of
the Alaska breeding range (Figure 1) of bar-tailed godwits. We did this
for two temporal periods: long-term trends encompassing the entire study
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

<!-- end list -->

``` r
## list all downloaded IMS ASCI files
fls.gz <- list.files("~/IMS/4km", pattern = ".gz", recursive = T,  full.names = T)

## get dates for ASCI files
dates  <- as.Date(as.POSIXct(unlist(lapply(strsplit(fls.gz, "ims"), function(x) {
  strsplit(x[[2]], "_4km")}))[c(TRUE, FALSE)], format = "%Y%j"))
```

2.  Initialize the file structure (projection, raster indices)

<!-- end list -->

``` r
prj <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"

## read random file to initialize raster
rastID  <- sample(1:length(fls.gz), 1)
asciDat <- readLines(fls.gz[rastID])
## delete non-data
tab     <- asciDat[-which(unlist(suppressWarnings(lapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x)))))))]
  
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

3.  Create matrix with each pixel that intersects the breeding range in
    rows and for each scene in columns. Change `mclapply` to `lapply`
    (and delete the `mc.cores` option) on systems running windows.

<!-- end list -->

``` r
snowM <- do.call("cbind", parallel::mclapply(1:length(dates), function(x) {  
    
    tmp <- suppressWarnings(readLines(fls.gz[x]))
    ind  <- tmp[-which(unlist(suppressWarnings(lapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x)))))))]
    
    z <- unlist(lapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]])))
    m <- matrix(z, ncol = 6144, nrow = 6144, byrow = T)[6144:1,]  
    
    raster(m)[][rInd]
    
}, mc.cores = parallel::detectCores()-1))
  
# snowRaw <- list(crds = coordinates(snowR)[rInd], dates = dates, snow = snowM)
save(snowRaw, file = "results/snowRaw_4km_2004_2020.RData")
```

4.  Define (1) asymmetric gaussian model, (2) log-likelihood function
    and (3) function that calculates the dates for specified thresholds
    in the model prediction.

<!-- end list -->

``` r
library(bbmle)
library(zoo)

### 1. assymetric gaussian model
gauss.curve <- function(parms, tab) {
  t <- 1:nrow(tab)
  parms <- as.list(parms)
  fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
  fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
  c(fit1, fit2[-1])
}

### 2. log-likelihood function with binomial error distribution.
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

### 3. Curve intersection
curveIntersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {

  curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
  curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)

  point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                    c(min(curve1$x), max(curve1$x)))$root

  point_y <- curve2_f(point_x)
  
  list(x = point_x, y = point_y)
} 
```

7.  Fit asymmetric gaussian model and extract date of 1/3 snow (model
    prediction intersects falls below 66.66) free per pixel and year. In
    addition dates will be extracted when the pixel was snow free for
    the first time in each year and the first date when the pixel
    remained snow free for at least 60 days.

### Noise-removed NDVI from NOAA STAR

### References

\[1\] IMS Daily Northern Hemisphere Snow and Ice Analysis at 1 km, 4 km,
and 24 km Resolutions, Version 1. <https://doi.org/10.7265/N52R3PMC>
