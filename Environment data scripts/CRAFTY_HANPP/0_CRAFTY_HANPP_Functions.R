
# Lon-Lat projection 
proj4.LL <- CRS("+proj=longlat +datum=WGS84")

# Proj4js.defs["EPSG:3035"] etrs89/etrs-laea
# Scope: Single CRS for all Europe. Used for statistical mapping at all scales and other purposes where true area representation is required.
# Reference: http://spatialreference.org/ref/epsg/3035/
proj4.etrs_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs";


# Cell ID and cooridnates 
# ctry.ids <- read.csv("~/Dropbox/KIT/CLIMSAVE/IAP/Cell_ID_LatLong.csv")
# saveRDS(ctry.ids, file = "../crafty_web/GISData/ctry.ids.Rds")
ctry.ids = readRDS("~/Dropbox/KIT_Modelling/CRAFTY/CRAFTY_WEB/GISData/ctry.ids.Rds")
x.lat.v = sort(unique(ctry.ids$Longitude))
y.lon.v = sort(unique(ctry.ids$Latitude))


getRaster <- function(fname, indicator_selected) {
    
    # Target outcome
    result.tmp = read.csv2(fname, sep = ",")
    result.tmp$lon = x.lat.v[result.tmp$X]
    result.tmp$lat = y.lon.v[result.tmp$Y]
    
    # Create a spatial pixels data frame using the lon-lat table (Cell_ID_LatLong.csv) and the input data 
    result.spdf <- SpatialPixelsDataFrame(points = SpatialPoints(cbind(result.tmp$lon, result.tmp$lat), proj4string = proj4.LL), data = data.frame(result.tmp), tolerance = 0.0011)
    result_r = stack(result.spdf)[[4:22]]
    
    names(result_r) <- indicator.names
    
    result_r = result_r[[which(names(result_r) %in% indicator_selected)]]
    
    return(result_r)
}

getCRAFTYRaster = function(paramset, sc_idx, y_idx, indicator_selected) { 
    
    scenfiles.m = sapply(targetyears, FUN = function(year) paste0("-", 0:(length(scens)-1), "-",input.offset, "-EU-Cell-", year, ".csv"))
    paths.m <- sapply(1:ncol(scenfiles.m), FUN = function(x) paste0(path.batchrun, paramset, "/", scens, "/", scens, scenfiles.m[,x]))
    
    out_r = getRaster(paths.m[sc_idx, y_idx],  indicator_selected)
    return(out_r)
}



# read crafty files 

crafty.raster.files <- sort(list.files("~/Dropbox/KIT_Modelling/CRAFTY/Calibration/CRAFTY model runs_TIF/ETRS_LAEA/", pattern = "Baseline-0-0-EU-Cell.*\\.tif$", full.names = T)) 
targetyears <- 2016:2100
n.targetyears <- length(targetyears)
 
input.offset = "99"
n.paramset = 5
 
# paramset = "Paramset3"
paramsets = paste0("Paramset", 1:n.paramset)






crafty.csv.colnames = str_split("Tick,X,Y,Service:Meat,Service:Crops,Service:Diversity,Service:Timber,Service:Carbon,Service:Urban,Service:Recreation,Capital:Crop.productivity,Capital:Forest.productivity,Capital:Grassland.productivity,Capital:Financial.capital,Capital:Human.capital,Capital:Social.capital,Capital:Manufactured.capital,Capital:Urban.capital,LandUse,LandUseIndex,Agent,Competitiveness", pattern = ",")[[1]]
# crafty.layer.names = flatten_chr(crafty.csv.colnames[4:22])
# regexpr(crafty.layer.names, "[A-Z]" ) 
# str_replace_all(crafty.layer.names, "[^[:alnum:]]", " ")
crafty.layer.names = str_replace_all(crafty.csv.colnames, "[[:punct:]]", ".")[4:22]


indicator.names =  c("Service.Meat","Service.Crops","Service.Diversity",
                     "Service.Timber","Service.Carbon","Service.Urban",
                     "Service.Recreation","Crop.productivity","Forest.productivity",
                     "Grassland.productivity","Financial.capital","Human.capital",
                     "Social.capital","Manufactured.capital","Urban.capital",
                     "LandUse_notuse", "LandUseIndex","Agent_notuse", "Competitiveness")

scens <-c( "Baseline", "RCP2_6-SSP1","RCP2_6-SSP4","RCP4_5-SSP1","RCP4_5-SSP3","RCP4_5-SSP4","RCP8_5-SSP3","RCP8_5-SSP5")

scen.colours = rich.colors(length(scens))

# LUI AFT 
# -1  ? Lazy FR
#  0  2 Extensive agro-forestry mosaic	Ext_AF          
#  1  3 Intensive arable farming	    IA             Wheat
#                                                      maize
# rice? 
#  2  4 Intesive agro-forestry farming	Int_AF         
#  3  5 Intensive farming	            Int_Fa        starchyRoots
#  4  6 Intensive pastoral farming	    IP            monogastrics  / ruminants
#  5  7 Managed forestry	            MF            energycrops

#  6  8 Minimal management	            Min_man
#  7  9 Mixed farming	                Mix_Fa        oilcrops
#  8 10 Mixed forest	                Mix_For 
#  9 11 Mixed pastoral farming      	Mix_P
# 10 12 Multifunctional	                Multifun      pulses
# 11 13 Peri-urban	                    P-Ur
# 12 14 Unmanaged land          	    UL
# 13 15 Unmanaged forest	            UMF
# 14 16 Urban	                        Ur
# 15 17 Very extensive pastoral farming	VEP            
# 16  1 Extensive pastoral farming  	EP            pasture 

serviceNames <- c("Meat","Crops", "Diversity", "Timber", "Carbon", "Urban", "Recreation")
# serviceColours <- c("Meat" = "coral1", "Crops" = "goldenrod1", "Diversity" = "red", "Timber" = "tan4", "Carbon" = "darkgreen", "Urban" = "grey", "Recreation" = "orange")
serviceColours = c("Meat" = "coral1", "Crops" = "goldenrod1", "Diversity"="turquoise", "Timber" = "tan4","Carbon"="black", "Urban" = "grey","Recreation"="dodgerblue2")


# aft.colors = rich.colors(17)
aft.shortnames.fromzero <- c( "Ext_AF", "IA", "Int_AF", "Int_Fa", "IP", "MF", "Min_man", "Mix_Fa", "Mix_For", "Mix_P", "Multifun", "P-Ur", "UL", "UMF", "Ur", "VEP", "EP")


# aftNames <- c("IA","Int_Fa","Mix_Fa","Int_AF","Ext_AF","MF", "Mix_For","UMF","IP","EP","Mix_P","VEP","Multifun","Min_man","UL","P-Ur","Ur", "Lazy FR")

aft.colors.fromzero <-  (c("Ext_AF" = "yellowgreen", "IA"  = "yellow1", "Int_AF" =  "darkolivegreen1", "Int_Fa" = "lightgoldenrod1",  "IP" = "red1", "MF" =  "green3", "Min_man" = "lightyellow3",  "Mix_Fa" = "darkgoldenrod",  "Mix_For" = "green4",   "Mix_P" = "violetred",  "Multifun" = "blueviolet", "P-Ur"="lightslategrey", "UL" = "grey", "UMF" = "darkgreen", "Ur" = "black", "VEP" = "red4", "EP" = "red3")) # , "Lazy FR" = "black")

aft.names.fromzero <- c("Ext. agro-forestry","Int. arable","Int. agro-forestry","Int. mixed farming","Int. pastoral","Managed forest","Minimal management",
                        "Ext. mixed farming","Mixed forest","Mixed pastoral","Multifunctional","Peri-Urban", "Unmanaged land","Umanaged forest","Urban", "Very ext. pastoral","Ext. pastoral")

aft.col.breaks.fromzero = seq(-0.5, 16.5, 1)
# plot(crafty.out.l[[10]]$LandUse, crafty.out.l[[10]]$Agent)
# plot(crafty.out.l[[10]]$LandUseIndex, crafty.out.l[[10]]$Agent)