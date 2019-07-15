### Input for the IPBES version of Madingley (https://github.com/mikeharfoot/C-sharp-version-of-Madingley/tree/IPBES_scenarios)
### Create primary and secondary land cover fractions and their losses and gains 

# Author: Alan B. Seo (bumsuk.seo@kit.edu)

library(rgdal)
library(raster)

library(stringr)
library(gplots)
library(RColorBrewer)

library(doMC)
library(lattice)
library(RNetCDF)

path_base = "~/Dropbox/KIT_Modelling/"
# path_base = "~/pd/"

path_wd = paste0(path_base, "Madingley/CRAFTY_HANPP/")
path_RData =paste0(path_base, "Madingley/CRAFTY_HANPP_DATA/")



setwd(path_wd)


source(paste0(path_wd, "RScripts/0_CRAFTY_HANPP_Functions.R"))


# Read the transition matrix 
# Values give the biomass correction (i.e. the proportion of initial biomass remaining after the transition) (From at top, to at bottom).
tm_df = readxl::read_excel("../CRAFTY_HANPP_DATA/AFT transition matrix for Madingley_FromCalum2019Feb.xlsx", 1, col_types = "numeric")
tm_df = data.frame((tm_df[,-1]))
# cbind(colnames(tm_df), aft.shortnames.fromzero[c(2, 5, 15, 17, 12, 10, 8, 4 , 1, 3, 11, 6, 9, 7,13, 16, 14 )])
rownames(tm_df) = colnames(tm_df) = aft.shortnames.fromzero[c(2, 5, 15, 17, 12, 10, 8, 4 , 1, 3, 11, 6, 9, 7,13, 16, 14 )]

as.matrix(tm_df)
tm_df[is.na(tm_df)] = 1 

tm_ord_df = tm_df[aft.shortnames.fromzero, aft.shortnames.fromzero]

levelplot(x =  as.matrix(tm_df)  ,  col.regions = rev(terrain.colors(30)),ylab = "From", xlab="To",scales=list(x=list(rot=90)), main = "Retained standing stock (%)")


eu_xlim = c(0,20)
eu_ylim = c(44, 58)

crafty_dummy_r = raster(paste0(path_RData, "Output_CRAFTY_HANPP/HANPP_Harvest_CRAFTY_RCP2_6-SSP1_2016.tif"))


### Reproject CRAFTY hanpp 
madingley_global_rs =   raster(paste0(path_RData, "Output/Madingley_Fprimary_2005.tif"))

crafty_dummy_global_r = projectRaster(crafty_dummy_r, madingley_global_rs, method = "bilinear")
# plot(!is.na(crafty_dummy_global_r))

crafty_dummy_global_05deg_r = aggregate(crafty_dummy_global_r, 2)

plot(crafty_dummy_global_05deg_r)







indicator_selected = indicator.names[c(1:15, 17)]
indicator_selected_idx = match(indicator_selected, indicator.names)


targetyears <- 2016:2096 # 2016~2096 = 2010-2090 ?


n.targetyears <- length(targetyears)
# crafty.out.l <-  vector("list", n.targetyears)

path.batchrun = "~/Dropbox/KIT_Modelling/CRAFTY/EUpaper/Batchruns/NewRunsMay2019/AnnualRuns/"


test = getCRAFTYRaster("Paramset1", 1, 1, indicator_selected)
plot(test)


crafty_dummy_global_05deg_poly = rasterToPolygons(crafty_dummy_global_05deg_r)

# names(crafty.out.l) <- paste0("Year", targetyears)


crafty.out.2010_LL = getCRAFTYRaster("Paramset1", 1, 1, indicator_selected)
dummy.crafty.r = projectRaster(crafty.out.2010_LL,crs =  proj4.etrs_laea, method = "ngb", res = 1E4)


crafty_aft_factor = setValues(dummy.crafty.r[[16]], factor(getValues(dummy.crafty.r[[16]]), levels = 0:16, labels = aft.shortnames.fromzero))

plot(crafty_aft_factor)



hyde_hanpp_filename = paste0("../Data_IPBES_repository/Madingley_IPBES_INPUT/Data/SSP/SSP1RCP26/SSP1-RCP2.6HANPP_HYDEpop2006-2070.nc")

nc.f = RNetCDF::open.nc(hyde_hanpp_filename)
RNetCDF::print.nc(nc.f)
lon = var.get.nc(nc.f, "longitude")
lat = var.get.nc(nc.f, "latitude")
years <- var.get.nc(nc.f, "years")
fillValue = att.get.nc(nc.f, variable = "HANPPharvest", "_FillValue")



# sc_idx = 2 
# y_idx = 2
# lui_idx = 1 

# Get remained standing stock (%)
remainedStock = function(prev_lui, curr_lui) {
    
    
    
    # Lazy FR to Unmanaged Land
    prev_lui_idx = ifelse(prev_lui == -1, yes = 13 + 1 , no = prev_lui + 1 )
    curr_lui_idx = ifelse(curr_lui == -1, yes = 13 + 1, no = curr_lui + 1 )
    
    prev_lui_idx[is.na(prev_lui_idx)] = ncol(tm_ord_df) +1 
    curr_lui_idx[is.na(curr_lui_idx)] = nrow(tm_ord_df) +1 
    
    # str(prev_lui_idx)
    # str(curr_lui_idx)
    
    tm_ord_tmp = as.matrix(rbind(cbind(tm_ord_df, NA), NA))
    
    
    tm_ord_tmp  [curr_lui_idx, prev_lui_idx]
}


# remainedStock(-1, 0) # LazyFR to Ext_AF 
# remainedStock(13, 0) # UMF to Ext_AF 
# remainedStock(13, 1) # UMF to IA 
# remainedStock(13, 15) # UMF to VEP 
# 
# # 
# sc_idx = 5

# a dummy 0.5 deg grid
remained_global_array= array(data = NA, dim = c(1440, 720, n.targetyears))

for (sc_idx in 7:length(scens)) { 
    
    beginCluster()
    registerDoMC()
    
    crafty_remained_m_l = foreach (y_idx = 1:length(targetyears)) %dopar% { 
        
        curr_idx = y_idx
        prev_idx = max(1, y_idx-1) # for the first year 
        
        year = targetyears[curr_idx]
        year_prev = targetyears[prev_idx]
        
        print(year)
        # crafty.out.rs.tmp = crafty.out.l[[y_idx]]
        crafty_current_rs_tmp_LL = getCRAFTYRaster("Paramset1", sc_idx, curr_idx, indicator_selected ) # RCP4_5-SSP3
        # crafty_rs_tmp = projectRaster(crafty_out_rs_tmp_LL,crs  =  proj4.etrs_laea, method = "ngb", res = 1E4)
        
        crafty_prev_rs_tmp_LL = getCRAFTYRaster("Paramset1", sc_idx, prev_idx, indicator_selected ) # RCP4_5-SSP3
        # crafty_prev_rs_tmp  = projectRaster(crafty_prev_rs_tmp_LL,crs  =  proj4.etrs_laea, method = "ngb", res = 1E4)
        
        a = getValues(crafty_prev_rs_tmp_LL[[16]])
        b = getValues(crafty_current_rs_tmp_LL[[16]])
        
        
        print(table(a[a!=b], b[a!=b]))
        
        system.time({
            remained_tmp  = mapply(remainedStock, a, b)
        }
        )
        
        # system.time({
        #     remained_tmp  = mcmapply(remainedStock, a, b, mc.cores = detectCores())
        # })
        
        
        # str(c)
        # summary(c)
        # table(c)
        print( summary(remained_tmp))
        
        # crafty_remained_rs_tmp_LL = setValues(crafty_current_rs_tmp_LL[[16]], remained_tmp)
        # names(crafty_remained_rs_tmp_LL) = "RemainedStock_perc"
        # crafty_remained_rs_tmp_LL[crafty_remained_rs_tmp_LL==1] = NA
        # plot(crafty_remained_rs_tmp_LL)
        
        
        return(remained_tmp)
        
    }
    
    remained_tmp_rs_LL = setValues( brick(crafty.out.2010_LL), do.call("cbind", crafty_remained_m_l))
    
    remained_global_r = projectRaster(remained_tmp_rs_LL, madingley_global_rs, method = "bilinear")
    
    
    remained_global_r[ is.na(remained_global_r)] = NaN
    
    # to visually confirm
    # > plot(hanpp_harv_annual_rs[[1]])
    # > plot(hanpp_harvest_crafty_rel_global_r, add=T)
    

    remained_global_array[,,] = aperm(as.array(remained_global_r), perm = c(2,1,3))
    
    
    
    # Create the netCDF to hold all of the data
    out_nc_name = paste0(scens[sc_idx], "_RemainedStock_CRAFTY_2006-2086")
    nc <- create.nc(paste0(path_RData, "tmp/", out_nc_name, ".nc"))
    
    # Create vectors of data corresponding to latitude, longitude and month dimensions
    
    # Define the dimensions in the netCDF
    dim.def.nc(nc,"longitude",1440)
    dim.def.nc(nc,"latitude",720)
    dim.def.nc(nc,"years",length(targetyears))
    
    # Define the variables in the netCDF
    var.def.nc(nc,"longitude","NC_DOUBLE",dimensions = "longitude")
    var.def.nc(nc,"latitude","NC_DOUBLE",dimensions = "latitude")
    var.def.nc(nc,"years","NC_DOUBLE","years")
    var.def.nc(nc,"standingstock_loss_perc","NC_DOUBLE",c("longitude","latitude","years"))
     
    # Add the missing-data attribute to the temperature variable in the netCDF
    att.put.nc(nc,"standingstock_loss_perc","_FillValue","NC_DOUBLE", NaN)
     
    # Add the dimension data to the corresponding netCDF variables
    var.put.nc(nc,"latitude",lat,1,length(lat))
    var.put.nc(nc,"longitude",lon,1,length(lon))
    var.put.nc(nc,"years",targetyears,1,length(targetyears))
    var.put.nc(nc,"standingstock_loss_perc",remained_global_array)
     
    var.put.nc(nc,"years",2006:2086,1,length(targetyears))
 
    print.nc(nc)
    
    # Sync and close the netCDF
    sync.nc(nc)
    close.nc(nc)
    gc()
    
}


 


