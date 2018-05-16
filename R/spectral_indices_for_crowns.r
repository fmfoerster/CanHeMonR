#' @title Calculate Spectral Indices For Individual Tree Crowns
#' @description Calculate one or more spectral indices for individual
#' tree crowns depicted in a polygon shapefile and a provided raster image and write the output optionally to a .shp
#' @param crown_SpatPoldf SpatialPolygonsDataFrame depicting crowns
#' @param extract_as_points Logical. Should the SpatialPolygons be converted to SpatialPoints before extracting pixel values?
#' This can be useful and save lots of time when the polygons are smaller than the cel size of the image. Default is FALSE.
#' @param r_fname Raster file of the image to calculate the spectral index from.
#' The bands of the raster should be named following the Quantalab convention *X700.000000.Nanometers*
#' @param index_names Character vector with the names of a spectral indices implemented as a function in CanHeMonR. e.g. "G_over_B".
#' If left to NULL (the default) and all_bands is T, then only the average band values, and no indices are returned. This is useful when
#' extracting values from a DEM image for example.
#' @param all_bands Logical. Should individual band values be returned? Default = T
#' @param bandnames Character. In case the bands aren't named according to wavelength and following csic convention, they
#' can be provided. Default is NULL in which cases bandnames are read from the image file and csic naming convention is assumed.
#' @param shp_outname Filename for a shapefile of the output. If NULL, no output file is written away. Default is NULL.
#' @return A point shapefile with attribute table, and the points placed at the centroid of each crown.
#' The attribute table of the shp file maintains the attribute table of the crown but adds a column for each spectral index.
#' @export
spectral_indices_for_crowns <- function(crown_SpatPoldf, r_fname, index_names = NULL,
                                        bandnames = NULL, all_bands = T,extract_as_points = F,
                                        shp_outname = NULL){

  raster::rasterOptions(progress = "text")

  #read in the crown shapefile
  #crown_pols <- raster::shapefile(crown_shp)
  #read in the image
  r <- raster::brick(r_fname)
  #assign the explicitly provided bandnames if necessary
  if (!is.null(bandnames)){
    if (raster::nlayers(r) != length(bandnames)){
      cat('You should provide as many bandnames as the image has layers\n!')
      browser()
    }else{
      names(r) <- bandnames
    }
  #read the bandnames contained in the image file, and if indices are requested
  }else{
    #check that the layernames of the raster reflect the wavelengths
    if (!is.null(index_names)){
      if (length(grep("Nano", names(r))) == 0){
        cat('spectral_indices_for_crowns expects that the image layernames indicate wavelengths \n
            This is not the case for ', r_fname,'\n')
        browser()
      }
    }
  }
  #extract individual band values for each crown
  if (extract_as_points){
    cat('Converting polygons to points before extracting image values.\n')
    #treating each crown as a point (this might save time when the image pixels are larger than crowns)
    crown_SpatPointdf <- rgeos::gCentroid(crown_SpatPoldf, byid = T)
    r_crown_vals <- raster::extract(r, crown_SpatPointdf)
  }else{
    #averaging the band values across each crown
    r_crown_vals <- raster::extract(r,crown_SpatPoldf,fun = mean, na.rm=T)
  }
#weights=T, small=T

  #get the additional attributes to append output to
  r_crown_vals_df <- slot(crown_SpatPoldf,"data")
  if (!is.null(index_names)){
    #remove any preexisting columns that coincide with the requested indices!
    pre_existing_indices <- index_names[index_names %in% colnames(r_crown_vals_df)]
    if (length(pre_existing_indices) > 0){
      cat('The following indices were already present in the attributes and will be overwritten!:\n')
      print(pre_existing_indices)
      r_crown_vals_df <- dplyr::select(r_crown_vals_df,-dplyr::one_of(index_names))
    }
  }

  #if requested, add the individual band values
  if (all_bands){
    r_crown_vals_df <- cbind.data.frame(r_crown_vals_df,as.data.frame(r_crown_vals))
  }


  if (!is.null(index_names)){
    #calculate and add a spectral index
    for (index_name in index_names){
      calculated_index <- get(index_name)(df = r_crown_vals)
      r_crown_vals_df <- cbind.data.frame(r_crown_vals_df,calculated_index)
      colnames(r_crown_vals_df)[ncol(r_crown_vals_df)] <- index_name
    }
  }
  slot(crown_SpatPoldf,"data") <- r_crown_vals_df
  #convert the polygons to points
  #crown_centers <- rgeos::gCentroid(crown_SpatPoldf, byid=TRUE)

  #write away the output
  if (!is.null(shp_outname)){
    raster::shapefile(x = crown_SpatPoldf, filename = shp_outname, overwrite=T)
    cat('Wrote away ', outname, '\n')
  }

  return(crown_SpatPoldf)
}
