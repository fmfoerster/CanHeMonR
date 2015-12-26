#' @title Create A Polygon .shp Of Raster File Extents
#' @description Mine through a directory and make a single polygon .shp
#' where the polygon shows the extent of an image and there's an attribute with the filename. So far only .tif files are considered !
#' @param dirname Name of the directory to mine through
#' @return A .shp file written to dirname and named images_overview.shp
#' @examples \dontrun{
#' create_overview_of_rasters('H:/FISE/forest/CanopyHealthMonitoring/PWN/imagery/PT_Orto/Ortos_DistrCasteloBranco/')
#' }
#' @export
create_overview_of_rasters <- function(dirname){
  #list all the files in the directory
  fnames <- list.files(dirname)
  fnames <- fnames[grep(".tif", fnames)]
  fnames <- fnames[substr(fnames, nchar(fnames)-3,nchar(fnames)) == ".tif"]
  fnames <- file.path(dirname,fnames)
  #get the extent
  rpols <- NULL
  for (fname in fnames){
    r <- raster::raster(fname)
    rpol <- as(raster::extent(r), 'SpatialPolygons')
    raster::projection(rpol) <- raster::projection(r)
    rpol <- sp::spChFIDs(rpol, basename(fname))

    if (is.null(rpols)){
      rpols <- rpol
      #establish the projection
      baseproj <- raster::projection(r)
    }else{
      rpols <- maptools::spRbind(rpols,rpol)
      #Check that this raster has the same extent as the first one
      if (raster::projection(r) != baseproj){
        cat('Error in create_overview_of_rasters.
            Rasters in this directory differ in projection.\n
            You will need to reproject before you can join the extents in a single polygon file\n')

      }
    }
  }

  attribs <- data.frame(fname = basename(fnames),fullname = fnames,row.names = basename(fnames))
  rpols.df <- sp::SpatialPolygonsDataFrame(Sr=rpols, data=attribs)
  raster::shapefile(rpols.df,filename = file.path(dirname, 'images_overview.shp'), overwrite = T)
  cat('Wrote away images_overview.shp in ',dirname,'\n')
  cat('The overview .shp contains ',length(rpols.df),' polygons\n')
  return()
}


