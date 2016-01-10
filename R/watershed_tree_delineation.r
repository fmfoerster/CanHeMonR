#' @title Watershed-Based Tree Detection
#' @description Detect potential crown using watershed based segmentation of an NDVI image
#' @param image_fname Filename of the image to run the tree detection on
#' @param extent Raster extent object. Default is NA, which means the entire image will be analyzed.
#' @param index_name Character (for a spectral index) or numeric (for an individual wavelength).
#' The spectral index or wavelength on which to perform the watershed segmentation. If a wavelength is provided, it should be expressed in nm.
#' Default is "NDVI".
#' @param bg_mask List. Uneven positions should provide index names or wavelength numbers. Even positions in the list
#' should provide the cut-off value of the spectral index, below which, pixels are considered background.
#' Default is list('NDVI',0,1).
#' @param neighbour_radius Numeric. In meters the radius to be considered for the detection of neigbouring objects in m.
#' Higher values causes greater smoothing in the detection of objects, and then trees. Default is 2 m.
#' @param watershed_tolerance . The tolerance setting in the EBImage::watershed operation. Lower values cause greater separation.
#' Defualt is 0.08, designed for NDVI.
#' @param plott Logical. Do you want to plot results? Default is False.
#' @param max_npix Integer. If set, image_name will be tiled to have fewer than max_npix pixels and processed
#' in tiles. Tiling, and setting max_npix needs to be set if parallel = T. Default is Inf.
#' @param parallel Logical. Would you like the tiles to be processed in parallel?  Default is False.
#' @param nWorkers If running the code in parallel, how many workers should be used? Default is 4.
#' @param rough_crowns_shp_fname Character. Filename for the output polygon shapefile. The default is to not write the output away.
#' @import foreach
#' @import maptools
#' @import rgeos
#' @note the raster::rasterToPolygons step is slow - look for ways to speed it up!
#' @return A SpatialPolygons object of suspected trees crowns.
#' @seealso polygons_to_seeds
#' @examples
#' \dontrun{
#' watershed_tree_delineation(
#'   image_fname = 'E:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq',
#'   extent = raster::extent(c(747400, 747490, 4463900, 4463990)),
#'   #extent = raster::extent(c(747400, 747420, 4463900, 4463920)),
#'   index_name = 'NDVI',
#'   bg_mask = list('NDVI',0.35,800, 1200),
#'   plott = T,
#'   neighbour_radius = 1,
#'   watershed_tolerance = 0.08,
#'   rough_crowns_shp_fname = 'C:/Users/pieterbeck/Documents/temp/test_treedetection_pols.shp')
#' }
#' @export
watershed_tree_detection <- function(image_fname, extent, index_name = "NDVI", bg_mask, neighbour_radius = 2, watershed_tolerance = .008,
                                     rough_crowns_shp_fname = '', plott = F, max_npix = Inf, parallel = F, nWorkers = 2){

  # Installing the EBImage package
  # source("http://bioconductor.org/biocLite.R")
  # biocLite()
  # biocLite("EBImage")

  # read in the image, and crop if requested
  input_image <- raster::brick(image_fname)
  if (!is.null(extent)){
    input_image <- raster::crop(input_image, extent)
  }

  # calculate the requested spectral index, or extract the requested wavelength
  index_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = index_name)
  if (plott){raster::plot(index_image, axes = F)}

  # filter out background areas, working through the list that is bg_mask
  for (i in seq(1,length(bg_mask),by=2)){
    mask_image <- CanHeMonR::get_index_or_wavelength_from_brick(br = input_image, index_name_or_wavelength = bg_mask[[i]])
    #()
    index_image[mask_image < bg_mask[[i+1]]] <- NA# -1
  }

  # smooth the spectral index image
  index_image <- raster::focal(index_image, w = matrix(1/9,nrow=3,ncol=3))

  # set up for tiled prcocessing if requested
  if (is.finite(max_npix)){
    tile_extents <- CanHeMonR::tile_raster_extent(r = index_image, max_pixs = max_npix)
  }else{
    tile_extents <- list(raster::extent(index_image))
  }

  #set up the cluster for parallel processing if requested
    try(parallel::stopCluster(cl), silent=T)
    # TO DO add a line that avoids allocating more workers than you have cores
  if (parallel){
    cl <- parallel::makeCluster(nWorkers)
    doParallel::registerDoParallel(cl)
  }

  #choose the appropriate operator for the foreach loop
  require(doParallel) #doesn't appear to work without
  `%op%` <- if (parallel) `%dopar%` else `%do%`
  crown_pols <- foreach::foreach(i = 1:length(tile_extents), .combine = maptools::spRbind,  .multicombine = F, .errorhandling='remove') %op% {
    require(rgeos)
    index_image_tile <- raster::crop(x = index_image, y = tile_extents[[i]])

    if(!is.na(raster::maxValue(index_image_tile))){
      # run the watershed on the smoothed NDVI in matrix form
      watershed_smoothness <- ceiling(neighbour_radius / raster::res(index_image_tile)[1])
      w <- EBImage::watershed(raster::as.matrix(index_image_tile), tolerance = watershed_tolerance, ext = watershed_smoothness)

      # create a template watershed image
      wshed <- index_image_tile

      # put the output in a template image
      wshed <- raster::setValues(wshed, w)

      # this ONLY works in parallel if the dissolve/aggregate step is executed seperately
      # and the by parameter  in raster::aggregate is set explicitly rather than by attribute name
      pol <- raster::rasterToPolygons(wshed, dissolve = F, digits = 1, na.rm = T, function(x){x >= 1})
      #browser()
      nn <- names(pol)
      require(sp) #code breaks if you remove this ! - unless you import(sp) in the namespace?
      pol <- raster::aggregate(pol, by = nn)

      pol <- sp::spChFIDs(pol, as.character(paste0(i,'_',1:length(pol))))
    }else{
      deliberate_error
    }
    pol
    #wshed
  }


 if (plott){raster::plot(crown_pols, col = NULL, add = T, border = 'red',lwd = 1)}

  #merge the crowns that got split because of the tiling
  #Where the tile edges intersect crowns, you end up with split crowns.
  #First, make a spatial lines object of the tile borders.
  tile_borders <- lapply(tile_extents, function(x){as(x, 'SpatialPolygons')})
  tile_borders <-foreach::foreach(i = 1:length(tile_extents), .combine = maptools::spRbind,.multicombine = F, .errorhandling='stop' ) %do% {
    x <- as(tile_extents[[i]],"SpatialLines")
    x <- sp::spChFIDs(x, as.character(i))
    x
  }

  tile_borders <- rgeos::gDifference(tile_borders,as(raster::extent(index_image),"SpatialLines"))
  raster::projection(tile_borders) <- raster::projection(crown_pols)
  #debug(merge_pols_intersected_by_lines)
  crown_pols2 <- merge_pols_intersected_by_lines(spat_pols = crown_pols, spat_lines = tile_borders)
    #this returns only the line segments that represent crown splits.


  #now intersect all the crowns with each of these segments.
  #this will each time identify the crowns that are split by this segment
  #flag these pairs for a merger operation.


  # somehow the pixels with NDVI set to -1 still make it into the polygons
  #pol <- pol[unlist(lapply(raster::extract(x = index_image, y = pol),function(x){mean(x, na.rm = T)} > 0)), ]
  if (plott){raster::plot(index_image, axes = F)}


  #save the polygons
  if (rough_crowns_shp_fname != ''){
    raster::shapefile(crown_pols2, filename = rough_crowns_shp_fname, overwrite = T)
  }

  return(crown_pols2)
}