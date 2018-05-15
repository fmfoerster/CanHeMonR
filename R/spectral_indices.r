#' @title Get Spectral Index or Chosen Wavelength From Brick
#' @description Extract from a rasterband a particular wavelength, or spectral index
#' @param br A multispectral raster brick. Layers should be named using the csic convention
#' @param index_name_or_wavelength Character or numeric. The name of a spectral index or the spectral wavelength,
#' in nanometers, to be extracted.
#' @return A raster layer
#' @export
get_index_or_wavelength_from_brick <- function(br, index_name_or_wavelength){
  require(CanHeMonR)
  if (is.character(index_name_or_wavelength)){
    index_image <- do.call(index_name_or_wavelength, list(df = br))
  }
  if (is.numeric(index_name_or_wavelength)){
    index_image <- CanHeMonR::get_band_of_wavelength(br, wavelengths_in_nm = index_name_or_wavelength)
  }
  return(index_image)
}

#' @title Get The Measurements Made In A Chosen Wavelength
#' @description Extract from a dataframe where columns are measurements in individual wavelengths,
#' and named following Quantalab's convention, the data that are closest to a chosen wavelegth.
#' @param spec_df Dataframe or RasterBrick of spectral measurements, with a column per band/spectral wavelength,
#' and columns (layers if spec_df is RasterBrick) named following Quantalab's convention.
#' @param wavelengths_in_nm Numeric vector. The wavelength(s) of interest, in nm. Defualt is NULL
#' @param wavelength_range Numeric vector of length 2. Instead of requesting observations (columns
#' or rasterlayers) for specific wavelengths  (through wavelengths_in_nm ), one can requests all
#' observations falling in a specific wavelength range to be returned in a given range. Default is NULL
#' @param band_txt Character. A piece of text that only occurs in columns/layers that have
#' radiance/reflectance values. Default is "X"
#' @param splitter character. The character that can be used in strpslit to isolate the wavelength value. Default is "[X.]"
#' @param i The position of the wavelength value after the stringsplit operation.
#' @return the column or rasterlayer with measurements in the chosen wavelength
#' @note TO DO: adapt to keep decimals in wavelength value written in layer name
#' @export
get_band_of_wavelength <- function(spec_df, wavelengths_in_nm = NULL, wavelength_range = NULL,
                                   band_txt = 'X', splitter="[X.]", i = 2){
  find_closest <- function(x,ys){
    closest <- NULL
    for (y in ys){
      closest <- c(closest, which.min(abs(x - y)))
    }
    return(closest)
  }
  if(class(spec_df) == "RasterLayer"){
    my_dat <- spec_df
  }else{
    #check inputs
    if (is.null(wavelengths_in_nm) + is.null(wavelength_range) != 1){
      cat('Error in get_band_of_wavelength. Please set the individual wavelengths you want OR the range')
    }
    #get available wavelengths for df or matrix input
    if((class(spec_df) == "data.frame")|(class(spec_df) == "matrix")){
      #keep only those attributes/columns that appear to be spectral measurements
      spec_df <- spec_df[,grep(band_txt, colnames(spec_df))]

      wavelengths <- as.numeric(unlist(lapply((strsplit(colnames(spec_df),split=splitter)),function(x){x[i]})))
    }
    #get available wavelengths from RasterBrick input
    if(class(spec_df) == "RasterBrick"){
      wavelengths <- as.numeric(unlist(lapply((strsplit(names(spec_df),split=splitter)),function(x){x[i]})))
    }

    #which columns/layers represent the requested wavelengths?
    #if individual wavelengths were requested
    if (!is.null(wavelengths_in_nm)){
      my_columns <- find_closest(wavelengths, wavelengths_in_nm)
    }
    #if a range of wavelengths was requested
    if (!is.null(wavelength_range)){
      #my_columns <- find_closest(wavelengths, wavelength_range[1]) : find_closest(wavelengths, wavelength_range[2])
      my_columns <- which((wavelengths > wavelength_range[1]) & ( wavelengths < wavelength_range[2]))
    }

    #extract the columns/layers representing the wavelenths
    if((class(spec_df) == "data.frame")|(class(spec_df) == "matrix")){
      my_dat <- spec_df[, my_columns]
    }
    if(class(spec_df) == "RasterBrick"){
      my_dat <- raster::subset(x = spec_df, subset = my_columns)
    }
  }
  return(my_dat)
}

#' @title Normalized Difference Vegetation Index simulating DMC sensor (DMC_NDVI)
#' @description Calculate DMC_NDVI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @references https://calval.cr.usgs.gov/JACIE_files/JACIE10/Presentations/ThurAM/Pagnutti_Mary_JACIE_2010_DigitalCameraCal.pdf
#' @return A vector or rasterlayer with the value of the index
#' @note TO DO: the averaging across the broad band could be weighted using spectral response curve
#' @export
DMC_NDVI <- function(df, outp_fname = NULL, ...){
  DMC_R620 <- get_band_of_wavelength(df, wavelength_range = c(600, 676), ...)
  DMC_R725 <- get_band_of_wavelength(df, wavelength_range = c(695, 831), ...)

  if((class(DMC_R620) == "data.frame")|(class(DMC_R620) == "matrix")){
    DMC_R620 <- rowMeans(DMC_R620, na.rm = T)
    DMC_R725 <- rowMeans(DMC_R725, na.rm = T)
  }
  if(class(DMC_R620) == "RasterBrick"){
    DMC_R620 <- raster::calc(DMC_R620, mean)
    DMC_R725 <- raster::calc(DMC_R725, mean)
  }

  outp <- (DMC_R725 - DMC_R620) / (DMC_R725 + DMC_R620)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Difference Vegetation Index simulating ADS100 sensor (ADS100_NDVI)
#' @description Calculate ADS100_NDVI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @references http://w3.leica-geosystems.com/downloads123/zz/airborne/ADS100/brochures-datasheet/Leica_ADS100_DS_en.pdf
#' @return A vector or rasterlayer with the value of the index
#' @note TO DO: the averaging across the broad band could be weighted using spectral response curve
#' @export
ADS100_NDVI <- function(df, outp_fname = NULL, ...){
  ADS100_R <- get_band_of_wavelength(df, wavelength_range = c(619, 651), ...)
  ADS100_NIR <- get_band_of_wavelength(df, wavelength_range = c(808, 882), ...)

  if((class(ADS100_R) == "data.frame")|(class(ADS100_R) == "matrix")){
    ADS100_R <- rowMeans(ADS100_R, na.rm = T)
    ADS100_NIR <- rowMeans(ADS100_NIR, na.rm = T)
  }
  if(class(ADS100_R) == "RasterBrick"){
    ADS100_R <- raster::calc(ADS100_R, mean)
    ADS100_NIR <- raster::calc(ADS100_NIR, mean)
  }

  outp <- (ADS100_NIR - ADS100_R) / (ADS100_NIR + ADS100_R)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Green_over_Red Index simulating DMC sensor (DMC_G2R)
#' @description Calculate DMC_G2R
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector or rasterlayer with the value of the index
#' @references https://calval.cr.usgs.gov/JACIE_files/JACIE10/Presentations/ThurAM/Pagnutti_Mary_JACIE_2010_DigitalCameraCal.pdf
#' @note TO DO: the averaging across the broad band could be weighted using spectral response curve
#' @export
DMC_G2R <- function(df, outp_fname = NULL, ...){
  DMC_R <- get_band_of_wavelength(df, wavelength_range = c(600, 676), ...)
  DMC_G <- get_band_of_wavelength(df, wavelength_range = c(514, 600), ...)

  if((class(DMC_R) == "data.frame")|(class(DMC_R) == "matrix")){
    DMC_R <- rowMeans(DMC_R, na.rm = T)
    DMC_G <- rowMeans(DMC_G, na.rm = T)
  }
  if(class(DMC_R) == "RasterBrick"){
    DMC_R <- raster::calc(DMC_R, mean)
    DMC_G <- raster::calc(DMC_G, mean)
  }

  outp <- (DMC_G) / (DMC_R)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Green_over_Red Index simulating ADS100 sensor (ADS100_G2R)
#' @description Calculate ADS100_G2R
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector or rasterlayer with the value of the index
#' @references http://w3.leica-geosystems.com/downloads123/zz/airborne/ADS100/brochures-datasheet/Leica_ADS100_DS_en.pdf
#' @note TO DO: the averaging across the broad band could be weighted using spectral response curve
#' @export
ADS100_G2R <- function(df, outp_fname = NULL, ...){
  ADS100_R <- get_band_of_wavelength(df, wavelength_range = c(619, 651), ...)
  ADS100_G <- get_band_of_wavelength(df, wavelength_range = c(525, 585), ...)

  if((class(ADS100_R) == "data.frame")|(class(ADS100_R) == "matrix")){
    ADS100_R <- rowMeans(ADS100_R, na.rm = T)
    ADS100_G <- rowMeans(ADS100_G, na.rm = T)
  }
  if(class(ADS100_R) == "RasterBrick"){
    ADS100_R <- raster::calc(ADS100_R, mean)
    ADS100_G <- raster::calc(ADS100_G, mean)
  }

  outp <- (ADS100_G) / (ADS100_R)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Difference Vegetation Index (NDVI)
#' @description Calculate NDVI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Rouse et al. 1974
#' @examples
#' \dontrun{
#' #Calculate the NDVI from a 6 band multispectral image
#' my_ms_data <- raster::brick('H:/FISE/forest/CanopyHealthMonitoring/PWN/flights_final/150727_mca/150727_mca.bsq')
#' my_NDVI <- CanHeMonR::NDVI(my_ms_data)
#' }
#' @export
NDVI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)

  outp <- (R800 - R670) / (R800 + R670 )

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Renormalized Difference Vegetation Index
#' @description Calculate Renormalized Difference Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength, and
#' columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename
#' to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Rouse et al. 1974
#' @export
RDVI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (R800 - R670) / sqrt(R800 + R670 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Simple Ratio
#' @description Calculate the simple ratio vegetation index
#' @param df A data frame where columns represent measurements in a single wavelength, and
#' columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename
#' to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Jordat 1969
#' @export
SR <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)

  outp <- R800 / R670

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Modified Simple Ratio
#' @description Calculate Modified Simple Ratio vegetation index
#' @param df A data frame where columns represent measurements in a single wavelength, and columns
#' are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Chen 1996
#' @export
mod_SR <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (R800 / R670 - 1) / (sqrt(R800 / R670 ) + 1)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Optimised Soil-Adjusted Vegetation Index
#' @description Calculate Renormalized Difference Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @note mostly structure-sensitive
#' @references Rondeaux et al. 1996
#' @export
OSAVI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (1 + 0.16) (R800 - R670) / (R800 + R670 + 0.16)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Modified Soil-Adjusted Vegetation Index
#' @description Calculate Modified Soil-Adjusted Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Qi et al. 1994
#' @export
MSAVI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (2 * R800 + 1 - sqrt((2 * R800 + 1)^2 - 8*(R800 - R670))) / 2
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Triangular Vegetation Index
#' @description Calculate Triangular Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Broge and Leblanc 2000
#' @export
TVI <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- 0.5 * (120 * (R750 - R550) - 200 * (R670 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Triangular Vegetation Index 1
#' @description Calculate Modified Triangular Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MTVI1 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- 1.2 * (1.2 * (R800 - R550) - 2.5 * (R670 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Triangular Vegetation Index 2
#' @description Calculate Modified Triangular Vegetation Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MTVI2 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (1.2 * (1.2 * (R800 - R550) - 2.5 * (R670 - R550)))/
    sqrt((2 * R800 + 1)^2 - (6 * R800 - 5*sqrt(R670) ) - 0.5)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}



#' @title Modified Chlorophyll Absorption Ratio Index 1
#' @description Calculate Modified Chlorophyll Absorption Ratio Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MCARI1 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- 1.2 * (2.5*(R800 - R670) - 1.3 * (R800 - R550))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Modified Chlorophyll Absorption Ratio Index 2
#' @description Calculate Modified Chlorophyll Absorption Ratio Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Haboudane et al 2004
#' @export
MCARI2 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- 1.2 * (2.5*(R800 - R670) - 1.3 * (R800 - R550))/
    sqrt((2 * R800 + 1)^2 - (6 * R800 - 5*sqrt(R670) ) - 0.5)


  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Enhanced Vegetation Index
#' @description Calculate Enhanced Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Liu and Huete 1995
#' @export
EVI <- function(df, outp_fname = NULL, ...){
  R400 <- get_band_of_wavelength(df, wavelengths_in_nm = 400, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- 2.5 * (R800 - R670) / (R800 + 6*R670 - 7.5*R400 + 1)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Lichtenthaler Index 1
#' @description Calculate Lichtenthaler Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Lichtenthaler et al 1996
#' @export
LIC1 <- function(df, outp_fname = NULL, ...){
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (R800 - R680) / (R800 + R680)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 1
#' @description Calculate Vogelmann Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG1 <- function(df, outp_fname = NULL, ...){
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)
  R740 <- get_band_of_wavelength(df, wavelengths_in_nm = 740, ...)
  outp <- R740 / R720
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 2
#' @description Calculate Vogelmann Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG2 <- function(df, outp_fname = NULL, ...){
  R715 <- get_band_of_wavelength(df, wavelengths_in_nm = 715, ...)
  R726 <- get_band_of_wavelength(df, wavelengths_in_nm = 726, ...)
  R734 <- get_band_of_wavelength(df, wavelengths_in_nm = 734, ...)
  R747 <- get_band_of_wavelength(df, wavelengths_in_nm = 747, ...)
  outp <- (R734 - R747 ) / (R715 + R726)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Vogelmann Index 3
#' @description Calculate Vogelmann Index 3
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Vogelmann et al 1993
#' @export
VOG3 <- function(df, outp_fname = NULL, ...){
  R715 <- get_band_of_wavelength(df, wavelengths_in_nm = 715, ...)
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)
  R734 <- get_band_of_wavelength(df, wavelengths_in_nm = 734, ...)
  R747 <- get_band_of_wavelength(df, wavelengths_in_nm = 747, ...)
  outp <- (R734 - R747 ) / (R715 + R720)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Gitelson & Merzlyak Index 1
#' @description CalculateGitelson & Merzlyak index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Gitelson & Merzlyak 1997
#' @export
GM1 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- R750 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Gitelson & Merzlyak Index 2
#' @description Calculate Gitelson & Merzlyak index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Gitelson & Merzlyak 1997
#' @export
GM2 <- function(df, outp_fname = NULL, ...){
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- R750 / R700
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Transformed Chlorophyll Absorption in Reflectance Index
#' @description Calculate Transformed Chlorophyll Absorption in Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @note TCARI is mostly pigment related but, at least in crops, is also sensitive to structure. Hence
#' TCARI/OSAVI estimates chlorohyll content better in eg mais crops over a range of LAI values
#' @references Haboudane et al. 2002
#' @export
TCARI <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  outp <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700/R670))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Transformed Chlorophyll Absorption in Reflectance Index Over Optimised Soil-Adjusted Vegetation Index
#' @description Calculate Transformed Chlorophyll Absorption in Reflectance Index Over Optimised Soil-Adjusted Vegetation Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Haboudane et al. 2002
#' @note TCARI is mostly pigment related but, at least in crops, is also sensitive to structure. Hence
#' TCARI/OSAVI estimates chlorohyll content better in eg mais crops over a range of LAI values
#' @export
TCARI_over_OSAVI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- TCARI(df, ...) / ((1 + 0.16) * (R800 - R670) / (R800 + R670 + 0.16))
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Chlorophyll Index Red Edge
#' @description Calculate Chlorophyll Index Red EdgeChlorophyll Index Red Edge
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CI <- function(df, outp_fname = NULL, ...){
  R710 <- get_band_of_wavelength(df, wavelengths_in_nm = 710, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- R750 / R710
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Simple Ratio Pigment Index
#' @description Calculate Simple Ratio Pigment Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
SRPI <- function(df, outp_fname = NULL, ...){
  R430 <- get_band_of_wavelength(df, wavelengths_in_nm = 430, ...)
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  outp <- R430 / R680
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Phaeophytinization Index
#' @description Calculate Normalized Phaeophytinization Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
NPQI <- function(df, outp_fname = NULL, ...){
  R415 <- get_band_of_wavelength(df, wavelengths_in_nm = 415, ...)
  R735 <- get_band_of_wavelength(df, wavelengths_in_nm = 735, ...)
  outp <- ( R415 - R735 ) / ( R415 + R735)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Normalized Pigments Index
#' @description Calculate Normalized Pigments Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
NPCI <- function(df, outp_fname = NULL, ...){
  R430 <- get_band_of_wavelength(df, wavelengths_in_nm = 430, ...)
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  outp <- ( R680 - R430 ) / ( R680 + R430)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carter Index 1
#' @description Calculate Carter Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CTRI1 <- function(df, outp_fname = NULL, ...){
  R420 <- get_band_of_wavelength(df, wavelengths_in_nm = 420, ...)
  R695 <- get_band_of_wavelength(df, wavelengths_in_nm = 695, ...)
  outp <- R695 / R420
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carter Index 2
#' @description Calculate Carter Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CTRI2 <- function(df, outp_fname = NULL, ...){
  R695 <- get_band_of_wavelength(df, wavelengths_in_nm = 695, ...)
  R760 <- get_band_of_wavelength(df, wavelengths_in_nm = 760, ...)
  outp <- R695 / R760
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Reflectance Band Ratio
#' @description Calculate Reflectance Band Ratio
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Datt et al. 1998
#' @export
datt_CabCx_c <- function(df, outp_fname = NULL, ...){
  R672 <- get_band_of_wavelength(df, wavelengths_in_nm = 672, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R708 <- get_band_of_wavelength(df, wavelengths_in_nm = 708, ...)
  outp <- R672 / (R550 * 3 * R708)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp )
}

#' @title Reflectance Band Ratio Using NIR
#' @description Calculate Reflectance Band Ratio Using NIR
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Datt et al. 1998
#' @export
datt_NIRCabCx_c <- function(df, outp_fname = NULL, ...){
  R860 <- get_band_of_wavelength(df, wavelengths_in_nm = 860, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R708 <- get_band_of_wavelength(df, wavelengths_in_nm = 708, ...)
  outp <- R860 / (R550 * R708)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Structure-Insensitive Pigment Index
#' @description Calculate Structure-Insensitive Pigment Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
SIPI <- function(df, outp_fname = NULL, ...){
  R445 <- get_band_of_wavelength(df, wavelengths_in_nm = 445, ...)
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- (R800 - R445) / (R800 + R680)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 550 nm
#' @description Calculate Carotenoid Reflectance Index 550 nm
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CRI550 <- function(df, outp_fname = NULL, ...){
  R515 <- get_band_of_wavelength(df, wavelengths_in_nm = 515, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  outp <- (1 / R515) - (1 / R550)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 700 nm
#' @description Calculate Carotenoid Reflectance Index 700 nm
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CRI700 <- function(df, outp_fname = NULL, ...){
  R515 <- get_band_of_wavelength(df, wavelengths_in_nm = 515, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  outp <- (1 / R515) - (1 / R700)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}


#' @title Carotenoid Reflectance Index 550 nm with NIR
#' @description Calculate Carotenoid Reflectance Index 550 nm with NIR
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
RNIR_CRI550 <- function(df, outp_fname = NULL, ...){
  R515 <- get_band_of_wavelength(df, wavelengths_in_nm = 515, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R770 <- get_band_of_wavelength(df, wavelengths_in_nm = 770, ...)
  outp <- (1 / R515) - (1 / R550) * R770
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Carotenoid Reflectance Index 700 nm with NIR
#' @description Calculate Carotenoid Reflectance Index 700 nm with NIR
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
RNIR_CRI700 <- function(df, outp_fname = NULL, ...){
  R515 <- get_band_of_wavelength(df, wavelengths_in_nm = 515, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R770 <- get_band_of_wavelength(df, wavelengths_in_nm = 770, ...)
  outp <- (1 / R515) - (1 / R700) * R770
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Plant Senescing Reflectance Index
#' @description Calculate Plant Senescing Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PSRI <- function(df, outp_fname = NULL, ...){
  R500 <- get_band_of_wavelength(df, wavelengths_in_nm = 500, ...)
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- (R680 - R500) / R750
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return( outp)
}

#' @title Lichtenhaler Index
#' @description Calculate Lichtenhaler index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
LIC3 <- function(df, outp_fname = NULL, ...){
  R440 <- get_band_of_wavelength(df, wavelengths_in_nm = 440, ...)
  R740 <- get_band_of_wavelength(df, wavelengths_in_nm = 740, ...)
  outp <- R440/R740
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Ratio Analysis Of Reflectance Spectra
#' @description Calculate Ratio Analysis Of Reflectance Spectra
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
RARS <- function(df, outp_fname = NULL, ...){
  R513 <- get_band_of_wavelength(df, wavelengths_in_nm = 513, ...)
  R746 <- get_band_of_wavelength(df, wavelengths_in_nm = 746, ...)
  outp <- R513/R746
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Chlorophyll a
#' @description Calculate Pigment Specific Simple Ratio Chlorophyll a
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRa <- function(df, outp_fname = NULL, ...){
  R675 <- get_band_of_wavelength(df, wavelengths_in_nm = 675, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- R800 / R675
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Chlorophyll b
#' @description Calculate Pigment Specific Simple Ratio Chlorophyll b
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRb <- function(df, outp_fname = NULL, ...){
  R650 <- get_band_of_wavelength(df, wavelengths_in_nm = 650, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- R800 / R650
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Simple Ratio Carotenoids
#' @description Calculate Pigment Specific Simple Ratio Carotenoids
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PSSRc <- function(df, outp_fname = NULL, ...){
  R500 <- get_band_of_wavelength(df, wavelengths_in_nm = 500, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- R800 / R500
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Pigment Specific Normalized Difference
#' @description Calculate Pigment Specific Normalized Difference
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PSNDc <- function(df, outp_fname = NULL, ...){
  R470 <- get_band_of_wavelength(df, wavelengths_in_nm = 470, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  outp <- ( R800 - R470 ) / ( R800 + R470 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 570
#' @description Calculate Photochemical Reflectance Index 570
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI570 <- function(df, outp_fname = NULL, ...){
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  R570 <- get_band_of_wavelength(df, wavelengths_in_nm = 570, ...)
  outp <- ( R570 - R531 ) / ( R570 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 515
#' @description Calculate Photochemical Reflectance Index 515
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI515 <- function(df, outp_fname = NULL, ...){
  R515 <- get_band_of_wavelength(df, wavelengths_in_nm = 515, ...)
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  outp <- ( R515 - R531 ) / ( R515 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 512
#' @description Calculate Photochemical Reflectance Index 512
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI512 <- function(df, outp_fname = NULL, ...){
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  R512 <- get_band_of_wavelength(df, wavelengths_in_nm = 512, ...)
  outp <- ( R512 - R531 ) / ( R512 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 600
#' @description Calculate Photochemical Reflectance Index 600
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI600 <- function(df, outp_fname = NULL, ...){
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  R600 <- get_band_of_wavelength(df, wavelengths_in_nm = 600, ...)
  outp <- ( R600 - R531 ) / ( R600 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 670
#' @description Calculate Photochemical Reflectance Index 670
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI670 <- function(df, outp_fname = NULL, ...){
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  outp <- ( R670 - R531 ) / ( R670 + R531 )
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Photochemical Reflectance Index 670 570
#' @description Calculate Photochemical Reflectance Index 670 570
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI670_570 <- function(df, outp_fname = NULL, ...){
  R531 <- get_band_of_wavelength(df, wavelengths_in_nm = 531, ...)
  R570 <- get_band_of_wavelength(df, wavelengths_in_nm = 570, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  outp <- ( R570 - R531 - R670) / ( R570 + R531 + R670)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Normalized Photochemical Reflectance Index
#' @description Calculate Normalized Photochemical Reflectance Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRIn <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  outp <- ( CanHeMonR::PRI570(df, ...)) / abs( CanHeMonR::RDVI(df, ...) * R700/R670)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Carotenoid/chlorophyll Ratio Index
#' @description Calculate Carotenoid/chlorophyll Ratio Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
PRI_CI <- function(df, outp_fname = NULL, ...){
  R500 <- get_band_of_wavelength(df, wavelengths_in_nm = 500, ...)
  R680 <- get_band_of_wavelength(df, wavelengths_in_nm = 680, ...)
  R750 <- get_band_of_wavelength(df, wavelengths_in_nm = 750, ...)
  outp <- (R680 - R500) / R750
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Reflectance Curvature Index
#' @description Calculate Reflectance Curvature Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
CUR <- function(df, outp_fname = NULL, ...){
  R675 <- get_band_of_wavelength(df, wavelengths_in_nm = 675, ...)
  R683 <- get_band_of_wavelength(df, wavelengths_in_nm = 683, ...)
  R690 <- get_band_of_wavelength(df, wavelengths_in_nm = 690, ...)
  outp <- (R675 - R690) / (R683 ^ 2)
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Redness Index
#' @description Calculate Redness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
Redness <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  outp <- R700 / R670
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Greenness Index
#' @description Calculate Greenness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
Greenness <- function(df, outp_fname = NULL, ...){
  R570 <- get_band_of_wavelength(df, wavelengths_in_nm = 570, ...)
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  outp <- R570 / R670
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Index
#' @description Calculate Greenness Index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
Blue_index <- function(df, outp_fname = NULL, ...){
  R450 <- get_band_of_wavelength(df, wavelengths_in_nm = 450, ...)
  R490 <- get_band_of_wavelength(df, wavelengths_in_nm = 490, ...)
  outp <- R450 / R490
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Green Index 1
#' @description Calculate Blue Green Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
BGI1 <- function(df, outp_fname = NULL, ...){
  R400 <- get_band_of_wavelength(df, wavelengths_in_nm = 400, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  outp <- R400 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Green Index 2
#' @description Calculate Blue Green Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
BGI2 <- function(df, outp_fname = NULL, ...){
  R450 <- get_band_of_wavelength(df, wavelengths_in_nm = 450, ...)
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  outp <- R450 / R550
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Red Index 1
#' @description Calculate Blue Red Index 1
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
BRI1 <- function(df, outp_fname = NULL, ...){
  R400 <- get_band_of_wavelength(df, wavelengths_in_nm = 400, ...)
  R690 <- get_band_of_wavelength(df, wavelengths_in_nm = 690, ...)
  outp <- R400 / R690
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Blue Red Index 2
#' @description Calculate Blue Red Index 2
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
BRI2 <- function(df, outp_fname = NULL, ...){
  R450 <- get_band_of_wavelength(df, wavelengths_in_nm = 450, ...)
  R690 <- get_band_of_wavelength(df, wavelengths_in_nm = 690, ...)
  outp <- R450 / R690
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Lichtenhaler (Blue Over Red) Index 2
#' @description Calculate Lichtenhaler red over blue index
#' @param df A data frame where columns represent measurements in a single wavelength,
#' and columns are named following Quantalab's conventions
#' @param outp_fname In case the input is raster data, this is the optional output filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references x
#' @export
LIC2 <- function(df, outp_fname = NULL, ...){
  R <- get_band_of_wavelength(df, wavelengths_in_nm = 690, ...)
  B <- get_band_of_wavelength(df, wavelengths_in_nm = 440, ...)
  outp <- B/R
  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

##################################################################
#
# Extended Indices mainly for SWIR
#
# and Nitrogen indices
#
##################################################################

#' @title Normalized Difference Red Edge (NDRE)
#' @description Calculate NDRE
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Barnes et al. 2000
#' @examples
#' @export
NDRE <- function(df, outp_fname = NULL, ...){
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)
  R790 <- get_band_of_wavelength(df, wavelengths_in_nm = 790, ...)

  outp <- (R790 - R720) / (R790 + R720 )

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Doube-peak Canopy N Index (DCNI)
#' @description Calculate DCNI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Chen et al. 2010
#' @examples
#' @export
DCNI <- function(df, outp_fname = NULL, ...){
  R670 <- get_band_of_wavelength(df, wavelengths_in_nm = 670, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)

  outp <- (R720 - R700)*(R700 - R670) / ((R720 - R670 ) + 0.3)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Modified Chlorophyll Absorption in Reflectance Index at 1510 nm (MCARI_1510)
#' @description Calculate MCARI_1510
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Herrmann et al. 2010
#' @examples
#' @export
MCARI_1510 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R1510 <- get_band_of_wavelength(df, wavelengths_in_nm = 1510, ...)

  outp <- ((R700 - R1510) - 0.2 * (R700 - R550)) * (R700 / R1510 )

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Transformed Chlorophyll Absorption in Reflectance Index at 1510 nm (TCARI_1510)
#' @description Calculate TCARI_1510
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Herrmann et al. 2010
#' @examples
#' @export
TCARI_1510 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R1510 <- get_band_of_wavelength(df, wavelengths_in_nm = 1510, ...)

  outp <- 3 * ((R700 - R1510) - 0.2 * (R700 - R550)) * (R700 / R1510 )

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title The TCARI and the Soil-Adjusted Vegetation Index at 1510 nm (TCARI_1510/OSAVI_1510)
#' @description Calculate TCARI_1510/OSAVI_1510
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Herrmann et al. 2010
#' @examples
#' @export
TOCARI_1510 <- function(df, outp_fname = NULL, ...){
  R550 <- get_band_of_wavelength(df, wavelengths_in_nm = 550, ...)
  R700 <- get_band_of_wavelength(df, wavelengths_in_nm = 700, ...)
  R800 <- get_band_of_wavelength(df, wavelengths_in_nm = 800, ...)
  R1510 <- get_band_of_wavelength(df, wavelengths_in_nm = 1510, ...)

  outp <- (3 * ((R700 - R1510) - 0.2 * (R700 - R550)) * (R700 / R1510 )) / ((1 + 0.16) * (R800 - R1510) / (R800 + R1510 + 0.16))

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Ratio Spectral Index (RSI)
#' @description Calculate RSI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Zhu et al. 2008
#' @examples
#' @export
RSI_990 <- function(df, outp_fname = NULL, ...){
  R990 <- get_band_of_wavelength(df, wavelengths_in_nm = 990, ...)
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)

  outp <- R990 / R720

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Combined Index (CI1)
#' @description Calculate CI1 (since the bands with 735nm and 736nm are used only a data-set with very narrow bandwidths can be used)
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Bao et al. 2013
#' @examples
#' @export
CI1 <- function(df, outp_fname = NULL, ...){
  R990 <- get_band_of_wavelength(df, wavelengths_in_nm = 990, ...)
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)
  R736 <- get_band_of_wavelength(df, wavelengths_in_nm = 736, ...)
  R735 <- get_band_of_wavelength(df, wavelengths_in_nm = 735, ...)

  outp <- (R736 - R735) * (R990 / R720)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Combined Index (CI2)
#' @description Calculate CI2 (since the bands with 735nm and 736nm are used only a data-set with very narrow bandwidths can be used)
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Bao et al. 2013
#' @examples
#' @export
CI2 <- function(df, outp_fname = NULL, ...){
  R900 <- get_band_of_wavelength(df, wavelengths_in_nm = 900, ...)
  R720 <- get_band_of_wavelength(df, wavelengths_in_nm = 720, ...)
  R736 <- get_band_of_wavelength(df, wavelengths_in_nm = 736, ...)
  R735 <- get_band_of_wavelength(df, wavelengths_in_nm = 735, ...)

  outp <- (R736 - R735) * (R900 / R720)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}


#' @title Normalized Difference Nitrogen Index (NDNI)
#' @description Calculate NDNI
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Serrano et al. 2002
#' @examples
#' @export
NDNI <- function(df, outp_fname = NULL, ...){
  R1510 <- get_band_of_wavelength(df, wavelengths_in_nm = 1510, ...)
  R1680 <- get_band_of_wavelength(df, wavelengths_in_nm = 1680, ...)


  outp <- (log(1/R1510)-log(1/R1680)) / (log(1/R1510)+log(1/R1680))

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title GnyLi
#' @description Calculate GnyLi
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Gnyp et al. 2014
#' @examples
#' @export
GnyLi <- function(df, outp_fname = NULL, ...){
  R900 <- get_band_of_wavelength(df, wavelengths_in_nm = 900, ...)
  R955 <- get_band_of_wavelength(df, wavelengths_in_nm = 955, ...)
  R1050 <- get_band_of_wavelength(df, wavelengths_in_nm = 1050, ...)
  R1220 <- get_band_of_wavelength(df, wavelengths_in_nm = 1220, ...)

  outp <- ((R900 * R1050) * (R955 * R1220)) / ((R900 * R1050) + (R955 * R1220))

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title N1645,1715
#' @description Calculate N1645,1715
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Pimstein et al. 2011
#' @examples
#' @export
N1645 <- function(df, outp_fname = NULL, ...){
  R1645 <- get_band_of_wavelength(df, wavelengths_in_nm = 1645, ...)
  R1715 <- get_band_of_wavelength(df, wavelengths_in_nm = 1715, ...)


  outp <- (R1645 - R1715) / (R1645 + R1715)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title N870,1450
#' @description Calculate N870,1450
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Pimstein et al. 2011
#' @examples
#' @export
N870 <- function(df, outp_fname = NULL, ...){
  R870 <- get_band_of_wavelength(df, wavelengths_in_nm = 870, ...)
  R1450 <- get_band_of_wavelength(df, wavelengths_in_nm = 1450, ...)


  outp <- (R870 - R1450) / (R870 + R1450)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title S1080,660
#' @description Calculate S1080,660
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Mahajan et al. 2014
#' @examples
#' @export
S1080 <- function(df, outp_fname = NULL, ...){
  R660 <- get_band_of_wavelength(df, wavelengths_in_nm = 660, ...)
  R1080 <- get_band_of_wavelength(df, wavelengths_in_nm = 1080, ...)


  outp <- (R1080 - R660) / (R1080 + R660)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title S1260,660
#' @description Calculate S1260,660
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references Mahajan et al. 2014
#' @examples
#' @export
S1260 <- function(df, outp_fname = NULL, ...){
  R660 <- get_band_of_wavelength(df, wavelengths_in_nm = 660, ...)
  R1260 <- get_band_of_wavelength(df, wavelengths_in_nm = 1260, ...)


  outp <- (R1260 - R660) / (R1260 + R660)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

##################################################################
#
# Extended Indices for exploring
#
# Florescense band (around 760nm)
#
##################################################################

#' @title Ratio before 760 nm
#' @description Calculates the ratio of the fluorescense absorption band and the band before
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references
#' @examples
#' @export
RB760 <- function(df, outp_fname = NULL, ...){
  R760 <- get_band_of_wavelength(df, wavelengths_in_nm = 760, ...)
  R755 <- get_band_of_wavelength(df, wavelengths_in_nm = 755, ...)


  outp <- (R755) / (R760)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Ratio after 760 nm
#' @description Calculates the ratio of the fluorescense absorption band and the band after
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references
#' @examples
#' @export
RA760 <- function(df, outp_fname = NULL, ...){
  R760 <- get_band_of_wavelength(df, wavelengths_in_nm = 760, ...)
  R765 <- get_band_of_wavelength(df, wavelengths_in_nm = 765, ...)


  outp <- (R765) / (R760)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}

#' @title Triangle Ratio around 760 nm
#' @description Calculates area of the fluorescense absorption band and the band before and after
#' @param df A data frame (or raster::brick) where columns (layers) represent reflectance measurements in a single wavelength, and
#' columns (layers) are named following Quantalab's conventions.
#' @param outp_fname In case the input is raster data, this is the optional output
#' filename to write the result to
#' @param ... Arguments to be passed to get_band_of_wavelength, particularly band_txt, splitter, and i.
#' @return A vector with the value of the index
#' @references
#' @examples
#' @export
FA760 <- function(df, outp_fname = NULL, ...){
  R755 <- get_band_of_wavelength(df, wavelengths_in_nm = 755, ...)
  R760 <- get_band_of_wavelength(df, wavelengths_in_nm = 760, ...)
  R765 <- get_band_of_wavelength(df, wavelengths_in_nm = 765, ...)


  outp <- ABS((756.6*(R760 - R765)) + (760.3*(R765 - R755)) + (763.9*(R755 - R760)) / 2)

  if ((!is.null(outp_fname)) & (class(outp) == "RasterLayer")){
    raster::writeRaster(outp, filename = outp_fname, overwrite = T)
    cat('Spectral index written away as raster to: ',outp_fname, '\n')
  }
  return(outp)
}
