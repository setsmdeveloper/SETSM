#ifndef _SETSMGEO_H_
#define _SETSMGEO_H_

#include "geotiffio.h"
#include "xtiffio.h"
#include "Typedefine.hpp"

int SetUpTIFFDirectory(TIFF *tif, size_t width, size_t height, double scale, double minX, double maxY, int data_type);
int SetUpGeoKeys(GTIF *gtif, int projection, int zone, int NS_hemisphere);
int WriteGeotiff(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type);
CSize ReadGeotiff_info(char *filename, double *minX, double *maxY, double *grid_size);

#endif

