#ifndef _SETSMGEO_H_
#define _SETSMGEO_H_

#include "geotiffio.h"
#include "xtiffio.h"
#include "Typedefine.hpp"

void SetUpTIFFDirectory(TIFF *tif, size_t width, size_t height, double scale, double minX, double maxY, int data_type);
void SetUpTIFFDirectory_round(TIFF *tif, size_t width, size_t height, double scale, double minX, double maxY, int data_type);
void SetUpGeoKeys(GTIF *gtif, int projection, int zone, int NS_hemisphere);
int WriteGeotiff(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type);
int WriteGeotiff_round(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type);
uint8 ReadGeotiff_bits(char *filename);
CSize ReadGeotiff_info(const char *filename, double *minX, double *maxY, double *grid_size);
CSize ReadGeotiff_info_dxy(char *filename, double *minX, double *maxY, double *grid_size_dx, double *grid_size_dy);

#endif
