#include "setsmgeo.hpp"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <cmath>

#define FLOAT 4
#define UCHAR 1
#define UINT16 12

#define N(a) (sizeof(a) / sizeof(a[0]))
#define GDAL_NODATA 42113
static const TIFFFieldInfo xtiffFieldInfo[] = {
    { GDAL_NODATA, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0, (char*)"GDAL_NODATA" }
};

#define STRIP_SIZE_DEFAULT 8192
#define max(a, b)  (((a) > (b)) ? (a) : (b))


int WriteGeotiff_real(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type, bool doRound)
{
    TIFF *tif;  /* TIFF-level descriptor */
    GTIF *gtif; /* GeoKey-level descriptor */
    size_t bytes = 0; /* Size of data type for computing buffer offset */

    switch (data_type)
    {
        case FLOAT:
            bytes = sizeof(float);
            break;
        case UCHAR:
            bytes = sizeof(unsigned char);
            break;
        case UINT16:
            bytes = sizeof(uint16);
            break;
        default:
            printf("unrecognized data type: %d\n", data_type);
            break;
    }

    if (bytes == 0)
    {
        // Unknown data type
        return -1;
    }

    long int int_mem = (long)(bytes*height*width);
    double tif_mem = (double)(int_mem/1024.0/1024.0/1024.0);
    printf("tif mem %f\n",tif_mem);

    tif = XTIFFOpen(filename, "w8");
    if (!tif)
    {
        printf("WriteGeotiff failed in XTIFFOpen\n");
        return -1;
    }

    gtif = GTIFNew(tif);
    if (!gtif)
    {
        printf("WriteGeotiff failed in GTIFNew\n");
        XTIFFClose(tif);
        return -1;
    }

    SetUpTIFFDirectory(tif, width, height, scale, minX, maxY, data_type);
    SetUpGeoKeys(gtif, projection, zone, NS_hemisphere);

    // scanline_buffer is needed because TIFFWriteScanline modifies the buffer if PREDICTOR != 1.
    void *scanline_buffer = malloc(bytes * width);
    if (!scanline_buffer) {
        printf("WriteGeotiff failed in malloc scanline_buffer\n");
        return -1;
    }

    for (int row=0; row<height; row++)
    {
        memcpy(scanline_buffer, ((char *)buffer) + (bytes * row * width), bytes * width);

        // If requested, round values to nearest 1/128
        if (doRound && data_type == FLOAT) {
            for(size_t i = 0 ; i < width; i++)
            {
                float &val = ((float*)scanline_buffer)[i];
                val = std::round(val*128) / 128;
            }
        }

        if (TIFFWriteScanline(tif, scanline_buffer, row, 0) == -1) // TODO: TIFFWriteScanline may return -1 on failure:
        {
            TIFFError("WriteGeotiff_DEM","failure in WriteScanline on row %d\n", row);
        }
    }

    free(scanline_buffer);
    
    GTIFWriteKeys(gtif);
    GTIFFree(gtif);
    XTIFFClose(tif);
    return 0;
}


int WriteGeotiff(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type)
{
	return WriteGeotiff_real(filename, buffer, width, height, scale, minX, maxY, projection, zone, NS_hemisphere, data_type, false);
}

int WriteGeotiff_round(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type)
{
	return WriteGeotiff_real(filename, buffer, width, height, scale, minX, maxY, projection, zone, NS_hemisphere, data_type, true);
}

uint8 ReadGeotiff_bits(char *filename)
{
    TIFF *tif;
    uint8 bits;
    
    tif = XTIFFOpen(filename, "r");
    if (tif)
    {
        size_t value = 0;
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &value);
        
        if(value == 16)
            bits = 12;
        else
            bits = value;
        
        XTIFFClose(tif);
    }
    return bits;
}

CSize ReadGeotiff_info(const char *filename, double *minX, double *maxY, double *grid_size)
{
    TIFF *tif;
    CSize image_size;

    tif = XTIFFOpen(filename, "r");
    if (tif)
    {
        uint16 count = 0;
        double *data = 0;

        TIFFGetField(tif, TIFFTAG_GEOTIEPOINTS, &count, &data);
        if (minX != NULL) *minX = data[3];
        if (maxY != NULL) *maxY = data[4];

        TIFFGetField(tif, TIFFTAG_GEOPIXELSCALE, &count, &data);
        if (grid_size != NULL) *grid_size = data[0];

        size_t value = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &value);
        image_size.width = value;
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &value);
        image_size.height = value;

        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &value);
        
        if(value == 16)
            value = 12;
       XTIFFClose(tif);
    }
    return image_size;
}

CSize ReadGeotiff_info_dxy(char *filename, double *minX, double *maxY, double *grid_size_dx, double *grid_size_dy)
{
    TIFF *tif;
    CSize image_size;
    
    tif = XTIFFOpen(filename, "r");
    if (tif)
    {
        uint16 count = 0;
        double *data = 0;
        
        TIFFGetField(tif, TIFFTAG_GEOTIEPOINTS, &count, &data);
        if (minX != NULL) *minX = data[3];
        if (maxY != NULL) *maxY = data[4];
        
        TIFFGetField(tif, TIFFTAG_GEOPIXELSCALE, &count, &data);
        if (grid_size_dx != NULL) *grid_size_dx = data[0];
        if (grid_size_dy != NULL) *grid_size_dy = data[1];
        
        size_t value = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &value);
        image_size.width = value;
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &value);
        image_size.height = value;
        
        XTIFFClose(tif);
    }
    return image_size;
}


void SetUpTIFFDirectory(TIFF *tif, size_t width, size_t height, double scale, double minX, double maxY, int data_type)
{
    double tiepoints[6] = {0};
    double pixscale[3] = {0};
    tiepoints[3] = minX;
    tiepoints[4] = maxY;
    pixscale[0] = scale;
    pixscale[1] = scale;

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_GEOTIEPOINTS, 6,tiepoints);
    TIFFSetField(tif, TIFFTAG_GEOPIXELSCALE, 3,pixscale);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    
    //TIFFSetField(tif, TIFFTAG_TILEWIDTH, 256);
    //TIFFSetField(tif, TIFFTAG_TILELENGTH,  256);
    
    switch (data_type)
    {
        case FLOAT:
            TIFFSetField(tif, TIFFTAG_PREDICTOR, 3);
            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
            TIFFMergeFieldInfo(tif, xtiffFieldInfo, N(xtiffFieldInfo));
            TIFFSetField(tif, GDAL_NODATA, "-9999");
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, max(1, (STRIP_SIZE_DEFAULT * 8) / (width * 32)));
            break;
        case UCHAR:
            TIFFSetField(tif, TIFFTAG_PREDICTOR, 1);
            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 2);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, max(1, (STRIP_SIZE_DEFAULT * 8) / (width * 8)));
            break;
        case UINT16:
            TIFFSetField(tif, TIFFTAG_PREDICTOR, 1);
            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, max(1, (STRIP_SIZE_DEFAULT * 8) / (width * 16)));
            break;
        default:
            break;
    }
}

void SetUpGeoKeys(GTIF *gtif, int projection, int zone, int NS_hemisphere)
{
    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelProjected);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
    GTIFKeySet(gtif, GeogCitationGeoKey, TYPE_ASCII, 0, "WGS 84");
    GTIFKeySet(gtif, GeogAngularUnitsGeoKey, TYPE_SHORT, 1, Angular_Degree);
    GTIFKeySet(gtif, ProjLinearUnitsGeoKey, TYPE_SHORT,  1, Linear_Meter);

    if (projection == 1)  // Polar Stereographic
    {
        GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, "unnamed");
        GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT,  1, GCS_WGS_84);
        GTIFKeySet(gtif, GeogSemiMajorAxisGeoKey, TYPE_DOUBLE, 1, 6378137.0);
        GTIFKeySet(gtif, GeogInvFlatteningGeoKey, TYPE_DOUBLE, 1, 298.257223563);
        GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, 32767);  // user-defined
        GTIFKeySet(gtif, ProjectionGeoKey, TYPE_SHORT, 1, 32767);  // user-defined
        GTIFKeySet(gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1, CT_PolarStereographic);
        GTIFKeySet(gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1, 0.0);
        GTIFKeySet(gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1, 0.0);
        GTIFKeySet(gtif, ProjScaleAtNatOriginGeoKey, TYPE_DOUBLE, 1, 1.0);

        if (NS_hemisphere != 0)
        {
            GTIFKeySet(gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1, 70.0);
            GTIFKeySet(gtif, ProjStraightVertPoleLongGeoKey, TYPE_DOUBLE, 1, -45.0);
        }
        else
        {
            GTIFKeySet(gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1, -71.0);
            GTIFKeySet(gtif, ProjStraightVertPoleLongGeoKey, TYPE_DOUBLE, 1, 0.0);
        }
    }
    else  // UTM
    {
        char GeoKeyStr[500];
        if (NS_hemisphere != 0)
        {
            sprintf(GeoKeyStr, "UTM Zone %d, Northern Hemisphere", zone);
            GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, GeoKeyStr);
            GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, PCS_WGS84_UTM_zone_1N - 1 + zone);
            printf("%s\n",GeoKeyStr);
        }
        else
        {
            sprintf(GeoKeyStr, "UTM Zone %d, Southern Hemisphere", zone);
            GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, GeoKeyStr);
            GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, PCS_WGS84_UTM_zone_1S - 1 + zone);
            printf("%s\n",GeoKeyStr);
        }
    }
}
