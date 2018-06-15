#include "setsmgeo.h"
#include <stdint.h>

#define FLOAT 4
#define UCHAR 1
#define UINT16 12

int WriteGeotiff_DEM(char *filename, void *buffer, size_t width, size_t height, double scale, double minX, double maxY, int projection, int zone, int NS_hemisphere, int data_type)
{
    TIFF *tif;  /* TIFF-level descriptor */
    GTIF *gtif; /* GeoKey-level descriptor */
    
    tif=XTIFFOpen(filename,"w");
    if (!tif)
    {
        printf("failed in XTIFFOpen\n");
        goto failure;
    }
    
    gtif = GTIFNew(tif);
    if (!gtif)
    {
        printf("failed in GTIFNew\n");
        goto failure;
    }
    
    SetUpTIFFDirectory(tif, width, height, scale, minX, maxY, data_type);
    SetUpGeoKeys(gtif, projection, zone, NS_hemisphere);

    for (int row=0; row<height; row++)
    {
	// Write from buffer after appropriate cast and collect status of write
        int status;
        switch (data_type)
        {
        case FLOAT:
            status = TIFFWriteScanline(tif, ((float *)buffer) + (row * width), row, 0);
        break;
        case UCHAR:
            status = TIFFWriteScanline(tif, ((unsigned char *)buffer) + (row * width), row, 0);
        break;
        case UINT16:
            status = TIFFWriteScanline(tif, ((uint16 *)buffer) + (row * width), row, 0);
        break;
        default:
            printf("unrecognized data type: %d\n", data_type);
            goto failure;
        break;
        }

        if (!status) // TODO: TIFFWriteScanline may return -1 on failure:
        {
            TIFFError("WriteGeotiff_DEM","failure in WriteScanline\n");
        }
    }
    
    GTIFWriteKeys(gtif);
    GTIFFree(gtif);
    XTIFFClose(tif);
    return 0;
    
failure:
    if (tif) TIFFClose(tif);
    if (gtif) GTIFFree(gtif);
    return -1;
}


int SetUpTIFFDirectory(TIFF *tif, size_t width, size_t height, double scale, double minX, double maxY, int data_type)
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
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(tif, TIFFTAG_GEOTIEPOINTS, 6,tiepoints);
    TIFFSetField(tif, TIFFTAG_GEOPIXELSCALE, 3,pixscale);

    switch (data_type)
    {
    case FLOAT:
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    break;
    case UCHAR:
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    break;
    case UINT16:
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    break;
    default:
        printf("unrecognized data type:  %d\n", data_type);
    break;
    }
}

int SetUpGeoKeys(GTIF *gtif, int projection, int zone, int NS_hemisphere)
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
        GTIFKeySet(gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1, 70.0);
        GTIFKeySet(gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1, 0.0);
        GTIFKeySet(gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1, 0.0);
        GTIFKeySet(gtif, ProjScaleAtNatOriginGeoKey, TYPE_DOUBLE, 1, 1.0);
        GTIFKeySet(gtif, ProjStraightVertPoleLongGeoKey, TYPE_DOUBLE, 1, -45.0);
    }
    else  // UTM
    {
    char GeoKeyStr[100];
    if (NS_hemisphere != 0)
    {
        sprintf(GeoKeyStr, "UTM Zone %d, Northern Hemisphere", zone);
        GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, GeoKeyStr);
        GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, PCS_WGS84_UTM_zone_1N - 1 + zone);
    }
    else
    {
        sprintf(GeoKeyStr, "UTM Zone %d, Southern Hemisphere", zone);
        GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, GeoKeyStr);
        GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, PCS_WGS84_UTM_zone_1S - 1 + zone);
    }
    }
}

