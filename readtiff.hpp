#include <tiff.h>
#include <tiffio.h>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <cerrno>
#include <type_traits>
#include <tuple>
#include"Typedefine.hpp"

constexpr ttag_t HACK_GDAL_NODATA = 42113;

/** Read image file from flat binary */
template <typename T>
T *_read_bin(const char *filename, CSize *Imagesize, long int *cols, long int *rows, CSize *data_size, T type) {
    (void)type;
    (void)Imagesize;
    FILE *bin  = fopen(filename,"rb");
    if(!bin) {
        printf("Error: failed to open %s for reading\n", filename);
        return NULL;
    }

    long r,a;
    data_size->width    = cols[1] - cols[0];
    data_size->height   = rows[1] - rows[0];
    
    long int data_length = data_size->height*data_size->width;
    
    T* out = (T*)malloc(sizeof(T)*data_length);
    
    T* t_data = (T*)malloc(sizeof(T)*data_size->width);

    for(r = rows[0]; r < rows[1] ; r++)
    {
        fseek(bin,sizeof(T)*(r*Imagesize->width + cols[0]),SEEK_SET);
        
        fread(t_data,sizeof(T),data_size->width,bin);
    
        for(a = cols[0];a<cols[1];a++)
        {
            long int pos = (r-rows[0])*data_size->width + (a-cols[0]);
            
            out[pos] = t_data[a-cols[0]];
        }
    }
    free(t_data);
    fclose(bin);

    return out;
}

/** returns true if the NODATA value should be converted to zero
 *   The second value in the pair is the value that should be converted
 *   to zero.
 */
template <typename T>
std::pair<bool, T> get_nodata_value(TIFF *tif, T type) {
    (void)type;

    // only consider conversion when type is uint16_t
    if(!std::is_same<T, uint16_t>::value) {
        return std::make_pair(false, 0);
    }

    bool needs_convert = false;
    T nodata_val = 0;

    uint32_t count = 0;
    char *nodata_string = NULL;

    if(TIFFGetField(tif, HACK_GDAL_NODATA, &count, &nodata_string) == 1) {
        errno = 0;
        if(std::is_floating_point<T>::value)
            nodata_val = strtof(nodata_string, NULL);
        else
            nodata_val = static_cast<T>(strtol(nodata_string, NULL, 10));
        if(errno)
        {
            printf("WARNING: failed to parse NODATA value\n");
            nodata_val = 0;
        } else if(nodata_val != 0) {
            printf("tiff file has a nonzero NODATA value (%s)."
                   "Will convert to zero when reading.\n", std::to_string(nodata_val).c_str());
            needs_convert = true;
        }
    }

    return std::make_pair(needs_convert, nodata_val);
}

/** Read scanline TIFF file */
template <typename T>
T *_read_scanline_tiff(TIFF *tif, CSize *Imagesize, long int *cols, long int *rows, CSize *data_size, int *num_samples, T type) {
    (void)type;
    (void)Imagesize; // not modified
    printf("tiff file is a scanline tiff\n");
    tsize_t scanline;
    tdata_t buf;
    uint16 s,nsamples;

    uint32_t row;
    
    int a;
    
    // scanline read
    data_size->width    = cols[1] - cols[0];
    data_size->height   = rows[1] - rows[0];
    long int data_length = (long)data_size->height*(long)data_size->width;
    T *out             = (T*)malloc(sizeof(T)*data_length);
    
    scanline        = TIFFScanlineSize(tif);
    
    buf             = _TIFFmalloc(scanline);
    
    TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);
    *num_samples = nsamples;

    // should we convert NODATA pixels to zero?
    T nodata_val = 0;
    bool req_nodata_conversion;
    std::tie(req_nodata_conversion, nodata_val) = get_nodata_value(tif, type);

    // check random access support
    // Random access is supported if no compression is used
    // or if there is one row per strip
    uint16 compression;
    uint32 rows_per_strip;

    bool supports_random_access = false;

    if((       1 == TIFFGetField(tif, TIFFTAG_COMPRESSION,  &compression))
            && 1 == TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rows_per_strip))
    {
        supports_random_access = (
                   (compression    == COMPRESSION_NONE)
                || (rows_per_strip == 1));
    }
    
    //TODO this needs to be updated to read in multi-sample images correctly.
    for(s =0;s< nsamples;s++)
    {
        if(!supports_random_access)
        {
            for (row=0;row<rows[0];row++)
            {
                TIFFReadScanline(tif,buf,row,s);
            }
        }
        for (row=rows[0];row<rows[1];row++)
        {
            T* t_data;
            TIFFReadScanline(tif,buf,row,s);
            t_data = (T*)buf;
            for(a = cols[0];a<cols[1];a++)
            {
                T v = t_data[a];
                if(req_nodata_conversion && v == nodata_val)
                    v = 0;
                out[(long)(row-rows[0])*(long)data_size->width + (long)(a-cols[0])] = v;
            }
        }
    }
    
    _TIFFfree(buf);
    return out;

}

/** Read tiled TIFF file */
template <typename T>
T *_read_tiled_tiff(TIFF *tif, CSize *Imagesize, long int *cols, long int *rows, CSize *data_size, int *num_samples, T type) {
    (void)type;

    T *out = NULL;
    // these are ultimately used to index an array,
    // so enure they are long enough if the array is
    // more than 2^32 elems long
    size_t i, j, k;

    // These need to be 32 bit unsigned per libtiff
    uint32_t tileW, tileL;

    // used to iterate through the TIFF tiles
    uint32_t row, col;
    printf("tiff file is a tile tiff\n");
    int count_W, count_L;
    uint32_t starttileL,starttileW;
    unsigned long start_row,start_col,end_row,end_col;
    tdata_t buf;
    T* t_data;

    uint16_t nsamples;
    
    // TILEWIDTH and TILELENGTH take pointers to uint32
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
    *num_samples = nsamples;
    
    if(nsamples > 1)
    {
        uint16_t planar_config = 999; // just some invalid number
        TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &planar_config);
        if(planar_config != PLANARCONFIG_CONTIG)
        {
            printf("Cannot read multi-band TIFFs with non-contiguous sample layout\n");
            exit(1);
        }
    }

    // should we convert NODATA pixels to zero?
    T nodata_val = 0;
    bool req_nodata_conversion;
    std::tie(req_nodata_conversion, nodata_val) = get_nodata_value(tif, type);


    starttileL      = (uint32_t)(rows[0]/tileL);
    start_row       = starttileL*tileL;
    end_row         = ((int)(rows[1]/tileL)+1)*tileL;
    printf("rows %ld\t%ld\ttileL %d\theight %d\n",rows[0],rows[1],tileL,Imagesize->height);
    if(end_row > Imagesize->height)
        end_row = Imagesize->height;
    
    starttileW      = (uint32_t)(cols[0]/tileW);
    start_col       = starttileW*tileW;
    end_col         = ((int)(cols[1]/tileW)+1)*tileW;
    printf("cols %ld\t%ld\ttileW %d\theight %d\n",cols[0],cols[1],tileW,Imagesize->width);
    if(end_col > Imagesize->width)
        end_col = Imagesize->width;
    
    printf("start %ld\t%ld\t end %ld\t%ld\n",start_col,start_row,end_col,end_row);
    cols[0]         = start_col;
    cols[1]         = end_col;
    rows[0]         = start_row;
    rows[1]         = end_row;
    
    data_size->width = end_col - start_col;
    data_size->height= end_row - start_row;
    
    long int data_length = (long int)data_size->height*(long int)data_size->width;
    
    out             = (T*)malloc(sizeof(T)*data_length * nsamples);
    
    buf             = _TIFFmalloc(TIFFTileSize(tif));
    
    count_L = ceil(data_size->height/(double)tileL);
    count_W = ceil(data_size->width/(double)tileW);
    
    int f_row_end = 0;
    int f_col_end = 0;
    
    if(count_L*tileL > data_size->height)
        f_row_end = tileL + data_size->height - (count_L*tileL);
    
    if(count_W*tileW > data_size->width)
        f_col_end = tileW + data_size->width - (count_W*tileW);
    
    printf("tile info %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",data_size->height, 
           data_size->width,starttileW,starttileL,count_W,count_L,tileW,tileL,f_col_end,f_row_end);
    
    for (row = 0; row < static_cast<unsigned>(count_L); row ++)
    {
        for (col = 0; col < static_cast<unsigned>(count_W); col ++)
        {
            int ret = TIFFReadTile(tif, buf, (col+starttileW)*tileW, (row+starttileL)*tileL, 0,0);
            if(ret < 0) {
                printf("ERROR: TIFFReadTile returned %d for row %d col %d\n", ret, row, col);
                exit(1);
            }
            t_data = (T*)buf;

            size_t end_row = tileL;
            if(f_row_end > 0 && row == static_cast<unsigned>(count_L) - 1)
            {
                end_row = f_row_end;
            }

            size_t end_col = tileW;
            if(f_col_end > 0 && col == static_cast<unsigned>(count_W) - 1)
            {
                end_col = f_col_end;
            }


            for(i = 0; i < end_row; i++)
            {
                for(j = 0; j < end_col; j++)
                {
                    for(k = 0; k < nsamples; k++)
                    {
                        size_t t_row = (row*tileL) + i;
                        size_t t_col = (col*tileW) + j;
                        if(t_row < data_size->height && t_col < data_size->width) {
                            T v = t_data[(i*tileW + j) * nsamples + k];
                            if(req_nodata_conversion && v == nodata_val)
                                v = 0;
                            out[(t_row*data_size->width + t_col) * nsamples + k] = v;
                        }
                    }

                }
            }
        }
    }
    _TIFFfree(buf);
    return out;
}

/** Check that the pixel sample type matches T */
template <typename T>
bool validate_format(TIFF *tif, T type) {
    (void)type;
    uint16_t sample_format, bps;
    int ret;

    // read format 
    ret = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sample_format);
    if(ret != 1) {
        printf("failed to read sample format from tiff file\n");
        return false;
    }

    // read bits per sample
    ret = TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    if(ret != 1) {
        printf("failed to get bits per sample for tiff file\n");
        return false;
    }

    // check that lengths match
    int bits_in_T = sizeof(T) * 8;
    if(bps != bits_in_T) {
        printf("mismatch! tiff has %d bits per sample, but T is %d bits\n", bps, bits_in_T);
        return false;
    }

    bool types_match = false;

    // check that types match
    switch(sample_format) {
    case SAMPLEFORMAT_UINT:
        types_match = std::is_integral<T>::value && std::is_unsigned<T>::value;
        break;
    case SAMPLEFORMAT_INT:
        types_match = std::is_integral<T>::value && !std::is_unsigned<T>::value;
        break;
    case SAMPLEFORMAT_IEEEFP:
        types_match = std::is_floating_point<T>::value;
        break;
    default:
        printf("Unsupported sample format %d\n", sample_format);
        types_match = false;
    };
    
    if(!types_match) {
        printf("sample format does not match template type\n");
    }
    return types_match;
}

/** Read in image data from tiff file */
template <typename T>
T *_read_tiff(const char *filename, CSize *Imagesize, long int *cols, long int *rows, CSize *data_size, int *num_samples,  T type) {

    T *out = NULL;
    
    TIFF *tif = TIFFOpen(filename, "r");
    if(!tif) {
        printf("failed to open tiff file %s\n", filename);
        return NULL;
    }

    if(!validate_format(tif, type)) {
        printf("Error: mismatch in type or size of tiff data and template parameter\n");
        TIFFClose(tif);
        return NULL;
    }

    uint32_t tileW;
    int ret = TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
    if(ret != 1) // TIFFGetField returns 1 on success
    {
        out = _read_scanline_tiff(tif, Imagesize, cols, rows, data_size, num_samples, type);
    }
    else
    {
        out = _read_tiled_tiff(tif, Imagesize, cols, rows, data_size, num_samples, type);
    }
    TIFFClose(tif);
    return out;
}

//** Read in image data from tiff or binary file */
template <typename T>
T *Readtiff_multi(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size, int *ns, T type)
{
    const char *ext = strrchr(filename,'.');
    if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
    {
        return _read_tiff(filename, Imagesize, cols, rows, data_size, ns, type);
    }
    else if(!strcmp("bin",ext+1))
    {
        // bin files don't have multiple samples per pixel
        *ns = 1;
        return _read_bin(filename, Imagesize, cols, rows, data_size, type);
    }
    printf("Error: unknown file extension for %s, cannot read\n", filename);

    return NULL;
}

/** Convert planet multi-band image to panchromatic */
template <typename T>
T *convert_to_pan(const T *out, const CSize *data_size, const int number_samples, const T type)
{
    (void)type;
    if(number_samples < 3)
    {
        printf("cannot convert to panchromatic. Only have %d samples."
               " Need at least 3.\n", number_samples);
    }
    T *pan = static_cast<T*>(calloc(data_size->width * data_size->height, sizeof(T)));
    if(!pan)
    {
        printf("failed to allocate space for panchromatic tiff image! w=%d h=%d\n",
            data_size->width, data_size->height);
        return NULL;
    }
    // assume the sample order is (Blue, Green, Red)
    // TODO we should add functionality to ensure this ordering is correct
    for(long i = 0; i < data_size->width * data_size->height; i++)
    {
        T blue = out[i * number_samples + 0];
        T green = out[i * number_samples + 1];
        T red = out[i * number_samples + 2];

        // important bit here
        pan[i] = round(0.2989 * red + 0.5870 * green + 0.1140 * blue);
    }
    return pan;
}


/** Convert WV multi-band image to NDVI image */
template <typename T>
T **convert_to_ndvi(const T *out, const CSize *data_size, const int number_samples, const T type, const BandInfo band)
{
    (void)type;
    if(number_samples < 8)
    {
        printf("cannot convert to ndvi. Only have %d samples."
               " Need at least 3.\n", number_samples);
    }
    T **ndvi = static_cast<T**>(calloc(3, sizeof(T))); //ndvi, red, nir
    for(int count = 0 ; count < 3 ; count++)
    {
        ndvi[count] = static_cast<T*>(calloc(data_size->width * data_size->height, sizeof(T)));
        if(!ndvi[count])
        {
            printf("failed to allocate space for ndvi tiff image! w=%d h=%d\n",
                data_size->width, data_size->height);
            return NULL;
        }
    }
    
    int red_band = 4;
    int nir_band = 6;
    
    printf("redband abscalfactor %10.9f\n",band.calibrated_abscal_multi[red_band]);
    printf("redband effbw %10.9f\n",band.calibrated_effbw_multi[red_band]);
    printf("redband tdi %lf\n",band.tdi_multi[red_band]);
    
    printf("nirband abscalfactor %10.9f\n",band.calibrated_abscal_multi[nir_band]);
    printf("nirband effbw %10.9f\n",band.calibrated_effbw_multi[nir_band]);
    printf("nirband tdi %lf\n",band.tdi_multi[nir_band]);
    
    double minrefl = 0.000001;

    // assume the sample order (8 bands) is (coastal blue, blue, green, yellow, red [4], red edge, nir1[6], nir2)
    // TODO we should add functionality to ensure this ordering is correct
    for(long i = 0; i < data_size->width * data_size->height; i++)
    {
        double red_ori = (double)(out[i * number_samples + red_band]);
        double nir_ori = (double)(out[i * number_samples + nir_band]);

        ndvi[1][i] = out[i * number_samples + red_band];
        ndvi[2][i] = out[i * number_samples + nir_band];

        if(red_ori >0 && nir_ori > 0)
        {
            double red = red_ori*(double)band.calibrated_abscal_multi[red_band] + (double)band.calibrated_effbw_multi[red_band];
            double nir = nir_ori*(double)band.calibrated_abscal_multi[nir_band] + (double)band.calibrated_effbw_multi[nir_band];
            
            if(red < 0)
            {
                printf("red less than 0 %f\n",red);
                red = minrefl;
            }
            if(red > 1)
            {
                printf("red more than 1 %f\n",red);
                red = 1.0;
            }
            if(nir < 0)
            {
                printf("nir less than 0 %f\n",nir);
                nir = minrefl;
            }
            if(nir > 1)
            {
                printf("nir more than 1 %f\n",nir);
                nir = 1.0;
            }
            //double red = (double)out[i * number_samples + red_band];
            //double nir = (double)out[i * number_samples + nir_band];
            // important bit here
            double real_ndvi = double(nir - red)/double(nir + red);
            if(real_ndvi < -1)
                printf("ndvi less than -1 %f\n", real_ndvi);
            if(real_ndvi > 1)
                printf("ndvi more than 1 %f\n", real_ndvi);
            //ndvi[i] = round(real_ndvi*(2047/2.0) + 2049/2);
            ndvi[0][i] = round(real_ndvi*1000.0 + 1000.0);

            //printf("red_ori %f\t nir_ori %f\n",red_ori,nir_ori);
            //printf("red %f\t nir %f\t real_ndvi %f\t ndvi[0][i] %d\n",red, nir,real_ndvi,ndvi[0][i]);
            //exit(1);
        }
        else
            ndvi[0][i] = 0;
    }
    return ndvi;
}

/** Read and return pointer to TIFF
 *
 * Arguments:
 *      filename - name of TIFF file to read
 *      Imagesize - Size of the image to read
 *      cols - (IN/OUT) array of length 2. First element is the column of
 *          the starting pixel. Second element is the column of the end
 *          pixel. These are updated if the start/end pixels shift based
 *          on TIFF tile boundaries. End is exlusive, start is inclusive.
 *      rows - (IN/OUT) same as cols, but for start/end pixel rows in images
 *      data_size - (OUT) size of the returned image. Length of returned
 *          image is width*height
 *      type - Return type for function. Variable value unused. I.e. to
 *          return data as floats, pass a float here
 *
 *  Returns a buffer holding image data in row-major order.
 */
template <typename T>
T *Readtiff_T(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size, T type) {
    int number_samples = 1;
    T *out = Readtiff_multi(filename, Imagesize, cols, rows, data_size, &number_samples, type);

    if(out != NULL && number_samples > 1) {
        printf("Read tiff file with %d samples, converting to panchromatic\n", number_samples);
        T *pan = convert_to_pan(out, data_size, number_samples, type);
        free(out);
        out = pan;
    }

    return out;
}

template <typename T>
T **Readtiff_T_NDVI(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size, T type, BandInfo band) {
    int number_samples = 1;
    T *out = Readtiff_multi(filename, Imagesize, cols, rows, data_size, &number_samples, type);

    if(out != NULL && number_samples > 1) {
        printf("Read tiff file with %d samples, converting to ndvi\n", number_samples);
        T **pan = convert_to_ndvi(out, data_size, number_samples, type, band);
        free(out);
        //out = pan;
        return pan;
    }
    else
        return NULL;

}


