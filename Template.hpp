//
//  Template.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef Template_h
#define Template_h

#include "Typedefine.hpp"

//declaration
template <typename T>
T* CreateImagePyramid(T* _input, CSize _img_size, int _filter_size, double _sigma);
template <typename T>
T BilinearResampling(T* input, const CSize img_size, D2DPOINT query_pt);
template <typename T>
T *Readtiff_T(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size, T type);



inline double SQRT(D2DPOINT a);
inline double SQRT(D2DPOINT a, D2DPOINT b);
inline double SQRT(D3DPOINT a, int dimension = 3);
inline double SQRT(D3DPOINT a, D3DPOINT b, int dimension = 3);

inline short DoubleToSignedChar_result(double val);
inline double SignedCharToDouble_result(short val);

inline short DoubleToSignedChar_grid(double val);
inline double SignedCharToDouble_grid(short val);

inline short DoubleToSignedChar_voxel(double val);
inline double SignedCharToDouble_voxel(short val);

inline signed char FloatToSignedChar(float val);
inline float SignedCharToFloat(signed char val);

//definition
template <typename T>
T* CreateImagePyramid(T* _input, CSize _img_size, int _filter_size, double _sigma)
{
    //_filter_size = 7, sigma = 1.6
    double sigma = _sigma;
    double temp,scale;
    double sum = 0;
    double** GaussianFilter;
    CSize result_size;
    T* result_img;

    GaussianFilter = (double**)malloc(sizeof(double*)*_filter_size);
    for(int i=0;i<_filter_size;i++)
        GaussianFilter[i] = (double*)malloc(sizeof(double)*_filter_size);

    
    result_size.width = _img_size.width/2;
    result_size.height = _img_size.height/2;
    scale=sqrt(2*PI)*sigma;
    
    int half_filter_size = (int)(_filter_size/2);
    
    result_img = (T*)malloc(sizeof(T)*result_size.height*result_size.width);

    for(int i=-half_filter_size;i<=half_filter_size;i++)
    {
        for(int j=-half_filter_size;j<=half_filter_size;j++)
        {
            temp = -1.0*(i*i+j*j)/(2*sigma*sigma);
            GaussianFilter[i+half_filter_size][j+half_filter_size]=exp(temp)/scale;
            sum += exp(temp)/scale;
        }
    }

#pragma omp parallel for schedule(guided)
    for(int i=-half_filter_size;i<=half_filter_size;i++)
    {
        for(int j=-half_filter_size;j<=half_filter_size;j++)
        {
            GaussianFilter[i+half_filter_size][j+half_filter_size]/=sum;
        }
    }

#pragma omp parallel for private(temp) schedule(guided)
    for(long int r=0;r<result_size.height;r++)
    {
        for(long int c=0;c<result_size.width;c++)
        {
            double temp_v = 0;
            int count = 0;
            for(int l=-half_filter_size;l<=half_filter_size;l++)
            {
                for(int k=-half_filter_size;k<=half_filter_size;k++)
                {
                    //r'->2r+m, c'->2c+n
                    if( (2*r + l) >= 0 && (2*c + k) >= 0 &&
                        (2*r + l) < _img_size.height && (2*c + k) < _img_size.width)
                    {
                        if(_input[(2*r + l)*_img_size.width +(2*c + k)] > Nodata)
                        {
                            temp_v += GaussianFilter[l + half_filter_size][k + half_filter_size]*_input[(2*r + l)*_img_size.width +(2*c + k)];
                            count ++;
                        }
                    }
                }
            }

            if(count == _filter_size*_filter_size)
                result_img[r*result_size.width + c] = (T)temp_v;
            else
                result_img[r*result_size.width + c] = _input[(2*r)*_img_size.width +(2*c)];
        }
    }
    
    for(int i=0;i<_filter_size;i++)
        if(GaussianFilter[i])
            free(GaussianFilter[i]);

    if(GaussianFilter)
        free(GaussianFilter);
    

    return result_img;
}

template <typename T>
T BilinearResampling(T* input, const CSize img_size, D2DPOINT query_pt)
{
    const long data_length = (long)img_size.width*(long)img_size.height;
    
    long int index1,index2,index3, index4;
    T value1, value2, value3, value4;
    double value;
    
    const long t_col_int   = (long int)(query_pt.m_X + 0.01);
    const long t_row_int   = (long int)(query_pt.m_Y + 0.01);
    
    const double dcol        = query_pt.m_X - t_col_int;
    const double drow        = query_pt.m_Y - t_row_int;
    
    
    index1  = (t_col_int    ) + (t_row_int    )*(long)img_size.width;
    index2  = (t_col_int + 1) + (t_row_int    )*(long)img_size.width;
    index3  = (t_col_int    ) + (t_row_int + 1)*(long)img_size.width;
    index4  = (t_col_int + 1) + (t_row_int + 1)*(long)img_size.width;
    
    if(index1 >= 0 && index1 < data_length && index2 >= 0 && index2 < data_length && index3 >= 0 && index3 < data_length && index4 >= 0 && index4 < data_length && t_col_int >= 0 && (t_col_int + 1) < img_size.width && t_row_int >= 0 && (t_row_int + 1) < img_size.height)
    {
        value1      = input[index1];
        value2      = input[index2];
        value3      = input[index3];
        value4      = input[index4];
    
        if(value1 > Nodata && value2 > Nodata && value3 > Nodata && value4 > Nodata)
            value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow) + value3*(1-dcol)*drow + value4*dcol*drow;
        else if(value1 > Nodata)
            value = input[index1];
        else if(value2 > Nodata)
            value = input[index2];
        else if(value3 > Nodata)
            value = input[index3];
        else if(value4 > Nodata)
            value = input[index4];
        else
            value = Nodata;
        
    }
    else
    {
        if(index1 >= 0 && index1 < data_length      && t_col_int >= 0       && (t_col_int) < img_size.width     && t_row_int >= 0 && (t_row_int) < img_size.height)
            value = input[index1];
        else if(index2 >= 0 && index2 < data_length && t_col_int + 1 >= 0   && (t_col_int + 1) < img_size.width && t_row_int >= 0 && (t_row_int) < img_size.height)
            value = input[index2];
        else if(index3 >= 0 && index3 < data_length && t_col_int >= 0       && (t_col_int) < img_size.width     && t_row_int + 1 >= 0 && (t_row_int + 1) < img_size.height)
            value = input[index3];
        else if(index4 >= 0 && index4 < data_length && t_col_int + 1 >= 0   && (t_col_int + 1) < img_size.width && t_row_int + 1 >= 0 && (t_row_int + 1) < img_size.height)
            value = input[index4];
        else
            value = Nodata;
    }
    
    return (T)value;
}

template <typename T>
T *Readtiff_T(const char *filename, CSize *Imagesize,long int *cols,long int *rows, CSize *data_size, T type)
{
    T *out;
    FILE *bin;
    int check_ftype = 1; // 1 = tif, 2 = bin
    TIFF *tif = NULL;
    const char *ext = strrchr(filename,'.');
    
    if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
    {
        tif  = TIFFOpen(filename,"r");
        check_ftype = 1;
    }
    else if(!strcmp("bin",ext+1))
    {
        bin  = fopen(filename,"rb");
        check_ftype = 2;
    }
    
    if(check_ftype == 1 && tif)
    {
        int i,j,row, col, tileW;
        
        tileW = -1;
        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
        if(tileW < 0)
        {
            printf("NO TILE\n");
            tsize_t scanline;
            tdata_t buf;
            uint16 s,nsamples;
            
            int a;
            
            // scanline read
            data_size->width    = cols[1] - cols[0];
            data_size->height   = rows[1] - rows[0];
            long int data_length = (long)data_size->height*(long)data_size->width;
            out             = (T*)malloc(sizeof(T)*data_length);
            
            scanline        = TIFFScanlineSize(tif);
            
            buf             = _TIFFmalloc(scanline);
            
            TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);
            
            for(s =0;s< nsamples;s++)
            {
                for (row=0;row<rows[0];row++)
                    TIFFReadScanline(tif,buf,row,s);
                for (row=rows[0];row<rows[1];row++)
                {
                    T* t_data;
                    TIFFReadScanline(tif,buf,row,s);
                    t_data = (T*)buf;
#pragma omp parallel for private(a) schedule(guided)
                    for(a = cols[0];a<cols[1];a++)
                    {
                        out[(long)(row-rows[0])*(long)data_size->width + (long)(a-cols[0])] = t_data[a];
                    }
                }
            }
            
            _TIFFfree(buf);
        }
        else
        {
            printf("tile\n");
            int tileL,count_W,count_L,starttileL,starttileW;
            uint16 start_row,start_col,end_row,end_col;
            tdata_t buf;
            T* t_data;
            
            TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
            TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);
            
            starttileL      = (int)(rows[0]/tileL);
            start_row       = starttileL*tileL;
            end_row         = ((int)(rows[1]/tileL)+1)*tileL;
            if(end_row > Imagesize->height)
                end_row = Imagesize->height;
            
            starttileW      = (int)(cols[0]/tileW);
            start_col       = starttileW*tileW;
            end_col         = ((int)(cols[1]/tileW)+1)*tileW;
            if(end_col > Imagesize->width)
                end_col = Imagesize->width;
            
            
            cols[0]         = start_col;
            cols[1]         = end_col;
            rows[0]         = start_row;
            rows[1]         = end_row;
            
            data_size->width = end_col - start_col;
            data_size->height= end_row - start_row;
            
            long int data_length = (long int)data_size->height*(long int)data_size->width;
            
            out             = (T*)malloc(sizeof(T)*data_length);
            
            buf             = _TIFFmalloc(TIFFTileSize(tif));
            
            count_L = ceil(data_size->height/(double)tileL);
            count_W = ceil(data_size->width/(double)tileW);
            
            int f_row_end = 0;
            int f_col_end = 0;
            
            if(count_L*tileL > data_size->height)
                f_row_end = tileL + data_size->height - (count_L*tileL);
            
            if(count_W*tileW > data_size->width)
                f_col_end = tileW + data_size->width - (count_W*tileW);
            
            printf("tile info %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",data_size->height,data_size->width,starttileW,starttileL,count_W,count_L,tileW,tileL,f_col_end,f_row_end);
            
            for (row = 0; row < count_L; row ++)
            {
                for (col = 0; col < count_W; col ++)
                {
                    TIFFReadTile(tif, buf, (col+starttileW)*tileW, (row+starttileL)*tileL, 0,0);
                    t_data = (T*)buf;
                    
                    if(f_row_end > 0 && f_col_end > 0)
                    {
                        if(row == count_L-1 && col == count_W -1)
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<f_row_end;i++)
                            {
                                for (j=0;j<f_col_end;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                        else if(row == count_L-1)
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<f_row_end;i++)
                            {
                                for (j=0;j<tileW;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                            
                        }
                        else if(col == count_W -1)
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<tileL;i++)
                            {
                                for (j=0;j<f_col_end;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                        else
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<tileL;i++)
                            {
                                for (j=0;j<tileW;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                    }
                    else if(f_row_end > 0)
                    {
                        if(row == count_L-1)
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<f_row_end;i++)
                            {
                                for (j=0;j<tileW;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                            
                        }
                        else
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<tileL;i++)
                            {
                                for (j=0;j<tileW;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                    }
                    else if(f_col_end > 0)
                    {
                        if(col == count_W -1)
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<tileL;i++)
                            {
                                for (j=0;j<f_col_end;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                        else
                        {
#pragma omp parallel for private(i,j) schedule(guided)
                            for (i=0;i<tileL;i++)
                            {
                                for (j=0;j<tileW;j++)
                                {
                                    int t_row = (row*tileL) + i;
                                    int t_col = (col*tileL) + j;
                                    if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                        out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                                }
                            }
                        }
                    }
                    else
                    {
                       
#pragma omp parallel for private(i,j) schedule(guided)
                        for (i=0;i<tileL;i++)
                        {
                            for (j=0;j<tileW;j++)
                            {
                                int t_row = (row*tileL) + i;
                                int t_col = (col*tileL) + j;
                                if(t_row >= 0 && t_row < data_size->height && t_col >= 0 && t_col < data_size->width)
                                    out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                            }
                        }
                    }
                }
            }
            _TIFFfree(buf);
        }
        TIFFClose(tif);
    }
    else if(check_ftype == 2 && bin)
    {
        long r,c,a;
        data_size->width    = cols[1] - cols[0];
        data_size->height   = rows[1] - rows[0];
        
        long int data_length = data_size->height*data_size->width;
        
        out             = (T*)malloc(sizeof(T)*data_length);
        
        for(r = rows[0]; r < rows[1] ; r++)
        {
            fseek(bin,sizeof(T)*(r*Imagesize->width + cols[0]),SEEK_SET);
            T* t_data = (T*)malloc(sizeof(T)*data_size->width);
            
            fread(t_data,sizeof(T),data_size->width,bin);
        
            for(a = cols[0];a<cols[1];a++)
            {
                long int pos = (r-rows[0])*data_size->width + (a-cols[0]);
                
                out[pos] = t_data[a-cols[0]];
            }
            free(t_data);
        }
        fclose(bin);
    }

    return out;
}

inline double SQRT(D2DPOINT a)
{
    return sqrt( SQ(a.m_X) + SQ(a.m_Y) );
}

inline double SQRT(D2DPOINT a, D2DPOINT b)
{
    return sqrt( SQ(a.m_X - b.m_X) + SQ(a.m_Y - b.m_Y) );
}

inline double SQRT(D3DPOINT a, int dimension)
{
    if(dimension == 2)
        return sqrt( SQ(a.m_X) + SQ(a.m_Y) );
    else
        return sqrt( SQ(a.m_X) + SQ(a.m_Y) + SQ(a.m_Z) );
}

inline double SQRT(D3DPOINT a, D3DPOINT b, int dimension)
{
    if(dimension == 2)
        return sqrt( SQ(a.m_X - b.m_X) + SQ(a.m_Y - b.m_Y) );
    else
        return sqrt( SQ(a.m_X - b.m_X) + SQ(a.m_Y - b.m_Y) + SQ(a.m_Z - b.m_Z) );
}

inline short DoubleToSignedChar_result(double val)
{
    return (short)round(val*1000.0);
}

inline double SignedCharToDouble_result(short val)
{
    return (double)(val)/1000.0;
}

inline short DoubleToSignedChar_grid(double val)
{
    return (short)round(val*1000.0);
}

inline double SignedCharToDouble_grid(short val)
{
    return (double)(val)/1000.0;
}

inline short DoubleToSignedChar_voxel(double val)
{
    return (short)(val*1000.0);
}

inline double SignedCharToDouble_voxel(short val)
{
    return (double)(val)/1000.0;
}


inline signed char FloatToSignedChar(float val)
{
    return (signed char)(val*100.0);
}

inline float SignedCharToFloat(signed char val)
{
    return (float)(val/100.0);
}



#endif /* Template_h */
