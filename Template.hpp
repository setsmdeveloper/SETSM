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
template <typename T>
void CoregParam_Image(ProInfo *proinfo, int ti, uint8 Pyramid_step, double *ImageAdjust, uint8 Template_size, T *Image_ref, CSize Imagesizes_ref, T *Image_tar, CSize Imagesizes_tar, double *Boundary_ref, double *Boundary_tar, D2DPOINT grid_dxy_ref, D2DPOINT grid_dxy_tar, int grid_space, double *over_Boundary, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, vector<D2DPOINT> &matched_MPs, vector<D2DPOINT> &matched_MPs_ref, vector<D2DPOINT> &MPs);
template <typename T>
bool postNCC_ortho(uint8 Pyramid_step, D2DPOINT Left, D2DPOINT Right, double subA[][6],double TsubA[][9],double InverseSubA[][6], uint8 Template_size, CSize leftsize, CSize rightsize, T* _leftimage, T* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, D2DPOINT *peak_pos);


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

inline void Set6by6Matrix(double subA[][6], double TsubA[][9], double InverseSubA[][6]);

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
    scale=1/(sqrt(2*PI)*sigma);
    
    int half_filter_size = (int)(_filter_size/2);
    
    result_img = (T*)malloc(sizeof(T)*result_size.height*result_size.width);
    double tmp = -1/(2*sigma*sigma);

#pragma omp parallel for schedule(guided) private(temp) collapse(2) reduction(+:sum)
    for(int i=-half_filter_size;i<=half_filter_size;i++)
    {
        for(int j=-half_filter_size;j<=half_filter_size;j++)
        {
            temp = (i*i+j*j)*tmp; //-1.0*(i*i+j*j)/(2*sigma*sigma);
            GaussianFilter[i+half_filter_size][j+half_filter_size]=exp(temp)*scale;
            sum += exp(temp)*scale;
        }
    }

#pragma omp parallel for schedule(guided) collapse(2)
    for(int i=-half_filter_size;i<=half_filter_size;i++)
    {
        for(int j=-half_filter_size;j<=half_filter_size;j++)
        {
            GaussianFilter[i+half_filter_size][j+half_filter_size]/=sum;
        }
    }

#pragma omp parallel for schedule(guided) collapse(2)
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
            long start_row,start_col,end_row,end_col;
            tdata_t buf;
            T* t_data;
            
            TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
            TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);
            
            starttileL      = (int)(rows[0]/tileL);
            start_row       = starttileL*tileL;
            end_row         = ((int)(rows[1]/tileL)+1)*tileL;
            printf("rows %d\t%d\ttileL %d\theight %d\n",rows[0],rows[1],tileL,Imagesize->height);
            if(end_row > Imagesize->height)
                end_row = Imagesize->height;
            
            starttileW      = (int)(cols[0]/tileW);
            start_col       = starttileW*tileW;
            end_col         = ((int)(cols[1]/tileW)+1)*tileW;
            printf("cols %d\t%d\ttileW %d\theight %d\n",cols[0],cols[1],tileW,Imagesize->width);
            if(end_col > Imagesize->width)
                end_col = Imagesize->width;
            
            printf("start %d\t%d\t end %d\t%d\n",start_col,start_row,end_col,end_row);
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

template <typename T>
void CoregParam_Image(ProInfo *proinfo, int ti, uint8 Pyramid_step, double *ImageAdjust, uint8 Template_size, T *Image_ref, CSize Imagesizes_ref, T *Image_tar, CSize Imagesizes_tar, double *Boundary_ref, double *Boundary_tar, D2DPOINT grid_dxy_ref, D2DPOINT grid_dxy_tar, int grid_space, double *over_Boundary, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, vector<D2DPOINT> &matched_MPs, vector<D2DPOINT> &matched_MPs_ref, vector<D2DPOINT> &MPs)
{
    double subA[9][6] = {0};
    double TsubA[6][9] = {0};
    double InverseSubA[6][6] = {0};
    
    Set6by6Matrix(subA, TsubA, InverseSubA);
    
    double GridSize_width = over_Boundary[2] - over_Boundary[0];
    double GridSize_height = over_Boundary[3] - over_Boundary[1];
    CSize grid_size(floor(GridSize_width/grid_space), floor(GridSize_height/grid_space));
    printf("Grid_size %d\t%d\n",grid_size.width,grid_size.height);
    
    char temp_path[500];
    //vector<D2DPOINT> MPs;
    if(MPs.size() == 0)
    {
        for(long row = 0 ; row < grid_size.height ; row ++)
        {
            for(long col = 0 ; col < grid_size.width ; col ++)
            {
                long index = row*(long)grid_size.width + col;
                D2DPOINT temp_pts(over_Boundary[0] + col*grid_space, over_Boundary[1] + row*grid_space);
                MPs.push_back(temp_pts);
            }
        }
        printf("no rock masked\n");
    }
    else
        printf("rock masked\n");
    
    long total_grid_counts = MPs.size();
    
    printf("total pts %d\n",total_grid_counts);
    
    sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",proinfo->save_filepath,ti,Pyramid_step);
    FILE *fid_pts = fopen(temp_path,"w");
    sprintf(temp_path,"%s/txt/CoregStat_Image_ID_%d_level_%d.txt",proinfo->save_filepath,ti,Pyramid_step);
    FILE *fid_stat = fopen(temp_path,"w");
    
    D3DPOINT* save_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),total_grid_counts);
    int* mps_index_save = (int*)calloc(sizeof(int),total_grid_counts);
    
    bool check_stop = false;
    const int max_iteration = 100;
    *iter_count = 1;
    while(!check_stop && *iter_count < max_iteration)
    {
        double sum_weight_X     = 0;
        double sum_weight_Y     = 0;
        double sum_max_roh      = 0;
        
        //calculation image coord from object coord by RFM in left and right image
        const double b_factor             = pow(2.0,(2-Pyramid_step))*2;
        const uint8 Half_template_size   = (int)(Template_size/2.0);
        
        int count_pts = 0;
#pragma omp parallel for reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
        for(long mps_index = 0; mps_index < total_grid_counts ; mps_index++)
        {
            D2DPOINT Startpos;
            
            D2DPOINT Left_Imagecoord, Right_Imagecoord,Left_Imagecoord_p, Right_Imagecoord_p;
            Left_Imagecoord_p.m_X = (MPs[mps_index].m_X - Boundary_ref[0])/grid_dxy_ref.m_X;
            Left_Imagecoord_p.m_Y = (Boundary_ref[3] - MPs[mps_index].m_Y)/grid_dxy_ref.m_Y;
            
            Right_Imagecoord_p.m_X = (MPs[mps_index].m_X - Boundary_tar[0])/grid_dxy_tar.m_X + ImageAdjust[1];
            Right_Imagecoord_p.m_Y = (Boundary_tar[3] - MPs[mps_index].m_Y)/grid_dxy_tar.m_Y + ImageAdjust[0];
            
            Left_Imagecoord = OriginalToPyramid_single(Left_Imagecoord_p,Startpos,Pyramid_step);
            Right_Imagecoord = OriginalToPyramid_single(Right_Imagecoord_p,Startpos,Pyramid_step);
            
            if(   Left_Imagecoord.m_Y  > Half_template_size*b_factor    + 10                    && Left_Imagecoord.m_X  > Half_template_size*b_factor + 10
               && Left_Imagecoord.m_Y  < Imagesizes_ref.height - Half_template_size*b_factor - 10    && Left_Imagecoord.m_X  < Imagesizes_ref.width - Half_template_size*b_factor - 10
               && Right_Imagecoord.m_Y > Half_template_size*b_factor + 10                    && Right_Imagecoord.m_X > Half_template_size*b_factor + 10
               && Right_Imagecoord.m_Y < Imagesizes_tar.height - Half_template_size*b_factor - 10    && Right_Imagecoord.m_X < Imagesizes_tar.width - Half_template_size*b_factor - 10)
            {
                D2DPOINT Left(Left_Imagecoord.m_X, Left_Imagecoord.m_Y);
                D2DPOINT Right(Right_Imagecoord.m_X, Right_Imagecoord.m_Y);
                long index_l = (long)Left.m_Y*(long)Imagesizes_ref.width + (long)Left.m_X;
                long index_r = (long)Right.m_Y*(long)Imagesizes_tar.width + (long)Right.m_X;
                if( (index_l > 0 && index_l < (long)Imagesizes_ref.height*(long)Imagesizes_ref.width) && (index_r > 0 && index_r < (long)Imagesizes_tar.height*(long)Imagesizes_tar.width) && Image_ref[index_l] > 0 && Image_tar[index_r] > 0 )
                {
                    D2DPOINT peak_pos;
                    double t_sum_weight_X       = 0;
                    double t_sum_weight_Y       = 0;
                    double t_sum_max_roh        = 0;
                    if(postNCC_ortho(Pyramid_step, Left, Right, subA, TsubA, InverseSubA, Template_size, Imagesizes_ref, Imagesizes_tar, Image_ref, Image_tar, &t_sum_weight_X, &t_sum_weight_Y, &t_sum_max_roh, &peak_pos))
                    {
                        sum_weight_X += t_sum_weight_X;
                        sum_weight_Y += t_sum_weight_Y;
                        sum_max_roh  += t_sum_max_roh;
                        
                        D3DPOINT save_pts_tmp;
                        save_pts_tmp.m_X = peak_pos.m_X*pwrtwo(Pyramid_step);
                        save_pts_tmp.m_Y = peak_pos.m_Y*pwrtwo(Pyramid_step);
                        save_pts_tmp.m_Z = 0;
                        save_pts_tmp.flag = 1;
                        
                        save_pts[count_pts] = save_pts_tmp;
                        mps_index_save[count_pts] = mps_index;
                        count_pts++;
                     }
                }
            }
        }
        
        if(count_pts > 10)
        {
            double shift_X             = sum_weight_X/sum_max_roh*pwrtwo(Pyramid_step);
            double shift_Y             = sum_weight_Y/sum_max_roh*pwrtwo(Pyramid_step);
            
            double sum_var_x = 0;
            double sum_var_y = 0;
            
            for(long c_i = 0 ; c_i < count_pts ; c_i++)
            {
                sum_var_x += (shift_X - save_pts[c_i].m_X)*(shift_X - save_pts[c_i].m_X);
                sum_var_y += (shift_Y - save_pts[c_i].m_Y)*(shift_Y - save_pts[c_i].m_Y);
            }
            
            adjust_std->m_X = sqrt(sum_var_x/count_pts);
            adjust_std->m_Y = sqrt(sum_var_y/count_pts);
            
            if(fabs(shift_Y) < 0.01 && fabs(shift_X) < 0.01)
                check_stop = true;
            
            fprintf(fid_stat,"%d\t%f\t%f\t%f\t%f\t%d\n",*iter_count,shift_X,shift_Y,ImageAdjust[1],ImageAdjust[0],count_pts);
            
            shift_X             += ImageAdjust[1];
            shift_Y             += ImageAdjust[0];
            
            ImageAdjust[1]      = shift_X;
            ImageAdjust[0]      = shift_Y;
        }
        else
            check_stop = true;
        
        (*iter_count)++;
        
        if(check_stop || *iter_count >= max_iteration)
        {
            for(long cc = 0 ; cc < count_pts ; cc++)
            {
                matched_MPs_ref.push_back(MPs[mps_index_save[cc]]);
                
                D2DPOINT temp_pts(MPs[mps_index_save[cc]].m_X + ImageAdjust[1]*grid_dxy_tar.m_X, MPs[mps_index_save[cc]].m_Y - ImageAdjust[0]*grid_dxy_tar.m_Y);
                matched_MPs.push_back(temp_pts);
                
                fprintf(fid_pts,"%8.2f\t%8.2f\t%8.2f\t%8.2f\n",matched_MPs[cc].m_X,matched_MPs[cc].m_Y,matched_MPs[cc].m_X,matched_MPs[cc].m_Y);
            }
            
            *avg_rho = sum_max_roh/(double)(count_pts);
            
            fclose(fid_pts);
            fclose(fid_stat);
        }
    }
    free(save_pts);
    free(mps_index_save);
    
    MPs.clear();
}

template <typename T>
bool postNCC_ortho(uint8 Pyramid_step, D2DPOINT Left, D2DPOINT Right, double subA[][6],double TsubA[][9],double InverseSubA[][6], uint8 Template_size, CSize leftsize, CSize rightsize, T* _leftimage, T* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, D2DPOINT *peak_pos)
{
    bool check_pt = false;
    
    const int Half_template_size  = (int)(Template_size/2);
    const int half_mask_size      = 1;
    
    double *result_rho  = (double*)calloc(9,sizeof(double));
    double *XX          = (double*)calloc(6,sizeof(double));
    double *ATLT        = (double*)calloc(6,sizeof(double));
    int i, j, k;
    uint8 cell_count = 0;
    
    for(j=0;j<9;j++)
        result_rho[j]       = -1.00;
    
    for(long mask_row = - half_mask_size ; mask_row <= half_mask_size ; mask_row++)
    {
        for(long mask_col = - half_mask_size ; mask_col <= half_mask_size ; mask_col++)
        {
            double Sum_LR = 0;
            double Sum_L = 0;
            double Sum_R = 0;
            double Sum_L2 = 0;
            double Sum_R2 = 0;
            double Sum_LR_2 = 0;
            double Sum_L_2 = 0;
            double Sum_R_2 = 0;
            double Sum_L2_2 = 0;
            double Sum_R2_2 = 0;
            double Sum_LR_3 = 0;
            double Sum_L_3 = 0;
            double Sum_R_3 = 0;
            double Sum_L2_3 = 0;
            double Sum_R2_3 = 0;
            int Count_N[3] = {0};
            
            for(long row = -Half_template_size; row <= Half_template_size ; row++)
            {
                for(long col = -Half_template_size; col <= Half_template_size ; col++)
                {
                    double radius  = sqrt((double)(row*row + col*col));
                    if(radius <= Half_template_size-1)
                    {
                        double pos_row_left      = Left.m_Y + row;
                        double pos_col_left      = Left.m_X + col;
                        
                        double pos_row_right     = Right.m_Y + row + mask_row;
                        double pos_col_right     = Right.m_X + col + mask_col;
                        
                        if(pos_row_right-3 >= 0 && pos_row_right+3 < rightsize.height && pos_col_right-3 >= 0 && pos_col_right+3 < rightsize.width &&
                           pos_row_left-3 >= 0 && pos_row_left+3 < leftsize.height && pos_col_left-3 >= 0 && pos_col_left+3 < leftsize.width)
                        {
                            //interpolate left_patch
                            double dx = pos_col_left - (int) (pos_col_left);
                            double dy = pos_row_left - (int) (pos_row_left);
                            double left_patch;
                            double right_patch;
                            double dxdy = dx * dy;
                            long position = (long) (pos_col_left) + (long) (pos_row_left) * leftsize.width;
                            
                            left_patch = (double) (_leftimage[position]) * (1 - dx - dy + dxdy) + (double) (_leftimage[position + 1]) * (dx - dxdy) +
                                (double) (_leftimage[position + leftsize.width]) * (dy - dxdy) + (double) (_leftimage[position + 1 + leftsize.width]) * (dxdy);
                           
                            //interpolate right_patch
                            dx = pos_col_right - (int) (pos_col_right);
                            dy = pos_row_right - (int) (pos_row_right);
                            dxdy = dx * dy;
                            position = (long) (pos_col_right) + (long) (pos_row_right) * rightsize.width;
                            
                            right_patch = (double) (_rightimage[position]) * (1 - dx - dy + dxdy) + (double) (_rightimage[position + 1]) * (dx - dxdy) +
                                (double) (_rightimage[position + rightsize.width]) * (dy - dxdy) + (double) (_rightimage[position + 1 + rightsize.width]) * (dxdy);
                            
                            if(left_patch > 1 && right_patch > 1)
                            {
                                Count_N[0]++;
                                
                                Sum_LR            = Sum_LR + left_patch*right_patch;
                                Sum_L             = Sum_L  + left_patch;
                                Sum_R             = Sum_R  + right_patch;
                                Sum_L2            = Sum_L2 + left_patch*left_patch;
                                Sum_R2            = Sum_R2 + right_patch*right_patch;
                                
                                int size_1, size_2;
                                size_1        = (int)(Half_template_size/2);
                                if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
                                {
                                    if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
                                    {
                                        Sum_LR_2  = Sum_LR_2 + left_patch*right_patch;
                                        Sum_L_2   = Sum_L_2  + left_patch;
                                        Sum_R_2   = Sum_R_2  + right_patch;
                                        Sum_L2_2  = Sum_L2_2 + left_patch*left_patch;
                                        Sum_R2_2  = Sum_R2_2 + right_patch*right_patch;
                                        Count_N[1]++;
                                    }
                                }
                                
                                size_2        = size_1 + (int)((size_1/2.0) + 0.5);
                                if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
                                {
                                    if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
                                    {
                                        Sum_LR_3  = Sum_LR_3 + left_patch*right_patch;
                                        Sum_L_3   = Sum_L_3  + left_patch;
                                        Sum_R_3   = Sum_R_3  + right_patch;
                                        Sum_L2_3  = Sum_L2_3 + left_patch*left_patch;
                                        Sum_R2_3  = Sum_R2_3 + right_patch*right_patch;
                                        Count_N[2]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            if(Count_N[0] > 0)
            {
                double N               = Count_N[0];
                double val1          = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
                double val2          = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
                double de            = sqrt(val1*val2);
                double de2           = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
                double ncc_1, ncc_2, ncc_3;
                
                if( val1*val2 > 0)
                    ncc_1           = de2/de;
                else
                    ncc_1           = -1.0;
                
                if(Count_N[1] > 0)
                {
                    N                   = Count_N[1];
                    val1                = (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
                    val2                = (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
                    de                  = sqrt(val1*val2);
                    de2                 = (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
                    if( val1*val2 > 0)
                        ncc_2         = de2/de;
                    else
                        ncc_2           = -1.0;
                }
                
                if(Count_N[2] > 0)
                {
                    N                   = Count_N[2];
                    val1                = (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
                    val2                = (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
                    de                  = sqrt(val1*val2);
                    de2                 = (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
                    if( val1*val2 > 0)
                        ncc_3         = de2/de;
                    else
                        ncc_3           = -1.0;
                }
                
                double temp_rho;
                if(Count_N[1] > 0 && Count_N[2] > 0)
                    temp_rho      = ((ncc_1 + ncc_2 + ncc_3)/3.0);
                else if(Count_N[1] > 0)
                    temp_rho      = ((ncc_1 + ncc_2)/2.0);
                else if(Count_N[2] > 0)
                    temp_rho      = ((ncc_1 + ncc_3)/2.0);
                else
                    temp_rho        = ncc_1;
                
                long grid_index           = (mask_row+1)*3 + (mask_col+1);
                if(grid_index < 9)
                    result_rho[grid_index] = temp_rho;
                cell_count++;
            }
        }
    }
    
    if(cell_count == 9)
    {
        for(i=0;i<6;i++)
        {
            for(j=0;j<1;j++)
            {
                double sum = 0.0;
                for(k=0;k<9;k++)
                    sum += TsubA[i][k]*result_rho[k*1 + j];
                ATLT[i*1 + j] = sum;
            }
        }
        
        for(i=0;i<6;i++)
        {
            for(j=0;j<1;j++)
            {
                double sum = 0.0;
                for(k=0;k<6;k++)
                    sum += InverseSubA[i][k]*ATLT[k*1 + j];
                XX[i*1 + j] = sum;
            }
        }
        
        double demnum      = -pow(XX[4],2.0) + 4*XX[3]*XX[5];
        if(demnum > 0 && XX[3] < 0)
        {
            double max_X = (- 2*XX[5]*XX[1] + XX[2]*XX[4])/demnum;
            double max_Y = (- 2*XX[2]*XX[3] + XX[1]*XX[4])/demnum;
            double max_roh =  XX[0]                + XX[1]*max_X           + XX[2]*max_Y
            + XX[3]*max_X*max_X + XX[4]*max_X*max_Y + XX[5]*max_Y*max_Y;
            
            bool find_index_1   = false;
            bool find_index_2   = false;
            bool find_index     = false;
            if(fabs(max_X) <= 1.0)
                find_index_1 = true;
            if(fabs(max_Y) <= 1.0)
                find_index_2 = true;
            if (Pyramid_step >= 2)
                find_index  = find_index_1 & find_index_2 & (max_roh > 0.80);
            else
                find_index  = find_index_1 & find_index_2 & (max_roh > 0.60);
            
            if(find_index)
            {
                *sum_weight_X  = max_X*max_roh;
                *sum_weight_Y  = max_Y*max_roh;
                *sum_max_roh   = max_roh;
                
                peak_pos->m_X = max_X;
                peak_pos->m_Y = max_Y;
                
                check_pt = true;
            }
        }
    }
    free(result_rho);
    free(ATLT);
    free(XX);
    
    if(!check_pt)
    {
        *sum_weight_X   = 0;
        *sum_weight_Y   = 0;
        *sum_max_roh    = 0;
    }
    
    return check_pt;
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

inline void Set6by6Matrix(double subA[][6], double TsubA[][9], double InverseSubA[][6])
{
    for(int ii=0;ii<9;ii++)
        subA[ii][0]   = 1.0;
    
    subA[0][1] = -1.0; subA[0][2] = -1.0; subA[0][3] =  1.0; subA[0][4] =  1.0; subA[0][5] =  1.0;
    subA[1][1] =  0.0; subA[1][2] = -1.0; subA[1][3] =  0.0; subA[1][4] =  0.0; subA[1][5] =  1.0;
    subA[2][1] =  1.0; subA[2][2] = -1.0; subA[2][3] =  1.0; subA[2][4] = -1.0; subA[2][5] =  1.0;
    subA[3][1] = -1.0; subA[3][2] =  0.0; subA[3][3] =  1.0; subA[3][4] =  0.0; subA[3][5] =  0.0;
    subA[4][1] =  0.0; subA[4][2] =  0.0; subA[4][3] =  0.0; subA[4][4] =  0.0; subA[4][5] =  0.0;
    subA[5][1] =  1.0; subA[5][2] =  0.0; subA[5][3] =  1.0; subA[5][4] =  0.0; subA[5][5] =  0.0;
    subA[6][1] = -1.0; subA[6][2] =  1.0; subA[6][3] =  1.0; subA[6][4] = -1.0; subA[6][5] =  1.0;
    subA[7][1] =  0.0; subA[7][2] =  1.0; subA[7][3] =  0.0; subA[7][4] =  0.0; subA[7][5] =  1.0;
    subA[8][1] =  1.0; subA[8][2] =  1.0; subA[8][3] =  1.0; subA[8][4] =  1.0; subA[8][5] =  1.0;
    
    for(int ii=0;ii<6;ii++)
        for(int kk=0;kk<9;kk++)
            TsubA[ii][kk]       = subA[kk][ii];
    
    InverseSubA[0][0] =  0.555556; InverseSubA[0][1] =  0.000000; InverseSubA[0][2] =  0.000000; InverseSubA[0][3] = -0.333333; InverseSubA[0][4] =  0.000000; InverseSubA[0][5] = -0.333333;
    InverseSubA[1][0] =  0.000000; InverseSubA[1][1] =  0.166667; InverseSubA[1][2] =  0.000000; InverseSubA[1][3] =  0.000000; InverseSubA[1][4] =  0.000000; InverseSubA[1][5] =  0.000000;
    InverseSubA[2][0] =  0.000000; InverseSubA[2][1] =  0.000000; InverseSubA[2][2] =  0.166667; InverseSubA[2][3] =  0.000000; InverseSubA[2][4] =  0.000000; InverseSubA[2][5] =  0.000000;
    InverseSubA[3][0] = -0.333333; InverseSubA[3][1] =  0.000000; InverseSubA[3][2] =  0.000000; InverseSubA[3][3] =  0.500000; InverseSubA[3][4] =  0.000000; InverseSubA[3][5] =  0.000000;
    InverseSubA[4][0] =  0.000000; InverseSubA[4][1] =  0.000000; InverseSubA[4][2] =  0.000000; InverseSubA[4][3] =  0.000000; InverseSubA[4][4] =  0.250000; InverseSubA[4][5] =  0.000000;
    InverseSubA[5][0] = -0.333333; InverseSubA[5][1] =  0.000000; InverseSubA[5][2] =  0.000000; InverseSubA[5][3] =  0.000000; InverseSubA[5][4] =  0.000000; InverseSubA[5][5] =  0.500000;
}

#endif /* Template_h */
