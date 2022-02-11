//
//  LSF.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "LSF.hpp"

//LSF smoothing
void LSFSmoothing_DEM(const char *savepath, const char* outputpath, const double MPP, const int divide)
{
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    total_ST = time(0);
    
    char str_DEMfile[1000];
    char str_smooth_file[500];
    char DEM_header[500];
    char DEM_GEOTIFF_filename[500];
    char DEM_GEOTIFF_filename_slope[500];
    char metafilename[500];
    char str_matchfile[1000];
    char str_matchfile_tif[1000];
    char result_file[500];
    
    if(divide == 0)
    {
        sprintf(str_DEMfile, "%s/%s_dem.tif", savepath,outputpath);
        sprintf(str_smooth_file,"%s/%s_smooth.raw",savepath,outputpath);
        sprintf(DEM_header, "%s/%s_smooth.hdr", savepath,outputpath);
        sprintf(DEM_GEOTIFF_filename, "%s/%s_dem_smooth.tif", savepath, outputpath);
        sprintf(DEM_GEOTIFF_filename_slope, "%s/%s_dem_smooth_slope.tif", savepath, outputpath);
        sprintf(metafilename,"%s/%s_meta.txt",savepath,outputpath);
        sprintf(str_matchfile,"%s/%s_matchtag.raw",savepath,outputpath);
        sprintf(str_matchfile_tif,"%s/%s_matchtag.tif",savepath,outputpath);
        sprintf(result_file,"%s/%s_smooth_result.txt",savepath,outputpath);
    }
    else
    {
        sprintf(str_DEMfile, "%s/%s_%d_dem.tif", savepath,outputpath,divide);
        sprintf(str_smooth_file,"%s/%s_%d_smooth.raw",savepath,outputpath,divide);
        sprintf(DEM_header, "%s/%s_%d_smooth.hdr", savepath,outputpath,divide);
        sprintf(DEM_GEOTIFF_filename, "%s/%s_%d_dem_smooth.tif", savepath, outputpath,divide);
        sprintf(metafilename,"%s/%s_meta.txt",savepath,outputpath);
        sprintf(str_matchfile,"%s/%s_%d_matchtag.raw",savepath,outputpath,divide);
        sprintf(str_matchfile_tif,"%s/%s_%d_matchtag.tif",savepath,outputpath,divide);
        sprintf(result_file,"%s/%s_%d_smooth_result.txt",savepath,outputpath,divide);
    }
    
    printf("dem file %s\n",str_DEMfile);
    printf("metafilename %s\n",metafilename);
    printf("matchfile %s\n",str_matchfile);
    
    printf("smooth DEM %s\n",str_smooth_file);
    printf("DEM_header %s\n",DEM_header);
    printf("DEM_GEOTIFF %s\n", DEM_GEOTIFF_filename);
    printf("result file %s\n",result_file);
    
    FILE *pFile_DEM = fopen(str_DEMfile,"r");
    printf("check exist %s %d\n",str_DEMfile,!!pFile_DEM);
    
    if(pFile_DEM)
    {
        printf("MPP_stereo_angle %f\n",MPP);
        
        double minX, maxY;
        double grid_size;
        CSize DEM_size;
        TransParam param;
        param.bHemisphere = 3; //no assigned
        
        DEM_size = GetDEMsize(str_DEMfile,metafilename,&param,&grid_size,&minX,&maxY);
        printf("DEM_size %d\t%d\n",DEM_size.width,DEM_size.height);
        
        float *seeddem = GetDEMValue(str_DEMfile,DEM_size);
        printf("Done seeddem\n");
        
        long DEM_length = (long)(DEM_size.height)*(long)(DEM_size.width);
        printf("Done matchtag memory allocation %ld\n",DEM_length);
        printf("%d\n",DEM_size.width);
        printf("%d\n",DEM_size.height);
        printf("%f\n",grid_size);
        
        LSFINFO *Grid_info = (LSFINFO*)calloc(sizeof(LSFINFO),DEM_length);
        
        for(long count_index = 0 ; count_index < DEM_length; count_index++)
            Grid_info[count_index].lsf_kernel = 2;
        
        float *smooth_DEM = (float*)calloc(sizeof(float),DEM_length);
        
        double max_std = -100000;
        int max_std_iter = -1;
        int min_std_iter = 100;
        double min_std = 100000;
        
        double sigma_avg = 10000;
        double sigma_std = 10000;
        const int max_iter_count = 5;
        int s_iter = 0;
        bool check_smooth_iter = true;
        
        while(check_smooth_iter && s_iter < max_iter_count)
        {
            printf("start LSF\n");
            if((sigma_avg < 0.5 && sigma_std < 1) || s_iter == max_iter_count-1)
            {
                if(s_iter > 2)
                {
                    printf("final local suface fitting\n");
                    check_smooth_iter = false;
                }
            }
           
            DEM_STDKenel_LSF(Grid_info, &sigma_avg,&sigma_std, seeddem,smooth_DEM, grid_size, s_iter, DEM_size,MPP);
            
            if(sigma_avg > max_std)
            {
                max_std = sigma_avg;
                max_std_iter = s_iter;
            }
            
            if(sigma_avg < min_std)
            {
                min_std = sigma_avg;
                min_std_iter = s_iter;
            }
            
            printf("End LSF %d\tsigma avg std %f\t%f\n",s_iter,sigma_avg,sigma_std);
            memcpy(seeddem,smooth_DEM,sizeof(float)*DEM_length);
            
            s_iter++;
        }
        free(smooth_DEM);
        free(Grid_info);
        
        WriteGeotiff(DEM_GEOTIFF_filename, seeddem, DEM_size.width, DEM_size.height, grid_size, minX, maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
        
        
        float *slope = (float*)malloc(sizeof(float)*DEM_size.width*DEM_size.height);
        MakeSlopeImage(DEM_size,seeddem,slope,grid_size);
        WriteGeotiff(DEM_GEOTIFF_filename_slope, slope, DEM_size.width, DEM_size.height, grid_size, minX, maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
        free(slope);
        
        FILE* presult = fopen(result_file,"w");
        fprintf(presult,"%d\t%f\t%d\t%f\n",max_std_iter,max_std,min_std_iter,min_std);
        fclose(presult);
        
        printf("%d\t%f\t%d\t%f\n",max_std_iter,max_std,min_std_iter,min_std);
        
        free(seeddem);
        
        fclose(pFile_DEM);
    }
    
    total_ET = time(0);
    
    total_gap = difftime(total_ET,total_ST);
    printf("LSF processing time[min] = %5.2f\n",total_gap/60.0);
}

CSize GetDEMsize(char *GIMP_path, char* metafilename,TransParam* param, double *grid_size, double* _minX, double* _maxY)
{
    CSize seeddem_size;
    
    int check_ftype = 1; // 1 = tif, 2 = raw
    char *ext = strrchr(GIMP_path,'.');
    
    if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
        check_ftype = 1;
    else if(!strcmp("raw",ext+1))
        check_ftype = 2;
    
    FILE *pFile_meta  = fopen(metafilename,"r");
    printf("meta file = %s\n",metafilename);
    printf("DEM file = %s\n",GIMP_path);
    bool check_open_header = false;
    
    double minX, maxX, minY,maxY;
    if(check_ftype == 1)
    {
        char *tmp_str = remove_ext(GIMP_path);
        char *hdr_path = (char*)calloc(strlen(tmp_str) + 1 + 4, 1); //null term + extlen
        sprintf(hdr_path,"%s.tfw",tmp_str);
        
        FILE *pfile = fopen(hdr_path,"r");
        if(pfile)
        {
            char dir[500];
            
            printf("tfw path %s \n",hdr_path);
            
            TFW_reader_LSFDEM(hdr_path, &minX, &maxY, grid_size, &param->utm_zone, dir);
            
            if(param->projection == 2)
            {
                sprintf(param->direction,"%s",dir);
                printf("projection %d\t%d\t%s\n",param->projection,param->utm_zone,param->direction);
            }
            
            GetImageSize(GIMP_path,&seeddem_size);
            fclose(pfile);
            check_open_header = true;
        }
        else
        {
            printf("geotiff info %s\n",GIMP_path);
            seeddem_size = ReadGeotiff_info(GIMP_path, &minX, &maxY, grid_size);
            SetTranParam_fromGeoTiff(param,GIMP_path);
            
            check_open_header = true;
        }

        free(tmp_str);
        free(hdr_path);
    }
    else if(check_ftype == 2)
    {
        char *tmp_str = remove_ext(GIMP_path);
        char *hdr_path = (char*)calloc(strlen(tmp_str) + 1 + 4, 1); //null term + extlen
        sprintf(hdr_path,"%s.hdr",tmp_str);
        
        printf("hdr path %s\n",hdr_path);
        FILE *phdr = fopen(hdr_path,"r");
        if(phdr)
        {
            seeddem_size  = Envihdr_reader_seedDEM(*param,hdr_path, &minX, &maxY, grid_size);
            fclose(phdr);
            check_open_header = true;
        }

        free(hdr_path);
        free(tmp_str);
    }
    
    if(pFile_meta && !check_open_header)
    {
        char bufstr[500];
        printf("open Boundary\n");
        while(!feof(pFile_meta))
        {
            fgets(bufstr,500,pFile_meta);
            if (strstr(bufstr,"Output Resolution=")!=NULL)
            {
                printf("%s\n",bufstr);
                sscanf(bufstr,"Output Resolution=%lf\n",grid_size);
            }
            else if (strstr(bufstr,"Output Projection=")!=NULL)
            {
                printf("%s\n",bufstr);
                char projection[500];
                sscanf(bufstr,"Output Projection='+proj=%s",projection);
                if(!strcmp(projection,"stere"))
                {
                    double lat;
                    char t_str1[500];
                    param->projection = 1;
                    sscanf(bufstr,"Output Projection='+proj=stere +lat_0=%lf +lat_ts=%s",&lat,t_str1);
                    if(lat >= 0)
                        param->bHemisphere = true;
                    else
                        param->bHemisphere = false;
                    
                    printf("projection %d\tlat %f\themisphere %d\n",param->projection,lat,param->bHemisphere);
                }
                else
                {
                    double lat;
                    char t_str1[500];
                    char hh[100];
                    char com[100] = "M";
                    
                    int t_int;
                    param->projection = 2;
                    sscanf(bufstr,"Output Projection='+proj=utm +zone=%d +north=%s +datum=%s",&t_int,hh,t_str1);
                    param->utm_zone = t_int;
                    if(strcmp(hh,com) > 0)
                        param->bHemisphere = true;
                    else
                        param->bHemisphere = false;
                    
                    printf("projection %d\themisphere %d\tzone %d\n",param->projection,param->bHemisphere,param->utm_zone);
                }
            }
            else if (strstr(bufstr,"Output dimensions=")!=NULL)
            {
                printf("%s\n",bufstr);
                sscanf(bufstr,"Output dimensions=%d\t%d\n",&seeddem_size.width,&seeddem_size.height);
            }
            else if (strstr(bufstr,"Upper left coordinates=")!=NULL)
            {
                printf("%s\n",bufstr);
                sscanf(bufstr,"Upper left coordinates=%lf\t%lf\n",&minX,&maxY);
            }
        }
        check_open_header = true;
        fclose(pFile_meta);
    }
    
    if(!check_open_header)
    {
        printf("No information about input tif!!\n");
        exit(1);
    }
    
    maxX    = minX + (*grid_size)*((double)seeddem_size.width);
    minY    = maxY - (*grid_size)*((double)seeddem_size.height);
    
    printf("%d\n",seeddem_size.width);
    printf("%d\n",seeddem_size.height);
    printf("%f\n",minX);
    printf("%f\n",minY);
    printf("%f\n",maxX);
    printf("%f\n",maxY);
    printf("%f\n",*grid_size);
    
    *_minX = minX;
    *_maxY = maxY;
    
    return seeddem_size;
}

void DEM_STDKenel_LSF(LSFINFO *Grid_info, double* sigma_average,double* sigma_std, float *seeddem, float *smooth_DEM, const double grid_size, const int smooth_iteration,const CSize seeddem_size, const double MPP_stereo_angle)
{
    long total_selected_points = 0;
    double sigma_sum = 0;
    double sigma2_sum = 0;
    long data_length = (long)seeddem_size.width*(long)seeddem_size.height;
    
#pragma omp parallel for schedule(guided) reduction(+:sigma_sum, total_selected_points, sigma2_sum)
    for(long iter_count = 0 ; iter_count < data_length ; iter_count++)
    {
        long pts_row = floor(iter_count/seeddem_size.width);
        long pts_col = iter_count % seeddem_size.width;
        double fitted_Z = seeddem[iter_count];
        double sigma;

        if(pts_col >= 0 && pts_col < seeddem_size.width && pts_row >= 0 && pts_row < seeddem_size.height)
        {
            if(seeddem[iter_count] > -50 )
            {
                long selected_count = 0;
                sigma = LocalSurfaceFitting_DEM(Grid_info, seeddem, selected_count, &fitted_Z, MPP_stereo_angle, smooth_iteration, seeddem_size.height, seeddem_size.width, grid_size, pts_col, pts_row);

                if(sigma < 20 && sigma > 0 && selected_count > 6)
                {
                    smooth_DEM[iter_count] = fitted_Z;
                    total_selected_points++;
                    sigma_sum += sigma;
                    sigma2_sum += (sigma*sigma);
                }
                else
                    smooth_DEM[iter_count] = seeddem[iter_count];
            }
            else
                smooth_DEM[iter_count] = seeddem[iter_count];
        }
    }
    printf("sigma %f\t%ld\n",sigma_sum,total_selected_points);
    
    *sigma_average = sigma_sum/(double)total_selected_points;
    *sigma_std = sqrt( sigma2_sum/(double)total_selected_points - (*sigma_average)*(*sigma_average) );
    printf("avg sigma %f\tstd sigma %f\ttotal_pts %ld\n",*sigma_average,*sigma_std,total_selected_points);
}

double LocalSurfaceFitting_DEM(LSFINFO *Grid_info, float *input, long &numpts, double *fitted_Z, const double MPP, const int smooth_iter, const long row_size, const long col_size, const double grid, const long X, const long Y)
{
    double sigma = 999999;
    
    long row,col;
    long final_interval = 5;
    
    const long data_length = row_size*col_size;
    const long t_index = Y*col_size + X;
    
    const double MPP_th = 5;
    long int mask_interval = 1;

    long int add_interval = 0;
    if(MPP > MPP_th)
        add_interval = 2;
    else if(MPP > MPP_th*2)
        add_interval = 3;

    if(smooth_iter > 0)
    {
        long pre_final_interval = Grid_info[t_index].lsf_kernel;

        if(add_interval > 0)
        {
            mask_interval = floor(pre_final_interval/5.0);
            if(mask_interval <= 0)
                mask_interval = 1;
        }
        else
        {
            if(grid < 1)
                mask_interval = 2;
        }
        final_interval = (long int)Grid_info[t_index].lsf_kernel;
    }
    else
    {
        int max_pts = 9;
        if(grid >= 8)
            max_pts = 7;
        
        long count1,count2,count3,count4;
        const long row_interval = 15;
        long interval = 2;
        
        bool check_stop = false;
        while (!check_stop)
        {
            numpts = 0;
            count1 = 0;
            count2 = 0;
            count3 = 0;
            count4 = 0;
            for(row = -interval;row <= interval;row++)
            {
                for(col = -interval;col <= interval ; col++)
                {
                    long int grid_pos = (long int)((Y+row)*(long int)col_size + (X+col));
                    if(grid_pos >= 0 && grid_pos < data_length && Y+row >= 0 && Y+row < row_size && X+col >= 0 && X+col < col_size)
                    {

                        if(input[grid_pos] > -100)
                        {
                            if(row >= 0 && row <=  interval && col >= 0 && col <=  interval)
                            {
                                count1++;
                                numpts++;
                            }
                            else if (row >= 0 && row <=  interval && col <= 0 && col >= -interval)
                            {
                                count2++;
                                numpts++;
                            }
                            else if (row < 0 && row >= -interval && col < 0 && col >= -interval)
                            {
                                count3++;
                                numpts++;
                            }
                            else if (row < 0 && row >= -interval && col > 0 && col <=  interval)
                            {
                                count4++;
                                numpts++;
                            }
                        }
                    }
                }
            }

            if (interval >= row_interval || (numpts > max_pts && count1 > 2 && count2 > 2 && count3 > 2 && count4 > 2))
            {
                check_stop = true;
                final_interval = interval;
                Grid_info[t_index].lsf_kernel = (unsigned char)final_interval;
            }
            else
                interval = interval + 1;
        }
    }

    double maxX_ptslists = -10000000;
    double maxY_ptslists = -10000000;
    double minX_ptslists =  10000000;
    double minY_ptslists =  10000000;
    double distX_ptslists = 0;
    double distY_ptslists = 0;
    const double Scale_ptslists = 1000;

    vector<D3DPOINT> XY_save;
    XY_save.clear();
    long count = 0;
    for(row = -final_interval;row <= final_interval;row+= mask_interval)
    {
        for(col = -final_interval;col <= final_interval ; col+= mask_interval)
        {
            long int grid_pos = (long int)((Y+row)*(long int)col_size + (X+col));
            long int grid_pos_col = (long int)((X+col));
            long int grid_pos_row = (long int)((Y+row));

            if(grid_pos >= 0 && grid_pos < data_length &&
                    grid_pos_row >= 0 && grid_pos_row < row_size && grid_pos_col >= 0 && grid_pos_col < col_size)
            {
                if(input[grid_pos] > - 100)
                {
                    D3DPOINT temp;
                    temp.m_X = grid_pos_col*grid;
                    temp.m_Y = grid_pos_row*grid;
                    temp.m_Z = input[grid_pos];
                    temp.flag = true;

                    XY_save.push_back(temp);

                    count++;
                    if(maxX_ptslists < grid_pos_col*grid)
                        maxX_ptslists = grid_pos_col*grid;
                    if(maxY_ptslists < grid_pos_row*grid)
                        maxY_ptslists = grid_pos_row*grid;
                    if(minX_ptslists > grid_pos_col*grid)
                        minX_ptslists = grid_pos_col*grid;
                    if(minY_ptslists > grid_pos_row*grid)
                        minY_ptslists = grid_pos_row*grid;
                }
            }
        }
    }

    distX_ptslists = maxX_ptslists - minX_ptslists;
    distY_ptslists = maxY_ptslists - minY_ptslists;

    double scale_factor_X = 1.0/distX_ptslists*Scale_ptslists;
    double scale_factor_Y = 1.0/distY_ptslists*Scale_ptslists;
    double X_scaled = ((double)X*grid - minX_ptslists)*scale_factor_X;
    double Y_scaled = ((double)Y*grid - minY_ptslists)*scale_factor_Y;
    double X_plane = ((double)X*grid - minX_ptslists);
    double Y_plane = ((double)Y*grid - minY_ptslists);
    numpts = 0;
    if(count > 15)
    {
        GMA_double *A_matrix = GMA_double_create(XY_save.size(), 3);
        GMA_double *L_matrix = GMA_double_create(XY_save.size(), 1);
        GMA_double *AT_matrix = GMA_double_create(3,XY_save.size());
        GMA_double *ATA_matrix = GMA_double_create(3,3);

        GMA_double *ATAI_matrix = GMA_double_create(3,3);
        GMA_double *ATL_matrix = GMA_double_create(3,1);

        GMA_double *X_matrix = GMA_double_create(3,1);
        GMA_double *AX_matrix = GMA_double_create(XY_save.size(),1);
        GMA_double *V_matrix = GMA_double_create(XY_save.size(),1);

        //plane fitting
        vector<D3DPOINT>::iterator it;
        count = 0;
        for(it = XY_save.begin() ; it != XY_save.end() ; ++it)
        {
            A_matrix->val[count][0] = it->m_X-minX_ptslists;
            A_matrix->val[count][1] = it->m_Y-minY_ptslists;
            A_matrix->val[count][2] = 1.0;

            L_matrix->val[count][0] = it->m_Z;
            count++;
        }

        GMA_double_Tran(A_matrix,AT_matrix);
        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
        GMA_double_inv(ATA_matrix,ATAI_matrix);
        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
        GMA_double_sub(AX_matrix,L_matrix,V_matrix);

        double sum_V = 0;
        for(row = 0; row < XY_save.size() ; row++)
            sum_V += V_matrix->val[row][0];

        if(!std::isnan(sum_V) && !std::isnan(X_matrix->val[0][0]) && !std::isnan(X_matrix->val[1][0]) && !std::isnan(X_matrix->val[2][0]))
        {
            double plane_Z = X_plane*X_matrix->val[0][0] + Y_plane*X_matrix->val[1][0] + X_matrix->val[2][0];
            if(plane_Z > -100 && plane_Z < 15000)
            {
                D3DPOINT N(X_matrix->val[0][0], X_matrix->val[1][0], 1.0 , 0);
                double norm  = SQRT(N);
                double angle = acos(fabs(N.m_Z)/norm)*RadToDeg;
                SetAngle(angle);

                long hist[20] = {0};
                for(row = 0; row < XY_save.size() ; row++)
                {
                    int hist_index = (int)(std::abs(V_matrix->val[row][0]));
                    if(hist_index > 19)
                        hist_index = 19;
                    if(hist_index >= 0 && hist_index <= 19)
                        hist[hist_index]++;
                }

                double hist_th = 0.8;
                if(grid >= 8)
                    hist_th = 0.9;
                else
                    hist_th = 0.8;
                
                int V_th = 20;
                int hist_sum = 0;
                double hist_rate;
                bool check_V = true;
                row = 0;
                while(check_V && row < 20)
                {
                    hist_sum += hist[row];
                    hist_rate = (double)hist_sum/(double)XY_save.size();
                    if(hist_rate > hist_th)
                    {
                        V_th = row;
                        check_V = false;
                    }
                    row++;
                }

                count = 0;
                vector<D3DPOINT> XY_selected;
                for(it = XY_save.begin() ; it != XY_save.end() ; ++it)
                {
                    if(std::abs(V_matrix->val[count][0]) < V_th+1)
                    {
                        D3DPOINT temp((it->m_X-minX_ptslists)*scale_factor_X, (it->m_Y-minY_ptslists)*scale_factor_Y, it->m_Z, it->flag);
                        XY_selected.push_back(temp);
                    }
                    count++;
                }

                GMA_double_destroy(A_matrix);
                GMA_double_destroy(L_matrix);
                GMA_double_destroy(AT_matrix);
                GMA_double_destroy(ATA_matrix);
                GMA_double_destroy(ATAI_matrix);
                GMA_double_destroy(ATL_matrix);
                GMA_double_destroy(X_matrix);
                GMA_double_destroy(AX_matrix);
                GMA_double_destroy(V_matrix);

                long selected_count = XY_selected.size();

                if(selected_count > 15)
                {
                    A_matrix = GMA_double_create(selected_count, 6);
                    L_matrix = GMA_double_create(selected_count, 1);
                    AT_matrix = GMA_double_create(6,selected_count);
                    ATA_matrix = GMA_double_create(6,6);

                    ATAI_matrix = GMA_double_create(6,6);
                    ATL_matrix = GMA_double_create(6,1);

                    X_matrix = GMA_double_create(6,1);
                    AX_matrix = GMA_double_create(selected_count,1);
                    V_matrix = GMA_double_create(selected_count,1);

                    count = 0;
                    for(it = XY_selected.begin() ; it != XY_selected.end() ; ++it)
                    {
                        A_matrix->val[count][0] = it->m_X*it->m_X;
                        A_matrix->val[count][1] = it->m_X*it->m_Y;
                        A_matrix->val[count][2] = it->m_Y*it->m_Y;
                        A_matrix->val[count][3] = it->m_X;
                        A_matrix->val[count][4] = it->m_Y;
                        A_matrix->val[count][5] = 1.0;

                        L_matrix->val[count][0] = it->m_Z;
                        count++;
                    }

                    GMA_double_Tran(A_matrix,AT_matrix);
                    GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
                    GMA_double_inv(ATA_matrix,ATAI_matrix);
                    GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
                    GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
                    GMA_double_mul(A_matrix,X_matrix,AX_matrix);
                    GMA_double_sub(AX_matrix,L_matrix,V_matrix);

                    double sum = 0;
                    count = 0;
                    vector<D3DPOINT>::iterator it_sel;
                    for(it_sel = XY_selected.begin() ; it_sel != XY_selected.end() ; ++it_sel)
                    {
                        sum += V_matrix->val[count][0] * V_matrix->val[count][0];

                        double temp_fitted_Z = X_matrix->val[0][0]*it_sel->m_X*it_sel->m_X + X_matrix->val[1][0]*it_sel->m_X*it_sel->m_Y +
                            X_matrix->val[2][0]*it_sel->m_Y*it_sel->m_Y + X_matrix->val[3][0]*it_sel->m_X + X_matrix->val[4][0]*it_sel->m_Y + X_matrix->val[5][0];

                        count++;
                    }

                    const double p = 1.5;
                    if(!std::isnan(sum) && sum > 0 )
                    {
                        sigma = sqrt(sum/(double)selected_count);

                        double A = X_matrix->val[0][0];
                        double B = X_matrix->val[2][0];
                        double C = X_matrix->val[3][0];
                        double D = X_matrix->val[4][0];
                        double E = X_matrix->val[1][0];

                        double det = 4*A*B - E*E;
                        double det1 = D*E - 2*C*B;
                        double det2 = 2*A*D - C*E;

                        bool check_clinder = false;
                        if(det == 0 && det1 == det2)
                            check_clinder = true;

                        if(!check_clinder)
                        {
                            *fitted_Z = X_matrix->val[0][0]*X_scaled*X_scaled + X_matrix->val[1][0]*X_scaled*Y_scaled + X_matrix->val[2][0]*Y_scaled*Y_scaled +X_matrix->val[3][0]*X_scaled + X_matrix->val[4][0]*Y_scaled + X_matrix->val[5][0];
                        }
                        else
                        {
                            double sum_weight = 0;
                            double sum_weigtdist = 0;
                            for(it_sel = XY_selected.begin() ; it_sel != XY_selected.end() ; ++it_sel)
                            {
                                double dist = sqrt((it_sel->m_X - X_scaled)*(it_sel->m_X - X_scaled) + (it_sel->m_Y - Y_scaled)*(it_sel->m_Y - Y_scaled));
                                sum_weight += (1.0/pow(dist,p));
                                sum_weigtdist += (it_sel->m_Z/pow(dist,p));
                            }

                            if(sum_weight > 0)
                            {
                                *fitted_Z = sum_weigtdist/sum_weight;
                                sigma = 1;
                            }
                            else
                                sigma = 999999;
                        }
                    }
                    else
                    {
                        double sum_weight = 0;
                        double sum_weigtdist = 0;
                        for(it_sel = XY_selected.begin() ; it_sel != XY_selected.end() ; ++it_sel)
                        {
                            double dist = sqrt((it_sel->m_X - X_scaled)*(it_sel->m_X - X_scaled) + (it_sel->m_Y - Y_scaled)*(it_sel->m_Y - Y_scaled));
                            sum_weight += (1.0/pow(dist,p));
                            sum_weigtdist += (it_sel->m_Z/pow(dist,p));
                        }

                        if(sum_weight > 0)
                        {
                            *fitted_Z = sum_weigtdist/sum_weight;
                            sigma = 1;
                        }
                        else
                            sigma = 999999;
                    }

                    if(grid > 2)
                    {
                        if(angle < 10)
                            Grid_info[t_index].lsf_kernel = 4;
                        else if(angle < 20)
                            Grid_info[t_index].lsf_kernel = 3;
                        else if(angle < 30)
                            Grid_info[t_index].lsf_kernel = 2;
                        else if(Grid_info[t_index].lsf_kernel < 2)
                            Grid_info[t_index].lsf_kernel = 2;
                    }
                    else if(grid == 2)
                    {
                        if(angle < 10)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(5 + add_interval);
                        else if(angle < 20)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(4 + add_interval);
                        else if(angle < 30)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(3 + add_interval);
                        else if(Grid_info[t_index].lsf_kernel < 2 + add_interval)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(2 + add_interval);
                    }
                    else if(grid == 1)
                    {
                        if(angle < 10)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(7 + add_interval*2);
                        else if(angle < 20)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(6 + add_interval*2);
                        else if(angle < 30)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(5 + add_interval*2);
                        else if(Grid_info[t_index].lsf_kernel < 4 + add_interval*2)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(4 + add_interval*2);
                    }
                    else
                    {
                        if(angle < 10)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(9 + add_interval*3);
                        else if(angle < 20)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(8 + add_interval*3);
                        else if(angle < 30)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(7 + add_interval*3);
                        else if(Grid_info[t_index].lsf_kernel < 6 + add_interval*3)
                            Grid_info[t_index].lsf_kernel = (unsigned char)(6 + add_interval*3);
                    }

                    GMA_double_destroy(A_matrix);
                    GMA_double_destroy(L_matrix);
                    GMA_double_destroy(AT_matrix);
                    GMA_double_destroy(ATA_matrix);

                    GMA_double_destroy(ATAI_matrix);
                    GMA_double_destroy(ATL_matrix);

                    GMA_double_destroy(X_matrix);
                    GMA_double_destroy(AX_matrix);
                    GMA_double_destroy(V_matrix);

                    numpts = selected_count;
                }
                else
                {
                    sigma = 999999;
                    numpts = 0;
                }
                XY_selected.clear();
                vector<D3DPOINT>().swap(XY_selected);
            }
            else
            {
                GMA_double_destroy(A_matrix);
                GMA_double_destroy(L_matrix);
                GMA_double_destroy(AT_matrix);
                GMA_double_destroy(ATA_matrix);
                GMA_double_destroy(ATAI_matrix);
                GMA_double_destroy(ATL_matrix);
                GMA_double_destroy(X_matrix);
                GMA_double_destroy(AX_matrix);
                GMA_double_destroy(V_matrix);

                sigma = 999999;
                numpts = 0;
            }
        }
    }
    else
    {
        sigma = 999999;
    }

    XY_save.clear();
    vector<D3DPOINT>().swap(XY_save);

    return sigma;
}


