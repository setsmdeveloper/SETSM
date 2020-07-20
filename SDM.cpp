//
//  SDM.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "SDM.hpp"

const char setsm_version[] = "4.2.0";

//SDM ortho
bool SDM_ortho(char* _filename, ARGINFO args, double** Coreg_param)
{
    ProInfo proinfo;
    if(OpenProject(_filename,&proinfo,args))
    {
        bool check_cal = false;
        char check_file[500];
        sprintf(check_file, "%s/%s_dmag.tif", proinfo.save_filepath, proinfo.Outputpath_name);
        FILE *pcheckFile = fopen(check_file,"r");
        if(!pcheckFile)
            check_cal = true;
        
        if(check_cal)
        {
            if(Maketmpfolders(&proinfo))
            {
                TransParam param;
                SetTranParam_fromGeoTiff(&param,proinfo.Imagefilename[0]);
                
                double Rimageparam[2] = {0.0};
                
                Rimageparam[0] = Coreg_param[1][0];
                Rimageparam[1] = Coreg_param[1][1];
                
                proinfo.SDM_SS = args.SDM_SS;
                proinfo.SDM_days = args.SDM_days;
                proinfo.SDM_AS = args.SDM_AS;
                
                time_t current_time;
                char*    c_time_string;
                
                current_time = time(NULL);
                c_time_string = ctime(&current_time);
                
                char metafilename[500];
                sprintf(metafilename, "%s/%s_meta.txt", proinfo.save_filepath, proinfo.Outputpath_name);
                
                FILE *pMetafile    = fopen(metafilename,"w");
                fprintf(pMetafile,"SETSM Version=%s\n",setsm_version);
                
                char temp_filepath[500];
                double Image1_gsd, Image2_gsd;
                ImageGSD GSD_image1, GSD_image2;
                double image1_minX, image1_maxY, image2_minX, image2_maxY, image1_grid, image2_grid;
                
                CSize Limagesize = ReadGeotiff_info(proinfo.Imagefilename[0], &image1_minX, &image1_maxY, &image1_grid);
                CSize Rimagesize = ReadGeotiff_info(proinfo.Imagefilename[1], &image2_minX, &image2_maxY, &image2_grid);
                
                GSD_image1.row_GSD = image1_grid;
                GSD_image1.col_GSD = image1_grid;
                GSD_image1.pro_GSD = image1_grid;
                GSD_image2.row_GSD = image2_grid;
                GSD_image2.col_GSD = image2_grid;
                GSD_image2.pro_GSD = image2_grid;
                printf("Limagesize %d\t%d\t Rimagesize %d\t%d\n",Limagesize.width,Limagesize.height,Rimagesize.width,Rimagesize.height);
                
                double Boundary[4], LBoundary[4], RBoundary[4];
                LBoundary[0] = image1_minX;
                LBoundary[1] = image1_maxY - Limagesize.height*image1_grid;
                LBoundary[2] = image1_minX + Limagesize.width*image1_grid;
                LBoundary[3] = image1_maxY;
                
                RBoundary[0] = image2_minX;
                RBoundary[1] = image2_maxY - Rimagesize.height*image2_grid;
                RBoundary[2] = image2_minX + Rimagesize.width*image2_grid;
                RBoundary[3] = image2_maxY;
                
                for(int i = 0 ; i < 4 ; i++)
                {
                    proinfo.LBoundary[i] = LBoundary[i];
                    proinfo.RBoundary[i] = RBoundary[i];
                }
                
                if(args.check_boundary)
                {
                    Boundary[0] = args.Min_X;
                    Boundary[1] = args.Min_Y;
                    Boundary[2] = args.Max_X;
                    Boundary[3] = args.Max_Y;
                    
                    printf("boundary = %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                }
                else
                {
                    for(int i=0;i<4;i++)
                    {
                        if(i<2)
                            Boundary[i] = ceil((max(LBoundary[i], RBoundary[i]) / 2.0)) * 2;
                        else
                            Boundary[i] = floor((min(LBoundary[i], RBoundary[i]) / 2.0)) * 2;
                    }
                }
                
                printf("boundary = %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                
                proinfo.resolution = image1_grid;
                
                printf("image resolution %f\n",proinfo.resolution);
                
                int end_level = floor(log10(proinfo.DEM_resolution/(proinfo.resolution*15))/log10(2));
                int th_grid = proinfo.resolution*pwrtwo(end_level);
                
                if(end_level < 0)
                    end_level = 0;
                if(end_level > 4)
                    end_level = 4;
                
                proinfo.end_level = end_level;
                proinfo.pyramid_level = args.pyramid_level;
                proinfo.number_of_images = 2;
                
                printf("pyramid level %d\tSDM_ss %d\tend_level = %d\t%d\n",proinfo.pyramid_level,proinfo.SDM_SS,end_level,th_grid);
                
                int sdm_kernal_size = floor( (double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution*pwrtwo(proinfo.pyramid_level)));
            printf("sdm_kernel size %f\t%f\t%f\t%f\t%d\n",sdm_kernal_size,proinfo.SDM_AS,proinfo.SDM_days,proinfo.resolution,proinfo.pyramid_level);
                if(proinfo.pre_DEMtif)
                {
                    sdm_kernal_size = 3;
                    proinfo.SDM_SS = sdm_kernal_size;
                    proinfo.pyramid_level = args.pyramid_level;
                    proinfo.end_level = 0;
                }
                else
                {
                    if(sdm_kernal_size > 3)
                        proinfo.SDM_SS = sdm_kernal_size;
                    else
                        proinfo.SDM_SS = 3;
                    
                    printf("ks %d\t sdm_ss %d\n",sdm_kernal_size,proinfo.SDM_SS);
                    
                    if(proinfo.SDM_SS > 5 && proinfo.pyramid_level >= 2)
                    {
                        bool check_while = false;
                        while(!check_while)
                        {
                            proinfo.pyramid_level++;
                            sdm_kernal_size = floor( (double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution*pwrtwo(proinfo.pyramid_level)));
                            if(sdm_kernal_size > 3)
                                proinfo.SDM_SS = sdm_kernal_size;
                            else
                                proinfo.SDM_SS = 3;
                            
                            if(proinfo.SDM_SS < 5)
                                check_while = true;
                        }
                    }
                    
                    if(proinfo.pyramid_level <= proinfo.end_level)
                    {
                        if(proinfo.pyramid_level >= 1)
                            proinfo.end_level = proinfo.pyramid_level - 1;
                        else
                            proinfo.end_level = proinfo.pyramid_level;
                    }
                    
                    if(proinfo.DEM_resolution < 50)
                        proinfo.end_level = 0;
                }
                
                printf("pyramid leve %d\tend level %d\n",proinfo.pyramid_level,proinfo.end_level);
                
                printf("ks %d\t sdm_ss %d\n",sdm_kernal_size,proinfo.SDM_SS);
                
                if (!args.check_DEM_space)
                    proinfo.DEM_resolution = proinfo.resolution;
                
                if(!proinfo.check_checktiff && !args.check_ortho)
                {
                    fprintf(pMetafile,"Creation Date=%s",c_time_string);
                    fprintf(pMetafile,"Image 1=%s\n",proinfo.Imagefilename[0]);
                    fprintf(pMetafile,"Image 2=%s\n",proinfo.Imagefilename[1]);
                    fprintf(pMetafile,"Output Resolution=%f\n",proinfo.DEM_resolution);
                }
                
                GetImageSize(proinfo.Imagefilename[0],&Limagesize);    GetImageSize(proinfo.Imagefilename[1],&Rimagesize);
                proinfo.image_bits = ReadGeotiff_bits(proinfo.Imagefilename[0]);
                printf("Limagesize %d\t%d\t Rimagesize %d\t%d\n",Limagesize.width,Limagesize.height,Rimagesize.width,Rimagesize.height);
                
                const uint8 pyramid_step = proinfo.pyramid_level;
                const uint8 Template_size    = 15;
                proinfo.pyramid_level = pyramid_step;
                printf("proinfo res %f\n",proinfo.resolution);
                int matching_number = 0;
                Matching_SETSM_SDM(proinfo, param, Template_size, Rimageparam, Limagesize, Rimagesize, Boundary, GSD_image1, GSD_image2, &matching_number);
                /*
                int final_iteration = 3;
                if(matching_number > 10)
                {
                    printf("Tile merging start final iteration %d!!\n",final_iteration);
                    int buffer_tile = 0;
                    double mt_grid_size = MergeTiles_SDM(proinfo,1,1,buffer_tile,final_iteration,param,pyramid_step);
                }
                 */
                fclose(pMetafile);
            }
        }
        else
        {
            printf("SDM has been already processed. *_dmag.tif file exists!!\n");
        }
    }

    return 0;
}

void Matching_SETSM_SDM(ProInfo proinfo, TransParam param, uint8 Template_size, double *Rimageparam, const CSize Limagesize, const CSize Rimagesize, const double *Boundary, const ImageGSD gsd_image1, const ImageGSD gsd_image2, int *matching_number)
{
    bool check_cal = false;
    char check_file[500];
    sprintf(check_file,"%s/txt/tin_xshift_level_1_1_0_iter_3.txt",proinfo.save_filepath);
    FILE *pcheckFile = fopen(check_file,"r");
    if(!pcheckFile)
        check_cal = true;
    
    if(check_cal)
    {
        printf("start cal tile\n");
        
        char save_file[500], Lsubsetfilename[500], Rsubsetfilename[500];
        sprintf(save_file,"%s/txt/echo_result_row_%d_col_%d.txt",proinfo.save_filepath,1,1);
        FILE *fid            = fopen(save_file,"w");
        fprintf(fid,"RA param X = %f\tY = %f\n",Rimageparam[0],Rimageparam[1]);
        printf("Coreg param X = %f\tY = %f\n",Rimageparam[0],Rimageparam[1]);
        sprintf(save_file,"%s/txt/headerinfo_row_%d_col_%d.txt",proinfo.save_filepath,1,1);
        FILE *fid_header    = fopen(save_file,"w");
    
        double subBoundary[4];
        subBoundary[0] =  (int)(floor(Boundary[0]/8))*8;
        subBoundary[1] =  (int)(floor(Boundary[1]/8))*8;
        subBoundary[2] =  (int)(floor(Boundary[2]/8))*8;
        subBoundary[3] =  (int)(floor(Boundary[3]/8))*8;
        
        printf("subBoundary = %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]);
        
        char *filename = GetFileName(proinfo.Imagefilename[0]);
        filename = remove_ext(filename);
        sprintf(Lsubsetfilename,"%s/%s_subset_%d_%d.raw",proinfo.tmpdir,filename,1,1);
        free(filename);

        filename = GetFileName(proinfo.Imagefilename[1]);
        filename = remove_ext(filename);
        sprintf(Rsubsetfilename,"%s/%s_subset_%d_%d.raw",proinfo.tmpdir,filename,1,1);
        free(filename);
        
        printf("subsetimage level %d\n",proinfo.pyramid_level);
        
        LevelInfo levelinfo = {NULL};
        levelinfo.Boundary = subBoundary;
        levelinfo.Template_size = &Template_size;
        levelinfo.param = &param;
        double **ImageAdjust = (double**)calloc(sizeof(double*),proinfo.number_of_images);
        for(int ti = 0 ; ti < proinfo.number_of_images ; ti++)
            ImageAdjust[ti] = (double*)calloc(sizeof(double*),2);
      
        ImageAdjust[1][0] = Rimageparam[0];
        ImageAdjust[1][1] = Rimageparam[1];
        levelinfo.ImageAdjust = ImageAdjust;
        
        CSize subsetsize[2];
        D2DPOINT startpos_ori[2];
        uint16 **SourceImages = (uint16 **)malloc(sizeof(uint16*)*2);
        if(subsetImage_SDM(proinfo, subBoundary, startpos_ori, subsetsize, SourceImages))
        {
            printf("Completion of subsetImage!!\n");
            const int pyramid_step = proinfo.pyramid_level;
            
            if( subsetsize[0].height > Template_size*pwrtwo(pyramid_step-1) && subsetsize[0].width > Template_size*pwrtwo(pyramid_step-1) &&
               subsetsize[1].height > Template_size*pwrtwo(pyramid_step-1) && subsetsize[1].width > Template_size*pwrtwo(pyramid_step-1) )
            {
                double pre_grid_resolution = 0;
                
                int level              = pyramid_step;
                bool lower_level_match    = true;
                bool flag_start            = false;
                
                CSize Size_Grid2D, pre_Size_Grid2D;
                Size_Grid2D.height    = 0;
                Size_Grid2D.width    = 0;
                
                CSize **data_size_lr = (CSize**)malloc(sizeof(CSize*)*proinfo.number_of_images);
                uint16 ***SubImages     = (uint16***)malloc(sizeof(uint16**)*(level+1));
                uint16 ***SubMagImages  = (uint16***)malloc(sizeof(uint16**)*(level+1));
                uint8 ***SubOriImages   = (uint8***)malloc(sizeof(uint8**)*(level+1));
                
                for(int iter_level = 0 ; iter_level < level + 1; iter_level++)
                {
                    SubImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo.number_of_images);
                    SubOriImages[iter_level] = (uint8**)malloc(sizeof(uint8*)*proinfo.number_of_images);
                    SubMagImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo.number_of_images);
                    
                    for(int image_index = 0 ; image_index < proinfo.number_of_images ; image_index++)
                    {
                        data_size_lr[image_index] = (CSize*)malloc(sizeof(CSize)*(level+1));
                        SetPySizes(data_size_lr[image_index], subsetsize[image_index], level);
                        
                        long data_length = (long)data_size_lr[image_index][iter_level].height*(long)data_size_lr[image_index][iter_level].width;
                        
                        if(iter_level == 0)
                            SubImages[iter_level][image_index] = SourceImages[image_index];
                        else
                            SubImages[iter_level][image_index] = NULL;
                        
                        SubOriImages[iter_level][image_index] = (uint8*)malloc(sizeof(uint8)*data_length);
                        SubMagImages[iter_level][image_index] = (uint16*)malloc(sizeof(uint16)*data_length);
                    }
                }
                //pyramid image generation
                printf("Preprocessing start!!\n");
                SetPyramidImages(&proinfo, level+1, data_size_lr, SubImages, SubMagImages, SubOriImages);
                printf("Preprocessing finish(time[m] )!!\n");
                
                printf("SDM generation start!!\n");
                
                UGRIDSDM *GridPT3 = NULL, *Pre_GridPT3 = NULL;
                D2DPOINT *GridPT = NULL;
                int final_level_iteration = 1;
                
                int total_matching_candidate_pts = 0;
                double matching_rate = 0;
                
                while(lower_level_match && level >= 0)
                {
                    int prc_level = level;
                    if(prc_level <= proinfo.end_level)
                        prc_level = proinfo.end_level;
                    
                    if(level > 4)
                    {
                        if(level <= 6)
                            Template_size = 11;
                        else if(level == 7)
                            Template_size = 7;
                        else
                            Template_size = 7;
                    }
                    
                    printf("level = %d\t final_level_iteration %d\n",level,final_level_iteration);
                    
                    D2DPOINT Startpos[2];
                    for(int image_index = 0 ; image_index < proinfo.number_of_images ; image_index++)
                    {
                        Startpos[image_index].m_X        = (double)(startpos_ori[image_index].m_X/pwrtwo(prc_level));
                        Startpos[image_index].m_Y        = (double)(startpos_ori[image_index].m_Y/pwrtwo(prc_level));
                    }
                    
                    levelinfo.py_Images = SubImages[level];
                    levelinfo.py_MagImages = SubMagImages[level];
                    levelinfo.py_OriImages = SubOriImages[level];
                    
                    levelinfo.py_Startpos = Startpos;
                    levelinfo.Pyramid_step = &prc_level;
                    levelinfo.py_Sizes = data_size_lr;
                    
                    double Th_roh, Th_roh_min, Th_roh_start, Th_roh_next;
                    SetThs_SDM(level,final_level_iteration, &Th_roh, &Th_roh_min, &Th_roh_next, &Th_roh_start,pyramid_step);
                    
                    double py_resolution = 0;
                    double grid_resolution = 0;
                    
                    GridPT    = SetGrids_SDM(proinfo,prc_level,level,pyramid_step, final_level_iteration, &Size_Grid2D, &py_resolution, &grid_resolution, subBoundary);
                    
                    const long int Grid_length = (long int)Size_Grid2D.width*(long int)Size_Grid2D.height;
                    
                    levelinfo.Size_Grid2D = &Size_Grid2D;
                    levelinfo.Grid_length = &Grid_length;
                    levelinfo.GridPts = GridPT;
                    levelinfo.grid_resolution = &grid_resolution;
                    
                    if(!flag_start)
                    {
                        printf("GridPT3 start\t seed flag %d\t filename %s\timage_resolution %f \n",proinfo.pre_DEMtif,proinfo.priori_DEM_tif,proinfo.resolution);
                        GridPT3 = SetGrid3PT_SDM(Size_Grid2D, Th_roh);
                    }
                    
                    if(flag_start)
                    {
                        printf("start ResizeGridPT3 pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                        GridPT3 = ResizeGirdPT3_SDM(pre_Size_Grid2D, Size_Grid2D, subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution);
                    }
                    
                    printf("end start ResizeGridPT3\n");
                    
                    pre_Size_Grid2D.width = Size_Grid2D.width;
                    pre_Size_Grid2D.height = Size_Grid2D.height;
                    pre_grid_resolution = grid_resolution;
                    
                    fprintf(fid,"level = %d, Completion of Gridinfo setup\t%d\t%d!!\n",level,Size_Grid2D.width,Size_Grid2D.height);
                    
                    fprintf(fid_header, "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", 1, 1, level, subBoundary[0], subBoundary[1], grid_resolution, Size_Grid2D.width,Size_Grid2D.height);
                    
                    double left_mag_var, left_mag_avg, right_mag_var, right_mag_avg;
                    
                    printf("load subimages\n");
                    
                    total_matching_candidate_pts = (long)Size_Grid2D.width*(long)Size_Grid2D.height;
                    
                    printf("End load subimages blunder_selected_level\n");
                    
                    int iteration        = 1;
                    if(level == 0)
                        iteration = final_level_iteration;
                    
                    int pre_matched_pts=10;
                    double matching_change_rate = 100;
                    double rate_th = 0.00999999;
                    int max_iteration = 10 - (pyramid_step - level)*2;
                    if(max_iteration < 3)
                        max_iteration = 3;
                    
                    NCCresultSDM *nccresult = (NCCresultSDM*)calloc(sizeof(NCCresultSDM),(long)Size_Grid2D.width*(long)Size_Grid2D.height);
                    
                    if(level == 0 &&  iteration == 3)
                        matching_change_rate = 0.001;
                    
                    bool check_pre_GridPT3 = false;
                    
                    while((Th_roh >= Th_roh_min || (matching_change_rate > rate_th)))
                    {
                        printf("%f \t %f\n",Th_roh,Th_roh_min);
                        
                        double Th_roh_update = 0;
                        
                        fprintf(fid,"Starting computation of NCC\n iteration = %u\tTh_roh = %f\tTh_roh_start = %f\tGrid size %d %d\n",
                                iteration, Th_roh,Th_roh_start,Size_Grid2D.width,Size_Grid2D.height);
                        
                        printf("template size =%d\n",Template_size);
                        
                        bool check_smooth = false;
                        if(grid_resolution >= 200 && level == 0)
                            check_smooth = true;
                        else if(level < 2 && grid_resolution < 300)
                            check_smooth = true;
                        
                        if(check_smooth)
                        {
                            bool check_while = false;
                            bool check_last_iter = false;
                            int check_size = 3;
                            if(grid_resolution >= 300)
                                check_size = 3;
                            else if(grid_resolution >= 200)
                                check_size = 4;
                            else if(grid_resolution >= 100)
                                check_size = 5;
                            else
                                check_size = 7;
                                
                            if(check_size < 3)
                                check_size = 3;
                            int total_size = (2*check_size+1)*(2*check_size+1);
                            int sm_iter = 0;
                            
                            uint8* t_ncc_array = (uint8*)calloc(sizeof(uint8),*levelinfo.Grid_length);
                            float* temp_col_shift = (float*)calloc(sizeof(float),*levelinfo.Grid_length);
                            float* temp_row_shift = (float*)calloc(sizeof(float),*levelinfo.Grid_length);
                            
                            while(check_while == 0 && sm_iter < 20)
                            {
                                #pragma omp parallel for schedule(guided)
                                for (long index = 0; index < *levelinfo.Grid_length; index++)
                                {
                                    long row = (floor(index/Size_Grid2D.width));
                                    long col = index%Size_Grid2D.width;
                                    long search_index = (long)row*(long)Size_Grid2D.width + (long)col;
                                    
                                    if(sm_iter > 0)
                                    {
                                        GridPT3[search_index].row_shift = temp_row_shift[search_index];
                                        GridPT3[search_index].col_shift = temp_col_shift[search_index];
                                    }
                                    
                                    if(check_last_iter)
                                    {
                                        t_ncc_array[search_index] = 1;
                                        check_while = 1;
                                    }
                                    else
                                    {
                                        if(t_ncc_array[search_index] != 2)
                                            t_ncc_array[search_index] = 1;
                                    }
                                }
                                
                                int count_null_cell = 0;
                                sm_iter++;
                                
                                #pragma omp parallel for schedule(guided) reduction(+:count_null_cell)
                                for (long index = 0; index < *levelinfo.Grid_length; index++)
                                {
                                    long row = (floor(index/Size_Grid2D.width));
                                    long col = index%Size_Grid2D.width;
                                    long search_index = row*(long)Size_Grid2D.width + col;
                                    const double p = 1.5;
                                    if (t_ncc_array[search_index] == 1)
                                    {
                                        double sum_row_shift=0;
                                        double sum_col_shift=0;
                                        double sum_roh = 0;
                                        
                                        long null_count_cell = 0;
                                        long count_cell = 0;
                                        
                                        for (long t_i = -check_size; t_i <= check_size;t_i++ )
                                        {
                                            for (long t_j = -check_size; t_j <= check_size; t_j++)
                                            {
                                                long index_row = row + t_i;
                                                long index_col = col + t_j;
                                                long t_index = index_row*(long)Size_Grid2D.width + index_col;
                                                
                                                if(index_row >= 0 && index_row < Size_Grid2D.height && index_col >= 0 && index_col < Size_Grid2D.width && t_i != 0 && t_j != 0)
                                                {
                                                    if(GridPT3[t_index].ortho_ncc > 0.0)
                                                    {
                                                        double ncc = 1 + GridPT3[t_index].ortho_ncc;
                                                        double weith_n = pow(ncc,10);
                                                        sum_row_shift += (GridPT3[t_index].row_shift*weith_n);
                                                        sum_col_shift += (GridPT3[t_index].col_shift*weith_n);
                                                        sum_roh += weith_n;
                                                        count_cell++;
                                                    }
                                                    
                                                    if(GridPT3[t_index].ortho_ncc < 0.3)
                                                        null_count_cell ++;
                                                }
                                            }
                                        }
                                        
                                        if (count_cell >= total_size*0.3 )
                                        {
                                            temp_col_shift[search_index] = sum_col_shift/sum_roh;
                                            temp_row_shift[search_index] = sum_row_shift/sum_roh;
                                            
                                            t_ncc_array[search_index] = 2;
                                            count_null_cell ++;
                                        }
                                        
                                        if (null_count_cell >= total_size*0.4 && level == 0 && final_level_iteration == 3)
                                        {
                                            temp_col_shift[search_index] = 0;
                                            temp_row_shift[search_index] = 0;
                                        }
                                    }
                                    else
                                    {
                                        temp_row_shift[search_index] = GridPT3[search_index].row_shift;
                                        temp_col_shift[search_index] = GridPT3[search_index].col_shift;
                                    }
                                }
                                
                                if(count_null_cell == 0)
                                    check_last_iter = true;
                                printf("sm iteration %d\t%d\n",sm_iter,count_null_cell);
                                
                                Update_ortho_NCC(proinfo, levelinfo, GridPT3, gsd_image1, gsd_image2, Rimageparam);
                            }
                            free(t_ncc_array);
                            free(temp_col_shift);
                            free(temp_row_shift);
                        }
                        
                        if(level < pyramid_step)
                            proinfo.SDM_SS = 3;
                        printf("kernel size %d\n",proinfo.SDM_SS);
                        
                        VerticalLineLocus_SDM(proinfo, levelinfo, nccresult, GridPT3, gsd_image1, gsd_image2, Rimageparam);
                        
                        vector<D3DPOINT> Matched_pts_col, Matched_pts_row;
                        long count_MPs = SelectMPs_SDM(proinfo, levelinfo, nccresult, GridPT3, Matched_pts_col, Matched_pts_row);
                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\t%d\n",1,1,level,iteration,count_MPs, Matched_pts_col.size());
                        
                        if(count_MPs > 100)
                            lower_level_match = true;
                        
                        if(lower_level_match)
                        {
                            if(level == 0 && iteration == 3)
                            {
                                echoprint_Gridinfo_SDM(proinfo, levelinfo, 1,1,level,iteration, GridPT3);
                                echoprint_adjustXYZ(proinfo, levelinfo, 1,1,level,iteration,GridPT3,1);
                            }
                            else
                            {
                                if(pre_matched_pts == 0)
                                    matching_change_rate = 0;
                                else
                                    matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
                                
                                printf("matching change rate pre curr %f\t%d\t%d\n",matching_change_rate,count_MPs,pre_matched_pts);
                                pre_matched_pts = count_MPs;
                                
                                if(level >= 4)
                                {
                                    if(level == pyramid_step)
                                    {
                                        if(iteration > 10)
                                            matching_change_rate = 0.001;
                                    }
                                    else if(iteration > 10 )
                                        matching_change_rate = 0.001;
                                }
                                else if(level == 3)
                                {
                                    if(iteration > 9 )
                                        matching_change_rate = 0.001;
                                }
                                else if(level <= 2)
                                {
                                    if(iteration > 2)
                                        matching_change_rate = 0.001;
                                    if(level == 0)
                                        matching_change_rate = 0.001;
                                }
                                
                                if(proinfo.pre_DEMtif && iteration >= 4)
                                {
                                    Th_roh = Th_roh_min - 1.0;
                                    matching_change_rate = 0.001;
                                }
                                
                                if(Th_roh >= Th_roh_min)
                                {
                                    if(level == 0)
                                        Th_roh_update        = (double)(Th_roh - 0.50);
                                    else if(level == 1)
                                        Th_roh_update        = (double)(Th_roh - 0.10);
                                    else if(level == 2)
                                        Th_roh_update        = (double)(Th_roh - 0.10);
                                    else if(level == 3)
                                        Th_roh_update        = (double)(Th_roh - 0.10);
                                    else
                                    {
                                        Th_roh_update        = (double)(Th_roh - 0.06);
                                    }
                                }
                                
                                if(level == 0)
                                {
                                    Pre_GridPT3        = CopyGridPT3(GridPT3, levelinfo);
                                    
                                    printf("update Pre_GridPT3\n");
                                }
                                else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
                                {
                                    Pre_GridPT3        = CopyGridPT3(GridPT3, levelinfo);
                                    printf("update Pre_GridPT3\n");
                                    
                                    check_pre_GridPT3 = true;
                                }
                                else
                                {
                                    GridPT3        = CopyGridPT3(GridPT3, levelinfo);
                                    printf("update GridPT3\n");
                                }
                                
                                //displacement computation
                                if(level <= pyramid_step)
                                {
                                    printf("Size_Grid2D %d\t%d\t%d\n",Size_Grid2D.width,Size_Grid2D.height,count_MPs);
                                    
                                    printf("TIN interpolation for col row shift %d\n",count_MPs);
                                   
                                    if(count_MPs > 5)
                                    {
                                        int count_tri;
                                        D3DPOINT * ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                        for(long i = 0 ; i < Matched_pts_col.size() ; i++)
                                            ptslists[i] = Matched_pts_col[i];
                                        
                                        Matched_pts_col.clear();
                                        vector<D3DPOINT>().swap(Matched_pts_col);
                                        
                                        double min_max[4] = {subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]};
                                        vector<UI3DPOINT> t_trilists;
                                        
                                        FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                        delete origTri;
                                        
                                        count_tri = t_trilists.size();
                                        UI3DPOINT *trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        
                                        for(long i = 0 ; i < t_trilists.size() ; i++)
                                            trilists[i] = t_trilists[i];
                                        
                                        t_trilists.clear();
                                        
                                        printf("end tingeneration %d\n",count_tri);
                                        
                                        if(level == 0)
                                            SetShiftFromTIN(count_MPs, count_tri, Pre_GridPT3, levelinfo, ptslists, trilists,0);
                                        else if(!check_pre_GridPT3)
                                            SetShiftFromTIN(count_MPs, count_tri, GridPT3,     levelinfo, ptslists, trilists,0);
                                        else
                                            SetShiftFromTIN(count_MPs, count_tri, Pre_GridPT3, levelinfo, ptslists, trilists,0);
                                        
                                        printf("end col computation %d\n",check_pre_GridPT3);
                                        
                                        free(trilists);
                                        free(ptslists);
                                        
                                        
                                        //row
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                        for(long i = 0 ; i < Matched_pts_row.size() ; i++)
                                            ptslists[i] = Matched_pts_row[i];
                                        
                                        origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                        delete origTri;
                                        
                                        count_tri = t_trilists.size();
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        
                                        for(long i = 0 ; i < t_trilists.size() ; i++)
                                            trilists[i] = t_trilists[i];
                                        
                                        t_trilists.clear();
                                        vector<UI3DPOINT>().swap(t_trilists);
                                        
                                        printf("end tingeneration %d\t%d\n",count_tri,count_MPs);
                                        
                                        double* ortho_ncc = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
                                        double* INCC = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
                                        
                                        if(level == 0)
                                        {
                                            SetShiftFromTIN(count_MPs, count_tri, Pre_GridPT3, levelinfo, ptslists,trilists,1);
                                            Update_ortho_NCC(proinfo, levelinfo, Pre_GridPT3, gsd_image1, gsd_image2, Rimageparam);
                                            
                                            if(iteration < 2)
                                                shift_filtering(proinfo, Pre_GridPT3, levelinfo);
                                            
                                            check_pre_GridPT3 = false;
                                        }
                                        else if(!check_pre_GridPT3)
                                        {
                                            SetShiftFromTIN(count_MPs, count_tri, GridPT3, levelinfo, ptslists,trilists,1);
                                            Update_ortho_NCC(proinfo, levelinfo, GridPT3, gsd_image1, gsd_image2, Rimageparam);
                                            
                                            if(level >= 2 || (level == 1 && iteration < 3))
                                                shift_filtering(proinfo, GridPT3, levelinfo);
                                        }
                                        else
                                        {
                                            SetShiftFromTIN(count_MPs, count_tri, Pre_GridPT3, levelinfo, ptslists,trilists,1);
                                            Update_ortho_NCC(proinfo, levelinfo, Pre_GridPT3, gsd_image1, gsd_image2, Rimageparam);
                                    
                                            if(level >= 2 || (level == 1 && iteration < 3))
                                                shift_filtering(proinfo, Pre_GridPT3, levelinfo);
                                            
                                            check_pre_GridPT3 = false;
                                        }
                                        printf("free(ortho_ncc)\n");
                                        free(ortho_ncc);
                                        printf("free(ortho_ncc)\n");
                                        free(INCC);
                                        printf("free(INCC)\n");
                                        free(trilists);
                                        printf("free(trilists)\n");
                                        free(ptslists);
                                        printf("free(ptslists)\n");
                                        
                                        printf("end compute\n");
                                    }
                                }
                            }
                            
                            printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SetHeightRange\n",1,1,level,iteration);
                            
                            
                            fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\n",
                                    1,1,level,iteration,count_MPs);
                            printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\n",
                                   1,1,level,iteration,count_MPs);
                            
                            printf("Size_Grid2D %d\t%d\n",Size_Grid2D.width,Size_Grid2D.height);
                        }
                        
                        if (level == 0 && iteration == 3)
                            *matching_number = count_MPs;
                        
                        if(lower_level_match)
                        {
                            flag_start            = true;
                            iteration++;
                        }
                        
                        if(level == 0)
                            Th_roh            = (double)(Th_roh - 0.50);
                        else if(level == 1)
                            Th_roh            = (double)(Th_roh - 0.10);
                        else if(level == 2)
                            Th_roh            = (double)(Th_roh - 0.10);
                        else if(level == 3)
                            Th_roh            = (double)(Th_roh - 0.10);
                        else
                            Th_roh            = (double)(Th_roh - 0.06);
                        
                        if(lower_level_match)
                        {
                            printf("th %f\t%f\t%f\t%f\n",Th_roh,Th_roh_min,matching_change_rate,rate_th);
                            if(Th_roh < Th_roh_min && matching_change_rate > rate_th)
                            {
                                if(level == 0)
                                    Th_roh            = (double)(Th_roh + 0.10);
                                else if(level == 1)
                                    Th_roh            = (double)(Th_roh + 0.10);
                                else if(level == 2)
                                    Th_roh            = (double)(Th_roh + 0.10);
                                else if(level == 3)
                                    Th_roh            = (double)(Th_roh + 0.10);
                                else
                                    Th_roh            = (double)(Th_roh + 0.06);
                                printf("th %f\t%f\n",Th_roh,Th_roh_min);
                            }
                        }
                        
                        if (!lower_level_match && Th_roh < Th_roh_min)
                        {
                            iteration++;
                            matching_change_rate = 0.001;
                            Th_roh_min = 0.4;
                        }
                        
                        if(level == 0)
                            final_level_iteration = iteration;
                        
                    }
                    
                    printf("DEM generation(%%) = %4.2f%% !!\n",(double)(pyramid_step+1 - level)/(double)(pyramid_step+1)*100);
                    
                    if(level > 0)
                        level    = level - 1;
                    
                    if(level == 0 && final_level_iteration == 4)
                        level = -1;
                    
                    
                    if(!lower_level_match && level > 1)
                    {
                        lower_level_match    = true;
                        flag_start            = false;
                    }
                    free(GridPT);
                    
                    printf("release Grid_wgs, nccresult\n");
                    free(nccresult);
                }
                printf("relese data size\n");
                
                for(int iter_level = 0 ; iter_level < level + 1; iter_level++)
                {
                    for(int image_index = 0 ; image_index < proinfo.number_of_images ; image_index++)
                    {
                        free(SubImages[iter_level][image_index]);
                        free(SubOriImages[iter_level][image_index]);
                        free(SubMagImages[iter_level][image_index]);
                    }
                    free(SubImages[iter_level]);
                    free(SubOriImages[iter_level]);
                    free(SubMagImages[iter_level]);
                }
                free(SubImages);
                free(SubOriImages);
                free(SubMagImages);
                
                for(int image_index = 0 ; image_index < proinfo.number_of_images ; image_index++)
                    free(data_size_lr[image_index]);
            
                free(data_size_lr);
                free(GridPT3);
                
                printf("release GridTP3\n");
                printf("DEM generation finish\n");
                
            }
        }
        fclose(fid);
        fclose(fid_header);
  
 //       RemoveFiles_SDM(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,proinfo.pyramid_level,0);
 //       printf("done remove tempfiles\n");
    }
}

bool subsetImage_SDM(ProInfo proinfo, double *subBoundary, D2DPOINT *startpos, CSize* subsetsize, uint16 **Sourceimage)
{
    bool ret = false;
    
    CSize LImagesize, RImagesize;
    
    if(GetImageSize(proinfo.Imagefilename[0],&LImagesize) && GetImageSize(proinfo.Imagefilename[1],&RImagesize))
    {
        long Lcols[2], Lrows[2], Rcols[2], Rrows[2];
        if(GetsubareaImage_GeoTiff(proinfo,proinfo.Imagefilename[0],LImagesize,proinfo.LBoundary,Lcols,Lrows) &&
           GetsubareaImage_GeoTiff(proinfo,proinfo.Imagefilename[1],RImagesize,proinfo.RBoundary,Rcols,Rrows) )
        {
            printf("image bits %d\n",proinfo.image_bits);
            Sourceimage[0] = SubsetImageFrombitsToUint16(proinfo.image_bits, proinfo.Imagefilename[0], Lcols, Lrows, &subsetsize[0]);
            Sourceimage[1] = SubsetImageFrombitsToUint16(proinfo.image_bits, proinfo.Imagefilename[1], Rcols, Rrows, &subsetsize[1]);
            
            startpos[0].m_X    = (double)(Lcols[0]);
            startpos[0].m_Y    = (double)(Lrows[0]);
            startpos[1].m_X    = (double)(Rcols[0]);
            startpos[1].m_Y    = (double)(Rrows[0]);
            
            printf("write subimage right\n");
    
            ret        = true;
        }
    }
    
    return ret;
}

void SetThs_SDM(int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, uint8 pyramid_step)
{
    if(level == 0)
    {
        *Th_roh             = (double)(0.30 - 0.10*(final_level_iteration-1));
        *Th_roh_min         = (double)(0.09);
    }
    else if(level == pyramid_step)
    {
        *Th_roh             = (double)(0.80);
        *Th_roh_min         = (double)(0.09);
    }
    else if(level == pyramid_step - 1)
    {
        *Th_roh             = (double)(0.70);
        *Th_roh_min         = (double)(0.09);
    }
    else if(level == pyramid_step - 2)
    {
        *Th_roh             = (double)(0.50);
        *Th_roh_min         = (double)(0.09);
    }
    else
    {
        *Th_roh             = (double)(0.30);
        *Th_roh_min         = (double)(0.09);
    }
    
    if(level >= 5)
        *Th_roh_next        = (double)(0.50);
    else if(level == 4)
        *Th_roh_next        = (double)(0.40);
    else if(level == 3)
        *Th_roh_next        = (double)(0.30);
    else if(level == 2)
        *Th_roh_next        = (double)(0.30);
    else if(level == 1)
        *Th_roh_next        = (double)(0.30);
    else
        *Th_roh_next        = (double)(0.10);
    
    *Th_roh_start        = (double)(*Th_roh);
}

D2DPOINT *SetGrids_SDM(ProInfo proinfo, const int prc_level, const int level, const int start_py, const int final_level_iteration, CSize *Size_Grid2D, double *py_resolution, double *grid_resolution, const double *subBoundary)
{
    if(level >= proinfo.end_level)
        *py_resolution     = (double)(proinfo.resolution*pwrtwo(prc_level + 1));
    else
        *py_resolution     = (double)(proinfo.resolution*pwrtwo(prc_level));
    
    *grid_resolution = *py_resolution;
    
    printf("pre resolution %f\t level %d\t final_level_iteration %d\n",*py_resolution,level,final_level_iteration);
    
    if(level == start_py || level > 2)
    {
        *py_resolution = (int)(*py_resolution/2.0);
        if(*py_resolution < proinfo.DEM_resolution/2.0)
            *py_resolution = (floor)(proinfo.DEM_resolution/2.0);
    }
    else
    {
        if(level > 1)
        {
            if((*py_resolution) <= 2 && proinfo.DEM_resolution <= 2) //low-res DEM more than 8m
                *py_resolution = 2;
            
            if(*py_resolution < proinfo.DEM_resolution/2.0)
                *py_resolution = (floor)(proinfo.DEM_resolution/2.0);
        }
        else if(level == 1)
        {
            if(*py_resolution > 2)
            {
                if((*py_resolution) >= proinfo.DEM_resolution/2.0)
                    *py_resolution = (floor)((*py_resolution)/2.0);
            }
            
            if(*py_resolution < proinfo.DEM_resolution/2.0)
                *py_resolution = (floor)(proinfo.DEM_resolution/2.0);
        }
        else if(level == 0)
        {
            if(final_level_iteration == 1)
            {
                if(*py_resolution < 2)
                {
                    if((*py_resolution)*4 > proinfo.DEM_resolution)
                        *py_resolution   = (*py_resolution)*4;
                    else
                    {
                        if((*py_resolution) >= proinfo.DEM_resolution/2.0)
                            *py_resolution = (floor)(proinfo.DEM_resolution/2.0);
                        else
                            *py_resolution = proinfo.DEM_resolution;
                    }
                }
                else
                    *py_resolution = proinfo.DEM_resolution;
            }
            else if(final_level_iteration == 2)
            {
                if(*py_resolution < 2)
                {
                    if((*py_resolution)*2 > proinfo.DEM_resolution)
                        *py_resolution   = (*py_resolution)*2;
                    else
                        *py_resolution   = proinfo.DEM_resolution;
                }
                else
                    *py_resolution = proinfo.DEM_resolution;
            }
            else
                *py_resolution   = proinfo.DEM_resolution;
        }
    }
    
    if(*py_resolution > 100)
        *py_resolution = (int)((*py_resolution)/50.0)*50.0;
    
    *grid_resolution = *py_resolution;
    
    Size_Grid2D->width    = (int)( ceil( ((double)(subBoundary[2] - subBoundary[0])/(*grid_resolution))/2.0)*2 );
    Size_Grid2D->height    = (int)( ceil( ((double)(subBoundary[3] - subBoundary[1])/(*grid_resolution))/2.0)*2 );

    
    printf("DEM resolution %f\tresolution %f\t size %d\t%d\n",proinfo.DEM_resolution,*py_resolution,Size_Grid2D->width,Size_Grid2D->height);
    D2DPOINT *GridPT = SetDEMGrid(subBoundary, *grid_resolution, *grid_resolution, Size_Grid2D);
    
    return GridPT;
}

UGRIDSDM *SetGrid3PT_SDM(const CSize Size_Grid2D, const double Th_roh)
{
    UGRIDSDM *GridPT3 = NULL;
    long total_grid_counts = (long)Size_Grid2D.height*(long)Size_Grid2D.width;
    
    GridPT3                    = (UGRIDSDM*)calloc(sizeof(UGRIDSDM),total_grid_counts);
    
#pragma omp parallel for
    for(long i=0;i<total_grid_counts;i++)
    {
        GridPT3[i].roh              = Th_roh;
        GridPT3[i].ortho_ncc        = 0;
        GridPT3[i].col_shift        = 0.0;
        GridPT3[i].row_shift        = 0.0;
    }

    return GridPT3;
}

UGRIDSDM* ResizeGirdPT3_SDM(const CSize preSize, const CSize resize_Size, const double* Boundary, const D2DPOINT *resize_Grid, UGRIDSDM *preGridPT3, const double pre_gridsize)
{
    
    UGRIDSDM *resize_GridPT3 = (UGRIDSDM *)calloc(sizeof(UGRIDSDM),(long)resize_Size.height*(long)resize_Size.width);
    
    for(long row=0;row<resize_Size.height;row++)
    {
        for(long col=0;col<resize_Size.width;col++)
        {
            long index = row*(long)resize_Size.width + col;
            double X = resize_Grid[index].m_X;
            double Y = resize_Grid[index].m_Y;
            int pos_c = (int)((X - Boundary[0])/pre_gridsize);
            int pos_r = (int)((Y - Boundary[1])/pre_gridsize);
            long pre_index = pos_r*(long)preSize.width + pos_c;
            if(pos_c >= 0 && pos_c < preSize.width && pos_r >= 0 && pos_r < preSize.height)
            {
                resize_GridPT3[index].roh            = preGridPT3[pre_index].roh;
                resize_GridPT3[index].ortho_ncc        = preGridPT3[pre_index].ortho_ncc;
                resize_GridPT3[index].col_shift        = preGridPT3[pre_index].col_shift;
                resize_GridPT3[index].row_shift        = preGridPT3[pre_index].row_shift;
            }
            else
            {
                resize_GridPT3[index].roh            = 0.0;
                resize_GridPT3[index].ortho_ncc        = 0;
                resize_GridPT3[index].col_shift     = 0;
                resize_GridPT3[index].row_shift     = 0;
            }
        }
    }
    
    printf("before release preGirdPT3\n");
    
    free(preGridPT3);
    
    printf("after release preGirdPT3\n");
    
    return resize_GridPT3;
}

bool VerticalLineLocus_SDM(ProInfo proinfo, LevelInfo &plevelinfo, NCCresultSDM* nccresult, UGRIDSDM *GridPT3, const ImageGSD gsd_image1, const ImageGSD gsd_image2, double* Coreg_param)
{
    const int Pyramid_step = *(plevelinfo.Pyramid_step);
    int Template_size = *(plevelinfo.Template_size);
    
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(proinfo.resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
        
        printf("VerticalLineLocus : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }
    
    int Half_template_size = (int)(Template_size/2);
    const double subBoundary[4] = {plevelinfo.Boundary[0], plevelinfo.Boundary[1], plevelinfo.Boundary[2], plevelinfo.Boundary[3]};
    
    const long numofpts = *(plevelinfo.Grid_length);
    
    printf("numofpts %d\t%d\t%d\tcoreg %f\t%f\n",numofpts,plevelinfo.Size_Grid2D->height,plevelinfo.Size_Grid2D->width,Coreg_param[0],Coreg_param[1]);
    const int reference_id = 0;
    const int target_id = 1;
#pragma omp parallel
    {
        SetKernel rsetkernel(reference_id,target_id,Half_template_size);
#pragma omp for schedule(guided)
        for(long iter_count = 0 ; iter_count < numofpts ; iter_count++)
        {
            long pts_row = floor(iter_count/plevelinfo.Size_Grid2D->width);
            long pts_col = iter_count % plevelinfo.Size_Grid2D->width;
            long pt_index = pts_row*(long)plevelinfo.Size_Grid2D->width + pts_col;
            
            if(pts_row >= 0 && pts_row < plevelinfo.Size_Grid2D->height && pts_col >= 0 && pts_col < plevelinfo.Size_Grid2D->width && pt_index >= 0 && pt_index < numofpts)
            {
                double max_1stroh = -1.0;
                
                int kernel_size = proinfo.SDM_SS;
                
                if(GridPT3[pt_index].ortho_ncc > 0.6)
                    kernel_size = kernel_size - 1;
                else if(GridPT3[pt_index].ortho_ncc > 0.2)
                    kernel_size = proinfo.SDM_SS;
                else
                    kernel_size = kernel_size + 1;
                
                if(GridPT3[pt_index].ortho_ncc < 0.8)
                {
                    bool check_false_h = false;
                    for(int kernel_row = -kernel_size ; kernel_row <= kernel_size ; kernel_row++)
                    {
                        for(int kernel_col = -kernel_size ; kernel_col <= kernel_size ; kernel_col++)
                        {
                            if(!check_false_h)
                            {
                                nccresult[pt_index].result0 = (-1.0);
                                check_false_h = true;
                            }
                            const CSize LImagesize(plevelinfo.py_Sizes[reference_id][Pyramid_step]);
                            const CSize RImagesize(plevelinfo.py_Sizes[target_id][Pyramid_step]);
                            
                            D2DPOINT temp_GP(plevelinfo.GridPts[pt_index]),temp_GP_R(plevelinfo.GridPts[pt_index]);
                            
                            D2DPOINT Left_Imagecoord        = GetObjectToImage_single(1,temp_GP,proinfo.LBoundary,proinfo.resolution);
                            D2DPOINT Right_Imagecoord        = GetObjectToImage_single(1,temp_GP_R,proinfo.RBoundary,proinfo.resolution);
                            
                            D2DPOINT Left_Imagecoord_py     = OriginalToPyramid_single(Left_Imagecoord,plevelinfo.py_Startpos[reference_id],Pyramid_step);
                            D2DPOINT Right_Imagecoord_py    = OriginalToPyramid_single(Right_Imagecoord,plevelinfo.py_Startpos[target_id],Pyramid_step);
                            
                            Right_Imagecoord_py.m_Y += (kernel_row + (GridPT3[pt_index].row_shift + Coreg_param[0])/pwrtwo(Pyramid_step));
                            Right_Imagecoord_py.m_X += (kernel_col + (GridPT3[pt_index].col_shift + Coreg_param[1])/pwrtwo(Pyramid_step));
                            
                            if( Right_Imagecoord_py.m_Y >= 0 && Right_Imagecoord_py.m_Y < RImagesize.height && Right_Imagecoord_py.m_X >= 0 && Right_Imagecoord_py.m_X < RImagesize.width && Left_Imagecoord_py.m_Y >= 0 && Left_Imagecoord_py.m_Y < LImagesize.height && Left_Imagecoord_py.m_X >= 0 && Left_Imagecoord_py.m_X < LImagesize.width)
                            {
                                int Count_N[3] = {0};
                                double total_NCC = 0;
                                double temp_INCC_roh = 0;
                                
                                double im_resolution_mask = (gsd_image1.pro_GSD + gsd_image2.pro_GSD)/2.0;
                                
                                for(int row = -Half_template_size; row <= Half_template_size ; row++)
                                {
                                    for(int col = -Half_template_size; col <= Half_template_size ; col++)
                                    {
                                        double row_distance = row*im_resolution_mask*pwrtwo(Pyramid_step);
                                        double col_distance = col*im_resolution_mask*pwrtwo(Pyramid_step);
                                        
                                        double row_pixel_left = row_distance/(gsd_image1.row_GSD*pwrtwo(Pyramid_step));
                                        double col_pixel_left = col_distance/(gsd_image1.col_GSD*pwrtwo(Pyramid_step));
                                        
                                        double row_pixel_right = row_distance/(gsd_image2.row_GSD*pwrtwo(Pyramid_step));
                                        double col_pixel_right = col_distance/(gsd_image2.col_GSD*pwrtwo(Pyramid_step));
                                        
                                        const int radius2  = (row*row + col*col);
                                        if(radius2 <= (Half_template_size+1)*(Half_template_size+1))
                                        {
                                            D2DPOINT pos_left(Left_Imagecoord_py.m_X + col_pixel_left, Left_Imagecoord_py.m_Y + row_pixel_left);
                                            D2DPOINT pos_right(Right_Imagecoord_py.m_X + col_pixel_right, Right_Imagecoord_py.m_Y + row_pixel_right);
                                            
                                            SetVecKernelValue(rsetkernel,plevelinfo.py_Sizes[rsetkernel.reference_id][*plevelinfo.Pyramid_step], plevelinfo.py_Sizes[rsetkernel.ti][*plevelinfo.Pyramid_step], plevelinfo.py_Images, plevelinfo.py_MagImages, row, col, pos_left,pos_right, radius2, Count_N);
                                        }
                                    }
                                }
                                
                                ComputeMultiNCC(rsetkernel, 0, Count_N, total_NCC, temp_INCC_roh);
                                
                                if(max_1stroh < temp_INCC_roh)
                                {
                                    max_1stroh = temp_INCC_roh;
                                    nccresult[pt_index].result0 = max_1stroh;
                                    
                                    nccresult[pt_index].result2.m_X = Left_Imagecoord_py.m_X;
                                    nccresult[pt_index].result2.m_Y = Left_Imagecoord_py.m_Y;
                                    nccresult[pt_index].result3.m_X = kernel_col + GridPT3[pt_index].col_shift/pwrtwo(Pyramid_step);
                                    nccresult[pt_index].result3.m_Y = kernel_row + GridPT3[pt_index].row_shift/pwrtwo(Pyramid_step);
                                }
                            }
                        }
                    }
                }
                else
                {
                    nccresult[pt_index].result0 = GridPT3[pt_index].ortho_ncc;
                    nccresult[pt_index].result3.m_X = GridPT3[pt_index].col_shift/pwrtwo(Pyramid_step);
                    nccresult[pt_index].result3.m_Y = GridPT3[pt_index].row_shift/pwrtwo(Pyramid_step);
                }
            }
        }
    }
    
    return true;
}

long SelectMPs_SDM(ProInfo proinfo, LevelInfo &rlevelinfo, NCCresultSDM* roh_height, UGRIDSDM *GridPT3, vector<D3DPOINT> &Matched_pts_col, vector<D3DPOINT> &Matched_pts_row)
{
    uint8 prc_level = *rlevelinfo.Pyramid_step;
    long count_MPs = 0;
    
    double minimum_Th = 0.3;
    printf("minimum TH %f\n",minimum_Th);
    double roh_th    = 0.10;
    
    for(long iter_index = 0 ; iter_index < *rlevelinfo.Grid_length ; iter_index++)
    {
        long row        = (floor(iter_index/rlevelinfo.Size_Grid2D->width));
        long col        = iter_index % rlevelinfo.Size_Grid2D->width;
        
        if(row >= 0 && row < rlevelinfo.Size_Grid2D->height && col >= 0 && col < rlevelinfo.Size_Grid2D->width)
        {
            long grid_index = row*(long)rlevelinfo.Size_Grid2D->width + col;
            GridPT3[grid_index].ortho_ncc = roh_height[grid_index].result0;
            
            bool roh_index        = false;
            
            if(roh_height[grid_index].result0 > GridPT3[grid_index].roh - roh_th && roh_height[grid_index].result0 > minimum_Th)
                roh_index = true;
            
            //Set the matched pts and information
            if(roh_index && rlevelinfo.GridPts[grid_index].m_X > rlevelinfo.Boundary[0] && rlevelinfo.GridPts[grid_index].m_X < rlevelinfo.Boundary[2] && rlevelinfo.GridPts[grid_index].m_Y > rlevelinfo.Boundary[1] && rlevelinfo.GridPts[grid_index].m_Y < rlevelinfo.Boundary[3])
            {
                count_MPs++;
                D3DPOINT temp_col(rlevelinfo.GridPts[grid_index].m_X, rlevelinfo.GridPts[grid_index].m_Y, roh_height[grid_index].result3.m_X*pwrtwo(prc_level), 1);
                D3DPOINT temp_row(rlevelinfo.GridPts[grid_index].m_X, rlevelinfo.GridPts[grid_index].m_Y, roh_height[grid_index].result3.m_Y*pwrtwo(prc_level), 1);
                Matched_pts_col.push_back(temp_col);
                Matched_pts_row.push_back(temp_row);
         
                // update max_roh value
                GridPT3[grid_index].roh        = roh_height[grid_index].result0;
                if(GridPT3[grid_index].roh < minimum_Th)
                    GridPT3[grid_index].roh = minimum_Th;
            }
        }
    }
    
    return count_MPs;
}

void echoprint_Gridinfo_SDM(ProInfo proinfo, LevelInfo &rlevelinfo, int row, int col, int level, int iteration, UGRIDSDM *GridPT3)
{
    uint8 prc_level = *rlevelinfo.Pyramid_step;
    CSize LImagesize = rlevelinfo.py_Sizes[0][prc_level];
    CSize RImagesize = rlevelinfo.py_Sizes[1][prc_level];
 
    printf("Size_Grid2D %d\t%d\n",rlevelinfo.Size_Grid2D->width,rlevelinfo.Size_Grid2D->height);
    
    float* Roh = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
    for(long k=0;k<rlevelinfo.Size_Grid2D->height;k++)
    {
        for(long j=0;j<rlevelinfo.Size_Grid2D->width;j++)
        {
            long grid_index    = (rlevelinfo.Size_Grid2D->height-1 - k)*(long)rlevelinfo.Size_Grid2D->width + j;
            long matlab_index    = k*(long)rlevelinfo.Size_Grid2D->width + j;
            double coord_x = rlevelinfo.Boundary[0] + j*proinfo.DEM_resolution;
            double coord_y = rlevelinfo.Boundary[1] + k*proinfo.DEM_resolution;
            
            int pos_lc = (int)((coord_x - proinfo.LBoundary[0])/(proinfo.resolution*pwrtwo(prc_level)));
            int pos_lr = (int)((proinfo.LBoundary[3] - coord_y)/(proinfo.resolution*pwrtwo(prc_level)));
            
            int pos_rc = (int)((coord_x - proinfo.RBoundary[0])/(proinfo.resolution*pwrtwo(prc_level)));
            int pos_rr = (int)((proinfo.RBoundary[3] - coord_y)/(proinfo.resolution*pwrtwo(prc_level)));
            
            long int l_index = pos_lr*LImagesize.width + pos_lc;
            long int r_index = pos_rr*RImagesize.width + pos_rc;
            
            if(pos_lc >= 0 && pos_lc < LImagesize.width && pos_lr >= 0 && pos_lr < LImagesize.height &&
               pos_rc >= 0 && pos_rc < RImagesize.width && pos_rr >= 0 && pos_rr < RImagesize.height)
            {
                int l_value = rlevelinfo.py_Images[0][l_index];
                int r_value = rlevelinfo.py_Images[1][r_index];
                
                if(l_value > 0 && r_value > 0)
                    Roh[grid_index] = GridPT3[matlab_index].ortho_ncc;
                else
                    Roh[grid_index] = -1.0;
            }
            else
                Roh[grid_index] = -1.0;
        }
    }
  
    char DEM_str[500];
    sprintf(DEM_str, "%s/%s_roh.tif", proinfo.save_filepath, proinfo.Outputpath_name);
    WriteGeotiff(DEM_str, Roh, rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height, proinfo.DEM_resolution, rlevelinfo.Boundary[0], rlevelinfo.Boundary[3], rlevelinfo.param->projection, rlevelinfo.param->utm_zone, rlevelinfo.param->bHemisphere, 4);
    
    free(Roh);
}

void echoprint_adjustXYZ(ProInfo proinfo, LevelInfo &rlevelinfo, int row, int col,int level, int iteration, UGRIDSDM *GridPT3, int d_date)
{
    uint8 prc_level = *rlevelinfo.Pyramid_step;
    CSize LImagesize = rlevelinfo.py_Sizes[0][prc_level];
    CSize RImagesize = rlevelinfo.py_Sizes[1][prc_level];
    
    FILE *outfile_Xshift, *outfile_Yshift;
    char t_str[500];
    sprintf(t_str,"%s/txt/tin_xshift_level_%d_%d_%d_iter_%d.txt",proinfo.save_filepath,row,col,level,iteration);
    outfile_Xshift    = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/tin_yshift_level_%d_%d_%d_iter_%d.txt",proinfo.save_filepath,row,col,level,iteration);
    outfile_Yshift    = fopen(t_str,"w");

    printf("Size_Grid2D %d\t%d\tres = %f\n",rlevelinfo.Size_Grid2D->width,rlevelinfo.Size_Grid2D->height,proinfo.resolution);
    
    float* Mag = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
    float* VxShift = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
    float* VyShift = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
    
    for(long k=0;k<rlevelinfo.Size_Grid2D->height;k++)
    {
        for(long j=0;j<rlevelinfo.Size_Grid2D->width;j++)
        {
            long grid_index    = (rlevelinfo.Size_Grid2D->height-1 - k)*(long)rlevelinfo.Size_Grid2D->width + j;
            long matlab_index    = k*(long)rlevelinfo.Size_Grid2D->width + j;
            
            double coord_x = rlevelinfo.Boundary[0] + j*proinfo.DEM_resolution;
            double coord_y = rlevelinfo.Boundary[1] + k*proinfo.DEM_resolution;
            
            int pos_lc = (int)((coord_x - proinfo.LBoundary[0])/(proinfo.resolution*pwrtwo(prc_level)));
            int pos_lr = (int)((proinfo.LBoundary[3] - coord_y)/(proinfo.resolution*pwrtwo(prc_level)));
            
            int pos_rc = (int)((coord_x - proinfo.RBoundary[0])/(proinfo.resolution*pwrtwo(prc_level)));
            int pos_rr = (int)((proinfo.RBoundary[3] - coord_y)/(proinfo.resolution*pwrtwo(prc_level)));
            
            long int l_index = pos_lr*LImagesize.width + pos_lc;
            long int r_index = pos_rr*RImagesize.width + pos_rc;
            
            double DX = GridPT3[matlab_index].col_shift*proinfo.resolution/proinfo.SDM_days;
            double DY = -GridPT3[matlab_index].row_shift*proinfo.resolution/proinfo.SDM_days;
      
            if(d_date == 2)
            {
                DX = -DX;
                DY = -DY;
            }
            
            if(pos_lc >= 0 && pos_lc < LImagesize.width && pos_lr >= 0 && pos_lr < LImagesize.height &&
               pos_rc >= 0 && pos_rc < RImagesize.width && pos_rr >= 0 && pos_rr < RImagesize.height)
            {
                int l_value = rlevelinfo.py_Images[0][l_index];
                int r_value = rlevelinfo.py_Images[1][r_index];
                
                if(l_value > 0 && r_value > 0)
                {
                    fprintf(outfile_Xshift,"%f\t", DX);
                    fprintf(outfile_Yshift,"%f\t", DY);
                    VxShift[grid_index] = DX;
                    VyShift[grid_index] = DY;
                    Mag[grid_index] = sqrt(DX*DX + DY*DY);
                }
                else
                {
                    fprintf(outfile_Xshift,"0\t");
                    fprintf(outfile_Yshift,"0\t");
                    VxShift[grid_index] = 0;
                    VyShift[grid_index] = 0;
                }
            }
            else
            {
                fprintf(outfile_Xshift,"0\t");
                fprintf(outfile_Yshift,"0\t");
                VxShift[grid_index] = 0;
                VyShift[grid_index] = 0;
            }
            fprintf(outfile_Xshift,"\n");
            fprintf(outfile_Yshift,"\n");
        }
    }
    fclose(outfile_Xshift);
    fclose(outfile_Yshift);
    
    char DEM_str[500];
    sprintf(DEM_str, "%s/%s_dx.tif", proinfo.save_filepath, proinfo.Outputpath_name);
    WriteGeotiff(DEM_str, VxShift, rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height, proinfo.DEM_resolution, rlevelinfo.Boundary[0], rlevelinfo.Boundary[3], rlevelinfo.param->projection, rlevelinfo.param->utm_zone, rlevelinfo.param->bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_dy.tif", proinfo.save_filepath, proinfo.Outputpath_name);
    WriteGeotiff(DEM_str, VyShift, rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height, proinfo.DEM_resolution, rlevelinfo.Boundary[0], rlevelinfo.Boundary[3], rlevelinfo.param->projection, rlevelinfo.param->utm_zone, rlevelinfo.param->bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_dmag.tif", proinfo.save_filepath, proinfo.Outputpath_name);
    WriteGeotiff(DEM_str, Mag, rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height, proinfo.DEM_resolution, rlevelinfo.Boundary[0], rlevelinfo.Boundary[3], rlevelinfo.param->projection, rlevelinfo.param->utm_zone, rlevelinfo.param->bHemisphere, 4);
    
    free(VxShift);
    free(VyShift);
    free(Mag);
}

bool Update_ortho_NCC(ProInfo proinfo, LevelInfo &rlevelinfo, UGRIDSDM *GridPT3, const ImageGSD gsd_image1, const ImageGSD gsd_image2, double* Coreg_param)
{
    int Pyramid_step = *rlevelinfo.Pyramid_step;
    int Template_size = *rlevelinfo.Template_size;
    
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(proinfo.resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
    }
    
    const int Half_template_size = (int)(Template_size/2);
    double subBoundary[4] = {rlevelinfo.Boundary[0], rlevelinfo.Boundary[1], rlevelinfo.Boundary[2], rlevelinfo.Boundary[3]};
    
    const long numofpts = *rlevelinfo.Grid_length;
    
    const int reference_id = 0;
    const int target_id = 1;
#pragma omp parallel
    {
        
        SetKernel rsetkernel(reference_id,target_id,Half_template_size);

#pragma omp for schedule(guided)
        for(long iter_count = 0 ; iter_count < numofpts ; iter_count++)
        {
            long pts_row = floor(iter_count/rlevelinfo.Size_Grid2D->width);
            long pts_col = iter_count % rlevelinfo.Size_Grid2D->width;
            long pt_index = pts_row*(long)rlevelinfo.Size_Grid2D->width + pts_col;
            
            if(pt_index < *rlevelinfo.Grid_length && pts_row < rlevelinfo.Size_Grid2D->height && pts_col < rlevelinfo.Size_Grid2D->width && pts_row >= 0 && pts_col >= 0)
            {
                const CSize LImagesize(rlevelinfo.py_Sizes[reference_id][Pyramid_step]);
                const CSize RImagesize(rlevelinfo.py_Sizes[target_id][Pyramid_step]);
                
                D2DPOINT temp_GP(rlevelinfo.GridPts[pt_index]),temp_GP_R(rlevelinfo.GridPts[pt_index]);
                
                D2DPOINT Left_Imagecoord        = GetObjectToImage_single(1,temp_GP,proinfo.LBoundary,proinfo.resolution);
                D2DPOINT Right_Imagecoord       = GetObjectToImage_single(1,temp_GP_R,proinfo.RBoundary,proinfo.resolution);
                
                D2DPOINT Left_Imagecoord_py     = OriginalToPyramid_single(Left_Imagecoord,rlevelinfo.py_Startpos[reference_id],Pyramid_step);
                D2DPOINT Right_Imagecoord_py    = OriginalToPyramid_single(Right_Imagecoord,rlevelinfo.py_Startpos[target_id],Pyramid_step);
                
                Right_Imagecoord_py.m_Y += (GridPT3[pt_index].row_shift + Coreg_param[0])/pwrtwo(Pyramid_step);
                Right_Imagecoord_py.m_X += (GridPT3[pt_index].col_shift + Coreg_param[1])/pwrtwo(Pyramid_step);
                
                if( Right_Imagecoord_py.m_Y >= 0 && Right_Imagecoord_py.m_Y < RImagesize.height && Right_Imagecoord_py.m_X >= 0 && Right_Imagecoord_py.m_X < RImagesize.width && Left_Imagecoord_py.m_Y >= 0 && Left_Imagecoord_py.m_Y < LImagesize.height && Left_Imagecoord_py.m_X >= 0 && Left_Imagecoord_py.m_X < LImagesize.width)
                {
                    int Count_N[3] = {0};
                    double total_NCC = 0;
                    double temp_INCC_roh = 0;
                    
                    double im_resolution_mask = (gsd_image1.pro_GSD + gsd_image2.pro_GSD)/2.0;
                    
                    for(int row = -Half_template_size; row <= Half_template_size ; row++)
                    {
                        for(int col = -Half_template_size; col <= Half_template_size ; col++)
                        {
                            double row_distance = row*im_resolution_mask*pwrtwo(Pyramid_step);
                            double col_distance = col*im_resolution_mask*pwrtwo(Pyramid_step);
                            
                            double row_pixel_left = row_distance/(gsd_image1.row_GSD*pwrtwo(Pyramid_step));
                            double col_pixel_left = col_distance/(gsd_image1.col_GSD*pwrtwo(Pyramid_step));
                            
                            double row_pixel_right = row_distance/(gsd_image2.row_GSD*pwrtwo(Pyramid_step));
                            double col_pixel_right = col_distance/(gsd_image2.col_GSD*pwrtwo(Pyramid_step));
                            
                            const int radius2  = (row*row + col*col);
                            if(radius2 <= (Half_template_size+1)*(Half_template_size+1))
                            {
                                D2DPOINT pos_left(Left_Imagecoord_py.m_X + col_pixel_left, Left_Imagecoord_py.m_Y + row_pixel_left);
                                D2DPOINT pos_right(Right_Imagecoord_py.m_X + col_pixel_right, Right_Imagecoord_py.m_Y + row_pixel_right);
                                
                                SetVecKernelValue(rsetkernel,rlevelinfo.py_Sizes[rsetkernel.reference_id][*rlevelinfo.Pyramid_step], rlevelinfo.py_Sizes[rsetkernel.ti][*rlevelinfo.Pyramid_step], rlevelinfo.py_Images, rlevelinfo.py_MagImages, row, col, pos_left,pos_right, radius2, Count_N);
                            }
                        }
                    }
                    ComputeMultiNCC(rsetkernel, 0, Count_N, total_NCC, temp_INCC_roh);
                    
                    GridPT3[pt_index].ortho_ncc = temp_INCC_roh;
                }
            }
        }
        
    }
    
    return true;
}

UGRIDSDM* CopyGridPT3(UGRIDSDM *GridPT3, LevelInfo &rlevelinfo)
{
    UGRIDSDM * result                    = (UGRIDSDM*)calloc(*rlevelinfo.Grid_length,sizeof(UGRIDSDM));
    
    for(long k=0;k<rlevelinfo.Size_Grid2D->height;k++)
    {
        for(long j=0;j<rlevelinfo.Size_Grid2D->width;j++)
        {
            long matlab_index    = k*(long)rlevelinfo.Size_Grid2D->width + j;
            
            result[matlab_index].col_shift                  = GridPT3[matlab_index].col_shift;
            result[matlab_index].row_shift                  = GridPT3[matlab_index].row_shift;
            
            result[matlab_index].ortho_ncc                  = GridPT3[matlab_index].ortho_ncc;
            result[matlab_index].roh                        = GridPT3[matlab_index].roh;
        }
    }
    free(GridPT3);

    return result;
}

void SetShiftFromTIN(const long numOfPts, const long numOfTri, UGRIDSDM *GridPT3, LevelInfo &rlevelinfo, D3DPOINT *pts, UI3DPOINT *tris, const int b_dir)
{
    double Total_Min_Z        =  100000;
    double Total_Max_Z        = -100000;
    
    uint8 *m_bHeight        = (uint8*)calloc(*rlevelinfo.Grid_length,sizeof(uint8));
    
    for(long tcnt=0;tcnt<numOfTri;tcnt++)
    {
        const UI3DPOINT t_tri(tris[tcnt]);
        const int pdex0 = t_tri.m_X;
        const int pdex1 = t_tri.m_Y;
        const int pdex2 = t_tri.m_Z;
        
        if(pdex0 < numOfPts && pdex1 < numOfPts && pdex2 < numOfPts)
        {
            D3DPOINT TriP1(pts[pdex0]);
            D3DPOINT TriP2(pts[pdex1]);
            D3DPOINT TriP3(pts[pdex2]);
            int PixelMinXY[2], PixelMaxXY[2];
            double temp_MinZ, temp_MaxZ;
            SetTinBoundary(rlevelinfo, TriP1, TriP2, TriP3, PixelMinXY, PixelMaxXY, Total_Min_Z, Total_Max_Z, temp_MinZ, temp_MaxZ);
            
            for (long Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
            {
                for (long Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
                {
                    long Index= (long)rlevelinfo.Size_Grid2D->width*Row + Col;
                    
                    if(!m_bHeight[Index])
                    {
                        float Z = -1000.0;
                        bool rtn = false;
                        
                        D3DPOINT CurGPXY((Col)*(*rlevelinfo.grid_resolution) + rlevelinfo.Boundary[0],(Row)*(*rlevelinfo.grid_resolution) + rlevelinfo.Boundary[1],0,0);
                        rtn = IsTinInside(CurGPXY, TriP1, TriP2, TriP3, Z);
                        
                        if (rtn)
                        {
                            const double diff1 = SQRT(CurGPXY, TriP1, 2);
                            const double diff2 = SQRT(CurGPXY, TriP2, 2);
                            const double diff3 = SQRT(CurGPXY, TriP3, 2);
                            
                            if(diff1 == 0)
                                Z    = TriP1.m_Z;
                            else if(diff2 == 0)
                                Z    = TriP2.m_Z;
                            else if(diff3 == 0)
                                Z    = TriP3.m_Z;
                     
                            m_bHeight[Index] = 1;
                            if(!b_dir)
                                GridPT3[Index].col_shift = Z;
                            else
                                GridPT3[Index].row_shift = Z;
                        }
                    }
                }
            }
        }
    }
    
    if(m_bHeight)
        free(m_bHeight);
}

void shift_filtering(ProInfo proinfo, UGRIDSDM *GridPT3, LevelInfo &rlevelinfo)
{
    const CSize gridsize(rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height);
    const int pyramid_step(*rlevelinfo.Pyramid_step);
    const double DEM_resolution(proinfo.DEM_resolution);
    
    const long data_length = *rlevelinfo.Grid_length;
    float *temp_col_shift = (float*)malloc(sizeof(float)*data_length);
    float *temp_row_shift = (float*)malloc(sizeof(float)*data_length);
    
    for(long r = 0 ; r < gridsize.height ; r ++)
    {
        for(long c = 0 ; c < gridsize.width; c++)
        {
            long ori_index = r*(long)gridsize.width + c;
            
            temp_col_shift[ori_index] = GridPT3[ori_index].col_shift;
            temp_row_shift[ori_index] = GridPT3[ori_index].row_shift;
        }
    }
    
    int da = 45;
    if(pyramid_step == 3)
        da = 30;
    else if(pyramid_step == 2)
        da = 20;
    
    const int slope_step = 360/da;
    
    printf("angle step %d\t%d\n",slope_step,data_length);
    
    const int shift_max_pixel = (int)(((double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution)) );
    
    int kernal_size = 1;
    if(DEM_resolution == 1)
        kernal_size = 3;
    if(pyramid_step >= 5)
        kernal_size = 9;
    else if(pyramid_step == 4)
        kernal_size = 9;
    else if(pyramid_step == 3)
        kernal_size = 7;
    else if(pyramid_step == 2)
        kernal_size = 5;
    
    printf("level %d\tshift_max_pixel %d\n",proinfo.pyramid_level,shift_max_pixel);
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < data_length ; iter_count ++)
    {
        long r = (floor(iter_count/gridsize.width));
        long c = iter_count % gridsize.width;
        long ori_index = r*(long)gridsize.width + c;
        if(ori_index >= 0 && r >= 0 && r < gridsize.height && c >= 0 && c < gridsize.width && ori_index < data_length)
        {
            vector<double> save_col, save_row, save_diff, save_slope, save_roh;
            double save_col_min = 10000.;
            double save_col_max = -10000.;
            double save_row_min = 10000.;
            double save_row_max = -10000.;
            
            int* hist_slope = (int*)calloc(sizeof(int),slope_step);
            
            int hist_col[6000] = {0};
            int hist_col_id[6000];
            int hist_row[6000] = {0};
            int hist_row_id[6000];
            
            for(int k=0;k<3000;k++)
            {
                hist_col_id[k] = k;
                hist_row_id[k] = k;
            }
            
            long save_count = 0;
            for(long k = -kernal_size ; k <= kernal_size ; k++)
            {
                for(long j=-kernal_size ; j <= kernal_size ; j++)
                {
                    long index = (r+k)*(long)gridsize.width + (c+j);
                    if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width && index >= 0 && index < data_length)
                    {
                        save_col.push_back(GridPT3[index].col_shift);
                        save_row.push_back(GridPT3[index].row_shift);
                        save_diff.push_back(sqrt(k*k + j*j));
                        save_slope.push_back(atan2(save_row[save_count],save_col[save_count])*RadToDeg);
                        save_roh.push_back(pow(3.0,GridPT3[index].ortho_ncc*10));
                        
                        if(save_col[save_count] < save_col_min)
                            save_col_min = save_col[save_count];
                        
                        if(save_col[save_count] > save_col_max)
                            save_col_max = save_col[save_count];
                        
                        if(save_col[save_count] < save_row_min)
                            save_row_min = save_col[save_count];
                        
                        if(save_col[save_count] > save_row_max)
                            save_row_max = save_col[save_count];
                        
                        if(save_slope[save_count] < 0)
                            save_slope[save_count] += 360;
                        
                        int t_hist = (int)(save_slope[save_count]/da);
                        if(t_hist < slope_step)
                            hist_slope[t_hist]++;
                        else
                            hist_slope[t_hist-1]++;
                        
                        t_hist = (int)(save_col[save_count]);
                        if(abs(t_hist) < shift_max_pixel)
                            hist_col[t_hist+shift_max_pixel]++;
                        
                        t_hist = (int)(save_row[save_count]);
                        if(abs(t_hist) < shift_max_pixel)
                            hist_row[t_hist+shift_max_pixel]++;
                        
                        save_count++;
                    }
                }
            }
   
            save_count = save_col.size();
            int th_count = (int)(((kernal_size*2+1)*(kernal_size*2+1))/10.0);
            if(th_count < 3 )
                th_count = 3;
            if(save_count > th_count)
            {
                int max_hist = -1;
                int max_hist_col = -1;
                int max_hist_row = -1;
                int min_hist_col = 10000;
                int min_hist_row = 10000;
                int max_hist_pos;
                int max_hist_col_pos, max_hist_row_pos;
                int min_hist_col_pos, min_hist_row_pos;
                
                for(int k=0;k<6000;k++)
                {
                    if(max_hist_col < hist_col[k])
                    {
                        max_hist_col = hist_col[k];
                        max_hist_col_pos = k;
                    }
                    
                    if(max_hist_row < hist_row[k])
                    {
                        max_hist_row = hist_row[k];
                        max_hist_row_pos = k;
                    }
                    
                    if(min_hist_col > hist_col[k])
                    {
                        min_hist_col = hist_col[k];
                        min_hist_col_pos = k;
                    }
                    
                    if(min_hist_row > hist_row[k])
                    {
                        min_hist_row = hist_row[k];
                        min_hist_row_pos = k;
                    }
                    
                    if(k<slope_step)
                    {
                        if(max_hist < hist_slope[k])
                        {
                            max_hist = hist_slope[k];
                            max_hist_pos = k;
                        }
                    }
                }
          
                vector<double> selected_col, selected_row, selected_diff_col, selected_diff_row, selected_roh;
                
                int total_seleted_count_col = 0;
                int total_seleted_count_row = 0;
                const double p = 1.5;
                
                for(long k= 0; k<save_count ; k++)
                {
                    int t_hist_col = (int)(save_col[k]);
                    int t_hist_row = (int)(save_row[k]);
                    int t_hist_slope = (int)(save_slope[k]/da);
                    
                    if(pyramid_step >= 0 )
                    {
                        //column
                        selected_col.push_back(save_col[k]);
                        selected_roh.push_back(save_roh[k]);
                        if( save_diff[k] == 0)
                            selected_diff_col.push_back(0.1);
                        else
                        {
                            if( t_hist_col+3000 == max_hist_col_pos)
                                selected_diff_col.push_back(0.1);
                            else
                            {
                                if(t_hist_slope == max_hist_pos)
                                    selected_diff_col.push_back(0.1);
                                else
                                    selected_diff_col.push_back(save_diff[k]*2.0);
                            }
                        }
                        total_seleted_count_col++;
                        
                        //row
                        selected_row.push_back(save_row[k]);
                        if( save_diff[k] == 0)
                            selected_diff_row.push_back(0.1);
                        else
                        {
                            if( t_hist_row+3000 == max_hist_row_pos)
                                selected_diff_row.push_back(0.1);
                            else
                            {
                                if(t_hist_slope == max_hist_pos)
                                    selected_diff_row.push_back(0.1);
                                else
                                    selected_diff_row.push_back(save_diff[k]*2.0);
                            }
                        }
                        total_seleted_count_row++;
                    }
                    else
                    {
                        if(min_hist_col_pos == max_hist_col_pos || min_hist_row_pos == max_hist_row_pos)
                        {
                            selected_col.push_back(save_col[k]);
                            selected_row.push_back(save_row[k]);
                            selected_roh.push_back(save_roh[k]);
                            
                            if( save_diff[k] == 0)
                            {
                                selected_diff_col.push_back(0.1);
                                selected_diff_row.push_back(0.1);
                            }
                            else
                            {
                                selected_diff_col.push_back(save_diff[k]);
                                selected_diff_row.push_back(save_diff[k]);
                            }
                            total_seleted_count_col ++;
                            total_seleted_count_row ++;
                        }
                        else if((t_hist_col+3000 != min_hist_col_pos || t_hist_row+3000 != min_hist_row_pos) )
                        {
                            selected_col.push_back(save_col[k]);
                            selected_row.push_back(save_row[k]);
                            selected_roh.push_back(save_roh[k]);
                            
                            if( save_diff[k] == 0)
                            {
                                selected_diff_col.push_back(0.1);
                                selected_diff_row.push_back(0.1);
                            }
                            else
                            {
                                selected_diff_col.push_back(save_diff[k]);
                                selected_diff_row.push_back(save_diff[k]);
                            }
                            total_seleted_count_col ++;
                            total_seleted_count_row ++;
                        }
                    }
                }
                
                //col
                if(total_seleted_count_col > 0)
                {
                    double sum1 = 0;
                    double sum2 = 0;
                    double sum1_r = 0;
                    double sum2_r = 0;
                    
                    if(pyramid_step >= 0)
                    {
                        if((max_hist_col > save_count*0.4))
                        {
                            //IDW
                            for(int k=0;k<total_seleted_count_col;k++)
                            {
                                sum1 += (selected_col[k]/pow(selected_diff_col[k],p))*selected_roh[k];
                                sum2 += (1.0/pow(selected_diff_col[k],p))*selected_roh[k];
                            }
                            temp_col_shift[ori_index] = sum1/sum2;
                        }
                        else
                        {
                            //Average
                            for(int k=0;k<total_seleted_count_col;k++)
                            {
                                sum1 += selected_col[k]*selected_roh[k];
                                sum2 += selected_roh[k];
                            }
                            temp_col_shift[ori_index] = sum1/sum2;
                        }
                    }
                    else
                    {
                        //Average
                        for(int k=0;k<total_seleted_count_col;k++)
                        {
                            sum1 += selected_col[k]*selected_roh[k];
                            sum2 += selected_roh[k];
                        }
                        temp_col_shift[ori_index] = sum1/sum2;
                    }
                }
                
                //row
                if(total_seleted_count_row > 0)
                {
                    double sum1, sum2, sum1_r, sum2_r;;
                    
                    sum1 = 0;
                    sum2 = 0;
                    sum1_r = 0;
                    sum2_r = 0;
                    
                    if((int)pyramid_step >= 0)
                    {
                        if((max_hist_row > save_count*0.4))
                        {
                            //IDW
                            for(int k=0;k<total_seleted_count_row;k++)
                            {
                                sum1_r += (selected_row[k]/pow(selected_diff_row[k],p))*selected_roh[k];
                                sum2_r += (1.0/pow(selected_diff_row[k],p))*selected_roh[k];
                            }
                            temp_row_shift[ori_index] = sum1_r/sum2_r;
                            
                        }
                        else
                        {
                            //Average
                            for(int k=0;k<total_seleted_count_row;k++)
                            {
                                sum1_r += selected_row[k]*selected_roh[k];
                                sum2_r += selected_roh[k];
                            }
                            temp_row_shift[ori_index] = sum1_r/sum2_r;
                        }
                    }
                    else
                    {
                        //Average
                        for(int k=0;k<total_seleted_count_row;k++)
                        {
                            sum1_r += selected_row[k]*selected_roh[k];
                            sum2_r += selected_roh[k];
                        }
                        temp_row_shift[ori_index] = sum1_r/sum2_r;
                    }
                }
                selected_col.clear();
                selected_row.clear();
                selected_diff_col.clear();
                selected_diff_row.clear();
                selected_roh.clear();
                
            }
            save_col.clear();
            save_row.clear();
            save_diff.clear();
            save_slope.clear();
            save_roh.clear();
            free(hist_slope);
        }
    }
    
    printf("start assign\n");
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < data_length ; iter_count ++)
    {
        long r = (floor(iter_count/gridsize.width));
        long c = iter_count % gridsize.width;
        long ori_index = r*(long)gridsize.width + c;
    
        double sum_c = 0;
        double sum_r = 0;
        int count_cr = 0;
    
        for(long k = -kernal_size ; k <= kernal_size ; k++)
        {
            for(long j=-kernal_size ; j <= kernal_size ; j++)
            {
                long index = (r+k)*(long)gridsize.width + (c+j);
                if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width  && index >= 0 && index < data_length)
                {
                    sum_c += temp_col_shift[index];
                    sum_r += temp_row_shift[index];
                    count_cr ++;
                }
            }
        }
        
        GridPT3[ori_index].col_shift = sum_c/count_cr;
        GridPT3[ori_index].row_shift = sum_r/count_cr;
    }
    
    printf("before free col_shift\n");
    free(temp_col_shift);
    printf("free col_shift\n");
    free(temp_row_shift);
    printf("free row_shift\n");
    printf("end updating grid set shift \n");
}


void echo_print_nccresults_SDM(char *save_path,int row,int col,int level, int iteration, NCCresultSDM *nccresult, CSize *Size_Grid2D, char *add_str)
{
    int k,j;
    FILE *outfile_min, *outfile_max, *outfile_h, *outfile_roh, *outfile_diff, *outfile_peak, *outINCC, *outGNCC,*outcount;
    CSize temp_S;
    char t_str[500];
    
    sprintf(t_str,"%s/txt/nccresult_roh1_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_min    = fopen(t_str,"w");
    
    temp_S.height    = Size_Grid2D->height;
    temp_S.width    = Size_Grid2D->width;
    
    for(k=0;k<temp_S.height;k++)
    {
        for(j=0;j<temp_S.width;j++)
        {
            int matlab_index    = k*temp_S.width + j;
            
            fprintf(outfile_min,"%f\t",SignedCharToDouble_result(nccresult[matlab_index].result0));
        }
        fprintf(outfile_min,"\n");
    }
    
    fclose(outfile_min);
}

bool average_filter_colrowshift(CSize Size_Grid2D, UGRIDSDM *GridPT3,uint8 Pyramid_step)
{
    int kernel_size = 40;
    if(Pyramid_step == 3)
        kernel_size = 30;
    else if(Pyramid_step == 2)
        kernel_size = 10;
    else if(Pyramid_step == 1)
        kernel_size = 5;
    else if(Pyramid_step == 0)
        kernel_size = 5;
    
    
#pragma omp parallel for schedule(guided)
    for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
    {
        int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        int pts_col = iter_count % Size_Grid2D.width;
        int pt_index;
        pt_index = pts_row*Size_Grid2D.width + pts_col;
        
        if(pt_index < Size_Grid2D.width * Size_Grid2D.height && pts_row < Size_Grid2D.height && pts_col < Size_Grid2D.width && pts_row >= 0 && pts_col >= 0)
        {
            double sum_c = 0;
            double sum_r = 0;
            double weight = 0;
            int count_cr=0;
            
            if(GridPT3[pt_index].ortho_ncc < 0.4)
            {
                for(int k = -kernel_size ; k <= kernel_size ; k++)
                {
                    for(int j=-kernel_size ; j <= kernel_size ; j++)
                    {
                        int index = (pts_row+k)*Size_Grid2D.width + (pts_col+j);
                        
                        if(pts_row+k >= 0 && pts_row+k < Size_Grid2D.height && pts_col+j >= 0 && pts_col+j < Size_Grid2D.width && index >= 0 && index < Size_Grid2D.height*Size_Grid2D.width)
                        {
                            if(GridPT3[index].ortho_ncc > 0.5)
                            {
                                sum_c += GridPT3[index].col_shift*GridPT3[index].ortho_ncc;
                                sum_r += GridPT3[index].row_shift*GridPT3[index].ortho_ncc;
                                count_cr ++;
                                weight += GridPT3[index].ortho_ncc;
                            }
                        }
                    }
                }
                
                if(count_cr > 1)
                {
                    GridPT3[pt_index].col_shift = sum_c/weight;
                    GridPT3[pt_index].row_shift = sum_r/weight;
                }
            }
            
        }
    }
    
    
    return true;
}

double MergeTiles_SDM(ProInfo info,int iter_row_end,int t_col_end, int buffer,int final_iteration,TransParam _param, uint8 pyramid_step)
{
    FILE *poutDEM;
    FILE *poutMatchtag;
    FILE *poutheader;
    FILE *poutrowshift;
    FILE *poutcolshift;
    FILE *poutvxshift;
    FILE *poutvyshift;
    
    int header_line = pyramid_step+3;
    int row,col;
    int row_end = iter_row_end;
    int col_end = t_col_end;
    int find_level = 0;
    int find_iter  = final_iteration;
    long int size;
    double grid_size;
    
    CSize DEM_size, DEMinter_size;
    
    double boundary[4];
    
    //double *DEM, *DEMinter;
    //float *ColShift, *RowShift, *VxShift, *VyShift, *Roh;
    float *Mag, *VxShift, *VyShift, *Roh;
    //bool *Matchtag;
    char DEM_str[500];
    
    bool check_gs = false;
    
    //find boundary of DEM
    boundary[0] = 10000000.0;
    boundary[1] = 10000000.0;
    boundary[2] = -10000000.0;
    boundary[3] = -10000000.0;
    
    for(long index_file = 0 ; index_file < row_end*col_end ; index_file++)
    {
        int row,col;
        
        row = (int)(floor(index_file/col_end));
        col = index_file%col_end;
        
        FILE *pfile;
        char t_str[500];
        
        sprintf(t_str,"%s/txt/tin_xshift_level_1_1_0_iter_3_final.txt",info.save_filepath);
        printf("%s\n",t_str);
        //sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info.save_filepath,row,col,find_level,find_iter);
        //sprintf(t_str,"%s/txt/col_shift_%d_%d_%d_%d.txt",info.save_filepath,find_level,find_iter,row,col);
        pfile    = fopen(t_str,"r");
        if(pfile)
        {
            fseek(pfile,0,SEEK_END);
            size = ftell(pfile);
            if(size > 0)
            {
                char h_t_str[500];
                FILE *p_hfile;
                //sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col,find_level,find_iter);
                sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col);
                
                p_hfile        = fopen(h_t_str,"r");
                if(p_hfile)
                {
                    printf("%s\n",h_t_str);
                    int t_row,t_col,t_level,t_col_size,t_row_size;
                    double t_grid_size;
                    double t_boundary[4];
                    while(!feof(p_hfile))
                    {
                        fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n", &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&t_col_size,&t_row_size);
                    }
                    
                    grid_size = t_grid_size;
                    t_boundary[2] = t_boundary[0] + t_grid_size*t_col_size;
                    t_boundary[3] = t_boundary[1] + t_grid_size*t_row_size;
                    
                    {
                        if(boundary[0] > t_boundary[0])
                            boundary[0]        = t_boundary[0];
                        if(boundary[1] > t_boundary[1])
                            boundary[1]        = t_boundary[1];
                        
                        if(boundary[2] < t_boundary[2])
                            boundary[2]        = t_boundary[2];
                        if(boundary[3] < t_boundary[3])
                            boundary[3]        = t_boundary[3];
                    }
                    
                    fclose(p_hfile);
                }
            }
            fclose(pfile);
        }
    }
    
    printf("boundary %f\t%f\t%f\t%f\n",boundary[0],boundary[1],boundary[2],boundary[3]);
    check_gs = false;
    
    buffer    = (int)(floor(buffer/grid_size));
    
    DEM_size.width        = (int)(ceil( (double)(boundary[2] - boundary[0]) /grid_size ));
    DEM_size.height        = (int)(ceil( (double)(boundary[3] - boundary[1]) /grid_size ));
    
    //ColShift = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    //RowShift = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    Mag = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    VxShift = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    VyShift = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    Roh = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    
    float *Vxsigma = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    float *Vysigma = (float*)calloc(DEM_size.height*DEM_size.width,sizeof(float));
    
    printf("dem size %d\t%d\n",DEM_size.width,DEM_size.height);
    
    //setting DEM value
    for(long index_file = 0 ; index_file < row_end*col_end ; index_file++)
    {
        int row,col;
        FILE *pfile;
        char t_str[500];
        
        row = (int)(floor(index_file/col_end));
        col = index_file%col_end;
        
        sprintf(t_str,"%s/txt/tin_xshift_level_1_1_0_iter_3_final.txt",info.save_filepath);
        printf("%s\n",t_str);
        //sprintf(t_str,"%s/txt/col_shift_%d_%d_%d_%d.txt",info.save_filepath,find_level,find_iter,row,col);
        pfile    = fopen(t_str,"r");
        if(pfile)
        {
            fseek(pfile,0,SEEK_END);
            size = ftell(pfile);
            fseek(pfile,0L,SEEK_SET);
            if(size > 0)
            {
                char h_t_str[500];
                FILE *p_hfile, *p_hvfile, *p_cshift, *p_rshift, *p_vxshift, *p_vyshift, *p_roh, *p_mag;
                
                sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col);
                p_hfile        = fopen(h_t_str,"r");
                printf("%s\n",h_t_str);
                //if(p_hfile)
                {
                    int iter;
                    char hv_t_str[500];
                    int row_size,col_size;
                    double t_boundary[4];
                    int t_row,t_col,t_level;
                    double t_grid_size;
                    
                    while(!feof(p_hfile))
                    {
                        fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n",  &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
                    }
                    
                    
                    sprintf(hv_t_str,"%s/txt/tin_xshift_level_%d_%d_%d_iter_3_final.txt",info.save_filepath,row,col,find_level);
                    p_vxshift    = fopen(hv_t_str,"r");
                    sprintf(hv_t_str,"%s/txt/tin_yshift_level_%d_%d_%d_iter_3_final.txt",info.save_filepath,row,col,find_level);
                    p_vyshift    = fopen(hv_t_str,"r");
                    /*
                    sprintf(hv_t_str,"%s/txt/tin_xsigma_level_%d_%d_%d_iter_3_final.txt",info.save_filepath,row,col,find_level);
                    p_cshift    = fopen(hv_t_str,"r");
                    sprintf(hv_t_str,"%s/txt/tin_ysigma_level_%d_%d_%d_iter_3_final.txt",info.save_filepath,row,col,find_level);
                    p_rshift    = fopen(hv_t_str,"r");
                    */
                    sprintf(hv_t_str,"%s/txt/tin_roh_level_%d_%d_%d_iter_3_final.txt",info.save_filepath,row,col,find_level);
                    p_roh    = fopen(hv_t_str,"r");
                    //if(p_hvfile)
                    {
                        int index_total;
                        for(index_total = 0; index_total < row_size*col_size ; index_total++)
                        {
                            int iter_row,iter_col;
                            iter_row = floor(index_total/col_size);
                            iter_col = index_total%col_size;
                            {
                                
                                double t_col = ( (double)(t_boundary[0] + grid_size*iter_col - boundary[0])  /grid_size);
                                double t_row = ( (double)(boundary[3] - (t_boundary[1] + grid_size*iter_row))/grid_size);
                                int index = (int)(t_row*DEM_size.width + t_col + 0.01);
                                
                                double DEM_value;
                                double Col_value, Row_value, Vx_value, Vy_value, Roh_value;
                                
                                //fscanf(p_cshift,"%lf\t",&Col_value);
                                //fscanf(p_rshift,"%lf\t",&Row_value);
                                
                                fscanf(p_vxshift,"%lf\t",&Vx_value);
                                fscanf(p_vyshift,"%lf\t",&Vy_value);
                                fscanf(p_roh,"%lf\t",&Roh_value);
                                
                                
                                if(index >= 0 && index < DEM_size.width*DEM_size.height &&
                                   iter_row > buffer && iter_row < row_size - buffer &&
                                   iter_col > buffer && iter_col < col_size - buffer)
                                {
                                    
                                    //if(fabs(Col_value) < 1000)
                                    //    Vxsigma[index] = Col_value;
                                    //if(fabs(Row_value) < 1000)
                                     //   Vysigma[index] = Row_value;
                                    
                                    //if(fabs(Vx_value) < 1000)
                                        VxShift[index] = Vx_value;
                                    //if(fabs(Vy_value) < 1000)
                                        VyShift[index] = Vy_value;
                                    Roh[index] = Roh_value;
                                    Mag[index] = sqrt(Vx_value*Vx_value + Vy_value*Vy_value);
                                }
                            }
                            //fscanf(p_cshift,"\n");
                            //fscanf(p_rshift,"\n");
                            fscanf(p_vxshift,"\n");
                            fscanf(p_vyshift,"\n");
                            fscanf(p_roh,"\n");
                        }
                 
                        //fclose(p_cshift);
                        //fclose(p_rshift);
                        fclose(p_vxshift);
                        fclose(p_vyshift);
                        fclose(p_roh);
                    }
                    
                    fclose(p_hfile);
                }
            }
            fclose(pfile);
        }
    }
    /*
    sprintf(DEM_str, "%s/%s_dx_sigma.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, Vxsigma, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_dy_sigma.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, Vysigma, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
     */
    sprintf(DEM_str, "%s/%s_dx.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, VxShift, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_dy.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, VyShift, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_dmag.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, Mag, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
    sprintf(DEM_str, "%s/%s_roh.tif", info.save_filepath, info.Outputpath_name);
    WriteGeotiff(DEM_str, Roh, DEM_size.width, DEM_size.height, grid_size, boundary[0], boundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 4);
    
    //free(Vxsigma);
    //free(Vysigma);
    
    free(VxShift);
    free(VyShift);
    free(Mag);
    free(Roh);
    
    return grid_size;
}
