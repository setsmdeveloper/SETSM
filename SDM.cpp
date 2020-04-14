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
                int th_grid = proinfo.resolution*pow(2.0,end_level);
                
                if(end_level < 0)
                    end_level = 0;
                if(end_level > 4)
                    end_level = 4;
                
                proinfo.end_level = end_level;
                proinfo.pyramid_level = args.pyramid_level;
                proinfo.number_of_images = 2;
                
                printf("pyramid level %d\tSDM_ss %d\tend_level = %d\t%d\n",proinfo.pyramid_level,proinfo.SDM_SS,end_level,th_grid);
                
                int sdm_kernal_size = floor( (double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution*pow(2.0,proinfo.pyramid_level)));
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
                            sdm_kernal_size = floor( (double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution*pow(2.0,proinfo.pyramid_level)));
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
                Matching_SETSM_SDM(proinfo, Template_size, Rimageparam, Limagesize, Rimagesize, Boundary, GSD_image1, GSD_image2, &matching_number);
                
                int final_iteration = 3;
                if(matching_number > 10)
                {
                    printf("Tile merging start final iteration %d!!\n",final_iteration);
                    int buffer_tile = 0;
                    double mt_grid_size = MergeTiles_SDM(proinfo,1,1,buffer_tile,final_iteration,param,pyramid_step);
                }
            
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



void Matching_SETSM_SDM(ProInfo proinfo, uint8 Template_size, double *Rimageparam, const CSize Limagesize, const CSize Rimagesize, const double *Boundary, const ImageGSD gsd_image1, const ImageGSD gsd_image2, int *matching_number)
{
    bool check_cal = false;
    char check_file[500];
    sprintf(check_file,"%s/txt/col_shift_0_3_%d_%d.txt",proinfo.save_filepath,1,1);
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
        filename = GetFileName(proinfo.Imagefilename[1]);
        filename = remove_ext(filename);
        sprintf(Rsubsetfilename,"%s/%s_subset_%d_%d.raw",proinfo.tmpdir,filename,1,1);
        
        printf("subsetimage level %d\n",proinfo.pyramid_level);
        
        LevelInfo levelinfo = {NULL};
        levelinfo.Boundary = subBoundary;
        levelinfo.Template_size = &Template_size;
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
            
            if( subsetsize[0].height > Template_size/2*pwrtwo(pyramid_step) && subsetsize[0].width > Template_size/2*pwrtwo(pyramid_step) &&
               subsetsize[1].height > Template_size/2*pwrtwo(pyramid_step) && subsetsize[1].width > Template_size/2*pwrtwo(pyramid_step) )
            {
                int count_tri;
                
                
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
                
                //Preprocessing_SDM(proinfo, proinfo.tmpdir, Sourceimage, level, &Lsubsetsize, &Rsubsetsize, data_size_l, data_size_r, fid);
                
                
                //PreST = time(0);
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
                    
                    if(proinfo.pre_DEMtif)
                        Template_size = 7;
                        
                    //Template_size = 15;
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
                    levelinfo.Pyramid_step = &level;
                    levelinfo.py_Sizes = data_size_lr;
                    
                    double Th_roh, Th_roh_min, Th_roh_start, Th_roh_next;
                    SetThs_SDM(level,final_level_iteration, &Th_roh, &Th_roh_min, &Th_roh_next, &Th_roh_start,pyramid_step);
                    
                    double py_resolution = 0;
                    double grid_resolution = 0;
                    
                    GridPT    = SetGrids_SDM(proinfo,prc_level,level,pyramid_step, final_level_iteration, &Size_Grid2D, &py_resolution, &grid_resolution, subBoundary);
                    
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
                        
                        BL blunder_param;
                        
                        char filename_mps[500];
                        
                        char filename_tri[500];
                        char v_temp_path[500];
                        
                        int count_results[2];
                        int count_results_anchor[2];
                        int count_MPs;
                        int count_blunder;
                        
                        fprintf(fid,"Starting computation of NCC\n iteration = %u\tTh_roh = %f\tTh_roh_start = %f\tGrid size %d %d\n",
                                iteration, Th_roh,Th_roh_start,Size_Grid2D.width,Size_Grid2D.height);
                        
                        //printf("sub size %d\t%d\t%d\t%d\n",data_size_l[level].width,data_size_l[level].height,data_size_r[level].width,data_size_r[level].height);
                        
                        sprintf(filename_mps,"%s/txt/matched_pts_%d_%d_%d_%d.txt",proinfo.save_filepath,1,1,level,iteration);
                        
                        sprintf(v_temp_path,"%s/tmp/vv_tmp",proinfo.save_filepath);
                        
                        printf("template size =%d\n",Template_size);
                        
                        //echoprint_shift(Lstartpos, Rstartpos,subBoundary,LRPCs, RRPCs, t_Rimageparam, param,GridPT,proinfo.save_filepath,row,col,level,iteration,update_flag,&Size_Grid2D,GridPT3,"final");
                        bool check_smooth = false;
                        if(grid_resolution >= 200 && level == 0)
                            check_smooth = true;
                        else if(level < 2 && grid_resolution < 300)
                            check_smooth = true;
                        
                        if(check_smooth)
                        {
                            bool check_while = false;
                            bool check_last_iter = false;
                            int check_size = 3;//*(int)(200/grid_resolution);
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
                            
                            int count_null_cell = 0;
                            
                            uint8* t_ncc_array = (uint8*)calloc(sizeof(uint8),(long)Size_Grid2D.width*(long)Size_Grid2D.height);
                            float* temp_col_shift = (float*)calloc(sizeof(float),(long)Size_Grid2D.width*(long)Size_Grid2D.height);
                            float* temp_row_shift = (float*)calloc(sizeof(float),(long)Size_Grid2D.width*(long)Size_Grid2D.height);
                            
                            while(check_while == 0 && sm_iter < 20)
                            {
                                
                                #pragma omp parallel for schedule(guided)
                                for (long index = 0; index < (long)Size_Grid2D.width*(long)Size_Grid2D.height; index++)
                                {
                                    int row, col;
                                    
                                    row = (int)(floor(index/Size_Grid2D.width));
                                    col = index%Size_Grid2D.width;
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
                                        //printf("sm last iteration %d\n",sm_iter);
                                    }
                                    else
                                    {
                                        if (/*GridPT3[search_index].ortho_ncc < 0.5 &&*/ t_ncc_array[search_index] != 2)
                                        {
                                            t_ncc_array[search_index] = 1;
                                        }
                                    }
                                }
                                
                                count_null_cell = 0;
                                sm_iter++;
                                
                                #pragma omp parallel for schedule(guided) reduction(+:count_null_cell)
                                for (long index = 0; index < (long)Size_Grid2D.width*(long)Size_Grid2D.height; index++)
                                {
                                    int row, col;
                                    
                                    long count_cell;
                                    int t_i, t_j;
                                    
                                    double sum_row_shift=0;
                                    double sum_col_shift=0;
                                    double sum_roh = 0;
                                    
                                    row = (int)(floor(index/Size_Grid2D.width));
                                    col = index%Size_Grid2D.width;
                                    long search_index = (long)row*(long)Size_Grid2D.width + (long)col;
                                    double p = 1.5;
                                    if (t_ncc_array[search_index] == 1)
                                    {
                                        long null_count_cell = 0;
                                        count_cell = 0;
                                        
                                        for (t_i = -check_size; t_i <= check_size;t_i++ )
                                        {
                                            for (t_j = -check_size; t_j <= check_size; t_j++)
                                            {
                                                int index_row = row + t_i;
                                                int index_col = col + t_j;
                                                long int t_index = (long)index_row*(long)Size_Grid2D.width + (long)index_col;
                                                
                                                if(index_row >= 0 && index_row < Size_Grid2D.height && index_col >= 0 && index_col < Size_Grid2D.width && t_i != 0 && t_j != 0)
                                                {
                                                    if(GridPT3[t_index].ortho_ncc > 0.0)
                                                    {
                                                        //double dist = sqrt((double)(t_i*t_i + t_j*t_j));
                                                        double ncc = 1 + GridPT3[t_index].ortho_ncc;
                                                        double weith_n = pow(ncc,10);
                                                        sum_row_shift += (GridPT3[t_index].row_shift*weith_n);
                                                        sum_col_shift += (GridPT3[t_index].col_shift*weith_n);
                                                        sum_roh += weith_n;
                                                        //sum_row_shift += (GridPT3[t_index].row_shift/pow(dist,p)*GridPT3[t_index].ortho_ncc);
                                                        //sum_col_shift += (GridPT3[t_index].col_shift/pow(dist,p)*GridPT3[t_index].ortho_ncc);
                                                        //sum_roh += ((1.0/pow(dist,p))*GridPT3[t_index].ortho_ncc);
                                                        count_cell++;
                                                    }
                                                    
                                                    if(GridPT3[t_index].ortho_ncc < 0.3)
                                                        null_count_cell ++;
                                                }
                                            }
                                        }
                                        
                                        if (count_cell >= total_size*0.3 )
                                        {
                                            //GridPT3[search_index].row_shift = sum_row_shift/(double)count_cell;
                                            //GridPT3[search_index].col_shift = sum_col_shift/(double)count_cell;
                                            //GridPT3[search_index].row_shift = sum_row_shift/sum_roh;
                                            //GridPT3[search_index].col_shift = sum_col_shift/sum_roh;
                                            
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
                                
                                Update_ortho_NCC(proinfo, SubMagImages_L,SubMagImages_R,grid_resolution, proinfo.resolution, Limagesize, data_size_l[prc_level], SubImages_L, Rimagesize,data_size_r[prc_level],SubImages_R, Template_size,
                                                 Size_Grid2D, GridPT, GridPT3, prc_level, Lstartpos, Rstartpos, iteration, subBoundary, gsd_image1, gsd_image2,Rimageparam,SubOriImages_L,SubOriImages_R);
                            }
                            free(t_ncc_array);
                            free(temp_col_shift);
                            free(temp_row_shift);
                        }
                        
                        if(level < pyramid_step)
                            proinfo.SDM_SS = 3;
                        printf("kernel size %d\n",proinfo.SDM_SS);
                        VerticalLineLocus_SDM(proinfo, nccresult, SubMagImages_L,SubMagImages_R,grid_resolution, proinfo.resolution, Limagesize, data_size_l[prc_level], SubImages_L, Rimagesize,data_size_r[prc_level],SubImages_R, Template_size,
                                          Size_Grid2D, GridPT, GridPT3, prc_level, pyramid_step, Lstartpos, Rstartpos, iteration, subBoundary, gsd_image1, gsd_image2,Rimageparam,SubOriImages_L,SubOriImages_R);
                        
                        count_MPs = SelectMPs_SDM(proinfo, nccresult,Size_Grid2D,GridPT,GridPT3,level,prc_level,
                                              iteration,filename_mps,1,1,subBoundary);
                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\n",1,1,level,iteration,count_MPs);
                        
                        //echo_print_nccresults(proinfo.save_filepath,row,col,level, iteration, nccresult, &Size_Grid2D, "final");
                        
                        if(count_MPs > 100)
                            lower_level_match = true;
                        
                        if(lower_level_match)
                        {
                            FILE* survey;
                            int i;
                            
                            survey    = fopen(filename_mps,"r");
                            
                            blunder_param.Boundary    = subBoundary;
                            blunder_param.gridspace    = grid_resolution;
                            blunder_param.height_check_flag = true;
                            blunder_param.iteration = iteration;
                            blunder_param.Pyramid_step = prc_level;
                            
                            blunder_param.Size_Grid2D.width = Size_Grid2D.width;
                            blunder_param.Size_Grid2D.height = Size_Grid2D.height;
                            
                            if(level == 0 && iteration == 3)
                            {
                                //echo_print_nccresults_SDM(proinfo.save_filepath,1,1,level, iteration, nccresult, &Size_Grid2D, "final");
                                const char *temp_str = {"final"};
                                echoprint_Gridinfo_SDM(prc_level,proinfo,proinfo.LBoundary,proinfo.RBoundary, data_size_l[prc_level], SubImages_L, data_size_r[prc_level], SubImages_R,subBoundary,proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,GridPT3,temp_str);
                                
                                echoprint_adjustXYZ(prc_level, proinfo.LBoundary, proinfo.RBoundary,data_size_l[prc_level], SubImages_L, data_size_r[prc_level], SubImages_R,subBoundary,proinfo, proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,GridPT3,temp_str,grid_resolution,1);
                            }
                            else
                            {
                                char bufstr[500];
                                
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
                                        if(proinfo.IsRA)
                                            Th_roh_update        = (double)(Th_roh - 0.10);
                                        else
                                            Th_roh_update        = (double)(Th_roh - 0.06);
                                    }
                                }
                                
                                if(level == 0)
                                {
                                    Pre_GridPT3        = SetHeightRange_SDM(GridPT3,blunder_param);
                                    
                                    printf("update Pre_GridPT3\n");
                                }
                                else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
                                {
                                    Pre_GridPT3        = SetHeightRange_SDM(GridPT3,blunder_param);
                                    printf("update Pre_GridPT3\n");
                                    
                                    
                                    check_pre_GridPT3 = true;
                                }
                                else
                                {
                                    GridPT3        = SetHeightRange_SDM(GridPT3,blunder_param);
                                    printf("update GridPT3\n");
                                }
                                
                                //displacement computation
                                if(level <= pyramid_step)
                                {
                                    int count_shift = 0;
                                    
                                    count_shift = count_MPs;
                                    
                                    printf("Size_Grid2D %d\t%d\t%d\n",Size_Grid2D.width,Size_Grid2D.height,count_shift);
                                    
                                    printf("TIN interpolation for col row shift %d\n",count_shift);
                                   
                                    if(count_shift > 5)
                                    {
                                        
                                        char bufstr[500];
                                        D3DPOINT *ptslists;
                                        FILE *fid_load;
                                        
                                        //column
                                        sprintf(bufstr,"%s/txt/col_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,1,1);
                                        fid_load = fopen(bufstr,"r");
                                        int i;
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_shift);
                                        
                                        i = 0;
                                        
                                        while( i < count_shift && (fscanf(fid_load,"%lf %lf %lf\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z)) != EOF )
                                        {
                                            i++;
                                        }
                                        
                                        fclose(fid_load);
                                        remove(bufstr);
                                        
                                        UI3DPOINT *trilists;
                                        double min_max[4] = {subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]};
                                        
                                        UI3DPOINT* t_trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_shift*4);
                                        
                                        sprintf(bufstr,"%s/txt/col_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,1,1);
                                        TINCreate(ptslists,count_shift,t_trilists,min_max,&count_tri,grid_resolution);
                                        printf("end tingeneration %d\n",count_tri);
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        i = 0;
                                        for(i=0;i<count_tri;i++)
                                        {
                                            trilists[i].m_X = t_trilists[i].m_X;
                                            trilists[i].m_Y = t_trilists[i].m_Y;
                                            trilists[i].m_Z = t_trilists[i].m_Z;
                                        }
                                        
                                        free(t_trilists);
                                        
                                        printf("end tingeneration %d\n",count_tri);
                                        
                                        if(level == 0)
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, Pre_GridPT3,blunder_param,ptslists,trilists,0);
                                            //echoprint_Gridinfo(proinfo.save_filepath,row,col,level,iteration,update_flag,&Size_Grid2D,Pre_GridPT3,"final");
                                        }
                                        else if(!check_pre_GridPT3)
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, GridPT3,blunder_param,ptslists,trilists,0);
                                            //echoprint_Gridinfo(proinfo.save_filepath,row,col,level,iteration,update_flag,&Size_Grid2D,GridPT3,"final");
                                        }
                                        else
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, Pre_GridPT3,blunder_param,ptslists,trilists,0);
                                            
                                            //echoprint_Gridinfo(proinfo.save_filepath,row,col,level,iteration,update_flag,&Size_Grid2D,Pre_GridPT3,"final");
                                        }
                                        
                                        printf("end col computation %d\n",check_pre_GridPT3);
                                        
                                        free(trilists);
                                        free(ptslists);
                                        
                                        
                                        //row
                                        sprintf(bufstr,"%s/txt/row_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,1,1);
                                        fid_load = fopen(bufstr,"r");
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_shift);
                                        
                                        i = 0;
                                        while( i < count_shift && (fscanf(fid_load,"%lf %lf %lf\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z)) != EOF )
                                        {
                                            i++;
                                        }
                                        
                                        fclose(fid_load);
                                        remove(bufstr);
                                        
                                        t_trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_shift*4);
                                        
                                        sprintf(bufstr,"%s/txt/row_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,1,1);
                                        TINCreate(ptslists,count_shift,t_trilists,min_max,&count_tri,grid_resolution);
                                        
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        i = 0;
                                        for(i=0;i<count_tri;i++)
                                        {
                                            trilists[i].m_X = t_trilists[i].m_X;
                                            trilists[i].m_Y = t_trilists[i].m_Y;
                                            trilists[i].m_Z = t_trilists[i].m_Z;
                                        }
                                        
                                        free(t_trilists);
                                        
                                        printf("end tingeneration %d\n",count_tri);
                                        
                                        double* ortho_ncc = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
                                        double* INCC = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
                                        
                                        if(level == 0)
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, Pre_GridPT3,blunder_param,ptslists,trilists,1);
                                            
                                            Update_ortho_NCC(proinfo, SubMagImages_L,SubMagImages_R,grid_resolution, proinfo.resolution, Limagesize, data_size_l[prc_level], SubImages_L, Rimagesize,data_size_r[prc_level],SubImages_R, Template_size,
                                                                  Size_Grid2D, GridPT, Pre_GridPT3, prc_level, Lstartpos, Rstartpos, iteration, subBoundary, gsd_image1, gsd_image2,Rimageparam,SubOriImages_L,SubOriImages_R);
                                            //if(level <= 4)
                                            {
                                                if(iteration < 2)
                                                {
                                                    //average_filter_colrowshift(Size_Grid2D, Pre_GridPT3,level);
                                                    blunder_param.Pyramid_step = level;
                                                    shift_filtering(proinfo, Pre_GridPT3, blunder_param, proinfo.DEM_resolution);
                                                }
                                            }
                                        //echoprint_Gridinfo_SDM(proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,Pre_GridPT3,"final");
                                            
                                            //echoprint_adjustXYZ(proinfo, proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,Pre_GridPT3,"final",grid_resolution,1);
                                            
                                            check_pre_GridPT3 = false;
                                        }
                                        else if(!check_pre_GridPT3)
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, GridPT3,blunder_param,ptslists,trilists,1);
                                            
                                            Update_ortho_NCC(proinfo, SubMagImages_L,SubMagImages_R,grid_resolution, proinfo.resolution, Limagesize, data_size_l[prc_level], SubImages_L, Rimagesize,data_size_r[prc_level],SubImages_R, Template_size,
                                                                  Size_Grid2D, GridPT, GridPT3, prc_level, Lstartpos, Rstartpos, iteration, subBoundary, gsd_image1, gsd_image2,Rimageparam,SubOriImages_L,SubOriImages_R);
                                            
                                            //if(level <= 4)
                                            {
                                                if(level >= 2 || (level == 1 && iteration < 3))
                                                {
                                                    //average_filter_colrowshift(Size_Grid2D, GridPT3,level);
                                                    blunder_param.Pyramid_step = level;
                                                    shift_filtering(proinfo, GridPT3, blunder_param, proinfo.DEM_resolution);
                                                }
                                            }
                                            
                                        //echoprint_Gridinfo_SDM(proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,GridPT3,"final");
                                            
                                            //echoprint_adjustXYZ(proinfo, proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,GridPT3,"final",grid_resolution,1);
                                        }
                                        else
                                        {
                                            SetHeightRange_shift(count_shift, count_tri, Pre_GridPT3,blunder_param,ptslists,trilists,1);
                                            
                                            Update_ortho_NCC(proinfo, SubMagImages_L,SubMagImages_R,grid_resolution, proinfo.resolution, Limagesize, data_size_l[prc_level], SubImages_L, Rimagesize,data_size_r[prc_level],SubImages_R, Template_size,
                                                                  Size_Grid2D, GridPT, Pre_GridPT3, prc_level, Lstartpos, Rstartpos, iteration, subBoundary, gsd_image1, gsd_image2,Rimageparam,SubOriImages_L,SubOriImages_R);
                                            //if(level <= 4)
                                            {
                                                if(level >= 2 || (level == 1 && iteration < 3))
                                                {
                                                    //average_filter_colrowshift(Size_Grid2D, Pre_GridPT3,level);
                                                    blunder_param.Pyramid_step = level;
                                                    shift_filtering(proinfo, Pre_GridPT3, blunder_param, proinfo.DEM_resolution);
                                                    
                                                }
                                            }
                                            
                                            
                                        //echoprint_Gridinfo_SDM(proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,Pre_GridPT3,"final");
                                            
                                            //echoprint_adjustXYZ(proinfo, proinfo.save_filepath,1,1,level,iteration,&Size_Grid2D,Pre_GridPT3,"final",grid_resolution,1);
                                            
                                            check_pre_GridPT3 = false;
                                        }
                                        
                                        free(ortho_ncc);
                                        free(INCC);
                                        
                                        free(trilists);
                                        free(ptslists);
                                        
                                        printf("end compute\n");
                                    }
                                    
                                    if(level == 0 && iteration == 2)
                                    {
                                        
                                    }
                                    else
                                    {
                                        /*
                                         char bufstr[500];
                                         sprintf(bufstr,"%s/txt/col_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,row,col);
                                         remove(bufstr);
                                         sprintf(bufstr,"%s/txt/row_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,level,iteration,row,col);
                                         remove(bufstr);
                                         */
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
                        {
                            *matching_number = count_MPs;
                        }
                        else
                        {
                            remove(filename_mps);
                        }
                        
                        
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
                        {
                            if(proinfo.IsRA)
                                Th_roh            = (double)(Th_roh - 0.10);
                            else
                                Th_roh            = (double)(Th_roh - 0.06);
                        }
                        
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
                                {
                                    if(proinfo.IsRA)
                                        Th_roh            = (double)(Th_roh + 0.10);
                                    else
                                        Th_roh            = (double)(Th_roh + 0.06);
                                }
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
                    /*
                    printf("release subImage L\n");
                    free(SubImages_L);
                    
                    printf("release subImage R\n");
                    free(SubImages_R);
                    
                    printf("release subimage Mag\n\n\n");
                    free(SubMagImages_L);
                    free(SubMagImages_R);
                    
                    free(SubOriImages_L);
                    free(SubOriImages_R);
                     */
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
  
        RemoveFiles_SDM(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,proinfo.pyramid_level,0);
        printf("done remove tempfiles\n");
    }
}

void SetHeightWithSeedDEM_SDM(ProInfo proinfo, UGRIDSDM *Grid, double *Boundary, CSize Grid_size, double Grid_set)
{
    double minX = 0, maxX = 0, minY = 0, maxY = 0, grid_size = 0, a_minX = 0, a_maxX = 0, a_minY = 0, a_maxY = 0;
    CSize seeddem_size;
    char* hdr_path;
    FILE *bin;
    TIFF *tif;
    char save_DEMfile[500];
    char *GIMP_path = proinfo.priori_DEM_tif;
    char row_shift_path[500];
    char col_shift_path[500];
    
    char* temp_outpathname = SetOutpathName(GIMP_path);
    sprintf(col_shift_path,"%s/%s_dx.tif",GIMP_path,temp_outpathname);
    sprintf(row_shift_path,"%s/%s_dy.tif",GIMP_path,temp_outpathname);
    
    printf("dx file %s\n",col_shift_path);
    printf("dy file %s\n",row_shift_path);
    int check_ftype = 1; // 1 = tif, 2 = raw
    
    seeddem_size = ReadGeotiff_info(col_shift_path, &minX, &maxY, &grid_size);
    
    maxX    = minX + grid_size*((double)seeddem_size.width);
    minY    = maxY - grid_size*((double)seeddem_size.height);
    
    printf("SDM_days %f\n",proinfo.SDM_days);
    printf("%d\n",seeddem_size.width);
    printf("%d\n",seeddem_size.height);
    printf("%f\n",minX);
    printf("%f\n",minY);
    printf("%f\n",maxX);
    printf("%f\n",maxY);
    printf("%f\n",grid_size);

    if(minX > (double)Boundary[0])
        a_minX = minX;
    else {
        a_minX = (double)Boundary[0];
    }
    if (maxX < (double)Boundary[2]) {
        a_maxX = maxX;
    }
    else {
        a_maxX = (double)Boundary[2];
    }

    if(minY > (double)Boundary[1])
        a_minY = minY;
    else {
        a_minY = (double)Boundary[1];
    }
    if (maxY < (double)Boundary[3]) {
        a_maxY = maxY;
    }
    else {
        a_maxY = (double)Boundary[3];
    }

    printf("%f %f %f %f\n",(double)Boundary[0],(double)Boundary[1],(double)Boundary[2],(double)Boundary[3]);
    printf("%f %f %f %f\n",a_minX, a_maxX, a_minY, a_maxY);
    if ( (a_minX < a_maxX) && (a_minY < a_maxY))
    {
        int row,col;
        
        float *seeddx = NULL;
        float *seeddy = NULL;
        long int cols[2], rows[2];
        CSize data_size;

        CSize *LImagesize = (CSize*)malloc(sizeof(CSize));
        LImagesize->width = seeddem_size.width;
        LImagesize->height = seeddem_size.height;
        
        cols[0] = 0;
        cols[1] = seeddem_size.width;
        
        rows[0] = 0;
        rows[1] = seeddem_size.height;
        
        float type(0);
        seeddx = Readtiff_T(col_shift_path,LImagesize,cols,rows,&data_size,type);
        seeddy = Readtiff_T(row_shift_path,LImagesize,cols,rows,&data_size,type);
        printf("Grid size %d\t%d\tcols rows %d\t%d\t%d\t%d\n",Grid_size.width,Grid_size.height,cols[0],cols[1],rows[0],rows[1]);
        
        for (row = 0; row < Grid_size.height; row ++) {
            for (col = 0; col < Grid_size.width; col ++) {
                int index_grid = row*Grid_size.width + col;
                double t_x,t_y;
                int index_seeddem;
                int row_seed, col_seed;
                
                t_x = Boundary[0] + col*Grid_set;
                t_y = Boundary[1] + row*Grid_set;
                col_seed = floor((t_x - minX)/grid_size);// - cols[0];
                row_seed = floor((maxY - t_y)/grid_size);// - rows[0];
                
                index_seeddem = row_seed*data_size.width + col_seed;
                if(index_seeddem >= 0 && index_seeddem < data_size.width*data_size.height)
                {
                    if(seeddx[index_seeddem] > -1000 && seeddy[index_seeddem] > -1000)
                    {
                        Grid[index_grid].col_shift = seeddx[index_seeddem]*proinfo.SDM_days/grid_size;
                        Grid[index_grid].row_shift = -seeddy[index_seeddem]*proinfo.SDM_days/grid_size;
                        
                        //printf("col row shift %f\t%f\tdx dy %f\t%f\n",Grid[index_grid].col_shift,Grid[index_grid].row_shift,seeddx[index_seeddem],seeddy[index_seeddem]);
                    }
                }
                
            }
        }
        
        free(seeddx);
        free(seeddy);
        printf("seeddem end\n");
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
/*
bool GetsubareaImage_SDM(ProInfo proinfo, char *ImageFilename, CSize *Imagesize, double *subBoundary, long int *cols,long int *rows)
{
    bool ret = false;
    
    if(GetImageSize(ImageFilename,Imagesize))
    {
        int i;
        
        D2DPOINT t_pts[8];
        D2DPOINT *ImageCoord;
        int buffer, null_buffer;
        
        double minX =  1000000;
        double maxX = -1000000;
        double minY =  1000000;
        double maxY = -1000000;
        
        t_pts[0].m_X    = subBoundary[0];
        t_pts[1].m_X    = subBoundary[2];
        t_pts[2].m_X    = subBoundary[0];
        t_pts[3].m_X    = subBoundary[2];
        t_pts[4].m_X    = subBoundary[0];
        t_pts[5].m_X    = subBoundary[2];
        t_pts[6].m_X    = subBoundary[0];
        t_pts[7].m_X    = subBoundary[2];
        
        t_pts[0].m_Y    = subBoundary[1];
        t_pts[1].m_Y    = subBoundary[3];
        t_pts[2].m_Y    = subBoundary[1];
        t_pts[3].m_Y    = subBoundary[3];
        t_pts[4].m_Y    = subBoundary[3];
        t_pts[5].m_Y    = subBoundary[1];
        t_pts[6].m_Y    = subBoundary[3];
        t_pts[7].m_Y    = subBoundary[1];
        
        ImageCoord        = GetObjectToImage(8, t_pts,subBoundary,proinfo.resolution);
        
        for(i=0;i<8;i++)
        {
            if(minX > ImageCoord[i].m_X)
                minX    = ImageCoord[i].m_X;
            if(maxX < ImageCoord[i].m_X)
                maxX    = ImageCoord[i].m_X;
            if(minY > ImageCoord[i].m_Y)
                minY    = ImageCoord[i].m_Y;
            if(maxY < ImageCoord[i].m_Y)
                maxY    = ImageCoord[i].m_Y;
            
            //printf("i %d\tImageCoord %f\t%f\n",i,ImageCoord[i].m_X,ImageCoord[i].m_Y);
        }
        
        buffer                = 0;
        cols[0]                = (int)(ceil(minX)-buffer);
        cols[1]                = (int)(ceil(maxX)+buffer);
        rows[0]                = (int)(ceil(minY)-buffer);
        rows[1]                = (int)(ceil(maxY)+buffer);
        
        null_buffer            = 0;
        // Null pixel value remove
        if(cols[0]            <= null_buffer)
            cols[0]            = null_buffer;
        if(rows[0]            <= null_buffer)
            rows[0]            = null_buffer;
        if(cols[0]            > Imagesize->width - null_buffer)
            cols[0]            = Imagesize->width - null_buffer;
        if(rows[0]            > Imagesize->height - null_buffer)
            rows[0]            = Imagesize->height - null_buffer;
        
        if(cols[1]            <= null_buffer)
            cols[1]            = null_buffer;
        if(rows[1]            <= null_buffer)
            rows[1]            = null_buffer;
        if(cols[1]            > Imagesize->width - null_buffer)
            cols[1]            = Imagesize->width - null_buffer;
        if(rows[1]            > Imagesize->height - null_buffer)
            rows[1]            = Imagesize->height - null_buffer;
        
        //printf("cols rows %d\t%d\t%d\t%d\n",cols[0],cols[1],rows[0],rows[1]);
        
        free(ImageCoord);
        
        ret    = true;
    }
    
    return ret;
}

D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres)
{
    D2DPOINT *IP;
    
    IP        = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
    
#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        IP[i].m_X = ( _GP[i].m_X - boundary[0])/imageres;
        IP[i].m_Y = (-_GP[i].m_Y + boundary[3])/imageres;
    }
    
    return IP;
}

void Preprocessing_SDM(ProInfo proinfo, char *save_path, char *Lsubsetfile, char *Rsubsetfile, uint8 py_level, CSize *Lsubsetsize, CSize *Rsubsetsize, CSize *data_size_l, CSize *data_size_r, FILE *fid)
{
    uint16 **pyimg;
    uint16 **magimg;
    int16  **dirimg;
    uint8  **oriimg;
    CSize *data_size;
    FILE *pFile_raw, *pFile_check_file;
    char* filename_py;
    
    uint8 count = 0;
    for(count = 0; count<2 ; count++)
    {
        char t_str[500];
        if(count == 0)
        {
            pFile_raw    = fopen(Lsubsetfile,"rb");
            filename_py        = GetFileName(Lsubsetfile);
            filename_py        = remove_ext(filename_py);
        }
        else
        {
            pFile_raw    = fopen(Rsubsetfile,"rb");
            filename_py        = GetFileName(Rsubsetfile);
            filename_py        = remove_ext(filename_py);
        }
        
        sprintf(t_str,"%s/%s_py_10.raw",save_path,filename_py);
        pFile_check_file    = fopen(t_str,"rb");
        if(pFile_raw && !pFile_check_file)
        {
            int i;
            
            FILE *pFile;
            
            pyimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
            magimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
            dirimg = (int16**)malloc(sizeof(int16*)*(py_level+1));
            oriimg = (uint8**)malloc(sizeof(uint8*)*(py_level+1));
            data_size = (CSize*)malloc(sizeof(CSize)*(py_level+1));
            
            if(count == 0)
            {
                fprintf(fid,"left\nlevel = 0\t width = %d\theight = %d\n",data_size_l[0].width,data_size_l[0].height);
                filename_py        = GetFileName(Lsubsetfile);
                filename_py        = remove_ext(filename_py);
                for(i=0;i<py_level+1;i++)
                    data_size[i]        = data_size_l[i];
                
            }
            else
            {
                fprintf(fid,"right\nlevel = 0\t width = %d\theight = %d\n",data_size_r[0].width,data_size_r[0].height);
                filename_py        = GetFileName(Rsubsetfile);
                filename_py        = remove_ext(filename_py);
                for(i=0;i<py_level+1;i++)
                    data_size[i]        = data_size_r[i];
            }
            
            long int data_length = (long int)data_size[0].height*(long int)data_size[0].width;
            pyimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
            magimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
            dirimg[0] = (int16*)malloc(sizeof(int16)*data_length);
            oriimg[0] = (uint8*)malloc(sizeof(uint8)*data_length);
            
            switch(proinfo.image_bits)
            {
                case 8:
                {
                    uint8* data8 = (uint8*)malloc(sizeof(uint8)*data_length);
                    fread(data8,sizeof(uint8),data_length,pFile_raw);
                    #pragma omp parallel for schedule(guided)
                    for(long int index = 0 ; index < data_length ; index++)
                        pyimg[0][index] = data8[index];
                    free(data8);
                }
                    break;
                case 12:
                {
                    uint16* data16 = (uint16*)malloc(sizeof(uint16)*data_length);
                    fread(data16,sizeof(uint16),data_length,pFile_raw);
#pragma omp parallel for schedule(guided)
                    for(long int index = 0 ; index < data_length ; index++)
                        pyimg[0][index] = data16[index];
                    free(data16);
                }
                    break;
            }
            
            MakeSobelMagnitudeImage(data_size[0],pyimg[0],magimg[0],dirimg[0]);
            Orientation(data_size[0],magimg[0],dirimg[0],15,oriimg[0]);
            
            sprintf(t_str,"%s/%s_py_0.raw",save_path,filename_py);
            pFile    = fopen(t_str,"wb");
            fwrite(pyimg[0],sizeof(uint16),data_length,pFile);
            fclose(pFile);
            
            sprintf(t_str,"%s/%s_py_0_mag.raw",save_path,filename_py);
            pFile    = fopen(t_str,"wb");
            fwrite(magimg[0],sizeof(uint16),data_length,pFile);
            fclose(pFile);
            free(magimg[0]);
            
            free(dirimg[0]);
            
            sprintf(t_str,"%s/%s_py_0_ori.raw",save_path,filename_py);
            pFile    = fopen(t_str,"wb");
            fwrite(oriimg[0],sizeof(uint8),data_length,pFile);
            fclose(pFile);
            free(oriimg[0]);
            
            for(i=0;i<py_level;i++)
            {
                pyimg[i+1] = CreateImagePyramid(pyimg[i],data_size[i],9,(double)(1.5));
                free(pyimg[i]);
                
                if(count == 0)
                {
                    data_size_l[i+1]    = data_size[i+1];
                    fprintf(fid,"level = %d\t width = %d\theight = %d\n",i+1,data_size_l[i+1].width,data_size_l[i+1].height);
                }
                else
                {
                    data_size_r[i+1]    = data_size[i+1];
                    fprintf(fid,"level = %d\t width = %d\theight = %d\n",i+1,data_size_r[i+1].width,data_size_r[i+1].height);
                }
                
                long int data_length_array[6];
                data_length_array[i] = (long int)data_size[i+1].height*(long int)data_size[i+1].width;
                magimg[i+1] = (uint16*)malloc(sizeof(uint16)*data_length_array[i]);
                dirimg[i+1] = (int16*)malloc(sizeof(int16)*data_length_array[i]);
                oriimg[i+1] = (uint8*)malloc(sizeof(uint8)*data_length_array[i]);
                MakeSobelMagnitudeImage(data_size[i+1],pyimg[i+1],magimg[i+1],dirimg[i+1]);
                Orientation(data_size[i+1],magimg[i+1],dirimg[i+1],15,oriimg[i+1]);
                
                sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,i+1);
                pFile    = fopen(t_str,"wb");
                fwrite(pyimg[i+1],sizeof(uint16),data_length_array[i],pFile);
                fclose(pFile);
                if(i == py_level-1)
                    free(pyimg[i+1]);
                
                sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,i+1);
                pFile    = fopen(t_str,"wb");
                fwrite(magimg[i+1],sizeof(uint16),data_length_array[i],pFile);
                fclose(pFile);
                free(magimg[i+1]);
                
                free(dirimg[i+1]);
                
                sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,i+1);
                pFile    = fopen(t_str,"wb");
                fwrite(oriimg[i+1],sizeof(uint8),data_length_array[i],pFile);
                fclose(pFile);
                free(oriimg[i+1]);
            }
            if(pyimg)
                free(pyimg);
            if(magimg)
                free(magimg);
            if(dirimg)
                free(dirimg);
            if(oriimg)
                free(oriimg);
            if(data_size)
                free(data_size);
        }
        
        if(pFile_raw)
            fclose(pFile_raw);
        if(pFile_check_file)
            fclose(pFile_check_file);
    }
}
*/

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

bool VerticalLineLocus_SDM(ProInfo proinfo, NCCresultSDM* nccresult, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, CSize Size_Grid2D, D2DPOINT* GridPts, UGRIDSDM *GridPT3, uint8 Pyramid_step, uint8 start_py, D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration, double* Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2, double* Coreg_param,uint8* SubOriImages_L,uint8* SubOriImages_R)
{
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(im_resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
        
        printf("VerticalLineLocus : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }
    
    int Half_template_size = (int)(Template_size/2);
    
    double subBoundary[4];
    
    D2DPOINT temp_p1, temp_p2;
    D3DPOINT temp_GrP;
    double temp_LIA[2] = {0,0};
    
    int numofpts;
    
    subBoundary[0]      = Boundary[0];
    subBoundary[1]      = Boundary[1];
    subBoundary[2]      = Boundary[2];
    subBoundary[3]      = Boundary[3];
    
    numofpts = Size_Grid2D.height*Size_Grid2D.width;
    
    int count_0_6 = 0;
    int count_0_3 = 0;
    int count_low = 0;
    printf("numofpts %d\t%d\t%d\tcoreg %f\t%f\n",numofpts,Size_Grid2D.height,Size_Grid2D.width,Coreg_param[0],Coreg_param[1]);
#pragma omp parallel for schedule(guided) //reduction(+:count_0_6,count_0_3,count_low)
    for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
    {
        int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        int pts_col = iter_count % Size_Grid2D.width;
        int pt_index = pts_row*Size_Grid2D.width + pts_col;
        
        if(pts_row >= 0 && pts_row < Size_Grid2D.height && pts_col >= 0 && pts_col < Size_Grid2D.width && pt_index >= 0 && pt_index < Size_Grid2D.height*Size_Grid2D.width)
        {
            double max_1stroh = -1.0;
            
            int i,j;
            
            //if(check_image_boundary(proinfo,Lstartpos,Rstartpos,GridPts[pt_index],LImagesize,RImagesize,Half_template_size,Pyramid_step))
            {
                char t_temp_path[500];
                bool check_blunder_cell = false;
                int kernel_size = proinfo.SDM_SS;
                
                if(GridPT3[pt_index].ortho_ncc > 0.6)
                {
                    kernel_size = kernel_size - 1;
                    //count_0_6++;
                }
                else if(GridPT3[pt_index].ortho_ncc > 0.2)
                {
                    //kernel_size = 3;
                    //count_0_3++;
                }
                else
                {
                    kernel_size = kernel_size + 1;
                    //count_low++;
                }
                
                if(GridPT3[pt_index].ortho_ncc < 0.8 /*|| (Pyramid_step == start_py )*/)
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
                            
                            D2DPOINT temp_GP,temp_GP_R;
                            D2DPOINT Left_Imagecoord;
                            D2DPOINT Right_Imagecoord;
                            D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
                            
                            temp_GP.m_X = GridPts[pt_index].m_X;
                            temp_GP.m_Y = GridPts[pt_index].m_Y;
                            
                            temp_GP_R.m_X = GridPts[pt_index].m_X;
                            temp_GP_R.m_Y = GridPts[pt_index].m_Y;
                            
                            Left_Imagecoord        = GetObjectToImage_single(1,temp_GP,proinfo.LBoundary,proinfo.resolution);
                            Right_Imagecoord        = GetObjectToImage_single(1,temp_GP_R,proinfo.RBoundary,proinfo.resolution);
                            
                            Left_Imagecoord_py    = OriginalToPyramid_single(Left_Imagecoord,Lstartpos,Pyramid_step);
                            Right_Imagecoord_py    = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,Pyramid_step);
                            
                            Right_Imagecoord_py.m_Y += (kernel_row + (GridPT3[pt_index].row_shift + Coreg_param[0])/pow(2.0,Pyramid_step));
                            Right_Imagecoord_py.m_X += (kernel_col + (GridPT3[pt_index].col_shift + Coreg_param[1])/pow(2.0,Pyramid_step));
                            
                            long int ori_left = (long int)Left_Imagecoord_py.m_Y * LImagesize.width + (long int)Left_Imagecoord_py.m_X;
                            long int ori_right = (long int)Right_Imagecoord_py.m_Y * RImagesize.width + (long int)Right_Imagecoord_py.m_X;
                            
                            if( Right_Imagecoord_py.m_Y >= 0 && Right_Imagecoord_py.m_Y < RImagesize.height && Right_Imagecoord_py.m_X    >= 0 && Right_Imagecoord_py.m_X    < RImagesize.width &&
                               Left_Imagecoord_py.m_Y >= 0 && Left_Imagecoord_py.m_Y      < LImagesize.height && Left_Imagecoord_py.m_X    >= 0 && Left_Imagecoord_py.m_X    < LImagesize.width)
                            {
                                int ori_diff = SubOriImages_L[ori_left] - SubOriImages_R[ori_right];
                                double bin_angle            = 360.0/18.0;
                                double Left_CR, Left_CC, Right_CR, Right_CC;
                                
                                Left_CR        = Left_Imagecoord_py.m_Y;
                                Left_CC        = Left_Imagecoord_py.m_X;
                                Right_CR    = Right_Imagecoord_py.m_Y;
                                Right_CC    = Right_Imagecoord_py.m_X;
                                
                                //if(check_orientation)
                                {
                                    double rot_theta = 0;//(double)((double)ori_diff*bin_angle*PI/180.0);
                                    
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
                                    
                                    double Sum_LR_mag = 0;
                                    double Sum_L_mag = 0;
                                    double Sum_R_mag = 0;
                                    double Sum_L2_mag = 0;
                                    double Sum_R2_mag = 0;
                                    double Sum_LR_2_mag = 0;
                                    double Sum_L_2_mag = 0;
                                    double Sum_R_2_mag = 0;
                                    double Sum_L2_2_mag = 0;
                                    double Sum_R2_2_mag = 0;
                                    double Sum_LR_3_mag = 0;
                                    double Sum_L_3_mag = 0;
                                    double Sum_R_3_mag = 0;
                                    double Sum_L2_3_mag = 0;
                                    double Sum_R2_3_mag = 0;
                                    
                                    int row, col;
                                    int N;
                                    
                                    double val1, val2, de, de2, ncc_1, ncc_2, ncc_3;
                                    double ncc_1_mag, ncc_2_mag, ncc_3_mag;
                                    
                                    double temp_rho;
                                    double temp_INCC_roh = 0;
                                    bool flag_value;
                                    int grid_index;
                                    double diff_rho;
                                    int t_direction;
                                    
                                    
                                    uint16 mag_center_l, mag_center_r;
                                    
                                    double im_resolution_mask = (gsd_image1.pro_GSD + gsd_image2.pro_GSD)/2.0;
                                    
                                    //printf("im resolution mask %f\n",im_resolution_mask);
                                    
                                    for(row = -Half_template_size; row <= Half_template_size ; row++)
                                    {
                                        for(col = -Half_template_size; col <= Half_template_size ; col++)
                                        {
                                            double row_distance = row*im_resolution_mask*pwrtwo(Pyramid_step);
                                            double col_distance = col*im_resolution_mask*pwrtwo(Pyramid_step);
                                            
                                            double row_pixel_left = row_distance/(gsd_image1.row_GSD*pwrtwo(Pyramid_step));
                                            double col_pixel_left = col_distance/(gsd_image1.col_GSD*pwrtwo(Pyramid_step));
                                            
                                            double row_pixel_right = row_distance/(gsd_image2.row_GSD*pwrtwo(Pyramid_step));
                                            double col_pixel_right = col_distance/(gsd_image2.col_GSD*pwrtwo(Pyramid_step));
                                            
                                            
                                            //printf("row col %d\t%d\tleft %f\t%f\tright %f\t%f\n",row, col, row_pixel_left,col_pixel_left,row_pixel_right,col_pixel_right);
                                            
                                            double radius  = sqrt((double)(row*row + col*col));
                                            if(radius <= Half_template_size+1)
                                            {
                                                //double pos_row_left         = (Left_CR + row);
                                                //double pos_col_left         = (Left_CC + col);
                                                
                                                double pos_row_left         = (Left_CR + row_pixel_left);
                                                double pos_col_left         = (Left_CC + col_pixel_left);
                                                
                                                //double temp_col           = (cos(-rot_theta)*col - sin(-rot_theta)*row);
                                                //double temp_row           = (sin(-rot_theta)*col + cos(-rot_theta)*row);
                                                
                                                double temp_col           = (cos(-rot_theta)*col_pixel_right - sin(-rot_theta)*row_pixel_right);
                                                double temp_row           = (sin(-rot_theta)*col_pixel_right + cos(-rot_theta)*row_pixel_right);
                                                
                                                double pos_row_right     = (Right_CR + temp_row);
                                                double pos_col_right     = (Right_CC + temp_col);
                                                
                                                
                                                if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right    >= 0 && pos_col_right+1    < RImagesize.width &&
                                                   pos_row_left >= 0 && pos_row_left+1      < LImagesize.height && pos_col_left    >= 0 && pos_col_left+1    < LImagesize.width)
                                                {
                                                    //interpolate left_patch
                                                    double dx           =  pos_col_left - (int)(pos_col_left);
                                                    double dy           =  pos_row_left - (int)(pos_row_left);
                                                    double dxdy = dx * dy;
                                                    double left_patch;
                                                    double right_patch;
                                                    double left_mag_patch;
                                                    double right_mag_patch;
                                                    long int position = (long int) pos_col_left + ((long int) pos_row_left) * LImagesize.width;
                                                    
                                                    left_patch = (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
                                                    (double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) +
                                                    (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
                                                    
                                                    left_mag_patch =
                                                    (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_L[position + 1]) * (dx - dxdy) +
                                                    (double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) +
                                                    (double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
                                                    
                                                    
                                                    //interpolate right_patch
                                                    dx            =  pos_col_right - (int)(pos_col_right);
                                                    dy            =  pos_row_right - (int)(pos_row_right);
                                                    dxdy = dx * dy;
                                                    position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
                                                    
                                                    right_patch =
                                                    (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
                                                    (double) (RightImage[position + RImagesize.width]) * (dy - dxdy) +
                                                    (double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
                                                    
                                                    right_mag_patch =
                                                    (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_R[position + 1]) * (dx - dxdy) +
                                                    (double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) +
                                                    (double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
                                                    
                                                    //end
                                                    Count_N[0]++;
                                                    
                                                    //*************
                                                    //Precomputing LR, L2 and R2 etc as they are used in different reductions. (Perhaps compiler handles that by itself)
                                                    double LR = left_patch * right_patch;
                                                    double L2 = left_patch * left_patch;
                                                    double R2 = right_patch * right_patch;
                                                    double LR_mag = left_mag_patch * right_mag_patch;
                                                    double L2_mag = left_mag_patch * left_mag_patch;
                                                    double R2_mag = right_mag_patch * right_mag_patch;
                                                    
                                                    Sum_LR              = Sum_LR + LR;
                                                    Sum_L              = Sum_L  + left_patch;
                                                    Sum_R              = Sum_R  + right_patch;
                                                    Sum_L2              = Sum_L2 + L2;
                                                    Sum_R2              = Sum_R2 + R2;
                                                    
                                                    Sum_LR_mag              = Sum_LR_mag + LR_mag;
                                                    Sum_L_mag              = Sum_L_mag  + left_mag_patch;
                                                    Sum_R_mag              = Sum_R_mag  + right_mag_patch;
                                                    Sum_L2_mag              = Sum_L2_mag + L2_mag;
                                                    Sum_R2_mag              = Sum_R2_mag + R2_mag;
                                                    
                                                    int size_1, size_2;
                                                    size_1          = (int)(Half_template_size/2);
                                                    if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
                                                    {
                                                        if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
                                                        {
                                                            Sum_LR_2  = Sum_LR_2 + LR;
                                                            Sum_L_2      = Sum_L_2     + left_patch;
                                                            Sum_R_2      = Sum_R_2     + right_patch;
                                                            Sum_L2_2  = Sum_L2_2 + L2;
                                                            Sum_R2_2  = Sum_R2_2 + R2;
                                                            Count_N[1]++;
                                                            
                                                            Sum_LR_2_mag  = Sum_LR_2_mag + LR_mag;
                                                            Sum_L_2_mag      = Sum_L_2_mag     + left_mag_patch;
                                                            Sum_R_2_mag      = Sum_R_2_mag     + right_mag_patch;
                                                            Sum_L2_2_mag  = Sum_L2_2_mag + L2_mag;
                                                            Sum_R2_2_mag  = Sum_R2_2_mag + R2_mag;
                                                            
                                                        }
                                                    }
                                                    
                                                    size_2          = size_1 + (int)((size_1/2.0) + 0.5);
                                                    if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
                                                    {
                                                        if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
                                                        {
                                                            Sum_LR_3  = Sum_LR_3 + LR;
                                                            Sum_L_3      = Sum_L_3     + left_patch;
                                                            Sum_R_3      = Sum_R_3     + right_patch;
                                                            Sum_L2_3  = Sum_L2_3 + L2;
                                                            Sum_R2_3  = Sum_R2_3 + R2;
                                                            Count_N[2]++;
                                                            
                                                            Sum_LR_3_mag  = Sum_LR_3_mag + LR_mag;
                                                            Sum_L_3_mag      = Sum_L_3_mag     + left_mag_patch;
                                                            Sum_R_3_mag      = Sum_R_3_mag     + right_mag_patch;
                                                            Sum_L2_3_mag  = Sum_L2_3_mag + L2_mag;
                                                            Sum_R2_3_mag  = Sum_R2_3_mag + R2_mag;
                                                            
                                                        }
                                                    }
                                                    
                                                    if(row == 0 && col == 0)
                                                    {
                                                        mag_center_l = left_patch;
                                                        mag_center_r = right_patch;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    N                = Count_N[0];
                                    val1          = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
                                    val2          = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de              = sqrt(val1*val2);
                                    de2              = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
                                    
                                    if( val1*val2 > 0)
                                        ncc_1            = de2/de;
                                    else
                                        ncc_1            = -1.0;
                                    
                                    val1          = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/N;
                                    val2          = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de              = sqrt(val1*val2);
                                    de2              = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/N;
                                    if( val1*val2 > 0)
                                        ncc_1_mag            = de2/de;
                                    else
                                        ncc_1_mag            = -1.0;
                                    
                                    N                    = Count_N[1];
                                    val1                = (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
                                    val2                = (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de                    = sqrt(val1*val2);
                                    de2                    = (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
                                    if( val1*val2 > 0)
                                        ncc_2          = de2/de;
                                    else
                                        ncc_2            = -1.0;
                                    
                                    val1                = (double)(Sum_L2_2_mag) - (double)(Sum_L_2_mag*Sum_L_2_mag)/N;
                                    val2                = (double)(Sum_R2_2_mag) - (double)(Sum_R_2_mag*Sum_R_2_mag)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de                    = sqrt(val1*val2);
                                    de2                    = (double)(Sum_LR_2_mag) - (double)(Sum_L_2_mag*Sum_R_2_mag)/N;
                                    if( val1*val2 > 0)
                                        ncc_2_mag          = de2/de;
                                    else
                                        ncc_2_mag            = -1.0;
                                    
                                    
                                    N                    = Count_N[2];
                                    val1                = (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
                                    val2                = (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de                    = sqrt(val1*val2);
                                    de2                    = (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
                                    if( val1*val2 > 0)
                                        ncc_3          = de2/de;
                                    else
                                        ncc_3            = -1.0;
                                    
                                    val1                = (double)(Sum_L2_3_mag) - (double)(Sum_L_3_mag*Sum_L_3_mag)/N;
                                    val2                = (double)(Sum_R2_3_mag) - (double)(Sum_R_3_mag*Sum_R_3_mag)/N;
                                    if(Pyramid_step <= 1)
                                    {
                                        if(val1 == 0)
                                            val1 = 0.00001;
                                        if(val2 == 0)
                                            val2 = 0.00001;
                                    }
                                    
                                    de                    = sqrt(val1*val2);
                                    de2                    = (double)(Sum_LR_3_mag) - (double)(Sum_L_3_mag*Sum_R_3_mag)/N;
                                    if( val1*val2 > 0)
                                        ncc_3_mag          = de2/de;
                                    else
                                        ncc_3_mag            = -1.0;
                                    
                                    flag_value        = true;
                                    
                                    temp_INCC_roh = (double)(ncc_1 + ncc_2 + ncc_3 + ncc_1_mag + ncc_2_mag + ncc_3_mag)/6.0;
                                    
                                    temp_rho = temp_INCC_roh;
                                    
                                    grid_index             = pt_index;
                                    
                                    
                                    if(max_1stroh < temp_rho)
                                    {
                                        max_1stroh = temp_rho;
                                        nccresult[grid_index].result0 = (max_1stroh);
                                        
                                        nccresult[grid_index].result2.m_X = Left_Imagecoord_py.m_X;
                                        nccresult[grid_index].result2.m_Y = Left_Imagecoord_py.m_Y;
                                        nccresult[grid_index].result3.m_X = kernel_col + GridPT3[pt_index].col_shift/pow(2.0,Pyramid_step);
                                        nccresult[grid_index].result3.m_Y = kernel_row + GridPT3[pt_index].row_shift/pow(2.0,Pyramid_step);
                                    }
                                    
                                }
                            }
                        }
                    }
                }
                else
                {
                /*
                    D2DPOINT temp_GP,temp_GP_R;
                    D2DPOINT Left_Imagecoord;
                    D2DPOINT Right_Imagecoord;
                    D2DPOINT Left_Imagecoord_py(0), Right_Imagecoord_py(0);
                    
                    temp_GP.m_X = GridPts[pt_index].m_X;
                    temp_GP.m_Y = GridPts[pt_index].m_Y;
                    
                    temp_GP_R.m_X = GridPts[pt_index].m_X;
                    temp_GP_R.m_Y = GridPts[pt_index].m_Y;
                    
                    Left_Imagecoord        = GetObjectToImage_single_SDM(1,temp_GP,proinfo.LBoundary,proinfo.resolution);
                    Right_Imagecoord        = GetObjectToImage_single_SDM(1,temp_GP_R,proinfo.RBoundary,proinfo.resolution);
                    */
                    nccresult[pt_index].result0 = GridPT3[pt_index].ortho_ncc;
                    //nccresult[pt_index].result2.m_X = Left_Imagecoord_py.m_X;
                    //nccresult[pt_index].result2.m_Y = Left_Imagecoord_py.m_Y;
                    nccresult[pt_index].result3.m_X = GridPT3[pt_index].col_shift/pow(2.0,Pyramid_step);
                    nccresult[pt_index].result3.m_Y = GridPT3[pt_index].row_shift/pow(2.0,Pyramid_step);
                }
            }
        }
    }
    
    //printf("count %d\t%d\t%d\n",count_0_6,count_0_3,count_low);
    
    return true;
}

int SelectMPs_SDM(ProInfo proinfo, NCCresultSDM* roh_height, CSize Size_Grid2D, D2DPOINT *GridPts_XY, UGRIDSDM *GridPT3, uint8 Pyramid_step, uint8 prc_level,uint8 iteration, char *filename_mps,uint8 tile_row,uint8 tile_col,double *Boundary)
{
    int count_MPs = 0;
    
    double roh_next;
    int row,col;
    double minimum_Th = 0.3;
    //if(Pyramid_step == 0)
    //    minimum_Th = 0.1;
    
    printf("minimum TH %f\n",minimum_Th);
    FILE* temp_fid;
    double h_divide;
    
    temp_fid            = fopen(filename_mps,"w");
    
    roh_next    = (double)(0.10);
    
    char shift_file[500];
    FILE *fid_col_shift, *fid_row_shift;
    sprintf(shift_file,"%s/txt/col_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,Pyramid_step,iteration,tile_row,tile_col);
    fid_col_shift            = fopen(shift_file,"w");
    sprintf(shift_file,"%s/txt/row_shift_%d_%d_%d_%d.txt",proinfo.save_filepath,Pyramid_step,iteration,tile_row,tile_col);
    fid_row_shift            = fopen(shift_file,"w");
    
    for(int iter_index = 0 ; iter_index < Size_Grid2D.height*Size_Grid2D.width ; iter_index++)
    {
        int row,col;
        row        = (int)(floor(iter_index/Size_Grid2D.width));
        col        = iter_index % Size_Grid2D.width;
        
        if(row >= 0 && row < Size_Grid2D.height && col >= 0 && col < Size_Grid2D.width)
        {
            bool index,index_1,index_2,index_3, roh_index;
            double temp, ROR, temp_h;
            double roh_th;
            D3DPOINT temp_mp;
            
            int grid_index = row*Size_Grid2D.width + col;
            
            GridPT3[grid_index].ortho_ncc = roh_height[grid_index].result0;
            
            roh_index        = false;
            
            roh_th        = roh_next;
            
            if(roh_height[grid_index].result0 > GridPT3[grid_index].roh - roh_th && roh_height[grid_index].result0 > minimum_Th)
                roh_index = true;
            
            {
                //Set the matched pts and information
                if(roh_index && GridPts_XY[grid_index].m_X > Boundary[0] && GridPts_XY[grid_index].m_X < Boundary[2] && GridPts_XY[grid_index].m_Y > Boundary[1] && GridPts_XY[grid_index].m_Y < Boundary[3])
                {
                    count_MPs++;
                    
                    fprintf(temp_fid,"%f %f %f %f %f %f\n",GridPts_XY[grid_index].m_X,GridPts_XY[grid_index].m_Y,roh_height[grid_index].result2.m_X,roh_height[grid_index].result2.m_Y,roh_height[grid_index].result3.m_X,roh_height[grid_index].result3.m_Y);
                    
                    fprintf(fid_col_shift,"%f\t%f\t%f\n",GridPts_XY[grid_index].m_X,GridPts_XY[grid_index].m_Y,roh_height[grid_index].result3.m_X*pow(2.0,prc_level));
                    fprintf(fid_row_shift,"%f\t%f\t%f\n",GridPts_XY[grid_index].m_X,GridPts_XY[grid_index].m_Y,roh_height[grid_index].result3.m_Y*pow(2.0,prc_level));
                    
                    // update max_roh value
                    GridPT3[grid_index].roh        = roh_height[grid_index].result0;
                    if(GridPT3[grid_index].roh < minimum_Th)
                    {
                        GridPT3[grid_index].roh = minimum_Th;
                    }
                    
                    
                }
            }
        }
    }
    fclose(temp_fid);
    
    fclose(fid_col_shift);
    fclose(fid_row_shift);
    
    return count_MPs;
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

void echoprint_Gridinfo_SDM(uint8 prc_level, ProInfo proinfo, double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary,char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3, const char *add_str)
{
    FILE *outfile_h,*outfile_hr, *outfile_min, *outfile_max,     *outfile_roh, *outfile_flag;
    CSize temp_S;
    char t_str[500];
    int k,j;
    
    /*if(strcmp(add_str,"final") == 0)
    {
        sprintf(t_str,"%s/txt/tin_colshift_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_min    = fopen(t_str,"w");
        sprintf(t_str,"%s/txt/tin_rowshift_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_max    = fopen(t_str,"w");
        sprintf(t_str,"%s/txt/tin_roh_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_roh    = fopen(t_str,"w");
    }
    else*/
    {
        if(level == 0 && iteration == 3)
        {
            sprintf(t_str,"%s/txt/tin_roh_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
            outfile_roh    = fopen(t_str,"w");
        }
    }
    
    printf("Size_Grid2D %d\t%d\n",Size_Grid2D->width,Size_Grid2D->height);

    temp_S.height    = Size_Grid2D->height;
    temp_S.width    = Size_Grid2D->width;
    
    for(k=0;k<temp_S.height;k++)
    {
        for(j=0;j<temp_S.width;j++)
        {
            int matlab_index    = k*temp_S.width + j;
            
            double coord_x = boundary[0] + j*proinfo.DEM_resolution;
            double coord_y = boundary[1] + k*proinfo.DEM_resolution;
            
            int pos_lc = (int)((coord_x - LBoundary[0])/(proinfo.resolution*pow(2.0,prc_level)));
            int pos_lr = (int)((LBoundary[3] - coord_y)/(proinfo.resolution*pow(2.0,prc_level)));
            
            int pos_rc = (int)((coord_x - RBoundary[0])/(proinfo.resolution*pow(2.0,prc_level)));
            int pos_rr = (int)((RBoundary[3] - coord_y)/(proinfo.resolution*pow(2.0,prc_level)));
            
            //printf("%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",j,k,pos_lc,pos_lr,pos_rc,pos_rr,proinfo.DEM_resolution,proinfo.resolution);
            long int l_index = pos_lr*LImagesize.width + pos_lc;
            long int r_index = pos_rr*RImagesize.width + pos_rc;
            
            if(pos_lc >= 0 && pos_lc < LImagesize.width && pos_lr >= 0 && pos_lr < LImagesize.height &&
               pos_rc >= 0 && pos_rc < RImagesize.width && pos_rr >= 0 && pos_rr < RImagesize.height)
            {
                int l_value = LeftImage[l_index];
                int r_value = RightImage[r_index];
                
                if(l_value > 0 && r_value > 0)
                {
                    fprintf(outfile_roh,"%f\t",GridPT3[matlab_index].ortho_ncc);
                }
                else
                    fprintf(outfile_roh,"-1.0\t");
            }
            else
                fprintf(outfile_roh,"-1.0\t");
        }
        fprintf(outfile_roh,"\n");
    }
    fclose(outfile_roh);

}

void echoprint_adjustXYZ(uint8 prc_level,double* LBoundary,double* RBoundary, CSize LImagesize, uint16* LeftImage, CSize RImagesize, uint16* RightImage,double* boundary, ProInfo proinfo, char *save_path,int row,int col,int level, int iteration, CSize *Size_Grid2D, UGRIDSDM *GridPT3, const char *add_str, double gridsize, int d_date)
{
    FILE *outfile_XYZ, *outfile_Xshift, *outfile_Yshift, *outfile_Xsigma, *outfile_Ysigma;
    CSize temp_S;
    char t_str[500];
    int k,j;
    
    if(strcmp(add_str,"final") == 0)
    {
        sprintf(t_str,"%s/txt/tin_xshift_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_Xshift    = fopen(t_str,"w");
        sprintf(t_str,"%s/txt/tin_yshift_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_Yshift    = fopen(t_str,"w");
        /*
        sprintf(t_str,"%s/txt/tin_xsigma_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_Xsigma    = fopen(t_str,"w");
        sprintf(t_str,"%s/txt/tin_ysigma_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
        outfile_Ysigma    = fopen(t_str,"w");
         */
    }
    
    UGRIDSDM* temp_GridPT3;
    long int data_length = (long int)(Size_Grid2D->width)*(long int)(Size_Grid2D->height);
    //temp_GridPT3 = (UGRIDSDM*)malloc(sizeof(UGRIDSDM)*data_length);
    //memcpy(temp_GridPT3,GridPT3,sizeof(UGRIDSDM)*data_length);
    
    printf("Size_Grid2D %d\t%d\tres = %f\n",Size_Grid2D->width,Size_Grid2D->height,proinfo.resolution);

    temp_S.height    = Size_Grid2D->height;
    temp_S.width    = Size_Grid2D->width;
    
    for(k=0;k<Size_Grid2D->height;k++)
    {
        for(j=0;j<Size_Grid2D->width;j++)
        {
            int matlab_index    = k*temp_S.width + j;
            double coord_x = boundary[0] + j*proinfo.DEM_resolution;
            double coord_y = boundary[1] + k*proinfo.DEM_resolution;
            
            int pos_lc = (int)((coord_x - LBoundary[0])/(proinfo.resolution*pow(2.0,prc_level)));
            int pos_lr = (int)((LBoundary[3] - coord_y)/(proinfo.resolution*pow(2.0,prc_level)));
            
            int pos_rc = (int)((coord_x - RBoundary[0])/(proinfo.resolution*pow(2.0,prc_level)));
            int pos_rr = (int)((RBoundary[3] - coord_y)/(proinfo.resolution*pow(2.0,prc_level)));
            
            //printf("%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",j,k,pos_lc,pos_lr,pos_rc,pos_rr,proinfo.DEM_resolution,proinfo.resolution);
            long int l_index = pos_lr*LImagesize.width + pos_lc;
            long int r_index = pos_rr*RImagesize.width + pos_rc;
            
            double DX = GridPT3[matlab_index].col_shift*proinfo.resolution/proinfo.SDM_days;
            double DY = -GridPT3[matlab_index].row_shift*proinfo.resolution/proinfo.SDM_days;
            
            //double Xsigma = GridPT3[matlab_index].col_sigma*proinfo.resolution;
            //double Ysigma = GridPT3[matlab_index].row_sigma*proinfo.resolution;
            
            if(d_date == 2)
            {
                DX = -DX;
                DY = -DY;
            }
            
            if(pos_lc >= 0 && pos_lc < LImagesize.width && pos_lr >= 0 && pos_lr < LImagesize.height &&
               pos_rc >= 0 && pos_rc < RImagesize.width && pos_rr >= 0 && pos_rr < RImagesize.height)
            {
                int l_value = LeftImage[l_index];
                int r_value = RightImage[r_index];
                
                if(l_value > 0 && r_value > 0)
                {
                    fprintf(outfile_Xshift,"%f\t", DX);
                    fprintf(outfile_Yshift,"%f\t", DY);
                    
                    //fprintf(outfile_Xsigma,"%f\t", Xsigma);
                    //fprintf(outfile_Ysigma,"%f\t", Ysigma);
                }
                else
                {
                    fprintf(outfile_Xshift,"0\t");
                    fprintf(outfile_Yshift,"0\t");
                    
                    //fprintf(outfile_Xsigma,"0\t");
                    //fprintf(outfile_Ysigma,"0\t");
                }
            }
            else
            {
                fprintf(outfile_Xshift,"0\t");
                fprintf(outfile_Yshift,"0\t");
                
                //fprintf(outfile_Xsigma,"0\t");
                //fprintf(outfile_Ysigma,"0\t");
            }
            
        }
        //if(strcmp(add_str,"final") == 0)
        {
            fprintf(outfile_Xshift,"\n");
            fprintf(outfile_Yshift,"\n");
            
            //fprintf(outfile_Xsigma,"\n");
            //fprintf(outfile_Ysigma,"\n");
        }
    }
    
    //if(strcmp(add_str,"final") == 0)
    {
        fclose(outfile_Xshift);
        fclose(outfile_Yshift);
        
        //fclose(outfile_Xsigma);
        //fclose(outfile_Ysigma);
    }
    
    //free(temp_GridPT3);
}

UGRIDSDM* SetHeightRange_SDM(UGRIDSDM *GridPT3, BL BL_param)
{
    UGRIDSDM *result;
    
    int i, tcnt;
    double gridspace;
    CSize gridsize;
    uint8 pyramid_step, iteration;
    uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
    
    pyramid_step    = BL_param.Pyramid_step;
    iteration        = BL_param.iteration;
    gridspace = BL_param.gridspace;
    const double *boundary    = BL_param.Boundary;
    gridsize.width    = BL_param.Size_Grid2D.width;
    gridsize.height    = BL_param.Size_Grid2D.height;
    TIN_Grid_Size_X= gridsize.width;
    TIN_Grid_Size_Y= gridsize.height;
    
    int k, j;
    
    result                    = (UGRIDSDM*)calloc(TIN_Grid_Size_Y*TIN_Grid_Size_X,sizeof(UGRIDSDM));
    
    //#pragma omp parallel for shared(TIN_Grid_Size_X,TIN_Grid_Size_Y,GridPT3,result,m_bHeight,Total_Min_Z,Total_Max_Z) private(k,j)
    for(k=0;k<(int)(TIN_Grid_Size_Y);k++)
    {
        for(j=0;j<(int)(TIN_Grid_Size_X);j++)
        {
            int matlab_index    = k*TIN_Grid_Size_X + j;
            
            result[matlab_index].col_shift                  = GridPT3[matlab_index].col_shift;
            result[matlab_index].row_shift                  = GridPT3[matlab_index].row_shift;
            
            result[matlab_index].ortho_ncc                  = GridPT3[matlab_index].ortho_ncc;
            result[matlab_index].roh                        = GridPT3[matlab_index].roh;
            
            //result[matlab_index].col_sigma                  = GridPT3[matlab_index].col_sigma;
            //result[matlab_index].row_sigma                  = GridPT3[matlab_index].row_sigma;
        }
    }
    
    
    printf("end updating grid set height!!\n");
    
    free(GridPT3);
    
    printf("end memory release updating grid set height!!\n");
    
    return result;
}

void SetHeightRange_shift(int numOfPts, int numOfTri, UGRIDSDM *GridPT3,BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int b_dir)
{
    uint32 num_triangles;
    int i, tcnt;
    double gridspace;
    CSize gridsize;
    uint8 pyramid_step, iteration;
    uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
    uint8* m_bHeight;
    
    double Total_Min_Z        =  100000;
    double Total_Max_Z        = -100000;
    
    int Col_C, Row_R;
    
    num_triangles            = numOfTri;
    pyramid_step    = BL_param.Pyramid_step;
    iteration        = BL_param.iteration;
    gridspace = BL_param.gridspace;
    const double *boundary    = BL_param.Boundary;
    gridsize.width    = BL_param.Size_Grid2D.width;
    gridsize.height    = BL_param.Size_Grid2D.height;
    TIN_Grid_Size_X= gridsize.width;
    TIN_Grid_Size_Y= gridsize.height;
    
    
    m_bHeight        = (uint8*)calloc(TIN_Grid_Size_Y*TIN_Grid_Size_X,sizeof(uint8));
    
    //    #pragma omp parallel shared(GridPT3,pts,tris,num_triangles,m_bHeight,pyramid_step,iteration,gridspace,boundary,gridsize,TIN_Grid_Size_X,TIN_Grid_Size_Y,DEM_error) private(tcnt)
    {
        //        #pragma omp for
        for(tcnt=0;tcnt<(int)(num_triangles);tcnt++)
        {
            UI3DPOINT t_tri;
            D3DPOINT pt0,pt1,pt2;
            int pdex0,pdex1,pdex2;
            double TriP1[3];
            double TriP2[3];
            double TriP3[3];
            
            
            double temp_MinZ, temp_MaxZ;
            double TriMinXY[2], TriMaxXY[2];
            int PixelMinXY[2]={0};
            int PixelMaxXY[2]={0};
            int Col, Row;
            int numPts = numOfPts;
            
            t_tri    = tris[tcnt];
            pdex0 = t_tri.m_X;
            pdex1 = t_tri.m_Y;
            pdex2 = t_tri.m_Z;
            
            if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
            {
                int node1_index, node2_index, node3_index;
                bool check_anchor_flag = false;
                pt0    = pts[pdex0];
                pt1    = pts[pdex1];
                pt2    = pts[pdex2];
                
                TriP1[0]        = pt0.m_X;
                TriP2[0]        = pt1.m_X;
                TriP3[0]        = pt2.m_X;
                
                TriP1[1]        = pt0.m_Y;
                TriP2[1]        = pt1.m_Y;
                TriP3[1]        = pt2.m_Y;
                
                TriP1[2]        = pt0.m_Z;
                TriP2[2]        = pt1.m_Z;
                TriP3[2]        = pt2.m_Z;
                
                temp_MinZ = min(min(TriP1[2],TriP2[2]),TriP3[2]);
                temp_MaxZ = max(max(TriP1[2],TriP2[2]),TriP3[2]);
                
                if(temp_MinZ < Total_Min_Z)
                    Total_Min_Z = temp_MinZ;
                if(temp_MaxZ > Total_Max_Z)
                    Total_Max_Z = temp_MaxZ;
                
                // calculation on BoundingBox(MinMax XY) of triangle
                TriMinXY[0]    = min(min(TriP1[0],TriP2[0]),TriP3[0]);
                TriMinXY[1]    = min(min(TriP1[1],TriP2[1]),TriP3[1]);
                TriMaxXY[0]    = max(max(TriP1[0],TriP2[0]),TriP3[0]);
                TriMaxXY[1]    = max(max(TriP1[1],TriP2[1]),TriP3[1]);
                
                PixelMinXY[0] = (int)((TriMinXY[0] - boundary[0])/gridspace + 0.5);
                PixelMinXY[1] = (int)((TriMinXY[1] - boundary[1])/gridspace + 0.5);
                PixelMaxXY[0] = (int)((TriMaxXY[0] - boundary[0])/gridspace + 0.5);
                PixelMaxXY[1] = (int)((TriMaxXY[1] - boundary[1])/gridspace + 0.5);
                
                PixelMinXY[0] -= 1;        PixelMinXY[1] -= 1;
                PixelMaxXY[0] += 1;        PixelMaxXY[1] += 1;
                if (PixelMaxXY[0] >= (int)(TIN_Grid_Size_X))
                    PixelMaxXY[0] =     (int)(TIN_Grid_Size_X-1);
                if (PixelMaxXY[1] >= (int)(TIN_Grid_Size_Y))
                    PixelMaxXY[1] =     (int)(TIN_Grid_Size_Y-1);
                if (PixelMinXY[0] < 0)
                    PixelMinXY[0] = 0;
                if (PixelMinXY[1] < 0)
                    PixelMinXY[1] = 0;
                
                for (Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
                {
                    for (Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
                    {
                        int Index= TIN_Grid_Size_X*Row + Col;
                        
                        if(!m_bHeight[Index])
                        {
                            double CurGPXY[2]={0.};
                            double Z = -1000.0;
                            bool rtn = false;
                            double _p1[3], _p2[3], _p3[3], v12[2], v1P[2];
                            double v23[2], v2P[2], v31[2], v3P[2];
                            int Sum;
                            int t_col_count = 0;
                            int t_row_count = 0;
                            
                            CurGPXY[0]    = (Col)*gridspace + boundary[0];
                            CurGPXY[1]    = (Row)*gridspace + boundary[1];
                            
                            _p1[0]        = TriP1[0];
                            _p1[1]        = TriP1[1];
                            _p1[2]        = TriP1[2];
                            
                            _p2[0]        = TriP2[0];
                            _p2[1]        = TriP2[1];
                            _p2[2]        = TriP2[2];
                            
                            _p3[0]        = TriP3[0];
                            _p3[1]        = TriP3[1];
                            _p3[2]        = TriP3[2];
                            
                            v12[0]        = _p2[0]-_p1[0];
                            v12[1]        = _p2[1]-_p1[1];
                            
                            v1P[0]        = CurGPXY[0]-_p1[0];
                            v1P[1]        = CurGPXY[1]-_p1[1];
                            
                            v23[0]        = _p3[0]-_p2[0];
                            v23[1]        = _p3[1]-_p2[1];
                            
                            v2P[0]        = CurGPXY[0]-_p2[0];
                            v2P[1]        = CurGPXY[1]-_p2[1];
                            
                            v31[0]        = _p1[0]-_p3[0];
                            v31[1]        = _p1[1]-_p3[1];
                            
                            v3P[0]        = CurGPXY[0]-_p3[0];
                            v3P[1]        = CurGPXY[1]-_p3[1];
                            
                            Sum = 3;
                            if (v12[0]*v1P[1]-v12[1]*v1P[0] <= 0)
                                Sum--;
                            if (v23[0]*v2P[1]-v23[1]*v2P[0] <= 0)
                                Sum--;
                            if (v31[0]*v3P[1]-v31[1]*v3P[0] <= 0)
                                Sum--;
                            
                            if (Sum==0 || Sum==3)
                            {
                                double v12[3], v13[3], Len, A, B, C, D;
                                double Normal[3]={0.};
                                
                                v12[0]    = _p2[0]-_p1[0];
                                v12[1]    = _p2[1]-_p1[1];
                                v12[2]    = _p2[2]-_p1[2];
                                
                                v13[0]    = _p3[0]-_p1[0];
                                v13[1]    = _p3[1]-_p1[1];
                                v13[2]    = _p3[2]-_p1[2];
                                
                                Normal[0]=v12[1]*v13[2] - v12[2]*v13[1];
                                Normal[1]=v12[2]*v13[0] - v12[0]*v13[2];
                                Normal[2]=v12[0]*v13[1] - v12[1]*v13[0];
                                
                                Len=sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
                                if(Len > 0 )
                                {
                                    Normal[0]/=Len;
                                    Normal[1]/=Len;
                                    Normal[2]/=Len;
                                    
                                    A = Normal[0];
                                    B = Normal[1];
                                    C = Normal[2];
                                    D = -(A*_p1[0]+B*_p1[1]+C*_p1[2]);
                                    
                                    if(C != 0)
                                    {
                                        Z = -1.0 * ((A * CurGPXY[0]) + (B * CurGPXY[1]) + D) / C;
                                        rtn = true;
                                    }
                                    else
                                        rtn = false;
                                }
                            }
                            
                            if (rtn)
                            {
                                //if(GridPT3[Index].pre_ncc < 1)
                                {
                                    double diff1, diff2, diff3, t1, t2;
                                    double pre_minHeight, pre_maxHeight;
                                    
                                    //GridPT3[Index].t_minHeight    = GridPT3[Index].minHeight;
                                    //GridPT3[Index].t_maxHeight    = GridPT3[Index].maxHeight;
                                    
                                    CurGPXY[0]    = (Col)*gridspace + boundary[0];
                                    CurGPXY[1]    = (Row)*gridspace + boundary[1];
                                    
                                    //// IDW
                                    diff1 = sqrt((CurGPXY[0] - TriP1[0])*(CurGPXY[0] - TriP1[0]) + (CurGPXY[1] - TriP1[1])*(CurGPXY[1] - TriP1[1]));
                                    diff2 = sqrt((CurGPXY[0] - TriP2[0])*(CurGPXY[0] - TriP2[0]) + (CurGPXY[1] - TriP2[1])*(CurGPXY[1] - TriP2[1]));
                                    diff3 = sqrt((CurGPXY[0] - TriP3[0])*(CurGPXY[0] - TriP3[0]) + (CurGPXY[1] - TriP3[1])*(CurGPXY[1] - TriP3[1]));
                                    
                                    if(diff1 == 0)
                                    {
                                        Z    = TriP1[2];
                                    }
                                    else if(diff2 == 0)
                                    {
                                        Z    = TriP2[2];
                                    }
                                    else if(diff3 == 0)
                                    {
                                        Z    = TriP3[2];
                                    }
                                    /*else
                                     {
                                     double p_z = 1.5;
                                     double sum_z   = TriP1[2]/pow(diff1,p_z) + TriP2[2]/pow(diff2,p_z) + TriP3[2]/pow(diff3,p_z);
                                     double sum_w   = 1.0/pow(diff1,p_z) + 1.0/pow(diff2,p_z) + 1.0/pow(diff3,p_z);
                                     Z = sum_z / sum_w;
                                     }
                                     */
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
        }
    }
    
    if(m_bHeight)
        free(m_bHeight);
}

bool Update_ortho_NCC(ProInfo proinfo, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, CSize Size_Grid2D, D2DPOINT* GridPts, UGRIDSDM *GridPT3, uint8 Pyramid_step, D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration, double* Boundary, ImageGSD gsd_image1, ImageGSD gsd_image2, double* Coreg_param,uint8* SubOriImages_L,uint8* SubOriImages_R)
{
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(im_resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
        
        //printf("Update ortho NCC : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }
    
    int Half_template_size = (int)(Template_size/2);
    
    double subBoundary[4];
    
    D2DPOINT temp_p1, temp_p2;
    D3DPOINT temp_GrP;
    double temp_LIA[2] = {0,0};
    
    int numofpts;
    
    subBoundary[0]      = Boundary[0];
    subBoundary[1]      = Boundary[1];
    subBoundary[2]      = Boundary[2];
    subBoundary[3]      = Boundary[3];
    
    numofpts = Size_Grid2D.height*Size_Grid2D.width;
    
    int count_0_6 = 0;
    int count_0_3 = 0;
    int count_low = 0;
    //printf("Update ortho NCC numofpts %d\t%d\t%d\tcoreg %f\t%f\n",numofpts,Size_Grid2D.height,Size_Grid2D.width,Coreg_param[0],Coreg_param[1]);
#pragma omp parallel for schedule(guided) //reduction(+:count_0_6,count_0_3,count_low)
    for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
    {
        int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        int pts_col = iter_count % Size_Grid2D.width;
        int pt_index = pts_row*Size_Grid2D.width + pts_col;
        
        if(pts_row >= 0 && pts_row < Size_Grid2D.height && pts_col >= 0 && pts_col < Size_Grid2D.width && pt_index >= 0 && pt_index < Size_Grid2D.height*Size_Grid2D.width)
        {
            double max_1stroh = -1.0;
            
            int i,j;
            
            char t_temp_path[500];
            bool check_blunder_cell = false;
            
            bool check_false_h = false;
            D2DPOINT temp_GP,temp_GP_R;
            D2DPOINT Left_Imagecoord;
            D2DPOINT Right_Imagecoord;
            D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
            
            temp_GP.m_X = GridPts[pt_index].m_X;
            temp_GP.m_Y = GridPts[pt_index].m_Y;
            
            temp_GP_R.m_X = GridPts[pt_index].m_X;
            temp_GP_R.m_Y = GridPts[pt_index].m_Y;
            
            Left_Imagecoord        = GetObjectToImage_single(1,temp_GP,proinfo.LBoundary,proinfo.resolution);
            Right_Imagecoord        = GetObjectToImage_single(1,temp_GP_R,proinfo.RBoundary,proinfo.resolution);
            
            Left_Imagecoord_py    = OriginalToPyramid_single(Left_Imagecoord,Lstartpos,Pyramid_step);
            Right_Imagecoord_py    = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,Pyramid_step);
            
            Right_Imagecoord_py.m_Y += (GridPT3[pt_index].row_shift + Coreg_param[0])/pow(2.0,Pyramid_step);
            Right_Imagecoord_py.m_X += (GridPT3[pt_index].col_shift + Coreg_param[1])/pow(2.0,Pyramid_step);
            
            long int ori_left = (long int)Left_Imagecoord_py.m_Y * LImagesize.width + (long int)Left_Imagecoord_py.m_X;
            long int ori_right = (long int)Right_Imagecoord_py.m_Y * RImagesize.width + (long int)Right_Imagecoord_py.m_X;
            
            if( Right_Imagecoord_py.m_Y >= 0 && Right_Imagecoord_py.m_Y < RImagesize.height && Right_Imagecoord_py.m_X    >= 0 && Right_Imagecoord_py.m_X    < RImagesize.width &&
               Left_Imagecoord_py.m_Y >= 0 && Left_Imagecoord_py.m_Y      < LImagesize.height && Left_Imagecoord_py.m_X    >= 0 && Left_Imagecoord_py.m_X    < LImagesize.width)
            {
                int ori_diff = SubOriImages_L[ori_left] - SubOriImages_R[ori_right];
                double bin_angle            = 360.0/18.0;
                double Left_CR, Left_CC, Right_CR, Right_CC;
                
                Left_CR        = Left_Imagecoord_py.m_Y;
                Left_CC        = Left_Imagecoord_py.m_X;
                Right_CR    = Right_Imagecoord_py.m_Y;
                Right_CC    = Right_Imagecoord_py.m_X;
                
                double rot_theta = 0;//(double)((double)ori_diff*bin_angle*PI/180.0);
                
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
                
                double Sum_LR_mag = 0;
                double Sum_L_mag = 0;
                double Sum_R_mag = 0;
                double Sum_L2_mag = 0;
                double Sum_R2_mag = 0;
                double Sum_LR_2_mag = 0;
                double Sum_L_2_mag = 0;
                double Sum_R_2_mag = 0;
                double Sum_L2_2_mag = 0;
                double Sum_R2_2_mag = 0;
                double Sum_LR_3_mag = 0;
                double Sum_L_3_mag = 0;
                double Sum_R_3_mag = 0;
                double Sum_L2_3_mag = 0;
                double Sum_R2_3_mag = 0;
                
                int row, col;
                int N;
                
                double val1, val2, de, de2, ncc_1, ncc_2, ncc_3;
                double ncc_1_mag, ncc_2_mag, ncc_3_mag;
                
                double temp_rho;
                double temp_INCC_roh = 0;
                bool flag_value;
                int grid_index;
                double diff_rho;
                int t_direction;
                
                
                uint16 mag_center_l, mag_center_r;
                
                double im_resolution_mask = (gsd_image1.pro_GSD + gsd_image2.pro_GSD)/2.0;
                
                //printf("im resolution mask %f\n",im_resolution_mask);
                
                for(row = -Half_template_size; row <= Half_template_size ; row++)
                {
                    for(col = -Half_template_size; col <= Half_template_size ; col++)
                    {
                        double row_distance = row*im_resolution_mask*pwrtwo(Pyramid_step);
                        double col_distance = col*im_resolution_mask*pwrtwo(Pyramid_step);
                        
                        double row_pixel_left = row_distance/(gsd_image1.row_GSD*pwrtwo(Pyramid_step));
                        double col_pixel_left = col_distance/(gsd_image1.col_GSD*pwrtwo(Pyramid_step));
                        
                        double row_pixel_right = row_distance/(gsd_image2.row_GSD*pwrtwo(Pyramid_step));
                        double col_pixel_right = col_distance/(gsd_image2.col_GSD*pwrtwo(Pyramid_step));
                        
                        
                        //printf("row col %d\t%d\tleft %f\t%f\tright %f\t%f\n",row, col, row_pixel_left,col_pixel_left,row_pixel_right,col_pixel_right);
                        
                        double radius  = sqrt((double)(row*row + col*col));
                        if(radius <= Half_template_size+1)
                        {
                            //double pos_row_left         = (Left_CR + row);
                            //double pos_col_left         = (Left_CC + col);
                            
                            double pos_row_left         = (Left_CR + row_pixel_left);
                            double pos_col_left         = (Left_CC + col_pixel_left);
                            
                            //double temp_col           = (cos(-rot_theta)*col - sin(-rot_theta)*row);
                            //double temp_row           = (sin(-rot_theta)*col + cos(-rot_theta)*row);
                            
                            double temp_col           = (cos(-rot_theta)*col_pixel_right - sin(-rot_theta)*row_pixel_right);
                            double temp_row           = (sin(-rot_theta)*col_pixel_right + cos(-rot_theta)*row_pixel_right);
                            
                            double pos_row_right     = (Right_CR + temp_row);
                            double pos_col_right     = (Right_CC + temp_col);
                            
                            
                            if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right    >= 0 && pos_col_right+1    < RImagesize.width &&
                               pos_row_left >= 0 && pos_row_left+1      < LImagesize.height && pos_col_left    >= 0 && pos_col_left+1    < LImagesize.width)
                            {
                                //interpolate left_patch
                                double dx           =  pos_col_left - (int)(pos_col_left);
                                double dy           =  pos_row_left - (int)(pos_row_left);
                                double dxdy = dx * dy;
                                double left_patch;
                                double right_patch;
                                double left_mag_patch;
                                double right_mag_patch;
                                long int position = (long int) pos_col_left + ((long int) pos_row_left) * LImagesize.width;
                                
                                left_patch = (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
                                (double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) +
                                (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
                                
                                left_mag_patch =
                                (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_L[position + 1]) * (dx - dxdy) +
                                (double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) +
                                (double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
                                
                                
                                //interpolate right_patch
                                dx            =  pos_col_right - (int)(pos_col_right);
                                dy            =  pos_row_right - (int)(pos_row_right);
                                dxdy = dx * dy;
                                position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
                                
                                right_patch =
                                (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
                                (double) (RightImage[position + RImagesize.width]) * (dy - dxdy) +
                                (double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
                                
                                right_mag_patch =
                                (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_R[position + 1]) * (dx - dxdy) +
                                (double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) +
                                (double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
                                
                                //end
                                Count_N[0]++;
                                
                                //*************
                                //Precomputing LR, L2 and R2 etc as they are used in different reductions. (Perhaps compiler handles that by itself)
                                double LR = left_patch * right_patch;
                                double L2 = left_patch * left_patch;
                                double R2 = right_patch * right_patch;
                                double LR_mag = left_mag_patch * right_mag_patch;
                                double L2_mag = left_mag_patch * left_mag_patch;
                                double R2_mag = right_mag_patch * right_mag_patch;
                                
                                Sum_LR              = Sum_LR + LR;
                                Sum_L              = Sum_L  + left_patch;
                                Sum_R              = Sum_R  + right_patch;
                                Sum_L2              = Sum_L2 + L2;
                                Sum_R2              = Sum_R2 + R2;
                                
                                Sum_LR_mag              = Sum_LR_mag + LR_mag;
                                Sum_L_mag              = Sum_L_mag  + left_mag_patch;
                                Sum_R_mag              = Sum_R_mag  + right_mag_patch;
                                Sum_L2_mag              = Sum_L2_mag + L2_mag;
                                Sum_R2_mag              = Sum_R2_mag + R2_mag;
                                
                                int size_1, size_2;
                                size_1          = (int)(Half_template_size/2);
                                if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
                                {
                                    if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
                                    {
                                        Sum_LR_2  = Sum_LR_2 + LR;
                                        Sum_L_2      = Sum_L_2     + left_patch;
                                        Sum_R_2      = Sum_R_2     + right_patch;
                                        Sum_L2_2  = Sum_L2_2 + L2;
                                        Sum_R2_2  = Sum_R2_2 + R2;
                                        Count_N[1]++;
                                        
                                        Sum_LR_2_mag  = Sum_LR_2_mag + LR_mag;
                                        Sum_L_2_mag      = Sum_L_2_mag     + left_mag_patch;
                                        Sum_R_2_mag      = Sum_R_2_mag     + right_mag_patch;
                                        Sum_L2_2_mag  = Sum_L2_2_mag + L2_mag;
                                        Sum_R2_2_mag  = Sum_R2_2_mag + R2_mag;
                                        
                                    }
                                }
                                
                                size_2          = size_1 + (int)((size_1/2.0) + 0.5);
                                if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
                                {
                                    if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
                                    {
                                        Sum_LR_3  = Sum_LR_3 + LR;
                                        Sum_L_3      = Sum_L_3     + left_patch;
                                        Sum_R_3      = Sum_R_3     + right_patch;
                                        Sum_L2_3  = Sum_L2_3 + L2;
                                        Sum_R2_3  = Sum_R2_3 + R2;
                                        Count_N[2]++;
                                        
                                        Sum_LR_3_mag  = Sum_LR_3_mag + LR_mag;
                                        Sum_L_3_mag      = Sum_L_3_mag     + left_mag_patch;
                                        Sum_R_3_mag      = Sum_R_3_mag     + right_mag_patch;
                                        Sum_L2_3_mag  = Sum_L2_3_mag + L2_mag;
                                        Sum_R2_3_mag  = Sum_R2_3_mag + R2_mag;
                                        
                                    }
                                }
                                
                                if(row == 0 && col == 0)
                                {
                                    mag_center_l = left_patch;
                                    mag_center_r = right_patch;
                                }
                            }
                        }
                    }
                }
                
                N                = Count_N[0];
                val1          = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
                val2          = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de              = sqrt(val1*val2);
                de2              = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
                
                if( val1*val2 > 0)
                    ncc_1            = de2/de;
                else
                    ncc_1            = -1.0;
                
                val1          = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/N;
                val2          = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de              = sqrt(val1*val2);
                de2              = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/N;
                if( val1*val2 > 0)
                    ncc_1_mag            = de2/de;
                else
                    ncc_1_mag            = -1.0;
                
                N                    = Count_N[1];
                val1                = (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
                val2                = (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de                    = sqrt(val1*val2);
                de2                    = (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
                if( val1*val2 > 0)
                    ncc_2          = de2/de;
                else
                    ncc_2            = -1.0;
                
                val1                = (double)(Sum_L2_2_mag) - (double)(Sum_L_2_mag*Sum_L_2_mag)/N;
                val2                = (double)(Sum_R2_2_mag) - (double)(Sum_R_2_mag*Sum_R_2_mag)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de                    = sqrt(val1*val2);
                de2                    = (double)(Sum_LR_2_mag) - (double)(Sum_L_2_mag*Sum_R_2_mag)/N;
                if( val1*val2 > 0)
                    ncc_2_mag          = de2/de;
                else
                    ncc_2_mag            = -1.0;
                
                
                N                    = Count_N[2];
                val1                = (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
                val2                = (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de                    = sqrt(val1*val2);
                de2                    = (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
                if( val1*val2 > 0)
                    ncc_3          = de2/de;
                else
                    ncc_3            = -1.0;
                
                val1                = (double)(Sum_L2_3_mag) - (double)(Sum_L_3_mag*Sum_L_3_mag)/N;
                val2                = (double)(Sum_R2_3_mag) - (double)(Sum_R_3_mag*Sum_R_3_mag)/N;
                if(Pyramid_step <= 1)
                {
                    if(val1 == 0)
                        val1 = 0.00001;
                    if(val2 == 0)
                        val2 = 0.00001;
                }
                
                de                    = sqrt(val1*val2);
                de2                    = (double)(Sum_LR_3_mag) - (double)(Sum_L_3_mag*Sum_R_3_mag)/N;
                if( val1*val2 > 0)
                    ncc_3_mag          = de2/de;
                else
                    ncc_3_mag            = -1.0;
                
                temp_INCC_roh = (double)(ncc_1 + ncc_2 + ncc_3 + ncc_1_mag + ncc_2_mag + ncc_3_mag)/6.0;
                
                GridPT3[pt_index].ortho_ncc = temp_INCC_roh;
            }
        }
    }
    
    //printf("count %d\t%d\t%d\n",count_0_6,count_0_3,count_low);
    
    return true;
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

void shift_filtering(ProInfo proinfo, UGRIDSDM *GridPT3, BL BL_param, double DEM_resolution)
{
    CSize gridsize;
    uint8 pyramid_step, iteration;
    pyramid_step    = BL_param.Pyramid_step;
    iteration        = BL_param.iteration;
    gridsize.width    = BL_param.Size_Grid2D.width;
    gridsize.height    = BL_param.Size_Grid2D.height;
    
    
    int interpolated_pts = 0;
    int non_iterpolated_pts = 0;
    int angle_pts = 0;
    long int data_length = (long int)(gridsize.width)*(long int)(gridsize.height);
    float *temp_col_shift = (float*)malloc(sizeof(float)*data_length);
    float *temp_row_shift = (float*)malloc(sizeof(float)*data_length);
    
    for(int r = 0 ; r < gridsize.height ; r ++)
    {
        for(int c = 0 ; c < gridsize.width; c++)
        {
            int ori_index = r*gridsize.width + c;
            
            temp_col_shift[ori_index] = GridPT3[ori_index].col_shift;
            temp_row_shift[ori_index] = GridPT3[ori_index].row_shift;
            
        }
    }
    
    int max_kernel_size = 40;
    double max_min_ortho_th = 0.5;
    /*if(pyramid_step >= 2)
        min_ortho_th = 0.0;
    else if(pyramid_step == 1)
        min_ortho_th = 0.05;
    
    max_kernel_size = 3;
    if(pyramid_step >= 5)
        max_kernel_size = 11;
    else if(pyramid_step == 4)
        max_kernel_size = 9;
    else if(pyramid_step == 3)
        max_kernel_size = 7;
    else if(pyramid_step == 2)
        max_kernel_size = 5;
    */
    //    if(b_dir/* && pyramid_step >= 5*/)
    {
        /*
        int kernal_size = 1;
        if(DEM_resolution < 2)
            kernal_size = 1;
        if(pyramid_step >= 5)
            kernal_size = 5;
        else if(pyramid_step == 4)
            kernal_size = 4;
        else if(pyramid_step == 3)
            kernal_size = 3;
        else if(pyramid_step == 2)
            kernal_size = 2;
        
        if(pyramid_step > 4)
         {
         if(pyramid_step <= 6)
         kernal_size = 7;
         else if(pyramid_step == 7)
         kernal_size = 7;
         else
         kernal_size = 7;
         }
         */
        //int total_kernal_size = kernal_size*2 + 1;
        int da = 45;
        if(pyramid_step == 3)
            da = 30;
        else if(pyramid_step == 2)
            da = 20;
        
        int slope_step = 360/da;
        
        printf("angle step %d\t%d\t%d\t%d\n",slope_step,gridsize.width,gridsize.height,data_length);
        
        //double* mean_colshift = (double*)calloc(sizeof(double), gridsize.height*gridsize.width);
        //double* mean_rowshift = (double*)calloc(sizeof(double), gridsize.height*gridsize.width);
        
        
        char t_str[500];
        double average_col = 0;
        double average_row = 0;
        double all_sum_col = 0;
        double all_sum_row = 0;
        double all_count = 0;
        double all_mean_slope = 0;
        double all_average_slope = 0;
        int mean_slope_count = 0;
        
        
        /*
         if(pyramid_step >=4 && iteration < 3)
         {
         for(int r = 0 ; r < gridsize.height ; r ++)
         {
         for(int c = 0 ; c < gridsize.width; c++)
         {
         int ori_index = r*gridsize.width + c;
         double col_save[200] = {0.};
         double row_save[200] = {0.};
         int save_count = 0;
         for(int k = -kernal_size ; k <= kernal_size ; k++)
         {
         for(int j=-kernal_size ; j <= kernal_size ; j++)
         {
         int index = (r+k)*gridsize.width + (c+j);
         if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width && index >= 0 && index < gridsize.height*gridsize.width)
         {
         double col_value = GridPT3[index].col_shift;
         double row_value = GridPT3[index].row_shift;
         
         
         {
         col_save[save_count] = col_value;
         row_save[save_count] = row_value;
         
         save_count++;
         }
         }
         }
         }
         //if(save_count > 0)
         //    printf("save count %d\n",save_count);
         if(save_count > 2)
         {
         bubble_sort(col_save,save_count);
         mean_colshift[ori_index] = col_save[(int)(save_count/2.0)];
         bubble_sort(row_save,save_count);
         mean_rowshift[ori_index] = row_save[(int)(save_count/2.0)];
         //printf("mean %f\t%f\n",mean_colshift[ori_index],mean_rowshift[ori_index]);
         }
         //fprintf(fid_col,"%f\t",GridPT3[ori_index].col_shift);
         //fprintf(fid_row,"%f\t",GridPT3[ori_index].row_shift);
         
         
         if(m_bHeight[ori_index] )
         {
         all_sum_col += mean_colshift[ori_index];
         all_sum_row += mean_rowshift[ori_index];
         all_count ++;
         
         if(GridPT3[ori_index].ortho_ncc > 0.2)
         {
         double temp_slope = atan2(mean_rowshift[ori_index],mean_colshift[ori_index])*180/3.141592;
         if(temp_slope < 0)
         temp_slope += 360;
         all_mean_slope += temp_slope;
         mean_slope_count++;
         }
         }
         
         }
         //fprintf(fid_col,"\n");
         //fprintf(fid_row,"\n");
         }
         //fclose(fid_col);
         //fclose(fid_row);
         
         average_col = all_sum_col/all_count;
         average_row = all_sum_row/all_count;
         all_average_slope = all_mean_slope/mean_slope_count;
         shift_angle_th = all_average_slope;
         
         printf("average %f\t%f\t%f\n",average_col,average_row,shift_angle_th);
         
         }
         */
        int shift_max_pixel = (int)(((double)(proinfo.SDM_AS * proinfo.SDM_days) / (proinfo.resolution)) );
        
        printf("level %d\tshift_max_pixel %d\n",proinfo.pyramid_level,shift_max_pixel);
#pragma omp parallel for schedule(guided)
        for(int iter_count = 0 ; iter_count < gridsize.height*gridsize.width ; iter_count ++)
        {
            //printf("id %d start hist allocation\n",iter_count);
            int r = (int)(floor(iter_count/gridsize.width));
            int c = iter_count % gridsize.width;
            int ori_index = r*gridsize.width + c;
            //printf("id %d %d start hist allocation\n",iter_count,ori_index);
            if(ori_index >= 0 && r >= 0 && r < gridsize.height && c >= 0 && c < gridsize.width && ori_index < gridsize.height*gridsize.width)
            {
                double min_ortho_th = max_min_ortho_th;
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
                /*
                bool check_while = false;
                int final_select_pixel = 0;
                while(kernal_size < max_kernel_size && !check_while)
                {
                    while(min_ortho_th > 0.0 && !check_while)
                    {
                        int total_select_pixel = 0;
                        for(int k = -kernal_size ; k <= kernal_size ; k++)
                        {
                            for(int j=-kernal_size ; j <= kernal_size ; j++)
                            {
                                int index = (r+k)*gridsize.width + (c+j);
                                if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width && index >= 0 && index < gridsize.height*gridsize.width)
                                {
                                    if(GridPT3[index].ortho_ncc > min_ortho_th)
                                        total_select_pixel++;
                                }
                            }
                        }
                        if(total_select_pixel > (2*kernal_size + 1)*(2*kernal_size + 1)*0.5)
                            check_while = true;
                        else
                            min_ortho_th = min_ortho_th - 0.02;
                    }
                    
                    if(!check_while)
                        kernal_size++;
                    //final_select_pixel = total_select_pixel;
                    //printf("total_select %d\t%f\tkernel_size %d\n",total_select_pixel,kernal_size,(2*kernal_size + 1)*(2*kernal_size + 1)*0.5);
                }
                */
                //if( (GridPT3[ori_index].ortho_ncc < 0.5 && pyramid_step <= 1) || pyramid_step > 1)
                {
                    double sum_col = 0.;
                    double sum_row = 0.;
                    int sum_count = 0;
                    
                    double save_col[10000] = {0.};
                    double save_row[10000] = {0.};
                    double save_diff[10000] = {0.};
                    double save_col_min = 10000.;
                    double save_col_max = -10000.;
                    double save_row_min = 10000.;
                    double save_row_max = -10000.;
                    
                    double save_slope[10000] = {0.};
                    double save_roh[10000] = {-1.};
                    
                    int* hist_slope = (int*)calloc(sizeof(int),slope_step);
                    
                    int save_count = 0;
                    
                    
                    
                    int hist_col[6000] = {0};
                    int hist_col_id[6000];
                    int hist_row[6000] = {0};
                    int hist_row_id[6000];
                    
                    for(int k=0;k<3000;k++)
                    {
                        hist_col_id[k] = k;
                        hist_row_id[k] = k;
                        
                        //if(k < 8)
                        //    hist_id[k] = k;
                    }
                    
                    
                    for(int k = -kernal_size ; k <= kernal_size ; k++)
                    {
                        for(int j=-kernal_size ; j <= kernal_size ; j++)
                        {
                            int index = (r+k)*gridsize.width + (c+j);
                            if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width && index >= 0 && index < gridsize.height*gridsize.width /*&&
                               GridPT3[index].ortho_ncc > min_ortho_th*/)
                            {
                                save_col[save_count] = GridPT3[index].col_shift;
                                save_row[save_count] = GridPT3[index].row_shift;
                                save_diff[save_count] = sqrt(k*k + j*j);
                                save_slope[save_count] = atan2(save_row[save_count],save_col[save_count])*180/3.141592;
                                
                                //if(pyramid_step <= 5)
                                {
                                    save_roh[save_count] = pow(3.0,GridPT3[index].ortho_ncc*10);
                                    /*if(GridPT3[index].ortho_ncc > 0)
                                        save_roh[save_count] = GridPT3[index].ortho_ncc;
                                    else
                                        save_roh[save_count] = 0.001;
                                     */
                                }
                                //else
                                //    save_roh[save_count] = 1.0;
                                
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
                                {
                                    hist_slope[t_hist-1]++;
                                    //printf("t_hist %f\t%f\t%d\n",save_slope[save_count],da,t_hist);
                                }
                                
                                
                                t_hist = (int)(save_col[save_count]);
                                if(abs(t_hist) < shift_max_pixel)
                                    hist_col[t_hist+shift_max_pixel]++;
                                //if(abs(t_hist) < 3000)
                                //    hist_col[t_hist+3000]++;
                                //else
                                //    printf("less than -2000 col shift\n");
                                
                                t_hist = (int)(save_row[save_count]);
                                if(abs(t_hist) < shift_max_pixel)
                                    hist_row[t_hist+shift_max_pixel]++;
                                //if(abs(t_hist) < 3000)
                                //    hist_row[t_hist+3000]++;
                                //else
                                //    printf("less than -2000 col shift\n");
                                
                                save_count++;
                                
                                //sum_col += mean_colshift[index];
                                //sum_row += mean_rowshift[index];
                                //sum_count++;
                            }
                        }
                    }
                    //printf("id %d end hist allocation\n",iter_count);
                    save_count--;
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
                        int max_hist_col_pos, max_hist_row_pos;
                        int min_hist_col_pos, min_hist_row_pos;
                        int min_hist = 10000;
                        int max_hist_pos;
                        int min_hist_k;
                        int while_iter = 0;
                        
                        //bubble_sort_int(hist_slope,hist_id,8);
                        //max_hist = hist_slope[0];
                        
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
                        /*
                         bubble_sort_int(hist_col,hist_col_id,100);
                         max_hist_col = hist_col[0];
                         bubble_sort_int(hist_row,hist_row_id,100);
                         max_hist_row = hist_row[0];
                         //printf("end sort\n");
                         */
                        double col_sum = 0;
                        double row_sum = 0;
                        double checked_count = 0;
                        
                        double selected_col[10000] = {0.};
                        double selected_row[10000] = {0.};
                        double selected_diff_col[10000] = {0.};
                        double selected_diff_row[10000] = {0.};
                        double selected_roh[10000] = {-1.0};
                        
                        int total_seleted_count_col = 0;
                        int total_seleted_count_row = 0;
                        double p = 1.5;
                        
                        for(int k= 0; k<save_count ; k++)
                        {
                            int t_hist_col = (int)(save_col[k]);
                            int t_hist_row = (int)(save_row[k]);
                            int t_hist_slope = (int)(save_slope[k]/da);
                            
                            if((int)pyramid_step >= 0 )
                            {
                                //column
                                selected_col[total_seleted_count_col] = save_col[k];
                                selected_roh[total_seleted_count_col] = save_roh[k];
                                if( save_diff[k] == 0)
                                    selected_diff_col[total_seleted_count_col] = 0.1;
                                else
                                {
                                    if( t_hist_col+3000 == max_hist_col_pos /*&& t_hist_row+1000 == max_hist_row_pos*/)
                                        selected_diff_col[total_seleted_count_col] = 0.1;//save_diff[k]/2.0;
                                    else
                                    {
                                        if(t_hist_slope == max_hist_pos)
                                        {
                                            selected_diff_col[total_seleted_count_col] = 0.1;//save_diff[k]/2.0;
                                            
                                            //angle_pts++;
                                            
                                        }
                                        else
                                        {
                                            selected_diff_col[total_seleted_count_col] = save_diff[k]*2.0;
                                        }
                                    }
                                }
                                total_seleted_count_col ++;
                                
                                //row
                                selected_row[total_seleted_count_row] = save_row[k];
                                if( save_diff[k] == 0)
                                    selected_diff_row[total_seleted_count_row] = 0.1;
                                else
                                {
                                    if( t_hist_row+3000 == max_hist_row_pos /*&& t_hist_col+1000 == max_hist_col_pos*/)
                                        selected_diff_row[total_seleted_count_row] = 0.1;//save_diff[k]/2.0;
                                    else
                                    {
                                        if(t_hist_slope == max_hist_pos)
                                        {
                                            selected_diff_row[total_seleted_count_row] = 0.1;//save_diff[k]/2.0;
                                            
                                            //angle_pts++;
                                            
                                        }
                                        else
                                        {
                                            selected_diff_row[total_seleted_count_row] = save_diff[k]*2.0;
                                        }
                                    }
                                }
                                total_seleted_count_row ++;
                            }
                            else
                            {
                                if(min_hist_col_pos == max_hist_col_pos || min_hist_row_pos == max_hist_row_pos)
                                {
                                    selected_col[total_seleted_count_col] = save_col[k];
                                    selected_row[total_seleted_count_row] = save_row[k];
                                    selected_roh[total_seleted_count_col] = save_roh[k];
                                    
                                    if( save_diff[k] == 0)
                                    {
                                        selected_diff_col[total_seleted_count_col] = 0.1;
                                        selected_diff_row[total_seleted_count_row] = 0.1;
                                    }
                                    else
                                    {
                                        /*if(t_hist_slope == max_hist_pos)
                                         {
                                         selected_diff_col[total_seleted_count_col] = save_diff[k]/2.0;
                                         selected_diff_row[total_seleted_count_row] = save_diff[k]/2.0;
                                         
                                         angle_pts++;
                                         
                                         }
                                         else*/
                                        {
                                            selected_diff_col[total_seleted_count_col] = save_diff[k];
                                            selected_diff_row[total_seleted_count_row] = save_diff[k];
                                        }
                                    }
                                    total_seleted_count_col ++;
                                    total_seleted_count_row ++;
                                }
                                else if((t_hist_col+3000 != min_hist_col_pos || t_hist_row+3000 != min_hist_row_pos) )
                                {
                                    selected_col[total_seleted_count_col] = save_col[k];
                                    selected_row[total_seleted_count_row] = save_row[k];
                                    selected_roh[total_seleted_count_col] = save_roh[k];
                                    
                                    if( save_diff[k] == 0)
                                    {
                                        selected_diff_col[total_seleted_count_col] = 0.1;
                                        selected_diff_row[total_seleted_count_row] = 0.1;
                                    }
                                    else
                                    {
                                        /*if(t_hist_slope == max_hist_pos)
                                         {
                                         selected_diff_col[total_seleted_count_col] = save_diff[k]/2.0;
                                         selected_diff_row[total_seleted_count_row] = save_diff[k]/2.0;
                                         
                                         angle_pts++;
                                         
                                         }
                                         else*/
                                        {
                                            selected_diff_col[total_seleted_count_col] = save_diff[k];
                                            selected_diff_row[total_seleted_count_row] = save_diff[k];
                                        }
                                    }
                                    total_seleted_count_col ++;
                                    total_seleted_count_row ++;
                                }
                            }
                        }
                        
                        //col
                        if(total_seleted_count_col > 0)
                        {
                            double sum1, sum2, sum1_r, sum2_r;;
                            
                            sum1 = 0;
                            sum2 = 0;
                            sum1_r = 0;
                            sum2_r = 0;
                            
                            if((int)pyramid_step >= 0)
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
                            
                            //interpolated_pts++;
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
                            
                            //interpolated_pts++;
                        }
                        /*else
                         {
                         printf("checked count is zero\n");
                         non_iterpolated_pts++;
                         checked_count = 0;
                         for(int k= 0; k<save_count ; k++)
                         {
                         int t_hist_col = (int)(save_col[k]);
                         int t_hist_row = (int)(save_row[k]);
                         int t_hist_slope = (int)(save_slope[k]/da);
                         
                         if((t_hist_col+50 == max_hist_col_pos) && (t_hist_row+50 == max_hist_row_pos) && (t_hist_slope == max_hist_pos))
                         {
                         col_sum += save_col[k];
                         row_sum += save_row[k];
                         
                         checked_count ++;
                         
                         printf("selected col id t_hist save_col col_sum %d\t%d\t%f\t%f\n",k,t_hist_col,save_col[k],col_sum);
                         }
                         }
                         }*/
                    }
                    //printf("end DWF\n");
                    //GridPT3[ori_index].col_shift = sum_col/(double)sum_count;
                    //GridPT3[ori_index].row_shift = sum_row/(double)sum_count;
                    
                    
                    /*if(!m_bHeight[ori_index] && pyramid_step >= 4 && iteration < 3)
                     {
                     GridPT3[ori_index].col_shift = average_col;
                     GridPT3[ori_index].row_shift = average_row;
                     }
                     */
                    free(hist_slope);
                    //printf("free hist_slope\n");
                }
            }
        }
        
        
        printf("start assign\n");
#pragma omp parallel for schedule(guided)
        for(int iter_count = 0 ; iter_count < gridsize.height*gridsize.width ; iter_count ++)
        {
            //printf("id %d start hist allocation\n",iter_count);
            int r = (int)(floor(iter_count/gridsize.width));
            int c = iter_count % gridsize.width;
            int ori_index = r*gridsize.width + c;
        
            double sum_c = 0;
            double sum_r = 0;
            int count_cr = 0;
        
            
            double min_ortho_th = max_min_ortho_th;
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
            /*
            bool check_while = false;
            int final_select_pixel = 0;
            while(kernal_size < max_kernel_size && !check_while)
            {
                while(min_ortho_th > 0.0 && !check_while)
                {
                    int total_select_pixel = 0;
                    for(int k = -kernal_size ; k <= kernal_size ; k++)
                    {
                        for(int j=-kernal_size ; j <= kernal_size ; j++)
                        {
                            int index = (r+k)*gridsize.width + (c+j);
                            if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width && index >= 0 && index < gridsize.height*gridsize.width)
                            {
                                if(GridPT3[index].ortho_ncc > min_ortho_th)
                                    total_select_pixel++;
                            }
                        }
                    }
                    if(total_select_pixel > (2*kernal_size + 1)*(2*kernal_size + 1)*0.5)
                        check_while = true;
                    else
                        min_ortho_th = min_ortho_th - 0.02;
                }
                
                if(!check_while)
                    kernal_size++;
                //final_select_pixel = total_select_pixel;
                //printf("total_select %d\t%f\tkernel_size %d\n",total_select_pixel,kernal_size,(2*kernal_size + 1)*(2*kernal_size + 1)*0.5);
            }
            */
            //if( (GridPT3[ori_index].ortho_ncc < 0.5 && pyramid_step <= 1) || pyramid_step > 1)
            {
                for(int k = -kernal_size ; k <= kernal_size ; k++)
                {
                    for(int j=-kernal_size ; j <= kernal_size ; j++)
                    {
                        int index = (r+k)*gridsize.width + (c+j);
                        if(r+k >= 0 && r+k < gridsize.height && c+j >= 0 && c+j < gridsize.width  && index >= 0 && index < (gridsize.height)*(gridsize.width))
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
            
        }
        //printf("end assign\n");
        //free(mean_colshift);
        //free(mean_rowshift);
    }
    
    printf("before free col_shift\n");
    free(temp_col_shift);
    printf("free col_shift\n");
    free(temp_row_shift);
    printf("free row_shift\n");
    
    
    printf("end updating grid set shift %d\t%d\t%d\n",interpolated_pts,non_iterpolated_pts,angle_pts);
    
    
    
    
}

void RemoveFiles_SDM(char *save_path, char *lfilename, char *rfilename, int py_level, bool flag)
{
    int status;
    int count;
    
    char *filename_py;
    char t_str[500];
    int start_lv, end_lv;
    
    start_lv    = py_level;
    
    if(flag)
        end_lv        = py_level - 1;
    else
        end_lv        = -2;
    
    if(py_level == 0)
    {
        start_lv    = 0;
        end_lv        = -2;
    }
    
    for(count = py_level ; count >= 0 ; count--)
    {
        bool bfile = false;
        
        filename_py        = GetFileName(lfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(rfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(lfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(rfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(lfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(rfilename);
        filename_py        = remove_ext(filename_py);
        sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,count);
        status = remove(t_str);
        
        filename_py        = GetFileName(lfilename);
        filename_py        = remove_ext(filename_py);
        
        sprintf(t_str,"%s/%s.raw",save_path,filename_py);
        
        status = remove(t_str);
        
        filename_py        = GetFileName(rfilename);
        filename_py        = remove_ext(filename_py);
        
        sprintf(t_str,"%s/%s.raw",save_path,filename_py);
        status = remove(t_str);
    }
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
    int index_file;
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
    
#pragma omp parallel for private(index_file) schedule(guided)
    for(index_file = 0 ; index_file < row_end*col_end ; index_file++)
    {
        int row,col;
        
        row = (int)(floor(index_file/col_end));
        col = index_file%col_end;
        
        FILE *pfile;
        char t_str[500];
        //sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info.save_filepath,row,col,find_level,find_iter);
        sprintf(t_str,"%s/txt/col_shift_%d_%d_%d_%d.txt",info.save_filepath,find_level,find_iter,row,col);
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
                    int iter;
                    for(iter=0;iter<header_line;iter++)
                    {
                        int t_row,t_col,t_level,t_col_size,t_row_size;
                        double t_grid_size;
                        double t_boundary[4];
                        fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n",
                               &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&t_col_size,&t_row_size);
                        if(iter == header_line-1)
                        {
                            grid_size = t_grid_size;
                            t_boundary[2] = t_boundary[0] + t_grid_size*t_col_size;
                            t_boundary[3] = t_boundary[1] + t_grid_size*t_row_size;
                            
#pragma omp critical
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
                        }
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
#pragma omp parallel for private(index_file) schedule(guided)
    for(index_file = 0 ; index_file < row_end*col_end ; index_file++)
    {
        int row,col;
        FILE *pfile;
        char t_str[500];
        
        row = (int)(floor(index_file/col_end));
        col = index_file%col_end;
        
        sprintf(t_str,"%s/txt/col_shift_%d_%d_%d_%d.txt",info.save_filepath,find_level,find_iter,row,col);
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
                //if(p_hfile)
                {
                    int iter;
                    char hv_t_str[500];
                    int row_size,col_size;
                    double t_boundary[4];
                    for(iter=0;iter<header_line;iter++)
                    {
                        int t_row,t_col,t_level;
                        double t_grid_size;
                        
                        fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n",
                               &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
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
