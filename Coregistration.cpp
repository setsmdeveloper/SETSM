//
//  Coregistration.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "Coregistration.hpp"

//Image Coregistration
double** ImageCoregistration(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt, D2DPOINT *adjust_std, bool* cal_check)
{
    ProInfo *proinfo = (ProInfo*)malloc(sizeof(ProInfo));
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    
    double **ImageAdjust_coreg = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
    
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    FILE *time_fid;
    
    total_ST = time(0);
    
    if(OpenProject(_filename,proinfo,args))
    {
        bool check_cal = false;
        char check_file[500];
        sprintf(check_file, "%s/%s_dmag.tif", proinfo->save_filepath, proinfo->Outputpath_name);
        printf("dmag file %s\n",check_file);
        FILE* pcheckFile;
        pcheckFile = fopen(check_file,"r");
        if(!pcheckFile)
            check_cal = true;
        
        if(!check_cal && args.check_sdm_ortho == 2)
        {
            printf("SDM has been already processed. *_dmag.tif file exists!!\n");
            *cal_check = false;
        }
        else
        {
            *cal_check = true;
            char temp_filepath[500];
            
            int check_folder = 1;
            
            sprintf(proinfo->save_filepath,"%s",args.Outputpath);
            
            int status;
            status = mkdir(proinfo->save_filepath,0777);
            if (opendir(proinfo->save_filepath) == NULL)
            {
                if (status == -1)
                {
                    printf("Outpath of '%s' cannot make, please check outpath!!\n",proinfo->save_filepath);
                    exit(1);
                }
            }
            
            sprintf(temp_filepath,"%s/txt",proinfo->save_filepath);
            mkdir(temp_filepath,0777);
            sprintf(temp_filepath,"%s/tif",proinfo->save_filepath);
            mkdir(temp_filepath,0777);
            sprintf(proinfo->tmpdir,"%s/tmp",proinfo->save_filepath);
            mkdir(proinfo->tmpdir,0777);
            
            uint8 *image_bits = (uint8*)malloc(sizeof(uint8)*proinfo->number_of_images);
            double *ortho_minX = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_maxY = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_dx = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_dy = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double **Boundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
            CSize tmp_datasize;
            long int cols[2], rows[2];
            double *GridSize_width = (double*)calloc(sizeof(double),proinfo->number_of_images);
            double *GridSize_height = (double*)calloc(sizeof(double),proinfo->number_of_images);
            //CSize *Grid_size = (CSize*)calloc(sizeof(CSize),proinfo->number_of_images);
            
            double Sum_grid = 0;
            NCCflag ncc_flag;
            uint16 **OriImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
            CSize *OriImagesizes = (CSize*)malloc(sizeof(CSize)*proinfo->number_of_images);
            double **ImageBoundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
            
            ncc_flag.rotate_flag     = 1;       ncc_flag.multi_flag      = 1;       ncc_flag.multi_flag_sum  = 1;       ncc_flag.inter_flag      = 1;
            ncc_flag.weight_flag     = 1;
            int py_level = proinfo->pyramid_level;
            
            if(args.check_sdm_ortho > 0)
            {
                if(py_level < 2)
                    py_level = 2;
            }
            
            for(int i=0;i<proinfo->number_of_images;i++)
            {
                Boundary[i] = (double*)calloc(sizeof(double),4);
                
                proinfo->check_selected_image[i] = true;
                
                ImageAdjust_coreg[i] = (double*)malloc(sizeof(double)*2);
                ImageAdjust_coreg[i][0] = 0.0;
                ImageAdjust_coreg[i][1] = 0.0;
                
                image_bits[i] = ReadGeotiff_bits(args.Image[i]);
                OriImagesizes[i] = ReadGeotiff_info_dxy(args.Image[i],&ortho_minX[i],&ortho_maxY[i],&ortho_dx[i],&ortho_dy[i]);
                /*if(i == 1)
                {
                    ortho_minX[i] += 20;
                    ortho_maxY[i] += 20;
                }
                 */
                ImageBoundary[i] = (double*)malloc(sizeof(double)*4);
                ImageBoundary[i][0] = ortho_minX[i];
                ImageBoundary[i][1] = ortho_maxY[i] - ortho_dy[i]*OriImagesizes[i].height;
                ImageBoundary[i][2] = ortho_minX[i] + ortho_dx[i]*OriImagesizes[i].width;
                ImageBoundary[i][3] = ortho_maxY[i];
                
                
                if(i == 0)
                {
                    Sum_grid = ((ortho_dx[i] + ortho_dy[i])/2.0);
                    
                    Boundary[i][0] = ImageBoundary[i][0];
                    Boundary[i][1] = ImageBoundary[i][1];
                    Boundary[i][2] = ImageBoundary[i][2];
                    Boundary[i][3] = ImageBoundary[i][3];
                }
                else
                {
                    if(Boundary[0][0] < ImageBoundary[i][0])
                        Boundary[i][0] = ImageBoundary[i][0];
                    else
                        Boundary[i][0] = Boundary[0][0];
                    
                    if(Boundary[0][1] < ImageBoundary[i][1])
                        Boundary[i][1] = ImageBoundary[i][1];
                    else
                        Boundary[i][1] = Boundary[0][1];
                    
                    if(Boundary[0][2] > ImageBoundary[i][2])
                        Boundary[i][2] = ImageBoundary[i][2];
                    else
                        Boundary[i][2] = Boundary[0][2];
                    
                    if(Boundary[0][3] > ImageBoundary[i][3])
                        Boundary[i][3] = ImageBoundary[i][3];
                    else
                        Boundary[i][3] = Boundary[0][3];
                }
                
                if(args.check_boundary)
                {
                    if(Boundary[i][0] < args.Min_X)
                        Boundary[i][0] = args.Min_X;
                    if(Boundary[i][1] < args.Min_Y)
                        Boundary[i][1] = args.Min_Y;
                    if(Boundary[i][2] > args.Max_X)
                        Boundary[i][2] = args.Max_X;
                    if(Boundary[i][3] > args.Max_Y)
                        Boundary[i][3] = args.Max_Y;
                    
                    Boundary[i][0] = args.Min_X;
                    Boundary[i][1] = args.Min_Y;
                    Boundary[i][2] = args.Max_X;
                    Boundary[i][3] = args.Max_Y;
                }
                
                GridSize_width[i] = Boundary[i][2] - Boundary[i][0];
                GridSize_height[i] = Boundary[i][3] - Boundary[i][1];
                
                printf("ID %d\tboundary = %f\t%f\t%f\t%f\t%f\t%f\n",i,Boundary[i][0],Boundary[i][1],Boundary[i][2],Boundary[i][3],GridSize_width[i],GridSize_height[i]);
                
                //printf("ID %d\t %s\t%d\t%d\t%f\t%f\t%f\t%f\n",i,args.Image[i],ortho_size[i].width,ortho_size[i].height,ortho_minX[i],ortho_maxY[i],ortho_dx[i],ortho_dy[i]);
                cols[0] = 0;
                cols[1] = OriImagesizes[i].width;
                
                rows[0] = 0;
                rows[1] = OriImagesizes[i].height;
                
                switch(image_bits[i])
                {
                    case 8:
                    {
                        uint8 type;
                        uint8* data8 = Readtiff_T(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],type);
                        long int data_size8 = (long int)OriImagesizes[i].width*(long int)OriImagesizes[i].height;
                        OriImages[i] = (uint16*)malloc(sizeof(uint16)*data_size8);
                        #pragma omp parallel for schedule(guided)
                        for(long int index = 0 ; index < data_size8 ; index++)
                            OriImages[i][index] = data8[index];
                        free(data8);
                    }
                        break;
                    case 12:
                    {
                        uint16 type;
                        uint16* data16 = Readtiff_T(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],type);
                        long int data_size16 = (long int)OriImagesizes[i].width*(long int)OriImagesizes[i].height;
                        OriImages[i] = (uint16*)malloc(sizeof(uint16)*data_size16);
                        #pragma omp parallel for schedule(guided)
                        for(long int index = 0 ; index < data_size16 ; index++)
                            OriImages[i][index] = data16[index];
                        free(data16);
                    }
                        
                        
                        //OriImages[i] = Readtiff(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],false);
                        break;
                }
                //OriImages[i] = Readtiff(args.Image[i], &OriImagesizes[i], cols, rows, &tmp_datasize, false);
                
                printf("ID %d\t %s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,args.Image[i],OriImagesizes[i].width,OriImagesizes[i].height,ortho_minX[i],ortho_maxY[i],ortho_dx[i],ortho_dy[i],ImageBoundary[i][0],ImageBoundary[i][1],ImageBoundary[i][2],ImageBoundary[i][3]);
            }
            
            
            /*
             ImageAdjust_coreg[1][0] = 2.33;
             ImageAdjust_coreg[1][1] = -3.647;
             
             ImageAdjust_coreg[2][0] = 15.0;
             ImageAdjust_coreg[2][1] = -15.0;
             */
            
            
            int *Grid_space = (int*)calloc(sizeof(int),py_level+1);
            for(int i = 0 ; i <= py_level ; i++)
            {
                //if(i == 0)
                //    Grid_space[i] = ceil(Sum_grid)*pwrtwo(i);
                //else
                    Grid_space[i] = ceil(Sum_grid)*pwrtwo(py_level+2);
                printf("image coregistration gridspace %d\t%d\n",i,Grid_space[i]);
            }
            /*
            for(int i=0;i<proinfo->number_of_images;i++)
            {
                Grid_size[i].width = floor(GridSize_width[i]/Grid_space[i]);
                Grid_size[i].height = floor(GridSize_height[i]/Grid_space[i]);
                
                printf("Grid_size %d\t%d\n",Grid_size[i].width,Grid_size[i].height);
            }
            */
            char **Subsetfilename;
            CSize **data_size_lr;
            char save_file[500];
            char out_file[500];
            char *filename;
            
            Subsetfilename = (char**)malloc(sizeof(char*)*proinfo->number_of_images);
            data_size_lr = (CSize**)malloc(sizeof(CSize*)*proinfo->number_of_images);
            
            for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
            {
                sprintf(proinfo->Imagefilename[ti],"%s",args.Image[ti]);
                //printf("image %d\t%s\n",ti,proinfo->Imagefilename[ti]);
                Subsetfilename[ti] = (char*)malloc(sizeof(char)*500);
                filename = GetFileName(args.Image[ti]);
                filename = remove_ext(filename);
                sprintf(Subsetfilename[ti],"%s/%s_subset.raw",proinfo->tmpdir,filename);
                //printf("subsetfilename %s\n",Subsetfilename[ti]);
                
                data_size_lr[ti] = (CSize*)malloc(sizeof(CSize)*(py_level+1));
                SetPySizes(data_size_lr[ti], OriImagesizes[ti], py_level);
                //for (int ttt = 0 ; ttt < py_level+1 ;ttt++)
                //    printf("data_size %d\t%d\n",data_size_lr[ti][ttt].width,data_size_lr[ti][ttt].height);
            }
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nbefor preprocessing %f\n\n",total_gap);
            total_ST = time(0);
            
            Preprocessing_Coreg(proinfo,proinfo->tmpdir,OriImages,Subsetfilename,py_level,OriImagesizes,data_size_lr);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\npreprocessing %f\n\n",total_gap);
            total_ST = time(0);
            
            uint16 **SubImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
            
            //printf("overlapped area %f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3],GridSize_width,GridSize_height,Grid_space,Grid_size.width,Grid_size.height,total_grid_counts);
            
            
            double *avg_roh = (double*)calloc(sizeof(double),proinfo->number_of_images);
            int *RA_iter_counts;
            
            for(int level = py_level ; level >= 0 ; level --)
            {
                
                printf("Processing level %d\n",level);
                printf("level\tImage ID\trow(pixel)\tcolumn(pixel)\tTy(meter)\tTx(meter)\tGCPS #\tavg_roh\t# of iteration\n");
                
                for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
                {
                    SubImages[ti]     = LoadPyramidImages(proinfo->tmpdir,Subsetfilename[ti],data_size_lr[ti][level],level);
                }
                
                
                
                int iter_counts;
                RA_iter_counts = CoregParam_Image(proinfo, level,py_level, ImageAdjust_coreg, ncc_flag,
                                                  15, SubImages, data_size_lr, ImageBoundary, ortho_dx, ortho_dy,
                                                  Grid_space,Boundary,proinfo->save_filepath,avg_roh,&iter_counts,adjust_std);
                
               
                
                for(int ti = 1 ; ti < proinfo->number_of_images ; ti ++)
                {
                    //printf("2 std %f\t%f\n",adjust_std[ti].m_X,adjust_std[ti].m_Y);
                    printf("%d\t%d\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%4.2f\t%d\n",level,ti,ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1],
                           -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],RA_iter_counts[ti],avg_roh[ti],iter_counts);
                    free(SubImages[ti]);
                }
                
                total_ET = time(0);
                total_gap = difftime(total_ET,total_ST);
                printf("\niter %d CoregParam_Image %f\n\n",level,total_gap);
                total_ST = time(0);
            }
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nCoregParam_Image %f\n\n",total_gap);
            total_ST = time(0);
            
            FILE* fid_out = NULL;
            sprintf(out_file,"%s/coreg_result.txt",proinfo->save_filepath);
            fid_out         = fopen(out_file,"w");
            fprintf(fid_out,"orthoimage name\tline(row) direction[pixel]\tsample(column) direction[pixel]\tTy[meter]\tTx[meter]\tavg_roh\n");
            for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
            {
                fprintf(fid_out,"%s\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%3.2f\n",proinfo->Imagefilename[ti],ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1],
                        -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],avg_roh[ti]);
            }
            fclose(fid_out);
            
            RemoveFiles(proinfo,proinfo->tmpdir,Subsetfilename,py_level,0);
            
            TransParam param;
            SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
            
            for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
            {
                FILE* p_GCP = NULL;
                sprintf(out_file,"%s/txt/GCPs_Image_ID_%d_level_0.txt",proinfo->save_filepath,ti);
                p_GCP = fopen(out_file,"r");
                if(p_GCP)
                {
                    D2DPOINT *ref_pts = (D2DPOINT*)malloc(sizeof(D2DPOINT)*RA_iter_counts[ti]);
                    D2DPOINT *tar_pts = (D2DPOINT*)malloc(sizeof(D2DPOINT)*RA_iter_counts[ti]);
                    int i = 0;
                    //printf("file %s\t%d\n",out_file,RA_iter_counts[ti]);
                    while( i < RA_iter_counts[ti] && (fscanf(p_GCP,"%f\t%f\t%f\t%f\n",&ref_pts[i].m_X,&ref_pts[i].m_Y,&tar_pts[i].m_X,&tar_pts[i].m_Y)) != EOF )
                    {
                        //printf("%f\t%f\t%f\t%f\n",ref_pts[i].m_X,ref_pts[i].m_Y,tar_pts[i].m_X,tar_pts[i].m_Y);
                        i++;
                    }
                    fclose(p_GCP);
                    
                    CSize GCP_size;
                    double GCP_grid;
                    if(gcp_opt == 1)
                        GCP_grid = (ortho_dx[ti] + ortho_dy[ti])/2.0*5;
                    else if(gcp_opt == 2)
                        GCP_grid = (ortho_dx[ti] + ortho_dy[ti])/2.0;
                    
                    GCP_size.width = floor(GridSize_width[ti]/GCP_grid);
                    GCP_size.height = floor(GridSize_height[ti]/GCP_grid);
                    
                    //printf("%f\t%d\t%d\t%d\n",GCP_grid,GCP_size.width ,GCP_size.height,i);
                    
                    long int data_length = (long int)(GCP_size.width)*(long int)(GCP_size.height);
                    unsigned char* GCP_value = (unsigned char*)calloc(sizeof(unsigned char),data_length);
                    for(int count = 0 ; count < data_length ; count ++)
                    {
                        GCP_value[count] = 255;
                    }
                    
                    int sel_count = 0;
                    for(int count = 0 ; count < i ;count++)
                    {
                        int pos_c = (int)((ref_pts[count].m_X - Boundary[ti][0])/GCP_grid);
                        int pos_r = (int)((Boundary[ti][3] - ref_pts[count].m_Y)/GCP_grid);
                        int index = pos_r*GCP_size.width + pos_c;
                        
                        //printf("id %d\t%d\t%d\n",count,pos_c,pos_r);
                        if(pos_c >= 0 && pos_c < GCP_size.width && pos_r >= 0 && pos_r < GCP_size.height)
                        {
                            GCP_value[index] = 0;
                            sel_count ++;
                        }
                    }
                    //printf("GCP_count %d\n",sel_count);
                    sprintf(out_file,"%s/GCPs_%d.tif",proinfo->save_filepath,ti);
                    WriteGeotiff(out_file, GCP_value, GCP_size.width, GCP_size.height, GCP_grid, Boundary[ti][0], Boundary[ti][3], param.projection, param.utm_zone, param.bHemisphere, 1);
             
                    char *Ifilename  = SetOutpathName(args.Image[ti]);
                    char *tmp_no_ext = remove_ext(Ifilename);
                    
                    sprintf(out_file,"%s/%s_coreg.tif",proinfo->save_filepath,tmp_no_ext);
                    
                    double left = - ImageAdjust_coreg[ti][1]*ortho_dx[ti];
                    double upper = ImageAdjust_coreg[ti][0]*ortho_dy[ti];
                    printf("coreg %s\t%d\t%d\t%f\t%f\n",out_file,OriImagesizes[ti].width,OriImagesizes[ti].height,left,upper);
                    
                    //WriteGeotiff(out_file, OriImages[ti], OriImagesizes[ti].width, OriImagesizes[ti].height, ortho_dx[ti], ImageBoundary[ti][0] +left , ImageBoundary[ti][3] + upper, param.projection, param.utm_zone, param.bHemisphere, 12);
                    
                    free(GCP_value);
                    free(ref_pts);
                    free(tar_pts);
                }
                free(OriImages[ti]);
                free(ImageBoundary[ti]);
            }
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nsave result %f\n\n",total_gap);
            
            
            free(Grid_space);
            free(ortho_dx);
            free(ortho_dy);
            free(avg_roh);
            free(RA_iter_counts);
            free(OriImages);
            free(OriImagesizes);
            
            free(ortho_minX);
            free(ortho_maxY);
            free(ImageBoundary);
            free(SubImages);
            
            free(Boundary);
            free(GridSize_width);
            free(GridSize_height);
            //free(Grid_size);
        }
    }
    
    return ImageAdjust_coreg;
}

void Preprocessing_Coreg(ProInfo *proinfo, char *save_path,uint16 **Oriimage, char **Subsetfile, uint8 py_level, CSize *Subsetsize, CSize **data_size_lr)
{
    char t_str[500];
    FILE *pFile_raw, *pFile_check_file;
    char* filename_py;
    
    uint8 count = 0;
    printf("start Preprocessing\n");
    for(count = 0; count<proinfo->number_of_images ; count++)
    {
        //pFile_raw   = fopen(proinfo->Imagefilename[count],"rb");
        filename_py     = GetFileName(Subsetfile[count]);
        filename_py     = remove_ext(filename_py);
        
        sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,py_level+1);
        pFile_check_file    = fopen(t_str,"rb");
        
        //if(pFile_raw && !pFile_check_file && proinfo->check_selected_image[count])
        {
            int i;
            
            FILE *pFile;
            uint16 **pyimg;
            CSize *data_size;
            
            pyimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
            data_size = (CSize*)malloc(sizeof(CSize)*(py_level+1));
            
            //printf("image %d\tlevel = 0\t width = %d\theight = %d\n",count,data_size_lr[count][0].width,data_size_lr[count][0].height);
            //fprintf(fid,"image %d\tlevel = 0\t width = %d\theight = %d\n",count,data_size_lr[count][0].width,data_size_lr[count][0].height);
            
            for(i=0;i<py_level+1;i++)
                data_size[i]        = data_size_lr[count][i];
            
            long int data_length = (long int)data_size[0].height*(long int)data_size[0].width;
            pyimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
            //memcpy(pyimg[0],Oriimage[count],sizeof(uint16)*(long)data_length);
            //fread(pyimg[0],sizeof(uint16),data_length,pFile_raw);
            
            sprintf(t_str,"%s/%s_py_0.raw",save_path,filename_py);
            pFile   = fopen(t_str,"wb");
            //fwrite(pyimg[0],sizeof(uint16),data_length,pFile);
            fwrite(Oriimage[count],sizeof(uint16),data_length,pFile);
            fclose(pFile);
            
            for(i=0;i<py_level;i++)
            {
                if(i==0)
                    pyimg[i+1] = CreateImagePyramid(Oriimage[count],data_size[i],9,(double)(1.5));
                else
                    pyimg[i+1] = CreateImagePyramid(pyimg[i],data_size[i],9,(double)(1.5));
                
                free(pyimg[i]);
                
                //printf("image %d\tlevel = %d\t width = %d\theight = %d\n",count,i+1,data_size_lr[count][i+1].width,data_size_lr[count][i+1].height);
                //fprintf(fid,"image %d\tlevel = %d\t width = %d\theight = %d\n",count,i+1,data_size_lr[count][i+1].width,data_size_lr[count][i+1].height);
                
                long int data_length_array[6];
                data_length_array[i] = (long int)data_size[i+1].height*(long int)data_size[i+1].width;
                
                sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,i+1);
                pFile   = fopen(t_str,"wb");
                fwrite(pyimg[i+1],sizeof(uint16),data_length_array[i],pFile);
                fclose(pFile);
                if(i == py_level-1)
                    free(pyimg[i+1]);
                
            }
            if(pyimg)
                free(pyimg);
            if(data_size)
                free(data_size);
        }
        
        //if(pFile_raw)
        //    fclose(pFile_raw);
        if(pFile_check_file)
            fclose(pFile_check_file);
    }
    
}

int* CoregParam_Image(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double **ImageAdjust, NCCflag _flag,
                      uint8 Template_size, uint16 **Images, CSize **Imagesizes, double **Boundary, double *grid_dx, double *grid_dy,
                      int* grid_space,double** over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std)
{
    int i;
    int *mp_iter_count = (int*)calloc(sizeof(int),proinfo->number_of_images);
    
    
    int reference_id = 0;
    CSize LImagesize;
    
    LImagesize.width  = Imagesizes[reference_id][Pyramid_step].width;
    LImagesize.height = Imagesizes[reference_id][Pyramid_step].height;
    
    double  left_IA[2];
    left_IA[0] = ImageAdjust[reference_id][0];
    left_IA[1] = ImageAdjust[reference_id][1];
    
    double **subA;
    double **TsubA;
    double **InverseSubA;
    
    FILE **fid_pts = (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    FILE **fid_stat= (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    
    FILE *fid_pts_pre;
    
    int ii,kk;
    
    subA    = (double**)malloc(9*sizeof(double*));
    TsubA   = (double**)malloc(6*sizeof(double*));
    InverseSubA = (double**)malloc(6*sizeof(double*));
    
    for(ii=0;ii<9;ii++)
    {
        subA[ii]    = (double*)malloc(6*sizeof(double));
        if(ii < 6)
        {
            TsubA[ii]       = (double*)malloc(9*sizeof(double));
            InverseSubA[ii] = (double*)malloc(6*sizeof(double));
        }
    }
    
    for(ii=0;ii<9;ii++)
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
    
    for(ii=0;ii<6;ii++)
        for(kk=0;kk<9;kk++)
            TsubA[ii][kk]       = subA[kk][ii];
    
    InverseSubA[0][0] =  0.555556; InverseSubA[0][1] =  0.000000; InverseSubA[0][2] =  0.000000; InverseSubA[0][3] = -0.333333; InverseSubA[0][4] =  0.000000; InverseSubA[0][5] = -0.333333;
    InverseSubA[1][0] =  0.000000; InverseSubA[1][1] =  0.166667; InverseSubA[1][2] =  0.000000; InverseSubA[1][3] =  0.000000; InverseSubA[1][4] =  0.000000; InverseSubA[1][5] =  0.000000;
    InverseSubA[2][0] =  0.000000; InverseSubA[2][1] =  0.000000; InverseSubA[2][2] =  0.166667; InverseSubA[2][3] =  0.000000; InverseSubA[2][4] =  0.000000; InverseSubA[2][5] =  0.000000;
    InverseSubA[3][0] = -0.333333; InverseSubA[3][1] =  0.000000; InverseSubA[3][2] =  0.000000; InverseSubA[3][3] =  0.500000; InverseSubA[3][4] =  0.000000; InverseSubA[3][5] =  0.000000;
    InverseSubA[4][0] =  0.000000; InverseSubA[4][1] =  0.000000; InverseSubA[4][2] =  0.000000; InverseSubA[4][3] =  0.000000; InverseSubA[4][4] =  0.250000; InverseSubA[4][5] =  0.000000;
    InverseSubA[5][0] = -0.333333; InverseSubA[5][1] =  0.000000; InverseSubA[5][2] =  0.000000; InverseSubA[5][3] =  0.000000; InverseSubA[5][4] =  0.000000; InverseSubA[5][5] =  0.500000;
    
    D2DPOINT **matched_MPs = (D2DPOINT**)malloc(sizeof(D2DPOINT*)*proinfo->number_of_images);
    uint8 **flag_MPs = (uint8**)malloc(sizeof(uint8*)*proinfo->number_of_images);
    
    CSize *grid_size = (CSize*)calloc(sizeof(CSize),proinfo->number_of_images);
    for(int i=0;i<proinfo->number_of_images;i++)
    {
        double GridSize_width = Boundary[i][2] - Boundary[i][0];
        double GridSize_height = Boundary[i][3] - Boundary[i][1];
        
        grid_size[i].width = floor(GridSize_width/grid_space[Pyramid_step]);
        grid_size[i].height = floor(GridSize_height/grid_space[Pyramid_step]);
        
        printf("Grid_size %d\t%d\n",grid_size[i].width,grid_size[i].height);
    }
    
    matched_MPs[reference_id] = (D2DPOINT*)malloc(sizeof(D2DPOINT)*grid_size[0].height*grid_size[0].width);
    flag_MPs[reference_id] = (uint8*)calloc(sizeof(uint8),grid_size[0].height*grid_size[0].width);
    
    
    bool check_top_level = false;
    //if(Pyramid_step < total_level)
    //    check_top_level = true;
    
    char temp_path[500];
    
    for(int ti = 1 ; ti < proinfo->number_of_images ; ti++)
    {
        long total_grid_counts = 0;
        D2DPOINT *MPs = NULL;
        if(check_top_level)
        {
            sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step+1);
            fid_pts_pre = fopen(temp_path,"r");
            float temp_v;
            while( fscanf(fid_pts_pre,"%f\t%f\t%f\t%f\n",&temp_v,&temp_v,&temp_v,&temp_v) != EOF )
                total_grid_counts++;
            
            fclose(fid_pts_pre);
            
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            int i = 0;
            fid_pts_pre = fopen(temp_path,"r");
            while( fscanf(fid_pts_pre,"%f\t%f\t%f\t%f\n",&MPs[i].m_X,&MPs[i].m_Y,&temp_v,&temp_v) != EOF )
            {
                //printf("mps id %d\t%f\t%f\n",i,MPs[i].m_X,MPs[i].m_Y);
                i++;
            }
            
            fclose(fid_pts_pre);
        }
        else
        {
            total_grid_counts = grid_size[ti].height*grid_size[ti].width;
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
    #pragma omp parallel for
            for(int row = 0 ; row < grid_size[ti].height ; row ++)
            {
                for(int col = 0 ; col < grid_size[ti].width ; col ++)
                {
                    int index = row*grid_size[ti].width + col;
                    MPs[index].m_X = over_Boundary[ti][0] + col*grid_space[Pyramid_step];
                    MPs[index].m_Y = over_Boundary[ti][1] + row*grid_space[Pyramid_step];
                }
            }
        }
        
        printf("ID %d\ttotal pts %d\n",ti,total_grid_counts);
        
        sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
        fid_pts[ti] = fopen(temp_path,"w");
        sprintf(temp_path,"%s/txt/CoregStat_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
        fid_stat[ti] = fopen(temp_path,"w");
        
        if(proinfo->check_selected_image[ti])
        {
            matched_MPs[ti] = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            flag_MPs[ti] = (uint8*)calloc(sizeof(uint8),total_grid_counts);
            
            bool check_stop = false;
            *iter_count = 1;
            while(!check_stop && *iter_count < 20)
            {
                uint8   Half_template_size;
                double b_factor;
                bool flag_boundary = false;
                int i;
                int count_pts = 0;
                double sum_weight_X     = 0;
                double sum_weight_Y     = 0;
                double sum_max_roh      = 0;
                double t_sum_weight_X       = 0;
                double t_sum_weight_Y       = 0;
                double t_sum_max_roh        = 0;
                double shift_X, shift_Y;
                
                //calculation image coord from object coord by RFM in left and right image
                b_factor             = pow(2.0,(2-Pyramid_step))*2;
                Half_template_size   = (int)(Template_size/2.0);
                
                F3DPOINT* save_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),total_grid_counts);
                
#pragma omp parallel for private(i,t_sum_weight_X,t_sum_weight_Y,t_sum_max_roh) reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
                for(i = 0; i<total_grid_counts ; i++)
                {
                    CSize RImagesize;
                    D2DPOINT Startpos;
                    Startpos.m_X = 0;
                    Startpos.m_Y = 0;
                    
                    RImagesize.width  = Imagesizes[ti][Pyramid_step].width;
                    RImagesize.height = Imagesizes[ti][Pyramid_step].height;
                    
                    D2DPOINT Left_Imagecoord, Right_Imagecoord,Left_Imagecoord_p, Right_Imagecoord_p;
                    Left_Imagecoord_p.m_X = (MPs[i].m_X - Boundary[reference_id][0])/grid_dx[reference_id];
                    Left_Imagecoord_p.m_Y = (Boundary[reference_id][3] - MPs[i].m_Y)/grid_dy[reference_id];
                    
                    Right_Imagecoord_p.m_X = (MPs[i].m_X - Boundary[ti][0])/grid_dx[ti] + ImageAdjust[ti][1];
                    Right_Imagecoord_p.m_Y = (Boundary[ti][3] - MPs[i].m_Y)/grid_dy[ti] + ImageAdjust[ti][0];
                    
                    Left_Imagecoord = OriginalToPyramid_single(Left_Imagecoord_p,Startpos,Pyramid_step);
                    Right_Imagecoord = OriginalToPyramid_single(Right_Imagecoord_p,Startpos,Pyramid_step);
                    
                    //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                    
                    if(   Left_Imagecoord.m_Y  > Half_template_size*b_factor    + 10                    && Left_Imagecoord.m_X  > Half_template_size*b_factor + 10
                       && Left_Imagecoord.m_Y  < LImagesize.height - Half_template_size*b_factor - 10    && Left_Imagecoord.m_X  < LImagesize.width - Half_template_size*b_factor - 10
                       && Right_Imagecoord.m_Y > Half_template_size*b_factor + 10                    && Right_Imagecoord.m_X > Half_template_size*b_factor + 10
                       && Right_Imagecoord.m_Y < RImagesize.height - Half_template_size*b_factor - 10    && Right_Imagecoord.m_X < RImagesize.width - Half_template_size*b_factor - 10)
                    {
                        double Left_X = Left_Imagecoord.m_X;
                        double Left_Y = Left_Imagecoord.m_Y;
                        double Right_X = Right_Imagecoord.m_X;
                        double Right_Y = Right_Imagecoord.m_Y;
                        int index_l = ((int)Left_Y)*LImagesize.width + (int)Left_X;
                        int index_r = ((int)Right_Y)*RImagesize.width + (int)Right_X;
                        double ori_diff;
                        if( (index_l > 0 && index_l < LImagesize.height*LImagesize.width) && (index_r > 0 && index_r < RImagesize.height*RImagesize.width) && Images[reference_id][index_l] > 0 && Images[ti][index_r] > 0 )
                        {
                            ori_diff = 0;
                            
                            if(postNCC_ortho(Pyramid_step, ori_diff, Left_Y,  Left_X, Right_Y, Right_X,
                                             subA,TsubA,InverseSubA,Template_size,_flag,1,LImagesize,RImagesize,Images[reference_id],Images[ti],&t_sum_weight_X,&t_sum_weight_Y,&t_sum_max_roh))
                            {
                                flag_MPs[reference_id][i] = true;
                                
                                flag_MPs[ti][i] = true;
                                
                                //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                                
                                
                                //printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n",ti,MPs[i].m_X,MPs[i].m_Y, matched_MPs[reference_id][count_pts].m_X,matched_MPs[reference_id][count_pts].m_Y,matched_MPs[ti][count_pts].m_X,matched_MPs[ti][count_pts].m_Y);
                                //#pragma omp critical
                                {
                                    sum_weight_X += t_sum_weight_X;
                                    sum_weight_Y += t_sum_weight_Y;
                                    sum_max_roh  += t_sum_max_roh;
                                    
                                    save_pts[i].m_X = t_sum_weight_X/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].m_Y = t_sum_weight_Y/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].flag = 1;
                                    
                                    //printf(" save pts %f\t%f\n",save_pts[i].m_X,save_pts[i].m_Y);
                                }
                                //#pragma omp atomic
                                count_pts++;
                            }
                        }
                    }
                }
                
                
                if(count_pts > 10)
                {
                    shift_X             = sum_weight_X/sum_max_roh*pwrtwo(Pyramid_step);
                    shift_Y             = sum_weight_Y/sum_max_roh*pwrtwo(Pyramid_step);
                    
                    double sum_var_x = 0;
                    double sum_var_y = 0;
                    
                    for(int c_i = 0 ; c_i < total_grid_counts ; c_i++)
                    {
                        if(save_pts[c_i].flag)
                        {
                            sum_var_x += (shift_X - save_pts[c_i].m_X)*(shift_X - save_pts[c_i].m_X);
                            sum_var_y += (shift_Y - save_pts[c_i].m_Y)*(shift_Y - save_pts[c_i].m_Y);
                            
                            //printf(" save pts %f\t%f\n",save_pts[c_i].m_X,save_pts[c_i].m_Y);
                        }
                    }
                    
                    adjust_std[ti].m_X = sqrt(sum_var_x/count_pts);
                    adjust_std[ti].m_Y = sqrt(sum_var_y/count_pts);
                    
                    //printf("1 std %f\t%f\t%d\n",adjust_std[ti].m_X,adjust_std[ti].m_Y,count_pts);
                    
                    if(fabs(shift_Y) < 0.01 && fabs(shift_X) < 0.01)
                        check_stop = true;
                    
                    //printf("ti %d\t%d\t%f\t%f\t%f\t%f\t%d\n",ti, *iter_count,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    fprintf(fid_stat[ti],"%d\t%d\t%f\t%f\t%f\t%f\t%d\n",*iter_count,ti,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    
                    shift_X             += ImageAdjust[ti][1];
                    shift_Y             += ImageAdjust[ti][0];
                    
                    ImageAdjust[ti][1]      = shift_X;
                    ImageAdjust[ti][0]      = shift_Y;
                    
                    
                    
                    
                    
                }
                else
                {
                    check_stop = true;
                }
                
                free(save_pts);
                
                (*iter_count)++;
                
                count_pts = 0;
                if(check_stop || *iter_count >= 20)
                {
                    for(int cc = 0 ; cc < total_grid_counts ; cc++)
                    {
                        if(flag_MPs[ti][cc])
                        {
                            matched_MPs[reference_id][cc].m_X = MPs[cc].m_X;
                            matched_MPs[reference_id][cc].m_Y = MPs[cc].m_Y;
                            
                            matched_MPs[ti][cc].m_X = MPs[cc].m_X + ImageAdjust[ti][1]*grid_dx[ti];
                            matched_MPs[ti][cc].m_Y = MPs[cc].m_Y - ImageAdjust[ti][0]*grid_dy[ti];
                            
                            //matched_MPs[ti][count_pts].m_X = Right_Imagecoord.m_X;
                            //matched_MPs[ti][count_pts].m_Y = Right_Imagecoord.m_Y;
                            //matched_MPs[reference_id][count_pts].m_X = Left_Imagecoord.m_X;
                            //matched_MPs[reference_id][count_pts].m_Y = Left_Imagecoord.m_Y;
                            
                            fprintf(fid_pts[ti],"%8.2f\t%8.2f\t%8.2f\t%8.2f\n",matched_MPs[reference_id][cc].m_X,matched_MPs[reference_id][cc].m_Y,matched_MPs[ti][cc].m_X,matched_MPs[ti][cc].m_Y);
                            count_pts++;
                        }
                    }
                    
                    avg_rho[ti] = sum_max_roh/(double)(count_pts);
                    
                    mp_iter_count[ti] = count_pts;
                    //printf("done %d\n",count_pts);
                    //*NumofPts = count_pts;
                    
                    fclose(fid_pts[ti]);
                    fclose(fid_stat[ti]);
                }
            }
        }
        
        free(MPs);
    }
    
    for(ii=0;ii<9;ii++)
    {
        free(subA[ii]);
        if(ii < 6)
        {
            free(TsubA[ii]);
            free(InverseSubA[ii]);
        }
    }
    free(grid_size);
    free(subA);
    free(TsubA);
    free(InverseSubA);
    
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(proinfo->check_selected_image[ti])
        {
            free(matched_MPs[ti]);
            free(flag_MPs[ti]);
        }
    }
    
    free(matched_MPs);
    free(flag_MPs);
    
    return mp_iter_count;
}

bool postNCC_ortho(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size,
                   NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh)
{
    
    // input order: NumOfPt, LeftImage, RightImage, Template_size, center_row_left, center_col_left, center_row_right, center_col_right, ncc_weight
    bool check_pt = false;
    int count_max_roh       = 0;
    
    uint8 rotate_flag,multi_flag,multi_flag_sum,inter_flag,weight_flag;
    
    //Image info input
    uint16 L_rowsize   = leftsize.height;
    uint16 L_colsize   = leftsize.width;
    
    uint16 R_rowsize   = rightsize.height;
    uint16 R_colsize   = rightsize.width;
    
    double t_weight_X   = 0;
    double t_weight_Y   = 0;
    double t_max_roh    = 0;
    double diff_theta;
    int mask_row, mask_col;
    
    int Half_template_size,half_mask_size;
    int count = 0;
    
    Half_template_size  = (int)(Template_size/2);
    half_mask_size      = 1;
    
    rotate_flag = _flag.rotate_flag;
    multi_flag  = _flag.multi_flag;
    multi_flag_sum  = _flag.multi_flag_sum;
    inter_flag  = _flag.inter_flag;
    weight_flag = _flag.weight_flag;
    
    
    diff_theta = Ori_diff;
    
    double *result_rho  = (double*)calloc(9,sizeof(double));
    double *XX          = (double*)calloc(6,sizeof(double));
    double *ATLT        = (double*)calloc(6,sizeof(double));
    int i, j, k;
    uint8 cell_count = 0;
    
    for(j=0;j<9;j++)
        result_rho[j]       = -1.00;
    
    for(mask_row = - half_mask_size ; mask_row <= half_mask_size ; mask_row++)
    {
        for(mask_col = - half_mask_size ; mask_col <= half_mask_size ; mask_col++)
        {
            double rot_theta = 0.0;
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
            int row, col;
            int N;
            double val1, val2, de, de2, ncc_1, ncc_2, ncc_3;
            double temp_rho;
            int grid_index;
            
            if(rotate_flag == 1)
                rot_theta = (double)(diff_theta*bin_angle*PI/180.0);
            
            
            for(row = -Half_template_size; row <= Half_template_size ; row++)
            {
                for(col = -Half_template_size; col <= Half_template_size ; col++)
                {
                    double radius  = sqrt((double)(row*row + col*col));
                    if(radius <= Half_template_size-1)
                    {
                        double pos_row_left      = (Left_CR + row);
                        double pos_col_left      = (Left_CC + col);
                        
                        double temp_col        = (cos(-rot_theta)*col - sin(-rot_theta)*row);
                        double temp_row        = (sin(-rot_theta)*col + cos(-rot_theta)*row);
                        double pos_row_right     = (Right_CR + temp_row + mask_row);
                        double pos_col_right     = (Right_CC + temp_col + mask_col);
                        
                        if(pos_row_right-3 >= 0 && pos_row_right+3 < R_rowsize && pos_col_right-3 >= 0 && pos_col_right+3 < R_colsize &&
                           pos_row_left-3 >= 0 && pos_row_left+3 < L_rowsize && pos_col_left-3 >= 0 && pos_col_left+3 < L_colsize)
                        {
                            //interpolate left_patch
                            double dx = pos_col_left - (int) (pos_col_left);
                            double dy = pos_row_left - (int) (pos_row_left);
                            double left_patch;
                            double right_patch;
                            double dxdy = dx * dy;
                            long int position = (long int) (pos_col_left) + (long int) (pos_row_left) * L_colsize;
                            
                            // Appears inter_flag is always == 1
                            if (inter_flag == 1) {
                                left_patch = (double) (_leftimage[position]) * (1 - dx - dy + dxdy) + (double) (_leftimage[position + 1]) * (dx - dxdy) +
                                (double) (_leftimage[position + L_colsize]) * (dy - dxdy) + (double) (_leftimage[position + 1 + L_colsize]) * (dxdy);
                            } else {
                                left_patch = (double) (_leftimage[position]);
                            }
                            //interpolate right_patch
                            dx = pos_col_right - (int) (pos_col_right);
                            dy = pos_row_right - (int) (pos_row_right);
                            dxdy = dx * dy;
                            position = (long int) (pos_col_right) + (long int) (pos_row_right) * R_colsize;
                            
                            // Appears inter_flag is always == 1
                            if (inter_flag == 1) {
                                right_patch = (double) (_rightimage[position]) * (1 - dx - dy + dxdy) + (double) (_rightimage[position + 1]) * (dx - dxdy) +
                                (double) (_rightimage[position + R_colsize]) * (dy - dxdy) + (double) (_rightimage[position + 1 + R_colsize]) * (dxdy);
                            } else {
                                right_patch = (double) (_rightimage[position]);
                            }
                            
                            if(left_patch > 1 && right_patch > 1)
                            {
                                //end
                                Count_N[0]++;
                                
                                Sum_LR            = Sum_LR + left_patch*right_patch;
                                Sum_L             = Sum_L  + left_patch;
                                Sum_R             = Sum_R  + right_patch;
                                Sum_L2            = Sum_L2 + left_patch*left_patch;
                                Sum_R2            = Sum_R2 + right_patch*right_patch;
                                
                                if(multi_flag == 1)
                                {
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
            }
            
            if(Count_N[0] > 0)
            {
                N               = Count_N[0];
                val1          = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
                val2          = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
                de            = sqrt(val1*val2);
                de2           = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
                if( val1*val2 > 0)
                    ncc_1           = de2/de;
                else
                    ncc_1           = -1.0;
                
                if(multi_flag == 1)
                {
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
                    
                }
                
                if(multi_flag == 1)
                {
                    if(Count_N[1] > 0 && Count_N[2] > 0)
                        temp_rho      = ((ncc_1 + ncc_2 + ncc_3)/3.0);
                    else if(Count_N[1] > 0)
                        temp_rho      = ((ncc_1 + ncc_2)/2.0);
                    else if(Count_N[2] > 0)
                        temp_rho      = ((ncc_1 + ncc_3)/2.0);
                    else
                        temp_rho        = ncc_1;
                }
                else
                {
                    temp_rho      = ncc_1;
                }
                
                
                grid_index           = (mask_row+1)*3 + (mask_col+1);
                if(grid_index < 9)
                    result_rho[grid_index] = temp_rho;
                cell_count++;
            }
            
            
        }
    }
    
    if(cell_count == 9)
    {
        double demnum;
        double max_X        = 100;
        double max_Y        = 100;
        double max_roh      = 0;
        bool find_index_1   = false;
        bool find_index_2   = false;
        bool find_index     = false;
        
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
        
        demnum      = -pow(XX[4],2.0) + 4*XX[3]*XX[5];
        if(demnum > 0 && XX[3] < 0)
        {
            max_X = (- 2*XX[5]*XX[1] + XX[2]*XX[4])/demnum;
            max_Y = (- 2*XX[2]*XX[3] + XX[1]*XX[4])/demnum;
            max_roh =  XX[0]                + XX[1]*max_X           + XX[2]*max_Y
            + XX[3]*max_X*max_X + XX[4]*max_X*max_Y + XX[5]*max_Y*max_Y;
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
                t_weight_X += max_X*max_roh;
                t_weight_Y += max_Y*max_roh;
                t_max_roh  += max_roh;
                
                check_pt = true;
                //printf("max XY %f\t%f\t%f\n",max_X,max_Y,max_roh);
            }
            
            
        }
    }
    free(result_rho);
    free(ATLT);
    free(XX);
    
    if(check_pt)
    {
        *sum_weight_X   = t_weight_X;
        *sum_weight_Y   = t_weight_Y;
        *sum_max_roh    = t_max_roh;
    }
    else
    {
        *sum_weight_X   = 0;
        *sum_weight_Y   = 0;
        *sum_max_roh    = 0;
    }
    
    return check_pt;
}

double *Readtiff_Coreg(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size)
{
    double *out;
    FILE *bin;
    int check_ftype = 1; // 1 = tif, 2 = bin
    TIFF *tif = NULL;
    char *ext;
    ext = strrchr(filename,'.');
    
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
            
            out             = (double*)malloc(sizeof(double)*data_size->height*data_size->width);
            
            scanline        = TIFFScanlineSize(tif);
            
            buf             = _TIFFmalloc(scanline);
            
            TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);
            
            for(s =0;s< nsamples;s++)
            {
                for (row=0;row<rows[0];row++)
                    TIFFReadScanline(tif,buf,row,s);
                for (row=rows[0];row<rows[1];row++)
                {
                    double* t_data;
                    TIFFReadScanline(tif,buf,row,s);
                    t_data = (double*)buf;
#pragma omp parallel for private(a) schedule(guided)
                    for(a = cols[0];a<cols[1];a++)
                    {
                        out[(row-rows[0])*data_size->width + (a-cols[0])] = t_data[a];
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
            double* t_data;
            
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
            
            data_size->width = Imagesize->width;
            data_size->height= Imagesize->height;
            
            long int data_length = (long int)data_size->height*(long int)data_size->width;
            
            out             = (double*)malloc(sizeof(double)*data_length);
            
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
                    t_data = (double*)buf;
                    
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
                    
                    /*
                     #pragma omp parallel for private(i,j) schedule(guided)
                     for (i=0;i<tileL;i++)
                     {
                     for (j=0;j<tileW;j++)
                     {
                     out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
                     }
                     }
                     */
                }
            }
            _TIFFfree(buf);
        }
        TIFFClose(tif);
    }
    
    
    return out;
}

double* LoadPyramidImages_double(char *save_path,char *subsetfile, CSize data_size, uint8 py_level)
{
    double *out = (double*)malloc(sizeof(double));
    FILE *pFile;
    char *filename_py;
    char t_str[500];
    
    filename_py     = GetFileName(subsetfile);
    filename_py     = remove_ext(filename_py);
    
    sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,py_level);
    pFile           = fopen(t_str,"rb");
    if(pFile)
    {
        out         = (double*)malloc(sizeof(double)*data_size.height*data_size.width);
        fread(out,sizeof(double),data_size.height*data_size.width,pFile);
    }
    fclose(pFile);
    return out;
}


//DEM coregistration with 2D Image Coregistration
double** DEM_ImageCoregistration(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt, D2DPOINT *adjust_std )
{
    ProInfo *proinfo = (ProInfo*)malloc(sizeof(ProInfo));
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    
    double **ImageAdjust_coreg = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
    
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    FILE *time_fid;
    
    total_ST = time(0);
    
    
    if(OpenProject(_filename,proinfo,args))
    {
        char temp_filepath[500];
        
        int check_folder = 1;
        
        sprintf(proinfo->save_filepath,"%s",args.Outputpath);
        
        int status;
        status = mkdir(proinfo->save_filepath,0777);
        if (opendir(proinfo->save_filepath) == NULL)
        {
            if (status == -1)
            {
                printf("Outpath of '%s' cannot make, please check outpath!!\n",proinfo->save_filepath);
                exit(1);
            }
        }
        
        sprintf(temp_filepath,"%s/txt",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(temp_filepath,"%s/tif",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(proinfo->tmpdir,"%s/tmp",proinfo->save_filepath);
        mkdir(proinfo->tmpdir,0777);
        
        CSize tmp_datasize;
        long int cols[2], rows[2];
        double Sum_grid = 0;
        NCCflag ncc_flag;
        int py_level = proinfo->pyramid_level;
        char save_file[500];
        char out_file[500];
        
        uint8 *image_bits = (uint8*)malloc(sizeof(uint8)*proinfo->number_of_images);
        double *ortho_minX = (double*)malloc(sizeof(double)*proinfo->number_of_images);
        double *ortho_maxY = (double*)malloc(sizeof(double)*proinfo->number_of_images);
        double *ortho_dx = (double*)malloc(sizeof(double)*proinfo->number_of_images);
        double *ortho_dy = (double*)malloc(sizeof(double)*proinfo->number_of_images);
        double *GridSize_width = (double*)calloc(sizeof(double),proinfo->number_of_images);
        double *GridSize_height = (double*)calloc(sizeof(double),proinfo->number_of_images);
        CSize *OriImagesizes = (CSize*)malloc(sizeof(CSize)*proinfo->number_of_images);
        double *avg_roh = (double*)calloc(sizeof(double),proinfo->number_of_images);
        int *Grid_space = (int*)calloc(sizeof(int),py_level+1);
        int *RA_iter_counts = (int*)calloc(sizeof(int),proinfo->number_of_images);
        D2DPOINT**MPs_2D = NULL;
        
        uint16 **OriImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
        double **ImageBoundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
        double **Boundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
        //char **Subsetfilename;
        CSize **data_size_lr = (CSize**)malloc(sizeof(CSize*)*proinfo->number_of_images);
        uint16 ***SubImages = (uint16***)malloc(sizeof(uint16**)*(py_level+1));
        
        ncc_flag.rotate_flag     = 1;       ncc_flag.multi_flag      = 0;       ncc_flag.multi_flag_sum  = 0;       ncc_flag.inter_flag      = 1;
        ncc_flag.weight_flag     = 0;
        
        //Subsetfilename = (char**)malloc(sizeof(char*)*proinfo->number_of_images);
        
        //read image and info
        for(int i=0;i<proinfo->number_of_images;i++)
        {
            //2d image
            Boundary[i] = (double*)calloc(sizeof(double),4);
            
            proinfo->check_selected_image[i] = true;
            
            ImageAdjust_coreg[i] = (double*)malloc(sizeof(double)*2);
            ImageAdjust_coreg[i][0] = 0.0;
            ImageAdjust_coreg[i][1] = 0.0;
            
            image_bits[i] = ReadGeotiff_bits(args.Image[i]);
            OriImagesizes[i] = ReadGeotiff_info_dxy(args.Image[i],&ortho_minX[i],&ortho_maxY[i],&ortho_dx[i],&ortho_dy[i]);
            /*if(i == 1)
             {
             ortho_minX[i] += 20;
             ortho_maxY[i] += 20;
             }
             */
            ImageBoundary[i] = (double*)malloc(sizeof(double)*4);
            ImageBoundary[i][0] = ortho_minX[i];
            ImageBoundary[i][1] = ortho_maxY[i] - ortho_dy[i]*OriImagesizes[i].height;
            ImageBoundary[i][2] = ortho_minX[i] + ortho_dx[i]*OriImagesizes[i].width;
            ImageBoundary[i][3] = ortho_maxY[i];
            
            
            if(i == 0)
            {
                Sum_grid = ((ortho_dx[i] + ortho_dy[i])/2.0);
                
                Boundary[i][0] = ImageBoundary[i][0];
                Boundary[i][1] = ImageBoundary[i][1];
                Boundary[i][2] = ImageBoundary[i][2];
                Boundary[i][3] = ImageBoundary[i][3];
            }
            else
            {
                if(Boundary[0][0] < ImageBoundary[i][0])
                    Boundary[i][0] = ImageBoundary[i][0];
                else
                    Boundary[i][0] = Boundary[0][0];
                
                if(Boundary[0][1] < ImageBoundary[i][1])
                    Boundary[i][1] = ImageBoundary[i][1];
                else
                    Boundary[i][1] = Boundary[0][1];
                
                if(Boundary[0][2] > ImageBoundary[i][2])
                    Boundary[i][2] = ImageBoundary[i][2];
                else
                    Boundary[i][2] = Boundary[0][2];
                
                if(Boundary[0][3] > ImageBoundary[i][3])
                    Boundary[i][3] = ImageBoundary[i][3];
                else
                    Boundary[i][3] = Boundary[0][3];
            }
            
            if(args.check_boundary)
            {
                if(Boundary[i][0] < args.Min_X)
                    Boundary[i][0] = args.Min_X;
                if(Boundary[i][1] < args.Min_Y)
                    Boundary[i][1] = args.Min_Y;
                if(Boundary[i][2] > args.Max_X)
                    Boundary[i][2] = args.Max_X;
                if(Boundary[i][3] > args.Max_Y)
                    Boundary[i][3] = args.Max_Y;
                
                Boundary[i][0] = args.Min_X;
                Boundary[i][1] = args.Min_Y;
                Boundary[i][2] = args.Max_X;
                Boundary[i][3] = args.Max_Y;
            }
            
            GridSize_width[i] = Boundary[i][2] - Boundary[i][0];
            GridSize_height[i] = Boundary[i][3] - Boundary[i][1];
            
            //printf("ID %d\tboundary = %f\t%f\t%f\t%f\t%f\t%f\n",i,Boundary[i][0],Boundary[i][1],Boundary[i][2],Boundary[i][3],GridSize_width[i],GridSize_height[i]);
            
            //printf("ID %d\t %s\t%d\t%d\t%f\t%f\t%f\t%f\n",i,args.Image[i],ortho_size[i].width,ortho_size[i].height,ortho_minX[i],ortho_maxY[i],ortho_dx[i],ortho_dy[i]);
            cols[0] = 0;
            cols[1] = OriImagesizes[i].width;
            
            rows[0] = 0;
            rows[1] = OriImagesizes[i].height;
            
            switch(image_bits[i])
            {
                case 8:
                {
                    uint8 type;
                    uint8 *data8 = Readtiff_T(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i], type);
                    long int data_size8 = (long int)OriImagesizes[i].width*(long int)OriImagesizes[i].height;
                    OriImages[i] = (uint16*)malloc(sizeof(uint16)*data_size8);
#pragma omp parallel for schedule(guided)
                    for(long int index = 0 ; index < data_size8 ; index++)
                        OriImages[i][index] = data8[index];
                    free(data8);
                }
                    break;
                case 12:
                {
                    uint16 type;
                    uint16* data16 = Readtiff_T(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],type);
                    long int data_size16 = (long int)OriImagesizes[i].width*(long int)OriImagesizes[i].height;
                    OriImages[i] = (uint16*)malloc(sizeof(uint16)*data_size16);
#pragma omp parallel for schedule(guided)
                    for(long int index = 0 ; index < data_size16 ; index++)
                        OriImages[i][index] = data16[index];
                    free(data16);
                }
                    
                    
                    //OriImages[i] = Readtiff(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],false);
                    break;
            }
            
            //OriImages[i] = Readtiff(args.Image[i], &OriImagesizes[i], cols, rows, &tmp_datasize, false);
            
            //printf("ID %d\t %s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,args.Image[i],OriImagesizes[i].width,OriImagesizes[i].height,ortho_minX[i],ortho_maxY[i],ortho_dx[i],ortho_dy[i],ImageBoundary[i][0],ImageBoundary[i][1],ImageBoundary[i][2],ImageBoundary[i][3]);
            
            /*
            char *filename;
            sprintf(proinfo->Imagefilename[i],"%s",args.Image[i]);
            //printf("image %d\t%s\n",ti,proinfo->Imagefilename[ti]);
            Subsetfilename[i] = (char*)malloc(sizeof(char)*500);
            filename = GetFileName(args.Image[i]);
            filename = remove_ext(filename);
            sprintf(Subsetfilename[i],"%s/%s_subset.raw",proinfo->tmpdir,filename);
            //printf("subsetfilename %s\n",Subsetfilename[ti]);
            */
            data_size_lr[i] = (CSize*)malloc(sizeof(CSize)*(py_level+1));
            SetPySizes(data_size_lr[i], OriImagesizes[i], py_level);
            //free(filename);
            //for (int ttt = 0 ; ttt < py_level+1 ;ttt++)
            //    printf("data_size %d\t%d\n",data_size_lr[i][ttt].width,data_size_lr[i][ttt].height);
        }
        
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nbefore preprocessing time %f\n\n",total_gap);
        total_ST = time(0);
        
        //Preprocessing_Coreg(proinfo,proinfo->tmpdir,OriImages,Subsetfilename,py_level,OriImagesizes,data_size_lr);
        
        printf("start Preprocessing\n");
        
        for(int level = 0 ; level <= py_level ; level++)
        {
            SubImages[level] = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
            
            for(uint8 ti = 0; ti<proinfo->number_of_images ; ti++)
            {
                CSize data_size;
                long int data_length;
 
                if(level == 0)
                {
                    data_size = data_size_lr[ti][level];
                    
                    data_length = (long int)data_size.height*(long int)data_size.width;
                    SubImages[level][ti] = (uint16*)malloc(sizeof(uint16)*data_length);
                    memcpy(SubImages[level][ti],OriImages[ti],sizeof(uint16)*(long)data_length);
                }
                else
                {
                    data_size = data_size_lr[ti][level-1];
                    
                    data_length = (long int)data_size.height*(long int)data_size.width;
                    SubImages[level][ti] = (uint16*)malloc(sizeof(uint16)*data_length);
                    
                    SubImages[level][ti] = CreateImagePyramid(SubImages[level-1][ti],data_size,5,(double)1.6);
                }
                //printf("preprocessing size level ID %d\t%d\t%d\t%d\n",level,ti,data_size.width,data_size.height);
            }
        }
        
        
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\npreprocessing time %f\n\n",total_gap);
        total_ST = time(0);
        
        
        
        //printf("overlapped area %f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3],GridSize_width,GridSize_height,Grid_space,Grid_size.width,Grid_size.height,total_grid_counts);
        for(int level = py_level ; level >= 0 ; level --)
        {
            Grid_space[level] = ceil(Sum_grid)*pwrtwo(py_level+2);
            if(Grid_space[level] < 50)
                Grid_space[level] = 50;
            //printf("image coregistration gridspace %d\t%d\n",level,Grid_space[level]);
            
            //printf("Processing level %d\n",level);
            printf("level\tImage ID\trow(pixel)\tcolumn(pixel)\tTy(meter)\tTx(meter)\tGCPS #\tavg_roh\t# of iteration\n");
            
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\niter %d LoadPyramidImages %f\n\n",level,total_gap);
            total_ST = time(0);
            */
            int iter_counts;
            MPs_2D = CoregParam_Image_MPs(proinfo, level,py_level, ImageAdjust_coreg, ncc_flag,
                                              5, SubImages[level], data_size_lr, ImageBoundary, ortho_dx, ortho_dy,
                                              Grid_space,Boundary,proinfo->save_filepath,avg_roh,&iter_counts,adjust_std,RA_iter_counts);
            
            for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
            {
                //printf("2 std %f\t%f\n",adjust_std[ti].m_X,adjust_std[ti].m_Y);
                printf("%d\t%d\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%4.2f\t%d\n",level,ti,ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1],
                       -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],RA_iter_counts[ti],avg_roh[ti],iter_counts);
                free(SubImages[level][ti]);
            }
            free(SubImages[level]);
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\niter %d CoregParam_Image %f\n\n",level,total_gap);
            total_ST = time(0);
             */
        }
        
        
        
        FILE* fid_out = NULL;
        sprintf(out_file,"%s/coreg_result.txt",proinfo->save_filepath);
        fid_out         = fopen(out_file,"w");
        fprintf(fid_out,"orthoimage name\tline(row) direction[pixel]\tsample(column) direction[pixel]\tTy[meter]\tTx[meter]\tavg_roh\n");
        for(int ti = 0 ; ti < proinfo->number_of_images ; ti ++)
        {
            fprintf(fid_out,"%s\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%3.2f\n",proinfo->Imagefilename[ti],ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1],
                    -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],avg_roh[ti]);
            
            free(OriImages[ti]);
            free(ImageBoundary[ti]);
            free(Boundary[ti]);
            //free(Subsetfilename[ti]);
            free(data_size_lr[ti]);
        }
        fclose(fid_out);
        /*
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nCoregParam_Image %f\n\n",total_gap);
        total_ST = time(0);
        //RemoveFiles(proinfo,proinfo->tmpdir,Subsetfilename,py_level,0);
        */
        free(SubImages);
        free(OriImages);
        free(ImageBoundary);
        free(Boundary);
        //free(Subsetfilename);
        free(data_size_lr);
        
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nCoregParam_Image time %f\n\n",total_gap);
        total_ST = time(0);
        
        
        //Tz and stat
        TransParam param;
        SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
        
        float *DEM = NULL;
        //dem image = ref
        char DEM_name_refoutfile[500];
        char ref_dem_name[500];
        char *ref_no_ext;
        ref_no_ext = remove_ext(args.Image[0]);
        int char_size = strlen(ref_no_ext);
        char* ref_fchar = (char*)malloc(sizeof(char)*(char_size-4));
        for(int c = 0 ; c < char_size - 5 ; c++)
        {
            ref_fchar[c] = ref_no_ext[c];
        }
        ref_fchar[char_size-5] = '\0';
        sprintf(ref_dem_name,"%sdem.tif",ref_fchar);
        printf("dem name = %s\n",ref_dem_name);
        free(ref_no_ext);
        free(ref_fchar);
        
        cols[0] = 0;
        cols[1] = OriImagesizes[0].width;
        
        rows[0] = 0;
        rows[1] = OriImagesizes[0].height;
        
        float type;
        DEM = Readtiff_T(ref_dem_name,&OriImagesizes[0],cols,rows,&OriImagesizes[0],type);
        
        char *reffilename  = SetOutpathName(args.Image[0]);
        ref_no_ext = remove_ext(reffilename);
        char_size = strlen(ref_no_ext);
        ref_fchar = (char*)malloc(sizeof(char)*(char_size-4));
        for(int c = 0 ; c < char_size - 5 ; c++)
        {
            ref_fchar[c] = ref_no_ext[c];
        }
        ref_fchar[char_size-5] = '\0';
        sprintf(DEM_name_refoutfile,"%sdem.tif",ref_fchar);
        free(ref_no_ext);
        free(ref_fchar);
        
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nRef DEM loading time %f\n\n",total_gap);
        total_ST = time(0);
        
        
        sprintf(out_file,"%s/DEM_coreg_result.txt",args.Outputpath);
        fid_out         = fopen(out_file,"w");
        fprintf(fid_out,"ref DEM name\t%s\n",DEM_name_refoutfile);
        fprintf(fid_out,"DEM name\t\t\t\t\t\t\tTx[meter]\tTy[meter]\tTz(average)\tTz(median)\tTx_std[meter]\tTy_std[meter]\tTz_std(average)\tTz_std(median)\tMean_all(avg)\tMedian_all(avg)\tdz_std(avg)\tMean_all(med)\tMedian_all(med)\tdz_std(med)\tNumberOfCPs\n");
        
        
        for(int ti = 1 ; ti < proinfo->number_of_images ; ti ++)
        {
            time_t total_ST_a = 0, total_ET_a = 0;
            total_ST_a = time(0);
            
            //dem image
            float *DEM_tar = NULL;
            char DEM_name_outfile[500];
            
            char ref_dem_name[500];
            char *ref_no_ext;
            ref_no_ext = remove_ext(args.Image[ti]);
            int char_size = strlen(ref_no_ext);
            char* ref_fchar = (char*)malloc(sizeof(char)*(char_size-4));
            for(int c = 0 ; c < char_size - 5 ; c++)
            {
                ref_fchar[c] = ref_no_ext[c];
            }
            ref_fchar[char_size-5] = '\0';
            sprintf(ref_dem_name,"%sdem.tif",ref_fchar);
            printf("dem name = %s\n",ref_dem_name);
            free(ref_no_ext);
            free(ref_fchar);
            
            cols[0] = 0;
            cols[1] = OriImagesizes[ti].width;
            
            rows[0] = 0;
            rows[1] = OriImagesizes[ti].height;
            
            float type;
            DEM_tar= Readtiff_T(ref_dem_name,&OriImagesizes[ti],cols,rows,&OriImagesizes[ti],type);
            
            char *reffilename  = SetOutpathName(args.Image[ti]);
            ref_no_ext = remove_ext(reffilename);
            char_size = strlen(ref_no_ext);
            ref_fchar = (char*)malloc(sizeof(char)*(char_size-4));
            for(int c = 0 ; c < char_size - 5 ; c++)
            {
                ref_fchar[c] = ref_no_ext[c];
            }
            ref_fchar[char_size-5] = '\0';
            sprintf(DEM_name_outfile,"%sdem.tif",ref_fchar);
            free(ref_no_ext);
            free(ref_fchar);
        
        
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nDEM loading time ID %d\t%f\n\n",ti,total_gap);
            total_ST = time(0);
        
            char co_dems_name[500];
            char copoly_dems_name[500];
            char co_dems_name_diff[500];
            char copoly_dems_name_diff[500];
            
            sprintf(co_dems_name,"%s/%s_coreg_average.tif",args.Outputpath,DEM_name_outfile);
            sprintf(copoly_dems_name,"%s/%s_coreg_median.tif",args.Outputpath,DEM_name_outfile);
            sprintf(co_dems_name_diff,"%s/%s_coreg_average_diff.tif",args.Outputpath,DEM_name_outfile);
            sprintf(copoly_dems_name_diff,"%s/%s_coreg_median_diff.tif",args.Outputpath,DEM_name_outfile);
            
            
            double Coreg_param[2];
            double ref_grid_size, tar_grid_size;
            F3DPOINT *save_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),RA_iter_counts[ti]);
            CSize ref_data_size = OriImagesizes[0];
            CSize tar_data_size = OriImagesizes[ti];
            double ref_minX = ortho_minX[0];
            double ref_maxY = ortho_maxY[0];
            double tar_minX = ortho_minX[ti];
            double tar_maxY = ortho_maxY[ti];
            
            Coreg_param[1] = ImageAdjust_coreg[ti][1]*ortho_dx[ti];
            Coreg_param[0] = - ImageAdjust_coreg[ti][0]*ortho_dx[ti];
            ref_grid_size = ortho_dx[0];
            tar_grid_size = ortho_dx[ti];
            
            //printf("image ID %d\tcoreg y=%f\tx=%f\t%d\n",ti,Coreg_param[0],Coreg_param[1],RA_iter_counts[ti]);
            //printf("ref size %d\t%d\ttar size %d\t%d\n",ref_data_size.width,ref_data_size.height,tar_data_size.width,tar_data_size.height);
            //printf("ref coord %f\t%f\ttar coord %f\t%f\n",ref_minX,ref_maxY,tar_minX,tar_maxY);
            int diff_count = 0;
            double sum_diff = 0;
            double W_th = 95;
            double dH_th = 30;
            
#pragma omp parallel for reduction(+:sum_diff,diff_count) schedule(guided)
            for(long int gcp_index = 0 ; gcp_index < RA_iter_counts[ti] ; gcp_index ++)
            {
                D2DPOINT gcp_coord,ref_img,tar_img;
                gcp_coord.m_X = MPs_2D[ti][gcp_index].m_X;
                gcp_coord.m_Y = MPs_2D[ti][gcp_index].m_Y;
                
                //printf("ID %d\tMPs %f\t%f\n",gcp_index,gcp_coord.m_X,gcp_coord.m_Y);
                
                ref_img.m_X = (gcp_coord.m_X - ortho_minX[0])/ref_grid_size;
                ref_img.m_Y = (ortho_maxY[0] - gcp_coord.m_Y)/ref_grid_size;
                
                tar_img.m_X = ( (gcp_coord.m_X + Coreg_param[1]) - ortho_minX[ti] )/tar_grid_size;
                tar_img.m_Y = ( ortho_maxY[ti] - (gcp_coord.m_Y + Coreg_param[0]) )/tar_grid_size;
                
                bool check_b = true;
                if(args.check_boundary)
                {
                    if(gcp_coord.m_X >= args.Min_X && gcp_coord.m_X <= args.Max_X && gcp_coord.m_Y >= args.Min_Y && gcp_coord.m_Y <= args.Max_Y)
                        check_b = true;
                    else
                        check_b = false;
                }
                
                if(ref_img.m_X - 2 >= 0 && ref_img.m_X + 2 < ref_data_size.width && ref_img.m_Y - 2 >= 0 && ref_img.m_Y + 2 < ref_data_size.height &&
                   tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_data_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_data_size.height && check_b)
                {
                    long ref_index = (int)(ref_img.m_Y)*ref_data_size.width + (int)(ref_img.m_X);
                    long tar_index = (int)(tar_img.m_Y)*tar_data_size.width + (int)(tar_img.m_X);
                    
                    if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                    {
                        double sum_zr = 0;
                        double sum_zt = 0;
                        int kernel_count = 0;
                        for(int kr = -2 ; kr <= 2 ; kr++ )
                        {
                            for(int kc = -2 ; kc <= 2 ; kc++ )
                            {
                                int k_refindex = (int)(ref_img.m_Y + kr)*ref_data_size.width + (int)(ref_img.m_X + kc);
                                int k_tarindex = (int)(tar_img.m_Y + kr)*tar_data_size.width + (int)(tar_img.m_X + kc);
                                
                                if(DEM[k_refindex] > -100 && DEM[k_refindex] < 10000 && DEM_tar[k_tarindex] > -100 && DEM_tar[k_tarindex] < 10000)
                                {
                                    kernel_count++;
                                    //printf("DEM value %f\t%f\n",DEM[0][k_refindex],DEM[ti][k_tarindex]);
                                }
                            }
                        }
                        
                        double* ref_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                        double* tar_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                        
                        GMA_double *A_matrix = GMA_double_create(kernel_count, 3);
                        GMA_double *L_matrix = GMA_double_create(kernel_count, 1);
                        GMA_double *AT_matrix = GMA_double_create(3,kernel_count);
                        GMA_double *ATA_matrix = GMA_double_create(3,3);
                        
                        GMA_double *ATAI_matrix = GMA_double_create(3,3);
                        GMA_double *ATL_matrix = GMA_double_create(3,1);
                        
                        GMA_double *X_matrix = GMA_double_create(3,1);
                        GMA_double *AX_matrix = GMA_double_create(kernel_count,1);
                        GMA_double *V_matrix = GMA_double_create(kernel_count,1);
                        
                        GMA_double *tA_matrix = GMA_double_create(kernel_count, 3);
                        GMA_double *tL_matrix = GMA_double_create(kernel_count, 1);
                        GMA_double *tAT_matrix = GMA_double_create(3,kernel_count);
                        GMA_double *tATA_matrix = GMA_double_create(3,3);
                        
                        GMA_double *tATAI_matrix = GMA_double_create(3,3);
                        GMA_double *tATL_matrix = GMA_double_create(3,1);
                        
                        GMA_double *tX_matrix = GMA_double_create(3,1);
                        GMA_double *tAX_matrix = GMA_double_create(kernel_count,1);
                        GMA_double *tV_matrix = GMA_double_create(kernel_count,1);
                        
                        kernel_count = 0;
                        
                        sum_zr = 0;
                        sum_zt = 0;
                        
                        for(int kr = -2 ; kr <= 2 ; kr++ )
                        {
                            for(int kc = -2 ; kc <= 2 ; kc++ )
                            {
                                int k_refindex = (int)(ref_img.m_Y + kr)*ref_data_size.width + (int)(ref_img.m_X + kc);
                                int k_tarindex = (int)(tar_img.m_Y + kr)*tar_data_size.width + (int)(tar_img.m_X + kc);
                                if(DEM[k_refindex] > -100 && DEM[k_refindex] < 10000 && DEM_tar[k_tarindex] > -100 && DEM_tar[k_tarindex] < 10000)
                                {
                                    ref_patch_vecs[kernel_count] = DEM[k_refindex];
                                    tar_patch_vecs[kernel_count] = DEM_tar[k_tarindex];
                                    
                                    A_matrix->val[kernel_count][0] = kc*ref_grid_size;
                                    A_matrix->val[kernel_count][1] = kr*ref_grid_size;
                                    A_matrix->val[kernel_count][2] = 1.0;
                                    L_matrix->val[kernel_count][0] = DEM[k_refindex];
                                    
                                    tA_matrix->val[kernel_count][0] = kc*tar_grid_size;
                                    tA_matrix->val[kernel_count][1] = kr*tar_grid_size;
                                    tA_matrix->val[kernel_count][2] = 1.0;
                                    tL_matrix->val[kernel_count][0] = DEM_tar[k_tarindex];
                                    
                                    kernel_count++;
                                }
                            }
                        }
                        //double var_diff = fabs(sqrt(sum_zr/(double)kernel_count) - sqrt(sum_zt/(double)kernel_count));
                        
                        double k_ncc = Correlate(ref_patch_vecs, tar_patch_vecs, kernel_count);
                        if (k_ncc != -99)
                            k_ncc = (k_ncc + 1)/2.0;
                        else
                            k_ncc = 0.0;
                        
                        GMA_double_Tran(A_matrix,AT_matrix);
                        GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
                        GMA_double_inv(ATA_matrix,ATAI_matrix);
                        GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
                        GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
                        GMA_double_mul(A_matrix,X_matrix,AX_matrix);
                        GMA_double_sub(AX_matrix,L_matrix,V_matrix);
                        
                        GMA_double_Tran(tA_matrix,tAT_matrix);
                        GMA_double_mul(tAT_matrix,tA_matrix,tATA_matrix);
                        GMA_double_inv(tATA_matrix,tATAI_matrix);
                        GMA_double_mul(tAT_matrix,tL_matrix,tATL_matrix);
                        GMA_double_mul(tATAI_matrix,tATL_matrix,tX_matrix);
                        GMA_double_mul(tA_matrix,tX_matrix,tAX_matrix);
                        GMA_double_sub(tAX_matrix,tL_matrix,tV_matrix);
                        
                        double N1 = X_matrix->val[0][0];
                        double N2 = X_matrix->val[1][0];
                        double N3 = 1.0;
                        
                        double norm  = sqrt(N1*N1 + N2*N2 + N3*N3);
                        double angle = acos(fabs(N3)/norm)*180/3.141592;
                        
                        if(angle <= 0 && angle >= -90)
                            angle = fabs(angle);
                        else if(angle <= -270 && angle >= -360)
                            angle = 360 + angle;
                        else if(angle >= 270 && angle <= 360)
                            angle = 360 - angle;
                        
                        double aspect = 90 - atan2(N2,N1)*180/3.141592;
                        
                        N1 = tX_matrix->val[0][0];
                        N2 = tX_matrix->val[1][0];
                        N3 = 1.0;
                        
                        norm  = sqrt(N1*N1 + N2*N2 + N3*N3);
                        double angle_tar = acos(fabs(N3)/norm)*180/3.141592;
                        
                        if(angle_tar <= 0 && angle_tar >= -90)
                            angle_tar = fabs(angle_tar);
                        else if(angle_tar <= -270 && angle_tar >= -360)
                            angle_tar = 360 + angle_tar;
                        else if(angle_tar >= 270 && angle_tar <= 360)
                            angle_tar = 360 - angle_tar;
                        
                        double aspect_tar = 90 - atan2(N2,N1)*180/3.141592;
                        
                        double SS = 1 - fabs(angle - angle_tar)/90;
                        double SA = 1 - fabs(aspect - aspect_tar)/360;
                        
                        double W = 40*k_ncc + 40*SS + 20*SA;
                        
                        if(W > W_th && fabs(DEM[ref_index] - DEM_tar[tar_index]) < dH_th)
                        {
                            long t_col_int   = (long)(tar_img.m_X + 0.01);
                            long t_row_int   = (long)(tar_img.m_Y + 0.01);
                            
                            double dcol        = tar_img.m_X - t_col_int;
                            double drow        = tar_img.m_Y - t_row_int;
                            
                            long index1,index2,index3, index4;
                            double value1, value2, value3, value4, value;
                            
                            index1  = (t_col_int   ) + (t_row_int   )*(long)tar_data_size.width;
                            index2  = (t_col_int +1) + (t_row_int   )*(long)tar_data_size.width;
                            index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_data_size.width;
                            index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_data_size.width;
                            
                            value1      = DEM_tar[index1];
                            value2      = DEM_tar[index2];
                            value3      = DEM_tar[index3];
                            value4      = DEM_tar[index4];
                            
                            if(value1 > -100 && value1 < 10000 && value2 > -100 && value2 < 10000 && value3 > -100 && value3 < 10000 && value4 > -100 && value4 < 10000)
                            {
                                value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                + value3*(1-dcol)*drow + value4*dcol*drow;
                                
                                if(fabs(DEM[ref_index] - value) < 30)
                                {
                                    sum_diff += (DEM[ref_index] - value);
                                    diff_count++;
                                    
                                    save_pts[gcp_index].m_X = gcp_coord.m_X + Coreg_param[1];
                                    save_pts[gcp_index].m_Y = gcp_coord.m_Y + Coreg_param[0];
                                    save_pts[gcp_index].m_Z = (DEM[ref_index] - value);
                                    save_pts[gcp_index].flag = 1;
                                    //printf("W kncc ss sa %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",W,k_ncc,SS,SA,angle,angle_tar,aspect,aspect_tar,save_pts[gcp_index].m_Z);
                                }
                            }
                        }
                        
                        //if(W > 90 && fabs(ref_dem[ref_index] - tar_dem[tar_index]) >= 30)
                        //    printf("W dz %f\t%f\t%f\t%f\t%f\n",W,k_ncc,SS,SA,ref_dem[ref_index] - tar_dem[tar_index]);
                        
                        GMA_double_destroy(A_matrix);
                        GMA_double_destroy(L_matrix);
                        GMA_double_destroy(AT_matrix);
                        GMA_double_destroy(ATA_matrix);
                        GMA_double_destroy(ATAI_matrix);
                        GMA_double_destroy(ATL_matrix);
                        GMA_double_destroy(X_matrix);
                        GMA_double_destroy(AX_matrix);
                        GMA_double_destroy(V_matrix);
                        
                        GMA_double_destroy(tA_matrix);
                        GMA_double_destroy(tL_matrix);
                        GMA_double_destroy(tAT_matrix);
                        GMA_double_destroy(tATA_matrix);
                        GMA_double_destroy(tATAI_matrix);
                        GMA_double_destroy(tATL_matrix);
                        GMA_double_destroy(tX_matrix);
                        GMA_double_destroy(tAX_matrix);
                        GMA_double_destroy(tV_matrix);
                        
                        free(ref_patch_vecs);
                        free(tar_patch_vecs);
                    }
                }
            }
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nselect vertical CPs time %f\n\n",total_gap);
            total_ST = time(0);
            
            double average = sum_diff/diff_count;
            
            float* save_Z = (float*)calloc(sizeof(float),diff_count);
            float sum_var = 0;
            float min_save_Z = 100000;
            float max_save_Z = -100000;
            long count = 0;
            
            char dem_gcp_filename[500];
            sprintf(dem_gcp_filename,"%s/DEM_gcps_%d.txt",args.Outputpath,ti);
            FILE *fid_dem_gcp = fopen(dem_gcp_filename,"w");
            
            for(long row = 0; row < RA_iter_counts[ti] ; row++)
            {
                if(save_pts[row].flag == 1)
                {
                    save_Z[count] = save_pts[row].m_Z;
                    
                    sum_var += (average - save_pts[row].m_Z)*(average - save_pts[row].m_Z);
                    fprintf(fid_dem_gcp,"%f\t%f\t%f\n",save_pts[row].m_X,save_pts[row].m_Y,save_pts[row].m_Z);
                    count++;
                }
            }
            fclose(fid_dem_gcp);
            
            float MED_z = quickselect(save_Z, diff_count, (int)(diff_count/2.0));
            //float MED_z = median(diff_count,save_Z,-(dH_th+5),dH_th+5);
            
            double SD_z = sqrt(sum_var/diff_count);
            
            double sum_var_med = 0;
            
#pragma omp parallel for reduction(+:sum_var_med) schedule(guided)
            for(int i=0;i<diff_count;i++)
            {
                sum_var_med += (MED_z -save_Z[i])*(MED_z -save_Z[i]);
            }
            
            double SD_z_med = sqrt(sum_var_med/diff_count);
            
            free(save_Z);
            free(save_pts);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nTz average and variation cal time %f\t%d\t%d\t%f\t%f\t%f\t%f\n\n",total_gap,count,diff_count,average,MED_z,SD_z,SD_z_med);
            total_ST = time(0);
            
            
            float* co_dem = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            float* copoly_dem = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            /*
            float* co_dem_diff = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            float* copoly_dem_diff = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            */
            float* save_dz = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            float* save_dz_med = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
            
            float dz_sum = 0;
            float dz_sum_med = 0;
            long dz_count = 0;
            min_save_Z = 100000;
            max_save_Z = -100000;
            for(long int co_index = 0 ; co_index < (long)(tar_data_size.width)*(long)(tar_data_size.height) ; co_index++)
            {
                int pts_row = (int)(floor(co_index/tar_data_size.width));
                int pts_col = co_index % tar_data_size.width;
                
                D2DPOINT gcp_coord,tar_img, ref_img,gcp_coord_ref;
                gcp_coord.m_X = pts_col*tar_grid_size + tar_minX + Coreg_param[1];
                gcp_coord.m_Y = tar_maxY - pts_row*tar_grid_size + Coreg_param[0];
                
                tar_img.m_X = ( gcp_coord.m_X - tar_minX )/tar_grid_size;
                tar_img.m_Y = ( tar_maxY - gcp_coord.m_Y )/tar_grid_size;
                
                gcp_coord_ref.m_X = pts_col*tar_grid_size + tar_minX;
                gcp_coord_ref.m_Y = tar_maxY - pts_row*tar_grid_size;
                
                ref_img.m_X = ( gcp_coord_ref.m_X - ref_minX )/ref_grid_size;
                ref_img.m_Y = ( ref_maxY - gcp_coord_ref.m_Y )/ref_grid_size;
                
                if(tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_data_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_data_size.height &&
                   ref_img.m_X >= 0 && ref_img.m_X < ref_data_size.width && ref_img.m_Y >= 0 && ref_img.m_Y < ref_data_size.height)
                {
                    long tar_index = (int)(tar_img.m_Y)*tar_data_size.width + (int)(tar_img.m_X);
                    long ref_index = (int)(ref_img.m_Y)*ref_data_size.width + (int)(ref_img.m_X);
                    
                    if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                    {
                        long t_col_int   = (long)(tar_img.m_X + 0.01);
                        long t_row_int   = (long)(tar_img.m_Y + 0.01);
                        
                        double dcol        = tar_img.m_X - t_col_int;
                        double drow        = tar_img.m_Y - t_row_int;
                        
                        long index1,index2,index3, index4;
                        double value1, value2, value3, value4, value;
                        
                        index1  = (t_col_int   ) + (t_row_int   )*(long)tar_data_size.width;
                        index2  = (t_col_int +1) + (t_row_int   )*(long)tar_data_size.width;
                        index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_data_size.width;
                        index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_data_size.width;
                        
                        value1      = DEM_tar[index1];
                        value2      = DEM_tar[index2];
                        value3      = DEM_tar[index3];
                        value4      = DEM_tar[index4];
                        
                        value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                        + value3*(1-dcol)*drow + value4*dcol*drow;
                        
                        //double value       = DEM[ti][tar_index];
                        co_dem[tar_index] = value + average;
                        
                        float copoly_dem_t = value + MED_z;
                        copoly_dem[tar_index] = copoly_dem_t;
                        
                        float co_dem_diff = (DEM[ref_index] - co_dem[tar_index]);
                        float copoly_dem_diff = (DEM[ref_index] - copoly_dem_t);
                        
                        //co_dem_diff[tar_index] = (DEM[0][ref_index] - co_dem[tar_index]);
                        //copoly_dem_diff[tar_index] = (DEM[0][ref_index] - copoly_dem[tar_index]);
                        
                        if(fabs(co_dem_diff) < dH_th && fabs(copoly_dem_diff) < dH_th)
                        {
                            dz_sum += co_dem_diff;
                            dz_sum_med += copoly_dem_diff;
                            
                            save_dz[dz_count] = co_dem_diff;
                            save_dz_med[dz_count] = copoly_dem_diff;
                            /*
                            if(save_dz[dz_count] < min_save_Z)
                                min_save_Z = save_dz[dz_count];
                            if(save_dz[dz_count] > max_save_Z)
                                max_save_Z = save_dz[dz_count];
                            
                            if(save_dz_med[dz_count] < min_save_Z)
                                min_save_Z = save_dz_med[dz_count];
                            if(save_dz_med[dz_count] > max_save_Z)
                                max_save_Z = save_dz_med[dz_count];
                            */
                            dz_count++;
                        }
                    }
                }
            }
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nWhole image stat %f\n\n",total_gap);
            total_ST = time(0);
            
            float all_average = dz_sum/dz_count;
            float all_med = quickselect(save_dz, dz_count, (int)(dz_count/2.0));
            //float all_med = binmedian(dz_count, save_dz);
            //float all_med = median(dz_count,save_dz,-(dH_th + 5),dH_th + 5);
            float all_sum_var = 0;
            
            float all_average_med = dz_sum_med/dz_count;
            float all_med_med = quickselect(save_dz_med, dz_count, (int)(dz_count/2.0));
            //float all_med_med = binmedian(dz_count,save_dz_med);
            //float all_med_med = median(dz_count,save_dz_med,-(dH_th + 5),dH_th + 5);
            float all_sum_var_med = 0;
            
#pragma omp parallel for reduction(+:all_sum_var,all_sum_var_med) schedule(guided)
            for(long index = 0 ; index < dz_count ; index++)
            {
                all_sum_var += (save_dz[index] - all_average)*(save_dz[index] - all_average);
                all_sum_var_med += (save_dz_med[index] - all_average_med)*(save_dz_med[index] - all_average_med);
            }
            //printf("done allsum\n");
            free(save_dz);
            free(save_dz_med);
            
            double all_std = sqrt(all_sum_var/dz_count);
            double all_std_med = sqrt(all_sum_var_med/dz_count);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nall median median time %f\n\n",total_gap);
            total_ST = time(0);
            
            printf("all stat %f\t%f\t%f\t%f\t%f\t%f\n",all_average,all_med,all_std,all_average_med,all_med_med,all_std_med);
            
            WriteGeotiff(co_dems_name, co_dem, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
            
            WriteGeotiff(copoly_dems_name, copoly_dem, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
            /*
            WriteGeotiff(co_dems_name_diff, co_dem_diff, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
            WriteGeotiff(copoly_dems_name_diff, copoly_dem_diff, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
            */
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\noutput save time %f\n\n",total_gap);
            total_ST = time(0);
            
            free(co_dem);
            free(copoly_dem);
            /*
            free(co_dem_diff);
            free(copoly_dem_diff);
            */
            fprintf(fid_out,"%s\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\n",DEM_name_outfile,Coreg_param[1],Coreg_param[0],average,MED_z,adjust_std[ti].m_Y,adjust_std[ti].m_X,SD_z,SD_z_med,
                all_average,all_med,all_std,all_average_med,all_med_med,all_std_med,diff_count);
           
            free(MPs_2D[ti]);
            free(DEM_tar);
            
            total_ET_a = time(0);
            total_gap = difftime(total_ET_a,total_ST_a);
            printf("\nTotal Tz computation time ID %d\t%f\n\n",ti,total_gap);
        }
        
        fclose(fid_out);
        
        free(image_bits);
        free(ortho_minX);
        free(ortho_maxY);
        free(ortho_dx);
        free(ortho_dy);
        free(GridSize_width);
        free(GridSize_height);
        free(OriImagesizes);
        free(avg_roh);
        free(Grid_space);
        free(RA_iter_counts);
        
        free(MPs_2D[0]);
        
        free(DEM);
        free(MPs_2D);
    }
    
    return ImageAdjust_coreg;
}

void DEM_ImageCoregistration_hillshade(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt)
{
    ProInfo *proinfo = (ProInfo*)malloc(sizeof(ProInfo));
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    FILE *time_fid;
    
    total_ST = time(0);
    int py_kernel_size = 3;
    
    bool check_start = true;
    if(args.check_txt_input == 0)
    {
        check_start = OpenProject(_filename,proinfo,args);
    }
    else if(args.check_txt_input == 1)
    {
        FILE *fid_read = fopen(args.DEM_input_file,"r");
        
        args.number_of_images = 0;
        while( (fscanf(fid_read,"%s\n",args.Image[args.number_of_images])) != EOF )
        {
            printf("input DEM %s\n",args.Image[args.number_of_images]);
            sprintf(proinfo->Imagefilename[args.number_of_images],"%s",args.Image[args.number_of_images]);
            args.number_of_images++;
        }
        fclose(fid_read);
        proinfo->number_of_images = args.number_of_images;
        printf("image num %d\n",args.number_of_images);
    }
    
    if(check_start)
    {
        char temp_filepath[500];
        
        int check_folder = 1;
        
        sprintf(proinfo->save_filepath,"%s",args.Outputpath);
        
        int status;
        status = mkdir(proinfo->save_filepath,0777);
        if (opendir(proinfo->save_filepath) == NULL)
        {
            if (status == -1)
            {
                printf("Outpath of '%s' cannot make, please check outpath!!\n",proinfo->save_filepath);
                exit(1);
            }
        }
        
        sprintf(temp_filepath,"%s/txt",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(temp_filepath,"%s/tif",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(proinfo->tmpdir,"%s/tmp",proinfo->save_filepath);
        mkdir(proinfo->tmpdir,0777);
        
     
        long int cols[2], rows[2];
        double Sum_grid = 0;
        NCCflag ncc_flag;
        int py_level = proinfo->pyramid_level;
        char save_file[500];
        char out_file[500];
        char out_file_bad[500];
        CSize ref_dem_size;
        double ref_minX,ref_maxY,ref_dx,ref_dy;
        
        ncc_flag.rotate_flag     = 1;       ncc_flag.multi_flag      = 0;       ncc_flag.multi_flag_sum  = 0;       ncc_flag.inter_flag      = 1;
        ncc_flag.weight_flag     = 0;
        
        TransParam param;
        SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
        
        //dem image = ref
        float *DEM = NULL;
        unsigned char *ref_hill = NULL;
        CSize *data_size_ref = (CSize*)malloc(sizeof(CSize)*(py_level+1));
        
        
        char DEM_name_refoutfile[500];
        char ref_hill_inputname[500];
        char ref_dem_name[500];
        
        sprintf(ref_dem_name,"%s",args.Image[0]);
        printf("dem name = %s\n",ref_dem_name);
        
        ref_dem_size = ReadGeotiff_info_dxy(ref_dem_name,&ref_minX,&ref_maxY,&ref_dx,&ref_dy);
        
        cols[0] = 0;
        cols[1] = ref_dem_size.width;
        
        rows[0] = 0;
        rows[1] = ref_dem_size.height;
        
        float type;
        DEM = Readtiff_T(ref_dem_name,&ref_dem_size,cols,rows,&ref_dem_size,type);
        
        char *reffilename  = SetOutpathName(args.Image[0]);
        char *ref_no_ext = remove_ext(reffilename);
        int char_size = strlen(ref_no_ext);
        char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
        for(int c = 0 ; c < char_size - 3 ; c++)
        {
            ref_fchar[c] = ref_no_ext[c];
        }
        ref_fchar[char_size-3] = '\0';
        sprintf(DEM_name_refoutfile,"%sdem.tif",ref_fchar);
        free(ref_no_ext);
        free(ref_fchar);
        free(reffilename);
        
        ref_no_ext = remove_ext(args.Image[0]);
        char_size = strlen(ref_no_ext);
        ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
        for(int c = 0 ; c < char_size - 3 ; c++)
        {
            ref_fchar[c] = ref_no_ext[c];
        }
        ref_fchar[char_size-3] = '\0';
        sprintf(ref_hill_inputname,"%sdem_hillshade.tif",ref_fchar);
        printf("ref hill name %s\n",ref_hill_inputname);
        free(ref_no_ext);
        free(ref_fchar);
        /*
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nRef DEM loading time %f\n\n",total_gap);
        total_ST = time(0);
        */
        //create hillshadeimage
        FILE *checktif = NULL;
        checktif = fopen(ref_hill_inputname,"r");
        if(checktif)
        {
            TIFF *refhilltif = NULL;
            refhilltif  = TIFFOpen(ref_hill_inputname,"r");
            uint8 type;
            ref_hill = Readtiff_T(ref_hill_inputname, &ref_dem_size, cols, rows, &ref_dem_size, type);
            TIFFClose(refhilltif);
            fclose(checktif);
            printf("loading existing hillshade image %s\n",ref_hill_inputname);
        }
        else
        {
            ref_hill = CreateHillshade(DEM,ref_dem_size,ref_dx);
        }
        //printf("hillshade %d\t%d\n",ref_dem_size.width,ref_dem_size.height);
        /*
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nhillshade time %f\n\n",total_gap);
        total_ST = time(0);
        */
        //create pyramidimage
        SetPySizes(data_size_ref, ref_dem_size, py_level);
        
        unsigned char **SubImages_ref = (unsigned char**)malloc(sizeof(unsigned char*)*(py_level+1));
        
        for(int level = 0 ; level <= py_level ; level++)
        {
            CSize data_size;
            long int data_length;
            
            if(level == 0)
            {
                data_size = data_size_ref[level];
                
                data_length = (long int)data_size.height*(long int)data_size.width;
                SubImages_ref[level] = (unsigned char*)malloc(sizeof(unsigned char)*data_length);
                memcpy(SubImages_ref[level],ref_hill,sizeof(unsigned char)*(long)data_length);
            }
            else
            {
                data_size = data_size_ref[level-1];
                
                data_length = (long int)data_size.height*(long int)data_size.width;
                //SubImages_ref[level]= CreateImagePyramid_BYTE(SubImages_ref[level-1],data_size,py_kernel_size,(double)1.6);
                SubImages_ref[level]= CreateImagePyramid(SubImages_ref[level-1],data_size,py_kernel_size,(double)1.6);
            }
            //printf("preprocessing size level ID %d\t%d\t%d\n",level,data_size.width,data_size.height);
        }
        /*
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nRef pyramid time %f\n\n",total_gap);
        total_ST = time(0);
        */
        double ImageBoundary_ref[4] = {0.0};
        ImageBoundary_ref[0] = ref_minX;
        ImageBoundary_ref[1] = ref_maxY - ref_dy*ref_dem_size.height;
        ImageBoundary_ref[2] = ref_minX + ref_dx*ref_dem_size.width;
        ImageBoundary_ref[3] = ref_maxY;
        
        sprintf(out_file,"%s/DEM_coreg_result.txt",args.Outputpath);
        FILE* fid_out         = fopen(out_file,"w");
        fprintf(fid_out,"ref DEM name\t%s\n",DEM_name_refoutfile);
        //printf("ref dem %s\n",DEM_name_refoutfile);
        fprintf(fid_out,"DEM name\t\t\t\t\t\t\tTx[meter]\tTy[meter]\tTz(average)\tTz(median)\tTx_std[meter]\tTy_std[meter]\tTz_std(average)\tTz_std(median)\tControls_rho\tMean_all(avg)\tMedian_all(avg)\tdz_std(avg)\tMean_all(med)\tMedian_all(med)\tdz_std(med)\tNumberOfCPs\tprocessing time\n");
        /*
        sprintf(out_file_bad,"%s/DEM_coreg_result_bad.txt",args.Outputpath);
        FILE* fid_out_bad         = fopen(out_file_bad,"w");
        fprintf(fid_out_bad,"ref DEM name\t%s\n",DEM_name_refoutfile);
        //printf("ref dem %s\n",DEM_name_refoutfile);
        fprintf(fid_out_bad,"DEM name\t\t\t\t\t\t\tTx[meter]\tTy[meter]\tTz(average)\tTz(median)\tTx_std[meter]\tTy_std[meter]\tTz_std(average)\tTz_std(median)\tControls_rho\tMean_all(avg)\tMedian_all(avg)\tdz_std(avg)\tMean_all(med)\tMedian_all(med)\tdz_std(med)\tNumberOfCPs\tprocessing time\n");
        */
        for(int ti = 1; ti < proinfo->number_of_images ; ti++)
        {
            time_t total_ST_iter = 0, total_ET_iter = 0;
            total_ST_iter = time(0);
            
            CSize tar_dem_size;
            double tar_minX,tar_maxY,tar_dx,tar_dy;
            
            unsigned char *tar_hill = NULL;
            CSize *data_size_tar = (CSize*)malloc(sizeof(CSize)*(py_level+1));
            
            //dem image
            float *DEM_tar = NULL;
            char DEM_name_outfile[500];
            char tar_hill_inputfile[500];
            
            char tar_dem_name[500];
            sprintf(tar_dem_name,"%s",args.Image[ti]);
            printf("dem name = %s\n",tar_dem_name);
            
            tar_dem_size = ReadGeotiff_info_dxy(tar_dem_name,&tar_minX,&tar_maxY,&tar_dx,&tar_dy);
            
            cols[0] = 0;
            cols[1] = tar_dem_size.width;
            
            rows[0] = 0;
            rows[1] = tar_dem_size.height;
            
            float type;
            DEM_tar= Readtiff_T(tar_dem_name,&tar_dem_size,cols,rows,&tar_dem_size,type);
            
            char *reffilename  = SetOutpathName(args.Image[ti]);
            char *ref_no_ext = remove_ext(reffilename);
            int char_size = strlen(ref_no_ext);
            char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
            {
                ref_fchar[c] = ref_no_ext[c];
            }
            ref_fchar[char_size-3] = '\0';
            sprintf(DEM_name_outfile,"%sdem",ref_fchar);
            free(ref_no_ext);
            free(ref_fchar);
            free(reffilename);
            
            ref_no_ext = remove_ext(args.Image[ti]);
            char_size = strlen(ref_no_ext);
            ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
            {
                ref_fchar[c] = ref_no_ext[c];
            }
            ref_fchar[char_size-3] = '\0';
            sprintf(tar_hill_inputfile,"%sdem_hillshade.tif",ref_fchar);
            printf("tar hill name %s\n",tar_hill_inputfile);
            free(ref_no_ext);
            free(ref_fchar);
            
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nDEM loading time ID %d\t%f\n\n",ti,total_gap);
            total_ST = time(0);
            */
            char co_dems_name[500];
            char copoly_dems_name[500];
            char co_dems_name_diff[500];
            char copoly_dems_name_diff[500];
            char ref_hill_name[500];
            char tar_hill_name[500];

            sprintf(co_dems_name,"%s/%s_coreg_avg.tif",args.Outputpath,DEM_name_outfile);
            sprintf(copoly_dems_name,"%s/%s_coreg_med.tif",args.Outputpath,DEM_name_outfile);
            sprintf(co_dems_name_diff,"%s/%s_coreg_avg_diff.tif",args.Outputpath,DEM_name_outfile);
            sprintf(copoly_dems_name_diff,"%s/%s_coreg_med_diff.tif",args.Outputpath,DEM_name_outfile);
            
            sprintf(ref_hill_name,"%s/%s_ref_hill.tif",args.Outputpath,DEM_name_outfile);
            sprintf(tar_hill_name,"%s/%s_tar_hill.tif",args.Outputpath,DEM_name_outfile);
            
            //create hillshadeimage
            FILE *checktif_tar = NULL;
            checktif_tar = fopen(tar_hill_inputfile,"r");
            if(checktif)
            {
                TIFF *tarhilltif = NULL;
                tarhilltif  = TIFFOpen(tar_hill_inputfile,"r");
                uint8 type;
                tar_hill = Readtiff_T(tar_hill_inputfile, &tar_dem_size, cols, rows, &tar_dem_size, type);
                TIFFClose(tarhilltif);
                fclose(checktif_tar);
                printf("loading existing hillshade image %s\n",tar_hill_inputfile);
            }
            else
                tar_hill = CreateHillshade(DEM_tar,tar_dem_size,tar_dx);
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nTar hillshade time %f\n\n",total_gap);
            total_ST = time(0);
             */
            
            //create pyramidimage
            SetPySizes(data_size_tar, tar_dem_size, py_level);
            
            unsigned char **SubImages_tar = (unsigned char**)malloc(sizeof(unsigned char*)*(py_level+1));
            
            for(int level = 0 ; level <= py_level ; level++)
            {
                CSize data_size;
                long int data_length;
                
                if(level == 0)
                {
                    data_size = data_size_tar[level];
                    
                    data_length = (long int)data_size.height*(long int)data_size.width;
                    SubImages_tar[level] = (unsigned char*)malloc(sizeof(unsigned char)*data_length);
                    memcpy(SubImages_tar[level],tar_hill,sizeof(unsigned char)*(long)data_length);
                }
                else
                {
                    data_size = data_size_tar[level-1];
                    
                    data_length = (long int)data_size.height*(long int)data_size.width;
                    //SubImages_tar[level]= CreateImagePyramid_BYTE(SubImages_tar[level-1],data_size,py_kernel_size,(double)1.6);
                    SubImages_tar[level]= CreateImagePyramid(SubImages_tar[level-1],data_size,py_kernel_size,(double)1.6);
                }
                //printf("preprocessing size level ID %d\t%d\t%d\t%d\n",level,ti,data_size.width,data_size.height);
            }
            free(tar_hill);
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nTar pyramid time %f\n\n",total_gap);
            total_ST = time(0);
            */
            double ImageBoundary_tar[4] = {0.0};
            double Boundary[4] = {0.0};
            double GridSize_width, GridSize_height;
            
            Boundary[0] = ImageBoundary_ref[0];
            Boundary[1] = ImageBoundary_ref[1];
            Boundary[2] = ImageBoundary_ref[2];
            Boundary[3] = ImageBoundary_ref[3];
            
            ImageBoundary_tar[0] = tar_minX;
            ImageBoundary_tar[1] = tar_maxY - tar_dy*tar_dem_size.height;
            ImageBoundary_tar[2] = tar_minX + tar_dx*tar_dem_size.width;
            ImageBoundary_tar[3] = tar_maxY;
            
            if(Boundary[0] < ImageBoundary_tar[0])
                Boundary[0] = ImageBoundary_tar[0];
            else
                Boundary[0] = Boundary[0];
            
            if(Boundary[1] < ImageBoundary_tar[1])
                Boundary[1] = ImageBoundary_tar[1];
            else
                Boundary[1] = Boundary[1];
            
            if(Boundary[2] > ImageBoundary_tar[2])
                Boundary[2] = ImageBoundary_tar[2];
            else
                Boundary[2] = Boundary[2];
            
            if(Boundary[3] > ImageBoundary_tar[3])
                Boundary[3] = ImageBoundary_tar[3];
            else
                Boundary[3] = Boundary[3];
            
            if(args.check_boundary)
            {
                if(Boundary[0] < args.Min_X)
                    Boundary[0] = args.Min_X;
                if(Boundary[1] < args.Min_Y)
                    Boundary[1] = args.Min_Y;
                if(Boundary[2] > args.Max_X)
                    Boundary[2] = args.Max_X;
                if(Boundary[3] > args.Max_Y)
                    Boundary[3] = args.Max_Y;
                
                Boundary[0] = args.Min_X;
                Boundary[1] = args.Min_Y;
                Boundary[2] = args.Max_X;
                Boundary[3] = args.Max_Y;
            }
            
            GridSize_width = Boundary[2] - Boundary[0];
            GridSize_height = Boundary[3] - Boundary[1];
            
            
            int iter_counts;
            int RA_iter_counts;
            D2DPOINT adjust_std;
            D2DPOINT *MPs_2D = NULL;
            double avg_roh;
            double ImageAdjust_coreg[2] = {0.0};
            
            for(int level = py_level ; level >= 0 ; level --)
            {
                double Grid_space = ceil((ref_dx + ref_dy)/2.0)*pwrtwo(py_level+2);
                //if(Grid_space < 50)
                    Grid_space = 50;
                //printf("image coregistration gridspace %d\t%d\n",level,Grid_space[level]);
                
                //printf("Processing level %d\n",level);
                //printf("level\t\trow(pixel)\tcolumn(pixel)\tTy(meter)\tTx(meter)\tGCPS #\tavg_roh\t# of iteration\n");
                
               
                
                MPs_2D = CoregParam_Image_MPs_stereo(proinfo, level,py_level, ImageAdjust_coreg, ncc_flag,
                                              5, SubImages_ref[level],SubImages_tar[level], data_size_ref,data_size_tar, ImageBoundary_ref, ImageBoundary_tar, ref_dx, ref_dy, tar_dx, tar_dy, Grid_space,Boundary,proinfo->save_filepath,&avg_roh,&iter_counts,&adjust_std,&RA_iter_counts);
                
                //printf("2 std %f\t%f\n",adjust_std[ti].m_X,adjust_std[ti].m_Y);
                printf("%d\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%4.2f\t%d\t%4.2f\t%4.2f\n",level,ImageAdjust_coreg[0], ImageAdjust_coreg[1],-ImageAdjust_coreg[0]*tar_dy, ImageAdjust_coreg[1]*tar_dx,RA_iter_counts,avg_roh,iter_counts,adjust_std.m_X,adjust_std.m_Y);
                free(SubImages_tar[level]);
            }
            free(SubImages_tar);
            free(data_size_tar);
            /*
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nCoregParam_Image time %f\n\n",total_gap);
            total_ST = time(0);
            */
            
            if(RA_iter_counts > 10)
            {
                time_t total_ST_a = 0, total_ET_a = 0;
                total_ST_a = time(0);
                
                double Coreg_param[2];
                
                CSize ref_data_size = ref_dem_size;
                CSize tar_data_size = tar_dem_size;
                
                Coreg_param[1] = ImageAdjust_coreg[1]*tar_dx;
                Coreg_param[0] = - ImageAdjust_coreg[0]*tar_dy;
                double ref_grid_size = ref_dx;
                double tar_grid_size = tar_dx;
                
                //printf("image ID %d\tcoreg y=%f\tx=%f\t%d\n",ti,Coreg_param[0],Coreg_param[1],RA_iter_counts[ti]);
                //printf("ref size %d\t%d\ttar size %d\t%d\n",ref_data_size.width,ref_data_size.height,tar_data_size.width,tar_data_size.height);
                //printf("ref coord %f\t%f\ttar coord %f\t%f\n",ref_minX,ref_maxY,tar_minX,tar_maxY);
                int diff_count = 0;
                
                double W_th = 80;
                double dH_th = 30;
                
                double average;
                float MED_z;
                double SD_z;
                double SD_z_med;
                
                bool check_stop = false;
                int while_iter = 0;
                double pre_std = 10000.0;
                double change_std_ratio = 1000;
                
                total_ST = time(0);
                bool check_cps = false;
                while(!check_stop && while_iter < 50 && !check_cps)
                {
                    double sum_diff = 0;
                    F3DPOINT *save_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),RA_iter_counts);
                    diff_count = 0;
    #pragma omp parallel for reduction(+:sum_diff,diff_count) schedule(guided)
                    for(long int gcp_index = 0 ; gcp_index < RA_iter_counts ; gcp_index ++)
                    {
                        D2DPOINT gcp_coord,ref_img,tar_img;
                        gcp_coord.m_X = MPs_2D[gcp_index].m_X;
                        gcp_coord.m_Y = MPs_2D[gcp_index].m_Y;
                        
                        //printf("ID %d\tMPs %f\t%f\n",gcp_index,gcp_coord.m_X,gcp_coord.m_Y);
                        
                        ref_img.m_X = (gcp_coord.m_X - ref_minX)/ref_grid_size;
                        ref_img.m_Y = (ref_maxY - gcp_coord.m_Y)/ref_grid_size;
                        
                        tar_img.m_X = ( (gcp_coord.m_X + Coreg_param[1]) - tar_minX )/tar_grid_size;
                        tar_img.m_Y = ( tar_maxY - (gcp_coord.m_Y + Coreg_param[0]) )/tar_grid_size;
                        
                        bool check_b = true;
                        if(args.check_boundary)
                        {
                            if(gcp_coord.m_X >= args.Min_X && gcp_coord.m_X <= args.Max_X && gcp_coord.m_Y >= args.Min_Y && gcp_coord.m_Y <= args.Max_Y)
                                check_b = true;
                            else
                                check_b = false;
                        }
                        
                        if(ref_img.m_X - 2 >= 0 && ref_img.m_X + 2 < ref_data_size.width && ref_img.m_Y - 2 >= 0 && ref_img.m_Y + 2 < ref_data_size.height &&
                           tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_data_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_data_size.height && check_b)
                        {
                            long ref_index = (int)(ref_img.m_Y)*ref_data_size.width + (int)(ref_img.m_X);
                            long tar_index = (int)(tar_img.m_Y)*tar_data_size.width + (int)(tar_img.m_X);
                            
                            if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                            {
                                double sum_zr = 0;
                                double sum_zt = 0;
                                int kernel_count = 0;
                                for(int kr = -2 ; kr <= 2 ; kr++ )
                                {
                                    for(int kc = -2 ; kc <= 2 ; kc++ )
                                    {
                                        int k_refindex = (int)(ref_img.m_Y + kr)*ref_data_size.width + (int)(ref_img.m_X + kc);
                                        int k_tarindex = (int)(tar_img.m_Y + kr)*tar_data_size.width + (int)(tar_img.m_X + kc);
                                        
                                        if(DEM[k_refindex] > -100 && DEM[k_refindex] < 10000 && DEM_tar[k_tarindex] > -100 && DEM_tar[k_tarindex] < 10000)
                                        {
                                            kernel_count++;
                                            //printf("DEM value %f\t%f\n",DEM[0][k_refindex],DEM[ti][k_tarindex]);
                                        }
                                    }
                                }
                                
                                double* ref_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                                double* tar_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                                
                                GMA_double *A_matrix = GMA_double_create(kernel_count, 3);
                                GMA_double *L_matrix = GMA_double_create(kernel_count, 1);
                                GMA_double *AT_matrix = GMA_double_create(3,kernel_count);
                                GMA_double *ATA_matrix = GMA_double_create(3,3);
                                
                                GMA_double *ATAI_matrix = GMA_double_create(3,3);
                                GMA_double *ATL_matrix = GMA_double_create(3,1);
                                
                                GMA_double *X_matrix = GMA_double_create(3,1);
                                GMA_double *AX_matrix = GMA_double_create(kernel_count,1);
                                GMA_double *V_matrix = GMA_double_create(kernel_count,1);
                                
                                GMA_double *tA_matrix = GMA_double_create(kernel_count, 3);
                                GMA_double *tL_matrix = GMA_double_create(kernel_count, 1);
                                GMA_double *tAT_matrix = GMA_double_create(3,kernel_count);
                                GMA_double *tATA_matrix = GMA_double_create(3,3);
                                
                                GMA_double *tATAI_matrix = GMA_double_create(3,3);
                                GMA_double *tATL_matrix = GMA_double_create(3,1);
                                
                                GMA_double *tX_matrix = GMA_double_create(3,1);
                                GMA_double *tAX_matrix = GMA_double_create(kernel_count,1);
                                GMA_double *tV_matrix = GMA_double_create(kernel_count,1);
                                
                                kernel_count = 0;
                                
                                sum_zr = 0;
                                sum_zt = 0;
                                
                                for(int kr = -2 ; kr <= 2 ; kr++ )
                                {
                                    for(int kc = -2 ; kc <= 2 ; kc++ )
                                    {
                                        int k_refindex = (int)(ref_img.m_Y + kr)*ref_data_size.width + (int)(ref_img.m_X + kc);
                                        int k_tarindex = (int)(tar_img.m_Y + kr)*tar_data_size.width + (int)(tar_img.m_X + kc);
                                        if(DEM[k_refindex] > -100 && DEM[k_refindex] < 10000 && DEM_tar[k_tarindex] > -100 && DEM_tar[k_tarindex] < 10000)
                                        {
                                            ref_patch_vecs[kernel_count] = DEM[k_refindex];
                                            tar_patch_vecs[kernel_count] = DEM_tar[k_tarindex];
                                            
                                            A_matrix->val[kernel_count][0] = kc*ref_grid_size;
                                            A_matrix->val[kernel_count][1] = kr*ref_grid_size;
                                            A_matrix->val[kernel_count][2] = 1.0;
                                            L_matrix->val[kernel_count][0] = DEM[k_refindex];
                                            
                                            tA_matrix->val[kernel_count][0] = kc*tar_grid_size;
                                            tA_matrix->val[kernel_count][1] = kr*tar_grid_size;
                                            tA_matrix->val[kernel_count][2] = 1.0;
                                            tL_matrix->val[kernel_count][0] = DEM_tar[k_tarindex];
                                            
                                            kernel_count++;
                                        }
                                    }
                                }
                                //double var_diff = fabs(sqrt(sum_zr/(double)kernel_count) - sqrt(sum_zt/(double)kernel_count));
                                
                                double k_ncc = Correlate(ref_patch_vecs, tar_patch_vecs, kernel_count);
                                if (k_ncc != -99)
                                    k_ncc = (k_ncc + 1)/2.0;
                                else
                                    k_ncc = 0.0;
                                
                                GMA_double_Tran(A_matrix,AT_matrix);
                                GMA_double_mul(AT_matrix,A_matrix,ATA_matrix);
                                GMA_double_inv(ATA_matrix,ATAI_matrix);
                                GMA_double_mul(AT_matrix,L_matrix,ATL_matrix);
                                GMA_double_mul(ATAI_matrix,ATL_matrix,X_matrix);
                                GMA_double_mul(A_matrix,X_matrix,AX_matrix);
                                GMA_double_sub(AX_matrix,L_matrix,V_matrix);
                                
                                GMA_double_Tran(tA_matrix,tAT_matrix);
                                GMA_double_mul(tAT_matrix,tA_matrix,tATA_matrix);
                                GMA_double_inv(tATA_matrix,tATAI_matrix);
                                GMA_double_mul(tAT_matrix,tL_matrix,tATL_matrix);
                                GMA_double_mul(tATAI_matrix,tATL_matrix,tX_matrix);
                                GMA_double_mul(tA_matrix,tX_matrix,tAX_matrix);
                                GMA_double_sub(tAX_matrix,tL_matrix,tV_matrix);
                                
                                double N1 = X_matrix->val[0][0];
                                double N2 = X_matrix->val[1][0];
                                double N3 = 1.0;
                                
                                double norm  = sqrt(N1*N1 + N2*N2 + N3*N3);
                                double angle = acos(fabs(N3)/norm)*180/3.141592;
                                
                                if(angle <= 0 && angle >= -90)
                                    angle = fabs(angle);
                                else if(angle <= -270 && angle >= -360)
                                    angle = 360 + angle;
                                else if(angle >= 270 && angle <= 360)
                                    angle = 360 - angle;
                                
                                double aspect = 90 - atan2(N2,N1)*180/3.141592;
                                
                                N1 = tX_matrix->val[0][0];
                                N2 = tX_matrix->val[1][0];
                                N3 = 1.0;
                                
                                norm  = sqrt(N1*N1 + N2*N2 + N3*N3);
                                double angle_tar = acos(fabs(N3)/norm)*180/3.141592;
                                
                                if(angle_tar <= 0 && angle_tar >= -90)
                                    angle_tar = fabs(angle_tar);
                                else if(angle_tar <= -270 && angle_tar >= -360)
                                    angle_tar = 360 + angle_tar;
                                else if(angle_tar >= 270 && angle_tar <= 360)
                                    angle_tar = 360 - angle_tar;
                                
                                double aspect_tar = 90 - atan2(N2,N1)*180/3.141592;
                                
                                double SS = 1 - fabs(angle - angle_tar)/90;
                                double SA = 1 - fabs(aspect - aspect_tar)/360;
                                
                                double W = 40*k_ncc + 40*SS + 20*SA;
                                
                                if(W > W_th && (DEM[ref_index] - DEM_tar[tar_index]) < dH_th + average && (DEM[ref_index] - DEM_tar[tar_index]) > - dH_th + average)
                                {
                                    long t_col_int   = (long)(tar_img.m_X + 0.01);
                                    long t_row_int   = (long)(tar_img.m_Y + 0.01);
                                    
                                    double dcol        = tar_img.m_X - t_col_int;
                                    double drow        = tar_img.m_Y - t_row_int;
                                    
                                    long index1,index2,index3, index4;
                                    double value1, value2, value3, value4, value;
                                    
                                    index1  = (t_col_int   ) + (t_row_int   )*(long)tar_data_size.width;
                                    index2  = (t_col_int +1) + (t_row_int   )*(long)tar_data_size.width;
                                    index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_data_size.width;
                                    index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_data_size.width;
                                    
                                    value1      = DEM_tar[index1];
                                    value2      = DEM_tar[index2];
                                    value3      = DEM_tar[index3];
                                    value4      = DEM_tar[index4];
                                    
                                    if(value1 > -100 && value1 < 10000 && value2 > -100 && value2 < 10000 && value3 > -100 && value3 < 10000 && value4 > -100 && value4 < 10000)
                                    {
                                        value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                        + value3*(1-dcol)*drow + value4*dcol*drow;
                                        
                                        if((DEM[ref_index] - value) < dH_th + average && (DEM[ref_index] - value) > -dH_th + average)
                                        {
                                            sum_diff += (DEM[ref_index] - value);
                                            diff_count++;
                                            
                                            save_pts[gcp_index].m_X = gcp_coord.m_X + Coreg_param[1];
                                            save_pts[gcp_index].m_Y = gcp_coord.m_Y + Coreg_param[0];
                                            save_pts[gcp_index].m_Z = (DEM[ref_index] - value);
                                            save_pts[gcp_index].flag = true;
                                            //printf("W kncc ss sa %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",W,k_ncc,SS,SA,angle,angle_tar,aspect,aspect_tar,save_pts[gcp_index].m_Z);
                                        }
                                        else
                                        {
                                            save_pts[gcp_index].flag = false;
                                        }
                                    }
                                }
                                
                                //if(W > 90 && fabs(ref_dem[ref_index] - tar_dem[tar_index]) >= 30)
                                //    printf("W dz %f\t%f\t%f\t%f\t%f\n",W,k_ncc,SS,SA,ref_dem[ref_index] - tar_dem[tar_index]);
                                
                                GMA_double_destroy(A_matrix);
                                GMA_double_destroy(L_matrix);
                                GMA_double_destroy(AT_matrix);
                                GMA_double_destroy(ATA_matrix);
                                GMA_double_destroy(ATAI_matrix);
                                GMA_double_destroy(ATL_matrix);
                                GMA_double_destroy(X_matrix);
                                GMA_double_destroy(AX_matrix);
                                GMA_double_destroy(V_matrix);
                                
                                GMA_double_destroy(tA_matrix);
                                GMA_double_destroy(tL_matrix);
                                GMA_double_destroy(tAT_matrix);
                                GMA_double_destroy(tATA_matrix);
                                GMA_double_destroy(tATAI_matrix);
                                GMA_double_destroy(tATL_matrix);
                                GMA_double_destroy(tX_matrix);
                                GMA_double_destroy(tAX_matrix);
                                GMA_double_destroy(tV_matrix);
                                
                                free(ref_patch_vecs);
                                free(tar_patch_vecs);
                            }
                        }
                    }
                    /*
                    total_ET = time(0);
                    total_gap = difftime(total_ET,total_ST);
                    printf("\nselect vertical CPs time %f\t%d\n\n",total_gap,diff_count);
                    //total_ST = time(0);
                    */
                    printf("iteration %d\tcps %d\n",while_iter,diff_count);
                    if(diff_count > 10)
                    {
                    
                        average = sum_diff/diff_count;
                        
                        float* save_Z = (float*)calloc(sizeof(float),diff_count);
                        float sum_var = 0;
                        
                        long count = 0;
                        
                        char dem_gcp_filename[500];
                        sprintf(dem_gcp_filename,"%s/DEM_gcps_%d.txt",args.Outputpath,ti);
                        FILE *fid_dem_gcp = fopen(dem_gcp_filename,"w");
                        
                        for(long row = 0; row < RA_iter_counts ; row++)
                        {
                            if(save_pts[row].flag == true)
                            {
                                save_Z[count] = save_pts[row].m_Z;
                                
                                sum_var += (average - save_pts[row].m_Z)*(average - save_pts[row].m_Z);
                                fprintf(fid_dem_gcp,"%f\t%f\t%f\n",save_pts[row].m_X,save_pts[row].m_Y,save_pts[row].m_Z);
                                count++;
                            }
                        }
                        fclose(fid_dem_gcp);
                    
                        MED_z = quickselect(save_Z, diff_count, (int)(diff_count/2.0));
                        //float MED_z = median(diff_count,save_Z,-(dH_th+5),dH_th+5);
                        
                        SD_z = sqrt(sum_var/diff_count);
                        
                        
                        float iter_dH_th = SD_z*3.29;
                        
                        change_std_ratio = fabs(pre_std - SD_z)/pre_std;
                        
                        if(change_std_ratio < 0.01 || dH_th < iter_dH_th)
                            check_stop = true;
                        
                        pre_std = SD_z;
                        dH_th = iter_dH_th;
                        
                        while_iter++;
                        
                        double sum_var_med = 0;
                        
            #pragma omp parallel for reduction(+:sum_var_med) schedule(guided)
                        for(int i=0;i<diff_count;i++)
                        {
                            sum_var_med += (MED_z -save_Z[i])*(MED_z -save_Z[i]);
                        }
                        
                        SD_z_med = sqrt(sum_var_med/diff_count);
                        
                        free(save_Z);
                    }
                    else
                        check_cps = true;
                    
                    free(save_pts);
                    
                    total_ET = time(0);
                    total_gap = difftime(total_ET,total_ST);
                    printf("\nTz average and variation cal time %d\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n\n",while_iter, total_gap,diff_count,average,MED_z,SD_z,SD_z_med,change_std_ratio,dH_th);
                    total_ST = time(0);
                }
                
                if(!check_cps)
                {
                    
                    float min_save_Z = 100000;
                    float max_save_Z = -100000;
                    
                    float* co_dem = NULL;
                    float* copoly_dem = NULL;
                    float* co_dem_diff = NULL;
                    float* copoly_dem_diff = NULL;
                    float* save_dz = NULL;
                    float* save_dz_med = NULL;
                    
                    float all_average = 0.0;
                    float all_med = 0.0;
                    float all_std = 0.0;
                    float all_average_med = 0.0;
                    float all_med_med = 0.0;
                    float all_std_med = 0.0;
                    
                    if(args.check_DEM_coreg_output == 2)
                    {
                        co_dem = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                        copoly_dem = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                        co_dem_diff = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                        copoly_dem_diff = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                        
                        for(long int co_index = 0 ; co_index < (long)(tar_data_size.width)*(long)(tar_data_size.height) ; co_index++)
                        {
                            co_dem[co_index] = Nodata;
                            copoly_dem[co_index] = Nodata;
                            co_dem_diff[co_index] = Nodata;
                            copoly_dem_diff[co_index] = Nodata;
                        }
                    }
                    
                    if(args.check_DEM_coreg_output == 1 || args.check_DEM_coreg_output == 2)
                    {
                        save_dz = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                        save_dz_med = (float*)calloc(sizeof(float),tar_data_size.width*tar_data_size.height);
                    
                        float dz_sum = 0;
                        float dz_sum_med = 0;
                        long dz_count = 0;
                        min_save_Z = 100000;
                        max_save_Z = -100000;
                        for(long int co_index = 0 ; co_index < (long)(tar_data_size.width)*(long)(tar_data_size.height) ; co_index++)
                        {
                            int pts_row = (int)(floor(co_index/tar_data_size.width));
                            int pts_col = co_index % tar_data_size.width;
                            
                            D2DPOINT gcp_coord,tar_img, ref_img,gcp_coord_ref;
                            gcp_coord.m_X = pts_col*tar_grid_size + tar_minX + Coreg_param[1];
                            gcp_coord.m_Y = tar_maxY - pts_row*tar_grid_size + Coreg_param[0];
                            
                            tar_img.m_X = ( gcp_coord.m_X - tar_minX )/tar_grid_size;
                            tar_img.m_Y = ( tar_maxY - gcp_coord.m_Y )/tar_grid_size;
                            
                            gcp_coord_ref.m_X = pts_col*tar_grid_size + tar_minX;
                            gcp_coord_ref.m_Y = tar_maxY - pts_row*tar_grid_size;
                            
                            ref_img.m_X = ( gcp_coord_ref.m_X - ref_minX )/ref_grid_size;
                            ref_img.m_Y = ( ref_maxY - gcp_coord_ref.m_Y )/ref_grid_size;
                            
                            if(tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_data_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_data_size.height &&
                               ref_img.m_X >= 0 && ref_img.m_X < ref_data_size.width && ref_img.m_Y >= 0 && ref_img.m_Y < ref_data_size.height)
                            {
                                long tar_index = (int)(tar_img.m_Y)*tar_data_size.width + (int)(tar_img.m_X);
                                long ref_index = (int)(ref_img.m_Y)*ref_data_size.width + (int)(ref_img.m_X);
                                
                                if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                                {
                                    long t_col_int   = (long)(tar_img.m_X + 0.01);
                                    long t_row_int   = (long)(tar_img.m_Y + 0.01);
                                    
                                    double dcol        = tar_img.m_X - t_col_int;
                                    double drow        = tar_img.m_Y - t_row_int;
                                    
                                    long index1,index2,index3, index4;
                                    double value1, value2, value3, value4, value;
                                    
                                    index1  = (t_col_int   ) + (t_row_int   )*(long)tar_data_size.width;
                                    index2  = (t_col_int +1) + (t_row_int   )*(long)tar_data_size.width;
                                    index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_data_size.width;
                                    index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_data_size.width;
                                    
                                    value1      = DEM_tar[index1];
                                    value2      = DEM_tar[index2];
                                    value3      = DEM_tar[index3];
                                    value4      = DEM_tar[index4];
                                    
                                    value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                    + value3*(1-dcol)*drow + value4*dcol*drow;
                                    
                                    //double value       = DEM[ti][tar_index];
                                    float co_dem_t = value + average;
                                    float copoly_dem_t = value + MED_z;
                                    float co_dem_diff_t = (DEM[ref_index] - co_dem_t);
                                    float copoly_dem_diff_t = (DEM[ref_index] - copoly_dem_t);
                                    
                                    if(args.check_DEM_coreg_output == 2)
                                    {
                                        co_dem[tar_index] = co_dem_t;
                                        copoly_dem[tar_index] = copoly_dem_t;
                                        
                                        co_dem_diff[tar_index] = co_dem_diff_t;
                                        copoly_dem_diff[tar_index] = copoly_dem_diff_t;
                                    }
                                    
                                    
                                    if(fabs(co_dem_diff_t) < dH_th && fabs(copoly_dem_diff_t) < dH_th)
                                    {
                                        dz_sum += co_dem_diff_t;
                                        dz_sum_med += copoly_dem_diff_t;
                                        
                                        save_dz[dz_count] = co_dem_diff_t;
                                        save_dz_med[dz_count] = copoly_dem_diff_t;
                                        
                                        dz_count++;
                                    }
                                }
                            }
                        }
                        /*
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\nWhole image stat %f\n\n",total_gap);
                        total_ST = time(0);
                        */
                        all_average = dz_sum/dz_count;
                        all_med = quickselect(save_dz, dz_count, (int)(dz_count/2.0));
                        //float all_med = binmedian(dz_count, save_dz);
                        //float all_med = median(dz_count,save_dz,-(dH_th + 5),dH_th + 5);
                        float all_sum_var = 0;
                        
                        all_average_med = dz_sum_med/dz_count;
                        all_med_med = quickselect(save_dz_med, dz_count, (int)(dz_count/2.0));
                        //float all_med_med = binmedian(dz_count,save_dz_med);
                        //float all_med_med = median(dz_count,save_dz_med,-(dH_th + 5),dH_th + 5);
                        float all_sum_var_med = 0;
                        
            #pragma omp parallel for reduction(+:all_sum_var,all_sum_var_med) schedule(guided)
                        for(long index = 0 ; index < dz_count ; index++)
                        {
                            all_sum_var += (save_dz[index] - all_average)*(save_dz[index] - all_average);
                            all_sum_var_med += (save_dz_med[index] - all_average_med)*(save_dz_med[index] - all_average_med);
                        }
                        //printf("done allsum\n");
                        free(save_dz);
                        free(save_dz_med);
                        
                        all_std = sqrt(all_sum_var/dz_count);
                        all_std_med = sqrt(all_sum_var_med/dz_count);
                        /*
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\nall median median time %f\n\n",total_gap);
                        total_ST = time(0);
                        */
                        printf("all stat %f\t%f\t%f\t%f\t%f\t%f\n",all_average,all_med,all_std,all_average_med,all_med_med,all_std_med);
                        
                        if(args.check_DEM_coreg_output == 2)
                        {
                            WriteGeotiff(co_dems_name, co_dem, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            WriteGeotiff(copoly_dems_name, copoly_dem, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            
                            WriteGeotiff(co_dems_name_diff, co_dem_diff, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            WriteGeotiff(copoly_dems_name_diff, copoly_dem_diff, tar_data_size.width, tar_data_size.height, tar_grid_size, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            
                            free(co_dem);
                            free(copoly_dem);
                            free(co_dem_diff);
                            free(copoly_dem_diff);
                        }
                        
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\noutput save time %f\n\n",total_gap);
                        total_ST = time(0);
                    }
                  
                    /*
                    total_ET_a = time(0);
                    total_gap = difftime(total_ET_a,total_ST_a);
                    printf("\nTotal Tz computation time ID %d\t%f\n\n",ti,total_gap);*/
                    total_ET_iter = time(0);
                    total_gap = difftime(total_ET_iter,total_ST_iter);
                    //printf("\nTotal computation time ID %d\t%f\n\n",ti,total_gap);
                    
                    if(avg_roh > 0.7 && SD_z < 5 && diff_count > 50)
                    { fprintf(fid_out,"%s.tif\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%f\n",DEM_name_outfile,Coreg_param[1],Coreg_param[0],average,MED_z,adjust_std.m_Y,adjust_std.m_X,SD_z,SD_z_med,avg_roh,
                            all_average,all_med,all_std,all_average_med,all_med_med,all_std_med,diff_count,total_gap);
                    }
                    else
                    {
                        fprintf(fid_out,"%s.tif\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%f\n",DEM_name_outfile,Coreg_param[1],Coreg_param[0],average,MED_z,adjust_std.m_Y,adjust_std.m_X,SD_z,SD_z_med,avg_roh,
                        all_average,all_med,all_std,all_average_med,all_med_med,all_std_med,diff_count,total_gap);
                    }
                }
                else
                {
                    fprintf(fid_out,"%s.tif\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t%f\n",DEM_name_outfile,total_gap);
                }
                    
            }
            else
            {
                fprintf(fid_out,"%s.tif\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t%f\n",DEM_name_outfile,total_gap);
            }
            
            if(MPs_2D)
                free(MPs_2D);
            free(DEM_tar);
        }
        fclose(fid_out);
        
        for(int level = py_level ; level >= 0 ; level --)
        {
            free(SubImages_ref[level]);
        }
        free(SubImages_ref);
        free(DEM);
        free(data_size_ref);
        free(ref_hill);
    }
}

void DEM_ImageCoregistration_GeomatricConstraint(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt)
{
    ProInfo *proinfo = (ProInfo*)malloc(sizeof(ProInfo));
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    FILE *time_fid;
    
    total_ST = time(0);
    int py_kernel_size = 3;
    
    bool check_start = true;
    if(args.check_txt_input == 0)
    {
        check_start = OpenProject(_filename,proinfo,args);
    }
    else if(args.check_txt_input == 1)
    {
        FILE *fid_read = fopen(args.DEM_input_file,"r");
        
        args.number_of_images = 0;
        while( (fscanf(fid_read,"%s\n",args.Image[args.number_of_images])) != EOF )
        {
            printf("input DEM %s\n",args.Image[args.number_of_images]);
            sprintf(proinfo->Imagefilename[args.number_of_images],"%s",args.Image[args.number_of_images]);
            args.number_of_images++;
        }
        fclose(fid_read);
        proinfo->number_of_images = args.number_of_images;
        printf("image num %d\n",args.number_of_images);
    }
    
    if(check_start)
    {
        char temp_filepath[500];
        
        int check_folder = 1;
        
        sprintf(proinfo->save_filepath,"%s",args.Outputpath);
        
        int status;
        status = mkdir(proinfo->save_filepath,0777);
        if (opendir(proinfo->save_filepath) == NULL)
        {
            if (status == -1)
            {
                printf("Outpath of '%s' cannot make, please check outpath!!\n",proinfo->save_filepath);
                exit(1);
            }
        }
        
        sprintf(temp_filepath,"%s/txt",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(temp_filepath,"%s/tif",proinfo->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(proinfo->tmpdir,"%s/tmp",proinfo->save_filepath);
        mkdir(proinfo->tmpdir,0777);
        
     
        long int cols[2], rows[2];
        double Sum_grid = 0;
        NCCflag ncc_flag;
        int py_level = proinfo->pyramid_level;
        char save_file[500];
        char out_file[500];
        char out_file_bad[500];
        CSize ref_dem_size;
        double ref_minX,ref_maxY,ref_dx,ref_dy;
        
        ncc_flag.rotate_flag     = 1;       ncc_flag.multi_flag      = 0;       ncc_flag.multi_flag_sum  = 0;       ncc_flag.inter_flag      = 1;
        ncc_flag.weight_flag     = 0;
        
        TransParam param;
        SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
        
        //dem image = ref
        float *DEM = NULL;
        CSize *data_size_ref = (CSize*)malloc(sizeof(CSize)*(py_level+1));
        
        
        char DEM_name_refoutfile[500];
        char ref_hill_inputname[500];
        char ref_dem_name[500];
        
        sprintf(ref_dem_name,"%s",args.Image[0]);
        printf("dem name = %s\n",ref_dem_name);
        
        ref_dem_size = ReadGeotiff_info_dxy(ref_dem_name,&ref_minX,&ref_maxY,&ref_dx,&ref_dy);
        
        double ImageBoundary_ref[4];
        ImageBoundary_ref[0] = ref_minX;
        ImageBoundary_ref[1] = ref_maxY - ref_dy*ref_dem_size.height;
        ImageBoundary_ref[2] = ref_minX + ref_dx*ref_dem_size.width;
        ImageBoundary_ref[3] = ref_maxY;
        
        cols[0] = 0;
        cols[1] = ref_dem_size.width;
        
        rows[0] = 0;
        rows[1] = ref_dem_size.height;
        
        float type;
        DEM = Readtiff_T(ref_dem_name,&ref_dem_size,cols,rows,&ref_dem_size,type);
        
        char *reffilename  = SetOutpathName(args.Image[0]);
        char *ref_no_ext = remove_ext(reffilename);
        int char_size = strlen(ref_no_ext);
        char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
        for(int c = 0 ; c < char_size - 3 ; c++)
        {
            ref_fchar[c] = ref_no_ext[c];
        }
        ref_fchar[char_size-3] = '\0';
        sprintf(DEM_name_refoutfile,"%s",ref_fchar);
        free(ref_no_ext);
        free(ref_fchar);
        free(reffilename);
        
        F3DPOINT *select_pts = NULL;
        long ref_tin_point_num = 0;
        
        CSize data_size;
        long int data_length;
        char pyimage[500];
        char pyslope[500];
        char pyascpect[500];
        
        bool check_slope_aspect = false;
        long tin_point_num = 0;
        long tri_count = 0;
        
        
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        printf("\nRef preparation time %f\n\n",total_gap);
        total_ST = time(0);
        
        sprintf(out_file,"%s/DEM_coreg_result.txt",args.Outputpath);
        FILE* fid_out         = fopen(out_file,"w");
        fprintf(fid_out,"ref DEM name\t%s\n",DEM_name_refoutfile);
        fprintf(fid_out,"DEM name\t\t\t\t\t\t\tDist_std[meter]\tTx[meter]\tTy[meter]\tTz[meter]\tTx_std[meter]\tTy_std[meter]\tTz_std[meter]\tdH_cp(mean)\tdH_cp_std(mean)\tdH_cp(med.)\tdH_cp_std(med.)\tdh_mean\t\tdh_med.\t\tdh_std\t\tNumberOfCPs\tprocessing time\n");
        /*
        sprintf(out_file_bad,"%s/DEM_coreg_result_bad.txt",args.Outputpath);
        FILE* fid_out_bad         = fopen(out_file_bad,"w");
        fprintf(fid_out_bad,"ref DEM name\t%s\n",DEM_name_refoutfile);
        fprintf(fid_out_bad,"DEM name\t\t\t\t\t\t\tDist_std[meter]\tTx[meter]\tTy[meter]\tTz[meter]\tTx_std[meter]\tTy_std[meter]\tTz_std[meter]\tdH_cp(mean)\tdH_cp_std(mean)\tdH_cp(med.)\tdH_cp_std(med.)\tdh_mean\t\tdh_med.\t\tdh_std\t\tNumberOfCPs\tprocessing time\n");
        */
        for(int ti = 1; ti < proinfo->number_of_images ; ti++)
        {
            time_t total_ST_iter = 0, total_ET_iter = 0;
            total_ST_iter = time(0);
            
            CSize tar_dem_size;
            double tar_minX,tar_maxY,tar_dx,tar_dy;
            
            CSize *data_size_tar = (CSize*)malloc(sizeof(CSize)*(py_level+1));
            
            //dem image
            float *DEM_tar = NULL;
            char DEM_name_outfile[500];
            char tar_hill_inputfile[500];
            
            char tar_dem_name[500];
            sprintf(tar_dem_name,"%s",args.Image[ti]);
            printf("dem name = %s\n",tar_dem_name);
            
            tar_dem_size = ReadGeotiff_info_dxy(tar_dem_name,&tar_minX,&tar_maxY,&tar_dx,&tar_dy);
            
            double ImageBoundary_tar[4] = {0.0};
            ImageBoundary_tar[0] = tar_minX;
            ImageBoundary_tar[1] = tar_maxY - tar_dy*tar_dem_size.height;
            ImageBoundary_tar[2] = tar_minX + tar_dx*tar_dem_size.width;
            ImageBoundary_tar[3] = tar_maxY;
            
            cols[0] = 0;
            cols[1] = tar_dem_size.width;
            
            rows[0] = 0;
            rows[1] = tar_dem_size.height;
            
            float type;
            DEM_tar= Readtiff_T(tar_dem_name,&tar_dem_size,cols,rows,&tar_dem_size,type);
            
            char *reffilename  = SetOutpathName(args.Image[ti]);
            char *ref_no_ext = remove_ext(reffilename);
            int char_size = strlen(ref_no_ext);
            char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
            {
                ref_fchar[c] = ref_no_ext[c];
            }
            ref_fchar[char_size-3] = '\0';
            sprintf(DEM_name_outfile,"%sdem",ref_fchar);
            free(ref_no_ext);
            free(ref_fchar);
            free(reffilename);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nDEM loading time ID %d\t%f\n\n",ti,total_gap);
            total_ST = time(0);
            
            char co_dems_name[500];
            char co_dems_name_diff[500];
            
            sprintf(co_dems_name,"%s/%s_coreg.tif",args.Outputpath,DEM_name_outfile);
            sprintf(co_dems_name_diff,"%s/%s_coreg_diff.tif",args.Outputpath,DEM_name_outfile);
            
            //overlapped boudnary
            double Boundary[4] = {0.0};
            CSize GridSize;
            Boundary[0] = ImageBoundary_ref[0];
            Boundary[1] = ImageBoundary_ref[1];
            Boundary[2] = ImageBoundary_ref[2];
            Boundary[3] = ImageBoundary_ref[3];
            
            printf("ref boundary %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
            printf("tar boundary %f\t%f\t%f\t%f\n",ImageBoundary_tar[0],ImageBoundary_tar[1],ImageBoundary_tar[2],ImageBoundary_tar[3]);
            if(Boundary[0] < ImageBoundary_tar[0])
                Boundary[0] = ImageBoundary_tar[0];
            else
                Boundary[0] = Boundary[0];
            
            if(Boundary[1] < ImageBoundary_tar[1])
                Boundary[1] = ImageBoundary_tar[1];
            else
                Boundary[1] = Boundary[1];
            
            if(Boundary[2] > ImageBoundary_tar[2])
                Boundary[2] = ImageBoundary_tar[2];
            else
                Boundary[2] = Boundary[2];
            
            if(Boundary[3] > ImageBoundary_tar[3])
                Boundary[3] = ImageBoundary_tar[3];
            else
                Boundary[3] = Boundary[3];
            
            if(args.check_boundary)
            {
                if(Boundary[0] < args.Min_X)
                    Boundary[0] = args.Min_X;
                if(Boundary[1] < args.Min_Y)
                    Boundary[1] = args.Min_Y;
                if(Boundary[2] > args.Max_X)
                    Boundary[2] = args.Max_X;
                if(Boundary[3] > args.Max_Y)
                    Boundary[3] = args.Max_Y;
                
                Boundary[0] = args.Min_X;
                Boundary[1] = args.Min_Y;
                Boundary[2] = args.Max_X;
                Boundary[3] = args.Max_Y;
            }
            printf("overlapped boundary %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
            
            //create pyramidimage
            SetPySizes(data_size_tar, tar_dem_size, py_level);
            
            F3DPOINT *select_pts_tar = NULL;// = (F3DPOINT*)malloc(sizeof(F3DPOINT));
            long tin_point_num = 0;
            
            //TINinfo tininfo_tar;
            float min_H = 100000000;
            float max_H = -100000000;
            
            CSize data_size;
            long int data_length;
            char pyimage[500];
            char pyslope[500];
            char pyascpect[500];
            char pyncc[500];
            
            select_pts_tar = SettingControls(DEM, DEM_tar, ref_dx, tar_dx, ImageBoundary_ref, ImageBoundary_tar, Boundary, ref_dem_size, tar_dem_size, &tin_point_num);
            
            printf("select point %d\n",tin_point_num);
            
            if(tin_point_num > 0)
            {
                
                total_ET = time(0);
                total_gap = difftime(total_ET,total_ST);
                printf("\nTar pyramid time %f\n\n",total_gap);
                total_ST = time(0);
                
                for(int index = 0 ; index < tin_point_num ; index ++)
                {
                    if(select_pts_tar[index].flag)
                    {
                        if(min_H > select_pts_tar[index].m_Z)
                            min_H = select_pts_tar[index].m_Z;
                        if(max_H < select_pts_tar[index].m_Z)
                            max_H = select_pts_tar[index].m_Z;
                        
                    }
                }
                /*
                float interval = 10.0;
                int max_pos;
                bool check_hist_stop = false;
                float sd_h;
                while(!check_hist_stop && interval < 500)
                {
                    int H_interval = floor(max_H/interval);
                    int *hist_H = (int*)calloc(sizeof(int*),H_interval);
                    
                    for(int index = 0 ; index < tin_point_num ; index ++)
                    {
                        int Height = floor(select_pts_tar[index].m_Z/interval);
                        if(Height >= 0 && Height <= floor(max_H))
                        {
                            hist_H[Height]++;
                        }
                    }
                    int max_hist = -10;
                    
                    int sum_hist = 0;
                    
                    for(int index = 0 ; index < H_interval ; index ++)
                    {
                        printf("hist id %d\t%d\n",index,hist_H[index]);
                        if(max_hist < hist_H[index])
                        {
                            max_hist = hist_H[index];
                            max_pos = index*(int)(interval);
                        }
                        sum_hist += hist_H[index]*index*interval;
                    }
                    
                    free(hist_H);
                    
                    if(max_hist > tin_point_num*0.5)
                    {
                        check_hist_stop = true;
                    }
                    else
                        interval += 10;
                    
                    float sum_var_h = 0;
                    for(int index = 0 ; index < tin_point_num ; index ++)
                    {
                        int Height = floor(select_pts_tar[index].m_Z);
                        sum_var_h = (max_pos - Height)*(max_pos - Height);
                    }
                    
                    sd_h = sqrt(sum_var_h/tin_point_num);
                    
                    printf("max hist %d\t%d\t%f\t%f\n",max_hist,max_pos,sd_h,interval);
                }
                
                min_H = max_pos*interval - sd_h*3.29;
                max_H = (max_pos+1)*interval + sd_h*3.29;
                */
                F3DPOINT coord_center, coord_scale;
                float scale_factor = 1;
                coord_center.m_X = (Boundary[2] - Boundary[0])/2.0 + Boundary[0];
                coord_center.m_Y = (Boundary[3] - Boundary[1])/2.0 + Boundary[1];
                coord_center.m_Z = (max_H - min_H)/2.0 + min_H;
                coord_scale.m_X = (Boundary[2] - Boundary[0])/2.0*scale_factor;
                coord_scale.m_Y = (Boundary[3] - Boundary[1])/2.0*scale_factor;
                coord_scale.m_Z = (max_H - min_H)/2.0*scale_factor;
                float ori_scale_Z = coord_scale.m_Z;
                
                float distance_scale = (sqrt(coord_scale.m_X*coord_scale.m_X + coord_scale.m_Y*coord_scale.m_Y + coord_scale.m_Z*coord_scale.m_Z));
                printf("minmaxH %f\t%f\t%f\t%f\t%f\t%f\n",min_H,max_H,coord_scale.m_X,coord_scale.m_Y,coord_scale.m_Z,distance_scale);
                
                F3DPOINT min_nor_coord,max_nor_coord;
                min_nor_coord.m_X = (Boundary[0] - coord_center.m_X)/coord_scale.m_X;
                min_nor_coord.m_Y = (Boundary[1] - coord_center.m_Y)/coord_scale.m_Y;
                min_nor_coord.m_Z = (min_H - coord_center.m_Z)/coord_scale.m_Z;
                
                max_nor_coord.m_X = (Boundary[2] - coord_center.m_X)/coord_scale.m_X;
                max_nor_coord.m_Y = (Boundary[3] - coord_center.m_Y)/coord_scale.m_Y;
                max_nor_coord.m_Z = (max_H - coord_center.m_Z)/coord_scale.m_Z;
                
                printf("total pts %d\n",tin_point_num);
                F3DPOINT *control_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                F3DPOINT *transformed_coord = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                float* SDistance = (float*)calloc(sizeof(float),tin_point_num);
                float* SDistance_ori = (float*)calloc(sizeof(float),tin_point_num);
                float* weight = (float*)calloc(sizeof(float),tin_point_num);
                F3DPOINT *tar_normal = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                F3DPOINT *tar_normal_ori = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                
                F3DPOINT* normalized_pt = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                F3DPOINT *ref_normal = (F3DPOINT*)calloc(sizeof(F3DPOINT),tin_point_num);
                
                float *ref_iter_height = (float*)calloc(sizeof(float),tin_point_num);
                double **ref_array = (double**)calloc(sizeof(double*),tin_point_num);
                for(int index = 0 ; index < tin_point_num ; index ++)
                {
                    if(select_pts_tar[index].flag)
                    {
                        ref_array[index] = (double*)calloc(sizeof(double),9);
                        transformed_coord[index].flag = true;
                    }
                }
                
                double level_ref_dx = ref_dx;
                double level_tar_dx = tar_dx;
                
                GridSize.width = ceil((Boundary[2] - Boundary[0])/level_ref_dx);
                GridSize.height = ceil((Boundary[3] - Boundary[1])/level_ref_dx);
                
                printf("GridSize %d\t%d\n",GridSize.width,GridSize.height);

                //select candidata pts and compute surface distance
                
                Conformalparam conparam;
                
                conparam.scale = 1.0;
                conparam.omega = 0.0;
                conparam.phi = 0.0;
                conparam.kappa = 0.0;
                conparam.Tx = 0.0;
                conparam.Ty = 0.0;
                conparam.Tz = 0.0;
                
                float* dX = NULL;
                float* sigmaX = (float*)calloc(sizeof(float),7);
                float sigma0;
                int while_iter = 0;
                bool check_stop = false;
                float W_th = 0.8;
                time_t total_ST_ad = 0, total_ET_ad = 0;
                total_ST_ad = time(0);
                
                float sum_distance = 0;
                long final_number_of_selected = 0;
                float average_distance, MED_distance;
                double vertical_average_distance;
                double SD_distance=0;
                float SD_distance_vertical=0;
                double SD_z_med = 0;
                double th_dH = 30.0/ori_scale_Z;
                double th_dH_vertical = 30.0;
                double pre_sigma = 10000.0;
                bool check_null_dem = false;
                
                float min_X = Boundary[0];
                float max_X = Boundary[2];
                
                float min_Y = Boundary[1];
                float max_Y = Boundary[3];
                
                while(!check_stop && while_iter < 50)
                {
                    /*coord_center.m_X = (max_X - min_X)/2.0 + min_X;
                    coord_center.m_Y = (max_Y - min_Y)/2.0 + min_Y;
                    coord_center.m_Z = (max_H - min_H)/2.0 + min_H;
                    coord_scale.m_X = (max_X - min_X)/2.0*scale_factor;
                    coord_scale.m_Y = (max_Y - min_Y)/2.0*scale_factor;*/
                    //coord_scale.m_Z = (max_H - min_H)/2.0*scale_factor;
                    /*
                    coord_scale.m_X = (max_X - min_X)/2.0*scale_factor;
                    coord_scale.m_Y = (max_Y - min_Y)/2.0*scale_factor;
                    coord_scale.m_Z = (max_H - min_H)/2.0*scale_factor;
                    
                    min_nor_coord.m_X = (min_X - coord_center.m_X)/coord_scale.m_X;
                    min_nor_coord.m_Y = (min_Y - coord_center.m_Y)/coord_scale.m_Y;
                    min_nor_coord.m_Z = (min_H - coord_center.m_Z)/coord_scale.m_Z;
                    
                    max_nor_coord.m_X = (max_X - coord_center.m_X)/coord_scale.m_X;
                    max_nor_coord.m_Y = (max_Y - coord_center.m_Y)/coord_scale.m_Y;
                    max_nor_coord.m_Z = (max_H - coord_center.m_Z)/coord_scale.m_Z;
                    */
                    printf("iter %d\tscale %f\t%f\t%f\t%f\n",while_iter,coord_scale.m_Z,min_H,max_H,W_th );
                    /*
                    min_H = 100000000;
                    max_H = -100000000;
                    
                    min_X = 100000000;
                    max_X = -100000000;
                    
                    min_Y = 100000000;
                    max_Y = -100000000;
                    */
                    
                    int* hist_W = (int*)calloc(sizeof(int),100);
                    
                    long number_of_selected = 0;
                    sum_distance = 0;
                    for(int index = 0 ; index < tin_point_num ; index ++)
                    {
                        //printf("pts %f\t%f\t%f\n",select_pts_tar[level][index].m_X,select_pts_tar[level][index].m_Y,select_pts_tar[level][index].m_Z);
                        if(select_pts_tar[index].flag && transformed_coord[index].flag)
                        {
                            long tininfo_col = (long)((select_pts_tar[index].m_X - Boundary[0])/level_ref_dx + 0.5);
                            long tininfo_row = (long)((Boundary[3] - select_pts_tar[index].m_Y)/level_ref_dx + 0.5);
                            long tininfo_pos = tininfo_row*GridSize.width + tininfo_col;
                            
                            long ref_col = (long)((select_pts_tar[index].m_X - ImageBoundary_ref[0])/level_ref_dx + 0.5);
                            long ref_row = (long)((ImageBoundary_ref[3] - select_pts_tar[index].m_Y)/level_ref_dx + 0.5);
                            long ref_pos = ref_row*ref_dem_size.width + ref_col;
                            
                            long tar_col = (long)((select_pts_tar[index].m_X - ImageBoundary_tar[0])/level_ref_dx + 0.5);
                            long tar_row = (long)((ImageBoundary_tar[3] - select_pts_tar[index].m_Y)/level_ref_dx + 0.5);
                            long tar_pos = tar_row*tar_dem_size.width + tar_col;
                            
                            if(ref_col >= 0 && ref_col < ref_dem_size.width && ref_row >= 0 && ref_row < ref_dem_size.height &&
                               tar_col >= 0 && tar_col < tar_dem_size.width && tar_row >= 0 && tar_row < tar_dem_size.height)
                            {
                                if(DEM[ref_pos] > -100 && DEM_tar[tar_pos] > -100 &&
                                   tininfo_col >= 0 && tininfo_col < GridSize.width && tininfo_row >= 0 && tininfo_row < GridSize.height)
                                {
                                    if(while_iter == 0)
                                        normalized_pt[index] = Normalize_coord(select_pts_tar[index],coord_center,coord_scale);
                                    
                                    transformed_coord[index] = ConformalTransform(normalized_pt[index],conparam);
                                    
                                    
                                    if(while_iter == 0)
                                    {
                                        F3DPOINT ref_normal_ori;
                                        ref_normal[index] = FindNormal(&ref_normal_ori,DEM, normalized_pt[index],coord_center,coord_scale, ImageBoundary_ref, conparam, ref_dem_size, level_ref_dx, ref_array[index], &ref_iter_height[index], 0);
                                    }
                                    if(ref_normal[index].m_X > -999)
                                    {
                                        double* tar_array = (double*)calloc(sizeof(double),9);
                                        float tar_iter_height;
                                        
                                        tar_normal[index] = FindNormal(&tar_normal_ori[index],DEM_tar, transformed_coord[index],coord_center,coord_scale, ImageBoundary_tar, conparam, tar_dem_size, level_ref_dx, tar_array, &tar_iter_height, 1);
                                        
                                        if(tar_normal[index].m_X > -999)
                                        {
                                            float ref_slope, ref_aspect, tar_slope, tar_aspect;
                                            SlopeAspect(ref_normal[index],coord_scale,&ref_slope,&ref_aspect);
                                            SlopeAspect(tar_normal[index],coord_scale,&tar_slope,&tar_aspect);
                                            double ncc = Correlate(ref_array[index],tar_array,9);
                                            
                                            if (ncc != -99)
                                                ncc = (ncc + 1)/2.0;
                                            else
                                                ncc = 0.0;
                                            
                                            float ration_slope = 1 - fabs(ref_slope - tar_slope)/90.0;
                                            float ration_aspect = 1 - fabs(ref_aspect - tar_aspect)/360.0;
                                            
                                            float W = (40*ncc + 40*ration_slope + 20*ration_aspect)/100.0;
                                            
                                            int t_W = floor(W*100);
                                            if(t_W >= 0 && t_W < 100)
                                                hist_W[t_W]++;
                                                       
                                            if(while_iter < 4)
                                                weight[index] = W;
                                            
                                            bool check_W_pos = false;
                                            if(W > W_th && ref_slope > 0 && tar_slope > 0 /*&& ref_aspect > 0 && tar_aspect > 0*/ &&
                                            tar_iter_height > min_nor_coord.m_Z && tar_iter_height < max_nor_coord.m_Z &&
                                            transformed_coord[index].m_X > min_nor_coord.m_X && transformed_coord[index].m_X < max_nor_coord.m_X &&
                                            transformed_coord[index].m_Y > min_nor_coord.m_Y && transformed_coord[index].m_Y < max_nor_coord.m_Y &&
                                            transformed_coord[index].m_Z > min_nor_coord.m_Z && transformed_coord[index].m_Z < max_nor_coord.m_Z)
                                                check_W_pos = true;
                                            
                                            //if(while_iter >= 4)
                                            //    check_W_pos = true;
                                            
                                            if(check_W_pos)
                                            {
                                                F3DPOINT t_coord = normalized_pt[index];
                                                t_coord.m_Z = tar_iter_height;
                                                
                                                F3DPOINT ref_pts = SurfaceDistance_ori(tar_normal_ori[index],DEM,tar_normal[index], t_coord, ImageBoundary_ref,ref_dem_size, level_ref_dx, conparam, coord_center,coord_scale, ref_iter_height[index]);
                                                
                                                float D = -(tar_normal[index].m_X*t_coord.m_X + tar_normal[index].m_Y*t_coord.m_Y + tar_normal[index].m_Z*t_coord.m_Z);
                                                float diff_distance = -(tar_normal[index].m_X*ref_pts.m_X + tar_normal[index].m_Y*ref_pts.m_Y + tar_normal[index].m_Z*ref_pts.m_Z + D);
                                                
                                                float dist_x = (t_coord.m_Z - ref_pts.m_Z)*coord_scale.m_Z;
                                                
                                                if(fabs(diff_distance) < th_dH /*&& fabs(dist_x) < th_dH_vertical*/)
                                                {
                                                    //printf("ID %d\tW %f\tnormal %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",index,W,ref_normal.m_X,ref_normal.m_Y,ref_normal.m_Z,tar_normal[index].m_X,tar_normal[index].m_Y,tar_normal[index].m_Z,ref_slope,ref_aspect,tar_slope,tar_aspect,ncc);
                                                    //printf("ref pts %f\t%f\t%f\tD %f\tdiff %f\t",ref_pts.m_X,ref_pts.m_Y,ref_pts.m_Z,D,diff_distance);
                                                    
                                                    SDistance[index] = diff_distance;
                                                    number_of_selected++;
                                                    
                                                    sum_distance += diff_distance;
                                                    control_pts[index] = ref_pts;
                                                    /*
                                                    if(min_H > ref_pts.m_Z * coord_scale.m_Z + coord_center.m_Z)
                                                        min_H = ref_pts.m_Z * coord_scale.m_Z + coord_center.m_Z;
                                                    if(max_H < ref_pts.m_Z * coord_scale.m_Z + coord_center.m_Z)
                                                        max_H = ref_pts.m_Z * coord_scale.m_Z + coord_center.m_Z;
                                                    
                                                    if(min_X > ref_pts.m_X * coord_scale.m_X + coord_center.m_X)
                                                        min_X = ref_pts.m_X * coord_scale.m_X + coord_center.m_X;
                                                    if(max_X < ref_pts.m_X * coord_scale.m_X + coord_center.m_X)
                                                        max_X = ref_pts.m_X * coord_scale.m_X + coord_center.m_X;
                                                    
                                                    if(min_Y > ref_pts.m_Y * coord_scale.m_Y + coord_center.m_Y)
                                                        min_Y = ref_pts.m_Y * coord_scale.m_Y + coord_center.m_Y;
                                                    if(max_Y < ref_pts.m_Y * coord_scale.m_Y + coord_center.m_Y)
                                                        max_Y = ref_pts.m_Y * coord_scale.m_Y + coord_center.m_Y;
                                                     */
                                                    //if(while_iter < 4)
                                                    {
                                                        select_pts_tar[index].flag = true;
                                                        transformed_coord[index].flag = true;
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                if(while_iter < 4)
                                                {
                                                    select_pts_tar[index].flag = false;
                                                    transformed_coord[index].flag = false;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            if(while_iter < 4)
                                            {
                                                select_pts_tar[index].flag = false;
                                                transformed_coord[index].flag = false;
                                            }
                                        }
                                        free(tar_array);
                                    }
                                    else
                                    {
                                        if(while_iter < 4)
                                        {
                                            select_pts_tar[index].flag = false;
                                            transformed_coord[index].flag = false;
                                        }
                                    }
                                }
                                else
                                {
                                    if(while_iter < 4)
                                    {
                                        select_pts_tar[index].flag = false;
                                        transformed_coord[index].flag = false;
                                    }
                                }
                            }
                            else
                            {
                                if(while_iter < 4)
                                {
                                    select_pts_tar[index].flag = false;
                                    transformed_coord[index].flag = false;
                                }
                            }
                        }
                        else
                        {
                            if(while_iter < 4)
                            {
                                select_pts_tar[index].flag = false;
                                transformed_coord[index].flag = false;
                            }
                        }
                    }
                    
                    double ration_pts = (double)number_of_selected/(double)tin_point_num*100.0;
                    
                    if(number_of_selected < 10)
                    {
                        check_stop = true;
                        check_null_dem = true;
                    }
                    final_number_of_selected = number_of_selected;
                    if(!check_null_dem)
                    {
                        //dX = AdjustmentConformal3D(tin_point_num, number_of_selected, normalized_pts, select_pts_tar, coord_center, coord_scale, SDistance, param, tininfo_tar, GridSize, Boundary, level_tar_dx, weight);
                        total_ET_ad = time(0);
                        total_gap = difftime(total_ET_ad,total_ST_ad);
                        printf("\npre matrix adjustment time %f\t%d\n\n",total_gap,number_of_selected);
                        total_ST_ad = time(0);
                        
                        dX = CoeffMatrix_25D(coord_center,coord_scale, tin_point_num, number_of_selected, transformed_coord, SDistance, conparam, tar_normal,weight,sigmaX,&sigma0);
                        
                        total_ET_ad = time(0);
                        total_gap = difftime(total_ET_ad,total_ST_ad);
                        printf("\nmatrix adjustment time %f\n\n",total_gap);
                        total_ST_ad = time(0);
                        
                        conparam.Tx += dX[4];
                        conparam.Ty += dX[5];
                        conparam.Tz += dX[6];
                        
                        double change_ratio_sigma = fabs(pre_sigma - sigma0*coord_scale.m_Z)/pre_sigma;
                        
                        if((fabs(dX[4]*coord_scale.m_X) < 0.01 && fabs(dX[5]*coord_scale.m_Y) < 0.01 && fabs(dX[6]*coord_scale.m_Z) < 0.01) && change_ratio_sigma < 0.01)
                            check_stop = true;
                        printf("pre sigma %f\t%f\t%f\n",pre_sigma,sigma0*coord_scale.m_Z,change_ratio_sigma);
                        pre_sigma = sigma0*coord_scale.m_Z;
                        //free(dX);
                        
                        average_distance = sum_distance/final_number_of_selected;
                        
                        float* save_Z = (float*)calloc(sizeof(float),final_number_of_selected);
                        float* save_Z_vertical = (float*)calloc(sizeof(float),final_number_of_selected);
                        float sum_var = 0;
                        long count = 0;
                        long count_vertical = 0;
                        char dem_gcp_filename[500];
                        sprintf(dem_gcp_filename,"%s/DEM_gcps_%d.txt",args.Outputpath,ti);
                        FILE *fid_dem_gcp = fopen(dem_gcp_filename,"w");
                        float sum_distance_vertical = 0;
                        
                        for(long i=0;i<tin_point_num;i++)
                        {
                            if(select_pts_tar[i].flag && transformed_coord[i].flag)
                            {
                                save_Z[count] = SDistance[i];
                                
                                sum_var += (average_distance - save_Z[count])*(average_distance - save_Z[count]);
                                fprintf(fid_dem_gcp,"%f\t%f\t%f\n",control_pts[i].m_X*coord_scale.m_X + coord_center.m_X,
                                        control_pts[i].m_Y*coord_scale.m_Y + coord_center.m_Y,
                                        control_pts[i].m_Z*coord_scale.m_Z + coord_center.m_Z);
                                count++;
                                
                                D2DPOINT gcp_coord,tar_img, ref_img,gcp_coord_ref;
                                gcp_coord.m_X = control_pts[i].m_X*coord_scale.m_X + coord_center.m_X + conparam.Tx*coord_scale.m_X;
                                gcp_coord.m_Y = control_pts[i].m_Y*coord_scale.m_Y + coord_center.m_Y + conparam.Ty*coord_scale.m_Y;
                                
                                tar_img.m_X = ( gcp_coord.m_X - tar_minX )/tar_dx;
                                tar_img.m_Y = ( tar_maxY - gcp_coord.m_Y )/tar_dx;
                                
                                gcp_coord_ref.m_X = control_pts[i].m_X*coord_scale.m_X + coord_center.m_X;
                                gcp_coord_ref.m_Y = control_pts[i].m_Y*coord_scale.m_Y + coord_center.m_Y;
                                
                                ref_img.m_X = ( gcp_coord_ref.m_X - ref_minX )/ref_dx;
                                ref_img.m_Y = ( ref_maxY - gcp_coord_ref.m_Y )/ref_dx;
                                
                                long tar_index = (int)(tar_img.m_Y)*tar_dem_size.width + (int)(tar_img.m_X);
                                long ref_index = (int)(ref_img.m_Y)*ref_dem_size.width + (int)(ref_img.m_X);
                                
                                if(tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_dem_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_dem_size.height &&
                                   ref_img.m_X >= 0 && ref_img.m_X < ref_dem_size.width && ref_img.m_Y >= 0 && ref_img.m_Y < ref_dem_size.height)
                                {
                                    if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                                    {
                                        long t_col_int   = (long)(tar_img.m_X + 0.01);
                                        long t_row_int   = (long)(tar_img.m_Y + 0.01);
                                        
                                        double dcol        = tar_img.m_X - t_col_int;
                                        double drow        = tar_img.m_Y - t_row_int;
                                        
                                        long index1,index2,index3, index4;
                                        double value1, value2, value3, value4, value;
                                        
                                        index1  = (t_col_int   ) + (t_row_int   )*(long)tar_dem_size.width;
                                        index2  = (t_col_int +1) + (t_row_int   )*(long)tar_dem_size.width;
                                        index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_dem_size.width;
                                        index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_dem_size.width;
                                        
                                        value1      = DEM_tar[index1];
                                        value2      = DEM_tar[index2];
                                        value3      = DEM_tar[index3];
                                        value4      = DEM_tar[index4];
                                        
                                        value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                        + value3*(1-dcol)*drow + value4*dcol*drow;
                                        
                                        //double value       = DEM[ti][tar_index];
                                        float co_dem_t = value + conparam.Tz*coord_scale.m_Z;
                                        float co_dem_diff = (DEM[ref_index] - co_dem_t);
                                        
                                        save_Z_vertical[count_vertical] = co_dem_diff;
                                        sum_distance_vertical += co_dem_diff;
                                        count_vertical++;
                                    }
                                }
                           }
                                
                        }
                        fclose(fid_dem_gcp);
                        
                        bool check_wth = false;
                        int wth_iter = 99;
                        int sum_hist_count = 0;
                        while(!check_wth && wth_iter >= 0)
                        {
                            sum_hist_count += hist_W[wth_iter];
                            
                            float hist_ratio = (float)sum_hist_count/(float)tin_point_num;
                            if(hist_ratio > 0.3)
                            {
                                check_wth = true;
                                W_th = wth_iter/100.0;
                            }
                            
                            wth_iter--;
                        }
                        if(W_th < 0.8)
                            W_th = 0.8;
                        
                        SD_distance = sqrt(sum_var/final_number_of_selected);
                        th_dH = SD_distance*3.29;
                        //th_dH = SD_distance*1.96;
                        /*if(th_dH < 5.0/ori_scale_Z)
                            th_dH = 5.0/ori_scale_Z;
                        */
                        vertical_average_distance = sum_distance_vertical/count_vertical;
                        float sum_var_vertical = 0;
                        for(long i=0;i<count_vertical;i++)
                        {
                            sum_var_vertical += (vertical_average_distance - save_Z_vertical[i])*(vertical_average_distance - save_Z_vertical[i]);
                        }
                        SD_distance_vertical = sqrt(sum_var_vertical/count_vertical);
                        
                        th_dH_vertical = SD_distance_vertical*3.29;
                        
                        MED_distance = quickselect(save_Z_vertical, count_vertical, (int)(count_vertical/2.0));
                        
                        double sum_var_med = 0;
                        
                        for(int i=0;i<count_vertical;i++)
                        {
                            sum_var_med += (MED_distance -save_Z_vertical[i])*(MED_distance -save_Z_vertical[i]);
                        }
                        printf("iter %d\tparam %f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",while_iter,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4]*sigma0*coord_scale.m_X,sigmaX[5]*sigma0*coord_scale.m_Y,sigmaX[6]*sigma0*coord_scale.m_Z, number_of_selected,sigma0*coord_scale.m_Z,change_ratio_sigma,th_dH);
                        while_iter++;
                        
                        SD_z_med = sqrt(sum_var_med/count_vertical);
                        printf("\nTz average and variation  %d\t%f\t%f\t%f\t%f\t%f\n\n",count,th_dH,vertical_average_distance,MED_distance,SD_distance_vertical,SD_z_med);
                        free(save_Z);
                        free(save_Z_vertical);
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\nadjustment time %f\n\n",total_gap);
                        total_ST = time(0);
                    }
                    
                    free(hist_W);
                }
                   
                if(!check_null_dem)
                {
                    total_ET_iter = time(0);
                    double total_gap_iter = difftime(total_ET_iter,total_ST_iter);
                    printf("\npair time %f\n\n",total_gap_iter);
                    total_ET_iter = time(0);
                        
                    float* co_dem = NULL;
                    float* copoly_dem = NULL;
                    float* save_dz = NULL;
                    
                    float all_average = 0.0;
                    float all_med = 0.0;
                    float all_std = 0.0;
                    double Coreg_param[2];
                    
                    Coreg_param[0] = conparam.Tx*coord_scale.m_X;
                    Coreg_param[1] = conparam.Ty*coord_scale.m_Y;
                    
                    if(args.check_DEM_coreg_output == 2)
                    {
                        co_dem = (float*)calloc(sizeof(float),tar_dem_size.width*tar_dem_size.height);
                        copoly_dem = (float*)calloc(sizeof(float),tar_dem_size.width*tar_dem_size.height);
                        for(long int co_index = 0 ; co_index < (long)(tar_dem_size.width)*(long)(tar_dem_size.height) ; co_index++)
                        {
                            co_dem[co_index] = Nodata;
                            copoly_dem[co_index] = Nodata;
                        }
                    }
                    
                    if(args.check_DEM_coreg_output == 1 || args.check_DEM_coreg_output == 2)
                    {
                        save_dz = (float*)calloc(sizeof(float),tar_dem_size.width*tar_dem_size.height);
                    
                        float dz_sum = 0;
                        float dz_sum_med = 0;
                        long dz_count = 0;
                
                        for(long int co_index = 0 ; co_index < (long)(tar_dem_size.width)*(long)(tar_dem_size.height) ; co_index++)
                        {
                            int pts_row = (int)(floor(co_index/tar_dem_size.width));
                            int pts_col = co_index % tar_dem_size.width;
                            
                            D2DPOINT gcp_coord,tar_img, ref_img,gcp_coord_ref;
                            gcp_coord.m_X = pts_col*tar_dx + tar_minX + Coreg_param[0];
                            gcp_coord.m_Y = tar_maxY - pts_row*tar_dx + Coreg_param[1];
                            
                            tar_img.m_X = ( gcp_coord.m_X - tar_minX )/tar_dx;
                            tar_img.m_Y = ( tar_maxY - gcp_coord.m_Y )/tar_dx;
                            
                            gcp_coord_ref.m_X = pts_col*tar_dx + tar_minX;
                            gcp_coord_ref.m_Y = tar_maxY - pts_row*tar_dx;
                            
                            ref_img.m_X = ( gcp_coord_ref.m_X - ref_minX )/ref_dx;
                            ref_img.m_Y = ( ref_maxY - gcp_coord_ref.m_Y )/ref_dx;
                            
                            long tar_index = (int)(tar_img.m_Y)*tar_dem_size.width + (int)(tar_img.m_X);
                            long ref_index = (int)(ref_img.m_Y)*ref_dem_size.width + (int)(ref_img.m_X);
                            
                            if(tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_dem_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_dem_size.height &&
                               ref_img.m_X >= 0 && ref_img.m_X < ref_dem_size.width && ref_img.m_Y >= 0 && ref_img.m_Y < ref_dem_size.height)
                            {
                                if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                                {
                                    long t_col_int   = (long)(tar_img.m_X + 0.01);
                                    long t_row_int   = (long)(tar_img.m_Y + 0.01);
                                    
                                    double dcol        = tar_img.m_X - t_col_int;
                                    double drow        = tar_img.m_Y - t_row_int;
                                    
                                    long index1,index2,index3, index4;
                                    double value1, value2, value3, value4, value;
                                    
                                    index1  = (t_col_int   ) + (t_row_int   )*(long)tar_dem_size.width;
                                    index2  = (t_col_int +1) + (t_row_int   )*(long)tar_dem_size.width;
                                    index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_dem_size.width;
                                    index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_dem_size.width;
                                    
                                    value1      = DEM_tar[index1];
                                    value2      = DEM_tar[index2];
                                    value3      = DEM_tar[index3];
                                    value4      = DEM_tar[index4];
                                    
                                    value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                    + value3*(1-dcol)*drow + value4*dcol*drow;
                                    
                                    //double value       = DEM[ti][tar_index];
                                    float co_dem_t = value + conparam.Tz*coord_scale.m_Z;
                                    float co_dem_diff = (DEM[ref_index] - co_dem_t);
                                    
                                    if(args.check_DEM_coreg_output == 2)
                                    {
                                        co_dem[tar_index] = co_dem_t;
                                        copoly_dem[tar_index] = co_dem_diff;
                                    }
                                    
                                    if(fabs(co_dem_diff) < 30)
                                    {
                                        dz_sum += co_dem_diff;
                                        save_dz[dz_count] = co_dem_diff;
                                        
                                        dz_count++;
                                    }
                                }
                            }
                            
                        }
                        
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\nWhole image stat %f\n\n",total_gap);
                        total_ST = time(0);
                        
                        all_average = dz_sum/dz_count;
                        all_med = quickselect(save_dz, dz_count, (int)(dz_count/2.0));
                        float all_sum_var = 0;
                        
            #pragma omp parallel for reduction(+:all_sum_var) schedule(guided)
                        for(long index = 0 ; index < dz_count ; index++)
                        {
                            all_sum_var += (save_dz[index] - all_average)*(save_dz[index] - all_average);
                        }
                        //printf("done allsum\n");
                        free(save_dz);
                        
                        all_std = sqrt(all_sum_var/dz_count);
                        
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\nall median median time %f\n\n",total_gap);
                        total_ST = time(0);
                        
                        printf("all stat %f\t%f\t%f\n",all_average,all_med,all_std);
                        
                        if(args.check_DEM_coreg_output == 2)
                        {
                            WriteGeotiff(co_dems_name, co_dem, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            WriteGeotiff(co_dems_name_diff, copoly_dem, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                            
                            free(co_dem);
                            free(copoly_dem);
                        }
                        
                        total_ET = time(0);
                        total_gap = difftime(total_ET,total_ST);
                        printf("\noutput save time %f\n\n",total_gap);
                        total_ST = time(0);
                        
                        
                    }
                    
                    if(sigma0*coord_scale.m_Z < 2 && vertical_average_distance < 1.0 && SD_distance_vertical < 5.0 && final_number_of_selected > 100 && sigmaX[4]*sigma0*coord_scale.m_X*(sqrt(distance_scale/coord_scale.m_X)) < tar_dx/2.0 && sigmaX[5]*sigma0*coord_scale.m_Y*(sqrt(distance_scale/coord_scale.m_Y)) < tar_dx/2.0)
                    {
                        
                           fprintf(fid_out,"%s\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t\t%4.2f\n",DEM_name_outfile,sigma0*coord_scale.m_Z,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4]*sigma0*coord_scale.m_X*(sqrt(distance_scale/coord_scale.m_X)),sigmaX[5]*sigma0*coord_scale.m_Y*(sqrt(distance_scale/coord_scale.m_Y)),sigmaX[6]*sigma0*coord_scale.m_Z*(sqrt(distance_scale/coord_scale.m_Z)),vertical_average_distance,SD_distance_vertical, MED_distance,SD_z_med,all_average,all_med,all_std,final_number_of_selected,total_gap_iter);
                        
                       /* fprintf(fid_out,"%s\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t\t%4.2f\n",DEM_name_outfile,sigma0*distance_scale,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4]*sigma0*coord_scale.m_X,sigmaX[5]*sigma0*coord_scale.m_Y,sigmaX[6]*sigma0*coord_scale.m_Z,vertical_average_distance,SD_distance_vertical, MED_distance,SD_z_med,all_average,all_med,all_std,final_number_of_selected,total_gap_iter);
                        
                        fprintf(fid_out,"%s\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t\t%4.2f\n",DEM_name_outfile,sigma0*coord_scale.m_Z,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4],sigmaX[5],sigmaX[6],vertical_average_distance,SD_distance_vertical, MED_distance,SD_z_med,all_average,all_med,all_std,final_number_of_selected,total_gap_iter);
                        */
                    }
                    else
                    {
                        fprintf(fid_out,"%s\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t\t%4.2f\n",DEM_name_outfile,sigma0*coord_scale.m_Z,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4]*sigma0*coord_scale.m_X*(sqrt(distance_scale/coord_scale.m_X)),sigmaX[5]*sigma0*coord_scale.m_Y*(sqrt(distance_scale/coord_scale.m_Y)),sigmaX[6]*sigma0*coord_scale.m_Z*(sqrt(distance_scale/coord_scale.m_Z)),vertical_average_distance,SD_distance_vertical, MED_distance,SD_z_med,all_average,all_med,all_std,final_number_of_selected,total_gap_iter);
                    }
                }
                else
                {
                    fprintf(fid_out,"%s\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\t%d\t\tNaN\n",DEM_name_outfile,final_number_of_selected);
                }
                free(dX);
                free(sigmaX);
                
                
                free(control_pts);
                free(transformed_coord);
                free(SDistance);
                free(SDistance_ori);
                free(weight);
                free(tar_normal);
                free(tar_normal_ori);
                
                free(normalized_pt);
                free(ref_normal);
                free(ref_iter_height);
                for(int index = 0 ; index < tin_point_num ; index ++)
                {
                    if(select_pts_tar[index].flag)
                    {
                        free(ref_array[index]);
                    }
                }
                free(ref_array);
                
                printf("done\n");
            }
            else
            {
                fprintf(fid_out,"%s\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\t%d\t\tNaN\n",DEM_name_outfile,tin_point_num);
            }
            free(DEM_tar);
            free(data_size_tar);
        }
        fclose(fid_out);
        //fclose(fid_out_bad);
        
        free(DEM);
        free(data_size_ref);
        
    }
}

F3DPOINT FindNormal(F3DPOINT *normal_ori, float* dem, F3DPOINT Pos,F3DPOINT Mean, F3DPOINT Scale, double* Boundary, Conformalparam X, CSize tinsize, double Gridspace, double *roh_array, float *Z, bool check_tar)
{
    F3DPOINT normal_vector;
    
    float geo_X = Pos.m_X*Scale.m_X + Mean.m_X;
    float geo_Y = Pos.m_Y*Scale.m_Y + Mean.m_Y;
    //printf("findnormal %f\t%f\n",geo_X,geo_Y);
    float col_float = (geo_X- Boundary[0])/Gridspace;
    float row_float = (Boundary[3]-geo_Y)/Gridspace;
    
    int col = floor(col_float);
    int row = floor(row_float);
    
    float X_diff = col_float - col;
    float Y_diff = row_float - row;
    
    //normalized height setting
    long index11 =  row   *tinsize.width + (col    );
    long index12 =  row   *tinsize.width + (col + 1);
    long index21 = (row+1)*tinsize.width + (col    );
    long index22 = (row+1)*tinsize.width + (col + 1);
    
    long index31 = (row-1)*tinsize.width + (col - 1);
    long index32 = (row-1)*tinsize.width + (col    );
    long index33 = (row-1)*tinsize.width + (col + 1);
    
    long index41 =  row   *tinsize.width + (col - 1);
    long index51 = (row+1)*tinsize.width + (col - 1);
    
    if(index11 >= 0 && index11 < tinsize.width*tinsize.height &&
       index12 >= 0 && index12 < tinsize.width*tinsize.height &&
       index21 >= 0 && index21 < tinsize.width*tinsize.height &&
       index22 >= 0 && index22 < tinsize.width*tinsize.height &&
       index31 >= 0 && index31 < tinsize.width*tinsize.height &&
       index32 >= 0 && index32 < tinsize.width*tinsize.height &&
       index33 >= 0 && index33 < tinsize.width*tinsize.height &&
       index41 >= 0 && index41 < tinsize.width*tinsize.height &&
       index51 >= 0 && index51 < tinsize.width*tinsize.height && row + 1 < tinsize.height && row - 1 >= 0 && col + 1 < tinsize.width && col -1 >= 0)
    {
        float Patch11,Patch12,Patch21,Patch22,Patch31,Patch32,Patch33,Patch41,Patch51;
        F3DPOINT P1,P2,P3,P4,P5,P6,P7,P8,P9;
        
        Patch11 = (dem[index11] - Mean.m_Z)/Scale.m_Z;
        Patch12 = (dem[index12] - Mean.m_Z)/Scale.m_Z;
        Patch21 = (dem[index21] - Mean.m_Z)/Scale.m_Z;
        Patch22 = (dem[index22] - Mean.m_Z)/Scale.m_Z;
        Patch31 = (dem[index31] - Mean.m_Z)/Scale.m_Z;
        Patch32 = (dem[index32] - Mean.m_Z)/Scale.m_Z;
        Patch33 = (dem[index33] - Mean.m_Z)/Scale.m_Z;
        Patch41 = (dem[index41] - Mean.m_Z)/Scale.m_Z;
        Patch51 = (dem[index51] - Mean.m_Z)/Scale.m_Z;
        
        if(dem[index11] > -100 && dem[index12] > -100 && dem[index21] > -100 && dem[index22] > -100 &&
           dem[index31] > -100 && dem[index32] > -100 && dem[index33] > -100 && dem[index41] > -100 && dem[index51] > -100)
        {
            P1.m_X =(    (col  )*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P1.m_Y =( - ((row  )*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P1.m_Z = Patch11;
             
            P2.m_X =(    (col+1)*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P2.m_Y =( - ((row  )*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P2.m_Z = Patch12;
             
            P3.m_X =(    (col+1)*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P3.m_Y =( - ((row+1)*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P3.m_Z = Patch22;
             
            P4.m_X =(    (col  )*Gridspace + Boundary[0]  - Mean.m_X)/Scale.m_X;
            P4.m_Y =( - ((row+1)*Gridspace - Boundary[3]) - Mean.m_Y)/Scale.m_Y;
            P4.m_Z = Patch21;
            
            P5.m_X =(    (col-1)*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P5.m_Y =( - ((row-1)*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P5.m_Z = Patch31;
        
            P6.m_X =(    (col  )*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P6.m_Y =( - ((row-1)*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P6.m_Z = Patch32;
        
            P8.m_X =(    (col+1)*Gridspace + Boundary[0]  - Mean.m_X)/Scale.m_X;
            P8.m_Y =( - ((row-1)*Gridspace - Boundary[3]) - Mean.m_Y)/Scale.m_Y;
            P8.m_Z = Patch33;
            
            P7.m_X =(    (col-1)*Gridspace + Boundary[0]   - Mean.m_X)/Scale.m_X;
            P7.m_Y =( - ((row  )*Gridspace - Boundary[3])  - Mean.m_Y)/Scale.m_Y;
            P7.m_Z = Patch41;
        
            P9.m_X =(    (col-1)*Gridspace + Boundary[0]  - Mean.m_X)/Scale.m_X;
            P9.m_Y =( - ((row+1)*Gridspace - Boundary[3]) - Mean.m_Y)/Scale.m_Y;
            P9.m_Z = Patch51;
            
            if(check_tar)
            {
                RM R = MakeRotationMatrix(X.omega, X.phi, X.kappa);
                
                float R11 = R.m11;
                float R12 = R.m12;
                float R13 = R.m13;
                
                float R21 = R.m21;
                float R22 = R.m22;
                float R23 = R.m23;
                
                float R31 = R.m31;
                float R32 = R.m32;
                float R33 = R.m33;
                
                P1.m_Z = X.scale*( R13*P1.m_X + R23*P1.m_Y + R33*P1.m_Z ) + X.Tz;
                P2.m_Z = X.scale*( R13*P2.m_X + R23*P2.m_Y + R33*P2.m_Z ) + X.Tz;
                P3.m_Z = X.scale*( R13*P3.m_X + R23*P3.m_Y + R33*P3.m_Z ) + X.Tz;
                P4.m_Z = X.scale*( R13*P4.m_X + R23*P4.m_Y + R33*P4.m_Z ) + X.Tz;
                P5.m_Z = X.scale*( R13*P5.m_X + R23*P5.m_Y + R33*P5.m_Z ) + X.Tz;
                P6.m_Z = X.scale*( R13*P6.m_X + R23*P6.m_Y + R33*P6.m_Z ) + X.Tz;
                P7.m_Z = X.scale*( R13*P7.m_X + R23*P7.m_Y + R33*P7.m_Z ) + X.Tz;
                P8.m_Z = X.scale*( R13*P8.m_X + R23*P8.m_Y + R33*P8.m_Z ) + X.Tz;
                P9.m_Z = X.scale*( R13*P9.m_X + R23*P9.m_Y + R33*P9.m_Z ) + X.Tz;
            }
            
            roh_array[0] = P1.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[1] = P2.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[2] = P3.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[3] = P4.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[4] = P5.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[5] = P6.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[6] = P7.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[7] = P8.m_Z*Scale.m_Z + Mean.m_Z;
            roh_array[8] = P9.m_Z*Scale.m_Z + Mean.m_Z;
            
            //for(int i = 0 ; i< 9 ; i++)
            //    printf("id %d\tdem val %f\n",i,roh_array[i]);
            
            //vector calculation
            F3DPOINT v1,v2,n1,n2;
            v1.m_X = P3.m_X - P1.m_X;
            v1.m_Y = P3.m_Y - P1.m_Y;
            v1.m_Z = P3.m_Z - P1.m_Z;
            
            v2.m_X = P2.m_X - P1.m_X;
            v2.m_Y = P2.m_Y - P1.m_Y;
            v2.m_Z = P2.m_Z - P1.m_Z;
            
            n1.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n1.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n1.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            v1.m_X = P4.m_X - P1.m_X;
            v1.m_Y = P4.m_Y - P1.m_Y;
            v1.m_Z = P4.m_Z - P1.m_Z;
            
            v2.m_X = P3.m_X - P1.m_X;
            v2.m_Y = P3.m_Y - P1.m_Y;
            v2.m_Z = P3.m_Z - P1.m_Z;
            
            n2.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n2.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n2.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            //plane vector decide
            F3DPOINT n;
            if( (X_diff == 0 && Y_diff == 0) || (X_diff != 0 & Y_diff == 0) || ( X_diff > 0 & Y_diff > 0 & (X_diff > Y_diff)) )
                n = n1;
            else
                n = n2;
               
            float mag = sqrt(n.m_X*n.m_X + n.m_Y*n.m_Y + n.m_Z*n.m_Z);
            n.m_X /= mag;
            n.m_Y /= mag;
            n.m_Z /= mag;
            
            float d = - P1.m_X*n.m_X - P1.m_Y*n.m_Y - P1.m_Z*n.m_Z;
            *Z = - (Pos.m_X*n.m_X + Pos.m_Y*n.m_Y + d)/n.m_Z;
            
            normal_vector = n;
            
            //denormailzed vector
            F3DPOINT P1_ori,P2_ori,P3_ori,P4_ori;
            P1_ori = Denormalize_coord(P1,Mean,Scale);
            P2_ori = Denormalize_coord(P2,Mean,Scale);
            P3_ori = Denormalize_coord(P3,Mean,Scale);
            P4_ori = Denormalize_coord(P4,Mean,Scale);
            
            P1 = P1_ori;
            P2 = P2_ori;
            P3 = P3_ori;
            P4 = P4_ori;
            
            v1.m_X = P3.m_X - P1.m_X;
            v1.m_Y = P3.m_Y - P1.m_Y;
            v1.m_Z = P3.m_Z - P1.m_Z;
            
            v2.m_X = P2.m_X - P1.m_X;
            v2.m_Y = P2.m_Y - P1.m_Y;
            v2.m_Z = P2.m_Z - P1.m_Z;
            
            n1.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n1.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n1.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            v1.m_X = P4.m_X - P1.m_X;
            v1.m_Y = P4.m_Y - P1.m_Y;
            v1.m_Z = P4.m_Z - P1.m_Z;
            
            v2.m_X = P3.m_X - P1.m_X;
            v2.m_Y = P3.m_Y - P1.m_Y;
            v2.m_Z = P3.m_Z - P1.m_Z;
            
            n2.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n2.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n2.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            //plane vector decide
            if( (X_diff == 0 && Y_diff == 0) || (X_diff != 0 & Y_diff == 0) || ( X_diff > 0 & Y_diff > 0 & (X_diff > Y_diff)) )
                n = n1;
            else
                n = n2;
               
            mag = sqrt(n.m_X*n.m_X + n.m_Y*n.m_Y + n.m_Z*n.m_Z);
            n.m_X /= mag;
            n.m_Y /= mag;
            n.m_Z /= mag;
            
            //float d = - P1.m_X*n.m_X - P1.m_Y*n.m_Y - P1.m_Z*n.m_Z;
            //*Z = - (Pos.m_X*n.m_X + Pos.m_Y*n.m_Y + d)/n.m_Z;
            
            normal_ori->m_X = n.m_X;
            normal_ori->m_Y = n.m_Y;
            normal_ori->m_Z = n.m_Z;
        }
        else
        {
            normal_vector.m_X = Nodata;
            normal_vector.m_Y = Nodata;
            normal_vector.m_Z = Nodata;
        }
    }
    else
    {
        normal_vector.m_X = Nodata;
        normal_vector.m_Y = Nodata;
        normal_vector.m_Z = Nodata;
    }
    
    return normal_vector;
}

F3DPOINT SurfaceDistance_ori(F3DPOINT tar_normal_ori, float* ref_dem, F3DPOINT tar_normal, F3DPOINT tar_pts, double *tin_boundary,CSize tinsize, double Gridspace, Conformalparam param, F3DPOINT Mean, F3DPOINT Scale, float p_ref_z)
{
    F3DPOINT ref_pts;
    
    //denomalizing normal vector
    F3DPOINT tar_normal_denor;
    tar_normal_denor.m_X    = tar_normal_ori.m_X;//tar_normal.m_X/Scale.m_X;
    tar_normal_denor.m_Y    = tar_normal_ori.m_Y;//tar_normal.m_Y/Scale.m_Y;
    tar_normal_denor.m_Z    = tar_normal_ori.m_Z;//tar_normal.m_Z/Scale.m_Z;
    
    float mag              = sqrt(tar_normal_denor.m_X*tar_normal_denor.m_X + tar_normal_denor.m_Y*tar_normal_denor.m_Y + tar_normal_denor.m_Z*tar_normal_denor.m_Z);
    tar_normal_denor.m_X = tar_normal_denor.m_X/mag;
    tar_normal_denor.m_Y = tar_normal_denor.m_Y/mag;
    tar_normal_denor.m_Z = tar_normal_denor.m_Z/mag;
    
    //vector slope about x and y directions
    float vector_slope_x      = atan2(fabs(tar_normal_denor.m_Z),fabs(tar_normal_denor.m_X));
    float vector_slope_y      = atan2(fabs(tar_normal_denor.m_Z),fabs(tar_normal_denor.m_Y));
        
    //set up for start points
    F3DPOINT line_node_pt = Denormalize_coord(tar_pts,Mean,Scale);
    F3DPOINT pre_pt_on_surface   = line_node_pt;
    pre_pt_on_surface.m_Z = p_ref_z*Scale.m_Z + Mean.m_Z;
    
    F3DPOINT after_pt_on_surface = pre_pt_on_surface;
    int max_iter            = 10;
    int count = 1;
    bool check_stop = false;
    while(count < max_iter && !check_stop)
    {
        F3DPOINT line_node_pt_after;
        line_node_pt_after.m_X = line_node_pt.m_X + (pre_pt_on_surface.m_Z - line_node_pt.m_Z)*(tar_normal_denor.m_X/tar_normal_denor.m_Z);
        line_node_pt_after.m_Y = line_node_pt.m_Y + (pre_pt_on_surface.m_Z - line_node_pt.m_Z)*(tar_normal_denor.m_Y/tar_normal_denor.m_Z);
        line_node_pt_after.m_Z = pre_pt_on_surface.m_Z;
        
        F3DPOINT find_pt    = Normalize_coord(line_node_pt_after,Mean,Scale);
        Conformalparam temp_X;
        temp_X.scale = 1.0;
        temp_X.omega = 0.0;
        temp_X.phi = 0.0;
        temp_X.kappa = 0.0;
        temp_X.Tx = 0.0;
        temp_X.Ty = 0.0;
        temp_X.Tz = 0.0;
        
        double temp_array[9];
        float Z_ref;
        F3DPOINT normal_ori;
        F3DPOINT temp = FindNormal(&normal_ori,ref_dem, find_pt, Mean, Scale, tin_boundary, temp_X, tinsize, Gridspace, temp_array, &Z_ref, false);
        
        Z_ref                           = Z_ref*Scale.m_Z + Mean.m_Z;
    
        F3DPOINT after_pt_on_surface;
        after_pt_on_surface.m_X  = line_node_pt_after.m_X;
        after_pt_on_surface.m_Y  = line_node_pt_after.m_Y;
        after_pt_on_surface.m_Z  = Z_ref;
        
        float dx                        = pre_pt_on_surface.m_X - after_pt_on_surface.m_X;
        float dy                        = pre_pt_on_surface.m_Y - after_pt_on_surface.m_Y;
        float dz                        = pre_pt_on_surface.m_Z - after_pt_on_surface.m_Z;
        float diff                      = sqrt(dx*dx + dy*dy + dz*dz);
        
        if(diff < 0.1)
        {
            ref_pts.m_X                 = (after_pt_on_surface.m_X - Mean.m_X) / Scale.m_X;
            ref_pts.m_Y                 = (after_pt_on_surface.m_Y - Mean.m_Y) / Scale.m_Y;
            ref_pts.m_Z                 = (after_pt_on_surface.m_Z - Mean.m_Z) / Scale.m_Z;
            check_stop = true;
        }
        else
        {
            float pre_angle_x               = atan2(fabs(dz),fabs(dx));
            float pre_angle_y               = atan2(fabs(dz),fabs(dy));
            bool divergent_index_x          = (pre_angle_x >= vector_slope_x & tar_normal_denor.m_X != 0 & dx != 0);
            bool divergent_index_y          = (pre_angle_y >= vector_slope_y & tar_normal_denor.m_Y != 0 & dy != 0);
            if( divergent_index_x || divergent_index_y)
            {
                F3DPOINT temp_pt;
                temp_pt.m_X                  = pre_pt_on_surface.m_X + (after_pt_on_surface.m_X - pre_pt_on_surface.m_X)/2.0;
                temp_pt.m_Y                  = pre_pt_on_surface.m_Y + (after_pt_on_surface.m_Y - pre_pt_on_surface.m_Y)/2.0;
                
                find_pt                     = Normalize_coord(temp_pt,Mean,Scale);
                float temp_Z = 0;
                temp = FindNormal(&normal_ori,ref_dem, find_pt, Mean, Scale, tin_boundary, temp_X, tinsize, Gridspace, temp_array, &temp_Z, false);
                temp_Z                   = temp_Z*Scale.m_Z + Mean.m_Z;
                after_pt_on_surface.m_X  = temp_pt.m_X;
                after_pt_on_surface.m_Y  = temp_pt.m_Y;
                after_pt_on_surface.m_Z  = temp_Z;
            }
        }
        
        pre_pt_on_surface               = after_pt_on_surface;
        
        count++;
    }
    
    if(!check_stop)
    {
        ref_pts.m_X = Nodata;
        ref_pts.m_Y = Nodata;
        ref_pts.m_Z = Nodata;
    }
    
    return ref_pts;
}

void SlopeAspect(F3DPOINT normal, F3DPOINT scale, float *slope, float *aspect)
{
    normal.m_X = normal.m_X/(scale.m_X);
    normal.m_Y = normal.m_Y/(scale.m_Y);
    normal.m_Z = normal.m_Z/(scale.m_Z);
    
    float denominator = sqrt(normal.m_X*normal.m_X + normal.m_Y*normal.m_Y + normal.m_Z*normal.m_Z);
    float value = normal.m_Z/denominator;
    
    *slope = acos(value)*(180.0/PI);
    
    *aspect = 90 - atan2(normal.m_Y,normal.m_X)*(180.0/PI);
    if(*aspect < 0)
        *aspect = *aspect + 360;
}

F3DPOINT ConformalTransform(F3DPOINT input, Conformalparam param)
{
    F3DPOINT out;
    
    RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);
       
    float S = param.scale;
    //3D conformal transformation
    out.m_X = S*(R.m11*input.m_X + R.m21*input.m_Y + R.m31*input.m_Z) + param.Tx;
    out.m_Y = S*(R.m12*input.m_X + R.m22*input.m_Y + R.m32*input.m_Z) + param.Ty;
    out.m_Z = S*(R.m13*input.m_X + R.m23*input.m_Y + R.m33*input.m_Z) + param.Tz;
    
    return out;
}

F3DPOINT Normalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale)
{
    F3DPOINT out;
    
    out.m_X = (input.m_X - Mean.m_X) / Scale.m_X;
    out.m_Y = (input.m_Y - Mean.m_Y) / Scale.m_Y;
    out.m_Z = (input.m_Z - Mean.m_Z) / Scale.m_Z;
    
    return out;
}

F3DPOINT Denormalize_coord(F3DPOINT input, F3DPOINT Mean, F3DPOINT Scale)
{
    F3DPOINT out;
    
    out.m_X = input.m_X * Scale.m_X + Mean.m_X;
    out.m_Y = input.m_Y * Scale.m_Y + Mean.m_Y;
    out.m_Z = input.m_Z * Scale.m_Z + Mean.m_Z;
    
    return out;
}

unsigned char* CreateHillshade(float* _input, CSize _img_size, double grid_size)
{
    CSize result_size;
    unsigned char* result_img;
    result_img = (unsigned char*)calloc(sizeof(unsigned char),(long)_img_size.height*(long)_img_size.width);
    
    double SobleFilter_X[3][3] = { {-1, 0 , 1} , {-2, 0 , 2}, {-1, 0 , 1} };
    double SobleFilter_Y[3][3] = { {1, 2 , 1} , {0, 0 , 0}, {-1, -2 , -1} };
    double a = 0;
    double b = 1.0/sqrt(2.0);
    double alpha = 45.0;
    double beta = 315.0;
    
    //printf("hillshade %d\t%d\t%f\n",_img_size.height,_img_size.width,grid_size);
     
#pragma omp parallel for schedule(guided)
    //for(long int ori_index=0 ; ori_index<(long)_img_size.height*(long)_img_size.width ; ori_index++)
    for(long int r=0;r<_img_size.height;r++)
    {
        for(long int c=0;c<_img_size.width;c++)
        {
            //long r = (int)(floor(ori_index/_img_size.width));
            //long c = ori_index % _img_size.width;
            long ori_index = r*_img_size.width + c;
            if(_input[ori_index] > -100)
            {
                float maxdx = 0;
                float maxdy = 0;
                bool check_null = false;
                for(int row = -1; row <= 1; row++)
                {
                    for(int col = -1 ; col <= 1 ;col++)
                    {
                        long row_index = r + row;
                        long col_index = c + col;
                        if(col_index >= 0 && col_index < _img_size.width && row_index >= 0 && row_index < _img_size.height)
                        {
                            long index = row_index*_img_size.width + col_index;
                            maxdx += SobleFilter_X[row+1][col+1]*_input[index];
                            maxdy += SobleFilter_Y[row+1][col+1]*_input[index];
                            
                            if(_input[index] < -100)
                                check_null = true;
                        }
                    }
                }
                if(!check_null)
                {
                    maxdx = maxdx/(8.0*grid_size);
                    maxdy = maxdy/(8.0*grid_size);
                    
                    double p0 = -cos(beta*DegToRad)*tan(alpha*DegToRad);
                    double q0 = -sin(beta*DegToRad)*tan(alpha*DegToRad);
                    
                    double p  = (p0*maxdx + q0*maxdy)/sqrt(p0*p0 + q0*q0);
                    int V  = (int)( ((0.5 + 0.5* ((p+a)/b)) ) * 255 );
                     
                    if(V > 255)
                        V = 255;
                    if(V < 0)
                        V = 0;
                    result_img[ori_index] = V;
                     
                }
                else
                    result_img[ori_index] = 0;
            }
            else
                result_img[ori_index] = 0;
        }
    }
    
    /*
    unsigned char* result_img_edge;

    result_img_edge = (unsigned char*)calloc(sizeof(unsigned char),(long)_img_size.height*(long)_img_size.width);
    double SharpFilter[3][3] = { {-1, -1 , -1} , {-1, 16 , -1}, {-1, -1 , -1} };

    
    for(long int r=0;r<_img_size.height;r++)
    {
        for(long int c=0;c<_img_size.width;c++)
        {
            long ori_index = r*_img_size.width + c;
            if(result_img[ori_index] > 0)
            {
                float sum = 0;
                bool check_null = false;
                for(int row = -1; row <= 1; row++)
                {
                    for(int col = -1 ; col <= 1 ;col++)
                    {
                        long row_index = r + row;
                        long col_index = c + col;
                        if(col_index >= 0 && col_index < _img_size.width && row_index >= 0 && row_index < _img_size.height)
                        {
                            long index = row_index*_img_size.width + col_index;
                            sum += SharpFilter[row+1][col+1]*(float)(result_img[index]);
                            //printf("row col %d\t%d\t%f\t%d\t%f\t%f\n",row,col,SharpFilter[row+1][col+1],result_img[index],SharpFilter[row+1][col+1]*result_img[index],sum);
                            if(result_img[index] == 0)
                                check_null = true;
                        }
                    }
                }
            
                //printf("sum %f\n",sum);
                if(!check_null)
                {
                    int V = int(sum/8.0);
                    if(V > 255)
                        V = 255;
                    if(V < 0)
                        V = 0;
                    
                    result_img_edge[ori_index] = V;
                }
            }
        }
    }
            
    free(result_img);
    */
    return result_img;
}

F3DPOINT* SettingControls(float* DEM_ref, float* DEM_tar, double grid_size_ref, double grid_size_tar, double *boundary_ref, double *boundary_tar, double* overlapped_br, CSize img_size_ref, CSize img_size_tar, long *tin_point_num)
{
    *tin_point_num = 0;
    F3DPOINT *coord = NULL;
    
    CSize coord_size;
    double control_spacing = grid_size_tar*20;
    control_spacing = 50;
    
    printf("control_spacing %f\n",control_spacing);
    
    coord_size.width = ceil((overlapped_br[2] - overlapped_br[0])/control_spacing);
    coord_size.height = ceil((overlapped_br[3] - overlapped_br[1])/control_spacing);
    long data_length = (long)coord_size.width*(long)coord_size.height;
    
    coord = (F3DPOINT*)calloc(sizeof(F3DPOINT),data_length);
    printf("coord_size %d\t%d\t%d\n",coord_size.width,coord_size.height,data_length);
    int kernel_size;
    
    if(grid_size_tar >= 2)
        kernel_size = 1;
    else if(grid_size_tar == 1)
        kernel_size = 2;
    else
        kernel_size = 4;
    
    long count_select = 0;
    for(long r=0;r<coord_size.height;r++)
    {
        for(long c=0;c<coord_size.width;c++)
        {
            long control_index = r*coord_size.width + c;
            
            double t_coord_x = overlapped_br[0] + c*control_spacing;
            double t_coord_y = overlapped_br[3] - r*control_spacing;
            
            long pos_c_ref = (t_coord_x - boundary_ref[0])/grid_size_ref;
            long pos_r_ref = (boundary_ref[3] - t_coord_y)/grid_size_ref;
            
            long pos_c_tar = (t_coord_x - boundary_tar[0])/grid_size_tar;
            long pos_r_tar = (boundary_tar[3] - t_coord_y)/grid_size_tar;
            
            if(pos_c_ref - kernel_size >= 0 && pos_c_ref + kernel_size < img_size_ref.width && pos_r_ref - kernel_size >= 0 && pos_r_ref + kernel_size < img_size_ref.height &&
               pos_c_tar - kernel_size >= 0 && pos_c_tar + kernel_size < img_size_tar.width && pos_r_tar - kernel_size >= 0 && pos_r_tar + kernel_size < img_size_tar.height)
            {
                long ref_index = pos_r_ref*img_size_ref.width + pos_c_ref;
                long tar_index = pos_r_tar*img_size_tar.width + pos_c_tar;
                
                double *ref_array = (double*)calloc(sizeof(double),(2*kernel_size + 1)*(2*kernel_size + 1));
                double *tar_array = (double*)calloc(sizeof(double),(2*kernel_size + 1)*(2*kernel_size + 1));
                
                for(int k = -kernel_size ; k <= kernel_size ; k++)
                {
                    for(int l = -kernel_size ; l <= kernel_size ; l++)
                    {
                        long array_index = (k+kernel_size)*(2*kernel_size + 1) + (l+kernel_size);
                        long k_ref_index = (pos_r_ref + k)*img_size_ref.width + pos_c_ref + l;
                        long k_tar_index = (pos_r_tar + k)*img_size_tar.width + pos_c_tar + l;
                        ref_array[array_index] = DEM_ref[k_ref_index];
                        tar_array[array_index] = DEM_tar[k_tar_index];
                        
                    }
                }
                
                double ncc = Correlate(ref_array,tar_array,(2*kernel_size + 1)*(2*kernel_size + 1));
                free(ref_array);
                free(tar_array);
                
                if(ncc > 0.7 && DEM_ref[ref_index] > -100 && DEM_ref[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                {
                    coord[control_index].m_X = t_coord_x;
                    coord[control_index].m_Y = t_coord_y;
                    coord[control_index].m_Z = DEM_tar[tar_index];
                    coord[control_index].flag = true;
                    
                    count_select++;
                    //printf("count_select %d\n",count_select);
                    //printf("inside ID %d\t%d\t sel pts %f\t%f\t%f\n",r,c,coord[control_index].m_X,coord[control_index].m_Y,coord[control_index].m_Z);
                }
                else
                    coord[control_index].flag = false;
            }
        }
    }
    
    *tin_point_num = count_select;
    printf("SettingControls %d\n",count_select);
    F3DPOINT *select_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),count_select);
    
    long t_count = 0;
    for(long r=0;r<coord_size.height;r++)
    {
        for(long c=0;c<coord_size.width;c++)
        {
            long ori_index = r*coord_size.width + c;
            if(coord[ori_index].flag)
            {
                select_pts[t_count] = coord[ori_index];
                
                //printf("ID %d sel pts %f\t%f\t%f\n",t_count,select_pts[t_count].m_X,select_pts[t_count].m_Y,select_pts[t_count].m_Z);
                t_count++;
                
                
            }
        }
    }
    
    free(coord);

    return select_pts;
}

F3DPOINT* CreateImagePyramid_DEM(float* _input, double grid_size, double *boundary, double* overlapped_br, uint8 pyramid_level, CSize _img_size, int _filter_size, double _sigma, float* result_img, long *tin_point_num, bool check_pts)
{
    double sigma = _sigma;
    double temp,scale;
    double sum = 0;
    double** GaussianFilter;
    CSize result_size;
    *tin_point_num = 0;
    //float* result_img;
    F3DPOINT *coord = NULL;
    
    if(pyramid_level > 0)
    {
        GaussianFilter = (double**)malloc(sizeof(double*)*_filter_size);
        for(int i=0;i<_filter_size;i++)
            GaussianFilter[i] = (double*)malloc(sizeof(double)*_filter_size);
        
        
        result_size.width = _img_size.width/2;
        result_size.height = _img_size.height/2;
        scale=sqrt(2*PI)*sigma;
        
        if(check_pts)
        {
            coord = (F3DPOINT*)calloc(sizeof(F3DPOINT),(long)result_size.height*(long)result_size.width);
            if (coord == NULL)
            {
                printf("ERROR: Out of memory - coord is NULL\n");
                exit(1);
            }
        }
        
        for(int i=-(int)(_filter_size/2);i<(int)(_filter_size/2)+1;i++)
            {
                for(int j=-(int)(_filter_size/2);j<(int)(_filter_size/2)+1;j++)
                {
                    temp = -1*(i*i+j*j)/(2*sigma*sigma);
                    GaussianFilter[i+(int)(_filter_size/2)][j+(int)(_filter_size/2)]=exp(temp)/scale;
                    sum += exp(temp)/scale;
                }
            }
            
        #pragma omp parallel for schedule(guided)
            for(int i=-(int)(_filter_size/2);i<(int)(_filter_size/2)+1;i++)
            {
                for(int j=-(int)(_filter_size/2);j<(int)(_filter_size/2)+1;j++)
                {
                    GaussianFilter[i+(int)(_filter_size/2)][j+(int)(_filter_size/2)]/=sum;
                    //int GI = (int)((GaussianFilter[i+(int)(_filter_size/2)][j+(int)(_filter_size/2)] + 0.001)*1000);
                    //GaussianFilter[i+(int)(_filter_size/2)][j+(int)(_filter_size/2)] = (double)(GI/1000.0);
                }
            }
    }
    else
    {
        result_size.width = _img_size.width;
        result_size.height = _img_size.height;
        long data_length = (long)result_size.height*(long)result_size.width;
        //result_img = (float*)malloc(sizeof(float)*data_length);
        memcpy(result_img,_input,sizeof(float)*data_length);
        
        if(check_pts)
        {
            coord = (F3DPOINT*)calloc(sizeof(F3DPOINT),data_length);
            printf("Done coord allocation %ld\n",data_length);
            
            if (coord == NULL)
            {
                printf("ERROR: Out of memory - coord is NULL\n");
                exit(1);
            }
        }
        
        printf("level 0 memory cp %d\n",data_length);
    }
    
    printf("py size %d\t%d\n",result_size.width,result_size.height);
//#pragma omp parallel for private(temp) schedule(guided)
    for(long r=0;r<result_size.height;r++)
    {
        for(long c=0;c<result_size.width;c++)
        {
            long data_length = (long)result_size.height*(long)result_size.width;
            temp = 0;
            long ori_index = r*(long)result_size.width + c;
            
            double t_coord_x = boundary[0] + c*grid_size;
            double t_coord_y = boundary[3] - r*grid_size;
            bool check_inside = false;
            
            if(t_coord_x >= overlapped_br[0] && t_coord_x <= overlapped_br[2] && t_coord_y >= overlapped_br[1] && t_coord_y <= overlapped_br[3])
                check_inside = true;
            
            if(check_pts && check_inside)
            {
                coord[ori_index].m_X = t_coord_x;
                coord[ori_index].m_Y = t_coord_y;
                
                //printf("coord %f\t%f\n",coord[ori_index].m_X ,coord[ori_index].m_Y);
            }
            
            if(r >= 0 && r < result_size.height && c >= 0 && c < result_size.width && ori_index >= 0 && ori_index < data_length)
            {
                if(pyramid_level > 0)
                {
                    long tt_index = 2*r*_img_size.width + 2*c;
                    if(2*r < _img_size.height && 2*r >= 0 && 2*c < _img_size.width && 2*c >=0)
                    {
                        if(_input[tt_index] > -100 && _input[tt_index] < 10000)
                        {
                            bool check_NULL = false;
                            for(long l=0;l<_filter_size;l++)
                            {
                                for(long k=0;k<_filter_size;k++)
                                {
                                    long temp_pos_r = (2*r + l-(int)(_filter_size/2));
                                    long temp_pos_c = (2*c + k-(int)(_filter_size/2));
                                    long temp_index = temp_pos_r*_img_size.width +temp_pos_c;
                                    //r'->2r+m, c'->2c+n
                                    if( temp_pos_r >= 0 && temp_pos_c >= 0 &&
                                       temp_pos_r < _img_size.height && temp_pos_c < _img_size.width)
                                    {
                                        if(_input[temp_index] > -100 && _input[temp_index] < 10000)
                                            temp += GaussianFilter[l][k]*_input[temp_index];
                                        else
                                            check_NULL = true;
                                    }
                                }
                            }
                            
                            
                            if(!check_NULL)
                            {
                                result_img[ori_index] = round(temp);
                                if(check_pts && check_inside)
                                {
                                    coord[ori_index].m_Z = result_img[ori_index];
                                    coord[ori_index].flag = true;
                                    (*tin_point_num)++;
                                }
                                else if(check_pts)
                                    coord[ori_index].flag = false;
                            }
                            else
                            {
                                result_img[ori_index] = -100;
                                if(check_pts)
                                    coord[ori_index].flag = false;
                            }
                        }
                        else
                        {
                            result_img[ori_index] = -100;
                            if(check_pts)
                                coord[ori_index].flag = false;
                        }
                        
                    }
                    else
                    {
                        result_img[ori_index] = -100;
                        if(check_pts)
                            coord[ori_index].flag = false;
                    }
                }
                else
                {
                    if(check_pts && check_inside)
                    {
                        //printf("ID %d\t%f\t%f\n",ori_index,coord[ori_index].m_X,coord[ori_index].m_Y);
                        if(result_img[ori_index] > -100 && result_img[ori_index] < 10000)
                        {
                            coord[ori_index].m_Z = result_img[ori_index];
                            coord[ori_index].flag = true;
                            (*tin_point_num)++;
                        }
                        else
                        {
                            coord[ori_index].flag = false;
                        }
                    }
                    else if(check_pts)
                        coord[ori_index].flag = false;
                }
            }
            
        }
    }
    
    
    printf("Done pyramid level %d\t%d\n",pyramid_level,*tin_point_num);
    
    F3DPOINT *select_pts = NULL;
    
    if(check_pts)
    {
        //char save[500] = "/home/noh.56/development/SETSM_osu_github/setsm/tin_points.txt";
        //FILE *fid = fopen(save,"w");
        select_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),*tin_point_num);
        printf("done select_pts allocation\n");
        
        long t_count = 0;
        for(long r=0;r<result_size.height;r++)
        {
            for(long c=0;c<result_size.width;c++)
            {
                long ori_index = r*result_size.width + c;
                if(coord[ori_index].flag)
                {
                    select_pts[t_count] = coord[ori_index];
                    //fprintf(fid,"%5.0f\t%5.0f\t%5.2f\n",select_pts[t_count].m_X,select_pts[t_count].m_Y,select_pts[t_count].m_Z);
                    //printf("ID %d sel pts %f\t%f\t%f\n",t_count,select_pts[t_count].m_X,select_pts[t_count].m_Y,select_pts[t_count].m_Z);
                    t_count++;
                }
            }
        }
        free(coord);
        //fclose(fid);
    }
    
    if(pyramid_level > 0)
    {
        for(int i=0;i<_filter_size;i++)
            if(GaussianFilter[i])
                free(GaussianFilter[i]);
        
        if(GaussianFilter)
            free(GaussianFilter);
    }
    
    return select_pts;
}

void SetHeightRange_slope_aspect(float* ref_img, double* ref_br, CSize ref_size, float* tar_img, double* tar_br, CSize tar_size, long numOfPts, long numOfTri, F3DPOINT *pts, UI3DPOINT *tris, double *boundary, CSize input_size, double input_grid,TINinfo* tininfo, bool check_tar)
{
    uint32 num_triangles;
    int i, tcnt;
    double gridspace = input_grid;
    CSize gridsize;
    uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
    int Col_C, Row_R;
    
    num_triangles           = numOfTri;
    boundary    = boundary;
    gridsize.width  = input_size.width;
    gridsize.height = input_size.height;
    TIN_Grid_Size_X= gridsize.width;
    TIN_Grid_Size_Y= gridsize.height;
    
//  #pragma omp parallel shared(GridPT3,pts,tris,num_triangles,m_bHeight,pyramid_step,iteration,gridspace,boundary,gridsize,TIN_Grid_Size_X,TIN_Grid_Size_Y,DEM_error) private(tcnt)
    {
//      #pragma omp for
        for(tcnt=0;tcnt<(int)(num_triangles);tcnt++)
        {
            UI3DPOINT t_tri;
            F3DPOINT pt0,pt1,pt2;
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
            
            t_tri   = tris[tcnt];
            pdex0 = t_tri.m_X;
            pdex1 = t_tri.m_Y;
            pdex2 = t_tri.m_Z;
            
            if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
            {
                int node1_index, node2_index, node3_index;
                bool check_anchor_flag = false;
                pt0 = pts[pdex0];
                pt1 = pts[pdex1];
                pt2 = pts[pdex2];
                
                TriP1[0]        = pt0.m_X;
                TriP2[0]        = pt1.m_X;
                TriP3[0]        = pt2.m_X;
                
                TriP1[1]        = pt0.m_Y;
                TriP2[1]        = pt1.m_Y;
                TriP3[1]        = pt2.m_Y;
                
                TriP1[2]        = pt0.m_Z;
                TriP2[2]        = pt1.m_Z;
                TriP3[2]        = pt2.m_Z;
                
                double U[3], V[3], N[3], normalized_N[3];
                double norm, angle, aspect;
                
                U[0] = pt1.m_X - pt0.m_X;
                V[0] = pt2.m_X - pt0.m_X;
                U[1] = pt1.m_Y - pt0.m_Y;
                V[1] = pt2.m_Y - pt0.m_Y;
                U[2] = pt1.m_Z - pt0.m_Z;
                V[2] = pt2.m_Z - pt0.m_Z;
                
                N[0]        =   U[1]*V[2] - V[1]*U[2];
                N[1]        = -(U[0]*V[2] - V[0]*U[2]);
                N[2]        =   U[0]*V[1] - V[0]*U[1];
                /*
                if(N[2] < 0)
                {
                    N[0] = -N[0];
                    N[1] = -N[1];
                    N[2] = -N[2];
                }
                */
                norm  = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
                
                normalized_N[0] = N[0]/norm;
                normalized_N[1] = N[1]/norm;
                normalized_N[2] = N[2]/norm;
                
                angle = acos(fabs(N[2])/norm)*180/3.141592;
                
                if(angle <= 0 && angle >= -90)
                    angle = fabs(angle);
                else if(angle <= -270 && angle >= -360)
                    angle = 360 + angle;
                else if(angle >= 270 && angle <= 360)
                    angle = 360 - angle;
                
                aspect = 90 - atan2(N[1],N[0])*180/3.141592;
                if(aspect < 0)
                    aspect = aspect + 360;
                
                // calculation on BoundingBox(MinMax XY) of triangle
                TriMinXY[0] = min(min(TriP1[0],TriP2[0]),TriP3[0]);
                TriMinXY[1] = min(min(TriP1[1],TriP2[1]),TriP3[1]);
                TriMaxXY[0] = max(max(TriP1[0],TriP2[0]),TriP3[0]);
                TriMaxXY[1] = max(max(TriP1[1],TriP2[1]),TriP3[1]);
                
                PixelMinXY[0] = (int)((TriMinXY[0] - boundary[0])/gridspace + 0.5);
                PixelMinXY[1] = (int)((boundary[3] - TriMinXY[1])/gridspace + 0.5);
                
                PixelMaxXY[0] = (int)((TriMaxXY[0] - boundary[0])/gridspace + 0.5);
                PixelMaxXY[1] = (int)((boundary[3] - TriMaxXY[1])/gridspace + 0.5);
                
                PixelMinXY[0] -= 1;     PixelMinXY[1] -= 1;
                PixelMaxXY[0] += 1;     PixelMaxXY[1] += 1;
                if (PixelMaxXY[0] >= (int)(TIN_Grid_Size_X))
                    PixelMaxXY[0] =  (int)(TIN_Grid_Size_X-1);
                if (PixelMaxXY[1] >= (int)(TIN_Grid_Size_Y))
                    PixelMaxXY[1] =  (int)(TIN_Grid_Size_Y-1);
                if (PixelMinXY[0] < 0)
                    PixelMinXY[0] = 0;
                if (PixelMinXY[1] < 0)
                    PixelMinXY[1] = 0;
                
                for (Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
                {
                    for (Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
                    {
                        long Index= TIN_Grid_Size_X*Row + Col;
                        D2DPOINT coord;
                        coord.m_X = Col * gridspace + boundary[0];
                        coord.m_Y = boundary[3] - (Row * gridspace);
                        
                        long ref_row = (long)((ref_br[3] - coord.m_Y)/gridspace + 0.5);
                        long ref_col = (long)((coord.m_X - ref_br[0])/gridspace + 0.5);
                        
                        long tar_row = (long)((tar_br[3] - coord.m_Y)/gridspace + 0.5);
                        long tar_col = (long)((coord.m_X - tar_br[0])/gridspace + 0.5);
                        
                        long ref_index = ref_row*ref_size.width + ref_col;
                        long tar_index = tar_row*tar_size.width + tar_col;
                        
                        if(ref_row >= 0 && ref_row < ref_size.height && ref_col >= 0 && ref_col < ref_size.width &&
                           tar_row >= 0 && tar_row < tar_size.height && tar_col >= 0 && tar_col < tar_size.width)
                        {
                            float ref_val = ref_img[ref_index];
                            float tar_val = tar_img[tar_index];
                            
                            if(Row >= 0 && Row < TIN_Grid_Size_Y && Col >= 0 && Col < TIN_Grid_Size_X && ref_val > -100 && ref_val < 10000 && tar_val > -100 && tar_val < 10000)
                            {
                                double CurGPXY[2]={0.};
                                float Z = -1000.0;
                                bool rtn = false;
                                double _p1[3], _p2[3], _p3[3], v12[2], v1P[2];
                                double v23[2], v2P[2], v31[2], v3P[2];
                                int Sum;
                                
                                CurGPXY[0]  = (Col)*gridspace + boundary[0];
                                CurGPXY[1]  = boundary[3] - (Row)*gridspace;
                                
                                _p1[0]      = TriP1[0];
                                _p1[1]      = TriP1[1];
                                _p1[2]      = TriP1[2];
                                
                                _p2[0]      = TriP2[0];
                                _p2[1]      = TriP2[1];
                                _p2[2]      = TriP2[2];
                                
                                _p3[0]      = TriP3[0];
                                _p3[1]      = TriP3[1];
                                _p3[2]      = TriP3[2];
                                
                                v12[0]      = _p2[0]-_p1[0];
                                v12[1]      = _p2[1]-_p1[1];
                                
                                v1P[0]      = CurGPXY[0]-_p1[0];
                                v1P[1]      = CurGPXY[1]-_p1[1];
                                
                                v23[0]      = _p3[0]-_p2[0];
                                v23[1]      = _p3[1]-_p2[1];
                                
                                v2P[0]      = CurGPXY[0]-_p2[0];
                                v2P[1]      = CurGPXY[1]-_p2[1];
                                
                                v31[0]      = _p1[0]-_p3[0];
                                v31[1]      = _p1[1]-_p3[1];
                                
                                v3P[0]      = CurGPXY[0]-_p3[0];
                                v3P[1]      = CurGPXY[1]-_p3[1];
                                
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
                                    
                                    v12[0]  = _p2[0]-_p1[0];
                                    v12[1]  = _p2[1]-_p1[1];
                                    v12[2]  = _p2[2]-_p1[2];
                                    
                                    v13[0]  = _p3[0]-_p1[0];
                                    v13[1]  = _p3[1]-_p1[1];
                                    v13[2]  = _p3[2]-_p1[2];
                                    
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
                                    //GridPT3[Index].angle = (double)angle;
                                }
                                
                                if (rtn)
                                {
                                    tininfo->dem[Index] = Z;
                                    /*tininfo->slope[Index] = angle;
                                    tininfo->aspect[Index] = aspect;
                                    tininfo->normal[Index].m_X = normalized_N[0];
                                    tininfo->normal[Index].m_Y = normalized_N[1];
                                    tininfo->normal[Index].m_Z = normalized_N[2];
                                    
                                    //float mag = sqrt(tininfo->normal[Index].m_X*tininfo->normal[Index].m_X + tininfo->normal[Index].m_Y*tininfo->normal[Index].m_Y + tininfo->normal[Index].m_Z*tininfo->normal[Index].m_Z);
                                    //printf("normal %f\t%f\t%f\t%f\n",tininfo->normal[Index].m_X,tininfo->normal[Index].m_Y,tininfo->normal[Index].m_Z,mag);
                                    
                                    if(check_tar)
                                    {
                                        int kernel_count = 0;
                                        
                                        for(int kernel_r = -2 ; kernel_r <= 2 ;kernel_r ++)
                                        {
                                            for(int kernel_c = -2 ; kernel_c <= 2 ;kernel_c ++)
                                            {
                                                long pos_ref_col = ref_col + kernel_c;
                                                long pos_ref_row = ref_row + kernel_r;
                                                
                                                long pos_tar_col = tar_col + kernel_c;
                                                long pos_tar_row = tar_row + kernel_r;
                                                
                                                long pos_ref = pos_ref_row*ref_size.width  +  pos_ref_col;
                                                long pos_tar = pos_tar_row*tar_size.width  +  pos_tar_col;
                                                
                                                if(pos_ref_col >= 0 && pos_ref_col < ref_size.width && pos_ref_row >=0 && pos_ref_row < ref_size.height &&
                                                   pos_tar_col >= 0 && pos_tar_col < tar_size.width && pos_tar_row >=0 && pos_tar_row < tar_size.height)
                                                {
                                                    if(ref_img[pos_ref] > -100 && ref_img[pos_ref] < 10000 && tar_img[pos_tar] > -100 && tar_img[pos_tar] < 10000)
                                                        kernel_count++;
                                                }
                                            }
                                        }
                                        
                                        if(kernel_count > 20)
                                        {
                                            double* ref_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                                            double* tar_patch_vecs = (double*)calloc(sizeof(double),kernel_count);
                                            
                                            kernel_count = 0;
                                            
                                            for(int kernel_r = -2 ; kernel_r <= 2 ;kernel_r ++)
                                            {
                                                for(int kernel_c = -2 ; kernel_c <= 2 ;kernel_c ++)
                                                {
                                                    long pos_ref_col = ref_col + kernel_c;
                                                    long pos_ref_row = ref_row + kernel_r;
                                                    
                                                    long pos_tar_col = tar_col + kernel_c;
                                                    long pos_tar_row = tar_row + kernel_r;
                                                    
                                                    long pos_ref = pos_ref_row*ref_size.width  +  pos_ref_col;
                                                    long pos_tar = pos_tar_row*tar_size.width  +  pos_tar_col;
                                                    
                                                    if(pos_ref_col >= 0 && pos_ref_col < ref_size.width && pos_ref_row >=0 && pos_ref_row < ref_size.height &&
                                                       pos_tar_col >= 0 && pos_tar_col < tar_size.width && pos_tar_row >=0 && pos_tar_row < tar_size.height)
                                                    {
                                                        if(ref_img[pos_ref] > -100 && ref_img[pos_ref] < 10000 && tar_img[pos_tar] > -100 && tar_img[pos_tar] < 10000)
                                                        {
                                                            ref_patch_vecs[kernel_count] = ref_img[pos_ref];
                                                            tar_patch_vecs[kernel_count] = tar_img[pos_tar];
                                                            kernel_count++;
                                                        }
                                                    }
                                                }
                                            }
                                            double ncc = Correlate(ref_patch_vecs,tar_patch_vecs,kernel_count);
                                            tininfo->ncc[Index] = ncc;
                                            free(ref_patch_vecs);
                                            free(tar_patch_vecs);
                                       }
                                    }
                                     */
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



D2DPOINT** CoregParam_Image_MPs(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double **ImageAdjust, NCCflag _flag,
                      uint8 Template_size, uint16 **Images, CSize **Imagesizes, double **Boundary, double *grid_dx, double *grid_dy,
                      int* grid_space,double** over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, int*mp_iter_count)
{
    int i;
    
    int reference_id = 0;
    CSize LImagesize;
    
    LImagesize.width  = Imagesizes[reference_id][Pyramid_step].width;
    LImagesize.height = Imagesizes[reference_id][Pyramid_step].height;
    
    double  left_IA[2];
    left_IA[0] = ImageAdjust[reference_id][0];
    left_IA[1] = ImageAdjust[reference_id][1];
    
    double **subA;
    double **TsubA;
    double **InverseSubA;
    
    //FILE **fid_pts = (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    //FILE **fid_stat= (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    
    FILE *fid_pts_pre;
    
    int ii,kk;
    
    subA    = (double**)malloc(9*sizeof(double*));
    TsubA   = (double**)malloc(6*sizeof(double*));
    InverseSubA = (double**)malloc(6*sizeof(double*));
    
    for(ii=0;ii<9;ii++)
    {
        subA[ii]    = (double*)malloc(6*sizeof(double));
        if(ii < 6)
        {
            TsubA[ii]       = (double*)malloc(9*sizeof(double));
            InverseSubA[ii] = (double*)malloc(6*sizeof(double));
        }
    }
    
    for(ii=0;ii<9;ii++)
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
    
    for(ii=0;ii<6;ii++)
        for(kk=0;kk<9;kk++)
            TsubA[ii][kk]       = subA[kk][ii];
    
    InverseSubA[0][0] =  0.555556; InverseSubA[0][1] =  0.000000; InverseSubA[0][2] =  0.000000; InverseSubA[0][3] = -0.333333; InverseSubA[0][4] =  0.000000; InverseSubA[0][5] = -0.333333;
    InverseSubA[1][0] =  0.000000; InverseSubA[1][1] =  0.166667; InverseSubA[1][2] =  0.000000; InverseSubA[1][3] =  0.000000; InverseSubA[1][4] =  0.000000; InverseSubA[1][5] =  0.000000;
    InverseSubA[2][0] =  0.000000; InverseSubA[2][1] =  0.000000; InverseSubA[2][2] =  0.166667; InverseSubA[2][3] =  0.000000; InverseSubA[2][4] =  0.000000; InverseSubA[2][5] =  0.000000;
    InverseSubA[3][0] = -0.333333; InverseSubA[3][1] =  0.000000; InverseSubA[3][2] =  0.000000; InverseSubA[3][3] =  0.500000; InverseSubA[3][4] =  0.000000; InverseSubA[3][5] =  0.000000;
    InverseSubA[4][0] =  0.000000; InverseSubA[4][1] =  0.000000; InverseSubA[4][2] =  0.000000; InverseSubA[4][3] =  0.000000; InverseSubA[4][4] =  0.250000; InverseSubA[4][5] =  0.000000;
    InverseSubA[5][0] = -0.333333; InverseSubA[5][1] =  0.000000; InverseSubA[5][2] =  0.000000; InverseSubA[5][3] =  0.000000; InverseSubA[5][4] =  0.000000; InverseSubA[5][5] =  0.500000;
    
    D2DPOINT **matched_MPs = (D2DPOINT**)malloc(sizeof(D2DPOINT*)*proinfo->number_of_images);
    uint8 **flag_MPs = (uint8**)malloc(sizeof(uint8*)*proinfo->number_of_images);
    
    CSize *grid_size = (CSize*)calloc(sizeof(CSize),proinfo->number_of_images);
    for(int i=0;i<proinfo->number_of_images;i++)
    {
        double GridSize_width = Boundary[i][2] - Boundary[i][0];
        double GridSize_height = Boundary[i][3] - Boundary[i][1];
        
        grid_size[i].width = floor(GridSize_width/grid_space[Pyramid_step]);
        grid_size[i].height = floor(GridSize_height/grid_space[Pyramid_step]);
        
        //printf("Grid_size %d\t%d\n",grid_size[i].width,grid_size[i].height);
    }
    
    matched_MPs[reference_id] = (D2DPOINT*)malloc(sizeof(D2DPOINT)*grid_size[0].height*grid_size[0].width);
    flag_MPs[reference_id] = (uint8*)calloc(sizeof(uint8),grid_size[0].height*grid_size[0].width);
    
    
    bool check_top_level = false;
    //if(Pyramid_step < total_level)
    //    check_top_level = true;
    
    char temp_path[500];
    long total_grid_counts = 0;
    int maximum_iter = 10;
    
    for(int ti = 1 ; ti < proinfo->number_of_images ; ti++)
    {
        
        D2DPOINT *MPs = NULL;
        if(check_top_level)
        {
            sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step+1);
            fid_pts_pre = fopen(temp_path,"r");
            float temp_v;
            while( fscanf(fid_pts_pre,"%f\t%f\t%f\t%f\n",&temp_v,&temp_v,&temp_v,&temp_v) != EOF )
                total_grid_counts++;
            
            fclose(fid_pts_pre);
            
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            int i = 0;
            fid_pts_pre = fopen(temp_path,"r");
            while( fscanf(fid_pts_pre,"%f\t%f\t%f\t%f\n",&MPs[i].m_X,&MPs[i].m_Y,&temp_v,&temp_v) != EOF )
            {
                //printf("mps id %d\t%f\t%f\n",i,MPs[i].m_X,MPs[i].m_Y);
                i++;
            }
            
            fclose(fid_pts_pre);
        }
        else
        {
            total_grid_counts = grid_size[ti].height*grid_size[ti].width;
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
#pragma omp parallel for
            for(int row = 0 ; row < grid_size[ti].height ; row ++)
            {
                for(int col = 0 ; col < grid_size[ti].width ; col ++)
                {
                    int index = row*grid_size[ti].width + col;
                    MPs[index].m_X = over_Boundary[ti][0] + col*grid_space[Pyramid_step];
                    MPs[index].m_Y = over_Boundary[ti][1] + row*grid_space[Pyramid_step];
                }
            }
        }
        
        //printf("ID %d\ttotal pts %d\n",ti,total_grid_counts);
        /*
        sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
        fid_pts[ti] = fopen(temp_path,"w");
        sprintf(temp_path,"%s/txt/CoregStat_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
        fid_stat[ti] = fopen(temp_path,"w");
        */
        if(proinfo->check_selected_image[ti])
        {
            matched_MPs[ti] = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            flag_MPs[ti] = (uint8*)calloc(sizeof(uint8),total_grid_counts);
            
            bool check_stop = false;
            *iter_count = 1;
            while(!check_stop && *iter_count < maximum_iter)
            {
                uint8   Half_template_size;
                double b_factor;
                bool flag_boundary = false;
                int i;
                int count_pts = 0;
                double sum_weight_X     = 0;
                double sum_weight_Y     = 0;
                double sum_max_roh      = 0;
                double t_sum_weight_X       = 0;
                double t_sum_weight_Y       = 0;
                double t_sum_max_roh        = 0;
                double shift_X, shift_Y;
                
                //calculation image coord from object coord by RFM in left and right image
                b_factor             = pow(2.0,(2-Pyramid_step))*2;
                Half_template_size   = (int)(Template_size/2.0);
                
                F3DPOINT* save_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),total_grid_counts);
                
#pragma omp parallel for private(i,t_sum_weight_X,t_sum_weight_Y,t_sum_max_roh) reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
                for(i = 0; i<total_grid_counts ; i++)
                {
                    CSize RImagesize;
                    D2DPOINT Startpos;
                    Startpos.m_X = 0;
                    Startpos.m_Y = 0;
                    
                    RImagesize.width  = Imagesizes[ti][Pyramid_step].width;
                    RImagesize.height = Imagesizes[ti][Pyramid_step].height;
                    
                    D2DPOINT Left_Imagecoord, Right_Imagecoord,Left_Imagecoord_p, Right_Imagecoord_p;
                    Left_Imagecoord_p.m_X = (MPs[i].m_X - Boundary[reference_id][0])/grid_dx[reference_id];
                    Left_Imagecoord_p.m_Y = (Boundary[reference_id][3] - MPs[i].m_Y)/grid_dy[reference_id];
                    
                    Right_Imagecoord_p.m_X = (MPs[i].m_X - Boundary[ti][0])/grid_dx[ti] + ImageAdjust[ti][1];
                    Right_Imagecoord_p.m_Y = (Boundary[ti][3] - MPs[i].m_Y)/grid_dy[ti] + ImageAdjust[ti][0];
                    
                    Left_Imagecoord = OriginalToPyramid_single(Left_Imagecoord_p,Startpos,Pyramid_step);
                    Right_Imagecoord = OriginalToPyramid_single(Right_Imagecoord_p,Startpos,Pyramid_step);
                    
                    //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                    
                    if(   Left_Imagecoord.m_Y  > Half_template_size*b_factor    + 10                    && Left_Imagecoord.m_X  > Half_template_size*b_factor + 10
                       && Left_Imagecoord.m_Y  < LImagesize.height - Half_template_size*b_factor - 10    && Left_Imagecoord.m_X  < LImagesize.width - Half_template_size*b_factor - 10
                       && Right_Imagecoord.m_Y > Half_template_size*b_factor + 10                    && Right_Imagecoord.m_X > Half_template_size*b_factor + 10
                       && Right_Imagecoord.m_Y < RImagesize.height - Half_template_size*b_factor - 10    && Right_Imagecoord.m_X < RImagesize.width - Half_template_size*b_factor - 10)
                    {
                        double Left_X = Left_Imagecoord.m_X;
                        double Left_Y = Left_Imagecoord.m_Y;
                        double Right_X = Right_Imagecoord.m_X;
                        double Right_Y = Right_Imagecoord.m_Y;
                        int index_l = ((int)Left_Y)*LImagesize.width + (int)Left_X;
                        int index_r = ((int)Right_Y)*RImagesize.width + (int)Right_X;
                        double ori_diff;
                        if( (index_l > 0 && index_l < LImagesize.height*LImagesize.width) && (index_r > 0 && index_r < RImagesize.height*RImagesize.width) && Images[reference_id][index_l] > 0 && Images[ti][index_r] > 0 )
                        {
                            ori_diff = 0;
                            
                            if(postNCC_ortho(Pyramid_step, ori_diff, Left_Y,  Left_X, Right_Y, Right_X,
                                             subA,TsubA,InverseSubA,Template_size,_flag,1,LImagesize,RImagesize,Images[reference_id],Images[ti],&t_sum_weight_X,&t_sum_weight_Y,&t_sum_max_roh))
                            {
                                flag_MPs[reference_id][i] = true;
                                
                                flag_MPs[ti][i] = true;
                                
                                //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                                
                                
                                //printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n",ti,MPs[i].m_X,MPs[i].m_Y, matched_MPs[reference_id][count_pts].m_X,matched_MPs[reference_id][count_pts].m_Y,matched_MPs[ti][count_pts].m_X,matched_MPs[ti][count_pts].m_Y);
                                //#pragma omp critical
                                {
                                    sum_weight_X += t_sum_weight_X;
                                    sum_weight_Y += t_sum_weight_Y;
                                    sum_max_roh  += t_sum_max_roh;
                                    
                                    save_pts[i].m_X = t_sum_weight_X/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].m_Y = t_sum_weight_Y/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].flag = 1;
                                    
                                    //printf(" save pts %f\t%f\n",save_pts[i].m_X,save_pts[i].m_Y);
                                }
                                //#pragma omp atomic
                                count_pts++;
                            }
                        }
                    }
                }
                
                
                if(count_pts > 10)
                {
                    shift_X             = sum_weight_X/sum_max_roh*pwrtwo(Pyramid_step);
                    shift_Y             = sum_weight_Y/sum_max_roh*pwrtwo(Pyramid_step);
                    
                    double sum_var_x = 0;
                    double sum_var_y = 0;
                    
                    for(int c_i = 0 ; c_i < total_grid_counts ; c_i++)
                    {
                        if(save_pts[c_i].flag)
                        {
                            sum_var_x += (shift_X - save_pts[c_i].m_X)*(shift_X - save_pts[c_i].m_X);
                            sum_var_y += (shift_Y - save_pts[c_i].m_Y)*(shift_Y - save_pts[c_i].m_Y);
                            
                            //printf(" save pts %f\t%f\n",save_pts[c_i].m_X,save_pts[c_i].m_Y);
                        }
                    }
                    
                    adjust_std[ti].m_X = sqrt(sum_var_x/count_pts);
                    adjust_std[ti].m_Y = sqrt(sum_var_y/count_pts);
                    
                    //printf("1 std %f\t%f\t%d\n",adjust_std[ti].m_X,adjust_std[ti].m_Y,count_pts);
                    
                    if(fabs(shift_Y) < 0.01 && fabs(shift_X) < 0.01)
                        check_stop = true;
                    
                    //printf("ti %d\t%d\t%f\t%f\t%f\t%f\t%d\n",ti, *iter_count,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    //fprintf(fid_stat[ti],"%d\t%d\t%f\t%f\t%f\t%f\t%d\n",*iter_count,ti,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    
                    shift_X             += ImageAdjust[ti][1];
                    shift_Y             += ImageAdjust[ti][0];
                    
                    ImageAdjust[ti][1]      = shift_X;
                    ImageAdjust[ti][0]      = shift_Y;
                    
                    
                    
                    
                    
                }
                else
                {
                    check_stop = true;
                }
                
                free(save_pts);
                
                (*iter_count)++;
                
                count_pts = 0;
                if(check_stop || *iter_count >= maximum_iter)
                {
                    for(int cc = 0 ; cc < total_grid_counts ; cc++)
                    {
                        if(flag_MPs[ti][cc])
                        {
                            matched_MPs[reference_id][cc].m_X = MPs[cc].m_X;
                            matched_MPs[reference_id][cc].m_Y = MPs[cc].m_Y;
                            
                            matched_MPs[ti][cc].m_X = MPs[cc].m_X;// + ImageAdjust[ti][1]*grid_dx[ti];
                            matched_MPs[ti][cc].m_Y = MPs[cc].m_Y;// - ImageAdjust[ti][0]*grid_dy[ti];
                            
                            //matched_MPs[ti][count_pts].m_X = Right_Imagecoord.m_X;
                            //matched_MPs[ti][count_pts].m_Y = Right_Imagecoord.m_Y;
                            //matched_MPs[reference_id][count_pts].m_X = Left_Imagecoord.m_X;
                            //matched_MPs[reference_id][count_pts].m_Y = Left_Imagecoord.m_Y;
                            
                            //fprintf(fid_pts[ti],"%8.2f\t%8.2f\t%8.2f\t%8.2f\n",matched_MPs[reference_id][cc].m_X,matched_MPs[reference_id][cc].m_Y,matched_MPs[ti][cc].m_X,matched_MPs[ti][cc].m_Y);
                            count_pts++;
                        }
                    }
                    
                    avg_rho[ti] = sum_max_roh/(double)(count_pts);
                    
                    mp_iter_count[ti] = count_pts;
                    //printf("done %d\t%d\n",count_pts,mp_iter_count[ti]);
                    //*NumofPts = count_pts;
                    
                    //fclose(fid_pts[ti]);
                    //fclose(fid_stat[ti]);
                }
            }
        }
        
        free(MPs);
    }
    
    for(ii=0;ii<9;ii++)
    {
        free(subA[ii]);
        if(ii < 6)
        {
            free(TsubA[ii]);
            free(InverseSubA[ii]);
        }
    }
    
    free(subA);
    free(TsubA);
    free(InverseSubA);
    
    
    D2DPOINT** Sel_matched_MPs = (D2DPOINT**)calloc(sizeof(D2DPOINT*),proinfo->number_of_images);
    
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        total_grid_counts = (long)grid_size[ti].height*(long)grid_size[ti].width;
        //printf("assing Mps starting %d\n",total_grid_counts);
        if(proinfo->check_selected_image[ti])
        {
            if(ti > 0)
            {
                
                Sel_matched_MPs[ti] = (D2DPOINT*)calloc(sizeof(D2DPOINT),mp_iter_count[ti]);
                
                int count = 0;
                for(int cc = 0 ; cc < total_grid_counts ; cc++)
                {
                    if(flag_MPs[ti][cc])
                    {
                        Sel_matched_MPs[ti][count] = matched_MPs[ti][cc];
                        count++;
                    }
                }
                //printf("assing Mps starting ID %d\t%d\t%d\t%d\n",ti,total_grid_counts,mp_iter_count[ti],count);
            }
            else
                Sel_matched_MPs[ti] = (D2DPOINT*)calloc(sizeof(D2DPOINT),1);
            
            free(matched_MPs[ti]);
            free(flag_MPs[ti]);
        }
    }
    //printf("assing Mps end\n");
    free(grid_size);
    free(matched_MPs);
    free(flag_MPs);
    
    return Sel_matched_MPs;
}


D2DPOINT* CoregParam_Image_MPs_stereo(ProInfo *proinfo,uint8 Pyramid_step, uint8 total_level, double *ImageAdjust, NCCflag _flag,
                                uint8 Template_size, unsigned char *Images_ref, unsigned char *Images_tar, CSize *Imagesizes_ref, CSize *Imagesizes_tar, double *Boundary_ref, double *Boundary_tar, double grid_dx_ref, double grid_dy_ref, double grid_dx_tar, double grid_dy_tar, int grid_space,double* over_Boundary, char* save_path, double* avg_rho, int* iter_count, D2DPOINT *adjust_std, int *mp_iter_count)
{
    int i;
    
    int reference_id = 0;
    CSize LImagesize;
    
    LImagesize.width  = Imagesizes_ref[Pyramid_step].width;
    LImagesize.height = Imagesizes_ref[Pyramid_step].height;
    
    double  left_IA[2];
    left_IA[0] = ImageAdjust[0];
    left_IA[1] = ImageAdjust[1];
    
    double **subA;
    double **TsubA;
    double **InverseSubA;
    
    //FILE **fid_pts = (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    //FILE **fid_stat= (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);
    
    FILE *fid_pts_pre;
    
    int ii,kk;
    
    subA    = (double**)malloc(9*sizeof(double*));
    TsubA   = (double**)malloc(6*sizeof(double*));
    InverseSubA = (double**)malloc(6*sizeof(double*));
    
    for(ii=0;ii<9;ii++)
    {
        subA[ii]    = (double*)malloc(6*sizeof(double));
        if(ii < 6)
        {
            TsubA[ii]       = (double*)malloc(9*sizeof(double));
            InverseSubA[ii] = (double*)malloc(6*sizeof(double));
        }
    }
    
    for(ii=0;ii<9;ii++)
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
    
    for(ii=0;ii<6;ii++)
        for(kk=0;kk<9;kk++)
            TsubA[ii][kk]       = subA[kk][ii];
    
    InverseSubA[0][0] =  0.555556; InverseSubA[0][1] =  0.000000; InverseSubA[0][2] =  0.000000; InverseSubA[0][3] = -0.333333; InverseSubA[0][4] =  0.000000; InverseSubA[0][5] = -0.333333;
    InverseSubA[1][0] =  0.000000; InverseSubA[1][1] =  0.166667; InverseSubA[1][2] =  0.000000; InverseSubA[1][3] =  0.000000; InverseSubA[1][4] =  0.000000; InverseSubA[1][5] =  0.000000;
    InverseSubA[2][0] =  0.000000; InverseSubA[2][1] =  0.000000; InverseSubA[2][2] =  0.166667; InverseSubA[2][3] =  0.000000; InverseSubA[2][4] =  0.000000; InverseSubA[2][5] =  0.000000;
    InverseSubA[3][0] = -0.333333; InverseSubA[3][1] =  0.000000; InverseSubA[3][2] =  0.000000; InverseSubA[3][3] =  0.500000; InverseSubA[3][4] =  0.000000; InverseSubA[3][5] =  0.000000;
    InverseSubA[4][0] =  0.000000; InverseSubA[4][1] =  0.000000; InverseSubA[4][2] =  0.000000; InverseSubA[4][3] =  0.000000; InverseSubA[4][4] =  0.250000; InverseSubA[4][5] =  0.000000;
    InverseSubA[5][0] = -0.333333; InverseSubA[5][1] =  0.000000; InverseSubA[5][2] =  0.000000; InverseSubA[5][3] =  0.000000; InverseSubA[5][4] =  0.000000; InverseSubA[5][5] =  0.500000;
    
    CSize grid_size;
    
    double GridSize_width = over_Boundary[2] - over_Boundary[0];
    double GridSize_height = over_Boundary[3] - over_Boundary[1];
    
    grid_size.width = floor(GridSize_width/grid_space);
    grid_size.height = floor(GridSize_height/grid_space);
    
    //printf("Grid_size %d\t%d\n",grid_size[i].width,grid_size[i].height);
    
    D2DPOINT *matched_MPs = NULL;
    uint8 *flag_MPs = NULL;
    
    
    bool check_top_level = false;
    //if(Pyramid_step < total_level)
    //    check_top_level = true;
    
    char temp_path[500];
    long total_grid_counts = 0;
    int maximum_iter = 20;
    
    //for(int ti = 1 ; ti < proinfo->number_of_images ; ti++)
    {
        
        D2DPOINT *MPs = NULL;
        /*if(check_top_level)
        {
            sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step+1);
            fid_pts_pre = fopen(temp_path,"r");
            float temp_v;
            while( fscanf(fid_pts_pre,"%f\t%f\t%f\t%f\n",&temp_v,&temp_v,&temp_v,&temp_v) != EOF )
                total_grid_counts++;
            
            fclose(fid_pts_pre);
            
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            int i = 0;
            fid_pts_pre = fopen(temp_path,"r");
            while( fscanf(fid_pts_pre,"%lf\t%lf\t%f\t%f\n",&MPs[i].m_X,&MPs[i].m_Y,&temp_v,&temp_v) != EOF )
            {
                //printf("mps id %d\t%f\t%f\n",i,MPs[i].m_X,MPs[i].m_Y);
                i++;
            }
            
            fclose(fid_pts_pre);
        }
        else*/
        {
            total_grid_counts = grid_size.height*grid_size.width;
            MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
#pragma omp parallel for
            for(int row = 0 ; row < grid_size.height ; row ++)
            {
                for(int col = 0 ; col < grid_size.width ; col ++)
                {
                    int index = row*grid_size.width + col;
                    MPs[index].m_X = over_Boundary[0] + col*grid_space;
                    MPs[index].m_Y = over_Boundary[1] + row*grid_space;
                }
            }
        }
        
        //printf("ID %d\ttotal pts %d\n",ti,total_grid_counts);
        /*
         sprintf(temp_path,"%s/txt/GCPs_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
         fid_pts[ti] = fopen(temp_path,"w");
         sprintf(temp_path,"%s/txt/CoregStat_Image_ID_%d_level_%d.txt",save_path,ti,Pyramid_step);
         fid_stat[ti] = fopen(temp_path,"w");
         */
        //if(proinfo->check_selected_image[ti])
        {
            matched_MPs = (D2DPOINT*)malloc(sizeof(D2DPOINT)*total_grid_counts);
            flag_MPs = (uint8*)calloc(sizeof(uint8),total_grid_counts);
            
            bool check_stop = false;
            *iter_count = 1;
            while(!check_stop && *iter_count < maximum_iter)
            {
                uint8   Half_template_size;
                double b_factor;
                bool flag_boundary = false;
                int i;
                int count_pts = 0;
                double sum_weight_X     = 0;
                double sum_weight_Y     = 0;
                double sum_max_roh      = 0;
                double t_sum_weight_X       = 0;
                double t_sum_weight_Y       = 0;
                double t_sum_max_roh        = 0;
                double shift_X, shift_Y;
                
                //calculation image coord from object coord by RFM in left and right image
                b_factor             = pow(2.0,(2-Pyramid_step))*2;
                Half_template_size   = (int)(Template_size/2.0);
                
                F3DPOINT* save_pts = (F3DPOINT*)calloc(sizeof(F3DPOINT),total_grid_counts);
                
#pragma omp parallel for private(i,t_sum_weight_X,t_sum_weight_Y,t_sum_max_roh) reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
                for(i = 0; i<total_grid_counts ; i++)
                {
                    CSize RImagesize;
                    D2DPOINT Startpos;
                    Startpos.m_X = 0;
                    Startpos.m_Y = 0;
                    
                    RImagesize.width  = Imagesizes_tar[Pyramid_step].width;
                    RImagesize.height = Imagesizes_tar[Pyramid_step].height;
                    
                    D2DPOINT Left_Imagecoord, Right_Imagecoord,Left_Imagecoord_p, Right_Imagecoord_p;
                    Left_Imagecoord_p.m_X = (MPs[i].m_X - Boundary_ref[0])/grid_dx_ref;
                    Left_Imagecoord_p.m_Y = (Boundary_ref[3] - MPs[i].m_Y)/grid_dy_ref;
                    
                    Right_Imagecoord_p.m_X = (MPs[i].m_X - Boundary_tar[0])/grid_dx_tar + ImageAdjust[1];
                    Right_Imagecoord_p.m_Y = (Boundary_tar[3] - MPs[i].m_Y)/grid_dy_tar + ImageAdjust[0];
                    
                    Left_Imagecoord = OriginalToPyramid_single(Left_Imagecoord_p,Startpos,Pyramid_step);
                    Right_Imagecoord = OriginalToPyramid_single(Right_Imagecoord_p,Startpos,Pyramid_step);
                    
                    //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                    
                    if(   Left_Imagecoord.m_Y  > Half_template_size*b_factor    + 10                    && Left_Imagecoord.m_X  > Half_template_size*b_factor + 10
                       && Left_Imagecoord.m_Y  < LImagesize.height - Half_template_size*b_factor - 10    && Left_Imagecoord.m_X  < LImagesize.width - Half_template_size*b_factor - 10
                       && Right_Imagecoord.m_Y > Half_template_size*b_factor + 10                    && Right_Imagecoord.m_X > Half_template_size*b_factor + 10
                       && Right_Imagecoord.m_Y < RImagesize.height - Half_template_size*b_factor - 10    && Right_Imagecoord.m_X < RImagesize.width - Half_template_size*b_factor - 10)
                    {
                        double Left_X = Left_Imagecoord.m_X;
                        double Left_Y = Left_Imagecoord.m_Y;
                        double Right_X = Right_Imagecoord.m_X;
                        double Right_Y = Right_Imagecoord.m_Y;
                        int index_l = ((int)Left_Y)*LImagesize.width + (int)Left_X;
                        int index_r = ((int)Right_Y)*RImagesize.width + (int)Right_X;
                        double ori_diff;
                        if( (index_l > 0 && index_l < LImagesize.height*LImagesize.width) && (index_r > 0 && index_r < RImagesize.height*RImagesize.width) && Images_ref[index_l] > 0 && Images_tar[index_r] > 0 )
                        {
                            ori_diff = 0;
                            F2DPOINT peak_pos;
                            if(postNCC_ortho_BYTE(Pyramid_step, ori_diff, Left_Y,  Left_X, Right_Y, Right_X,
                                             subA,TsubA,InverseSubA,Template_size,_flag,1,LImagesize,RImagesize,Images_ref,Images_tar,&t_sum_weight_X,&t_sum_weight_Y,&t_sum_max_roh,&peak_pos))
                            {
                                flag_MPs[i] = true;
                                
                                //printf("point %f\t%f\t%f\t%f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,Left_Imagecoord.m_X,Left_Imagecoord.m_Y,Right_Imagecoord.m_X,Right_Imagecoord.m_Y);
                                
                                
                                //printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n",ti,MPs[i].m_X,MPs[i].m_Y, matched_MPs[reference_id][count_pts].m_X,matched_MPs[reference_id][count_pts].m_Y,matched_MPs[ti][count_pts].m_X,matched_MPs[ti][count_pts].m_Y);
                                //#pragma omp critical
                                {
                                    sum_weight_X += t_sum_weight_X;
                                    sum_weight_Y += t_sum_weight_Y;
                                    sum_max_roh  += t_sum_max_roh;
                                    
                                    save_pts[i].m_X = peak_pos.m_X*pwrtwo(Pyramid_step);//t_sum_weight_X/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].m_Y = peak_pos.m_Y*pwrtwo(Pyramid_step);//t_sum_weight_Y/t_sum_max_roh*pwrtwo(Pyramid_step);
                                    save_pts[i].flag = 1;
                                    
                                    //printf(" save pts %f\t%f\n",save_pts[i].m_X,save_pts[i].m_Y);
                                }
                                //#pragma omp atomic
                                count_pts++;
                            }
                        }
                    }
                }
                
                
                if(count_pts > 10)
                {
                    shift_X             = sum_weight_X/sum_max_roh*pwrtwo(Pyramid_step);
                    shift_Y             = sum_weight_Y/sum_max_roh*pwrtwo(Pyramid_step);
                    
                    double sum_var_x = 0;
                    double sum_var_y = 0;
                    
                    for(int c_i = 0 ; c_i < total_grid_counts ; c_i++)
                    {
                        if(save_pts[c_i].flag)
                        {
                            sum_var_x += (shift_X - save_pts[c_i].m_X)*(shift_X - save_pts[c_i].m_X);
                            sum_var_y += (shift_Y - save_pts[c_i].m_Y)*(shift_Y - save_pts[c_i].m_Y);
                            
                            //printf(" save pts %f\t%f\n",save_pts[c_i].m_X,save_pts[c_i].m_Y);
                        }
                    }
                    
                    adjust_std->m_X = sqrt(sum_var_x/count_pts);
                    adjust_std->m_Y = sqrt(sum_var_y/count_pts);
                    
                    //printf("1 std %f\t%f\t%d\n",adjust_std[ti].m_X,adjust_std[ti].m_Y,count_pts);
                    
                    if(fabs(shift_Y) < 0.01 && fabs(shift_X) < 0.01)
                        check_stop = true;
                    
                    //printf("ti %d\t%d\t%f\t%f\t%f\t%f\t%d\n",ti, *iter_count,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    //fprintf(fid_stat[ti],"%d\t%d\t%f\t%f\t%f\t%f\t%d\n",*iter_count,ti,shift_X,shift_Y,ImageAdjust[ti][1],ImageAdjust[ti][0],count_pts);
                    
                    shift_X             += ImageAdjust[1];
                    shift_Y             += ImageAdjust[0];
                    
                    ImageAdjust[1]      = shift_X;
                    ImageAdjust[0]      = shift_Y;
                    
                    
                    
                    
                    
                }
                else
                {
                    check_stop = true;
                }
                
                free(save_pts);
                
                (*iter_count)++;
                
                count_pts = 0;
                if(check_stop || *iter_count >= maximum_iter)
                {
                    for(int cc = 0 ; cc < total_grid_counts ; cc++)
                    {
                        if(flag_MPs[cc])
                        {
                            matched_MPs[cc].m_X = MPs[cc].m_X;
                            matched_MPs[cc].m_Y = MPs[cc].m_Y;
                            
                            matched_MPs[cc].m_X = MPs[cc].m_X;// + ImageAdjust[ti][1]*grid_dx[ti];
                            matched_MPs[cc].m_Y = MPs[cc].m_Y;// - ImageAdjust[ti][0]*grid_dy[ti];
                            
                            //matched_MPs[ti][count_pts].m_X = Right_Imagecoord.m_X;
                            //matched_MPs[ti][count_pts].m_Y = Right_Imagecoord.m_Y;
                            //matched_MPs[reference_id][count_pts].m_X = Left_Imagecoord.m_X;
                            //matched_MPs[reference_id][count_pts].m_Y = Left_Imagecoord.m_Y;
                            
                            //fprintf(fid_pts[ti],"%8.2f\t%8.2f\t%8.2f\t%8.2f\n",matched_MPs[reference_id][cc].m_X,matched_MPs[reference_id][cc].m_Y,matched_MPs[ti][cc].m_X,matched_MPs[ti][cc].m_Y);
                            count_pts++;
                        }
                    }
                    
                    *avg_rho = sum_max_roh/(double)(count_pts);
                    
                    *mp_iter_count = count_pts;
                    //printf("done %d\t%d\n",count_pts,mp_iter_count[ti]);
                    //*NumofPts = count_pts;
                    
                    //fclose(fid_pts[ti]);
                    //fclose(fid_stat[ti]);
                }
            }
        }
        
        free(MPs);
    }
    
    for(ii=0;ii<9;ii++)
    {
        free(subA[ii]);
        if(ii < 6)
        {
            free(TsubA[ii]);
            free(InverseSubA[ii]);
        }
    }
    
    free(subA);
    free(TsubA);
    free(InverseSubA);
    
    
    D2DPOINT* Sel_matched_MPs = NULL;
    
    //for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        total_grid_counts = (long)grid_size.height*(long)grid_size.width;
        //printf("assing Mps starting %d\n",total_grid_counts);
        //if(proinfo->check_selected_image[ti])
        {
            //if(ti > 0)
            {
                
                Sel_matched_MPs = (D2DPOINT*)calloc(sizeof(D2DPOINT),(*mp_iter_count));
                
                int count = 0;
                for(int cc = 0 ; cc < total_grid_counts ; cc++)
                {
                    if(flag_MPs[cc])
                    {
                        Sel_matched_MPs[count] = matched_MPs[cc];
                        count++;
                    }
                }
                //printf("assing Mps starting ID %d\t%d\t%d\t%d\n",ti,total_grid_counts,mp_iter_count[ti],count);
            }
        }
    }
    //printf("assing Mps end\n");
    free(matched_MPs);
    free(flag_MPs);
    
    return Sel_matched_MPs;
}


bool postNCC_ortho_BYTE(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size,
                   NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, unsigned char* _leftimage, unsigned char* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, F2DPOINT *peak_pos)
{
    
    // input order: NumOfPt, LeftImage, RightImage, Template_size, center_row_left, center_col_left, center_row_right, center_col_right, ncc_weight
    bool check_pt = false;
    int count_max_roh       = 0;
    
    uint8 rotate_flag,multi_flag,multi_flag_sum,inter_flag,weight_flag;
    
    //Image info input
    uint16 L_rowsize   = leftsize.height;
    uint16 L_colsize   = leftsize.width;
    
    uint16 R_rowsize   = rightsize.height;
    uint16 R_colsize   = rightsize.width;
    
    double t_weight_X   = 0;
    double t_weight_Y   = 0;
    double t_max_roh    = 0;
    double diff_theta;
    int mask_row, mask_col;
    
    int Half_template_size,half_mask_size;
    int count = 0;
    
    Half_template_size  = (int)(Template_size/2);
    half_mask_size      = 1;
    
    rotate_flag = _flag.rotate_flag;
    multi_flag  = _flag.multi_flag;
    multi_flag_sum  = _flag.multi_flag_sum;
    inter_flag  = _flag.inter_flag;
    weight_flag = _flag.weight_flag;
    
    
    diff_theta = Ori_diff;
    
    double *result_rho  = (double*)calloc(9,sizeof(double));
    double *XX          = (double*)calloc(6,sizeof(double));
    double *ATLT        = (double*)calloc(6,sizeof(double));
    int i, j, k;
    uint8 cell_count = 0;
    
    for(j=0;j<9;j++)
        result_rho[j]       = -1.00;
    
    for(mask_row = - half_mask_size ; mask_row <= half_mask_size ; mask_row++)
    {
        for(mask_col = - half_mask_size ; mask_col <= half_mask_size ; mask_col++)
        {
            double rot_theta = 0.0;
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
            int row, col;
            int N;
            double val1, val2, de, de2, ncc_1, ncc_2, ncc_3;
            double temp_rho;
            int grid_index;
            
            if(rotate_flag == 1)
                rot_theta = (double)(diff_theta*bin_angle*PI/180.0);
            
            
            for(row = -Half_template_size; row <= Half_template_size ; row++)
            {
                for(col = -Half_template_size; col <= Half_template_size ; col++)
                {
                    double radius  = sqrt((double)(row*row + col*col));
                    if(radius <= Half_template_size-1)
                    {
                        double pos_row_left      = (Left_CR + row);
                        double pos_col_left      = (Left_CC + col);
                        
                        double temp_col        = (cos(-rot_theta)*col - sin(-rot_theta)*row);
                        double temp_row        = (sin(-rot_theta)*col + cos(-rot_theta)*row);
                        double pos_row_right     = (Right_CR + temp_row + mask_row);
                        double pos_col_right     = (Right_CC + temp_col + mask_col);
                        
                        if(pos_row_right-3 >= 0 && pos_row_right+3 < R_rowsize && pos_col_right-3 >= 0 && pos_col_right+3 < R_colsize &&
                           pos_row_left-3 >= 0 && pos_row_left+3 < L_rowsize && pos_col_left-3 >= 0 && pos_col_left+3 < L_colsize)
                        {
                            //interpolate left_patch
                            double dx = pos_col_left - (int) (pos_col_left);
                            double dy = pos_row_left - (int) (pos_row_left);
                            double left_patch;
                            double right_patch;
                            double dxdy = dx * dy;
                            long int position = (long int) (pos_col_left) + (long int) (pos_row_left) * L_colsize;
                            
                            // Appears inter_flag is always == 1
                            if (inter_flag == 1) {
                                left_patch = (double) (_leftimage[position]) * (1 - dx - dy + dxdy) + (double) (_leftimage[position + 1]) * (dx - dxdy) +
                                (double) (_leftimage[position + L_colsize]) * (dy - dxdy) + (double) (_leftimage[position + 1 + L_colsize]) * (dxdy);
                            } else {
                                left_patch = (double) (_leftimage[position]);
                            }
                            //interpolate right_patch
                            dx = pos_col_right - (int) (pos_col_right);
                            dy = pos_row_right - (int) (pos_row_right);
                            dxdy = dx * dy;
                            position = (long int) (pos_col_right) + (long int) (pos_row_right) * R_colsize;
                            
                            // Appears inter_flag is always == 1
                            if (inter_flag == 1) {
                                right_patch = (double) (_rightimage[position]) * (1 - dx - dy + dxdy) + (double) (_rightimage[position + 1]) * (dx - dxdy) +
                                (double) (_rightimage[position + R_colsize]) * (dy - dxdy) + (double) (_rightimage[position + 1 + R_colsize]) * (dxdy);
                            } else {
                                right_patch = (double) (_rightimage[position]);
                            }
                            
                            if(left_patch > 1 && right_patch > 1)
                            {
                                //end
                                Count_N[0]++;
                                
                                Sum_LR            = Sum_LR + left_patch*right_patch;
                                Sum_L             = Sum_L  + left_patch;
                                Sum_R             = Sum_R  + right_patch;
                                Sum_L2            = Sum_L2 + left_patch*left_patch;
                                Sum_R2            = Sum_R2 + right_patch*right_patch;
                                
                                if(multi_flag == 1)
                                {
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
            }
            
            if(Count_N[0] > 0)
            {
                N               = Count_N[0];
                val1          = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
                val2          = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
                de            = sqrt(val1*val2);
                de2           = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
                if( val1*val2 > 0)
                    ncc_1           = de2/de;
                else
                    ncc_1           = -1.0;
                
                if(multi_flag == 1)
                {
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
                    
                }
                
                if(multi_flag == 1)
                {
                    if(Count_N[1] > 0 && Count_N[2] > 0)
                        temp_rho      = ((ncc_1 + ncc_2 + ncc_3)/3.0);
                    else if(Count_N[1] > 0)
                        temp_rho      = ((ncc_1 + ncc_2)/2.0);
                    else if(Count_N[2] > 0)
                        temp_rho      = ((ncc_1 + ncc_3)/2.0);
                    else
                        temp_rho        = ncc_1;
                }
                else
                {
                    temp_rho      = ncc_1;
                }
                
                
                grid_index           = (mask_row+1)*3 + (mask_col+1);
                if(grid_index < 9)
                    result_rho[grid_index] = temp_rho;
                cell_count++;
            }
            
            
        }
    }
    
    if(cell_count == 9)
    {
        double demnum;
        double max_X        = 100;
        double max_Y        = 100;
        double max_roh      = 0;
        bool find_index_1   = false;
        bool find_index_2   = false;
        bool find_index     = false;
        
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
        
        demnum      = -pow(XX[4],2.0) + 4*XX[3]*XX[5];
        if(demnum > 0 && XX[3] < 0)
        {
            max_X = (- 2*XX[5]*XX[1] + XX[2]*XX[4])/demnum;
            max_Y = (- 2*XX[2]*XX[3] + XX[1]*XX[4])/demnum;
            max_roh =  XX[0]                + XX[1]*max_X           + XX[2]*max_Y
            + XX[3]*max_X*max_X + XX[4]*max_X*max_Y + XX[5]*max_Y*max_Y;
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
                t_weight_X += max_X*max_roh;
                t_weight_Y += max_Y*max_roh;
                t_max_roh  += max_roh;
                
                peak_pos->m_X = max_X;
                peak_pos->m_Y = max_Y;
                
                check_pt = true;
                //printf("max XY %f\t%f\t%f\n",max_X,max_Y,max_roh);
            }
            
            
        }
    }
    free(result_rho);
    free(ATLT);
    free(XX);
    
    if(check_pt)
    {
        *sum_weight_X   = t_weight_X;
        *sum_weight_Y   = t_weight_Y;
        *sum_max_roh    = t_max_roh;
    }
    else
    {
        *sum_weight_X   = 0;
        *sum_weight_Y   = 0;
        *sum_max_roh    = 0;
    }
    
    return check_pt;
}


float SurfaceDistance(TINinfo tininfo_ref, TINinfo tininfo_tar, F3DPOINT tar_pts, F3DPOINT tar_normal, CSize tinsize, double *tin_boundary, double gridsize,Conformalparam param, int *f_iter)
{
    float SurfaceDist;
    
    float th = 0.10;
    bool check_stop = false;
    
    float dh = Nodata;
    
    RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);

    //tar_pts.m_Z += param.Tz;
    
    //3D conformal transformation
    //float P1 = param.scale*(R.m11*tar_pts.m_X + R.m21*tar_pts.m_Y + R.m31*tar_pts.m_Z) + param.Tx;
    //float P2 = param.scale*(R.m12*tar_pts.m_X + R.m22*tar_pts.m_Y + R.m32*tar_pts.m_Z) + param.Ty;
    //float P3 = param.scale*(R.m13*tar_pts.m_X + R.m23*tar_pts.m_Y + R.m33*tar_pts.m_Z) + param.Tz;

    tar_pts.m_X = tar_pts.m_X + param.Tx;
    tar_pts.m_Y = tar_pts.m_Y + param.Ty;
    
    
    long ref_col = (long)((tar_pts.m_X - tin_boundary[0])/gridsize + 0.5);
    long ref_row = (long)((tin_boundary[3] - tar_pts.m_Y)/gridsize + 0.5);
    long ref_index = ref_row*tinsize.width + ref_col;
    
    tar_normal = tininfo_tar.normal[ref_index];
    tar_pts.m_Z = tininfo_tar.dem[ref_index] + param.Tz;
    
    float Znext = tininfo_ref.dem[ref_index];
    float angleX = atan2(fabs(tar_normal.m_Z),fabs(tar_normal.m_X));
    float angleY = atan2(fabs(tar_normal.m_Z),fabs(tar_normal.m_Y));
    
    int max_iter = 20;
    int iteration = 0;
    
    //printf("pts %f\t%f\t%f\tnormal %f\t%f\t%f\n",tar_pts.m_X,tar_pts.m_Y,tar_pts.m_Z,tar_normal.m_X,tar_normal.m_Y,tar_normal.m_Z);
    
    float Xpre = tar_pts.m_X;
    float Ypre = tar_pts.m_Y;
    while(!check_stop && iteration < max_iter)
    {
        //printf("iteration %d\n",iteration);
        float Xnext = tar_pts.m_X + (Znext - tar_pts.m_Z)*tar_normal.m_X/tar_normal.m_Z;
        float Ynext = tar_pts.m_Y + (Znext - tar_pts.m_Z)*tar_normal.m_Y/tar_normal.m_Z;
        
        long col = (long)((Xnext - tin_boundary[0])/gridsize + 0.5);
        long row = (long)((tin_boundary[3] - Ynext)/gridsize + 0.5);
        long index = row*tinsize.width + col;
        
        if(col >= 0 && col < tinsize.width && row >= 0 && row < tinsize.height)
        {
            float Znew = tininfo_ref.dem[index];
            //printf("next pts %f\t%f\t%f\n",Xnext,Ynext,Znew);
            float deltadx = Xpre - Xnext;
            float deltady = Ypre - Ynext;
            float deltadz = Znext - Znew;
            float deltaL  = sqrt(deltadx*deltadx + deltady*deltady + deltadz*deltadz);
            
            float delta_angleX = atan2(fabs(deltadz),fabs(deltadx));
            float delta_angleY = atan2(fabs(deltadz),fabs(deltady));
            
            if(deltaL < th)
            {
                check_stop = true;
                double d = -(tar_normal.m_X*tar_pts.m_X + tar_normal.m_Y*tar_pts.m_Y + tar_normal.m_Z*tar_pts.m_Z);
                double normal_mag = tar_normal.m_X*tar_normal.m_X + tar_normal.m_Y*tar_normal.m_Y + tar_normal.m_Z*tar_normal.m_Z;
                dh = (tar_normal.m_X*Xnext + tar_normal.m_Y*Ynext + tar_normal.m_Z*Znew + d) / (sqrt(normal_mag));
                
                //printf("done deltaL %f\t dh %f\n",deltaL,dh);
            }
            else if( (fabs(deltadx) > 0 && delta_angleX >= angleX) || (fabs(deltady) > 0 && delta_angleY >= angleY) )
            {
                Xnext = tar_pts.m_X + (Xnext - tar_pts.m_X)/2.0;
                Ynext = tar_pts.m_Y + (Ynext - tar_pts.m_Y)/2.0;
                
                col = (long)((Xnext - tin_boundary[0])/gridsize + 0.5);
                row = (long)((tin_boundary[3] - Ynext)/gridsize + 0.5);
                index = row*tinsize.width + col;
                
                //printf("angle compare %f\t%f\t%f\t%f\n",angleX,angleY,delta_angleX,delta_angleY);
                if(col >= 0 && col < tinsize.width && row >= 0 && row < tinsize.height)
                {
                    Znew = tininfo_ref.dem[index];
                    
                    //printf("divergenc next pts %f\t%f\t%f\n",Xnext,Ynext,Znew);
                }
                else
                    check_stop = true;
            }
                    
            Znext = Znew;
            Xpre = Xnext;
            Ypre = Ynext;
        }
        else
            check_stop = true;
        
        iteration++;
    }
    
    *f_iter = iteration;
    
    return dh;
}

float* CoeffMatrix_25D(F3DPOINT coord_center, F3DPOINT coord_scale,long pts_nums, long selected_pts, F3DPOINT* transformed_coord, float* dH, Conformalparam param, F3DPOINT* tar_normal,float *weight,float* sigmaX, float *sigma0)
{
    float* C = (float*)calloc(sizeof(float),21);
    
    RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);
    
    float A1 = param.omega*DegToRad;
    float A2 = param.phi*DegToRad;
    float A3 = param.kappa*DegToRad;
    float S = param.scale;
    
    float R11 = R.m11;
    float R12 = R.m12;
    float R13 = R.m13;
    
    float R21 = R.m21;
    float R22 = R.m22;
    float R23 = R.m23;
    
    float R31 = R.m31;
    float R32 = R.m32;
    float R33 = R.m33;
    
    GMA_double *A_matrix = GMA_double_create(selected_pts, 7);
    GMA_double *AW_matrix = GMA_double_create(selected_pts, 7);
    GMA_double *AWT_matrix = GMA_double_create(7,selected_pts);
    GMA_double *AWTA_matrix = GMA_double_create(7,7);
    GMA_double *AWTAWb_matrix = GMA_double_create(7,7);
    
    GMA_double *Wb_matrix = GMA_double_create(7,7);
    GMA_double *Qxx_matrix = GMA_double_create(7,7);
    
    GMA_double *L_matrix = GMA_double_create(selected_pts, 1);
    GMA_double *AWTL_matrix = GMA_double_create(7,1);
    
    GMA_double *Lb_matrix = GMA_double_create(7,1);
    GMA_double *WbLb_matrix = GMA_double_create(7,1);
    GMA_double *AWTLWbLb_matrix = GMA_double_create(7,1);
    
    GMA_double *X_matrix = GMA_double_create(7,1);
    GMA_double *AX_matrix = GMA_double_create(selected_pts,1);
    GMA_double *V_matrix = GMA_double_create(selected_pts,1);
    GMA_double *VT_matrix = GMA_double_create(1,selected_pts);
    GMA_double *VTV_matrix = GMA_double_create(1,1);
    
    /*
    GMA_double *A_matrix_ori = GMA_double_create(selected_pts, 7);
    GMA_double *AW_matrix_ori = GMA_double_create(selected_pts, 7);
    GMA_double *AWT_matrix_ori = GMA_double_create(7,selected_pts);
    GMA_double *AWTA_matrix_ori = GMA_double_create(7,7);
    GMA_double *Qxx_matrix_ori = GMA_double_create(7,7);
    */
    for(int i=0;i<7;i++)
    {
        for(int j=0;j<7;j++)
        {
            Wb_matrix->val[i][j] = 0.0;
        }
        Lb_matrix->val[i][0] = 0.000000001;
    }
    
    Wb_matrix->val[0][0] = 100000000000;
    Wb_matrix->val[1][1] = 100000000000;
    Wb_matrix->val[2][2] = 100000000000;
    Wb_matrix->val[3][3] = 100000000000;
    Wb_matrix->val[4][4] = 1.0;
    Wb_matrix->val[5][5] = 1.0;
    Wb_matrix->val[6][6] = 1.0;
    
    
    long count = 0;
    
    for(long i = 0 ; i < pts_nums ; i++)
    {
        if(transformed_coord[i].flag)
        {
            float P1 = transformed_coord[i].m_X;
            float P2 = transformed_coord[i].m_Y;
            float P3 = transformed_coord[i].m_Z;
            
            float Gx = tar_normal[i].m_X;///coord_scale.m_X;
            float Gy = tar_normal[i].m_Y;///coord_scale.m_Y;
            float Gz = -tar_normal[i].m_Z;///coord_scale.m_Z;
            
            //printf("pt.z %f\t%f\t%f\tG %f\t%f\t%f\tDh %f\n",P1,P2,P3,Gx,Gy,Gz,dH[i]);
            //X equation
            C[0] = R11*P1 + R21*P2 + R31*P3; //d_scale
            C[1] = 0.0; //d_omega
            C[2] = ( -sin(A2)*cos(A3)*P1 + sin(A2)*sin(A3)*P2 + cos(A2)*P3 )*S; //d_phi
            C[3] = ( R21*P1 - R11*P2 )*S; //d_kappa
            C[4] = 1.0;
            C[5] = 0.0;
            C[6] = 0.0;
            
            //Y equation
            C[7] = R12*P1 + R22*P2 + R32*P3;
            C[8] = (-R13*P1 - R23*P2 - R33*P3)*S;
            C[9] =( sin(A1)*cos(A2)*cos(A3)*P1 - sin(A1)*cos(A2)*sin(A3)*P2 + sin(A1)*sin(A2)*P3 )*S;
            C[10] = ( R22*P1 - R12*P2 )*S;
            C[11] = 0.0;
            C[12] = 1.0;
            C[13] = 0.0;
            
            //Z equation
            C[14] = R13*P1 + R23*P2 + R33*P3;
            C[15] = ( R12*P1 + R22*P2 + R32*P3 )*S;
            C[16] = ( -cos(A1)*cos(A2)*cos(A3)*P1 + cos(A1)*cos(A2)*sin(A3)*P2 - cos(A1)*sin(A2)*P3 )*S;
            C[17] = ( R23*P1 - R13*P2 )*S;
            C[18] = 0.0;
            C[19] = 0.0;
            C[20] = 1.0;
            
            A_matrix->val[count][0] = (Gx*C[0] + Gy*C[7]  + Gz*C[14]); //Scale
            A_matrix->val[count][1] = (Gx*C[1] + Gy*C[8]  + Gz*C[15]); //omega
            A_matrix->val[count][2] = (Gx*C[2] + Gy*C[9]  + Gz*C[16]); //phi
            A_matrix->val[count][3] = (Gx*C[3] + Gy*C[10] + Gz*C[17]); //kappa
            A_matrix->val[count][4] = (Gx*C[4] + Gy*C[11] + Gz*C[18]); //tx
            A_matrix->val[count][5] = (Gx*C[5] + Gy*C[12] + Gz*C[19]); //ty
            A_matrix->val[count][6] = (Gx*C[6] + Gy*C[13] + Gz*C[20]); //tz
            
            AW_matrix->val[count][0] = A_matrix->val[count][0]*weight[i];
            AW_matrix->val[count][1] = A_matrix->val[count][1]*weight[i];
            AW_matrix->val[count][2] = A_matrix->val[count][2]*weight[i];
            AW_matrix->val[count][3] = A_matrix->val[count][3]*weight[i];
            AW_matrix->val[count][4] = A_matrix->val[count][4]*weight[i];
            AW_matrix->val[count][5] = A_matrix->val[count][5]*weight[i];
            AW_matrix->val[count][6] = A_matrix->val[count][6]*weight[i];
            
            L_matrix->val[count][0] = dH[i];
            
            /*
            //sigmaX cal with original coords
            P1 = transformed_coord[i].m_X*coord_scale.m_X;// + coord_center.m_X;
            P2 = transformed_coord[i].m_Y*coord_scale.m_Y;// + coord_center.m_Y;
            P3 = transformed_coord[i].m_Z*coord_scale.m_Z;// + coord_center.m_Z;
            
            Gx = tar_normal[i].m_X*coord_scale.m_X;///coord_scale.m_X;
            Gy = tar_normal[i].m_Y*coord_scale.m_Y;///coord_scale.m_Y;
            Gz = -tar_normal[i].m_Z*coord_scale.m_Z;///coord_scale.m_Z;
            
            float mag = sqrt(Gx*Gx + Gy*Gy + Gz*Gz);
            Gx /= mag;
            Gy /= mag;
            Gz /= mag;
            //printf("pt.z %f\t%f\t%f\tG %f\t%f\t%f\tDh %f\n",P1,P2,P3,Gx,Gy,Gz,dH[i]);
            //X equation
            C[0] = R11*P1 + R21*P2 + R31*P3; //d_scale
            C[1] = 0.0; //d_omega
            C[2] = ( -sin(A2)*cos(A3)*P1 + sin(A2)*sin(A3)*P2 + cos(A2)*P3 )*S; //d_phi
            C[3] = ( R21*P1 - R11*P2 )*S; //d_kappa
            C[4] = 1.0;
            C[5] = 0.0;
            C[6] = 0.0;
            
            //Y equation
            C[7] = R12*P1 + R22*P2 + R32*P3;
            C[8] = (-R13*P1 - R23*P2 - R33*P3)*S;
            C[9] =( sin(A1)*cos(A2)*cos(A3)*P1 - sin(A1)*cos(A2)*sin(A3)*P2 + sin(A1)*sin(A2)*P3 )*S;
            C[10] = ( R22*P1 - R12*P2 )*S;
            C[11] = 0.0;
            C[12] = 1.0;
            C[13] = 0.0;
            
            //Z equation
            C[14] = R13*P1 + R23*P2 + R33*P3;
            C[15] = ( R12*P1 + R22*P2 + R32*P3 )*S;
            C[16] = ( -cos(A1)*cos(A2)*cos(A3)*P1 + cos(A1)*cos(A2)*sin(A3)*P2 - cos(A1)*sin(A2)*P3 )*S;
            C[17] = ( R23*P1 - R13*P2 )*S;
            C[18] = 0.0;
            C[19] = 0.0;
            C[20] = 1.0;
            
            A_matrix_ori->val[count][0] = (Gx*C[0] + Gy*C[7]  + Gz*C[14]); //Scale
            A_matrix_ori->val[count][1] = (Gx*C[1] + Gy*C[8]  + Gz*C[15]); //omega
            A_matrix_ori->val[count][2] = (Gx*C[2] + Gy*C[9]  + Gz*C[16]); //phi
            A_matrix_ori->val[count][3] = (Gx*C[3] + Gy*C[10] + Gz*C[17]); //kappa
            A_matrix_ori->val[count][4] = (Gx*C[4] + Gy*C[11] + Gz*C[18]); //tx
            A_matrix_ori->val[count][5] = (Gx*C[5] + Gy*C[12] + Gz*C[19]); //ty
            A_matrix_ori->val[count][6] = (Gx*C[6] + Gy*C[13] + Gz*C[20]); //tz
            
            AW_matrix_ori->val[count][0] = A_matrix->val[count][0]*weight[i];
            AW_matrix_ori->val[count][1] = A_matrix->val[count][1]*weight[i];
            AW_matrix_ori->val[count][2] = A_matrix->val[count][2]*weight[i];
            AW_matrix_ori->val[count][3] = A_matrix->val[count][3]*weight[i];
            AW_matrix_ori->val[count][4] = A_matrix->val[count][4]*weight[i];
            AW_matrix_ori->val[count][5] = A_matrix->val[count][5]*weight[i];
            AW_matrix_ori->val[count][6] = A_matrix->val[count][6]*weight[i];
            */
            
            count++;
        }
    }
    
    free(C);
    
    GMA_double_Tran(AW_matrix,AWT_matrix);
    GMA_double_mul(AWT_matrix,A_matrix,AWTA_matrix);
    GMA_double_sum(AWTA_matrix,Wb_matrix,AWTAWb_matrix);
        
    GMA_double_inv(AWTAWb_matrix,Qxx_matrix);
    
    GMA_double_mul(AWT_matrix,L_matrix,AWTL_matrix);
    GMA_double_mul(Wb_matrix,Lb_matrix,WbLb_matrix);
    
    GMA_double_sum(AWTL_matrix,WbLb_matrix,AWTLWbLb_matrix);
        
    
    GMA_double_mul(Qxx_matrix,AWTLWbLb_matrix,X_matrix);
    GMA_double_mul(A_matrix,X_matrix,AX_matrix);
    
    GMA_double_sub(AX_matrix,L_matrix,V_matrix);
    GMA_double_Tran(V_matrix,VT_matrix);
    GMA_double_mul(VT_matrix,V_matrix,VTV_matrix);
    
    
    /*
    GMA_double_Tran(AW_matrix_ori,AWT_matrix_ori);
    GMA_double_mul(AWT_matrix_ori,A_matrix_ori,AWTA_matrix_ori);
    GMA_double_inv(AWTA_matrix_ori,Qxx_matrix_ori);
    */
    *sigma0 = sqrt((VTV_matrix->val[0][0])/(selected_pts - 7));
    float* dx = (float*)calloc(sizeof(float),7);
    
    for(int i = 0 ; i < 7 ; i++)
    {
        sigmaX[i] = sqrt(Qxx_matrix->val[i][i]);
        dx[i] = X_matrix->val[i][0];
        
        printf("dx %d\t%f\tsigmaX %f\n",i,dx[i],sigmaX[i]);
    }
    /*
    printf("dx %f\tsigmaX %f\n",dx[4]*coord_scale.m_X,sigmaX[4]*coord_scale.m_X);
    printf("dy %f\tsigmaX %f\n",dx[5]*coord_scale.m_Y,sigmaX[5]*coord_scale.m_Y);
    printf("dz %f\tsigmaX %f\n",dx[6]*coord_scale.m_Z,sigmaX[6]*coord_scale.m_Z);
    */
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(AW_matrix);
    GMA_double_destroy(AWT_matrix);
    GMA_double_destroy(AWTA_matrix);
    GMA_double_destroy(AWTAWb_matrix);
    
    GMA_double_destroy(Wb_matrix);
    GMA_double_destroy(Qxx_matrix);
    
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AWTL_matrix);
    
    GMA_double_destroy(Lb_matrix);
    GMA_double_destroy(WbLb_matrix);
    GMA_double_destroy(AWTLWbLb_matrix);
    
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
    GMA_double_destroy(VT_matrix);
    GMA_double_destroy(VTV_matrix);
    /*
    GMA_double_destroy(A_matrix_ori);
    GMA_double_destroy(AW_matrix_ori);
    GMA_double_destroy(AWT_matrix_ori);
    GMA_double_destroy(AWTA_matrix_ori);
    GMA_double_destroy(Qxx_matrix_ori);
    */
    return dx;
}

float* AdjustmentConformal3D(long pts_nums, long selected_pts, F3DPOINT* normalized_input, F3DPOINT* input,F3DPOINT coord_center, F3DPOINT coord_scale, float* dH, Conformalparam param, TINinfo tininfo_tar, CSize tinsize, double *tin_boundary, double gridsize, float *weight)
{
    float* C = (float*)calloc(sizeof(float),21);
    
    RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);
    
    float A1 = param.omega*DegToRad;
    float A2 = param.phi*DegToRad;
    float A3 = param.kappa*DegToRad;
    float S = param.scale;
    
    float R11 = R.m11;
    float R12 = R.m12;
    float R13 = R.m13;
    
    float R21 = R.m21;
    float R22 = R.m22;
    float R23 = R.m23;
    
    float R31 = R.m31;
    float R32 = R.m32;
    float R33 = R.m33;
    
    GMA_double *A_matrix = GMA_double_create(selected_pts, 7);
    GMA_double *AW_matrix = GMA_double_create(selected_pts, 7);
    GMA_double *AWT_matrix = GMA_double_create(7,selected_pts);
    GMA_double *AWTA_matrix = GMA_double_create(7,7);
    GMA_double *AWTAWb_matrix = GMA_double_create(7,7);
    
    GMA_double *Wb_matrix = GMA_double_create(7,7);
    GMA_double *Qxx_matrix = GMA_double_create(7,7);
    
    GMA_double *L_matrix = GMA_double_create(selected_pts, 1);
    GMA_double *AWTL_matrix = GMA_double_create(7,1);
    
    GMA_double *Lb_matrix = GMA_double_create(7,1);
    GMA_double *WbLb_matrix = GMA_double_create(7,1);
    GMA_double *AWTLWbLb_matrix = GMA_double_create(7,1);
    
    GMA_double *X_matrix = GMA_double_create(7,1);
    GMA_double *AX_matrix = GMA_double_create(selected_pts,1);
    GMA_double *V_matrix = GMA_double_create(selected_pts,1);
    GMA_double *VT_matrix = GMA_double_create(1,selected_pts);
    GMA_double *VTV_matrix = GMA_double_create(1,1);
    
    for(int i=0;i<7;i++)
    {
        for(int j=0;j<7;j++)
        {
            Wb_matrix->val[i][j] = 0.0;
        }
        Lb_matrix->val[i][0] = 0.000000001;
    }
    
    Wb_matrix->val[0][0] = 100000000000;
    Wb_matrix->val[1][1] = 100000000000;
    Wb_matrix->val[2][2] = 100000000000;
    Wb_matrix->val[3][3] = 100000000000;
    Wb_matrix->val[4][4] = 1.0;
    Wb_matrix->val[5][5] = 1.0;
    Wb_matrix->val[6][6] = 1.0;
    
    
    long count = 0;
    
    for(long i = 0 ; i < pts_nums ; i++)
    {
        if(input[i].flag)
        {
            //3D conformal transformation
            //float P1 = S*(R.m11*normalized_input[i].m_X + R.m21*normalized_input[i].m_Y + R.m31*normalized_input[i].m_Z) + param.Tx;
            //float P2 = S*(R.m12*normalized_input[i].m_X + R.m22*normalized_input[i].m_Y + R.m32*normalized_input[i].m_Z) + param.Ty;
            //float P3 = S*(R.m13*normalized_input[i].m_X + R.m23*normalized_input[i].m_Y + R.m33*normalized_input[i].m_Z) + param.Tz;
            
            float P1 = (input[i].m_X + param.Tx - coord_center.m_X);
            float P2 = (input[i].m_Y + param.Ty - coord_center.m_Y);
            
            long col = (long)(((input[i].m_X + param.Tx) - tin_boundary[0])/gridsize + 0.5);
            long row = (long)((tin_boundary[3] - (input[i].m_Y + param.Ty))/gridsize + 0.5);
            long index = row*tinsize.width + col;
            
            float Gx = tininfo_tar.normal[index].m_X;///coord_scale.m_X;
            float Gy = tininfo_tar.normal[index].m_Y;///coord_scale.m_Y;
            float Gz = tininfo_tar.normal[index].m_Z;///coord_scale.m_Z;
            
            col = (long)(((input[i].m_X) - tin_boundary[0])/gridsize + 0.5);
            row = (long)((tin_boundary[3] - (input[i].m_Y))/gridsize + 0.5);
            index = row*tinsize.width + col;
            
            float P3 = tininfo_tar.dem[index] + param.Tz - coord_center.m_Z;
            
            //printf("pt.z %f\t%f\n",input[i].m_Z, tininfo_tar.dem[index]);
            //X equation
            C[0] = R11*P1 + R21*P2 + R31*P3; //d_scale
            C[1] = 0.0; //d_omega
            C[2] = ( -sin(A2)*cos(A3)*P1 + sin(A2)*sin(A3)*P2 + cos(A2)*P3 )*S; //d_phi
            C[3] = ( R21*P1 - R11*P2 )*S; //d_kappa
            C[4] = 1.0;
            C[5] = 0.0;
            C[6] = 0.0;
            
            //Y equation
            C[7] = R12*P1 + R22*P2 + R32*P3;
            C[8] = (-R13*P1 - R23*P2 - R33*P3)*S;
            C[9] =( sin(A1)*cos(A2)*cos(A3)*P1 - sin(A1)*cos(A2)*sin(A3)*P2 + sin(A1)*sin(A2)*P3 )*S;
            C[10] = ( R22*P1 - R12*P2 )*S;
            C[11] = 0.0;
            C[12] = 1.0;
            C[13] = 0.0;
            
            //Z equation
            C[14] = R13*P1 + R23*P2 + R33*P3;
            C[15] = ( R12*P1 + R22*P2 + R32*P3 )*S;
            C[16] = ( -cos(A1)*cos(A2)*cos(A3)*P1 + cos(A1)*cos(A2)*sin(A3)*P2 - cos(A1)*sin(A2)*P3 )*S;
            C[17] = ( R23*P1 - R13*P2 )*S;
            C[18] = 0.0;
            C[19] = 0.0;
            C[20] = 1.0;
            
            
            
            
            float mag = sqrt(Gx*Gx + Gy*Gy + Gz*Gz);
            Gx /= mag;
            Gy /= mag;
            Gz /= mag;
            
            A_matrix->val[count][0] = (Gx*C[0] + Gy*C[7]  + Gz*C[14]); //Scale
            A_matrix->val[count][1] = (Gx*C[1] + Gy*C[8]  + Gz*C[15]); //omega
            A_matrix->val[count][2] = (Gx*C[2] + Gy*C[9]  + Gz*C[16]); //phi
            A_matrix->val[count][3] = (Gx*C[3] + Gy*C[10] + Gz*C[17]); //kappa
            A_matrix->val[count][4] = (Gx*C[4] + Gy*C[11] + Gz*C[18]); //tx
            A_matrix->val[count][5] = (Gx*C[5] + Gy*C[12] + Gz*C[19]); //ty
            A_matrix->val[count][6] = (Gx*C[6] + Gy*C[13] + Gz*C[20]); //tz
            
            AW_matrix->val[count][0] = A_matrix->val[count][0]*weight[i];
            AW_matrix->val[count][1] = A_matrix->val[count][1]*weight[i];
            AW_matrix->val[count][2] = A_matrix->val[count][2]*weight[i];
            AW_matrix->val[count][3] = A_matrix->val[count][3]*weight[i];
            AW_matrix->val[count][4] = A_matrix->val[count][4]*weight[i];
            AW_matrix->val[count][5] = A_matrix->val[count][5]*weight[i];
            AW_matrix->val[count][6] = A_matrix->val[count][6]*weight[i];
            
            L_matrix->val[count][0] = dH[i];
            
            
            //printf("index %d\t pts %f\t%f\t%f\t%f\t%f\t%f\n",count,P1,P2,P3,normalized_input[i].m_X,normalized_input[i].m_Y,normalized_input[i].m_Z);
            //printf("Gxyz %f\t%f\t%f\t%f\t dH %f\n",Gx,Gy,Gz,mag,dH[i]);
            count++;
        }
    }
    
    free(C);
    
    GMA_double_Tran(AW_matrix,AWT_matrix);
    GMA_double_mul(AWT_matrix,A_matrix,AWTA_matrix);
    //GMA_double_sum(AWTA_matrix,Wb_matrix,AWTAWb_matrix);
        
    GMA_double_inv(AWTA_matrix,Qxx_matrix);
    //GMA_double_inv(AWTAWb_matrix,Qxx_matrix);
    
    GMA_double_mul(AWT_matrix,L_matrix,AWTL_matrix);
    //GMA_double_mul(Wb_matrix,Lb_matrix,WbLb_matrix);
    
    //GMA_double_sum(AWTL_matrix,WbLb_matrix,AWTLWbLb_matrix);
        
    
    GMA_double_mul(Qxx_matrix,AWTL_matrix,X_matrix);
    //GMA_double_mul(Qxx_matrix,AWTLWbLb_matrix,X_matrix);
    GMA_double_mul(A_matrix,X_matrix,AX_matrix);
    
    GMA_double_sub(AX_matrix,L_matrix,V_matrix);
    GMA_double_Tran(V_matrix,VT_matrix);
    GMA_double_mul(VT_matrix,V_matrix,VTV_matrix);
    
    float sigma0 = sqrt(VTV_matrix->val[0][0])/(selected_pts - 7);
    float sigmaX[7] = {0.};
    float* dx = (float*)calloc(sizeof(float),7);
    
    for(int i = 0 ; i < 7 ; i++)
    {
        sigmaX[i] = sqrt(Qxx_matrix->val[i][i]);
        dx[i] = X_matrix->val[i][0];
        
        printf("dx %d\t%f\tsigmaX %f\n",i,dx[i],sigmaX[i]);
    }
    /*
    printf("dx %f\tsigmaX %f\n",dx[4]*coord_scale.m_X,sigmaX[4]*coord_scale.m_X);
    printf("dy %f\tsigmaX %f\n",dx[5]*coord_scale.m_Y,sigmaX[5]*coord_scale.m_Y);
    printf("dz %f\tsigmaX %f\n",dx[6]*coord_scale.m_Z,sigmaX[6]*coord_scale.m_Z);
    */
    GMA_double_destroy(A_matrix);
    GMA_double_destroy(AW_matrix);
    GMA_double_destroy(AWT_matrix);
    GMA_double_destroy(AWTA_matrix);
    GMA_double_destroy(AWTAWb_matrix);
    
    GMA_double_destroy(Wb_matrix);
    GMA_double_destroy(Qxx_matrix);
    
    GMA_double_destroy(L_matrix);
    GMA_double_destroy(AWTL_matrix);
    
    GMA_double_destroy(Lb_matrix);
    GMA_double_destroy(WbLb_matrix);
    GMA_double_destroy(AWTLWbLb_matrix);
    
    GMA_double_destroy(X_matrix);
    GMA_double_destroy(AX_matrix);
    GMA_double_destroy(V_matrix);
    GMA_double_destroy(VT_matrix);
    GMA_double_destroy(VTV_matrix);
    
    return dx;
}


//End Image Coregistration
