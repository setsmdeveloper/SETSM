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
    if(args.GCP_spacing > 0)
        proinfo->GCP_spacing = args.GCP_spacing;
    else
        proinfo->GCP_spacing = -9;
    
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
            
            int status = Maketmpfolders(proinfo);
            
            uint8 *image_bits = (uint8*)malloc(sizeof(uint8)*proinfo->number_of_images);
            double *ortho_minX = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_maxY = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_dx = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double *ortho_dy = (double*)malloc(sizeof(double)*proinfo->number_of_images);
            double **Boundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
     
            long int cols[2], rows[2];
            double *GridSize_width = (double*)calloc(sizeof(double),proinfo->number_of_images);
            double *GridSize_height = (double*)calloc(sizeof(double),proinfo->number_of_images);
            
            double Sum_grid = 0;
            uint16 **OriImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
            CSize *OriImagesizes = (CSize*)malloc(sizeof(CSize)*proinfo->number_of_images);
            double **ImageBoundary = (double**)malloc(sizeof(double*)*proinfo->number_of_images);
            
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
                
                cols[0] = 0;
                cols[1] = OriImagesizes[i].width;
                
                rows[0] = 0;
                rows[1] = OriImagesizes[i].height;
                
                OriImages[i] = SubsetImageFrombitsToUint16(image_bits[i], args.Image[i], cols, rows, &OriImagesizes[i]);
                /*
                switch(image_bits[i])
                {
                    case 8:
                    {
                        uint8 type(0);
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
                        uint16 type(0);
                        uint16* data16 = Readtiff_T(args.Image[i], &OriImagesizes[i], cols, rows, &OriImagesizes[i],type);
                        long int data_size16 = (long int)OriImagesizes[i].width*(long int)OriImagesizes[i].height;
                        OriImages[i] = (uint16*)malloc(sizeof(uint16)*data_size16);
                        #pragma omp parallel for schedule(guided)
                        for(long int index = 0 ; index < data_size16 ; index++)
                            OriImages[i][index] = data16[index];
                        free(data16);
                    }
                        break;
                }
                 */
                printf("ID %d\t %s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,args.Image[i],OriImagesizes[i].width,OriImagesizes[i].height,ortho_minX[i],ortho_maxY[i],ortho_dx[i],ortho_dy[i],ImageBoundary[i][0],ImageBoundary[i][1],ImageBoundary[i][2],ImageBoundary[i][3]);
            }
            
            int *Grid_space = (int*)calloc(sizeof(int),py_level+1);
            for(int i = 0 ; i <= py_level ; i++)
            {
                Grid_space[i] = ceil(Sum_grid)*pwrtwo(py_level+2);
                if(proinfo->GCP_spacing > 0)
                    Grid_space[i] = proinfo->GCP_spacing;
                printf("image coregistration gridspace %d\t%d\n",i,Grid_space[i]);
            }
        
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
                Subsetfilename[ti] = (char*)malloc(sizeof(char)*500);
                filename = GetFileName(args.Image[ti]);
                filename = remove_ext(filename);
                sprintf(Subsetfilename[ti],"%s/%s_subset.raw",proinfo->tmpdir,filename);
                
                data_size_lr[ti] = (CSize*)malloc(sizeof(CSize)*(py_level+1));
                SetPySizes(data_size_lr[ti], OriImagesizes[ti], py_level);

                free(filename);
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
            
            sprintf(out_file,"%s/coreg_result.txt",proinfo->save_filepath);
            FILE *fid_out         = fopen(out_file,"w");
            fprintf(fid_out,"orthoimage name\tline(row) direction[pixel]\tsample(column) direction[pixel]\tTy[meter]\tTx[meter]\tavg_roh\n");
            
            TransParam param;
            SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
            int reference_id = 0;
            for(int ti = 1 ; ti < proinfo->number_of_images ; ti ++)
            {
                double avg_roh = 0;
                vector<D2DPOINT> matched_MPs;
                vector<D2DPOINT> matched_MPs_ref;
                for(int level = py_level ; level >= 0 ; level --)
                {
                    printf("Processing level %d\n",level);
                    printf("level\tImage ID\trow(pixel)\tcolumn(pixel)\tTy(meter)\tTx(meter)\tGCPS #\tavg_roh\t# of iteration\n");
                    
                    uint16 *SubImages_ref = LoadPyramidImages(proinfo->tmpdir,Subsetfilename[reference_id],data_size_lr[reference_id][level],level);
                    uint16 *SubImages_tar = LoadPyramidImages(proinfo->tmpdir,Subsetfilename[ti],data_size_lr[ti][level],level);
                    
                    int iter_counts;
                    D2DPOINT grid_dxy_ref(ortho_dx[reference_id], ortho_dy[reference_id]);
                    D2DPOINT grid_dxy_tar(ortho_dx[ti], ortho_dy[ti]);
                    
                    CoregParam_Image(proinfo, ti, level, ImageAdjust_coreg[ti], 15, SubImages_ref, data_size_lr[reference_id][level], SubImages_tar, data_size_lr[ti][level], ImageBoundary[reference_id], ImageBoundary[ti], grid_dxy_ref, grid_dxy_tar, Grid_space[level], Boundary[ti], &avg_roh, &iter_counts, &adjust_std[ti], matched_MPs, matched_MPs_ref);
                     
                    printf("%d\t%d\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%4.2f\t%d\n",level,ti,ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1],
                           -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],matched_MPs.size(),avg_roh,iter_counts);
                    free(SubImages_ref);
                    free(SubImages_tar);
                    
                    total_ET = time(0);
                    total_gap = difftime(total_ET,total_ST);
                    printf("\niter %d CoregParam_Image %f\n\n",level,total_gap);
                    total_ST = time(0);
                }
                
                total_ET = time(0);
                total_gap = difftime(total_ET,total_ST);
                printf("\nCoregParam_Image %f\n\n",total_gap);
                total_ST = time(0);
                
                
                fprintf(fid_out,"%s\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%3.2f\n",proinfo->Imagefilename[ti],ImageAdjust_coreg[ti][0], ImageAdjust_coreg[ti][1], -ImageAdjust_coreg[ti][0]*ortho_dy[ti], ImageAdjust_coreg[ti][1]*ortho_dx[ti],avg_roh);
                
                
                FILE* p_GCP = NULL;
                sprintf(out_file,"%s/txt/GCPs_Image_ID_%d_level_0.txt",proinfo->save_filepath,ti);
                p_GCP = fopen(out_file,"r");
                if(p_GCP)
                {
                    double GCP_grid;
                    if(gcp_opt == 1)
                        GCP_grid = (ortho_dx[ti] + ortho_dy[ti])/2.0*5;
                    else if(gcp_opt == 2)
                        GCP_grid = (ortho_dx[ti] + ortho_dy[ti])/2.0;
                    
                    CSize GCP_size;
                    GCP_size.width = floor(GridSize_width[ti]/GCP_grid);
                    GCP_size.height = floor(GridSize_height[ti]/GCP_grid);
                    
                    long data_length = (long )(GCP_size.width)*(long )(GCP_size.height);
                    unsigned char* GCP_value = (unsigned char*)calloc(sizeof(unsigned char),data_length);
                    for(long count = 0 ; count < data_length ; count ++)
                        GCP_value[count] = 255;
                    
                    for(long count = 0 ; count < matched_MPs.size() ;count++)
                    {
                        long pos_c = (long)((matched_MPs_ref[count].m_X - Boundary[ti][0])/GCP_grid);
                        long pos_r = (long)((Boundary[ti][3] - matched_MPs_ref[count].m_Y)/GCP_grid);
                        long index = pos_r*(long)GCP_size.width + pos_c;
                        
                        if(pos_c >= 0 && pos_c < GCP_size.width && pos_r >= 0 && pos_r < GCP_size.height)
                            GCP_value[index] = 0;
                    }
                    
                    sprintf(out_file,"%s/GCPs_%d.tif",proinfo->save_filepath,ti);
                    WriteGeotiff(out_file, GCP_value, GCP_size.width, GCP_size.height, GCP_grid, Boundary[ti][0], Boundary[ti][3], param.projection, param.utm_zone, param.bHemisphere, 1);
                    
                    char *Ifilename  = SetOutpathName(args.Image[ti]);
                    char *tmp_no_ext = remove_ext(Ifilename);
                    
                    sprintf(out_file,"%s/%s_coreg.tif",proinfo->save_filepath,tmp_no_ext);
                    
                    double left = - ImageAdjust_coreg[ti][1]*ortho_dx[ti];
                    double upper = ImageAdjust_coreg[ti][0]*ortho_dy[ti];
                    printf("coreg %s\t%d\t%d\t%f\t%f\n",out_file,OriImagesizes[ti].width,OriImagesizes[ti].height,left,upper);
                    
                    free(tmp_no_ext);
                    free(GCP_value);
                }
                
                matched_MPs_ref.clear();
                matched_MPs.clear();
                free(OriImages[ti]);
                free(ImageBoundary[ti]);
            }
            fclose(fid_out);
            RemoveFiles(proinfo,proinfo->tmpdir,Subsetfilename,py_level,0);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nsave result %f\n\n",total_gap);
            
            free(Grid_space);
            free(ortho_dx);
            free(ortho_dy);
            free(OriImages);
            free(OriImagesizes);
            
            free(ortho_minX);
            free(ortho_maxY);
            free(ImageBoundary);
            
            free(Boundary);
            free(GridSize_width);
            free(GridSize_height);
        }
    }
    
    return ImageAdjust_coreg;
}

void Preprocessing_Coreg(ProInfo *proinfo, char *save_path,uint16 **Oriimage, char **Subsetfile, uint8 py_level, CSize *Subsetsize, CSize **data_size_lr)
{
    printf("start Preprocessing\n");
    for(int count = 0; count<proinfo->number_of_images ; count++)
    {
        char* filename_py     = GetFileName(Subsetfile[count]);
        filename_py     = remove_ext(filename_py);
        
        char t_str[500];
        sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,py_level+1);
        FILE *pFile_check_file    = fopen(t_str,"rb");
        
        uint16 **pyimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
        CSize *data_size = (CSize*)malloc(sizeof(CSize)*(py_level+1));
        
        for(int i=0;i<py_level+1;i++)
            data_size[i]        = data_size_lr[count][i];
        
        long int data_length = (long int)data_size[0].height*(long int)data_size[0].width;
        pyimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);

        sprintf(t_str,"%s/%s_py_0.raw",save_path,filename_py);
        FILE *pFile   = fopen(t_str,"wb");

        fwrite(Oriimage[count],sizeof(uint16),data_length,pFile);
        fclose(pFile);
        
        for(int i=0;i<py_level;i++)
        {
            if(i==0)
                pyimg[i+1] = CreateImagePyramid(Oriimage[count],data_size[i],9,(double)(1.5));
            else
                pyimg[i+1] = CreateImagePyramid(pyimg[i],data_size[i],9,(double)(1.5));
            
            free(pyimg[i]);
            
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
        if(pFile_check_file)
            fclose(pFile_check_file);
        free(filename_py);
    }
}

void DEM_ImageCoregistration_hillshade(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt)
{
    ProInfo *proinfo = (ProInfo*)malloc(sizeof(ProInfo));
    
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    sprintf(proinfo->save_filepath,"%s",args.Outputpath);
    if(args.GCP_spacing > 0)
        proinfo->GCP_spacing = args.GCP_spacing;
    else
        proinfo->GCP_spacing = -9;
    
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
        if(Maketmpfolders(proinfo))
        {
            int py_level = proinfo->pyramid_level;
            TransParam param;
            SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
            
            char ref_dem_name[500];
            sprintf(ref_dem_name,"%s",args.Image[0]);
            printf("dem name = %s\n",ref_dem_name);
            
            double ref_minX,ref_maxY,ref_dx,ref_dy;
            CSize ref_dem_size = ReadGeotiff_info_dxy(ref_dem_name,&ref_minX,&ref_maxY,&ref_dx,&ref_dy);
            
            long cols[2] = {0, ref_dem_size.width};
            long rows[2] = {0, ref_dem_size.height};
            
            //dem image = ref
            float type(0);
            float *DEM = Readtiff_T(ref_dem_name,&ref_dem_size,cols,rows,&ref_dem_size,type);
            
            char *reffilename  = SetOutpathName(args.Image[0]);
            char *ref_no_ext = remove_ext(reffilename);
            int char_size = strlen(ref_no_ext);
            char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
                ref_fchar[c] = ref_no_ext[c];
      
            ref_fchar[char_size-3] = '\0';
            char DEM_name_refoutfile[500];
            sprintf(DEM_name_refoutfile,"%sdem.tif",ref_fchar);
            free(ref_no_ext);
            free(ref_fchar);
            free(reffilename);
            
            ref_no_ext = remove_ext(args.Image[0]);
            char_size = strlen(ref_no_ext);
            ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
                ref_fchar[c] = ref_no_ext[c];
            
            ref_fchar[char_size-3] = '\0';
            char ref_hill_inputname[500];
            sprintf(ref_hill_inputname,"%sdem_hillshade.tif",ref_fchar);
            printf("ref hill name %s\n",ref_hill_inputname);
            free(ref_no_ext);
            free(ref_fchar);
            
            //create hillshadeimage
            unsigned char *ref_hill = NULL;
            FILE *checktif = fopen(ref_hill_inputname,"r");
            if(checktif)
            {
                TIFF *refhilltif  = TIFFOpen(ref_hill_inputname,"r");
                uint8 type(0);
                ref_hill = Readtiff_T(ref_hill_inputname, &ref_dem_size, cols, rows, &ref_dem_size, type);
                TIFFClose(refhilltif);
                fclose(checktif);
                printf("loading existing hillshade image %s\n",ref_hill_inputname);
            }
            else
                ref_hill = CreateHillshade(DEM,ref_dem_size,ref_dx);
   
            //create pyramidimage
            CSize *data_size_ref = (CSize*)malloc(sizeof(CSize)*(py_level+1));
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
                    SubImages_ref[level]= CreateImagePyramid(SubImages_ref[level-1],data_size,py_kernel_size,(double)1.6);
                }
            }
            free(ref_hill);
            
            double ImageBoundary_ref[4] = {0.0};
            ImageBoundary_ref[0] = ref_minX;
            ImageBoundary_ref[1] = ref_maxY - ref_dy*ref_dem_size.height;
            ImageBoundary_ref[2] = ref_minX + ref_dx*ref_dem_size.width;
            ImageBoundary_ref[3] = ref_maxY;
            
            char out_file[500];
            sprintf(out_file,"%s/DEM_coreg_result.txt",args.Outputpath);
            FILE* fid_out         = fopen(out_file,"w");
            fprintf(fid_out,"ref DEM name\t%s\n",DEM_name_refoutfile);
            fprintf(fid_out,"DEM name\t\t\t\t\t\t\tTx[meter]\tTy[meter]\tTz(average)\tTz(median)\tTx_std[meter]\tTy_std[meter]\tTz_std(average)\tTz_std(median)\tControls_rho\tMean_all(avg)\tMedian_all(avg)\tdz_std(avg)\tMean_all(med)\tMedian_all(med)\tdz_std(med)\tNumberOfCPs\tprocessing time\n");
            
            double Grid_space = ref_dx*15.0;
            if(proinfo->GCP_spacing > 0)
                Grid_space = proinfo->GCP_spacing;
            printf("Grid space %f\n", Grid_space);
            for(int ti = 1; ti < proinfo->number_of_images ; ti++)
            {
                time_t total_ST_iter = 0, total_ET_iter = 0;
                total_ST_iter = time(0);
                
                char tar_dem_name[500];
                sprintf(tar_dem_name,"%s",args.Image[ti]);
                printf("dem name = %s\n",tar_dem_name);
                
                double tar_minX,tar_maxY,tar_dx,tar_dy;
                CSize tar_dem_size = ReadGeotiff_info_dxy(tar_dem_name,&tar_minX,&tar_maxY,&tar_dx,&tar_dy);
                
                cols[0] = 0;
                cols[1] = tar_dem_size.width;
                
                rows[0] = 0;
                rows[1] = tar_dem_size.height;
                
                //dem image
                float type(0);
                float *DEM_tar= Readtiff_T(tar_dem_name,&tar_dem_size,cols,rows,&tar_dem_size,type);
                
                char *reffilename  = SetOutpathName(args.Image[ti]);
                char *ref_no_ext = remove_ext(reffilename);
                int char_size = strlen(ref_no_ext);
                char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
                for(int c = 0 ; c < char_size - 3 ; c++)
                    ref_fchar[c] = ref_no_ext[c];
        
                ref_fchar[char_size-3] = '\0';
                char DEM_name_outfile[500];
                sprintf(DEM_name_outfile,"%sdem",ref_fchar);
                free(ref_no_ext);
                free(ref_fchar);
                free(reffilename);
                
                ref_no_ext = remove_ext(args.Image[ti]);
                char_size = strlen(ref_no_ext);
                ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
                for(int c = 0 ; c < char_size - 3 ; c++)
                    ref_fchar[c] = ref_no_ext[c];
         
                ref_fchar[char_size-3] = '\0';
                char tar_hill_inputfile[500];
                sprintf(tar_hill_inputfile,"%sdem_hillshade.tif",ref_fchar);
                printf("tar hill name %s\n",tar_hill_inputfile);
                free(ref_no_ext);
                free(ref_fchar);
                
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
                unsigned char *tar_hill = NULL;
                FILE *checktif_tar = fopen(tar_hill_inputfile,"r");
                if(checktif)
                {
                    TIFF *tarhilltif = NULL;
                    tarhilltif  = TIFFOpen(tar_hill_inputfile,"r");
                    uint8 type(0);
                    tar_hill = Readtiff_T(tar_hill_inputfile, &tar_dem_size, cols, rows, &tar_dem_size, type);
                    TIFFClose(tarhilltif);
                    fclose(checktif_tar);
                    printf("loading existing hillshade image %s\n",tar_hill_inputfile);
                }
                else
                    tar_hill = CreateHillshade(DEM_tar,tar_dem_size,tar_dx);
                
                //create pyramidimage
                CSize *data_size_tar = (CSize*)malloc(sizeof(CSize)*(py_level+1));
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
                        SubImages_tar[level]= CreateImagePyramid(SubImages_tar[level-1],data_size,py_kernel_size,(double)1.6);
                    }
                }
                free(tar_hill);
                
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
                
                D2DPOINT adjust_std;
                vector<D2DPOINT> MPs_2D;
                vector<D2DPOINT> MPs_2D_mat;
                double avg_roh;
                double ImageAdjust_coreg[2] = {0};
                
                for(int level = py_level ; level >= 0 ; level --)
                {
                    int iter_counts;
                    D2DPOINT grid_dxy_ref(ref_dx, ref_dy);
                    D2DPOINT grid_dxy_tar(tar_dx, tar_dy);
                    
                    CoregParam_Image(proinfo, ti, level, ImageAdjust_coreg, 15, SubImages_ref[level], data_size_ref[level], SubImages_tar[level],  data_size_tar[level], ImageBoundary_ref, ImageBoundary_tar, grid_dxy_ref, grid_dxy_tar, Grid_space, Boundary, &avg_roh, &iter_counts, &adjust_std, MPs_2D_mat, MPs_2D);
                    
                    printf("%d\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%4.2f\t\t%d\t%4.2f\t%d\t%4.2f\t%4.2f\n",level,ImageAdjust_coreg[0], ImageAdjust_coreg[1],-ImageAdjust_coreg[0]*tar_dy, ImageAdjust_coreg[1]*tar_dx,MPs_2D.size(),avg_roh,iter_counts,adjust_std.m_X,adjust_std.m_Y);
                    free(SubImages_tar[level]);
                }
                free(SubImages_tar);
                free(data_size_tar);
        
                int RA_iter_counts = MPs_2D.size();
                if(RA_iter_counts > 10)
                {
                    time_t total_ST_a = 0, total_ET_a = 0;
                    total_ST_a = time(0);
                    
                    const double Coreg_param[2] = {- ImageAdjust_coreg[0]*tar_dy, ImageAdjust_coreg[1]*tar_dx};
                    const double W_th = 80;
                    
                    int diff_count = 0;
                    double dH_th = 30;
                    double SD_z = 0;
                    double SD_z_med = 0;
                    double MED_z = 0;
                    double average = 0;
                    
                    bool check_stop = false;
                    int while_iter = 0;
                    double pre_std = 10000.0;
                    double change_std_ratio = 1000;
                    
                    total_ST = time(0);
                    bool check_cps = false;
                    while(!check_stop && while_iter < 50 && !check_cps)
                    {
                        double sum_diff = 0;
                        D3DPOINT *save_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),RA_iter_counts);
                        diff_count = 0;
#pragma omp parallel for reduction(+:sum_diff,diff_count) schedule(guided)
                        for(long int gcp_index = 0 ; gcp_index < RA_iter_counts ; gcp_index ++)
                        {
                            D2DPOINT gcp_coord,ref_img,tar_img;
                            gcp_coord.m_X = MPs_2D[gcp_index].m_X;
                            gcp_coord.m_Y = MPs_2D[gcp_index].m_Y;
                            
                            ref_img.m_X = (gcp_coord.m_X - ref_minX)/ref_dx;
                            ref_img.m_Y = (ref_maxY - gcp_coord.m_Y)/ref_dx;
                            
                            tar_img.m_X = ( (gcp_coord.m_X + Coreg_param[1]) - tar_minX )/tar_dx;
                            tar_img.m_Y = ( tar_maxY - (gcp_coord.m_Y + Coreg_param[0]) )/tar_dx;
                            
                            bool check_b = true;
                            if(args.check_boundary)
                            {
                                if(gcp_coord.m_X >= args.Min_X && gcp_coord.m_X <= args.Max_X && gcp_coord.m_Y >= args.Min_Y && gcp_coord.m_Y <= args.Max_Y)
                                    check_b = true;
                                else
                                    check_b = false;
                            }
                            
                            if(ref_img.m_X - 2 >= 0 && ref_img.m_X + 2 < ref_dem_size.width && ref_img.m_Y - 2 >= 0 && ref_img.m_Y + 2 < ref_dem_size.height &&
                               tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_dem_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_dem_size.height && check_b)
                            {
                                long ref_index = (long)(ref_img.m_Y)*(long)ref_dem_size.width + (long)(ref_img.m_X);
                                long tar_index = (long)(tar_img.m_Y)*(long)tar_dem_size.width + (long)(tar_img.m_X);
                                
                                if(DEM[ref_index] > -100 && DEM[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                                {
                                    vector<double> ref_patch_vecs, tar_patch_vecs;
                                    vector<long> kc_pos, kr_pos;
                                    int kernel_count = 0;
                                    for(long kr = -2 ; kr <= 2 ; kr++ )
                                    {
                                        for(long kc = -2 ; kc <= 2 ; kc++ )
                                        {
                                            long k_refindex = (long)(ref_img.m_Y + kr)*(long)ref_dem_size.width + (long)(ref_img.m_X + kc);
                                            long k_tarindex = (long)(tar_img.m_Y + kr)*(long)tar_dem_size.width + (long)(tar_img.m_X + kc);
                                            
                                            if(DEM[k_refindex] > -100 && DEM[k_refindex] < 10000 && DEM_tar[k_tarindex] > -100 && DEM_tar[k_tarindex] < 10000)
                                            {
                                                ref_patch_vecs.push_back(DEM[k_refindex]);
                                                tar_patch_vecs.push_back(DEM_tar[k_tarindex]);
                                                kc_pos.push_back(kc);
                                                kr_pos.push_back(kr);
                                            }
                                        }
                                    }
                                    
                                    kernel_count = ref_patch_vecs.size();
                                    
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
                                    
                                    for(long k_index = 0 ; k_index < kernel_count ; k_index++)
                                    {
                                        A_matrix->val[k_index][0] = kc_pos[k_index]*ref_dx;
                                        A_matrix->val[k_index][1] = kr_pos[k_index]*ref_dx;
                                        A_matrix->val[k_index][2] = 1.0;
                                        
                                        L_matrix->val[k_index][0] = ref_patch_vecs[k_index];
                                        
                                        tA_matrix->val[k_index][0] = kc_pos[k_index]*tar_dx;
                                        tA_matrix->val[k_index][1] = kr_pos[k_index]*tar_dx;
                                        tA_matrix->val[k_index][2] = 1.0;
                                        tL_matrix->val[k_index][0] = tar_patch_vecs[k_index];
                                    }
                                    
                                    //DEM height correlation
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
                                    
                                    ref_patch_vecs.clear();
                                    tar_patch_vecs.clear();
                                    kc_pos.clear();
                                    kr_pos.clear();
                                    
                                    D3DPOINT ref_normal(X_matrix->val[0][0], X_matrix->val[1][0], 1.0, 0);
                                    D3DPOINT tar_normal(tX_matrix->val[0][0], tX_matrix->val[1][0], 1.0, 0);
                                    D3DPOINT scale(1,1,1,0);
                                    double ref_angle, ref_aspect, tar_angle, tar_aspect;
                                    
                                    SlopeAspect(ref_normal, scale, &ref_angle, &ref_aspect);
                                    SlopeAspect(tar_normal, scale, &tar_angle, &tar_aspect);
                                    
                                    double SS = 1 - fabs(ref_angle - tar_angle)/90.0;
                                    double SA = 1 - fabs(ref_aspect - tar_aspect)/360.0;
                                    
                                    double W = 40*k_ncc + 40*SS + 20*SA;
                                    
                                    if(W > W_th && (DEM[ref_index] - DEM_tar[tar_index]) < dH_th + average && (DEM[ref_index] - DEM_tar[tar_index]) > - dH_th + average)
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
                                            }
                                            else
                                            {
                                                save_pts[gcp_index].flag = false;
                                            }
                                        }
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
                                    
                                    GMA_double_destroy(tA_matrix);
                                    GMA_double_destroy(tL_matrix);
                                    GMA_double_destroy(tAT_matrix);
                                    GMA_double_destroy(tATA_matrix);
                                    GMA_double_destroy(tATAI_matrix);
                                    GMA_double_destroy(tATL_matrix);
                                    GMA_double_destroy(tX_matrix);
                                    GMA_double_destroy(tAX_matrix);
                                    GMA_double_destroy(tV_matrix);
                                }
                            }
                        }
      
                        printf("iteration %d\tcps %d\n",while_iter,diff_count);
                        if(diff_count > 10)
                        {
                            average = sum_diff/diff_count;
                            
                            vector<double> save_Z;
                            double sum_var = 0;
                            
                            char dem_gcp_filename[500];
                            sprintf(dem_gcp_filename,"%s/DEM_gcps_%d.txt",args.Outputpath,ti);
                            FILE *fid_dem_gcp = fopen(dem_gcp_filename,"w");
                            
                            for(long row = 0; row < RA_iter_counts ; row++)
                            {
                                if(save_pts[row].flag == true)
                                {
                                    save_Z.push_back(save_pts[row].m_Z);
                                    
                                    sum_var += (average - save_pts[row].m_Z)*(average - save_pts[row].m_Z);
                                    fprintf(fid_dem_gcp,"%f\t%f\t%f\n",save_pts[row].m_X,save_pts[row].m_Y,save_pts[row].m_Z);
                                }
                            }
                            fclose(fid_dem_gcp);
                            
                            diff_count = save_Z.size();
                            
                            MED_z = quickselect(save_Z, diff_count, (int)(diff_count/2.0));
                            SD_z = sqrt(sum_var/diff_count);
                            
                            double iter_dH_th = SD_z*3.29;
                            change_std_ratio = fabs(pre_std - SD_z)/pre_std;
                            
                            if(change_std_ratio < 0.01 || dH_th < iter_dH_th)
                                check_stop = true;
                            
                            pre_std = SD_z;
                            dH_th = iter_dH_th;
                            
                            while_iter++;
                            
                            double sum_var_med = 0;
#pragma omp parallel for reduction(+:sum_var_med) schedule(guided)
                            for(int i=0;i<diff_count;i++)
                                sum_var_med += (MED_z -save_Z[i])*(MED_z -save_Z[i]);
                            
                            SD_z_med = sqrt(sum_var_med/save_Z.size());
                            
                            save_Z.clear();
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
                        double min_save_Z = 100000;
                        double max_save_Z = -100000;
                        
                        float* co_dem = NULL;
                        float* copoly_dem = NULL;
                        float* co_dem_diff = NULL;
                        float* copoly_dem_diff = NULL;
                        
                        double all_average = 0.0;
                        double all_med = 0.0;
                        double all_std = 0.0;
                        double all_average_med = 0.0;
                        double all_med_med = 0.0;
                        double all_std_med = 0.0;
                        
                        long tar_data_length = (long)tar_dem_size.width*(long)tar_dem_size.height;
                        if(args.check_DEM_coreg_output == 2)
                        {
                            co_dem = (float*)malloc(sizeof(float)*tar_data_length);
                            copoly_dem = (float*)malloc(sizeof(float)*tar_data_length);
                            co_dem_diff = (float*)malloc(sizeof(float)*tar_data_length);
                            copoly_dem_diff = (float*)malloc(sizeof(float)*tar_data_length);
                            
                            for(long int co_index = 0 ; co_index < tar_data_length ; co_index++)
                            {
                                co_dem[co_index] = Nodata;
                                copoly_dem[co_index] = Nodata;
                                co_dem_diff[co_index] = Nodata;
                                copoly_dem_diff[co_index] = Nodata;
                            }
                        }
                        
                        if(args.check_DEM_coreg_output == 1 || args.check_DEM_coreg_output == 2)
                        {
                            vector<double> save_dz, save_dz_med;
                            
                            double dz_sum = 0;
                            double dz_sum_med = 0;
                            min_save_Z = 100000;
                            max_save_Z = -100000;
                            for(long co_index = 0 ; co_index < tar_data_length ; co_index++)
                            {
                                long pts_row = floor(co_index/tar_dem_size.width);
                                long pts_col = co_index % tar_dem_size.width;
                                
                                D2DPOINT gcp_coord,tar_img, ref_img,gcp_coord_ref;
                                gcp_coord.m_X = pts_col*tar_dx + tar_minX + Coreg_param[1];
                                gcp_coord.m_Y = tar_maxY - pts_row*tar_dx + Coreg_param[0];
                                
                                tar_img.m_X = ( gcp_coord.m_X - tar_minX )/tar_dx;
                                tar_img.m_Y = ( tar_maxY - gcp_coord.m_Y )/tar_dx;
                                
                                gcp_coord_ref.m_X = pts_col*tar_dx + tar_minX;
                                gcp_coord_ref.m_Y = tar_maxY - pts_row*tar_dx;
                                
                                ref_img.m_X = ( gcp_coord_ref.m_X - ref_minX )/ref_dx;
                                ref_img.m_Y = ( ref_maxY - gcp_coord_ref.m_Y )/ref_dx;
                                
                                if(tar_img.m_X - 2 >= 0 && tar_img.m_X + 2 < tar_dem_size.width && tar_img.m_Y - 2 >= 0 && tar_img.m_Y + 2 < tar_dem_size.height &&
                                   ref_img.m_X >= 0 && ref_img.m_X < ref_dem_size.width && ref_img.m_Y >= 0 && ref_img.m_Y < ref_dem_size.height)
                                {
                                    long tar_index = (long)(tar_img.m_Y)*(long)tar_dem_size.width + (long)(tar_img.m_X);
                                    long ref_index = (long)(ref_img.m_Y)*(long)ref_dem_size.width + (long)(ref_img.m_X);
                                    
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
                                        
                                        double co_dem_t = value + average;
                                        double copoly_dem_t = value + MED_z;
                                        double co_dem_diff_t = (DEM[ref_index] - co_dem_t);
                                        double copoly_dem_diff_t = (DEM[ref_index] - copoly_dem_t);
                                        
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
                                            
                                            save_dz.push_back(co_dem_diff_t);
                                            save_dz_med.push_back(copoly_dem_diff_t);
                                        }
                                    }
                                }
                            }
                           
                            long dz_count = save_dz.size();
                            all_average = dz_sum/dz_count;
                            all_med = quickselect(save_dz, dz_count, (int)(dz_count/2.0));
                            
                            all_average_med = dz_sum_med/dz_count;
                            all_med_med = quickselect(save_dz_med, dz_count, (int)(dz_count/2.0));
                   
                            double all_sum_var = 0;
                            double all_sum_var_med = 0;
#pragma omp parallel for reduction(+:all_sum_var,all_sum_var_med) schedule(guided)
                            for(long index = 0 ; index < dz_count ; index++)
                            {
                                all_sum_var += (save_dz[index] - all_average)*(save_dz[index] - all_average);
                                all_sum_var_med += (save_dz_med[index] - all_average_med)*(save_dz_med[index] - all_average_med);
                            }
                            save_dz.clear();
                            save_dz_med.clear();
                            
                            all_std = sqrt(all_sum_var/dz_count);
                            all_std_med = sqrt(all_sum_var_med/dz_count);
                            
                            printf("all stat %f\t%f\t%f\t%f\t%f\t%f\n",all_average,all_med,all_std,all_average_med,all_med_med,all_std_med);
                            
                            if(args.check_DEM_coreg_output == 2)
                            {
                                WriteGeotiff(co_dems_name, co_dem, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                                WriteGeotiff(copoly_dems_name, copoly_dem, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                                
                                WriteGeotiff(co_dems_name_diff, co_dem_diff, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                                WriteGeotiff(copoly_dems_name_diff, copoly_dem_diff, tar_dem_size.width, tar_dem_size.height, tar_dx, tar_minX, tar_maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                                
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
                        
                        total_ET_iter = time(0);
                        total_gap = difftime(total_ET_iter,total_ST_iter);
                        
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
                
                MPs_2D.clear();
                free(DEM_tar);
            }
            fclose(fid_out);
            
            for(int level = py_level ; level >= 0 ; level --)
                free(SubImages_ref[level]);
    
            free(SubImages_ref);
            free(DEM);
            free(data_size_ref);
        }
    }
    free(proinfo);
}

void DEM_ImageCoregistration_GeomatricConstraint(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath, int gcp_opt)
{
    ProInfo *proinfo = new ProInfo;
    
    proinfo->number_of_images = args.number_of_images;
    proinfo->pyramid_level = args.pyramid_level;
    sprintf(proinfo->save_filepath,"%s",args.Outputpath);
    if(args.GCP_spacing > 0)
        proinfo->GCP_spacing = args.GCP_spacing;
    else
        proinfo->GCP_spacing = -9;
    
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    
    total_ST = time(0);
    
    bool check_start = true;
    if(args.check_txt_input == 0)
        check_start = OpenProject(_filename,proinfo,args);
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
        if(Maketmpfolders(proinfo))
        {
            int py_level = proinfo->pyramid_level;
            TransParam param;
            SetTranParam_fromGeoTiff(&param,proinfo->Imagefilename[0]);
            
            char ref_dem_name[500];
            sprintf(ref_dem_name,"%s",args.Image[0]);
            printf("dem name = %s\n",ref_dem_name);
            
            double ref_minX,ref_maxY,ref_dx,ref_dy;
            CSize ref_dem_size = ReadGeotiff_info_dxy(ref_dem_name,&ref_minX,&ref_maxY,&ref_dx,&ref_dy);
            
            double ImageBoundary_ref[4] = {ref_minX, ref_maxY - ref_dy*ref_dem_size.height, ref_minX + ref_dx*ref_dem_size.width, ref_maxY};
            
            long cols[2] = {0, ref_dem_size.width};
            long rows[2] = {0, ref_dem_size.height};
            
            //dem image = ref
            float type(0);
            float *DEM = Readtiff_T(ref_dem_name,&ref_dem_size,cols,rows,&ref_dem_size,type);
            
            char *reffilename  = SetOutpathName(args.Image[0]);
            char *ref_no_ext = remove_ext(reffilename);
            int char_size = strlen(ref_no_ext);
            char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
            for(int c = 0 ; c < char_size - 3 ; c++)
                ref_fchar[c] = ref_no_ext[c];
   
            ref_fchar[char_size-3] = '\0';
            char DEM_name_refoutfile[500];
            sprintf(DEM_name_refoutfile,"%s",ref_fchar);
            free(ref_no_ext);
            free(ref_fchar);
            free(reffilename);
            
            total_ET = time(0);
            total_gap = difftime(total_ET,total_ST);
            printf("\nRef preparation time %f\n\n",total_gap);
            total_ST = time(0);
            
            char out_file[500];
            sprintf(out_file,"%s/DEM_coreg_result.txt",args.Outputpath);
            FILE* fid_out         = fopen(out_file,"w");
            fprintf(fid_out,"ref DEM name\t%s\n",DEM_name_refoutfile);
            fprintf(fid_out,"DEM name\t\t\t\t\t\t\tDist_std[meter]\tTx[meter]\tTy[meter]\tTz[meter]\tTx_std[meter]\tTy_std[meter]\tTz_std[meter]\tdH_cp(mean)\tdH_cp_std(mean)\tdH_cp(med.)\tdH_cp_std(med.)\tdh_mean\t\tdh_med.\t\tdh_std\t\tNumberOfCPs\tprocessing time\n");
        
            for(int ti = 1; ti < proinfo->number_of_images ; ti++)
            {
                time_t total_ST_iter = 0, total_ET_iter = 0;
                total_ST_iter = time(0);
                
                char tar_dem_name[500];
                sprintf(tar_dem_name,"%s",args.Image[ti]);
                printf("dem name = %s\n",tar_dem_name);
                
                double tar_minX,tar_maxY,tar_dx,tar_dy;
                CSize tar_dem_size = ReadGeotiff_info_dxy(tar_dem_name,&tar_minX,&tar_maxY,&tar_dx,&tar_dy);
                
                double ImageBoundary_tar[4] = {tar_minX, tar_maxY - tar_dy*tar_dem_size.height, tar_minX + tar_dx*tar_dem_size.width, tar_maxY};
            
                cols[0] = 0;
                cols[1] = tar_dem_size.width;
                
                rows[0] = 0;
                rows[1] = tar_dem_size.height;
                
                //dem image
                float type(0);
                float *DEM_tar= Readtiff_T(tar_dem_name,&tar_dem_size,cols,rows,&tar_dem_size,type);
                
                char *reffilename  = SetOutpathName(args.Image[ti]);
                char *ref_no_ext = remove_ext(reffilename);
                int char_size = strlen(ref_no_ext);
                char *ref_fchar = (char*)malloc(sizeof(char)*(char_size-2));
                for(int c = 0 ; c < char_size - 3 ; c++)
                {
                    ref_fchar[c] = ref_no_ext[c];
                }
                ref_fchar[char_size-3] = '\0';
                char DEM_name_outfile[500];
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
                double Boundary[4] = {ImageBoundary_ref[0], ImageBoundary_ref[1], ImageBoundary_ref[2], ImageBoundary_ref[3]};
                
                
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
                CSize *data_size_tar = (CSize*)malloc(sizeof(CSize)*(py_level+1));
                SetPySizes(data_size_tar, tar_dem_size, py_level);
                
                long tin_point_num = 0;
                vector<D3DPOINT> select_pts_tar;
                SettingControls(proinfo, DEM, DEM_tar, ref_dx, tar_dx, ImageBoundary_ref, ImageBoundary_tar, Boundary, ref_dem_size, tar_dem_size, tin_point_num, select_pts_tar);
                
                printf("select point %d\t%d\n",tin_point_num, select_pts_tar.size());
                
                if(tin_point_num > 0)
                {
                    total_ET = time(0);
                    total_gap = difftime(total_ET,total_ST);
                    printf("\nTar pyramid time %f\n\n",total_gap);
                    total_ST = time(0);
                    
                    double min_H = 100000000;
                    double max_H = -100000000;
                    vector<D3DPOINT>::iterator it;
                    for(it = select_pts_tar.begin() ; it != select_pts_tar.end() ; it ++)
                    {
                        if(min_H > it->m_Z)
                            min_H = it->m_Z;
                        if(max_H < it->m_Z)
                            max_H = it->m_Z;
                    }
                    
                    double scale_factor = 1;
                    D3DPOINT coord_center((Boundary[2] - Boundary[0])/2.0 + Boundary[0], (Boundary[3] - Boundary[1])/2.0 + Boundary[1], (max_H - min_H)/2.0 + min_H, 0);
                    D3DPOINT coord_scale((Boundary[2] - Boundary[0])/2.0*scale_factor, (Boundary[3] - Boundary[1])/2.0*scale_factor, (max_H - min_H)/2.0*scale_factor, 0);
                    
                    double ori_scale_Z = coord_scale.m_Z;
                    double distance_scale = (sqrt(coord_scale.m_X*coord_scale.m_X + coord_scale.m_Y*coord_scale.m_Y + coord_scale.m_Z*coord_scale.m_Z));
                    printf("minmaxH %f\t%f\t%f\t%f\t%f\t%f\n",min_H,max_H,coord_scale.m_X,coord_scale.m_Y,coord_scale.m_Z,distance_scale);
                    
                    D3DPOINT min_nor_coord, max_nor_coord;
                    min_nor_coord.m_X = (Boundary[0] - coord_center.m_X)/coord_scale.m_X;
                    min_nor_coord.m_Y = (Boundary[1] - coord_center.m_Y)/coord_scale.m_Y;
                    min_nor_coord.m_Z = (min_H - coord_center.m_Z)/coord_scale.m_Z;
                    
                    max_nor_coord.m_X = (Boundary[2] - coord_center.m_X)/coord_scale.m_X;
                    max_nor_coord.m_Y = (Boundary[3] - coord_center.m_Y)/coord_scale.m_Y;
                    max_nor_coord.m_Z = (max_H - coord_center.m_Z)/coord_scale.m_Z;
                    
                    printf("total pts %d\n",tin_point_num);
                    
                    vector<double> weight(0.0);
                    vector<double>::iterator w_it;
                    for(it = select_pts_tar.begin() ; it != select_pts_tar.end() ; it ++)
                        weight.push_back(0.0);
                    
                    double level_ref_dx = ref_dx;
                    double level_tar_dx = tar_dx;
                    
                    CSize GridSize;
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
                    
                    double* dX = NULL;
                    double* sigmaX = (double*)calloc(sizeof(double),7);
                    double sigma0;
                    int while_iter = 0;
                    bool check_stop = false;
                    double W_th = 0.8;
                    time_t total_ST_ad = 0, total_ET_ad = 0;
                    total_ST_ad = time(0);
                    
                    double sum_distance = 0;
                    long final_number_of_selected = 0;
                    double average_distance, MED_distance;
                    double vertical_average_distance;
                    double SD_distance=0;
                    double SD_distance_vertical=0;
                    double SD_z_med = 0;
                    const double th_dH_nonscale = 30;
                    double th_dH = th_dH_nonscale/ori_scale_Z;
                    double pre_sigma = 10000.0;
                    bool check_null_dem = false;
                    
                    double min_X = Boundary[0];
                    double max_X = Boundary[2];
                    double min_Y = Boundary[1];
                    double max_Y = Boundary[3];
                    
                    int max_W_update = 4;
                    while(!check_stop && while_iter < 50)
                    {
                        vector<D3DPOINT> control_pts;
                        vector<D3DPOINT> transformed_coord;
                        vector<D3DPOINT> tar_normal_save;
                        vector<double> SDistance;
                        vector<double> weight_save;
                        
                        int* hist_W = (int*)calloc(sizeof(int),100);
                        
                        long number_of_selected = 0;
                        sum_distance = 0;
                        
                        for(long index = 0 ; index < select_pts_tar.size() ; index++)
                        {
                            if(select_pts_tar[index].flag)
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
                                        D3DPOINT normalized_pt = Normalize_coord(select_pts_tar[index],coord_center,coord_scale);
                                        D3DPOINT transformed_coord_pt = ConformalTransform(normalized_pt,conparam);
                                        
                                        D3DPOINT ref_normal_ori;
                                        double ref_iter_height;
                                        double* ref_array = (double*)calloc(sizeof(double),9);
                                        D3DPOINT ref_normal = FindNormal(&ref_normal_ori,DEM, normalized_pt,coord_center,coord_scale, ImageBoundary_ref, conparam, ref_dem_size, level_ref_dx, ref_array, &ref_iter_height, 0);
                                        
                                        if(ref_normal.m_X > -999)
                                        {
                                            double* tar_array = (double*)calloc(sizeof(double),9);
                                            double tar_iter_height;
                                            
                                            D3DPOINT tar_normal_ori;
                                            D3DPOINT tar_normal = FindNormal(&tar_normal_ori,DEM_tar, transformed_coord_pt,coord_center,coord_scale, ImageBoundary_tar, conparam, tar_dem_size, level_ref_dx, tar_array, &tar_iter_height, 1);
                                            
                                            if(tar_normal.m_X > -999)
                                            {
                                                double ref_slope, ref_aspect, tar_slope, tar_aspect;
                                                SlopeAspect(ref_normal,coord_scale,&ref_slope,&ref_aspect);
                                                SlopeAspect(tar_normal,coord_scale,&tar_slope,&tar_aspect);
                                                
                                                double ncc = Correlate(ref_array,tar_array,9);
                                                if (ncc != -99)
                                                    ncc = (ncc + 1)/2.0;
                                                else
                                                    ncc = 0.0;
                                                
                                                double ration_slope = 1 - fabs(ref_slope - tar_slope)/90.0;
                                                double ration_aspect = 1 - fabs(ref_aspect - tar_aspect)/360.0;
                                                
                                                double W = (40*ncc + 40*ration_slope + 20*ration_aspect)/100.0;
                                                
                                                int t_W = floor(W*100);
                                                if(t_W >= 0 && t_W < 100)
                                                    hist_W[t_W]++;
                                                           
                                                if(while_iter < max_W_update)
                                                {
                                                    w_it = weight.begin() + index;
                                                    *w_it = W;
                                                }
                                                
                                                bool check_W_pos = false;
                                                if(W > W_th && ref_slope > 0 && tar_slope > 0 /*&& ref_aspect > 0 && tar_aspect > 0*/ &&
                                                tar_iter_height > min_nor_coord.m_Z && tar_iter_height < max_nor_coord.m_Z &&
                                                transformed_coord_pt.m_X > min_nor_coord.m_X && transformed_coord_pt.m_X < max_nor_coord.m_X &&
                                                transformed_coord_pt.m_Y > min_nor_coord.m_Y && transformed_coord_pt.m_Y < max_nor_coord.m_Y &&
                                                transformed_coord_pt.m_Z > min_nor_coord.m_Z && transformed_coord_pt.m_Z < max_nor_coord.m_Z)
                                                    check_W_pos = true;
                                                
                                                if(check_W_pos)
                                                {
                                                    D3DPOINT t_coord = normalized_pt;
                                                    t_coord.m_Z = tar_iter_height;
                                                    
                                                    D3DPOINT ref_pts = SurfaceDistance_ori(tar_normal_ori,DEM,tar_normal, t_coord, ImageBoundary_ref,ref_dem_size, level_ref_dx, conparam, coord_center,coord_scale, ref_iter_height);
                                                    
                                                    double D = -(tar_normal.m_X*t_coord.m_X + tar_normal.m_Y*t_coord.m_Y + tar_normal.m_Z*t_coord.m_Z);
                                                    double diff_distance = -(tar_normal.m_X*ref_pts.m_X + tar_normal.m_Y*ref_pts.m_Y + tar_normal.m_Z*ref_pts.m_Z + D);
                                                    
                                                    double dist_x = (t_coord.m_Z - ref_pts.m_Z)*coord_scale.m_Z;
                                                    
                                                    if(fabs(diff_distance) < th_dH )
                                                    {
                                                        transformed_coord.push_back(transformed_coord_pt);
                                                        tar_normal_save.push_back(tar_normal);
                                                        SDistance.push_back(diff_distance);
                                                        weight_save.push_back(weight[index]);
                                                        
                                                        
                                                        sum_distance += diff_distance;
                                                        control_pts.push_back(ref_pts);
                                                        
                                                        number_of_selected++;
                                                        
                                                        select_pts_tar[index].flag = true;
                                                    }
                                                }
                                                else
                                                {
                                                    if(while_iter < max_W_update)
                                                        select_pts_tar[index].flag = false;
                                                }
                                            }
                                            else
                                            {
                                                if(while_iter < max_W_update)
                                                    select_pts_tar[index].flag = false;
                                            }
                                        }
                                        else
                                        {
                                            if(while_iter < max_W_update)
                                                select_pts_tar[index].flag = false;
                                        }
                                    }
                                    else
                                    {
                                        if(while_iter < max_W_update)
                                            select_pts_tar[index].flag = false;
                                    }
                                }
                                else
                                {
                                    if(while_iter < max_W_update)
                                        select_pts_tar[index].flag = false;
                                }
                            }
                            else
                            {
                                if(while_iter < max_W_update)
                                    select_pts_tar[index].flag = false;
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
                            total_ET_ad = time(0);
                            total_gap = difftime(total_ET_ad,total_ST_ad);
                            printf("\npre matrix adjustment time %f\t%d\n\n",total_gap,number_of_selected);
                            total_ST_ad = time(0);
                            
                            dX = CoeffMatrix_25D(coord_center,coord_scale, number_of_selected, transformed_coord, SDistance, conparam, tar_normal_save,weight,sigmaX,&sigma0);
                            
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
                            
                            average_distance = sum_distance/final_number_of_selected;
                            
                            vector<double> save_Z_vertical;
                            double sum_var = 0;
                            char dem_gcp_filename[500];
                            sprintf(dem_gcp_filename,"%s/DEM_gcps_%d.txt",args.Outputpath,ti);
                            FILE *fid_dem_gcp = fopen(dem_gcp_filename,"w");
                            double sum_distance_vertical = 0;
                            
                            for(long i = 0 ; i < transformed_coord.size() ; i++)
                            {
                                sum_var += (average_distance - SDistance[i])*(average_distance - SDistance[i]);
                                fprintf(fid_dem_gcp,"%f\t%f\t%f\n",control_pts[i].m_X*coord_scale.m_X + coord_center.m_X,
                                        control_pts[i].m_Y*coord_scale.m_Y + coord_center.m_Y,
                                        control_pts[i].m_Z*coord_scale.m_Z + coord_center.m_Z);
                                
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
                                        
                                        double co_dem_t = value + conparam.Tz*coord_scale.m_Z;
                                        double co_dem_diff = (DEM[ref_index] - co_dem_t);
                                        save_Z_vertical.push_back(co_dem_diff);
                                        sum_distance_vertical += co_dem_diff;
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
                                
                                double hist_ratio = (double)sum_hist_count/(double)tin_point_num;
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
                            long count_vertical = save_Z_vertical.size();
                            vertical_average_distance = sum_distance_vertical/count_vertical;
                            double sum_var_vertical = 0;
                            for(long i=0;i<count_vertical;i++)
                                sum_var_vertical += (vertical_average_distance - save_Z_vertical[i])*(vertical_average_distance - save_Z_vertical[i]);
                    
                            SD_distance_vertical = sqrt(sum_var_vertical/count_vertical);
                            
                            MED_distance = quickselect(save_Z_vertical, count_vertical, (int)(count_vertical/2.0));
                            
                            double sum_var_med = 0;
                            
                            for(int i=0;i<count_vertical;i++)
                                sum_var_med += (MED_distance -save_Z_vertical[i])*(MED_distance -save_Z_vertical[i]);
                   
                            printf("iter %d\tparam %f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",while_iter,conparam.Tx*coord_scale.m_X,conparam.Ty*coord_scale.m_Y,conparam.Tz*coord_scale.m_Z,sigmaX[4]*sigma0*coord_scale.m_X,sigmaX[5]*sigma0*coord_scale.m_Y,sigmaX[6]*sigma0*coord_scale.m_Z, number_of_selected,sigma0*coord_scale.m_Z,change_ratio_sigma,th_dH*ori_scale_Z);
                            while_iter++;
                            
                            SD_z_med = sqrt(sum_var_med/count_vertical);
                            printf("\nTz average and variation  %d\t%f\t%f\t%f\t%f\t%f\n\n",final_number_of_selected,th_dH*ori_scale_Z,vertical_average_distance,MED_distance,SD_distance_vertical,SD_z_med);
                            save_Z_vertical.clear();
                            total_ET = time(0);
                            total_gap = difftime(total_ET,total_ST);
                            printf("\nadjustment time %f\n\n",total_gap);
                            total_ST = time(0);
                        }
                        
                        free(hist_W);
                        
                        control_pts.clear();
                        transformed_coord.clear();
                        tar_normal_save.clear();
                        SDistance.clear();
                        weight_save.clear();
                    }
                       
                    if(!check_null_dem)
                    {
                        total_ET_iter = time(0);
                        double total_gap_iter = difftime(total_ET_iter,total_ST_iter);
                        printf("\npair time %f\n\n",total_gap_iter);
                        total_ET_iter = time(0);
                            
                        float* co_dem = NULL;
                        float* copoly_dem = NULL;
                        vector<double> save_dz;
                        
                        double all_average = 0.0;
                        double all_med = 0.0;
                        double all_std = 0.0;
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
                            double dz_sum = 0;
                            
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
                                        
                                        long index1  = (t_col_int   ) + (t_row_int   )*(long)tar_dem_size.width;
                                        long index2  = (t_col_int +1) + (t_row_int   )*(long)tar_dem_size.width;
                                        long index3  = (t_col_int   ) + (t_row_int +1)*(long)tar_dem_size.width;
                                        long index4  = (t_col_int +1) + (t_row_int +1)*(long)tar_dem_size.width;
                                        
                                        double value1      = DEM_tar[index1];
                                        double value2      = DEM_tar[index2];
                                        double value3      = DEM_tar[index3];
                                        double value4      = DEM_tar[index4];
                                        
                                        double value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                        + value3*(1-dcol)*drow + value4*dcol*drow;
                                        
                                        float co_dem_t = value + conparam.Tz*coord_scale.m_Z;
                                        float co_dem_diff = (DEM[ref_index] - co_dem_t);
                                        
                                        if(args.check_DEM_coreg_output == 2)
                                        {
                                            co_dem[tar_index] = co_dem_t;
                                            copoly_dem[tar_index] = co_dem_diff;
                                        }
                                        
                                        if(fabs(co_dem_diff) < th_dH_nonscale)
                                        {
                                            dz_sum += co_dem_diff;
                                            save_dz.push_back(co_dem_diff);
                                        }
                                    }
                                }
                            }
                            
                            total_ET = time(0);
                            total_gap = difftime(total_ET,total_ST);
                            printf("\nWhole image stat %f\n\n",total_gap);
                            total_ST = time(0);
                            
                            long dz_count = save_dz.size();
                            all_average = dz_sum/dz_count;
                            all_med = quickselect(save_dz, dz_count, (int)(dz_count/2.0));
                            double all_sum_var = 0;
                            
                #pragma omp parallel for reduction(+:all_sum_var) schedule(guided)
                            for(long index = 0 ; index < dz_count ; index++)
                                all_sum_var += (save_dz[index] - all_average)*(save_dz[index] - all_average);
                        
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
                    
                    
                    weight.clear();
                    printf("done\n");
                }
                else
                {
                    fprintf(fid_out,"%s\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\tNaN\t\t%d\t\tNaN\n",DEM_name_outfile,tin_point_num);
                }
                free(DEM_tar);
                free(data_size_tar);
                select_pts_tar.clear();
            }
            fclose(fid_out);
             
            free(DEM);
        }
    }
    delete proinfo;
}

D3DPOINT FindNormal(D3DPOINT *normal_ori, const float* dem, const D3DPOINT Pos, const D3DPOINT Mean, const D3DPOINT Scale, const double* Boundary, const Conformalparam X, const CSize tinsize, const double Gridspace, double *roh_array, double *Z, const bool check_tar)
{
    D3DPOINT normal_vector;
    
    const double geo_X = Pos.m_X*Scale.m_X + Mean.m_X;
    const double geo_Y = Pos.m_Y*Scale.m_Y + Mean.m_Y;
    const double col_float = (geo_X- Boundary[0])/Gridspace;
    const double row_float = (Boundary[3]-geo_Y)/Gridspace;
    
    const int col = floor(col_float);
    const int row = floor(row_float);
    
    double X_diff = col_float - col;
    double Y_diff = row_float - row;
    
    //normalized height setting
    const long index11 =  row   *tinsize.width + (col    );
    const long index12 =  row   *tinsize.width + (col + 1);
    const long index21 = (row+1)*tinsize.width + (col    );
    const long index22 = (row+1)*tinsize.width + (col + 1);
    
    const long index31 = (row-1)*tinsize.width + (col - 1);
    const long index32 = (row-1)*tinsize.width + (col    );
    const long index33 = (row-1)*tinsize.width + (col + 1);
    
    const long index41 =  row   *tinsize.width + (col - 1);
    const long index51 = (row+1)*tinsize.width + (col - 1);
    
    const long data_length = (long)tinsize.width*(long)tinsize.height;
    
    if(index11 >= 0 && index11 < data_length &&
       index12 >= 0 && index12 < data_length &&
       index21 >= 0 && index21 < data_length &&
       index22 >= 0 && index22 < data_length &&
       index31 >= 0 && index31 < data_length &&
       index32 >= 0 && index32 < data_length &&
       index33 >= 0 && index33 < data_length &&
       index41 >= 0 && index41 < data_length &&
       index51 >= 0 && index51 < data_length && row + 1 < tinsize.height && row - 1 >= 0 && col + 1 < tinsize.width && col -1 >= 0)
    {
        double Patch11,Patch12,Patch21,Patch22,Patch31,Patch32,Patch33,Patch41,Patch51;
        D3DPOINT P1,P2,P3,P4,P5,P6,P7,P8,P9;
        
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
                
                const double R11 = R.m11;
                const double R12 = R.m12;
                const double R13 = R.m13;
                
                const double R21 = R.m21;
                const double R22 = R.m22;
                const double R23 = R.m23;
                
                const double R31 = R.m31;
                const double R32 = R.m32;
                const double R33 = R.m33;
                
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
            
            //vector calculation
            D3DPOINT v1 = P3 - P1;
            D3DPOINT v2 = P2 - P1;
            D3DPOINT n1, n2;
            
            n1.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n1.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n1.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            v1 = P4 - P1;
            v2 = P3 - P1;
            
            n2.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n2.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n2.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            //plane vector decide
            D3DPOINT n;
            if( (X_diff == 0 && Y_diff == 0) || (X_diff != 0 & Y_diff == 0) || ( X_diff > 0 & Y_diff > 0 & (X_diff > Y_diff)) )
                n = n1;
            else
                n = n2;
               
            double mag = SQRT(n);
            n.m_X /= mag;
            n.m_Y /= mag;
            n.m_Z /= mag;
            
            double d = - P1.m_X*n.m_X - P1.m_Y*n.m_Y - P1.m_Z*n.m_Z;
            *Z = - (Pos.m_X*n.m_X + Pos.m_Y*n.m_Y + d)/n.m_Z;
            
            normal_vector = n;
            
            //denormailzed vector
            D3DPOINT P1_ori,P2_ori,P3_ori,P4_ori;
            P1_ori = Denormalize_coord(P1,Mean,Scale);
            P2_ori = Denormalize_coord(P2,Mean,Scale);
            P3_ori = Denormalize_coord(P3,Mean,Scale);
            P4_ori = Denormalize_coord(P4,Mean,Scale);
            
            P1 = P1_ori;
            P2 = P2_ori;
            P3 = P3_ori;
            P4 = P4_ori;
            
            v1 = P3 - P1;
            v2 = P2 - P1;
            
            n1.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n1.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n1.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            v1 = P4 - P1;
            v2 = P3 - P1;
            
            n2.m_X =    v1.m_Y*v2.m_Z  - v2.m_Y*v1.m_Z;
            n2.m_Y = - (v1.m_X*v2.m_Z  - v2.m_X*v1.m_Z);
            n2.m_Z =    v1.m_X*v2.m_Y  - v2.m_X*v1.m_Y;
            
            //plane vector decide
            if( (X_diff == 0 && Y_diff == 0) || (X_diff != 0 & Y_diff == 0) || ( X_diff > 0 & Y_diff > 0 & (X_diff > Y_diff)) )
                n = n1;
            else
                n = n2;
               
            mag = SQRT(n);
            n.m_X /= mag;
            n.m_Y /= mag;
            n.m_Z /= mag;
            
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

D3DPOINT SurfaceDistance_ori(const D3DPOINT tar_normal_ori, const float* ref_dem, const D3DPOINT tar_normal, const D3DPOINT tar_pts, const double *tin_boundary,const CSize tinsize, const double Gridspace, const Conformalparam param, const D3DPOINT Mean, const D3DPOINT Scale, const double p_ref_z)
{
    D3DPOINT ref_pts;
    
    //denomalizing normal vector
    D3DPOINT tar_normal_denor(tar_normal_ori);
    
    const double mag              = SQRT(tar_normal_denor);
    tar_normal_denor.m_X = tar_normal_denor.m_X/mag;
    tar_normal_denor.m_Y = tar_normal_denor.m_Y/mag;
    tar_normal_denor.m_Z = tar_normal_denor.m_Z/mag;
    
    //vector slope about x and y directions
    const double vector_slope_x      = atan2(fabs(tar_normal_denor.m_Z),fabs(tar_normal_denor.m_X));
    const double vector_slope_y      = atan2(fabs(tar_normal_denor.m_Z),fabs(tar_normal_denor.m_Y));
        
    //set up for start points
    const D3DPOINT line_node_pt(Denormalize_coord(tar_pts,Mean,Scale));
    D3DPOINT pre_pt_on_surface(line_node_pt);
    pre_pt_on_surface.m_Z = p_ref_z*Scale.m_Z + Mean.m_Z;
    
    int max_iter            = 10;
    int count = 1;
    bool check_stop = false;
    Conformalparam temp_X;
    temp_X.scale = 1.0;
    temp_X.omega = 0.0;
    temp_X.phi = 0.0;
    temp_X.kappa = 0.0;
    temp_X.Tx = 0.0;
    temp_X.Ty = 0.0;
    temp_X.Tz = 0.0;
    while(count < max_iter && !check_stop)
    {
        D3DPOINT line_node_pt_after;
        line_node_pt_after.m_X = line_node_pt.m_X + (pre_pt_on_surface.m_Z - line_node_pt.m_Z)*(tar_normal_denor.m_X/tar_normal_denor.m_Z);
        line_node_pt_after.m_Y = line_node_pt.m_Y + (pre_pt_on_surface.m_Z - line_node_pt.m_Z)*(tar_normal_denor.m_Y/tar_normal_denor.m_Z);
        line_node_pt_after.m_Z = pre_pt_on_surface.m_Z;
        
        D3DPOINT find_pt    = Normalize_coord(line_node_pt_after,Mean,Scale);
        
        double temp_array[9];
        double Z_ref;
        D3DPOINT normal_ori;
        D3DPOINT temp = FindNormal(&normal_ori,ref_dem, find_pt, Mean, Scale, tin_boundary, temp_X, tinsize, Gridspace, temp_array, &Z_ref, false);
        
        Z_ref                           = Z_ref*Scale.m_Z + Mean.m_Z;
    
        D3DPOINT after_pt_on_surface;
        after_pt_on_surface.m_X  = line_node_pt_after.m_X;
        after_pt_on_surface.m_Y  = line_node_pt_after.m_Y;
        after_pt_on_surface.m_Z  = Z_ref;
        
        const double dx                        = pre_pt_on_surface.m_X - after_pt_on_surface.m_X;
        const double dy                        = pre_pt_on_surface.m_Y - after_pt_on_surface.m_Y;
        const double dz                        = pre_pt_on_surface.m_Z - after_pt_on_surface.m_Z;
        const double diff                      = sqrt(dx*dx + dy*dy + dz*dz);
        
        if(diff < 0.1)
        {
            ref_pts                     = Normalize_coord(after_pt_on_surface,Mean,Scale);
            check_stop = true;
        }
        else
        {
            const double pre_angle_x               = atan2(fabs(dz),fabs(dx));
            const double pre_angle_y               = atan2(fabs(dz),fabs(dy));
            const bool divergent_index_x          = (pre_angle_x >= vector_slope_x & tar_normal_denor.m_X != 0 & dx != 0);
            const bool divergent_index_y          = (pre_angle_y >= vector_slope_y & tar_normal_denor.m_Y != 0 & dy != 0);
            if( divergent_index_x || divergent_index_y)
            {
                D3DPOINT temp_pt;
                temp_pt.m_X                  = pre_pt_on_surface.m_X + (after_pt_on_surface.m_X - pre_pt_on_surface.m_X)/2.0;
                temp_pt.m_Y                  = pre_pt_on_surface.m_Y + (after_pt_on_surface.m_Y - pre_pt_on_surface.m_Y)/2.0;
                
                find_pt                     = Normalize_coord(temp_pt,Mean,Scale);
                double temp_Z = 0;
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

void SlopeAspect(D3DPOINT normal, const D3DPOINT scale, double *slope, double *aspect)
{
    normal.m_X = normal.m_X/(scale.m_X);
    normal.m_Y = normal.m_Y/(scale.m_Y);
    normal.m_Z = normal.m_Z/(scale.m_Z);
    
    const double denominator = SQRT(normal);
    const double value = normal.m_Z/denominator;
    
    *slope = acos(value)*RadToDeg;
    
    if(*slope <= 0 && *slope >= -90)
        *slope = fabs(*slope);
    else if(*slope <= -270 && *slope >= -360)
        *slope = 360 + *slope;
    else if(*slope >= 270 && *slope <= 360)
        *slope = 360 - *slope;
    
    *aspect = 90 - atan2(normal.m_Y,normal.m_X)*(180.0/PI);
    if(*aspect < 0)
        *aspect = *aspect + 360;
}

D3DPOINT ConformalTransform(const D3DPOINT input, const Conformalparam param)
{
    D3DPOINT out;
    
    const RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);
       
    const double S = param.scale;
    //3D conformal transformation
    out.m_X = S*(R.m11*input.m_X + R.m21*input.m_Y + R.m31*input.m_Z) + param.Tx;
    out.m_Y = S*(R.m12*input.m_X + R.m22*input.m_Y + R.m32*input.m_Z) + param.Ty;
    out.m_Z = S*(R.m13*input.m_X + R.m23*input.m_Y + R.m33*input.m_Z) + param.Tz;
    
    return out;
}

D3DPOINT Normalize_coord(const D3DPOINT input, const D3DPOINT Mean, const D3DPOINT Scale)
{
    D3DPOINT out;
    
    out.m_X = (input.m_X - Mean.m_X) / Scale.m_X;
    out.m_Y = (input.m_Y - Mean.m_Y) / Scale.m_Y;
    out.m_Z = (input.m_Z - Mean.m_Z) / Scale.m_Z;
    
    return out;
}

D3DPOINT Denormalize_coord(const D3DPOINT input, const D3DPOINT Mean, const D3DPOINT Scale)
{
    D3DPOINT out;
    
    out.m_X = input.m_X * Scale.m_X + Mean.m_X;
    out.m_Y = input.m_Y * Scale.m_Y + Mean.m_Y;
    out.m_Z = input.m_Z * Scale.m_Z + Mean.m_Z;
    
    return out;
}

unsigned char* CreateHillshade(const float* _input, const CSize _img_size, const double grid_size)
{
    unsigned char* result_img;
    const long data_length = (long)_img_size.height*(long)_img_size.width;
    result_img = (unsigned char*)calloc(sizeof(unsigned char),data_length);
    
    const double SobleFilter_X[3][3] = { {-1, 0 , 1} , {-2, 0 , 2}, {-1, 0 , 1} };
    const double SobleFilter_Y[3][3] = { {1, 2 , 1} , {0, 0 , 0}, {-1, -2 , -1} };
    const double a = 0;
    const double b = 1.0/sqrt(2.0);
    const double alpha = 45.0;
    const double beta = 315.0;
    
#pragma omp parallel for schedule(guided)
    for(long r=0;r<_img_size.height;r++)
    {
        for(long c=0;c<_img_size.width;c++)
        {
            long ori_index = r*(long)_img_size.width + c;
            if(_input[ori_index] > -100)
            {
                double maxdx = 0;
                double maxdy = 0;
                bool check_null = false;
                for(long row = -1; row <= 1; row++)
                {
                    for(long col = -1 ; col <= 1 ;col++)
                    {
                        long row_index = r + row;
                        long col_index = c + col;
                        if(col_index >= 0 && col_index < _img_size.width && row_index >= 0 && row_index < _img_size.height)
                        {
                            long index = row_index*(long)_img_size.width + col_index;
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
    return result_img;
}

void SettingControls(const ProInfo *proinfo, const float* DEM_ref, const float* DEM_tar, const double grid_size_ref, const double grid_size_tar, const double *boundary_ref, const double *boundary_tar, const double* overlapped_br, const CSize img_size_ref, const CSize img_size_tar, long &tin_point_num, vector<D3DPOINT> &select_pts)
{
    tin_point_num = 0;
    
    double control_spacing = grid_size_tar*20;
    control_spacing = 50;
    if(proinfo->GCP_spacing > 0)
        control_spacing = proinfo->GCP_spacing;
    
    printf("control_spacing %f\n",control_spacing);
    
    CSize coord_size;
    coord_size.width = ceil((overlapped_br[2] - overlapped_br[0])/control_spacing);
    coord_size.height = ceil((overlapped_br[3] - overlapped_br[1])/control_spacing);
    
    printf("coord_size %d\t%d\n",coord_size.width,coord_size.height);
    
    int kernel_size;
    if(grid_size_tar >= 2)
        kernel_size = 1;
    else if(grid_size_tar == 1)
        kernel_size = 2;
    else
        kernel_size = 4;
    
    const long patch_size = (2*kernel_size + 1)*(2*kernel_size + 1);
    
    for(long r=0;r<coord_size.height;r++)
    {
        for(long c=0;c<coord_size.width;c++)
        {
            long control_index = r*(long)coord_size.width + c;
            
            double t_coord_x = overlapped_br[0] + c*control_spacing;
            double t_coord_y = overlapped_br[3] - r*control_spacing;
            
            long pos_c_ref = (t_coord_x - boundary_ref[0])/grid_size_ref;
            long pos_r_ref = (boundary_ref[3] - t_coord_y)/grid_size_ref;
            
            long pos_c_tar = (t_coord_x - boundary_tar[0])/grid_size_tar;
            long pos_r_tar = (boundary_tar[3] - t_coord_y)/grid_size_tar;
            
            if(pos_c_ref - kernel_size >= 0 && pos_c_ref + kernel_size < img_size_ref.width && pos_r_ref - kernel_size >= 0 && pos_r_ref + kernel_size < img_size_ref.height &&
               pos_c_tar - kernel_size >= 0 && pos_c_tar + kernel_size < img_size_tar.width && pos_r_tar - kernel_size >= 0 && pos_r_tar + kernel_size < img_size_tar.height)
            {
                long ref_index = pos_r_ref*(long)img_size_ref.width + pos_c_ref;
                long tar_index = pos_r_tar*(long)img_size_tar.width + pos_c_tar;
                
                double *ref_array = (double*)calloc(sizeof(double),(2*kernel_size + 1)*(2*kernel_size + 1));
                double *tar_array = (double*)calloc(sizeof(double),(2*kernel_size + 1)*(2*kernel_size + 1));
                
                for(long k = -kernel_size ; k <= kernel_size ; k++)
                {
                    for(long l = -kernel_size ; l <= kernel_size ; l++)
                    {
                        long array_index = (k+kernel_size)*(2*kernel_size + 1) + (l+kernel_size);
                        long k_ref_index = (pos_r_ref + k)*(long)img_size_ref.width + pos_c_ref + l;
                        long k_tar_index = (pos_r_tar + k)*(long)img_size_tar.width + pos_c_tar + l;
                        ref_array[array_index] = DEM_ref[k_ref_index];
                        tar_array[array_index] = DEM_tar[k_tar_index];
                    }
                }
                
                double ncc = Correlate(ref_array,tar_array,patch_size);
                free(ref_array);
                free(tar_array);
                
                if(ncc > 0.7 && DEM_ref[ref_index] > -100 && DEM_ref[ref_index] < 10000 && DEM_tar[tar_index] > -100 && DEM_tar[tar_index] < 10000)
                {
                    D3DPOINT temp;
                    
                    temp.m_X = t_coord_x;
                    temp.m_Y = t_coord_y;
                    temp.m_Z = DEM_tar[tar_index];
                    temp.flag = true;
                    
                    select_pts.push_back(temp);
                    
                    tin_point_num++;
                }
            }
        }
    }

    printf("SettingControls %d\n",tin_point_num);
}

double* CoeffMatrix_25D(const D3DPOINT coord_center, const D3DPOINT coord_scale, const long selected_pts, const vector<D3DPOINT> &transformed_coord, const vector<double> &dH, const Conformalparam param, const vector<D3DPOINT> &tar_normal, const vector<double> &weight, double* sigmaX, double *sigma0)
{
    double C[21] = {0.0};
    
    const RM R = MakeRotationMatrix(param.omega, param.phi, param.kappa);
    
    const double A1 = param.omega*DegToRad;
    const double A2 = param.phi*DegToRad;
    const double A3 = param.kappa*DegToRad;
    const double S = param.scale;
    
    const double R11 = R.m11;
    const double R12 = R.m12;
    const double R13 = R.m13;
    
    const double R21 = R.m21;
    const double R22 = R.m22;
    const double R23 = R.m23;
    
    const double R31 = R.m31;
    const double R32 = R.m32;
    const double R33 = R.m33;
    
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
            Wb_matrix->val[i][j] = 0.0;
 
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
    
    for(long count = 0 ; count < selected_pts ; count++)
    {
        const double P1 = transformed_coord[count].m_X;
        const double P2 = transformed_coord[count].m_Y;
        const double P3 = transformed_coord[count].m_Z;
        
        const double Gx = tar_normal[count].m_X;///coord_scale.m_X;
        const double Gy = tar_normal[count].m_Y;///coord_scale.m_Y;
        const double Gz = -tar_normal[count].m_Z;///coord_scale.m_Z;
        
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
        
        AW_matrix->val[count][0] = A_matrix->val[count][0]*weight[count];
        AW_matrix->val[count][1] = A_matrix->val[count][1]*weight[count];
        AW_matrix->val[count][2] = A_matrix->val[count][2]*weight[count];
        AW_matrix->val[count][3] = A_matrix->val[count][3]*weight[count];
        AW_matrix->val[count][4] = A_matrix->val[count][4]*weight[count];
        AW_matrix->val[count][5] = A_matrix->val[count][5]*weight[count];
        AW_matrix->val[count][6] = A_matrix->val[count][6]*weight[count];
        
        L_matrix->val[count][0] = dH[count];
    }
    
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
    
    *sigma0 = sqrt((VTV_matrix->val[0][0])/(selected_pts - 7));
    double* dx = (double*)calloc(sizeof(double),7);
    
    for(int i = 0 ; i < 7 ; i++)
    {
        sigmaX[i] = sqrt(Qxx_matrix->val[i][i]);
        dx[i] = X_matrix->val[i][0];
        
        printf("dx %d\t%f\tsigmaX %f\n",i,dx[i],sigmaX[i]);
    }
 
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
