/*
 * Copyright 2017 Myoung-Jong Noh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "setsm_code.hpp"

#ifdef BUILDMPI
#include "mpi.h"
#endif

const char setsm_version[] = "4.3.0";

int main(int argc,char *argv[])
{
    setlogmask (LOG_UPTO (LOG_NOTICE));

    openlog ("setsm", LOG_CONS | LOG_PID | LOG_NDELAY, LOG_NOTICE);

    syslog (LOG_NOTICE, "git_description: %s user: %s\n", GIT_DESCRIPTION, getlogin());
    #ifdef DEBUG_SYSLOG
      printf ("setsm git_description: %s user: %s\n", GIT_DESCRIPTION, getlogin());
    #endif

    closelog ();
    
    setbuf(stdout, NULL);
    TIFFSetWarningHandler(NULL);
    char* projectfilename   = (char*)"default.txt";
    int i=0;
    char* image1_name = NULL;
    char* image2_name = NULL;
    char* output_directory_name = NULL;
    ARGINFO args;
    
    args.sensor_type = SB; //1 = satellite, 2 = Aerial Photo
    args.pyramid_level = 4;
    args.SDM_SS = 3;
    args.SDM_AS = 20.0;
    args.SDM_days = 1;
    args.number_of_images = 2;
    args.check_arg = 0;
    args.check_DEM_space = false;
    args.check_Threads_num = false;
    args.check_seeddem = false;
    args.check_minH = false;
    args.check_maxH = false;
    args.check_tiles_SR = false;
    args.check_tiles_ER = false;
    args.check_tiles_SC = false;
    args.check_tiles_EC = false;
    args.check_RA_line  = false;
    args.check_RA_sample= false;
    args.check_gridonly = false;
    args.check_RA_tileR = false;
    args.check_RA_tileC = false;
    args.check_tilesize = false;
    args.check_boundary = false;
    args.check_checktiff = false;
    args.check_ortho = false;
    args.check_RA_only = false;
    args.check_LSF   = false;
    args.check_LSF_DEM = false;
    args.check_LSFDEMpath = false;
    args.check_LSF2  = 0;
    args.check_Matchtag = false;
    args.check_EO = false;
    args.check_fl = false;
    args.check_ccd = false;
    args.check_full_cal = false;
    args.check_coreg = 0;     //image coreg = 1, DEM coreg = 2, image + DEM = 3
    args.check_sdm_ortho = 0; //no coreg = 1 , with coreg = 2
    args.check_DEM_coreg_output = false;
    args.check_txt_input = 0; //no txt input = 0, DEM coregistration txt input = 1;
    args.check_downsample = false;
    args.check_DS_txy = false;
    
    args.number_of_images = 2;
    args.projection = 3;//PS = 1, UTM = 2
    args.sensor_provider = DG; //DG = 1, Pleiades = 2, Planet = 3
    args.check_imageresolution = false;
    args.utm_zone = -99;
    args.ortho_count = 1;
    args.overlap_length = 50;
    args.RA_only = 0;
    args.focal_length = 120;
    args.CCD_size = 12;
    args.SGM_py = 1;
    args.DS_sigma = 1.6;
    args.DS_kernel = 9;
    args.GCP_spacing = -9;
    
    TransParam param;
    param.bHemisphere = 1;
    
    int image_count = 0;
    args.RA_line_count = 0;
    args.RA_sample_count = 0;
    
    args.ra_line[0] = 0.0;
    args.ra_sample[0] = 0.0;
    args.System_memory = getSystemMemory(); //default system memory
    
    int DEM_divide = 0;
    double **Imageparams = NULL;
    
    if(argc == 1)
    {
        char save_filepath[500];
        char LeftImagefilename[500];
        
        args.check_arg = 0;
        
        Imageparams = (double**)calloc(MaxNCC, sizeof(double*));
        for(int ti = 0 ; ti < MaxNCC ; ti++)
        {
            Imageparams[ti] = (double*)calloc(sizeof(double),2);
            Imageparams[ti][0] = 0.0;
            Imageparams[ti][1] = 0.0;
        }
        
        DEM_divide = SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);
        
        char DEMFilename[500];
        char Outputpath[500];
        sprintf(DEMFilename, "%s/%s_dem.tif", save_filepath,args.Outputpath_name);
        sprintf(Outputpath, "%s", save_filepath);
            
        printf("%s\n",LeftImagefilename);
        printf("%s\n",DEMFilename);
        printf("%s\n",Outputpath);
        if(DEM_divide == 0)
            orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1,DEM_divide,Imageparams);
        else
            orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1,DEM_divide,Imageparams);
    }
    else if(argc == 2)
    {
        if (strcmp("-help",argv[1]) == 0) 
        {
            printf("Usage:./setsm [-help] : general explanation\n");
            printf("\texecution 1 : ./setsm : execute setsm with default.txt\n");
            printf("\t\texample usage : ./setsm\n");
            printf("\texecution 2 : ./setsm [image1] [image2] [output_directory]\n");
            printf("\t\t(execute setsm with image1, image2, and output directory for saving the results)\n");
            printf("\t\texample usage : ./setsm /home/image1.tif /home/image2.tif /home/output\n");
            printf("\t\t\t or ./setsm /home/image1.bin /home/image2.bin /home/output\n");
            printf("\texecution 3 : ./setsm [image1] [image2] [output_directory] [-options]\n");
            printf("\t\t(execute setsm with image1, image2 and output directory for saving the results with user-defined options\n");
            printf("\t\texample usage : ./setsm /home/image1.tif /home/image2.tif /home/output -outres 10 -threads 12 -seed /home/seed_dem.bin 50\n\n");
            
            printf("setsm version : %s\n", setsm_version);
            printf("supported image format : tif with xml, and binary with envi header file\n");
            printf("options\n");
            printf("\t[-outres value]\t: Output grid spacing[m] of Digital Elevation Model(DEM)\n");
            printf("\t[-seed filepath sigma]\t: Seed DEM(tif or binary format with envi header file for reducing computation loads with a height accuracy[m] of seed dem\n");
            printf("\t\t(if this option is set, setsm defines serach-spaces between + 1*sigma and - 1*sigma for calculating a height of grid position.\n");
            printf("\t\tThis option is not recommended on very changeable area. setsm can reconstruct 3D surface information without any seed dem\n");
            printf("\t[-boundary_min_X value1 -boundary_min_Y value2 -boundary_max_X value3 -boundary_max_Y value4]\t: Define specific DEM area to generate. The X and Y coordinate values should be Polar Stereographic or UTM\n");
            printf("\t[-tilesize value]\t: Set a tilesize for one time processing. Default is 8,000 pixels\n");
            printf("\t[-projection value]\t: Set planemetric coordinate projection. The value is 'ps' or 'utm'. Default projection is automatically defined by latitude of the input information of xml (between 60N and 60S is utm, and other is ps\n");
            printf("\t[-threads value]\t : Total number of threads for utilizing openmp parallel codes\n");
            printf("\t\t(if you don't know about this value, input '0'. Openmp can automatically detect a best value of your system)\n");
            printf("\t[-RAonly value]\t: If set to 1 (true), program will exit after RA calculation. Default = 0 (false)\n");
        }
    }
    else if(argc == 3)
    {
        if(strcmp("-gridonly",argv[1]) == 0)
        {
            char save_filepath[500];
            char LeftImagefilename[500];
            
            args.check_gridonly = true;
            
            sprintf(args.Outputpath,"%s",argv[2]);
            
            char *Outputpath_name  = SetOutpathName(args.Outputpath);
            sprintf(args.Outputpath_name,"%s",Outputpath_name);
            printf("after pathname %s\n",args.Outputpath_name);
            
            printf("%s\n",args.Outputpath);
            printf("%s\n", args.Outputpath_name);

            free(Outputpath_name);
            
            Imageparams = (double**)calloc(MaxNCC, sizeof(double*));
            for(int ti = 0 ; ti < MaxNCC ; ti++)
            {
                Imageparams[ti] = (double*)calloc(sizeof(double),2);
                Imageparams[ti][0] = 0.0;
                Imageparams[ti][1] = 0.0;
            }
            
            DEM_divide = SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);
        }
    }
    else if(argc == 4)
    {
        args.check_arg = 1;
        sprintf(args.Image[0],"%s",argv[1]);
        sprintf(args.Image[1],"%s",argv[2]);
        sprintf(args.Outputpath,"%s",argv[3]);
        
        char *Outputpath_name  = SetOutpathName(args.Outputpath);
        sprintf(args.Outputpath_name,"%s",Outputpath_name);
        printf("after pathname %s\n",args.Outputpath_name);
        
        printf("%s\n",args.Image[0]);
        printf("%s\n",args.Image[1]);
        printf("%s\n",args.Outputpath);
        printf("%s\n", args.Outputpath_name);

        free(Outputpath_name);
        
        char save_filepath[500];
        char LeftImagefilename[500];
        
        if( strcmp(args.Image[0],args.Image[1]) != 0)
        {
            
            Imageparams = (double**)calloc(MaxNCC, sizeof(double*));
            for(int ti = 0 ; ti < MaxNCC ; ti++)
            {
                Imageparams[ti] = (double*)calloc(sizeof(double),2);
                Imageparams[ti][0] = 0.0;
                Imageparams[ti][1] = 0.0;
            }
            
            DEM_divide = SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);
            
            char DEMFilename[500];
            char Outputpath[500];
            
            sprintf(Outputpath, "%s", save_filepath);
            if(DEM_divide == 0)
            {
                sprintf(DEMFilename, "%s/%s_dem.tif", save_filepath,args.Outputpath_name);
                orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1,DEM_divide,Imageparams);
            }
            else
            {
                for(int iter = 1 ; iter <= DEM_divide ; iter++)
                {
                    sprintf(DEMFilename, "%s/%s_%d_dem.tif", save_filepath,args.Outputpath_name,iter);
                    orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1,iter,Imageparams);
                }
            }
        }
        else
            printf("Please check input 1 and input 2. Both is same\n");
            
    }
    else if(argc > 4)
    {
        bool cal_flag = true;
        
        for (i=0; i<argc; i++)
        {
            if (strcmp("-LSFDEM",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input '1' for Yes or '0' for No\n");
                    cal_flag = false;
                }
                else
                {
                    args.check_LSF_DEM = atoi(argv[i+1]);
                    printf("LSFDEM %d\n",args.check_LSF_DEM);
                }
            }
            
            if (strcmp("-LSFDEM_path",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input DEM path for LSF\n");
                    cal_flag = false;
                }
                else
                {
                    sprintf(args.Outputpath,"%s",argv[i+1]);
                    args.check_LSFDEMpath = true;
                    printf("LSFDEMpath %s\n",args.Outputpath);
                }
            }
            
            if (strcmp("-coreg",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input '1' for Orthoimage or '2' for DEM\n");
                    cal_flag = false;
                }
                else
                {
                    args.check_coreg = atoi(argv[i+1]);
                    printf("Coregistration %d\n",args.check_coreg);
                }
            }
            
            if (strcmp("-txt_input",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input image txt file\n");
                    cal_flag = false;
                }
                else
                {
                    args.check_txt_input = 1;
                    printf("txt input %d\n",args.check_txt_input);
                    if(args.check_txt_input == 1)
                    {
                        sprintf(args.DEM_input_file,"%s",argv[i+1]);
                        printf("txt input %s\n",args.DEM_input_file);
                    }
                }
            }
            
            if (strcmp("-coreg_output",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input '2' for generation both coregistered image and whole image diff stat or '1' for generation whole image diff stat or '0' for no generation\n");
                    cal_flag = false;
                }
                else
                {
                    if(args.check_coreg == 2 || args.check_coreg == 3)
                    {
                        args.check_DEM_coreg_output = atoi(argv[i+1]);
                        printf("Coregistration output tif gneneration %d\n",args.check_DEM_coreg_output);
                    }
                }
            }
            
            if (strcmp("-gcp_space",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input the GCP space value\n");
                    cal_flag = false;
                }
                else
                {
                    args.GCP_spacing = atof(argv[i+1]);
                    printf("%f\n",args.GCP_spacing);
                }
            }
            
            if (strcmp("-SDM",argv[i]) == 0 || strcmp("-sdm",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input '1' for Orthoimage only without coreg, '2' for orthoimage with coreg, '3' for DEM with original image\n");
                    cal_flag = false;
                }
                else
                {
                    args.check_sdm_ortho = atoi(argv[i+1]);
                    printf("SDM with orthoimage %d\n",args.check_sdm_ortho);
                }
            }
            
            if (strcmp("-outres",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input the outres value\n");
                    cal_flag = false;
                }
                else
                {
                    args.DEM_space = atof(argv[i+1]);
                    printf("%f\n",args.DEM_space);
                    args.check_DEM_space = true;
                }
            }
            
            if (strcmp("-projection", argv[i]) == 0) {
                if (argc == i + 1) {
                    printf("Please input Projection info\n");
                    cal_flag = false;
                } else {
                    if(strcmp("utm",argv[i+1]) == 0 || strcmp("UTM",argv[i+1]) == 0)
                    {
                        args.projection = 2;
                        printf("UTM projection \n");
                        
                    }
                    else
                    {
                        args.projection = 1;
                        printf("PS projection\n");
                    }
                }
            }
            
            if (strcmp("-North",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input North value\n");
                    cal_flag = false;
                }
                else
                {
                    int dir = atoi(argv[i+1]);
                    if(dir == 1)
                    {
                        param.bHemisphere = 1;
                        param.pm = 1;
                    }
                    else
                    {
                        param.bHemisphere = 2;
                        param.pm = -1;
                    }
                }
            }
            
            if (strcmp("-mem",argv[i]) == 0)
            {
                if (argc == i+1) {
                    printf("Please input System memory\n");
                    cal_flag = false;
                }
                else
                {
                    args.System_memory = atof(argv[i+1]);
                    printf("System memory %f\n",args.System_memory);
                }
            }
        }
        
        if(args.check_LSF_DEM)
        {
            char str_DEMfile[1000];
            char str_matchfile[1000];
            char str_matchfile_tif[1000];
            
            char str_smooth_file[500];
            char DEM_header[500];
            char smooth_GEOTIFF_filename[500];
            char result_file[500];
            char metafilename[500];
            
            if(args.projection == 3)
                param.projection = 1;
            else
                param.projection = args.projection;
            
            char *ext = strrchr(args.Outputpath,'.');
            if(ext)
            {
                if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1) || !strcmp("raw",ext+1))
                {
                    printf("DEM file open\n");
                    
                    sprintf(str_DEMfile, "%s", args.Outputpath);
                    char *tmp_chr = remove_ext(args.Outputpath);
                    
                    sprintf(str_smooth_file, "%s_smooth.raw", tmp_chr);
                    sprintf(DEM_header, "%s_smooth.hdr", tmp_chr);
                    sprintf(smooth_GEOTIFF_filename, "%s_smooth.tif", tmp_chr);
                    
                    int full_size = strlen(tmp_chr);
                    char *t_name = (char*)malloc(sizeof(char)*(full_size-4));
                    for(int kk = 0; kk < full_size-4 ; kk++)
                        t_name[kk] = tmp_chr[kk];
           
                    sprintf(metafilename,"%s_meta.txt",t_name);
                    sprintf(str_matchfile,"%s_matchtag.raw",t_name);
                    sprintf(str_matchfile_tif,"%s_matchtag.tif",t_name);
                    sprintf(result_file,"%s_smooth_result.txt",t_name);

                    free(tmp_chr);
                }
                else
                {
                    printf("No DEM(tif or raw) exists!!\n");
                    exit(1);
                }
            }
            else
            {
                char *Outputpath_name  = SetOutpathName(args.Outputpath);
                
                sprintf(str_DEMfile, "%s/%s_dem.tif", args.Outputpath,Outputpath_name);
                sprintf(str_smooth_file,"%s/%s_dem_smooth.raw",args.Outputpath,Outputpath_name);
                sprintf(DEM_header, "%s/%s_dem_smooth.hdr", args.Outputpath,Outputpath_name);
                sprintf(smooth_GEOTIFF_filename, "%s/%s_dem_smooth.tif", args.Outputpath, Outputpath_name);
                sprintf(metafilename,"%s/%s_meta.txt",args.Outputpath,Outputpath_name);
                sprintf(str_matchfile,"%s/%s_matchtag.raw",args.Outputpath,Outputpath_name);
                sprintf(str_matchfile_tif,"%s/%s_matchtag.tif",args.Outputpath,Outputpath_name);
                sprintf(result_file,"%s/%sdem_smooth_result.txt",args.Outputpath,Outputpath_name);

                free(Outputpath_name);
            }
            
            printf("dem file %s\n",str_DEMfile);
            printf("metafilename %s\n",metafilename);
            printf("matchfile %s\n",str_matchfile);
            
            printf("smooth DEM %s\n",str_smooth_file);
            printf("DEM_header %s\n",DEM_header);
            printf("smooth_GEOTIFF_filename %s\n", smooth_GEOTIFF_filename);
            printf("result file %s\n",result_file);
            
            FILE *pFile_DEM = fopen(str_DEMfile,"r");
            printf("check exist %s %d\n",str_DEMfile,!!pFile_DEM);
            
            if(pFile_DEM)
            {
                double MPP_stereo_angle = 1;
                FILE *pMetafile   = fopen(metafilename,"r");
                if(pMetafile)
                {
                    printf("meta file exist!!");
                    
                    char linestr[1000];
                    char linestr1[1000];
                    char* pos1;
                    char* token = NULL;
                    bool check_mpp = false;
                    
                    while(!feof(pMetafile) && !check_mpp)
                    {
                        fgets(linestr,sizeof(linestr),pMetafile);
                        strcpy(linestr1,linestr);
                        token = strtok(linestr,"=");
                        
                        if(strcmp(token,"Setereo_pair_expected_height_accuracy") == 0)
                        {
                            pos1 = strstr(linestr1,"=")+1;
                            MPP_stereo_angle            = atof(pos1);
                            
                            check_mpp = true;
                        }
                    }
                    fclose(pMetafile);
                }
                else
                    printf("meta file doesn't exist!!\n");
                
                printf("MPP_stereo_angle %f\n",MPP_stereo_angle);
                
                double minX, maxY;
                double grid_size = args.DEM_space;
                CSize DEM_size = GetDEMsize(str_DEMfile,metafilename,&param,&grid_size,&minX,&maxY);
                float *seeddem = GetDEMValue(str_DEMfile,DEM_size);
                
                long data_length = (long)DEM_size.width*(long)DEM_size.height;
                printf("data_length %ld\n",data_length);
                printf("projection %d\tHemisphere %d\n",param.projection,param.bHemisphere);
                
                printf("%d\n",DEM_size.width);
                printf("%d\n",DEM_size.height);
                printf("%f\n",grid_size);
                
                LSFINFO *Grid_info = (LSFINFO*)malloc(sizeof(LSFINFO)*data_length);
                if(Grid_info == NULL)
                {
                    printf("Insufficient memory available\n");
                    exit(1);
                }
                
                for(long count_index = 0 ; count_index < data_length; count_index++)
                    Grid_info[count_index].lsf_kernel = 2;
                
                float *smooth_DEM = (float*)malloc(sizeof(float)*data_length);
                if(smooth_DEM == NULL)
                {
                    printf("Insufficient memory available\n");
                    exit(1);
                }
                
                double max_std = -100000;
                int max_std_iter = -1;
                int min_std_iter = 100;
                double min_std = 100000;
                int max_iter_th;
                if(grid_size >= 2)
                    max_iter_th = 2;
                else
                    max_iter_th = 2;
                
                const int max_iter_count = 5;
                int s_iter = 0;
                bool check_smooth_iter = true;
                
                double sigma_avg = 10000;
                double sigma_std = 10000;
                
                while(check_smooth_iter && s_iter < max_iter_count)
                {
                    printf("start LSF\n");
                    int selected_numpts;
                    
                    if((sigma_avg < 0.5 && sigma_std < 1) || s_iter == max_iter_count-1)
                    {
                        if(s_iter > max_iter_th)
                        {
                            printf("final local suface fitting\n");
                            check_smooth_iter = false;
                        }
                    }
                    
                    DEM_STDKenel_LSF(Grid_info, &sigma_avg, &sigma_std,seeddem,smooth_DEM,grid_size,s_iter,DEM_size, MPP_stereo_angle);
                    
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
                    memcpy(seeddem,smooth_DEM,sizeof(float)*data_length);
                    
                    s_iter++;
                }
                free(smooth_DEM);
                free(Grid_info);
                
                WriteGeotiff(smooth_GEOTIFF_filename, seeddem, DEM_size.width, DEM_size.height, grid_size, minX, maxY, param.projection, param.utm_zone, param.bHemisphere, 4);
                
                FILE *presult = fopen(result_file,"w");
                fprintf(presult,"%d\t%f\t%d\t%f\n",max_std_iter,max_std,min_std_iter,min_std);
                fclose(presult);
                
                printf("%d\t%f\t%d\t%f\n",max_std_iter,max_std,min_std_iter,min_std);
                
                free(seeddem);
      
                fclose(pFile_DEM);
            }
        }
        else
        {
            args.check_arg = 1;
            
            bool bminx  = false;
            bool bmaxx  = false;
            bool bminy  = false;
            bool bmaxy  = false;
            
            for (i=0; i<argc; i++)
            {
                if (strcmp("-SGMpy",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input sgm pyramid level steps (default is -1)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.SGM_py = atoi(argv[i+1]);
                        printf("Steps of sgm pyramid level %d\n",args.SGM_py);
                    }
                }
                
                if (strcmp("-PL",argv[i]) == 0 || strcmp("-pl",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input pyramid level steps (default is 4)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.pyramid_level = atoi(argv[i+1]);
                        printf("Steps of pyramid level %d\n",args.pyramid_level);
                    }
                }
                
                if (strcmp("-SDM_SS",argv[i]) == 0 || strcmp("-sdm_ss",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input SDM search size (default is 3)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.SDM_SS = atoi(argv[i+1]);
                        printf("SDM search size %d\n",args.SDM_SS);
                    }
                }
                
                if (strcmp("-SDM_DAYS",argv[i]) == 0 || strcmp("-sdm_days",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input time difference(day) between images (default is 1)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.SDM_days = atof(argv[i+1]);
                        printf("SDM time gap %f(days)\n",args.SDM_days);
                    }
                }
                
                if (strcmp("-SDM_AS",argv[i]) == 0 || strcmp("-sdm_as",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input average velocity per day(m/day) for SDM\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.SDM_AS = atof(argv[i+1]);
                        printf("Average velocity %f(m/day)\n",args.SDM_AS);
                    }
                }
                
                if (strcmp("-PC",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please select option of full resolution DEM grid (default is 0)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.check_full_cal = atoi(argv[i+1]);
                        printf("Option of full resolution DEM grid %d\n",args.check_full_cal);
                    }
                }
                
                if (strcmp("-Sensor",argv[i]) == 0 || strcmp("-sensor",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input '1' for RFM(default) or '0' for Collinear Equation(Frame)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        int temp;
                        temp = atoi(argv[i+1]);
                        if(temp == 1)
                        {
                            args.sensor_type = SB;
                            printf("Spaceborne RFM Sensor %d\n",args.sensor_type);
                        }
                        else
                        {
                            args.sensor_type = AB;
                            printf("Airborne frame Sensor %d\n",args.sensor_type);
                        }
                    }
                }
                
                if (strcmp("-NImages",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input Number of images (default is 2 for stereo)\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.number_of_images = atoi(argv[i+1]);
                        printf("Number of Images %d\n",args.number_of_images);
                    }
                }
                
                if (strcmp("-FL",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input focal length[mm]\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.focal_length = atof(argv[i+1]);
                        printf("focal length %f\n",args.focal_length);
                        args.check_fl = true;
                    }
                }
                
                if (strcmp("-CCDsize",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input CCD size[um]\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.CCD_size = atof(argv[i+1]);
                        printf("CCD size %f\n",args.CCD_size);
                        args.check_ccd = true;
                    }
                }
                
                if (strcmp("-EO",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input Exterior Orientation informations\n");
                        cal_flag = false;
                    }
                    else
                    {
                        sprintf(args.EO_Path,"%s",argv[i+1]);
                        args.check_EO = true;
                        printf("EO_path %s\n",args.EO_Path);
                    }
                }
                
                if (strcmp("-image",argv[i]) == 0 || strcmp("-Image",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input image path\n");
                        cal_flag = false;
                    }
                    else
                    {
                        sprintf(args.Image[image_count],"%s",argv[i+1]);
                        printf("image%d %s\n",image_count,args.Image[image_count]);
                        
                        image_count++;
                    }
                }
                
                if (strcmp("-outpath",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input outpath\n");
                        cal_flag = false;
                    }
                    else
                    {
                        sprintf(args.Outputpath,"%s",argv[i+1]);
                        printf("Out path %s\n",args.Outputpath);
                    }
                }
                
                if (strcmp("-provider", argv[i]) == 0) {
                    if (argc == i + 1) {
                        printf("Please input Provider info\n");
                        cal_flag = false;
                    } else {
                        if(strcmp("DG",argv[i+1]) == 0 || strcmp("dg",argv[i+1]) == 0)
                        {
                            args.sensor_provider = DG;
                            printf("Image Provider : Digital Globe\n");
                            
                        }
                        else if(strcmp("Pleiades",argv[i+1]) == 0 || strcmp("pleiades",argv[i+1]) == 0)
                        {
                            args.sensor_provider = PL;
                            printf("Image Provider : Pleiades\n");
                        }
                        else if(strcmp("Planet",argv[i+1]) == 0 || strcmp("planet",argv[i+1]) == 0)
                        {
                            args.sensor_provider = PT;
                            printf("Image Provider : Planet\n");
                        }
                        else
                        {
                            args.projection = 1;    
                            printf("Not suppoted Provider. Please use Digital Globe (Worldview, QuickBird, GeoEye) or Pleiades\n");
                            exit(1);
                        }
                    }
                }
                
                if (strcmp("-GSD",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input the image resolution (GSD) value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.image_resolution = atof(argv[i+1]);
                        printf("%f\n",args.image_resolution);
                        args.check_imageresolution = true;
                    }
                }
                
                if (strcmp("-LSF",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input '1 ' for applying LSF smoothing\n ");
                        cal_flag = false;
                    }
                    else
                    {
                        int input_number = atoi(argv[i+1]);
                        if(input_number == 1 || input_number == 2)
                        {
                            args.check_LSF2 = input_number;
                            printf("LSF %d\n",args.check_LSF2);
                        }
                        else
                            printf("LSF is not applied!!\n");
                    }
                }
                
                if (strcmp("-MT",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input '1 ' for applying v331 matchtag info\n ");
                        cal_flag = false;
                    }
                    else
                    {
                        int input_number = atoi(argv[i+1]);
                        if(input_number == 1)
                        {
                            args.check_Matchtag = true;
                            printf("MT %d\n",args.check_Matchtag);
                        }
                        else
                            printf("MT is not applied!!\n");
                    }
                }
                
                if (strcmp("-outres",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the outres value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.DEM_space = atof(argv[i+1]);
                        printf("%f\n",args.DEM_space);
                        args.check_DEM_space = true;
                    }
                }
                
                if (strcmp("-threads",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the threads value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.Threads_num = atoi(argv[i+1]);
                        printf("%d\n",args.Threads_num);
                        args.check_Threads_num = true;
                    }
                }
     
                if (strcmp("-downsample",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input the threads value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        sprintf(args.seedDEMfilename,"%s",argv[i+1]);
                        printf("source %s\n",args.seedDEMfilename);
                        sprintf(args.Outputpath_name,"%s",argv[i+2]);
                        printf("target %s\n",args.Outputpath_name);
                        args.check_downsample = true;
                    }
                }
                
                if (strcmp("-txy",argv[i]) == 0)
                {
                    args.DS_tx = atof(argv[i+1]);
                    printf("Target origin X %f\n",args.DS_tx);
                    args.DS_ty = atof(argv[i+2]);
                    printf("Target origin Y %f\n",args.DS_ty);
                    args.check_DS_txy = true;
                }
                
                if (strcmp("-kernel",argv[i]) == 0)
                {
                    args.DS_kernel = atoi(argv[i+1]);
                    printf("kernel size %d\n",args.DS_kernel);
                }
                
                if (strcmp("-sigma",argv[i]) == 0)
                {
                    args.DS_sigma = atof(argv[i+1]);
                    printf("sigma value %f\n",args.DS_sigma);
                }
                
                if (strcmp("-seed",argv[i]) == 0)
                {
                    if(args.check_sdm_ortho == 0)
                    {
                        TIFF *tif = NULL;
                        
                        if (argc == i+1) {
                            printf("Please input the seed filepath\n");
                            cal_flag = false;
                        }
                        else
                        {
                            int str_size = strlen(argv[i+1]);
                            int kk=0;
                            
                            for(kk=0;kk<str_size;kk++)
                            {
                                args.seedDEMfilename[kk] = argv[i+1][kk];
                            }
                            args.seedDEMfilename[str_size] = '\0';
                            
                            printf("%s\n",args.seedDEMfilename);
                            
                            printf("DEM check\n");
                            tif  = TIFFOpen(args.seedDEMfilename,"r");
                            printf("first check\n");

                            if(tif)
                            {
                                printf("DEM is exist!! \n");
                                
                                printf("%s\n",args.seedDEMfilename);
                                args.check_seeddem = true;
                            }
                            else {
                                printf("second check\n");
                                char* temp_path = remove_ext(args.seedDEMfilename);
                                
                                sprintf(args.seedDEMfilename,"%s.tif",temp_path);
                                
                                tif  = TIFFOpen(args.seedDEMfilename,"r");
                                if(tif)
                                {
                                    printf("%s\n",args.seedDEMfilename);
                                    args.check_seeddem = true;
                                }
                                else
                                {
                                    sprintf(args.seedDEMfilename,"%s.raw",temp_path);
                                    printf("%s\n",args.seedDEMfilename);
                                    FILE *pfile = fopen(args.seedDEMfilename,"r");
                                    if(pfile)
                                    {
                                        printf("DEM is exist!! \n");
                                        
                                        printf("%s\n",args.seedDEMfilename);
                                        args.check_seeddem = true;
                                    }
                                    else
                                    {
                                        printf("DEM doesn't exist!! Please check DEM path\n");
                                        args.check_seeddem = false;
                                        exit(0);
                                    }
                                }

                                free(temp_path);
                            }

                            if(args.check_seeddem)
                            {
                                char* temp_path = remove_ext(args.seedDEMfilename);
                                printf("seedem %s\n",temp_path);
                                int full_size;
                                full_size       = strlen(args.seedDEMfilename);
                                char* Metafile1 = (char*)malloc(sizeof(char)*(full_size-8+1));
                                char Metafile[500];
                                int i;
                                for(i=0;i<full_size-8;i++)
                                    Metafile1[i] = args.seedDEMfilename[i];
                                Metafile1[full_size-8] = '\0';
                                
                                for(i=0;i<full_size-8;i++)
                                    Metafile[i] = args.seedDEMfilename[i];
                                Metafile[full_size-8] = '\0';
                                
                                printf("%s\n",Metafile);
                                char *str = (char*)"meta.txt";
                                
                                sprintf(args.metafilename,"%s_%s",Metafile,str);
                                printf("Meta file %s\n",args.metafilename);
                                
                                FILE* pFile_meta;
                                pFile_meta  = fopen(args.metafilename,"r");
                                if(pFile_meta)
                                {
                                    printf("meta file loading successful\n");
                                    printf("Meta file %s\n",args.metafilename);
                                    fclose(pFile_meta);
                                }
                                else
                                {
                                    printf("%s\n",Metafile1);
                                    
                                    sprintf(args.metafilename,"%s_%s",Metafile1,str);
                                    FILE* pFile_meta1;
                                    pFile_meta1 = fopen(args.metafilename,"r");
                                    if(pFile_meta1)
                                    {
                                        printf("meta file loading successful\n");
                                        printf("Meta file %s\n",args.metafilename);
                                        fclose(pFile_meta1);
                                    }
                                    else
                                    {
                                        printf("meta file loading failed\n");
                                        args.check_seeddem = false;
                                    }
                                }
                                
                                free(Metafile1);
                                free(temp_path);
                            }
                            
                        }
                        
                        if (argc == i+2) {
                            printf("Please input the seed sigma value\n");
                            cal_flag = false;
                        }
                        else
                        {
                            args.seedDEMsigma = atof(argv[i+2]);
                            printf("%f\n",args.seedDEMsigma);
                            if(tif)
                                args.check_seeddem = true;
                        }
                    }
                    else
                    {
                        int str_size = strlen(argv[i+1]);
                        int kk=0;
                        
                        for(kk=0;kk<str_size;kk++)
                        {
                            args.seedDEMfilename[kk] = argv[i+1][kk];
                        }
                        args.seedDEMfilename[str_size] = '\0';
                        
                        printf("%s\n",args.seedDEMfilename);
                        args.check_seeddem = true;
                    }
                }
                
                if (strcmp("-tilesSR",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the start row tile value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.start_row = atoi(argv[i+1]);
                        printf("%d\n",args.start_row);
                        args.check_tiles_SR = true;
                    }
                }
                
                if (strcmp("-tilesER",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the end row tile value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        
                        args.end_row = atoi(argv[i+1]);
                        printf("%d\n",args.end_row);
                        args.check_tiles_ER = true;
                    }
                }
                
                if (strcmp("-tilesSC",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the start col tile value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.start_col = atoi(argv[i+1]);
                        printf("%d\n",args.start_col);
                        args.check_tiles_SC = true;
                    }
                }
                
                if (strcmp("-tilesEC",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the end col tile value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.end_col = atoi(argv[i+1]);
                        printf("%d\n",args.end_col);
                        args.check_tiles_EC = true;
                    }
                }
                
                if (strcmp("-minH",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the min Height value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.minHeight = atof(argv[i+1]);
                        printf("%f\n",args.minHeight);
                        args.check_minH = true;
                    }
                }
                
                if (strcmp("-maxH",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the max Height value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.maxHeight = atof(argv[i+1]);
                        printf("%f\n",args.maxHeight);
                        args.check_maxH = true;
                    }
                }
                
                if (strcmp("-RAline",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the RA line value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.ra_line[args.RA_line_count] = atof(argv[i+1]);
                        printf("%f\n",args.ra_line[args.RA_line_count]);
                        args.check_RA_line = true;
                        
                        args.RA_line_count++;
                    }
                }
                
                if (strcmp("-RAsample",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the RA sample value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.ra_sample[args.RA_sample_count] = atof(argv[i+1]);
                        printf("%f\n",args.ra_sample[args.RA_sample_count]);
                        args.check_RA_sample = true;
                        
                        args.RA_sample_count ++;
                    }
                }
                
                if (strcmp("-RAtileR",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the RA tileR value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.RA_row = atof(argv[i+1]);
                        printf("%d\n",args.RA_row);
                        args.check_RA_tileR = true;
                    }
                }
                
                if (strcmp("-RAtileC",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input the RA tileC value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.RA_col = atof(argv[i+1]);
                        printf("%d\n",args.RA_col);
                        args.check_RA_tileC = true;
                    }
                }

                if (strcmp("-RAonly",argv[i]) == 0)
                {
                    if (argc == i+1 || (atoi(argv[i+1]) != 0 && atoi(argv[i+1]) != 1)) {
                        printf("Please input the RA_only 0 or 1\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.RA_only = atoi(argv[i+1]);
                        printf("%d\n",args.RA_only);
                        args.check_RA_only = true;
                    }
                }

                if (strcmp("-tilesize",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input tile size value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.tilesize = atoi(argv[i+1]);
                        if(args.tilesize > 20000)
                            args.tilesize = 100000;
                        
                        printf("%d\n",args.tilesize);
                        args.check_tilesize = true;
                    }
                }
                
                if (strcmp("-LOO",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input length of overlapped area value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.overlap_length = atof(argv[i+1]);
                        printf("%f\n",args.overlap_length);
                    }
                }
                
                
                if (strcmp("-boundary_min_X",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input boundary min_X value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.Min_X = atof(argv[i+1]);
                        printf("%f\n",args.Min_X);
                        
                        bminx   = true;
                    }
                }
                if (strcmp("-boundary_min_Y",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input boundary min_Y value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.Min_Y = atof(argv[i+1]);
                        printf("%f\n",args.Min_Y);
                        
                        bminy   = true;
                    }
                }
                if (strcmp("-boundary_max_X",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input boundary max_X value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.Max_X = atof(argv[i+1]);
                        printf("%f\n",args.Max_X);
                        
                        bmaxx   = true;
                    }
                }
                if (strcmp("-boundary_max_Y",argv[i]) == 0) 
                {
                    if (argc == i+1) {
                        printf("Please input boundary max_Y value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.Max_Y = atof(argv[i+1]);
                        printf("%f\n",args.Max_Y);
                        
                        bmaxy   = true;
                    }
                }
                
                if (bminx && bmaxx && bminy && bmaxy)
                    args.check_boundary = true;
                
                if(strcmp("-checktiff",argv[i]) == 0)
                {
                    args.check_checktiff = true;
                }
                
                if (strcmp("-projection", argv[i]) == 0) {
                    if (argc == i + 1) {
                        printf("Please input Projection info\n");
                        cal_flag = false;
                    } else {
                        if(strcmp("utm",argv[i+1]) == 0 || strcmp("UTM",argv[i+1]) == 0)
                        {
                            args.projection = 2;
                            printf("UTM projection \n");
                            
                        }
                        else
                        {
                            args.projection = 1;    
                            printf("PS projection\n");
                        }
                    }
                }
                
                if (strcmp("-utm_zone",argv[i]) == 0)
                {
                    if (argc == i+1) {
                        printf("Please input utm_zome value\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.utm_zone = atoi(argv[i+1]);
                        printf("%d\n",args.utm_zone);
                    }
                }
                
                
                if (strcmp("-ortho",argv[i]) == 0)
                {
                    args.check_ortho = true;
                    printf("ortho only\n");
                    
                    if (argc == i + 1) {
                        printf("Please input ortho # 1 or 2\n");
                        cal_flag = false;
                    }
                    else
                    {
                        args.ortho_count = atoi(argv[i+1]);
                        printf("%d\n",args.ortho_count);
                    }
                }
            }
            
            if(cal_flag)
            {
                char save_filepath[500];
                char LeftImagefilename[500];
                
                bool check_frame_info = true;
                
                if(args.check_txt_input == 0)
                {
                    if(args.check_downsample)
                    {
                        DownSample(args);
                    }
                    else if(image_count > 1 || args.check_ortho)
                    {
                        args.number_of_images = image_count;
                        
                        char *Outputpath_name  = SetOutpathName(args.Outputpath);
                        sprintf(args.Outputpath_name,"%s",Outputpath_name);
                        printf("after pathname %s\n",args.Outputpath_name);
                        
                        printf("%s\n",args.Outputpath);
                        printf("%s\n", args.Outputpath_name);

                        free(Outputpath_name);

                        
                        Imageparams = (double**)calloc(MaxNCC, sizeof(double*));
                        for(int ti = 0 ; ti < MaxNCC ; ti++)
                        {
                            Imageparams[ti] = (double*)calloc(sizeof(double),2);
                            Imageparams[ti][0] = 0.0;
                            Imageparams[ti][1] = 0.0;
                        }
                        
                        if(args.check_checktiff)
                        {
                            DEM_divide = SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);
                        }
                        else if( strcmp(args.Image[0],args.Image[1]) != 0)
                        {
                            DEM_divide = SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);

                            char DEMFilename[500];
                            char Outputpath[500];
                            
                            sprintf(Outputpath, "%s", save_filepath);
                            
                            printf("param %s %d %d\n", param.direction,param.utm_zone,param.projection);
                            if(args.projection != 3)
                                param.projection = args.projection;
                            
                            printf("imageparams %f\t%f\t%f\t%f\n",Imageparams[0][0],Imageparams[0][1],Imageparams[1][0],Imageparams[1][1]);
                            
                            if(DEM_divide == 0)
                            {
                                //if(!args.check_ortho)
                                    sprintf(DEMFilename, "%s/%s_dem.tif", save_filepath,args.Outputpath_name);
                                //else
                                //    sprintf(DEMFilename, "%s", args.seedDEMfilename);
                                
                                if(!args.check_Matchtag)
                                {
                                    orthogeneration(param,args,args.Image[0], DEMFilename, Outputpath,1,DEM_divide,Imageparams);
                                    orthogeneration(param,args,args.Image[1], DEMFilename, Outputpath,2,DEM_divide,Imageparams);
                                }
                                //else if(args.ortho_count == 2)
                                //    orthogeneration(param,args,args.Image[1], DEMFilename, Outputpath,2,DEM_divide,Imageparams);
                                
                                if(args.check_LSF2 == 2)
                                    remove(DEMFilename);
                            }
                            else
                            {
                                for(int iter = 1 ; iter <= DEM_divide ; iter++)
                                {
                                    sprintf(DEMFilename, "%s/%s_%d_dem.tif", save_filepath,args.Outputpath_name,iter);
                                    if(!args.check_Matchtag)
                                    {
                                        orthogeneration(param,args,args.Image[0], DEMFilename, Outputpath,1,iter,Imageparams);
                                        orthogeneration(param,args,args.Image[1], DEMFilename, Outputpath,2,iter,Imageparams);
                                    }
                                    //else if(args.ortho_count == 2)
                                    //    orthogeneration(param,args,args.Image[1], DEMFilename, Outputpath,2,iter,Imageparams);
                                    
                                    if(args.check_LSF2 == 2)
                                        remove(DEMFilename);
                                }
                            }
                             
                        }
                        else
                            printf("Please check input 1 and input 2. Both is same\n");
                    }
                    else
                    {
                        printf("Plese check input images\n");
                    }
                }
                else if(args.check_txt_input == 1)
                {
                    char *Outputpath_name  = SetOutpathName(args.Outputpath);
                    sprintf(args.Outputpath_name,"%s",Outputpath_name);
                    printf("after pathname %s\n",args.Outputpath_name);
                    
                    printf("%s\n",args.Outputpath);
                    printf("%s\n", args.Outputpath_name);

                    free(Outputpath_name);
                    
                    SETSMmainfunction(&param,projectfilename,args,save_filepath,Imageparams);
                }
                    
            }
        }
    }
    
    printf("# of allocated threads = %d\n",omp_get_max_threads());
            
    
    if(Imageparams)
    {
        for(int ti = 0 ; ti < MaxNCC ; ti++)
        {
            free(Imageparams[ti]);
        }
    }
    free(Imageparams);
    
    return 0;
}

void DownSample(ARGINFO &args)
{
    long cols[2];
    long rows[2];
    CSize data_size;
    double minX, maxY, grid_size;
    
    CSize seeddem_size = ReadGeotiff_info(args.seedDEMfilename, &minX, &maxY, &grid_size);
    int downsample_step = ceil(log2(args.DEM_space/grid_size));
    if(downsample_step < 1)
    {
        printf("No downsampling necessary. Please check source(%3.1f) and target resolution(%3.1f)!!\n",grid_size,args.DEM_space);
        exit(1);
    }
    else
    {
        TransParam param;
        SetTranParam_fromGeoTiff(&param,args.seedDEMfilename);
        
        
        CSize *Imagesize = (CSize*)malloc(sizeof(CSize));
        Imagesize->width = seeddem_size.width;
        Imagesize->height = seeddem_size.height;
        
        cols[0] = 0;
        cols[1] = seeddem_size.width;
        
        rows[0] = 0;
        rows[1] = seeddem_size.height;
        
        float type(0);
        float *seeddem = Readtiff_T(args.seedDEMfilename,Imagesize,cols,rows,&data_size,type);
        
        float **pyimg = (float**)malloc(sizeof(float*)*downsample_step);
        
        CSize out_size;
        
        for(int level = 0 ; level < downsample_step ; level ++)
        {
            if(level == 0)
            {
                pyimg[0] = CreateImagePyramid(seeddem,seeddem_size,args.DS_kernel,args.DS_sigma);
                free(seeddem);
                out_size.width = seeddem_size.width/2;
                out_size.height = seeddem_size.height/2;
            }
            else
            {
                pyimg[level] = CreateImagePyramid(pyimg[level-1],out_size,args.DS_kernel,args.DS_sigma);
                out_size.width = out_size.width/2;
                out_size.height = out_size.height/2;
            }
        }
        /*
        for(int i=0 ; i< downsample_step-1 ; i++)
            free(pyimg[i]);
        
        WriteGeotiff(args.Outputpath_name, pyimg[downsample_step-1], out_size.width, out_size.height, args.DEM_space, minX, maxY, param.projection, param.zone, param.bHemisphere, 4);
        
        free(pyimg[downsample_step-1]);
        free(pyimg);
        */
        printf("Done Gaussian processing\n");
        
        D2DPOINT target;
        if(!args.check_DS_txy)
        {
            target.m_X = (int)(minX/2.0)*2;
            target.m_Y = (int)(maxY/2.0)*2;
        }
        else
        {
            target.m_X = (int)(args.DS_tx/2.0)*2;
            target.m_Y = (int)(args.DS_ty/2.0)*2;
        }
        
        D2DPOINT Dxy(minX - (int)target.m_X, maxY - (int)target.m_Y);
        D2DPOINT Dgrid(-Dxy.m_X/args.DEM_space,Dxy.m_Y/args.DEM_space);
        long data_length = (long)out_size.width*(long)out_size.height;
        
        float *outimg = (float*)malloc(sizeof(float)*data_length);
#pragma omp parallel for schedule(guided)
        for(long int iter_count = 0 ; iter_count < data_length ; iter_count++)
        {
            long int pts_row = (long int)(floor(iter_count/out_size.width));
            long int pts_col = iter_count % out_size.width;
            long int pt_index = pts_row*(long int)out_size.width + pts_col;
            
            D2DPOINT query_pt(pts_col + Dgrid.m_X,pts_row + Dgrid.m_Y);
            outimg[iter_count] = BilinearResampling(pyimg[downsample_step-1],out_size,query_pt);
        }
        printf("Done resampling\n");
        free(pyimg[downsample_step-1]);
        free(pyimg);
        
        WriteGeotiff(args.Outputpath_name, outimg, out_size.width, out_size.height, args.DEM_space, target.m_X, target.m_Y, param.projection, param.utm_zone, param.bHemisphere, 4);

        free(outimg);
         
    }
}

int SETSMmainfunction(TransParam *return_param, char* _filename, ARGINFO args, char *_save_filepath,double **Imageparams)
{
#ifdef BUILDMPI
    char a;
    char *pa = &a;
    char **ppa = &pa;
    int argc = 0;
    int provided = 1;
    MPI_Init_thread(&argc, &ppa, MPI_THREAD_FUNNELED, &provided);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      printf("MPI: Number of processes: %d\n", size);
    }
#endif

    int DEM_divide = 0;
    char computation_file[500];
    time_t total_ST = 0, total_ET = 0;
    double total_gap;
    FILE *time_fid;
    
    total_ST = time(0);
    
    bool cal_check;
    
    ProInfo *proinfo = new ProInfo;
    proinfo->number_of_images = args.number_of_images;
    proinfo->sensor_type = args.sensor_type;
    proinfo->sensor_provider = args.sensor_provider;
    proinfo->System_memory = args.System_memory;
    proinfo->pyramid_level = args.pyramid_level;
    proinfo->check_full_cal = args.check_full_cal;
    proinfo->SGM_py = args.SGM_py;
    sprintf(proinfo->save_filepath,"%s",args.Outputpath);
    printf("sgm level %d\t system memory %f\n",proinfo->SGM_py,proinfo->System_memory);
    
    if(args.check_ortho)
    {
        printf("number of images %d\n",proinfo->number_of_images);
        if(OpenProject(_filename,proinfo,args))
        {
            if(Maketmpfolders(proinfo))
            {
                TransParam param;
                char Outputpath[500];
                sprintf(Outputpath, "%s", args.Outputpath);
                double ***RPCs = (double***)calloc(proinfo->number_of_images, sizeof(double**));
                
                for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                {
                    if(args.sensor_provider == DG)
                    {
                        double row_grid_size, col_grid_size, product_grid_size;
                        BandInfo band;
                        RPCs[ti]       = OpenXMLFile(proinfo, ti, &row_grid_size, &col_grid_size, &product_grid_size, &band);
                    }
                    else if(args.sensor_provider == PL)
                    {
                        RPCs[ti]       = OpenXMLFile_Pleiades(proinfo->RPCfilename[ti]);
                    }
                    else if(args.sensor_provider == PT)
                    {
                        RPCs[ti]       = OpenXMLFile_Planet(proinfo->RPCfilename[ti]);
                    }
                }
                param.projection = args.projection;
                param.utm_zone   = args.utm_zone;
                
                double minLat, minLon;
                if(proinfo->sensor_type == SB)
                {
                    minLat      = RPCs[0][0][3];
                    minLon      = RPCs[0][0][2];
                    SetTransParam(minLat, minLon, &param);
                }
               
                printf("param projection %d\tzone %d\n",param.projection,param.utm_zone);
                char DEMFilename[500];
                sprintf(DEMFilename, "%s", args.seedDEMfilename);
                for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    orthogeneration(param,args,args.Image[ti], DEMFilename, Outputpath,1,0,Imageparams);

                for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    RPCsFree(RPCs[ti]);
                free(RPCs);
            }
        }
    }
    else if(args.check_coreg == 1)
    {
        D2DPOINT *adjust_sdt = (D2DPOINT*)calloc(sizeof(D2DPOINT),proinfo->number_of_images);
        double **Coreg_param = ImageCoregistration(return_param, _filename, args, _save_filepath, 1,adjust_sdt,&cal_check);
        if(cal_check)
            free(Coreg_param);
    }
    else if(args.check_coreg == 2)
        DEM_ImageCoregistration_GeomatricConstraint(return_param, _filename, args, _save_filepath, 2);
    else if(args.check_coreg == 3)
        DEM_ImageCoregistration_hillshade(return_param, _filename, args, _save_filepath, 2);
    else if(args.check_sdm_ortho == 1)
    {
        double** Coreg_param = (double**)calloc(sizeof(double*),2);
        Coreg_param[0] = (double*)calloc(sizeof(double),2);
        Coreg_param[1] = (double*)calloc(sizeof(double),2);
        SDM_ortho(_filename, args, Coreg_param);
        free(Coreg_param);
    }
    else if(args.check_sdm_ortho == 2)
    {
        D2DPOINT *adjust_sdt = (D2DPOINT*)calloc(sizeof(D2DPOINT),proinfo->number_of_images);
        double **Coreg_param = ImageCoregistration(return_param, _filename, args, _save_filepath, 1,adjust_sdt,&cal_check);
        printf("cal_check %d\n",cal_check);
        if(cal_check)
        {
            SDM_ortho(_filename, args, Coreg_param);
            free(Coreg_param);
        }
        free(adjust_sdt);
    }
    else
    {
        if(OpenProject(_filename,proinfo,args))
        {
            if(Maketmpfolders(proinfo))
            {
                const uint8 NumOfIAparam  = 2;
                
                ImageInfo *image_info = (ImageInfo*)malloc(sizeof(ImageInfo)*proinfo->number_of_images);
                CSize *Limagesize = (CSize*)malloc(sizeof(CSize)*proinfo->number_of_images); //original imagesize
                double ***RPCs = (double***)calloc(proinfo->number_of_images, sizeof(double**));
                
                sprintf(_save_filepath,"%s",proinfo->save_filepath);
                
                printf("Completion of loading project file!!\n");
                printf("# of detected threads by openmp = %d\n",omp_get_max_threads());
                printf("# of allocated threads = %d\tinput image counts = %d\n",omp_get_max_threads(),proinfo->number_of_images);
                
                char metafilename[500];
                
                FILE *pMetafile = NULL;
                sprintf(metafilename, "%s/%s_meta.txt", proinfo->save_filepath, proinfo->Outputpath_name);
                if(args.check_Matchtag)
                    sprintf(metafilename, "%s/%s_new_matchtag_meta.txt", proinfo->save_filepath, proinfo->Outputpath_name);
                
                if(!proinfo->check_checktiff && !args.check_ortho)
                {
                    pMetafile   = fopen(metafilename,"w");
                    
                    fprintf(pMetafile,"SETSM Version=%s\n", setsm_version);
                }
                
                time_t current_time;
                char*   c_time_string;
                
                current_time = time(NULL);
                c_time_string = ctime(&current_time);
                
                char temp_filepath[500];
                double *Image_gsd_r = (double*)calloc(sizeof(double),proinfo->number_of_images);
                double *Image_gsd_c = (double*)calloc(sizeof(double),proinfo->number_of_images);
                double *Image_gsd = (double*)calloc(sizeof(double),proinfo->number_of_images);
                
                ImageGSD *GSD_image = (ImageGSD*)calloc(sizeof(ImageGSD),proinfo->number_of_images);
                BandInfo *leftright_band = (BandInfo*)calloc(sizeof(BandInfo),proinfo->number_of_images);
                
                ImageGSD GSD_image1;
                GSD_image1.row_GSD = 0;
                GSD_image1.col_GSD = 0;
                GSD_image1.pro_GSD = 0;
                
                printf("sensor_provider %d\t%d\n",args.sensor_type,args.sensor_provider);
                
                double mean_product_res;
                double convergence_angle;
                if(args.sensor_type == AB)
                {
                    for(int ti=0;ti<proinfo->number_of_images ;ti++)
                    {
                        RPCs[ti]       = OpenXMLFile(proinfo,ti,&Image_gsd_r[ti],&Image_gsd_c[ti],&Image_gsd[ti],&leftright_band[ti]);
                        
                        image_info[ti].GSD.row_GSD = Image_gsd_r[ti];
                        image_info[ti].GSD.col_GSD = Image_gsd_c[ti];
                        image_info[ti].GSD.pro_GSD = Image_gsd[ti];
                        
                        GSD_image1.row_GSD += Image_gsd_r[ti];
                        GSD_image1.col_GSD += Image_gsd_c[ti];
                        GSD_image1.pro_GSD += Image_gsd[ti];
                        mean_product_res = GSD_image1.pro_GSD/proinfo->number_of_images;
                        
                        GetImageSize(proinfo->Imagefilename[ti],&Limagesize[ti]);
                        proinfo->frameinfo.m_Camera.m_ImageSize.width = Limagesize[ti].width;
                        proinfo->frameinfo.m_Camera.m_ImageSize.height = Limagesize[ti].height;
                        
                        printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
                               proinfo->frameinfo.Photoinfo[ti].path,
                               proinfo->frameinfo.Photoinfo[ti].m_Xl,proinfo->frameinfo.Photoinfo[ti].m_Yl,proinfo->frameinfo.Photoinfo[ti].m_Zl,
                               proinfo->frameinfo.Photoinfo[ti].m_Wl,proinfo->frameinfo.Photoinfo[ti].m_Pl,proinfo->frameinfo.Photoinfo[ti].m_Kl,
                               proinfo->frameinfo.m_Camera.m_ImageSize.width,proinfo->frameinfo.m_Camera.m_ImageSize.height);
                    }
                    convergence_angle = 40;
                    
                    printf("Load EO\n");
                }
                else
                {
                    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    {
                        if(args.sensor_provider == DG)
                        {
                            RPCs[ti]       = OpenXMLFile(proinfo,ti,&Image_gsd_r[ti],&Image_gsd_c[ti],&Image_gsd[ti],&leftright_band[ti]);
                            
                            image_info[ti].GSD.row_GSD = Image_gsd_r[ti];
                            image_info[ti].GSD.col_GSD = Image_gsd_c[ti];
                            image_info[ti].GSD.pro_GSD = Image_gsd[ti];
                            
                            GSD_image1.row_GSD += Image_gsd_r[ti];
                            GSD_image1.col_GSD += Image_gsd_c[ti];
                            GSD_image1.pro_GSD += Image_gsd[ti];
                            mean_product_res = GSD_image1.pro_GSD/proinfo->number_of_images;
                            
                            GetImageSize(proinfo->Imagefilename[ti],&Limagesize[ti]);
                            
                            OpenXMLFile_orientation(proinfo->RPCfilename[ti],&image_info[ti]);
                            
                            image_info[0].convergence_angle = acos(sin(image_info[0].Mean_sat_elevation*DegToRad)*sin(image_info[1].Mean_sat_elevation*DegToRad) + cos(image_info[0].Mean_sat_elevation*DegToRad)*cos(image_info[1].Mean_sat_elevation*DegToRad)*cos( (image_info[0].Mean_sat_azimuth_angle - image_info[1].Mean_sat_azimuth_angle)*DegToRad))*RadToDeg;
                            
                            convergence_angle = image_info[0].convergence_angle;
                            
                            printf("%d_image info\nSatID = %s\nAcquisition_time = %s\nMean_row_GSD = %f\nMean_col_GSD = %f\nMean_GSD = %f\nMean_sun_azimuth_angle = %f\nMean_sun_elevation = %f\nMean_sat_azimuth_angle = %f\nMean_sat_elevation = %f\nIntrack_angle = %f\nCrosstrack_angle = %f\nOffnadir_angle = %f\ntdi = %d\neffbw = %f\nabscalfact = %f\nconvergence_angle = %f\n",ti+1,image_info[ti].SatID,image_info[ti].imagetime,Image_gsd_r[ti],Image_gsd_c[ti],Image_gsd[ti],image_info[ti].Mean_sun_azimuth_angle,image_info[ti].Mean_sun_elevation,image_info[ti].Mean_sat_azimuth_angle,image_info[ti].Mean_sat_elevation,image_info[ti].Intrack_angle,image_info[ti].Crosstrack_angle,image_info[ti].Offnadir_angle,(int)leftright_band[ti].tdi,leftright_band[ti].effbw,leftright_band[ti].abscalfactor,image_info[0].convergence_angle);
                        }
                        else if(args.sensor_provider == PL)
                        {
                            RPCs[ti]       = OpenXMLFile_Pleiades(proinfo->RPCfilename[ti]);
                            GetImageSize(proinfo->Imagefilename[ti],&Limagesize[ti]);
                            convergence_angle = 40;
                            mean_product_res = 0.5;
                        }
                        else if(args.sensor_provider == PT)
                        {
                            printf("planet header loading\n");
                            RPCs[ti]       = OpenXMLFile_Planet(proinfo->RPCfilename[ti]);
                            GetImageSize(proinfo->Imagefilename[ti],&Limagesize[ti]);
                            convergence_angle = 40;
                        }
                    }
                }
                
                if(!args.check_imageresolution)
                {
                    if(proinfo->sensor_type == AB)
                    {
                        proinfo->resolution = (int)((mean_product_res)*10 + 0.5)/10.0;
                    }
                    else if(args.sensor_provider == DG)
                    {
                        proinfo->resolution = (int)((mean_product_res)*10 + 0.5)/10.0;
                        
                        if (proinfo->resolution < 0.75)
                        {
                            proinfo->resolution = 0.5;
                        }
                        else if(proinfo->resolution < 2.0)
                            proinfo->resolution = 1.0;
                        else
                            proinfo->resolution = floor(proinfo->resolution);
                    }
                    else
                        proinfo->resolution = 0.5;
                }
                else
                {
                    proinfo->resolution  = args.image_resolution;
                    mean_product_res = args.image_resolution;
                }
                
                double MPP_stereo_angle = 1;
                CalMPP_pair(convergence_angle,mean_product_res, proinfo->resolution, &MPP_stereo_angle);
                
                printf("image resolution %f\t%f\n",proinfo->resolution,mean_product_res);
                
                if(args.check_Matchtag)
                    proinfo->check_Matchtag = args.check_Matchtag;
                
                if(!proinfo->check_checktiff && !args.check_ortho)
                {
                    fprintf(pMetafile,"Creation Date=%s",c_time_string);
                    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                        fprintf(pMetafile,"Image %d=%s\n",ti+1,proinfo->Imagefilename[ti]);
                    fprintf(pMetafile,"Output Resolution=%f\n",proinfo->DEM_resolution);
                }
                
                double Image_res[2] = {proinfo->resolution, proinfo->resolution};
                double Res[2] = {proinfo->resolution, proinfo->DEM_resolution};
                
                TransParam param;
                param.projection = args.projection;
                param.utm_zone   = args.utm_zone;
                
                double minLat, minLon;
                if(proinfo->sensor_type == SB)
                {
                    minLat      = RPCs[0][0][3];
                    minLon      = RPCs[0][0][2];
                    SetTransParam(minLat,minLon,&param);
                    
                    printf("minLat Lon %f\t%f\n",minLat, minLon);
                }
                
                printf("param projection %d\tzone %d\n",param.projection,param.utm_zone);
                *return_param = param;
                
                double Boundary[4], LBoundary[4],RBoundary[4],LminmaxHeight[2],RminmaxHeight[2],ori_minmaxHeight[2];
                double LHinterval, Hinterval;
                double lonlatboundary[4] = {0.0};
                if(args.check_boundary)
                {
                    Boundary[0] = args.Min_X;
                    Boundary[1] = args.Min_Y;
                    Boundary[2] = args.Max_X;
                    Boundary[3] = args.Max_Y;
                    
                    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    {
                        Imageparams[ti][0]  = proinfo->RA_param[ti][0];
                        Imageparams[ti][1]  = proinfo->RA_param[ti][1];
                        
                        if(proinfo->sensor_type == SB)
                            SetDEMBoundary(proinfo, RPCs[ti],Image_res,param,LBoundary,LminmaxHeight,&LHinterval);
                        else
                            SetDEMBoundary_photo(proinfo->frameinfo.Photoinfo[ti], proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm, LBoundary,LminmaxHeight,&LHinterval);
                        
                        if(ti == 0)
                        {
                            Hinterval   = LHinterval;
                            
                            ori_minmaxHeight[0] = LminmaxHeight[0];
                            ori_minmaxHeight[1] = LminmaxHeight[1];
                        }
                        else
                        {
                            if(LHinterval > Hinterval)
                                Hinterval   = LHinterval;
                            
                            ori_minmaxHeight[0] = min(LminmaxHeight[0],ori_minmaxHeight[0]);
                            ori_minmaxHeight[1] = max(LminmaxHeight[1],ori_minmaxHeight[1]);
                        }
                    }
                }
                else
                {
                    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    {
                        Imageparams[ti][0]  = proinfo->RA_param[ti][0];
                        Imageparams[ti][1]  = proinfo->RA_param[ti][1];
                        
                        if(proinfo->sensor_type == SB)
                            SetDEMBoundary(proinfo,RPCs[ti],Image_res,param,LBoundary,LminmaxHeight,&LHinterval);
                        else
                            SetDEMBoundary_photo(proinfo->frameinfo.Photoinfo[ti], proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm, LBoundary,LminmaxHeight,&LHinterval);
                        
                        if(ti == 0)
                        {
                            for(int i=0;i<4;i++)
                                lonlatboundary[i] = LBoundary[i];
                            
                            Hinterval   = LHinterval;
                            
                            ori_minmaxHeight[0] = LminmaxHeight[0];
                            ori_minmaxHeight[1] = LminmaxHeight[1];
                        }
                        else
                        {
                            for(int i=0;i<4;i++)
                            {
                                if(i<2)
                                    lonlatboundary[i] = max(LBoundary[i], lonlatboundary[i]);
                                else
                                    lonlatboundary[i] = min(LBoundary[i], lonlatboundary[i]);
                            }
                            
                            if(LHinterval > Hinterval)
                                Hinterval   = LHinterval;
                            
                            ori_minmaxHeight[0] = min(LminmaxHeight[0],ori_minmaxHeight[0]);
                            ori_minmaxHeight[1] = max(LminmaxHeight[1],ori_minmaxHeight[1]);
                        }
                    }
                    
                    findOverlappArea(proinfo,param,RPCs,Image_res,lonlatboundary);
                    
                    if(proinfo->sensor_type == SB)
                    {
                        printf("lonlatboundary = %f\t%f\t%f\t%f\n",lonlatboundary[0],lonlatboundary[1],lonlatboundary[2],lonlatboundary[3]);
                        
                        D2DPOINT *lonlat = (D2DPOINT *) malloc(sizeof(D2DPOINT) * 4);
                        lonlat[0].m_X = lonlatboundary[0];
                        lonlat[0].m_Y = lonlatboundary[1];
                        lonlat[1].m_X = lonlatboundary[0];
                        lonlat[1].m_Y = lonlatboundary[3];
                        lonlat[2].m_X = lonlatboundary[2];
                        lonlat[2].m_Y = lonlatboundary[3];
                        lonlat[3].m_X = lonlatboundary[2];
                        lonlat[3].m_Y = lonlatboundary[1];
                        
                        D2DPOINT *XY = wgs2ps(param, 4, lonlat);
                        
                        printf("XY X %f\t%f\t%f\t%f\n",XY[0].m_X,XY[1].m_X,XY[2].m_X,XY[3].m_X);
                        printf("XY Y %f\t%f\t%f\t%f\n",XY[0].m_Y,XY[1].m_Y,XY[2].m_Y,XY[3].m_Y);
                        
                        if( lonlatboundary[1] < 0 && lonlatboundary[3] > 0)
                        {
                            double below_eq = 10000000 - XY[0].m_Y;
                            double above_eq = XY[1].m_Y;
                            if(below_eq > above_eq)
                            {
                                XY[1].m_Y = 10000000;
                                XY[2].m_Y = 10000000;
                            }
                            else
                            {
                                XY[0].m_Y = 0;
                                XY[3].m_Y = 0;
                            }
                        }
                        
                        printf("XY X %f\t%f\t%f\t%f\n",XY[0].m_X,XY[1].m_X,XY[2].m_X,XY[3].m_X);
                        printf("XY Y %f\t%f\t%f\t%f\n",XY[0].m_Y,XY[1].m_Y,XY[2].m_Y,XY[3].m_Y);
                        
                        double minX = min(XY[3].m_X, min(XY[2].m_X, min(XY[0].m_X, XY[1].m_X)));
                        double maxX = max(XY[3].m_X, max(XY[2].m_X, max(XY[0].m_X, XY[1].m_X)));
                        
                        double minY = min(XY[3].m_Y, min(XY[2].m_Y, min(XY[0].m_Y, XY[1].m_Y)));
                        double maxY = max(XY[3].m_Y, max(XY[2].m_Y, max(XY[0].m_Y, XY[1].m_Y)));
                        
                        free(lonlat);
                        free(XY);
                        
                        Boundary[0] = ceil(minX/2.0)*2;
                        Boundary[1] = ceil(minY/2.0)*2;
                        Boundary[2] = floor(maxX/2.0)*2;
                        Boundary[3] = floor(maxY/2.0)*2;
                    }
                    else
                    {
                        Boundary[0] = lonlatboundary[0];
                        Boundary[1] = lonlatboundary[1];
                        Boundary[2] = lonlatboundary[2];
                        Boundary[3] = lonlatboundary[3];
                    }
                    
                    
                }
                printf("boundary = %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                
                CSize Boundary_size(ceil(Boundary[2] - Boundary[0]), ceil(Boundary[3] - Boundary[1]));
                
                if(Boundary_size.height/1000.0 > args.overlap_length || Boundary_size.width/1000.0 > args.overlap_length)
                {
                    double new_height = (10000000 - Boundary[3] + Boundary[1])/1000.0;
                    if( new_height > args.overlap_length)
                    {
                        printf("Overlapped area between stereo pair is very long along the strip(height=%3.2f(km), width=%3.2f(km)), so that the assumption of RPC bias computation (less than 50km) is not satisfied,\nso relative RPC bias can be not accurately compensated. \nPlease process after split the overlapped area in strip direction into several small area less than 30 km\nBounary(minX, minY, maxX, maxY[m]) = %f %f %f %f\n",Boundary_size.height/1000.0,Boundary_size.width/1000.0,Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                        exit(1);
                    }
                    else
                    {
                        if(10000000 - Boundary[3] > Boundary[1])
                        {
                            Boundary[1] = Boundary[3];
                            Boundary[3] = 10000000;
                        }
                        else
                        {
                            Boundary[3] = Boundary[1];
                            Boundary[1] = 0;
                        }
                        
                        Boundary_size.width     = Boundary[2] - Boundary[0];
                        Boundary_size.height    = Boundary[3] - Boundary[1];
                        printf("cross equator boundary_size %f\t%d\t%d\n",new_height,Boundary_size.width,Boundary_size.height);
                        printf("boundary = %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                    }
                }
                
                CSize Matchtag_seeddem_size;
                if(proinfo->check_Matchtag && !args.check_tiles_SR && !args.check_tiles_SC && !args.check_tiles_ER && !args.check_tiles_EC)
                {
                    double t_minX,t_maxY, t_grid_size;
                    char *seed_path = proinfo->priori_DEM_tif;
                    
                    Matchtag_seeddem_size = ReadGeotiff_info(seed_path, &t_minX, &t_maxY, &t_grid_size);
                    
                    Boundary_size = Matchtag_seeddem_size;
                    
                    Boundary[0] = t_minX;
                    Boundary[1] = t_maxY - t_grid_size*Boundary_size.height;
                    Boundary[2] = t_minX + t_grid_size*Boundary_size.width;
                    Boundary[3] = t_maxY;
                    
                    printf("size %d\t%d\t boundary %f\t%f\t%f\t%f\n",Boundary_size.width,Boundary_size.height, Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                }
                
                if (args.check_minH) {
                    ori_minmaxHeight[0] = (int)args.minHeight;
                    printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
                }
                if (args.check_maxH) {
                    ori_minmaxHeight[1] = (int)args.maxHeight;
                    printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
                }
                
                if (proinfo->check_minH) {
                    ori_minmaxHeight[0] = (int)proinfo->minHeight;
                    printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
                }
                if (proinfo->check_maxH) {
                    ori_minmaxHeight[1] = (int)proinfo->maxHeight;
                    printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
                }
                
                if(proinfo->pre_DEMtif)
                {
                    ori_minmaxHeight[1] += 1000;
                    ori_minmaxHeight[0] -= 1000;
                    if(ori_minmaxHeight[0] < -100)
                        ori_minmaxHeight[0] = -100;
                }
                
                if (ori_minmaxHeight[1] >9000) {
                    ori_minmaxHeight[1] = 9000;
                }
                if (ori_minmaxHeight[0] < -100)
                    ori_minmaxHeight[0] = -100;
                
                if(!args.check_ortho)
                {
                    printf("minmaxH = %f\t%f\n",ori_minmaxHeight[0],ori_minmaxHeight[1]);
                    printf("seed fff %d\n",proinfo->pre_DEMtif);
                    
                    PairInfo pairinfo;
                    pairinfo.pairs = (UI2DPOINT*)malloc(sizeof(UI2DPOINT)*MaxNCC);
                    pairinfo.BHratio = (float*)malloc(sizeof(float)*MaxNCC);
                    
                    if(SetupParam(proinfo,&proinfo->pre_DEMtif))
                    {
                        const double bin_angle = 360.0/18.0;
                        const uint8 Template_size = 15;
                        uint16 buffer_area;
                        uint8 pyramid_step = proinfo->pyramid_level;
                        
                        double seedDEM_gridsize;
                        
                        uint8 iter_row_start, iter_row_end, t_col_start, t_col_end;
                        int tile_size = 0;
                        double subX, subY;
                        
                        time_t ST = 0, ET = 0;
                        double gap;
                        
                        printf("IsRA = %d\n",proinfo->IsRA);
                        if(proinfo->IsRA)
                        {
                            uint8 RA_row_iter = 1;
                            uint8 RA_col_iter = 1;
                            
                            tile_size           = 8000;
                            if(args.check_tilesize)
                                tile_size       = args.tilesize*2;
                            printf("tileszie %d\n",tile_size);
                            
                            bool check_RA_1000  = false;
                            
                            if (args.check_RA_tileR)
                            {
                                if(args.RA_row == 100)
                                {
                                    tile_size           = 50000;
                                    check_RA_1000       = true;
                                    args.RA_row         = 1;
                                }
                            }
                            
                            SetTiles_RA(proinfo, Boundary, tile_size, &pyramid_step, &iter_row_start, &iter_row_end, &RA_row_iter, &t_col_start, &t_col_end, &RA_col_iter, &subX, &subY);
                            
                            printf("before RA row:col => row = %d\t%d\t;col = %d\t%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end);
                            
                            if (args.check_RA_tileR)
                            {
                                iter_row_start  = args.RA_row;
                                iter_row_end    = iter_row_start+1;
                            }
                            
                            if (args.check_RA_tileC)
                            {
                                t_col_start     = args.RA_col;
                                t_col_end       = t_col_start + 1;
                            }
                            
                            printf("RA row:col = row = %d\t%d\t;col = %d\t%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end);
                            
                            //loading existing RA parameters
                            bool check_load_RA = false;
                            char bufstr[500];
                            
                            //check RA info in seeddem folder
                            printf("check seedem %d\n",args.check_seeddem);
                            if(args.check_seeddem)
                            {
                                //meta file check in seeddem folder
                                printf("Meta file %s\n",proinfo->metafilename);
                                FILE* pFile_meta  = fopen(proinfo->metafilename,"r");
                                if(pFile_meta)
                                {
                                    printf("meta file exist!!\n");
                                    int ti = 0;
                                    while(!feof(pFile_meta))
                                    {
                                        fgets(bufstr,500,pFile_meta);
                                        
                                        bool check_load_RA_each = false;
                                        if (strstr(bufstr,"RA Params=")!=NULL)
                                        {
                                            printf("%s\n",bufstr);
                                            int temp_ref,temp_tar;
                                            double bh_ratio;
                                            sscanf(bufstr,"RA Params=%d\t%d\t%lf\t%lf\t%lf\n",&temp_ref,&temp_tar,&Imageparams[ti][0],&Imageparams[ti][1],&bh_ratio);
                                            pairinfo.pairs[ti].m_X = temp_ref;
                                            pairinfo.pairs[ti].m_Y = temp_tar;
                                            pairinfo.BHratio[ti] = bh_ratio;
                                            /*if(ti == 0)
                                                proinfo->check_selected_image[ti] = true;
                                            else*/
                                            {
                                                if(Imageparams[ti][0] == 0 || Imageparams[ti][1] == 0)
                                                    check_load_RA_each = true;
                                                //else
                                                    
                                                    //proinfo->check_selected_image[ti] = false;
                                            }
                                            
                                            printf("meta RA %d\t%d\t%d\t%f\t%f\n",ti,temp_ref,temp_tar,Imageparams[ti][0],Imageparams[ti][1]);
                                            /*
                                            if(proinfo->number_of_images == 2 && ti == 0 && Imageparams[ti][0] != 0 && Imageparams[ti][1] != 0)
                                            {
                                                check_load_RA = true;
                                                Imageparams[1][0] = Imageparams[ti][0];
                                                Imageparams[1][1] = Imageparams[ti][1];

                                                Imageparams[0][0] = 0;
                                                Imageparams[0][1] = 0;
                                                proinfo->check_selected_image[1] = true;
                                            }
                                            */
                                            ti++;
                                            
                                            pairinfo.NumberOfPairs = ti;
                                        }
                                        else if(strstr(bufstr,"SETSM Version=")!=NULL)
                                        {
                                            printf("%s\n",bufstr);
                                            double version;
                                            sscanf(bufstr,"SETSM Version=%lf\n",&version);
                                            printf("version %f\n",version);
                                            
                                            if(!args.check_Matchtag)
                                            {
                                                if (version > 2.0128) {
                                                    if(args.seedDEMsigma < 10)
                                                        proinfo->seedDEMsigma = 10;
                                                    else
                                                        proinfo->seedDEMsigma = args.seedDEMsigma;
                                                }
                                                else {
                                                    proinfo->seedDEMsigma = 100;
                                                }
                                            }
                                            printf("sigma %f\n",proinfo->seedDEMsigma);
                                        }
                                        else if(strstr(bufstr,"Output Resolution=")!=NULL)
                                        {
                                            printf("%s\n",bufstr);
                                            double version;
                                            sscanf(bufstr,"Output Resolution=%lf\n",&seedDEM_gridsize);
                                            printf("seed DEM gridsize %f\n",seedDEM_gridsize);
                                        }
                                        
                                        if(!check_load_RA_each)
                                            check_load_RA = true;
                                    }
                                    fclose(pFile_meta);
                                }
                                else
                                    printf("meta file doesn't exist!! \n");
                                
                                //RAinfo file check in seeddem folder
                                if(!check_load_RA)
                                {
                                    int full_size       = strlen(args.seedDEMfilename);
                                    printf("full_size %d\n",full_size);
                                    char RAfile_raw[500];
                                    for (int i = 0; i < full_size - 14; i++)
                                        RAfile_raw[i] = args.seedDEMfilename[i];
                                    
                                    char str_rafile_2[500];
                                    sprintf(str_rafile_2, "%s/txt/RAinfo.txt", RAfile_raw);
                                    printf("RA file %s\n", str_rafile_2);
                                    
                                    check_load_RA = GetRAinfo(proinfo, str_rafile_2, Imageparams,pairinfo);
                                }
                                
                                //echo_result file check in seeddem folder
                                if(!check_load_RA)
                                {
                                    char *lastSlash   = (char*)malloc(strlen(args.seedDEMfilename) + 1);
                                    strcpy(lastSlash, args.seedDEMfilename);
                                    char *fullseeddir = dirname(lastSlash);
                                    int dir_size        = strlen(fullseeddir);
                                    
                                    printf("fullseeddir %s\tdir_size %d\n",fullseeddir,dir_size);
                                    printf("lastSlash %s\n",lastSlash);
                                    
                                    char seeddir[500];
                                    for(int i=0;i<dir_size-4;i++)
                                        seeddir[i] = fullseeddir[i];
                                    
                                    char str_echofile[500];
                                    sprintf(str_echofile,"%s/txt/echo_result_row_1_col_1.txt",fullseeddir);
                                    
                                    check_load_RA = GetRAinfoFromEcho(proinfo, str_echofile, Imageparams,pairinfo);
                                }
                            }
                            else
                            {
                                //RA_echo and echo_result file check in current working folder
                                char str_echofile[500];
                                sprintf(str_echofile,"%s/txt/echo_result_row_1_col_1.txt",proinfo->save_filepath);
                                
                                printf("echo %s\n",str_echofile);
                                
                                check_load_RA = GetRAinfoFromEcho(proinfo, str_echofile, Imageparams,pairinfo);
                                
                                //RAinfo file check in current working folder
                                if(!check_load_RA)
                                {
                                    char str_rafile[500];
                                    FILE* pFile_info;
                                    sprintf(str_rafile,"%s/txt/RAinfo.txt",proinfo->save_filepath);
                                    printf("RAinfo %s\n",str_rafile);
                                    
                                    check_load_RA = GetRAinfo(proinfo, str_rafile, Imageparams,pairinfo);
                                }
                            }
                            if(check_load_RA)
                            {
                                for(int ti = 0; ti < pairinfo.NumberOfPairs ; ti++)
                                    printf("check load RA %d\t%f\t%f\t%d\t%d\t%f\n",check_load_RA,Imageparams[ti][0],Imageparams[ti][1], pairinfo.pairs[ti].m_X,pairinfo.pairs[ti].m_Y,pairinfo.BHratio[ti]);
                            }
                            
                            if(!check_load_RA)
                            {
                                tile_size           = 40000;
                                
                                args.RA_row         = 1;
                                
                                SetTiles_RA(proinfo, Boundary,  tile_size, &pyramid_step, &iter_row_start, &iter_row_end, &RA_row_iter, &t_col_start, &t_col_end, &RA_col_iter, &subX, &subY);
                                
                                double temp_DEM_resolution = proinfo->DEM_resolution;
                                proinfo->DEM_resolution = Image_res[0]*pwrtwo(pyramid_step+1);
                                buffer_area  = 400;
                                
                                Matching_SETSM(proinfo,image_info,pyramid_step, Template_size, buffer_area,1,2,1,2,subX,subY,bin_angle,Hinterval,Image_res, Imageparams, RPCs, NumOfIAparam, Limagesize,param, ori_minmaxHeight,Boundary,convergence_angle,mean_product_res,&MPP_stereo_angle,pairinfo);
                                proinfo->DEM_resolution = temp_DEM_resolution;
                            }
                        }
                        
                        if (args.check_RA_line && args.check_RA_sample)
                        {
                            if(args.RA_line_count != args.RA_sample_count)
                            {
                                printf("please check RPC bias input in line and sample direction!");
                                exit(1);
                            }
                            else
                            {
                                pairinfo.NumberOfPairs = args.RA_line_count;
                                printf("RA_line_count %d\n",args.RA_line_count);
                                for(int ti = 0; ti < pairinfo.NumberOfPairs ; ti++)
                                {
                                    Imageparams[ti][0] = args.ra_line[ti];
                                    Imageparams[ti][1] = args.ra_sample[ti];
                                    
                                    //need to handle manually input RA pairs
                                    pairinfo.pairs[ti].m_X = 0;
                                    pairinfo.pairs[ti].m_Y = 1;
                                    
                                    printf("loaded imageprarm %f\t%f\n",Imageparams[ti][0],Imageparams[ti][1]);
                                }
                            }
                        }
                        
                        if(!proinfo->check_checktiff)
                        {
                            for(int ti = 0; ti < pairinfo.NumberOfPairs ; ti++)
                                fprintf(pMetafile,"RA Params=%d\t%d\t%f\t%f\t%f\n",pairinfo.pairs[ti].m_X,pairinfo.pairs[ti].m_Y,Imageparams[ti][0],Imageparams[ti][1],pairinfo.BHratio[ti]);
                   
                            fprintf(pMetafile,"RA tilesize=%d\n",tile_size);
                        }
                        
                        proinfo->IsRA        = false;
                        
                        if(!args.RA_only)
                        {
                            tile_size           = 4000;
                            
                            if(Boundary_size.width < tile_size && Boundary_size.height < tile_size)
                            {
                                if(Boundary_size.width > Boundary_size.height)
                                    tile_size = Boundary_size.width;
                                else
                                    tile_size = Boundary_size.height;
                            }
                            
                            if(args.check_tilesize)
                                tile_size       = args.tilesize;
                            printf("tilesize %d\n",tile_size);
                            
                            if(!proinfo->check_checktiff)
                            {
                                fprintf(pMetafile,"tilesize=%d\n",tile_size);
                                if(proinfo->pre_DEMtif)
                                    fprintf(pMetafile,"Seed DEM=%s\n",proinfo->priori_DEM_tif);
                                else
                                    fprintf(pMetafile,"Seed DEM=\n");
                                
                                if(param.projection == 1)
                                {
                                    if (param.bHemisphere)
                                        fprintf(pMetafile, "Output Projection='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45, +k=1 +x_0=0 +y_0=0 +datum=WGS84 +unit=m +no_defs'\n");
                                    else
                                        fprintf(pMetafile, "Output Projection='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0, +k=1 +x_0=0 +y_0=0 +datum=WGS84 +unit=m +no_defs'\n");
                                }
                                else
                                {
                                    if (param.bHemisphere)
                                        fprintf(pMetafile, "Output Projection='+proj=utm +zone=%d +north=%s +datum=WGS84 +unit=m +no_defs'\n", param.utm_zone, param.direction);
                                    else
                                        fprintf(pMetafile, "Output Projection='+proj=utm +zone=%d +south=%s +datum=WGS84 +unit=m +no_defs'\n", param.utm_zone, param.direction);
                                }
                            }
                            
                            buffer_area  = tile_size*0.1;
                            SetTiles(proinfo, Boundary, tile_size, &pyramid_step, &buffer_area, &iter_row_start, &iter_row_end, &t_col_start, &t_col_end, &subX, &subY);
                            
                            if (args.check_tiles_SR)
                                iter_row_start    = args.start_row;
                            
                            if (args.check_tiles_ER)
                                iter_row_end      = args.end_row;
                            
                            if (args.check_tiles_SC)
                                t_col_start       = args.start_col;
                            
                            if (args.check_tiles_EC)
                                t_col_end         = args.end_col;
                            
                            if (proinfo->check_tiles_SR)
                                iter_row_start    = proinfo->start_row;
                            
                            if (proinfo->check_tiles_ER)
                                iter_row_end      = proinfo->end_row;
                            
                            if (proinfo->check_tiles_SC)
                                t_col_start       = proinfo->start_col;
                            
                            if (proinfo->check_tiles_EC)
                                t_col_end         = proinfo->end_col;
                            
                            for(int ti = 0; ti < pairinfo.NumberOfPairs ; ti++)
                                printf("RA param = %f\t%f\n",Imageparams[ti][0],Imageparams[ti][1]);
                            
                            printf("Tiles row:col = row = %d\t%d\t;col = %d\t%d\tseed flag =%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end,proinfo->pre_DEMtif);
                            
                            if(!args.check_gridonly)
                            {
                                Matching_SETSM(proinfo,image_info,pyramid_step, Template_size, buffer_area,iter_row_start, iter_row_end,t_col_start,t_col_end,subX,subY,bin_angle,Hinterval,Image_res,Imageparams,RPCs, NumOfIAparam, Limagesize,param,ori_minmaxHeight,Boundary,convergence_angle,mean_product_res,&MPP_stereo_angle,pairinfo);
                            }
#ifdef BUILDMPI
                            MPI_Barrier(MPI_COMM_WORLD);
                            MPI_Finalize();
                            if(rank != 0)
                                exit(0);
#endif
                            if(!args.check_ortho)
                            {
                                char check_file[500];
                                int max_row = 0;
                                int max_col = 0;
                                for(int row = 1; row < 100 ; row++)
                                {
                                    for(int col = 1; col < 100 ; col++)
                                    {
                                        sprintf(check_file,"%s/txt/matched_pts_%d_%d_0_3.txt",proinfo->save_filepath,row,col);
                                        FILE *pcheckFile = fopen(check_file,"r");
                                        if(pcheckFile)
                                        {
                                            if(max_row < row)
                                                max_row = row;
                                            if(max_col < col)
                                                max_col = col;
                                            
                                            fclose(pcheckFile);
                                        }
                                    }
                                }
                                
                                if (args.check_tiles_SR)
                                    iter_row_start    = args.start_row;
                                
                                if (args.check_tiles_ER)
                                    iter_row_end      = args.end_row;
                                else
                                {
                                    if (proinfo->check_tiles_ER)
                                        iter_row_end      = proinfo->end_row;
                                    else
                                        iter_row_end = max_row + 1;
                                }
                                
                                if (args.check_tiles_SC)
                                    t_col_start       = args.start_col;
                                
                                if (args.check_tiles_EC)
                                    t_col_end         = args.end_col;
                                else
                                {
                                    if (proinfo->check_tiles_EC)
                                        t_col_end         = proinfo->end_col;
                                    else
                                        t_col_end    = max_col + 1;
                                }
                                
                                if (proinfo->check_tiles_SR)
                                    iter_row_start    = proinfo->start_row;
                                
                                if (proinfo->check_tiles_SC)
                                    t_col_start       = proinfo->start_col;
                            }
                            
                            printf("Tiles row:col = row = %d\t%d\t;col = %d\t%d\tseed flag =%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end,proinfo->pre_DEMtif);
                            
                            if(iter_row_end < 2 && t_col_end < 2)
                            {
                                printf("No matching results. Please check overlapped area of stereo pair, or image textures\n");
                                exit(1);
                            }
                            
                            char str_DEMfile[500];
                            float *DEM_values = NULL;
                            signed char *Ortho_values = NULL;
                            float *H_value = NULL;
                            unsigned char* MT_value = NULL;
                            
                            int buffer_tile = buffer_area + 20;
                            const int final_iteration = 3;
                            
                            double FinalDEM_boundary[4] = {0.};
                            CSize Final_DEMsize;
                            
                            int total_row_num = iter_row_end - iter_row_start;
                            int total_col_num = t_col_end - t_col_start;
                            double DEM_width = Boundary[2] - Boundary[0];
                            double DEM_height = Boundary[3] - Boundary[1];
                            
                            if((total_col_num == 1 && total_row_num == 1) || (DEM_width < buffer_tile || DEM_height < buffer_tile))
                                buffer_tile = 0;
                            
                            if(proinfo->check_Matchtag && !args.check_tiles_SR && !args.check_tiles_SC && !args.check_tiles_ER && !args.check_tiles_EC)
                            {
                                buffer_tile = 0;
                                
                                FinalDEM_boundary[0] = Boundary[0];
                                FinalDEM_boundary[1] = Boundary[1];
                                FinalDEM_boundary[2] = Boundary[2];
                                FinalDEM_boundary[3] = Boundary[3];
                                
                                Final_DEMsize = Matchtag_seeddem_size;
                            }
                            else
                                Final_DEMsize = DEM_final_Size(proinfo->save_filepath, iter_row_start,t_col_start, iter_row_end,t_col_end,proinfo->DEM_resolution,FinalDEM_boundary);
                            
                            double total_memory = CalMemorySize_Post(Final_DEMsize,Final_DEMsize);
                            
                            if(args.check_LSF2 == 1 || args.check_LSF2 == 2)
                                total_memory = CalMemorySize_Post_LSF(Final_DEMsize,Final_DEMsize);
                            
                            printf("total tile memory %f\t%f\t%d\t%d\n",proinfo->System_memory,total_memory,Final_DEMsize.width,Final_DEMsize.height);
                            
                            if(total_memory > proinfo->System_memory - 5)
                            {
                                double define_row = 2.0;
                                int tile_row_step;
                                const int tile_row_half = ceil(iter_row_end/2.0);
                                double row_tilememory;
                                
                                char outsefile[500];
                                sprintf(outsefile, "%s/DEM_splited.txt", proinfo->save_filepath);
                                FILE *p_sefile = fopen(outsefile,"w");
                                
                                bool check_tile = true;
                                while(check_tile && define_row < 10.0)
                                {
                                    tile_row_step = ceil((iter_row_end-1)/define_row);
                                    
                                    int tile_row_step_half = floor(tile_row_step/2.0);
                                    int tile_row_start = tile_row_half - tile_row_step_half;
                                    int tile_row_end = tile_row_start + tile_row_step - 1;
                                    
                                    printf("%d\t%d\t%d\t%d\t%d\n",tile_row_half,tile_row_step,tile_row_step_half,tile_row_start,tile_row_end);
                                    
                                    CSize tile_Final_DEMsize = DEM_final_Size(proinfo->save_filepath, tile_row_start,t_col_start, tile_row_end,t_col_end,proinfo->DEM_resolution,FinalDEM_boundary);
                                    
                                    row_tilememory = CalMemorySize_Post(tile_Final_DEMsize,tile_Final_DEMsize);
                                    if(args.check_LSF2 == 1 || args.check_LSF2 == 2)
                                        row_tilememory = CalMemorySize_Post_LSF(tile_Final_DEMsize,tile_Final_DEMsize);
                                    
                                    printf("tile_Final_DEMsize %d\t%d\t%d\t%d\t%d\t%d\n",tile_Final_DEMsize.width,tile_Final_DEMsize.height,tile_row_start,t_col_start, tile_row_end,t_col_end);
                                    printf("while tile_row_step %d\tdivide %f\tmemory %f\n",tile_row_step,define_row,row_tilememory);
                                    
                                    if(row_tilememory < proinfo->System_memory - 5)
                                        check_tile = false;
                                    else
                                        define_row++;
                                }
                                
                                printf("define_row tile_row_step %f\t%d\t%f\n",define_row,tile_row_step,row_tilememory);
                                DEM_divide = define_row;
                                fprintf(p_sefile,"total_splited_DEMs %d\n",(int)define_row);
                                
                                for(int tile_row = 1 ; tile_row <= define_row ; tile_row++)
                                {
                                    iter_row_start  = 1 + tile_row_step*(tile_row-1);
                                    iter_row_end    = tile_row_step*tile_row;
                                    
                                    printf("iter %d\tTiles row:col = row = %d\t%d\t;col = %d\t%d\n",tile_row,iter_row_start,iter_row_end,t_col_start,t_col_end);
                                    
                                    sprintf(str_DEMfile, "%s/%s_%d_dem.tif", proinfo->save_filepath,proinfo->Outputpath_name,tile_row);
                                    
                                    FILE *pFile_DEM = fopen(str_DEMfile,"r");
                                    
                                    ST = time(0);
                                    printf("Tile merging start final iteration %d!!\n",final_iteration);
                                    
                                    CSize tile_Final_DEMsize =
                                    DEM_final_Size(proinfo->save_filepath, iter_row_start,t_col_start, iter_row_end,t_col_end,proinfo->DEM_resolution,FinalDEM_boundary);
                                    
                                    printf("DEM_size %d\t%d\t%d\t%d\t%d\t%d\n",tile_Final_DEMsize.width,tile_Final_DEMsize.height,iter_row_start,iter_row_end,t_col_start,t_col_end);
                                    
                                    if(tile_row == 1)
                                    {
                                        FinalDEM_boundary[3] += 200;
                                        iter_row_end += 1;
                                    }
                                    else if(tile_row == define_row)
                                        FinalDEM_boundary[1] -= 200;
                                    else
                                    {
                                        FinalDEM_boundary[1] -= 200;
                                        FinalDEM_boundary[3] += 200;
                                        iter_row_end += 1;
                                        iter_row_start -= 1;
                                    }
                                    printf("re boundary %f\t%f\t%f\t%f\n",FinalDEM_boundary[0],FinalDEM_boundary[1],FinalDEM_boundary[2],FinalDEM_boundary[3]);
                                    tile_Final_DEMsize.width = (int)((FinalDEM_boundary[2] - FinalDEM_boundary[0])/proinfo->DEM_resolution) + 1;
                                    tile_Final_DEMsize.height = (int)((FinalDEM_boundary[3] - FinalDEM_boundary[1])/proinfo->DEM_resolution) + 1;
                                    printf("re DEM_size %d\t%d\t%d\t%d\t%d\t%d\n",tile_Final_DEMsize.width,tile_Final_DEMsize.height,iter_row_start,iter_row_end,t_col_start,t_col_end);
                                    
                                    fprintf(p_sefile,"DEM %d\n",tile_row);
                                    fprintf(p_sefile,"Output dimensions=%d\t%d\n",tile_Final_DEMsize.width,tile_Final_DEMsize.height);
                                    fprintf(p_sefile,"Upper left coordinates=%f\t%f\n",FinalDEM_boundary[0],FinalDEM_boundary[3]);
                                    
                                    DEM_values = (float*)malloc(sizeof(float)*tile_Final_DEMsize.width*tile_Final_DEMsize.height);
                                    MergeTiles(proinfo,iter_row_start,t_col_start,iter_row_end,t_col_end,buffer_tile,final_iteration,DEM_values,tile_Final_DEMsize,FinalDEM_boundary);
                                    
                                    printf("Interpolation start!!\n");
                                    printf("%f %f\n",proinfo->DEM_resolution,proinfo->DEM_resolution);
                                    
                                    H_value = (float*)malloc(sizeof(float)*(long)tile_Final_DEMsize.width*(long)tile_Final_DEMsize.height);
                                    MT_value = (unsigned char*)calloc(sizeof(unsigned char),(long)tile_Final_DEMsize.width*(long)tile_Final_DEMsize.height);
                                    NNA_M(proinfo, param, iter_row_start, t_col_start, iter_row_end, t_col_end, buffer_tile, final_iteration, tile_row, tile_Final_DEMsize, DEM_values, H_value, MT_value, FinalDEM_boundary);
                                    
                                    ET = time(0);
                                    gap = difftime(ET,ST);
                                    printf("DEM finish(time[m] = %5.2f)!!\n",gap/60.0);
                                    
                                    double MT_memory = CalMemorySize_Post_MT(tile_Final_DEMsize,tile_Final_DEMsize);
                                    if(MT_memory > proinfo->System_memory - 5)
                                    {
                                        ST = time(0);
                                        
                                        printf("not enough memory for matchtag filtering[%f]!!\nMake original matchtag file!!\n",MT_memory);
                                        char GEOTIFF_matchtag_filename[500];
                                        sprintf(GEOTIFF_matchtag_filename, "%s/%s_%d_matchtag.tif", proinfo->save_filepath, proinfo->Outputpath_name,tile_row);
                                        
                                        WriteGeotiff(GEOTIFF_matchtag_filename, MT_value, tile_Final_DEMsize.width, tile_Final_DEMsize.height, proinfo->DEM_resolution, FinalDEM_boundary[0], FinalDEM_boundary[3], param.projection, param.utm_zone, param.bHemisphere, 1);
                                        
                                        free(H_value);
                                        free(MT_value);
                                        
                                        ET = time(0);
                                        gap = difftime(ET,ST);
                                        printf("Matched PT finish(time[m] = %5.2f)!!\n",gap/60.0);
                                    }
                                    else
                                    {
                                        ST = time(0);
                                        
                                        Ortho_values = (signed char*)malloc(sizeof(signed char)*tile_Final_DEMsize.width*tile_Final_DEMsize.height);
                                        MergeTiles_Ortho(proinfo,iter_row_start,t_col_start,iter_row_end,t_col_end,buffer_tile,final_iteration,Ortho_values,tile_Final_DEMsize,FinalDEM_boundary);
                                        
                                        NNA_M_MT(proinfo, param, iter_row_start,t_col_start, iter_row_end, t_col_end, buffer_tile, final_iteration, tile_row, Ortho_values, H_value, MT_value, tile_Final_DEMsize, FinalDEM_boundary);
                                        
                                        ET = time(0);
                                        gap = difftime(ET,ST);
                                        printf("Matched PT finish(time[m] = %5.2f)!!\n",gap/60.0);
                                    }
                                    
                                    if(args.check_LSF2 == 1 || args.check_LSF2 == 2)
                                    {
                                        ST = time(0);
                                        
                                        LSFSmoothing_DEM(proinfo->save_filepath,proinfo->Outputpath_name,MPP_stereo_angle,tile_row);
                                        
                                        ET = time(0);
                                        gap = difftime(ET,ST);
                                        printf("LSF finish(time[m] = %5.2f)!!\n",gap/60.0);
                                    }
                                    
                                    CSize seeddem_size;
                                    double tminX, tmaxY;
                                    
                                    char tiff_path[500];
                                    sprintf(tiff_path, "%s/%s_%d_dem.tif", proinfo->save_filepath, proinfo->Outputpath_name,tile_row);
                                    seeddem_size = ReadGeotiff_info(tiff_path, &tminX, &tmaxY, NULL);
                                    
                                    fprintf(pMetafile,"DEM %d\n",tile_row);
                                    fprintf(pMetafile,"Output dimensions=%d\t%d\n",seeddem_size.width,seeddem_size.height);
                                    fprintf(pMetafile,"Upper left coordinates=%f\t%f\n",tminX,tmaxY);
                                }
                                fclose(p_sefile);
                            }
                            else
                            {
                                sprintf(str_DEMfile, "%s/%s_dem.tif", proinfo->save_filepath,proinfo->Outputpath_name);
                                
                                FILE* pFile_DEM = NULL;
                                
                                pFile_DEM = fopen(str_DEMfile,"r");
                                printf("check exist %s %d\n",str_DEMfile,!!pFile_DEM);
                                
                                ST = time(0);
                                printf("Tile merging start final iteration %d!!\n",final_iteration);
                                
                                DEM_values = (float*)malloc(sizeof(float)*Final_DEMsize.width*Final_DEMsize.height);
                                MergeTiles(proinfo,iter_row_start,t_col_start,iter_row_end,t_col_end,buffer_tile,final_iteration,DEM_values,Final_DEMsize,FinalDEM_boundary);
                                printf("%f\t",DEM_values[10]);
                                
                                printf("Interpolation start!!\n");
                                printf("%f %f\n",proinfo->DEM_resolution,proinfo->DEM_resolution);
                                
                                H_value = (float*)malloc(sizeof(float)*(long)Final_DEMsize.width*(long)Final_DEMsize.height);
                                MT_value = (unsigned char*)calloc(sizeof(unsigned char),(long)Final_DEMsize.width*(long)Final_DEMsize.height);
                                NNA_M(proinfo, param, iter_row_start, t_col_start, iter_row_end, t_col_end, buffer_tile, final_iteration, 0, Final_DEMsize, DEM_values, H_value, MT_value, FinalDEM_boundary);
                                
                                ET = time(0);
                                gap = difftime(ET,ST);
                                printf("DEM finish(time[m] = %5.2f)!!\n",gap/60.0);
                                
                                double MT_memory = CalMemorySize_Post_MT(Final_DEMsize,Final_DEMsize);
                                //MT_memory = 100;
                                if(MT_memory > proinfo->System_memory - 5)
                                {
                                    ST = time(0);
                                    
                                    printf("not enough memory for matchtag filtering[%f]!!\nMake original matchtag file!!\n",MT_memory);
                                    char GEOTIFF_matchtag_filename[500];
                                    sprintf(GEOTIFF_matchtag_filename, "%s/%s_matchtag.tif", proinfo->save_filepath, proinfo->Outputpath_name);
                                    
                                    WriteGeotiff(GEOTIFF_matchtag_filename, MT_value, Final_DEMsize.width, Final_DEMsize.height, proinfo->DEM_resolution, FinalDEM_boundary[0], FinalDEM_boundary[3], param.projection, param.utm_zone, param.bHemisphere, 1);
                                    
                                    free(H_value);
                                    free(MT_value);
                                    
                                    ET = time(0);
                                    gap = difftime(ET,ST);
                                    printf("Matched PT finish(time[m] = %5.2f)!!\n",gap/60.0);
                                }
                                else
                                {
                                    ST = time(0);
                                    
                                    Ortho_values = (signed char*)malloc(sizeof(signed char)*Final_DEMsize.width*Final_DEMsize.height);
                                    MergeTiles_Ortho(proinfo,iter_row_start,t_col_start,iter_row_end,t_col_end,buffer_tile,final_iteration,Ortho_values,Final_DEMsize,FinalDEM_boundary);
                                    
                                    NNA_M_MT(proinfo, param, iter_row_start, t_col_start, iter_row_end, t_col_end, buffer_tile, final_iteration, 0, Ortho_values, H_value, MT_value, Final_DEMsize, FinalDEM_boundary);
                                    
                                    ET = time(0);
                                    gap = difftime(ET,ST);
                                    printf("Matched PT finish(time[m] = %5.2f)!!\n",gap/60.0);
                                }
                                
                                CSize seeddem_size;
                                double tminX, tmaxY;
                                
                                char tiff_path[500];
                                if(proinfo->check_Matchtag)
                                    sprintf(tiff_path, "%s", proinfo->priori_DEM_tif);
                                else
                                    sprintf(tiff_path, "%s/%s_dem.tif", proinfo->save_filepath, proinfo->Outputpath_name);
                                seeddem_size = ReadGeotiff_info(tiff_path, &tminX, &tmaxY, NULL);
                                fprintf(pMetafile,"Output dimensions=%d\t%d\n",seeddem_size.width,seeddem_size.height);
                                fprintf(pMetafile,"Upper left coordinates=%f\t%f\n",tminX,tmaxY);
                                
                                if (args.check_minH)
                                    fprintf(pMetafile,"user_defined_minH=%f\n",args.minHeight);
                                if (args.check_maxH)
                                    fprintf(pMetafile,"user_defined_maxH=%f\n",args.maxHeight);
                                
                                if(args.check_LSF2 == 1 || args.check_LSF2 == 2)
                                {
                                    ST = time(0);
                                    
                                    LSFSmoothing_DEM(proinfo->save_filepath,proinfo->Outputpath_name,MPP_stereo_angle,0);
                                    
                                    ET = time(0);
                                    gap = difftime(ET,ST);
                                    printf("LSF finish(time[m] = %5.2f)!!\n",gap/60.0);
                                }
                            }
                            
                            for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                            {
                                fprintf(pMetafile,"Image %d info\nImage_%d_satID=%s\nImage_%d_Acquisition_time=%s\nImage_%d_Mean_row_GSD=%f\nImage_%d_Mean_col_GSD=%f\nImage_%d_Mean_GSD=%f\nImage_%d_Mean_sun_azimuth_angle=%f\nImage_%d_Mean_sun_elevation=%f\nImage_%d_Mean_sat_azimuth_angle=%f\nImage_%d_Mean_sat_elevation=%f\nImage_%d_Intrack_angle=%f\nImage_%d_Crosstrack_angle=%f\nImage_%d_Offnadir_angle=%f\nImage_%d_tdi=%d\nImage_%d_effbw=%f\nImage_%d_abscalfact=%f\n",ti+1,ti+1,image_info[ti].SatID,ti+1,image_info[ti].imagetime,ti+1,Image_gsd_r[ti],ti+1,Image_gsd_c[ti],ti+1,Image_gsd[ti],ti+1,image_info[ti].Mean_sun_azimuth_angle,ti+1,image_info[ti].Mean_sun_elevation,ti+1,image_info[ti].Mean_sat_azimuth_angle,ti+1,image_info[ti].Mean_sat_elevation,ti+1,image_info[ti].Intrack_angle,ti+1,image_info[ti].Crosstrack_angle,ti+1,image_info[ti].Offnadir_angle,ti+1,(int)leftright_band[ti].tdi,ti+1,leftright_band[ti].effbw,ti+1,leftright_band[ti].abscalfactor);
                            }
                            
                            fprintf(pMetafile,"Stereo_pair_convergence_angle=%f\n",convergence_angle);
                            fprintf(pMetafile,"Stereo_pair_expected_height_accuracy=%f\n",MPP_stereo_angle);
                            fclose(pMetafile);
                        } // if (!RA_only)
                    }
                    else
                        printf("out of boundary!! please check boundary infomation!!\n");
                    
                    sprintf(temp_filepath,"%s/tmp",proinfo->save_filepath);
                }

                if(proinfo->sensor_type != AB)
                {
                    for(int i = 0; i < proinfo->number_of_images; i++)
                        RPCsFree(RPCs[i]);
                    free(RPCs);
                }
                free(leftright_band);
                free(Image_gsd_r);
                free(Image_gsd_c);
                free(Image_gsd);
                free(Limagesize);
                free(image_info);
            }
            else
                printf("Check output directory path!!\n");
        }
    }
    
    if(!cal_check && args.check_sdm_ortho == 2)
    {
        printf("No start computation time\n");
    }
    else
    {
        printf("start computation time\n");
        total_ET = time(0);
        total_gap = difftime(total_ET,total_ST);
        
        sprintf(computation_file,"%s/txt/computation_time.txt",proinfo->save_filepath);
        time_fid            = fopen(computation_file,"w");
        fprintf(time_fid,"Computation_time[m] = %5.2f\n",total_gap/60.0);
        printf("Computation_time[m] = %5.2f\n",total_gap/60.0);
        
        fclose(time_fid);
    }
    
    delete proinfo;

    if(args.check_coreg == 1)
        exit(1);

    if(args.check_coreg == 2 || args.check_coreg == 3)
        exit(1);
    else if(args.check_sdm_ortho == 1)
        exit(1);
    else if(args.check_sdm_ortho == 2)
        exit(1);
    else if(args.check_ortho)
        exit(1);
    
#ifdef BUILDMPI
    // Make sure to finalize
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
    {
        MPI_Finalize();
    }
#endif
    
    return DEM_divide;
}

#ifdef BUILDMPI
// Reorder tiles for static load balancing with MPI
int reorder_list_of_tiles(int iterations[], int length, int col_length, int row_length)
{
    int i,j;

    int *temp = (int*)malloc(col_length*row_length*2*sizeof(int));
    int midrow = ceil(row_length / 2.0);
    int midcol = ceil(col_length / 2.0);

    for (i = 0; i < length*2; i += 2)
    {
        int closest = 0;
        int closest_dist = midrow + midcol;
        for (j = 0; j < length*2; j += 2)
        {
            int new_dist = abs(midrow - iterations[j]) + abs(midcol - iterations[j+1]);
            if (new_dist <= closest_dist)
            {
                closest_dist = new_dist;
                closest = j;
            }
        }
        temp[i] = iterations[closest];
        temp[i+1] = iterations[closest+1];
        iterations[closest] = -1;
        iterations[closest+1] = -1;
    }

    for (i = 0; i < length*2; i++)
    {
        iterations[i] = temp[i];
    }
    free(temp);
    return length;
}
#endif

void SetPairs(ProInfo *proinfo, PairInfo &pairinfo, const ImageInfo *image_info)
{
    int pair_number = 0;
    float minBH = 100;
    float maxBH = -100;
    for(int ref_ti = 0 ; ref_ti < proinfo->number_of_images - 1 ; ref_ti++)
    {
        for(int ti = ref_ti + 1 ; ti < proinfo->number_of_images ; ti ++)
            pair_number++;
    }
    pairinfo.NumberOfPairs = pair_number;
    
    pairinfo.pairs = (UI2DPOINT*)malloc(sizeof(UI2DPOINT)*pair_number);
    pairinfo.BHratio = (float*)malloc(sizeof(float)*pair_number);
    pairinfo.ConvergenceAngle = (float*)malloc(sizeof(float)*pair_number);
    
    pair_number = 0;
    for(int ref_ti = 0 ; ref_ti < proinfo->number_of_images - 1 ; ref_ti++)
    {
        for(int ti = ref_ti + 1 ; ti < proinfo->number_of_images ; ti ++)
        {
            pairinfo.pairs[pair_number].m_X = ref_ti;
            pairinfo.pairs[pair_number].m_Y = ti;
            
            if(proinfo->sensor_type == AB)
            {
                D3DPOINT image1_PL(proinfo->frameinfo.Photoinfo[ref_ti].m_Xl,proinfo->frameinfo.Photoinfo[ref_ti].m_Yl,proinfo->frameinfo.Photoinfo[ref_ti].m_Zl);
                D3DPOINT image2_PL(proinfo->frameinfo.Photoinfo[ti].m_Xl,proinfo->frameinfo.Photoinfo[ti].m_Yl,proinfo->frameinfo.Photoinfo[ti].m_Zl);
                double PC_distance = SQRT(image1_PL,image2_PL,3);
                pairinfo.BHratio[pair_number] = PC_distance/((proinfo->frameinfo.Photoinfo[ref_ti].m_Zl + proinfo->frameinfo.Photoinfo[ti].m_Zl)/2.0);
            }
            else
            {
                if(proinfo->sensor_provider == DG)
                {
                    double convergence_angle = acos(sin(image_info[ref_ti].Mean_sat_elevation*DegToRad)*sin(image_info[ti].Mean_sat_elevation*DegToRad) + cos(image_info[ref_ti].Mean_sat_elevation*DegToRad)*cos(image_info[ti].Mean_sat_elevation*DegToRad)*cos( (image_info[ref_ti].Mean_sat_azimuth_angle - image_info[ti].Mean_sat_azimuth_angle)*DegToRad))*RadToDeg;
                    
                    pairinfo.ConvergenceAngle[pair_number] = convergence_angle;
                    pairinfo.BHratio[pair_number] = 2.0*tan(convergence_angle*DegToRad*0.5);
                }
                else if(proinfo->sensor_provider = PT)
                {
                    pairinfo.BHratio[pair_number] = 0.5;
                }
            }
            
            if(minBH > pairinfo.BHratio[pair_number])
                minBH = pairinfo.BHratio[pair_number];
            
            if(maxBH < pairinfo.BHratio[pair_number])
                maxBH = pairinfo.BHratio[pair_number];
            
            printf("pairnumber %d\timage %d\t%d\tBHratio %f\t%f\t%f\n",pair_number,pairinfo.pairs[pair_number].m_X,pairinfo.pairs[pair_number].m_Y,pairinfo.BHratio[pair_number],minBH,maxBH);
            
            pair_number++;
        }
    }
    /*
     for(int count=0;count<pairinfo.NumberOfPairs;count++)
     {
     pairinfo.BHratio[count] =  (pairinfo.BHratio[count] - minBH)/(maxBH - minBH)*0.5 + 0.5;
     printf("pairnumber %d\timage %d\t%d\tBHratio %f\t%f\t%f\n",count,pairinfo.pairs[count].m_X,pairinfo.pairs[count].m_Y,pairinfo.BHratio[count],minBH,maxBH);
     }
     */
}

void actual_pair(const ProInfo *proinfo, LevelInfo &plevelinfo, double *minmaxHeight, vector<unsigned char> &save_pair, PairInfo &pairinfo)
{
    save_pair.clear();
    
    for(long int iter_count = 0 ; iter_count < (*plevelinfo.Grid_length) ; iter_count++)
    {
        long int pt_index = iter_count;
        const int start_H     = minmaxHeight[0];
        const int end_H       = minmaxHeight[1];
        int select_pair = select_referenceimage(pt_index, proinfo, plevelinfo, start_H, end_H);

        if(save_pair.size() > 0)
        {
            bool check_stop = false;
            int count = 0;
            while(count < save_pair.size() && !check_stop)
            {
                if(save_pair[count] == select_pair)
                    check_stop = true;
                count++;
            }
            if(!check_stop)
            {
                save_pair.push_back(select_pair);
                //for(int i = 0 ; i < save_pair.size() ; i++)
                //    printf("save_pair.push_back(select_pair) %d\t%d\t%d\n",save_pair.size(),save_pair[i],select_pair);
            }
        }
        else
        {
            save_pair.push_back(select_pair);
            //printf("save_pair.push_back(select_pair) %d\n",select_pair);
        }
    }
    
    vector<unsigned char> actual_pair_save;
    for(int pair_number = 0 ; pair_number < plevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        const int reference_id = plevelinfo.pairinfo->pairs[pair_number].m_X;
        const int ti = plevelinfo.pairinfo->pairs[pair_number].m_Y;
        
        bool check_stop = false;
        int count = 0;
        while(count < save_pair.size() && !check_stop)
        {
            if(reference_id == save_pair[count] || ti == save_pair[count])
                check_stop = true;
            count++;
        }
        if(check_stop)
        {
            actual_pair_save.push_back(pair_number);
            
            //printf("pair number %d\t%d\t%d\n",pair_number,reference_id,ti);
        }
    }
    
    plevelinfo.pairinfo->NumberOfPairs = actual_pair_save.size();
    
    pairinfo.NumberOfPairs = actual_pair_save.size();
    pairinfo.pairs = (UI2DPOINT*)malloc(sizeof(UI2DPOINT)*pairinfo.NumberOfPairs);
    pairinfo.BHratio = (float*)malloc(sizeof(float)*(pairinfo.NumberOfPairs));
    pairinfo.ConvergenceAngle = (float*)malloc(sizeof(float)*(pairinfo.NumberOfPairs));
    
    UI2DPOINT *temp_pairs = (UI2DPOINT*)calloc(sizeof(UI2DPOINT),plevelinfo.pairinfo->NumberOfPairs);
    float *temp_BHratio = (float*)calloc(sizeof(float),plevelinfo.pairinfo->NumberOfPairs);
    float *temp_CA = (float*)calloc(sizeof(float),plevelinfo.pairinfo->NumberOfPairs);
    
    for(int count = 0 ; count < plevelinfo.pairinfo->NumberOfPairs ; count++)
    {
        int pair_number = actual_pair_save[count];
        temp_pairs[count].m_X = plevelinfo.pairinfo->pairs[pair_number].m_X;
        temp_pairs[count].m_Y = plevelinfo.pairinfo->pairs[pair_number].m_Y;
        temp_BHratio[count] = plevelinfo.pairinfo->BHratio[pair_number];
        temp_CA[count] = plevelinfo.pairinfo->ConvergenceAngle[pair_number];
    }
    
    for(int count = 0 ; count < plevelinfo.pairinfo->NumberOfPairs ; count++)
    {
        plevelinfo.pairinfo->pairs[count].m_X = temp_pairs[count].m_X;
        plevelinfo.pairinfo->pairs[count].m_Y = temp_pairs[count].m_Y;
        plevelinfo.pairinfo->BHratio[count] = temp_BHratio[count];
        plevelinfo.pairinfo->ConvergenceAngle[count] = temp_CA[count];
        
        pairinfo.pairs[count].m_X = temp_pairs[count].m_X;
        pairinfo.pairs[count].m_Y = temp_pairs[count].m_Y;
        pairinfo.BHratio[count] = temp_BHratio[count];
        pairinfo.ConvergenceAngle[count] = temp_CA[count];
        //printf("count %d\t%d\t%d\t%f\n",count, plevelinfo.pairinfo->pairs[count].m_X,plevelinfo.pairinfo->pairs[count].m_Y,plevelinfo.pairinfo->BHratio[count]);
    }
    
    free(temp_pairs);
    free(temp_BHratio);
    free(temp_CA);
}

void findOverlappArea(ProInfo *proinfo, TransParam param, double*** RPCs, double *Image_res, double Boundary[])
{
    Boundary[0] = 10000000;
    Boundary[1] = 10000000;
    Boundary[2] = -10000000;
    Boundary[3] = -10000000;
    
    double LBoundary[4],LminmaxHeight[2];
    double LHinterval;
    
    
    for(int ref_ti = 0 ; ref_ti < proinfo->number_of_images - 1 ; ref_ti++)
    {
        double lonlatboundary_ref[4] = {0.0};
        
        if(proinfo->sensor_type == SB)
            SetDEMBoundary(proinfo,RPCs[ref_ti],Image_res,param,LBoundary,LminmaxHeight,&LHinterval);
        else
            SetDEMBoundary_photo(proinfo->frameinfo.Photoinfo[ref_ti], proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ref_ti].m_Rm, LBoundary,LminmaxHeight,&LHinterval);
        
        for(int i=0;i<4;i++)
            lonlatboundary_ref[i] = LBoundary[i];
        
        for(int ti = ref_ti + 1 ; ti < proinfo->number_of_images ; ti ++)
        {
            double lonlatboundary[4] = {0.0};
            if(proinfo->sensor_type == SB)
                SetDEMBoundary(proinfo,RPCs[ti],Image_res,param,LBoundary,LminmaxHeight,&LHinterval);
            else
                SetDEMBoundary_photo(proinfo->frameinfo.Photoinfo[ti], proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm, LBoundary,LminmaxHeight,&LHinterval);
            
            printf("tar %d\tboundary %f\t%f\t%f\t%f\n",ti,LBoundary[0],LBoundary[1],LBoundary[2],LBoundary[3]);
            
            for(int i=0;i<4;i++)
            {
                if(i<2)
                    lonlatboundary[i] = max(LBoundary[i], lonlatboundary_ref[i]);
                else
                    lonlatboundary[i] = min(LBoundary[i], lonlatboundary_ref[i]);
            }
            
            printf("ref %d tar %d\tboundary %f\t%f\t%f\t%f\n",ref_ti,ti,lonlatboundary[0],lonlatboundary[1],lonlatboundary[2],lonlatboundary[3]);
            
            for(int i=0;i<4;i++)
            {
                if(i<2)
                    Boundary[i] = min(Boundary[i], lonlatboundary[i]);
                else
                    Boundary[i] = max(Boundary[i], lonlatboundary[i]);
            }
            
            printf("all boundary %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
        }
    }
}


int Matching_SETSM(ProInfo *proinfo,const ImageInfo *image_info, const uint8 pyramid_step, const uint8 Template_size, const uint16 buffer_area, const uint8 iter_row_start, const uint8 iter_row_end, const uint8 t_col_start, const uint8 t_col_end, const double subX,const double subY,const double bin_angle,const double Hinterval,const double *Image_res, double **Imageparams, const double *const*const*RPCs, const uint8 NumOfIAparam, const CSize *Imagesizes,const TransParam param, double *ori_minmaxHeight,const double *Boundary, const double CA,const double mean_product_res, double *stereo_angle_accuracy, PairInfo &pairinfo_return)
{
#ifdef BUILDMPI
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    int final_iteration = -1;
    int row,col;
    int* RA_count = (int*)calloc(sizeof(int),MaxNCC);
    
    int row_length = iter_row_end-iter_row_start;
    int col_length = t_col_end-t_col_start;
    int *iterations = (int*)malloc(col_length*row_length*2*sizeof(int));
    int length = 0;
    int count_tri;
    
    for(row = iter_row_start; row < iter_row_end ; row++)
    {
        for(col = t_col_start ; col < t_col_end ; col++)
        {
            iterations[2*length] = row;
            iterations[2*length+1] = col;
            length+=1;
        }
    }
    
#ifdef BUILDMPI
    //Reorder list of tiles for static load balancing
    if (length > 1) {
        reorder_list_of_tiles(iterations, length, col_length, row_length);
    }
#endif
    
    int tile_iter, i;
    for(tile_iter = 0; tile_iter < length; tile_iter += 1)
    {
        row = iterations[2*tile_iter];
        col = iterations[2*tile_iter+1];
        
#ifdef BUILDMPI
        // Skip this tile if it belongs to a different MPI rank
        if (tile_iter % size != rank)
            continue;
        printf("MPI: Rank %d is analyzing row %d, col %d\n", rank, row, col);
#endif
        bool check_cal = false;
        bool check_cal_2 = false;
        if(proinfo->IsRA)
            check_cal = true;
        else 
        {
            char check_file[500];
            sprintf(check_file,"%s/txt/matched_pts_%d_%d_0_3.txt",proinfo->save_filepath,row,col);
            FILE* pcheckFile = fopen(check_file,"r");
            if(!pcheckFile)
                check_cal = true;
            else
                fclose(pcheckFile);
            
            sprintf(check_file,"%s/txt/matched_BR_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
            FILE* pcheckFile_2 = fopen(check_file,"r");
            if(!pcheckFile_2)
                check_cal_2 = true;
            else
                fclose(pcheckFile_2);
            
            printf("check existing tile reuslt %d\t%d\t%d\t%d\n",row,col,check_cal,check_cal_2);
        }
        
        if(check_cal || check_cal_2)
        {
            printf("start cal tile\n");
            
            bool temp_asc_fprint = false;
            
            FILE *fid = NULL;
            FILE *fid_header = NULL;
            
            D2DPOINT *Startpos_ori = (D2DPOINT*)calloc(sizeof(D2DPOINT),proinfo->number_of_images);
            CSize *Subsetsize = (CSize*)calloc(sizeof(CSize),proinfo->number_of_images);
            
            double **t_Imageparams = (double**)calloc(sizeof(double*),MaxNCC);
            for(int ti = 0 ; ti < MaxNCC ; ti++)
                t_Imageparams[ti] = (double*)calloc(sizeof(double),2);
            
            char save_file[500];
            if(proinfo->IsRA)
            {
                if(!proinfo->check_checktiff)
                {
                    sprintf(save_file,"%s/txt/RA_echo_result_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    fid         = fopen(save_file,"w");
                    for(int ti = 0 ; ti < MaxNCC ; ti++)
                        fprintf(fid,"RA param X = %f\tY = %f\t0\t0\t0.0\n",t_Imageparams[ti][0],t_Imageparams[ti][1]);
                    
                    sprintf(save_file,"%s/txt/RA_headerinfo_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    fid_header  = fopen(save_file,"w");
                    
                    for(int ti = 0 ; ti < MaxNCC ; ti++)
                    {
                        Imageparams[ti][0] = 0;
                        Imageparams[ti][1] = 0;
                    }
                }
            }
            else
            {
                for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
                {
                    t_Imageparams[ti][0]    = Imageparams[ti][0];
                    t_Imageparams[ti][1]    = Imageparams[ti][1];
                }
                
                if(!proinfo->check_checktiff)
                {
                    sprintf(save_file,"%s/txt/echo_result_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    fid         = fopen(save_file,"w");
                    for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
                        fprintf(fid,"RA param X = %f\tY = %f\t%d\t%d\t%f\n",t_Imageparams[ti][0],t_Imageparams[ti][1],pairinfo_return.pairs[ti].m_X,pairinfo_return.pairs[ti].m_Y,pairinfo_return.BHratio[ti]);
                    
                    sprintf(save_file,"%s/txt/headerinfo_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    fid_header  = fopen(save_file,"w");
                }
            }
            
            double subBoundary[4];
            double minmaxHeight[2] = {ori_minmaxHeight[0], ori_minmaxHeight[1]};
            printf("minmaxH = %f\t%f\n",minmaxHeight[0],minmaxHeight[1]);
            SetSubBoundary(Boundary,subX,subY,buffer_area,col,row,subBoundary);
            
            printf("subBoundary = %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]);
            
            printf("subsetimage\n");
            
            if(proinfo->IsRA)
            {
                double RA_memory = 0;
                for(int ti = 0; ti < MaxNCC ; ti++)
                {
                    printf("RA Params=%f\t%f\t\n",Imageparams[ti][0],Imageparams[ti][1]);
                }
                
                for(int ti = 0; ti < proinfo->number_of_images ; ti++)
                {
                    CSize t_Imagesize((subBoundary[2] - subBoundary[0])/0.7, (subBoundary[3] - subBoundary[1])/0.7);
                    long int data_length =(long int)t_Imagesize.width*(long int)t_Imagesize.height;
                    RA_memory += (sizeof(uint16)*data_length)*2;
                    RA_memory += (sizeof(int16)*data_length);
                    RA_memory += (sizeof(uint8)*data_length);
                    
                    proinfo->check_selected_image[ti] = true;
                }
                
                RA_memory = RA_memory/1024.0/1024.0/1024.0;
            
                printf("RPC bias calculation required Memory : System %f\t SETSM required %f\n",proinfo->System_memory,RA_memory);
                if(RA_memory > proinfo->System_memory - 2)
                {
                    printf("System memory is not enough to run a relative RPC bias computation module of SETSM. Please reduce RA tilesize or assign more physical memory!!\n");
                    exit(1);
                }
            }
            else
            {
                for(int ti = 0; ti < pairinfo_return.NumberOfPairs ; ti++)
                {
                    printf("RA Params=%f\t%f\t\n",Imageparams[ti][0],Imageparams[ti][1]);
                }
            }
            
            time_t PreST = 0, PreET = 0;
            double Pregab = 0;
            
            LevelInfo levelinfo = {NULL};
            levelinfo.RPCs = RPCs;
            levelinfo.Boundary = subBoundary;
            levelinfo.Template_size = &Template_size;
            levelinfo.ImageAdjust = t_Imageparams;
            levelinfo.param = &param;
            levelinfo.NumOfIAparam = &NumOfIAparam;
            levelinfo.Imagesize_ori = Imagesizes;
            levelinfo.imageinfo = image_info;
            
            uint16 **SourceImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
            int count_available_images = 0;
            for(int index_image = 0 ; index_image < proinfo->number_of_images ; index_image++)
            {
                SourceImages[index_image] = SetsubsetImage(proinfo, levelinfo, index_image,param,NumOfIAparam,RPCs,t_Imageparams,subBoundary,minmaxHeight,Startpos_ori,Subsetsize);
                if(proinfo->check_selected_image[index_image])
                    count_available_images++;
            }
            
            if(count_available_images >= 2)
            {
                printf("Completion of subsetImage!!\n");
                
                if( check_kernel_size(proinfo, Subsetsize, Template_size, pyramid_step))
                {
                    PairInfo pairinfo;
                    SetPairs(proinfo,pairinfo,image_info);
                    
                    levelinfo.pairinfo = &pairinfo;
                    
                    for(int kk = 0 ; kk < levelinfo.pairinfo->NumberOfPairs ; kk++)
                    {
                        printf("ref\t%d\ttar\t%d\n",levelinfo.pairinfo->pairs[kk].m_X,levelinfo.pairinfo->pairs[kk].m_Y);
                    }
                    
                    double total_memory = 0.0;
                    double py_resolution = 0;
                    double grid_resolution = 0;
                    double pre_grid_resolution = 0;
                    bool lower_level_match = true;
                    
                    bool flag_start = false;
                    
                    int level             = pyramid_step;
                    if(proinfo->check_Matchtag)
                    {
                        level   = 0;
                        printf("reprocessing Matchtag\n");
                    }
                    
                    CSize Size_Grid2D(0,0), pre_Size_Grid2D(0,0);
                    CSize **data_size_lr = (CSize**)malloc(sizeof(CSize*)*proinfo->number_of_images);
                    UGRID *GridPT3 = NULL, *Pre_GridPT3 = NULL;
                     
                    if(!proinfo->check_Matchtag)
                    {
                        for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                        {
                            if(proinfo->check_selected_image[ti])
                            {
                                data_size_lr[ti] = (CSize*)malloc(sizeof(CSize)*(level+1));
                                SetPySizes(data_size_lr[ti], Subsetsize[ti], level);
                                
                                for (int ttt = 0 ; ttt < level+1 ;ttt++)
                                    printf("data_size %d\t%d\n",data_size_lr[ti][ttt].width,data_size_lr[ti][ttt].height);
                            }
                        }
                    }
                    else
                    {
                        for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                        {
                            if(proinfo->check_selected_image[ti])
                            {
                                data_size_lr[ti] = (CSize*)malloc(sizeof(CSize)*(level+2));
                                SetPySizes(data_size_lr[ti], Subsetsize[ti], level+1);
                                for (int ttt = 0 ; ttt < level+2 ;ttt++)
                                    printf("data_size %d\t%d\n",data_size_lr[ti][ttt].width,data_size_lr[ti][ttt].height);
                            }
                        }
                    }
                    
                    PreST = time(0);
                    printf("row = %d/%d\tcol = %d/%d\tPreprocessing start!!\n",row,iter_row_end,col,t_col_end);
                    
                    //Set Pyramid Images memory
                    uint8 py_level_set;
                    if(proinfo->check_Matchtag)
                        py_level_set = 2;
                    else
                        py_level_set = level + 1;
                    
                    uint16 ***SubImages = (uint16***)malloc(sizeof(uint16**)*py_level_set);
                    uint8  ***SubOriImages = (uint8***)malloc(sizeof(uint8**)*py_level_set);
                    uint16 ***SubMagImages = (uint16***)malloc(sizeof(uint16**)*py_level_set);
                    
                    for(int iter_level = 0 ; iter_level < py_level_set; iter_level++)
                    {
                        SubImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
                        SubOriImages[iter_level] = (uint8**)malloc(sizeof(uint8*)*proinfo->number_of_images);
                        SubMagImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
                        
                        for(int image_index = 0 ; image_index < proinfo->number_of_images ; image_index++)
                        {
                            if(proinfo->check_selected_image[image_index])
                            {
                                long int data_length = (long int)data_size_lr[image_index][iter_level].height*(long int)data_size_lr[image_index][iter_level].width;
                                
                                if(iter_level == 0)
                                    SubImages[iter_level][image_index] = SourceImages[image_index];
                                else
                                    SubImages[iter_level][image_index] = NULL;//(uint16*)malloc(sizeof(uint16)*data_length);
                                
                                SubOriImages[iter_level][image_index] = (uint8*)malloc(sizeof(uint8)*data_length);
                                SubMagImages[iter_level][image_index] = (uint16*)malloc(sizeof(uint16)*data_length);
                            }
                        }
                    }
                    //pyramid image generation
                    SetPyramidImages(proinfo, py_level_set, data_size_lr, SubImages, SubMagImages, SubOriImages);
                    
                    PreET = time(0);
                    Pregab = difftime(PreET,PreST);
                    printf("row = %d/%d\tcol = %d/%d\tPreprocessing finish(time[m] = %5.2f)!!\n",row,iter_row_end,col,t_col_end,Pregab/60.0);
                    
                    PreST = time(0);
                    printf("row = %d/%d\tcol = %d/%d\tDEM generation start!!\n",row,iter_row_end,col,t_col_end);
                    
                    int blunder_selected_level;
                    
                    int final_level_iteration = 1;
                    if(proinfo->check_Matchtag)
                        final_level_iteration = 2;
                    
                    double matching_rate = 0;
                    bool check_matching_rate = false;
                    const double th_mr = 0.05;
                    
                    bool check_RA_divide = false;
                    const double tilesize_RA = 20000;
                    const double lengthOfX = subBoundary[2] - subBoundary[0];
                    const double lengthOfY = subBoundary[3] - subBoundary[1];
                    int division_X = 0, division_Y = 0;
                    int total_tile = 0;
                    double new_subBoundary_RA[4];
                    bool check_new_subBoundary_RA = false;
                    double preBoundary[4] = {0};
                    const double coverage = lengthOfY*lengthOfX/1000000.0;
                    if(coverage > tilesize_RA*tilesize_RA/1000000.0)
                    {
                        if(lengthOfX < lengthOfY)
                        {
                            if(lengthOfY > tilesize_RA)
                            {
                                check_RA_divide = true;
                                division_X = (int) (ceil(lengthOfX / (tilesize_RA)));
                                division_Y = (int) (ceil(lengthOfY / (tilesize_RA)));
                                total_tile = division_X*division_Y;
                            }
                        }
                        else
                        {
                            if(lengthOfX > tilesize_RA)
                            {
                                check_RA_divide = true;
                                division_X = (int) (ceil(lengthOfX / (tilesize_RA)));
                                division_Y = (int) (ceil(lengthOfY / (tilesize_RA)));
                                total_tile = division_X*division_Y;
                            }
                        }
                    }
                    
                    printf("length %f\t%f\tdivision %d\t%d\ntotal_tile %d\tcheck_RA_divide %d\n", lengthOfX, lengthOfY, division_X, division_Y, total_tile, check_RA_divide);
                    
                    const int Py_combined_level = 0;
                    int RA_resize_level = 0;
                    vector<unsigned char> compute_pair;
                    
                    while(lower_level_match && level >= 0)
                    {
                        printf("level = %d\t final_level_iteration %d\n",level,final_level_iteration);
                        
                        for(int ti = 0; ti < pairinfo_return.NumberOfPairs ; ti++)
                        {
                            printf("RA Params=%f\t%f\t%f\t%f\n",t_Imageparams[ti][0],t_Imageparams[ti][1],levelinfo.ImageAdjust[ti][0],levelinfo.ImageAdjust[ti][1]);
                        }
                        
                        if(proinfo->IsRA && check_new_subBoundary_RA)
                        {
                            printf("Resize RA tile\n");
                            //delete pre-asigned image memory
                            for(int t_level = 0 ; t_level < py_level_set ; t_level++)
                            {
                                for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                                {
                                    if(proinfo->check_selected_image[ti])
                                    {
                                        free(SubImages[t_level][ti]);
                                        free(SubMagImages[t_level][ti]);
                                        free(SubOriImages[t_level][ti]);
                                    }
                                }
                                free(SubImages[t_level]);
                                free(SubMagImages[t_level]);
                                free(SubOriImages[t_level]);
                            }
                            free(SubImages);
                            free(SubMagImages);
                            free(SubOriImages);
                            free(SourceImages);
                            
                            preBoundary[0] = subBoundary[0];
                            preBoundary[1] = subBoundary[1];
                            preBoundary[2] = subBoundary[2];
                            preBoundary[3] = subBoundary[3];
                            
                            subBoundary[0] = new_subBoundary_RA[0];
                            subBoundary[1] = new_subBoundary_RA[1];
                            subBoundary[2] = new_subBoundary_RA[2];
                            subBoundary[3] = new_subBoundary_RA[3];
                            
                            SourceImages = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
                            
                            for(int image_index = 0 ; image_index < proinfo->number_of_images ; image_index++)
                            {
                                if(proinfo->check_selected_image[image_index])
                                {
                                    SourceImages[image_index] = SetsubsetImage(proinfo, levelinfo, image_index,param,NumOfIAparam,RPCs,t_Imageparams,subBoundary,minmaxHeight,Startpos_ori,Subsetsize);
                                    
                                    SetPySizes(data_size_lr[image_index], Subsetsize[image_index], pyramid_step);
                                }
                            }
                            
                            //new memory allocate
                            SubImages = (uint16***)malloc(sizeof(uint16**)*(pyramid_step+1));
                            SubOriImages = (uint8***)malloc(sizeof(uint8**)*(pyramid_step+1));
                            SubMagImages = (uint16***)malloc(sizeof(uint16**)*(pyramid_step+1));
                            
                            for(int iter_level = 0 ; iter_level < pyramid_step+1; iter_level++)
                            {
                                SubImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
                                SubOriImages[iter_level] = (uint8**)malloc(sizeof(uint8*)*proinfo->number_of_images);
                                SubMagImages[iter_level] = (uint16**)malloc(sizeof(uint16*)*proinfo->number_of_images);
                                
                                for(int image_index = 0 ; image_index < proinfo->number_of_images ; image_index++)
                                {
                                    if(proinfo->check_selected_image[image_index])
                                    {
                                        long int data_length = (long int)data_size_lr[image_index][iter_level].height*(long int)data_size_lr[image_index][iter_level].width;
                                        
                                        if(iter_level == 0)
                                            SubImages[iter_level][image_index] = SourceImages[image_index];
                                        else
                                            SubImages[iter_level][image_index] = NULL;//(uint16*)malloc(sizeof(uint16)*data_length);
                                        
                                        SubOriImages[iter_level][image_index] = (uint8*)malloc(sizeof(uint8)*data_length);
                                        SubMagImages[iter_level][image_index] = (uint16*)malloc(sizeof(uint16)*data_length);
                                    }
                                }
                            }
                            //pyramid image generation
                            SetPyramidImages(proinfo, pyramid_step+1, data_size_lr, SubImages, SubMagImages, SubOriImages);
                            
                            printf("Resize RA tile end\n");
                        }
                        
                        printf("subBoundary %f\t%f\t%f\t%f\t preBoundary %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3], preBoundary[0], preBoundary[1], preBoundary[2], preBoundary[3]);
                        
                        double Th_roh, Th_roh_min, Th_roh_start, Th_roh_next;
                        double minH_mps, maxH_mps;
                        double minH_grid, maxH_grid;
                        double MPP;
                        double MPP_simgle_image;
                        double MPP_stereo_angle = 1;
                        
                        if(level >= pyramid_step)
                            blunder_selected_level = level;
                        else if(level >= 2)
                            blunder_selected_level = level;
                        else
                            blunder_selected_level = level + 1;
                        
                        //blunder_selected_level = level;
                        printf("selected_bl %d\n",blunder_selected_level);
                        
                        uint8 iteration;
                        
                        D2DPOINT *Startpos = (D2DPOINT*)malloc(sizeof(D2DPOINT)*proinfo->number_of_images);
                        D2DPOINT *BStartpos= (D2DPOINT*)malloc(sizeof(D2DPOINT)*proinfo->number_of_images);
                        
                        D2DPOINT *Startpos_next = NULL;
                        if(level > Py_combined_level)
                        {
                            Startpos_next = (D2DPOINT*)malloc(sizeof(D2DPOINT)*proinfo->number_of_images);
                        }
                        
                        for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                        {
                            Startpos[ti].m_X       = (double)(Startpos_ori[ti].m_X/pwrtwo(level));
                            Startpos[ti].m_Y       = (double)(Startpos_ori[ti].m_Y/pwrtwo(level));
                            
                            BStartpos[ti].m_X       = (double)(Startpos_ori[ti].m_X/pwrtwo(blunder_selected_level));
                            BStartpos[ti].m_Y       = (double)(Startpos_ori[ti].m_Y/pwrtwo(blunder_selected_level));
                            
                            printf("Startpos %f\t%f\t%f\t%f\t%f\t%f\n",Startpos_ori[ti].m_X,Startpos_ori[ti].m_Y,Startpos[ti].m_X,Startpos[ti].m_Y,BStartpos[ti].m_X,BStartpos[ti].m_Y);
                            
                            if(level > Py_combined_level)
                            {
                                Startpos_next[ti].m_X       = (double)(Startpos_ori[ti].m_X/pwrtwo(level-1));
                                Startpos_next[ti].m_Y       = (double)(Startpos_ori[ti].m_Y/pwrtwo(level-1));
                            }
                        }
                        
                        levelinfo.py_Images = SubImages[level];
                        levelinfo.py_MagImages = SubMagImages[level];
                        levelinfo.py_OriImages = SubOriImages[level];
                        
                        levelinfo.py_BImages = SubImages[blunder_selected_level];
                        levelinfo.py_BMagImages = SubMagImages[blunder_selected_level];
                        
                        if(level > Py_combined_level)
                        {
                            levelinfo.py_Images_next = SubImages[level-1];
                            levelinfo.py_OriImages_next = SubOriImages[level-1];
                            levelinfo.py_MagImages_next = SubMagImages[level-1];
                        }
                        
                        levelinfo.py_Startpos = Startpos;
                        levelinfo.py_BStartpos = BStartpos;
                        levelinfo.py_Startpos_next = Startpos_next;
                        
                        levelinfo.RPCs = RPCs;
                        levelinfo.Boundary = subBoundary;
                        levelinfo.Pyramid_step = &level;
                        levelinfo.py_Sizes = data_size_lr;
                        levelinfo.Template_size = &Template_size;
                        levelinfo.ImageAdjust = t_Imageparams;
                        levelinfo.param = &param;
                        levelinfo.NumOfIAparam = &NumOfIAparam;
                        levelinfo.bin_angle = &bin_angle;
                        levelinfo.blunder_selected_level = &blunder_selected_level;
                        levelinfo.Py_combined_level = &Py_combined_level;
                        levelinfo.check_SGM = proinfo->SGM_py >= level;
                        
                        printf("levelinfo %d\t%d\n",*levelinfo.Pyramid_step,*levelinfo.blunder_selected_level);
                        
                        SetThs(proinfo,level,final_level_iteration, &Th_roh, &Th_roh_min, &Th_roh_next, &Th_roh_start);
                        
                        D2DPOINT *GridPT = NULL;
                        if(proinfo->IsRA)
                        {
                            py_resolution           = Image_res[0]*pwrtwo(pyramid_step+1);
                            grid_resolution         = Image_res[0]*pwrtwo(pyramid_step+1);
                            
                            GridPT                  = SetDEMGrid(subBoundary, grid_resolution, grid_resolution,&Size_Grid2D);
                            
                            printf("RA grid size %f\tSize_Grid2D %d\t%d\n", py_resolution, Size_Grid2D.width, Size_Grid2D.height);
                        }
                        else
                            GridPT  = SetGrids(proinfo, level, final_level_iteration, proinfo->resolution, &Size_Grid2D, proinfo->DEM_resolution, &py_resolution, &grid_resolution, subBoundary);
                        
                        const long int Grid_length = (long int)Size_Grid2D.width*(long int)Size_Grid2D.height;
                        
                        levelinfo.Size_Grid2D = &Size_Grid2D;
                        levelinfo.Grid_length = &Grid_length;
                        levelinfo.GridPts = GridPT;
                        levelinfo.grid_resolution = &grid_resolution;
                        levelinfo.Hinterval = &Hinterval;
                        
                        if(!flag_start)
                        {
                            printf("GridPT3 start\t seed flag %d\t filename %s\timage_resolution %f minmax %f %f\n", proinfo->pre_DEMtif, proinfo->priori_DEM_tif, Image_res[0], minmaxHeight[0], minmaxHeight[1]);
                            if (GridPT3)
                                free(GridPT3);
                            GridPT3 = SetGrid3PT(proinfo, levelinfo, Th_roh, minmaxHeight);
                        }
                        
                        if(flag_start)
                        {
                            if(proinfo->IsRA)
                            {
                                if(check_new_subBoundary_RA)
                                {
                                    GridPT3 = ResizeGirdPT3_RA(proinfo,levelinfo, pre_Size_Grid2D, Size_Grid2D, preBoundary,subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                                    
                                    check_new_subBoundary_RA = false;
                                    
                                    printf("start ResizeGridPT3 with newBoundry pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                }
                                else
                                {
                                    printf("start ResizeGridPT3 pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                    GridPT3 = ResizeGirdPT3(proinfo, levelinfo,pre_Size_Grid2D, Size_Grid2D, subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                                }
                            }
                            else
                            {
                                printf("start ResizeGridPT3 pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                GridPT3 = ResizeGirdPT3(proinfo, levelinfo, pre_Size_Grid2D, Size_Grid2D, subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                            }
                        }
                        
                        printf("end start ResizeGridPT3 minmax height %f\t%f\n",minmaxHeight[0],minmaxHeight[1]);
                        
                        pre_Size_Grid2D.width = Size_Grid2D.width;
                        pre_Size_Grid2D.height = Size_Grid2D.height;
                        pre_grid_resolution = grid_resolution;
                        
                        fprintf(fid,"level = %d, Completion of Gridinfo setup\t%d\t%d!!\n",level,Size_Grid2D.width,Size_Grid2D.height);
                        
                        fprintf(fid_header, "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", row, col, level, subBoundary[0], subBoundary[1], grid_resolution, Size_Grid2D.width,Size_Grid2D.height);
                        
                        D2DPOINT *Grid_wgs = ps2wgs(param,Grid_length,GridPT);
                        
                        levelinfo.Grid_wgs = Grid_wgs;
                        levelinfo.reference_id = 0;
                        
                        double height_step = GetHeightStep(level,Image_res[0]);
                        
                        if(proinfo->sensor_provider == PT)
                            height_step = GetHeightStep_Planet(level,Image_res[0]);
                        
                        levelinfo.height_step = &height_step;
                        
                        
                        if(!flag_start)
                            actual_pair(proinfo, levelinfo, minmaxHeight, compute_pair, pairinfo_return);
                    
                        
                        printf("done\n");
                        
                        for(int count = 0 ; count < levelinfo.pairinfo->NumberOfPairs ; count++)
                        {
                            printf("count %d\t%d\t%d\t%f\n",count, levelinfo.pairinfo->pairs[count].m_X,levelinfo.pairinfo->pairs[count].m_Y,levelinfo.pairinfo->BHratio[count]);
                        }
                        
                        printf("Height step %f\t%f\t%f\n",minmaxHeight[1],minmaxHeight[0],height_step);
                        
                        //exit(1);
                        
                        if(proinfo->pre_DEMtif && !flag_start)
                            VerticalLineLocus_seeddem(proinfo,levelinfo, GridPT3, minmaxHeight);
                        
                        iteration       = 1;
                        if(level == 0)
                            iteration = final_level_iteration;
                        
                        int pre_matched_pts=10;
                        double matching_change_rate = 100;
                        const double rate_th = 0.00999999;
                        const int max_iteration = 9;
                        float mem_th = 5;
                        if(level <= 1)
                            mem_th = 10;
                        
                        double minimum_memory;
                        const double level_total_memory =  CalMemorySize(proinfo,levelinfo,GridPT3,&minimum_memory,iteration,minmaxHeight);
                        
                        if(total_memory < level_total_memory)
                            total_memory = level_total_memory;
                        
                        printf("Memory : System %f\t SETSM required %f\tminimum %f\tcheck_matching_rate %d\n",proinfo->System_memory, total_memory,minimum_memory,check_matching_rate);
                        if(minimum_memory > proinfo->System_memory -2)
                        {
                            printf("System memory is not enough to run SETSM. Please reduce tilesize or assign more physical memory!!\n");
                            exit(1);
                        }
                        
                        //check_matching_rate = true;
                        
                        if(!check_matching_rate)
                        {
                            if(proinfo->IsRA || level <= 2 && (total_memory > proinfo->System_memory - mem_th  ))
                                check_matching_rate = true;
                        }
                        
                        if(proinfo->check_Matchtag && proinfo->DEM_resolution < 2)
                            check_matching_rate = true;
                        
                        //if(level == 0 && proinfo->DEM_resolution < 2)
                        //    check_matching_rate = true;
                        
                        printf("Creating GridVoxel with num pairs: %d\n", pairinfo_return.NumberOfPairs);
                        auto grid_voxel_size = check_matching_rate ? 0 : Size_Grid2D.width*Size_Grid2D.height;
                        GridVoxel grid_voxel = GridVoxel(grid_voxel_size, pairinfo_return.NumberOfPairs);
                        
                        if(proinfo->sensor_type == SB)
                        {
                            double sum_MPP_simgle_image = 0;
                            double sum_MPP_stereo_angle = 0;
                            for(int pair_number = 0 ; pair_number < levelinfo.pairinfo->NumberOfPairs ; pair_number++)
                            {
                                if(proinfo->IsRA)
                                    CalMPP_8(proinfo, levelinfo, minmaxHeight, CA, mean_product_res, &MPP_simgle_image, &MPP_stereo_angle, pair_number);
                                else
                                {
                                    if(proinfo->DEM_resolution <= 4)
                                        CalMPP(proinfo, levelinfo, minmaxHeight, CA, mean_product_res, &MPP_simgle_image, &MPP_stereo_angle, pair_number);
                                    else
                                        CalMPP_8(proinfo, levelinfo, minmaxHeight, CA, mean_product_res, &MPP_simgle_image, &MPP_stereo_angle, pair_number);
                                }
                                
                                if(sum_MPP_simgle_image < MPP_simgle_image)
                                    sum_MPP_simgle_image = MPP_simgle_image;
                                if(sum_MPP_stereo_angle < MPP_stereo_angle)
                                    sum_MPP_stereo_angle = MPP_stereo_angle;
                                
                                //sum_MPP_simgle_image += MPP_simgle_image;
                                //sum_MPP_stereo_angle += MPP_stereo_angle;
                            }
                            MPP_simgle_image = sum_MPP_simgle_image;
                            MPP_stereo_angle = sum_MPP_stereo_angle;
                            
                            //MPP_simgle_image = sum_MPP_simgle_image/levelinfo.pairinfo->NumberOfPairs;
                            //MPP_stereo_angle = sum_MPP_stereo_angle/levelinfo.pairinfo->NumberOfPairs;
                        }
                        else
                        {
                            MPP_simgle_image = proinfo->resolution*1.5;
                            MPP_stereo_angle = MPP_simgle_image;
                        }
                        
                        *stereo_angle_accuracy = MPP_stereo_angle;
                        
                        printf("final MPP %f\t%f\n",MPP_simgle_image,MPP_stereo_angle);
                        
                        NCCresult *nccresult = (NCCresult*)calloc(sizeof(NCCresult),Grid_length);
                        
                        
                        levelinfo.check_matching_rate = &check_matching_rate;
                        
                        bool level_check_matching_rate = false;
                        
                        if(levelinfo.pairinfo->NumberOfPairs < 2 || proinfo->IsRA)
                        {
                            while((Th_roh >= Th_roh_min || (matching_change_rate > rate_th)) )
                            {
                                levelinfo.ImageAdjust = t_Imageparams;
                                
                                levelinfo.iteration = &iteration;
                                
                                if(level == 0 &&  iteration == 3)
                                    matching_change_rate = 0.001;
                                
                                printf("%f \t %f\n",Th_roh,Th_roh_min);
                                
                                double Th_roh_update = 0;
                                
                                char filename_mps[500];
                                char filename_mps_asc[500];
                                char filename_mps_anchor[500];
                                char filename_mps_blunder[500];
                                char filename_mps_all[500];
                                
                                int count_results[2];
                                int count_results_anchor[2];
                                int count_MPs;
                                bool check_ortho_cal = false;
                                
                                uint8 ortho_level = 2;
                                if(proinfo->DEM_resolution >= 8)
                                    ortho_level = 3;
                                
                                if(level >= ortho_level)
                                {
                                    check_ortho_cal = true;
                                    
                                    if(level == 2 && iteration > 5)
                                        check_ortho_cal = false;
                                }
                                else
                                    check_ortho_cal = false;
                                
                                printf("ortho level = %d\n",ortho_level);
                                
                                fprintf(fid,"Starting computation of NCC\n iteration = %u\tTh_roh = %f\tTh_roh_start = %f\tGrid size %d %d\n",
                                        iteration, Th_roh,Th_roh_start,Size_Grid2D.width,Size_Grid2D.height);
                                
                                if(proinfo->IsRA)
                                {
                                    sprintf(filename_mps,"%s/txt/RA_matched_pts_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_asc,"%s/txt/RA_matched_pts_%d_%d_%d_%d_asc.txt",proinfo->save_filepath,row,col,level,iteration);
                                    
                                    sprintf(filename_mps_anchor,"%s/txt/RA_matched_pts_%d_%d_%d_%d_anchor.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_blunder,"%s/txt/RA_matched_pts_%d_%d_%d_%d_blunder.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_all,"%s/txt/RA_matched_pts_%d_%d_%d_%d_all.txt",proinfo->save_filepath,row,col,level,iteration);
                                }
                                else
                                {
                                    sprintf(filename_mps,"%s/txt/matched_pts_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_asc,"%s/txt/matched_pts_%d_%d_%d_%d_asc.txt",proinfo->save_filepath,row,col,level,iteration);
                                    
                                    sprintf(filename_mps_anchor,"%s/txt/matched_pts_%d_%d_%d_%d_anchor.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_blunder,"%s/txt/matched_pts_%d_%d_%d_%d_blunder.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_all,"%s/txt/matched_pts_%d_%d_%d_%d_all.txt",proinfo->save_filepath,row,col,level,iteration);
                                }
                                
                                printf("template size =%d\n",Template_size);
                                
                                if(!check_matching_rate)
                                    InitializeVoxel(proinfo,grid_voxel,levelinfo,GridPT3, nccresult,iteration,minmaxHeight);
                      
                                const long int Accessable_grid = VerticalLineLocus(grid_voxel,proinfo,image_info,nccresult,levelinfo,GridPT3,iteration,minmaxHeight);
                                
                                printf("Done VerticalLineLocus\tgrid %d\n",Accessable_grid);
                                
                                if(!check_matching_rate)
                                {
                                    if(level == 0 && iteration > 1 && proinfo->DEM_resolution >= 2)
                                    {
                                        fprintf(fid_header, "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", row, col, level, subBoundary[0], subBoundary[1], grid_resolution, Size_Grid2D.width,Size_Grid2D.height);
                                        
                                        printf("header file write %d\t%d\n",level,iteration);
                                    }
                                    
                                    const int MaxNumberofHeightVoxel = (int)((minmaxHeight[1] - minmaxHeight[0])/height_step);
                                    
                                    if(levelinfo.check_SGM)
                                        AWNCC_SGM(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,0);
                                    else
                                    {
                                        if(levelinfo.pairinfo->NumberOfPairs < 2)
                                            AWNCC_single(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,0);
                                        else
                                            AWNCC_multi(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight);
                                    }
                                    printf("Done AWNCC\n");
                                }
                                
                                printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd computation of NCC!! minmax %f %f\n",row,col,level,iteration,minmaxHeight[0], minmaxHeight[1]);
                                
                                minH_mps = 1000000;
                                maxH_mps = -1000000;
                                
                                if(level >= 5)
                                    MPP = MPP_simgle_image;
                                else if(MPP_stereo_angle > 5)
                                    MPP = MPP_stereo_angle;
                                else
                                    MPP = MPP_simgle_image;
                                
                                vector<D3DPOINT> MatchedPts_list;
                                
                                count_MPs = SelectMPs(proinfo, levelinfo, nccresult, GridPT3, Th_roh, Th_roh_min, Th_roh_start, Th_roh_next, iteration, MPP, final_level_iteration, MPP_stereo_angle, &MatchedPts_list);
                                
                                printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\t%d\n",row,col,level,iteration,count_MPs,MatchedPts_list.size());
                                
                                D3DPOINT *ptslists = NULL;
                                
                                vector<D3DPOINT>::iterator it;
                                vector<D3DPOINT> MatchedPts_list_mps;
                                vector<D3DPOINT> MatchedPts_list_blunder;
                                vector<D3DPOINT> MatchedPts_list_anchor;
                                
                                FILE *pMT_all = NULL;
                                if(count_MPs > 10)
                                {
                                    if (check_ortho_cal && proinfo->IsRA != 1)
                                    {
                                        printf("blunder detection for anchor points\n");
                                        //anchor points
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*MatchedPts_list.size());
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_all,"w");
                                        for(long count_pt = 0 ; count_pt < MatchedPts_list.size() ; count_pt ++)
                                        {
                                            ptslists[count_pt] = MatchedPts_list[count_pt];
                                            if(temp_asc_fprint)
                                                fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[count_pt].m_X,ptslists[count_pt].m_Y,ptslists[count_pt].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        DecisionMPs(proinfo, levelinfo, false,count_MPs,GridPT3, iteration, Hinterval,count_results_anchor, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_anchor,"w");
                                        long tcnt;
                                        count_results_anchor[0] = 0;
                                        for(tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_anchor.push_back(ptslists[tcnt]);
                                                count_results_anchor[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                            }
                                            ptslists[tcnt].flag = 0;
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd anchor points\n",row,col,level,iteration);
                                        
                                        printf("blunder detection for all points\n");
                                        //blunder detection
                                        DecisionMPs(proinfo, levelinfo, true,count_MPs,GridPT3, iteration, Hinterval,count_results, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_blunder,"w");
                                        count_results[0] = 0;
                                        for(tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_blunder.push_back(ptslists[tcnt]);
                                                count_results[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                                
                                                if(minH_mps > ptslists[tcnt].m_Z)
                                                    minH_mps        = ptslists[tcnt].m_Z;
                                                if(maxH_mps < ptslists[tcnt].m_Z)
                                                    maxH_mps       = ptslists[tcnt].m_Z;
                                            }
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        free(ptslists);
                                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
                                    }
                                    else
                                    {
                                        printf("blunder detection for all points\n");
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_all,"w");
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*MatchedPts_list.size());
                                        for(long count_pt = 0 ; count_pt < MatchedPts_list.size() ; count_pt ++)
                                        {
                                            ptslists[count_pt] = MatchedPts_list[count_pt];
                                            if(temp_asc_fprint)
                                                fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[count_pt].m_X,ptslists[count_pt].m_Y,ptslists[count_pt].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        DecisionMPs(proinfo, levelinfo, true,count_MPs,GridPT3, iteration, Hinterval,count_results, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_blunder,"w");
                                        
                                        count_results[0] = 0;
                                        for(int tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_mps.push_back(ptslists[tcnt]);
                                                count_results[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                                
                                                if(minH_mps > ptslists[tcnt].m_Z)
                                                    minH_mps        = ptslists[tcnt].m_Z;
                                                if(maxH_mps < ptslists[tcnt].m_Z)
                                                    maxH_mps       = ptslists[tcnt].m_Z;
                                            }
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        count_MPs       = count_results[0];
                                        
                                        free(ptslists);
                                        printf("RA row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
                                    }
                                }
                                else
                                    lower_level_match = false;
                                
                                MatchedPts_list.clear();
                                vector<D3DPOINT>().swap(MatchedPts_list);
                                
                                fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tcheck = %d(%d)\tEnd blunder detection\n",row,col,level,iteration,lower_level_match,count_MPs);
                                
                                printf("End computation of blunder!! Mps = %d\tTris = %d\tminz Mp = %f\tmaxz Mp = %f minmax %f %f \n",
                                       count_results[0],count_results[1],minH_mps,maxH_mps,minmaxHeight[0],minmaxHeight[1]);
                                
                                if(lower_level_match)
                                {
                                    if(check_ortho_cal && proinfo->IsRA != 1)
                                    {
                                        printf("settingflag %d\t%d\n",MatchedPts_list_anchor.size(),MatchedPts_list_blunder.size());
                                        count_MPs = SetttingFlagOfGrid(levelinfo, GridPT3, MatchedPts_list_anchor, MatchedPts_list_blunder, &MatchedPts_list_mps);
                                        MatchedPts_list_anchor.clear();
                                        vector<D3DPOINT>().swap(MatchedPts_list_anchor);
                                        MatchedPts_list_blunder.clear();
                                        vector<D3DPOINT>().swap(MatchedPts_list_blunder);
                                    }
                                    
                                    printf("count_MPs %d\t%d\n",count_MPs,MatchedPts_list_mps.size());
                                    count_MPs = MatchedPts_list_mps.size();
                                    
                                    vector<UI3DPOINT> t_trilists;
                                    vector<UI3DPOINT>::iterator it_tri;
                                    
                                    if(level == 0 && iteration == 3)
                                    {
                                        D3DPOINTSAVE *ptslists_save = (D3DPOINTSAVE*)calloc(count_MPs,sizeof(D3DPOINTSAVE));
                                        double minmaxBR[6] = {10000000, 10000000, -10000000, -10000000, 100000, -100000};
                                        
                                        int i = 0;
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_asc,"w");
                                        for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                        {
                                            ptslists_save[i].m_X = MatchedPts_list_mps[i].m_X;
                                            ptslists_save[i].m_Y = MatchedPts_list_mps[i].m_Y;
                                            ptslists_save[i].m_Z = MatchedPts_list_mps[i].m_Z;
                                            
                                            if(minmaxBR[0] > ptslists_save[i].m_X)
                                                minmaxBR[0]     = ptslists_save[i].m_X;
                                            if(minmaxBR[1] > ptslists_save[i].m_Y)
                                                minmaxBR[1]     = ptslists_save[i].m_Y;
                                            
                                            if(minmaxBR[2] < ptslists_save[i].m_X)
                                                minmaxBR[2]     = ptslists_save[i].m_X;
                                            if(minmaxBR[3] < ptslists_save[i].m_Y)
                                                minmaxBR[3]     = ptslists_save[i].m_Y;
                                            if(minmaxBR[4] > ptslists_save[i].m_Z)
                                                minmaxBR[4] = ptslists_save[i].m_Z;
                                            if(minmaxBR[5] < ptslists_save[i].m_Z)
                                                minmaxBR[5] = ptslists_save[i].m_Z;
                                            
                                            if(temp_asc_fprint)
                                            fprintf(pMT_all,"%f\t%f\t%f\n",ptslists_save[i].m_X,ptslists_save[i].m_Y,ptslists_save[i].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        MatchedPts_list_mps.clear();
                                        vector<D3DPOINT>().swap(MatchedPts_list_mps);
                                        
                                        FILE *pFile = fopen(filename_mps,"wb");
                                        fwrite(ptslists_save,sizeof(D3DPOINTSAVE),count_MPs,pFile);
                                        fclose(pFile);
                                        
                                        if(!proinfo->IsRA)
                                        {
                                            FILE *fid_BR;
                                            FILE *fid_count;
                                            
                                            char save_file[500];
                                            sprintf(save_file,"%s/txt/matched_BR_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                                            fid_BR  = fopen(save_file,"w");
                                            fprintf(fid_BR,"%f\t%f\t%f\t%f\t%f\t%f\n",minmaxBR[0],minmaxBR[1],minmaxBR[2],minmaxBR[3],minmaxBR[4],minmaxBR[5]);
                                            fclose(fid_BR);
                                            
                                            sprintf(save_file,"%s/txt/count_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                                            fid_count  = fopen(save_file,"w");
                                            
                                            fprintf(fid_count,"%d\n",count_MPs);
                                            fclose(fid_count);
                                            
                                            if(temp_asc_fprint)
                                            {
                                                echoprint_Gridinfo_asc(proinfo,levelinfo, nccresult,row,col,level,iteration,Size_Grid2D,GridPT3);
                                            }
                                            
                                            echoprint_Gridinfo(proinfo,nccresult,row,col,level,iteration,0,&Size_Grid2D,GridPT3,(char*)"final");
                                        }
                                        free(ptslists_save);
                                        
                                        matching_change_rate = 0.001;
                                    }
                                    else
                                    {
                                        if(check_ortho_cal && proinfo->IsRA != 1)
                                        {
                                            double maxX_ptslists = -100000000;
                                            double maxY_ptslists = -100000000;
                                            double minX_ptslists =  100000000;
                                            double minY_ptslists =  100000000;
                                            
                                            int i;
                                            
                                            ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                            
                                            for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                            {
                                                ptslists[i] = MatchedPts_list_mps[i];
                                                //if(level >= 3)
                                                //    ptslists[i].flag = 1; //temporary blunders flag for ortho blunder
                                            }
                                            
                                            MatchedPts_list_mps.clear();
                                            vector<D3DPOINT>().swap(MatchedPts_list_mps);
                           
                                            double min_max[4] = {subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]};
                                            UI3DPOINT *trilists;
                                            
                                            FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                            delete origTri;
                                            
                                            count_tri = t_trilists.size();
                                            trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                            
                                            for(i = 0 ; i < t_trilists.size() ; i++)
                                                trilists[i] = t_trilists[i];
                                            
                                            t_trilists.clear();
                                            vector<UI3DPOINT>().swap(t_trilists);
                                            
                                            fprintf(fid,"level = %d\tMatching Pts = %d\n",level,count_results[0]);
                                            
                                            printf("ortho minmax %f %f pts anchor blunder %d %d \n",minmaxHeight[0],minmaxHeight[1],count_MPs,count_tri);
                                            
                                            count_results[0] = Ortho_blunder(proinfo, levelinfo, MPP_simgle_image, ptslists, count_MPs, trilists,count_tri, GridPT3);
                                            free(trilists);
                                            
                                            printf("end ortho_blunder %d\n",count_results[0]);
                                            
                                            if(level >= 3)
                                            {
                                                int matched_pts = 0;
                                                vector<D3DPOINT> ortho_list;
                                                for(i=0;i<count_MPs;i++)
                                                {
                                                    if(ptslists[i].flag != 1)
                                                    {
                                                        ortho_list.push_back(ptslists[i]);
                                                        matched_pts++;
                                                    }
                                                }
                                                
                                                free(ptslists);
                                                count_MPs = matched_pts;
                                                ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                                
                                                for( i = 0 ; i < ortho_list.size() ; i++)
                                                {
                                                    ptslists[i] = ortho_list[i];
                                                }
                                                ortho_list.clear();
                                                vector<D3DPOINT>().swap(ortho_list);
                                            }
                                            
                                            printf("load ortho_blunder pts %d\n",count_MPs);
                                            
                                            //Save triangulation and delete it since we will not use it
                                            FullTriangulation *origTri_2 = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                            delete origTri_2;
                                            
                                            count_tri = t_trilists.size();
                                            trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                            
                                            for(i = 0 ; i < t_trilists.size() ; i++)
                                                trilists[i] = t_trilists[i];
                                            t_trilists.clear();
                                            vector<UI3DPOINT>().swap(t_trilists);
                              
                                            DecisionMPs_setheight(proinfo, levelinfo, count_MPs,GridPT3, iteration, Hinterval, minmaxHeight, ptslists,trilists,count_tri);
                                            
                                            printf("end decision_setheight\n");
                                            
                                            if(pre_matched_pts == 0)
                                                matching_change_rate = 0;
                                            else
                                                matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
                                            
                                            matching_rate = count_MPs/(double)(Accessable_grid);
                                            
                                            printf("matching change rate pre curr %f\t%d\t%d\tmatching rate %f\t%d\n",matching_change_rate,count_MPs,pre_matched_pts,matching_rate,Accessable_grid);
                                            pre_matched_pts = count_MPs;
                                            
                                            //if(level <= 2 && matching_rate < th_mr && proinfo->DEM_resolution <= 4)
                                            //    level_check_matching_rate = true;
                                            
                                            if(iteration > max_iteration)
                                                matching_change_rate = 0.001;
                                            
                                            if(level == 0)
                                            {
                                                if(proinfo->DEM_resolution < 2)
                                                    matching_change_rate = 0.001;
                                                if(iteration > 2)
                                                    matching_change_rate = 0.001;
                                            }
                                            
                                            if(level <= 1)
                                            {
                                                if(iteration >= ceil(max_iteration/2.0))
                                                    matching_change_rate = 0.001;
                                            }
                                            
                                            if(proinfo->IsRA)
                                                matching_change_rate = 0.001;
                                            
                                            if(proinfo->DEM_resolution >= 8)
                                                matching_change_rate = 0.001;
                                 
                                            if(Th_roh >= Th_roh_min)
                                            {
                                                if(level == 0)
                                                {
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                    if(iteration >= 2)
                                                        Th_roh_update          = (double)(Th_roh - 0.50);
                                                }
                                                else if(level == 1)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else if(level == 2)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else if(level == 3)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else
                                                {
                                                    if(proinfo->IsRA)
                                                        Th_roh_update       = (double)(Th_roh - 0.10);
                                                    else
                                                        Th_roh_update       = (double)(Th_roh - 0.06);
                                                }
                                            }
                                            
                                            printf("matching change rate pre curr %f\t%d\t%d\tTh_roh %f\t%f\n",matching_change_rate,count_MPs,pre_matched_pts,Th_roh,Th_roh_min);
                                            
                                            if(temp_asc_fprint)
                                                echoprint_Gridinfo_asc(proinfo, levelinfo, nccresult,row,col,level,iteration,Size_Grid2D,GridPT3);
                                            
                                            if(level == 0)
                                            {
                                                if(MPP_stereo_angle > 5)
                                                    MPP = MPP_stereo_angle;
                                                else
                                                    MPP = MPP_simgle_image;
                                                
                                                if(proinfo->DEM_resolution < 2)
                                                    Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                                else
                                                    GridPT3         = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
                                            {
                                                if(MPP_stereo_angle > 5)
                                                {
                                                    if(level > 2)
                                                        MPP = MPP_stereo_angle;
                                                    else
                                                        MPP = MPP_simgle_image;
                                                }
                                                else
                                                    MPP = MPP_simgle_image;
                                                
                                                Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            else
                                            {
                                                if(MPP_stereo_angle > 5)
                                                {
                                                    if(level <= 1)
                                                        MPP = MPP_stereo_angle;
                                                    else
                                                        MPP = MPP_simgle_image;
                                                }
                                                else
                                                    MPP = MPP_simgle_image;
                                                
                                                GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            
                                            if(temp_asc_fprint)
                                            {
                                                pMT_all = fopen(filename_mps_asc,"w");
                                                for( i = 0 ; i < count_MPs ; i++)
                                                {
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[i].m_X,ptslists[i].m_Y,ptslists[i].m_Z);
                                                }
                                                fclose(pMT_all);
                                            }
                                            
                                            free(trilists);
                                            free(ptslists);
                                            
                                            final_iteration = iteration;
                                        }
                                        else
                                        {
                                            double maxX_ptslists = -100000000;
                                            double maxY_ptslists = -100000000;
                                            double minX_ptslists =  100000000;
                                            double minY_ptslists =  100000000;
                                            int i;
                                            
                                            ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                            
                                            for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                                ptslists[i] = MatchedPts_list_mps[i];
                                    
                                            MatchedPts_list_mps.clear();
                                            vector<D3DPOINT>().swap(MatchedPts_list_mps);
                                            
                                            UI3DPOINT *trilists;
                                            
                                            double min_max[4] = {subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]};
                                            
                                            //Save triangulation and delete it since we will not use it
                                            FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                            delete origTri;
                                            
                                            count_tri = t_trilists.size();
                                            trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                            
                                            for(i = 0 ; i < t_trilists.size() ; i++)
                                                trilists[i] = t_trilists[i];
                                 
                                            t_trilists.clear();
                                            vector<UI3DPOINT>().swap(t_trilists);
                                            
                                            if(pre_matched_pts == 0)
                                                matching_change_rate = 0;
                                            else
                                                matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
                                            
                                            printf("matching change rate pre curr %f\t%d\t%d\n",matching_change_rate,count_MPs,pre_matched_pts);
                                            pre_matched_pts = count_results[0];
                                            
                                            if(iteration > max_iteration)
                                                matching_change_rate = 0.001;
                                            
                                            if(level == 0)
                                            {
                                                if(proinfo->DEM_resolution < 2)
                                                    matching_change_rate = 0.001;
                                                if(iteration > 2)
                                                    matching_change_rate = 0.001;
                                            }
                                            
                                            if(level <= 1)
                                            {
                                                if(iteration >= ceil(max_iteration/2.0))
                                                    matching_change_rate = 0.001;
                                            }
                                            
                                            if(proinfo->IsRA)
                                                matching_change_rate = 0.001;
                                            
                                            if(proinfo->DEM_resolution >= 8)
                                                matching_change_rate = 0.001;
                                   
                                            if(Th_roh >= Th_roh_min)
                                            {
                                                if(level == 0)
                                                {
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                    if(iteration >= 2)
                                                        Th_roh_update          = (double)(Th_roh - 0.50);
                                                }
                                                else if(level == 1)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else if(level == 2)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else if(level == 3)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else
                                                {
                                                    if(proinfo->IsRA)
                                                        Th_roh_update       = (double)(Th_roh - 0.10);
                                                    else
                                                        Th_roh_update       = (double)(Th_roh - 0.06);
                                                }
                                            }
                                            
                                            matching_rate = count_MPs/(double)(Accessable_grid);
                                            printf("matching change rate pre curr %f\t%d\t%d\tTh_roh %f\t%f\tmatching rate %f\t%d\n",matching_change_rate,count_MPs,pre_matched_pts,Th_roh,Th_roh_min,matching_rate,Accessable_grid);
                                            
                                            //if(level <= 2 && matching_rate < th_mr && proinfo->DEM_resolution <= 4)
                                            //    level_check_matching_rate = true;
                                            
                                            bool check_level_end = false;
                                            
                                            if(level != 0)
                                            {
                                                if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th)
                                                    check_level_end = true;
                                            }
                                            
                                            if(temp_asc_fprint)
                                                echoprint_Gridinfo_asc(proinfo, levelinfo, nccresult,row,col,level,iteration,Size_Grid2D,GridPT3);
                                            
                                            if(level == 0)
                                            {
                                                if(MPP_stereo_angle > 5)
                                                    MPP = MPP_stereo_angle;
                                                else
                                                    MPP = MPP_simgle_image;
                                                if(proinfo->DEM_resolution < 2)
                                                {
                                                    Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                                }
                                                else
                                                {
                                                    GridPT3         = SetHeightRange(proinfo,levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                                }
                                            }
                                            else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
                                            {
                                                if(MPP_stereo_angle > 5)
                                                {
                                                    if(level > 2)
                                                        MPP = MPP_stereo_angle;
                                                    else
                                                        MPP = MPP_simgle_image;
                                                }
                                                else
                                                    MPP = MPP_simgle_image;
                                                
                                                Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            else
                                            {
                                                if(MPP_stereo_angle > 5)
                                                {
                                                    if(level <= 1)
                                                        MPP = MPP_stereo_angle;
                                                    else
                                                        MPP = MPP_simgle_image;
                                                }
                                                else
                                                    MPP = MPP_simgle_image;
                                                
                                                GridPT3     = SetHeightRange(proinfo,levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            
                                            free(trilists);
                                            
                                            if(proinfo->IsRA && level <= 3)
                                            {
                                                int RA_iter_counts = 0;
                                                RA_iter_counts = AdjustParam(proinfo, levelinfo, count_MPs, t_Imageparams, pyramid_step, ptslists);
                                                for(int ti = 0 ; ti < levelinfo.pairinfo->NumberOfPairs ; ti++)
                                                {
                                                    //if(proinfo->check_selected_image[ti])
                                                    {
                                                        fprintf(fid,"RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Imageparams[ti][0],t_Imageparams[ti][1]);
                                                        printf("RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Imageparams[ti][0],t_Imageparams[ti][1]);
                                                    }
                                                }
                                                
                                                if (level <= 1)
                                                {
                                                    char save_file[500];
                                                    sprintf(save_file,"%s/txt/RAinfo.txt",proinfo->save_filepath);
                                                    FILE *fid_RAinfo  = fopen(save_file,"w");
                                                    for(int ti = 0 ; ti < levelinfo.pairinfo->NumberOfPairs ; ti++)
                                                        fprintf(fid_RAinfo,"%d\t%d\t%f\t%f\t%f\n",levelinfo.pairinfo->pairs[ti].m_X,levelinfo.pairinfo->pairs[ti].m_Y,t_Imageparams[ti][0],t_Imageparams[ti][1],levelinfo.pairinfo->BHratio[ti]);
                                                    fclose(fid_RAinfo);
                                                }
                                            }
                                            
                                            if(proinfo->IsRA)
                                            {
                                                if(check_RA_divide)
                                                {
                                                    if(level <= 3 && iteration > 2)
                                                    {
                                                        int* t_count = (int*)calloc(sizeof(int),total_tile);
                                                        
                                                        for(int k=0;k<count_MPs;k++)
                                                        {
                                                            int t_col = floor((ptslists[k].m_X - subBoundary[0])/(double)tilesize_RA);
                                                            int t_row = floor((ptslists[k].m_Y - subBoundary[1])/(double)tilesize_RA);
                                                            
                                                            t_count[t_col+division_X*t_row] ++;
                                                        }
                                                        
                                                        int saved_count = 0;
                                                        int selected_X = 0;
                                                        int selected_Y = 0;
                                                        int total_count = 0;
                                                        for(int k=0;k<total_tile;k++)
                                                        {
                                                            if(t_count[k] > saved_count)
                                                            {
                                                                saved_count = t_count[k];
                                                                selected_Y  = floor(k/division_X);
                                                                selected_X  = k % division_X;
                                                            }
                                                            total_count += t_count[k];
                                                        }
                                                        printf("total_count %d\tsaved_count %d\tselected_X %d\tselected_Y %d\n",total_count,saved_count,selected_X,selected_Y);
                                                        printf("selected br %f\t%f\t%f\t%f\n",subBoundary[0] + selected_X*tilesize_RA,subBoundary[1] + selected_Y*tilesize_RA,
                                                               subBoundary[0] + (selected_X+1)*tilesize_RA,subBoundary[1] + (selected_Y+1)*tilesize_RA);
                                                        
                                                        fprintf(fid,"total_count %d\tsaved_count %d\tselected_X %d\tselected_Y %d\n",total_count,saved_count,selected_X,selected_Y);
                                                        
                                                        new_subBoundary_RA[0] = subBoundary[0] + selected_X*tilesize_RA;
                                                        new_subBoundary_RA[1] = subBoundary[1] + selected_Y*tilesize_RA;
                                                        new_subBoundary_RA[2] = subBoundary[0] + (selected_X+1)*tilesize_RA;
                                                        new_subBoundary_RA[3] = subBoundary[1] + (selected_Y+1)*tilesize_RA;
                                                        
                                                        check_new_subBoundary_RA = true;
                                                        check_RA_divide = false;
                                                        
                                                        free(t_count);
                                                    }
                                                }
                                            }
                                            
                                            if(temp_asc_fprint)
                                            {
                                                pMT_all = fopen(filename_mps_asc,"w");
                                                for( i = 0 ; i < count_MPs ; i++)
                                                {
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[i].m_X,ptslists[i].m_Y,ptslists[i].m_Z);
                                                }
                                                fclose(pMT_all);
                                            }
                                            
                                            free(ptslists);
                                        }
                                    }
                                    
                                    fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\tmatching_rate = %f\n",
                                            row,col,level,iteration,count_MPs,minmaxHeight[0],minmaxHeight[1],matching_rate);
                                    printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\n",
                                           row,col,level,iteration,count_MPs,minmaxHeight[0],minmaxHeight[1]);
                                }
                                
                                if(lower_level_match)
                                {
                                    flag_start          = true;
                                    iteration++;
                                }
                                else if(level == 0)
                                    iteration++;
                                
                                if(level == 0)
                                {
                                    if(proinfo->DEM_resolution >= 2)
                                        Th_roh          = (double)(Th_roh - 0.10);
                                    else
                                        Th_roh          = (double)(Th_roh - 0.50);
                                    
                                    if(iteration > 3)
                                        Th_roh          = (double)(Th_roh - 0.50);
                                }
                                else if(level == 1)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else if(level == 2)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else if(level == 3)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else
                                {
                                    if(proinfo->IsRA)
                                        Th_roh          = (double)(Th_roh - 0.10);
                                    else
                                        Th_roh          = (double)(Th_roh - 0.06);
                                }
                                
                                if(lower_level_match)
                                {
                                    if(Th_roh < Th_roh_min && matching_change_rate > rate_th)
                                    {
                                        if(level == 0)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 1)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 2)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 3)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else
                                        {
                                            if(proinfo->IsRA)
                                                Th_roh          = (double)(Th_roh + 0.10);
                                            else
                                                Th_roh          = (double)(Th_roh + 0.06);
                                        }
                                    }
                                }
                                
                                if (!lower_level_match && Th_roh < Th_roh_min)
                                {
                                    if(level > 0)
                                    {
                                        iteration++;
                                        matching_change_rate = 0.001;
                                        Th_roh_min = 0.4;
                                    }
                                }
                                
                                if(level == 0)
                                    final_level_iteration = iteration;
                                
                                printf("Memory : System %f\t SETSM required %f\n",proinfo->System_memory, total_memory);
                            }
                        }
                        else if(!check_matching_rate)
                        {
                            while((Th_roh >= Th_roh_min || (matching_change_rate > rate_th)) )
                            {
                                levelinfo.ImageAdjust = t_Imageparams;
                                levelinfo.iteration = &iteration;
                                
                                if(level == 0 &&  iteration == 3)
                                    matching_change_rate = 0.001;
                                
                                printf("%f \t %f\n",Th_roh,Th_roh_min);
                                
                                double Th_roh_update = 0;
                                
                                char filename_mps[500];
                                char filename_mps_asc[500];
                                char filename_mps_anchor[500];
                                char filename_mps_blunder[500];
                                char filename_mps_all[500];
                                
                                int count_results[2];
                                int count_results_anchor[2];
                                int count_MPs;
                                bool check_ortho_cal = false;
                                
                                uint8 ortho_level = 2;
                                if(proinfo->DEM_resolution >= 8)
                                    ortho_level = 3;
                                
                                if(level >= ortho_level)
                                {
                                    check_ortho_cal = true;
                                    
                                    if(level == 2 && iteration > 5)
                                        check_ortho_cal = false;
                                }
                                else
                                    check_ortho_cal = false;
                                
                                printf("ortho level = %d\n",ortho_level);
                                
                                fprintf(fid,"Starting computation of NCC\n iteration = %u\tTh_roh = %f\tTh_roh_start = %f\tGrid size %d %d\n",
                                        iteration, Th_roh,Th_roh_start,Size_Grid2D.width,Size_Grid2D.height);
                                
                                minH_mps = 1000000;
                                maxH_mps = -1000000;
                                
                                if(level >= 5)
                                    MPP = MPP_simgle_image;
                                else if(MPP_stereo_angle > 5)
                                    MPP = MPP_stereo_angle;
                                else
                                    MPP = MPP_simgle_image;
                                
                                vector<D3DPOINT> MatchedPts_list;
                                
                                if(proinfo->IsRA)
                                {
                                    sprintf(filename_mps,"%s/txt/RA_matched_pts_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_asc,"%s/txt/RA_matched_pts_%d_%d_%d_%d_asc.txt",proinfo->save_filepath,row,col,level,iteration);
                                    
                                    sprintf(filename_mps_anchor,"%s/txt/RA_matched_pts_%d_%d_%d_%d_anchor.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_blunder,"%s/txt/RA_matched_pts_%d_%d_%d_%d_blunder.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_all,"%s/txt/RA_matched_pts_%d_%d_%d_%d_all.txt",proinfo->save_filepath,row,col,level,iteration);
                                }
                                else
                                {
                                    sprintf(filename_mps,"%s/txt/matched_pts_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_asc,"%s/txt/matched_pts_%d_%d_%d_%d_asc.txt",proinfo->save_filepath,row,col,level,iteration);
                                    
                                    sprintf(filename_mps_anchor,"%s/txt/matched_pts_%d_%d_%d_%d_anchor.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_blunder,"%s/txt/matched_pts_%d_%d_%d_%d_blunder.txt",proinfo->save_filepath,row,col,level,iteration);
                                    sprintf(filename_mps_all,"%s/txt/matched_pts_%d_%d_%d_%d_all.txt",proinfo->save_filepath,row,col,level,iteration);
                                }
                                
                                printf("template size =%d\n",Template_size);
                                
                                if(!check_matching_rate)
                                    InitializeVoxel(proinfo,grid_voxel,levelinfo,GridPT3, nccresult,iteration,minmaxHeight);
                                
                                const long int Accessable_grid = VerticalLineLocus(grid_voxel,proinfo,image_info,nccresult,levelinfo,GridPT3,iteration,minmaxHeight);
                                
                                printf("Done VerticalLineLocus\tgrid %d\n",Accessable_grid);
                                
                                vector<D3DPOINT>::iterator it;
                                vector<D3DPOINT> MatchedPts_list_mps;
                                vector<D3DPOINT> MatchedPts_list_blunder;
                                vector<D3DPOINT> MatchedPts_list_anchor;
                                
                                D3DPOINT *ptslists = NULL;
                                FILE *pMT_all = NULL;
                                
                                vector<UI3DPOINT> t_trilists;
                                vector<UI3DPOINT>::iterator it_tri;
                                
                                UI3DPOINT *trilists;
                                
                                double maxX_ptslists = -100000000;
                                double maxY_ptslists = -100000000;
                                double minX_ptslists =  100000000;
                                double minY_ptslists =  100000000;
                                
                                double min_max[4] = {subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]};
                                
                                if(level == 0 && iteration > 1 && proinfo->DEM_resolution >= 2)
                                {
                                    fprintf(fid_header, "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", row, col, level, subBoundary[0], subBoundary[1], grid_resolution, Size_Grid2D.width,Size_Grid2D.height);
                                    
                                    printf("header file write %d\t%d\n",level,iteration);
                                }
                                
                                const int MaxNumberofHeightVoxel = (int)((minmaxHeight[1] - minmaxHeight[0])/height_step);
                                
                                lower_level_match = false;
                                MultiMPs **multimps = (MultiMPs**)calloc(sizeof(MultiMPs*),Grid_length);
                                for(long count = 0 ; count < Grid_length ; count++)
                                {
                                    multimps[count] = (MultiMPs*)calloc(sizeof(MultiMPs),levelinfo.pairinfo->NumberOfPairs+1);
                                    for(int temp_i = 0 ; temp_i < levelinfo.pairinfo->NumberOfPairs+1 ; temp_i++)
                                    {
                                        multimps[count][temp_i].peak_height = Nodata;
                                        multimps[count][temp_i].check_matched = false;
                                    }
                                }
                                
                                long max_count_MPs = 0;
                                long total_count_MPs = 0 ;
                                float min_pair_H = 9999;
                                float max_pair_H = -9999;
                                //if(!levelinfo.check_SGM) //no SGM
                                {
                                    AWNCC_AWNCC(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight);
                                    
                                    count_MPs = SelectMPs(proinfo, levelinfo, nccresult, GridPT3, Th_roh, Th_roh_min, Th_roh_start, Th_roh_next, iteration, MPP, final_level_iteration, MPP_stereo_angle, &MatchedPts_list);
                                    
                                    for(long tcnt=0;tcnt<MatchedPts_list.size();tcnt++)
                                    {
                                        long t_col         = (long)((MatchedPts_list[tcnt].m_X - subBoundary[0])/grid_resolution + 0.5);
                                        long t_row         = (long)((MatchedPts_list[tcnt].m_Y - subBoundary[1])/grid_resolution + 0.5);
                                        const long ref_index((long)Size_Grid2D.width*t_row + t_col);
                                        
                                        if(t_col >= 0 && t_col < Size_Grid2D.width && t_row >=0 && t_row < Size_Grid2D.height)
                                        {
                                            multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_roh = MatchedPts_list[tcnt].m_roh;
                                            multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_height = MatchedPts_list[tcnt].m_Z;
                                            multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].check_matched = true;
                                            total_count_MPs++;
                                            //printf("pos %d\t%d\t%d\t%f\n",t_col,t_row,multimps[ref_index][pair_number].peak_roh,multimps[ref_index][pair_number].peak_height);
                                            
                                            //fprintf(pfile_pair,"%f\t%f\t%f\n",MatchedPts_list[tcnt].m_X,MatchedPts_list[tcnt].m_Y,MatchedPts_list[tcnt].m_Z);
                                            if(min_pair_H > multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_height)
                                                min_pair_H = multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_height;
                                            if(max_pair_H < multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_height)
                                                max_pair_H = multimps[ref_index][levelinfo.pairinfo->NumberOfPairs].peak_height;
                                            
                                        }
                                    }
                                    
                                    printf("AWNCC row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\t%d\t%f\t%f\n",row,col,level,iteration,count_MPs,MatchedPts_list.size(),min_pair_H,max_pair_H);
                                    
                                    MatchedPts_list.clear();
                                }
                                
                                char filename_pair[500];
                                FILE *pfile_pair = NULL;
                                for(int pair_number = 0 ; pair_number < levelinfo.pairinfo->NumberOfPairs ; pair_number++)
                                {
                                    min_pair_H = 9999;
                                    max_pair_H = -9999;
                                    
                                    if(levelinfo.check_SGM) //SGM
                                        AWNCC_SGM(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,pair_number);
                                    else
                                        AWNCC_single(proinfo,grid_voxel,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,pair_number);
                                    
                                    printf("pairnumber = %d\trow = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd computation of NCC!! minmax %f %f\n",pair_number,row,col,level,iteration,minmaxHeight[0], minmaxHeight[1]);
                                    
                                    
                                    count_MPs = SelectMPs(proinfo, levelinfo, nccresult, GridPT3, Th_roh, Th_roh_min, Th_roh_start, Th_roh_next, iteration, MPP, final_level_iteration, MPP_stereo_angle, &MatchedPts_list);
                                    
                                    if(max_count_MPs < count_MPs)
                                        max_count_MPs = count_MPs;
                                    
                                    if(temp_asc_fprint)
                                    {
                                        sprintf(filename_pair,"%s/txt/matched_pts_%d_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration,pair_number);
                                        pfile_pair = fopen(filename_pair,"w");
                                    }
                                    for(long tcnt=0;tcnt<MatchedPts_list.size();tcnt++)
                                    {
                                        long t_col         = (long)((MatchedPts_list[tcnt].m_X - subBoundary[0])/grid_resolution + 0.5);
                                        long t_row         = (long)((MatchedPts_list[tcnt].m_Y - subBoundary[1])/grid_resolution + 0.5);
                                        const long ref_index((long)Size_Grid2D.width*t_row + t_col);
                                        
                                        if(t_col >= 0 && t_col < Size_Grid2D.width && t_row >=0 && t_row < Size_Grid2D.height)
                                        {
                                            multimps[ref_index][pair_number].peak_roh = MatchedPts_list[tcnt].m_roh;
                                            multimps[ref_index][pair_number].peak_height = MatchedPts_list[tcnt].m_Z;
                                            multimps[ref_index][pair_number].check_matched = true;
                                            total_count_MPs++;
                                            //printf("pos %d\t%d\t%d\t%f\n",t_col,t_row,multimps[ref_index][pair_number].peak_roh,multimps[ref_index][pair_number].peak_height);
                                            if(temp_asc_fprint)
                                                fprintf(pfile_pair,"%f\t%f\t%f\n",MatchedPts_list[tcnt].m_X,MatchedPts_list[tcnt].m_Y,MatchedPts_list[tcnt].m_Z);
                                            
                                            if(min_pair_H > multimps[ref_index][pair_number].peak_height)
                                                min_pair_H = multimps[ref_index][pair_number].peak_height;
                                            if(max_pair_H < multimps[ref_index][pair_number].peak_height)
                                                max_pair_H = multimps[ref_index][pair_number].peak_height;
                                        }
                                    }
                                    
                                    printf("pair %d\trow = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\t%d\t%f\t%f\n",pair_number,row,col,level,iteration,count_MPs,MatchedPts_list.size(),min_pair_H,max_pair_H);
                                    
                                    //if(count_MPs < 1)
                                    //    exit(1);
                                    
                                    MatchedPts_list.clear();
                                    
                                    if(temp_asc_fprint)
                                        fclose(pfile_pair);
                                }
                                
                                
                                //MatchedPts_list.clear();
                                printf("start AWNCC_MPs %d\t%d\n",MatchedPts_list.size(),total_count_MPs);
                                
                                /*if(level <= proinfo->SGM_py) //SGM
                                    AWNCC_MPs_SGM(proinfo,levelinfo,Size_Grid2D, GridPT3,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,multimps,MatchedPts_list);
                                else*/
                                    AWNCC_MPs(proinfo,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,multimps,MatchedPts_list);
                                
                                count_MPs = MatchedPts_list.size();
                                printf("Done AWNCC_MPs %d\t%d\n",count_MPs,max_count_MPs);
                                for(long count = 0 ; count < Grid_length ; count++)
                                    free(multimps[count]);
                                free(multimps);
                                printf("free multimps\n");
                                
                                min_pair_H = 9999;
                                max_pair_H = -9999;
                                
                                if(temp_asc_fprint)
                                {
                                    sprintf(filename_pair,"%s/txt/matched_pts_%d_%d_%d_%d_%d.txt",proinfo->save_filepath,row,col,level,iteration,levelinfo.pairinfo->NumberOfPairs);
                                    pfile_pair = fopen(filename_pair,"w");
                                }
                                
                                for(long tcnt=0;tcnt<MatchedPts_list.size();tcnt++)
                                {
                                    if(temp_asc_fprint)
                                        fprintf(pfile_pair,"%f\t%f\t%f\n",MatchedPts_list[tcnt].m_X,MatchedPts_list[tcnt].m_Y,MatchedPts_list[tcnt].m_Z);
                                    
                                    if(min_pair_H > MatchedPts_list[tcnt].m_Z)
                                        min_pair_H = MatchedPts_list[tcnt].m_Z;
                                    if(max_pair_H < MatchedPts_list[tcnt].m_Z)
                                        max_pair_H = MatchedPts_list[tcnt].m_Z;
                                }
                                if(temp_asc_fprint)
                                    fclose(pfile_pair);
                                
                                printf("merge row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\t%d\t%f\t%f\n",row,col,level,iteration,count_MPs,MatchedPts_list.size(),min_pair_H,max_pair_H);
                     
                                //exit(1);
                                if(count_MPs > 10)
                                {
                                    if (check_ortho_cal && proinfo->IsRA != 1)
                                    {
                                        printf("blunder detection for anchor points\n");
                                        //anchor points
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*MatchedPts_list.size());
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_all,"w");
                                        for(long count_pt = 0 ; count_pt < MatchedPts_list.size() ; count_pt ++)
                                        {
                                            ptslists[count_pt] = MatchedPts_list[count_pt];
                                            if(temp_asc_fprint)
                                                fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[count_pt].m_X,ptslists[count_pt].m_Y,ptslists[count_pt].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        MatchedPts_list.clear();
                                        vector<D3DPOINT>().swap(MatchedPts_list);
                                        
                                        DecisionMPs(proinfo, levelinfo, false,count_MPs,GridPT3, iteration, Hinterval,count_results_anchor, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_anchor,"w");
                                        long tcnt;
                                        count_results_anchor[0] = 0;
                                        for(tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_anchor.push_back(ptslists[tcnt]);
                                                count_results_anchor[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                            }
                                            ptslists[tcnt].flag = 0;
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd anchor points\n",row,col,level,iteration);
                                        
                                        printf("blunder detection for all points\n");
                                        //blunder detection
                                        DecisionMPs(proinfo, levelinfo, true,count_MPs,GridPT3, iteration, Hinterval,count_results, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_blunder,"w");
                                        count_results[0] = 0;
                                        for(tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_blunder.push_back(ptslists[tcnt]);
                                                count_results[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                                
                                                if(minH_mps > ptslists[tcnt].m_Z)
                                                    minH_mps        = ptslists[tcnt].m_Z;
                                                if(maxH_mps < ptslists[tcnt].m_Z)
                                                    maxH_mps       = ptslists[tcnt].m_Z;
                                            }
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        free(ptslists);
                                        printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
                                        
                                        
                                        if(check_ortho_cal && proinfo->IsRA != 1)
                                        {
                                            printf("settingflag %d\t%d\n",MatchedPts_list_anchor.size(),MatchedPts_list_blunder.size());
                                            count_MPs = SetttingFlagOfGrid(levelinfo, GridPT3, MatchedPts_list_anchor, MatchedPts_list_blunder, &MatchedPts_list_mps);
                                            MatchedPts_list_anchor.clear();
                                            vector<D3DPOINT>().swap(MatchedPts_list_anchor);
                                            MatchedPts_list_blunder.clear();
                                            vector<D3DPOINT>().swap(MatchedPts_list_blunder);
                                        }
                                        printf("count_MPs %d\t%d\n",count_MPs,MatchedPts_list_mps.size());
                                        count_MPs = MatchedPts_list_mps.size();
                                        
                                        
                                        
                                        long i;
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                        
                                        for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                        {
                                            ptslists[i] = MatchedPts_list_mps[i];
                                            //if(level >= 3)
                                            //    ptslists[i].flag = 1; //temporary blunders flag for ortho blunder
                                        }
                                        
                                        //MatchedPts_list_mps.clear();
                                        
                                        //vector<D3DPOINT>().swap(MatchedPts_list_mps);
                                        
                                        FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                        delete origTri;
                                        
                                        count_tri = t_trilists.size();
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        
                                        for(i = 0 ; i < t_trilists.size() ; i++)
                                            trilists[i] = t_trilists[i];
                                        
                                        t_trilists.clear();
                                        vector<UI3DPOINT>().swap(t_trilists);
                                        
                                        fprintf(fid,"level = %d\tMatching Pts = %d\n",level,count_results[0]);
                                        
                                        printf("ortho minmax %f %f pts anchor blunder %d %d \n",minmaxHeight[0],minmaxHeight[1],count_MPs,count_tri);
                                        
                                        count_results[0] = Ortho_blunder(proinfo, levelinfo, MPP_simgle_image, ptslists, count_MPs, trilists,count_tri, GridPT3);
                                        free(trilists);
                                        
                                        printf("end ortho_blunder %d\n",count_results[0]);
                                        
                                        if(level >= 3)
                                        {
                                            int matched_pts = 0;
                                            vector<D3DPOINT> ortho_list;
                                            for(i=0;i<count_MPs;i++)
                                            {
                                                if(ptslists[i].flag != 1)
                                                {
                                                    ortho_list.push_back(ptslists[i]);
                                                    matched_pts++;
                                                }
                                            }
                                            
                                            free(ptslists);
                                            count_MPs = matched_pts;
                                            ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                            
                                            MatchedPts_list_mps.clear();
                                            for( i = 0 ; i < ortho_list.size() ; i++)
                                            {
                                                ptslists[i] = ortho_list[i];
                                                MatchedPts_list_mps.push_back(ptslists[i]);
                                            }
                                            ortho_list.clear();
                                            vector<D3DPOINT>().swap(ortho_list);
                                        }
                                        
                                        printf("load ortho_blunder pts %d\n",count_MPs);
                                        
                                        //Save triangulation and delete it since we will not use it
                                        FullTriangulation *origTri_2 = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                        delete origTri_2;
                                        
                                        count_tri = t_trilists.size();
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        
                                        for(i = 0 ; i < t_trilists.size() ; i++)
                                            trilists[i] = t_trilists[i];
                                        
                                        t_trilists.clear();
                                        vector<UI3DPOINT>().swap(t_trilists);
                                        
                                        DecisionMPs_setheight(proinfo, levelinfo, count_MPs,GridPT3, iteration, Hinterval, minmaxHeight, ptslists,trilists,count_tri);
                                        
                                        free(ptslists);
                                        free(trilists);
                                        
                                        printf("end decision_setheight\n");
                                        
                                        if(pre_matched_pts == 0)
                                            matching_change_rate = 0;
                                        else
                                            matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
                                        
                                        matching_rate = count_MPs/(double)(Accessable_grid);
                                        
                                        printf("matching change rate pre curr %f\t%d\t%d\tmatching rate %f\t%d\n",matching_change_rate,count_MPs,pre_matched_pts,matching_rate,Accessable_grid);
                                        pre_matched_pts = count_MPs;
                                        
                                        //if(level <= 2 && matching_rate < th_mr && proinfo->DEM_resolution <= 4)
                                        //    level_check_matching_rate = true;
                                        
                                        //if(max_count_MPs < count_MPs)
                                        //    max_count_MPs = count_MPs;
                                        
                                        printf("matching change rate pre curr %f\t%d\t%d\tTh_roh %f\t%f\n",matching_change_rate,count_MPs,pre_matched_pts,Th_roh,Th_roh_min);
                                        
                                        final_iteration = iteration;
                                        
                                    }
                                    else
                                    {
                                        printf("blunder detection for all points\n");
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_all,"w");
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*MatchedPts_list.size());
                                        for(long count_pt = 0 ; count_pt < MatchedPts_list.size() ; count_pt ++)
                                        {
                                            ptslists[count_pt] = MatchedPts_list[count_pt];
                                            if(temp_asc_fprint)
                                                fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[count_pt].m_X,ptslists[count_pt].m_Y,ptslists[count_pt].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        MatchedPts_list.clear();
                                        vector<D3DPOINT>().swap(MatchedPts_list);
                                        
                                        DecisionMPs(proinfo, levelinfo, true,count_MPs,GridPT3, iteration, Hinterval,count_results, &minH_mps,&maxH_mps,minmaxHeight, ptslists);
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_blunder,"w");
                                        
                                        count_results[0] = 0;
                                        for(int tcnt=0;tcnt<count_MPs;tcnt++)
                                        {
                                            if(ptslists[tcnt].flag != 1 && ptslists[tcnt].m_X >= subBoundary[0] && ptslists[tcnt].m_X <= subBoundary[2] && ptslists[tcnt].m_Y >= subBoundary[1] && ptslists[tcnt].m_Y <= subBoundary[3])
                                            {
                                                MatchedPts_list_mps.push_back(ptslists[tcnt]);
                                                count_results[0]++;
                                                
                                                if(temp_asc_fprint)
                                                    fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
                                                
                                                if(minH_mps > ptslists[tcnt].m_Z)
                                                    minH_mps        = ptslists[tcnt].m_Z;
                                                if(maxH_mps < ptslists[tcnt].m_Z)
                                                    maxH_mps       = ptslists[tcnt].m_Z;
                                            }
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        count_MPs       = count_results[0];
                                        
                                        free(ptslists);
                                        printf("RA row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
                                        
                                        printf("count_MPs %d\t%d\n",count_MPs,MatchedPts_list_mps.size());
                                        count_MPs = MatchedPts_list_mps.size();
                                        
                                        /*
                                        long i;
                                        
                                        ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                        
                                        for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                            ptslists[i] = MatchedPts_list_mps[i];
                                        
                                        //MatchedPts_list_mps.clear();
                                        
                                        //vector<D3DPOINT>().swap(MatchedPts_list_mps);
                                        /*
                                        //Save triangulation and delete it since we will not use it
                                        FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                        delete origTri;
                                        
                                        count_tri = t_trilists.size();
                                        trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                        
                                        for(i = 0 ; i < t_trilists.size() ; i++)
                                            trilists[i] = t_trilists[i];
                                        
                                        t_trilists.clear();
                                        vector<UI3DPOINT>().swap(t_trilists);
                                        */
                                        if(pre_matched_pts == 0)
                                            matching_change_rate = 0;
                                        else
                                            matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
                                        
                                        printf("matching change rate pre curr %f\t%d\t%d\n",matching_change_rate,count_MPs,pre_matched_pts);
                                        pre_matched_pts = count_results[0];
                                        
                                        matching_rate = count_MPs/(double)(Accessable_grid);
                                        printf("matching change rate pre curr %f\t%d\t%d\tTh_roh %f\t%f\tmatching rate %f\t%d\n",matching_change_rate,count_MPs,pre_matched_pts,Th_roh,Th_roh_min,matching_rate,Accessable_grid);
                                        
                                        //if(level <= 2 && matching_rate < th_mr && proinfo->DEM_resolution <= 4)
                                        //    level_check_matching_rate = true;
                                        //if(max_count_MPs < count_MPs)
                                        //    max_count_MPs = count_MPs;
                                        
                                        bool check_level_end = false;
                                        
                                        if(level != 0)
                                        {
                                            if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th)
                                                check_level_end = true;
                                        }
                                    }
                                    /*
                                    for(long tcnt=0;tcnt<count_MPs;tcnt++)
                                    {
                                        long t_col         = (long)((ptslists[tcnt].m_X - subBoundary[0])/grid_resolution + 0.5);
                                        long t_row         = (long)((ptslists[tcnt].m_Y - subBoundary[1])/grid_resolution + 0.5);
                                        const long ref_index((long)Size_Grid2D.width*t_row + t_col);
                                        
                                        if(t_col >= 0 && t_col < Size_Grid2D.width && t_row >=0 && t_row < Size_Grid2D.height)
                                        {
                                            multimps[ref_index][pair_number].peak_roh = ptslists[tcnt].m_roh;
                                            multimps[ref_index][pair_number].peak_height = ptslists[tcnt].m_Z;
                                            multimps[ref_index][pair_number].check_matched = true;
                                            
                                            //printf("pos %d\t%d\t%d\t%f\n",t_col,t_row,multimps[ref_index][pair_number].peak_roh,multimps[ref_index][pair_number].peak_height);
                                        }
                                    }
                                     */
                                    //exit(1);
                                    
                                    //free(ptslists);
                                    //free(trilists);
                                    
                                    lower_level_match = true;
                                }
                                //else
                                //    lower_level_match = false;
                                
                                //MatchedPts_list.clear();
                                //vector<D3DPOINT>().swap(MatchedPts_list);
                                
                                fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tcheck = %d(%d)\tEnd blunder detection\n",row,col,level,iteration,lower_level_match,count_MPs);
                                
                                printf("End computation of blunder!! Mps = %d\tTris = %d\tminz Mp = %f\tmaxz Mp = %f minmax %f %f \n",
                                       count_results[0],count_results[1],minH_mps,maxH_mps,minmaxHeight[0],minmaxHeight[1]);
                                    
                                //exit(1);
                                
                                if(lower_level_match)
                                {
                                    /*
                                    MatchedPts_list_mps.clear();
                                    printf("start AWNCC_MPs %d\n",MatchedPts_list_mps.size());
                                    AWNCC_MPs(proinfo,levelinfo,Size_Grid2D, GridPT3,nccresult,height_step,level,iteration,MaxNumberofHeightVoxel,minmaxHeight,multimps,MatchedPts_list_mps);
                                    
                                    printf("Done AWNCC_MPs %d\t%d\n",MatchedPts_list_mps.size(),max_count_MPs);
                                    for(long count = 0 ; count < Grid_length ; count++)
                                        free(multimps[count]);
                                    free(multimps);
                                    printf("free multimps\n");
                                     */
                                    //exit(1);
                                    count_MPs = MatchedPts_list_mps.size();
                                    ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
                                    
                                    for(long i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                        ptslists[i] = MatchedPts_list_mps[i];
                                    
                                    if(temp_asc_fprint)
                                    {
                                        pMT_all = fopen(filename_mps_asc,"w");
                                        for(long i = 0 ; i < count_MPs ; i++)
                                        {
                                            fprintf(pMT_all,"%f\t%f\t%f\n",ptslists[i].m_X,ptslists[i].m_Y,ptslists[i].m_Z);
                                        }
                                        fclose(pMT_all);
                                    }
                                    printf("print results\n");
                                    //Save triangulation and delete it since we will not use it
                                    FullTriangulation *origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, grid_resolution);
                                    delete origTri;
                                    
                                    count_tri = t_trilists.size();
                                    trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
                                    
                                    for(int i = 0 ; i < t_trilists.size() ; i++)
                                        trilists[i] = t_trilists[i];
                                    
                                    t_trilists.clear();
                                    vector<UI3DPOINT>().swap(t_trilists);
                                }
                
                                printf("Done AWNCC\n");
                                
                                if(lower_level_match)
                                {
                                    count_MPs = MatchedPts_list_mps.size();
                                    if(level == 0 && iteration == 3)
                                    {
                                        D3DPOINTSAVE *ptslists_save = (D3DPOINTSAVE*)calloc(count_MPs,sizeof(D3DPOINTSAVE));
                                        double minmaxBR[6] = {10000000, 10000000, -10000000, -10000000, 100000, -100000};
                                        
                                        int i = 0;
                                        
                                        if(temp_asc_fprint)
                                            pMT_all = fopen(filename_mps_asc,"w");
                                        for( i = 0 ; i < MatchedPts_list_mps.size() ; i++)
                                        {
                                            ptslists_save[i].m_X = MatchedPts_list_mps[i].m_X;
                                            ptslists_save[i].m_Y = MatchedPts_list_mps[i].m_Y;
                                            ptslists_save[i].m_Z = MatchedPts_list_mps[i].m_Z;
                                            
                                            if(minmaxBR[0] > ptslists_save[i].m_X)
                                                minmaxBR[0]     = ptslists_save[i].m_X;
                                            if(minmaxBR[1] > ptslists_save[i].m_Y)
                                                minmaxBR[1]     = ptslists_save[i].m_Y;
                                            
                                            if(minmaxBR[2] < ptslists_save[i].m_X)
                                                minmaxBR[2]     = ptslists_save[i].m_X;
                                            if(minmaxBR[3] < ptslists_save[i].m_Y)
                                                minmaxBR[3]     = ptslists_save[i].m_Y;
                                            if(minmaxBR[4] > ptslists_save[i].m_Z)
                                                minmaxBR[4] = ptslists_save[i].m_Z;
                                            if(minmaxBR[5] < ptslists_save[i].m_Z)
                                                minmaxBR[5] = ptslists_save[i].m_Z;
                                            
                                            if(temp_asc_fprint)
                                                fprintf(pMT_all,"%f\t%f\t%f\n",ptslists_save[i].m_X,ptslists_save[i].m_Y,ptslists_save[i].m_Z);
                                        }
                                        if(temp_asc_fprint)
                                            fclose(pMT_all);
                                        
                                        FILE *pFile = fopen(filename_mps,"wb");
                                        fwrite(ptslists_save,sizeof(D3DPOINTSAVE),count_MPs,pFile);
                                        fclose(pFile);
                                        
                                        if(!proinfo->IsRA)
                                        {
                                            FILE *fid_BR;
                                            FILE *fid_count;
                                            
                                            char save_file[500];
                                            sprintf(save_file,"%s/txt/matched_BR_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                                            fid_BR  = fopen(save_file,"w");
                                            fprintf(fid_BR,"%f\t%f\t%f\t%f\t%f\t%f\n",minmaxBR[0],minmaxBR[1],minmaxBR[2],minmaxBR[3],minmaxBR[4],minmaxBR[5]);
                                            fclose(fid_BR);
                                            
                                            sprintf(save_file,"%s/txt/count_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                                            fid_count  = fopen(save_file,"w");
                                            
                                            fprintf(fid_count,"%d\n",count_MPs);
                                            fclose(fid_count);
                                            
                                            if(temp_asc_fprint)
                                            {
                                                echoprint_Gridinfo_asc(proinfo,levelinfo, nccresult,row,col,level,iteration,Size_Grid2D,GridPT3);
                                            }
                                            
                                            echoprint_Gridinfo(proinfo,nccresult,row,col,level,iteration,0,&Size_Grid2D,GridPT3,(char*)"final");
                                        }
                                        free(ptslists_save);
                                        
                                        matching_change_rate = 0.001;
                                    }
                                    else
                                    {
                                        if(iteration > max_iteration)
                                            matching_change_rate = 0.001;
                                        
                                        if(level == 0)
                                        {
                                            if(proinfo->DEM_resolution < 2)
                                                matching_change_rate = 0.001;
                                            if(iteration > 2)
                                                matching_change_rate = 0.001;
                                        }
                                        
                                        if(level <= 1)
                                        {
                                            if(iteration >= ceil(max_iteration/2.0))
                                                matching_change_rate = 0.001;
                                        }
                                        
                                        if(proinfo->IsRA)
                                            matching_change_rate = 0.001;
                                        
                                        if(proinfo->DEM_resolution >= 8)
                                            matching_change_rate = 0.001;
                                        
                                        if(Th_roh >= Th_roh_min)
                                        {
                                            if(level == 0)
                                            {
                                                Th_roh_update       = (double)(Th_roh - 0.10);
                                                if(iteration >= 2)
                                                    Th_roh_update          = (double)(Th_roh - 0.50);
                                            }
                                            else if(level == 1)
                                                Th_roh_update       = (double)(Th_roh - 0.10);
                                            else if(level == 2)
                                                Th_roh_update       = (double)(Th_roh - 0.10);
                                            else if(level == 3)
                                                Th_roh_update       = (double)(Th_roh - 0.10);
                                            else
                                            {
                                                if(proinfo->IsRA)
                                                    Th_roh_update       = (double)(Th_roh - 0.10);
                                                else
                                                    Th_roh_update       = (double)(Th_roh - 0.06);
                                            }
                                        }
                                        
                                        if(temp_asc_fprint)
                                            echoprint_Gridinfo_asc(proinfo, levelinfo, nccresult,row,col,level,iteration,Size_Grid2D,GridPT3);
                                        
                                        if(level == 0)
                                        {
                                            if(MPP_stereo_angle > 5)
                                                MPP = MPP_stereo_angle;
                                            else
                                                MPP = MPP_simgle_image;
                                            if(proinfo->DEM_resolution < 2)
                                            {
                                                Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                            else
                                            {
                                                GridPT3         = SetHeightRange(proinfo,levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                            }
                                        }
                                        else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
                                        {
                                            if(MPP_stereo_angle > 5)
                                            {
                                                if(level > 2)
                                                    MPP = MPP_stereo_angle;
                                                else
                                                    MPP = MPP_simgle_image;
                                            }
                                            else
                                                MPP = MPP_simgle_image;
                                            
                                            Pre_GridPT3     = SetHeightRange(proinfo, levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                        }
                                        else
                                        {
                                            if(MPP_stereo_angle > 5)
                                            {
                                                if(level <= 1)
                                                    MPP = MPP_stereo_angle;
                                                else
                                                    MPP = MPP_simgle_image;
                                            }
                                            else
                                                MPP = MPP_simgle_image;
                                            
                                            GridPT3     = SetHeightRange(proinfo,levelinfo, count_MPs, count_tri, GridPT3, iteration, &minH_grid, &maxH_grid, ptslists, trilists, MPP, check_matching_rate);
                                        }
                                        
                                        
                                        /*
                                        if(proinfo->IsRA)
                                        {
                                            if(level <= 3)
                                            {
                                                int RA_iter_counts = 0;
                                                RA_iter_counts = AdjustParam(proinfo, levelinfo, count_MPs, t_Imageparams, pyramid_step, ptslists);
                                                for(int ti = 0 ; ti < levelinfo.pairinfo->NumberOfPairs ; ti++)
                                                {
                                                    //if(proinfo->check_selected_image[ti])
                                                    {
                                                        fprintf(fid,"RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Imageparams[ti][0],t_Imageparams[ti][1]);
                                                        printf("RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Imageparams[ti][0],t_Imageparams[ti][1]);
                                                    }
                                                }
                                                
                                                if (level <= 1)
                                                {
                                                    char save_file[500];
                                                    sprintf(save_file,"%s/txt/RAinfo.txt",proinfo->save_filepath);
                                                    FILE *fid_RAinfo  = fopen(save_file,"w");
                                                    for(int ti = 0 ; ti < levelinfo.pairinfo->NumberOfPairs ; ti++)
                                                        fprintf(fid_RAinfo,"%d\t%d\t%f\t%f\t%f\n",levelinfo.pairinfo->pairs[ti].m_X,levelinfo.pairinfo->pairs[ti].m_Y,t_Imageparams[ti][0],t_Imageparams[ti][1],levelinfo.pairinfo->BHratio[ti]);
                                                    fclose(fid_RAinfo);
                                                }
                                            }
                                            
                                            if(check_RA_divide)
                                            {
                                                if(level <= 3 && iteration > 2)
                                                {
                                                    int* t_count = (int*)calloc(sizeof(int),total_tile);
                                                    
                                                    for(int k=0;k<count_MPs;k++)
                                                    {
                                                        int t_col = floor((ptslists[k].m_X - subBoundary[0])/(double)tilesize_RA);
                                                        int t_row = floor((ptslists[k].m_Y - subBoundary[1])/(double)tilesize_RA);
                                                        
                                                        t_count[t_col+division_X*t_row] ++;
                                                    }
                                                    
                                                    int saved_count = 0;
                                                    int selected_X = 0;
                                                    int selected_Y = 0;
                                                    int total_count = 0;
                                                    for(int k=0;k<total_tile;k++)
                                                    {
                                                        if(t_count[k] > saved_count)
                                                        {
                                                            saved_count = t_count[k];
                                                            selected_Y  = floor(k/division_X);
                                                            selected_X  = k % division_X;
                                                        }
                                                        total_count += t_count[k];
                                                    }
                                                    printf("total_count %d\tsaved_count %d\tselected_X %d\tselected_Y %d\n",total_count,saved_count,selected_X,selected_Y);
                                                    printf("selected br %f\t%f\t%f\t%f\n",subBoundary[0] + selected_X*tilesize_RA,subBoundary[1] + selected_Y*tilesize_RA,
                                                           subBoundary[0] + (selected_X+1)*tilesize_RA,subBoundary[1] + (selected_Y+1)*tilesize_RA);
                                                    
                                                    fprintf(fid,"total_count %d\tsaved_count %d\tselected_X %d\tselected_Y %d\n",total_count,saved_count,selected_X,selected_Y);
                                                    
                                                    new_subBoundary_RA[0] = subBoundary[0] + selected_X*tilesize_RA;
                                                    new_subBoundary_RA[1] = subBoundary[1] + selected_Y*tilesize_RA;
                                                    new_subBoundary_RA[2] = subBoundary[0] + (selected_X+1)*tilesize_RA;
                                                    new_subBoundary_RA[3] = subBoundary[1] + (selected_Y+1)*tilesize_RA;
                                                    
                                                    check_new_subBoundary_RA = true;
                                                    check_RA_divide = false;
                                                    
                                                    free(t_count);
                                                }
                                            }
                                        }
                                        */
                                    }
                                    
                                    MatchedPts_list_mps.clear();
                                    vector<D3DPOINT>().swap(MatchedPts_list_mps);
                                    
                                    free(trilists);
                                    free(ptslists);
                                    
                                    fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\tmatching_rate = %f\n",
                                            row,col,level,iteration,count_MPs,minmaxHeight[0],minmaxHeight[1],matching_rate);
                                    printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\n",
                                           row,col,level,iteration,count_MPs,minmaxHeight[0],minmaxHeight[1]);
                                }
                                    
                                //exit(1);
                                    
                                if(lower_level_match)
                                {
                                    flag_start          = true;
                                    iteration++;
                                }
                                else if(level == 0)
                                    iteration++;
                                
                                if(level == 0)
                                {
                                    if(proinfo->DEM_resolution >= 2)
                                        Th_roh          = (double)(Th_roh - 0.10);
                                    else
                                        Th_roh          = (double)(Th_roh - 0.50);
                                    
                                    if(iteration > 3)
                                        Th_roh          = (double)(Th_roh - 0.50);
                                }
                                else if(level == 1)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else if(level == 2)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else if(level == 3)
                                    Th_roh          = (double)(Th_roh - 0.10);
                                else
                                {
                                    if(proinfo->IsRA)
                                        Th_roh          = (double)(Th_roh - 0.10);
                                    else
                                        Th_roh          = (double)(Th_roh - 0.06);
                                }
                                
                                if(lower_level_match)
                                {
                                    if(Th_roh < Th_roh_min && matching_change_rate > rate_th)
                                    {
                                        if(level == 0)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 1)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 2)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else if(level == 3)
                                            Th_roh          = (double)(Th_roh + 0.10);
                                        else
                                        {
                                            if(proinfo->IsRA)
                                                Th_roh          = (double)(Th_roh + 0.10);
                                            else
                                                Th_roh          = (double)(Th_roh + 0.06);
                                        }
                                    }
                                }
                                
                                if (!lower_level_match && Th_roh < Th_roh_min)
                                {
                                    if(level > 0)
                                    {
                                        iteration++;
                                        matching_change_rate = 0.001;
                                        Th_roh_min = 0.4;
                                    }
                                }
                                
                                if(level == 0)
                                    final_level_iteration = iteration;
                                
                                printf("Memory : System %f\t SETSM required %f\n",proinfo->System_memory, total_memory);
                            }
                        }
                        
                        if(flag_start)
                        {
                            double min_after   = (double)(minH_mps - pwrtwo(level)*10*MPP);
                            double max_after   = (double)(maxH_mps + pwrtwo(level)*10*MPP);
                            if(level <= 2)
                            {
                                if(minmaxHeight[0] < min_after)
                                    minmaxHeight[0]     = (double)(floor(min_after));
                                if(minmaxHeight[1] > max_after)
                                    minmaxHeight[1]     = (double)(ceil(max_after));
                            }
                            else
                            {
                                if(min_after > minH_grid)
                                    min_after = minH_grid;
                                if(max_after < maxH_grid)
                                    max_after = maxH_grid;
                                
                                minmaxHeight[0]     = min_after;
                                minmaxHeight[1]     = max_after;
                            }
                            printf("minmax %f\t%f\t\n", minmaxHeight[0],minmaxHeight[1]);
                            fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd of level processing!! minmaxHeight = [%f \t%f]\n",
                                    row,col,level,iteration,minmaxHeight[0],minmaxHeight[1]);
                        }
                        
                        printf("\trow = %d/%d\tcol = %d/%d\tDEM generation(%%) = %4.2f%% !!\n",row,iter_row_end,col,t_col_end,(double)(pyramid_step+1 - level)/(double)(pyramid_step+1)*100);
                        
                        if(proinfo->IsRA)
                        {
                            if(!lower_level_match)
                            {
                                lower_level_match   = true;
                                flag_start          = false;
                            }
                        }
                        else
                        {
                            if(!lower_level_match && level > 2)
                            {
                                lower_level_match   = true;
                                flag_start          = false;
                            }
                        }
                        
                        printf("release Grid_wgs, nccresult\n");
                        free(GridPT);
                        free(Grid_wgs);
                        
                        if(!check_matching_rate)
                        {
                            grid_voxel.clearall();
                            printf("free grid_voxel\n");
                        }
                        
                        if(!check_matching_rate)
                            check_matching_rate = level_check_matching_rate;
                        
                        free(nccresult);
                        free(Startpos);
                        free(BStartpos);
                        
                        if(level > Py_combined_level)
                            free(Startpos_next);
                        
                        if(level > 0)
                            level   = level - 1;
                        
                        if(level == 0 && final_level_iteration == 4)
                            level = -1;
                        
                    }
                    printf("relese data size\n");
                    
                    for(int t_level = 0 ; t_level < py_level_set ; t_level++)
                    {
                        for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                        {
                            if(proinfo->check_selected_image[ti])
                            {
                                free(SubImages[t_level][ti]);
                                free(SubMagImages[t_level][ti]);
                                free(SubOriImages[t_level][ti]);
                            }
                        }
                        free(SubImages[t_level]);
                        free(SubMagImages[t_level]);
                        free(SubOriImages[t_level]);
                    }
                    
                    free(SubImages);
                    free(SubMagImages);
                    free(SubOriImages);
                    
                    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                    {
                        if(proinfo->check_selected_image[ti])
                             free(data_size_lr[ti]);
                    }
                    free(data_size_lr);
                    
                    for(long int t_i = 0 ; t_i < (long)Size_Grid2D.width*(long)Size_Grid2D.height ; t_i++)
                        free(GridPT3[t_i].ortho_ncc);
                    free(GridPT3);
                    
                    
                    free(pairinfo.pairs);
                    free(pairinfo.BHratio);
                    
                    printf("release GridTP3\n");
                    PreET = time(0);
                    Pregab = difftime(PreET,PreST);
                    printf("row = %d/%d \tcol = %d/%d\tDEM generation finish(time[m] = %5.2f)!!\n",row,iter_row_end,col,t_col_end,Pregab/60.0);
                }
            }
            
            free(SourceImages);
            
            fclose(fid);
            fclose(fid_header);
            printf("proinfo -ISRA\n");
            if(proinfo->IsRA)
            {
                for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
                {
                    if(t_Imageparams[ti][0] != 0 && t_Imageparams[ti][1] != 0)
                    {
                        RA_count[ti]++;
                        Imageparams[ti][0] += t_Imageparams[ti][0];
                        Imageparams[ti][1] += t_Imageparams[ti][1];
                    }
                }
            }
            
            printf("free t_Imageparams\n");
            for(int ti = 0 ; ti < MaxNCC ; ti++)
            {
                //if(proinfo->check_selected_image[ti])
                    free(t_Imageparams[ti]);
            }
            free(t_Imageparams);
             
            free(Startpos_ori);
            free(Subsetsize);
            printf("free Subsetsize\n");
        }
    }
    
    free(iterations);
    if(proinfo->IsRA)
    {
        for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
        {
            if(Imageparams[ti][0] != 0 && Imageparams[ti][1] != 0 && RA_count[ti] > 0)
            {
                Imageparams[ti][0] /= RA_count[ti];
                Imageparams[ti][1] /= RA_count[ti];
                
                printf("RPC bias %d\t%f\t%f\t%d\t%d\t%f\n",ti,Imageparams[ti][0],Imageparams[ti][1],pairinfo_return.pairs[ti].m_X,pairinfo_return.pairs[ti].m_Y,pairinfo_return.BHratio[ti]);
            }
        }
    }
#ifdef BUILDMPI
    for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
        MPI_Bcast(Imageparams[ti], 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    for(int ti = 0 ; ti < pairinfo_return.NumberOfPairs ; ti++)
        printf("Num of RAs = %d\tRA param = %f\t%f\n",RA_count[ti],Imageparams[ti][0],Imageparams[ti][1]);
    free(RA_count);
    
    return final_iteration;
}

double CalMemorySize(const ProInfo *info,LevelInfo &plevelinfo,const UGRID *GridPT3, double *minimum_memory, const uint8 iteration,const double *minmaxHeight)
{
    double Memory = 0;
    
    double image = 0;
    double oriimage = 0;
    double magimage = 0;
    
    const int level = *(plevelinfo.Pyramid_step);
    const int blunder_selected_level = *(plevelinfo.blunder_selected_level);
    printf("%d\t%d\n",level,blunder_selected_level);
    
    for(int iter_level = 0 ; iter_level <= level ; iter_level++)
    {
        for(int image_index = 0 ; image_index < info->number_of_images ; image_index++)
        {
            long int data_length = (long int)plevelinfo.py_Sizes[image_index][iter_level].height*(long int)plevelinfo.py_Sizes[image_index][iter_level].width;
                
            image += (double)(sizeof(uint16)*data_length)*2.0;
            image += (double)(sizeof(uint8)*data_length);
        }
    }
    
    //printf("%f\t%f\t%f\t%f\n",subBoundary[0],subBoundary[1],subBoundary[2],subBoundary[3]);
    Memory += (image + oriimage + magimage);
    //printf("memory 1 %f\t%f\t%f\t%f\n",Memory,image,oriimage,magimage);
    long int GridPT  = (double)(sizeof(D2DPOINT)*(*(plevelinfo.Grid_length)));
    Memory += (GridPT)*3;
    //printf("memory 2 %f\n",Memory);
    long int GridPT3_size = (double)(sizeof(UGRID)*(*(plevelinfo.Grid_length)));
    Memory += (GridPT3_size);
    //printf("memory 3 %f\n",Memory);
    long int nccresult_size = (double)(sizeof(NCCresult)*(*(plevelinfo.Grid_length)));
    Memory += (nccresult_size);
    
    if(plevelinfo.pairinfo->NumberOfPairs > 1)
    {
        long int multimps_size = (double)(sizeof(MultiMPs)*(*(plevelinfo.Grid_length))*plevelinfo.pairinfo->NumberOfPairs);
        Memory += (multimps_size);
    }
    //printf("memory 4 %f\n",Memory);
    
    bool check_ortho = true;
    if((level == 4 && iteration == 1) || info->IsRA == true)
        check_ortho = false;
    if(check_ortho)
    {
        int sub_imagesize_w, sub_imagesize_h;
        int sub_imagesize_w_next, sub_imagesize_h_next;
        double all_im_cd = 0;
        double all_im_cd_next = 0;
        double im_resolution = info->resolution*pwrtwo(level);
        double im_resolution_next;
        long int sub_imagesize_total_next;
        
        if(level > 0)
            im_resolution_next = im_resolution*pwrtwo(level-1);
        
        sub_imagesize_w = (int)((plevelinfo.Boundary[2] - plevelinfo.Boundary[0])/im_resolution)+1;
        sub_imagesize_h = (int)((plevelinfo.Boundary[3] - plevelinfo.Boundary[1])/im_resolution)+1;
        
        if(level > 0)
        {
            sub_imagesize_w_next = (int)((plevelinfo.Boundary[2] - plevelinfo.Boundary[0])/im_resolution_next)+1;
            sub_imagesize_h_next = (int)((plevelinfo.Boundary[3] - plevelinfo.Boundary[1])/im_resolution_next)+1;
            sub_imagesize_total_next = (long int)sub_imagesize_w_next * (long int)sub_imagesize_h_next;
        }
        
        long int sub_imagesize_total = (long int)sub_imagesize_w * (long int)sub_imagesize_h;
        
        for(int ti = 0 ; ti < info->number_of_images ; ti++)
        {
            all_im_cd += (double)(sizeof(D2DPOINT)*sub_imagesize_total);
            if(level > 0)
                all_im_cd_next += (double)(sizeof(D2DPOINT)*sub_imagesize_total_next);
        }
        Memory += (all_im_cd + all_im_cd_next);
    }
    *minimum_memory = (double)(Memory/1024.0/1024.0/1024.0);
    
    if(!info->IsRA)
    {
        long int grid_voxel = (double)(sizeof(VOXEL)*(*(plevelinfo.Grid_length)));
        grid_voxel += (double)(sizeof(float)*(*(plevelinfo.Grid_length)));
        for(long int t_i = 0 ; t_i < (*(plevelinfo.Grid_length)); t_i++)
        {
            int selected_images;
            if(check_image_boundary(info,plevelinfo,plevelinfo.GridPts[t_i],plevelinfo.Grid_wgs[t_i],GridPT3[t_i].minHeight,GridPT3[t_i].maxHeight,7,selected_images))
            {
                bool check_blunder_cell = true;
                double th_height = 1000;
                
                if ( level >= 2)
                    check_blunder_cell = false;
                else if( GridPT3[t_i].Matched_flag != 0)
                {
                    if(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight > 0)
                    {
                        if(level <= 1)
                        {
                            if(GridPT3[t_i].maxHeight <= minmaxHeight[1] && GridPT3[t_i].minHeight >= minmaxHeight[0])
                                check_blunder_cell = false;
                        }
                        
                        if(level == 1)
                        {
                            if((abs(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight) < 1000))
                                check_blunder_cell = false;
                        }
                        
                        if(level == 0)
                        {
                            if(abs(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight) < th_height)
                                check_blunder_cell = false;
                        }
                    }
                }
                
                if(!check_blunder_cell)
                {
                    const int NumberofHeightVoxel = (int)((GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight)/(*(plevelinfo.height_step)));
                    
                    if(NumberofHeightVoxel > 0 )
                    {
                        grid_voxel += (double)(NumberofHeightVoxel*(sizeof(VOXEL)));
                        //grid_voxel += (double)(NumberofHeightVoxel*(sizeof(short)*plevelinfo.pairinfo->NumberOfPairs));
                    }
                }
            }
        }
        Memory += (grid_voxel);
    }
    //printf("memory 5 %f\n",Memory);
    
    double result = (double)(Memory/1024.0/1024.0/1024.0);
    
    return result;
}


bool SetupParam(ProInfo *info, bool *pre_DEMtif)
{
    bool cal_check      = true;
    
    printf("seed flag %d\n",info->pre_DEMtif);
    
    if (info->pre_DEMtif)
    {
        FILE *pFile_DEM;
        printf("dem filename %s\n",info->priori_DEM_tif);
        pFile_DEM   = fopen(info->priori_DEM_tif,"r");
        if(pFile_DEM)
            *pre_DEMtif = true;
        else
            *pre_DEMtif = false;
        fclose(pFile_DEM);
    }
    
    printf("seed flag %d\n",*pre_DEMtif);

    return cal_check;
}


void SetTiles(const ProInfo *info, const double *Boundary, const int tile_size, uint8 *pyramid_step, uint16 *buffer_area, uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, double *subX, double *subY)
{
    int division_X, division_Y;
    SetsubsetBR(info, Boundary, tile_size, subX, subY, division_X, division_Y);
    printf("dx = %d\tdy = %d\t%f\t%f\n", division_X, division_Y, *subX, *subY);
    
    if(info->pyramid_level != 4)
        *pyramid_step = info->pyramid_level;
    else
    {
        if(info->sensor_provider == PT)
            *pyramid_step   = 2;
        else
            *pyramid_step   = 4;
    }
    
    
    if(info->pre_DEMtif)
    {
        if(info->seedDEMsigma <= 15)
            *pyramid_step   = 2;
        else if(info->seedDEMsigma <= 30)
            *pyramid_step   = 3;
    }
    
    if(info->DEM_resolution  >= 10)
        *buffer_area    = (uint16)(*buffer_area * 1.5);

    *iter_row_start = 1;
    *iter_row_end   = division_Y+1;
    
    *t_col_start    = 1;
    *t_col_end      = division_X+1;
    
    if(*iter_row_start == *iter_row_end)
        *iter_row_end += 1;
    if(*t_col_start == *t_col_end)
        *t_col_end += 1;
}

void SetTiles_RA(const ProInfo *info, const double *Boundary, const int tile_size, uint8 *pyramid_step, uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, double *subX, double *subY)
{
    int division_X, division_Y;
    SetsubsetBR(info, Boundary, tile_size, subX, subY, division_X, division_Y);
    
    printf("dx = %d\tdy = %d\t%f\t%f\n", division_X, division_Y, *subX, *subY);
    
    *pyramid_step = 4;
    *RA_col_iter = 2;
    *RA_row_iter = 2;
    
    const uint8 ceilDivisionX = (uint8) ceil(division_X / 2.0);
    *RA_row_start = ceilDivisionX;
    *RA_row_end = (uint8) (ceilDivisionX + 1);
    *t_col_start = ceilDivisionX;
    *t_col_end = (uint8) (ceilDivisionX + 1);
}

void SetsubsetBR(const ProInfo *info, const double *Boundary, const int tile_size, double *subX, double *subY, int &division_X, int &division_Y)
{
    int lengthOfX = Boundary[2] - Boundary[0];
    int lengthOfY = Boundary[3] - Boundary[1];
    division_X = (int) (ceil(lengthOfX / (double)(tile_size)));
    division_Y = (int) (ceil(lengthOfY / (double)(tile_size)));
    *subX = floor((ceil(ceil(lengthOfX / division_X) / info->DEM_resolution) * info->DEM_resolution) / 2) * 2;
    *subY = floor((ceil(ceil(lengthOfY / division_Y) / info->DEM_resolution) * info->DEM_resolution) / 2) * 2;
}

void SetThs(const ProInfo *proinfo,const int level, const int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start)
{
    if(proinfo->IsRA)
    {
        if(level >= 4)
        {
            *Th_roh          = (double)(0.80);
            *Th_roh_min      = (double)(0.39);
            
            if(proinfo->IsRA)
            {
                *Th_roh          = (double)(0.80);
                *Th_roh_min      = (double)(0.49);
            }
        }
        else if(level >= 3)
        {
            *Th_roh          = (double)(0.70);
            *Th_roh_min      = (double)(0.39);
        }
        else if(level == 2)
        {
            *Th_roh          = (double)(0.60);
            *Th_roh_min      = (double)(0.49);
        }
        else if(level == 1)
        {
            *Th_roh          = (double)(0.50);
            *Th_roh_min      = (double)(0.39);
        }
        else if(level == 0)
        {
            *Th_roh          = (double)(0.40 - 0.10*(final_level_iteration-1));
            *Th_roh_min      = (double)(0.19);
        }
        
        
        if(level == 5)
            *Th_roh_next        = (double)(0.80);
        else if(level == 4)
            *Th_roh_next        = (double)(0.70);
        else if(level == 3)
            *Th_roh_next        = (double)(0.60);
        else if(level == 2)
            *Th_roh_next        = (double)(0.50);
        else if(level == 1)
            *Th_roh_next        = (double)(0.40);
        else
            *Th_roh_next        = (double)(0.20);
        
        *Th_roh_start       = (double)(*Th_roh);
    }
    else
    {
        if(proinfo->DEM_resolution >= 8)
        {
            if(level >= 4)
            {
                *Th_roh          = (double)(0.80);
                *Th_roh_min      = (double)(0.39);
                
                if(proinfo->IsRA)
                {
                    *Th_roh          = (double)(0.80);
                    *Th_roh_min      = (double)(0.49);
                }
            }
            else if(level >= 3)
            {
                *Th_roh          = (double)(0.70);
                *Th_roh_min      = (double)(0.39);
            }
            else if(level == 2)
            {
                *Th_roh          = (double)(0.60);
                *Th_roh_min      = (double)(0.49);
            }
            else if(level == 1)
            {
                *Th_roh          = (double)(0.50);
                *Th_roh_min      = (double)(0.39);
            }
            else if(level == 0)
            {
                *Th_roh          = (double)(0.40 - 0.10*(final_level_iteration-1));
                *Th_roh_min      = (double)(0.19);
            }

            
            if(level == 5)
                *Th_roh_next        = (double)(0.80);
            else if(level == 4)
                *Th_roh_next        = (double)(0.70);
            else if(level == 3)
                *Th_roh_next        = (double)(0.60);
            else if(level == 2)
                *Th_roh_next        = (double)(0.50);
            else if(level == 1)
                *Th_roh_next        = (double)(0.40);
            else
                *Th_roh_next        = (double)(0.20);
       
            *Th_roh_start       = (double)(*Th_roh);
        }
        else
        {
            if(level >= 4)
            {
                *Th_roh          = (double)0.60;//(0.80);
                *Th_roh_min      = (double)0.19;//(0.39);
            }
            else if(level >= 3)
            {
                *Th_roh          = (double)0.50;//(0.70);
                *Th_roh_min      = (double)0.19;//(0.39);
            }
            else if(level == 2)
            {
                *Th_roh          = (double)(0.40);
                *Th_roh_min      = (double)(0.09);
            }
            else if(level == 1)
            {
                *Th_roh          = (double)(0.30);
                *Th_roh_min      = (double)(0.09);
            }
            else if(level == 0)
            {
                *Th_roh          = (double)(0.30 - 0.10*(final_level_iteration-1));
                *Th_roh_min      = (double)(0.09);
            }
            
            
            if(level == 5)
                *Th_roh_next        = (double)0.60;//(0.80);
            else if(level == 4)
                *Th_roh_next        = (double)0.50;//(0.70);
            else if(level == 3)
                *Th_roh_next        = (double)(0.40);
            else if(level == 2)
                *Th_roh_next        = (double)(0.30);
            else if(level == 1)
                *Th_roh_next        = (double)(0.30);
            else
                *Th_roh_next        = (double)(0.20);
            
            *Th_roh_start       = (double)(*Th_roh);
        }
        
        if(proinfo->check_Matchtag)
        {
            *Th_roh          = (double)(0.35 - 0.15*(final_level_iteration-2));
            *Th_roh_min      = (double)(0.19);
            
            *Th_roh_start       = (double)(*Th_roh);
        }
    }
}

D2DPOINT *SetGrids(const ProInfo *info, const int level, const int final_level_iteration, const double resolution, CSize *Size_Grid2D, const double DEM_resolution, double *py_resolution, double *grid_resolution, const double *subBoundary)
{
    *py_resolution   = (double)(resolution*pwrtwo(level));
    *grid_resolution = *py_resolution;
    
    printf("pre resolution %f\t level %d\t final_level_iteration %d\n",*py_resolution,level,final_level_iteration);
    
    if(*py_resolution > 32)//low-res original imagery
    {
        *py_resolution = (int)(*py_resolution/4.0);
        *grid_resolution = *py_resolution;
    }
    else
    {
        if(resolution >= 0.4)
        {
            if(level > 0)
            {
                if((*py_resolution)*3 > DEM_resolution) //low-res DEM more than 8m
                {
                    if(DEM_resolution > 8)
                    {
                        *py_resolution = DEM_resolution;
                    }
                    else
                    {
                        if((*py_resolution)*3 > 8)
                            *py_resolution = 8;
                        else
                            *py_resolution   = (*py_resolution)*3;
                    }
                }
                else
                    *py_resolution = DEM_resolution;
            }
            else if(level == 0)
            {
                if(DEM_resolution >= 2)
                    *py_resolution = DEM_resolution;
                else
                {
                    if(final_level_iteration == 1)
                    {
                        if(*py_resolution < 2)
                        {
                            if((*py_resolution)*4 > DEM_resolution)
                                *py_resolution   = (*py_resolution)*4;
                            else
                                *py_resolution   = DEM_resolution;
                        }
                        else
                        {
                            *py_resolution = DEM_resolution;
                        }
                    }
                    else if(final_level_iteration == 2)
                    {
                        if(*py_resolution < 2)
                        {
                            if((*py_resolution)*2 > DEM_resolution)
                                *py_resolution   = (*py_resolution)*2;
                            else
                                *py_resolution   = DEM_resolution;
                        }
                        else
                        {
                            *py_resolution = DEM_resolution;
                        }
                    }
                    else
                        *py_resolution   = DEM_resolution;
                }
                
            }
        }
        else
        {
            if(DEM_resolution >= 2)
                *py_resolution = DEM_resolution;
            else
            {
                if(level == 0)
                {
                    if(final_level_iteration == 1)
                    {
                        if((*py_resolution)*4 > DEM_resolution)
                            *py_resolution   = (*py_resolution)*4;
                        else
                            *py_resolution   = DEM_resolution;
                    }
                    else if(final_level_iteration == 2)
                    {
                        if((*py_resolution)*2 > DEM_resolution)
                            *py_resolution   = (*py_resolution)*2;
                        else
                            *py_resolution   = DEM_resolution;
                    }
                    else
                        *py_resolution   = DEM_resolution;
                }
            }
        }
    }
    
    //full computation
    if(info->check_full_cal)
    {
        *py_resolution   = (double)(resolution*pwrtwo(level));
        *grid_resolution = *py_resolution;
    }
    
    if(*py_resolution < DEM_resolution)
        *py_resolution = DEM_resolution;
    
    *grid_resolution = *py_resolution;
    
    
    Size_Grid2D->width  = (int)(ceil((double)(subBoundary[2] - subBoundary[0])/(*grid_resolution) ));
    Size_Grid2D->height = (int)(ceil((double)(subBoundary[3] - subBoundary[1])/(*grid_resolution) ));
 
    printf("DEM resolution %f\tresolution %f\t size %d\t%d\n",DEM_resolution,*py_resolution,Size_Grid2D->width,Size_Grid2D->height);
    D2DPOINT *GridPT = SetDEMGrid(subBoundary, *grid_resolution, *grid_resolution, Size_Grid2D);

    return GridPT;
}

UGRID *SetGrid3PT(const ProInfo *proinfo, LevelInfo &rlevelinfo, const double Th_roh, double *minmaxHeight)
{
    UGRID *GridPT3 = NULL;
    long int total_grid_counts = *rlevelinfo.Grid_length;;

    GridPT3                 = (UGRID*)calloc(sizeof(UGRID),total_grid_counts);
    
#pragma omp parallel for
    for(long i=0;i<total_grid_counts;i++)
    {
        GridPT3[i].Matched_flag     = 0;
        GridPT3[i].roh              = DoubleToSignedChar_grid(Th_roh);
        GridPT3[i].anchor_flag      = 0;
        GridPT3[i].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
        GridPT3[i].Mean_ortho_ncc   = 0;

        GridPT3[i].minHeight        = floor(minmaxHeight[0] - 0.5);
        GridPT3[i].maxHeight        = ceil(minmaxHeight[1] + 0.5);
        GridPT3[i].Height           = -1000.0;
        GridPT3[i].selected_pair    = 100;
        GridPT3[i].ncc_seleceted_pair    = -1;
        GridPT3[i].total_images     = 0;
    }
    if(proinfo->pre_DEMtif)
    {
        printf("seedem load\n");
        SetHeightWithSeedDEM(proinfo, rlevelinfo, GridPT3, minmaxHeight);
    }

    return GridPT3;
}

void SetSubBoundary(const double *Boundary, const double subX, const double subY, const double buffer_area, const int col, const int row, double *subBoundary)
{
    subBoundary[0]      = Boundary[0] + subX*(col - 1) - buffer_area;
    subBoundary[1]      = Boundary[1] + subY*(row - 1) - buffer_area;
    subBoundary[2]      = Boundary[0] + subX*(col    ) + buffer_area;
    subBoundary[3]      = Boundary[1] + subY*(row    ) + buffer_area;
    
    subBoundary[0] =  (int)(floor(subBoundary[0]/8))*8 - 40;
    subBoundary[1] =  (int)(floor(subBoundary[1]/8))*8 - 40;
    subBoundary[2] =  (int)(floor(subBoundary[2]/8))*8 + 40;
    subBoundary[3] =  (int)(floor(subBoundary[3]/8))*8 + 40;
    
    if(subBoundary[0] < Boundary[0])
        subBoundary[0] = Boundary[0];
    if(subBoundary[1] < Boundary[1])
        subBoundary[1] = Boundary[1];
    if(subBoundary[2] > Boundary[2])
        subBoundary[2] = Boundary[2];
    if(subBoundary[3] > Boundary[3])
        subBoundary[3] = Boundary[3];
}

void SetHeightWithSeedDEM(const ProInfo *proinfo, LevelInfo &rlevelinfo, UGRID *Grid, double *minmaxHeight)
{
    const char *GIMP_path = proinfo->priori_DEM_tif;
    int IsRA = proinfo->IsRA;
    const char* metafilename = proinfo->metafilename;

    int check_ftype = 1; // 1 = tif, 2 = raw
    const char *ext = strrchr(GIMP_path,'.');

    if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
    {
        check_ftype = 1;
    }
    else if(!strcmp("raw",ext+1))
    {
        check_ftype = 2;
    }

    const double oriminH = minmaxHeight[0];
    const double orimaxH = minmaxHeight[1];

    printf("oriminmaxH %f\t%f\n",oriminH,orimaxH);

    double seedDEM_sigma = proinfo->seedDEMsigma;
    if(!proinfo->check_Matchtag)
    {
        if(seedDEM_sigma < 10)
            seedDEM_sigma = 10;

        if(IsRA == 1)
        {
            if(seedDEM_sigma < 50)
                seedDEM_sigma = 50;
        }
    }

    printf("ttt1 %f\n",seedDEM_sigma);

    FILE *pFile_meta  = fopen(metafilename,"r");
    printf("meta file = %s\n",metafilename);


    double minX = 0, maxX = 0, minY = 0, maxY = 0, grid_size = 0;
    CSize seeddem_size;
    if(pFile_meta)
    {
        char bufstr[500];
        printf("open Boundary\n");
        while(!feof(pFile_meta))
        {
            fgets(bufstr,500,pFile_meta);
            if (strstr(bufstr,"Output Resolution=")!=NULL)
            {
                printf("%s\n",bufstr);
                sscanf(bufstr,"Output Resolution=%lf\n",&grid_size);
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
        fclose(pFile_meta);
    }
    if(check_ftype == 1)
    {
        seeddem_size = ReadGeotiff_info(GIMP_path, &minX, &maxY, &grid_size);
    }
    else if(check_ftype == 2)
    {
        char* hdr_path = remove_ext(GIMP_path);
        sprintf(hdr_path,"%s.hdr",hdr_path);

        printf("hdr path %s\n",hdr_path);
        seeddem_size  = Envihdr_reader_seedDEM(*rlevelinfo.param,hdr_path, &minX, &maxY, &grid_size);
        free(hdr_path);
    }

    maxX    = minX + grid_size*((double)seeddem_size.width);
    minY    = maxY - grid_size*((double)seeddem_size.height);

    printf("%d\n",seeddem_size.width);
    printf("%d\n",seeddem_size.height);
    printf("%f\n",minX);
    printf("%f\n",minY);
    printf("%f\n",maxX);
    printf("%f\n",maxY);
    printf("%f\n",grid_size);

    double a_minX = 0, a_maxX = 0, a_minY = 0, a_maxY = 0;
    if(minX > rlevelinfo.Boundary[0])
        a_minX = minX;
    else 
        a_minX = rlevelinfo.Boundary[0];

    if (maxX < rlevelinfo.Boundary[2])
        a_maxX = maxX;
    else 
        a_maxX = rlevelinfo.Boundary[2];

    if(minY > rlevelinfo.Boundary[1])
        a_minY = minY;
    else 
        a_minY = rlevelinfo.Boundary[1];

    if (maxY < rlevelinfo.Boundary[3])
        a_maxY = maxY;
    else 
        a_maxY = rlevelinfo.Boundary[3];

    double total_minH = 999999;
    double total_maxH = -999999;
    printf("%f %f %f %f\n",rlevelinfo.Boundary[0],rlevelinfo.Boundary[1],rlevelinfo.Boundary[2],rlevelinfo.Boundary[3]);
    printf("%f %f %f %f\n",a_minX, a_maxX, a_minY, a_maxY);
    if ( (a_minX < a_maxX) && (a_minY < a_maxY))
    {
        printf("seeddem cal %f\n",seedDEM_sigma);

        if(check_ftype == 2)
        {
            FILE *bin = fopen(GIMP_path,"rb");
            float *seeddem = (float*)malloc(sizeof(float)*(long int)seeddem_size.width*(long int)seeddem_size.height);
            fread(seeddem,sizeof(float),(long int)seeddem_size.width*(long int)seeddem_size.height,bin);

            SetGridHeightFromSeed(rlevelinfo, Grid, seeddem, seeddem_size, grid_size, minX, maxY, seedDEM_sigma, minmaxHeight);
            
            free(seeddem);
        }
        else
        {
            long int cols[2];
            long int rows[2];
            CSize data_size;

            CSize *LImagesize = (CSize*)malloc(sizeof(CSize));
            LImagesize->width = seeddem_size.width;
            LImagesize->height = seeddem_size.height;

            cols[0] = 0;
            cols[1] = seeddem_size.width;

            rows[0] = 0;
            rows[1] = seeddem_size.height;

            float type(0);
            float *seeddem = Readtiff_T(GIMP_path,LImagesize,cols,rows,&data_size, type);
            printf("Grid size %d\t%d\tcols rows %d\t%d\t%d\t%d\n",rlevelinfo.Size_Grid2D->width,rlevelinfo.Size_Grid2D->height,cols[0],cols[1],rows[0],rows[1]);

            SetGridHeightFromSeed(rlevelinfo, Grid, seeddem, seeddem_size, grid_size, minX, maxY, seedDEM_sigma, minmaxHeight);
            
            free(seeddem);
        }
        printf("seeddem end\n");
        printf("%f %f\n",minmaxHeight[0],minmaxHeight[1]);
    }
}

void SetGridHeightFromSeed(LevelInfo &rlevelinfo, UGRID *Grid, float *seeddem, CSize seeddem_size, double seed_grid, double minX, double maxY, double seedDEM_sigma, double *minmaxHeight)
{
    double total_minH = 999999;
    double total_maxH = -999999;
    
    for (long row = 0; row < rlevelinfo.Size_Grid2D->height; row ++)
    {
        for (long col = 0; col < rlevelinfo.Size_Grid2D->width; col ++)
        {
            long int index_grid = row*(long int)rlevelinfo.Size_Grid2D->width + col;
            double t_x = rlevelinfo.Boundary[0] + col*(*rlevelinfo.grid_resolution);
            double t_y = rlevelinfo.Boundary[1] + row*(*rlevelinfo.grid_resolution);
            
            long int col_seed = floor((t_x - minX)/seed_grid);
            long int row_seed = floor((maxY - t_y)/seed_grid);
            long int index_seeddem = row_seed*(long)seeddem_size.width + col_seed;
            if(index_seeddem >= 0 && index_seeddem < (long int)seeddem_size.width*(long int)seeddem_size.height)
            {
                if(seeddem[index_seeddem] > -1000)
                {
                    if(seeddem[index_seeddem] >= minmaxHeight[0] && seeddem[index_seeddem] <= minmaxHeight[1])
                    {
                        if(minmaxHeight[0] > (int)(seeddem[index_seeddem] - seedDEM_sigma - 0.5))
                            Grid[index_grid].minHeight = floor(minmaxHeight[0]);
                        else
                            Grid[index_grid].minHeight = floor(seeddem[index_seeddem] - seedDEM_sigma - 0.5);

                        if(minmaxHeight[1] < (int)(seeddem[index_seeddem] + seedDEM_sigma + 0.5))
                            Grid[index_grid].maxHeight = ceil(minmaxHeight[1]);
                        else
                            Grid[index_grid].maxHeight = ceil(seeddem[index_seeddem] + seedDEM_sigma + 0.5);

                        Grid[index_grid].Height        = seeddem[index_seeddem];
                    }
                    else
                    {
                        Grid[index_grid].minHeight      = floor(minmaxHeight[0] - 0.5);
                        Grid[index_grid].maxHeight      = ceil(minmaxHeight[1] + 0.5);
                        Grid[index_grid].Height         = -1000.0;
                    }

                    if(seeddem[index_seeddem] >= minmaxHeight[0] && seeddem[index_seeddem] <= minmaxHeight[1])
                    {
                        if(total_minH > seeddem[index_seeddem] -seedDEM_sigma)
                            total_minH = seeddem[index_seeddem] -seedDEM_sigma;

                        if(total_maxH < seeddem[index_seeddem] +seedDEM_sigma)
                            total_maxH = seeddem[index_seeddem] +seedDEM_sigma;
                    }
                }
                else
                {
                    Grid[index_grid].minHeight = Nodata;
                    Grid[index_grid].maxHeight = Nodata;
                }
            }
        }
    }
    minmaxHeight[0] = total_minH;
    minmaxHeight[1] = total_maxH;
}

void SetDEMBoundary(const ProInfo *info, double** _rpcs, double* _res,TransParam _param, double* _boundary, double* _minmaxheight, double* _Hinterval)
{
    double minLon = (double) (-1.2 * _rpcs[1][2] + _rpcs[0][2]);
    double maxLon = (double) (1.2 * _rpcs[1][2] + _rpcs[0][2]);
    double minLat = (double) (-1.2 * _rpcs[1][3] + _rpcs[0][3]);
    double maxLat = (double) (1.2 * _rpcs[1][3] + _rpcs[0][3]);
    
    if(info->sensor_provider == DG)
    {
        _minmaxheight[0] =  floor((-1.5 * _rpcs[1][4] + _rpcs[0][4])/10.0)*10;
        _minmaxheight[1] =  ceil((1.5 * _rpcs[1][4] + _rpcs[0][4])/10.0)*10;
    }
    else if(info->sensor_provider == PT)
    {
        printf("planet height\n");
        _minmaxheight[0] =  floor((-0.5 * _rpcs[1][4] + _rpcs[0][4])/10.0)*10;
        _minmaxheight[1] =  ceil((0.5 * _rpcs[1][4] + _rpcs[0][4])/10.0)*10;
    }
    
    printf("minmaxheight %f\t%f\n",_minmaxheight[0],_minmaxheight[1]);
    int oriminmaxH[2];
    oriminmaxH[0] = floor(-1.0 * _rpcs[1][4] + _rpcs[0][4]);
    oriminmaxH[1] = ceil(1.0 * _rpcs[1][4] + _rpcs[0][4]);
    
    
    if (oriminmaxH[0] < -100) {
        oriminmaxH[0] = -100;
    }
    
    _Hinterval[0] = _minmaxheight[1] - _minmaxheight[0];
    
    if (_minmaxheight[0] < -100) {
        _minmaxheight[0] = -100;
    }
    
    D2DPOINT lonlat[4];
    lonlat[0].m_X = minLon;
    lonlat[0].m_Y = minLat;
    lonlat[1].m_X = minLon;
    lonlat[1].m_Y = maxLat;
    lonlat[2].m_X = maxLon;
    lonlat[2].m_Y = maxLat;
    lonlat[3].m_X = maxLon;
    lonlat[3].m_Y = minLat;
    
    double minX = min(lonlat[3].m_X, min(lonlat[2].m_X, min(lonlat[0].m_X, lonlat[1].m_X)));
    double maxX = max(lonlat[3].m_X, max(lonlat[2].m_X, max(lonlat[0].m_X, lonlat[1].m_X)));
    
    double minY = min(lonlat[3].m_Y, min(lonlat[2].m_Y, min(lonlat[0].m_Y, lonlat[1].m_Y)));
    double maxY = max(lonlat[3].m_Y, max(lonlat[2].m_Y, max(lonlat[0].m_Y, lonlat[1].m_Y)));
    
    printf("lonlat X %f\t%f\t%f\t%f\n",lonlat[0].m_X,lonlat[1].m_X,lonlat[2].m_X,lonlat[3].m_X);
    printf("lonlat Y %f\t%f\t%f\t%f\n",lonlat[0].m_Y,lonlat[1].m_Y,lonlat[2].m_Y,lonlat[3].m_Y);
    
    printf("minmaxXY %f\t%f\t%f\t%f\n",minX,minY,maxX,maxY);
    _boundary[0] =  (minX);
    _boundary[1] =  (minY);
    _boundary[2] =  (maxX);
    _boundary[3] =  (maxY);
    
    //_imagesize->height = (unsigned int) (ceil((_boundary[3] - _boundary[1]) / _res[1]));
    //_imagesize->width = (unsigned int) (ceil((_boundary[2] - _boundary[0]) / _res[0]));
}

void SetDEMBoundary_photo(EO Photo, CAMERA_INFO m_Camera, RM M, double* _boundary, double* _minmaxheight, double* _Hinterval)
{
    // manual setup of minmaxheight for test
    double MSL = 0;
    _minmaxheight[0] =  0;
    _minmaxheight[1] =  200;
    
    int oriminmaxH[2];
    oriminmaxH[0] = 0;
    oriminmaxH[1] = 150;
    
    
    if (oriminmaxH[0] < -100) {
        oriminmaxH[0] = -100;
    }
    
    _Hinterval[0] = _minmaxheight[1] - _minmaxheight[0];
    
    if (_minmaxheight[0] < -100) {
        _minmaxheight[0] = -100;
    }
    
    D3DPOINT top_left_3D,top_right_3D,bottom_right_3D,bottom_left_3D;
    D2DPOINT top_left, top_right, bottom_right,bottom_left;
    
    top_left.m_X = 0.0;
    top_left.m_Y = 0.0;
    top_right.m_X = m_Camera.m_ImageSize.width;
    top_right.m_Y = 0.0;
    bottom_right.m_X = m_Camera.m_ImageSize.width;
    bottom_right.m_Y = m_Camera.m_ImageSize.height;
    bottom_left.m_X = 0.0;
    bottom_left.m_Y = m_Camera.m_ImageSize.height;
    
    top_left = ImageToPhoto_single(top_left,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    top_right = ImageToPhoto_single(top_right,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    bottom_right = ImageToPhoto_single(bottom_right,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    bottom_left = ImageToPhoto_single(bottom_left,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    
    top_left_3D     = GetObjectCoordinate_single(top_left,MSL,Photo, m_Camera, M);
    top_right_3D    = GetObjectCoordinate_single(top_right,MSL,Photo, m_Camera, M);
    bottom_right_3D = GetObjectCoordinate_single(bottom_right,MSL,Photo, m_Camera, M);
    bottom_left_3D  = GetObjectCoordinate_single(bottom_left,MSL,Photo, m_Camera, M);
    
    double minX = (top_left_3D.m_X < bottom_left_3D.m_X) ? top_left_3D.m_X : bottom_left_3D.m_X;
    double minY = (bottom_left_3D.m_Y < bottom_right_3D.m_Y) ? bottom_left_3D.m_Y : bottom_right_3D.m_Y;
    double maxX = (top_right_3D.m_X > bottom_right_3D.m_X) ? top_right_3D.m_X : bottom_right_3D.m_X;
    double maxY = (top_left_3D.m_Y > top_right_3D.m_Y) ? top_left_3D.m_Y : top_right_3D.m_Y;
    
    _boundary[0] =  floor(minX);
    _boundary[1] =  floor(minY);
    _boundary[2] =  ceil(maxX);
    _boundary[3] =  ceil(maxY);
}

uint16 *SetsubsetImage(ProInfo *proinfo, LevelInfo &rlevelinfo, const int index_image, const TransParam transparam, const uint8 NumofIAparam, const double * const * const *RPCs, const double * const * ImageParams, const double *subBoundary, const double *minmaxHeight, D2DPOINT *Startpos, CSize *Subsetsize)
{
    bool ret = false;

    CSize Imagesize;
    uint16 *outimage = NULL;
    
    if(GetImageSize(proinfo->Imagefilename[index_image],&Imagesize))
    {
        long int Lcols[2], Lrows[2];
        if(GetsubareaImage(proinfo->sensor_type, proinfo->frameinfo, index_image, *rlevelinfo.param, rlevelinfo.ImageAdjust[index_image],  rlevelinfo.RPCs[index_image], proinfo->Imagefilename[index_image], Imagesize, rlevelinfo.Boundary, minmaxHeight, Lcols, Lrows))
        {
            printf("read image %d\n", index_image);
            uint16 type(0);
            outimage   = Readtiff_T(proinfo->Imagefilename[index_image],&Imagesize,Lcols,Lrows,&Subsetsize[index_image],type);
            if(proinfo->check_checktiff)
                exit(1);
            
            Startpos[index_image].m_X  = (double)(Lcols[0]);
            Startpos[index_image].m_Y  = (double)(Lrows[0]);
            
            if(!outimage)
                proinfo->check_selected_image[index_image] = false;
        }
    }
    else
        proinfo->check_selected_image[index_image] = false;
  
    return outimage;
}

// temporary, 1st and 2nd image
void CalMPP_pair(double CA,double mean_product_res, double im_resolution, double *MPP_stereo_angle)
{
    double ccdsize = 0.00001;
    double scale = ccdsize/mean_product_res;
    double convergence_mpp = 1.0/tan(CA*DegToRad*0.5)*mean_product_res;
    double BH_ratio = mean_product_res*2/convergence_mpp;
    double sigmaZ = 1.414*ccdsize/BH_ratio/scale;
    
    *MPP_stereo_angle = sigmaZ;
    
    if(*MPP_stereo_angle < im_resolution*1.5)
        *MPP_stereo_angle = im_resolution*1.5;
    
    printf("mpp = %f\n",*MPP_stereo_angle);
}

void CalMPP(ProInfo *proinfo, LevelInfo &rlevelinfo, const double* minmaxHeight, double CA,const double mean_product_res, double *MPP_simgle_image, double *MPP_stereo_angle, const int pair_number)
{
    const int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
    const int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
    
    D2DPOINT temp_p1, temp_p2;
    double temp_LIA[2] = {0.0, 0.0};
    
    const long int numofpts = *rlevelinfo.Grid_length;
    D3DPOINT temp_GrP(rlevelinfo.Grid_wgs[(int)(numofpts/2.0)].m_X, rlevelinfo.Grid_wgs[(int)(numofpts/2.0)].m_Y, minmaxHeight[0], 0);
    temp_p1     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[reference_id],*rlevelinfo.NumOfIAparam,temp_LIA,temp_GrP);
    
    temp_GrP.m_Z = minmaxHeight[1];
    temp_p2     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[reference_id],*rlevelinfo.NumOfIAparam,temp_LIA,temp_GrP);
    
    const double left_mpp = (minmaxHeight[1] - minmaxHeight[0]) / sqrt( pow(temp_p1.m_X - temp_p2.m_X,2.0) + pow(temp_p1.m_Y - temp_p2.m_Y,2.0));
    
    temp_GrP.m_Z = minmaxHeight[0];
    temp_GrP.flag = 0;
    temp_p1     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[ti],*rlevelinfo.NumOfIAparam,rlevelinfo.ImageAdjust[pair_number],temp_GrP);
    
    temp_GrP.m_Z = minmaxHeight[1];
    temp_p2     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[ti],*rlevelinfo.NumOfIAparam,rlevelinfo.ImageAdjust[pair_number],temp_GrP);
    
    const double right_mpp = (minmaxHeight[1] - minmaxHeight[0]) / sqrt( pow(temp_p1.m_X - temp_p2.m_X,2.0) + pow(temp_p1.m_Y - temp_p2.m_Y,2.0));
    
    printf("left right mpp %f\t%f\n",left_mpp,right_mpp);
    
    if(left_mpp > 5*proinfo->resolution)
        *MPP_simgle_image = right_mpp;
    else if(right_mpp > 5*proinfo->resolution)
        *MPP_simgle_image = left_mpp;
    else
    {
        if(left_mpp > right_mpp)
            *MPP_simgle_image = left_mpp;
        else
            *MPP_simgle_image = right_mpp;
    }
    
    CA = rlevelinfo.pairinfo->ConvergenceAngle[pair_number];
    if(proinfo->sensor_provider == DG)
    {
        const double ccdsize = 0.00001;
        const double scale = ccdsize/mean_product_res;
        const double convergence_mpp = 1.0/tan(CA*DegToRad*0.5)*mean_product_res;
        const double BH_ratio = rlevelinfo.pairinfo->BHratio[pair_number];// mean_product_res*2/convergence_mpp;
        const double sigmaZ = 1.414*ccdsize/BH_ratio/scale;
        
        printf("scale = %f\tConvergence_Angle = %f\tconvergnece_mpp = %f\tBH_ratio = %f\t sigmaZ = %f\n",scale,CA,convergence_mpp,BH_ratio,sigmaZ);
        if(*rlevelinfo.Pyramid_step > 2)
            *MPP_stereo_angle = (*MPP_simgle_image)*sigmaZ;
        else
        {
            *MPP_stereo_angle = sigmaZ;
            *MPP_simgle_image = sigmaZ;
        }
    }
    else if(proinfo->sensor_provider == PT)
    {
        if(*rlevelinfo.Pyramid_step > 2)
            *MPP_stereo_angle = (*MPP_simgle_image)*proinfo->resolution*1.5;
        else
        {
            *MPP_stereo_angle = proinfo->resolution*1.5;
            *MPP_simgle_image = proinfo->resolution*1.5;
        }
    }
        
    if(*MPP_stereo_angle < proinfo->resolution*1.5)
        *MPP_stereo_angle = proinfo->resolution*1.5;
    
    printf("pair_number = %d\tmpp = %f\t mpr = %f\n",pair_number,*MPP_simgle_image,*MPP_stereo_angle);
}

void CalMPP_8(ProInfo *proinfo, LevelInfo &rlevelinfo, const double* minmaxHeight, double CA,const double mean_product_res, double *MPP_simgle_image, double *MPP_stereo_angle, const int pair_number)
{
    const int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
    const int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
    
    D2DPOINT temp_p1, temp_p2;
    double temp_LIA[2] = {0.0, 0.0};
    
    const long int numofpts = *rlevelinfo.Grid_length;
    D3DPOINT temp_GrP(rlevelinfo.Grid_wgs[(int)(numofpts/2.0)].m_X, rlevelinfo.Grid_wgs[(int)(numofpts/2.0)].m_Y, minmaxHeight[0], 0);
    temp_p1     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[reference_id],*rlevelinfo.NumOfIAparam,temp_LIA,temp_GrP);
    
    temp_GrP.m_Z = minmaxHeight[1];
    temp_p2     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[reference_id],*rlevelinfo.NumOfIAparam,temp_LIA,temp_GrP);
    
    const double left_mpp = (minmaxHeight[1] - minmaxHeight[0]) / sqrt( pow(temp_p1.m_X - temp_p2.m_X,2.0) + pow(temp_p1.m_Y - temp_p2.m_Y,2.0));
    
    temp_GrP.m_Z = minmaxHeight[0];
    temp_GrP.flag = 0;
    temp_p1     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[ti],*rlevelinfo.NumOfIAparam,rlevelinfo.ImageAdjust[pair_number],temp_GrP);
    
    temp_GrP.m_Z = minmaxHeight[1];
    temp_p2     = GetObjectToImageRPC_single_mpp(rlevelinfo.RPCs[ti],*rlevelinfo.NumOfIAparam,rlevelinfo.ImageAdjust[pair_number],temp_GrP);
    
    const double right_mpp = (minmaxHeight[1] - minmaxHeight[0]) / sqrt( pow(temp_p1.m_X - temp_p2.m_X,2.0) + pow(temp_p1.m_Y - temp_p2.m_Y,2.0));
    
    printf("left right mpp %f\t%f\n",left_mpp,right_mpp);
    
    if(left_mpp > 5*proinfo->resolution)
        *MPP_simgle_image = right_mpp;
    else if(right_mpp > 5*proinfo->resolution)
        *MPP_simgle_image = left_mpp;
    else
    {
        if(left_mpp > right_mpp)
            *MPP_simgle_image = left_mpp;
        else
            *MPP_simgle_image = right_mpp;
    }
    
    CA = rlevelinfo.pairinfo->ConvergenceAngle[pair_number];
    if(proinfo->sensor_provider == DG)
    {
        const double ccdsize = 0.00001;
        const double scale = ccdsize/mean_product_res;
        const double convergence_mpp = 1.0/tan(CA*DegToRad*0.5)*mean_product_res;
        const double BH_ratio = rlevelinfo.pairinfo->BHratio[pair_number];//mean_product_res*2.0/convergence_mpp;
        const double sigmaZ = 1.414*ccdsize/BH_ratio/scale;
        
        printf("scale = %f\tConvergence_Angle = %f\tconvergnece_mpp = %f\tBH_ratio = %f\t sigmaZ = %f\n",scale,CA,convergence_mpp,BH_ratio,sigmaZ);
        if(sigmaZ > 1)
            *MPP_stereo_angle = (*MPP_simgle_image)*sigmaZ;
    }
    else if(proinfo->sensor_provider == PT)
    {
        if(*rlevelinfo.Pyramid_step > 2)
            *MPP_stereo_angle = (*MPP_simgle_image)*proinfo->resolution*2.0;
        else
        {
            *MPP_stereo_angle = proinfo->resolution*2.0;
            *MPP_simgle_image = proinfo->resolution*2.0;
        }
    }
    
    if(*MPP_stereo_angle < proinfo->resolution*2.0)
        *MPP_stereo_angle = proinfo->resolution*2.0;
    
    printf("pair_number = %d\tmpp = %f\t mpr = %f\n",pair_number,*MPP_simgle_image,*MPP_stereo_angle);
}

double GetHeightStep(int Pyramid_step, double im_resolution)
{
    const double h_divide = 2;
    
    im_resolution = im_resolution*pwrtwo(Pyramid_step);
    
    double HS = (double)(im_resolution/h_divide);
    
    double &&tt1 = HS*1000.0;
    double &&tt2 = floor(tt1 + 0.1);
    HS = tt2/1000.0;
    
    if(HS > 3)
        HS = 3;
    
    return HS;
}

double GetHeightStep_Planet(int Pyramid_step, double im_resolution)
{
    const double h_divide = 4;
    
    im_resolution = im_resolution*pwrtwo(Pyramid_step);
    
    double HS = (double)(im_resolution/h_divide);
    
    double &&tt1 = HS*1000.0;
    double &&tt2 = floor(tt1 + 0.1);
    HS = tt2/1000.0;
    
    if(HS > 3)
        HS = 3;
    
    return HS;
}

void InitializeVoxel(const ProInfo *proinfo, GridVoxel &grid_voxel,LevelInfo &plevelinfo, UGRID *GridPT3, NCCresult* nccresult,const int iteration, const double *minmaxHeight)
{
    const double height_step = *plevelinfo.height_step;
    const uint8 pyramid_step = *plevelinfo.Pyramid_step;
    
#pragma omp parallel for schedule(guided)
    for(long int t_i = 0 ; t_i < *plevelinfo.Grid_length; t_i++)
    {
        int selected_images;
        if(check_image_boundary(proinfo,plevelinfo,plevelinfo.GridPts[t_i],plevelinfo.Grid_wgs[t_i],GridPT3[t_i].minHeight,GridPT3[t_i].maxHeight,7,selected_images))
        {
            //if(check_image_boundary_each(proinfo,plevelinfo,plevelinfo.GridPts[t_i],plevelinfo.Grid_wgs[t_i],GridPT3[t_i].minHeight,GridPT3[t_i].maxHeight,7,selected_images))
            //    GridPT3[t_i].total_images = selected_images;
            
            int change_step_min = 0;
            int change_step_max = 0;
            bool check_blunder_cell = true;
            double th_height = 1000;
            if(proinfo->DEM_resolution <= 4)
                th_height = 500;
            
            if(proinfo->check_Matchtag || iteration == 1)
            {
                nccresult[t_i].minHeight = GridPT3[t_i].minHeight;
                nccresult[t_i].maxHeight = GridPT3[t_i].maxHeight;
                //nccresult[t_i].GNCC = DoubleToSignedChar_result(-1.0);
                nccresult[t_i].check_height_change = true;
                check_blunder_cell = false;
            }
            else
            {
                if ( pyramid_step >= 2 )
                    check_blunder_cell = false;
                else if( GridPT3[t_i].Matched_flag != 0)
                {
                    if(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight > 0)
                    {
                        if(pyramid_step <= 1)
                        {
                            if(iteration > 1)
                            {
                                if(GridPT3[t_i].maxHeight <= minmaxHeight[1] && GridPT3[t_i].minHeight >= minmaxHeight[0] && nccresult[t_i].maxHeight <= minmaxHeight[1] && nccresult[t_i].minHeight >= minmaxHeight[0])
                                    check_blunder_cell = false;
                                else
                                    check_blunder_cell = true;
                            }
                            else
                            {
                                if(GridPT3[t_i].maxHeight <= minmaxHeight[1] && GridPT3[t_i].minHeight >= minmaxHeight[0])
                                    check_blunder_cell = false;
                                else
                                    check_blunder_cell = true;
                            }
                        }
                        
                        if(pyramid_step == 1)
                        {
                            if((abs(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight) < 1000))
                                check_blunder_cell = false;
                            else
                                check_blunder_cell = true;
                        }
                        
                        if(pyramid_step == 0)
                        {
                            if(iteration > 1)
                            {
                                if((abs(nccresult[t_i].maxHeight - nccresult[t_i].minHeight) < th_height) && (abs(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight) < th_height))
                                    check_blunder_cell = false;
                                else
                                    check_blunder_cell = true;
                            }
                            else
                            {
                                if((abs(GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight) < th_height))
                                    check_blunder_cell = false;
                                else
                                    nccresult[t_i].NumOfHeight = 0;
                            }
                        }
                    }
                    else
                    {
                        if(nccresult[t_i].NumOfHeight > 0)
                        {
                            grid_voxel[t_i].clear();
                        }
                        
                        nccresult[t_i].NumOfHeight = 0;
                        check_blunder_cell = true;
                    }
                }
                else
                {
                    if(nccresult[t_i].NumOfHeight > 0)
                    {
                        grid_voxel[t_i].clear();
                    }
                    
                    nccresult[t_i].NumOfHeight = 0;
                    check_blunder_cell = true;
                }
            }
            
            if(!check_blunder_cell)
            {
                
                if(iteration == 1 || (pyramid_step == 0 && proinfo->DEM_resolution < 2))
                {
                    nccresult[t_i].minHeight = GridPT3[t_i].minHeight;
                    nccresult[t_i].maxHeight = GridPT3[t_i].maxHeight;
                    //nccresult[t_i].GNCC = DoubleToSignedChar_result(-1.0);
                    nccresult[t_i].check_height_change = true;
                    
                    if(pyramid_step == 1)
                    {
                        if((abs((int)GridPT3[t_i].maxHeight - (int)GridPT3[t_i].minHeight) > 1000))
                        {
                            GridPT3[t_i].maxHeight = -100;
                            GridPT3[t_i].minHeight = -100;
                            nccresult[t_i].check_height_change = false;
                        }
                    }
                    else if(pyramid_step == 0)
                    {
                        if((abs((int)GridPT3[t_i].maxHeight - (int)GridPT3[t_i].minHeight) > th_height))
                        {
                            GridPT3[t_i].maxHeight = -100;
                            GridPT3[t_i].minHeight = -100;
                            nccresult[t_i].check_height_change = false;
                        }
                    }
                }
                else
                {
                    const int NumberofHeightVoxel = (int)((GridPT3[t_i].maxHeight - GridPT3[t_i].minHeight)/height_step);
                    //new search height
                    if(nccresult[t_i].NumOfHeight == 0 && NumberofHeightVoxel > 0)
                    {
                        change_step_min = 0;
                        change_step_max = 0;
                        nccresult[t_i].check_height_change = true;
                    }
                    else
                    {
                        //extension search height than GridPT3 (compare before and after search height)
                        if(nccresult[t_i].minHeight > GridPT3[t_i].minHeight || nccresult[t_i].maxHeight < GridPT3[t_i].maxHeight)
                        {
                            if(nccresult[t_i].minHeight > GridPT3[t_i].minHeight)
                                change_step_min = (int)((nccresult[t_i].minHeight - GridPT3[t_i].minHeight)/height_step + 0.5);
                            
                            if(nccresult[t_i].maxHeight < GridPT3[t_i].maxHeight)
                                change_step_max = (int)((GridPT3[t_i].maxHeight - nccresult[t_i].maxHeight)/height_step + 0.5);
                            
                            nccresult[t_i].minHeight = floor(nccresult[t_i].minHeight - change_step_min*height_step);
                            nccresult[t_i].maxHeight = ceil(nccresult[t_i].maxHeight + change_step_max*height_step);
                            
                            nccresult[t_i].check_height_change = true;
                            
                            if(pyramid_step == 0 && abs(nccresult[t_i].maxHeight - nccresult[t_i].minHeight) > th_height)
                                nccresult[t_i].check_height_change = false;
                        }
                        else
                        {
                            nccresult[t_i].check_height_change = false;
                        }
                    }
                }
                
                if(nccresult[t_i].check_height_change)
                {
                    nccresult[t_i].minHeight = floor(nccresult[t_i].minHeight - change_step_min*height_step);
                    nccresult[t_i].maxHeight = ceil(nccresult[t_i].maxHeight + change_step_max*height_step);
                    
                    const int NumberofHeightVoxel = (int)((float)(nccresult[t_i].maxHeight - nccresult[t_i].minHeight)/height_step);
                    
                    if(NumberofHeightVoxel > 0 )
                    {
                        if(nccresult[t_i].NumOfHeight > 0)
                        {
                            grid_voxel[t_i].clear();
                        }
                        
                        nccresult[t_i].NumOfHeight = NumberofHeightVoxel;
                        
                        grid_voxel[t_i].allocate(NumberofHeightVoxel);
                        
                    }
                    else
                    {
                        if(nccresult[t_i].NumOfHeight > 0)
                        {
                            grid_voxel[t_i].clear();
                        }
                        
                        nccresult[t_i].NumOfHeight = 0;
                        nccresult[t_i].check_height_change = false;
                    }
                }
            }
            else
            {
                if(nccresult[t_i].NumOfHeight > 0)
                {
                    grid_voxel[t_i].clear();
                }
                
                nccresult[t_i].NumOfHeight = 0;
                nccresult[t_i].check_height_change = false;
            }
        }
        else
        {
            if(nccresult[t_i].NumOfHeight > 0)
            {
                grid_voxel[t_i].clear();
            }
            
            nccresult[t_i].NumOfHeight = 0;
            nccresult[t_i].check_height_change = false;
        }
        
        if((nccresult[t_i].minHeight == 0 || nccresult[t_i].maxHeight == 0))
        {
            if(nccresult[t_i].NumOfHeight > 0)
            {
                grid_voxel[t_i].clear();
            }
            
            nccresult[t_i].NumOfHeight = 0;
            nccresult[t_i].check_height_change = false;
        }
        
        if(pyramid_step == 0 && nccresult[t_i].NumOfHeight > 1000)
        {
            if(nccresult[t_i].NumOfHeight > 0)
            {
                grid_voxel[t_i].clear();
            }
            
            nccresult[t_i].NumOfHeight = 0;
            nccresult[t_i].check_height_change = false;
        }
    }

}

double SetNCC_alpha(const int Pyramid_step, const int iteration, bool IsRA)
{
    double ncc_alpha;
    
    if(Pyramid_step >= 4)
    {
        ncc_alpha = 1.0 - ((4-Pyramid_step)*0.2 + (iteration-1)*0.05);
        if(ncc_alpha < 0.8)
            ncc_alpha = 0.8;
    }
    else if(Pyramid_step >= 3)
    {
        ncc_alpha = 1.0 - ((4-Pyramid_step)*0.2 + (iteration-1)*0.05);
        if(ncc_alpha < 0.6)
            ncc_alpha = 0.6;
    }
    else if(Pyramid_step == 2)
    {
        ncc_alpha = 1.0 - (0.4 + (iteration-1)*0.05);
        if(ncc_alpha < 0.4)
            ncc_alpha = 0.4;
    }
    else
        ncc_alpha = 1.0 - ((4-Pyramid_step)*0.2 + (iteration-1)*0.05);
    
    if(ncc_alpha < 0.1)
        ncc_alpha = 0.1;
    
    if(IsRA == 1)
        ncc_alpha = 1.0;
    
    return ncc_alpha;
    
}
/*
double SetGnccWeight(int Pyramid_step, double GNCC, double INCC, double matched_height, double tar_height, double step_height)
{
    double gncc_weight = 1.0;
    
    if(INCC > 0)
    {
        if(Pyramid_step < 2)
        {
            gncc_weight = 1.0;
        }
        else
        {
            double height_diff_th = 5;
            
            double weight = 0.5;
            double height_step_diff = fabs(matched_height - tar_height)/step_height;
         
            if(INCC >= 0 && GNCC >= 0)
            {
                if(height_step_diff <= height_diff_th)
                {
                    if(height_step_diff == 0)
                        gncc_weight = 1.0;
                    else
                        gncc_weight = 1.0 - weight*(1.0 - GNCC)*(1.0/(1.0+(height_diff_th - height_step_diff)*(height_diff_th - height_step_diff)));
                }
                else
                    gncc_weight = 1.0 - weight*(1.0 - GNCC);
            }
            else
            {
                gncc_weight = 1.0 - weight;
            }
        }
    }
    else
        gncc_weight = 1.0;
    return gncc_weight;
}
*/
double SetGnccWeight(int Pyramid_step, double GNCC, double INCC, double matched_height, double tar_height, double step_height)
{
    double gncc_weight = 1.0;
    
    if(INCC > 0)
    {
        if(Pyramid_step < 2)
        {
            gncc_weight = 1.0;
        }
        else
        {
            double height_diff_th = 5;
            
            double weight = 0.2;
            double height_step_diff = fabs(matched_height - tar_height)/step_height;
         
            if(GNCC >= 0)
            {
                if(height_step_diff <= height_diff_th)
                {
                    if(height_step_diff == 0)
                        gncc_weight = 1.0 + weight*GNCC;
                    else
                        gncc_weight = 1.0 + weight*GNCC*(1.0/(1.0+(height_step_diff)*(height_step_diff)));
                }
                else
                    gncc_weight = 1.0;
            }
            else
            {
                gncc_weight = 1.0;
            }
        }
    }
    else
        gncc_weight = 1.0;
    gncc_weight = 1.0;
    
    return gncc_weight;
}

int select_referenceimage(const long pt_index, const ProInfo *proinfo, LevelInfo &plevelinfo, double start_H, double end_H)
{
    int selected_ref = 0;
    double diff_min = 1000000;
    
    double min_off_nadir = 90;
    
    for(int pair_number = 0 ; pair_number < plevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        const int reference_id = plevelinfo.pairinfo->pairs[pair_number].m_X;
        const int ti = plevelinfo.pairinfo->pairs[pair_number].m_Y;
        
         // Image point setting
        const double temp_LIA[2] = {0.0, 0.0};
        
        D2DPOINT Ref_Imagecoord[1];
        D2DPOINT Tar_Imagecoord[1];
        D3DPOINT temp_GP[1];
        D2DPOINT photo;
        temp_GP[0].m_Z = (double)(end_H - start_H)/2.0;
        if(proinfo->sensor_type == SB)
        {
            if(min_off_nadir > plevelinfo.imageinfo[reference_id].Offnadir_angle)
            {
                min_off_nadir = plevelinfo.imageinfo[reference_id].Offnadir_angle;
                selected_ref = reference_id;
            }
            if(min_off_nadir > plevelinfo.imageinfo[ti].Offnadir_angle)
            {
                min_off_nadir > plevelinfo.imageinfo[ti].Offnadir_angle;
                selected_ref = ti;
            }
            /*
            temp_GP[0] = plevelinfo.Grid_wgs[pt_index];
            
            Ref_Imagecoord[0]      = GetObjectToImageRPC_single(plevelinfo.RPCs[reference_id],*plevelinfo.NumOfIAparam,temp_LIA,temp_GP[0]);
            
            Tar_Imagecoord[0]     = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pair_number],temp_GP[0]);
             */
        }
        else
        {
            temp_GP[0] = plevelinfo.GridPts[pt_index];
            
            photo = GetPhotoCoordinate_single(temp_GP[0],proinfo->frameinfo.Photoinfo[reference_id],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[reference_id].m_Rm);
            Ref_Imagecoord[0] = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
            
            photo = GetPhotoCoordinate_single(temp_GP[0],proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[ti].m_Rm);
            Tar_Imagecoord[0] = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
        
            D2DPOINT diff;
            diff.m_X = Ref_Imagecoord[0].m_X - plevelinfo.Imagesize_ori[reference_id].width/2.0;
            diff.m_Y = Ref_Imagecoord[0].m_Y - plevelinfo.Imagesize_ori[reference_id].height/2.0;
            
            double diff_x = sqrt(diff.m_X*diff.m_X + diff.m_Y*diff.m_Y);
            diff.m_X = Tar_Imagecoord[0].m_X - plevelinfo.Imagesize_ori[ti].width/2.0;
            diff.m_Y = Tar_Imagecoord[0].m_Y - plevelinfo.Imagesize_ori[ti].height/2.0;
            double diff_y = sqrt(diff.m_X*diff.m_X + diff.m_Y*diff.m_Y);
            
            if(diff_x < diff_y)
            {
                if(diff_x < diff_min)
                {
                    diff_min = diff_x;
                    selected_ref = reference_id;
                }
            }
            else
            {
                if(diff_y < diff_min)
                {
                    diff_min = diff_y;
                    selected_ref = ti;
                }
            }
        }
    }
    
    return selected_ref;
}

int VerticalLineLocus(GridVoxel &grid_voxel,const ProInfo *proinfo, const ImageInfo *image_info, NCCresult* nccresult, LevelInfo &plevelinfo, UGRID *GridPT3, const uint8 iteration, const double *minmaxHeight)
{
    const bool check_matchtag = proinfo->check_Matchtag;
    const char* save_filepath = proinfo->save_filepath;
    const bool pre_DEMtif = proinfo->pre_DEMtif;
    const bool IsRA = proinfo->IsRA;
    long int Accessable_grid = 0;
    
    const int Pyramid_step = *(plevelinfo.Pyramid_step);
    int Template_size = *(plevelinfo.Template_size);
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(proinfo->resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
        if(*(plevelinfo.Template_size) < t_Template_size)
            Template_size = t_Template_size;
        
        printf("VerticalLineLocus : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }
    
    int Half_template_size = (int)(Template_size/2);
    const double subBoundary[4] = {plevelinfo.Boundary[0], plevelinfo.Boundary[1], plevelinfo.Boundary[2], plevelinfo.Boundary[3]};
    
    bool check_ortho = false;
    const int TH_N = 0;
    if(Pyramid_step == 3)
        Half_template_size = Half_template_size - 1;
    else if(Pyramid_step < 3)
        Half_template_size = Half_template_size - 2;
    
    const long int numofpts = *(plevelinfo.Grid_length);
    
    bool check_combined_WNCC = false;
    bool check_combined_WNCC_INCC = false;
    
    if(!IsRA && ((Pyramid_step > (*plevelinfo.Py_combined_level))))
        check_combined_WNCC = true;
    
    if(check_combined_WNCC && ( (Pyramid_step >= 1 && iteration%2 == 0) ) )
        check_combined_WNCC_INCC = true;
    
    printf("check_combined_WNCC %d\tcheck_combined_WNCC_INCC %d\tnumber of pairs %d\n",check_combined_WNCC,check_combined_WNCC_INCC,plevelinfo.pairinfo->NumberOfPairs);
    
    printf("*plevelinfo.check_matching_rate %d\n",*plevelinfo.check_matching_rate);
    
    double im_resolution_next = proinfo->resolution;
    double im_resolution = proinfo->resolution;
    if(check_combined_WNCC)
        im_resolution_next = proinfo->resolution*pwrtwo(Pyramid_step-1);
    
    im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    
    if((Pyramid_step == 4 && iteration == 1) || IsRA == true)
        check_ortho = false;
    else
        check_ortho = true;
    
    if(pre_DEMtif && !(Pyramid_step == 4 && iteration == 1))
        check_ortho = true;
    
    if(check_matchtag)
        check_ortho = true;
    
    //if(levelinfo.check_SGM)
    //    check_ortho = false;
    
    const double ncc_alpha = SetNCC_alpha(Pyramid_step,iteration, IsRA);
    const double ncc_beta = 1.0 - ncc_alpha;
        
    //int sum_data2 = 0;
    //int sum_data = 0;
    const double ortho_th = 0.7 - (4 - Pyramid_step)*0.10;
    
    
    printf("start omp\n");
#pragma omp parallel
    {
        SetKernel rsetkernel(0,1,Half_template_size);
        SetKernel rsetkernel_next(0,1,Half_template_size);
        
#pragma omp for schedule(dynamic,1) reduction(+:Accessable_grid/*,sum_data2, sum_data*/)
        for(long int iter_count = 0 ; iter_count < numofpts ; iter_count++)
        {
            long int pts_row = (long int)(floor(iter_count/plevelinfo.Size_Grid2D->width));
            long int pts_col = iter_count % plevelinfo.Size_Grid2D->width;
            long int pt_index = iter_count;//pts_row*(long int)plevelinfo.Size_Grid2D->width + pts_col;
            
            GridPT3[pt_index].total_images = 0;
            /*
            if(nccresult[pt_index].check_height_change)
                sum_data2++;
            else
                sum_data++;
            */
            const int start_H     = GridPT3[pt_index].minHeight;
            const int end_H       = GridPT3[pt_index].maxHeight;
     
            int selected_images;
            if(check_image_boundary(proinfo,plevelinfo,plevelinfo.GridPts[pt_index],plevelinfo.Grid_wgs[pt_index],start_H,end_H,Half_template_size,selected_images) &&
               pts_row >= 0 && pts_row < plevelinfo.Size_Grid2D->height && pts_col >= 0 && pts_col < plevelinfo.Size_Grid2D->width)
            {
                bool check_blunder_cell = false;
                double th_height = 1000;
                if(proinfo->DEM_resolution <= 4)
                    th_height = 500;
                
                if ( Pyramid_step >= 2)
                    check_blunder_cell = false;
                else
                {
                    if(GridPT3[pt_index].Matched_flag == 0)
                        check_blunder_cell = true;
                    else if( GridPT3[pt_index].Matched_flag != 0)
                    {
                        if(Pyramid_step == 1)
                        {
                            if((abs(GridPT3[pt_index].maxHeight - GridPT3[pt_index].minHeight) < th_height))
                                check_blunder_cell = false;
                            else
                                check_blunder_cell = true;
                        }
                        else if((abs(GridPT3[pt_index].maxHeight - GridPT3[pt_index].minHeight) < th_height))
                            check_blunder_cell = false;
                        else
                            check_blunder_cell = true;
                    }
                    else
                        check_blunder_cell = true;
                }
                
                if(Pyramid_step <= 1)
                {
                    if(GridPT3[pt_index].maxHeight > minmaxHeight[1] || GridPT3[pt_index].minHeight < minmaxHeight[0])
                        check_blunder_cell = true;
                }
                
                if(check_matchtag && iteration  <= 2)
                {
                    check_blunder_cell = false;
                    if(GridPT3[pt_index].maxHeight > minmaxHeight[1] || GridPT3[pt_index].minHeight < minmaxHeight[0])
                        check_blunder_cell = true;
                    else if( (abs(GridPT3[pt_index].maxHeight - GridPT3[pt_index].minHeight) > th_height) )
                        check_blunder_cell = true;
                }
                
                const int NumOfHeights = (int)((end_H -  start_H)/(*plevelinfo.height_step));
                if(*plevelinfo.check_matching_rate)
                {
                    nccresult[pt_index].check_height_change = true;
                    nccresult[pt_index].NumOfHeight = NumOfHeights;
                    nccresult[pt_index].minHeight = start_H;
                    nccresult[pt_index].maxHeight = end_H;
                }
                
                if(!check_blunder_cell && nccresult[pt_index].NumOfHeight > 0)
                {
                    vector<double> *WNCC_save = NULL;
                    vector<unsigned char> *WNCC_save_pair_ID = NULL;
                    
                    if(*plevelinfo.check_matching_rate)
                    {
                        WNCC_save = (vector<double>*)calloc(sizeof(vector<double>),nccresult[pt_index].NumOfHeight);
                        WNCC_save_pair_ID = (vector<unsigned char>*)calloc(sizeof(vector<unsigned char>),nccresult[pt_index].NumOfHeight);
                    }

                    GridPT3[pt_index].selected_pair = select_referenceimage(pt_index, proinfo, plevelinfo, nccresult[pt_index].minHeight, nccresult[pt_index].maxHeight);
                    
                    double max_WNCC = -1;
                    
                    //printf("start pair corrleation\n");
                    if( *plevelinfo.check_matching_rate || nccresult[pt_index].check_height_change)
                    {
                        Accessable_grid ++;
                        
                        nccresult[pt_index].result0 = DoubleToSignedChar_result(-1.0);
                        nccresult[pt_index].result2 = -1000;
                        nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
                        nccresult[pt_index].result3 = -1000;
                        nccresult[pt_index].result4 = 0;
                        
                        for(int grid_voxel_hindex = 0 ; grid_voxel_hindex < nccresult[pt_index].NumOfHeight ; grid_voxel_hindex++)
                        {
                            float iter_height = nccresult[pt_index].minHeight + grid_voxel_hindex*(*plevelinfo.height_step);
                            
                            if(iter_height >= start_H && iter_height <= end_H)
                            {
                                for(int pair_number = 0 ; pair_number < plevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                                {
                                    const int reference_id = plevelinfo.pairinfo->pairs[pair_number].m_X;
                                    const int ti = plevelinfo.pairinfo->pairs[pair_number].m_Y;
                                    
                                    double ref_gsd = plevelinfo.imageinfo[reference_id].GSD.pro_GSD;
                                    double ti_gsd = plevelinfo.imageinfo[ti].GSD.pro_GSD;
                                    double gsd_ratio = ref_gsd/ti_gsd;
                                    
                                    bool check_compute = false;
                                    /*if(GridPT3[pt_index].Matched_flag == 1 || GridPT3[pt_index].Matched_flag == 2)
                                    {
                                        if(check_image_boundary_each(proinfo,plevelinfo,plevelinfo.GridPts[pt_index],plevelinfo.Grid_wgs[pt_index],iter_height,iter_height,Half_template_size,reference_id) && check_image_boundary_each(proinfo,plevelinfo,plevelinfo.GridPts[pt_index],plevelinfo.Grid_wgs[pt_index],iter_height,iter_height,Half_template_size,ti))
                                            check_compute = true;
                                    }
                                    else*/ if(check_image_boundary_each(proinfo,plevelinfo,plevelinfo.GridPts[pt_index],plevelinfo.Grid_wgs[pt_index],nccresult[pt_index].minHeight,nccresult[pt_index].maxHeight,Half_template_size,reference_id,pair_number, true) && check_image_boundary_each(proinfo,plevelinfo,plevelinfo.GridPts[pt_index],plevelinfo.Grid_wgs[pt_index],nccresult[pt_index].minHeight,nccresult[pt_index].maxHeight,Half_template_size,ti,pair_number, false))
                                        check_compute = true;
                                    
                                    if(check_compute)
                                    {
                                        if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                                        {
                                            if(grid_voxel_hindex == floor(nccresult[pt_index].NumOfHeight/2.0))
                                                GridPT3[pt_index].total_images = 1;
                                            
                                            
                                            rsetkernel.reference_id = reference_id;
                                            rsetkernel_next.reference_id = reference_id;
                                            
                                            rsetkernel.ti = ti;
                                            rsetkernel_next.ti = ti;
                                            
                                            KernelPatchArg patch{
                                                rsetkernel,
                                                plevelinfo.py_Sizes[rsetkernel.reference_id][*plevelinfo.Pyramid_step],
                                                plevelinfo.py_Sizes[rsetkernel.ti][*plevelinfo.Pyramid_step],
                                                plevelinfo.py_Images[rsetkernel.reference_id],
                                                plevelinfo.py_MagImages[rsetkernel.reference_id],
                                                plevelinfo.py_Images[rsetkernel.ti],
                                                plevelinfo.py_MagImages[rsetkernel.ti]};
                                            
                                            KernelPatchArg patch_next = check_combined_WNCC ?
                                            KernelPatchArg{
                                                rsetkernel_next,
                                                plevelinfo.py_Sizes[rsetkernel_next.reference_id][*plevelinfo.Pyramid_step-1],
                                                plevelinfo.py_Sizes[rsetkernel_next.ti][*plevelinfo.Pyramid_step-1],
                                                plevelinfo.py_Images_next[rsetkernel_next.reference_id],
                                                plevelinfo.py_MagImages_next[rsetkernel_next.reference_id],
                                                plevelinfo.py_Images_next[rsetkernel_next.ti],
                                                plevelinfo.py_MagImages_next[rsetkernel_next.ti]}
                                            : KernelPatchArg{rsetkernel_next,}; // unused, bit need to init ref type
                                            
                                            
                                            //GNCC computation
                                            double temp_rho = -1.0;
                                            double db_INCC(-1);
                                            double db_GNCC = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[pair_number]);
                                            //nccresult[pt_index].GNCC = GridPT3[pt_index].ortho_ncc[pair_number];
                                            //printf("db_GNCC\n");
                                            //INCC computation
                                            bool check_height_orientation = false;
                                            
                                            //printf("before grid_voxel_hindex\n");
                                            
                                            const CSize LImagesize(plevelinfo.py_Sizes[reference_id][Pyramid_step]);
                                            const CSize RImagesize(plevelinfo.py_Sizes[ti][Pyramid_step]);
                                            
                                            // Image point setting
                                            const double temp_LIA[2] = {0.0, 0.0};
                                            
                                            D2DPOINT Ref_Imagecoord[1], Ref_Imagecoord_py[1];
                                            D2DPOINT Tar_Imagecoord[1], Tar_Imagecoord_py[1];
                                            D3DPOINT temp_GP[1];
                                            D2DPOINT photo;
                                            temp_GP[0].m_Z = (double)iter_height;
                                            if(proinfo->sensor_type == SB)
                                            {
                                                temp_GP[0] = plevelinfo.Grid_wgs[pt_index];
                                                
                                                Ref_Imagecoord[0]      = GetObjectToImageRPC_single(plevelinfo.RPCs[reference_id],*plevelinfo.NumOfIAparam,temp_LIA,temp_GP[0]);
                                                
                                                Tar_Imagecoord[0]     = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pair_number],temp_GP[0]);
                                            }
                                            else
                                            {
                                                temp_GP[0] = plevelinfo.GridPts[pt_index];
                                                
                                                photo = GetPhotoCoordinate_single(temp_GP[0],proinfo->frameinfo.Photoinfo[reference_id],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[reference_id].m_Rm);
                                                Ref_Imagecoord[0] = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
                                                
                                                photo = GetPhotoCoordinate_single(temp_GP[0],proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[ti].m_Rm);
                                                Tar_Imagecoord[0] = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
                                            }
                                            
                                            //D2DPOINT center_Ref(plevelinfo.Imagesize_ori[reference_id].width/2.0,plevelinfo.Imagesize_ori[reference_id].height/2.0);
                                            //D2DPOINT center_Tar(plevelinfo.Imagesize_ori[ti].width/2.0,plevelinfo.Imagesize_ori[ti].height/2.0);
                                            
                                            Ref_Imagecoord_py[0]   = OriginalToPyramid_single(Ref_Imagecoord[0],plevelinfo.py_Startpos[reference_id],Pyramid_step);
                                            Tar_Imagecoord_py[0]  = OriginalToPyramid_single(Tar_Imagecoord[0],plevelinfo.py_Startpos[ti],Pyramid_step);
                                            
                                            const bool check_py_image_pt = (int)Ref_Imagecoord_py[0].m_Y >= 0 && (int)Ref_Imagecoord_py[0].m_Y + 1 < LImagesize.height && (int)Ref_Imagecoord_py[0].m_X >= 0 && (int)Ref_Imagecoord_py[0].m_X + 1 < LImagesize.width && (int)Tar_Imagecoord_py[0].m_Y >= 0 && (int)Tar_Imagecoord_py[0].m_Y + 1 < RImagesize.height && (int)Tar_Imagecoord_py[0].m_X >= 0 && (int)Tar_Imagecoord_py[0].m_X + 1< RImagesize.width;
                                            
                                            // image point setting for next level
                                            bool check_py_image_pt_next = true;
                                            D2DPOINT Ref_Imagecoord_py_next[1];
                                            CSize LImagesize_next, RImagesize_next;
                                            D2DPOINT Tar_Imagecoord_py_next[1];
                                            if(check_combined_WNCC_INCC)
                                            {
                                                LImagesize_next = plevelinfo.py_Sizes[reference_id][Pyramid_step-1];
                                                RImagesize_next = plevelinfo.py_Sizes[ti][Pyramid_step-1];
                                                
                                                Ref_Imagecoord_py_next[0]   = OriginalToPyramid_single(Ref_Imagecoord[0],plevelinfo.py_Startpos_next[reference_id],Pyramid_step-1);
                                                
                                                Tar_Imagecoord_py_next[0]  = OriginalToPyramid_single(Tar_Imagecoord[0],plevelinfo.py_Startpos_next[ti],Pyramid_step-1);
                                                
                                                check_py_image_pt_next = (int)Ref_Imagecoord_py_next[0].m_Y >= 0 && (int)Ref_Imagecoord_py_next[0].m_Y + 1 < LImagesize_next.height && (int)Ref_Imagecoord_py_next[0].m_X >= 0 && (int)Ref_Imagecoord_py_next[0].m_X + 1 < LImagesize_next.width && (int)Tar_Imagecoord_py_next[0].m_Y >= 0 && (int)Tar_Imagecoord_py_next[0].m_Y + 1 < RImagesize_next.height && (int)Tar_Imagecoord_py_next[0].m_X >= 0 && (int)Tar_Imagecoord_py_next[0].m_X + 1< RImagesize_next.width;
                                            }
                                            
                                            if(check_py_image_pt && check_py_image_pt_next)
                                            {
                                                //orientation diff setting
                                                const int ori_diff = plevelinfo.py_OriImages[reference_id][(int)Ref_Imagecoord_py[0].m_Y*LImagesize.width + (int)Ref_Imagecoord_py[0].m_X] - plevelinfo.py_OriImages[ti][(int)Tar_Imagecoord_py[0].m_Y*RImagesize.width + (int)Tar_Imagecoord_py[0].m_X];
                                                
                                                int ori_diff_next = 0;
                                                if(check_combined_WNCC_INCC)
                                                    ori_diff_next = plevelinfo.py_OriImages_next[reference_id][(int)Ref_Imagecoord_py_next[0].m_Y*LImagesize_next.width + (int)Ref_Imagecoord_py_next[0].m_X] - plevelinfo.py_OriImages_next[ti][(int)Tar_Imagecoord_py_next[0].m_Y*RImagesize_next.width + (int)Tar_Imagecoord_py_next[0].m_X];
                                                
                                                double diff_theta  = (double)(ori_diff);
                                                if(check_combined_WNCC_INCC)
                                                    diff_theta = (diff_theta + (double)(ori_diff_next)) / 2.0;
                                                
                                                bool check_orientation = false;
                                                
                                                if(IsRA != 1)
                                                {
                                                    if(Pyramid_step >= 3 && diff_theta*(*plevelinfo.bin_angle) < 90 && diff_theta*(*plevelinfo.bin_angle) > - 90)
                                                        check_orientation = true;
                                                    else if(Pyramid_step <= 2)
                                                        check_orientation = true;
                                                }
                                                else
                                                    check_orientation = true;
                                                
                                                //printf("before check_orientation\n");
                                                //INCC computation
                                                if(check_orientation)
                                                {
                                                    check_height_orientation = true;
                                                    
                                                    double sum_INCC_multi = 0;
                                                    double count_INCC = 0;
                                                    int Count_N[3] = {0};
                                                    int Count_N_next[3] = {0};
                                                    
                                                    const double rot_theta = (double)(diff_theta*(*plevelinfo.bin_angle)*PI/180.0);
                                                    const double cos0 = cos(-rot_theta);
                                                    const double sin0 = sin(-rot_theta);
                                                    
                                                    for(int row = -Half_template_size; row <= Half_template_size ; row++)
                                                    {
                                                        for(int col = -Half_template_size; col <= Half_template_size ; col++)
                                                        {
                                                            int radius2  =  row*row + col*col;
                                                            double right_col = col*gsd_ratio;
                                                            double right_row = row*gsd_ratio;
                                                            if(radius2 <= (Half_template_size + 1)*(Half_template_size + 1))
                                                            {
                                                                D2DPOINT pos_left(Ref_Imagecoord_py[0].m_X + col,Ref_Imagecoord_py[0].m_Y + row);
                                                                D2DPOINT temp_pos(cos0*right_col - sin0*right_row, sin0*right_col + cos0*right_row);
                                                                D2DPOINT pos_right(Tar_Imagecoord_py[0].m_X + temp_pos.m_X,Tar_Imagecoord_py[0].m_Y + temp_pos.m_Y);
                                                                
                                                                SetVecKernelValue(patch, row, col, pos_left,pos_right, radius2, Count_N);
                                                                
                                                                if(check_combined_WNCC_INCC)
                                                                {
                                                                    D2DPOINT pos_left_next(Ref_Imagecoord_py_next[0].m_X + col,Ref_Imagecoord_py_next[0].m_Y + row);
                                                                    D2DPOINT temp_pos_next(cos0*right_col - sin0*right_row, sin0*right_col + cos0*right_row);
                                                                    D2DPOINT pos_right_next(Tar_Imagecoord_py_next[0].m_X + temp_pos_next.m_X,Tar_Imagecoord_py_next[0].m_Y + temp_pos_next.m_Y);
                                                                    
                                                                    SetVecKernelValue(patch_next, row, col, pos_left_next,pos_right_next, radius2, Count_N_next);
                                                                    
                                                                }  // if(check_combined_WNCC_INCC)
                                                            }
                                                        }  // end col loop
                                                    }  // end row loop
                                                    
                                                    // Compute correlations
                                                    ComputeMultiNCC(rsetkernel, TH_N, Count_N, count_INCC,  sum_INCC_multi);
                                                    
                                                    if(check_combined_WNCC_INCC)
                                                        ComputeMultiNCC(rsetkernel_next, TH_N, Count_N_next, count_INCC,  sum_INCC_multi);
                                                    
                                                    //printf("sum_INCC_multi %f\n",sum_INCC_multi);
                                                    if(count_INCC > 0)
                                                    {
                                                        db_INCC = sum_INCC_multi/count_INCC;
                                                        /*if(db_INCC > 0 && plevelinfo.pairinfo->NumberOfPairs > 2)
                                                        {
                                                            double DistToCenter_Ref = SQRT(Ref_Imagecoord[0],center_Ref);
                                                            double DistToCenter_Tar = SQRT(Tar_Imagecoord[0],center_Tar);
                                                            double DistToCenter_avg = (DistToCenter_Ref+DistToCenter_Tar)/2.0;
                                                            double WeightDist = 1.0 + (20 - (DistToCenter_avg/SQRT(center_Ref))*20.0)/100.0;
                                                            
                                                            db_INCC *= WeightDist;
                                                        }*/
                                                        
                                                        double gncc_weight = SetGnccWeight(Pyramid_step, db_GNCC, db_INCC, GridPT3[pt_index].Height, iter_height, *plevelinfo.height_step);
                                                        
                                                        if((Pyramid_step == 4 && iteration == 1) || IsRA)
                                                            gncc_weight = 1.0;
                                                        
                                                        if(check_ortho) // GNCC check
                                                        {
                                                            if(db_GNCC > ortho_th)
                                                                temp_rho = db_INCC*ncc_alpha + db_GNCC*ncc_beta;
                                                            else
                                                                temp_rho = db_INCC*gncc_weight;
                                                        }
                                                        else
                                                        {
                                                            temp_rho = db_INCC*gncc_weight;
                                                        }
                                                        
                                                        if(max_WNCC < db_INCC)
                                                            max_WNCC = db_INCC;
                                                    }
                                                    else
                                                    {
                                                        temp_rho = -1.0;
                                                        db_INCC = -1.0;
                                                    }
                                                }
                                                //printf("after check_orientation\n");
                                                /*                       else
                                                 {
                                                 if(!(*plevelinfo.check_matching_rate) && !IsRA)
                                                 grid_voxel[pt_index][grid_voxel_hindex].flag_cal = false;
                                                 }*/
                                            }
                                            
                                            if(!(*plevelinfo.check_matching_rate) && nccresult[pt_index].check_height_change)
                                            {
                                                if(check_height_orientation && temp_rho > -1)
                                                {
                                                    grid_voxel[pt_index].flag_cal(grid_voxel_hindex, pair_number) = true;
                                                    grid_voxel[pt_index].INCC(grid_voxel_hindex, pair_number) = DoubleToSignedChar_voxel(db_INCC);
                                                }
                                                else
                                                {
                                                    grid_voxel[pt_index].flag_cal(grid_voxel_hindex, pair_number) = false;
                                                    grid_voxel[pt_index].INCC(grid_voxel_hindex, pair_number) = DoubleToSignedChar_voxel(-1);
                                                }
                                            }
                                            
                                            //find peak position
                                            if(*plevelinfo.check_matching_rate)
                                            {
                                                if(check_height_orientation && temp_rho > -1)
                                                {
                                                    WNCC_save[grid_voxel_hindex].push_back(temp_rho);
                                                    WNCC_save_pair_ID[grid_voxel_hindex].push_back(pair_number);
                                                    //WNCC_save_count[grid_voxel_hindex].push_back(pair_number);
                                                }
                                            }
                                            
                                        }//end iter_height loop
                                    }// end grid_voxel_hindex loop
                                    //printf("after grid_voxel_hindex\n");
                                    //if(nccresult[pt_index].check_height_change)
                                    //    nccresult[pt_index].max_WNCC = DoubleToSignedChar_result(max_WNCC);
                                }//end INCC computation
                            }
                        }
                    }  // end ti loop
                    
                    //average of WNCCs
                    if( *plevelinfo.check_matching_rate || nccresult[pt_index].check_height_change )
                    {
                        //printf("before peak\n");
                        if(*plevelinfo.check_matching_rate) //no 3D Voxel structure
                        {
                            //find peak
                            double pre_rho  = -1.0;
                            float pre_height= 0.0;
                            int direction   = 0;
                            bool check_rho  = false;
                            double temp_max_WNCC = -100;
                            double pre_rho_WNCC = -1.0;
                            
                            double temp_nccresult = -100.0;
                            double temp_nccresult_sec = -100.0;
                            
                            double max_WNCC_save = -100;
                            int max_WNCC_pair_ID;
                            for(int grid_voxel_hindex = 0 ; grid_voxel_hindex < nccresult[pt_index].NumOfHeight ; grid_voxel_hindex++)
                            {
                                float iter_height = nccresult[pt_index].minHeight + grid_voxel_hindex*(*plevelinfo.height_step);
                                
                                double temp_rho(0);
                                if(WNCC_save[grid_voxel_hindex].size() > 0)
                                {
                                    for(int pair = 0 ; pair < WNCC_save[grid_voxel_hindex].size() ; pair++)
                                    {
                                        temp_rho += WNCC_save[grid_voxel_hindex][pair];
                                        if(max_WNCC_save < WNCC_save[grid_voxel_hindex][pair])
                                        {
                                            max_WNCC_save = WNCC_save[grid_voxel_hindex][pair];
                                            max_WNCC_pair_ID = WNCC_save_pair_ID[grid_voxel_hindex][pair];
                                        }
                                            
                                    }
                                    
                                    temp_rho /= WNCC_save[grid_voxel_hindex].size();
                                    double WNCC_temp_rho = temp_rho;
                                    //find peak position
                                    FindPeakNcc2(Pyramid_step, iteration, temp_rho, iter_height, check_rho, pre_rho, pre_rho_WNCC, WNCC_temp_rho, pre_height, direction, temp_max_WNCC, nccresult[pt_index], temp_nccresult, temp_nccresult_sec);
                                    //FindPeakNcc(*plevelinfo.Pyramid_step, iteration, pt_index, temp_rho, iter_height, check_rho, pre_rho, pre_height, direction, max_WNCC, nccresult);
                                }
                            }
                            
                            if(temp_nccresult > 1.0)
                            {
                                temp_nccresult = 1.0;
                                temp_nccresult_sec = -1.0;
                            }
                            
                            nccresult[pt_index].result0 = DoubleToSignedChar_result(temp_nccresult);
                            nccresult[pt_index].result1 = DoubleToSignedChar_result(temp_nccresult_sec);
                            
                            if(plevelinfo.pairinfo->NumberOfPairs == 1)
                                GridPT3[pt_index].ncc_seleceted_pair = 0;
                            else
                                GridPT3[pt_index].ncc_seleceted_pair = max_WNCC_pair_ID;
                            
                            if(GridPT3[pt_index].ncc_seleceted_pair >= plevelinfo.pairinfo->NumberOfPairs)
                            {
                                printf("VerticalLineLocus GridPT3[pt_index].ncc_seleceted_pair %d\n",GridPT3[pt_index].ncc_seleceted_pair);
                                exit(1);
                            }
                            //printf("temp_nccresult %d\t%d\t%f\t%f\t%d\n",nccresult[pt_index].result0,nccresult[pt_index].result1,nccresult[pt_index].result2,nccresult[pt_index].result3,nccresult[pt_index].result4);
                        }
                        
                        if(nccresult[pt_index].check_height_change)
                            nccresult[pt_index].max_WNCC = DoubleToSignedChar_result(max_WNCC);
                        //printf("after peak Pyramid_step %d\t%d\n",Pyramid_step,iteration);
                    }
                    
                    if(*plevelinfo.check_matching_rate)
                    {
                        for(int grid_voxel_hindex = 0 ; grid_voxel_hindex < nccresult[pt_index].NumOfHeight ; grid_voxel_hindex++)
                        {
                            WNCC_save[grid_voxel_hindex].clear();
                            vector<double>().swap(WNCC_save[grid_voxel_hindex]);
                            
                            WNCC_save_pair_ID[grid_voxel_hindex].clear();
                            vector<unsigned char>().swap(WNCC_save_pair_ID[grid_voxel_hindex]);
                        }
                        
                        //for(int grid_voxel_hindex = 0 ; grid_voxel_hindex < nccresult[pt_index].NumOfHeight ; grid_voxel_hindex++)
                        //    WNCC_save_count[grid_voxel_hindex].clear();
                            
                        free(WNCC_save);
                        free(WNCC_save_pair_ID);
                        //free(WNCC_save_count);
                    }
                    
                }//end check blunder cell
            }//end check_image_boundary
        } // end omp for
    } // end omp parallel
    
    //printf("check height cell %d\t%d\t%d\t%f\t%f\n",sum_data2,sum_data,sum_data2+sum_data,(double)sum_data2/(double)(sum_data2+sum_data)*100,(double)sum_data/(double)(sum_data2+sum_data)*100);
    
    return Accessable_grid;
}  // end VerticalLineLocus

void SetOrthoImageCoord(const ProInfo *proinfo, LevelInfo &plevelinfo, const UGRID *GridPT3, const bool check_combined_WNCC, enum PyImageSelect check_pyimage, const double im_resolution, const double im_resolution_next, long int &sub_imagesize_w, long int &sub_imagesize_h, long int &sub_imagesize_w_next, long int &sub_imagesize_h_next, D2DPOINT **all_im_cd, D2DPOINT **all_im_cd_next)
{
    sub_imagesize_w = (long int)((plevelinfo.Boundary[2] - plevelinfo.Boundary[0])/im_resolution)+1;
    sub_imagesize_h = (long int)((plevelinfo.Boundary[3] - plevelinfo.Boundary[1])/im_resolution)+1;
    
    if(check_combined_WNCC)
    {
        sub_imagesize_w_next = (long int)((plevelinfo.Boundary[2] - plevelinfo.Boundary[0])/im_resolution_next)+1;
        sub_imagesize_h_next = (long int)((plevelinfo.Boundary[3] - plevelinfo.Boundary[1])/im_resolution_next)+1;
    }
    
    const long int sub_imagesize_total = (long int)sub_imagesize_w * (long int)sub_imagesize_h;
    long int sub_imagesize_total_next;
    if(check_combined_WNCC)
        sub_imagesize_total_next = (long int)sub_imagesize_w_next * (long int)sub_imagesize_h_next;
    
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(proinfo->check_selected_image[ti])
        {
            all_im_cd[ti] = (D2DPOINT*)calloc(sizeof(D2DPOINT),sub_imagesize_total);
            if (all_im_cd[ti] == NULL)
            {
                printf("ERROR: Out of memory - all_left_im_cd is NULL\n");
                exit(1);
            }
            
            if(check_combined_WNCC)
            {
                all_im_cd_next[ti] = (D2DPOINT*)calloc(sizeof(D2DPOINT),sub_imagesize_total_next);
                if (all_im_cd_next[ti] == NULL)
                {
                    printf("ERROR: Out of memory - all_left_im_cd_next is NULL\n");
                    exit(1);
                }
            }
        }
    }
    
#pragma omp parallel for schedule(guided)
    for(long int iter_count = 0 ; iter_count < sub_imagesize_total ; iter_count++)
    {
        long int pts_row = (long int)(floor(iter_count/sub_imagesize_w));
        long int pts_col = iter_count % sub_imagesize_w;
        
        double t_X     = plevelinfo.Boundary[0] + pts_col*im_resolution;
        double t_Y     = plevelinfo.Boundary[1] + pts_row*im_resolution;
        
        long int t_col   = (long int)((t_X - plevelinfo.Boundary[0])/(*plevelinfo.grid_resolution));
        long int t_row   = (long int)((t_Y - plevelinfo.Boundary[1])/(*plevelinfo.grid_resolution));
        
        long int pt_index    = t_row*(long int)plevelinfo.Size_Grid2D->width + t_col;
        long int pt_index_im = pts_row*(long int)sub_imagesize_w + pts_col;
        
        
        long int t_col_next, t_row_next;
        long int pt_index_im_next;
        if(check_combined_WNCC)
        {
            t_col_next   = (long int)((t_X - plevelinfo.Boundary[0])/im_resolution_next);
            t_row_next   = (long int)((t_Y - plevelinfo.Boundary[1])/im_resolution_next);
            
            pt_index_im_next = t_row_next*(long int)sub_imagesize_w_next + t_col_next;
        }
        
        if(pt_index < *plevelinfo.Grid_length && t_col < plevelinfo.Size_Grid2D->width && t_row < plevelinfo.Size_Grid2D->height)
        {
            if(GridPT3[pt_index].Height != -1000)
            {
                D3DPOINT temp_GP;
                D2DPOINT temp_GP_p;
                
                if(proinfo->sensor_type == SB)
                {
                    temp_GP_p.m_X = t_X;
                    temp_GP_p.m_Y = t_Y;
                    
                    temp_GP     = ps2wgs_single(*plevelinfo.param,temp_GP_p);
                }
                else
                {
                    temp_GP.m_X   = t_X;
                    temp_GP.m_Y   = t_Y;
                }
                temp_GP.m_Z   = (double)GridPT3[pt_index].Height;
                
                for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
                {
                    if(proinfo->check_selected_image[ti])
                    {
                        D2DPOINT Imagecoord;
                        D2DPOINT Imagecoord_py;
                        if(proinfo->sensor_type == SB)
                            Imagecoord      = GetObjectToImageRPC_single_noBias(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,temp_GP);
                        else
                        {
                            D2DPOINT photo  = GetPhotoCoordinate_single(temp_GP,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[ti].m_Rm);
                            Imagecoord      = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
                        }
                        
                        if(check_pyimage == OR)
                            Imagecoord_py  = OriginalToPyramid_single(Imagecoord,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);
                        else
                            Imagecoord_py  = OriginalToPyramid_single(Imagecoord,plevelinfo.py_BStartpos[ti],*plevelinfo.blunder_selected_level);
                        
                        all_im_cd[ti][pt_index_im] = Imagecoord_py;
                        
                        if(check_combined_WNCC)
                        {
                            Imagecoord_py  = OriginalToPyramid_single(Imagecoord,plevelinfo.py_Startpos_next[ti],*plevelinfo.Pyramid_step - 1);
                            if(pt_index_im_next < sub_imagesize_w_next*sub_imagesize_h_next && t_col_next < sub_imagesize_w_next && t_row_next < sub_imagesize_h_next)
                            {
                                all_im_cd_next[ti][pt_index_im_next] = Imagecoord_py;
                            }
                        }
                    }
                }
            }
        }
    }
}

void FindPeakNcc(const int Pyramid_step, const int iteration, const long int grid_index, const double temp_rho, const float iter_height, bool &check_rho, double &pre_rho, float &pre_height, int &direction, double &max_WNCC, NCCresult *nccresult)
{
    double diff_rho;
    int t_direction;
    if(temp_rho >= 0 && pre_rho >= 0)
        diff_rho = temp_rho - pre_rho;
    else if(temp_rho < 0 && pre_rho < 0)
        diff_rho = -(fabs(temp_rho) - fabs(pre_rho));
    else if(temp_rho >= 0 && pre_rho < 0)
        diff_rho = 1;
    else
        diff_rho = -1;
    
    if(diff_rho > 0)
        t_direction = 1;
    else if(diff_rho < 0)
        t_direction = -1;
    else
        t_direction = 0;
    
    if(!check_rho)
    {
        check_rho   = true;
        t_direction = -1;
    }
    else if(pre_rho != -1)
    {
        if((iteration <= 3 && Pyramid_step == 4) || (iteration <= 2 && Pyramid_step == 3) || (iteration <= 1 && Pyramid_step == 2)) //Max AWNCC
        {
            if(nccresult[grid_index].result0 < DoubleToSignedChar_result(temp_rho))
            {
                if(nccresult[grid_index].result4 < 10)
                    nccresult[grid_index].result4 += 1;
                
                nccresult[grid_index].result0 = SignedCharToDouble_result(temp_rho);
                nccresult[grid_index].result2 = iter_height;
                
                nccresult[grid_index].result1 = SignedCharToDouble_result(-1.0);
                nccresult[grid_index].result2 = Nodata;
            }
        }
        else
        {
            bool check_peak = false;
            
            if (direction > 0 && t_direction < 0)
                check_peak = true;
            
            if( check_peak )
            {
                if(nccresult[grid_index].result4 < 10)
                    nccresult[grid_index].result4 += 1;
                if(nccresult[grid_index].result0 < DoubleToSignedChar_result(pre_rho))
                {
                    double temp_1 = SignedCharToDouble_result(nccresult[grid_index].result0);
                    nccresult[grid_index].result0 = DoubleToSignedChar_result(pre_rho);
                    
                    float temp_2 = nccresult[grid_index].result2;
                    nccresult[grid_index].result2 = pre_height;
                    
                    if(nccresult[grid_index].result1 < DoubleToSignedChar_result(temp_1))
                    {
                        nccresult[grid_index].result1 = DoubleToSignedChar_result(temp_1);
                        nccresult[grid_index].result3 = temp_2;
                    }
                }
                else
                {
                    if(nccresult[grid_index].result1 < DoubleToSignedChar_result(pre_rho))
                    {
                        nccresult[grid_index].result1 = DoubleToSignedChar_result(pre_rho);
                        nccresult[grid_index].result3 = pre_height;
                    }
                }
            }
        }
    }
    
    pre_rho                = temp_rho;
    pre_height             = iter_height;
    direction              = t_direction;
    
    if(max_WNCC < temp_rho)
        max_WNCC = temp_rho;
}

void FindPeakNcc2(const int Pyramid_step, const int iteration, const double temp_rho, const float iter_height, bool &check_rho, double &pre_rho, double &pre_rho_WNCC, double WNCC_temp_rho, float &pre_height, int &direction, double &max_roh, NCCresult &nccresult, double &temp_nccresult, double &temp_nccresult_sec)
{
    double diff_rho;
    int t_direction;
    if(temp_rho >= 0 && pre_rho >= 0)
        diff_rho = temp_rho - pre_rho;
    else if(temp_rho < 0 && pre_rho < 0)
        diff_rho = -(fabs(temp_rho) - fabs(pre_rho));
    else if(temp_rho >= 0 && pre_rho < 0)
        diff_rho = 1;
    else
        diff_rho = -1;
    
    if(diff_rho > 0)
        t_direction = 1;
    else if(diff_rho < 0)
        t_direction = -1;
    else
        t_direction = 0;
    
    if(!check_rho)
    {
        check_rho   = true;
        t_direction = -1;
    }
    else if(pre_rho != -1)
    {
        if((iteration <= 3 && Pyramid_step == 4) || (iteration <= 2 && Pyramid_step == 3) || (iteration <= 1 && Pyramid_step == 2)) //Max AWNCC
        {
            if(temp_nccresult < temp_rho)
            {
                if(nccresult.result4 < 10)
                    nccresult.result4 += 1;
                
                temp_nccresult = temp_rho;
                nccresult.result2 = iter_height;
                temp_nccresult_sec = -1.0;
                nccresult.result3 = Nodata;

                max_roh = WNCC_temp_rho;
                
            }
        }
        else
        {
            bool check_peak = false;
            
            if (direction > 0 && t_direction < 0)
                check_peak = true;
            
            if( check_peak )
            {
                if(nccresult.result4 < 10)
                    nccresult.result4 += 1;
                if(temp_nccresult < pre_rho)
                {
                    double temp_1;
                    float temp_2;
                    temp_1 = temp_nccresult;
                    temp_nccresult = pre_rho;
                    
                    temp_2 = nccresult.result2;
                    
                    nccresult.result2 = pre_height;
                    
                    max_roh = pre_rho_WNCC;
                    
                    if(temp_nccresult_sec < temp_1)
                    {
                        temp_nccresult_sec = temp_1;
                        nccresult.result3 = temp_2;
                    }
                }
                else
                {
                    if(temp_nccresult_sec < pre_rho)
                    {
                        temp_nccresult_sec = pre_rho;
                        nccresult.result3 = pre_height;
                    }
                }
            }
        }
    }
    
    pre_rho                = temp_rho;
    pre_rho_WNCC           = WNCC_temp_rho; //for max_roh < 1.0
    pre_height             = iter_height;
    
    direction              = t_direction;
}

void FindPeakNcc_SGM(const int Pyramid_step, const int iteration, const double temp_rho, const float iter_height, bool &check_rho, double &pre_rho, double &pre_rho_WNCC, double WNCC_temp_rho, float &pre_height, int &direction, double &max_roh, double &max_roh_sec, NCCresult &nccresult, double &temp_nccresult, double &temp_nccresult_sec)
{
    double diff_rho;
    int t_direction;
    if(temp_rho >= 0 && pre_rho >= 0)
        diff_rho = temp_rho - pre_rho;
    else if(temp_rho < 0 && pre_rho < 0)
        diff_rho = -(fabs(temp_rho) - fabs(pre_rho));
    else if(temp_rho >= 0 && pre_rho < 0)
        diff_rho = 1;
    else
        diff_rho = -1;
    
    if(diff_rho > 0)
        t_direction = 1;
    else if(diff_rho < 0)
        t_direction = -1;
    else
        t_direction = 0;
    
    if(!check_rho)
    {
        check_rho   = true;
        t_direction = -1;
    }
    else if(pre_rho != -1)
    {
        if((iteration <= 3 && Pyramid_step == 4) || (iteration <= 2 && Pyramid_step == 3) || (iteration <= 1 && Pyramid_step == 2)) //Max AWNCC
        {
            if(temp_nccresult < temp_rho)
            {
                if(nccresult.result4 < 10)
                    nccresult.result4 += 1;
                
                temp_nccresult = temp_rho;
                nccresult.result2 = iter_height;
                temp_nccresult_sec = -1.0;
                nccresult.result3 = Nodata;

                max_roh = WNCC_temp_rho;
                
            }
        }
        else
        {
            bool check_peak = false;
            
            if (direction > 0 && t_direction < 0)
                check_peak = true;
            
            if( check_peak )
            {
                if(nccresult.result4 < 10)
                    nccresult.result4 += 1;
                if(temp_nccresult < pre_rho)
                {
                    double temp_1, temp_max_roh;
                    float temp_2;
                    temp_1 = temp_nccresult;
                    temp_nccresult = pre_rho;
                    
                    temp_max_roh = max_roh;
                    
                    temp_2 = nccresult.result2;
                    
                    nccresult.result2 = pre_height;
                    
                    max_roh = pre_rho_WNCC;
                    
                    if(temp_nccresult_sec < temp_1)
                    {
                        max_roh_sec = temp_max_roh;
                        
                        temp_nccresult_sec = temp_1;
                        nccresult.result3 = temp_2;
                    }
                }
                else
                {
                    if(temp_nccresult_sec < pre_rho)
                    {
                        max_roh_sec = pre_rho_WNCC;
                        
                        temp_nccresult_sec = pre_rho;
                        nccresult.result3 = pre_height;
                    }
                }
            }
        }
    }
    
    pre_rho                = temp_rho;
    pre_rho_WNCC           = WNCC_temp_rho; //for max_roh < 1.0
    pre_height             = iter_height;
    
    direction              = t_direction;
}


void SGM_start_pos(NCCresult *nccresult, GridVoxel &grid_voxel, LevelInfo &rlevelinfo, UGRID *GridPT3, long pt_index, float* LHcost_pre,float **SumCost, double height_step_interval, const int pairnumber)
{
    for(int height_step = 0 ; height_step < nccresult[pt_index].NumOfHeight ; height_step++)
    {
        float iter_height = nccresult[pt_index].minHeight + height_step*height_step_interval;//
        if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
        {
            double WNCC_sum = 0;
            int pair_count = 0;
            
            int count = pairnumber;
            //for(int count = 0 ; count < rlevelinfo.pairinfo->NumberOfPairs ; count++)
            {
                double bhratio = 1.0;
                //if(count > 1)
                //    bhratio = 1.0+rlevelinfo.pairinfo->BHratio[count];
                
                int reference_id = rlevelinfo.pairinfo->pairs[count].m_X;
                int ti = rlevelinfo.pairinfo->pairs[count].m_Y;
                if((reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& GridPT3[pt_index].ncc_seleceted_pair == count*/)
                {
                    if(grid_voxel[pt_index].flag_cal(height_step, count))
                    {
                        WNCC_sum += SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, count))*bhratio;
                        pair_count++;
                    }
                }
            }
            if(pair_count > 0)
                WNCC_sum /= (double)pair_count;
           
            LHcost_pre[height_step] = WNCC_sum;
            SumCost[pt_index][height_step] += LHcost_pre[height_step];
        }
    }
}

void SGM_con_pos(int pts_col, int pts_row, CSize Size_Grid2D, int direction_iter, double step_height, int P_HS_step, int *u_col, int *v_row, NCCresult *nccresult, GridVoxel &grid_voxel,UGRID *GridPT3, LevelInfo &rlevelinfo, long pt_index, double P1, double P2, float* LHcost_pre, float* LHcost_curr, float **SumCost, const int pairnumber)
{
    for(int height_step = 0 ; height_step < nccresult[pt_index].NumOfHeight ; height_step++)
    {
        const float iter_height = nccresult[pt_index].minHeight + height_step*step_height;
        if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
        {
            double WNCC_sum = 0;
            int pair_count = 0;
            double sum_LHcost_curr = 0;
            
            int count = pairnumber;
            //for(int count = 0 ; count < rlevelinfo.pairinfo->NumberOfPairs ; count++)
            {
                double bhratio = 1.0;
                //if(count > 1)
                //    bhratio = 1.0+rlevelinfo.pairinfo->BHratio[count];
                
                int reference_id = rlevelinfo.pairinfo->pairs[count].m_X;
                int ti = rlevelinfo.pairinfo->pairs[count].m_Y;
                if((reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& GridPT3[pt_index].ncc_seleceted_pair == count*/)
                {
                    if(grid_voxel[pt_index].flag_cal(height_step, count))
                    {
                        WNCC_sum = SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, count));
                        
                        const int t_col = pts_col + u_col[direction_iter];
                        const int t_row = pts_row + v_row[direction_iter];
                        const int t_index = t_row*Size_Grid2D.width + t_col;
                        
                        float maxWNCC = SignedCharToDouble_result(nccresult[t_index].max_WNCC);
                        
                        if(t_col >= 0 && t_col < Size_Grid2D.width && t_row >= 0 && t_row < Size_Grid2D.height)
                        {
                            const int t_index_h_index   = (int)((iter_height - nccresult[t_index].minHeight)/step_height);
                            const int t_index_h_index_1 = (int)((iter_height - nccresult[t_index].minHeight)/step_height) - P_HS_step;
                            const int t_index_h_index_2 = (int)((iter_height - nccresult[t_index].minHeight)/step_height) + P_HS_step;
                            double V1 = 0;
                            double V2 = 0;
                            double V3 = 0;
                            double V4 = 0;
                            
                            if(t_index_h_index == 0 && t_index_h_index_2 >= 0 && t_index_h_index_2 < nccresult[t_index].NumOfHeight)
                            {
                                if(grid_voxel[t_index].flag_cal(t_index_h_index, count) && grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value23 = V2 > V3 ? V2 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index, count))
                                {
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V3 > V4 ? V3 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V2 > V4 ? V2 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else
                                {
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                            }
                            else if(t_index_h_index == nccresult[t_index].NumOfHeight - 1 && t_index_h_index_1 >= 0 && t_index_h_index_1 < nccresult[t_index].NumOfHeight)
                            {
                                if(grid_voxel[t_index].flag_cal(t_index_h_index, count) && grid_voxel[t_index].flag_cal(t_index_h_index_1, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value23 = V1 > V3 ? V1 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index, count))
                                {
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V3 > V4 ? V3 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index_1, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V1 > V4 ? V1 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else
                                {
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                            }
                            else if(t_index_h_index   > 0 && t_index_h_index   < nccresult[t_index].NumOfHeight - 1 &&
                                    t_index_h_index_1 >= 0 && t_index_h_index_1 < nccresult[t_index].NumOfHeight &&
                                    t_index_h_index_2 >= 0 && t_index_h_index_2 < nccresult[t_index].NumOfHeight )
                            {
                                
                                if(grid_voxel[t_index].flag_cal(t_index_h_index, count) && grid_voxel[t_index].flag_cal(t_index_h_index_1, count)
                                   && grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value12 = V1 > V2 ? V1 : V2;
                                    double max_value23 = max_value12 > V3 ? max_value12 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index, count) && grid_voxel[t_index].flag_cal(t_index_h_index_1, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value12 = V1;
                                    double max_value23 = max_value12 > V3 ? max_value12 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index_1, count) && grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value12 = V1 > V2 ? V1 : V2;
                                    double max_value23 = max_value12 > V3 ? max_value12 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index, count) && grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value12 = V2;
                                    double max_value23 = max_value12 > V3 ? max_value12 : V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index, count))
                                {
                                    V3 = LHcost_pre[t_index_h_index];
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value23 = V3;
                                    double max_value = max_value23 > V4 ? max_value23 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index_1, count))
                                {
                                    V1 = LHcost_pre[t_index_h_index_1] - P1;
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V1 > V4 ? V1 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else if(grid_voxel[t_index].flag_cal(t_index_h_index_2, count))
                                {
                                    V2 = LHcost_pre[t_index_h_index_2] - P1;
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V2 > V4 ? V2 : V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                                else
                                {
                                    V4 = maxWNCC - P2;
                                    
                                    double max_value = V4;
                                    double t_WNCC_sum = max_value - maxWNCC;
                                    sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                    pair_count++;
                                    //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                    //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                                }
                            }
                            else
                            {
                                V4 = maxWNCC - P2;
                                
                                double max_value = V4;
                                double t_WNCC_sum = max_value - maxWNCC;
                                
                                sum_LHcost_curr += (WNCC_sum + t_WNCC_sum)*bhratio;
                                pair_count++;
                                //LHcost_curr[height_step] = WNCC_sum + t_WNCC_sum;
                                //SumCost[pt_index][height_step] += LHcost_curr[height_step];
                            }
                        }
                        
                    }
                }
            }
    
            if(pair_count > 0)
                LHcost_curr[height_step] = sum_LHcost_curr/(double)pair_count;
            else
                LHcost_curr[height_step] = sum_LHcost_curr;
            
            SumCost[pt_index][height_step] += LHcost_curr[height_step];
            
            //if(pair_count > 0)
            //    WNCC_sum /= (double)pair_count;
            
            
        }
    }
}


void AWNCC_single(ProInfo *proinfo, GridVoxel &grid_voxel,LevelInfo &rlevelinfo,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel, double *minmaxHeight, int pair_number)
{
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    printf("AWNCC single find peak pairs %d\t%d\t%f\n",rlevelinfo.pairinfo->NumberOfPairs,pair_number,im_resolution);
    
    bool pre_DEMtif = proinfo->pre_DEMtif;
    bool check_ortho = false;
    
    if(Pyramid_step == 4 && iteration == 1)
        check_ortho = false;
    else
        check_ortho = true;
    
    if(pre_DEMtif && !(Pyramid_step == 4 && iteration == 1))
        check_ortho = true;
    
    const double ncc_alpha = SetNCC_alpha(Pyramid_step,iteration, proinfo->IsRA);
    const double ncc_beta = 1.0 - ncc_alpha;
    const double ortho_th = 0.7 - (4 - Pyramid_step)*0.10;
    
    int height_buffer = 2;
    if(Pyramid_step <= 1)
        height_buffer = 0;
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        
        double pre_rho  = -1.0;
        double pre_rho_WNCC = -1.0;
        float pre_height= 0.0;
        int pre_height_step = 0;
        int direction   = 0;
        bool check_rho  = false;
        
        double max_roh = 0;
        
        double temp_nccresult = -100.0;
        double temp_nccresult_sec = -100.0;
        
        nccresult[pt_index].result0 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result2 = Nodata;
        nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result3 = Nodata;
        nccresult[pt_index].result4 = 0;
    
        double db_GNCC, db_INCC;
        
        for(long height_step = 0 + height_buffer ; height_step < nccresult[pt_index].NumOfHeight - height_buffer ; height_step++)
        {
            float iter_height = nccresult[pt_index].minHeight + height_step*step_height;
            
            if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
            {
                //int pair_number = 0 ;
                int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                {
                    db_GNCC = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[pair_number]);
                    db_INCC = SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, pair_number));
                    
                    double temp_rho = 0;
                    double WNCC_temp_rho = 0;
                    WNCC_temp_rho = db_INCC;
                    
                    if(grid_voxel[pt_index].flag_cal(height_step, pair_number) && db_INCC > -1)
                    {
                        double gncc_weight = SetGnccWeight(Pyramid_step, db_GNCC, db_INCC, GridPT3[pt_index].Height, iter_height, step_height);
                        
                        if((Pyramid_step == 4 && iteration == 1))
                            gncc_weight = 1.0;
                        
                        if(check_ortho) // GNCC check
                        {
                            if(db_GNCC > ortho_th)
                                temp_rho = db_INCC*ncc_alpha + db_GNCC*ncc_beta;
                            else
                                temp_rho = db_INCC*gncc_weight;
                        }
                        else
                        {
                            temp_rho = db_INCC*gncc_weight;
                        }
                        
                        FindPeakNcc2(Pyramid_step, iteration, temp_rho, iter_height, check_rho, pre_rho, pre_rho_WNCC, WNCC_temp_rho, pre_height, direction, max_roh, nccresult[pt_index], temp_nccresult, temp_nccresult_sec);
                    }
                }
            }
        }
        
        GridPT3[pt_index].total_images = 1;
        
        if(rlevelinfo.pairinfo->NumberOfPairs == 1)
            GridPT3[pt_index].ncc_seleceted_pair = pair_number;
        
        /*if(temp_nccresult > 1.0)
        {
            temp_nccresult = 1.0;
            temp_nccresult_sec = -1.0;
        }
        */
        nccresult[pt_index].result0 = DoubleToSignedChar_result(temp_nccresult);
        nccresult[pt_index].result1 = DoubleToSignedChar_result(temp_nccresult_sec);
        
        //printf("nccresult %d\t%d\t%f\t%f\t%d\n",nccresult[pt_index].result0,nccresult[pt_index].result1,nccresult[pt_index].result2,nccresult[pt_index].result3,nccresult[pt_index].result4);
    }
}

void AWNCC_AWNCC(ProInfo *proinfo, GridVoxel &grid_voxel,LevelInfo &rlevelinfo,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel, double *minmaxHeight)
{
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    printf("AWNCC single find peak pairs %d\t%d\t%f\n",rlevelinfo.pairinfo->NumberOfPairs,rlevelinfo.pairinfo->NumberOfPairs,im_resolution);
    
    bool pre_DEMtif = proinfo->pre_DEMtif;
    bool check_ortho = false;
    
    if(Pyramid_step == 4 && iteration == 1)
        check_ortho = false;
    else
        check_ortho = true;
    
    if(pre_DEMtif && !(Pyramid_step == 4 && iteration == 1))
        check_ortho = true;
    
    const double ncc_alpha = SetNCC_alpha(Pyramid_step,iteration, proinfo->IsRA);
    const double ncc_beta = 1.0 - ncc_alpha;
    const double ortho_th = 0.7 - (4 - Pyramid_step)*0.10;
    
    int AWNCC_id = rlevelinfo.pairinfo->NumberOfPairs;
    
    int height_buffer = 2;
    if(Pyramid_step <= 1)
        height_buffer = 0;
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        
        double pre_rho  = -1.0;
        double pre_rho_WNCC = -1.0;
        float pre_height= 0.0;
        int pre_height_step = 0;
        int direction   = 0;
        bool check_rho  = false;
        
        double max_roh = 0;
        
        double temp_nccresult = -100.0;
        double temp_nccresult_sec = -100.0;
        
        nccresult[pt_index].result0 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result2 = Nodata;
        nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result3 = Nodata;
        nccresult[pt_index].result4 = 0;
    
        double db_GNCC, db_INCC;
        
        for(long height_step = 0 + height_buffer ; height_step < nccresult[pt_index].NumOfHeight - height_buffer ; height_step++)
        {
            float iter_height = nccresult[pt_index].minHeight + height_step*step_height;
            
            if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
            {
                double sum_INCC = 0;
                double sum_WNCC = 0;
                double sum_GNCC = 0;
                int INCC_count = 0;
                
                for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                {
                    int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                    int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                    if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                    {
                        db_GNCC = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[pair_number]);
                        db_INCC = SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, pair_number));
                        
                        double temp_rho = 0;
                        
                        if(grid_voxel[pt_index].flag_cal(height_step, pair_number) && db_INCC > -1)
                        {
                            double gncc_weight = SetGnccWeight(Pyramid_step, db_GNCC, db_INCC, GridPT3[pt_index].Height, iter_height, step_height);
                            
                            if((Pyramid_step == 4 && iteration == 1))
                                gncc_weight = 1.0;
                            
                            if(check_ortho) // GNCC check
                            {
                                if(db_GNCC > ortho_th)
                                    temp_rho = db_INCC*ncc_alpha + db_GNCC*ncc_beta;
                                else
                                    temp_rho = db_INCC*gncc_weight;
                            }
                            else
                            {
                                temp_rho = db_INCC*gncc_weight;
                            }
                            
                            sum_WNCC += temp_rho;
                            sum_INCC += db_INCC;
                            sum_GNCC += db_GNCC;
                            INCC_count++;
                        }
                    }
                }
                
                double temp_rho = -1;
                double WNCC_temp_rho = -1;
                if(INCC_count > 0)
                    WNCC_temp_rho = sum_INCC / (double)INCC_count;
                else
                    WNCC_temp_rho = -1;
                
                if(INCC_count > 0)
                {
                    double Multiweight = 1.0;
                    //if(INCC_count > 1)
                    //    Multiweight = SetMultiWeight(INCC_count, save_roh_positive);
                    
                    double avg_WNCC = sum_WNCC /(double)INCC_count;
                    double avg_INCC = sum_INCC /(double)INCC_count;
                    double avg_GNCC = sum_GNCC /(double)INCC_count;
                    if(avg_WNCC >= 0)
                        temp_rho = avg_WNCC*Multiweight;
                    else
                        temp_rho = avg_WNCC;
                    
                    /*
                     if(save_roh_positive.size() > 1)
                     {
                     double sum_save_po = 0;
                     double sum_variance = 0;
                     double variance = 0;
                     
                     for(int i = 0 ; i < save_roh_positive.size() ; i++)
                     sum_save_po += save_roh_positive[i];
                     double avg_save_po = sum_save_po/save_roh_positive.size();
                     
                     for(int i = 0 ; i < save_roh_positive.size() ; i++)
                     sum_variance += (avg_save_po - save_roh_positive[i])*(avg_save_po - save_roh_positive[i]);
                     
                     variance = sqrt(sum_variance/save_roh_positive.size());
                     
                     double var_weight = 0.1/(1.0+variance);
                     
                     var_weight = 0.0;
                     
                     temp_rho = temp_rho*(1.0 + var_weight);
                     }
                     */
                    //temp_rho = avg_INCC;
                    
                    /*
                     if(temp_nccresult_multi[AWNCC_id] < temp_rho)
                     {
                     temp_nccresult_multi[AWNCC_id] = temp_rho;
                     temp_nccresult_sec_multi[AWNCC_id] = -1.0;
                     nccresult_pairs[AWNCC_id].result2 = iter_height;
                     }*/
                    FindPeakNcc2(Pyramid_step, iteration, temp_rho, iter_height, check_rho, pre_rho, pre_rho_WNCC, WNCC_temp_rho, pre_height, direction, max_roh, nccresult[pt_index], temp_nccresult, temp_nccresult_sec);
                }
                
            }
        }
        
        GridPT3[pt_index].total_images = 1;
        
        GridPT3[pt_index].ncc_seleceted_pair = AWNCC_id;
        
        /*if(temp_nccresult > 1.0)
        {
            temp_nccresult = 1.0;
            temp_nccresult_sec = -1.0;
        }
        */
        nccresult[pt_index].result0 = DoubleToSignedChar_result(temp_nccresult);
        nccresult[pt_index].result1 = DoubleToSignedChar_result(temp_nccresult_sec);
    }
}

void AWNCC_MPs(ProInfo *proinfo, LevelInfo &rlevelinfo,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel, double *minmaxHeight, MultiMPs **multimps, vector<D3DPOINT> &MatchedPts_list_mps)
{
    long count_MPs = 0;
    
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    printf("MPs multi find peak pairs %d\t%f\n",rlevelinfo.pairinfo->NumberOfPairs,im_resolution);
    
    double min_bhratio = 1000000;
    double max_bhratio = -1000000;
    
    vector<double> weight_bhratio;
    for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        double temp_ratio = pow(1.0/rlevelinfo.pairinfo->BHratio[pair_number],1.0);
        if(min_bhratio > temp_ratio)
            min_bhratio = temp_ratio;
        
        if(max_bhratio < temp_ratio)
            max_bhratio = temp_ratio;
        
        weight_bhratio.push_back(temp_ratio);
        //printf("before weight_bhratio %d\t%f\n",pair_number,weight_bhratio[pair_number]);
    }
    
    double bhratio_interval = max_bhratio - min_bhratio;
    double bhratio_norm = 100.0;
    
    for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        weight_bhratio[pair_number] = (weight_bhratio[pair_number] - min_bhratio)/bhratio_interval*bhratio_norm;
        
        //printf("after weight_bhratio %d\t%f\t%f\n",pair_number,weight_bhratio[pair_number],bhratio_interval);
    }
    //exit(1);
    
    int AWNCC_id = rlevelinfo.pairinfo->NumberOfPairs;
    
    D3DPOINT* temp_points = (D3DPOINT*)calloc(sizeof(D3DPOINT),(long)Size_Grid2D.height*(long)Size_Grid2D.width);
    
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        
        D3DPOINT I_point;
        I_point.m_X = pts_col*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[0];
        I_point.m_Y = pts_row*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[1];
        I_point.m_Z = Nodata;
        I_point.flag = true;
        temp_points[pt_index] = I_point;
        
//        GridPT3[pt_index].height_counts = 0;
        GridPT3[pt_index].total_images = 0;
        
        double height_interval = (*rlevelinfo.grid_resolution)*10;
        
        int selected_pair = -1;
        
        if(GridPT3[pt_index].ncc_seleceted_pair == AWNCC_id)
            GridPT3[pt_index].ncc_seleceted_pair = selected_pair;
        
        double max_wncc = -10;
        double min_wncc = 10;
        double height_diff;
        double final_height;
        
        /*
        //count peaks
        int count_peaks = 0;
        int peak_pair = -1;
        for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
        {
            int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
            int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
            if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
            {
                if(multimps[pt_index][pair_number].check_matched)
                {
                    count_peaks++;
                    peak_pair = pair_number;
                }
            }
        }
        
        if(count_peaks == 0)
        {
            selected_pair = -1;
        }
        else if(count_peaks == 1)
        {
            selected_pair = peak_pair;
            
            D3DPOINT point;
            point.m_X = pts_col*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[0];
            point.m_Y = pts_row*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[1];
            point.m_Z = multimps[pt_index][selected_pair].peak_height;
            point.m_roh = multimps[pt_index][selected_pair].peak_roh;
            point.flag = false;
            
            temp_points[pt_index] = point;
            
            GridPT3[pt_index].Mean_ortho_ncc = GridPT3[pt_index].ortho_ncc[selected_pair];
            GridPT3[pt_index].ncc_seleceted_pair = selected_pair;
            
            printf("sole peak %d\t%d\t%f\n",selected_pair,multimps[pt_index][selected_pair].peak_roh,multimps[pt_index][selected_pair].peak_height);
        }
        else*/
        {
            //select max_ncc pair and selected_pair
            for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
            {
                int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                {
                    double pair_peak_roh = SignedCharToDouble_result(multimps[pt_index][pair_number].peak_roh)*weight_bhratio[pair_number];
                    if(pair_peak_roh > -1 && multimps[pt_index][pair_number].check_matched)
                    {
                        //printf("pair_peak_roh %d\t%f\t%f\n",pt_index,pair_peak_roh,multimps[pt_index][pair_number].peak_height);
                        GridPT3[pt_index].total_images++;
                        
                        if(Pyramid_step == 4 && iteration == 1)
                        {
                            if(pair_peak_roh > 0.0)
                            {
                                if(max_wncc < pair_peak_roh)
                                {
                                    max_wncc = pair_peak_roh;
                                    selected_pair = pair_number;
                                }
                                
                                if(min_wncc > pair_peak_roh)
                                {
                                    min_wncc = pair_peak_roh;
                                }
                            }
                        }
                        else
                        {
                            if(pair_peak_roh > 0.0)
                            {
                                if(max_wncc < pair_peak_roh)
                                {
                                    max_wncc = pair_peak_roh;
                                    selected_pair = pair_number;
                                }
                                
                                if(min_wncc > pair_peak_roh)
                                {
                                    min_wncc = pair_peak_roh;
                                }
                            }
                        }
                        //printf("pair_number %d\tpair_peak_roh %f\n",pair_number,pair_peak_roh);
                    }
                }
            }
        
            double wncc_interval = max_wncc - min_wncc;
            //select matched height from average WNCC and single WNCC, set minmax height with weighted method
            
            
            if(selected_pair > -1)
            {
                double min_Height = 9999;
                double max_Height = -9999;
                //printf("interval %f\t%f\tminmax %f\t%f\t%f\t%f\n",bhratio_interval,wncc_interval,min_bhratio,max_bhratio,min_wncc,max_wncc);
                if(wncc_interval < 0)
                    exit(1);
                
                final_height = multimps[pt_index][selected_pair].peak_height;
                
                vector<unsigned char>* save_pair = (vector<unsigned char>*)calloc(sizeof(vector<unsigned char>),(rlevelinfo.pairinfo->NumberOfPairs+1));
                vector<double>* save_height = (vector<double>*)calloc(sizeof(vector<double>),(rlevelinfo.pairinfo->NumberOfPairs+1));
                double* mid_H = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
                double* weight_height = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
                double mean_bhratio;
                
                //printf("mem allocate\n");
                double sum_weight_height = 0;
                double sum_weight = 0;
                int max_weight = 0.0;
                int final_selected_pair = selected_pair;
                int end_query_pair = rlevelinfo.pairinfo->NumberOfPairs;
                if(rlevelinfo.pairinfo->NumberOfPairs == 1)
                    end_query_pair = 1;
                
                for(int query_pair = 0 ; query_pair < end_query_pair ; query_pair++)
                {
                    //int query_pair = selected_pair;
                    weight_height[query_pair] = Nodata;
                    
                    save_pair[query_pair].clear();
                    save_height[query_pair].clear();
                    
                    if(multimps[pt_index][query_pair].check_matched && multimps[pt_index][selected_pair].check_matched)
                    {
                        double query_peak_roh = SignedCharToDouble_result(multimps[pt_index][query_pair].peak_roh);
                        
                        height_diff = fabs(multimps[pt_index][selected_pair].peak_height - multimps[pt_index][query_pair].peak_height);
                        if(height_diff < height_interval && query_peak_roh > 0.0)
                        {
                            double mean_bhratio = min_bhratio;
                            int count = 0;
                            bool check_query_pair = false;
                            
                            if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)
                            {
                                int reference_id_q = rlevelinfo.pairinfo->pairs[query_pair].m_X;
                                int ti_q = rlevelinfo.pairinfo->pairs[query_pair].m_Y;
                                check_query_pair = (reference_id_q == GridPT3[pt_index].selected_pair || ti_q == GridPT3[pt_index].selected_pair);
                            }
                            else
                                check_query_pair = true;
                                
                            if(check_query_pair)
                            {
                                for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                                {
                                    int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                                    int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                                    bool check_select_pair = (reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& (query_pair != pair_number)*/;
                                    
                                    if(check_select_pair)
                                    {
                                        if(multimps[pt_index][pair_number].check_matched)//Single peak
                                        {
                                            double pair_peak_roh = SignedCharToDouble_result(multimps[pt_index][pair_number].peak_roh);
                                            height_diff = fabs(multimps[pt_index][query_pair].peak_height - multimps[pt_index][pair_number].peak_height);
                                            if(height_diff < height_interval && pair_peak_roh > 0.0)
                                            {
                                                save_pair[query_pair].push_back(pair_number);
                                                save_height[query_pair].push_back(multimps[pt_index][pair_number].peak_height);
                                            }
                                        }
                                    }
                                }
                                
                                if(!rlevelinfo.check_SGM)
                                {
                                    //AWNCC add check
                                    if(multimps[pt_index][AWNCC_id].check_matched)
                                    {
                                        height_diff = fabs(multimps[pt_index][query_pair].peak_height - multimps[pt_index][AWNCC_id].peak_height);
                                        double awncc_peak_roh = SignedCharToDouble_result(multimps[pt_index][AWNCC_id].peak_roh);
                                        if(height_diff < height_interval && awncc_peak_roh > 0.0)
                                        {
                                            //AWNCC bhratio setting
                                            double sum_bhratio = 0;
                                            int count_pair = 0;
                                            for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                                            {
                                                int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                                                int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                                                bool check_select_pair = (reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& (query_pair != pair_number)*/;
                                                
                                                if(check_select_pair)
                                                {
                                                    if(multimps[pt_index][pair_number].check_matched)//Single peak
                                                    {
                                                        double pair_peak_roh = SignedCharToDouble_result(multimps[pt_index][pair_number].peak_roh);
                                                        height_diff = fabs(multimps[pt_index][AWNCC_id].peak_height - multimps[pt_index][pair_number].peak_height);
                                                        if(height_diff < height_interval && pair_peak_roh > 0.0)
                                                        {
                                                            sum_bhratio += pow(1.0/rlevelinfo.pairinfo->BHratio[pair_number],1.0);
                                                            count_pair++;
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            if(count_pair > 0)
                                            {
                                                save_pair[query_pair].push_back(AWNCC_id);
                                                save_height[query_pair].push_back(multimps[pt_index][AWNCC_id].peak_height);
                                                mean_bhratio = sum_bhratio/(double)count_pair;
                                            }
                                        }
                                    }
                                }
                            }
                            
                            if(save_pair[query_pair].size() > 0)
                            {
                                if(save_height[query_pair].size() > 1)
                                {
                                    mid_H[query_pair] = quickselect(save_height[query_pair],save_height[query_pair].size(),(int)(save_height[query_pair].size()/2.0));
                                }
                                else
                                    mid_H[query_pair] = save_height[query_pair][0];
                  
                                int mid_H_pair_number = save_pair[query_pair][(int)(save_height[query_pair].size()/2.0)];
                                double mid_H_wheight = mid_H[query_pair]*weight_bhratio[mid_H_pair_number];
                                
                                mid_H[query_pair] = multimps[pt_index][query_pair].peak_height;
                                mid_H_wheight = mid_H[query_pair];///(rlevelinfo.pairinfo->BHratio[query_pair]*100.0);
                                
                                double wheight_idw = 0;
                                double weight_idw = 0;
                                double wheight_bh = 0;
                                double weight_bh = 0;
                                double wheight_wncc = 0;
                                double weight_wncc = 0;
                                double wheight_awncc = 0;
                                double weight_awncc = 0;
                                
                                double total_weight = 0.0;
                                
                                for(int count = 0 ; count < save_pair[query_pair].size() ; count++)
                                {
                                    int pair_number = save_pair[query_pair][count];
                                    double pair_peak_roh = SignedCharToDouble_result(multimps[pt_index][pair_number].peak_roh);
                                    
                                    if(min_Height > save_height[query_pair][count])
                                        min_Height = save_height[query_pair][count];
                                    
                                    if(max_Height < save_height[query_pair][count])
                                        max_Height = save_height[query_pair][count];
                                    
                                    //height_diff = fabs(save_height[query_pair][count]*weight_bhratio[pair_number] - mid_H_wheight);
                                    
                                    height_diff = fabs(save_height[query_pair][count] - mid_H_wheight);
                                    
                                    double IDW_w;
                                    if(height_diff < 1.0)
                                        IDW_w = bhratio_norm;
                                    else
                                        IDW_w = bhratio_norm/pow(height_diff,1.0);
                                    
                                    double w_bhratio,w_ncc, weightAWNCC(1.0);
                                    
                                    
                                    if(bhratio_interval > 0)
                                    {
                                        if(pair_number < rlevelinfo.pairinfo->NumberOfPairs)//Single peak
                                            w_bhratio = weight_bhratio[pair_number];
                                        else//AWNCC peak
                                            w_bhratio = (mean_bhratio-min_bhratio)/bhratio_interval*bhratio_norm;
                                        
                                        if(w_bhratio == 0)
                                            w_bhratio = 0.001;
                                    }
                                    else
                                        w_bhratio = 1.0;
                                    
                                    double ncc_diff = (fabs(pair_peak_roh - query_peak_roh))*10.0;
                                    
                                    /*if(rlevelinfo.check_SGM)
                                    {
                                        weightAWNCC = 1.0;
                                        //w_ncc = 1.0;
                                    }
                                    else*/
                                    {
                                        /*
                                        if(ncc_diff < 1)
                                            w_ncc = 1.0;
                                        else
                                            w_ncc = 1.0/pow(ncc_diff,0.5);
                                        */
                                        if(pair_number == AWNCC_id)
                                            weightAWNCC = 1.0 + 0.1*save_pair[AWNCC_id].size();
                                        else
                                            weightAWNCC = 1.0;
                                    }
                                    
                                    if(ncc_diff < 1)
                                        w_ncc = bhratio_norm;
                                    else
                                        w_ncc = bhratio_norm/pow(ncc_diff,1.0);
                                    
                                    
                                    wheight_idw += multimps[pt_index][pair_number].peak_height*IDW_w*weightAWNCC;
                                    weight_idw += IDW_w*weightAWNCC;
                                    
                                    wheight_bh += multimps[pt_index][pair_number].peak_height*w_bhratio*weightAWNCC;
                                    weight_bh += w_bhratio*weightAWNCC;
                                    
                                    wheight_wncc += multimps[pt_index][pair_number].peak_height*w_ncc*weightAWNCC;
                                    weight_wncc += w_ncc*weightAWNCC;
                                    
                                    //wheight_awncc += multimps[pt_index][pair_number].peak_height*weightAWNCC;
                                    //weight_awncc += weightAWNCC;
                                    
                                    if(IDW_w > bhratio_norm || IDW_w < 0 || w_bhratio > bhratio_norm || w_bhratio < 0.0 || w_ncc > bhratio_norm || w_ncc < 0.0)
                                    {
                                        printf("weight %d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",query_pair,pair_number,IDW_w,w_bhratio,w_ncc,weightAWNCC,min_bhratio,max_bhratio,weight_bhratio[pair_number],mean_bhratio);
                                        exit(1);
                                    }
                                }
                                
                                total_weight = wheight_bh+weight_idw+weight_wncc;
                                total_weight = wheight_bh;
                                if(weight_bh > 0)
                                {
                                    if(rlevelinfo.check_SGM)
                                        weight_height[query_pair] = 0.40*wheight_idw/weight_idw + 0.60*wheight_bh/weight_bh;
                                    else
                                        weight_height[query_pair] = 0.30*wheight_idw/weight_idw + 0.40*wheight_bh/weight_bh + 0.30*wheight_wncc/weight_wncc;
                                }
                                else
                                    weight_height[query_pair] = 0.50*wheight_idw/weight_idw + 0.50*wheight_wncc/weight_wncc;
                                
                                double final_weight = 1.0;//query_peak_roh*total_weight;//*save_pair[query_pair].size();
                                if(max_weight < total_weight && query_pair < AWNCC_id)
                                {
                                    max_weight = total_weight;
                                    selected_pair = query_pair;
                                }
                                //if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)//Single peak
                                {
                                    sum_weight_height += weight_height[query_pair]*final_weight;
                                    sum_weight += final_weight;
                                }
                                
                                save_pair[query_pair].clear();
                                vector<unsigned char>().swap(save_pair[query_pair]);
                                save_height[query_pair].clear();
                                vector<double>().swap(save_height[query_pair]);
                             }
                        }
                    }
                }
                
                free(save_height);
                free(save_pair);
                free(mid_H);
                
                
                if(sum_weight > 0)
                {
                    final_height = sum_weight_height/sum_weight;
                    //final_height = weight_height[selected_pair];
                }
                free(weight_height);
                //if( final_height >= min_Height && final_height <= max_Height)
                {
                    D3DPOINT point;
                    point.m_X = pts_col*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[0];
                    point.m_Y = pts_row*(*rlevelinfo.grid_resolution)+rlevelinfo.Boundary[1];
                    point.m_Z = final_height;
                    point.m_roh = multimps[pt_index][selected_pair].peak_roh;
                    point.flag = false;
                    
                    temp_points[pt_index] = point;
                    
                    GridPT3[pt_index].Mean_ortho_ncc = GridPT3[pt_index].ortho_ncc[selected_pair];
                    GridPT3[pt_index].ncc_seleceted_pair = selected_pair;
                }
             }
        }
    }
    
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        if(!temp_points[pt_index].flag && temp_points[pt_index].m_Z > Nodata)
            MatchedPts_list_mps.push_back(temp_points[pt_index]);
    }
    
    free(temp_points);
}

void AWNCC_multi(ProInfo *proinfo, GridVoxel &grid_voxel,LevelInfo &rlevelinfo,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel, double *minmaxHeight)
{
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    printf("AWNCC multi find peak pairs %d\t%f\n",rlevelinfo.pairinfo->NumberOfPairs,im_resolution);
    
    bool pre_DEMtif = proinfo->pre_DEMtif;
    bool check_ortho = false;
    
    if(Pyramid_step == 4 && iteration == 1)
        check_ortho = false;
    else
        check_ortho = true;
    
    if(pre_DEMtif && !(Pyramid_step == 4 && iteration == 1))
        check_ortho = true;
    
    const double ncc_alpha = SetNCC_alpha(Pyramid_step,iteration, proinfo->IsRA);
    const double ncc_beta = 1.0 - ncc_alpha;
    const double ortho_th = 0.7 - (4 - Pyramid_step)*0.10;
    
    double min_bhratio = 100;
    double max_bhratio = -100;
    for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        if(min_bhratio > rlevelinfo.pairinfo->BHratio[pair_number])
            min_bhratio = rlevelinfo.pairinfo->BHratio[pair_number];
        
        if(max_bhratio < rlevelinfo.pairinfo->BHratio[pair_number])
            max_bhratio = rlevelinfo.pairinfo->BHratio[pair_number];
        
    }
    double bhratio_interval = max_bhratio - min_bhratio;
    double bhratio_norm = 1.0;
    double ncc_norm = 0.5;
    
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        
        nccresult[pt_index].result0 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result2 = Nodata;
        nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result3 = Nodata;
        nccresult[pt_index].result4 = 0;
    
//        GridPT3[pt_index].height_counts = 0;
        GridPT3[pt_index].total_images = 0;
        
        double db_GNCC, db_INCC;
        
        double *pre_rho_multi = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        double *pre_rho_WNCC_multi = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        double *max_roh_multi = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        double *temp_nccresult_multi = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        double *temp_nccresult_sec_multi = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        
        float *pre_height_multi = (float*)malloc(sizeof(float)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        int *direction_multi = (int*)malloc(sizeof(int)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        bool *check_rho_multi = (bool*)malloc(sizeof(bool)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        
        NCCresult *nccresult_pairs = (NCCresult*)malloc(sizeof(NCCresult)*(rlevelinfo.pairinfo->NumberOfPairs+1));
        
        for(int i = 0 ; i <= rlevelinfo.pairinfo->NumberOfPairs ; i++)
        {
            pre_rho_multi[i] = -1.0;
            pre_rho_WNCC_multi[i] = -1.0;
            
            max_roh_multi[i] = 0.0;
            temp_nccresult_multi[i] = -100.0;
            temp_nccresult_sec_multi[i] = -100.0;
            
            pre_height_multi[i] = 0.0;
            direction_multi[i] = 0;
            check_rho_multi[i] = false;
            
            nccresult_pairs[i].result0 = DoubleToSignedChar_result(-1.0);
            nccresult_pairs[i].result2 = Nodata;
            nccresult_pairs[i].result1 = DoubleToSignedChar_result(-1.0);
            nccresult_pairs[i].result3 = Nodata;
            nccresult_pairs[i].result4 = 0;
        }
        
        int AWNCC_id = rlevelinfo.pairinfo->NumberOfPairs;
        
        for(long height_step = 0 ; height_step < nccresult[pt_index].NumOfHeight ; height_step++)
        {
            float iter_height = nccresult[pt_index].minHeight + height_step*step_height;
            
            if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
            {
                double sum_INCC = 0;
                double sum_WNCC = 0;
                double sum_GNCC = 0;
                int INCC_count = 0;
                
                //averaged WNCC with multi pair
                vector<double> save_roh_positive;
                for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                {
                    int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                    int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                    if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                    {
                        db_GNCC = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[pair_number]);
                        db_INCC = SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, pair_number));
                        
                        double temp_rho = 0;
                        double WNCC_temp_rho = 0;
                        WNCC_temp_rho = db_INCC;
                        
                        if(grid_voxel[pt_index].flag_cal(height_step, pair_number) && db_INCC > -1)
                        {
                            double gncc_weight = SetGnccWeight(Pyramid_step, db_GNCC, db_INCC, GridPT3[pt_index].Height, iter_height, step_height);
                            
                            if((Pyramid_step == 4 && iteration == 1))
                                gncc_weight = 1.0;
                            
                            if(check_ortho) // GNCC check
                            {
                                if(db_GNCC > ortho_th)
                                    temp_rho = db_INCC*ncc_alpha + db_GNCC*ncc_beta;
                                else
                                    temp_rho = db_INCC*gncc_weight;
                            }
                            else
                            {
                                temp_rho = db_INCC*gncc_weight;
                            }
                            
                            if(temp_rho >= 0)
                                save_roh_positive.push_back(temp_rho);
                            
                            FindPeakNcc2(Pyramid_step, iteration, temp_rho, iter_height, check_rho_multi[pair_number], pre_rho_multi[pair_number], pre_rho_WNCC_multi[pair_number], WNCC_temp_rho, pre_height_multi[pair_number], direction_multi[pair_number], max_roh_multi[pair_number], nccresult_pairs[pair_number], temp_nccresult_multi[pair_number], temp_nccresult_sec_multi[pair_number]);
                            
                            if(bhratio_interval > 0)
                                sum_WNCC += temp_rho;//*(1.0 + (rlevelinfo.pairinfo->BHratio[pair_number] - min_bhratio)/bhratio_interval*0.2);
                            else
                                sum_WNCC += temp_rho;
                            
                            sum_INCC += db_INCC;
                            sum_GNCC += db_GNCC;
                            INCC_count++;
                        }
                    }
                }
                
                double temp_rho = -1;
                double WNCC_temp_rho = -1;
                if(INCC_count > 0)
                    WNCC_temp_rho = sum_INCC / (double)INCC_count;
                else
                    WNCC_temp_rho = -1;
                
                if(INCC_count > 0)
                {
                    double Multiweight = 1.0;
                    //if(INCC_count > 1)
                    //    Multiweight = SetMultiWeight(INCC_count, save_roh_positive);
                    
                    double avg_WNCC = sum_WNCC /(double)INCC_count;
                    double avg_INCC = sum_INCC /(double)INCC_count;
                    double avg_GNCC = sum_GNCC /(double)INCC_count;
                    if(avg_WNCC >= 0)
                        temp_rho = avg_WNCC*Multiweight;
                    else
                        temp_rho = avg_WNCC;
                    
                    /*
                    if(save_roh_positive.size() > 1)
                    {
                        double sum_save_po = 0;
                        double sum_variance = 0;
                        double variance = 0;
                        
                        for(int i = 0 ; i < save_roh_positive.size() ; i++)
                            sum_save_po += save_roh_positive[i];
                        double avg_save_po = sum_save_po/save_roh_positive.size();
                        
                        for(int i = 0 ; i < save_roh_positive.size() ; i++)
                            sum_variance += (avg_save_po - save_roh_positive[i])*(avg_save_po - save_roh_positive[i]);
                        
                        variance = sqrt(sum_variance/save_roh_positive.size());
                        
                        double var_weight = 0.1/(1.0+variance);
                        
                        var_weight = 0.0;
                        
                        temp_rho = temp_rho*(1.0 + var_weight);
                    }
                    */
                    //temp_rho = avg_INCC;
           
                    /*
                    if(temp_nccresult_multi[AWNCC_id] < temp_rho)
                    {
                        temp_nccresult_multi[AWNCC_id] = temp_rho;
                        temp_nccresult_sec_multi[AWNCC_id] = -1.0;
                        nccresult_pairs[AWNCC_id].result2 = iter_height;
                    }*/
                    FindPeakNcc2(Pyramid_step, iteration, temp_rho, iter_height, check_rho_multi[AWNCC_id], pre_rho_multi[AWNCC_id], pre_rho_WNCC_multi[AWNCC_id], WNCC_temp_rho, pre_height_multi[AWNCC_id], direction_multi[AWNCC_id], max_roh_multi[AWNCC_id], nccresult_pairs[AWNCC_id], temp_nccresult_multi[AWNCC_id], temp_nccresult_sec_multi[AWNCC_id]);
                }
                save_roh_positive.clear();
            }
        }
        
        double height_interval = (*rlevelinfo.grid_resolution)*10;
        
        int selected_pair = -1;
        int selected_ncc_pair = selected_pair;
        
        double max_wncc = -10;
        double min_wncc = 10;
        double max_gncc = -10;
        double height_diff;
        
        //select max_ncc pair and selected_pair
        for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
        {
            int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
            int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
            if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
            {
                if(temp_nccresult_multi[pair_number] > -1)
                {
                    GridPT3[pt_index].total_images++;
                    
                    if(Pyramid_step == 4 && iteration == 1)
                    {
                        if(temp_nccresult_multi[pair_number] > 0.0)
                        {
                            if(max_wncc < temp_nccresult_multi[pair_number])
                            {
                                max_wncc = temp_nccresult_multi[pair_number];
                                selected_pair = pair_number;
                            }
                            
                            if(min_wncc > temp_nccresult_multi[pair_number])
                            {
                                min_wncc = temp_nccresult_multi[pair_number];
                            }
                        }
                    }
                    else
                    {
                        if(temp_nccresult_multi[pair_number] > 0.0)
                        {
                            if(max_wncc < temp_nccresult_multi[pair_number])
                            {
                                max_wncc = temp_nccresult_multi[pair_number];
                                max_gncc = db_GNCC;
                                selected_pair = pair_number;
                            }
                            
                            if(min_wncc > temp_nccresult_multi[pair_number])
                            {
                                min_wncc = temp_nccresult_multi[pair_number];
                            }
                        }
                    }
                }
                //if(GridPT3[pt_index].selected_pair == 1 || GridPT3[pt_index].selected_pair == 4)
                //    printf("pair_number %d\t%d\t%d\n",pair_number,reference_id,ti);
            }
        }
        //if(GridPT3[pt_index].selected_pair == 1 || GridPT3[pt_index].selected_pair == 4)
        //    exit(1);
        if(temp_nccresult_multi[AWNCC_id] > 0.0)
        {
            if(temp_nccresult_multi[AWNCC_id] < min_wncc)
                min_wncc = temp_nccresult_multi[AWNCC_id];
            
            if(temp_nccresult_multi[AWNCC_id] > max_wncc)
                max_wncc = temp_nccresult_multi[AWNCC_id];
        }
        
        double wncc_interval = max_wncc - min_wncc;
        
        //select matched height from average WNCC and single WNCC, set minmax height with weighted method
        //printf("interval %f\t%f\tminmax %f\t%f\t%f\t%f\n",bhratio_interval,wncc_interval,min_bhratio,max_bhratio,min_wncc,max_wncc);
        if(selected_pair > -1)
        {
            vector<unsigned char>* save_pair = (vector<unsigned char>*)calloc(sizeof(vector<unsigned char>),(rlevelinfo.pairinfo->NumberOfPairs+1));
            vector<double>* save_height = (vector<double>*)calloc(sizeof(vector<double>),(rlevelinfo.pairinfo->NumberOfPairs+1));
            double* mid_H = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
            double* weight_height = (double*)malloc(sizeof(double)*(rlevelinfo.pairinfo->NumberOfPairs+1));
            double mean_bhratio;
            
            //printf("mem allocate\n");
            double sum_weight_height = 0;
            double sum_weight = 0;
            int max_weight = 0.0;
            int final_selected_pair = selected_pair;
            
            //for(int query_pair = 0 ; query_pair < rlevelinfo.pairinfo->NumberOfPairs+1 ; query_pair++)
            {
                int query_pair = selected_pair;
                weight_height[query_pair] = Nodata;
                
                save_pair[query_pair].clear();
                save_height[query_pair].clear();
                
                height_diff = fabs(nccresult_pairs[selected_pair].result2 - nccresult_pairs[query_pair].result2);
                if(height_diff < height_interval*2 && temp_nccresult_multi[query_pair] > 0.0)
                {
                    double sum_bhratio = 0;
                    int count = 0;
                    bool check_single_pair = false;
                    bool check_multi_pair = false;
                    
                    if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)
                    {
                        int reference_id_q = rlevelinfo.pairinfo->pairs[query_pair].m_X;
                        int ti_q = rlevelinfo.pairinfo->pairs[query_pair].m_Y;
                        if((reference_id_q == GridPT3[pt_index].selected_pair || ti_q == GridPT3[pt_index].selected_pair))
                        {
                            
                            for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                            {
                                int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                                int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                                if((reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) && (query_pair != pair_number))
                                {
                                    //printf("query_pair pair_number %d\t%d\tselected_pair %d\n",query_pair,pair_number,selected_pair);
                                    
                                    if(nccresult_pairs[query_pair].result2 > Nodata && nccresult_pairs[pair_number].result2 > Nodata)//Single peak
                                    {
                                        //printf("hd before\n");
                                        height_diff = fabs(nccresult_pairs[query_pair].result2 - nccresult_pairs[pair_number].result2);
                                        //printf("hd after %f\t%f\t%f\n",nccresult_pairs[query_pair].result2,nccresult_pairs[pair_number].result2,height_diff);
                                        if(height_diff < height_interval && temp_nccresult_multi[pair_number] > 0.0)
                                        {
                                            //printf("add before\n");
                                            save_pair[query_pair].push_back(pair_number);
                                            save_height[query_pair].push_back(nccresult_pairs[pair_number].result2);
                                            //printf("after before %f\t%f\n",nccresult_pairs[query_pair].result2,nccresult_pairs[pair_number].result2);
                                            
                                            check_single_pair = true;
                                            
                                            //printf("pair number %d\tsave pair %d\tsave_height %f\n",pair_number,save_pair[query_pair][count],save_height[query_pair][count]);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                        {
                            int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                            int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                            if((reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& (query_pair != pair_number)*/)
                            {
                                if(nccresult_pairs[AWNCC_id].result2 > Nodata && nccresult_pairs[pair_number].result2 > Nodata)//AWNCC peak
                                {
                                    height_diff = fabs(nccresult_pairs[AWNCC_id].result2 - nccresult_pairs[pair_number].result2);
                                    if(height_diff < height_interval && temp_nccresult_multi[pair_number] > 0.0)
                                    {
                                        save_pair[query_pair].push_back(pair_number);
                                        save_height[query_pair].push_back(nccresult_pairs[pair_number].result2);
                                        sum_bhratio += rlevelinfo.pairinfo->BHratio[pair_number];
                                        
                                        check_multi_pair = true;
                                        //printf("pair number %d\tsave pair %d\tsave_height %f\n",pair_number,save_pair[query_pair][count],save_height[query_pair][count]);
                                        count++;
                                    }
                                }
                                
                            }
                        }
                    }
                    
                    if(save_pair[query_pair].size() > 0)
                    {
                        /*
                         if(query_pair == rlevelinfo.pairinfo->NumberOfPairs)
                         {
                         printf("save_height[query_pair].size() %d\t%d\tmean_bhratio %f\n",count,save_height[query_pair].size(),mean_bhratio);
                         exit(1);
                         }
                         */
                        
                        /*if(save_height[query_pair].size() > 1)
                        {
                            
                            mid_H[query_pair] = quickselect(save_height[query_pair],save_height[query_pair].size(),(int)(save_height[query_pair].size()/2.0));
                        }
                        else
                            mid_H[query_pair] = save_height[query_pair][0];
                        */
                        mid_H[query_pair] = nccresult_pairs[query_pair].result2;
                        
                        save_height[query_pair].clear();
                        
                        double wheight_idw = 0;
                        double weight_idw = 0;
                        double wheight_bh = 0;
                        double weight_bh = 0;
                        double wheight_wncc = 0;
                        double weight_wncc = 0;
                        
                        double minH = 9999;
                        double maxH = -9999;
                        double total_weight = 0.0;
                        double sum_of_HD = 0;
                        double sum_of_ND = 0;
                        double sum_of_BH = 0;
                        
                        vector<double> save_height_diff, save_ncc_diff;
                        for(int count = 0 ; count < save_pair[query_pair].size() ; count++)
                        {
                            double IDW_w;
                            int pair_number = save_pair[query_pair][count];
                            
                            if(minH > nccresult_pairs[pair_number].result2)
                                minH = nccresult_pairs[pair_number].result2;
                            
                            if(maxH < nccresult_pairs[pair_number].result2)
                                maxH = nccresult_pairs[pair_number].result2;
                            
                            height_diff = fabs(nccresult_pairs[pair_number].result2 - mid_H[query_pair]);
                            sum_of_HD += height_diff;
                            save_height_diff.push_back(height_diff);
                            
                            
                            if(height_diff < 1.0)
                                IDW_w = 1.0;
                            else
                                IDW_w = 1.0/pow(height_diff,1.5);
                            
                            
                            double w_idw,w_bhratio,w_ncc;
                            
                            if(bhratio_interval > 0)
                            {
                                if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)//Single peak
                                    w_bhratio = (rlevelinfo.pairinfo->BHratio[pair_number]-min_bhratio)/bhratio_interval*bhratio_norm;
                                else//AWNCC peak
                                {
                                    mean_bhratio = sum_bhratio/save_pair[query_pair].size();
                                    w_bhratio = (mean_bhratio-min_bhratio)/bhratio_interval*bhratio_norm;
                                }
                            }
                            else
                                w_bhratio = 1.0;
                            
                            //w_bhratio = rlevelinfo.pairinfo->BHratio[pair_number];
                            sum_of_BH += w_bhratio;
                            
                            double ncc_diff = (fabs(temp_nccresult_multi[pair_number] - temp_nccresult_multi[query_pair]))*10.0;
                            sum_of_ND += ncc_diff;
                            save_ncc_diff.push_back(ncc_diff);
                            
                            //printf("query_pair %d\tpair_number %d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",query_pair,pair_number,w_bhratio,nccresult_pairs[pair_number].result2,mid_H[query_pair],save_height_diff[count],temp_nccresult_multi[pair_number],temp_nccresult_multi[query_pair],save_ncc_diff[count]);
                            
                            if(ncc_diff > 1)
                                w_ncc = 1.0/pow(ncc_diff,1.5);
                            else
                                w_ncc = 1.0;
                            /*
                            if(wncc_interval > 0)
                                w_ncc = 1.0 + (temp_nccresult_multi[pair_number]-min_wncc)/wncc_interval*ncc_norm;
                            else
                                w_ncc = 1.0;
                            */
                            //total_weight += IDW_w*w_bhratio*w_ncc;
                            //total_weight += (0.5*IDW_w + 0.2*w_bhratio + 0.3*w_ncc);
                            //if(save_pair.size() > 2)
                            //    printf("single pt_index %d\tpairindex %d\t%f\tBHratio %f\t%f\t%f\t%f\t%f\t%f\t%f\n",pt_index,pair_number,sum_bhratio,rlevelinfo.pairinfo->BHratio[pair_number],mid_H,IDW_w,nccresult_pairs[pair_number].result2,w_bhratio,w_ncc,total_weight);
                            /*if(w_bhratio < 1 || w_ncc < 1)
                            {
                                printf("%d\tIDW_w %f\t%f\t%f\t%f\n",pt_index,total_weight,IDW_w,w_bhratio,w_ncc);
                                exit(1);
                            }*/
                            
                            
                            wheight_idw += nccresult_pairs[pair_number].result2*IDW_w;
                            weight_idw += IDW_w;
                            
                            wheight_bh += nccresult_pairs[pair_number].result2*w_bhratio;
                            weight_bh += w_bhratio;
                            
                            wheight_wncc += nccresult_pairs[pair_number].result2*w_ncc;
                            weight_wncc += w_ncc;
                            
                            
                            //wheight_awncc += nccresult_pairs[pair_number].result2*weightAWNCC;
                            //weight_awncc += weightAWNCC;
                        }
                        /*
                        double av_HD = sum_of_HD/save_pair[query_pair].size();
                        double av_ND = sum_of_ND/save_pair[query_pair].size();
                        //printf("average %f\t%f\n",av_HD,av_ND);
                        
                        sum_of_HD = 0;
                        sum_of_ND = 0;
                        for(int count = 0 ; count < save_pair[query_pair].size() ; count++)
                        {
                            int pair_number = save_pair[query_pair][count];
                            //printf("pair_number %d\t%f\t%f\t%f\t%f\n",pair_number,save_height_diff[count],save_ncc_diff[count],sum_of_HD,sum_of_ND);
                            sum_of_HD += ( save_height_diff[count] - av_HD)*( save_height_diff[count] - av_HD);
                            sum_of_ND += ( save_ncc_diff[count] - av_ND)*( save_ncc_diff[count] - av_ND);
                        }
                        
                        wheight_bh = sum_of_BH/save_pair[query_pair].size();
                        double var = sqrt(sum_of_HD/save_pair[query_pair].size());
                        if(var > 1)
                            weight_idw = 1.0/pow(var,0.2);
                        else
                            weight_idw = 1.0;
                        
                        var = sqrt(sum_of_ND/save_pair[query_pair].size());
                        if(var > 1)
                            weight_wncc = 1.0/var;
                        else
                            weight_wncc = 1.0;
                        
                        total_weight = wheight_bh*weight_idw*weight_wncc;//*save_pair[query_pair].size();
                        */
                         
                        //if(query_pair == AWNCC_id)
                        //    total_weight = total_weight*2.0;
                        //printf("total_weight %f\t%f\t%f\t%f\t%d\n",total_weight,wheight_bh,weight_idw,weight_wncc,save_pair[query_pair].size());
                        //if(save_pair[query_pair].size() > 3)
                        //    exit(1);
                        weight_height[query_pair] = 0.3*wheight_idw/weight_idw + 0.3*wheight_bh/weight_bh + 0.4*wheight_wncc/weight_wncc;// + 0.4*wheight_awncc/weight_awncc;
                        
                        //weight_height[query_pair] = nccresult_pairs[query_pair].result2;
                        
                        /*
                         if(weight_height >= minH -1 && weight_height <= maxH +1)
                         nccresult_pairs[selected_pair].result2 = weight_height;
                         else
                         {
                         
                         for(int count = 0 ; count < save_pair.size() ; count++)
                         {
                         int pair_number = save_pair[count];
                         printf("pt_index %d\tinterval %f\t%f\tminmax %f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n",pt_index,bhratio_interval,wncc_interval,min_bhratio,max_bhratio,min_wncc,max_wncc,count,pair_number,nccresult_pairs[pair_number].result2,minH,maxH,weight_height,mid_H);
                         }
                         //nccresult_pairs[selected_pair].result2 = mid_H;
                         exit(1);
                         }
                         */
                        //GridPT3[pt_index].height_counts = save_pair.size();
                        
                        //if(save_pair.size() > 2)
                        //    exit(1);
                        /*
                         if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)//Single peak
                         {
                         sum_weight_height += weight_height[query_pair]*pow(save_pair[query_pair].size(),2.0)*temp_nccresult_multi[query_pair];
                         sum_weight += pow(save_pair[query_pair].size(),2.0)*temp_nccresult_multi[query_pair];
                         }
                         else
                         {
                         sum_weight_height += weight_height[query_pair]*pow(save_pair[query_pair].size(),4.0)*temp_nccresult_multi[AWNCC_id];
                         sum_weight += pow(save_pair[query_pair].size(),4.0)*temp_nccresult_multi[AWNCC_id];
                         }
                         */
                        double final_weight = 1.0;//temp_nccresult_multi[query_pair]*total_weight*save_pair[query_pair].size();
                        /*if(max_weight < final_weight)
                         {
                         max_weight = final_weight;
                         final_selected_pair = query_pair;
                         }*/
                        //if(query_pair < rlevelinfo.pairinfo->NumberOfPairs)//Single peak
                        {
                            sum_weight_height += weight_height[query_pair]*final_weight;
                            sum_weight += final_weight;
                        }
                        /*else
                         {
                         sum_weight_height += weight_height[query_pair]*final_weight;
                         sum_weight += final_weight;
                         }
                         */
                        //printf("pt_index %d\tquery_pair %d\t%f\t%f\t%f\t%d\n",pt_index,query_pair,sum_weight_height,sum_weight,weight_height[query_pair],save_pair[query_pair].size());
                        save_pair[query_pair].clear();
                    }
                }
            }
            
            free(save_height);
            free(save_pair);
            free(mid_H);
            free(weight_height);
            
            //if(final_selected_pair != AWNCC_id)
            //    selected_pair = final_selected_pair;
            if(sum_weight > 0)
                nccresult_pairs[selected_pair].result2 = sum_weight_height/sum_weight;
            //else
            //    printf("sum_weight %f\t%f\n",sum_weight_height,sum_weight);
            //printf("pt_index %d\t%f\n",pt_index,nccresult_pairs[selected_pair].result2);
            //exit(1);
            /*
            if(weight_height_awncc > -100 && weight_height_single > -100)
            {
                //printf("pt_index %d\tweightH %f\t%f\n",pt_index,weight_height_awncc,weight_height_single);
                nccresult_pairs[selected_pair].result2 = weight_height_awncc*0.5 + weight_height_single*0.5;
                //exit(1);
            }
            else if(weight_height_awncc > -100)
                nccresult_pairs[selected_pair].result2 = weight_height_awncc;
            else if(weight_height_single > -100)
                nccresult_pairs[selected_pair].result2 = weight_height_single;
             */
        }
        
        
        GridPT3[pt_index].ncc_seleceted_pair = selected_pair;
        
        if(selected_pair > -1)
        {
            double height_sum_1 = 0;
            double height_sum_2 = 0;
            double db_gncc_sum = 0;
            double db_gncc_sum_2 = 0;
            
            short minH = 10000;
            short maxH = -1000;
            
            if(temp_nccresult_multi[selected_pair] > 1.0)
            {
                temp_nccresult_multi[selected_pair] = 1.0;
                temp_nccresult_sec_multi[selected_pair] = -1.0;
            }
            
            nccresult[pt_index].result0 = DoubleToSignedChar_result(temp_nccresult_multi[selected_pair]);
            nccresult[pt_index].result1 = DoubleToSignedChar_result(temp_nccresult_sec_multi[selected_pair]);
            
            nccresult[pt_index].result4 = nccresult_pairs[selected_pair].result4;
            nccresult[pt_index].result2 = nccresult_pairs[selected_pair].result2;
            nccresult[pt_index].result3 = nccresult_pairs[selected_pair].result3;
        }
        else
        {
            nccresult[pt_index].result0 = DoubleToSignedChar_result(-1);
            nccresult[pt_index].result1 = DoubleToSignedChar_result(-1);
            
            nccresult[pt_index].result2 = 0;
            nccresult[pt_index].result3 = 0;
            
            nccresult[pt_index].result4 = 0;
        }
        
        free(pre_rho_multi);
        free(pre_rho_WNCC_multi);
        free(max_roh_multi);
        free(temp_nccresult_multi);
        free(temp_nccresult_sec_multi);
        
        free(pre_height_multi);
        free(direction_multi);
        free(check_rho_multi);
        
        free(nccresult_pairs);
    }
}


void AWNCC_SGM(ProInfo *proinfo, GridVoxel &grid_voxel,LevelInfo &rlevelinfo,CSize Size_Grid2D, UGRID *GridPT3, NCCresult *nccresult, double step_height, uint8 Pyramid_step, uint8 iteration,int MaxNumberofHeightVoxel, double *minmaxHeight, const int pairnumber)
{
    // P2 >= P1
    const double P1 = 0.3;
    const double P2 = 0.6;
    
    const int P_HS_step = 1;
    float **SumCost = NULL;
    
    bool check_diagonal = true;
    
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    
    long total_grid_size = (long)Size_Grid2D.width*(long)Size_Grid2D.height;
    SumCost = (float**)calloc(sizeof(float*),total_grid_size);
    for(long i=0;i<Size_Grid2D.height;i++)
    {
        for(long j=0;j<Size_Grid2D.width;j++)
        {
            long t_index = i*(long)Size_Grid2D.width + j;
            if(Pyramid_step == 0 && iteration == 3 && nccresult[t_index].NumOfHeight > 2000)
                printf("gridsize %d\t%d\t pos %d\t%d\t numofheight %d\t%d\t%d\n",Size_Grid2D.width,Size_Grid2D.height,j,i,nccresult[t_index].NumOfHeight,nccresult[t_index].maxHeight,nccresult[t_index].minHeight);
            if(nccresult[t_index].NumOfHeight > 0)
            {
                SumCost[t_index] = (float*)calloc(sizeof(float),nccresult[t_index].NumOfHeight);
            }
        }
    }
    
    //left , right, top, bottom, upper left, upper right, bottom left, bottom right
    int v_row[8]    = { 0, 0, -1, 1, -1, -1 ,  1, 1};
    int u_col[8]    = {-1, 1,  0, 0, -1,  1 , -1, 1};
    
    int row_iter[8] = { 1, 1, 1, -1,  1,  1, -1, -1};
    int col_iter[8] = { 1,-1, 1,  1,  1, -1,  1, -1};
    
    int start_row[8]    = {                         0,                        0,                        0, (int)Size_Grid2D.height-1,                        0,                        0, (int) Size_Grid2D.height-1 , (int)Size_Grid2D.height-1};
    int end_row[8]      = {(int)Size_Grid2D.height   , (int)Size_Grid2D.height , (int)Size_Grid2D.height ,                         0, (int)Size_Grid2D.height , (int)Size_Grid2D.height ,                          0 ,                         0};
    int start_col[8]    = {                         0, (int)Size_Grid2D.width-1,                        0,                         0,                        0, (int)Size_Grid2D.width-1,                          0 , (int)Size_Grid2D.width-1 };
    int end_col[8]      = {(int)Size_Grid2D.width    ,                        0, (int)Size_Grid2D.width  , (int)Size_Grid2D.width   , (int)Size_Grid2D.width  ,                        0, (int) Size_Grid2D.width    ,                         0};

    int direction_iter = 0;

    short maxHeight = 0;
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
      if(maxHeight < nccresult[iter_count].NumOfHeight)
        maxHeight = nccresult[iter_count].NumOfHeight;
    
    //4 directional
    {
#pragma omp parallel
      {
        
        float *LHcost_pre = (float*)calloc(sizeof(float), maxHeight);
        float *LHcost_curr = (float*)calloc(sizeof(float), maxHeight);
        float *temp;

        //printf("left to right\n");
#pragma omp for
        for(long pts_row = start_row[direction_iter] ; pts_row < end_row[direction_iter] ; pts_row = pts_row + row_iter[direction_iter])
        {
            for(long pts_col = start_col[direction_iter] ; pts_col < end_col[direction_iter] ; pts_col = pts_col + col_iter[direction_iter])
            {
                long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
                
                if(pts_col == start_col[direction_iter])
                {
                    memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_start_pos(nccresult, grid_voxel,rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                }
                else
                {
                    memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                    SWAP(LHcost_pre, LHcost_curr);
                    
                }
            }
        }
        
        //printf("right to left\n");
        direction_iter = 1;
#pragma omp for
        for(long pts_row = start_row[direction_iter] ; pts_row < end_row[direction_iter] ; pts_row = pts_row + row_iter[direction_iter])
        {
            for(long pts_col = start_col[direction_iter] ; pts_col >= end_col[direction_iter] ; pts_col = pts_col + col_iter[direction_iter])
            {
                long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
                
                if(pts_col == start_col[direction_iter])
                {
                    memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                }
                else
                {
                    memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                    SWAP(LHcost_pre, LHcost_curr);
                    
                }
            }
        }
        
        //printf("top to bottom\n");
        direction_iter = 2;
        
#pragma omp for
        for(long pts_col = start_col[direction_iter] ; pts_col < end_col[direction_iter] ; pts_col = pts_col + col_iter[direction_iter])
        {
            for(long pts_row = start_row[direction_iter] ; pts_row < end_row[direction_iter] ; pts_row = pts_row + row_iter[direction_iter])
            {
                long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
                
                if(pts_row == start_row[direction_iter])
                {
                    memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                }
                else
                {
                    memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                    SWAP(LHcost_pre, LHcost_curr);
                }
            }
        }
        
        //printf("bottom to top\n");
        direction_iter = 3;
#pragma omp for
        for(long pts_col = start_col[direction_iter] ; pts_col < end_col[direction_iter] ; pts_col = pts_col + col_iter[direction_iter])
        {
            for(long pts_row = start_row[direction_iter] ; pts_row >= end_row[direction_iter] ; pts_row = pts_row + row_iter[direction_iter])
            {
                long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
                
                if(pts_row == start_row[direction_iter])
                {
                    memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                }
                else
                {
                    memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                    SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                    SWAP(LHcost_pre, LHcost_curr);
                }
            }
        }
      }
    }
    
    //8 directional
    if(check_diagonal)
    {
#pragma omp parallel
      {
        float *LHcost_pre = (float*)calloc(sizeof(float), maxHeight);
        float *LHcost_curr = (float*)calloc(sizeof(float), maxHeight);
        float *temp;
        long ref_iter;
        long pts_row, pts_col;
        {
            //printf("upper left to right\n");
            direction_iter = 4;
            
#pragma omp for
            for(ref_iter = start_row[direction_iter] ; ref_iter < end_row[direction_iter] ; ref_iter = ref_iter + row_iter[direction_iter]) //left wall row direction
            {
                pts_row = ref_iter;
                bool check_end = false;
                pts_col = start_col[direction_iter];
         
                while(pts_col < end_col[direction_iter] && !check_end)
                {
                    double WNCC_sum = 0;
         
                    if(pts_col == start_col[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_row = pts_row - v_row[direction_iter];
         
                        if(pts_row >= 0 && pts_row < Size_Grid2D.height)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_col = pts_col - u_col[direction_iter];
               }
            }
         
#pragma omp for
            for(ref_iter = start_col[direction_iter] ; ref_iter < end_col[direction_iter] ; ref_iter = ref_iter + col_iter[direction_iter]) //top wall col direction
            {
                pts_col = ref_iter;
                bool check_end = false;
                pts_row = start_row[direction_iter];
         
                while(pts_row < end_row[direction_iter] && !check_end)
                {
                    double WNCC_sum = 0;
         
                    if(pts_row == start_row[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_col = pts_col - u_col[direction_iter];
         
                        if(pts_col >= 0 && pts_col < Size_Grid2D.width)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_row = pts_row - v_row[direction_iter];
                }
            }
        }
        
        {
            //printf("upper right to left\n");
            direction_iter = 5;
#pragma omp for
            for(ref_iter = start_row[direction_iter] ; ref_iter < end_row[direction_iter] ; ref_iter = ref_iter + row_iter[direction_iter]) //left wall row direction
            {
                pts_row = ref_iter;
                bool check_end = false;
                pts_col = start_col[direction_iter];
         
                while(pts_col >= end_col[direction_iter] && !check_end)
                {
                    double WNCC_sum = 0;
         
                    if(pts_col == start_col[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_row = pts_row - v_row[direction_iter];
         
                        if(pts_row >= 0 && pts_row < Size_Grid2D.height)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_col = pts_col - u_col[direction_iter];
                }
            }
         
#pragma omp for
            for(ref_iter = start_col[direction_iter] ; ref_iter >= end_col[direction_iter] ; ref_iter = ref_iter + col_iter[direction_iter]) //top wall col direction
            {
                pts_col = ref_iter;
                bool check_end = false;
                pts_row = start_row[direction_iter];
         
                while(pts_row < end_row[direction_iter] && !check_end)
                {
                    double WNCC_sum = 0;
         
                    if(pts_row == start_row[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_col = pts_col - u_col[direction_iter];
         
                        if(pts_col >= 0 && pts_col < Size_Grid2D.width)
                        {
                            long pt_index = pts_row*Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_row = pts_row - v_row[direction_iter];
                }
            }
         
            //printf("bottom left to right\n");
            direction_iter = 6;
#pragma omp for
            for(ref_iter = start_row[direction_iter] ; ref_iter >= end_row[direction_iter] ; ref_iter = ref_iter + row_iter[direction_iter]) //left wall row direction
            {
                pts_row = ref_iter;
                bool check_end = false;
                pts_col = start_col[direction_iter];
         
                while(pts_col < end_col[direction_iter] && !check_end)
                {
                    if(pts_col == start_col[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_row = pts_row - v_row[direction_iter];
         
                        if(pts_row >= 0 && pts_row < Size_Grid2D.height)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_col = pts_col - u_col[direction_iter];
                }
            }
         
#pragma omp for
            for(ref_iter = start_col[direction_iter] ; ref_iter < end_col[direction_iter] ; ref_iter = ref_iter + col_iter[direction_iter]) //top wall col direction
            {
                pts_col = ref_iter;
                bool check_end = false;
                pts_row = start_row[direction_iter];
         
                while(pts_row >= end_row[direction_iter] && !check_end)
                {
                    if(pts_row == start_row[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_col = pts_col - u_col[direction_iter];
         
                        if(pts_col >= 0 && pts_col < Size_Grid2D.width)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_row = pts_row - v_row[direction_iter];
                }
            }
         
            //printf("bottom right to left\n");
            direction_iter = 7;
#pragma omp for
            for(ref_iter = start_row[direction_iter] ; ref_iter >= end_row[direction_iter] ; ref_iter = ref_iter + row_iter[direction_iter]) //left wall row direction
            {
                pts_row = ref_iter;
                bool check_end = false;
                pts_col = start_col[direction_iter];
         
                while(pts_col >= end_col[direction_iter] && !check_end)
                {
                    if(pts_col == start_col[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_row = pts_row - v_row[direction_iter];
         
                        if(pts_row >= 0 && pts_row < Size_Grid2D.height)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_col = pts_col - u_col[direction_iter];
                }
            }
         
#pragma omp for
            for(ref_iter = start_col[direction_iter] ; ref_iter >= end_col[direction_iter] ; ref_iter = ref_iter + col_iter[direction_iter]) //top wall col direction
            {
                pts_col = ref_iter;
                bool check_end = false;
                pts_row = start_row[direction_iter];
         
                while(pts_row >= end_row[direction_iter] && !check_end)
                {
                    if(pts_row == start_row[direction_iter])
                    {
                        long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                        memset(LHcost_pre, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                        SGM_start_pos(nccresult, grid_voxel, rlevelinfo, GridPT3, pt_index, LHcost_pre, SumCost, step_height, pairnumber);
                    }
                    else
                    {
                        pts_col = pts_col - u_col[direction_iter];
         
                        if(pts_col >= 0 && pts_col < Size_Grid2D.width)
                        {
                            long pt_index = pts_row*(long)Size_Grid2D.width + pts_col;
         
                            memset(LHcost_curr, 0, nccresult[pt_index].NumOfHeight*sizeof(float));
                            SGM_con_pos(pts_col, pts_row, Size_Grid2D, direction_iter, step_height, P_HS_step, u_col, v_row, nccresult, grid_voxel, GridPT3, rlevelinfo, pt_index, P1, P2, LHcost_pre, LHcost_curr, SumCost, pairnumber);
                            SWAP(LHcost_pre, LHcost_curr);
                        }
                        else
                            check_end = true;
                    }
                    pts_row = pts_row - v_row[direction_iter];
                }
            }
        }
      }
    }
   
    printf("SGM find peak pairs %d\t%f\n",pairnumber,im_resolution);
    
    
    bool check_ortho = false;
    
    if(Pyramid_step == 4 && iteration == 1)
        check_ortho = false;
    else
        check_ortho = true;
    
    if(proinfo->pre_DEMtif && !(Pyramid_step == 4 && iteration == 1))
        check_ortho = true;
    
    const double ncc_alpha = SetNCC_alpha(Pyramid_step,iteration, proinfo->IsRA);
    const double ncc_beta = 1.0 - ncc_alpha;
    const double ortho_th = 0.7 - (4 - Pyramid_step)*0.10;
#pragma omp parallel for schedule(guided)
    for(long iter_count = 0 ; iter_count < (long)Size_Grid2D.height*(long)Size_Grid2D.width ; iter_count++)
    {
        long pts_row = (int)(floor(iter_count/Size_Grid2D.width));
        long pts_col = iter_count % Size_Grid2D.width;
        long pt_index = iter_count;
        
        double pre_rho  = -1.0;
        double pre_rho_WNCC = -1.0;
        float pre_height= 0.0;
        int pre_height_step = 0;
        int direction   = 0;
        bool check_rho  = false;
        
        double max_roh = 0;
        double max_roh_sec = 0;
        
        double temp_nccresult = -100.0;
        double temp_nccresult_sec = -100.0;
        
        nccresult[pt_index].result0 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result2 = -1000;
        nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
        nccresult[pt_index].result3 = -1000;
        nccresult[pt_index].result4 = 0;
    
        double db_INCC, db_GNCC;
        
        for(long height_step = 0 ; height_step < nccresult[pt_index].NumOfHeight ; height_step++)
        {
            float iter_height = nccresult[pt_index].minHeight + height_step*step_height;
            
            if(iter_height >= GridPT3[pt_index].minHeight && iter_height <= GridPT3[pt_index].maxHeight)
            {
                double sum_INCC = 0;
                int INCC_count = 0;
                
                int reference_id = rlevelinfo.pairinfo->pairs[pairnumber].m_X;
                int ti = rlevelinfo.pairinfo->pairs[pairnumber].m_Y;
                if(reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair)
                {
                    db_GNCC = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[pairnumber]);
                    db_INCC = SignedCharToDouble_voxel(grid_voxel[pt_index].INCC(height_step, pairnumber));
                    
                    if(grid_voxel[pt_index].flag_cal(height_step, pairnumber) && db_INCC > -1)
                    {
                        //double gncc_weight = SetGnccWeight(Pyramid_step, db_GNCC, db_INCC, GridPT3[pt_index].Height, iter_height, step_height);
                        
                        //if((Pyramid_step == 4 && iteration == 1))
                        //    gncc_weight = 1.0;
                        
                        if(check_ortho) // GNCC check
                        {
                            if(db_GNCC > ortho_th)
                                sum_INCC = db_INCC*ncc_alpha + db_GNCC*ncc_beta;
                            else
                                sum_INCC = db_INCC;//*gncc_weight;
                            INCC_count++;
                        }
                        else
                        {
                            sum_INCC = db_INCC;//*gncc_weight;
                            INCC_count++;
                        }
                        sum_INCC = db_INCC;//*gncc_weight;
                        INCC_count++;
                    }
                }
                
                double temp_rho = -1;
                double WNCC_temp_rho = -1;
                if(INCC_count > 0)
                    WNCC_temp_rho = sum_INCC / (double)INCC_count;
                else
                    WNCC_temp_rho = -1;
                
                if(INCC_count > 0)
                {
                    temp_rho = SumCost[pt_index][height_step];

                    FindPeakNcc_SGM(Pyramid_step, iteration, temp_rho, iter_height, check_rho, pre_rho, pre_rho_WNCC, WNCC_temp_rho, pre_height, direction, max_roh, max_roh_sec, nccresult[pt_index], temp_nccresult, temp_nccresult_sec);
                }
            }
        }
        
        if(max_roh > 0 && temp_nccresult > -100)
        {
            if(Pyramid_step == 0 && iteration >= 2)
            {
                nccresult[pt_index].result0 = DoubleToSignedChar_result(max_roh);
                nccresult[pt_index].result1 = DoubleToSignedChar_result(-1.0);
                nccresult[pt_index].result3 = Nodata;
            }
            else
            {
                if(temp_nccresult < temp_nccresult_sec)
                    printf("ori peak 1 2 %f\t%f\n",temp_nccresult,temp_nccresult_sec);
                
                if(fabs(temp_nccresult) > 30.0)
                {
                    if(temp_nccresult > 30)
                        temp_nccresult = 30.0;
                    else if(temp_nccresult < -30.0)
                        temp_nccresult = -30.0;
                }
                
                if(fabs(temp_nccresult_sec) > 30.0)
                {
                    if(temp_nccresult_sec > 30)
                        temp_nccresult_sec = 30.0;
                    else if(temp_nccresult_sec < -30)
                        temp_nccresult_sec = -30.0;
                }
                
                nccresult[pt_index].result0 = DoubleToSignedChar_result(temp_nccresult);
                nccresult[pt_index].result1 = DoubleToSignedChar_result(temp_nccresult_sec);
                
                if(nccresult[pt_index].result0 < nccresult[pt_index].result1)
                    printf("after peak 1 2 %d\t%d\n",nccresult[pt_index].result0,nccresult[pt_index].result1);
                
            }
        }
    }
    
    for(int i=0;i<Size_Grid2D.height;i++)
    {
        for(int j=0;j<Size_Grid2D.width;j++)
        {
            long t_index = (long)(i*Size_Grid2D.width) + (long)j;
            if(nccresult[t_index].NumOfHeight > 0)
                free(SumCost[t_index]);
        }
    }
    free(SumCost);
    printf("done SumCost free\n");
    //if(pairnumber == 2)
    //    exit(1);
}

void VerticalLineLocus_seeddem(const ProInfo *proinfo,LevelInfo &rlevelinfo, UGRID *GridPT3, const double* minmaxHeight)
{
    const int Pyramid_step = *(rlevelinfo.Pyramid_step);
    int Template_size = *(rlevelinfo.Template_size);
    double template_area = 5.0;
    int t_Template_size = (int)((template_area/(proinfo->resolution*pwrtwo(Pyramid_step)))/2.0)*2+1;
    if(Template_size < t_Template_size)
        Template_size = t_Template_size;

    int Half_template_size = (int)(Template_size/2.0);
    const double subBoundary[4] = {rlevelinfo.Boundary[0], rlevelinfo.Boundary[1], rlevelinfo.Boundary[2], rlevelinfo.Boundary[3]};
    
    double im_resolution = proinfo->resolution*pwrtwo(Pyramid_step);
    double im_resolution_next = proinfo->resolution;
    D2DPOINT **all_im_cd = NULL;
    D2DPOINT **all_im_cd_next = NULL;
    
    long int sub_imagesize_w, sub_imagesize_h;
    long int sub_imagesize_w_next, sub_imagesize_h_next;
    
    all_im_cd = (D2DPOINT**)malloc(sizeof(D2DPOINT*)*proinfo->number_of_images);
    if(proinfo->check_Matchtag)
        SetOrthoImageCoord(proinfo, rlevelinfo, GridPT3, 0, BD, im_resolution, im_resolution_next, sub_imagesize_w, sub_imagesize_h, sub_imagesize_w_next, sub_imagesize_h_next, all_im_cd, all_im_cd_next);
    else
        SetOrthoImageCoord(proinfo, rlevelinfo, GridPT3, 0, OR, im_resolution, im_resolution_next, sub_imagesize_w, sub_imagesize_h, sub_imagesize_w_next, sub_imagesize_h_next, all_im_cd, all_im_cd_next);
    
    
    const int reference_id = rlevelinfo.reference_id;
    //int count_total = 0;
    //int count_low = 0;
#pragma omp parallel
    {
        SetKernel rsetkernel(reference_id,1,Half_template_size);
        
#pragma omp for /*reduction(+:count_low, count_total)*/ schedule(guided)
        for(long int iter_count = 0 ; iter_count < *rlevelinfo.Grid_length ; iter_count++)
        {
            long int pts_row = (long int)(floor(iter_count/rlevelinfo.Size_Grid2D->width));
            long int pts_col = iter_count % rlevelinfo.Size_Grid2D->width;
            long int pt_index = iter_count;//pts_row*(long int)rlevelinfo.Size_Grid2D->width + pts_col;

            if(pt_index < *rlevelinfo.Grid_length && pts_row < rlevelinfo.Size_Grid2D->height && pts_col < rlevelinfo.Size_Grid2D->width && pts_row >= 0 && pts_col >= 0)
            {
                for(int ti = 1 ; ti < proinfo->number_of_images ; ti++)
                {
                    if(proinfo->check_selected_image[ti])
                    {
                        rsetkernel.ti = ti;
                        int Count_N[3] = {0};
                        double count_GNCC = 0;
                        double nccresult = 0.0;

                        KernelPatchArg patch = proinfo->check_Matchtag ?
                            KernelPatchArg{
                                rsetkernel,
                                rlevelinfo.py_Sizes[rsetkernel.reference_id][*rlevelinfo.blunder_selected_level],
                                rlevelinfo.py_Sizes[rsetkernel.ti][*rlevelinfo.blunder_selected_level],
                                rlevelinfo.py_BImages[rsetkernel.reference_id],
                                rlevelinfo.py_BMagImages[rsetkernel.reference_id],
                                rlevelinfo.py_BImages[rsetkernel.ti],
                                rlevelinfo.py_BMagImages[rsetkernel.ti]}
                          : KernelPatchArg{
                                rsetkernel,
                                rlevelinfo.py_Sizes[rsetkernel.reference_id][*rlevelinfo.Pyramid_step],
                                 rlevelinfo.py_Sizes[rsetkernel.ti][*rlevelinfo.Pyramid_step],
                                 rlevelinfo.py_Images[rsetkernel.reference_id],
                                 rlevelinfo.py_MagImages[rsetkernel.reference_id],
                                 rlevelinfo.py_Images[rsetkernel.ti],
                                 rlevelinfo.py_MagImages[rsetkernel.ti]};
                        
                        //if(GridPT3[pt_index].Height != -1000)
                        //    count_total+=1;

                        for(int row = -Half_template_size; row <= Half_template_size ; row++)
                        {
                            for(int col = -Half_template_size; col <= Half_template_size ; col++)
                            {
                                const int radius2  =  row*row + col*col;
                                if(radius2 <= (Half_template_size+1)*(Half_template_size+1))
                                {
                                    const double t_X     = rlevelinfo.GridPts[pt_index].m_X + col*im_resolution;
                                    const double t_Y     = rlevelinfo.GridPts[pt_index].m_Y + row*im_resolution;
                                    
                                    long int t_col   = (long int)((t_X - subBoundary[0])/im_resolution);
                                    long int t_row   = (long int)((t_Y - subBoundary[1])/im_resolution);
                                    long int pt_index_temp = t_row*sub_imagesize_w + t_col;
                                    
                                    const long int tt_col  = (long int)((t_X - subBoundary[0])/(*rlevelinfo.grid_resolution));
                                    const long int tt_row  = (long int)((t_Y - subBoundary[1])/(*rlevelinfo.grid_resolution));
                                    const long int pt_index_dem  = tt_row*(long int)rlevelinfo.Size_Grid2D->width + tt_col;
                                    
                                    if(pt_index_temp >= 0 && pt_index_temp < sub_imagesize_w *sub_imagesize_h && t_col >= 0 && t_col < sub_imagesize_w && t_row >=0 && t_row < sub_imagesize_h && pt_index_dem >= 0 && pt_index_dem < *rlevelinfo.Grid_length && tt_col >= 0 && tt_col < rlevelinfo.Size_Grid2D->width && tt_row >=0 && tt_row < rlevelinfo.Size_Grid2D->height && all_im_cd[reference_id] != NULL && all_im_cd[ti] != NULL)
                                    {
                                        if(GridPT3[pt_index_dem].Height != -1000)
                                        {
                                            D2DPOINT pos_left(all_im_cd[reference_id][pt_index_temp]);
                                            D2DPOINT pos_right(all_im_cd[ti][pt_index_temp]);
                                            
                                            SetVecKernelValue(patch, row, col, pos_left,pos_right, radius2, Count_N);
                                        }
                                    }
                                }  // if(radius2 <= ...
                            }  // end col loop
                        }
                        
                        // Compute correlations
                        ComputeMultiNCC(rsetkernel, 0, Count_N, count_GNCC,  nccresult);
                        if(Count_N[0] == 0)
                            nccresult = -10;
                        
                        if(nccresult < 0.3)
                        {
                            if(!proinfo->check_Matchtag)
                            {
                                GridPT3[pt_index].minHeight     -= 100;
                                if(GridPT3[pt_index].minHeight < minmaxHeight[0])
                                    GridPT3[pt_index].minHeight     = floor(minmaxHeight[0]);
                                GridPT3[pt_index].maxHeight     += 100;
                                if(GridPT3[pt_index].maxHeight > minmaxHeight[1])
                                    GridPT3[pt_index].maxHeight     = ceil(minmaxHeight[1]);

                                if(nccresult < 0.1)
                                {
                                    GridPT3[pt_index].minHeight  = Nodata;
                                    GridPT3[pt_index].maxHeight  = Nodata;
                                }
                            }

                            //if(nccresult > -1)
                            //    count_low += 1;
                        }

                        GridPT3[pt_index].ortho_ncc[ti] = DoubleToSignedChar_grid(nccresult);


                        if(proinfo->check_Matchtag)
                        {
                            if(nccresult >= 0.2)
                                GridPT3[pt_index].roh   = DoubleToSignedChar_grid(nccresult);
                            else
                                GridPT3[pt_index].roh   = DoubleToSignedChar_grid(-0.2);
                        }
                    }
                } // end ti loop
            }
        } // end omp for

    } // end omp parallel
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(proinfo->check_selected_image[ti])
        {
            if(all_im_cd[ti])
                free(all_im_cd[ti]);
        }
    }
    if (all_im_cd)
        free(all_im_cd);

    //printf("%d %d\n",count_low,count_total);
}

bool VerticalLineLocus_blunder(const ProInfo *proinfo,LevelInfo &rlevelinfo, float* nccresult, UGRID *GridPT3, uint8 iteration, bool bblunder)
{
    double im_resolution = proinfo->resolution;
    uint8 Template_size = *rlevelinfo.Template_size;
    int Pyramid_step = *rlevelinfo.Pyramid_step;
    
    double template_area = 5.0;
    int t_Template_size = (int)((template_area/(im_resolution*pwrtwo(*rlevelinfo.blunder_selected_level)))/2.0)*2+1;
    if(Template_size < t_Template_size)
        Template_size = t_Template_size;
    
    int Half_template_size = (int)(Template_size/2.0);
    
    /*if(bblunder)
    {
        if(*rlevelinfo.Pyramid_step == 4)
        {
            Half_template_size = Template_size - iteration*2;
            if(Half_template_size < (int)(Template_size/2.0))
                Half_template_size = (int)(Template_size/2.0);
        }
        else if(*rlevelinfo.Pyramid_step == 3)
        {
            Half_template_size = Template_size - 2  - (iteration)*2;
            if(Half_template_size < (int)(Template_size/2.0))
                Half_template_size = (int)(Template_size/2.0);
        }
        else
        {
            Half_template_size = (int)(Template_size/2.0);
        }
    }
    */
    
    const double subBoundary[4] = {rlevelinfo.Boundary[0], rlevelinfo.Boundary[1], rlevelinfo.Boundary[2], rlevelinfo.Boundary[3]};
    const long numofpts = *rlevelinfo.Grid_length;
    im_resolution = proinfo->resolution*pwrtwo(*rlevelinfo.blunder_selected_level);
    
    D2DPOINT **all_im_cd = (D2DPOINT**)malloc(sizeof(D2DPOINT*)*proinfo->number_of_images);
    D2DPOINT **all_im_cd_next = NULL;
    
    long sub_imagesize_w, sub_imagesize_h;
    long sub_imagesize_w_next, sub_imagesize_h_next;
    
    SetOrthoImageCoord(proinfo, rlevelinfo, GridPT3, false, BD, im_resolution, 0, sub_imagesize_w, sub_imagesize_h, sub_imagesize_w_next, sub_imagesize_h_next, all_im_cd, all_im_cd_next);
    
    //printf("VerticalLineLocus_blunder omp\n");
    
#pragma omp parallel
    {
        SetKernel rsetkernel(0,1,Half_template_size);
        // Make patch vectors thread private rather than private to each loop iteration
#pragma omp for schedule(guided)
        for(long iter_count = 0 ; iter_count < numofpts ; iter_count++)
        {
            long pts_row = (int)(floor(iter_count/rlevelinfo.Size_Grid2D->width));
            long pts_col = iter_count % rlevelinfo.Size_Grid2D->width;
            long pt_index = iter_count;//pts_row*(long)rlevelinfo.Size_Grid2D->width + pts_col;
            /*
            if(GridPT3[pt_index].ncc_seleceted_pair < 0 || GridPT3[pt_index].ncc_seleceted_pair >= rlevelinfo.pairinfo->NumberOfPairs)
            {
                printf("VerticalLineLocus_blunder check ncc_seleceted_pair %d\n",GridPT3[pt_index].ncc_seleceted_pair);
                exit(1);
            }*/
            //printf("GridPT3[pt_index].ncc_seleceted_pair %d\n",GridPT3[pt_index].ncc_seleceted_pair);
            if(pt_index < *rlevelinfo.Grid_length && pts_row < rlevelinfo.Size_Grid2D->height && pts_col < rlevelinfo.Size_Grid2D->width && pts_row >= 0 && pts_col >= 0 /*&& GridPT3[pt_index].ncc_seleceted_pair > -1*/)
            {
                nccresult[pt_index] = -1.0;
                double max_ncc = 100;
                double sum_ncc = 0;
                int total_ncc = 0;
                
                bool check_AWNCC = false;
                if(GridPT3[pt_index].ncc_seleceted_pair == rlevelinfo.pairinfo->NumberOfPairs)
                    check_AWNCC = true;
                
                for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
                {
                    const int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                    const int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                    
                    if((reference_id == GridPT3[pt_index].selected_pair || ti == GridPT3[pt_index].selected_pair) /*&& GridPT3[pt_index].ncc_seleceted_pair == pair_number*/)
                    {
                        rsetkernel.reference_id = reference_id;
                        rsetkernel.ti = ti;
                        
                        KernelPatchArg patch{
                            rsetkernel,
                            rlevelinfo.py_Sizes[rsetkernel.reference_id][*rlevelinfo.blunder_selected_level],
                            rlevelinfo.py_Sizes[rsetkernel.ti][*rlevelinfo.blunder_selected_level],
                            rlevelinfo.py_BImages[rsetkernel.reference_id],
                            rlevelinfo.py_BMagImages[rsetkernel.reference_id],
                            rlevelinfo.py_BImages[rsetkernel.ti],
                            rlevelinfo.py_BMagImages[rsetkernel.ti]};
                        
                        
                        int Count_N[3] = {0};
                        double count_GNCC = 0;
                        double t_nccresult = 0;
                        
                        for(int row = -Half_template_size; row <= Half_template_size ; row++)
                        {
                            for(int col = -Half_template_size; col <= Half_template_size ; col++)
                            {
                                const int radius2  = row*row + col*col;
                                if(radius2 <= (Half_template_size-1)*(Half_template_size-1))
                                {
                                    const double t_X     = rlevelinfo.GridPts[pt_index].m_X + col*im_resolution;
                                    const double t_Y     = rlevelinfo.GridPts[pt_index].m_Y + row*im_resolution;
                                    
                                    long int t_col   = (long int)((t_X - subBoundary[0])/im_resolution);
                                    long int t_row   = (long int)((t_Y - subBoundary[1])/im_resolution);
                                    long int pt_index_temp = t_row*sub_imagesize_w + t_col;
                                    
                                    const long int tt_col  = (long int)((t_X - subBoundary[0])/(*rlevelinfo.grid_resolution));
                                    const long int tt_row  = (long int)((t_Y - subBoundary[1])/(*rlevelinfo.grid_resolution));
                                    const long int pt_index_dem  = tt_row*(long int)rlevelinfo.Size_Grid2D->width + tt_col;
                                    
                                    if(pt_index_temp >= 0 && pt_index_temp < sub_imagesize_w *sub_imagesize_h && t_col >= 0 && t_col < sub_imagesize_w && t_row >=0 && t_row < sub_imagesize_h && pt_index_dem >= 0 && pt_index_dem < numofpts && tt_col >= 0 && tt_col < rlevelinfo.Size_Grid2D->width && tt_row >=0 && tt_row < rlevelinfo.Size_Grid2D->height && all_im_cd[reference_id] != NULL && all_im_cd[ti] != NULL)
                                    {
                                        if(GridPT3[pt_index_dem].Height != -1000)
                                        {
                                            D2DPOINT pos_left(all_im_cd[reference_id][pt_index_temp]);
                                            D2DPOINT pos_right(all_im_cd[ti][pt_index_temp]);
                                            D2DPOINT pos_right_before(pos_right);
                                            
                                            pos_right.m_X = pos_right.m_X + rlevelinfo.ImageAdjust[pair_number][1]/pwrtwo(*rlevelinfo.blunder_selected_level);
                                            pos_right.m_Y = pos_right.m_Y + rlevelinfo.ImageAdjust[pair_number][0]/pwrtwo(*rlevelinfo.blunder_selected_level);
                                            
                                            //printf("rpc bias %f\t%f\t before %f\t%f\t after %f\t%f\n",rlevelinfo.ImageAdjust[pair_number][0],rlevelinfo.ImageAdjust[pair_number][1],pos_right_before.m_X,pos_right_before.m_Y,pos_right.m_X,pos_right.m_Y);
                                            //exit(1);
                                            SetVecKernelValue(patch, row, col, pos_left,pos_right, radius2, Count_N);
                                        }
                                    }
                                } // if(radius <= Half_template_size-1)
                            } // end col loop
                        } // end row loop
                        
                        // Compute collelations
                        ComputeMultiNCC(rsetkernel, 0, Count_N, count_GNCC,  t_nccresult);
                        
                        if(t_nccresult > -1.0)
                        {
                            if(max_ncc > t_nccresult)
                                max_ncc = t_nccresult;
                            
                            sum_ncc += t_nccresult;
                            total_ncc++;
                            
                            GridPT3[pt_index].ortho_ncc[pair_number] = DoubleToSignedChar_grid(t_nccresult);
                        }
                    }
                } // end ti loop
                
                //if(!check_AWNCC && GridPT3[pt_index].ncc_seleceted_pair > -1)
                {
                    if(Pyramid_step >= 0)
                    {
                        if(check_AWNCC)
                        {
                            GridPT3[pt_index].Mean_ortho_ncc = DoubleToSignedChar_grid(max_ncc); //WNCC weight
                            nccresult[pt_index] = SignedCharToDouble_grid(GridPT3[pt_index].Mean_ortho_ncc); //blunder detection
                        }
                        else
                        {
                            //GridPT3[pt_index].Max_ortho_ncc = DoubleToSignedChar_grid(0); //WNCC weight
                            GridPT3[pt_index].Mean_ortho_ncc = GridPT3[pt_index].ortho_ncc[GridPT3[pt_index].ncc_seleceted_pair]; //WNCC weight
                            //nccresult[pt_index] = 0; //blunder detection
                            nccresult[pt_index] = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[GridPT3[pt_index].ncc_seleceted_pair]);
                            
                            //GridPT3[pt_index].Mean_ortho_ncc = DoubleToSignedChar_grid(max_ncc); //WNCC weight
                            //nccresult[pt_index] = /*(sum_ncc/(double)total_ncc);//*/max_ncc; //blunder detection
                        }
                    }
                    else
                    {
                        GridPT3[pt_index].Mean_ortho_ncc = GridPT3[pt_index].ortho_ncc[GridPT3[pt_index].ncc_seleceted_pair]; //WNCC weight
                        nccresult[pt_index] = SignedCharToDouble_grid(GridPT3[pt_index].ortho_ncc[GridPT3[pt_index].ncc_seleceted_pair]);//*/weight_ncc; //blunder detection
                    }
                }
                
                //GridPT3[pt_index].Mean_ortho_ncc = DoubleToSignedChar_grid(max_ncc);
                
                
                //nccresult[pt_index] = max_ncc;
            }
        } // end omp for
    } // end omp parallel
 
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(proinfo->check_selected_image[ti])
        {
            if(all_im_cd[ti])
                free(all_im_cd[ti]);
        }
    }
    if (all_im_cd)
        free(all_im_cd);
    
    return true;
}

double SetMultiWeight(int pairnumbers, vector<double> &save_roh_positive)
{
    double additional_weight = 0.2;//0.05*(pairnumbers - count_nagative);
    /*if(additional_weight < 0)
        additional_weight = 0;
    if(additional_weight > 0.3)
        additional_weight = 0.3;
     */
    //double ncc_plus_rate = 1.0 + additional_weight*(pairnumbers - count_nagative);
    
    double ncc_plus_rate = 1.0 + additional_weight*(double)save_roh_positive.size()/(double)pairnumbers;
    if(ncc_plus_rate > 1.0 + additional_weight)
        ncc_plus_rate = 1.0 + additional_weight;
    //if(count_nagative == 0)
    //    ncc_plus_rate = 1.0;
    
    
    return ncc_plus_rate;
}

int Ortho_blunder(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3)
{
    const int max_count = 200;
    int while_count = 0;
    
    const long num_triangles = numOfTri;
    const int pyramid_step    = *rlevelinfo.Pyramid_step;
    const double gridspace = *rlevelinfo.grid_resolution;
    const double *boundary    = rlevelinfo.Boundary;
 
    uint8* tris_check = (uint8*)calloc(num_triangles,sizeof(uint8));
    bool check_stop_TIN = false;
    
    printf("tric_check done\n");
    while(!check_stop_TIN && while_count < max_count)
    {
        bool check_ortho_cal = false;
      
        while_count++;
      
        check_stop_TIN = true;
      
        double *updated_height = (double*)malloc(sizeof(double)*num_triangles);
        int *selected_index = (int*)malloc(sizeof(int)*num_triangles);
        bool *updated_check = (bool*)calloc(sizeof(bool),num_triangles);
        double *selected_count = (double*)calloc(sizeof(double),num_triangles);
        double *FNCC = (double*)calloc(sizeof(double),num_triangles);
        int *selected_target_index = (int*)malloc(sizeof(int)*num_triangles);
        double* com_count = (double*)calloc(sizeof(double),numOfPts);
        double* com_FNCC = (double*)calloc(sizeof(double),numOfPts);
        
#pragma omp parallel for schedule(dynamic, 1)
        for(int tcnt=0;tcnt<(int)(num_triangles);tcnt++)
        {
            if(tris_check[tcnt] == 0)
            {
                const UI3DPOINT t_tri   = tris[tcnt];
                const long pdex0 = t_tri.m_X;
                const long pdex1 = t_tri.m_Y;
                const long pdex2 = t_tri.m_Z;
                
                if(pdex0 < numOfPts && pdex1 < numOfPts && pdex2 < numOfPts)
                {
                    const D3DPOINT pt0(pts[pdex0]);
                    const D3DPOINT pt1(pts[pdex1]);
                    const D3DPOINT pt2(pts[pdex2]);
                    
                    const long node1_index = (long)((pt0.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt0.m_X - boundary[0])/gridspace + 0.5);
                    const long node2_index = (long)((pt1.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt1.m_X - boundary[0])/gridspace + 0.5);
                    const long node3_index = (long)((pt2.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt2.m_X - boundary[0])/gridspace + 0.5);
                    
                    const uint8 node1_F     =  GridPT3[node1_index].anchor_flag;
                    const uint8 node2_F     =  GridPT3[node2_index].anchor_flag;
                    const uint8 node3_F     =  GridPT3[node3_index].anchor_flag;
                    
                    int node_f_count = 0;
                    
                    if(node1_F == 1 || node1_F == 3)
                        node_f_count++;
                    if(node2_F == 1 || node2_F == 3)
                        node_f_count++;
                    if(node3_F == 1 || node3_F == 3)
                        node_f_count++;
                    
                    if(node_f_count == 3)
                        tris_check[tcnt] = 1;
                    else if(node_f_count == 2)
                    {
                        //ortho matching process
                        D3DPOINT ref1_pt, ref2_pt, target_pt;
                        long target_index,ref1_index, ref2_index;
                        long target_pt_index;
                        
                        if(node1_F == 2)
                        {
                            ref1_pt     = pt1;
                            ref2_pt     = pt2;
                            target_pt   = pt0;
                            target_index= node1_index;
                            ref1_index  = node2_index;
                            ref2_index  = node3_index;
                            target_pt_index = pdex0;
                        }
                        else if(node2_F == 2)
                        {
                            ref1_pt     = pt0;
                            ref2_pt     = pt2;
                            target_pt   = pt1;
                            target_index= node2_index;
                            ref1_index  = node1_index;
                            ref2_index  = node3_index;
                            target_pt_index = pdex1;
                        }
                        else
                        {
                            ref1_pt     = pt0;
                            ref2_pt     = pt1;
                            target_pt   = pt2;
                            target_index= node3_index;
                            ref1_index  = node1_index;
                            ref2_index  = node2_index;
                            target_pt_index = pdex2;
                        }
                        
                        double F_SNCC, F_height;
                        double t_selected_count = VerticalLineLocus_Ortho(proinfo,rlevelinfo, MPP, &F_height, ref1_pt,ref2_pt,target_pt, GridPT3,target_index,&F_SNCC);
             
                        if(F_height != Nodata )
                        {
                            updated_height[tcnt] = F_height;
                            selected_index[tcnt] = target_pt_index;
                            updated_check[tcnt] = true;
                            selected_count[tcnt] = t_selected_count;
                            FNCC[tcnt] = F_SNCC;
                            selected_target_index[tcnt] = target_index;
                        }
                    }
                }
            }
        }
        
        for(int tcnt=0;tcnt<(int)(num_triangles);tcnt++)
        {
            if(tris_check[tcnt] == 0)
            {
                const UI3DPOINT t_tri   = tris[tcnt];
                const long pdex0 = t_tri.m_X;
                const long pdex1 = t_tri.m_Y;
                const long pdex2 = t_tri.m_Z;
                
                if(pdex0 < numOfPts && pdex1 < numOfPts && pdex2 < numOfPts)
                {
                    const D3DPOINT pt0(pts[pdex0]);
                    const D3DPOINT pt1(pts[pdex1]);
                    const D3DPOINT pt2(pts[pdex2]);
                    
                    const long node1_index = (long)((pt0.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt0.m_X - boundary[0])/gridspace + 0.5);
                    const long node2_index = (long)((pt1.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt1.m_X - boundary[0])/gridspace + 0.5);
                    const long node3_index = (long)((pt2.m_Y - boundary[1])/gridspace + 0.5)*rlevelinfo.Size_Grid2D->width + (long)((pt2.m_X - boundary[0])/gridspace + 0.5);
                    
                    if(updated_check[tcnt])
                    {
                        const int target_pt_index = selected_index[tcnt];
                        if(com_count[target_pt_index] < selected_count[tcnt])
                        {
                            com_count[target_pt_index] = selected_count[tcnt];
                            com_FNCC[target_pt_index] = FNCC[tcnt];
                            
                            const int target_index = selected_target_index[tcnt];
                            GridPT3[target_index].anchor_flag = 3;
                            pts[target_pt_index].m_Z = updated_height[tcnt];
                            pts[target_pt_index].flag = 0;
                            
                            check_stop_TIN = false;
                            check_ortho_cal = true;
                            tris_check[tcnt] = 1;
                        }
                    }
                    
                    const uint8 node1_F     =  GridPT3[node1_index].anchor_flag;
                    const uint8 node2_F     =  GridPT3[node2_index].anchor_flag;
                    const uint8 node3_F     =  GridPT3[node3_index].anchor_flag;
                    
                    int node_f_count = 0;
                    if(node1_F == 1 || node1_F == 3)
                    {
                        node_f_count++;
                        pts[pdex0].flag = 0;
                    }
                    if(node2_F == 1 || node2_F == 3)
                    {
                        node_f_count++;
                        pts[pdex1].flag = 0;
                    }
                    if(node3_F == 1 || node3_F == 3)
                    {
                        node_f_count++;
                        pts[pdex2].flag = 0;
                    }
                    
                    if(node_f_count == 3)
                        tris_check[tcnt] = 1;
                }
            }
        }
        
        free(updated_height);
        free(selected_index);
        free(updated_check);
        free(selected_count);
        free(FNCC);
        free(selected_target_index);
        free(com_count);
        free(com_FNCC);
        
        if(check_ortho_cal == false)
            check_stop_TIN = true;
    }
    
    printf("ortho bluncer iteration %d\n",while_count);
    
    free(tris_check);
    
    return numOfPts;
}

int VerticalLineLocus_Ortho(ProInfo *proinfo, LevelInfo &rlevelinfo, double MPP, double *F_Height, D3DPOINT ref1_pt, D3DPOINT ref2_pt, D3DPOINT target_pt, UGRID *GridPT3, int target_index, double *F_sncc)
{
    double F_NCC = -1.0;
    bool check_ncc = false;
    const double gridspace = *rlevelinfo.grid_resolution;
    const int Pyramid_step = *rlevelinfo.Pyramid_step;
    
    *F_Height = Nodata;
    
    const uint32 TIN_Grid_Size_X = rlevelinfo.Size_Grid2D->width;
    const uint32 TIN_Grid_Size_Y = rlevelinfo.Size_Grid2D->height;
    
    if(MPP > 3)
        MPP = 3;
    if(MPP < 1)
        MPP = 1;
    
    double height_step;
    //if(*rlevelinfo.iteration <= 2)
        height_step = proinfo->resolution*pwrtwo(Pyramid_step)*MPP;
    //else
    //    height_step = pwrtwo(Pyramid_step)*MPP;
    
    int start_H     = GridPT3[target_index].minHeight;
    int end_H       = GridPT3[target_index].maxHeight;
    
    const int NumOfHeights = (int)((end_H -  start_H)/height_step) + 1;
    
    const int th_count = 10;
    int min_th_count = 2;
    double min_th_sncc = 0.2;
    if(*rlevelinfo.Pyramid_step < 4)
    {
        min_th_count = 0;
        min_th_sncc = 0.5;
    }
    
    int reference_id = rlevelinfo.reference_id;
    double selected_count = 0;
    
    for(long count_height = 0 ; count_height < NumOfHeights ; count_height++)
    {
        float iter_height      = start_H + count_height*height_step;
        
        D3DPOINT TriP1(ref1_pt.m_X, ref1_pt.m_Y, ref1_pt.m_Z, 0);
        D3DPOINT TriP2(ref2_pt.m_X, ref2_pt.m_Y, ref2_pt.m_Z, 0);
        D3DPOINT TriP3(target_pt.m_X, target_pt.m_Y, (double)iter_height, 0);
        
        // calculation on BoundingBox(MinMax XY) of triangle
        const double TriMinXY[2] = {min(min(TriP1.m_X,TriP2.m_X),TriP3.m_X), min(min(TriP1.m_Y,TriP2.m_Y),TriP3.m_Y)};
        const double TriMaxXY[2] = {max(max(TriP1.m_X,TriP2.m_X),TriP3.m_X), max(max(TriP1.m_Y,TriP2.m_Y),TriP3.m_Y)};
        
        int PixelMinXY[2] = {(int)((TriMinXY[0] - rlevelinfo.Boundary[0])/gridspace + 0.5), (int)((TriMinXY[1] - rlevelinfo.Boundary[1])/gridspace + 0.5)};
        int PixelMaxXY[2] = {(int)((TriMaxXY[0] - rlevelinfo.Boundary[0])/gridspace + 0.5), (int)((TriMaxXY[1] - rlevelinfo.Boundary[1])/gridspace + 0.5)};
        
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
        
        if(PixelMaxXY[1]-PixelMinXY[1] > 0 && PixelMaxXY[0]-PixelMinXY[0] > 0)
        {
            int patch_size = (PixelMaxXY[1]-PixelMinXY[1]+1) * (PixelMaxXY[0]-PixelMinXY[0]+1);
            double *left_patch_vec = (double *)malloc(sizeof(double)*patch_size);
            double *right_patch_vec = (double *)malloc(sizeof(double)*patch_size);
            double *left_mag_patch_vec = (double *)malloc(sizeof(double)*patch_size);
            double *right_mag_patch_vec = (double *)malloc(sizeof(double)*patch_size);
            
            double sum_NCC = -1.0;
            bool count_NCC = false;
            int sum_count_N = 0;
            int total_ncc_count = 0;
            
            for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
            {
                reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
                const int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
                
                if((reference_id == GridPT3[target_index].selected_pair || ti == GridPT3[target_index].selected_pair) && GridPT3[target_index].ncc_seleceted_pair == pair_number)
                {
                    int Count_N = 0;

                    CSize LImagesize, RImagesize;
                    
                    LImagesize.width  = rlevelinfo.py_Sizes[reference_id][Pyramid_step].width;
                    LImagesize.height = rlevelinfo.py_Sizes[reference_id][Pyramid_step].height;
                    RImagesize.width  = rlevelinfo.py_Sizes[ti][Pyramid_step].width;
                    RImagesize.height = rlevelinfo.py_Sizes[ti][Pyramid_step].height;
                    
                    for (long Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
                    {
                        for (long Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
                        {
                            bool rtn = false;
                            double Z = 0.0;
                            D3DPOINT CurGPXY( Col*gridspace + rlevelinfo.Boundary[0] , Row*gridspace + rlevelinfo.Boundary[1], 0, 0);
                            const D3DPOINT v12(TriP2 - TriP1);
                            const D3DPOINT v1P(CurGPXY - TriP1);
                            const D3DPOINT v23(TriP3 - TriP2);
                            const D3DPOINT v2P(CurGPXY - TriP2);
                            const D3DPOINT v31(TriP1 - TriP3);
                            const D3DPOINT v3P(CurGPXY - TriP3);
                            
                            int Sum = 3;
                            if (v12.m_X*v1P.m_Y-v12.m_Y*v1P.m_X <= 0)
                                Sum--;
                            if (v23.m_X*v2P.m_Y-v23.m_Y*v2P.m_X <= 0)
                                Sum--;
                            if (v31.m_X*v3P.m_Y-v31.m_Y*v3P.m_X <= 0)
                                Sum--;
                            
                            if (Sum==0 || Sum==3)
                            {
                                const D3DPOINT v12(TriP2 - TriP1);
                                const D3DPOINT v13(TriP3 - TriP1);
                                D3DPOINT Normal(v12.m_Y*v13.m_Z - v12.m_Z*v13.m_Y, v12.m_Z*v13.m_X - v12.m_X*v13.m_Z, v12.m_X*v13.m_Y - v12.m_Y*v13.m_X, 0);
                                const double Len = SQRT(Normal);
                                
                                if(Len > 0)
                                {
                                    Normal.m_X/=Len; Normal.m_Y/=Len; Normal.m_Z/=Len;
                                    const double A = Normal.m_X;
                                    const double B = Normal.m_Y;
                                    const double C = Normal.m_Z;
                                    const double D = -(A*TriP1.m_X+B*TriP1.m_Y+C*TriP1.m_Z);
                                    
                                    if(C != 0)
                                    {
                                        Z = -1.0 * ((A * CurGPXY.m_X) + (B * CurGPXY.m_Y) + D) / C;
                                        rtn = true;
                                    }
                                    else
                                        rtn = false;
                                }
                            }
                            
                            if (rtn)
                            {
                                D2DPOINT Ref_Imagecoord, temp_GP_p;
                                D3DPOINT temp_GP;
                                
                                const double temp_LIA[2] = {0.0, 0.0};
                                
                                temp_GP.m_Z = Z;
                                if(proinfo->sensor_type == SB)
                                {
                                    temp_GP_p.m_X = CurGPXY.m_X;
                                    temp_GP_p.m_Y = CurGPXY.m_Y;
                                    
                                    temp_GP     = ps2wgs_single(*rlevelinfo.param,temp_GP_p);
                                    
                                    Ref_Imagecoord     = GetObjectToImageRPC_single(rlevelinfo.RPCs[reference_id],*rlevelinfo.NumOfIAparam,temp_LIA,temp_GP);
                                }
                                else
                                {
                                    temp_GP.m_X = CurGPXY.m_X;
                                    temp_GP.m_Y = CurGPXY.m_Y;
                                    
                                    D2DPOINT photo = GetPhotoCoordinate_single(temp_GP,proinfo->frameinfo.Photoinfo[reference_id],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[reference_id].m_Rm);
                                    
                                    Ref_Imagecoord = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
                                }
                                const D2DPOINT Ref_Imagecoord_py  = OriginalToPyramid_single(Ref_Imagecoord,rlevelinfo.py_Startpos[reference_id],Pyramid_step);
                                
                                D2DPOINT Tar_Imagecoord;
                                double pos_row_left;
                                double pos_col_left;
                                double pos_row_right;
                                double pos_col_right;
                                
                                if(proinfo->sensor_type == SB)
                                {
                                    temp_GP_p.m_X = CurGPXY.m_X;
                                    temp_GP_p.m_Y = CurGPXY.m_Y;
                                    
                                    temp_GP     = ps2wgs_single(*rlevelinfo.param,temp_GP_p);
                                    
                                    Tar_Imagecoord     = GetObjectToImageRPC_single(rlevelinfo.RPCs[ti],*rlevelinfo.NumOfIAparam,rlevelinfo.ImageAdjust[pair_number],temp_GP);
                                }
                                else
                                {
                                    temp_GP.m_X = CurGPXY.m_X;
                                    temp_GP.m_Y = CurGPXY.m_Y;
                                    
                                    D2DPOINT photo = GetPhotoCoordinate_single(temp_GP,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera,proinfo->frameinfo.Photoinfo[ti].m_Rm);
                                    
                                    Tar_Imagecoord = PhotoToImage_single(photo,proinfo->frameinfo.m_Camera.m_CCDSize,proinfo->frameinfo.m_Camera.m_ImageSize);
                                }
                                const D2DPOINT Tar_Imagecoord_py = OriginalToPyramid_single(Tar_Imagecoord,rlevelinfo.py_Startpos[ti],Pyramid_step);
                                
                                pos_row_left      = Ref_Imagecoord_py.m_Y;
                                pos_col_left      = Ref_Imagecoord_py.m_X;
                                
                                pos_row_right     = Tar_Imagecoord_py.m_Y;
                                pos_col_right     = Tar_Imagecoord_py.m_X;
                                
                                if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right  >= 0 && pos_col_right+1 < RImagesize.width &&
                                   pos_row_left >= 0 && pos_row_left+1   < LImagesize.height && pos_col_left   >= 0 && pos_col_left+1  < LImagesize.width)
                                {
                                    //interpolate left_patch
                                    double dx = pos_col_left - (int) (pos_col_left);
                                    double dy = pos_row_left - (int) (pos_row_left);
                                    long position = (long int) (pos_col_left) + (long int) (pos_row_left) * LImagesize.width;
                                    double left_patch = InterpolatePatch(rlevelinfo.py_Images[reference_id], position, LImagesize, dx, dy);
                                    double left_mag_patch = InterpolatePatch(rlevelinfo.py_MagImages[reference_id], position, LImagesize, dx, dy);
                                    
                                    //interpolate right_patch
                                    dx = pos_col_right - (int) (pos_col_right);
                                    dy = pos_row_right - (int) (pos_row_right);
                                    position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
                                    double right_patch = InterpolatePatch(rlevelinfo.py_Images[ti], position, RImagesize, dx, dy);
                                    double right_mag_patch = InterpolatePatch(rlevelinfo.py_MagImages[ti], position, RImagesize, dx, dy);
                                    
                                    //end
                                    left_patch_vec[Count_N] = left_patch;
                                    left_mag_patch_vec[Count_N] = left_mag_patch;
                                    right_patch_vec[Count_N] = right_patch;
                                    right_mag_patch_vec[Count_N] = right_mag_patch;
                                    Count_N++;
                                }
                            }  // if (rtn)
                        }  //  // end Col=PixelMinXY[0] loop
                    }  // end Row=PixelMinXY[1] loop
            
                    // Compute correlations
                    if(Count_N >= 1)
                    {
                        double temp_roh = 0;
                        double count_roh = 0;
                        double ncc_1 = Correlate(left_patch_vec, right_patch_vec, Count_N);
                        if (ncc_1 != -99)
                        {
                            count_roh++;
                            temp_roh += ncc_1;
                        }
                            
                        double ncc_2 = Correlate(left_mag_patch_vec, right_mag_patch_vec, Count_N);
                        if (ncc_2 != -99)
                        {
                            count_roh++;
                            temp_roh += ncc_2;
                        }
                        
                        double ncc;
                        if (count_roh > 0)
                        {
                            ncc = temp_roh/count_roh;
                            
                            if(sum_NCC < ncc)//maximum
                            {
                                sum_count_N = Count_N;
                                sum_NCC = ncc;
                                count_NCC = true;
                                total_ncc_count = 1;
                            }
                        }
                        else
                            ncc = -1;

                        //sum_count_N += Count_N;
                        //sum_NCC += ncc;
                        //count_NCC++;
                    }
                }
            }  // end ti loop
            
            free(left_patch_vec);
            free(right_patch_vec);
            free(left_mag_patch_vec);
            free(right_mag_patch_vec);
        
            if(count_NCC)
            {
                double sncc = sum_NCC;
                double avg_ncc_count = (double)sum_count_N;
                if(avg_ncc_count >= th_count)
                {
                    if(sncc > min_th_sncc)
                    {
                        if(F_NCC < sncc)
                        {
                            F_NCC = sncc;
                            *F_Height = (double)iter_height;
                            selected_count =  avg_ncc_count;
                            *F_sncc = F_NCC;
                        }
                    }
                }
                else if(avg_ncc_count > min_th_count)
                {
                    if(sncc > min_th_sncc  + 0.2)
                    {
                        if(F_NCC < sncc)
                        {
                            F_NCC = sncc;
                            *F_Height = (ref1_pt.m_Z + ref2_pt.m_Z)/2.0;
                            selected_count =  avg_ncc_count;
                            *F_sncc = F_NCC;
                        }
                    }
                }
            }
            
            /*
            if(count_NCC > 0)
            {
                if(sum_count_N/count_NCC >= th_count)
                {
                    double sncc = sum_NCC/count_NCC;
                    if(sncc > min_th_sncc)
                    {
                        check_ncc = true;
                        if(F_NCC < sncc)
                        {
                            F_NCC = sncc;
                            *F_Height = (double)iter_height;
                            selected_count =  ((double)(sum_count_N)/(double)(count_NCC));
                            *F_sncc = F_NCC;
                        }
                    }
                }
                else if(sum_count_N/count_NCC > min_th_count)
                {
                    double sncc = sum_NCC/count_NCC;
                    if(sncc > min_th_sncc + 0.2)
                    {
                        check_ncc = true;
                        if(F_NCC < sncc)
                        {
                            F_NCC = sncc;
                            *F_Height = (ref1_pt.m_Z + ref2_pt.m_Z)/2.0;
                            selected_count =  ((double)(sum_count_N)/(double)(count_NCC));
                            *F_sncc = F_NCC;
                        }
                    }
                }
            }
             */
        }
    }
    
    return selected_count;
}

long SelectMPs(const ProInfo *proinfo,LevelInfo &rlevelinfo, const NCCresult* roh_height, UGRID *GridPT3, const double Th_roh, const double Th_roh_min, const double Th_roh_start, const double Th_roh_next, const int iteration, const double MPP, const int final_level_iteration,const double MPP_stereo_angle, vector<D3DPOINT> *linkedlist)
{
    long int count_MPs = 0;

    double minimum_Th = 0.2;
    double minGrid_th = 0.2;
    
    int SGM_th_py = proinfo->SGM_py;
    const int Pyramid_step = *rlevelinfo.Pyramid_step;
    if(proinfo->IsRA)
    {
        minimum_Th = 0.2;
    }
    else
    {
        /*
        if(Pyramid_step > 0)
        {
            const int rohconvert = 100;
            long int* hist = (long int*)calloc(sizeof(long int),rohconvert);
            long int total_roh = 0;
            
            for(long int iter_index = 0 ; iter_index < *rlevelinfo.Grid_length ; iter_index++)
            {
                long row     = (long)(floor(iter_index/rlevelinfo.Size_Grid2D->width));
                long col     = iter_index % rlevelinfo.Size_Grid2D->width;
                long grid_index = iter_index;//row*(long)rlevelinfo.Size_Grid2D->width + col;
                
                if(row >= 0 && row < rlevelinfo.Size_Grid2D->height && col >= 0 && col < rlevelinfo.Size_Grid2D->width && roh_height[grid_index].NumOfHeight > 2)
                {
                    int roh_int = ceil(SignedCharToDouble_result(roh_height[grid_index].result0)*rohconvert);
                    //printf("roh_int %d\t%f\n", roh_int,SignedCharToDouble_result(roh_height[grid_index].result0));
                    if(roh_int > rohconvert - 1)
                        roh_int = rohconvert - 1;
                    if(roh_int >= 0)
                        hist[roh_int]++;
                    
                    total_roh++;
                }
            }
            
            if(total_roh > 0)
            {
                double min_roh_th;
                
                if(Pyramid_step == 4)
                    min_roh_th = 0.20 + (iteration-1)*0.02;
                else if(Pyramid_step == 3)
                    min_roh_th = 0.25 + (iteration-1)*0.03;
                else if(Pyramid_step == 2)
                    min_roh_th = 0.60;
                else if(Pyramid_step == 1)
                    min_roh_th = 0.80;
                
                if(Pyramid_step <= SGM_th_py )
                {
                    min_roh_th = 0.95 - Pyramid_step*0.1;
                    if(min_roh_th < 0.5)
                        min_roh_th = 0.5;
                }
                
                int roh_iter = rohconvert - 1;
                bool check_stop = false;
                minimum_Th = 0.2;
                int sum_roh_count = 0;
                double sum_roh_rate = 0;
                while(!check_stop && roh_iter > 0)
                {
                    sum_roh_count += hist[roh_iter];
                    sum_roh_rate = sum_roh_count/(double)total_roh;
                    //printf("roh_iter %d\t%d\t%d\t%f\n",roh_iter,hist[roh_iter],total_roh,sum_roh_rate);
                    if(sum_roh_rate > min_roh_th)
                    {
                        check_stop = true;
                        minimum_Th = (double)roh_iter/(double)rohconvert;
                    }
                    roh_iter--;
                }
                if(minimum_Th > 0.95)
                    minimum_Th = 0.95;
            }
            else
                minimum_Th = 0.2;
            
            free(hist);
        }
        else
            minimum_Th = 0.2;
        */
        
        
        if(Pyramid_step > 0)
        {
            const int rohconvert = 1000;
            long int* hist = (long int*)calloc(sizeof(long int),rohconvert);
            long int total_roh = 0;
            
            for(long int iter_index = 0 ; iter_index < *rlevelinfo.Grid_length ; iter_index++)
            {
                long row     = (long)(floor(iter_index/rlevelinfo.Size_Grid2D->width));
                long col     = iter_index % rlevelinfo.Size_Grid2D->width;
                long grid_index = iter_index;//row*(long)rlevelinfo.Size_Grid2D->width + col;
                
                if(row >= 0 && row < rlevelinfo.Size_Grid2D->height && col >= 0 && col < rlevelinfo.Size_Grid2D->width && roh_height[grid_index].NumOfHeight > 2)
                {
                    int roh_int = ceil(SignedCharToDouble_result(roh_height[grid_index].result0)*rohconvert);
                    //printf("roh_int %d\t%f\n", roh_int,SignedCharToDouble_result(roh_height[grid_index].result0));
                    if(roh_int > rohconvert - 1)
                        roh_int = rohconvert - 1;
                    if(roh_int >= 0)
                        hist[roh_int]++;
                    
                    total_roh++;
                }
            }
            
            if(total_roh > 0)
            {
                double min_roh_th;
                double min_Gridroh_th;
                
                if(Pyramid_step == 4)
                    min_roh_th = 0.20 + (iteration-1)*0.01; //0.3 //0.05
                else if(Pyramid_step == 3)
                    min_roh_th = 0.30 + (iteration-1)*0.01; //0.5 //0.15
                else if(Pyramid_step == 2)
                    min_roh_th = 0.70;//0.60;
                else if(Pyramid_step == 1)
                    min_roh_th = 0.90;//0.80;
                
                min_Gridroh_th = min_roh_th - 0.2;
                //if(min_Gridroh_th > 0.9)
                //    min_Gridroh_th = 0.9;
                
                int roh_iter = rohconvert - 1;
                bool check_stop = false;
                
                minimum_Th = 0.2;
                minGrid_th = 0.2;
                
                int sum_roh_count = 0;
                double sum_roh_rate = 0;
                while(!check_stop && roh_iter > 0)
                {
                    sum_roh_count += hist[roh_iter];
                    sum_roh_rate = sum_roh_count/(double)total_roh;
                    //printf("roh_iter %d\t%d\t%d\t%f\n",roh_iter,hist[roh_iter],total_roh,sum_roh_rate);
                    if(sum_roh_rate > min_roh_th)
                    {
                        check_stop = true;
                        minimum_Th = (double)roh_iter/(double)rohconvert;
                    }
                    roh_iter--;
                }
                if(minimum_Th > 0.80)
                    minimum_Th = 0.80;
                
                roh_iter = rohconvert - 1;
                check_stop = false;
                sum_roh_count = 0;
                sum_roh_rate = 0;
                while(!check_stop && roh_iter > 0)
                {
                    sum_roh_count += hist[roh_iter];
                    sum_roh_rate = sum_roh_count/(double)total_roh;
                    //printf("roh_iter %d\t%d\t%d\t%f\n",roh_iter,hist[roh_iter],total_roh,sum_roh_rate);
                    if(sum_roh_rate > min_Gridroh_th)
                    {
                        check_stop = true;
                        minGrid_th = (double)roh_iter/(double)rohconvert;
                    }
                    roh_iter--;
                }
                
                if(minGrid_th > 0.90)
                    minGrid_th = 0.90;
            }
            else
                minimum_Th = 0.2;
            
            free(hist);
        }
        else
            minimum_Th = 0.2;
    }
    
    //minGrid_th = 0.1;
    printf("minimum TH %f\tminGrid_th %f\n",minimum_Th,minGrid_th);
    //exit(1);
    bool check_iter_end = false;
    
    double roh_th = 0.05;
    if(iteration == 1)
        roh_th = 0.10;
    
    if(Pyramid_step == 0)
    {
        if(Th_roh - 0.10 < Th_roh_min)
            check_iter_end  = true;
    }
    else if(Pyramid_step <= 1)
    {
        if(Th_roh - 0.10 < Th_roh_min)
            check_iter_end  = true;
    }
    else if(Pyramid_step == 2)
    {
        if(Th_roh - 0.10 < Th_roh_min)
            check_iter_end  = true;
    }
    else if(Pyramid_step == 3) 
    {
        if(Th_roh - 0.10 < Th_roh_min)
            check_iter_end  = true;
    }
    else 
    {
        if(proinfo->IsRA)
        {
            if(Th_roh - 0.10 < Th_roh_min)
                check_iter_end  = true;
        }
        else
        {
            if(Th_roh - 0.06 < Th_roh_min)
                check_iter_end  = true;
        }
    }

    double height_step = GetHeightStep(Pyramid_step, proinfo->resolution);
    
    if(proinfo->sensor_provider == PT)
        height_step = GetHeightStep_Planet(Pyramid_step, proinfo->resolution);
    
    for(long int iter_index = 0 ; iter_index < *rlevelinfo.Grid_length ; iter_index++)
    {
        long row     = (long)(floor(iter_index/rlevelinfo.Size_Grid2D->width));
        long col     = iter_index % rlevelinfo.Size_Grid2D->width;
        long grid_index = iter_index;//row*(long)rlevelinfo.Size_Grid2D->width + col;
        
        if(row >= 0 && row < rlevelinfo.Size_Grid2D->height && col >= 0 && col < rlevelinfo.Size_Grid2D->width && roh_height[grid_index].NumOfHeight > 2 && GridPT3[grid_index].ncc_seleceted_pair > -1)
        {
            bool index,index_1,index_2,index_3, roh_index;
            double ROR;
            
            if(Pyramid_step == 4 && iteration == 1)
                GridPT3[grid_index].Height          = -1000.0;
            
            index           = false;
            index_1         = false;
            index_2         = false;
            index_3         = false;
            roh_index       = false;
     
            if(SignedCharToDouble_grid(GridPT3[grid_index].roh) > minGrid_th)
                GridPT3[grid_index].roh = DoubleToSignedChar_grid(minGrid_th);
            
            //ratio of 1st peak roh / 2nd peak roh
            //double peak_step_diff = fabs(roh_height[grid_index].result2 - roh_height[grid_index].result3)/height_step;
            bool check_ROR = (iteration <= 2 && Pyramid_step >= 3) || (iteration <= 1 && Pyramid_step == 2);
            if(proinfo->sensor_provider == PT)
                check_ROR = (iteration <= 2 && Pyramid_step >= 2) || (iteration <= 1 && Pyramid_step == 1);
            
            if(check_ROR)
                ROR = 1.0;
            else
            {
                if(fabs(SignedCharToDouble_result(roh_height[grid_index].result0)) > 0)
                    ROR         = (SignedCharToDouble_result(roh_height[grid_index].result0) - SignedCharToDouble_result(roh_height[grid_index].result1))/SignedCharToDouble_result(roh_height[grid_index].result0);
                else
                    ROR         = 0;
            }
            
            if(Pyramid_step == 0)
            {
                if(ROR >= (0.1 - 0.03*(3 - Pyramid_step)) && SignedCharToDouble_result(roh_height[grid_index].result0) > minimum_Th)
                    index_2 = true;
            }
            else
            {
                if(ROR >= 0.1 && SignedCharToDouble_result(roh_height[grid_index].result0) > minimum_Th)
                    index_2 = true;
            }
            
            //if(!index_2 && peak_step_diff < 10 && SignedCharToDouble_result(roh_height[grid_index].result0) > minimum_Th)
            //    index_2 = true;
            
            // threshold of 1st peak roh
            if(SignedCharToDouble_result(roh_height[grid_index].result0) > SignedCharToDouble_grid(GridPT3[grid_index].roh) - roh_th)
                index_3     = true;

            if(roh_height[grid_index].result4 > 0)
                index       = true;
            
            if(MPP_stereo_angle > 5)
            {
                if(Pyramid_step >= 2)
                    roh_index   = index & index_2 & index_3;
                else if(Pyramid_step >= 1)
                {
                    if( index_3  && SignedCharToDouble_result(roh_height[grid_index].result0) > minimum_Th)
                        index_1 = true;
                    
                    roh_index   = ((index_3 & index_2) | index_1);
                }
                else if(Pyramid_step == 0 && final_level_iteration < 3)
                {
                    roh_index   = index & index_2 & index_3;
                }
                else
                {
                    if( index_3)
                        index_1 = true;
                    
                    roh_index   = (index_2 & index_1);
                }
            }
            else
            {
                if(Pyramid_step >= 1)
                    roh_index   = index & index_2 & index_3;
                else if(Pyramid_step == 0 && final_level_iteration < 3)
                {
                    roh_index   = index & index_2 & index_3;
                }
                else
                {
                    if( index_3)
                        index_1 = true;
                    
                    roh_index   = (index_2 & index_1);
                }
            }
            
            if(GridPT3[grid_index].Matched_flag != 0)
            {
                index           = false;
                index_1         = false;
                index_2         = false;
                
                if(Pyramid_step <= 0)
                {
                    if( (ROR < (0.1 - 0.03*(3 - Pyramid_step)) ) && (SignedCharToDouble_result(roh_height[grid_index].result1) > SignedCharToDouble_grid(GridPT3[grid_index].roh) - roh_th) )
                        index   = true;
                }
                else
                {
                    if( (ROR < 0.1) && (SignedCharToDouble_result(roh_height[grid_index].result1) > SignedCharToDouble_grid(GridPT3[grid_index].roh) - roh_th) )
                        index   = true;
                }
                
                if(roh_height[grid_index].result2 > roh_height[grid_index].result3)
                    index_1     = true;
                index_2         = !index_1;
                index_1         = index & index_1;
                index_2         = index & index_2;
                
                if(roh_height[grid_index].result3 > Nodata && roh_height[grid_index].result2 > Nodata)
                {
                    if(index_1)
                    {
                        GridPT3[grid_index].minHeight = floor(roh_height[grid_index].result3 - 0.5);
                        GridPT3[grid_index].maxHeight = ceil(roh_height[grid_index].result2 + 0.5);
                        GridPT3[grid_index].Matched_flag = 4;
                        
                    }

                    if(index_2)
                    {
                        GridPT3[grid_index].minHeight = floor(roh_height[grid_index].result2 - 0.5);
                        GridPT3[grid_index].maxHeight = ceil(roh_height[grid_index].result3 + 0.5);
                        GridPT3[grid_index].Matched_flag = 4;
                    }
                }
            }
        
            D3DPOINT temp_mp;
            {
                //Set the matched pts and information
                if(roh_index)
                {
                    double pre_H, pre_range;
                    uint8 pre_Match;
                    bool index_1, index_2, index_3, index_4, index_5, index_41, index_42, index_6, index_7, index;
                    index           = false;
                    index_1         = false;
                    index_2         = false;
                    index_3         = false;
                    index_4         = false;
                    index_41        = false;
                    index_42        = false;
                    index_5         = false;
                    index_6         = false;
                    index_7         = false;

                    temp_mp.m_Z = (double)roh_height[grid_index].result2;
                    temp_mp.m_X = rlevelinfo.GridPts[grid_index].m_X;
                    temp_mp.m_Y = rlevelinfo.GridPts[grid_index].m_Y;
                    temp_mp.m_roh = roh_height[grid_index].result0;
                    
                    //update MPs by previsous level results
                    if(Pyramid_step < 2)
                    {
                        pre_H           = (double)GridPT3[grid_index].Height;
                        pre_range       = (double)(GridPT3[grid_index].maxHeight - GridPT3[grid_index].minHeight);
                        pre_Match       = GridPT3[grid_index].Matched_flag;
                        
                        if(pre_Match == 0) //non matched grid
                            index_1     = true;
                        else if(pre_Match == 2) //matched grid
                            index_2     = true;
                        else if(pre_Match == 1 || pre_Match == 4) //inside triangle and min,max change status by matched grid
                            index_3     = true;
                        
                        
                        if(fabs(pre_H - temp_mp.m_Z) <= MPP*2*pwrtwo(Pyramid_step))
                            index_4     = true;
                        if(fabs(pre_H - temp_mp.m_Z) <= pre_range/2.0)
                            index_5     = true;
                        if(fabs(temp_mp.m_Z - GridPT3[grid_index].minHeight) <= MPP*3*pwrtwo(Pyramid_step))
                            index_41    = true;
                        if(fabs(temp_mp.m_Z - GridPT3[grid_index].maxHeight) <= MPP*3*pwrtwo(Pyramid_step))
                            index_42    = true;
                        
                        index_6         = (index_4 | index_5 | index_41 | index_42) & index_3; // inside triangle by matched grid
                        index_7         = index_2 & index_4; // matched grid
                        
                        index           = index_1 | index_7 | index_6;
                        
                        if(index && temp_mp.m_Z > -100 && temp_mp.m_Z < 10000)
                        {
                            count_MPs++;
                            
                            temp_mp.flag = 0;
                            linkedlist->push_back(temp_mp);
       
                            // update max_roh value
                            GridPT3[grid_index].roh     = roh_height[grid_index].result0;
                            if(GridPT3[grid_index].roh < DoubleToSignedChar_grid(minimum_Th))
                                GridPT3[grid_index].roh = DoubleToSignedChar_grid(minimum_Th);
                            
                        }
                        else
                        {
                            GridPT3[grid_index].roh = DoubleToSignedChar_grid(Th_roh);
                            if(check_iter_end && Th_roh_start == SignedCharToDouble_grid(GridPT3[grid_index].roh))
                                GridPT3[grid_index].roh     = DoubleToSignedChar_grid(Th_roh_next);
                        }
                    }
                    else
                    {
                        if(temp_mp.m_Z > -100 && temp_mp.m_Z < 10000)
                        {
                            count_MPs++;
                            
                            temp_mp.flag = 0;
                            linkedlist->push_back(temp_mp);
                        }
                        // update max_roh value
                        GridPT3[grid_index].roh     = roh_height[grid_index].result0;
                        if(GridPT3[grid_index].roh < DoubleToSignedChar_grid(minimum_Th))
                            GridPT3[grid_index].roh = DoubleToSignedChar_grid(minimum_Th);
                    }
                }
                else
                {
                    GridPT3[grid_index].roh = DoubleToSignedChar_grid(Th_roh);
                    // update max_roh value
                    if(check_iter_end && Th_roh_start == SignedCharToDouble_grid(GridPT3[grid_index].roh))
                        GridPT3[grid_index].roh     = DoubleToSignedChar_grid(Th_roh_next);
                }
            }
            
            if(GridPT3[grid_index].minHeight < Nodata)
                exit(1);
        }
    }

    return count_MPs;
}

void DecisionMPs(const ProInfo *proinfo, LevelInfo &rlevelinfo, const bool flag_blunder,const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, int *count_Results, double *minz_mp, double *maxz_mp, const double *minmaxHeight, D3DPOINT *ptslists)
{
    
    *minz_mp = 100000;
    *maxz_mp = -100000;

    if((*rlevelinfo.Pyramid_step == 0 && iteration == 3))
    {
        return;
    }

    const long int count_MPs       = count_MPs_input;
    const int Th_pts = 1;
 
    // Determine max count
    uint8 max_count         = 30;
    if(!flag_blunder) //anchor points
        max_count = 10;
    else //blunder points
    {
        if(*rlevelinfo.Pyramid_step >= 4)
            max_count  = 40;
        else if(*rlevelinfo.Pyramid_step >= 3)
            max_count  = 20;
        else
            max_count = 10;
    }
    
    long count = 0;
    double min_max[4] = {rlevelinfo.Boundary[0], rlevelinfo.Boundary[1], rlevelinfo.Boundary[2], rlevelinfo.Boundary[3]};
    
    UI3DPOINT *trilists = NULL;
    FullTriangulation *origTri = NULL;
    vector<UI3DPOINT> t_trilists;
    int count_tri;
    
    //Save triangulation for later use as we will remove blunders directly from this triangulation
    origTri = TINCreate_list(ptslists,count_MPs,&t_trilists,min_max,&count_tri, *rlevelinfo.grid_resolution);
    
    count_tri = t_trilists.size();
    
    //TODO why do we do this copy here?
    trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
    for(long i = 0 ; i < t_trilists.size() ; i++)
    {
        trilists[i] = t_trilists[i];
    }
    t_trilists.clear();
    vector<UI3DPOINT>().swap(t_trilists);
    
    
    bool flag               = true;
    double blunder_dh = 0;
    long pre_count_blunder = -100;
    long pre_count_MPs = -100;
    
    
    // runs between 10 and 30 times (by max_count at least)
    while(flag == true && count < max_count && count_tri > 20)
    {
        long int count_blunders = 0;
        long blunder_count[2] = {0};
        
        float* ortho_ncc = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
        
        bool *detBlunders = (bool*)calloc(sizeof(bool), count_MPs);
        
        if(proinfo->IsRA != 1)
        {
            //printf("DecisionMP : start VerticalLineLocus_blunder 1 iteration %d\n",count);
            SetHeightRange_blunder(rlevelinfo,ptslists, count_MPs, trilists,count_tri, GridPT3);
            VerticalLineLocus_blunder(proinfo, rlevelinfo, ortho_ncc, GridPT3, iteration, true);
            //printf("DecisionMP : end VerticalLineLocus_blunder 1\n");
        }
        
        //printf("blunder_detection_TIN start\n");
        blunder_detection_TIN(proinfo, rlevelinfo, iteration, ortho_ncc, flag_blunder, count, ptslists, detBlunders, count_MPs, trilists, count_tri, GridPT3, blunder_count,minz_mp,maxz_mp);
        //printf("blunder_detection_TIN end\n");
        
        free(ortho_ncc);
        
        if(count > 0)
            count_blunders = abs(int(blunder_count[1]) - pre_count_blunder);
        else
            count_blunders = blunder_count[1];
        
        if(count_blunders < Th_pts)
            flag = false;
        
        if(blunder_count[0] < Th_pts)
            flag = false;
        
        if(pre_count_blunder == blunder_count[1] && pre_count_MPs == blunder_count[0])
            flag = false;
        else
        {
            pre_count_MPs = blunder_count[0];
            pre_count_blunder = blunder_count[1];
        }
        
        //printf("start TIN\n");
        D3DPOINT *input_blunder_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),count_blunders);
        D3DPOINT *input_tri_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),blunder_count[0]);
        uint32 *check_id        = (uint32*)calloc(sizeof(uint32),blunder_count[0]);
        
        long int new_blunder_cnt = 0;
        long int t_tri_counts = 0;
        for(long int i=0;i<count_MPs;i++)
        {
            //Check for newly found blunders, and save them to be removed from triangulation
            if (detBlunders[i])
            {
                input_blunder_pts[new_blunder_cnt].m_X = ptslists[i].m_X;
                input_blunder_pts[new_blunder_cnt].m_Y = ptslists[i].m_Y;
                new_blunder_cnt++;
            }
            if(flag)
            {
                if(ptslists[i].flag != 1 && ptslists[i].flag != 2)
                {
                    input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
                    input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
                    check_id[t_tri_counts]          = i;
                    t_tri_counts++;
                }
            }
            else
            {
                if(ptslists[i].flag != 1)
                {
                    input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
                    input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
                    check_id[t_tri_counts]          = i;
                    t_tri_counts++;
                }
            }
        }
        
        if (new_blunder_cnt > 0)
        {
            free(trilists);
            vector<UI3DPOINT> tt_trilists;
            //If we have many blunders compared to points in triangulations, almost certainly faster to ditch old triangulation
            //Change the threshold by adding a scale factor to either t_tri_counts or t_blunder_counts
            if (TINUPD_THRSHLD*new_blunder_cnt > t_tri_counts)
            {
                //Must delete old triangulation and create new one, should be faster
                delete origTri;
                origTri = TINCreate_list(input_tri_pts, t_tri_counts, &tt_trilists, min_max, &count_tri, *rlevelinfo.grid_resolution);
            }
            else
            {
                //Rather than recreating entire triangulation, edit saved triangulation and only remove new blunders
                TINUpdate_list(input_tri_pts, t_tri_counts, &tt_trilists, min_max, &count_tri, *rlevelinfo.grid_resolution, origTri, input_blunder_pts, new_blunder_cnt);
            }
            
            count_tri = tt_trilists.size();
            
            trilists    = (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
            
            vector<UI3DPOINT>::iterator itt;
            for(long i = 0 ; i < tt_trilists.size() ; i++)
            //for(itt = tt_trilists.begin(); itt != tt_trilists.end() ; ++itt)
            {
                trilists[i].m_X = check_id[tt_trilists[i].m_X];
                trilists[i].m_Y = check_id[tt_trilists[i].m_Y];
                trilists[i].m_Z = check_id[tt_trilists[i].m_Z];
                //i++;
            }
            tt_trilists.clear();
            vector<UI3DPOINT>().swap(tt_trilists);
        }
        
        free(input_blunder_pts);
        free(input_tri_pts);
        free(detBlunders);
        free(check_id);
        
        //printf("end TIN\n");
        
        count_Results[0]    = count_MPs;
        count_Results[1]    = count_tri;
        
        count++;
        if(count >= max_count || count_tri <= 20)
            flag = false;
        
        //printf("iter = %d\tGridsize = %f\tMPs = %d\tBlunder = %d\tcount_tri = %d\tflag = %d\n",count,*rlevelinfo.grid_resolution,blunder_count[0],count_blunders,count_tri,flag);
    }

    SetHeightRange_blunder(rlevelinfo,ptslists, count_Results[0], trilists,count_tri, GridPT3);
   
    free(trilists);
    
    float* ortho_ncc = (float*)calloc(sizeof(float),*rlevelinfo.Grid_length);
    
    //printf("DecisionMP : start VerticalLineLocus_blunder 2\n");
    VerticalLineLocus_blunder(proinfo, rlevelinfo, ortho_ncc, GridPT3, iteration, false);
    //printf("DecisionMP : end VerticalLineLocus_blunder 2\n");
    
    free(ortho_ncc);
    delete origTri;
}

void DecisionMPs_setheight(const ProInfo *proinfo, LevelInfo &rlevelinfo, const long int count_MPs_input,UGRID *GridPT3, const uint8 iteration, const double Hinterval, const double *minmaxHeight, D3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri)
{
    const int count_MPs       = count_MPs_input;
    const double gridspace            = *rlevelinfo.grid_resolution;
    
    BL blunder_param;
    blunder_param.Boundary  = rlevelinfo.Boundary;
    blunder_param.gridspace = *rlevelinfo.grid_resolution;
    blunder_param.height_check_flag = true;
    blunder_param.Hinterval = Hinterval;
    blunder_param.iteration = iteration;
    blunder_param.Pyramid_step = *rlevelinfo.Pyramid_step;
    
    blunder_param.Size_Grid2D.width = rlevelinfo.Size_Grid2D->width;
    blunder_param.Size_Grid2D.height = rlevelinfo.Size_Grid2D->height;
    
    float* ortho_ncc = (float*)calloc(*rlevelinfo.Grid_length,sizeof(float));
    
    SetHeightRange_blunder(rlevelinfo,ptslists, count_MPs, trilists,numoftri, GridPT3);
    
    VerticalLineLocus_blunder(proinfo, rlevelinfo, ortho_ncc, GridPT3, iteration, false);

    free(ortho_ncc);

}

void set_blunder(long int index, uint8_t val, D3DPOINT *pts, bool *detectedBlunders) {
    // Only update the detectedBlunders value
    // if we know that this is a new blunder.
    // if it was previously zero, we know it is
    // a new blunder. If it was not previously zero,
    // there are two possible cases. First, it
    // was already a blunder. In that case, we don't
    // want to set detectedBlunders. Alternatively,
    // it could be a new blunder but another thread
    // changed it from zero to 1 or 3. In that case,
    // the other thread/iteration would have
    // updated detectedBlunders, so we don't need to.
    uint8_t prev;

#pragma omp atomic read
    prev = pts[index].flag;
    
    // if the flag is already 1 and it's an "old blunder",
    // then we don't want to update it to 3. This can cause
    // nondeterminism in the loop. However, the race here is
    // okay. If it's an "old blunder", we'll always see 1 here.
    // If it's a new blunder, we may see any of 0, 1 or 3 here.
    // If it's a new blunder and 1, then detectedBlunders was
    // already set, so we can exit early. Same with 3. If it's
    // a new blunder and zero, but there's a race, that's fine
    // too.
    if(prev != 0)
        return;

#pragma omp atomic capture
    {
        prev = pts[index].flag;
        pts[index].flag = val;
    }
    if(prev == 0)
    {
#pragma omp atomic write
        detectedBlunders[index] = true;
    }
}

bool blunder_detection_TIN(const ProInfo *proinfo, LevelInfo &rlevelinfo, const int iteration, float* ortho_ncc, bool flag_blunder, uint16 count_bl, D3DPOINT *pts, bool *detectedBlunders, long int num_points, UI3DPOINT *tris, long int num_triangles, UGRID *Gridpts, long *blunder_count,double *minz_mp, double *maxz_mp)
{
    int IsRA(proinfo->IsRA);
    const uint8 pyramid_step(*rlevelinfo.Pyramid_step);
    double gridspace = *rlevelinfo.grid_resolution;
    
    if(!(pyramid_step == 0 && iteration == 3))
    {
        uint32 hdiffcount = (uint32)(*rlevelinfo.Hinterval);
        const uint8 max_nodes       = 30;
        //savenode will now be a pointer to a uint
        //make savenode a giant array with row size max_nodes
        uint32 *savenode    = (uint32*)malloc((long)sizeof(uint32)*(long)max_nodes*(long)num_points);
        uint32 *nodecount   = (uint32*)calloc(num_points,sizeof(uint32));
        uint32 *hdiffbin    = (uint32*)calloc(hdiffcount+1,sizeof(uint32));
        
        const double *boundary    = rlevelinfo.Boundary;
        const CSize gridsize(rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height);
        double sum_oncc = 0;
        double sum2_oncc = 0;
        long int total_oncc_count = 0;
        double sum_data2    = 0.0;
        double sum_data = 0.0;
        long int dh_count = 0;
        
        //compute dh statistics
        const double ortho_ncc_thread = -1;
        for(long tcnt=0;tcnt<num_triangles;tcnt++)
        {
            bool check_pt_index = true;
            
            D3DPOINT pt0,pt1,pt2;
            int t_col,t_row;
            
            if(tris[tcnt].m_X < num_points)
            {
                if(nodecount[tris[tcnt].m_X] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_X*(long)max_nodes+(long)nodecount[tris[tcnt].m_X]] = tcnt;
                    nodecount[tris[tcnt].m_X]++;
                }
                pt0     = pts[tris[tcnt].m_X];
                
                t_col         = (int)((pt0.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt0.m_Y - boundary[1])/gridspace + 0.5);
                int pt0_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt0_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt0_index];
                    sum2_oncc += ortho_ncc[pt0_index]* ortho_ncc[pt0_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            if(tris[tcnt].m_Y < num_points)
            {
                if(nodecount[tris[tcnt].m_Y] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_Y*(long)max_nodes+(long)nodecount[tris[tcnt].m_Y]] = tcnt;
                    nodecount[tris[tcnt].m_Y]++;
                }
                pt1     = pts[tris[tcnt].m_Y];
                
                t_col         = (int)((pt1.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt1.m_Y - boundary[1])/gridspace + 0.5);
                int pt1_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt1_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt1_index];
                    sum2_oncc += ortho_ncc[pt1_index]* ortho_ncc[pt1_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            if(tris[tcnt].m_Z < num_points)
            {
                if(nodecount[tris[tcnt].m_Z] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_Z*(long)max_nodes+(long)nodecount[tris[tcnt].m_Z]] = tcnt;
                    nodecount[tris[tcnt].m_Z]++;
                }
                pt2     = pts[tris[tcnt].m_Z];
                
                t_col         = (int)((pt2.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt2.m_Y - boundary[1])/gridspace + 0.5);
                int pt2_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt2_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt2_index];
                    sum2_oncc += ortho_ncc[pt2_index]* ortho_ncc[pt2_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            double dh1,dh2,dh3;
            if(check_pt_index)
            {
                dh1     = (double)fabs(pt0.m_Z - pt1.m_Z);
                sum_data += dh1;
                sum_data2 += dh1*dh1;
                dh_count++;
                
                dh2     = (double)fabs(pt0.m_Z - pt2.m_Z);
                sum_data += dh2;
                sum_data2 += dh2*dh2;
                dh_count++;
                
                dh3     = (double)fabs(pt1.m_Z - pt2.m_Z);
                sum_data += dh3;
                sum_data2 += dh3*dh3;
                dh_count++;
                
                dh1   = (uint32)dh1;
                if (dh1 < hdiffcount)
                    hdiffbin[(uint32)dh1]++;
                
                dh2   = (uint32)dh2;
                if (dh2 < hdiffcount)
                    hdiffbin[(uint32)dh2]++;
                
                dh3   = (uint32)dh3;
                if (dh3 < hdiffcount)
                    hdiffbin[(uint32)dh3]++;
            }
        }
        
        //set ortho_ncc
        const double oncc_mean    =  sum_oncc/total_oncc_count;
        const double oncc_std =   sqrt(fabs(sum2_oncc - (sum_oncc)*(sum_oncc)/total_oncc_count)/total_oncc_count);
        
        double ortho_ncc_th;// = 0.7 + (iteration-1)*0.02;
        
        if(pyramid_step == 4 )
            ortho_ncc_th = 0.6;// - (iteration - 1)*0.01;
        else if(pyramid_step >= 3)
            ortho_ncc_th = 0.5;// - (iteration - 1)*0.01;
        else if(pyramid_step == 2)
            ortho_ncc_th = 0.4;// - (iteration - 1)*0.01;
        else if(pyramid_step == 1)
            ortho_ncc_th = 0.3 ;
        else
            ortho_ncc_th = 0.2 ;
        
        double temp_oncc_th = oncc_mean - 1.5*oncc_std;
        if(temp_oncc_th < 0)
            temp_oncc_th = oncc_mean;
        if(temp_oncc_th > 0)
        {
            if(temp_oncc_th < ortho_ncc_th)
                ortho_ncc_th = temp_oncc_th;
            else
            {
                temp_oncc_th = oncc_mean - 3*oncc_std;
                if(temp_oncc_th < 0)
                    temp_oncc_th = oncc_mean;
                if(temp_oncc_th < ortho_ncc_th)
                    ortho_ncc_th = temp_oncc_th;
            }
        }
        if(ortho_ncc_th < 0.2)
            ortho_ncc_th = 0.2;
         
         double ortho_ancc_th = oncc_mean + 1.0*oncc_std;//0.9 - 0.1*(4-pyramid_step);//100.;
         if(ortho_ancc_th > 0.95)
             ortho_ancc_th = 0.95;
         
         double th_ref_ncc = 0.1;// + (iteration-1)*0.05;
         if(th_ref_ncc > ortho_ncc_th)
             th_ref_ncc = ortho_ncc_th;
        
        /*
        double ortho_ncc_th = 0.5 + (iteration-1)*0.02;
        double temp_oncc_th = oncc_mean - oncc_std;
        if(temp_oncc_th < ortho_ncc_th)
            ortho_ncc_th = temp_oncc_th;
        else
        {
            temp_oncc_th = oncc_mean - 2*oncc_std;
            if(temp_oncc_th < ortho_ncc_th)
                ortho_ncc_th = temp_oncc_th;
        }
        
        if(pyramid_step == 4 )
            ortho_ncc_th = 0.6 - (iteration - 1)*0.01;
        else if(pyramid_step >= 3)
            ortho_ncc_th = 0.5 - (iteration - 1)*0.01;
        else if(pyramid_step == 2)
            ortho_ncc_th = 0.4 - (iteration - 1)*0.01;
        else if(pyramid_step == 1)
            ortho_ncc_th = 0.3 ;
        else
            ortho_ncc_th = 0.2 ;
 
        const double ortho_ancc_th = 100.;
        double th_ref_ncc = 0.1 + (iteration-1)*0.05;
        if(th_ref_ncc > ortho_ncc_th)
            th_ref_ncc = ortho_ncc_th;
        */
        
        
        
        
        //set height_th
        const int blunder_pyramid_step = 3;
        bool check_dh       = false;
        if(pyramid_step <= 4 && pyramid_step >= 3)
            check_dh        = true;
        else if(pyramid_step == 2)
        {
            if(iteration <= 4)
                check_dh        = true;
        }
        else if(pyramid_step <= 1)
        {
            if(iteration <= 2)
                check_dh        = true;
        }
        
        double height_th;
        double lw_3sigma, up_3sigma;
        bool check_total_dh = false;
        double GSD;
        double mean,std;
        uint32 th80 = 1000;
        
        if(pyramid_step >= 1)
        {
            GSD = gridspace;
            if(pyramid_step >= 3 && gridspace < proinfo->resolution*pwrtwo(pyramid_step))
                GSD = proinfo->resolution*pwrtwo(pyramid_step);
        }
        else
            GSD = proinfo->resolution*pwrtwo(pyramid_step);
        
        const double height_th_1(*rlevelinfo.Hinterval);
        double height_th_2(GSD*10);
        
        if(IsRA == 1)
            height_th_2= GSD*30;
        
        if(height_th_1 < height_th_2)
            height_th       = height_th_1;
        else
            height_th       = height_th_2;
        
        if(pyramid_step >= blunder_pyramid_step)
        {
            if(check_dh)
            {
                mean    =   sum_data/dh_count;
                std =   sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
                lw_3sigma       = mean - 3*std;
                up_3sigma       = mean + 3*std;
            }
            else
            {
                lw_3sigma       = -100;
                up_3sigma       = height_th;
            }
            
            if(up_3sigma < height_th)
            {
                if(up_3sigma < height_th/2.0)
                    height_th = height_th/2.0;
                else
                    height_th = up_3sigma;
            }
        }
        else
        {
            uint32 total_dh = 0;
            long tcnt = 0;
            while(tcnt<hdiffcount && !check_total_dh)
            {
                double per;
                total_dh += hdiffbin[tcnt];
                per = (double)total_dh/(double)dh_count*100.0;
                if (per > 99)
                {
                    th80 = tcnt;
                    check_total_dh = true;
                }
                tcnt++;
            }
            if (!check_total_dh)
                th80 = tcnt -1;
            
            mean    =   sum_data/dh_count;
            std =   sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
            lw_3sigma       = mean - 3*std;
            up_3sigma       = mean + 3*std;
            
            if(up_3sigma < th80)
                height_th = up_3sigma;
            else {
                height_th = th80;
            }
            
            if(height_th > 50)
                height_th = 50;
        }
        
        free(hdiffbin);

        // be very careful modifying this code. The thread safety
        // here is a bit complicated. This loop touches three
        // shared arrays: pts, detectedBlunders, and GridPts.
        // GridPts and detecetedBlunders are write-only, while
        // pts is read/write. GridPts is the least problematic.
        // It is written is a "safe" way, with one-to-one mapping
        // between iterations and the index.
        //
        // That's not the case for pts and detectedBlunders.
        // pts is only read from the iteration index. But, it is
        // written from other indices. detectedBlunders tracks
        // whether a given pts index was updated from zero to
        // one or three. It is also written from other threads.
        //
        // openmp atomics are used to coordinate this. Take
        // care when modifying anything in here.
        //
        // detectedBlunders is updated when for an index when
        // the pts value is changed from 0 to 1 or from 0 to 3.
#pragma omp parallel for schedule(guided)
        for(long index=0;index<num_points;index++)
        {
            if(pts[index].flag != 1)
            {
                // use this flag instead of setting the point directly
                bool pt_is_blunder = false;

                int count_th_positive   = 0;
                int count_th_negative   = 0;
                int count = 0;
                int max_iter;
                
                bool check_neigh = false;
                
                if(nodecount[index] < max_nodes)
                    max_iter    = nodecount[index];
                else
                    max_iter    = max_nodes;
                
                const D3DPOINT ref_index_pt(pts[index]);
                long t_col         = (long)((ref_index_pt.m_X - boundary[0])/gridspace + 0.5);
                long t_row         = (long)((ref_index_pt.m_Y - boundary[1])/gridspace + 0.5);
                const long ref_index((long)gridsize.width*t_row + t_col);
                
                // assume that the mapping between index
                // and ref_index is one-to-one, so this
                // is okay.
                Gridpts[ref_index].anchor_flag = 0;
                
                if(!IsRA)
                {
                    if(ortho_ncc[ref_index] < th_ref_ncc && pyramid_step >= 2)
                        pt_is_blunder = true;
                }
                /*
                if(pyramid_step >= 3 && iteration >= 4)
                {
                    if(ortho_ncc[ref_index] >= 0.95 && flag_blunder)
                    {
                        // Assume this is okay
                        Gridpts[ref_index].anchor_flag = 1;
                    }
                }
                */
                for(long iter = 0 ; iter < max_iter ; iter++)
                {
                    if(savenode[(long)index*(long)max_nodes+(long)iter] < num_triangles && savenode[(long)index*(long)max_nodes+(long)iter] > 0)
                    {
                        long reference_index = 0;
                        long target_index_0 = 0, target_index_1 = 0;
                        const uint32 temp_tri[3] = {tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_X, tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Y, tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Z};
                        
                        bool check_index = false;
                        for(int kk=0;kk<3;kk++)
                        {
                            if(temp_tri[kk] == index)
                                reference_index = index;
                            else
                            {
                                if(!check_index)
                                {
                                    target_index_0 = temp_tri[kk];
                                    check_index    = true;
                                }
                                else
                                    target_index_1 = temp_tri[kk];
                            }
                        }
                        
                        if(reference_index < num_points && target_index_0 < num_points && target_index_1 < num_points &&
                           target_index_0 >= 0 && target_index_1 >= 0)
                        {
                            //plane normal angle
                            const D3DPOINT pt0(pts[reference_index]);
                            const D3DPOINT pt1(pts[target_index_0]);
                            const D3DPOINT pt2(pts[target_index_1]);
                            
                            double angle = SetNormalAngle(pt0, pt1, pt2);
                            
                            const double dh1 = pt0.m_Z - pt1.m_Z;
                            const double dh2 = pt0.m_Z - pt2.m_Z;
                            const double dh3 = pt1.m_Z - pt2.m_Z;
                            
                            if(dh1 >= 0 && dh2 >= 0 && angle > 30)
                                count_th_positive ++;
                            if(dh1 <  0 && dh2 < 0  && angle > 30)
                                count_th_negative ++;
                            
                            bool check_match = false;
                            if(pyramid_step  > 1)
                            {
                                check_match     = true;
                            }
                            else if(pyramid_step == 1)
                            {
                                if(iteration == 1 && count_bl <= 5)
                                    check_match = true;
                            }
                            else if(pyramid_step == 0)
                            {
                                if(iteration == 1 && count_bl <= 5)
                                    check_match = true;
                                else if(iteration == 2 && count_bl <= 5)
                                    check_match = false;
                            }
                            else
                                check_match = false;
                            
                            if(check_match)
                            {
                                const double ddh_1     = fabs(dh1) - fabs(dh2);
                                const double ddh_2     = fabs(dh1) - fabs(dh3);
                                const double ddh_3     = fabs(dh2) - fabs(dh3);
                                
                                if((fabs(ddh_1) > height_th || fabs(ddh_2) > height_th || fabs(ddh_3) > height_th))
                                {
                                    // all node check with height difference.
                                    double h1,h2,dh;
                                    int t_o_min,t_o_max,t_o_mid;
                                    long order[3]    = {reference_index, target_index_0, target_index_1};
                                    double height[3] = {pt0.m_Z, pt1.m_Z, pt2.m_Z};
                                    
                                    long Col         = (long)((pt0.m_X - boundary[0])/gridspace + 0.5);
                                    long Row         = (long)((pt0.m_Y- boundary[1])/gridspace + 0.5);
                                    long Index[3];
                                    Index[0]     = (long)gridsize.width*Row + Col;
                                    Col             = (long)((pt1.m_X - boundary[0])/gridspace + 0.5);
                                    Row             = (long)((pt1.m_Y - boundary[1])/gridspace + 0.5);
                                    Index[1]     = (long)gridsize.width*Row + Col;
                           
                                    Col             = (long)((pt2.m_X - boundary[0])/gridspace + 0.5);
                                    Row             = (long)((pt2.m_Y - boundary[1])/gridspace + 0.5);
                                    Index[2]     = (long)gridsize.width*Row + Col;
                                    
                                    check_neigh = true;
                                    
                                    if(height[0] > height[1])
                                    {
                                        if(height[0] > height[2])
                                        {
                                            t_o_max = 0;
                                            if(height[1] > height[2])
                                            {
                                                t_o_min = 2;
                                                t_o_mid = 1;
                                            }
                                            else
                                            {
                                                t_o_min = 1;
                                                t_o_mid = 2;
                                            }
                                        }
                                        else
                                        {
                                            t_o_max = 2;
                                            t_o_min = 1;
                                            t_o_mid = 0;
                                        }
                                    }
                                    else
                                    {
                                        if(height[1] > height[2])
                                        {
                                            t_o_max = 1;
                                            if(height[0] > height[2])
                                            {
                                                t_o_min = 2;
                                                t_o_mid = 0;
                                            }
                                            else
                                            {
                                                t_o_min = 0;
                                                t_o_mid = 2;
                                            }
                                        }
                                        else
                                        {
                                            t_o_max = 2;
                                            t_o_mid = 1;
                                            t_o_min = 0;
                                        }
                                    }
                                    
                                    h1        = height[t_o_mid] - height[t_o_min];
                                    h2        = height[t_o_max] - height[t_o_mid];
                                    dh        = h1 - h2;

                                    long int blunder_neighbor_index = -1;
                                    if((dh > 0 && dh > height_th))
                                    {
                                        if(IsRA == 1)
                                        {
                                            blunder_neighbor_index = order[t_o_min];
                                        }
                                        else
                                        {
                                            if(flag_blunder)
                                            {
                                                if(ortho_ncc[Index[t_o_min]] < ortho_ncc_th)
                                                    blunder_neighbor_index = order[t_o_min];
                                            }
                                            else
                                                if(ortho_ncc[Index[t_o_min]] < ortho_ancc_th)
                                                    blunder_neighbor_index = order[t_o_min];
                                            
                                        }
                                    }
                                    else if((dh < 0 && fabs(dh) > height_th))
                                    {
                                        if(IsRA == 1)
                                        {
                                            blunder_neighbor_index = order[t_o_max];
                                        }
                                        else
                                        {
                                            if(flag_blunder)
                                            {
                                                if(ortho_ncc[Index[t_o_max]] < ortho_ncc_th)
                                                    blunder_neighbor_index = order[t_o_max];
                                            }
                                            else
                                                if(ortho_ncc[Index[t_o_max]] < ortho_ancc_th)
                                                    blunder_neighbor_index = order[t_o_max];
                                        }
                                    }

                                    if(blunder_neighbor_index >= 0)
                                        set_blunder(blunder_neighbor_index, 3, pts, detectedBlunders);
                                }
                            }
                            count++;
                        }
                    }
                }
                
                if(IsRA == 1)
                {
                    if(check_neigh == true)
                    {
                        if(!flag_blunder)
                            pt_is_blunder = true;
                    }
      
                    if((count_th_positive >= (int)(count*0.7 + 0.5) || count_th_negative >= (int)(count*0.7 + 0.5)) && count > 1)
                        pt_is_blunder = true;
                }
                else
                {
                    if(check_neigh == true)
                    {
                        if(!flag_blunder)
                        {
                            if(pyramid_step >= 3 && iteration == 4)
                            {
                                if(ortho_ncc[ref_index] < ortho_ancc_th)
                                    pt_is_blunder = true;
                            }
                            else
                                if(ortho_ncc[ref_index] < ortho_ancc_th)
                                    pt_is_blunder = true;
                        }
                    }
                    
                    double count_threshold;
                    count_threshold = count*0.5;
                    if(pyramid_step >= 2)
                        count_threshold = count*0.5;
                    else if(pyramid_step == 1)
                        count_threshold = (int)(count*0.6 + 0.5);
                    else
                        count_threshold = (int)(count*0.7 + 0.5);
                    
                    if((count_th_positive >= count_threshold || count_th_negative >= count_threshold) && count > 1)
                    {
                        if(flag_blunder)
                        {
                            double tmp_th = 0.6 - (4-pyramid_step)*0.1;
                            if(pyramid_step == 1)
                                tmp_th = 0.6;
                            
                            if(pyramid_step >= 1)
                            {
                                if(ortho_ncc[ref_index] < tmp_th)
                                    pt_is_blunder = true;
                            }
                            else if(pyramid_step == 0)
                            {
                                if(iteration <= 1)
                                {
                                    if(ortho_ncc[ref_index] < 0.9)
                                        pt_is_blunder = true;
                                }
                                else if(iteration == 2 )
                                    pt_is_blunder = true;
                                else {
                                    if(ortho_ncc[ref_index] < ortho_ncc_th)
                                        pt_is_blunder = true;
                                }
                            }
                        }
                        else
                        {
                            if(ortho_ncc[ref_index] < ortho_ancc_th)
                                pt_is_blunder = true;
                        }
                    }
                }
                // only set detectedBlunders[index] to true for the index associated
                // with this loop iteration. The neighbor points set detectedBlunders
                // above as appropriate. We have to use the private variable instead
                // of just reading the array value because other threads may set
                // the value to 3 as we're processing.
                if(pt_is_blunder)
                    set_blunder(index, 1, pts, detectedBlunders);
            }
        }
        
        for(long tcnt=0;tcnt<num_points;tcnt++)
        {
            if(pts[tcnt].flag == 1)
                blunder_count[1]++;
            else if(pts[tcnt].flag == 3)
            {
                blunder_count[1]++;
                pts[tcnt].flag = 1;
            }
            else
            {
                blunder_count[0]++;
                
                if(*minz_mp > pts[tcnt].m_Z)
                    *minz_mp        = pts[tcnt].m_Z;
                if(*maxz_mp < pts[tcnt].m_Z)
                    *maxz_mp        = pts[tcnt].m_Z;
            }
        }
        free(savenode);
        free(nodecount);
    }
  
    return true;
}

/*
bool blunder_detection_TIN(const ProInfo *proinfo, LevelInfo &rlevelinfo, const int iteration, float* ortho_ncc, bool flag_blunder, uint16 count_bl, D3DPOINT *pts, bool *detectedBlunders, long int num_points, UI3DPOINT *tris, long int num_triangles, UGRID *Gridpts, long *blunder_count,double *minz_mp, double *maxz_mp)
{
    int IsRA(proinfo->IsRA);
    const uint8 pyramid_step(*rlevelinfo.Pyramid_step);
    double gridspace = *rlevelinfo.grid_resolution;
    
    if(!(pyramid_step == 0 && iteration == 3))
    {
        uint32 hdiffcount = (uint32)(*rlevelinfo.Hinterval);
        const uint8 max_nodes       = 30;
        //savenode will now be a pointer to a uint
        //make savenode a giant array with row size max_nodes
        uint32 *savenode    = (uint32*)malloc((long)sizeof(uint32)*(long)max_nodes*(long)num_points);
        uint32 *nodecount   = (uint32*)calloc(num_points,sizeof(uint32));
        uint32 *hdiffbin    = (uint32*)calloc(hdiffcount+1,sizeof(uint32));
        
        const double *boundary    = rlevelinfo.Boundary;
        const CSize gridsize(rlevelinfo.Size_Grid2D->width, rlevelinfo.Size_Grid2D->height);
        double sum_oncc = 0;
        double sum2_oncc = 0;
        long int total_oncc_count = 0;
        double sum_data2    = 0.0;
        double sum_data = 0.0;
        long int dh_count = 0;
        
        //compute dh statistics
        const double ortho_ncc_thread = 0;
        for(long tcnt=0;tcnt<num_triangles;tcnt++)
        {
            bool check_pt_index = true;
            
            D3DPOINT pt0,pt1,pt2;
            int t_col,t_row;
            
            if(tris[tcnt].m_X < num_points)
            {
                if(nodecount[tris[tcnt].m_X] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_X*(long)max_nodes+(long)nodecount[tris[tcnt].m_X]] = tcnt;
                    nodecount[tris[tcnt].m_X]++;
                }
                pt0     = pts[tris[tcnt].m_X];
                
                t_col         = (int)((pt0.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt0.m_Y - boundary[1])/gridspace + 0.5);
                int pt0_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt0_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt0_index];
                    sum2_oncc += ortho_ncc[pt0_index]* ortho_ncc[pt0_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            if(tris[tcnt].m_Y < num_points)
            {
                if(nodecount[tris[tcnt].m_Y] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_Y*(long)max_nodes+(long)nodecount[tris[tcnt].m_Y]] = tcnt;
                    nodecount[tris[tcnt].m_Y]++;
                }
                pt1     = pts[tris[tcnt].m_Y];
                
                t_col         = (int)((pt1.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt1.m_Y - boundary[1])/gridspace + 0.5);
                int pt1_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt1_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt1_index];
                    sum2_oncc += ortho_ncc[pt1_index]* ortho_ncc[pt1_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            if(tris[tcnt].m_Z < num_points)
            {
                if(nodecount[tris[tcnt].m_Z] < max_nodes)
                {
                    savenode[(long)tris[tcnt].m_Z*(long)max_nodes+(long)nodecount[tris[tcnt].m_Z]] = tcnt;
                    nodecount[tris[tcnt].m_Z]++;
                }
                pt2     = pts[tris[tcnt].m_Z];
                
                t_col         = (int)((pt2.m_X - boundary[0])/gridspace + 0.5);
                t_row         = (int)((pt2.m_Y - boundary[1])/gridspace + 0.5);
                int pt2_index     = gridsize.width*t_row + t_col;
                if( ortho_ncc[pt2_index] > ortho_ncc_thread)
                {
                    sum_oncc += ortho_ncc[pt2_index];
                    sum2_oncc += ortho_ncc[pt2_index]* ortho_ncc[pt2_index];
                    total_oncc_count++;
                }
            }
            else
                check_pt_index = false;
            
            double dh1,dh2,dh3;
            if(check_pt_index)
            {
                dh1     = (double)fabs(pt0.m_Z - pt1.m_Z);
                sum_data += dh1;
                sum_data2 += dh1*dh1;
                dh_count++;
                
                dh2     = (double)fabs(pt0.m_Z - pt2.m_Z);
                sum_data += dh2;
                sum_data2 += dh2*dh2;
                dh_count++;
                
                dh3     = (double)fabs(pt1.m_Z - pt2.m_Z);
                sum_data += dh3;
                sum_data2 += dh3*dh3;
                dh_count++;
                
                dh1   = (uint32)dh1;
                if (dh1 < hdiffcount)
                    hdiffbin[(uint32)dh1]++;
                
                dh2   = (uint32)dh2;
                if (dh2 < hdiffcount)
                    hdiffbin[(uint32)dh2]++;
                
                dh3   = (uint32)dh3;
                if (dh3 < hdiffcount)
                    hdiffbin[(uint32)dh3]++;
            }
        }
        
        //set ortho_ncc
        const double oncc_mean    =  sum_oncc/total_oncc_count;
        const double oncc_std =   sqrt(fabs(sum2_oncc - (sum_oncc)*(sum_oncc)/total_oncc_count)/total_oncc_count);
        
        //printf("oncc mean %f\tstd%f\n",oncc_mean,oncc_std);
        
        double ortho_ncc_th;// = 0.7 + (iteration-1)*0.02;
        
        if(pyramid_step == 4 )
            ortho_ncc_th = 0.6;// - (iteration - 1)*0.01;
        else if(pyramid_step >= 3)
            ortho_ncc_th = 0.5;// - (iteration - 1)*0.01;
        else if(pyramid_step == 2)
            ortho_ncc_th = 0.4;// - (iteration - 1)*0.01;
        else if(pyramid_step == 1)
            ortho_ncc_th = 0.3 ;
        else
            ortho_ncc_th = 0.2 ;
        
        double temp_oncc_th = oncc_mean - 1.5*oncc_std;
        if(temp_oncc_th < 0)
            temp_oncc_th = oncc_mean;
        if(temp_oncc_th > 0)
        {
            if(temp_oncc_th < ortho_ncc_th)
                ortho_ncc_th = temp_oncc_th;
            else
            {
                temp_oncc_th = oncc_mean - 3*oncc_std;
                if(temp_oncc_th < 0)
                    temp_oncc_th = oncc_mean;
                if(temp_oncc_th < ortho_ncc_th)
                    ortho_ncc_th = temp_oncc_th;
            }
        }
        if(ortho_ncc_th < 0.2)
            ortho_ncc_th = 0.2;
        
        double ortho_ancc_th = oncc_mean + 1.0*oncc_std;//0.9 - 0.1*(4-pyramid_step);//100.;
        if(ortho_ancc_th > 0.95)
            ortho_ancc_th = 0.95;
        
        double th_ref_ncc = 0.1;// + (iteration-1)*0.05;
        if(th_ref_ncc > ortho_ncc_th)
            th_ref_ncc = ortho_ncc_th;
        
        //printf("ortho_ncc_th %f\t ortho_ancc_th %f\t th_ref_ncc %f\n",ortho_ncc_th,ortho_ancc_th,temp_oncc_th);
        
        //set height_th
        const int blunder_pyramid_step = 3;
        bool check_dh       = false;
        if(pyramid_step <= 4 && pyramid_step >= 3)
            check_dh        = true;
        else if(pyramid_step == 2)
        {
            if(iteration <= 4)
                check_dh        = true;
        }
        else if(pyramid_step <= 1)
        {
            if(iteration <= 2)
                check_dh        = true;
        }
        
        double height_th;
        double lw_3sigma, up_3sigma;
        bool check_total_dh = false;
        double GSD;
        double mean,std;
        uint32 th80 = 1000;
        
        if(pyramid_step >= 1)
        {
            GSD = gridspace;
            if(pyramid_step >= 3 && gridspace < proinfo->resolution*pow(2.0,pyramid_step))
                GSD = proinfo->resolution*pow(2.0,pyramid_step);
        }
        else
            GSD = proinfo->resolution*pow(2.0,pyramid_step);
        
        const double height_th_1(*rlevelinfo.Hinterval);
        double height_th_2(GSD*10);
        
        if(IsRA == 1)
            height_th_2= GSD*30;
        
        if(height_th_1 < height_th_2)
            height_th       = height_th_1;
        else
            height_th       = height_th_2;
        
        if(pyramid_step >= blunder_pyramid_step)
        {
            if(check_dh)
            {
                mean    =   sum_data/dh_count;
                std =   sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
                lw_3sigma       = mean - 3*std;
                up_3sigma       = mean + 3*std;
            }
            else
            {
                lw_3sigma       = -100;
                up_3sigma       = height_th;
            }
            
            if(up_3sigma < height_th)
            {
                if(up_3sigma < height_th/2.0)
                    height_th = height_th/2.0;
                else
                    height_th = up_3sigma;
            }
        }
        else
        {
            uint32 total_dh = 0;
            long tcnt = 0;
            while(tcnt<hdiffcount && !check_total_dh)
            {
                double per;
                total_dh += hdiffbin[tcnt];
                per = (double)total_dh/(double)dh_count*100.0;
                if (per > 99)
                {
                    th80 = tcnt;
                    check_total_dh = true;
                }
                tcnt++;
            }
            if (!check_total_dh)
                th80 = tcnt -1;
            
            mean    =   sum_data/dh_count;
            std =   sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
            lw_3sigma       = mean - 3*std;
            up_3sigma       = mean + 3*std;
            
            if(up_3sigma < th80)
                height_th = up_3sigma;
            else {
                height_th = th80;
            }
            
            if(height_th > 50)
                height_th = 50;
        }
        
        free(hdiffbin);
        
        //printf("height_th %f\titeration %d\n",height_th,iteration);
        
#pragma omp parallel for schedule(guided)
        for(long index=0;index<num_points;index++)
        {
            if(pts[index].flag != 1)
            {
                int count_th_positive   = 0;
                int count_th_negative   = 0;
                int count = 0;
                int max_iter;
                
                bool check_neigh = false;
                
                if(nodecount[index] < max_nodes)
                    max_iter    = nodecount[index];
                else
                    max_iter    = max_nodes;
                
                const D3DPOINT ref_index_pt(pts[index]);
                long t_col         = (long)((ref_index_pt.m_X - boundary[0])/gridspace + 0.5);
                long t_row         = (long)((ref_index_pt.m_Y - boundary[1])/gridspace + 0.5);
                const long ref_index((long)gridsize.width*t_row + t_col);
                
                Gridpts[ref_index].anchor_flag = 0;
                
                if(!IsRA)
                {
                    if(ortho_ncc[ref_index] < th_ref_ncc && pyramid_step >= 2)
                        pts[index].flag = 1;
                }
                
                for(long iter = 0 ; iter < max_iter ; iter++)
                {
                    if(savenode[(long)index*(long)max_nodes+(long)iter] < num_triangles && savenode[(long)index*(long)max_nodes+(long)iter] > 0)
                    {
                        long reference_index = 0;
                        long target_index_0 = 0, target_index_1 = 0;
                        const uint32 temp_tri[3] = {tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_X, tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Y, tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Z};
                        
                        bool check_index = false;
                        for(int kk=0;kk<3;kk++)
                        {
                            if(temp_tri[kk] == index)
                                reference_index = index;
                            else
                            {
                                if(!check_index)
                                {
                                    target_index_0 = temp_tri[kk];
                                    check_index    = true;
                                }
                                else
                                    target_index_1 = temp_tri[kk];
                            }
                        }
                        
                        if(reference_index < num_points && target_index_0 < num_points && target_index_1 < num_points &&
                           target_index_0 >= 0 && target_index_1 >= 0)
                        {
                            //plane normal angle
                            const D3DPOINT pt0(pts[reference_index]);
                            const D3DPOINT pt1(pts[target_index_0]);
                            const D3DPOINT pt2(pts[target_index_1]);
                            
                            double angle = SetNormalAngle(pt0, pt1, pt2);
                            
                            const double dh1 = pt0.m_Z - pt1.m_Z;
                            const double dh2 = pt0.m_Z - pt2.m_Z;
                            const double dh3 = pt1.m_Z - pt2.m_Z;
                            
                            if(dh1 >= 0 && dh2 >= 0 && angle > 30)
                                count_th_positive ++;
                            if(dh1 <  0 && dh2 < 0  && angle > 30)
                                count_th_negative ++;
                            
                            bool check_match = false;
                            if(pyramid_step  > 1)
                            {
                                check_match     = true;
                            }
                            else if(pyramid_step == 1)
                            {
                                if(iteration == 1 && count_bl <= 5)
                                    check_match = true;
                            }
                            else if(pyramid_step == 0)
                            {
                                if(iteration == 1 && count_bl <= 5)
                                    check_match = true;
                                else if(iteration == 2 && count_bl <= 5)
                                    check_match = false;
                            }
                            else
                                check_match = false;
                            
                            if(check_match)
                            {
                                const double ddh_1     = fabs(dh1) - fabs(dh2);
                                const double ddh_2     = fabs(dh1) - fabs(dh3);
                                const double ddh_3     = fabs(dh2) - fabs(dh3);
                                
                                if((fabs(ddh_1) > height_th || fabs(ddh_2) > height_th || fabs(ddh_3) > height_th))
                                {
                                    // all node check with height difference.
                                    double h1,h2,dh;
                                    int t_o_min,t_o_max,t_o_mid;
                                    long order[3]    = {reference_index, target_index_0, target_index_1};
                                    double height[3] = {pt0.m_Z, pt1.m_Z, pt2.m_Z};
                                    
                                    long Col         = (long)((pt0.m_X - boundary[0])/gridspace + 0.5);
                                    long Row         = (long)((pt0.m_Y- boundary[1])/gridspace + 0.5);
                                    long Index[3];
                                    Index[0]     = (long)gridsize.width*Row + Col;
                                    Col             = (long)((pt1.m_X - boundary[0])/gridspace + 0.5);
                                    Row             = (long)((pt1.m_Y - boundary[1])/gridspace + 0.5);
                                    Index[1]     = (long)gridsize.width*Row + Col;
                           
                                    Col             = (long)((pt2.m_X - boundary[0])/gridspace + 0.5);
                                    Row             = (long)((pt2.m_Y - boundary[1])/gridspace + 0.5);
                                    Index[2]     = (long)gridsize.width*Row + Col;
                                    
                                    check_neigh = true;
                                    
                                    if(height[0] > height[1])
                                    {
                                        if(height[0] > height[2])
                                        {
                                            t_o_max = 0;
                                            if(height[1] > height[2])
                                            {
                                                t_o_min = 2;
                                                t_o_mid = 1;
                                            }
                                            else
                                            {
                                                t_o_min = 1;
                                                t_o_mid = 2;
                                            }
                                        }
                                        else
                                        {
                                            t_o_max = 2;
                                            t_o_min = 1;
                                            t_o_mid = 0;
                                        }
                                    }
                                    else
                                    {
                                        if(height[1] > height[2])
                                        {
                                            t_o_max = 1;
                                            if(height[0] > height[2])
                                            {
                                                t_o_min = 2;
                                                t_o_mid = 0;
                                            }
                                            else
                                            {
                                                t_o_min = 0;
                                                t_o_mid = 2;
                                            }
                                        }
                                        else
                                        {
                                            t_o_max = 2;
                                            t_o_mid = 1;
                                            t_o_min = 0;
                                        }
                                    }
                                    
                                    h1        = height[t_o_mid] - height[t_o_min];
                                    h2        = height[t_o_max] - height[t_o_mid];
                                    dh        = h1 - h2;
                                    if((dh > 0 && dh > height_th))
                                    {
                                        if(IsRA == 1)
                                        {
                                            pts[order[t_o_min]].flag = 3;
                                        }
                                        else
                                        {
                                            if(flag_blunder)
                                            {
                                                if(ortho_ncc[Index[t_o_min]] < ortho_ncc_th)
                                                    pts[order[t_o_min]].flag = 3;
                                            }
                                            else
                                                if(ortho_ncc[Index[t_o_min]] < ortho_ancc_th)
                                                    pts[order[t_o_min]].flag = 3;
                                            
                                        }
                                    }
                                    else if((dh < 0 && fabs(dh) > height_th))
                                    {
                                        if(IsRA == 1)
                                        {
                                            pts[order[t_o_max]].flag = 3;
                                        }
                                        else
                                        {
                                            if(flag_blunder)
                                            {
                                                if(ortho_ncc[Index[t_o_max]] < ortho_ncc_th)
                                                    pts[order[t_o_max]].flag = 3;
                                            }
                                            else
                                                if(ortho_ncc[Index[t_o_max]] < ortho_ancc_th)
                                                    pts[order[t_o_max]].flag = 3;
                                        }
                                    }
                                }
                            }
                            count++;
                        }
                    }
                }
                
                if(IsRA == 1)
                {
                    if(check_neigh == true)
                    {
                        if(!flag_blunder)
                            pts[index].flag = 1;
                    }
      
                    if((count_th_positive >= (int)(count*0.7 + 0.5) || count_th_negative >= (int)(count*0.7 + 0.5)) && count > 1)
                        pts[index].flag = 1;
                }
                else
                {
                    if(check_neigh == true)
                    {
                        if(!flag_blunder)
                        {
                            if(pyramid_step >= 3 && iteration == 4)
                            {
                                if(ortho_ncc[ref_index] < ortho_ancc_th)
                                    pts[index].flag = 1;
                            }
                            else
                                if(ortho_ncc[ref_index] < ortho_ancc_th)
                                    pts[index].flag = 1;
                        }
                    }
                    
                    double count_threshold;
                    count_threshold = count*0.5;
                    if(pyramid_step >= 2)
                        count_threshold = count*0.5;
                    else if(pyramid_step == 1)
                        count_threshold = (int)(count*0.6 + 0.5);
                    else
                        count_threshold = (int)(count*0.7 + 0.5);
                    
                    if((count_th_positive >= count_threshold || count_th_negative >= count_threshold) && count > 1)
                    {
                        if(flag_blunder)
                        {
                            double tmp_th = 0.6 - (4-pyramid_step)*0.1;
                            if(pyramid_step == 1)
                                tmp_th = 0.6;
                            
                            if(pyramid_step >= 1)
                            {
                                if(ortho_ncc[ref_index] < tmp_th)
                                    pts[index].flag = 1;
                            }
                            else if(pyramid_step == 0)
                            {
                                if(iteration <= 1)
                                {
                                    if(ortho_ncc[ref_index] < 0.9)
                                        pts[index].flag = 1;
                                }
                                else if(iteration == 2 )
                                    pts[index].flag = 1;
                                else {
                                    if(ortho_ncc[ref_index] < ortho_ncc_th)
                                        pts[index].flag = 1;
                                }
                            }
                        }
                        else
                        {
                            if(ortho_ncc[ref_index] < ortho_ancc_th)
                                pts[index].flag = 1;
                        }
                    }
                }
                //If a point was found to be a blunder in the current blunder detection, mark detectedBlunders index as true
                detectedBlunders[index] = (pts[index].flag==1 || pts[index].flag==3);
            }
        }
        
        for(long tcnt=0;tcnt<num_points;tcnt++)
        {
            if(pts[tcnt].flag == 1)
                blunder_count[1]++;
            else if(pts[tcnt].flag == 3)
            {
                blunder_count[1]++;
                pts[tcnt].flag = 1;
            }
            else
            {
                blunder_count[0]++;
                
                if(*minz_mp > pts[tcnt].m_Z)
                    *minz_mp        = pts[tcnt].m_Z;
                if(*maxz_mp < pts[tcnt].m_Z)
                    *maxz_mp        = pts[tcnt].m_Z;
            }
        }
        free(savenode);
        free(nodecount);
    }
  
    return true;
}
*/
int SetttingFlagOfGrid(LevelInfo &rlevelinfo, UGRID *GridPT3, vector<D3DPOINT> MatchedPts_list_anchor,vector<D3DPOINT> MatchedPts_list_blunder,vector<D3DPOINT> *MatchedPts_list_mps)
{
    int total_count = 0;
    double X,Y,Z;
    short roh;
    int t_flag;
    long i = 0;
    long grid_index;
    long t_col, t_row;
    vector<D3DPOINT>::iterator it;
    
    i = 0;
    for( i = 0 ; i < MatchedPts_list_anchor.size() ; i++)
    //for(it = MatchedPts_list_anchor.begin(); it != MatchedPts_list_anchor.end() ; ++it)
    {
        X = MatchedPts_list_anchor[i].m_X;
        Y = MatchedPts_list_anchor[i].m_Y;
        Z = MatchedPts_list_anchor[i].m_Z;
        roh = MatchedPts_list_anchor[i].m_roh;
        t_flag = MatchedPts_list_anchor[i].flag;
        
        t_col         = (int)((X - rlevelinfo.Boundary[0])/(*rlevelinfo.grid_resolution) + 0.5);
        t_row         = (int)((Y - rlevelinfo.Boundary[1])/(*rlevelinfo.grid_resolution) + 0.5);
        grid_index     = (long)rlevelinfo.Size_Grid2D->width*t_row + t_col;
        if(grid_index >= 0 && grid_index < *rlevelinfo.Grid_length && t_col >=0 && t_col < rlevelinfo.Size_Grid2D->width && t_row >= 0 && t_row < rlevelinfo.Size_Grid2D->height)
        {
            GridPT3[grid_index].anchor_flag = 1;
        }
        
   //     i++;
    }
 
    i=0;
    for( i = 0 ; i < MatchedPts_list_blunder.size() ; i++)
    //for(it = MatchedPts_list_blunder.begin(); it != MatchedPts_list_blunder.end() ; ++it)
    {
        X = MatchedPts_list_blunder[i].m_X;
        Y = MatchedPts_list_blunder[i].m_Y;
        Z = MatchedPts_list_blunder[i].m_Z;
        roh = MatchedPts_list_blunder[i].m_roh;
        t_flag = MatchedPts_list_blunder[i].flag;
        
        t_col         = (int)((X - rlevelinfo.Boundary[0])/(*rlevelinfo.grid_resolution) + 0.5);
        t_row         = (int)((Y - rlevelinfo.Boundary[1])/(*rlevelinfo.grid_resolution) + 0.5);
        
        grid_index     = (long)rlevelinfo.Size_Grid2D->width*t_row + t_col;
        if(grid_index >= 0 && grid_index < *rlevelinfo.Grid_length && t_col >=0 && t_col < rlevelinfo.Size_Grid2D->width && t_row >= 0 && t_row < rlevelinfo.Size_Grid2D->height)
        {
            D3DPOINT temp_pts;
            temp_pts.m_X = X;
            temp_pts.m_Y = Y;
            temp_pts.m_Z = Z;
            temp_pts.m_roh = roh;
            temp_pts.flag = t_flag;
            
            MatchedPts_list_mps->push_back(temp_pts);
            
            total_count++;
            if(GridPT3[grid_index].anchor_flag != 1)
            {
                GridPT3[grid_index].anchor_flag = 2;
                
            }
        }
        
        //i++;
    }

    return total_count;
}

UGRID* SetHeightRange(ProInfo *proinfo, LevelInfo &rlevelinfo, const int numOfPts, const int num_triangles, UGRID *GridPT3, const int iteration, double *minH_grid, double *maxH_grid, D3DPOINT *pts, const UI3DPOINT *tris, const double MPP, const bool level_check_matching_rate)
{
    UGRID *result = NULL;
    
    double Total_Min_Z      =  100000;
    double Total_Max_Z      = -100000;
    
    printf("MPP of setheightrange = %f\n",MPP);
    
    const int pyramid_step      = *rlevelinfo.Pyramid_step;
    const double gridspace      = *rlevelinfo.grid_resolution;
    const double *boundary      = rlevelinfo.Boundary;
    
    double BufferOfHeight   = MPP*4.0*pwrtwo(pyramid_step);
    if (pyramid_step == 1)
    {
        if(iteration >= 2)
            BufferOfHeight = MPP*2;
        else
            BufferOfHeight = MPP*3;
    }
    else if(pyramid_step == 0)
    {
        if(iteration == 1)
            BufferOfHeight = MPP*2;
        else
            BufferOfHeight = MPP;
        
        if (BufferOfHeight < 0.5)
            BufferOfHeight = 0.5;
    }
    
    if(proinfo->pre_DEMtif)
    {
        if(BufferOfHeight > proinfo->seedDEMsigma)
            BufferOfHeight = proinfo->seedDEMsigma;
        
        printf("buff %f seed %f \n",BufferOfHeight,proinfo->seedDEMsigma);
    }
    
    BufferOfHeight = ceil(BufferOfHeight);
    
    if(BufferOfHeight > 100)
        BufferOfHeight = 100;
    
    printf("BufferOfHeight = %f\n",BufferOfHeight);
    
    double th_HG = 100000;
    
    if(proinfo->DEM_resolution >= 4)
    {
        if(pyramid_step <= 1)
            th_HG = 1000;
    }
    else
    {
        if(!level_check_matching_rate)
        {
            if(pyramid_step == 2 && iteration >= 5)
                th_HG = 1000;
            else if(pyramid_step == 1)
                th_HG = 500;
            else if(pyramid_step == 0)
                th_HG = 100;
        }
        else
        {
            if(pyramid_step <= 2)
                th_HG = 100;
        }
        
        if(pyramid_step == 0 && proinfo->DEM_resolution < 2 && iteration > 1)
            th_HG = 50;
    }
    
    printf("BufferOfHeight = %f\tth_HG = %f\n",BufferOfHeight,th_HG);
    
    uint8 *m_bHeight       = (uint8*)calloc(*rlevelinfo.Grid_length,sizeof(uint8));
    float *NewHeight = (float*)malloc(sizeof(float)*(*rlevelinfo.Grid_length));
#pragma omp parallel for schedule(static)
    for (long counter = 0; counter < *rlevelinfo.Grid_length; counter++)
        NewHeight[counter]          = -1000.0;
    
    for(long tcnt = 0;tcnt < num_triangles;tcnt++)
    {
        const UI3DPOINT &t_tri = (tris[tcnt]);
        const int pdex0 = t_tri.m_X;
        const int pdex1 = t_tri.m_Y;
        const int pdex2 = t_tri.m_Z;
        
        if(pdex0 < numOfPts && pdex1 < numOfPts && pdex2 < numOfPts)
        {
            const D3DPOINT &TriP1(pts[pdex0]);
            const D3DPOINT &TriP2(pts[pdex1]);
            const D3DPOINT &TriP3(pts[pdex2]);
            int PixelMinXY[2], PixelMaxXY[2];
            double temp_MinZ, temp_MaxZ;
            SetTinBoundary(rlevelinfo, TriP1, TriP2, TriP3, PixelMinXY, PixelMaxXY, Total_Min_Z, Total_Max_Z, temp_MinZ, temp_MaxZ);
            
            double angle = SetNormalAngle(TriP1, TriP2, TriP3);
            
            double angle_weight;
            if(angle <= 30)
                angle_weight = 1.0 + 0.05*(angle)/5.0;
            else
                angle_weight = 1.3;
            
            double BF;
            if(pyramid_step >= 3)
                BF = ceil(BufferOfHeight*angle_weight);
            else
                BF = ceil(BufferOfHeight);
            
            for (long Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
            {
                for (long Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
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
                            if(pyramid_step < 2)
                            {
                                // IDW
                                const double diff1 = SQRT(CurGPXY, TriP1, 2);
                                const double diff2 = SQRT(CurGPXY, TriP2, 2);
                                const double diff3 = SQRT(CurGPXY, TriP3, 2);
                                
                                if(diff1 == 0)
                                {
                                    Z   = TriP1.m_Z;
                                    GridPT3[Index].Matched_flag = 2;
                                }
                                else if(diff2 == 0)
                                {
                                    Z   = TriP2.m_Z;
                                    GridPT3[Index].Matched_flag = 2;
                                }
                                else if(diff3 == 0)
                                {
                                    Z   = TriP3.m_Z;
                                    GridPT3[Index].Matched_flag = 2;
                                }
                            }
                            
                            m_bHeight[Index] = 1;
                            GridPT3[Index].Height = (float)Z;
                            NewHeight[Index] = GridPT3[Index].Height;
                            
                            double ortho_ncc_th;
                            if(pyramid_step == 4 )
                                ortho_ncc_th = 0.7;
                            else if(pyramid_step == 3)
                                ortho_ncc_th = 0.6 - (iteration - 1)*0.02;
                            else if(pyramid_step == 2)
                                ortho_ncc_th = 0.5 - (iteration - 1)*0.02;
                            else if(pyramid_step == 1)
                                ortho_ncc_th = 0.4;
                            else
                                ortho_ncc_th = 0.3;
                            
                            if(ortho_ncc_th < 0.3)
                                ortho_ncc_th = 0.3;
                            
                            //min, max height setting
                            double t1, t2;
                            t1       = min(temp_MinZ, Z);
                            if(GridPT3[Index].Matched_flag == 4) //extension minHeight
                            {
                                if(t1 - BF <= GridPT3[Index].minHeight)
                                    GridPT3[Index].minHeight   = floor(t1 - BF);
                            }
                            else
                                GridPT3[Index].minHeight   = floor(t1 - BF);
                            
                            t2       = max(temp_MaxZ, Z);
                            if(GridPT3[Index].Matched_flag == 4)
                            {
                                if(t2 + BF >= GridPT3[Index].maxHeight)
                                    GridPT3[Index].maxHeight   = ceil(t2 + BF);
                            }
                            else
                                GridPT3[Index].maxHeight   = ceil(t2 + BF);
                            
                            if(GridPT3[Index].minHeight > GridPT3[Index].maxHeight)
                            {
                                int temp_H = GridPT3[Index].minHeight;
                                GridPT3[Index].minHeight = GridPT3[Index].maxHeight;
                                GridPT3[Index].maxHeight = temp_H;
                            }
                            
                            if(GridPT3[Index].Matched_flag == 2)
                            {
                                int mintemp, maxtemp;
                                mintemp = floor(Z - BF);
                                maxtemp = ceil(Z + BF);
                                
                                GridPT3[Index].minHeight = mintemp;
                                GridPT3[Index].maxHeight = maxtemp;
                            }
                            else
                            {
                                GridPT3[Index].Matched_flag = 1;
                                if(GridPT3[Index].Mean_ortho_ncc > DoubleToSignedChar_grid(ortho_ncc_th))
                                {
                                    GridPT3[Index].minHeight = floor(Z - BF);
                                    GridPT3[Index].maxHeight = ceil(Z + BF);
                                }
                                else if(pyramid_step >= 2)
                                {
                                    GridPT3[Index].minHeight = floor(t1 - BF);
                                    GridPT3[Index].maxHeight = ceil(t2 + BF);
                                }
                            }
                            
                            GridPT3[Index].minHeight = floor(GridPT3[Index].minHeight);
                            GridPT3[Index].maxHeight = ceil(GridPT3[Index].maxHeight);
                            
                            int HeightGap = ((GridPT3[Index].maxHeight - GridPT3[Index].minHeight));
                            if(HeightGap > th_HG)
                            {
                                GridPT3[Index].minHeight = -100;
                                GridPT3[Index].maxHeight = -100;
                            }
                        }
                        
                        if(GridPT3[Index].minHeight < Nodata)
                        {
                            printf("minHeight %d\n",GridPT3[Index].minHeight);
                            exit(1);
                        }
                    }
                }
            }
        }
    }
    
    double minH_temp = 100000.0;
    double maxH_temp = Nodata;
    
    UGRID *GridPT3_temp = (UGRID*)malloc(sizeof(UGRID)*(*rlevelinfo.Grid_length));
#pragma omp parallel for schedule(static)
    for (long counter = 0; counter < *rlevelinfo.Grid_length; counter++)
    {
        GridPT3[counter].Height = NewHeight[counter];
        GridPT3_temp[counter] = GridPT3[counter];
    }
    
    free(NewHeight);
    
    // minmaxheight setup for no matched grids
#pragma omp parallel for reduction(max:maxH_temp) reduction(min:minH_temp)
    for (long Row_R=0; Row_R < (long)rlevelinfo.Size_Grid2D ->height; Row_R++)
    {
        for (long Col_C=0; Col_C < (long)rlevelinfo.Size_Grid2D->width; Col_C++)
        {
            long Index = rlevelinfo.Size_Grid2D->width*Row_R + Col_C;
            if(GridPT3[Index].Matched_flag == 0)
            {
                int min_H = 9999;
                int max_H = -9999;
                bool  t_flag_min = false;
                bool  t_flag_max = false;
                for(long t_r = -1 ; t_r <= 1 ; t_r++)
                {
                    for(long t_c = -1 ; t_c <= 1 ; t_c++)
                    {
                        if(Col_C + t_c >= 0 && Col_C + t_c < (rlevelinfo.Size_Grid2D->width) && Row_R + t_r >= 0 && Row_R + t_r < (rlevelinfo.Size_Grid2D->height))
                        {
                            long t_index = (Col_C + t_c) + (Row_R + t_r)*(long)rlevelinfo.Size_Grid2D->width;
                            if(GridPT3_temp[t_index].Matched_flag != 0)
                            {
                                if(min_H > GridPT3_temp[t_index].minHeight)
                                {
                                    min_H   = GridPT3_temp[t_index].minHeight;
                                    t_flag_min = true;
                                }
                                if(max_H < GridPT3_temp[t_index].maxHeight)
                                {
                                    max_H   = GridPT3_temp[t_index].maxHeight;
                                    t_flag_max = true;
                                }
                            }
                        }
                    }
                }
                
                if(t_flag_min || t_flag_max)
                {
                    if(t_flag_min)
                        GridPT3[Index].minHeight = min_H;
                    if(t_flag_max)
                        GridPT3[Index].maxHeight = max_H;
                    m_bHeight[Index] = 2;
                }
                else
                {
                    if(pyramid_step <= 1)
                    {
                        double diff_H = fabs(Total_Max_Z - Total_Min_Z)/3.0;
                        double BF;
                        if(BufferOfHeight > diff_H)
                            BF = diff_H;
                        else
                            BF = BufferOfHeight;
                        
                        if(BF < pwrtwo(pyramid_step)*proinfo->resolution)
                            BF = pwrtwo(pyramid_step)*proinfo->resolution;
                        
                        if(GridPT3[Index].minHeight < Total_Min_Z - BF)
                            GridPT3[Index].minHeight   =  floor(Total_Min_Z - BF);
                        if(GridPT3[Index].maxHeight > Total_Max_Z + BF)
                            GridPT3[Index].maxHeight   =  ceil(Total_Max_Z + BF);
                    }
                }
            }
            else
            {
                if(m_bHeight[Index] == 0)
                    GridPT3[Index].Matched_flag = 0;
            }
            
            if(minH_temp > GridPT3[Index].minHeight)
                minH_temp   = (double)GridPT3[Index].minHeight;
            if(maxH_temp < GridPT3[Index].maxHeight)
                maxH_temp   = (double)GridPT3[Index].maxHeight;
            
            GridPT3[Index].minHeight = floor(GridPT3[Index].minHeight);
            GridPT3[Index].maxHeight = ceil(GridPT3[Index].maxHeight);
            
            int HeightGap = ((GridPT3[Index].maxHeight - GridPT3[Index].minHeight));
            if(HeightGap > th_HG)
            {
                GridPT3[Index].minHeight = -100;
                GridPT3[Index].maxHeight = -100;
            }
        }
    }
    
    *minH_grid = minH_temp;
    *maxH_grid = maxH_temp;
    
    free(GridPT3_temp);
    
    printf("end grid set height!!\t%f\t%f\t%f\t%f\n",*minH_grid,*maxH_grid,Total_Min_Z,Total_Max_Z);
    
    result                  = (UGRID*)calloc(sizeof(UGRID),*rlevelinfo.Grid_length);
    
#pragma omp parallel for
    for(long k=0;k<rlevelinfo.Size_Grid2D->height;k++)
    {
        for(long j=0;j<rlevelinfo.Size_Grid2D->width;j++)
        {
            long matlab_index    = k*(long)rlevelinfo.Size_Grid2D->width + j;
            
            result[matlab_index].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
            
            if(GridPT3[matlab_index].minHeight > GridPT3[matlab_index].maxHeight)
            {
                int temp = GridPT3[matlab_index].maxHeight;
                GridPT3[matlab_index].maxHeight = GridPT3[matlab_index].minHeight;
                GridPT3[matlab_index].minHeight = temp;
            }
            
            result[matlab_index].Matched_flag               = GridPT3[matlab_index].Matched_flag;
            if(m_bHeight[matlab_index] == 2 && GridPT3[matlab_index].Matched_flag == 0)
                result[matlab_index].Matched_flag           = 1;
            
            result[matlab_index].roh                        = GridPT3[matlab_index].roh;
            
            result[matlab_index].Height                     = GridPT3[matlab_index].Height;
            
            result[matlab_index].anchor_flag                = 0;
            
            for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs ; ti++)
                result[matlab_index].ortho_ncc[ti]      = GridPT3[matlab_index].ortho_ncc[ti];
       
            result[matlab_index].Mean_ortho_ncc             = GridPT3[matlab_index].Mean_ortho_ncc;
            result[matlab_index].minHeight                  = GridPT3[matlab_index].minHeight;
            
            if(result[matlab_index].minHeight < -100)
                result[matlab_index].minHeight          = -100;
            
            result[matlab_index].maxHeight                  = GridPT3[matlab_index].maxHeight;
            
            result[matlab_index].selected_pair              = GridPT3[matlab_index].selected_pair;
            result[matlab_index].ncc_seleceted_pair         = GridPT3[matlab_index].ncc_seleceted_pair;
            
            result[matlab_index].total_images               = GridPT3[matlab_index].total_images;
            //if(m_bHeight[matlab_index] == 0)
            //    result[matlab_index].Height = -1000.0;
            
            if(result[matlab_index].minHeight == 0)
                result[matlab_index].minHeight = -1;
            if(result[matlab_index].maxHeight == 0)
                result[matlab_index].maxHeight = 1;
            
            free(GridPT3[matlab_index].ortho_ncc);
        }
    }
    
    if(m_bHeight)
        free(m_bHeight);
    
    printf("end updating grid set height!!\n");
    
    
    free(GridPT3);
    
    printf("end memory release updating grid set height!!\n");
    
    return result;
}

UGRID* ResizeGirdPT3(ProInfo *proinfo, LevelInfo &rlevelinfo, CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight)
{
    if(resize_Size.height > 8000 || resize_Size.width > 8000)
        printf("resize memory allocation start\n");
    
    UGRID *resize_GridPT3 = (UGRID *)calloc(sizeof(UGRID),(long)resize_Size.height*(long)resize_Size.width);
    
    if(resize_Size.height > 8000 || resize_Size.width > 8000)
        printf("resize memory allocation start %ld\n",sizeof(UGRID)*(long)resize_Size.height*(long)resize_Size.width);
    
    printf("preresize memory allocation start %ld\n",sizeof(UGRID)*(long)preSize.height*(long)preSize.width);
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
            if(pos_c >= 0 && pos_c < preSize.width && pos_r >= 0 && pos_r < preSize.height && pre_index >= 0 && pre_index < (long)preSize.width*(long)preSize.height)
            {
                resize_GridPT3[index].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
                
                resize_GridPT3[index].minHeight     = preGridPT3[pre_index].minHeight;
                resize_GridPT3[index].maxHeight     = preGridPT3[pre_index].maxHeight;
                resize_GridPT3[index].Height        = preGridPT3[pre_index].Height;
                resize_GridPT3[index].Matched_flag  = preGridPT3[pre_index].Matched_flag;
                resize_GridPT3[index].roh           = preGridPT3[pre_index].roh;
                resize_GridPT3[index].anchor_flag   = preGridPT3[pre_index].anchor_flag;
                
                resize_GridPT3[index].Mean_ortho_ncc=preGridPT3[pre_index].Mean_ortho_ncc;
                
                resize_GridPT3[index].selected_pair = preGridPT3[pre_index].selected_pair;
                resize_GridPT3[index].ncc_seleceted_pair = preGridPT3[pre_index].ncc_seleceted_pair;
                
                resize_GridPT3[index].total_images  = preGridPT3[pre_index].total_images;
                
                for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs ; ti++)
                {
                    //if(proinfo->check_selected_image[ti])
                        resize_GridPT3[index].ortho_ncc[ti]     = preGridPT3[pre_index].ortho_ncc[ti];
                }
            }
            else
            {
                resize_GridPT3[index].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
                
                resize_GridPT3[index].minHeight     = floor(minmaxheight[0] - 0.5);
                resize_GridPT3[index].maxHeight     = ceil(minmaxheight[1] + 0.5);
                resize_GridPT3[index].Height        = -1000;
                resize_GridPT3[index].Matched_flag  = 0;
                resize_GridPT3[index].roh           = 0.0;
                resize_GridPT3[index].anchor_flag   = 0;
                for(int ti = 0 ; ti <= rlevelinfo.pairinfo->NumberOfPairs ; ti++)
                {
                    //if(proinfo->check_selected_image[ti])
                        resize_GridPT3[index].ortho_ncc[ti]     = 0;
                }
                resize_GridPT3[index].Mean_ortho_ncc= 0;
                resize_GridPT3[index].selected_pair = 100;
                resize_GridPT3[index].ncc_seleceted_pair = -1;
                resize_GridPT3[index].total_images  = 0;
            }
        }
    }
    
    printf("before release preGirdPT3\n");

    for(long row=0;row<preSize.height;row++)
    {
        for(long col=0;col<preSize.width;col++)
        {
            long index = row*(long)preSize.width + col;
            free(preGridPT3[index].ortho_ncc);
        }
    }
    free(preGridPT3);
    
    printf("after release preGirdPT3\n");

    return resize_GridPT3;
}

UGRID* ResizeGirdPT3_RA(const ProInfo *proinfo,LevelInfo &rlevelinfo, const CSize preSize,const CSize resize_Size,const double* preBoundary,const double* Boundary, const D2DPOINT *resize_Grid, UGRID *preGridPT3, const double pre_gridsize,const double* minmaxheight)
{
    
    UGRID *resize_GridPT3 = (UGRID *)calloc(sizeof(UGRID),resize_Size.height*resize_Size.width);
    
    for(long int row=0;row<resize_Size.height;row++)
    {
        for(long int col=0;col<resize_Size.width;col++)
        {
            long int index = row*(long int)resize_Size.width + col;
            double X = resize_Grid[index].m_X;
            double Y = resize_Grid[index].m_Y;
            long int pos_c = (long int)((X - preBoundary[0])/pre_gridsize);
            long int pos_r = (long int)((Y - preBoundary[1])/pre_gridsize);
            long int pre_index = pos_r*(long int)preSize.width + pos_c;
            if(pos_c >= 0 && pos_c < preSize.width && pos_r >= 0 && pos_r < preSize.height && pre_index >= 0 && pre_index < (long)preSize.width*(long)preSize.height)
            {
                resize_GridPT3[index].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
                
                resize_GridPT3[index].minHeight     = preGridPT3[pre_index].minHeight;
                resize_GridPT3[index].maxHeight     = preGridPT3[pre_index].maxHeight;
                resize_GridPT3[index].Height        = preGridPT3[pre_index].Height;
                resize_GridPT3[index].Matched_flag  = preGridPT3[pre_index].Matched_flag;
                resize_GridPT3[index].roh           = preGridPT3[pre_index].roh;
                resize_GridPT3[index].anchor_flag   = preGridPT3[pre_index].anchor_flag;
                for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs ; ti++)
                {
                    //if(proinfo->check_selected_image[ti])
                        resize_GridPT3[index].ortho_ncc[ti]     = preGridPT3[pre_index].ortho_ncc[ti];
                }
                resize_GridPT3[index].Mean_ortho_ncc=preGridPT3[pre_index].Mean_ortho_ncc;
                resize_GridPT3[index].selected_pair = preGridPT3[pre_index].selected_pair;
            }
            else
            {
                resize_GridPT3[index].ortho_ncc = (short*)calloc(sizeof(short),rlevelinfo.pairinfo->NumberOfPairs);
                
                resize_GridPT3[index].minHeight     = floor(minmaxheight[0] - 0.5);
                resize_GridPT3[index].maxHeight     = ceil(minmaxheight[1] + 0.5);
                resize_GridPT3[index].Height        = -1000;
                resize_GridPT3[index].Matched_flag  = 0;
                resize_GridPT3[index].roh           = 0.0;
                resize_GridPT3[index].anchor_flag   = 0;
                for(int ti = 0 ; ti <= rlevelinfo.pairinfo->NumberOfPairs ; ti++)
                {
                    //if(proinfo->check_selected_image[ti])
                        resize_GridPT3[index].ortho_ncc[ti]     = 0;
                }
                resize_GridPT3[index].selected_pair = 100;
            }
        }
    }
    
    printf("before release preGirdPT3\n");
    
    for(long row=0;row<preSize.height;row++)
    {
        for(long col=0;col<preSize.width;col++)
        {
            long index = row*(long)preSize.width + col;
            free(preGridPT3[index].ortho_ncc);
        }
    }
    
    free(preGridPT3);
    
    printf("after release preGirdPT3\n");
    
    return resize_GridPT3;
}

bool SetHeightRange_blunder(LevelInfo &rlevelinfo, const D3DPOINT *pts, const int numPts, UI3DPOINT *tris,const long num_triangles, UGRID *GridPT3)
{
    
#pragma omp parallel for schedule(guided)
    for(long tcnt=0;tcnt<num_triangles;tcnt++)
    {
        double Total_Min_Z      =  100000;
        double Total_Max_Z      = -100000;

        const UI3DPOINT &t_tri = (tris[tcnt]);
        const int pdex0 = t_tri.m_X;
        const int pdex1 = t_tri.m_Y;
        const int pdex2 = t_tri.m_Z;
        
        if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
        {
            const D3DPOINT &TriP1 = (pts[pdex0]);
            const D3DPOINT &TriP2 = (pts[pdex1]);
            const D3DPOINT &TriP3 = (pts[pdex2]);
            int PixelMinXY[2], PixelMaxXY[2];
            double temp_MinZ, temp_MaxZ;
            SetTinBoundary(
                rlevelinfo,
                TriP1,
                TriP2,
                TriP3,
                PixelMinXY,
                PixelMaxXY,
                Total_Min_Z,
                Total_Max_Z,
                temp_MinZ,
                temp_MaxZ);
            
            for (long Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
            {
                for (long Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
                {
                    const int Index= (long)rlevelinfo.Size_Grid2D->width*Row + Col;

                    float Z = -1000.0;
                    bool rtn = false;
                     
                    D3DPOINT CurGPXY(
                        (Col)*(*rlevelinfo.grid_resolution) + rlevelinfo.Boundary[0],
                        (Row)*(*rlevelinfo.grid_resolution) + rlevelinfo.Boundary[1],
                        0,
                        0);

                    rtn = IsTinInside(CurGPXY, TriP1, TriP2, TriP3, Z);
                    
                    if (rtn)
                    {
#pragma omp atomic write
                        GridPT3[Index].Height = (float)Z;
                    }
                }
            }
        }
    }
    return true;
}

void echoprint_Gridinfo(ProInfo *proinfo,NCCresult* roh_height,int row,int col,int level, int iteration, double update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str)
{
    FILE *outfile_h,*outfile_min, *outfile_max,   *outfile_flag, *outMean_ortho, *outMean_ortho_asc;
    CSize temp_S;
    char t_str[500];
    int k,j;
    /*FILE **outfile_roh;
    outfile_roh = (FILE**)malloc(sizeof(FILE*)*proinfo->number_of_images);*/
    
    //sprintf(t_str,"%s/txt/tin_min_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    //outfile_min   = fopen(t_str,"w");
    //sprintf(t_str,"%s/txt/tin_max_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    //outfile_max   = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/tin_h_level_%d_%d_%d_iter_%d_%s.txt",proinfo->save_filepath,row,col,level,iteration,add_str);
    outfile_h   = fopen(t_str,"wb");
    sprintf(t_str,"%s/txt/tin_ortho_ncc_level_%d_%d_%d_iter_%d_%s.txt",proinfo->save_filepath,row,col,level,iteration,add_str);
    outMean_ortho = fopen(t_str,"wb");
    
    
    //sprintf(t_str,"%s/txt/tin_ortho_ncc_level_%d_%d_%d_iter_%d_%s_asc.txt",save_path,row,col,level,iteration,add_str);
    //outMean_ortho_asc = fopen(t_str,"w");
    /*for(int ti = 0 ; ti < proinfo->number_of_images; ti++)
    {
        sprintf(t_str,"%s/txt/tin_ortho_ncc_level_%d_%d_%d_iter_%d_%s_%d.txt",save_path,row,col,level,iteration,add_str,ti);
        outfile_roh[ti] = fopen(t_str,"w");
    }*/
    /*sprintf(t_str,"%s/txt/tin_flag_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
      outfile_flag  = fopen(t_str,"w");*/
    
//  if(outfile_min && outfile_max && outfile_h && outfile_roh && outfile_flag)
    {
        if(update_flag)
        {
            temp_S.height   = Size_Grid2D->height*2;
            temp_S.width    = Size_Grid2D->width*2;
        }
        else
        {
            temp_S.height   = Size_Grid2D->height;
            temp_S.width    = Size_Grid2D->width;
        }
        
        float* temp_height = (float*)malloc(sizeof(float)*temp_S.height*temp_S.width);
        float* temp_ncc = (float*)malloc(sizeof(float)*temp_S.height*temp_S.width);
        
        for(k=0;k<temp_S.height;k++)
        {
            for(j=0;j<temp_S.width;j++)
            {
                int matlab_index    = k*temp_S.width + j;

                //fprintf(outfile_min,"%f\t",GridPT3[matlab_index].minHeight);
                //fprintf(outfile_max,"%f\t",GridPT3[matlab_index].maxHeight - GridPT3[matlab_index].minHeight);
                //if(GridPT3[matlab_index].Matched_flag != 0)
                //if(roh_height[matlab_index].NumOfHeight > 0)
                {
                    //fprintf(outfile_h,"%f\t",GridPT3[matlab_index].Height);
                    //fprintf(outMean_ortho_asc,"%f\t",SignedCharToDouble_grid(GridPT3[matlab_index].ortho_ncc[1]));
                    //printf("%f\t",GridPT3[matlab_index].Height);
                    temp_height[matlab_index] = GridPT3[matlab_index].Height;
                    temp_ncc[matlab_index] = SignedCharToDouble_grid(GridPT3[matlab_index].Mean_ortho_ncc);
                    /*for(int ti = 0 ; ti < proinfo->number_of_images; ti++)
                        fprintf(outfile_roh[ti],"%f\t",GridPT3[matlab_index].ortho_ncc[ti]);*/
                    /*fprintf(outfile_flag,"%d\t",GridPT3[matlab_index].Matched_flag);*/
                    
                }/*
                else
                {
                    fprintf(outfile_h,"-1000.0\t");
                    fprintf(outMean_ortho,"-1.0\t");
                    
                    for(int ti = 0 ; ti < proinfo->number_of_images; ti++)
                        fprintf(outfile_roh[ti],"-1.0\t");
                }
                  */
            }
            //fprintf(outfile_min,"\n");
            //fprintf(outfile_max,"\n");
            //fprintf(outfile_h,"\n");
            //fprintf(outMean_ortho_asc,"\n");
            /*for(int ti = 0 ; ti < proinfo->number_of_images; ti++)
                fprintf(outfile_roh[ti],"\n");*/
            /*fprintf(outfile_flag,"\n");*/
        }

        fwrite(temp_height,sizeof(float),temp_S.height*temp_S.width,outfile_h);
        fwrite(temp_ncc,sizeof(float),temp_S.height*temp_S.width,outMean_ortho);
        free(temp_height);
        free(temp_ncc);
        //fclose(outfile_min);
        //fclose(outfile_max);
        fclose(outfile_h);
        fclose(outMean_ortho);
        //fclose(outMean_ortho_asc);
        /*for(int ti = 0 ; ti < proinfo->number_of_images; ti++)
            fclose(outfile_roh[ti]);*/
        /*fclose(outfile_flag);*/
    }
}

void echoprint_Gridinfo_asc(ProInfo *proinfo,LevelInfo &rlevelinfo, NCCresult* roh_height,int row,int col,int level, int iteration, CSize Size_Grid2D, UGRID *GridPT3)
{
    CSize temp_S;
    char t_str[500];
    int k,j;
    FILE **outfile_roh;
    FILE *outfile_sp, *outfile_tp, * outfile_ttp;
    
    outfile_roh = (FILE**)malloc(sizeof(FILE*)*rlevelinfo.pairinfo->NumberOfPairs);
    sprintf(t_str,"%s/txt/selected_pair_level_%d_%d_%d_iter_%d.txt",proinfo->save_filepath,row,col,level,iteration);
    outfile_sp = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/total_pair_level_%d_%d_%d_iter_%d.txt",proinfo->save_filepath,row,col,level,iteration);
    outfile_tp = fopen(t_str,"w");
    
//    sprintf(t_str,"%s/txt/bin_pair_level_%d_%d_%d_iter_%d.txt",proinfo->save_filepath,row,col,level,iteration);
//    outfile_ttp = fopen(t_str,"w");
    for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs; ti++)
    {
        sprintf(t_str,"%s/txt/tin_ortho_ncc_level_%d_%d_%d_iter_%d_%d_asc.txt",proinfo->save_filepath,row,col,level,iteration,ti);
        outfile_roh[ti] = fopen(t_str,"w");
    }

    temp_S.height   = Size_Grid2D.height;
    temp_S.width    = Size_Grid2D.width;
    
    for(k=0;k<temp_S.height;k++)
    {
        for(j=0;j<temp_S.width;j++)
        {
            int matlab_index    = k*temp_S.width + j;

            fprintf(outfile_sp,"%d\t",GridPT3[matlab_index].ncc_seleceted_pair);
            fprintf(outfile_tp,"%d\t",GridPT3[matlab_index].total_images);
 //           fprintf(outfile_ttp,"%3.2f\t",GridPT3[matlab_index].height_counts);
            
            for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
            {
                fprintf(outfile_roh[pair_number],"%f\t",SignedCharToDouble_grid(GridPT3[matlab_index].ortho_ncc[pair_number]));
            }
        }
        fprintf(outfile_sp,"\n");
        fprintf(outfile_tp,"\n");
 //       fprintf(outfile_ttp,"\n");
        for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs; ti++)
            fprintf(outfile_roh[ti],"\n");
    }

    fclose(outfile_sp);
    fclose(outfile_tp);
//    fclose(outfile_ttp);
    for(int ti = 0 ; ti < rlevelinfo.pairinfo->NumberOfPairs; ti++)
        fclose(outfile_roh[ti]);
    
}

void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str)
{
    int k,j;
    FILE *outfile_min, *outfile_max, *outfile_h, *outfile_roh, *outfile_diff, *outfile_peak, *outINCC, *outGNCC,*outcount;
    CSize temp_S;
    char t_str[500];
    
    sprintf(t_str,"%s/txt/nccresult_roh1_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_min = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/nccresult_roh2_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_max = fopen(t_str,"w");
    
    sprintf(t_str,"%s/txt/nccresult_max_NCC_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_h   = fopen(t_str,"w");
    
    sprintf(t_str,"%s/txt/nccresult_max_NCC_pos_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_roh = fopen(t_str,"w");
    /*
    sprintf(t_str,"%s/txt/nccresult_peaks_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_peak    = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/nccresult_diff_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
    outfile_diff    = fopen(t_str,"w");\*/
    /*sprintf(t_str,"%s/txt/INCC_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
    outINCC = fopen(t_str,"w");
    sprintf(t_str,"%s/txt/GNCC_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
    outGNCC = fopen(t_str,"w");
    
    sprintf(t_str,"%s/txt/rohcount_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
    outcount    = fopen(t_str,"w");
    */
    temp_S.height   = Size_Grid2D->height;
    temp_S.width    = Size_Grid2D->width;
        
    for(k=0;k<temp_S.height;k++)
    {
        for(j=0;j<temp_S.width;j++)
        {
            int matlab_index    = k*temp_S.width + j;
           
        fprintf(outfile_min,"%f\t",SignedCharToDouble_result(nccresult[matlab_index].result0));
            fprintf(outfile_max,"%f\t",SignedCharToDouble_result(nccresult[matlab_index].result1));
            fprintf(outfile_h,"%d\t",nccresult[matlab_index].NumOfHeight); 
            //fprintf(outfile_roh,"%d\t",nccresult[matlab_index].max_WNCC_pos);
            //fprintf(outfile_peak,"%f\t",nccresult[matlab_index].result4);
            
            
            //fprintf(outINCC,"%f\t",nccresult[matlab_index].INCC);
            //fprintf(outGNCC,"%f\t",nccresult[matlab_index].GNCC);
            //fprintf(outcount,"%d\t",nccresult[matlab_index].roh_count);
        }
        fprintf(outfile_min,"\n");
        fprintf(outfile_max,"\n");
        fprintf(outfile_h,"\n");
        //fprintf(outfile_roh,"\n");
        //fprintf(outfile_peak,"\n");
        //fprintf(outfile_diff,"\n");
        
        //fprintf(outINCC,"\n");
        //fprintf(outGNCC,"\n");
        //fprintf(outcount,"\n");
    }

    fclose(outfile_min);
    fclose(outfile_max);
    fclose(outfile_h);
    //fclose(outfile_roh);
    //fclose(outfile_peak);
    //fclose(outfile_diff);
    
    //fclose(outINCC);
    //fclose(outGNCC);
    //fclose(outcount);
}

int AdjustParam(ProInfo *proinfo, LevelInfo &rlevelinfo, int NumofPts, double **ImageAdjust, uint8 total_pyramid, D3DPOINT* ptslists)
{
    
    const int Pyramid_step = *rlevelinfo.Pyramid_step;
    
    double subA[9][6] = {0};
    double TsubA[6][9] = {0};
    double InverseSubA[6][6] = {0};

    Set6by6Matrix(subA,TsubA,InverseSubA);

    D3DPOINT *Coord           = ps2wgs_3D(*rlevelinfo.param,NumofPts,ptslists);

    int iter_count = 0;
    for(int pair_number = 0 ; pair_number < rlevelinfo.pairinfo->NumberOfPairs ; pair_number++)
    {
        const int reference_id = rlevelinfo.pairinfo->pairs[pair_number].m_X;
        const int ti = rlevelinfo.pairinfo->pairs[pair_number].m_Y;
        
        CSize LImagesize(rlevelinfo.py_Sizes[reference_id][Pyramid_step].width, rlevelinfo.py_Sizes[reference_id][Pyramid_step].height);
        const double left_IA[2] = {0.0, 0.0};
        
        //if(proinfo->check_selected_image[ti])
        {
            bool check_stop = false;
            iter_count = 1;
            while(!check_stop && iter_count < 10)
            {
                bool flag_boundary = false;
                int count_pts = 0;
                double sum_weight_X     = 0;
                double sum_weight_Y     = 0;
                double sum_max_roh      = 0;
                double t_sum_weight_X       = 0;
                double t_sum_weight_Y       = 0;
                double t_sum_max_roh        = 0;
                
                const double b_factor             = pwrtwo(total_pyramid-Pyramid_step+1);
                const int Half_template_size   = (int)(*rlevelinfo.Template_size/2.0);
                int patch_size = (2*Half_template_size+1) * (2*Half_template_size+1);

                std::array<double*, 3>* left_patch_vecs_array;
                std::array<double*, 3>* right_patch_vecs_array;

#pragma omp parallel private(t_sum_weight_X,t_sum_weight_Y,t_sum_max_roh) reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
                {
                    
#pragma omp single
                    {
                        // Make patch vectors thread private rather than private to each loop iteration
                        // These are used by postNCC but allocated here for efficiency
                        left_patch_vecs_array = new std::array<double*, 3>[omp_get_num_threads()];
                        right_patch_vecs_array = new std::array<double*, 3>[omp_get_num_threads()];
                        for (int th=0; th<omp_get_num_threads(); th++) {
                            for (int k=0; k<3; k++)
                            {
                                left_patch_vecs_array[th][k] = (double *)malloc(sizeof(double)*patch_size);
                                right_patch_vecs_array[th][k] = (double *)malloc(sizeof(double)*patch_size);
                            }
                        }
                    }
                    
#pragma omp for schedule(guided)
                    for(long i = 0; i<NumofPts ; i++)
                    {
                        double** left_patch_vecs = left_patch_vecs_array[omp_get_thread_num()].data();
                        double** right_patch_vecs = right_patch_vecs_array[omp_get_thread_num()].data();
                        
                        //calculation image coord from object coord by RFM in left and right image
                        D2DPOINT Left_Imagecoord_p   = GetObjectToImageRPC_single(rlevelinfo.RPCs[reference_id],2,left_IA,Coord[i]);
                        D2DPOINT Left_Imagecoord     = OriginalToPyramid_single(Left_Imagecoord_p,rlevelinfo.py_Startpos[reference_id],Pyramid_step);
                        
                        CSize RImagesize(rlevelinfo.py_Sizes[ti][Pyramid_step].width, rlevelinfo.py_Sizes[ti][Pyramid_step].height);
                        
                        D2DPOINT Right_Imagecoord_p  = GetObjectToImageRPC_single(rlevelinfo.RPCs[ti],2,ImageAdjust[pair_number],Coord[i]);
                        D2DPOINT Right_Imagecoord    = OriginalToPyramid_single(Right_Imagecoord_p,rlevelinfo.py_Startpos[ti],Pyramid_step);
                        
                        if(Left_Imagecoord.m_Y  > Half_template_size*b_factor + 10 && Left_Imagecoord.m_X  > Half_template_size*b_factor + 10 && Left_Imagecoord.m_Y  < LImagesize.height - Half_template_size*b_factor - 10 && Left_Imagecoord.m_X  < LImagesize.width - Half_template_size*b_factor - 10 && Right_Imagecoord.m_Y > Half_template_size*b_factor + 10 && Right_Imagecoord.m_X > Half_template_size*b_factor + 10 && Right_Imagecoord.m_Y < RImagesize.height - Half_template_size*b_factor - 10 && Right_Imagecoord.m_X < RImagesize.width - Half_template_size*b_factor - 10)
                        {
                            long index_l = (long)Left_Imagecoord.m_Y*(long)LImagesize.width + (long)Left_Imagecoord.m_X;
                            long index_r = (long)Right_Imagecoord.m_Y*(long)RImagesize.width + (long)Right_Imagecoord.m_X;
                            if( (index_l > 0 && index_l < (long)LImagesize.height*(long)LImagesize.width) && (index_r > 0 && index_r < (long)RImagesize.height*(long)RImagesize.width) )
                            {
                                double ori_diff = rlevelinfo.py_OriImages[reference_id][index_l] - rlevelinfo.py_OriImages[ti][index_r];
                                
                                if(postNCC(rlevelinfo, ori_diff, Left_Imagecoord, Right_Imagecoord, subA, TsubA, InverseSubA, Half_template_size, reference_id, ti, &t_sum_weight_X, &t_sum_weight_Y, &t_sum_max_roh, left_patch_vecs, right_patch_vecs))
                                {
                                    sum_weight_X += t_sum_weight_X;
                                    sum_weight_Y += t_sum_weight_Y;
                                    sum_max_roh  += t_sum_max_roh;
                                    count_pts++;
                                }
                            }
                        }
                    } // end omp for
                    
#pragma omp single
                    {
                        // free thread-private vectors
                        for (int th=0; th<omp_get_num_threads(); th++) {
                            for (int k=0; k<3; k++)
                            {
                                free(left_patch_vecs_array[th][k]);
                                free(right_patch_vecs_array[th][k]);
                            }
                        }
                        delete[] left_patch_vecs_array;
                        delete[] right_patch_vecs_array;
                    }
                    
                } // end omp parallel

                if(count_pts > 10)
                {
                    double shift_X             = sum_weight_X/sum_max_roh*pwrtwo(Pyramid_step);
                    double shift_Y             = sum_weight_Y/sum_max_roh*pwrtwo(Pyramid_step);
                    if(fabs(shift_Y) < 0.1 && fabs(shift_X) < 0.1)
                        check_stop = true;
         
                    printf("pair_number %d\t%d\t%f\t%f\t%f\t%f\n",pair_number, iter_count,shift_X,shift_Y,ImageAdjust[pair_number][1],ImageAdjust[pair_number][0]);
                    
                    shift_X             += ImageAdjust[pair_number][1];
                    shift_Y             += ImageAdjust[pair_number][0];

                    ImageAdjust[pair_number][1]      = shift_X;
                    ImageAdjust[pair_number][0]      = shift_Y;
                }
                else
                    check_stop = true;

                iter_count++;
            }
        }
    }
    
    free(Coord);
    
    
    
    return iter_count;
}


bool postNCC(LevelInfo &rlevelinfo, const double Ori_diff, const D2DPOINT left_pt, const D2DPOINT right_pt, double subA[][6], double TsubA[][9], double InverseSubA[][6], uint8 Half_template_size, const int reference_ID, const int target_ID, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh, double **left_patch_vecs, double **right_patch_vecs)
{
    bool check_pt = false;
 
    CSize leftsize = rlevelinfo.py_Sizes[reference_ID][*rlevelinfo.Pyramid_step];
    CSize rightsize = rlevelinfo.py_Sizes[target_ID][*rlevelinfo.Pyramid_step];
 
    const int half_mask_size = 1;
    double result_rho[9]  = {};
    
    uint8 cell_count = 0;

    for(int j=0;j<9;j++)
        result_rho[j]       = -1.00;

    for(long mask_row = - half_mask_size ; mask_row <= half_mask_size ; mask_row++)
    {
        for(long mask_col = - half_mask_size ; mask_col <= half_mask_size ; mask_col++)
        {
            int Count_N[3] = {0};
            double rot_theta = (double)(Ori_diff*(*rlevelinfo.bin_angle)*PI/180.0);
        
            for(long row = -Half_template_size; row <= Half_template_size ; row++)
            {
                for(long col = -Half_template_size; col <= Half_template_size ; col++)
                {
                    double radius2  = (double)(row*row + col*col);
                    if(radius2 <= (Half_template_size-1)*(Half_template_size-1))
                    {
                        double pos_row_left      = (left_pt.m_Y + row);
                        double pos_col_left      = (left_pt.m_X + col);

                        double temp_col        = (cos(-rot_theta)*col - sin(-rot_theta)*row);
                        double temp_row        = (sin(-rot_theta)*col + cos(-rot_theta)*row);
                        double pos_row_right     = (right_pt.m_Y + temp_row + mask_row);
                        double pos_col_right     = (right_pt.m_X + temp_col + mask_col);

                        if(pos_row_right-3 >= 0 && pos_row_right+3 < rightsize.height && pos_col_right-3 >= 0 && pos_col_right+3 < rightsize.width &&
                           pos_row_left-3 >= 0 && pos_row_left+3 < leftsize.height && pos_col_left-3 >= 0 && pos_col_left+3 < leftsize.width)
                        {
                            //interpolate left_patch
                            double dx = pos_col_left - (int) (pos_col_left);
                            double dy = pos_row_left - (int) (pos_row_left);
                            long position = (long int) (pos_col_left) + (long int) (pos_row_left) * (long)leftsize.width;
                            
                            double left_patch = InterpolatePatch(rlevelinfo.py_Images[reference_ID], position, leftsize, dx, dy);
                            left_patch_vecs[0][Count_N[0]] = left_patch;

                            //interpolate right_patch
                            dx = pos_col_right - (int) (pos_col_right);
                            dy = pos_row_right - (int) (pos_row_right);
                            position = (long int) (pos_col_right) + (long int) (pos_row_right) * (long)rightsize.width;
                            
                            double right_patch = InterpolatePatch(rlevelinfo.py_Images[target_ID], position, rightsize, dx, dy);
                            right_patch_vecs[0][Count_N[0]] = right_patch;
                            
                            //end
                            Count_N[0]++;

                            int size_1        = (int)(Half_template_size/2);
                            if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
                            {
                                if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
                                {
                                    left_patch_vecs[1][Count_N[1]] = left_patch;
                                    right_patch_vecs[1][Count_N[1]] = right_patch;
                                    Count_N[1]++;
                                }
                            }

                            int size_2        = size_1 + (int)((size_1/2.0) + 0.5);
                            if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
                            {
                                if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
                                {
                                    left_patch_vecs[2][Count_N[2]] = left_patch;
                                    right_patch_vecs[2][Count_N[2]] = right_patch;
                                    Count_N[2]++;
                                }
                            }
                        }
                    }
                }  // end col loop
            }  // end row loop

            if(Count_N[0] > 0 && Count_N[1] && Count_N[2])
            {
                double count_rho = 0;
                double temp_rho = 0;
                for (int k=0; k<3; k++)
                {
                    if (Count_N[k] > 0)
                    {
                        double ncc = Correlate(left_patch_vecs[k],right_patch_vecs[k],Count_N[k]);
                        if (ncc != -99)
                        {
                            count_rho++;
                            temp_rho += ncc;
                        }
                    }
                }

                if (count_rho > 0)
                    temp_rho = temp_rho/count_rho;
                else
                    temp_rho = -1;

                long grid_index           = (mask_row+1)*3 + (mask_col+1);
                if(grid_index < 9)
                    result_rho[grid_index] = temp_rho;
                cell_count++;
            }

        }  // end mask_col loop
    }  // end mask_row loop

    double t_weight_X   = 0;
    double t_weight_Y   = 0;
    double t_max_roh    = 0;
    double XX[6]          = {};
    double ATLT[6]        = {};
    if(cell_count == 9)
    {
        double demnum;
        double max_X        = 100;
        double max_Y        = 100;
        double max_roh      = 0;
        bool find_index_1   = false;
        bool find_index_2   = false;
        bool find_index     = false;

        for(int i=0;i<6;i++)
        {
            for(int j=0;j<1;j++)
            {
                double sum = 0.0;
                for(int k=0;k<9;k++)
                    sum += TsubA[i][k]*result_rho[k*1 + j];
                ATLT[i*1 + j] = sum;
            }
        }

        for(int i=0;i<6;i++)
        {
            for(int j=0;j<1;j++)
            {
                double sum = 0.0;
                for(int k=0;k<6;k++)
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
            if (*rlevelinfo.Pyramid_step >= 2)
                find_index  = find_index_1 & find_index_2 & (max_roh > 0.80);
            else
                find_index  = find_index_1 & find_index_2 & (max_roh > 0.90);

            if(find_index)
            {
                t_weight_X += max_X*max_roh;
                t_weight_Y += max_Y*max_roh;
                t_max_roh  += max_roh;

                check_pt = true;
            }
        }
    }

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

bool check_kernel_size(ProInfo *proinfo, const CSize *Subsetsize, const int Template_size, const int pyramid_step)
{
    bool ret = false;
    int count_image = 0;
    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(Subsetsize[ti].height > Template_size/2*pwrtwo(pyramid_step) && Subsetsize[ti].width > Template_size/2*pwrtwo(pyramid_step))
        {
            proinfo->check_selected_image[ti] = true;
            count_image++;
        }
        else
            proinfo->check_selected_image[ti] = false;
    }
    
    if(!proinfo->check_selected_image[0])
    {
        printf("SubsetImage : Subset of reference image is too small to cover the kernel size!!\n");
        ret = false;
    }
    else if(count_image < 2)
    {
        printf("SubsetImage : not enough image taken to matching (less than 2)!!\n");
        ret = false;
    }
    else
        ret = true;
    return ret;
}

bool check_image_boundary(const ProInfo *proinfo,LevelInfo &plevelinfo, const D2DPOINT pos_xy_m,const D2DPOINT pos_xy,const double minH,const double maxH,const int H_template_size, int &selected_images)
{
    bool check = true;
    selected_images = 0;
    const int buff_pixel = 1;
    bool check_reference = true;
    
    bool stop_check = false;
    int pairs = 0;
    while(!stop_check && pairs < plevelinfo.pairinfo->NumberOfPairs)
    //for(int pairs = 0 ; pairs < plevelinfo.pairinfo->NumberOfPairs ; pairs++)
    {
        for(int count = 0 ; count < 2 ; count ++)
        {
            int ti;
            if(count == 0)
                ti = plevelinfo.pairinfo->pairs[pairs].m_X;
            else
                ti = plevelinfo.pairinfo->pairs[pairs].m_Y;
            
            bool bleft_s, bright_s, bleft_e, bright_e;
            D3DPOINT temp_gp;
            D2DPOINT temp;
            D2DPOINT photo;
            
            //min height
            temp_gp.m_Z = minH;
            if(proinfo->sensor_type == SB)
            {
                temp_gp.m_X = pos_xy.m_X;
                temp_gp.m_Y = pos_xy.m_Y;
                if(count == 0)
                {
                    double temp_LIA[2] = {0.0, 0.0};
                    temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,temp_LIA,temp_gp);
                }
                else
                    temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pairs],temp_gp);
            }
            else
            {
                temp_gp.m_X = pos_xy_m.m_X;
                temp_gp.m_Y = pos_xy_m.m_Y;

                photo  = GetPhotoCoordinate_single(temp_gp,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm);
                temp            = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
            }
            temp        = OriginalToPyramid_single(temp,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);

            if(temp.m_X > H_template_size +buff_pixel && temp.m_X < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].width - H_template_size -buff_pixel && temp.m_Y > H_template_size +buff_pixel && temp.m_Y < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].height- H_template_size -buff_pixel)
                bleft_s     = true;
            else
                bleft_s     = false;

            //max height
            temp_gp.m_Z = maxH;
            if(proinfo->sensor_type == SB)
            {
                temp_gp.m_X = pos_xy.m_X;
                temp_gp.m_Y = pos_xy.m_Y;
                if(count == 0)
                {
                    double temp_LIA[2] = {0.0, 0.0};
                    temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,temp_LIA,temp_gp);
                }
                else
                    temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pairs],temp_gp);
            }
            else
            {
                temp_gp.m_X = pos_xy_m.m_X;
                temp_gp.m_Y = pos_xy_m.m_Y;

                photo  = GetPhotoCoordinate_single(temp_gp,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm);
                temp            = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
            }
            temp        = OriginalToPyramid_single(temp,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);

            if(temp.m_X > H_template_size +buff_pixel && temp.m_X < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].width - H_template_size -buff_pixel && temp.m_Y > H_template_size +buff_pixel && temp.m_Y < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].height- H_template_size -buff_pixel)
                bleft_e     = true;
            else
                bleft_e     = false;
            
            if( bleft_s && bleft_e)
                selected_images++;
            //else if(ti == 0)
            //    check_reference = false;
        }
        
        if(selected_images >= 2)
            check   = true;
        else
            check   = false;
        
        stop_check = check;
        pairs++;
    }
    /*
    if(!check_reference)
        check = false;
    else
    {
        if(selected_images >= 2)
            check   = true;
        else
            check   = false;
    }
     */
    return check;

}

bool check_image_boundary_each(const ProInfo *proinfo,LevelInfo &plevelinfo, const D2DPOINT pos_xy_m,const D2DPOINT pos_xy,const double minH,const double maxH,const int H_template_size, const int image_number, const int pair_number, bool check_ref)
{
    bool check = false;
    const int buff_pixel = 1;
    
    //while(!stop_check && pairs < plevelinfo.pairinfo->NumberOfPairs)
    //for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        int ti = image_number;
        bool bleft_s, bright_s, bleft_e, bright_e;
        D3DPOINT temp_gp;
        D2DPOINT temp;
        D2DPOINT photo;
        
        //min height
        temp_gp.m_Z = minH;
        if(proinfo->sensor_type == SB)
        {
            temp_gp.m_X = pos_xy.m_X;
            temp_gp.m_Y = pos_xy.m_Y;
            if(check_ref)
            {
                double temp_LIA[2] = {0.0, 0.0};
                temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,temp_LIA,temp_gp);
            }
            else
                temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pair_number],temp_gp);
        }
        else
        {
            temp_gp.m_X = pos_xy_m.m_X;
            temp_gp.m_Y = pos_xy_m.m_Y;

            photo  = GetPhotoCoordinate_single(temp_gp,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm);
            temp            = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
        }
        temp        = OriginalToPyramid_single(temp,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);

        if(temp.m_X > H_template_size +buff_pixel && temp.m_X < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].width - H_template_size -buff_pixel && temp.m_Y > H_template_size +buff_pixel && temp.m_Y < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].height- H_template_size -buff_pixel)
            bleft_s     = true;
        else
            bleft_s     = false;

        //max height
        temp_gp.m_Z = maxH;
        if(proinfo->sensor_type == SB)
        {
            temp_gp.m_X = pos_xy.m_X;
            temp_gp.m_Y = pos_xy.m_Y;
            
            if(check_ref)
            {
                double temp_LIA[2] = {0.0, 0.0};
                temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,temp_LIA,temp_gp);
            }
            else
                temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[pair_number],temp_gp);
        }
        else
        {
            temp_gp.m_X = pos_xy_m.m_X;
            temp_gp.m_Y = pos_xy_m.m_Y;

            photo  = GetPhotoCoordinate_single(temp_gp,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm);
            temp            = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
        }
        temp        = OriginalToPyramid_single(temp,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);

        if(temp.m_X > H_template_size +buff_pixel && temp.m_X < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].width - H_template_size -buff_pixel && temp.m_Y > H_template_size +buff_pixel && temp.m_Y < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].height- H_template_size -buff_pixel)
            bleft_e     = true;
        else
            bleft_e     = false;
        
        if( bleft_s && bleft_e)
            check = true;
        else
            check = false;
    }
    
    return check;
}

bool check_image_boundary_height(const ProInfo *proinfo,LevelInfo &plevelinfo, const D2DPOINT pos_xy_m,const D2DPOINT pos_xy,const double Height, const int H_template_size, const int pair_number)
{
    bool check = true;
    const int buff_pixel = 1;
    
    //for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        int ti = pair_number;
        D3DPOINT temp_gp;
        D2DPOINT temp;
        D2DPOINT photo;
        
        temp_gp.m_Z = Height;
        if(proinfo->sensor_type == SB)
        {
            temp_gp.m_X = pos_xy.m_X;
            temp_gp.m_Y = pos_xy.m_Y;
            temp        = GetObjectToImageRPC_single(plevelinfo.RPCs[ti],*plevelinfo.NumOfIAparam,plevelinfo.ImageAdjust[ti],temp_gp);
        }
        else
        {
            temp_gp.m_X = pos_xy_m.m_X;
            temp_gp.m_Y = pos_xy_m.m_Y;

            photo  = GetPhotoCoordinate_single(temp_gp,proinfo->frameinfo.Photoinfo[ti],proinfo->frameinfo.m_Camera, proinfo->frameinfo.Photoinfo[ti].m_Rm);
            temp            = PhotoToImage_single(photo, proinfo->frameinfo.m_Camera.m_CCDSize, proinfo->frameinfo.m_Camera.m_ImageSize);
        }
        temp        = OriginalToPyramid_single(temp,plevelinfo.py_Startpos[ti],*plevelinfo.Pyramid_step);

        if(temp.m_X > H_template_size +buff_pixel && temp.m_X < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].width - H_template_size -buff_pixel && temp.m_Y > H_template_size +buff_pixel && temp.m_Y < plevelinfo.py_Sizes[ti][*plevelinfo.Pyramid_step].height- H_template_size -buff_pixel)
            check     = true;
        else
            check     = false;
    }
    
    return check;
}

double CalMemorySize_Post(CSize DEM_size, CSize Final_DEMsize)
{
    long int Memory = 0;
    
    long int DEM_values = (long)(sizeof(float)*DEM_size.width*DEM_size.height);
    
    printf("DEM_value %f\t%d\t%d\n",DEM_values/1024.0/1024.0/1024.0,DEM_size.width,DEM_size.height);
    
    long int value = (long)(sizeof(float)*(long)Final_DEMsize.width*(long)Final_DEMsize.height);
    printf("value %f\t%d\t%d\n",value/1024.0/1024.0/1024.0,Final_DEMsize.width,Final_DEMsize.height);
    
    long int value_pt = (long)(sizeof(unsigned char)*(long)Final_DEMsize.width*(long)Final_DEMsize.height);
    printf("value_pt %f\n",value_pt/1024.0/1024.0/1024.0);
    
    Memory = DEM_values + value + value_pt;
    
    double result = (double)(Memory/1024.0/1024.0/1024.0);
    printf("DEM total memory %f\n",result);
    
    return result;
}

double CalMemorySize_Post_MT(CSize DEM_size, CSize Final_DEMsize)
{
    long int Memory = 0;
    
    long int DEM_values = (long)(sizeof(float)*DEM_size.width*DEM_size.height);
    
    printf("DEM_value %f\t%d\t%d\n",DEM_values/1024.0/1024.0/1024.0,DEM_size.width,DEM_size.height);
    
    long int value = (long)(sizeof(signed char)*(long)Final_DEMsize.width*(long)Final_DEMsize.height)*2;
    printf("value %f\t%d\t%d\n",value/1024.0/1024.0/1024.0,Final_DEMsize.width,Final_DEMsize.height);
    
    long int value_pt = (long)(sizeof(unsigned char)*(long)Final_DEMsize.width*(long)Final_DEMsize.height)*2;
    printf("value_pt %f\n",value_pt/1024.0/1024.0/1024.0);
    
    Memory = DEM_values + value + value_pt;
    
    double result = (double)(Memory/1024.0/1024.0/1024.0);
    printf("MT total memory %f\n",result);
    
    return result;
}

double CalMemorySize_Post_LSF(CSize DEM_size, CSize Final_DEMsize)
{
    long int Memory = 0;
    
    long int DEM_values = (long)(sizeof(float)*DEM_size.width*DEM_size.height);
    
    printf("DEM_value %f\t%d\t%d\n",DEM_values/1024.0/1024.0/1024.0,DEM_size.width,DEM_size.height);
    
    long int value = (long)(sizeof(float)*(long)Final_DEMsize.width*(long)Final_DEMsize.height);
    printf("value %f\t%d\t%d\n",value/1024.0/1024.0/1024.0,Final_DEMsize.width,Final_DEMsize.height);
    
    long int value_pt = (long)(sizeof(LSFINFO)*(long)Final_DEMsize.width*(long)Final_DEMsize.height);
    printf("value_pt %f\n",value_pt/1024.0/1024.0/1024.0);
    
    Memory = DEM_values + value + value_pt;
    
    double result = (double)(Memory/1024.0/1024.0/1024.0) + 5.0;
    printf("LSF total memory %f\n",result);
    
    return result;
}

void MergeTiles(const ProInfo *info, const int iter_row_start, const int t_col_start, const int iter_row_end, const int t_col_end, int buffer, const int final_iteration, float *DEM, const CSize Final_DEMsize, double *FinalDEM_boundary)
{
    const int find_level = 0;
    const double grid_size = info->DEM_resolution;

    printf("MergeTile boundary %f\t%f\t%f\t%f\n",FinalDEM_boundary[0],FinalDEM_boundary[1],FinalDEM_boundary[2],FinalDEM_boundary[3]);
    
    buffer  = floor(buffer/grid_size);

    printf("dem size %d\t%d\t%d\n",Final_DEMsize.width,Final_DEMsize.height,buffer);
    
    long DEM_data_size = (long)Final_DEMsize.height*(long)Final_DEMsize.width;
#pragma omp parallel for schedule(guided)
    for(long index = 0 ; index < DEM_data_size ; index++)
        DEM[index] = Nodata;
    
    //setting DEM value
    for(int row = iter_row_start ; row <= iter_row_end ; row ++)
    {
        for(int col = t_col_start ; col <= t_col_end ; col++)
        {
            char t_str[500];
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info->save_filepath,row,col,find_level,final_iteration);
            FILE *pfile   = fopen(t_str,"r");
            if(pfile)
            {
                printf("matched tiles %s\n",t_str);
                fseek(pfile,0,SEEK_END);
                long int size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    char h_t_str[500];
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info->save_filepath,row,col);
                    FILE *p_hfile     = fopen(h_t_str,"r");
                    if(p_hfile)
                    {
                        long row_size,col_size;
                        double t_boundary[4];
                        while(!feof(p_hfile))
                        {
                            int t_row,t_col,t_level;
                            double t_grid_size;
                            
                            fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%ld\t%ld\n",
                                   &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
                        }
                        
                        printf("header %f\t%f\t%d\t%d\t%f\t%d\n",t_boundary[0],t_boundary[1],col_size,row_size,grid_size,buffer);
                        
                        char hv_t_str[500];
                        sprintf(hv_t_str,"%s/txt/tin_h_level_%d_%d_%d_iter_%d_final.txt",info->save_filepath,row,col,find_level,final_iteration);
                        FILE *p_hvfile    = fopen(hv_t_str,"rb");
                        
                        if(p_hvfile)
                        {
                            float* temp_height = (float*)malloc(sizeof(float)*col_size*row_size);
                            fread(temp_height,sizeof(float),col_size*row_size,p_hvfile);
                            
                            #pragma omp parallel for schedule(guided)
                            for(long iter_row = 0 ; iter_row < row_size ; iter_row ++)
                            {
                                for(long iter_col = 0 ; iter_col < col_size ; iter_col++)
                                {
                                    long t_col = (long)( (t_boundary[0] + grid_size*iter_col - FinalDEM_boundary[0])  /grid_size);
                                    long t_row = (long)( (FinalDEM_boundary[3] - (t_boundary[1] + grid_size*iter_row))/grid_size);
                                    long index = t_row*(long)Final_DEMsize.width + t_col;
                                    
                                    if(t_col >= 0 && t_col < Final_DEMsize.width && t_row >= 0 && t_row < Final_DEMsize.height &&
                                       iter_row > buffer && iter_row < row_size - buffer &&
                                       iter_col > buffer && iter_col < col_size - buffer)
                                    {
                                        float DEM_value = temp_height[iter_row*col_size + iter_col];
                                        if(DEM_value > -1000)
                                            DEM[index] = DEM_value;
                                    }
                                }
                            }
                            free(temp_height);
                            fclose(p_hvfile);
                        }
                        fclose(p_hfile);
                    }
                }
                fclose(pfile);
            }
        }
    }
    
    char DEM_str[500];
    sprintf(DEM_str, "%s/%s_dem_header_tin.txt", info->save_filepath, info->Outputpath_name);
    
    printf("name %s\t%f\t%f\t%f\t%d\t%d\n",DEM_str,FinalDEM_boundary[0],FinalDEM_boundary[3],grid_size,Final_DEMsize.width,Final_DEMsize.height);
    
    FILE *poutheader = fopen(DEM_str,"w");
    fprintf(poutheader,"%f\t%f\t%f\t%d\t%d\n",FinalDEM_boundary[0],FinalDEM_boundary[3],grid_size,Final_DEMsize.width,Final_DEMsize.height);
    fclose(poutheader);
}

CSize DEM_final_Size(const char *save_path, const int row_start, const int col_start,const int row_end, const int col_end, const double grid_resolution, double *boundary)
{
    CSize final_size;
    
    double minX, maxX, minY, maxY, minHeight, maxHeight;
    minX = 10000000;
    minY = 10000000;
    maxX = -10000000;
    maxY = -10000000;
    
    minHeight = 100000;
    maxHeight = -100000;
    
    printf("DEM grid = %lf\n",grid_resolution);
    int count_matched_files = 0;
    
    for(int index_file = 0 ; index_file <= row_end*col_end ; index_file++)
    {
        const int row = (int)(floor(index_file/col_end)) + 1;
        const int col = index_file%col_end+1;
        
        if(row >= row_start && row <= row_end && col >= col_start &&  col <= col_end)
        {
            char t_str[500];
            double minmaxBR[6];
            sprintf(t_str,"%s/txt/matched_BR_row_%d_col_%d.txt",save_path,row,col);
            FILE *pfile = fopen(t_str,"r");
            if(pfile)
            {
                sprintf(t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",save_path,row,col);
                FILE *p_hfile     = fopen(t_str,"r");
                if(p_hfile)
                {
                    int row_size,col_size;
                    double t_boundary[4] = {0.0};
                    double t_grid_size;
                    
                    while(!feof(p_hfile))
                    {
                        int t_row,t_col,t_level;
                        fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n",
                        &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
                        
                        printf("header file %s\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\n",t_str,t_row,t_col,t_level,t_boundary[0],t_boundary[1],t_grid_size,col_size,row_size);
                    }
                    
                    t_boundary[2] = t_boundary[0] + col_size*t_grid_size;
                    t_boundary[3] = t_boundary[1] + row_size*t_grid_size;
                    
                    count_matched_files++;
                    fscanf(pfile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&minmaxBR[0],&minmaxBR[1],&minmaxBR[2],&minmaxBR[3],&minmaxBR[4],&minmaxBR[5]);
             
                    printf("BR file %s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t_str,minmaxBR[0],minmaxBR[1],minmaxBR[2],minmaxBR[3],minmaxBR[4],minmaxBR[5]);
                    printf("t_bounary %f\t%f\t%f\t%f\n",minmaxBR[0] - t_boundary[0],t_boundary[2] - minmaxBR[2],minmaxBR[1] - t_boundary[1],t_boundary[3] - minmaxBR[3]);
                    if(minmaxBR[0] >= t_boundary[0] - t_grid_size  && minmaxBR[2] <= t_boundary[2] + t_grid_size && minmaxBR[1] >= t_boundary[1] - t_grid_size && minmaxBR[3] <= t_boundary[3] + t_grid_size)
                    {
                        if(minX > minmaxBR[0])
                            minX     = minmaxBR[0];
                        if(minY > minmaxBR[1])
                            minY     = minmaxBR[1];
                        
                        if(maxX < minmaxBR[2])
                            maxX     = minmaxBR[2];
                        if(maxY < minmaxBR[3])
                            maxY     = minmaxBR[3];
                        
                        if(minHeight > minmaxBR[4])
                            minHeight = minmaxBR[4];
                        if(maxHeight < minmaxBR[5])
                            maxHeight = minmaxBR[5];
                        
                        printf("minmax %f\t%f\t%f\t%f\n",minX,maxX,minY,maxY);
                    }
                    fclose(p_hfile);
                }
                fclose(pfile);
            }
        }
    }
    
    if(count_matched_files < 1)
    {
        printf("No matched tiles. Please check overlapped area or image quality!!\n");
        exit(1);
    }
    
    boundary[0] = minX - 10*grid_resolution;
    boundary[1] = minY - 10*grid_resolution;
    boundary[2] = maxX + 10*grid_resolution;
    boundary[3] = maxY + 10*grid_resolution;
    
    final_size.width = (int)((maxX - minX)/grid_resolution) + 1;
    final_size.height = (int)((maxY - minY)/grid_resolution) + 1;
    
    printf("DEM_final_Size %d\t%d\n",final_size.width,final_size.height);
    printf("BR %f\t%f\t%f\t%f\t%f\t%f\n",boundary[0],boundary[1],boundary[2],boundary[3],minHeight,maxHeight);
    
    return final_size;
}

void NNA_M(const ProInfo *proinfo, const TransParam _param, const int row_start, const int col_start, const int row_end, const int col_end, int buffer_clip, const int final_iteration, const int divide, const CSize Final_DEMsize, float* DEM_values, float* value, unsigned char* value_pt, const double *FinalDEM_boundary)
{
    time_t total_ST = 0, total_ET = 0;
    double total_gap;

    printf("DEM grid = %lf\n",proinfo->DEM_resolution);
    
    double minX,maxX, minY, maxY,DEM_minX,DEM_maxY;
    minX = FinalDEM_boundary[0];
    minY = FinalDEM_boundary[1];
    maxX = FinalDEM_boundary[2];
    maxY = FinalDEM_boundary[3];
  
    long DEM_rows, DEM_cols;
    long col_count = Final_DEMsize.width;
    long row_count = Final_DEMsize.height;
    
    char DEM_header[500];
    sprintf(DEM_header, "%s/%s_dem_header_tin.txt", proinfo->save_filepath,proinfo->Outputpath_name);
    FILE *fheader = fopen(DEM_header,"r");
    double dummy;
    fscanf(fheader,"%lf\t%lf\t%lf\t%ld\t%ld",&DEM_minX,&DEM_maxY,&dummy,&DEM_cols,&DEM_rows);
    printf("%f\t%f\t%d\t%d\n",DEM_minX,DEM_maxY,DEM_cols,DEM_rows);
    fclose(fheader);
    
    remove(DEM_header);
    
    total_ST = time(0);
    
    for(long i=0;i<row_count;i++)
    {
        for(long j=0;j<col_count;j++)
            value[i*col_count + j] = Nodata;
    }
    
    long total_search_count = 0;
    long DEM_data_size = DEM_cols*DEM_rows;
    #pragma omp parallel for schedule(guided) reduction(+:total_search_count)
    for(long ix=0;ix<DEM_data_size;ix++)
    {
        long row = floor(ix/DEM_cols);
        long col = ix%DEM_cols;
        
        double t_x = DEM_minX + col*proinfo->DEM_resolution;
        double t_y = DEM_maxY - row*proinfo->DEM_resolution;
        double t_z = DEM_values[row*DEM_cols + col];
        
        if(t_x >= minX && t_x < maxX && t_y >= minY && t_y < maxY)
        {
            long d_row = (long)((maxY - t_y)/proinfo->DEM_resolution);
            long d_col = (long)((t_x - minX)/proinfo->DEM_resolution);
            
            long t_index = d_row*col_count + d_col;
            
            if(t_index < col_count*row_count && d_row < row_count && d_col < col_count)
            {
                value[t_index] = t_z;
                value_pt[t_index] = 0;
                if(t_z > -1000)
                {
                    value_pt[t_index] = 2;
                    total_search_count++;
                }
            }
        }
    }
    free(DEM_values);
    
    printf("Done tfile write\n");
   
    buffer_clip  = floor(buffer_clip/proinfo->DEM_resolution);
    
    for(int index_file = 0 ; index_file <= row_end*col_end ; index_file++)
    {
        long row = floor(index_file/col_end) + 1;
        long col = index_file%col_end + 1;
        
        if(row >= row_start && row <= row_end && col >= col_start &&  col <= col_end)
        {
            char t_str[500];
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_0_%d.txt",proinfo->save_filepath,row,col,final_iteration);
            FILE *pfile   = fopen(t_str,"rb");
            if(pfile)
            {
                long int size;
                fseek(pfile,0,SEEK_END);
                size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    char h_t_str[500];
                    char c_t_str[500];
                    
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    FILE *p_hfile     = fopen(h_t_str,"r");
                    
                    sprintf(c_t_str,"%s/txt/count_row_%d_col_%d.txt",proinfo->save_filepath,row,col);
                    FILE *c_hfile     = fopen(c_t_str,"r");
                    
                    if(p_hfile)
                    {
                        char hv_t_str[500];
                        long row_size, col_size;
                        double t_boundary[4];
                        while(!feof(p_hfile))
                        {
                            int t_row,t_col,t_level;
                            double t_grid_size;
                            
                            fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%ld\t%ld\n",
                                   &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
                        }
                        
                        long count_MPs;
                        fscanf(c_hfile,"%ld\n",&count_MPs);
                        D3DPOINTSAVE *temp_pts = (D3DPOINTSAVE*)malloc(sizeof(D3DPOINTSAVE)*count_MPs);
                        fread(temp_pts,sizeof(D3DPOINTSAVE),count_MPs,pfile);
                        printf("read count_MPs %d\t%d\n",count_MPs,buffer_clip);
                        
                        long count_read = 0;
                        long count_out = 0;
#pragma omp parallel for schedule(guided) reduction(+:count_read,count_out)
                        for(long int t_i = 0 ; t_i < count_MPs ; t_i++)
                        {
                            //if(temp_pts[t_i].flag != 1)
                            {
                                long pos_col = (long)((temp_pts[t_i].m_X - minX)/proinfo->DEM_resolution);
                                long pos_row = (long)((maxY - temp_pts[t_i].m_Y)/proinfo->DEM_resolution);
                                
                                long clip_pos_col = (long)((temp_pts[t_i].m_X - t_boundary[0])/proinfo->DEM_resolution);
                                long clip_pos_row = (long)((temp_pts[t_i].m_Y - t_boundary[1])/proinfo->DEM_resolution);
                                
                                if(pos_row >= 0 && pos_row < row_count && pos_col >= 0 && pos_col < col_count &&
                                   clip_pos_col > buffer_clip && clip_pos_col < col_size - buffer_clip &&
                                   clip_pos_row > buffer_clip && clip_pos_row < row_size - buffer_clip)
                                {
                                    long t_index = pos_row*col_count + pos_col;
                                    value[t_index] = temp_pts[t_i].m_Z;
                                    value_pt[t_index] = 1;
                                    
                                    count_read++;
                                }
                            }
                        }
                        printf("read_done count_MPs %ld\t%ld\n",count_read,count_out);
                        
                        free(temp_pts);
                        fclose(p_hfile);
                    }
                    else
                    {
                        printf("No header file exist : %s\n",h_t_str);
                        printf("Removed %s\n",t_str);
                        printf("Please reprocess!!\n");
                        remove(t_str);
                        exit(1);
                    }
                }
                fclose(pfile);
            }
        }
    }
    
    printf("Done read matched pts\n");
    total_ET = time(0);
    total_gap = difftime(total_ET,total_ST);
    printf("time[m] = %5.2f\n",total_gap/60.0);
    
    //IDW Interpolation
    total_ST = time(0);
    long total_interpolated = 0;
    const int row_interval    = 50;
    if(total_search_count > 0 && !proinfo->check_Matchtag)
    {
#pragma omp parallel for schedule(guided) reduction(+:total_interpolated)
        for(long count = 0;count < row_count*col_count;count++)
        {
            long pos_row = floor(count/col_count);
            long pos_col = count%col_count;
            
            int check = 1;
            
            double query[2] = {pos_col*proinfo->DEM_resolution + minX, maxY - pos_row*proinfo->DEM_resolution};
            long pos_index = pos_row*col_count + pos_col;
            
            if (pos_col >= 0 && pos_col < col_count && pos_row >= 0 && pos_row < row_count)
            {
                if(value_pt[pos_index] == 2 && value[pos_index] > - 1000)
                {
                    check = 0;
                    value_pt[pos_index] = 0;
                }
            }
            
            if(!check)
            {
                //IDW Interpolation
                double height = FindNebPts_F_M_IDW(value, value_pt, row_count, col_count, proinfo->DEM_resolution, minX, maxY, query[0], query[1],row_interval);
                
                if(height > Nodata)
                {
                    value[pos_index] = height;
                    total_interpolated++;
                }
            }
            
            if(value_pt[pos_index] == 2)
                value_pt[pos_index] = 0;
        }
    }

    if(proinfo->check_Matchtag)
    {
    #pragma omp parallel for schedule(guided)
        for(long count = 0;count < row_count*col_count;count++)
        {
            if(value_pt[count] == 2)
                value_pt[count] = 0;
        }
    }
    
    printf("end interpolation\t%ld\t%ld\n",total_search_count,total_interpolated);
    total_ET = time(0);
    total_gap = difftime(total_ET,total_ST);
    printf("time[m] = %5.2f\n",total_gap/60.0);
    
    total_ST = time(0);
    
    //smoothing
    float *value_sm = (float*)malloc(sizeof(float)*col_count*row_count);
    
#pragma omp parallel for schedule(guided)
    for (long index = 0; index < col_count*row_count; index++)
    {
        long row = floor(index/col_count);
        long col = index%col_count;
        long count_cell = 0;
        long null_count_cell = 0;
        double sum_h = 0;
        double null_sum_h = 0;
        
        for (long t_i = -1; t_i <= 1;t_i++ )
        {
            for (long t_j = -1; t_j <= 1; t_j++)
            {
                long index_row = row + t_i;
                long index_col = col + t_j;
                long int t_index = index_row*col_count + index_col;
                
                if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
                {
                    if(value[t_index] != Nodata)
                    {
                        count_cell++;
                        sum_h += value[t_index];
                    }
                    
                    if (value[row*col_count + col] == Nodata)
                    {
                        if(value[t_index] != Nodata && t_i != 0 && t_j != 0)
                        {
                            null_count_cell++;
                            null_sum_h += value[index_row*(long)col_count + index_col];
                        }
                    }
                }
            }
        }
        
        if(count_cell > 0)
            value_sm[row*col_count + col] = sum_h / count_cell;
        else
            value_sm[row*col_count + col] = value[row*col_count + col];
        
        if(null_count_cell > 6)
            value_sm[row*col_count + col] = null_sum_h / null_count_cell;
    }
    printf("end smoothing\n");
    total_ET = time(0);
    total_gap = difftime(total_ET,total_ST);
    printf("time[m] = %5.2f\n",total_gap/60.0);
    
    if(!proinfo->check_Matchtag)
    {
        if(divide == 0)
        {
            char GEOTIFF_dem_filename[500];
            sprintf(GEOTIFF_dem_filename, "%s/%s_dem.tif", proinfo->save_filepath, proinfo->Outputpath_name);
            WriteGeotiff(GEOTIFF_dem_filename, value_sm, col_count, row_count, proinfo->DEM_resolution, minX, maxY, _param.projection, _param.utm_zone, _param.bHemisphere, 4);
        }
        else
        {
            char GEOTIFF_dem_filename[500];
            sprintf(GEOTIFF_dem_filename, "%s/%s_%d_dem.tif", proinfo->save_filepath, proinfo->Outputpath_name,divide);
            WriteGeotiff(GEOTIFF_dem_filename, value_sm, col_count, row_count, proinfo->DEM_resolution, minX, maxY, _param.projection, _param.utm_zone, _param.bHemisphere, 4);
        }
    }
    printf("Done writing DEM tif\n");
    free(value_sm);
}

void MergeTiles_Ortho(const ProInfo *info, const int iter_row_start, const int t_col_start, const int iter_row_end,const int t_col_end, int buffer,const int final_iteration, signed char *DEM_ortho, const CSize Final_DEMsize, const double *FinalDEM_boundary)
{
    const int find_level = 0;
    double grid_size = info->DEM_resolution;
    
    printf("MergeOrtho boundary %f\t%f\t%f\t%f\n",FinalDEM_boundary[0],FinalDEM_boundary[1],FinalDEM_boundary[2],FinalDEM_boundary[3]);
    
    buffer  = floor(buffer/grid_size);
    
    printf("dem size %d\t%d\t%d\t%f\n",Final_DEMsize.width,Final_DEMsize.height,buffer,grid_size);
    
    long DEM_data_size = (long)Final_DEMsize.height*(long)Final_DEMsize.width;
#pragma omp parallel for schedule(guided)
    for(long index_a = 0 ; index_a < DEM_data_size ; index_a++)
        DEM_ortho[index_a] = FloatToSignedChar(-1.0);
    
    //setting Ortho NCC value
    for(long row = iter_row_start ; row <= iter_row_end ; row ++)
    {
        for(long col = t_col_start ; col <= t_col_end ; col++)
        {
            char t_str[500];
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info->save_filepath,row,col,find_level,final_iteration);
            FILE *pfile   = fopen(t_str,"r");
            if(pfile)
            {
                fseek(pfile,0,SEEK_END);
                long int size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    char h_t_str[500];
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info->save_filepath,row,col);
                    FILE *p_hfile     = fopen(h_t_str,"r");
                    if(p_hfile)
                    {
                        long row_size,col_size;
                        double t_boundary[4];
                        while(!feof(p_hfile))
                        {
                            int t_row,t_col,t_level;
                            double t_grid_size;
                            
                            fscanf(p_hfile,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%ld\t%ld\n",
                                   &t_row,&t_col,&t_level,&t_boundary[0],&t_boundary[1],&t_grid_size,&col_size,&row_size);
                        }
                        
                        char ortho_str[500];
                        sprintf(ortho_str,"%s/txt/tin_ortho_ncc_level_%d_%d_%d_iter_%d_final.txt",info->save_filepath,row,col,find_level,final_iteration);
                        FILE *p_orthofile = fopen(ortho_str,"rb");
                        
                        if(p_orthofile)
                        {
                            float* temp_ncc = (float*)malloc(sizeof(float)*col_size*row_size);
                            fread(temp_ncc,sizeof(float),col_size*row_size,p_orthofile);
                            
                            #pragma omp parallel for schedule(guided)
                            for(long iter_row = 0 ; iter_row < row_size ; iter_row ++)
                            {
                                for(long iter_col = 0 ; iter_col < col_size ; iter_col++)
                                {
                                    long t_col = (long)((t_boundary[0] + grid_size*iter_col - FinalDEM_boundary[0])  /grid_size);
                                    long t_row = (long)((FinalDEM_boundary[3] - (t_boundary[1] + grid_size*iter_row))/grid_size);
                                    long index = t_row*(long)Final_DEMsize.width + t_col;
                                    
                                    if(t_col >= 0 && t_col < Final_DEMsize.width && t_row >= 0 && t_row < Final_DEMsize.height &&
                                       iter_row > buffer && iter_row < row_size - buffer &&
                                       iter_col > buffer && iter_col < col_size - buffer)
                                    {
                                        float ortho_value = temp_ncc[iter_row*col_size + iter_col];
                                        if(ortho_value > -1.0)
                                            DEM_ortho[index] = FloatToSignedChar(ortho_value);
                                    }
                                }
                            }
                            free(temp_ncc);
                            fclose(p_orthofile);
                        }
                        
                        fclose(p_hfile);
                    }
                }
                fclose(pfile);
            }
        }
    }
}

void NNA_M_MT(const ProInfo *proinfo, const TransParam _param, const int row_start, const int col_start,const int row_end, const int col_end, int buffer_clip, const int final_iteration, const int divide, signed char* Ortho_values, float* value, unsigned char* value_pt, const CSize Final_DEMsize, const double *FinalDEM_boundary)
{
    double minX = FinalDEM_boundary[0];
    double minY = FinalDEM_boundary[1];
    double maxX = FinalDEM_boundary[2];
    double maxY = FinalDEM_boundary[3];
    const int max_mt_iteration = 10;
    long col_count = (long)Final_DEMsize.width;
    long row_count = (long)Final_DEMsize.height;
    long data_size = row_count*col_count;
    printf("start null\n");
    long count_null_cell = 0;
    int check_while = 0;
    
    signed char *t_value_orthoncc = (signed char*)malloc(sizeof(signed char)*data_size);
    unsigned char *t_value_pt = (unsigned char*)malloc(sizeof(unsigned char)*data_size);
    
    memcpy(t_value_orthoncc,Ortho_values,sizeof(signed char)*data_size);
    memcpy(t_value_pt,value_pt,sizeof(unsigned char)*data_size);
    
    int iter_check = 0;
    int check_size = 0;
    int max_total_iteration = 3;
    if(proinfo->DEM_resolution > 2)
        max_total_iteration = 3;
    int total_iteration;
    const double th_rate = 0.6;
    int grid_rate = (int)(8.0/(double)proinfo->DEM_resolution);
    if(grid_rate < 1)
        grid_rate = 1;
    
    printf("first grid rate %d\n",grid_rate);
    
    unsigned char *matchtag_check = (unsigned char*)calloc(sizeof(unsigned char),data_size);
    
    for(total_iteration = 1  ; total_iteration <= max_total_iteration ; total_iteration++)
    {
        check_while = 0;
        check_size = total_iteration*grid_rate;
        iter_check = 0;
        while(check_while == 0)
        {
            iter_check ++;
            count_null_cell = 0;
#pragma omp parallel for schedule(guided) reduction(+:count_null_cell)
            for (long index = 0; index < data_size; index++)
            {
                long count_cell = 0;
                long count_matched_cell = 0;
                long null_cell = 0;
                
                double sum_ortho = 0;
                
                long row = floor(index/col_count);
                long col = index%col_count;
                long grid_pos = row*col_count + col;
                
                if (value_pt[grid_pos] == 0 && value[grid_pos] > -100 && !matchtag_check[grid_pos])
                {
                    long total_size = 0;
                    
                    for (long t_i = -check_size; t_i <= check_size;t_i+=grid_rate )
                    {
                        for (long t_j = -check_size; t_j <= check_size; t_j+=grid_rate)
                        {
                            long index_row = row + t_i;
                            long index_col = col + t_j;
                            long t_index = index_row*col_count + index_col;
                            
                            if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count && value[t_index] > -100)
                            {
                                if(value_pt[t_index] == 1 && SignedCharToFloat(Ortho_values[t_index]) > 0.0)
                                {
                                    sum_ortho += SignedCharToFloat(Ortho_values[t_index]);
                                    count_cell++;
                                }
                                
                                if(value_pt[t_index] == 0)
                                    count_matched_cell++;
                                
                                total_size++;
                            }
                        }
                    }
                    
                    if(count_matched_cell == total_size)
                        matchtag_check[grid_pos] = true;
                    
                    if (count_cell >= 7 && count_cell >= total_size*th_rate )
                    {
                        t_value_orthoncc[grid_pos] = FloatToSignedChar(sum_ortho/(double)count_cell);
                        t_value_pt[grid_pos] = 1;
                        count_null_cell ++;
                    }
                }
            }
            
            memcpy(Ortho_values,t_value_orthoncc,sizeof(signed char)*data_size);
            memcpy(value_pt,t_value_pt,sizeof(unsigned char)*data_size);
            
            if(count_null_cell == 0 || iter_check > max_mt_iteration)
                check_while = 1;
        }
    }
    
    memset(matchtag_check,0,sizeof(unsigned char)*data_size);
    
    check_size = 10*grid_rate;
    total_iteration = 0;
    printf("second null\n");
    check_while = 0;
    iter_check = 0;
    while(check_while == 0)
    {
        printf("2nd null iteration %d\t%d\n",check_size,iter_check);
        iter_check ++;
        count_null_cell = 0;
        int count_highnull_cell = 0;
#pragma omp parallel for schedule(guided) reduction(+:count_null_cell,count_highnull_cell)
        for (long index = 0; index < data_size ; index++)
        {
            int count_low_cell = 0;
            int count_high_cell = 0;
            
            int count_matched_cell = 0;
            int count_unmatched_cell = 0;
            
            int total_size = 0;
            
            long row = floor(index/col_count);
            long col = index%col_count;
            
            long grid_pos = row*col_count + col;
            
            if (value_pt[grid_pos] == 1 && value[grid_pos] > -100 && !matchtag_check[grid_pos])
            {
                for (long t_i = -check_size; t_i <= check_size;t_i+=grid_rate )
                {
                    for (long t_j = -check_size; t_j <= check_size; t_j+=grid_rate)
                    {
                        long index_row = row + t_i;
                        long index_col = col + t_j;
                        long t_index = index_row*col_count + index_col;
                        
                        if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
                        {
                            if(value_pt[t_index] == 0)
                            {
                                if(SignedCharToFloat(Ortho_values[t_index]) < 0.0)
                                    count_low_cell++;
                                else if(SignedCharToFloat(Ortho_values[t_index]) > 0.0)
                                    count_high_cell++;
                            }
                            
                            if(value_pt[t_index] == 1)
                                count_matched_cell++;
                            
                            total_size++;
                        }
                    }
                }
                
                if(count_matched_cell == total_size)
                    matchtag_check[grid_pos] = true;
                else
                {
                    if(count_low_cell >= total_size*0.2)
                    {
                        t_value_pt[grid_pos] = 0;
                        count_null_cell ++;
                    }
                    
                    if(count_high_cell >= total_size*0.7)
                    {
                        t_value_pt[grid_pos] = 2;
                        count_highnull_cell++;
                    }
                }
            }
        }
        printf("%ld\t%d\n",count_null_cell,count_highnull_cell);
        if(count_null_cell == 0 || iter_check > max_mt_iteration)
            check_while = 1;
        
        memcpy(value_pt,t_value_pt,sizeof(unsigned char)*data_size);
    }
    
#pragma omp parallel for schedule(guided)
    for (long index = 0; index < data_size ; index++)
    {
        if(t_value_pt[index] > 0)
        {
            t_value_pt[index] = 1;
            value_pt[index] = 1;
        }
    }
    
    memset(matchtag_check,0,sizeof(unsigned char)*data_size);
    
    if(grid_rate > 1)
    {
        printf("last iteration\n");
        
        check_size = 0;
        max_total_iteration = 3;
        total_iteration;
        
        for(total_iteration = 1  ; total_iteration <= max_total_iteration ; total_iteration++)
        {
            check_while = 0;
            check_size = total_iteration;
            iter_check = 0;
            
            while(check_while == 0)
            {
                iter_check ++;
                count_null_cell = 0;
                
#pragma omp parallel for schedule(guided) reduction(+:count_null_cell)
                for (long index = 0; index < data_size ; index++)
                {
                    long count_cell = 0;
                    long count_all_cell = 0;
                    double sum_h = 0;
                    double sum_ortho = 0;
                    
                    int total_size = 0;
                    
                    long row = floor(index/col_count);
                    long col = index%col_count;
                    
                    long grid_pos = row*col_count + col;
                    if (value_pt[grid_pos] == 0 && value[grid_pos] > -100 && !matchtag_check[grid_pos])
                    {
                        for (long t_i = -check_size; t_i <= check_size;t_i++ )
                        {
                            for (long t_j = -check_size; t_j <= check_size; t_j++)
                            {
                                long index_row = row + t_i;
                                long index_col = col + t_j;
                                long t_index = index_row*col_count + index_col;
                                
                                if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
                                {
                                    if(value_pt[t_index] > 0 && SignedCharToFloat(Ortho_values[t_index]) > 0.0 && value[t_index] > -100)
                                    {
                                        sum_ortho += SignedCharToFloat(Ortho_values[t_index]);
                                        count_cell++;
                                    }
                                    
                                    if(value_pt[t_index] == 0)
                                        count_all_cell++;
                                    
                                    total_size++;
                                }
                            }
                        }
                        
                        if(count_all_cell == total_size)
                            matchtag_check[grid_pos] = true;
                        
                        if (count_cell >= total_size*th_rate )
                        {
                            t_value_orthoncc[grid_pos] = FloatToSignedChar(sum_ortho/(double)count_cell);
                            
                            t_value_pt[grid_pos] = 1;
                            
                            count_null_cell ++;
                        }
                    }
                }
                memcpy(Ortho_values,t_value_orthoncc,sizeof(signed char)*data_size);
                memcpy(value_pt,t_value_pt,sizeof(unsigned char)*data_size);
                
                if(count_null_cell == 0 || iter_check > max_mt_iteration)
                    check_while = 1;
            }
        }
    }
    free(matchtag_check);
    
    free(t_value_orthoncc);
    free(Ortho_values);
    
    printf("end last iteration\n");
#pragma omp parallel for schedule(guided)
    for (long index = 0; index < data_size ; index++)
    {
        int count_null_grid = 0;
        int count_match_grid = 0;
        
        long row = floor(index/col_count);
        long col = index%col_count;
        long count_cell = 0;
        long null_count_cell = 0;
        
        long p_index = row*col_count + col;
        
        if(value[p_index] > Nodata)
        {
            for (long t_i = -1; t_i <= 1;t_i++ )
            {
                for (long t_j = -1; t_j <= 1; t_j++)
                {
                    long index_row = row + t_i;
                    long index_col = col + t_j;
                    long t_index = index_row*col_count + index_col;
                    
                    if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
                    {
                        if (value_pt[p_index] == 0 && value_pt[t_index] > 0)
                            count_match_grid++;
                        
                        if (value_pt[p_index] == 1 && value_pt[t_index] == 0)
                            count_null_grid++;
                    }
                }
            }
        }
        
        if(count_match_grid >= 8)
            t_value_pt[p_index] = 1;
        if(count_null_grid >= 8)
            t_value_pt[p_index] = 0;
    }
    
    free(value);
    free(value_pt);
    
    printf("end smoothing\n");
    printf("end null\n");
    
    char GEOTIFF_matchtag_filename[500];
    if(divide == 0)
        sprintf(GEOTIFF_matchtag_filename, "%s/%s_matchtag.tif", proinfo->save_filepath, proinfo->Outputpath_name);
    else
        sprintf(GEOTIFF_matchtag_filename, "%s/%s_%d_matchtag.tif", proinfo->save_filepath, proinfo->Outputpath_name,divide);
    
    WriteGeotiff(GEOTIFF_matchtag_filename, t_value_pt, col_count, row_count, proinfo->DEM_resolution, minX, maxY, _param.projection, _param.utm_zone, _param.bHemisphere, 1);
    
    free(t_value_pt);
}

double FindNebPts_F_M_IDW(const float *input, const unsigned char *matching_flag, const long row_size, const long col_size, const double grid, const double minX, const double maxY, const double X, const double Y, const int row_interval)
{
    double result;
    
    long row,col;
    int check_stop = 0;

    typedef struct tagKV
    {
        double diff;
        double height;
        
        tagKV(double diff, double height):diff(diff),height(height)
        {
        }
    }KV;
    
    vector<KV> Kernel;
    
    double col_pos = ((X - minX)/grid);
    double row_pos = ((maxY - Y)/grid);
    
    long interval = 1;
    int numpts = 0;
    int count1(0), count2(0), count3(0), count4(0);
    while (check_stop == 0)
    {
        //top
        row = interval;
        for(col = -interval ; col <= interval ; col++)
        {
            long grid_pos = ((row_pos+row)*col_size + (col_pos+col));
            if(grid_pos >= 0 && grid_pos < row_size*col_size &&
               row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size && col != 0 && row != 0)
            {
                if(input[grid_pos] != Nodata && matching_flag[grid_pos] == 1)
                {
                    numpts++;
                    if (row >= 0 && row <=interval && col >= 0 && col <= interval)
                        count1++;
                    
                    if (row >= 0 && row <=interval && col < 0 && col >= -interval)
                        count2++;
                    
                    if (row < 0 && row >= -interval && col < 0 && col >= -interval)
                        count3++;
                    
                    if (row < 0 && row >= -interval && col >= 0 && col <= interval)
                        count4++;
                    
                    double dis_X = (col*grid);
                    double dis_Y = (row*grid);
                    KV t_kv(sqrt(dis_X*dis_X + dis_Y*dis_Y), input[grid_pos]);
                    Kernel.push_back(t_kv);
                }
            }
        }
        
        //bottom
        row = -interval;
        for(col = -interval ; col <= interval ; col++)
        {
            long grid_pos = ((row_pos+row)*col_size + (col_pos+col));
            if(grid_pos >= 0 && grid_pos < row_size*col_size &&
               row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size && col != 0 && row != 0)
            {
                if(input[grid_pos] != Nodata && matching_flag[grid_pos] == 1)
                {
                    numpts++;
                    
                    if (row >= 0 && row <=interval && col >= 0 && col <= interval)
                        count1++;
                    
                    if (row >= 0 && row <=interval && col < 0 && col >= -interval)
                        count2++;
                    
                    if (row < 0 && row >= -interval && col < 0 && col >= -interval)
                        count3++;
                    
                    if (row < 0 && row >= -interval && col >= 0 && col <= interval)
                        count4++;
                    
                    double dis_X = (col*grid);
                    double dis_Y = (row*grid);
                    KV t_kv(sqrt(dis_X*dis_X + dis_Y*dis_Y), input[grid_pos]);
                    Kernel.push_back(t_kv);
                }
            }
        }
        
        //right
        col = interval;
        for(row = -interval+1 ; row <= interval-1 ; row++)
        {
            int grid_pos = (int)((row_pos+row)*col_size + (col_pos+col));
            if(grid_pos >= 0 && grid_pos < row_size*col_size &&
               row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size && col != 0 && row != 0)
            {
                if(input[grid_pos] != Nodata && matching_flag[grid_pos] == 1)
                {
                    numpts++;
                    
                    if (row >= 0 && row <=interval && col >= 0 && col <= interval)
                        count1++;
                    
                    if (row >= 0 && row <=interval && col < 0 && col >= -interval)
                        count2++;
                    
                    if (row < 0 && row >= -interval && col < 0 && col >= -interval)
                        count3++;
                    
                    if (row < 0 && row >= -interval && col >= 0 && col <= interval)
                        count4++;
                   
                    double dis_X = (col*grid);
                    double dis_Y = (row*grid);
                    KV t_kv(sqrt(dis_X*dis_X + dis_Y*dis_Y), input[grid_pos]);
                    Kernel.push_back(t_kv);
                }
            }
        }
        
        //left
        col = -interval;
        for(row = -interval+1 ; row <= interval-1 ; row++)
        {
            int grid_pos = (int)((row_pos+row)*col_size + (col_pos+col));
            if(grid_pos >= 0 && grid_pos < row_size*col_size &&
               row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size && col != 0 && row != 0)
            {
                if(input[grid_pos] != Nodata && matching_flag[grid_pos] == 1)
                {
                    numpts++;
                    
                    if (row >= 0 && row <=interval && col >= 0 && col <= interval)
                        count1++;
                    
                    if (row >= 0 && row <=interval && col < 0 && col >= -interval)
                        count2++;
                    
                    if (row < 0 && row >= -interval && col < 0 && col >= -interval)
                        count3++;
                    
                    if (row < 0 && row >= -interval && col >= 0 && col <= interval)
                        count4++;
                    
                    double dis_X = (col*grid);
                    double dis_Y = (row*grid);
                    KV t_kv(sqrt(dis_X*dis_X + dis_Y*dis_Y), input[grid_pos]);
                    Kernel.push_back(t_kv);
                }
            }
        }
        
        if (interval >= row_interval || ((numpts) >= 10 && count1 >= 2 && count2 >= 2 && count3 >= 2 && count4 >= 2))
            check_stop = 1;
        else
            interval = interval + 1;
    }
    
    double sum1(0), sum2(0);
    double p = 1.5;

    vector<KV>::iterator it;
    
    for(it = Kernel.begin(); it != Kernel.end() ; ++it)
    {
        double height = it->height;
        double diff = it->diff;
        sum1 += (height/pow(diff,p));
        sum2 += (1.0/pow(diff,p));
    }
    
    if(sum2 > 0)
        result = sum1/sum2;
    else
        result = Nodata;
   
    Kernel.clear();
    vector<KV>().swap(Kernel);
 
    return result;
}


// Find the interpolated value of a patch given the nominal position and the X and Y offsets 
// along with the image itself

















