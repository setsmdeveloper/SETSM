/*
 * Copyright 2017 Myoung-Jong Noh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *	   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Includes code derived from the voronoi algorithm by Steven Fortune
 * (http://ect.bell-labs.com/who/sjf/)
 * as modified by Derek Bradley
 * (http://zurich.disneyresearch.com/derekbradley/voronoi.html)
 *
 * Reference: Steve J. Fortune (1987) A Sweepline Algorithm for Voronoi Diagrams,
 * Algorithmica 2, 153-174.
 */

#define REAL double

#include "setsm_code.h"
#ifdef TRILIBRARY
#include "triangle.h"
#else
#define VLIBRARY
#include "voronoi_setsm.h"
#endif
#include "math.h"
#include <omp.h>
#include <time.h>
#include <dirent.h>
#include <libgen.h>
#include <sys/stat.h>
#ifdef BUILDMPI
#include "mpi.h"
#endif

char *dirname(char *path);

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

int main(int argc,char *argv[])
{
	setbuf(stdout, NULL);
	TIFFSetWarningHandler(NULL);
	char* projectfilename	= "default.txt";
	int i=0;
	char* image1_name = NULL;
	char* image2_name = NULL;
	char* output_directory_name = NULL;
	ARGINFO args;
	
	args.check_arg = 0;
	args.check_DEM_space = false;
	args.check_Threads_num = false;
	args.check_seeddem = false;
	args.check_minH = false;
	args.check_maxH = false;
	args.check_tiles_SR	= false;
	args.check_tiles_ER	= false;
	args.check_tiles_SC	= false;
	args.check_tiles_EC	= false;
	args.check_RA_line	= false;
	args.check_RA_sample= false;
	args.check_gridonly = false;
	args.check_RA_tileR = false;
	args.check_RA_tileC = false;
	args.check_tilesize = false;
	args.check_boundary = false;
	args.check_checktiff = false;
	args.check_ortho = false;
	args.check_RA_only = false;
    
	args.projection = 3;//PS = 1, UTM = 2
    args.sensor_provider = 1; //DG = 1, Pleiades = 2
    args.check_imageresolution = false;
    args.utm_zone = -99;
    args.ortho_count = 1;
    args.overlap_length = 50;
    args.RA_only = 0;
    
	TransParam param;
	
	if(argc == 1)
	{
		char save_filepath[500];
		char LeftImagefilename[500];
		
		args.check_arg = 0;
		SETSMmainfunction(&param,projectfilename,args,LeftImagefilename,save_filepath);
		
		char DEMFilename[500];
		char Outputpath[500];
		sprintf(DEMFilename,"%s/tif/DEM_nn.raw",save_filepath);
		sprintf(Outputpath,"%s/tif",save_filepath);
			
		printf("%s\n",LeftImagefilename);
		printf("%s\n",DEMFilename);
		printf("%s\n",Outputpath);
		orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1);
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
			
			printf("setsm version : 3.2.8\n");
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
			printf("\t[-RAonly value]\t: If set to 1 (true), program will exit after RA comprehension. Default = 0 (false)\n");
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
			
			SETSMmainfunction(&param,projectfilename,args,LeftImagefilename,save_filepath);
		}
	}
	else if(argc == 4)
	{
		args.check_arg = 1;
		sprintf(args.Image1,"%s",argv[1]);
		sprintf(args.Image2,"%s",argv[2]);
		sprintf(args.Outputpath,"%s",argv[3]);
		
		char *Outputpath_name  = SetOutpathName(args.Outputpath);
		sprintf(args.Outputpath_name,"%s",Outputpath_name);
		printf("after pathname %s\n",args.Outputpath_name);
		
		printf("%s\n",args.Image1);
		printf("%s\n",args.Image2);
		printf("%s\n",args.Outputpath);
		printf("%s\n", args.Outputpath_name);
		
		char save_filepath[500];
		char LeftImagefilename[500];
		
		if( strcmp(args.Image1,args.Image2) != 0)
		{
			SETSMmainfunction(&param,projectfilename,args,LeftImagefilename,save_filepath);
			
			char DEMFilename[500];
			char Outputpath[500];
			sprintf(DEMFilename, "%s/%s_dem.raw", save_filepath,args.Outputpath_name);
			sprintf(Outputpath, "%s", save_filepath);
			orthogeneration(param,args,LeftImagefilename, DEMFilename, Outputpath,1);
		}
		else
			printf("Please check input 1 and input 2. Both is same\n");
			
	}
	else if(argc > 4)
	{
		bool cal_flag = true;
		args.check_arg = 1;
		sprintf(args.Image1,"%s",argv[1]);
		sprintf(args.Image2,"%s",argv[2]);
		sprintf(args.Outputpath,"%s",argv[3]);
		
		char *Outputpath_name  = SetOutpathName(args.Outputpath);
		sprintf(args.Outputpath_name,"%s",Outputpath_name);
		printf("after pathname %s\n",args.Outputpath_name);
		
		printf("%s\n",args.Image1);
		printf("%s\n",args.Image2);
		printf("%s\n",args.Outputpath);
		printf("%s\n", args.Outputpath_name);
		
		
		bool bminx	= false;
		bool bmaxx	= false;
		bool bminy	= false;
		bool bmaxy	= false;
		
		for (i=4; i<argc; i++)
		{
            if (strcmp("-provider", argv[i]) == 0) {
                if (argc == i + 1) {
                    printf("Please input Provider info\n");
                    cal_flag = false;
                } else {
                    if(strcmp("DG",argv[i+1]) == 0 || strcmp("dg",argv[i+1]) == 0)
                    {
                        args.sensor_provider = 1;
                        printf("Image Provider : Digital Globe\n");
                        
                    }
                    else if(strcmp("Pleiades",argv[i+1]) == 0 || strcmp("pleiades",argv[i+1]) == 0)
                    {
                        args.sensor_provider = 2;
                        printf("Image Provider : Pleiades\n");
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
                    
                    /*if(args.DEM_space < 1.0)
                    {
                        args.DEM_space = 1.0;
                        printf("Minimum size of DEM grid is 1m. outres set 1.0\n");
                    }
                     */
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
 
			if (strcmp("-seed",argv[i]) == 0)
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
					tif	 = TIFFOpen(args.seedDEMfilename,"r");
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
						
						tif	 = TIFFOpen(args.seedDEMfilename,"r");
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
					}

					if(args.check_seeddem)
					{
						char* temp_path = remove_ext(args.seedDEMfilename);
						printf("seedem %s\n",temp_path);
						int full_size;
						full_size		= strlen(args.seedDEMfilename);
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
						char *str = "meta.txt";
						
						sprintf(args.metafilename,"%s_%s",Metafile,str);
						printf("Meta file %s\n",args.metafilename);
						
						FILE* pFile_meta;
						pFile_meta	= fopen(args.metafilename,"r");
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
							pFile_meta1	= fopen(args.metafilename,"r");
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
					args.ra_line = atof(argv[i+1]);
					printf("%f\n",args.ra_line);
					args.check_RA_line = true;
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
					args.ra_sample = atof(argv[i+1]);
					printf("%f\n",args.ra_sample);
					args.check_RA_sample = true;
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
					
					bminx	= true;
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
					
					bminy	= true;
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
					
					bmaxx	= true;
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
					
					bmaxy	= true;
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
			
            
            if(args.check_checktiff)
			{
				SETSMmainfunction(&param,projectfilename,args,LeftImagefilename,save_filepath);
			}
			else if( strcmp(args.Image1,args.Image2) != 0)
			{
				SETSMmainfunction(&param,projectfilename,args,LeftImagefilename,save_filepath);

				char DEMFilename[500];
				char Outputpath[500];
				sprintf(DEMFilename, "%s/%s_dem.raw", save_filepath,args.Outputpath_name);
				sprintf(Outputpath, "%s", save_filepath);
				
				printf("param %s %d\n", param.direction,param.zone);
				param.projection = args.projection;
		orthogeneration(param,args,args.Image1, DEMFilename, Outputpath,1);
		if(!args.check_ortho)
		    orthogeneration(param,args,args.Image2, DEMFilename, Outputpath,2);
		else if(args.ortho_count == 2)
		    orthogeneration(param,args,args.Image2, DEMFilename, Outputpath,2);
			}
			else
				printf("Please check input 1 and input 2. Both is same\n");
		}
	}
	
	printf("# of allocated threads = %d\n",omp_get_max_threads());
			
	
	return 0;
}

char* SetOutpathName(char *_path)
{
	char *t_name;
	char lastSlash[500];
	strcpy(lastSlash, _path);
	printf("path %s\n",lastSlash);
	char *fullseeddir = dirname(lastSlash);
	printf("path %s\n",lastSlash);
	int dir_size = strlen(fullseeddir);
	int full_size = strlen(_path);
	printf("fullseeddir %s\tdir_size %d\n", fullseeddir, dir_size);
	printf("lastSlash %s\tfull_size %d\n", lastSlash, full_size);
	
	if(dir_size > 1)
	{
		int start = dir_size+1;
		int end = full_size;
		int lenth = end	 - start+1;
		t_name = (char*)(malloc(sizeof(char)*(lenth+1)));
		int path_size = strlen(t_name);
		for (int i = 0; i < lenth; i++) {
			t_name[i] = fullseeddir[i+start];
		}
		t_name[lenth] = '\0';
	}
	else {
		t_name = (char*)(malloc(sizeof(char)*(full_size+1)));
		strcpy(t_name,_path);
	}
	
	printf("Outputpath_name %s\n",t_name);
	
	return t_name;
	
}

void SETSMmainfunction(TransParam *return_param, char* _filename, ARGINFO args, char *_LeftImagefilename, char *_save_filepath)
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

	char computation_file[500];
	time_t total_ST = 0, total_ET = 0;
	double total_gap;
	FILE *time_fid;
	
	total_ST = time(0);
	
	
	ProInfo proinfo;
	if(OpenProject(_filename,&proinfo,args))
	{
		double Boundary[4]	= {0.0};
		double LBoundary[4],RBoundary[4],LminmaxHeight[2],RminmaxHeight[2],ori_minmaxHeight[2];
		double LHinterval[1], RHinterval[1], Hinterval;
		double Image_res[2]	= {0.0};
		double Res[2]		= {0.0};
		double Limageparam[2] = {0.0};
		double Rimageparam[2] = {0.0};
		int final_iteration;
		
		double **LRPCs, **RRPCs, minLat, minLon;
        ImageInfo leftimage_info;
        ImageInfo rightimage_info;
        
		uint8 pre_DEM_level = 0;
		uint8 DEM_level		= 0;
		uint8 NumOfIAparam	= 2;

		int i;

		bool check_tile_array = false;
		bool Hemisphere;
		bool* tile_array = NULL;

		CSize Limagesize, Rimagesize;//original imagesize
		CSize LBRsize, RBRsize;//DEM boudary size

		TransParam param;

		sprintf(_LeftImagefilename,"%s",proinfo.LeftImagefilename);
		sprintf(_save_filepath,"%s",proinfo.save_filepath);
		

		printf("Completion of loading project file!!\n");

		printf("# of detected threads by openmp = %d\n",omp_get_max_threads());
		
		if (proinfo.threads_num != 0) {
			omp_set_num_threads(proinfo.threads_num);
		}
		
		printf("# of allocated threads = %d\n",omp_get_max_threads());
		
		if(Maketmpfolders(proinfo))
		{
			char metafilename[500];
			
			FILE *pMetafile = NULL;
			sprintf(metafilename, "%s/%s_meta.txt", proinfo.save_filepath, proinfo.Outputpath_name);
			
			if(!proinfo.check_checktiff && !args.check_ortho)
			{
				pMetafile	= fopen(metafilename,"w");
			
				fprintf(pMetafile,"SETSM Version=3.2.8\n");
			}
			
			time_t current_time;
			char*	c_time_string;
			
			current_time = time(NULL);
			c_time_string = ctime(&current_time);
			
            char temp_filepath[500];
            double Image1_gsd_r,Image1_gsd_c,Image2_gsd_r,Image2_gsd_c, Image1_gsd, Image2_gsd;
            ImageGSD GSD_image1, GSD_image2;
            BandInfo left_band, right_band;
            
            if(args.sensor_provider == 1)
            {
                LRPCs		= OpenXMLFile(proinfo.LeftRPCfilename,&Image1_gsd_r,&Image1_gsd_c,&Image1_gsd,&left_band);
                RRPCs		= OpenXMLFile(proinfo.RightRPCfilename,&Image2_gsd_r,&Image2_gsd_c,&Image2_gsd,&right_band);
                
                GSD_image1.row_GSD = Image1_gsd_r;
                GSD_image1.col_GSD = Image1_gsd_c;
                GSD_image1.pro_GSD = Image1_gsd;
                GSD_image2.row_GSD = Image2_gsd_r;
                GSD_image2.col_GSD = Image2_gsd_c;
                GSD_image2.pro_GSD = Image2_gsd;
                
                OpenXMLFile_orientation(proinfo.LeftRPCfilename,&leftimage_info);
                OpenXMLFile_orientation(proinfo.RightRPCfilename,&rightimage_info);
                
                printf("leftimage info\nSatID = %s\nAcquisition_time = %s\nMean_row_GSD = %f\nMean_col_GSD = %f\nMean_GSD = %f\nMean_sun_azimuth_angle = %f\nMean_sun_elevation = %f\nMean_sat_azimuth_angle = %f\nMean_sat_elevation = %f\nIntrack_angle = %f\nCrosstrack_angle = %f\nOffnadir_angle = %f\ntdi = %d\neffbw = %f\nabscalfact = %f\n",leftimage_info.SatID,leftimage_info.imagetime,Image1_gsd_r,Image1_gsd_c,Image1_gsd,leftimage_info.Mean_sun_azimuth_angle,leftimage_info.Mean_sun_elevation,leftimage_info.Mean_sat_azimuth_angle,leftimage_info.Mean_sat_elevation,leftimage_info.Intrack_angle,leftimage_info.Crosstrack_angle,leftimage_info.Offnadir_angle,(int)left_band.tdi,left_band.effbw,left_band.abscalfactor);
                printf("rightimage info\nSatID = %s\nAcquisition_time = %s\nMean_row_GSD = %f\nMean_col_GSD = %f\nMean_GSD = %f\nMean_sun_azimuth_angle = %f\nMean_sun_elevation = %f\nMean_sat_azimuth_angle = %f\nMean_sat_elevation = %f\nIntrack_angle = %f\nCrosstrack_angle = %f\nOffnadir_angle = %f\ntdi = %d\neffbw = %f\nabscalfact = %f\n",rightimage_info.SatID,rightimage_info.imagetime,Image2_gsd_r,Image2_gsd_c,Image2_gsd,rightimage_info.Mean_sun_azimuth_angle,rightimage_info.Mean_sun_elevation,rightimage_info.Mean_sat_azimuth_angle,rightimage_info.Mean_sat_elevation,rightimage_info.Intrack_angle,rightimage_info.Crosstrack_angle,rightimage_info.Offnadir_angle,(int)right_band.tdi,right_band.effbw,right_band.abscalfactor);
            }
            else
            {
                LRPCs		= OpenXMLFile_Pleiades(proinfo.LeftRPCfilename);
                RRPCs		= OpenXMLFile_Pleiades(proinfo.RightRPCfilename);
            }
            

            if(!args.check_imageresolution)
            {
                if(args.sensor_provider == 1)
                {
                    proinfo.resolution = (int)(((Image1_gsd_r + Image1_gsd_c + Image2_gsd_r + Image2_gsd_c)/4.0)*10 + 0.5)/10.0;
                    
                    if (proinfo.resolution < 0.75)
                    {
                        proinfo.resolution = 0.5;
                    }
                    else if(proinfo.resolution < 1.25)
                        proinfo.resolution = 1.0;
                    else
                        proinfo.resolution = floor(proinfo.resolution);
                }
                else
                    proinfo.resolution = 0.5;
            }
            else
            {
                proinfo.resolution  = args.image_resolution;
            }
            
            printf("image resolution %f\n",proinfo.resolution);
            
            if (!args.check_DEM_space)
            {
                proinfo.DEM_resolution = proinfo.resolution;
            }
            
            //if(proinfo.DEM_resolution < 1.0)
            //    proinfo.DEM_resolution = 1.0;
            
            if(!proinfo.check_checktiff && !args.check_ortho)
            {
                fprintf(pMetafile,"Creation Date=%s",c_time_string);
                fprintf(pMetafile,"Image 1=%s\n",proinfo.LeftImagefilename);
                fprintf(pMetafile,"Image 2=%s\n",proinfo.RightImagefilename);
                fprintf(pMetafile,"Output Resolution=%f\n",proinfo.DEM_resolution);
            }
            
			Res[0]		= proinfo.resolution;						Res[1]		= proinfo.DEM_resolution;
			Image_res[0]= proinfo.resolution;						Image_res[1]= proinfo.resolution;

			Rimageparam[0]	= proinfo.RA_param[0];					Rimageparam[1]	= proinfo.RA_param[1];

			minLat		= LRPCs[0][3];
			minLon = (double) LRPCs[0][2];
			param.projection = args.projection;
            param.utm_zone   = args.utm_zone;
            
			SetTransParam((double)(minLat),(double)(minLon),&Hemisphere, &param);

			printf("param projection %d\tzone %d\n",param.projection,param.utm_zone);
			*return_param = param;
			
			SetDEMBoundary(LRPCs,Image_res,param,Hemisphere,LBoundary,LminmaxHeight,&LBRsize,LHinterval);
			SetDEMBoundary(RRPCs,Image_res,param,Hemisphere,RBoundary,RminmaxHeight,&RBRsize,RHinterval);
		
			
			if(LHinterval[0] > RHinterval[0])
				Hinterval	= LHinterval[0];
			else
				Hinterval	= RHinterval[0];
		
			GetImageSize(proinfo.LeftImagefilename,&Limagesize);	GetImageSize(proinfo.RightImagefilename,&Rimagesize);

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
				for(i=0;i<4;i++)
				{
					if(i<2)
						Boundary[i] = ceil((max(LBoundary[i], RBoundary[i]) / 2.0)) * 2;
					else
						Boundary[i] = floor((min(LBoundary[i], RBoundary[i]) / 2.0)) * 2;
				}
			}
			
			printf("boundary = %f\t%f\t%f\t%f\n",Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
			
            
            CSize Boundary_size;
            Boundary_size.width     = Boundary[2] - Boundary[0];
            Boundary_size.height    = Boundary[3] - Boundary[1];
            
            if(Boundary_size.height/1000.0 > args.overlap_length || Boundary_size.width/1000.0 > args.overlap_length)
            {
                printf("Overlapped area between stereo pair is very long along the strip(height=%3.2f(km), width=%3.2f(km)), so that the assumption of RPC bias computation (less than 50km) is not satisfied,\nso relative RPC bias can be not accurately compensated. \nPlease process after split the overlapped area in strip direction into several small area less than 30 km\nBounary(minX, minY, maxX, maxY[m]) = %f %f %f %f\n",Boundary_size.height/1000.0,Boundary_size.width/1000.0,Boundary[0],Boundary[1],Boundary[2],Boundary[3]);
                exit(1);
            }
            
            
			ori_minmaxHeight[0] = min(LminmaxHeight[0],RminmaxHeight[0]);
			ori_minmaxHeight[1] = max(LminmaxHeight[1],RminmaxHeight[1]);

			if (args.check_minH) {
				ori_minmaxHeight[0] = (int)args.minHeight;
				printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
			}
			if (args.check_maxH) {
				ori_minmaxHeight[1] = (int)args.maxHeight;
				printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
			}
			
			if (proinfo.check_minH) {
				ori_minmaxHeight[0] = (int)proinfo.minHeight;
				printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
			}
			if (proinfo.check_maxH) {
				ori_minmaxHeight[1] = (int)proinfo.maxHeight;
				printf("minmaxH = %f\t%f\n", ori_minmaxHeight[0], ori_minmaxHeight[1]);
			}
			
			if (ori_minmaxHeight[1] >9000) {
				ori_minmaxHeight[1] = 9000;
			}
			
            
            if(!args.check_ortho)
            {
                printf("minmaxH = %f\t%f\n",ori_minmaxHeight[0],ori_minmaxHeight[1]);
                printf("seed fff %d\n",proinfo.pre_DEMtif);
                if(SetupParam(proinfo,&NumOfIAparam, &pre_DEM_level, &DEM_level,&proinfo.pre_DEMtif,&check_tile_array ))
                {
                    uint8 pyramid_step;
                    uint8 Template_size	= 15;
                    uint16 buffer_area	= 400;
                    
                    uint8 iter_row_start, iter_row_end, t_col_start, t_col_end;
                    int	 total_count = 0, tile_size = 0;
                    double subX, subY;
                    double bin_angle;
                    double mt_grid_size;
                    char temp_c[500];
                    double temp_br_x, temp_br_y;
                    CSize temp_size;
                    time_t ST = 0, ET = 0;
                    double gap;
                
                    printf("IsRA = %d\n",proinfo.IsRA);
                    
                    if(proinfo.IsRA)
                    {
                        uint8 RA_row_iter = 1;
                        uint8 RA_col_iter = 1;

                        tile_size			= 8000;
                        if(args.check_tilesize)
                            tile_size		= args.tilesize*2;
                        printf("tileszie %d\n",tile_size);

                        bin_angle			= 360.0/18.0;
                        bool check_RA_1000	= false;
                        
                        if (args.check_RA_tileR)
                        {
                            if(args.RA_row == 100)
                            {
                                tile_size			= 50000;
                                check_RA_1000		= true;
                                args.RA_row			= 1;
                            }
                        }
                        
                        

                        SetTiles_RA(proinfo,proinfo.IsSP,proinfo.IsRR, Boundary, Res, tile_size, proinfo.pre_DEMtif, &pyramid_step, &buffer_area, 
                                    &iter_row_start, &iter_row_end,&RA_row_iter, &t_col_start, &t_col_end, &RA_col_iter, &subX, &subY);	

                        
                        printf("before RA row:col => row = %d\t%d\t;col = %d\t%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end);
                        
                        if (args.check_RA_tileR)
                        {
                            iter_row_start	= args.RA_row;
                            iter_row_end	= iter_row_start+1;
                        }
                        
                        if (args.check_RA_tileC)
                        {
                            t_col_start		= args.RA_col;
                            t_col_end		= t_col_start + 1;
                        }
                        
                        total_count			= 0;
                        printf("RA row:col = row = %d\t%d\t;col = %d\t%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end);
                        
                        
                        //loading existing RA parameters
                        bool check_load_RA = false;
                        char str_rafile[500];
                        char bufstr[500];
                        sprintf(str_rafile,"%s/txt/RA_echo_result_row_%d_col_%d.txt",proinfo.save_filepath,iter_row_start,t_col_start);
                        printf("RA file %s\n",str_rafile);
                        FILE* pFile;
                        pFile		= fopen(str_rafile,"r");
                        if(pFile)
                        {
                            char str_echofile[500];
                            sprintf(str_echofile,"%s/txt/echo_result_row_1_col_1.txt",proinfo.save_filepath);
                            FILE* pFile_echo;
                            
                            printf("echo %s\n",str_echofile);
                            
                            pFile_echo	= fopen(str_echofile,"r");
                            if(pFile_echo)
                            {
                                printf("open RA\n");
                                fgets(bufstr,500,pFile_echo);
                                if (strstr(bufstr,"RA param X")!=NULL)
                                {
                                    sscanf(bufstr,"RA param X = %lf Y = %lf\n",&Rimageparam[0],&Rimageparam[1]);
                                    if(Rimageparam[0] != 0 && Rimageparam[1] != 0)
                                        check_load_RA = true;
                                }
                                fclose(pFile_echo);
                            }
                            
                            fclose(pFile);
                        }
                        
                        printf("check seedem %d\n",args.check_seeddem);
                        if(args.check_seeddem)
                        {
                            int dir_size;
                            int full_size;
                            char str_rafile_1[500];
                            char *fullseeddir	= NULL;
                            char *lastSlash = NULL;
                            lastSlash	= malloc(strlen(args.seedDEMfilename) + 1);
                            strcpy(lastSlash, args.seedDEMfilename);
                            fullseeddir = dirname(lastSlash);
                            
                            dir_size		= strlen(fullseeddir);
                            full_size		= strlen(args.seedDEMfilename);
                            printf("fullseeddir %s\tdir_size %d\n",fullseeddir,dir_size);
                            printf("lastSlash %s\tfull_size %d\n",lastSlash,full_size);
                            char seeddir[500];
                            for(i=0;i<dir_size-4;i++)
                                seeddir[i] = fullseeddir[i];
                            sprintf(str_rafile_1,"%s/txt/RA_echo_result_row_1_col_1.txt",fullseeddir);
                            printf("RA file %s\n",str_rafile_1);

                            char str_echofile[500];
                            sprintf(str_echofile,"%s/txt/echo_result_row_1_col_1.txt",fullseeddir);
                            FILE* pFile_echo;
                                
                            pFile_echo	= fopen(str_echofile,"r");
                            if(pFile_echo)
                            {
                                fgets(bufstr,500,pFile_echo);
                                if (strstr(bufstr,"RA param X")!=NULL)
                                {
                                    sscanf(bufstr,"RA param X = %lf Y = %lf\n",&Rimageparam[0],&Rimageparam[1]);
                                    if(Rimageparam[0] != 0 && Rimageparam[1] != 0)
                                        check_load_RA = true;
                                }
                                fclose(pFile_echo);
                            }
                            
                            char str_rafile_2[500];
                            char RAfile[500];
                            for (int i = 0; i < full_size - 7; i++) {
                                RAfile[i] = args.seedDEMfilename[i];
                            }
                            char RAfile_raw[500];
                            for (int i = 0; i < full_size - 14; i++) {
                                RAfile_raw[i] = args.seedDEMfilename[i];
                            }

                            sprintf(str_rafile_2, "%s/txt/RAinfo.txt", RAfile_raw);
                            printf("Meta file %s\n", args.metafilename);
                            printf("Meta file %s\n", proinfo.metafilename);
                            printf("RA file %s\n", str_rafile_2);

                            if(!check_load_RA)
                            {
                                FILE* pFile_echo;
                                pFile_echo	= fopen(str_rafile,"r");
                                if(pFile_echo)
                                {
                                    fscanf(pFile_echo,"%lf %lf",&Rimageparam[0],&Rimageparam[1]);
                                    if(Rimageparam[0] != 0 && Rimageparam[1] != 0)
                                        check_load_RA = true;
                                    fclose(pFile_echo);
                                }
                                else
                                {
                                    printf("Meta file %s\n",proinfo.metafilename);
                                    
                                    FILE* pFile_meta;
                                    pFile_meta	= fopen(proinfo.metafilename,"r");
                                    if(pFile_meta)
                                    {
                                        printf("meta file exist!!\n");
                                        while(!feof(pFile_meta))
                                        {
                                            
                                            fgets(bufstr,500,pFile_meta);
                                            if (strstr(bufstr,"RA Params=")!=NULL)
                                            {
                                                printf("%s\n",bufstr);
                                                sscanf(bufstr,"RA Params=%lf\t%lf\n",&Rimageparam[0],&Rimageparam[1]);
                                            }
                                            else if(strstr(bufstr,"SETSM Version=")!=NULL)
                                            {
                                                printf("%s\n",bufstr);
                                                double version;
                                                sscanf(bufstr,"SETSM Version=%lf\n",&version);
                                                printf("version %f\n",version);
                                                
                                                if (version > 2.0128) {
                                                    proinfo.seedDEMsigma = 20;
                                                }
                                                else {
                                                    proinfo.seedDEMsigma = 100;
                                                }
                                                printf("sigma %f\n",proinfo.seedDEMsigma);
                                            }
                                            
                                        }
                                                
                                        if(Rimageparam[0] != 0 && Rimageparam[1] != 0)
                                            check_load_RA = true;
                                        fclose(pFile_meta);
                                    }
                                    else {
                                        printf("meta file doesn't exist!! \n");
                                    }

                                }
                            }
                        }
                        
                        if(!check_load_RA)
                        {
                            FILE* pFile_info;
                            sprintf(str_rafile,"%s/txt/RAinfo.txt",proinfo.save_filepath);
                            printf("RAinfo %s\n",str_rafile);
                            pFile_info		= fopen(str_rafile,"r");
                            if(pFile_info)
                            {
                                printf("open RA\n");
                                fscanf(pFile_info,"%lf\t%lf",&Rimageparam[0],&Rimageparam[1]);
                                if(Rimageparam[0] != 0 && Rimageparam[1] != 0)
                                    check_load_RA = true;
                                fclose(pFile_info);
                            }
                        }
                        
                        printf("check load RA %d %f %f\n",check_load_RA,Rimageparam[0],Rimageparam[1]);
                        
                        if(!check_load_RA)
                        {
                            tile_size			= 50000;
                            
                            args.RA_row			= 1;
                    
                            SetTiles_RA(proinfo,proinfo.IsSP,proinfo.IsRR, Boundary, Res, tile_size, proinfo.pre_DEMtif, &pyramid_step, &buffer_area, 
                                        &iter_row_start, &iter_row_end,&RA_row_iter, &t_col_start, &t_col_end, &RA_col_iter, &subX, &subY);	
                            
                            iter_row_start	= 1;
                            iter_row_end	= iter_row_start+1;
                            t_col_start		= 1;
                            t_col_end		= t_col_start + 1;
                            
                            final_iteration = Matching_SETSM(proinfo,pyramid_step, Template_size, buffer_area,iter_row_start, iter_row_end,t_col_start,t_col_end,
                                                             subX,subY,bin_angle,Hinterval,Image_res,Res, Limageparam, Rimageparam,
                                                             LRPCs, RRPCs, pre_DEM_level, DEM_level,	NumOfIAparam, check_tile_array,Hemisphere,tile_array,
                                                             Limagesize,Rimagesize,LBRsize,RBRsize,param,total_count,ori_minmaxHeight,Boundary,RA_row_iter,RA_col_iter);
                        }
                    }
                    
                    if (args.check_RA_line)
                        Rimageparam[0] = args.ra_line;
                    if (args.check_RA_sample)
                        Rimageparam[1] = args.ra_sample;
                    
                    if(!proinfo.check_checktiff)
                    {
                        fprintf(pMetafile,"RA Params=%f\t%f\t\n",Rimageparam[0],Rimageparam[1]);
                        fprintf(pMetafile,"RA tilesize=%d\n",tile_size);
                    }


		    proinfo.IsRA		= false;

					if (!args.RA_only)
					{
		    tile_size			= 4000;
		    
		    if(Boundary_size.width < tile_size && Boundary_size.height < tile_size)
		    {
			if(Boundary_size.width > Boundary_size.height)
			    tile_size = Boundary_size.width;
			else
			    tile_size = Boundary_size.height;
		    }
		    
		    if(args.check_tilesize)
			tile_size		= args.tilesize;
		    printf("tilesize %d\n",tile_size);

		    if(!proinfo.check_checktiff)
		    {
			fprintf(pMetafile,"tilesize=%d\n",tile_size);
			if(proinfo.pre_DEMtif)
			    fprintf(pMetafile,"Seed DEM=%s\n",proinfo.priori_DEM_tif);
			else
			    fprintf(pMetafile,"Seed DEM=\n");
		    
			if(param.projection == 1)
			{
			    if (Hemisphere) {
				fprintf(pMetafile, "Output Projection='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45, +k=1 +x_0=0 +y_0=0 +datum=WGS84 +unit=m +no_defs'\n");
			    } else {
				fprintf(pMetafile, "Output Projection='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0, +k=1 +x_0=0 +y_0=0 +datum=WGS84 +unit=m +no_defs'\n");
			    }
			}
			else {
			    if (Hemisphere) {
				fprintf(pMetafile, "Output Projection='+proj=utm +zone=%d +north=%s +datum=WGS84 +unit=m +no_defs'\n",param.zone,param.direction);
			    } else {
				fprintf(pMetafile, "Output Projection='+proj=utm +zone=%d +south=%s +datum=WGS84 +unit=m +no_defs'\n",param.zone,param.direction);
			    }
			}
		    }
		    
		    
		    bin_angle			= 360.0/18.0;

		    SetTiles(proinfo,proinfo.IsSP,proinfo.IsRR, Boundary, Res, tile_size, proinfo.pre_DEMtif, &pyramid_step, &buffer_area, 
			     &iter_row_start, &iter_row_end, &t_col_start, &t_col_end, &subX, &subY);	

                    total_count			= 0;

                    if (args.check_tiles_SR)
                    {
                        iter_row_start	  = args.start_row;
                    }
                    
                    if (args.check_tiles_ER)
                    {
                        iter_row_end	  = args.end_row;
                    }
                    
                    if (args.check_tiles_SC)
                    {	
                        t_col_start		  = args.start_col;
                    }
                    
                    if (args.check_tiles_EC)
                    {	
                        t_col_end		  = args.end_col;
                    }
                    
                    if (proinfo.check_tiles_SR)
                    {
                        iter_row_start	  = proinfo.start_row;
                    }
                    
                    if (proinfo.check_tiles_ER)
                    {
                        iter_row_end	  = proinfo.end_row;
                    }
                    
                    if (proinfo.check_tiles_SC)
                    {	
                        t_col_start		  = proinfo.start_col;
                    }
                    
                    if (proinfo.check_tiles_EC)
                    {	
                        t_col_end		  = proinfo.end_col;
                    }
                
                    printf("RA param = %f\t%f\n",Rimageparam[0],Rimageparam[1]);
                    
                    printf("Tiles row:col = row = %d\t%d\t;col = %d\t%d\tseed flag =%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end,proinfo.pre_DEMtif);
                    
                    if(!args.check_gridonly)
                    {
                        final_iteration = Matching_SETSM(proinfo,pyramid_step, Template_size, buffer_area,iter_row_start, iter_row_end,t_col_start,t_col_end,
                                                         subX,subY,bin_angle,Hinterval,Image_res,Res, Limageparam, Rimageparam,
                                                         LRPCs, RRPCs, pre_DEM_level, DEM_level,	NumOfIAparam, check_tile_array,Hemisphere,tile_array,
                                                         Limagesize,Rimagesize,LBRsize,RBRsize,param,total_count,ori_minmaxHeight,Boundary,1,1);
                    }
#ifdef BUILDMPI
					MPI_Barrier(MPI_COMM_WORLD);
					MPI_Finalize();
					if(rank != 0)
					{
						exit(0);
					}
#endif
		    if(!args.check_ortho)
		    {
			char check_file[500];
			FILE* pcheckFile;
			int max_row = 0;
			int max_col = 0;
			int row,col;
			for(row = 1; row < 100 ; row++)
			{
			    for(col = 1; col < 100 ; col++)
			    {
				sprintf(check_file,"%s/txt/matched_pts_%d_%d_0_3.txt",proinfo.save_filepath,row,col);
				pcheckFile = fopen(check_file,"r");
				if(pcheckFile)
				{
				    if(max_row < row)
					max_row = row;
				    if(max_col < col)
					max_col = col;
				}
			    }
			}
			
			if (args.check_tiles_SR)
			    iter_row_start	  = args.start_row;
			
			if (args.check_tiles_ER)
			    iter_row_end	  = args.end_row;
			else
			{
			    if (proinfo.check_tiles_ER)
				iter_row_end	  = proinfo.end_row;
			    else
			    {
				iter_row_end = max_row + 1;
			    }
			}
			
			if (args.check_tiles_SC)
			    t_col_start		  = args.start_col;
			
			if (args.check_tiles_EC)
			    t_col_end		  = args.end_col;
			else
			{
			    if (proinfo.check_tiles_EC)
				t_col_end		  = proinfo.end_col;
			    else
			    {
				t_col_end	 = max_col + 1;
			    }
			}
			
			if (proinfo.check_tiles_SR)
			    iter_row_start	  = proinfo.start_row;
			
			
			if (proinfo.check_tiles_SC)
			    t_col_start		  = proinfo.start_col;
			
		    }

		    printf("Tiles row:col = row = %d\t%d\t;col = %d\t%d\tseed flag =%d\n",iter_row_start,iter_row_end,t_col_start,t_col_end,proinfo.pre_DEMtif);
		    
		    
		    if(iter_row_end < 2 && t_col_end < 2)
		    {
			printf("No matching results. Please check overlapped area of stereo pair, or image textures\n");
			exit(1);
		    }
		    
		    char str_DEMfile[500];
		    sprintf(str_DEMfile, "%s/%s_dem.raw", proinfo.save_filepath,proinfo.Outputpath_name);
		    
		    FILE* pFile_DEM = NULL;
		    
		    pFile_DEM = fopen(str_DEMfile,"r");
		    printf("check exist %s %d\n",str_DEMfile,pFile_DEM);
		    final_iteration = 3;
		    if(!pFile_DEM)
		    {
			ST = time(0);
			printf("Tile merging start final iteration %d!!\n",final_iteration);
			int buffer_tile = 420;
			mt_grid_size = MergeTiles(proinfo,iter_row_start,t_col_start,iter_row_end,t_col_end,buffer_tile,final_iteration);
			
			mt_grid_size = proinfo.DEM_resolution;
			
			ET = time(0);
			gap = difftime(ET,ST);
			printf("Tile merging finish(time[m] = %5.2f)!!\n",gap/60.0);
			
			ST = time(0);
			printf("Interpolation start!!\n");
			sprintf(temp_c, "%s/%s_dem_tin.txt", proinfo.save_filepath,proinfo.Outputpath_name);
			printf("%f %f\n",proinfo.DEM_resolution,mt_grid_size);
			
			NNA_M(param,proinfo.save_filepath, proinfo.Outputpath_name,temp_c,iter_row_start,t_col_start, iter_row_end,t_col_end,proinfo.DEM_resolution,mt_grid_size,buffer_tile,Hemisphere,final_iteration);
			ET = time(0);
			gap = difftime(ET,ST);
			printf("Interpolation finish(time[m] = %5.2f)!!\n",gap/60.0);
		    }
		    
		    char hdr_path[500];
		    CSize seeddem_size;
		    double tminX, tmaxY, tgrid_size;
		    sprintf(hdr_path, "%s/%s_dem.hdr", proinfo.save_filepath, proinfo.Outputpath_name);
		    seeddem_size  = Envihdr_reader_seedDEM(param,hdr_path, &tminX, &tmaxY, &tgrid_size);
		    
		    fprintf(pMetafile,"Output dimensions=%d\t%d\n",seeddem_size.width,seeddem_size.height);
		    fprintf(pMetafile,"Upper left coordinates=%f\t%f\n",tminX,tmaxY);
		    
		    fprintf(pMetafile,"Image 1 info\nImage_1_satID=%s\nImage_1_Acquisition_time=%s\nImage_1_Mean_row_GSD=%f\nImage_1_Mean_col_GSD=%f\nImage_1_Mean_GSD=%f\nImage_1_Mean_sun_azimuth_angle=%f\nImage_1_Mean_sun_elevation=%f\nImage_1_Mean_sat_azimuth_angle=%f\nImage_1_Mean_sat_elevation=%f\nImage_1_Intrack_angle=%f\nImage_1_Crosstrack_angle=%f\nImage_1_Offnadir_angle=%f\nImage_1_tdi=%d\nImage_1_effbw=%f\nImage_1_abscalfact=%f\n",leftimage_info.SatID,leftimage_info.imagetime,Image1_gsd_r,Image1_gsd_c,Image1_gsd,leftimage_info.Mean_sun_azimuth_angle,leftimage_info.Mean_sun_elevation,leftimage_info.Mean_sat_azimuth_angle,leftimage_info.Mean_sat_elevation,leftimage_info.Intrack_angle,leftimage_info.Crosstrack_angle,leftimage_info.Offnadir_angle,(int)left_band.tdi,left_band.effbw,left_band.abscalfactor);
		    fprintf(pMetafile,"Image 2 info\nImage_2_satID=%s\nImage_2_Acquisition_time=%s\nImage_2_Mean_row_GSD=%f\nImage_2_Mean_col_GSD=%f\nImage_2_Mean_GSD=%f\nImage_2_Mean_sun_azimuth_angle=%f\nImage_2_Mean_sun_elevation=%f\nImage_2_Mean_sat_azimuth_angle=%f\nImage_2_Mean_sat_elevation=%f\nImage_2_Intrack_angle=%f\nImage_2_Crosstrack_angle=%f\nImage_2_Offnadir_angle=%f\nImage_2_tdi=%d\nImage_2_effbw=%f\nImage_2_abscalfact=%f\n",rightimage_info.SatID,rightimage_info.imagetime,Image2_gsd_r,Image2_gsd_c,Image2_gsd,rightimage_info.Mean_sun_azimuth_angle,rightimage_info.Mean_sun_elevation,rightimage_info.Mean_sat_azimuth_angle,rightimage_info.Mean_sat_elevation,rightimage_info.Intrack_angle,rightimage_info.Crosstrack_angle,rightimage_info.Offnadir_angle,(int)right_band.tdi,right_band.effbw,right_band.abscalfactor);
					} // if (!RA_only)
		    fclose(pMetafile);
                }
                else
                    printf("out of boundary!! please check boundary infomation!!\n");
                
                sprintf(temp_filepath,"%s/tmp",proinfo.save_filepath);
            }
		}
		else
		{
			printf("Check output directory path!!\n");
		}
	}
	else {
	}
	
	total_ET = time(0);
	total_gap = difftime(total_ET,total_ST);
	
	sprintf(computation_file,"%s/txt/computation_time.txt",proinfo.save_filepath);
	time_fid			= fopen(computation_file,"w");
	fprintf(time_fid,"Computation_time[m] = %5.2f\n",total_gap/60.0);
	fclose(time_fid);
#ifdef BUILDMPI
	// Make sure to finalize
	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized)
	{
		MPI_Finalize();
	}
#endif

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

int Matching_SETSM(ProInfo proinfo,uint8 pyramid_step, uint8 Template_size, uint16 buffer_area,uint8 iter_row_start, uint8 iter_row_end,uint8 t_col_start,uint8 t_col_end,
				   double subX,double subY,double bin_angle,double Hinterval,double *Image_res,double *Res, double *Limageparam, double *Rimageparam,
				   double **LRPCs, double **RRPCs, uint8 pre_DEM_level, uint8 DEM_level,	uint8 NumOfIAparam, bool check_tile_array,bool Hemisphere,bool* tile_array,
				   CSize Limagesize,CSize Rimagesize,CSize LBRsize,CSize RBRsize,TransParam param,int total_count,double *ori_minmaxHeight,double *Boundary, int row_iter, int col_iter)
{
#ifdef BUILDMPI
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
	int final_iteration = -1;
	bool lower_level_match;
	int row,col;
	int RA_count		= 0;

	int row_length = iter_row_end-iter_row_start;
	int col_length = t_col_end-t_col_start;
	int *iterations = (int*)malloc(col_length*row_length*2*sizeof(int));
	int length = 0;
	int count_tri;

	for(row = iter_row_start; row < iter_row_end ; row+=row_iter)
	{
		for(col = t_col_start ; col < t_col_end ; col+= col_iter)
		{
			iterations[2*length] = row;
			iterations[2*length+1] = col;
			length+=1;
		}
	}

#ifdef BUILDMPI
	//Reorder list of tiles for static load balancing
	if (length > 1) {
		reorder_list_of_tiles(proinfo, iterations, length, col_length, row_length);
	}
#endif

	int i;
	for(i = 0; i < length; i += 1)
	{
		row = iterations[2*i];
		col = iterations[2*i+1];

#ifdef BUILDMPI
		// Skip this tile if it belongs to a different MPI rank
		if (i % size != rank)
			continue;	
		printf("MPI: Rank %d is analyzing row %d, col %d\n", rank, row, col);
#endif

		char save_file[500], Lsubsetfilename[500], Rsubsetfilename[500];
		char *filename;
			
		FILE *fid = NULL;
		FILE *fid_header = NULL;
		FILE *fid_RAinfo = NULL;

		double minmaxHeight[2];
		double subBoundary[4];
		double t_Rimageparam[2] = {0.0};

		D2DPOINT Lstartpos_ori, Rstartpos_ori;

		CSize Lsubsetsize, Rsubsetsize;

		bool check_cal = false;
		if(proinfo.IsRA)
			check_cal = true;
		else {
			
			char check_file[500];
			sprintf(check_file,"%s/txt/matched_pts_%d_%d_0_3.txt",proinfo.save_filepath,row,col);
			FILE* pcheckFile;
			pcheckFile = fopen(check_file,"r");
			if(!pcheckFile)
				check_cal = true;
		}

		if(check_cal)
		{
			printf("start cal tile\n");
			total_count += 1;
			
			if(proinfo.IsRA)
			{
				if(!proinfo.check_checktiff)
				{
					sprintf(save_file,"%s/txt/RA_echo_result_row_%d_col_%d.txt",proinfo.save_filepath,row,col);
					fid			= fopen(save_file,"w");
					fprintf(fid,"RA param X = %f\tY = %f\n",t_Rimageparam[0],t_Rimageparam[1]);

					sprintf(save_file,"%s/txt/RA_headerinfo_row_%d_col_%d.txt",proinfo.save_filepath,row,col);
					fid_header	= fopen(save_file,"w");
				}
			}
			else
			{
				t_Rimageparam[0]	= Rimageparam[0];
				t_Rimageparam[1]	= Rimageparam[1];
				
				if(!proinfo.check_checktiff)
				{
					sprintf(save_file,"%s/txt/echo_result_row_%d_col_%d.txt",proinfo.save_filepath,row,col);
					fid			= fopen(save_file,"w");
					fprintf(fid,"RA param X = %f\tY = %f\n",t_Rimageparam[0],t_Rimageparam[1]);

					sprintf(save_file,"%s/txt/headerinfo_row_%d_col_%d.txt",proinfo.save_filepath,row,col);
					fid_header	= fopen(save_file,"w");
				}
			}

			minmaxHeight[0]		= ori_minmaxHeight[0];
			minmaxHeight[1]		= ori_minmaxHeight[1];
			//minmaxHeight[0]		-= 100;
			if(minmaxHeight[0] < -100)
				minmaxHeight[0] = -100;

			printf("minmaxH = %f\t%f\n",minmaxHeight[0],minmaxHeight[1]);
			SetSubBoundary(Boundary,subX,subY,buffer_area,col,row,subBoundary);

			printf("subBoundary = %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]);
			filename = GetFileName(proinfo.LeftImagefilename);
			filename = remove_ext(filename);
			sprintf(Lsubsetfilename,"%s/%s_subset_%d_%d.raw",proinfo.tmpdir,filename,row,col);
			filename = GetFileName(proinfo.RightImagefilename);
			filename = remove_ext(filename);
			sprintf(Rsubsetfilename,"%s/%s_subset_%d_%d.raw",proinfo.tmpdir,filename,row,col);

			printf("subsetimage\n");
			
			if(subsetImage(param,NumOfIAparam,LRPCs,Limageparam,proinfo.LeftImagefilename,RRPCs,t_Rimageparam,proinfo.RightImagefilename,subBoundary,minmaxHeight,
						   &Lstartpos_ori,&Rstartpos_ori,Lsubsetfilename,Rsubsetfilename,&Lsubsetsize, &Rsubsetsize, fid,proinfo.check_checktiff))
			{
				bool check_tile_pr = true;
				if(check_tile_array)
				{
					if(proinfo.NumOfTile_row >= row && proinfo.NumOfTile_col >= col)
					{
						if(tile_array[proinfo.NumOfTile_col*row + col])
							check_tile_pr	= true;
						else
							check_tile_pr	= false;
					}
					else
						check_tile_pr		= false;
				}

				printf("Completion of subsetImage!!\n");

				if( Lsubsetsize.height > Template_size/2*pow(2,pyramid_step) && Lsubsetsize.width > Template_size/2*pow(2,pyramid_step) &&
					Rsubsetsize.height > Template_size/2*pow(2,pyramid_step) && Rsubsetsize.width > Template_size/2*pow(2,pyramid_step) &&
					check_tile_pr)
				{
					double pre_3sigma= 1000;
					double pre_mean	= 1000;
					double ratio = 0;
					double py_resolution = 0;
					double grid_resolution = 0;
					double pre_grid_resolution = 0;
					
					bool flag_start, dem_update_flag;

					int level			  = pyramid_step;
					uint8 sub_total_count	= 0;
					
					CSize Size_Grid2D, pre_Size_Grid2D;

					CSize *data_size_l, *data_size_r;
					UGRID *GridPT3 = NULL, *Pre_GridPT3 = NULL;
					
					
					time_t PreST = 0, PreET = 0;
					double Pregab;
					
					lower_level_match	= true;
					flag_start			= false;
					if(proinfo.IsRA)
						dem_update_flag		= false;
					else
						dem_update_flag		= true;

					
					Size_Grid2D.height	= 0;
					Size_Grid2D.width	= 0;
					data_size_l = (CSize*)malloc(sizeof(CSize)*(level+1));
					data_size_r = (CSize*)malloc(sizeof(CSize)*(level+1));

					sub_total_count = CalTotalIteration(DEM_level,level);
					SetPySizes(data_size_l, data_size_r, Lsubsetsize, Rsubsetsize, level);

					PreST = time(0);
					printf("row = %d/%d\tcol = %d/%d\tPreprocessing start!!\n",row,iter_row_end,col,t_col_end);
					Preprocessing(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,level,&Lsubsetsize, &Rsubsetsize,data_size_l,data_size_r, fid);
					PreET = time(0);
					Pregab = difftime(PreET,PreST);
					printf("row = %d/%d\tcol = %d/%d\tPreprocessing finish(time[m] = %5.2f)!!\n",row,iter_row_end,col,t_col_end,Pregab/60.0);
					
					PreST = time(0);
					printf("row = %d/%d\tcol = %d/%d\tDEM generation start!!\n",row,iter_row_end,col,t_col_end);
					
					int blunder_selected_level;
					int pre_blunder_selected_level;
					if(proinfo.DEM_resolution == 8)
						pre_blunder_selected_level = 4;
					else if(proinfo.DEM_resolution == 4)
						pre_blunder_selected_level = 3;
					else if(proinfo.DEM_resolution == 2)
						pre_blunder_selected_level = 2;
					else if(proinfo.DEM_resolution == 1)
						pre_blunder_selected_level = 1;
					else
						pre_blunder_selected_level = 0;
					
					D2DPOINT *GridPT = NULL;
                    int final_level_iteration = 1;
                        
                    int total_matching_candidate_pts = 0;
                    double matching_rate = 0;
                        
                    bool check_RA_divide = false;
                    double tilesize_RA = 20000;
                    double lengthOfX = subBoundary[2] - subBoundary[0];
                    double lengthOfY = subBoundary[3] - subBoundary[1];
                    int division_X = 0, division_Y = 0;
                    int total_tile = 0;
                    double new_subBoundary_RA[4];
                    bool check_new_subBoundary_RA = false;
                    double preBoundary[4] = {0};
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
                 
                    printf("length %f\t%f\tdivision %d\t%d\ntotal_tile %d\tcheck_RA_divide %d\n",lengthOfX,lengthOfY,division_X,division_Y,total_tile,check_RA_divide);
                    
                  	while(lower_level_match && level >= DEM_level)
					{
                    	printf("level = %d\t final_level_iteration %d\n",level,final_level_iteration);
                           
                        if(proinfo.IsRA && check_new_subBoundary_RA)
                        {
                            preBoundary[0] = subBoundary[0];
                            preBoundary[1] = subBoundary[1];
                            preBoundary[2] = subBoundary[2];
                            preBoundary[3] = subBoundary[3];
                            
                            subBoundary[0] = new_subBoundary_RA[0];
                            subBoundary[1] = new_subBoundary_RA[1];
                            subBoundary[2] = new_subBoundary_RA[2];
                            subBoundary[3] = new_subBoundary_RA[3];
                            
                            RemoveFiles(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,0,0);
                            
                            subsetImage(param,NumOfIAparam,LRPCs,Limageparam,proinfo.LeftImagefilename,RRPCs,t_Rimageparam,proinfo.RightImagefilename,subBoundary,minmaxHeight,
                                        &Lstartpos_ori,&Rstartpos_ori,Lsubsetfilename,Rsubsetfilename,&Lsubsetsize, &Rsubsetsize, fid,proinfo.check_checktiff);
                            
                            SetPySizes(data_size_l, data_size_r, Lsubsetsize, Rsubsetsize, level+1);
                            
                            Preprocessing(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,level+1,&Lsubsetsize, &Rsubsetsize,data_size_l,data_size_r, fid);
                        }
                        
                        printf("subBoundary %f\t%f\t%f\t%f\t preBoundary %f\t%f\t%f\t%f\n",subBoundary[0],subBoundary[1],subBoundary[2],subBoundary[3],preBoundary[0],preBoundary[1],preBoundary[2],preBoundary[3]);
                        
						double Th_roh, Th_roh_min, Th_roh_start, Th_roh_next;
						double minH_mps, maxH_mps;
						double minH_grid, maxH_grid;
						double MPP;
						
						uint8 iteration;

						D2DPOINT Lstartpos, Rstartpos;
						D2DPOINT BLstartpos, BRstartpos;
						
						uint16 *SubImages_L, *SubImages_R;
						uint8  *SubOriImages_L, *SubOriImages_R;
						uint16 *SubMagImages_L, *SubMagImages_R;
						uint16 *SubImages_BL, *SubImages_BR;
						uint16 *SubMagImages_BL, *SubMagImages_BR;
						
						D2DPOINT *Grid_wgs;
						
						if(level == 4)
							blunder_selected_level = level;
						else
							blunder_selected_level = level + 1;
							
						printf("selected_bl %d\n",blunder_selected_level);
						
						SetThs(level,final_level_iteration, &Th_roh, &Th_roh_min, &Th_roh_next, &Th_roh_start,proinfo.pre_DEMtif,proinfo.IsRA,proinfo.seedDEMsigma,proinfo.DEM_resolution);
							
						Lstartpos.m_X		= (double)(Lstartpos_ori.m_X/pow(2,level));		Lstartpos.m_Y		= (double)(Lstartpos_ori.m_Y/pow(2,level));
						Rstartpos.m_X		= (double)(Rstartpos_ori.m_X/pow(2,level));		Rstartpos.m_Y		= (double)(Rstartpos_ori.m_Y/pow(2,level));

						BLstartpos.m_X		 = (double)(Lstartpos_ori.m_X/pow(2,blunder_selected_level));		BLstartpos.m_Y		 = (double)(Lstartpos_ori.m_Y/pow(2,blunder_selected_level));
						BRstartpos.m_X		 = (double)(Rstartpos_ori.m_X/pow(2,blunder_selected_level));		BRstartpos.m_Y		 = (double)(Rstartpos_ori.m_Y/pow(2,blunder_selected_level));
						
						if(proinfo.IsRA)
						{
							dem_update_flag			= false;
							py_resolution			= Image_res[0]*pow(2,pyramid_step+1);
							grid_resolution			= Image_res[0]*pow(2,pyramid_step+1);
                               
                            printf("RA grid size %f\n",py_resolution);
                               
							//if(!flag_start)
							{
								GridPT					= SetDEMGrid(subBoundary, grid_resolution, grid_resolution,&Size_Grid2D);
                                   
                                printf("Size_Grid2D %d\t%d\n",Size_Grid2D.width,Size_Grid2D.height);
							}
						}
						else
							GridPT	= SetGrids(&dem_update_flag, flag_start, level, final_level_iteration, proinfo.resolution, &Size_Grid2D, proinfo.pre_DEMtif,proinfo.priori_DEM_tif, proinfo.DEM_resolution, minmaxHeight,&py_resolution, &grid_resolution, subBoundary);

						if(!flag_start)
						{
							printf("GridPT3 start\t seed flag %d\t filename %s\timage_resolution %f minmax %f %f\n",proinfo.pre_DEMtif,proinfo.priori_DEM_tif,Image_res[0],minmaxHeight[0],minmaxHeight[1]);
							GridPT3 = SetGrid3PT(param, dem_update_flag, flag_start, Size_Grid2D, Th_roh, level, minmaxHeight,subBoundary,grid_resolution,proinfo.priori_DEM_tif,proinfo.pre_DEMtif,proinfo.seedDEMsigma,proinfo.IsRA,proinfo.metafilename);
						}
						
						if(flag_start)
						{
                            if(proinfo.IsRA)
                            {
                                if(check_new_subBoundary_RA)
                                {
                                    GridPT3 = ResizeGirdPT3_RA(pre_Size_Grid2D, Size_Grid2D, preBoundary,subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                                
                                    check_new_subBoundary_RA = false;
                                    
                                    printf("start ResizeGridPT3 with newBoundry pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                }
                                else
                                {
                                    printf("start ResizeGridPT3 pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                    GridPT3 = ResizeGirdPT3(pre_Size_Grid2D, Size_Grid2D, subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                                }
                            }
                            else
                            {
                                printf("start ResizeGridPT3 pre size %d %d size %d %d pre_resol %f\n",pre_Size_Grid2D.width,pre_Size_Grid2D.height,Size_Grid2D.width,Size_Grid2D.height,pre_grid_resolution);
                                GridPT3 = ResizeGirdPT3(pre_Size_Grid2D, Size_Grid2D, subBoundary, GridPT, Pre_GridPT3, pre_grid_resolution,minmaxHeight);
                            }
                		}
			
                        printf("end start ResizeGridPT3\n");
                        
						pre_Size_Grid2D.width = Size_Grid2D.width;
						pre_Size_Grid2D.height = Size_Grid2D.height;
						pre_grid_resolution = grid_resolution;
						
						fprintf(fid,"level = %d, Completion of Gridinfo setup\t%d\t%d!!\n",level,Size_Grid2D.width,Size_Grid2D.height);
						
						fprintf(fid_header, "%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n", row, col, level, subBoundary[0], subBoundary[1], grid_resolution, Size_Grid2D.width,Size_Grid2D.height);

						double left_mag_var, left_mag_avg, right_mag_var, right_mag_avg;
						double left_mag_var_B, left_mag_avg_B, right_mag_var_B, right_mag_avg_B;
						
                        printf("load subimages\n");
                        
						SubImages_L		= LoadPyramidImages(proinfo.tmpdir,Lsubsetfilename,data_size_l[level],level);
						SubImages_R		= LoadPyramidImages(proinfo.tmpdir,Rsubsetfilename,data_size_r[level],level);
						SubOriImages_L	= LoadPyramidOriImages(proinfo.tmpdir,Lsubsetfilename,data_size_l[level],level);
						SubOriImages_R	= LoadPyramidOriImages(proinfo.tmpdir,Rsubsetfilename,data_size_r[level],level);
						SubMagImages_L	= LoadPyramidMagImages(proinfo.tmpdir,Lsubsetfilename,data_size_l[level],level,&left_mag_var,&left_mag_avg);
						SubMagImages_R	= LoadPyramidMagImages(proinfo.tmpdir,Rsubsetfilename,data_size_r[level],level,&right_mag_var,&right_mag_avg);

						printf("mag var avg %f\t%f\t%f\t%f\n",left_mag_var,left_mag_avg,right_mag_var,right_mag_avg);
						
						double mag_avg = (left_mag_avg + right_mag_avg)/2.0;
						double mag_var = (left_mag_var + right_mag_var)/2.0;
						
						SubImages_BL		= LoadPyramidImages(proinfo.tmpdir,Lsubsetfilename,data_size_l[blunder_selected_level],blunder_selected_level);
						SubImages_BR		= LoadPyramidImages(proinfo.tmpdir,Rsubsetfilename,data_size_r[blunder_selected_level],blunder_selected_level);
						SubMagImages_BL		= LoadPyramidMagImages(proinfo.tmpdir,Lsubsetfilename,data_size_l[blunder_selected_level],blunder_selected_level,&left_mag_var_B,&left_mag_avg_B);
						SubMagImages_BR		= LoadPyramidMagImages(proinfo.tmpdir,Rsubsetfilename,data_size_r[blunder_selected_level],blunder_selected_level,&right_mag_var_B,&right_mag_avg_B);

						printf("load subimages blunder_selected_level\n");
                           
                        total_matching_candidate_pts = Size_Grid2D.width*Size_Grid2D.height;
                           
						Grid_wgs = ps2wgs(param,Size_Grid2D.width*Size_Grid2D.height,GridPT);
						
						if(proinfo.pre_DEMtif && !flag_start)
						{
							ratio = VerticalLineLocus_seeddem(SubMagImages_L,SubMagImages_R,grid_resolution, Image_res[0], LRPCs, RRPCs,
															  Limagesize, data_size_l[level], SubImages_L, Rimagesize, data_size_r[level], SubImages_R, Template_size, 
															  Size_Grid2D, param, GridPT, Grid_wgs, GridPT3,
															  NumOfIAparam, t_Rimageparam, level, Lstartpos, Rstartpos, 
															  proinfo.save_filepath, row, col, 1,1,subBoundary,minmaxHeight,proinfo.seedDEMsigma);
							printf("ratio %f\n",ratio);
						}
							
						if(ratio > 70)
						{
							SetThs_ratio(level, &Th_roh, &Th_roh_min, &Th_roh_next, &Th_roh_start,proinfo.pre_DEMtif,proinfo.IsRA,proinfo.DEM_resolution);
							
							printf("%f \t %f\n",Th_roh,Th_roh_min);
						}
						
						iteration		= 1;
                        if(level == 0)
                            iteration = final_level_iteration;
                        
						int pre_matched_pts=10;
						double matching_change_rate = 100;
						double rate_th = 0.00999999;
						int max_iteration = 10 - (4 - level)*2;
						if(max_iteration < 3)
							max_iteration = 3;
						
						NCCresult *nccresult;
						nccresult = (NCCresult*)calloc(sizeof(NCCresult),Size_Grid2D.width*Size_Grid2D.height);
						
                        if(level == 0 &&  iteration == 3)
                            matching_change_rate = 0.001;
                        
                         while((Th_roh >= Th_roh_min || (matching_change_rate > rate_th)) )
						{
							printf("%f \t %f\n",Th_roh,Th_roh_min);
							
							double pre_3sigma, pre_mean;
							double Th_roh_update = 0;

							NCCflag flag;

							BL blunder_param;

							char filename_mps[500];
							char filename_mps_mid[500];
							char filename_mps_fin[500];
							char filename_mps_pre[500];
							char filename_mps_anchor[500];
							char filename_mps_aft[500];
								
							char filename_tri[500];
							char v_temp_path[500];
								
							int count_results[2];
							int count_results_anchor[2];
							int count_MPs;
							int count_blunder;
							bool update_flag = false;
							bool check_ortho_cal = false;
								
							uint8 ortho_level = 2;
							if(proinfo.DEM_resolution >= 8)
								ortho_level = 3;
							
							if(level >= ortho_level)
							{
								check_ortho_cal = true;
							}
							else
								check_ortho_cal = false;
								
							printf("ortho level = %d\n",ortho_level);
							
							pre_3sigma	= 0;	pre_mean	= 0;
							
							if(level > 1)
							{
								flag.rotate_flag	 = 1;		flag.multi_flag		 = 1;		flag.multi_flag_sum	 = 1;		flag.inter_flag		 = 1;
							}
							else
							{
								flag.rotate_flag	 = 1;		flag.multi_flag		 = 1;		flag.multi_flag_sum	 = 1;		flag.inter_flag		 = 1;		  
							}

							if(iteration > 1)
								flag.weight_flag	 = 1;
							else
								flag.weight_flag	 = 0;

							fprintf(fid,"Starting computation of NCC\n iteration = %u\tTh_roh = %f\tTh_roh_start = %f\tGrid size %d %d\n",
									iteration, Th_roh,Th_roh_start,Size_Grid2D.width,Size_Grid2D.height);

							printf("sub size %d\t%d\t%d\t%d\n",data_size_l[level].width,data_size_l[level].height,data_size_r[level].width,data_size_r[level].height);
								
							if(proinfo.IsRA)
							{
								sprintf(filename_mps,"%s/txt/RA_matched_pts_%d_%d_%d_%d.txt",proinfo.save_filepath,row,col,level,iteration);
								sprintf(filename_mps_pre,"%s/txt/matched_pts_%d_%d_%d_%d_pre.txt",proinfo.save_filepath,row,col,level,iteration);
							}
							else
							{
								sprintf(filename_mps,"%s/txt/matched_pts_%d_%d_%d_%d.txt",proinfo.save_filepath,row,col,level,iteration);
								sprintf(filename_mps_fin,"%s/txt/matched_pts_%d_%d_%d_%d_fin.txt",proinfo.save_filepath,row,col,level,iteration);
								sprintf(filename_mps_pre,"%s/txt/matched_pts_%d_%d_%d_%d_pre.txt",proinfo.save_filepath,row,col,level,iteration);
								sprintf(filename_mps_aft,"%s/txt/matched_pts_%d_%d_%d_%d_aft.txt",proinfo.save_filepath,row,col,level,iteration);
								sprintf(filename_mps_anchor,"%s/txt/matched_pts_%d_%d_%d_%d_anchor.txt",proinfo.save_filepath,row,col,level,iteration);
							}
							
							sprintf(v_temp_path,"%s/tmp/vv_tmp",proinfo.save_filepath);
							
							printf("template size =%d\n",Template_size);
							
							VerticalLineLocus(nccresult,SubMagImages_L,SubMagImages_R,grid_resolution, Image_res[0],LRPCs,RRPCs,Limagesize,data_size_l[level],SubImages_L,Rimagesize,data_size_r[level],SubImages_R,Template_size,Size_Grid2D,
											  param,GridPT,Grid_wgs,GridPT3,flag,NumOfIAparam,t_Rimageparam,minmaxHeight,level,Lstartpos,Rstartpos,iteration,SubOriImages_L,SubOriImages_R,bin_angle,1,0,fid,true,Hemisphere,
											  proinfo.save_filepath,row,col,subBoundary,proinfo.pre_DEMtif,v_temp_path,&MPP,proinfo.IsRA,mag_avg,mag_var);
							
							printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd computation of NCC!! minmax %f %f\n",row,col,level,iteration,minmaxHeight[0], minmaxHeight[1]);

							minH_mps = 1000000;
							maxH_mps = -1000000;

							count_MPs = SelectMPs(nccresult,Size_Grid2D,GridPT,GridPT3,Th_roh,Th_roh_min,Th_roh_start,Th_roh_next,level,pyramid_step,
												  iteration,0,filename_mps_pre,proinfo.pre_DEMtif,proinfo.IsRA,MPP,proinfo.DEM_resolution,Image_res[0],final_level_iteration);
							printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SelectMPs\tcount_mps = %d\n",row,col,level,iteration,count_MPs);
								
							if (check_ortho_cal && proinfo.IsRA != 1)
							{
								//anchor points
								count_blunder = DecisionMPs(false,count_MPs,subBoundary,GridPT3,level,grid_resolution,iteration,Size_Grid2D,filename_mps_pre,filename_mps_anchor,proinfo.save_filepath,
															Hinterval,&lower_level_match,&pre_3sigma,&pre_mean,count_results_anchor,&minH_mps,&maxH_mps,minmaxHeight,
															SubMagImages_BL,SubMagImages_BR,grid_resolution, Image_res[0],LRPCs,RRPCs,
															Limagesize,data_size_l[blunder_selected_level],SubImages_BL,Rimagesize,data_size_r[blunder_selected_level],SubImages_BR,Template_size,
															param, Grid_wgs,GridPT,
															NumOfIAparam, t_Rimageparam, BLstartpos, BRstartpos,
															proinfo.save_filepath,row,col,proinfo.pre_DEMtif,proinfo.IsRA,SubOriImages_L,SubOriImages_R,blunder_selected_level,proinfo.seedDEMsigma,proinfo.DEM_resolution);
								printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd anchor points\n",row,col,level,iteration);

								//blunder detection
								count_blunder = DecisionMPs(true,count_MPs,subBoundary,GridPT3,level,grid_resolution,iteration,Size_Grid2D,filename_mps_pre,filename_mps_aft,proinfo.save_filepath,
															Hinterval,&lower_level_match,&pre_3sigma,&pre_mean,count_results,&minH_mps,&maxH_mps,minmaxHeight,
															SubMagImages_BL,SubMagImages_BR,grid_resolution, Image_res[0],LRPCs,RRPCs,
															Limagesize,data_size_l[blunder_selected_level],SubImages_BL,Rimagesize,data_size_r[blunder_selected_level],SubImages_BR,Template_size,
															param, Grid_wgs,GridPT,
															NumOfIAparam, t_Rimageparam, BLstartpos, BRstartpos,
															proinfo.save_filepath,row,col,proinfo.pre_DEMtif,proinfo.IsRA,SubOriImages_L,SubOriImages_R,blunder_selected_level,proinfo.seedDEMsigma,proinfo.DEM_resolution);
								
								printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
							}
							else
							{
								count_blunder = DecisionMPs(true,count_MPs,subBoundary,GridPT3,level,grid_resolution,iteration,Size_Grid2D,filename_mps_pre,filename_mps,proinfo.save_filepath,
															Hinterval,&lower_level_match,&pre_3sigma,&pre_mean,count_results,&minH_mps,&maxH_mps,minmaxHeight,
															SubMagImages_BL,SubMagImages_BR,grid_resolution, Image_res[0],LRPCs,RRPCs,
															Limagesize,data_size_l[blunder_selected_level],SubImages_BL,Rimagesize,data_size_r[blunder_selected_level],SubImages_BR,Template_size,
															param, Grid_wgs,GridPT,
															NumOfIAparam, t_Rimageparam, BLstartpos, BRstartpos,
															proinfo.save_filepath,row,col,proinfo.pre_DEMtif,proinfo.IsRA,SubOriImages_L,SubOriImages_R,blunder_selected_level,proinfo.seedDEMsigma,proinfo.DEM_resolution);
								count_MPs		= count_results[0];
								printf("RA row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd blunder points\n",row,col,level,iteration);
							}
							
							fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tcheck = %d(%d,%d)\tEnd blunder detection\n",row,col,level,iteration,lower_level_match,count_MPs,count_blunder);

							printf("End computation of blunder!! Mps = %d\tTris = %d\tminz Mp = %f\tmaxz Mp = %f minmax %f %f \n",
								   count_results[0],count_results[1],minH_mps,maxH_mps,minmaxHeight[0],minmaxHeight[1]);
								
							if(lower_level_match)
							{
								FILE* survey;
								int i;
								
								if(check_ortho_cal && proinfo.IsRA != 1)
								{
									count_MPs = SetttingFlagOfGrid(subBoundary,GridPT3,level,grid_resolution,iteration,Size_Grid2D,filename_mps_anchor,filename_mps_aft,count_results_anchor[0],count_results[0],filename_mps_fin);
									survey	= fopen(filename_mps_fin,"r");
								}
								else
									survey	= fopen(filename_mps,"r");
								
								blunder_param.Boundary	= subBoundary;
								blunder_param.gridspace	= grid_resolution;
								blunder_param.height_check_flag = true;
								blunder_param.Hinterval = Hinterval;
								blunder_param.iteration = iteration;
								blunder_param.Pyramid_step = level;
								
								blunder_param.Size_Grid2D.width = Size_Grid2D.width;
								blunder_param.Size_Grid2D.height = Size_Grid2D.height;
																	
								int TIN_split_level = 0;
								if (grid_resolution <= 8)
								{
									if(grid_resolution <= 8 && grid_resolution > 4)
										TIN_split_level = 0;
									else if(grid_resolution <= 4)
										TIN_split_level = 2;
								}
								
								if(proinfo.IsRA)
								{
									TIN_split_level = 4;
								}
								
								if(level == 0 && iteration == 3)
								{
									D3DPOINT *ptslists;
									ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
									
									i = 0;
									while( i < count_MPs && (fscanf(survey,"%lf %lf %lf %hhd\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z,&ptslists[i].flag)) != EOF )
									{
										i++;
									}
									fclose(survey);
									
									FILE *pFile = fopen(filename_mps,"w");
									for(i=0;i<count_MPs;i++)
									{
										if(ptslists[i].flag != 1)
										{
											if(ptslists[i].m_X >= subBoundary[0] && ptslists[i].m_X <= subBoundary[2] && ptslists[i].m_Y >= subBoundary[1] && ptslists[i].m_Y <= subBoundary[3])
												fprintf(pFile,"%f %f %f\n",ptslists[i].m_X,ptslists[i].m_Y,ptslists[i].m_Z);
										}
										//i++;
									}
									fclose(pFile);
									
									if(!proinfo.IsRA)
										echoprint_Gridinfo(proinfo.save_filepath,row,col,level,iteration,update_flag,&Size_Grid2D,GridPT3,"final");
									
									free(ptslists);
								}
								else
								{
									if(check_ortho_cal && proinfo.IsRA != 1)
									{
										char bufstr[500];
										D3DPOINT *ptslists;
										
										FILE *pTri;
										double maxX_ptslists = -100000000;
										double maxY_ptslists = -100000000;
										double minX_ptslists =	100000000;
										double minY_ptslists =	100000000;
										int i;
										
										ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
										
										i = 0;
										while( i < count_MPs && (fscanf(survey,"%lf %lf %lf %hhd\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z,&ptslists[i].flag)) != EOF )
										{
											if(maxX_ptslists < ptslists[i].m_X)
												maxX_ptslists = ptslists[i].m_X;
											if(maxY_ptslists < ptslists[i].m_Y)
												maxY_ptslists = ptslists[i].m_Y;
											if(minX_ptslists > ptslists[i].m_X)
												minX_ptslists = ptslists[i].m_X;
											if(minY_ptslists > ptslists[i].m_Y)
												minY_ptslists = ptslists[i].m_Y;
											i++;
										}
										
										fclose(survey);
										double min_max[4] = {minX_ptslists, minY_ptslists, maxX_ptslists, maxY_ptslists};
										UI3DPOINT *trilists;
										
										if(level >= TIN_split_level || count_MPs < 10000)
										{
											UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_MPs*4);
											
											sprintf(bufstr,"%s/txt/tri_ortho.txt",proinfo.save_filepath);
											TINCreate(ptslists,bufstr,count_MPs,t_trilists,min_max,&count_tri);
											
											trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
											i = 0;
											for(i=0;i<count_tri;i++)
											{
												trilists[i].m_X = t_trilists[i].m_X;
												trilists[i].m_Y = t_trilists[i].m_Y;
												trilists[i].m_Z = t_trilists[i].m_Z;
											}
											
											free(t_trilists);
										
										}
										else
										{
											int iter_row,iter_col;
											
											trilists = TINgeneration(true,proinfo.save_filepath, level, Size_Grid2D, Image_res[0], grid_resolution,
																	 min_max,
																	 subBoundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
										}
										
									
										fprintf(fid,"level = %d\tMatching Pts = %d\n",level,count_results[0]);
										
										printf("ortho minmax %f %f pts anchor blunder %d %d \n",minmaxHeight[0],minmaxHeight[1],count_MPs,count_tri);
										
										count_results[0] = Ortho_blunder(ptslists, count_MPs, trilists,count_tri, update_flag,&minH_grid,&maxH_grid,blunder_param,
																		 SubMagImages_L,SubMagImages_R,SubImages_L,SubImages_R,
																		 grid_resolution, Image_res[0],LRPCs,RRPCs,
																		 data_size_l[level],data_size_r[level],Size_Grid2D,param,NumOfIAparam,
																		 t_Rimageparam,minmaxHeight,level,MPP,
																		 Lstartpos,Rstartpos,iteration,GridPT3,filename_mps,proinfo.save_filepath);
										free(trilists);
										
										printf("end ortho_blunder %d\n",count_results[0]);
										
										int matched_pts = 0;
										i = 0;
										{
											FILE *fid_all = fopen(filename_mps,"w");
											for(i=0;i<count_MPs;i++)
											{
												int t_col, t_row,grid_index;
												t_col		  = (int)((ptslists[i].m_X - subBoundary[0])/grid_resolution + 0.5);
												t_row		  = (int)((ptslists[i].m_Y - subBoundary[1])/grid_resolution + 0.5);
												grid_index	   = Size_Grid2D.width*t_row + t_col;
												
												fprintf(fid_all,"%f %f %f\n",ptslists[i].m_X,ptslists[i].m_Y,ptslists[i].m_Z);
												i++;
												matched_pts++;
											}
											fclose(fid_all);
										}
										
										free(ptslists);
											
										printf("load ortho_blunder pts %d\n",matched_pts);
										count_MPs = matched_pts;
										
										survey	= fopen(filename_mps,"r");
										ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
										
										i = 0;
										while( i < count_MPs && (fscanf(survey,"%lf %lf %lf\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z)) != EOF )
										{
											if(maxX_ptslists < ptslists[i].m_X)
												maxX_ptslists = ptslists[i].m_X;
											if(maxY_ptslists < ptslists[i].m_Y)
												maxY_ptslists = ptslists[i].m_Y;
											if(minX_ptslists > ptslists[i].m_X)
												minX_ptslists = ptslists[i].m_X;
											if(minY_ptslists > ptslists[i].m_Y)
												minY_ptslists = ptslists[i].m_Y;
											
											i++;
										}
										fclose(survey);
										double min_max2[4] = {minX_ptslists,minY_ptslists,maxX_ptslists,maxY_ptslists};
										
										if(level >= TIN_split_level || count_MPs < 10000)
										{
											UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_MPs*4);
											
											sprintf(bufstr,"%s/txt/tri_ortho.txt",proinfo.save_filepath);
											TINCreate(ptslists,bufstr,count_MPs,t_trilists,min_max2,&count_tri);
											
											trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
											i = 0;
											for(i=0;i<count_tri;i++)
											{
												trilists[i].m_X = t_trilists[i].m_X;
												trilists[i].m_Y = t_trilists[i].m_Y;
												trilists[i].m_Z = t_trilists[i].m_Z;
											}
											
											free(t_trilists);
												
										}
										else
										{
											int iter_row,iter_col;
											
											trilists = TINgeneration(true,proinfo.save_filepath, level, Size_Grid2D, Image_res[0], grid_resolution,
																	 min_max2,
																	 subBoundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
										}
										
										
										count_blunder = DecisionMPs_setheight(true,count_MPs,subBoundary,GridPT3,level,grid_resolution,iteration,Size_Grid2D,filename_mps,proinfo.save_filepath,
																			  Hinterval,&lower_level_match,&pre_3sigma,&pre_mean,count_results,&minH_mps,&maxH_mps,minmaxHeight,
																			  SubMagImages_BL,SubMagImages_BR,grid_resolution, Image_res[0],LRPCs,RRPCs,
																			  Limagesize,data_size_l[blunder_selected_level],SubImages_BL,Rimagesize,data_size_r[blunder_selected_level],SubImages_BR,Template_size,
																			  param, Grid_wgs,GridPT,
																			  NumOfIAparam, t_Rimageparam, BLstartpos, BRstartpos,
																			  proinfo.save_filepath,row,col,ptslists,trilists,count_tri,SubOriImages_L,SubOriImages_R,blunder_selected_level);

										printf("end decision_setheight\n");
										
                                        if(pre_matched_pts == 0)
                                            matching_change_rate = 0;
                                        else
                                            matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
											
										printf("matching change rate pre curr %f\t%d\t%d\n",matching_change_rate,count_MPs,pre_matched_pts);
										pre_matched_pts = count_MPs;
										
										if(iteration > 9 || level == 0 )
											matching_change_rate = 0.001;
										
										if(level <= 1)
                                        {
                                            if(iteration > 5)
                                                matching_change_rate = 0.001;
                                        }
                                        
										if(proinfo.IsRA)
											matching_change_rate = 0.001;
										
										if(proinfo.DEM_resolution >= 8)
										{
											if(level <= 2)
												matching_change_rate = 0.001;
											//if(level == 3 && iteration > 4)
											//	matching_change_rate = 0.001;
											
											matching_change_rate = 0.001;
										}
										
										if(proinfo.pre_DEMtif)
										{
											if(level >= 4)
												matching_change_rate = 0.001;
										}
										
										if(Th_roh >= Th_roh_min)
										{
											if(level == 0)
												Th_roh_update		= (double)(Th_roh - 0.50);
											else if(level == 1)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else if(level == 2)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else if(level == 3)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else
											{
												if(proinfo.IsRA)
													Th_roh_update		= (double)(Th_roh - 0.10);
												else
													Th_roh_update		= (double)(Th_roh - 0.06);
												
											}
										}
										
										bool check_level_end = false;
										if(level != 0)
										{
											if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th)
											{
												check_level_end = true;
												
												if(dem_update_flag)
													update_flag		= true;
												else
													update_flag		= false;
											}
											else
												update_flag			= false;
										}
										else
											update_flag				= false;
										
										if(level == 0)
                                        {
                                            Pre_GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
                                            printf("update GridPT3\n");
                                        }
                                        else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
										{
											Pre_GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
											printf("update GridPT3\n");
										}
										else
											GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
										
										free(trilists);
										free(ptslists);
										
										final_iteration = iteration;
									}
									else
									{
										char bufstr[500];
										D3DPOINT *ptslists;
										
										FILE *pTri;
										double maxX_ptslists = -100000000;
										double maxY_ptslists = -100000000;
										double minX_ptslists =	100000000;
										double minY_ptslists =	100000000;
										int i;
										
										ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
										
										i = 0;
										while( i < count_MPs && (fscanf(survey,"%lf %lf %lf %hhd\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z,&ptslists[i].flag)) != EOF )
										{
											if(maxX_ptslists < ptslists[i].m_X)
												maxX_ptslists = ptslists[i].m_X;
											if(maxY_ptslists < ptslists[i].m_Y)
												maxY_ptslists = ptslists[i].m_Y;
											if(minX_ptslists > ptslists[i].m_X)
												minX_ptslists = ptslists[i].m_X;
											if(minY_ptslists > ptslists[i].m_Y)
												minY_ptslists = ptslists[i].m_Y;
											
											i++;
										}
										fclose(survey);
										
										UI3DPOINT *trilists;
										
										double min_max[4] = {minX_ptslists, minY_ptslists, maxX_ptslists, maxY_ptslists};
										
										if(level >= TIN_split_level || count_MPs < 10000)
										{
											UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_MPs*4);
											
											sprintf(bufstr,"%s/txt/tri_ortho.txt",proinfo.save_filepath);
											TINCreate(ptslists,bufstr,count_MPs,t_trilists,min_max,&count_tri);
											
											trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
											i = 0;
											for(i=0;i<count_tri;i++)
											{
												trilists[i].m_X = t_trilists[i].m_X;
												trilists[i].m_Y = t_trilists[i].m_Y;
												trilists[i].m_Z = t_trilists[i].m_Z;
											}
											
											free(t_trilists);
										}
										else
										{
											int iter_row,iter_col;
											
											trilists = TINgeneration(true,proinfo.save_filepath, level, Size_Grid2D, Image_res[0], grid_resolution,
																	 min_max,
																	 subBoundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
										}
											
                                           
                                        if(pre_matched_pts == 0)
                                             matching_change_rate = 0;
                                         else
                                             matching_change_rate = fabs( (double)pre_matched_pts - (double)count_MPs ) /(double)pre_matched_pts;
										
										printf("matching change rate pre curr %f\t%d\t%d\n",matching_change_rate,count_MPs,pre_matched_pts);
										pre_matched_pts = count_results[0];
										
										if(iteration > 9 || level == 0 )
											matching_change_rate = 0.001;
										
                                        if(level <= 1)
                                        {
                                            if(iteration > 5)
                                                matching_change_rate = 0.001;
                                        }
									
										if(proinfo.IsRA)
											matching_change_rate = 0.001;
										
										if(proinfo.DEM_resolution >= 8)
										{
											if(level <= 2)
												matching_change_rate = 0.001;
											//if(level == 3 && iteration > 4)
											//	matching_change_rate = 0.001;
											matching_change_rate = 0.001;
										}
										
										if(proinfo.pre_DEMtif)
										{
											if(level >= 4)
												matching_change_rate = 0.001;
										}
										
										if(Th_roh >= Th_roh_min)
										{
											if(level == 0)
												Th_roh_update		= (double)(Th_roh - 0.50);
											else if(level == 1)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else if(level == 2)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else if(level == 3)
												Th_roh_update		= (double)(Th_roh - 0.10);
											else
											{
												if(proinfo.IsRA)
													Th_roh_update		= (double)(Th_roh - 0.10);
												else
													Th_roh_update		= (double)(Th_roh - 0.06);
											}
										}
										
										bool check_level_end = false;
										
										if(level != 0)
										{
											if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th)
											{
												check_level_end = true;
												
												if(dem_update_flag)
													update_flag		= true;
												else
													update_flag		= false;
											}
											else
												update_flag			= false;
										}
										else
											update_flag				= false;

                                        if(level == 0)
                                        {
                                            Pre_GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
                                            printf("update GridPT3\n");
                                        }
										else if(Th_roh_update < Th_roh_min && matching_change_rate < rate_th && level > 0)
										{
											Pre_GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
											printf("update GridPT3\n");
										}
										else
											GridPT3		= SetHeightRange(proinfo.pre_DEMtif,minmaxHeight,count_MPs, count_tri, GridPT3,update_flag,&minH_grid,&maxH_grid,blunder_param,ptslists,trilists,proinfo.IsRA,MPP,proinfo.save_filepath,row,col,check_level_end,proinfo.seedDEMsigma);
										
										free(trilists);
										free(ptslists);
                                    }
								}
								printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd SetHeightRange\n",row,col,level,iteration);
								
								fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\n",
										row,col,level,iteration,count_results[0],minH_grid,maxH_grid);
								printf("row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd iterpolation of Grids!! Mps = %d\tminH = %f\tmaxH = %f\n",
									   row,col,level,iteration,count_results[0],minH_grid,maxH_grid);

								if(proinfo.IsRA && level <= 3)
								{
									int RA_iter_counts = 0;
									RA_iter_counts = AdjustParam(level, count_MPs, filename_mps, Lstartpos, Rstartpos, LRPCs, RRPCs, t_Rimageparam, flag,
																 Template_size, SubImages_L, data_size_l[level], SubImages_R, data_size_r[level], SubOriImages_L, SubOriImages_R, param,
																 bin_angle, pyramid_step, Hemisphere, proinfo.save_filepath, proinfo.tmpdir);
									fprintf(fid,"RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Rimageparam[0],t_Rimageparam[1]);
									printf("RA iter = %d\tRA Line = %f\tSamp = %f\n",RA_iter_counts,t_Rimageparam[0],t_Rimageparam[1]);
									
									if (level <= 1)
									{
										//sprintf(save_file,"%s/txt/RAinfo.txt",proinfo.save_filepath,row,col);
										sprintf(save_file,"%s/txt/RAinfo.txt",proinfo.save_filepath);
										fid_RAinfo	= fopen(save_file,"w");
										fprintf(fid_RAinfo,"%f\t%f\n",t_Rimageparam[0],t_Rimageparam[1]);
										fclose(fid_RAinfo);
									}
								}
                                    
                                if(proinfo.IsRA)
                                {
                                    if(check_RA_divide)
                                    {
                                        if(level <= 3 && iteration > 2)
                                        {
                                            FILE *fid_pts		= fopen(filename_mps,"r");
                                            double t_X,t_Y,t_Z;
                                            bool t_flag;
                                            int* t_count = (int*)calloc(sizeof(int),total_tile);
                                            
                                            for(int k=0;k<count_results[0];k++)
                                            {
                                                fscanf(fid_pts,"%lf %lf %lf %hhd\n",&t_X,&t_Y,&t_Z,&t_flag);
                                                
                                                int t_col = floor((t_X - subBoundary[0])/(double)tilesize_RA);
                                                int t_row = floor((t_Y - subBoundary[1])/(double)tilesize_RA);
                                                
                                                t_count[t_col+division_X*t_row] ++;
                                            }
                                            fclose(fid_pts);
                                            
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
                                                
                                                //printf("k %d\tt_count %d\tsaved_count %d\tselected_X %d\tselected_Y %d\n",k,t_count[k],saved_count,selected_X,selected_Y);
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
                                    
                                //matching_rate = (double)count_results[0]/(double)total_matching_candidate_pts;
                                
                                //printf("total_matching_candidate_pts = %d\tMPs = %d\tmatching_rate = %f\n",total_matching_candidate_pts,count_results[0],matching_rate);
							}
							if (level == 0 && iteration == 3)
							{
								remove(filename_mps_pre);
								if(level >= ortho_level && proinfo.IsRA != 1)
								{
									remove(filename_mps_aft);
									remove(filename_mps_fin);
									remove(filename_mps_anchor);
								
								}
							}
							else 
							{
								remove(filename_mps);
								remove(filename_mps_pre);
								if(level >= ortho_level && proinfo.IsRA != 1)
								{
									remove(filename_mps_aft);
									remove(filename_mps_fin);
									remove(filename_mps_anchor);
								}
							}
							
							if(lower_level_match)
							{
								flag_start			= true;
								iteration++;
							}

							if(level == 0)
								Th_roh			= (double)(Th_roh - 0.50);
							else if(level == 1)
								Th_roh			= (double)(Th_roh - 0.10);
							else if(level == 2)
								Th_roh			= (double)(Th_roh - 0.10);
							else if(level == 3)
								Th_roh			= (double)(Th_roh - 0.10);
							else
							{
								if(proinfo.IsRA)
									Th_roh			= (double)(Th_roh - 0.10);
								else
									Th_roh			= (double)(Th_roh - 0.06);
							}

							if(lower_level_match)
							{
								if(Th_roh < Th_roh_min && matching_change_rate > rate_th)
								{
									if(level == 0)
										Th_roh			= (double)(Th_roh + 0.10);
									else if(level == 1)
										Th_roh			= (double)(Th_roh + 0.10);
									else if(level == 2)
										Th_roh			= (double)(Th_roh + 0.10);
									else if(level == 3)
										Th_roh			= (double)(Th_roh + 0.10);
									else
									{
										if(proinfo.IsRA)
											Th_roh			= (double)(Th_roh + 0.10);
										else
											Th_roh			= (double)(Th_roh + 0.06);
									}
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

						if(flag_start)
						{
                            double min_after, max_after;
                            min_after	= (double)(minH_mps - pow(2, level)*2*MPP);		max_after	= (double)(maxH_mps + pow(2, level)*2*MPP);
                            printf("minmax MP %f\t%f\n",min_after, max_after);
							if(level <= 2)
							{
                                if(minmaxHeight[0] < min_after)
                                    minmaxHeight[0]		= (double)(floor(min_after));
                                if(minmaxHeight[1] > max_after)
                                    minmaxHeight[1]		= (double)(ceil(max_after));
							}
                            else
                            {
                                if(min_after > minH_grid)
                                    min_after = minH_grid;
                                if(max_after < maxH_grid)
                                    max_after = maxH_grid;
                                
                                minmaxHeight[0]		= min_after;
                                minmaxHeight[1]		= max_after;
                            }
							printf("minmax %f\t%f\t\n", minmaxHeight[0],minmaxHeight[1]);
							fprintf(fid,"row = %d\tcol = %d\tlevel = %d\titeration = %d\tEnd of level processing!! minmaxHeight = [%f \t%f]\n",
									row,col,level,iteration,minmaxHeight[0],minmaxHeight[1]);
						}
						
						printf("\trow = %d/%d\tcol = %d/%d\tDEM generation(%%) = %4.2f%% !!\n",row,iter_row_end,col,t_col_end,(double)(pyramid_step+1 - level)/(double)(pyramid_step+1)*100);
						
                        if(level > 0)
                            level	= level - 1;

                        if(level == 0 && final_level_iteration == 4)
                            level = -1;
			
						if(proinfo.IsRA)
						{
							//if(!lower_level_match || level < DEM_level)
								free(GridPT);
							
							if(!lower_level_match)
							{
								lower_level_match	= true;
								flag_start			= false;
							}
						}
						else
						{
							if(!lower_level_match && level > 1)
							{
								lower_level_match	= true;
								flag_start			= false;
							}
							free(GridPT);
						}

                        printf("release Grid_wgs, nccresult\n");
						free(Grid_wgs);
						free(nccresult);
						
                        printf("release subImage L\n");
						free(SubImages_L);
						free(SubImages_BL);
						free(SubOriImages_L);
						
                        printf("release subImage R\n");
						free(SubImages_R);
						free(SubImages_BR);
						free(SubOriImages_R);

                        printf("release subimage Mag\n\n\n");
						free(SubMagImages_L);
						free(SubMagImages_R);
						free(SubMagImages_BL);
						free(SubMagImages_BR);
					}
                    printf("relese data size\n");
					free(data_size_l);
					free(data_size_r);
					free(GridPT3);
					
                    printf("release GridTP3\n");
					PreET = time(0);
					Pregab = difftime(PreET,PreST);
					printf("row = %d/%d \tcol = %d/%d\tDEM generation finish(time[m] = %5.2f)!!\n",row,iter_row_end,col,t_col_end,Pregab/60.0);
					
				}
			}
			fclose(fid);
			fclose(fid_header);

			RemoveFiles(proinfo.tmpdir,Lsubsetfilename,Rsubsetfilename,0,0);
			
			if(proinfo.IsRA)
			{
				if(t_Rimageparam[0] != 0 && t_Rimageparam[1] != 0)
				{
					RA_count++;
					Rimageparam[0] += t_Rimageparam[0];
					Rimageparam[1] += t_Rimageparam[1];
				}
			}
		}
	}
	free(iterations);
	if(proinfo.IsRA && RA_count > 0)
	{
		Rimageparam[0] /= RA_count;
		Rimageparam[1] /= RA_count;
	}
#ifdef BUILDMPI
	MPI_Bcast(Rimageparam, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	printf("Num of RAs = %d\tRA param = %f\t%f\n",RA_count,Rimageparam[0],Rimageparam[1]);
	return final_iteration;
}

bool OpenProject(char* _filename, ProInfo *info, ARGINFO args)
{
	bool bopened = false;
	FILE *pFile;
	char bufstr[500];

	info->IsRA			= false;
	info->IsRR			= false;
	info->IsSP			= false;
	info->IsSaveStep	= false;
	info->pre_DEMtif	= false;
	info->check_tiles_SR= false;
	info->check_tiles_ER= false;
	info->check_tiles_SC= false;
	info->check_tiles_EC= false;
	info->check_minH	= false;
	info->check_maxH	= false;
	info->check_checktiff	= false;
	info->tile_info[0]		= '\0';
	info->priori_DEM_tif[0]	= '\0';
	info->metafilename[0]	= '\0';

	pFile		= fopen(_filename,"r");
	if( pFile == NULL)
	{
		printf("'default.txt' file doesn't exist in SETSM execution folder, please check it!!\n");
        exit(1);
	}
	else
	{
		FILE *limage, *rimage;
		bopened = true;
		char* tmp_chr;
		
				
		while(!feof(pFile))
		{
			fgets(bufstr,500,pFile);
			if (strstr(bufstr,"Image_GSD")!=NULL)
				sscanf(bufstr,"Image_GSD %lf\n",&info->resolution);
			else if (!args.check_DEM_space && strstr(bufstr,"DEM_space")!=NULL)
				sscanf(bufstr,"DEM_space %lf\n",&info->DEM_resolution);
			else if (strstr(bufstr,"Image_RA")!=NULL)
			{
				if(!args.check_RA_line || !args.check_RA_sample)
				{
					sscanf(bufstr,"Image_RA %hhd %lf %lf\n",&info->IsRA,&info->RA_param[0],&info->RA_param[1]);
					if (info->IsRA == 1) {
						info->RA_param[0] = 0.0;
						info->RA_param[1] = 0.0;
					}
				}
			}
			else if (!args.check_Threads_num && strstr(bufstr,"Threads_num")!=NULL)
				sscanf(bufstr,"Threads_num %d\n",&info->threads_num);
			else if (args.check_arg == 0 && strstr(bufstr,"Image1")!=NULL)
			{
				sscanf(bufstr,"Image1 %s\n",info->LeftImagefilename);
                
                FILE *temp_ptif = fopen(info->LeftImagefilename,"r");
                if(temp_ptif)
                    printf("image1 load completed!\n");
                else
                {
                    printf("image1 load faied. Please check filename!!\n");
                    exit(1);
                }
                
				tmp_chr = remove_ext(info->LeftImagefilename);
				sprintf(info->LeftRPCfilename,"%s.xml",tmp_chr);
				
				FILE *temp_pFile;
				temp_pFile			 = fopen(info->LeftRPCfilename,"r");
				//printf("xml file %s\n",info->LeftRPCfilename);
				if(temp_pFile)
					printf("image1 xml load completed!\n");
				else
				{
					sprintf(info->LeftRPCfilename,"%s.XML",tmp_chr);
					//printf("xml file %s\n",info->LeftRPCfilename);
					temp_pFile			 = fopen(info->LeftRPCfilename,"r");
					if(temp_pFile)
						printf("image1 XML load completed!\n");
					else
					{
						printf("image1 xml/XML load failed!\n");
						exit(0);
					}
				}
			}
			else if (args.check_arg == 0 && strstr(bufstr,"Image2")!=NULL)
			{
				sscanf(bufstr,"Image2 %s\n",info->RightImagefilename);
                
                FILE *temp_ptif = fopen(info->RightImagefilename,"r");
                if(temp_ptif)
                    printf("image2 load completed!\n");
                else
                {
                    printf("image2 load faied. Please check filename!!\n");
                    exit(1);
                }
                
				tmp_chr = remove_ext(info->RightImagefilename);
				sprintf(info->RightRPCfilename,"%s.xml",tmp_chr);
				
				FILE *temp_pFile;
				temp_pFile			 = fopen(info->RightRPCfilename,"r");
				//printf("xml file %s\n",info->RightRPCfilename);
				if(temp_pFile)
					printf("image2 xml load completed!\n");
				else
				{
					sprintf(info->RightRPCfilename,"%s.XML",tmp_chr);
					//printf("xml file %s\n",info->RightRPCfilename);
					temp_pFile			 = fopen(info->RightRPCfilename,"r");
					if(temp_pFile)
						printf("image2 XML load completed!\n");
					else
					{
						printf("image2 xml/XML load failed!\n");
						exit(0);
					}
				}
			}
			else if (args.check_arg == 0 && strstr(bufstr,"outputpath")!=NULL)
			{
				sscanf(bufstr,"outputpath %s\n",info->save_filepath);
				
				
				char *Outputpath_name  = SetOutpathName(info->save_filepath);
				sprintf(info->Outputpath_name,"%s",Outputpath_name);
				printf("after pathname %s\n",info->Outputpath_name);
				
				sprintf(info->tmpdir, "%s/tmp", info->save_filepath);
			}
			else if (args.check_seeddem == 0 && strstr(bufstr,"seeddempath")!=NULL)
			{
				sscanf(bufstr,"seeddempath %hhd ",&info->pre_DEMtif);
				if(info->pre_DEMtif)
				{
					sscanf(bufstr,"seeddempath %hhd %s %lf\n",&info->pre_DEMtif,info->priori_DEM_tif,&info->seedDEMsigma);
				}
				printf("seeddempath %d %s %lf\n",info->pre_DEMtif,info->priori_DEM_tif,info->seedDEMsigma);
			}
			else if (args.check_tiles_SR == 0 && strstr(bufstr,"tile_SR")!=NULL)
			{
				sscanf(bufstr,"tile_SR %d\n",&info->start_row);
				info->check_tiles_SR = true;
				printf("%d\n",info->start_row);
			}
			else if (args.check_tiles_ER == 0 && strstr(bufstr,"tile_ER")!=NULL)
			{
				sscanf(bufstr,"tile_ER %d\n",&info->end_row);
				info->check_tiles_ER = true;
				printf("%d\n",info->end_row);
			}
			else if (args.check_tiles_SC == 0 && strstr(bufstr,"tile_SC")!=NULL)
			{
				sscanf(bufstr,"tile_SC %d\n",&info->start_col);
				info->check_tiles_SC = true;
				printf("%d\n",info->start_col);
			}
			else if (args.check_tiles_EC == 0 && strstr(bufstr,"tile_EC")!=NULL)
			{
				sscanf(bufstr,"tile_EC %d\n",&info->end_col);
				info->check_tiles_EC = true;
				printf("%d\n",info->end_col);
			}
			else if (args.check_minH == 0 && strstr(bufstr,"minH")!=NULL)
			{
				sscanf(bufstr,"minH %lf\n",&info->minHeight);
				info->check_minH = true;
				printf("%f\n",info->minHeight);
			}
			else if (args.check_maxH == 0 && strstr(bufstr,"maxH")!=NULL)
			{
				sscanf(bufstr,"maxH %lf\n",&info->maxHeight);
				info->check_maxH = true;
				printf("%f\n",info->maxHeight);
			}
		}
				
		if (args.check_arg != 0) {
			sprintf(info->tmpdir,"%s/tmp",args.Outputpath);
			sprintf(info->LeftImagefilename,"%s", args.Image1);
			sprintf(info->RightImagefilename,"%s",args.Image2);
			
            FILE *temp_ptif = fopen(info->LeftImagefilename,"r");
            if(temp_ptif)
                printf("image1 load completed!\n");
            else
            {
                printf("image1 load faied. Please check filename!!\n");
                exit(1);
            }
            
			tmp_chr = remove_ext(args.Image1);
			sprintf(info->LeftRPCfilename,"%s.xml",tmp_chr);
			
			FILE *temp_pFile;
			temp_pFile			 = fopen(info->LeftRPCfilename,"r");
			//printf("xml file %s\n",info->LeftRPCfilename);
			if(temp_pFile)
				printf("image1 xml load completed!\n");
			else
			{
				sprintf(info->LeftRPCfilename,"%s.XML",tmp_chr);
				//printf("xml file %s\n",info->LeftRPCfilename);
				temp_pFile			 = fopen(info->LeftRPCfilename,"r");
				if(temp_pFile)
					printf("image1 XML load completed!\n");
				else
				{
					printf("image1 xml/XML load failed!\n");
					exit(0);
				}
			}
			
            temp_ptif = fopen(info->RightImagefilename,"r");
            if(temp_ptif)
                printf("image2 load completed!\n");
            else
            {
                printf("image2 load faied. Please check filename!!\n");
                exit(1);
            }
            
			tmp_chr = remove_ext(args.Image2);
			sprintf(info->RightRPCfilename,"%s.xml",tmp_chr);
			
			temp_pFile			 = fopen(info->RightRPCfilename,"r");
			//printf("xml file %s\n",info->RightRPCfilename);
			if(temp_pFile)
				printf("image2 xml load completed!\n");
			else
			{
				sprintf(info->RightRPCfilename,"%s.XML",tmp_chr);
				//printf("xml file %s\n",info->RightRPCfilename);
				temp_pFile			 = fopen(info->RightRPCfilename,"r");
				if(temp_pFile)
					printf("image2 XML load completed!\n");
				else
				{
					printf("image2 xml/XML load failed!\n");
					exit(0);
				}
			}
			
			sprintf(info->save_filepath,"%s",args.Outputpath);
			sprintf(info->Outputpath_name, "%s", args.Outputpath_name);
		}
		
		if (args.check_gridonly) {
			info->check_gridonly = true;
			sprintf(info->save_filepath,"%s",args.Outputpath);
		}
		
		if (args.check_checktiff)
		{
			info->check_checktiff = true;
		}
		
		if (args.check_DEM_space) {
			info->DEM_resolution = args.DEM_space;
		}
		if (args.check_Threads_num) {
			info->threads_num = args.Threads_num;
		}
		if (args.check_seeddem) {
			sprintf(info->priori_DEM_tif,"%s",args.seedDEMfilename);
			sprintf(info->metafilename,"%s",args.metafilename);
			info->pre_DEMtif = args.check_seeddem;
			info->seedDEMsigma = args.seedDEMsigma;
		}
		
		printf("left image = %s\n",info->LeftImagefilename);
		printf("right image = %s\n",info->RightImagefilename);
		printf("left rpc = %s\n",info->LeftRPCfilename);
		printf("right rpc = %s\n",info->RightRPCfilename);
		printf("save = %s\n",info->save_filepath);
		printf("save_name = %s\n", info->Outputpath_name);
		printf("DEM space = %f\n",info->DEM_resolution);
		printf("Threads = %d\n",info->threads_num);
		printf("seeddem flag = %d\n",info->pre_DEMtif);
		if(info->pre_DEMtif)
		{
			printf("seeddem = %s seedDEMsigma = %f\n",info->priori_DEM_tif, info->seedDEMsigma);
			printf("meta file = %s\n",info->metafilename);
			
		}
		limage = fopen(info->LeftImagefilename,"r");
		rimage = fopen(info->RightImagefilename,"r");
		if(!limage || !rimage)
		{
			bopened		= false;
			printf("Check input image filename!!\n");
		}
		if(limage)
			fclose(limage);
		if(rimage)
			fclose(rimage);
	}

	
	return bopened;
}


int Maketmpfolders(ProInfo info)
{
	char temp_filepath[500];

	int check_folder = 1;
	
	if(!info.check_checktiff)
	{
        int status;
		status = mkdir(info.save_filepath,0777);
        if (opendir(info.save_filepath) == NULL)
        {
            if (status == -1)
            {
                printf("Outpath of '%s' cannot make, please check outpath!!\n",info.save_filepath);
                exit(1);
            }
        }
        sprintf(temp_filepath,"%s/txt",info.save_filepath);
		mkdir(temp_filepath,0777);
		sprintf(temp_filepath,"%s/tif",info.save_filepath);
		mkdir(temp_filepath,0777);
		sprintf(temp_filepath,"%s/tmp",info.save_filepath);
		mkdir(temp_filepath,0777);
	}
	return check_folder;
}

bool SetupParam(ProInfo info,uint8 *NumOfIAparam, uint8 *pre_DEM_level, uint8 *DEM_level,  bool *pre_DEMtif, bool *check_tile_array)
{
	bool cal_check		= true;
	
	
	*NumOfIAparam= 2;
	
	printf("seed flag %d\n",info.pre_DEMtif);
	
	if (info.pre_DEMtif) 
	{
		FILE *pFile_DEM;
		printf("dem filename %s\n",info.priori_DEM_tif);
		pFile_DEM	= fopen(info.priori_DEM_tif,"r");
		if(pFile_DEM)
			*pre_DEMtif	= true;
		else
			*pre_DEMtif	= false;
		fclose(pFile_DEM);
	}
	
	printf("seed flag %d\n",*pre_DEMtif);

	*pre_DEM_level = 0;
	*DEM_level	= 0;

	
	*check_tile_array	= false;
	*check_tile_array	= false;

	return cal_check;
}


void SetTransParam(double minLat, double minLon, bool *Hemisphere, TransParam *param)
{
	if(minLat > 0)
	{
		*Hemisphere	= true;

		param->a = (double)(6378137.0);
		param->e = (double)(0.08181919);
		param->phi_c = (double)(70.0*DegToRad);
		param->lambda_0 = (double)(-45.0*DegToRad);
		param->bHemisphere = *Hemisphere;
	}
	else
	{
		*Hemisphere	= false;

		param->a = (double)(6378137.0);
		param->e = (double)(0.08181919);
		param->phi_c = (double)(-71.0*DegToRad);
		param->lambda_0 = 0;
		param->bHemisphere = *Hemisphere;
	}
	
	if(param->phi_c < 0)
	{
		param->pm		= -1;
		param->phi_c	= -param->phi_c;
		param->lambda_0 = -param->lambda_0;
	}
	else
		param->pm		= 1;
	
	param->t_c	= tan(PI/4.0 - param->phi_c/2.0)/(pow((1.0-param->e*sin(param->phi_c)) / (1.0+param->e*sin(param->phi_c)) ,(param->e/2.0)));
	param->m_c	= cos(param->phi_c)/sqrt(1.0-pow(param->e,2)*pow(sin(param->phi_c),2));
	
	//UTM param
	param->sa	= 6378137.000000;
	param->sb	= 6356752.314245;
	param->e2	= ( sqrt( ( param->sa*param->sa ) - ( param->sb*param->sb ) ) ) / param->sb;
	param->e2cuadrada = param->e2*param->e2;
	param->c	= (param->sa*param->sa)/param->sb;
	
	
	printf("%f %f %f %f %f\n",param->sa,param->sb,param->e2,param->e2cuadrada,param->c);
	double Lat = minLat;
	char direction[10];
	
	if(param->projection == 3)
	{
		if(Lat > -60 && Lat < 60)
			param->projection = 2;
		else
			param->projection = 1;
	}
		
	if(minLat < -72)
		sprintf(direction,"%s","C");
	else if(Lat<-64 && Lat>=-72)
		sprintf(direction,"%s","D");
	else if(Lat<-56 && Lat>=-64)
		sprintf(direction,"%s","E");
	else if(Lat<-48 && Lat>=-56)
		sprintf(direction,"%s","F");
	else if(Lat<-40 && Lat>=-48)
		sprintf(direction,"%s","G");
	else if(Lat<-32 && Lat>=-40)
		sprintf(direction,"%s","H");
	else if(Lat<-24 && Lat>=-32)
		sprintf(direction,"%s","J");
	else if(Lat<-16 && Lat>=-24)
		sprintf(direction,"%s","K");
	else if(Lat<-8 && Lat>=-16)
		sprintf(direction,"%s","L");
	else if(Lat<0 && Lat>=-8)
		sprintf(direction,"%s","M");
	else if(Lat<8 && Lat>=0)
		sprintf(direction,"%s","N");
	else if(Lat<16 && Lat>=8)
		sprintf(direction,"%s","P");
	else if(Lat<24 && Lat>=16)
		sprintf(direction,"%s","Q");
	else if(Lat<32 && Lat>=24)
		sprintf(direction,"%s","R");
	else if(Lat<40 && Lat>=32)
		sprintf(direction,"%s","S");
	else if(Lat<48 && Lat>=40)
		sprintf(direction,"%s","T");
	else if(Lat<56 && Lat>=48)
		sprintf(direction,"%s","U");
	else if(Lat<64 && Lat>=56)
		sprintf(direction,"%s","V");
	else if(Lat<72 && Lat>=64)
		sprintf(direction,"%s","W");
	else
		sprintf(direction,"%s","X");
	
	
	sprintf(param->direction,"%s",direction);
	
    if(param->utm_zone < -1)
		param->zone = (int)( ( minLon / 6 ) + 31);
    else
        param->zone = param->utm_zone;
    
}
void SetTiles(ProInfo info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
			  uint8 *iter_row_start, uint8 *iter_row_end, uint8 *t_col_start, uint8 *t_col_end, double *subX, double *subY)
{
	int lengthOfX = Boundary[2] - Boundary[0];
	int lengthOfY = Boundary[3] - Boundary[1];
	int division_X = (int) (ceil(lengthOfX / (double)(tile_size)));
	int division_Y = (int) (ceil(lengthOfY / (double)(tile_size)));
	*subX = floor((ceil(ceil(lengthOfX / division_X) / info.DEM_resolution) * info.DEM_resolution) / 2) * 2;
	*subY = floor((ceil(ceil(lengthOfY / division_Y) / info.DEM_resolution) * info.DEM_resolution) / 2) * 2;
	
	*pyramid_step	= 4;

	if(info.DEM_resolution	>= 10)
		*buffer_area	= (uint16)(*buffer_area * 1.5);

	if(!IsSP)
	{
		*iter_row_start	= 1;
		*iter_row_end	= division_Y+1;
	}
	else
	{
		int iter_row;
		iter_row		= (int)(ceil(division_Y/info.SPnumber[1]));
		*iter_row_start	= iter_row*(info.SPnumber[0] - 1) + 1;
		*iter_row_end	= iter_row*info.SPnumber[0];

		if(info.SPnumber[0] == info.SPnumber[1])
			*iter_row_end= division_Y+1;
	}
	
	*t_col_start			= 1;
	*t_col_end			= division_X+1;
	if(IsRR)
		*iter_row_start	= info.start_row;
	
	if(*iter_row_start == *iter_row_end)
		*iter_row_end += 1;
	if(*t_col_start == *t_col_end)
		*t_col_end += 1;
}


void SetTiles_RA(ProInfo info, bool IsSP, bool IsRR, double *Boundary, double *Res, int tile_size, bool pre_DEMtif, uint8 *pyramid_step, uint16 *buffer_area,
				 uint8 *RA_row_start, uint8 *RA_row_end, uint8 * RA_row_iter, uint8 *t_col_start, uint8 *t_col_end, uint8 *RA_col_iter, double *subX, double *subY)
{
	int lengthOfX = Boundary[2] - Boundary[0];
	int lengthOfY = Boundary[3] - Boundary[1];
	int division_X = (int) (ceil(lengthOfX / (double)(tile_size)));
	int division_Y = (int) (ceil(lengthOfY / (double)(tile_size)));
	*subX = floor((ceil(ceil(lengthOfX / division_X) / info.DEM_resolution) * info.DEM_resolution) / 2) * 2;
	*subY = floor((ceil(ceil(lengthOfY / division_Y) / info.DEM_resolution) * info.DEM_resolution) / 2) * 2;
	
	printf("%d %d tile size %d\n", lengthOfX, lengthOfY, tile_size);
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


void SetPySizes(CSize *data_size_l, CSize *data_size_r, CSize Lsubsetsize, CSize Rsubsetsize, int level)
{
	int i;

	data_size_l[0].height		= Lsubsetsize.height;
	data_size_r[0].height		= Rsubsetsize.height;
	data_size_l[0].width		= Lsubsetsize.width;
	data_size_r[0].width		= Rsubsetsize.width;
	for(i=0;i<level;i++)
	{
		data_size_l[i+1].width	= data_size_l[i].width/2;
		data_size_l[i+1].height = data_size_l[i].height/2;
		data_size_r[i+1].width	= data_size_r[i].width/2;
		data_size_r[i+1].height = data_size_r[i].height/2;
	}
}

void SetThs(int level, int final_level_iteration, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, int pre_DEMtif, int IsRA, double seedDEMsigma, double f_demsize)
{
	if(f_demsize >= 8)
	{
		if(level >= 4)
		{
			*Th_roh			 = (double)(0.80);
			*Th_roh_min		 = (double)(0.39);
			
			if(IsRA)
			{
				*Th_roh			 = (double)(0.80);
				*Th_roh_min		 = (double)(0.49);
			}
		}
		else if(level >= 3)
		{
			*Th_roh			 = (double)(0.70);
			*Th_roh_min		 = (double)(0.39);
		}
		else if(level == 2)
		{
			*Th_roh			 = (double)(0.60);
			*Th_roh_min		 = (double)(0.49);
		}
		else if(level == 1)
		{
			*Th_roh			 = (double)(0.50);
			*Th_roh_min		 = (double)(0.39);
		}
		else if(level == 0)
		{
			*Th_roh			 = (double)(0.40 - 0.10*(final_level_iteration-1));
			*Th_roh_min		 = (double)(0.19);
		}

		
		if(level == 5)
			*Th_roh_next		= (double)(0.80);
		else if(level == 4)
			*Th_roh_next		= (double)(0.70);
		else if(level == 3)
			*Th_roh_next		= (double)(0.60);
		else if(level == 2)
			*Th_roh_next		= (double)(0.50);
		else if(level == 1)
			*Th_roh_next		= (double)(0.40);
		else
			*Th_roh_next		= (double)(0.20);

		
		
		*Th_roh_start		= (double)(*Th_roh);
		
		
		if(pre_DEMtif && !IsRA)
		{
			if(level >= 4)
			{
				*Th_roh			 = (double)(0.80);
				*Th_roh_min		 = (double)(0.39);
			}
			else if(level >= 3)
			{
				*Th_roh			 = (double)(0.70);
				*Th_roh_min		 = (double)(0.39);
			}
			else if(level == 2)
			{
				*Th_roh			 = (double)(0.60);
				*Th_roh_min		 = (double)(0.49);
			}
			else if(level == 1)
			{
				*Th_roh			 = (double)(0.50);
				*Th_roh_min		 = (double)(0.39);
			}
			else if(level == 0)
			{
				*Th_roh			 = (double)(0.40 - 0.10*(final_level_iteration-1));
				*Th_roh_min		 = (double)(0.19);
			}
			
			
			if(level == 5)
				*Th_roh_next		= (double)(0.80);
			else if(level == 4)
				*Th_roh_next		= (double)(0.70);
			else if(level == 3)
				*Th_roh_next		= (double)(0.60);
			else if(level == 2)
				*Th_roh_next		= (double)(0.50);
			else if(level == 1)
				*Th_roh_next		= (double)(0.40);
			else
				*Th_roh_next		= (double)(0.20);
			
			
			if(seedDEMsigma <= 20)
			{
				if(level == 2)
				{
					*Th_roh			 = (double)(0.50);
					*Th_roh_min		 = (double)(0.29);
				}
				else if(level == 1)
				{
					*Th_roh			 = (double)(0.40);
					*Th_roh_min		 = (double)(0.19);
				}
				else if(level == 0)
				{
					*Th_roh			 = (double)(0.40 - 0.10*(final_level_iteration-1));
					*Th_roh_min		 = (double)(0.19);
				}
				
				
				
				if(level == 2)
					*Th_roh_next		= (double)(0.40);
				else if(level == 1)
					*Th_roh_next		= (double)(0.40);
				else
					*Th_roh_next		= (double)(0.20);
			}
			

			*Th_roh_start		= (double)(*Th_roh);
			
			
		}
	}
	else
	{
		if(level >= 4)
		{
			*Th_roh			 = (double)(0.80);
			*Th_roh_min		 = (double)(0.39);
			
			if(IsRA)
			{
				*Th_roh			 = (double)(0.80);
				*Th_roh_min		 = (double)(0.49);
			}
		}
		else if(level >= 3)
		{
			*Th_roh			 = (double)(0.70);
			*Th_roh_min		 = (double)(0.39);
		}
		else if(level == 2)
		{
			*Th_roh			 = (double)(0.50);
			*Th_roh_min		 = (double)(0.19);
		}
		else if(level == 1)
		{
			*Th_roh			 = (double)(0.40);
			*Th_roh_min		 = (double)(0.19);
		}
		else if(level == 0)
		{
			*Th_roh			 = (double)(0.40 - 0.10*(final_level_iteration-1));
			*Th_roh_min		 = (double)(0.19);
		}
		
		
		if(level == 5)
			*Th_roh_next		= (double)(0.80);
		else if(level == 4)
			*Th_roh_next		= (double)(0.70);
		else if(level == 3)
			*Th_roh_next		= (double)(0.50);
		else if(level == 2)
			*Th_roh_next		= (double)(0.40);
		else if(level == 1)
			*Th_roh_next		= (double)(0.40);
		else
			*Th_roh_next		= (double)(0.20);
		
		
		
		*Th_roh_start		= (double)(*Th_roh);
		
		
		if(pre_DEMtif && !IsRA)
		{
			if(level >= 4)
			{
				*Th_roh			 = (double)(0.80);
				*Th_roh_min		 = (double)(0.39);
			}
			else if(level >= 3)
			{
				*Th_roh			 = (double)(0.70);
				*Th_roh_min		 = (double)(0.39);
			}
			else if(level == 2)
			{
				*Th_roh			 = (double)(0.50);
				*Th_roh_min		 = (double)(0.19);
			}
			else if(level == 1)
			{
				*Th_roh			 = (double)(0.40);
				*Th_roh_min		 = (double)(0.19);
			}
			else if(level == 0)
			{
				*Th_roh			 = (double)(0.40 - 0.10*(final_level_iteration-1));
				*Th_roh_min		 = (double)(0.19);
			}
			
			
			if(level == 5)
				*Th_roh_next		= (double)(0.80);
			else if(level == 4)
				*Th_roh_next		= (double)(0.70);
			else if(level == 3)
				*Th_roh_next		= (double)(0.50);
			else if(level == 2)
				*Th_roh_next		= (double)(0.40);
			else if(level == 1)
				*Th_roh_next		= (double)(0.40);
			else
				*Th_roh_next		= (double)(0.20);
			
			*Th_roh_start		= (double)(*Th_roh);
		}
	}
}

void SetThs_ratio(int level, double *Th_roh, double *Th_roh_min, double *Th_roh_next, double *Th_roh_start, int pre_DEMtif, int IsRA, double f_demsize)
{
	if(pre_DEMtif && !IsRA)
	{
		if(level >= 4)
		{
			*Th_roh			 = (double)(0.85);
			*Th_roh_min		 = (double)(0.59);
		}
		else if(level >= 3)
		{
			*Th_roh			 = (double)(0.80);
			*Th_roh_min		 = (double)(0.49);
		}
		else if(level == 2)
		{
			*Th_roh			 = (double)(0.70);
			*Th_roh_min		 = (double)(0.49);
		}
		else if(level == 1)
		{
			*Th_roh			 = (double)(0.60);
			*Th_roh_min		 = (double)(0.399);
		}
		else if(level == 0)
		{
			*Th_roh			 = (double)(0.30);
			*Th_roh_min		 = (double)(0.099);
		}
		
		
		if(level == 5)
			*Th_roh_next		= (double)(0.85);
		else if(level == 4)
			*Th_roh_next		= (double)(0.80);
		else if(level == 3)
			*Th_roh_next		= (double)(0.70);
		else if(level == 2)
			*Th_roh_next		= (double)(0.60);
		else if(level == 1)
			*Th_roh_next		= (double)(0.30);
		else
			*Th_roh_next		= (double)(0.20);
		
		*Th_roh_start		= (double)(*Th_roh);
	}
}

D2DPOINT *SetGrids(bool *dem_update_flag, bool flag_start, int level, int final_level_iteration, double resolution, CSize *Size_Grid2D, bool pre_DEMtif, char *priori_DEM_tif, double DEM_resolution, double *minmaxHeight,
				   double *py_resolution, double *grid_resolution, double *subBoundary)
{
    *py_resolution	 = (double)(resolution*pow(2,level));
    *grid_resolution = *py_resolution;
    
    double R_resolution = 4*(*py_resolution);
    
    printf("pre resolution %f\t level %d\t final_level_iteration %d\n",*py_resolution,level,final_level_iteration);
    
    
    
	if(*dem_update_flag)
	{
        if(*py_resolution > 32)//low-res original imagery
        {
            *py_resolution = (int)(*py_resolution/4.0);
            *grid_resolution = *py_resolution;
        }
        else
        {
            if(resolution >= 0.4)
            {
                /*if(level >= 3)
                {
                    
                }
                else*/ if(level > 0)
                {
                    /*if(DEM_resolution >= 2)
                    {
                        if(R_resolution > DEM_resolution)
                        {
                            if(R_resolution > 8)
                                *py_resolution = 8;
                            else
                                *py_resolution	 = R_resolution;
                        }
                        else
                            *py_resolution = DEM_resolution;
                    }
                    else
                    {*/
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
                                    *py_resolution	 = (*py_resolution)*3;
                            }
                        }
                        else
                            *py_resolution = DEM_resolution;
                    //}
                }
                else if(level == 0)
                {
                    if(final_level_iteration == 1)
                    {
                        if(*py_resolution < 2)
                        {
                            if((*py_resolution)*4 > DEM_resolution)
                                *py_resolution	 = (*py_resolution)*4;
                            else
                                *py_resolution	 = DEM_resolution;
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
                                *py_resolution	 = (*py_resolution)*2;
                            else
                                *py_resolution	 = DEM_resolution;
                        }
                        else
                        {
                            *py_resolution = DEM_resolution;
                        }
                    }
                    else
                        *py_resolution	 = DEM_resolution;
                }
            }
            else
            {
                /*if(level > 0)
                {
                    if(*py_resolution < DEM_resolution)
                        *py_resolution = DEM_resolution;
                }
                else*/ if(level == 0)
                {
                    if(final_level_iteration == 1)
                    {
                        if((*py_resolution)*4 > DEM_resolution)
                            *py_resolution	 = (*py_resolution)*4;
                        else
                            *py_resolution	 = DEM_resolution;
                    }
                    else if(final_level_iteration == 2)
                    {
                        if((*py_resolution)*2 > DEM_resolution)
                            *py_resolution	 = (*py_resolution)*2;
                        else
                            *py_resolution	 = DEM_resolution;
                    }
                    else
                        *py_resolution	 = DEM_resolution;
                }
            }
        }
        
        *grid_resolution = *py_resolution;
        
        
        /*if(resolution < 0.5)
        {
            if(*py_resolution < DEM_resolution)
            {
                *py_resolution	 = DEM_resolution;
                *grid_resolution = DEM_resolution;
            }
        }
        else
        {
            if(DEM_resolution >= 8)
            {
                if(*py_resolution < DEM_resolution)
                {
                    *py_resolution	 = DEM_resolution;
                    *grid_resolution = DEM_resolution;
                }
                
                if(level == 0)
                {
                    *py_resolution	 = DEM_resolution;
                    *grid_resolution = DEM_resolution;
                }
            }
            else if(DEM_resolution >= 4 && DEM_resolution < 8)
            {
                if(level == 0)
                {
                    *py_resolution	 = DEM_resolution;
                    *grid_resolution = DEM_resolution;
                }
                else if(level <= 3)
                {
                    if(*py_resolution <= 4)
                    {
                        *py_resolution	 = DEM_resolution;
                        *grid_resolution = DEM_resolution;
                    }
                }
                else
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
            }
            else if(DEM_resolution >= 2 && DEM_resolution < 4)
            {
                if(level == 1)
                {
                    if(*py_resolution < 4)
                    {
                        *py_resolution	 = 4;
                        *grid_resolution = 4;
                    }
                    
                }
                else if(level >= 2)
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
                else if(level == 0)
                {
                    *py_resolution	 = DEM_resolution;
                    *grid_resolution = DEM_resolution;
                }
                else
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
            }
            else if(DEM_resolution < 2)
            {
                
                if(level == 1)
                {
                    if(*py_resolution < 4)
                    {
                        *py_resolution	 = 4;
                        *grid_resolution = 4;
                    }
                    
                }
                else if(level == 2)
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
                else if(level >= 3)
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
                else if(level == 0)
                {
                    if(final_level_iteration == 1)
                    {
                        if((*py_resolution)*4 > DEM_resolution)
                        {
                            *py_resolution	 = (*py_resolution)*4;
                            *grid_resolution = *py_resolution;
                        }
                        else
                        {
                            *py_resolution	 = DEM_resolution;
                            *grid_resolution = DEM_resolution;
                        }
                    }
                    else if(final_level_iteration == 2)
                    {
                        if((*py_resolution)*2 > DEM_resolution)
                        {
                            *py_resolution	 = (*py_resolution)*2;
                            *grid_resolution = *py_resolution;
                        }
                        else
                        {
                            *py_resolution	 = DEM_resolution;
                            *grid_resolution = DEM_resolution;
                        }
                    }
                    else
                    {
                        *py_resolution	 = DEM_resolution;
                        *grid_resolution = DEM_resolution;
                    }
                }
                else
                {
                    if(*py_resolution < 8)
                    {
                        *py_resolution	 = 8;
                        *grid_resolution = 8;
                    }
                }
                
            }
         
            if(pre_DEMtif && level == 2)
            {
                if(DEM_resolution == 2)
                {
                    if(level == 1)
                    {
                        if(*py_resolution < 4)
                        {
                            *py_resolution	 = 4;
                            *grid_resolution = 4;
                        }
                        
                    }
                    else if(level == 0)
                    {
                        *py_resolution	 = DEM_resolution;
                        *grid_resolution = DEM_resolution;
                }
                    else
                    {
                        if(*py_resolution < 8)
                        {
                            *py_resolution	 = 8;
                            *grid_resolution = 8;
                        }
                    }
                }
                else if(DEM_resolution < 2)
                {
                    if(level == 1)
                    {
                        if(*py_resolution < 2)
                        {
                            *py_resolution	 = 2;
                            *grid_resolution = 2;
                        }
                        
                    }
                    else if(level == 0)
                    {
                        *py_resolution	 = DEM_resolution;
                        *grid_resolution = DEM_resolution;
                    }
                    else
                    {
                        if(*py_resolution < 4)
                        {
                            *py_resolution	 = 4;
                            *grid_resolution = 4;
                        }
                    }
                }
                
            }
         
        }
        */
		Size_Grid2D->width	= (int)(ceil((double)(subBoundary[2] - subBoundary[0])/(*grid_resolution) ));
		Size_Grid2D->height	= (int)(ceil((double)(subBoundary[3] - subBoundary[1])/(*grid_resolution) ));
	}
 
	printf("DEM resolution %f\tresolution %f\t size %d\t%d\n",DEM_resolution,*py_resolution,Size_Grid2D->width,Size_Grid2D->height);
	D2DPOINT *GridPT = SetDEMGrid(subBoundary, *grid_resolution, *grid_resolution, Size_Grid2D);

	return GridPT;
}

UGRID *SetGrid3PT(TransParam param, bool dem_update_flag, bool flag_start, CSize Size_Grid2D, double Th_roh, int level, double *minmaxHeight,double *subBoundary,double py_resolution,
				  char* priori_DEM_tif,bool pre_DEMtif, double seedDEMsigma, int IsRA,char* metafilename)
{
	UGRID *GridPT3;
	int total_grid_counts, i;

	if(!flag_start)
	{
		
		total_grid_counts		= Size_Grid2D.height*Size_Grid2D.width;
		GridPT3					= (UGRID*)calloc(sizeof(UGRID),total_grid_counts);

#pragma omp parallel for shared(total_grid_counts,GridPT3,minmaxHeight,Th_roh) private(i)
		for(i=0;i<total_grid_counts;i++)
		{
			GridPT3[i].Matched_flag		= 0;
			GridPT3[i].roh				= Th_roh;
			GridPT3[i].anchor_flag		= 0;
			GridPT3[i].Matched_height	= -1000.0;
			GridPT3[i].ortho_ncc		= 0;
			GridPT3[i].angle			= 0;
//			GridPT3[i].false_h_count	= 0;

			GridPT3[i].minHeight		= (double)(minmaxHeight[0] - 0.5);
			GridPT3[i].maxHeight		= (double)(minmaxHeight[1] + 0.5);
			GridPT3[i].Height			= -1000.0;
		}
		if(pre_DEMtif)
		{
			printf("seedem load\n");
			SetHeightWithSeedDEM(param, GridPT3,subBoundary,Size_Grid2D,py_resolution, priori_DEM_tif,minmaxHeight,seedDEMsigma,IsRA,metafilename);
		}
	}

	return GridPT3;
}

int	 CalTotalIteration(uint8 DEM_level,int level)
{
	int sub_total_count;
	//compute a total iterations. (need to modify)
	if(DEM_level == 0)
	{
		if(level >= 3)
			sub_total_count = 4*(level - 2) + 2*3;
		else
			sub_total_count = 2*3;
	}
	else if(DEM_level == 1)
	{
		if(level >= 3)
			sub_total_count = 4*(level - 2) + 2*2;
		else
			sub_total_count = 2*2;
	}
	else if(DEM_level == 2)
	{
		if(level >= 3)
			sub_total_count = 4*(level - 2) + 2*1;
		else
			sub_total_count = 2*1;
	}
	else if(DEM_level == 3)
		sub_total_count = 4*(level - 2);
	else if(DEM_level == 4)
		sub_total_count = 4*(level - 3);
	else
		sub_total_count = 1;

	return sub_total_count;
}

char* remove_ext(char* mystr)
{
	char *retstr;
	char *lastdot;
	if (mystr == NULL)
		return NULL;
	if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
		return NULL;
	strcpy (retstr, mystr);
	lastdot = strrchr (retstr, '.');
	if (lastdot != NULL)
		*lastdot = '\0';
	return retstr;
}

char* GetFileName(char file_path[])
{
 
	char *file_name = NULL;
 
	while(*file_path)
	{
		if(*file_path == '/' && (file_path +1) != NULL )
		{
			file_name = file_path+1;
		}
		else
		{ }
			
		file_path++;
			   
	} 
	return file_name;
}

char* GetFileDir(char file_path[],int *size)
{
	
	char *file_name = NULL;
	*size = 0;
	while(*file_path)
	{
		if(*file_path == '/' && (file_path +1) != NULL )
		{
			
		}
		else
		{
			file_name = file_path+1;
			printf("str %s\n",file_name);
		}
		
		file_path++;
		*size = *size + 1;
	} 
	return file_name;
}


bool GetImageSize(char *filename, CSize *Imagesize)
{
	bool ret = false;
	
	char *ext;
	ext = strrchr(filename,'.');
	
	if(!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		TIFF *tif  = TIFFOpen(filename,"r");

		if(tif)
		{
			TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&Imagesize->width);
			TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&Imagesize->height);
			TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&Imagesize->width);
			TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&Imagesize->height);
			ret = true;
		}

		TIFFClose(tif);
	}
	else if(!strcmp("bin",ext+1))
	{
		char *tmp;
		tmp = remove_ext(filename);
		sprintf(tmp,"%s.hdr",tmp);
		*Imagesize = Envihdr_reader(tmp);
		
		ret = true;
		
	}
	return ret;
}


bool GetsubareaImage(TransParam transparam, uint8 NumofIAparam, double **RPCs, double *ImageParam, char *ImageFilename, CSize *Imagesize,
					 double *subBoundary, double *minmaxHeight, int *cols, int *rows)
{
	bool ret = false;

	if(GetImageSize(ImageFilename,Imagesize))
	{
		int i;
		
		D3DPOINT t_pts[8];
		D3DPOINT *t_pts1;
		D2DPOINT *ImageCoord;
		int buffer, null_buffer;

		double minX =  1000000;
		double maxX = -1000000;
		double minY =  1000000;
		double maxY = -1000000;

		t_pts[0].m_X	= subBoundary[0];
		t_pts[1].m_X	= subBoundary[2];
		t_pts[2].m_X	= subBoundary[0];
		t_pts[3].m_X	= subBoundary[2];
		t_pts[4].m_X	= subBoundary[0];
		t_pts[5].m_X	= subBoundary[2];
		t_pts[6].m_X	= subBoundary[0];
		t_pts[7].m_X	= subBoundary[2];
		
		t_pts[0].m_Y	= subBoundary[1];
		t_pts[1].m_Y	= subBoundary[3];
		t_pts[2].m_Y	= subBoundary[1];
		t_pts[3].m_Y	= subBoundary[3];
		t_pts[4].m_Y	= subBoundary[3];
		t_pts[5].m_Y	= subBoundary[1];
		t_pts[6].m_Y	= subBoundary[3];
		t_pts[7].m_Y	= subBoundary[1];

		t_pts[0].m_Z	= minmaxHeight[0];
		t_pts[1].m_Z	= minmaxHeight[0];
		t_pts[2].m_Z	= minmaxHeight[1];
		t_pts[3].m_Z	= minmaxHeight[1];
		t_pts[4].m_Z	= minmaxHeight[0];
		t_pts[5].m_Z	= minmaxHeight[0];
		t_pts[6].m_Z	= minmaxHeight[1];
		t_pts[7].m_Z	= minmaxHeight[1];

        t_pts1			= ps2wgs_3D(transparam,8,t_pts);

		ImageCoord		= GetObjectToImageRPC(RPCs, NumofIAparam, ImageParam, 8, t_pts1);
        
		for(i=0;i<8;i++)
		{
			if(minX > ImageCoord[i].m_X)
				minX	= ImageCoord[i].m_X;
			if(maxX < ImageCoord[i].m_X)
				maxX	= ImageCoord[i].m_X;
			if(minY > ImageCoord[i].m_Y)
				minY	= ImageCoord[i].m_Y;
			if(maxY < ImageCoord[i].m_Y)
				maxY	= ImageCoord[i].m_Y;
            
            //printf("i %d\tImageCoord %f\t%f\n",i,ImageCoord[i].m_X,ImageCoord[i].m_Y);
      	}

		buffer				= 200;
		cols[0]				= (int)(ceil(minX)-buffer);
		cols[1]				= (int)(ceil(maxX)+buffer);
		rows[0]				= (int)(ceil(minY)-buffer);
		rows[1]				= (int)(ceil(maxY)+buffer);

        null_buffer			= 1;
		// Null pixel value remove
		if(cols[0]			<= null_buffer)
			cols[0]			= null_buffer;
		if(rows[0]			<= null_buffer)
			rows[0]			= null_buffer;
		if(cols[0]			> Imagesize->width - null_buffer)
			cols[0]			= Imagesize->width - null_buffer;
		if(rows[0]			> Imagesize->height - null_buffer)
			rows[0]			= Imagesize->height - null_buffer;

		if(cols[1]			<= null_buffer)
			cols[1]			= null_buffer;
		if(rows[1]			<= null_buffer)
			rows[1]			= null_buffer;
		if(cols[1]			> Imagesize->width - null_buffer)
			cols[1]			= Imagesize->width - null_buffer;
		if(rows[1]			> Imagesize->height - null_buffer)
			rows[1]			= Imagesize->height - null_buffer;

        //printf("cols rows %d\t%d\t%d\t%d\n",cols[0],cols[1],rows[0],rows[1]);
        
		free(t_pts1);
		free(ImageCoord);

		ret	= true;
	}

	return ret;
}
		

uint16 *Readtiff(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size, bool check_checktiff)
{
	uint16 *out;
	FILE *bin;
	int check_ftype = 1; // 1 = tif, 2 = bin 
	TIFF *tif = NULL;
	char *ext;
	ext = strrchr(filename,'.');
	
	if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		printf("tif open\n");
		tif	 = TIFFOpen(filename,"r");
		check_ftype = 1;
		printf("tif open end\n");
	}
	else if(!strcmp("bin",ext+1))
	{
		bin	 = fopen(filename,"rb");
		check_ftype = 2;
	}
	
	if(check_ftype == 1 && tif)
	{
		int i,j,row, col, tileW;

		tileW = -1;
		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
		
		printf("tileW = %d\n",tileW);
		
		if(tileW < 0)
		{
			tsize_t scanline;
			tdata_t buf;
			uint16 s,nsamples;
			
			int a;

			// scanline read
			data_size->width	= cols[1] - cols[0];
			data_size->height	= rows[1] - rows[0];

			long int data_length = (long int)data_size->height*(long int)data_size->width;
			
			printf("memory allocation %d\t%d\n",data_size->height,data_size->width);
			
			out				= (uint16*)malloc(sizeof(uint16)*data_length);

			printf("memory allocation %d\t%d\n",data_size->height,data_size->width);
			scanline		= TIFFScanlineSize(tif);

			buf				= _TIFFmalloc(scanline);

			TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);

			for(s =0;s< nsamples;s++)
			{
				for (row=0;row<rows[0];row++)
				{
					if(TIFFReadScanline(tif,buf,row,s) != -1)
					{
						
					}
					else {
						printf("tif error : can't read tiff information\n");
						printf("Error reading image file: %s\n",filename);
						if(!check_checktiff)
							exit(1);
						else
							goto nextmove;
					}

				}
				
				for (row=rows[0];row<rows[1];row++)
				{
					uint16* t_data;
					if(TIFFReadScanline(tif,buf,row,s) != -1)
					{
						if(!check_checktiff)
						{
							t_data = (uint16*)buf;
#pragma omp parallel for shared(cols,rows,out,data_size,t_data) private(a)
							for(a = cols[0];a<cols[1];a++)
							{
								long int pos = (row-rows[0])*data_size->width + (a-cols[0]);
								out[pos] = t_data[a];
							}
						}
					}
					else {
						printf("tif error : can't read tiff information\n");
						printf("Error reading image file: %s\n",filename);
						if(!check_checktiff)
							exit(1);
						else
							goto nextmove;
					}
				}
				
			}
		nextmove:
			_TIFFfree(buf);
		}
		else
		{
            int tileL,count_W,count_L,starttileL,starttileW;
			int start_row,start_col,end_row,end_col;
			tdata_t buf;
			uint16* t_data;

			TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
			TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);

			starttileL		= (int)(rows[0]/tileL);
			start_row		= starttileL*tileL;
			end_row			= ((int)(rows[1]/tileL)+1)*tileL;
			if(end_row > Imagesize->height)
				end_row = Imagesize->height;

			starttileW		= (int)(cols[0]/tileW);
			start_col		= starttileW*tileW;
			end_col			= ((int)(cols[1]/tileW)+1)*tileW;
			if(end_col > Imagesize->width)
				end_col = Imagesize->width;


			cols[0]			= start_col;
			cols[1]			= end_col;
			rows[0]			= start_row;
			rows[1]			= end_row;
			data_size->width = end_col - start_col;
			data_size->height= end_row - start_row;

			long int data_length = (long int)data_size->height*(long int)data_size->width;
			printf("memory allocation %d\t%d\t%li\n",data_size->height,data_size->width,data_length);
			
			out				= (uint16*)malloc(sizeof(uint16)*data_length);

			printf("memory allocation end %d\t%d\n",data_size->height,data_size->width);
			
			buf				= _TIFFmalloc(TIFFTileSize(tif));

			count_L = (int)(data_size->height/tileL);
			count_W = (int)(data_size->width/tileW);
			
			for (row = 0; row < count_L; row ++)
			{
				for (col = 0; col < count_W; col ++)
				{
					if(TIFFReadTile(tif, buf, (col+starttileW)*tileW, (row+starttileL)*tileL, 0,0) != -1)
					{
						if(!check_checktiff)
						{
							t_data = (uint16*)buf;
#pragma omp parallel for shared(tileW,tileL,out,data_size,t_data) private(i,j)
							for (i=0;i<tileL;i++)
							{
								for (j=0;j<tileW;j++)
								{
									long int pos = ((row*tileL) + i)*data_size->width + ((col*tileL) + j);
									out[pos] = t_data[i*tileW + j];
								}
							}
						}
					}
					else {
						printf("tiff error : can't read tiff information\n");
						printf("Error reading image file: %s\n",filename);
						if(!check_checktiff)
							exit(1);
						else
							goto nextmove2;
					}
				}
			}
		nextmove2:
			_TIFFfree(buf);
		}
		TIFFClose(tif);
	}
	else if(check_ftype == 2 && bin)
	{
		int r,c,a;
		data_size->width	= cols[1] - cols[0];
		data_size->height	= rows[1] - rows[0];
		
		long int data_length = data_size->height*data_size->width;
		
		out				= (uint16*)malloc(sizeof(uint16)*data_length);
		
		for(r = rows[0]; r < rows[1] ; r++)
		{
			fseek(bin,sizeof(uint16)*(r*Imagesize->width + cols[0]),SEEK_SET);
			uint16* t_data = (uint16*)malloc(sizeof(uint16)*data_size->width);
			
			fread(t_data,sizeof(uint16),data_size->width,bin);
		
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


float *Readtiff_DEM(char *filename, CSize *Imagesize, int *cols, int *rows, CSize *data_size)
{
	float *out;
	FILE *bin;
	int check_ftype = 1; // 1 = tif, 2 = bin 
	TIFF *tif = NULL;
	char *ext;
	ext = strrchr(filename,'.');
	
	if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		tif	 = TIFFOpen(filename,"r");
		check_ftype = 1;
	}
	else if(!strcmp("bin",ext+1))
	{
		bin	 = fopen(filename,"rb");
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
			data_size->width	= cols[1] - cols[0];
			data_size->height	= rows[1] - rows[0];
			
			out				= (float*)malloc(sizeof(float)*data_size->height*data_size->width);
			
			scanline		= TIFFScanlineSize(tif);
			
			buf				= _TIFFmalloc(scanline);
			
			TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);
			
			for(s =0;s< nsamples;s++)
			{
				for (row=0;row<rows[0];row++)
					TIFFReadScanline(tif,buf,row,s);
				for (row=rows[0];row<rows[1];row++)
				{
					float* t_data;
					TIFFReadScanline(tif,buf,row,s);
					t_data = (float*)buf;
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
			float* t_data;
			
			TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
			TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);
			
			starttileL		= (int)(rows[0]/tileL);
			start_row		= starttileL*tileL;
			end_row			= ((int)(rows[1]/tileL)+1)*tileL;
			if(end_row > Imagesize->height)
				end_row = Imagesize->height;
			
			starttileW		= (int)(cols[0]/tileW);
			start_col		= starttileW*tileW;
			end_col			= ((int)(cols[1]/tileW)+1)*tileW;
			if(end_col > Imagesize->width)
				end_col = Imagesize->width;
			
			
			cols[0]			= start_col;
			cols[1]			= end_col;
			rows[0]			= start_row;
			rows[1]			= end_row;
			
			data_size->width = end_col - start_col;
			data_size->height= end_row - start_row;
			
			out				= (float*)malloc(sizeof(float)*data_size->height*data_size->width);
			
			buf				= _TIFFmalloc(TIFFTileSize(tif));
			
			count_L = (int)(data_size->height/tileL);
			count_W = (int)(data_size->width/tileW);
			
			for (row = 0; row < count_L; row ++)
			{
				for (col = 0; col < count_W; col ++)
				{
					TIFFReadTile(tif, buf, (col+starttileW)*tileW, (row+starttileL)*tileL, 0,0);
					t_data = (float*)buf;
#pragma omp parallel for private(i,j) schedule(guided)
					for (i=0;i<tileL;i++)
					{
						for (j=0;j<tileW;j++)
						{
							out[((row*tileL) + i)*data_size->width + ((col*tileL) + j)] = t_data[i*tileW + j];
						}
					}
				}
			}
			_TIFFfree(buf);
		}
		TIFFClose(tif);
	}
	
	
	return out;
}


void SetSubBoundary(double *Boundary, double subX, double subY, double buffer_area, int col, int row, double *subBoundary)
{
	subBoundary[0]		= Boundary[0] + subX*(col - 1) - buffer_area;
	subBoundary[1]		= Boundary[1] + subY*(row - 1) - buffer_area;
	subBoundary[2]		= Boundary[0] + subX*(col	 ) + buffer_area;
	subBoundary[3]		= Boundary[1] + subY*(row	 ) + buffer_area;
	
	subBoundary[0] =  (int)(floor(subBoundary[0]/8))*8 - 40;
	subBoundary[1] =  (int)(floor(subBoundary[1]/8))*8 - 40;
	subBoundary[2] =  (int)(floor(subBoundary[2]/8))*8 + 40;
	subBoundary[3] =  (int)(floor(subBoundary[3]/8))*8 + 40;
}

D2DPOINT *SetDEMGrid(double *Boundary, double Grid_x, double Grid_y, CSize *Size_2D)
{
	D2DPOINT *GridPT;

	//if(Size_2D->height == 0)
	{
		Size_2D->width	= (int)(ceil((double)(Boundary[2] - Boundary[0])/Grid_x));
		Size_2D->height	= (int)(ceil((double)(Boundary[3] - Boundary[1])/Grid_y));
	}
	
	GridPT	= (D2DPOINT*)malloc(sizeof(D2DPOINT)*Size_2D->height*Size_2D->width);
#pragma omp parallel for schedule(guided)
	for(int row = 0 ; row < Size_2D->height ; row++)
		for(int col = 0; col < Size_2D->width ; col++)
		{
			int index = row*Size_2D->width + col;
			GridPT[index].m_X = Boundary[0] + col*Grid_x;
			GridPT[index].m_Y = Boundary[1] + row*Grid_y;
		}

	return GridPT;
}


void SetHeightWithSeedDEM(TransParam param, UGRID *Grid, double *Boundary, CSize Grid_size, double Grid_set, char *GIMP_path, double *minmaxHeight, double seedDEM_sigma, int IsRA,char* metafilename)
{
	double minX, maxX, minY,maxY,grid_size,a_minX,a_maxX,a_minY,a_maxY;
	CSize seeddem_size;
	char* hdr_path;
	FILE *bin;
	TIFF *tif;
	char save_DEMfile[500];
	
	int check_ftype = 1; // 1 = tif, 2 = raw
	char *ext;
	ext = strrchr(GIMP_path,'.');
	
	if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		check_ftype = 1;
	}
	else if(!strcmp("raw",ext+1))
	{
		check_ftype = 2;
	}
	
	double oriminH = minmaxHeight[0];
	double orimaxH = minmaxHeight[1];
	
	printf("oriminmaxH %f\t%f\n",oriminH,orimaxH);
	
	if(seedDEM_sigma < 10)
		seedDEM_sigma = 10;

	if(IsRA == 1)
	{
		if(seedDEM_sigma < 50)
			seedDEM_sigma = 50;
	}
	
	printf("ttt1 %f\n",seedDEM_sigma);
	
	FILE* pFile_meta;
	pFile_meta	= fopen(metafilename,"r");
	printf("meta file = %s\n",metafilename);
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
	else
	{
		if(check_ftype == 1)
		{
			hdr_path = remove_ext(GIMP_path);
			sprintf(hdr_path,"%s.tfw",hdr_path);
			
			
			FILE *pfile = fopen(hdr_path,"r");
			if(pfile)
			{
				printf("tfw path %s \n",hdr_path);
				TFW_reader_seedDEM(hdr_path, &minX, &maxY, &grid_size);
				GetImageSize(GIMP_path,&seeddem_size);
				fclose(pfile);
			}
			else
			{
				hdr_path = remove_ext(GIMP_path);
				sprintf(hdr_path,"%s.hdr",hdr_path);
				printf("hdr path %s\n",hdr_path);
				seeddem_size  = Envihdr_reader_seedDEM(param,hdr_path, &minX, &maxY, &grid_size);
				fclose(pfile);
			}
		}
		else if(check_ftype == 2)
		{
			hdr_path = remove_ext(GIMP_path);
			sprintf(hdr_path,"%s.hdr",hdr_path);

			printf("hdr path %s\n",hdr_path);
			seeddem_size  = Envihdr_reader_seedDEM(param,hdr_path, &minX, &maxY, &grid_size);
		}
	}
	
	maxX	= minX + grid_size*((double)seeddem_size.width);
	minY	= maxY - grid_size*((double)seeddem_size.height);
	
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

	double total_minH = 999999;
	double total_maxH = -999999;
	printf("%f %f %f %f\n",(double)Boundary[0],(double)Boundary[1],(double)Boundary[2],(double)Boundary[3]);
	printf("%f %f %f %f\n",a_minX, a_maxX, a_minY, a_maxY);
	if ( (a_minX < a_maxX) && (a_minY < a_maxY))
	{
		printf("seeddem cal %f\n",seedDEM_sigma);
		int row,col;
		
		if(check_ftype == 2)
		{
			float *seeddem = NULL;
			bin = fopen(GIMP_path,"rb");
			seeddem = (float*)malloc(sizeof(float)*seeddem_size.width*seeddem_size.height);
			fread(seeddem,sizeof(float),seeddem_size.width*seeddem_size.height,bin);
			
			for (row = 0; row < Grid_size.height; row ++) {
				for (col = 0; col < Grid_size.width; col ++) {
					int index_grid = row*Grid_size.width + col;
					double t_x,t_y;
					int index_seeddem;
					int row_seed, col_seed;
					
					t_x = Boundary[0] + col*Grid_set;
					t_y = Boundary[1] + row*Grid_set;
					col_seed = floor((t_x - minX)/grid_size);
					row_seed = floor((maxY - t_y)/grid_size);
					index_seeddem = row_seed*seeddem_size.width + col_seed;
					if(index_seeddem >= 0 && index_seeddem < seeddem_size.width*seeddem_size.height)
					{
						if(seeddem[index_seeddem] > -1000)
						{
							if(seeddem[index_seeddem] >= oriminH && seeddem[index_seeddem] <= orimaxH)
							{
								if(minmaxHeight[0] > (int)(seeddem[index_seeddem] - seedDEM_sigma - 0.5))
									Grid[index_grid].minHeight = minmaxHeight[0];
								else
									Grid[index_grid].minHeight = (int)(seeddem[index_seeddem] - seedDEM_sigma - 0.5);
								
								if(minmaxHeight[1] < (int)(seeddem[index_seeddem] + seedDEM_sigma + 0.5))
									Grid[index_grid].maxHeight = minmaxHeight[1];
								else
									Grid[index_grid].maxHeight = (int)(seeddem[index_seeddem] + seedDEM_sigma + 0.5);
								
								Grid[index_grid].Height		   = seeddem[index_seeddem];
								
							}
							else
							{
								Grid[index_grid].minHeight		= (int)(minmaxHeight[0] - 0.5);
								Grid[index_grid].maxHeight		= (int)(minmaxHeight[1] + 0.5);
								Grid[index_grid].Height			= -1000.0;
								
								printf("null seed %f\n",seeddem[index_seeddem]);
							}
							
							if(seeddem[index_seeddem] >= oriminH && seeddem[index_seeddem] <= orimaxH)
							{
								if(total_minH > seeddem[index_seeddem] -seedDEM_sigma)
									total_minH = seeddem[index_seeddem] -seedDEM_sigma;
								
								if(total_maxH < seeddem[index_seeddem] +seedDEM_sigma)
									total_maxH = seeddem[index_seeddem] +seedDEM_sigma;
								
								
							}
						}
						else {
							Grid[index_grid].minHeight = -9999;
							Grid[index_grid].maxHeight = -9999;
						}
					}
					else
					{
						//printf("out of seeddem\n");
					}
					
				}
			}
			
			minmaxHeight[0] = total_minH;
			minmaxHeight[1] = total_maxH;
			
			free(seeddem);
		}
		else
		{
			float *seeddem = NULL;
			int cols[2];
			int rows[2];
			CSize data_size;

			CSize *LImagesize = (CSize*)malloc(sizeof(CSize));
			LImagesize->width = seeddem_size.width;
			LImagesize->height = seeddem_size.height;
			
			
			cols[0] = (int)((a_minX - minX)/grid_size + 0.5);
			cols[1] = (int)((a_maxX - minX)/grid_size + 0.5);
			
			rows[0] = (int)((maxY - a_maxY)/grid_size + 0.5);
			rows[1] = (int)((maxY - a_minY)/grid_size + 0.5);
			
			seeddem = Readtiff_DEM(GIMP_path,LImagesize,cols,rows,&data_size);
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
						if(seeddem[index_seeddem] > -1000)
						{
							if(seeddem[index_seeddem] >= oriminH && seeddem[index_seeddem] <= orimaxH)
							{
								if(minmaxHeight[0] > (int)(seeddem[index_seeddem] - seedDEM_sigma - 0.5)) 
									Grid[index_grid].minHeight = minmaxHeight[0];
								else
									Grid[index_grid].minHeight = (int)(seeddem[index_seeddem] - seedDEM_sigma - 0.5);
								
								if(minmaxHeight[1] < (int)(seeddem[index_seeddem] + seedDEM_sigma + 0.5))
									Grid[index_grid].maxHeight = minmaxHeight[1];
								else
									Grid[index_grid].maxHeight = (int)(seeddem[index_seeddem] + seedDEM_sigma + 0.5);
								
								Grid[index_grid].Height		   = seeddem[index_seeddem];
							}
							else
							{
								Grid[index_grid].minHeight		= (int)(minmaxHeight[0] - 0.5);
								Grid[index_grid].maxHeight		= (int)(minmaxHeight[1] + 0.5);
								Grid[index_grid].Height			= -1000.0;
							}
							
							if(seeddem[index_seeddem] >= oriminH && seeddem[index_seeddem] <= orimaxH)
							{
								if(total_minH > seeddem[index_seeddem] -seedDEM_sigma)
									total_minH = seeddem[index_seeddem] -seedDEM_sigma;
								
								if(total_maxH < seeddem[index_seeddem] +seedDEM_sigma)
									total_maxH = seeddem[index_seeddem] +seedDEM_sigma;
							}
						}
						else {
							Grid[index_grid].minHeight = -9999;
							Grid[index_grid].maxHeight = -9999;
						}

					}
					
				}
			}
			
			minmaxHeight[0] = total_minH;
			minmaxHeight[1] = total_maxH;
			
			free(seeddem);
		}
		printf("seeddem end\n");
		printf("%f %f\n",minmaxHeight[0],minmaxHeight[1]);
	}
}


double** OpenXMLFile(char* _filename, double* gsd_r, double* gsd_c, double* gsd, BandInfo* band)
{
    double** out = NULL;
    
    FILE *pFile;
    char temp_str[1000];
    char linestr[1000];
    char linestr1[1000];
    int i;
    char* pos1;
    char* pos2;
    char* token = NULL;
    char* token1 = NULL;
    char* token2 = NULL;
    
    double aa;
    bool band_check = false;
    
    pFile			= fopen(_filename,"r");
    if(pFile)
    {
        out = (double**)malloc(sizeof(double*)*7);
        out[0] = (double*)malloc(sizeof(double)*5);
        out[1] = (double*)malloc(sizeof(double)*5);
        out[6] = (double*)malloc(sizeof(double)*2);
        
        while(!feof(pFile))
        {
            fgets(linestr,sizeof(linestr),pFile);
            strcpy(linestr1,linestr);
            token1 = strstr(linestr,"<");
            token = strtok(token1,">");
            if(strcmp(token,"<ABSCALFACTOR") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                band->abscalfactor			= atof(pos2);
                printf("abscalfactor %f\n",band->abscalfactor);
            }
            
            if(strcmp(token,"<EFFECTIVEBANDWIDTH") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                band->effbw			= atof(pos2);
                printf("effbw %f\n",band->effbw);
            }
            
            if(strcmp(token,"<TDILEVEL") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                band->tdi			= atof(pos2);
                printf("tdi %f\n",band->tdi);
            }
            
            if(strcmp(token,"<MEANCOLLECTEDROWGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd_r			= atof(pos2);
                printf("collect row %f\n",*gsd_r);
            }
            if(strcmp(token,"<MEANCOLLECTEDCOLGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd_c			= atof(pos2);
                printf("collect col %f\n",*gsd_c);
            }
            if(strcmp(token,"<MEANCOLLECTEDGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd			= atof(pos2);
                printf("collect gsd %f\n",*gsd);
            }
            
            if(strcmp(token,"<MEANPRODUCTROWGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd_r			= atof(pos2);
                printf("product row %f\n",*gsd_r);
            }
            if(strcmp(token,"<MEANPRODUCTCOLGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd_c			= atof(pos2);
                printf("product col %f\n",*gsd_c);
            }
            if(strcmp(token,"<MEANPRODUCTGSD") == 0)
            {
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                *gsd			= atof(pos2);
                printf("product gsd %f\n",*gsd);
            }
            
            
            if(strcmp(token,"<ERRBIAS") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[6][0]			= atof(pos2);
                printf("ERRBIAS %f\n",out[6][0]);
            }
            if(strcmp(token,"<ERRRAND") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[6][1]			= atof(pos2);
            }
            if(strcmp(token,"<LINEOFFSET") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][0]			= atof(pos2);
            }
            if(strcmp(token,"<SAMPOFFSET") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][1]			= atof(pos2);
            }
            if(strcmp(token,"<LATOFFSET") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][2]			= atof(pos2);
            }
            if(strcmp(token,"<LONGOFFSET") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][3]			= atof(pos2);
            }
            if(strcmp(token,"<HEIGHTOFFSET") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][4]			= atof(pos2);
            }
            
            
            if(strcmp(token,"<LINESCALE") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][0]			= atof(pos2);
            }
            if(strcmp(token,"<SAMPSCALE") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][1]			= atof(pos2);
            }
            if(strcmp(token,"<LATSCALE") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][2]			= atof(pos2);
            }
            if(strcmp(token,"<LONGSCALE") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][3]			= atof(pos2);
            }
            if(strcmp(token,"<HEIGHTSCALE") == 0)
            {
                //fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr1,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][4]			= atof(pos2);
            }
         
            
            if(strcmp(token,"<LINENUMCOEFList") == 0)
            {
                out[2] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token2 = strtok(pos2," ");
                
                
                i=0;
                while(token2 != NULL && i<20)
                {
                    out[2][i]			= atof(token2);
                    token2 = strtok(NULL," ");
                    
                    //printf("out[2][i] %f\n",out[2][i]);
                    
                    i++;
                    
                    
                }
            }
            if(strcmp(token,"<LINEDENCOEFList") == 0)
            {
                out[3] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token2 = strtok(pos2," ");
                i=0;
                while(token != NULL && i<20)
                {
                    out[3][i]			= atof(token2);
                    token2 = strtok(NULL," ");
                    
                    //printf("out[3][i] %f\n",out[3][i]);
                    
                    i++;
                }
            }
            if(strcmp(token,"<SAMPNUMCOEFList") == 0)
            {
                out[4] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token2 = strtok(pos2," ");
                
                printf("token2 %s\n",token2);
                
                i=0;
                while(token != NULL && i<20)
                {
                    out[4][i]			= atof(token2);
                    token2 = strtok(NULL," ");
                    
                    //printf("out[4][i] %f\n",out[4][i]);
                    
                    i++;
                }
            }
            if(strcmp(token,"<SAMPDENCOEFList") == 0)
            {
                out[5] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token2 = strtok(pos2," ");
                i=0;
                while(token != NULL && i<20)
                {
                    out[5][i]			= atof(token2);
                    token2 = strtok(NULL," ");
                    
                    //printf("out[5][i] %f\n",out[5][i]);
                    
                    i++;
                }
            }
            /*
            fscanf(pFile,"%s",temp_str);
            if(strcmp(temp_str,"<BAND_P>") == 0 && !band_check)
            {
                fgets(temp_str,sizeof(temp_str),pFile);
                bool check_end = false;
                int t_count = 0;
                while(!check_end && t_count < 20)
                {
                    t_count++;
                    fgets(linestr,sizeof(linestr),pFile);
                    strcpy(linestr1,linestr);
                    token1 = strstr(linestr,"<");
                    token = strtok(token1,">");
                    if(strcmp(token,"<ABSCALFACTOR") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        band->abscalfactor			= atof(pos2);
                    }
                    
                    if(strcmp(token,"<EFFECTIVEBANDWIDTH") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        band->effbw			= atof(pos2);
                    }
                    
                    if(strcmp(token,"<TDILEVEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        band->tdi			= atof(pos2);
                        check_end = true;
                    }
                }
                
                if(check_end)
                    band_check = true;
            }
            
            if(strcmp(temp_str,"<IMAGE>") == 0)
            {
                fgets(temp_str,sizeof(temp_str),pFile);
                
                bool check_end = false;
                int t_count = 0;
                while(!check_end && t_count < 30)
                {
                    t_count++;
                    fgets(linestr,sizeof(linestr),pFile);
                    strcpy(linestr1,linestr);
                    token1 = strstr(linestr,"<");
                    token = strtok(token1,">");
                    bool b_rowgsd = false;
                    bool b_colgsd = false;
                    if(strcmp(token,"<MEANCOLLECTEDROWGSD") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd_r			= atof(pos2);
                        b_rowgsd = true;
                    }
                    if(strcmp(token,"<MEANCOLLECTEDCOLGSD") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd_c			= atof(pos2);
                        b_colgsd = true;
                    }
                    if(strcmp(token,"<MEANCOLLECTEDGSD") == 0)
                    {
                        fgets(linestr,sizeof(linestr),pFile);
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd			= atof(pos2);
                    }

                    if(strcmp(token,"<MEANPRODUCTROWGSD") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd_r			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANPRODUCTCOLGSD") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd_c			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANPRODUCTGSD") == 0)
                    {
                        fgets(linestr,sizeof(linestr),pFile);
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        *gsd			= atof(pos2);
                        check_end = true;
                    }
                    
                    if(b_rowgsd && b_colgsd)
                        check_end = true;
                }
            }
            
            if(strcmp(linestr,"<RPB>")==0)
            {
                for(i=0;i<5;i++)
                    fgets(temp_str,sizeof(temp_str),pFile);
                
                out[6] = (double*)malloc(sizeof(double)*2);
                for(i=0;i<2;i++)
                {
                    
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[6][i]			= atof(pos2);
                    printf("Err %f\n",out[6][i]);
                }
                
                out[0] = (double*)malloc(sizeof(double)*5);
                for(i=0;i<5;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][i]			= atof(pos2);
                    
                    printf("param %f\n",out[0][i]);
                }
                aa						= out[0][2];
                out[0][2]				= out[0][3];
                out[0][3]				= aa;
                
                out[1] = (double*)malloc(sizeof(double)*5);
                for(i=0;i<5;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][i]			= atof(pos2);
                }
                
                aa						= out[1][2];
                out[1][2]				= out[1][3];
                out[1][3]				= aa;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                out[2] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token = strtok(pos2," ");
                i=0;
                while(token != NULL)
                {
                    out[2][i]			= atof(token);
                    token = strtok(NULL," ");
                    i++;
                }
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                out[3] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token = strtok(pos2," ");
                i=0;
                while(token != NULL)
                {
                    out[3][i]			= atof(token);
                    token = strtok(NULL," ");
                    i++;
                }
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                out[4] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token = strtok(pos2," ");
                i=0;
                while(token != NULL)
                {
                    out[4][i]			= atof(token);
                    token = strtok(NULL," ");
                    i++;
                }
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                out[5] = (double*)malloc(sizeof(double)*20);
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                token = strtok(pos2," ");
                i=0;
                while(token != NULL)
                {
                    out[5][i]			= atof(token);
                    token = strtok(NULL," ");
                    i++;
                }
            }
             */
        }
        
        aa						= out[0][2];
        out[0][2]				= out[0][3];
        out[0][3]				= aa;
        
        aa						= out[1][2];
        out[1][2]				= out[1][3];
        out[1][3]				= aa;
        
        fclose(pFile);
        
        if(pos1)
            pos1 = NULL;
        if(pos2)
            pos2 = NULL;
        if(token)
            token = NULL;
    }
    return out;
}


double** OpenXMLFile_Pleiades(char* _filename)
{
    double** out;
    
    FILE *pFile;
    char temp_str[1000];
    char linestr[1000];
    int i;
    char* pos1;
    char* pos2;
    char* token = NULL;
    
    double aa;
    
    pFile			= fopen(_filename,"r");
    if(pFile)
    {
        out = (double**)malloc(sizeof(double*)*7);
        while(!feof(pFile))
        {
            fscanf(pFile,"%s",temp_str);
            if(strcmp(temp_str,"<Inverse_Model>") == 0)
            {
                fgets(linestr,sizeof(linestr),pFile);
                
                printf("sample num\n");
                out[4] = (double*)malloc(sizeof(double)*20);
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[4][i]			= atof(pos2);
                    
                    printf("%f\n",out[4][i]);
                }
                
                printf("sample den\n");
                out[5] = (double*)malloc(sizeof(double)*20);
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[5][i]			= atof(pos2);
                    
                    printf("%f\n",out[5][i]);
                }
                
                printf("line num\n");
                out[2] = (double*)malloc(sizeof(double)*20);
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[2][i]			= atof(pos2);
                    
                    printf("%f\n",out[2][i]);
                }
                
                printf("line den\n");
                out[3] = (double*)malloc(sizeof(double)*20);
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[3][i]			= atof(pos2);
                    
                    printf("%f\n",out[3][i]);
                }
                
                out[6] = (double*)malloc(sizeof(double)*2);
                for(i=0;i<2;i++)
                {
                    
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[6][i]			= atof(pos2);
                    
                    printf("%f\n",out[6][i]);
                }
                
                for(i=0;i<14;i++)
                    fgets(linestr,sizeof(linestr),pFile);
                
                out[0] = (double*)malloc(sizeof(double)*5); //offset
                out[1] = (double*)malloc(sizeof(double)*5); //scale
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][2]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][2]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][3]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][3]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][4]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][4]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][1]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][1]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][0]			= atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][0]			= atof(pos2);
                
                for(i=0;i<2;i++)
                {
                    for(int j=0;j<5;j++)
                        printf("%f\n",out[i][j]);
                }
            }
        }
        fclose(pFile);
        
        if(pos1)
            pos1 = NULL;
        if(pos2)
            pos2 = NULL;
    }
    return out;
}

void OpenXMLFile_orientation(char* _filename, ImageInfo *Iinfo)
{
    
    FILE *pFile;
    char temp_str[1000];
    char linestr[1000];
    char linestr1[1000];
    char imagetime[1000];
    char SatID[1000];
    
    int i;
    char* pos1;
    char* pos2;
    char* token = NULL;
    char* token1 = NULL;
    char direction[100];
    
    double dx, dy;
    double MSUNAz, MSUNEl, MSATAz, MSATEl, MIntrackangle, MCrosstrackangle, MOffnadirangle, Cloud;
    double UL[3], UR[3], LR[3], LL[3];
    double angle;
    
    //printf("%s\n",_filename);
    
    pFile           = fopen(_filename,"r");
    if(pFile)
    {
        bool check_br = false;
        bool check_d = false;
        
        while(!feof(pFile) && (!check_br || !check_d))
        {
            fscanf(pFile,"%s",temp_str);
            if(strcmp(temp_str,"<BAND_P>") == 0 && !check_br)
            {
                check_br = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    UL[i]			= atof(pos2);
                    Iinfo->UL[i]    = UL[i];
                    
                }
                //printf("UL %f %f \n",UL[0],UL[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    UR[i]			= atof(pos2);
                    Iinfo->UR[i]    = UR[i];
                }
                //printf("UR %f %f \n",UR[0],UR[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LR[i]			= atof(pos2);
                    Iinfo->LR[i]    = LR[i];
                }
                //printf("LR %f %f \n",LR[0],LR[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LL[i]			= atof(pos2);
                    Iinfo->LL[i]    = LL[i];
                }
                //printf("LL %f %f \n",LL[0],LL[1]);
            }
            
            if(strcmp(temp_str,"<IMAGE>") == 0 && !check_d)
            {
                check_d = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                bool check_end = false;
                int t_count = 0;
                while(!check_end && t_count < 60)
                {
                    t_count++;
                    fgets(linestr,sizeof(linestr),pFile);
                    strcpy(linestr1,linestr);
                    token1 = strstr(linestr,"<");
                    token = strtok(token1,">");
                    
                    if(strcmp(token,"<SATID") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        sprintf(SatID,"%s",pos2);
                    }
                    if(strcmp(token,"<FIRSTLINETIME") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        sprintf(imagetime,"%s",pos2);
                    }
                    if(strcmp(token,"<MEANSUNAZ") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSUNAz			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANSUNEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSUNEl			= atof(pos2);
                    }
                    
                    if(strcmp(token,"<MEANSATAZ") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATAz			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANSATEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATEl			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANINTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MIntrackangle			= atof(pos2);
                    }
                    
                    if(strcmp(token,"<MEANCROSSTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MCrosstrackangle			= atof(pos2);
                    }
                    if(strcmp(token,"<MEANOFFNADIRVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MOffnadirangle			= atof(pos2);
                    }
                    if(strcmp(token,"<CLOUDCOVER") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        Cloud			= atof(pos2);
                        
                        check_end = true;
                    }
                }
            }
        }
        
        fclose(pFile);
        
        Iinfo->Mean_sun_azimuth_angle   = MSUNAz;
        Iinfo->Mean_sun_elevation       = MSUNEl;
        Iinfo->Mean_sat_azimuth_angle   = MSATAz;
        Iinfo->Mean_sat_elevation       = MSATEl;
        Iinfo->Intrack_angle            = MIntrackangle;
        Iinfo->Crosstrack_angle         = MCrosstrackangle;
        Iinfo->Offnadir_angle           = MOffnadirangle;
        Iinfo->cloud                    = Cloud;
        sprintf(Iinfo->imagetime,"%s",imagetime);
        sprintf(Iinfo->SatID,"%s",SatID);
        
        if(pos1)
            pos1 = NULL;
        if(pos2)
            pos2 = NULL;
        if(token)
            token = NULL;
    }
}


void SetDEMBoundary(double** _rpcs, double* _res,TransParam _param, bool _hemisphere, double* _boundary, double* _minmaxheight, CSize* _imagesize, double* _Hinterval)
{
	double minLon = (double) (-1.2 * _rpcs[1][2] + _rpcs[0][2]);
	double maxLon = (double) (1.2 * _rpcs[1][2] + _rpcs[0][2]);
	double minLat = (double) (-1.2 * _rpcs[1][3] + _rpcs[0][3]);
	double maxLat = (double) (1.2 * _rpcs[1][3] + _rpcs[0][3]);
	
	_minmaxheight[0] =	floor(-1.5 * _rpcs[1][4] + _rpcs[0][4]);
	_minmaxheight[1] =	ceil(1.5 * _rpcs[1][4] + _rpcs[0][4]);
	
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
	
	D2DPOINT *lonlat = (D2DPOINT *) malloc(sizeof(D2DPOINT) * 4);
	lonlat[0].m_X = minLon;
	lonlat[0].m_Y = minLat;
	lonlat[1].m_X = minLon;
	lonlat[1].m_Y = maxLat;
	lonlat[2].m_X = maxLon;
	lonlat[2].m_Y = maxLat;
	lonlat[3].m_X = maxLon;
	lonlat[3].m_Y = minLat;
	D2DPOINT *XY = wgs2ps(_param, 4, lonlat);
	
	double minX = min(XY[3].m_X, min(XY[2].m_X, min(XY[0].m_X, XY[1].m_X)));
	double maxX = max(XY[3].m_X, max(XY[2].m_X, max(XY[0].m_X, XY[1].m_X)));
	
	double minY = min(XY[3].m_Y, min(XY[2].m_Y, min(XY[0].m_Y, XY[1].m_Y)));
	double maxY = max(XY[3].m_Y, max(XY[2].m_Y, max(XY[0].m_Y, XY[1].m_Y)));
	
	_boundary[0] =	floor(minX);
	_boundary[1] =	floor(minY);
	_boundary[2] =	ceil(maxX);
	_boundary[3] =	ceil(maxY);
	
	_imagesize->height = (unsigned int) (ceil((_boundary[3] - _boundary[1]) / _res[1]));
	_imagesize->width = (unsigned int) (ceil((_boundary[2] - _boundary[0]) / _res[0]));
	
	free(lonlat);
	free(XY);
}


bool subsetImage(TransParam transparam, uint8 NumofIAparam, double **LRPCs, double *LImageParam, char *LImageFilename, double **RRPCs, double *RImageParam, char *RImageFilename, 
				 double *subBoundary, double *minmaxHeight, D2DPOINT *Lstartpos, D2DPOINT *Rstartpos, char *LsubsetImage, char *RsubsetImage, CSize* Lsubsetsize, CSize* Rsubsetsize, FILE *fid,bool check_checktiff)
{
	bool ret = false;

	CSize *LImagesize = (CSize*)malloc(sizeof(CSize));
	CSize *RImagesize = (CSize*)malloc(sizeof(CSize));

	if(GetImageSize(LImageFilename,LImagesize) && GetImageSize(RImageFilename,RImagesize))
	{
		int Lcols[2], Lrows[2], Rcols[2], Rrows[2];
		if(GetsubareaImage(transparam,NumofIAparam,LRPCs,LImageParam,LImageFilename,LImagesize,subBoundary,minmaxHeight,Lcols,Lrows) &&
		   GetsubareaImage(transparam,NumofIAparam,RRPCs,RImageParam,RImageFilename,RImagesize,subBoundary,minmaxHeight,Rcols,Rrows) )
		{
			uint16 *leftimage, *rightimage;

			printf("read leftimage\n");
			leftimage	= Readtiff(LImageFilename,LImagesize,Lcols,Lrows,Lsubsetsize,check_checktiff);
			if(check_checktiff)
				exit(1);
			
			printf("read rightimage\n");
			rightimage	= Readtiff(RImageFilename,RImagesize,Rcols,Rrows,Rsubsetsize,check_checktiff);
			
			
			
			Lstartpos->m_X	= (double)(Lcols[0]);
			Lstartpos->m_Y	= (double)(Lrows[0]);
			Rstartpos->m_X	= (double)(Rcols[0]);
			Rstartpos->m_Y	= (double)(Rrows[0]);
			

			printf("write subimage left\n");
			if(leftimage)
			{
				FILE *pFile = fopen(LsubsetImage,"wb");

				CSize result_size;

				result_size.width	= Lsubsetsize->width;
				result_size.height	= Lsubsetsize->height;
					
				long int data_length = result_size.height*result_size.width;
				fwrite(leftimage,sizeof(uint16),data_length,pFile);
				fclose(pFile);
				
				free(leftimage);
			}

			printf("write subimage right\n");
			if(rightimage)
			{
				FILE *pFile = fopen(RsubsetImage,"wb");

				CSize result_size;
				

				result_size.width  = Rsubsetsize->width;
				result_size.height = Rsubsetsize->height;

				long int data_length = result_size.height*result_size.width;
				fwrite(rightimage,sizeof(uint16),data_length,pFile);
				fclose(pFile);

				free(rightimage);
			}

			ret		= true;
		}
	}
	free(LImagesize);
	free(RImagesize);

	return ret;
}

D2DPOINT *wgs2ps(TransParam _param, int _numofpts, D2DPOINT *_wgs) {
	int m_NumOfPts = _numofpts;
	
	if(_param.projection == 1)
	{
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
		D2DPOINT *m_sPS;
		
		if (m_NumOfPts > 0) {
			m_sPS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
			double lambda = (double) (_wgs[i].m_X * pm * DegToRad);
			double phi= (double) (_wgs[i].m_Y * pm * DegToRad);
		        
			double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
			double rho = a * m_c * t / t_c;
			
			double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
			m_sPS[i].m_X = (double) (pm * rho * sin(lambda - lambda_0));
			m_sPS[i].m_Y = (double) (-pm * rho * cos(lambda - lambda_0));
		}
	        
		return m_sPS;
	}
	else
	{
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
		int Huso = _param.zone;
		
		D2DPOINT *m_sPS;
		if (m_NumOfPts > 0) {
			m_sPS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
			double Lat = _wgs[i].m_Y;
			double Lon = _wgs[i].m_X;
			double lon = Lon * DegToRad;
			double lat = Lat * DegToRad;
			
			int S = ( ( Huso * 6 ) - 183 );
			double deltaS = lon - ( (double)S*DegToRad) ;
			
			double a = cos(lat)*sin(deltaS);
			double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
			double nu = atan( tan(lat)/cos(deltaS) ) - lat;
			double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
			double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
			double a1 = sin( 2* lat );
			double a2 = a1* ( cos(lat)*cos(lat) );
			double j2 = lat + ( a1 / 2.0 );
			double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
			double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
			double alfa = ( 3 / 4.0 ) * e2cuadrada;
			double beta = ( 5 / 3.0 ) * alfa * alfa;
			double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
			double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
			double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
			double yy = nu * v * ( 1 + ta ) + Bm;
			
			if (yy < 0)
				yy = 9999999 + yy;
			
			m_sPS[i].m_X = xx;
			m_sPS[i].m_Y = yy;
		}
		
		return m_sPS;
	}
}

D2DPOINT wgs2ps_single(TransParam _param, D2DPOINT _wgs) {
	if(_param.projection == 1)
	{
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
		
		D2DPOINT m_sWGS = _wgs;
		D2DPOINT m_sPS;
		
		{
			m_sWGS.m_X = m_sWGS.m_X * pm * DegToRad;
			m_sWGS.m_Y = m_sWGS.m_Y * pm * DegToRad;
			double lambda = m_sWGS.m_X;
			double phi = m_sWGS.m_Y;
			
			double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
			double rho = a * m_c * t / t_c;
			
			double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
			m_sPS.m_X = pm * rho * sin(lambda - lambda_0);
			m_sPS.m_Y = -pm * rho * cos(lambda - lambda_0);
		}
		
		return m_sPS;
	}
	else
	{
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
		int Huso = _param.zone;
		
		D2DPOINT m_sPS;
		
		double Lat = _wgs.m_Y;
		double Lon = _wgs.m_X;
		double lon = Lon * DegToRad;
		double lat = Lat * DegToRad;
		
		int S = ( ( Huso * 6 ) - 183 );
		double deltaS = lon - ( (double)S*DegToRad) ;
		
		double a = cos(lat)*sin(deltaS);
		double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
		double nu = atan( tan(lat)/cos(deltaS) ) - lat;
		double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
		double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
		double a1 = sin( 2* lat );
		double a2 = a1* ( cos(lat)*cos(lat) );
		double j2 = lat + ( a1 / 2.0 );
		double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
		double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
		double alfa = ( 3 / 4.0 ) * e2cuadrada;
		double beta = ( 5 / 3.0 ) * alfa * alfa;
		double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
		double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
		double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
		double yy = nu * v * ( 1 + ta ) + Bm;
		
		
		if (yy < 0)
			yy = 9999999 + yy;
		
		m_sPS.m_X = xx;
		m_sPS.m_Y = yy;
		
		return m_sPS;
	}
	
}

D3DPOINT *wgs2ps_3D(TransParam _param, int _numofpts, D3DPOINT *_wgs) {
	int m_NumOfPts = _numofpts;
	
	if(_param.projection == 1)
	{
		
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
		D3DPOINT *m_sPS;
		
		if (m_NumOfPts > 0) {
			m_sPS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
		  double lambda = _wgs[i].m_X * pm * DegToRad;
		  double phi = _wgs[i].m_Y * pm * DegToRad;
			
			double t = tan(PI / 4.0 - phi / 2.0) / pow((1.0 - e * sin(phi)) / (1.0 + e * sin(phi)), e / 2.0);
			double rho = a * m_c * t / t_c;
			
			double m = cos(phi) / sqrt(1.0 - pow(e, 2) * pow(sin(phi), 2));
			m_sPS[i].m_X = pm * rho * sin(lambda - lambda_0);
			m_sPS[i].m_Y = -pm * rho * cos(lambda - lambda_0);
			m_sPS[i].m_Z = _wgs[i].m_Z;
		}
	        
		return m_sPS;
	}
	else
	{
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
		int Huso = _param.zone;
		
		D3DPOINT *m_sPS;
		
		if (m_NumOfPts > 0) {
			m_sPS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
			double Lat = _wgs[i].m_Y;
			double Lon = _wgs[i].m_X;
			double lon = Lon * DegToRad;
			double lat = Lat * DegToRad;
			
			int S = ( ( Huso * 6 ) - 183 );
			double deltaS = lon - ( (double)S*DegToRad) ;
			
			double a = cos(lat)*sin(deltaS);
			double epsilon = 0.5 * log( ( 1 +  a)/ ( 1 - a ) );
			double nu = atan( tan(lat)/cos(deltaS) ) - lat;
			double v = ( c /  sqrt( 1 + ( e2cuadrada* cos(lat)*cos(lat)) ) )* 0.9996;
			double ta = ( e2cuadrada/ 2.0 ) * (epsilon*epsilon)* ( cos(lat)*cos(lat) );
			double a1 = sin( 2* lat );
			double a2 = a1* ( cos(lat)*cos(lat) );
			double j2 = lat + ( a1 / 2.0 );
			double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
			double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
			double alfa = ( 3 / 4.0 ) * e2cuadrada;
			double beta = ( 5 / 3.0 ) * alfa * alfa;
			double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
			double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
			double xx = epsilon * v * ( 1 + ( ta / 3.0 ) ) + 500000;
			double yy = nu * v * ( 1 + ta ) + Bm;
			
			
			if (yy < 0)
				yy = 9999999 + yy;
			
			m_sPS[i].m_X = xx;
			m_sPS[i].m_Y = yy;
			m_sPS[i].m_Z = _wgs[i].m_Z;
		}
		
		return m_sPS;
		
	}
	
}

D2DPOINT *ps2wgs(TransParam _param, int _numofpts, D2DPOINT *_ps) {
	int m_NumOfPts = _numofpts;
	
	if(_param.projection == 1)
	{
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
	        
		D2DPOINT *m_sWGS;
		
		if (m_NumOfPts > 0) {
			m_sWGS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
		double e2 = e * e;
		double e4 = e2 * e2;
		double e6 = e2 * e4;
		double e8 = e4 * e4;
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
		  double x = _ps[i].m_X * pm;
		  double y = _ps[i].m_Y * pm;
			
		  double rho = sqrt(pow(x, 2) + pow(y, 2));
			double t = rho * t_c / (a * m_c);
			
			double chi = PI / 2 - 2 * atan(t);
			double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
				(7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
			
			double lambda = lambda_0 + atan2(x, -y);
			phi = pm * phi;
			lambda = pm * lambda;
			
			if (lambda > PI) {
				lambda = lambda - 2 * PI;
			} else if (lambda < -PI) {
				lambda = lambda + 2 * PI;
			}
			
			m_sWGS[i].m_Y = RadToDeg * phi;
			m_sWGS[i].m_X = RadToDeg * lambda;
		}
	        
		return m_sWGS;
	}
	else
	{
		int hemis = _param.pm; //1 = north, -1 = south
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
	        
		D2DPOINT *m_sWGS;
		
		if (m_NumOfPts > 0) {
			m_sWGS = (D2DPOINT *) malloc(sizeof(D2DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
			
			double x = _ps[i].m_X;
			double y = _ps[i].m_Y;
			
			double X = x - 500000;
			double Y = y;
			if (hemis < 0)
				Y = Y - 10000000;
			
			int S = ( ( _param.zone* 6 ) - 183 );
			double lat =  Y / ( 6366197.724 * 0.9996 );
			double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
			double a = X / v;
			
			
			double a1 = sin( 2* lat );
			double a2 = a1* ( cos(lat)*cos(lat) );
			double j2 = lat + ( a1 / 2.0 );
			double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
			double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
			double alfa = ( 3 / 4.0 ) * e2cuadrada;
			double beta = ( 5 / 3.0 ) * alfa * alfa;
			double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
			double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
			double b = ( Y - Bm ) / v;
			double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
			double Eps = a * ( 1 - ( Epsi / 3.0 ) );
			double nab = ( b * ( 1 - Epsi ) ) + lat;
			double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
			double Delt = atan(senoheps / (cos(nab) ) );
			double TaO = atan(cos(Delt) * tan(nab));
			double longitude = (Delt * (180/PI) ) + S;
			
			double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
			
			m_sWGS[i].m_Y = latitude;
			m_sWGS[i].m_X = longitude;
		}
		return m_sWGS;
		
	}
	
}

D2DPOINT ps2wgs_single(TransParam _param, D2DPOINT _ps) {
	if(_param.projection == 1)
	{
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
		
		D2DPOINT m_sWGS;
		D2DPOINT m_sPS = _ps;
		
		{
			double e2 = e * e;
			double e4 = e2 * e2;
			double e6 = e2 * e4;
			double e8 = e4 * e4;
			
			m_sPS.m_X = m_sPS.m_X * pm;
			m_sPS.m_Y = m_sPS.m_Y * pm;
			
			double rho = sqrt(pow(m_sPS.m_X, 2) + pow(m_sPS.m_Y, 2));
			double t = rho * t_c / (a * m_c);
			
			double chi = PI / 2 - 2 * atan(t);
			double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
				(7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
			
			double lambda = lambda_0 + atan2(m_sPS.m_X, -m_sPS.m_Y);
			phi = pm * phi;
			lambda = pm * lambda;
			
			if (lambda > PI) {
				lambda = lambda - 2 * PI;
			} else if (lambda < -PI) {
				lambda = lambda + 2 * PI;
			}
			
			m_sWGS.m_Y = RadToDeg * phi;
			m_sWGS.m_X = RadToDeg * lambda;
		}
		
		return m_sWGS;
	}
	else
	{
		int hemis = _param.pm; //1 = north, -1 = south
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
		
		D2DPOINT m_sWGS;
		D2DPOINT m_sPS = _ps;
		
		double x = m_sPS.m_X;
		double y = m_sPS.m_Y;
		
		double X = x - 500000;
		double Y = y;
		if (hemis < 0)
			Y = Y - 10000000;
		
		int S = ( ( _param.zone* 6 ) - 183 );
		double lat =  Y / ( 6366197.724 * 0.9996 );
		double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
		double a = X / v;
		
		
		double a1 = sin( 2* lat );
		double a2 = a1* ( cos(lat)*cos(lat) );
		double j2 = lat + ( a1 / 2.0 );
		double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
		double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
		double alfa = ( 3 / 4.0 ) * e2cuadrada;
		double beta = ( 5 / 3.0 ) * alfa * alfa;
		double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
		double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
		double b = ( Y - Bm ) / v;
		double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
		double Eps = a * ( 1 - ( Epsi / 3.0 ) );
		double nab = ( b * ( 1 - Epsi ) ) + lat;
		double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
		double Delt = atan(senoheps / (cos(nab) ) );
		double TaO = atan(cos(Delt) * tan(nab));
		double longitude = (Delt * (180/PI) ) + S;
		
		double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
		
		m_sWGS.m_Y = latitude;
		m_sWGS.m_X = longitude;
		
		return m_sWGS;
	}
	
}

D3DPOINT *ps2wgs_3D(TransParam _param, int _numofpts, D3DPOINT *_ps) {
	int m_NumOfPts = _numofpts;
	if(_param.projection == 1)
	{
		bool m_bHemisphere = _param.bHemisphere;
		double a = _param.a;
		double e = _param.e;
		double phi_c = _param.phi_c;
		double lambda_0 = _param.lambda_0;
		
		int pm = _param.pm;
		double t_c = _param.t_c;
		double m_c = _param.m_c;
	        
		D3DPOINT *m_sWGS;
		
		if (m_NumOfPts > 0) {
			m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
		double e2 = e * e;
		double e4 = e2 * e2;
		double e6 = e2 * e4;
		double e8 = e4 * e4;
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
		  double x = _ps[i].m_X * pm;
		  double y = _ps[i].m_Y * pm;
			
			double rho = sqrt(pow(x, 2) + pow(y, 2));
			double t = rho * t_c / (a * m_c);
			
			double chi = PI / 2 - 2 * atan(t);
			double phi = chi + (e2 / 2 + 5 * e4 / 24 + e6 / 12 + 13 * e8 / 360) * sin(2 * chi) + (7 * e4 / 48 + 29 * e6 / 240 + 811 * e8 / 11520) * sin(4 * chi) +
				(7 * e6 / 120 + 81 * e8 / 1120) * sin(6 * chi) + (4279 * e8 / 161280) * sin(8 * chi);
			
			double lambda = lambda_0 + atan2(x, -y);
			phi = pm * phi;
			lambda = pm * lambda;
			if (lambda > PI) {
				lambda = lambda - 2 * PI;
			} else if (lambda < -PI) {
				lambda = lambda + 2 * PI;
			}
			
			m_sWGS[i].m_Y = RadToDeg * phi;
			m_sWGS[i].m_X = RadToDeg * lambda;
			m_sWGS[i].m_Z = _ps[i].m_Z;
		}
	        
		return m_sWGS;
	}
	else
	{
		int hemis = _param.pm; //1 = north, -1 = south
		double sa = _param.sa;
		double sb = _param.sb;
		double e2 = _param.e2;
		double e2cuadrada = _param.e2cuadrada;
		double c = _param.c;
	        
		D3DPOINT *m_sWGS;
		
		if (m_NumOfPts > 0) {
			m_sWGS = (D3DPOINT *) malloc(sizeof(D3DPOINT) * m_NumOfPts);
		} else {
			return false;
		}
		
		
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < m_NumOfPts; i++) {
			
			double x = _ps[i].m_X;
			double y = _ps[i].m_Y;
			
			double X = x - 500000;
			double Y = y;
			if (hemis < 0)
				Y = Y - 10000000;
			
			int S = ( ( _param.zone* 6 ) - 183 ); 
			double lat =  Y / ( 6366197.724 * 0.9996 );									   
			double v = ( c / ( sqrt( 1.0 + ( e2cuadrada * ( cos(lat)*cos(lat) ) ) ) ) ) * 0.9996;
			double a = X / v;
			
			
			double a1 = sin( 2* lat );
			double a2 = a1* ( cos(lat)*cos(lat) );
			double j2 = lat + ( a1 / 2.0 );
			double j4 = ( ( 3 * j2 ) + a2 ) / 4.0;
			double j6 = ( ( 5 * j4 ) + ( a2 * ( cos(lat)*cos(lat) )) ) / 3;
			double alfa = ( 3 / 4.0 ) * e2cuadrada;
			double beta = ( 5 / 3.0 ) * alfa * alfa;
			double gama = ( 35 / 27.0 ) * alfa * alfa * alfa;
			double Bm = 0.9996 * c * ( lat - alfa * j2 + beta * j4 - gama * j6 );
			double b = ( Y - Bm ) / v;
			double Epsi = ( ( e2cuadrada * a *a ) / 2.0 ) * ( cos(lat)*cos(lat) );
			double Eps = a * ( 1 - ( Epsi / 3.0 ) );
			double nab = ( b * ( 1 - Epsi ) ) + lat;
			double senoheps = ( exp(Eps) - exp(-Eps) ) / 2;
			double Delt = atan(senoheps / (cos(nab) ) );
			double TaO = atan(cos(Delt) * tan(nab));
			double longitude = (Delt * (180/PI) ) + S;
			
			double latitude = ( lat + ( 1 + e2cuadrada * (cos(lat)*cos(lat)) - ( 3/2.0 )* e2cuadrada * sin(lat) * cos(lat) * ( TaO - lat ) )* ( TaO - lat ) ) * (180/PI);
			
			m_sWGS[i].m_Y = latitude;
			m_sWGS[i].m_X = longitude;
			m_sWGS[i].m_Z = _ps[i].m_Z;
		}
		return m_sWGS;
	}
	
}

D2DPOINT* GetObjectToImageRPC(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP)
{
	D2DPOINT *IP;
	
	IP		= (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);
	
#pragma omp parallel for schedule(guided)
	for(int i=0;i<_numofpts;i++)
	{
		double L, P, H, Line, Samp, deltaP, deltaR;
		double Coeff[4];

		L		= (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
		P		= (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
		H		= (_GP[i].m_Z - _rpc[0][4])/_rpc[1][4];
        
        //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        
        if(L < -10.0 || L > 10.0)
        {
            if(_GP[i].m_X > 0)
                _GP[i].m_X = _GP[i].m_X - 360;
            else
                _GP[i].m_X = _GP[i].m_X + 360;
            
            L		= (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
            
        }
        if(P < -10.0 || P > 10.0)
        {
            if(_GP[i].m_Y > 0)
                _GP[i].m_Y = _GP[i].m_Y - 360;
            else
                _GP[i].m_Y = _GP[i].m_Y + 360;
            
            P		= (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
            
            //printf("original %f\t%f\t%f\tL P H %f\t%f\t%f\n",_GP[i].m_X,_GP[i].m_Y,_GP[i].m_Z,L,P,H);
        }
        
        //printf("L P H %f\t%f\t%f\n",L,P,H);
        
		for(int j=0;j<4;j++)
		{
			Coeff[j]	= _rpc[j+2][0]*1.0			+ _rpc[j+2][1]*L			+ _rpc[j+2][2]*P
				+ _rpc[j+2][3]*H			+ _rpc[j+2][4]*L*P			+ _rpc[j+2][5]*L*H
				+ _rpc[j+2][6]*P*H			+ _rpc[j+2][7]*L*L			+ _rpc[j+2][8]*P*P
				+ _rpc[j+2][9]*H*H			+ _rpc[j+2][10]*(P*L)*H		+ _rpc[j+2][11]*(L*L)*L
				+ _rpc[j+2][12]*(L*P)*P		+ _rpc[j+2][13]*(L*H)*H		+ _rpc[j+2][14]*(L*L)*P
				+ _rpc[j+2][15]*(P*P)*P		+ _rpc[j+2][16]*(P*H)*H		+ _rpc[j+2][17]*(L*L)*H
				+ _rpc[j+2][18]*(P*P)*H		+ _rpc[j+2][19]*(H*H)*H;
		}

		Line	 = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
		Samp	 = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
		deltaP		= _imageparam[0];
		deltaR		= _imageparam[1];

		IP[i].m_Y		= deltaP + Line;
		IP[i].m_X		= deltaR + Samp;
        
        if(IP[i].m_Y < 0)
            IP[i].m_Y = 0;
        
        if(IP[i].m_Y > _rpc[0][0] + _rpc[1][0]*1.2)
            IP[i].m_Y = _rpc[0][0] + _rpc[1][0]*1.2;
        
        if(IP[i].m_X < 0)
            IP[i].m_X = 0;
        
        if(IP[i].m_X > _rpc[0][1] + _rpc[1][1]*1.2)
            IP[i].m_X = _rpc[0][1] + _rpc[1][1]*1.2;
	}

	return IP;
}

D2DPOINT GetObjectToImageRPC_single(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP)
{
	D2DPOINT IP;
	
	int j;
	

	double L, P, H, Line, Samp, deltaP, deltaR;
	double Coeff[4];

	L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
	P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
	H		= (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
    if(L < -10.0 || L > 10.0)
    {
        if(_GP.m_X > 0)
            _GP.m_X = _GP.m_X - 360;
        else
            _GP.m_X = _GP.m_X + 360;
        
        L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    }
    
    if(P < -10.0 || P > 10.0)
    {
        if(_GP.m_Y > 0)
            _GP.m_Y = _GP.m_Y - 360;
        else
            _GP.m_Y = _GP.m_Y + 360;
        
        P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    }
    
	for(j=0;j<4;j++)
	{
		Coeff[j]	= _rpc[j+2][0]*1.0			+ _rpc[j+2][1]*L			+ _rpc[j+2][2]*P
			+ _rpc[j+2][3]*H			+ _rpc[j+2][4]*L*P			+ _rpc[j+2][5]*L*H
			+ _rpc[j+2][6]*P*H			+ _rpc[j+2][7]*L*L			+ _rpc[j+2][8]*P*P
			+ _rpc[j+2][9]*H*H			+ _rpc[j+2][10]*(P*L)*H		+ _rpc[j+2][11]*(L*L)*L
			+ _rpc[j+2][12]*(L*P)*P		+ _rpc[j+2][13]*(L*H)*H		+ _rpc[j+2][14]*(L*L)*P
			+ _rpc[j+2][15]*(P*P)*P		+ _rpc[j+2][16]*(P*H)*H		+ _rpc[j+2][17]*(L*L)*H
			+ _rpc[j+2][18]*(P*P)*H		+ _rpc[j+2][19]*(H*H)*H;
	}

	Line	 = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
	Samp	 = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample

	deltaP		= _imageparam[0];
	deltaR		= _imageparam[1];

	IP.m_Y		= deltaP + Line;
	IP.m_X		= deltaR + Samp;

    if(IP.m_Y < 0)
        IP.m_Y = 0;
    
    if(IP.m_Y > _rpc[0][0] + _rpc[1][0]*1.2)
        IP.m_Y = _rpc[0][0] + _rpc[1][0]*1.2;
    
    if(IP.m_X < 0)
        IP.m_X = 0;
    
    if(IP.m_X > _rpc[0][1] + _rpc[1][1]*1.2)
        IP.m_X = _rpc[0][1] + _rpc[1][1]*1.2;
    
	return IP;
}

D2DPOINT GetObjectToImageRPC_single_mpp(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP)
{
    D2DPOINT IP;
    
    int j;
    
    
    double L, P, H, Line, Samp, deltaP, deltaR;
    double Coeff[4];
    
    L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    H		= (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
    if(L < -10.0 || L > 10.0)
    {
        if(_GP.m_X > 0)
            _GP.m_X = _GP.m_X - 360;
        else
            _GP.m_X = _GP.m_X + 360;
        
        L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    }
    
    if(P < -10.0 || P > 10.0)
    {
        if(_GP.m_Y > 0)
            _GP.m_Y = _GP.m_Y - 360;
        else
            _GP.m_Y = _GP.m_Y + 360;
        
        P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    }
    
    for(j=0;j<4;j++)
    {
        Coeff[j]	= _rpc[j+2][0]*1.0			+ _rpc[j+2][1]*L			+ _rpc[j+2][2]*P
        + _rpc[j+2][3]*H			+ _rpc[j+2][4]*L*P			+ _rpc[j+2][5]*L*H
        + _rpc[j+2][6]*P*H			+ _rpc[j+2][7]*L*L			+ _rpc[j+2][8]*P*P
        + _rpc[j+2][9]*H*H			+ _rpc[j+2][10]*(P*L)*H		+ _rpc[j+2][11]*(L*L)*L
        + _rpc[j+2][12]*(L*P)*P		+ _rpc[j+2][13]*(L*H)*H		+ _rpc[j+2][14]*(L*L)*P
        + _rpc[j+2][15]*(P*P)*P		+ _rpc[j+2][16]*(P*H)*H		+ _rpc[j+2][17]*(L*L)*H
        + _rpc[j+2][18]*(P*P)*H		+ _rpc[j+2][19]*(H*H)*H;
    }
    
    Line	 = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
    Samp	 = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
    
    deltaP		= _imageparam[0];
    deltaR		= _imageparam[1];
    
    IP.m_Y		= deltaP + Line;
    IP.m_X		= deltaR + Samp;
    
    return IP;
}

void Preprocessing(char *save_path,char *Lsubsetfile, char *Rsubsetfile, uint8 py_level, CSize *Lsubsetsize, CSize *Rsubsetsize, CSize *data_size_l, CSize *data_size_r, FILE *fid)
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
			pFile_raw	= fopen(Lsubsetfile,"rb");
			filename_py		= GetFileName(Lsubsetfile);
			filename_py		= remove_ext(filename_py);
		}
		else
		{
			pFile_raw	= fopen(Rsubsetfile,"rb");
			filename_py		= GetFileName(Rsubsetfile);
			filename_py		= remove_ext(filename_py);
		}

		sprintf(t_str,"%s/%s_py_5.raw",save_path,filename_py);
		pFile_check_file	= fopen(t_str,"rb");
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
				filename_py		= GetFileName(Lsubsetfile);
				filename_py		= remove_ext(filename_py);
				for(i=0;i<py_level+1;i++)
					data_size[i]		= data_size_l[i];
				
			}
			else
			{
				fprintf(fid,"right\nlevel = 0\t width = %d\theight = %d\n",data_size_r[0].width,data_size_r[0].height);
				filename_py		= GetFileName(Rsubsetfile);
				filename_py		= remove_ext(filename_py);
				for(i=0;i<py_level+1;i++)
					data_size[i]		= data_size_r[i];
			}

			long int data_length = (long int)data_size[0].height*(long int)data_size[0].width;
			pyimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
			magimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
			dirimg[0] = (int16*)malloc(sizeof(int16)*data_length);
			oriimg[0] = (uint8*)malloc(sizeof(uint8)*data_length);

			fread(pyimg[0],sizeof(uint16),data_length,pFile_raw);

			MakeSobelMagnitudeImage(data_size[0],pyimg[0],magimg[0],dirimg[0]);
			Orientation(data_size[0],magimg[0],dirimg[0],15,oriimg[0]);
			
			sprintf(t_str,"%s/%s_py_0.raw",save_path,filename_py);
			pFile	= fopen(t_str,"wb");
			fwrite(pyimg[0],sizeof(uint16),data_length,pFile);
			fclose(pFile);
			
			sprintf(t_str,"%s/%s_py_0_mag.raw",save_path,filename_py);
			pFile	= fopen(t_str,"wb");
			fwrite(magimg[0],sizeof(uint16),data_length,pFile);
			fclose(pFile);
			free(magimg[0]);

			free(dirimg[0]);

			sprintf(t_str,"%s/%s_py_0_ori.raw",save_path,filename_py);
			pFile	= fopen(t_str,"wb");
			fwrite(oriimg[0],sizeof(uint8),data_length,pFile);
			fclose(pFile);
			free(oriimg[0]);

			for(i=0;i<py_level;i++)
			{
				pyimg[i+1] = CreateImagePyramid(pyimg[i],data_size[i],9,(double)(1.5));
				free(pyimg[i]);

				if(count == 0)
				{
					data_size_l[i+1]	= data_size[i+1];
					fprintf(fid,"level = %d\t width = %d\theight = %d\n",i+1,data_size_l[i+1].width,data_size_l[i+1].height);
				}
				else
				{
					data_size_r[i+1]	= data_size[i+1];
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
				pFile	= fopen(t_str,"wb");
				fwrite(pyimg[i+1],sizeof(uint16),data_length_array[i],pFile);
				fclose(pFile);
				if(i == py_level-1)
					free(pyimg[i+1]);

				sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,i+1);
				pFile	= fopen(t_str,"wb");
				fwrite(magimg[i+1],sizeof(uint16),data_length_array[i],pFile);
				fclose(pFile);
				free(magimg[i+1]);

				free(dirimg[i+1]);

				sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,i+1);
				pFile	= fopen(t_str,"wb");
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

uint16* LoadPyramidImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level)
{
	uint16 *out = (uint16*)malloc(sizeof(uint16));
	FILE *pFile;
	char *filename_py;
	char t_str[500];

	filename_py		= GetFileName(subsetfile);
	filename_py		= remove_ext(filename_py);

	sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,py_level);
	pFile			= fopen(t_str,"rb");
	if(pFile)
	{
		out			= (uint16*)malloc(sizeof(uint16)*data_size.height*data_size.width);
		fread(out,sizeof(uint16),data_size.height*data_size.width,pFile);
	}
	fclose(pFile);
	return out;
}

uint8* LoadPyramidOriImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level)
{
	uint8 *out = (uint8*)malloc(sizeof(uint8));
	FILE *pFile;
	char *filename_py;
	char t_str[500];

	filename_py		= GetFileName(subsetfile);
	filename_py		= remove_ext(filename_py);
	sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,py_level);
	pFile			= fopen(t_str,"rb");
	if(pFile)
	{
		out		= (uint8*)malloc(sizeof(uint8)*data_size.height*data_size.width);
		fread(out,sizeof(uint8),data_size.height*data_size.width,pFile);
	}
	fclose(pFile);

	return out;
}


uint16* LoadPyramidMagImages(char *save_path,char *subsetfile, CSize data_size, uint8 py_level, double *val, double *avg)
{
	uint16 *out = (uint16*)malloc(sizeof(uint16));
	FILE *pFile;
	char *filename_py;
	char t_str[500];
	
	filename_py		= GetFileName(subsetfile);
	filename_py		= remove_ext(filename_py);
	sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,py_level);
	pFile			= fopen(t_str,"rb");
	if(pFile)
	{
		out		= (uint16*)malloc(sizeof(uint16)*data_size.height*data_size.width);
		fread(out,sizeof(uint16),data_size.height*data_size.width,pFile);
	}
	fclose(pFile);
	
	
	double sum2 = 0;
	double sum = 0;
	int count = 0;
	double residual = 0;
	for(int i=0;i<data_size.height;i++)
	{
		for(int j=0;j<data_size.width;j++)
		{
			sum += out[i*data_size.width + j];
			sum2 += out[i*data_size.width + j]*out[i*data_size.width + j];
			count++;
		}
	}
	
	for(int i=0;i<data_size.height;i++)
	{
		for(int j=0;j<data_size.width;j++)
		{
			residual += (out[i*data_size.width + j] - *avg)*(out[i*data_size.width + j] - *avg);
		}
	}
	
	*val		  = sqrt((double)residual/count);
	*avg		  = sum/count;
	
	return out;
}


uint16* CreateImagePyramid(uint16* _input, CSize _img_size, int _filter_size, double _sigma)
{
	//_filter_size = 7, sigma = 1.6
	double sigma = _sigma;
	double temp,scale;
	double sum = 0;
	double** GaussianFilter;
	CSize result_size;
	uint16* result_img;

	GaussianFilter = (double**)malloc(sizeof(double*)*_filter_size);
	for(int i=0;i<_filter_size;i++)
		GaussianFilter[i] = (double*)malloc(sizeof(double)*_filter_size);

	
	result_size.width = _img_size.width/2;
	result_size.height = _img_size.height/2;
	scale=sqrt(2*PI)*sigma;
	
	result_img = (uint16*)malloc(sizeof(uint16)*result_size.height*result_size.width);

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
		}
	}

#pragma omp parallel for private(temp) schedule(guided)
	for(long int r=0;r<result_size.height;r++)
	{
		for(long int c=0;c<result_size.width;c++)
		{
			temp = 0;

			for(int l=0;l<_filter_size;l++)
			{
				for(int k=0;k<_filter_size;k++)
				{
					//r'->2r+m, c'->2c+n
					if( (2*r + l-(int)(_filter_size/2)) >= 0 && (2*c + k-(int)(_filter_size/2)) >= 0 && 
						(2*r + l-(int)(_filter_size/2)) < _img_size.height && (2*c + k-(int)(_filter_size/2)) < _img_size.width)
					{
						temp += GaussianFilter[l][k]*_input[(2*r + l-(int)(_filter_size/2))*_img_size.width +(2*c + k-(int)(_filter_size/2))];
						
					}
				}
			}

			result_img[r*result_size.width + c] = (uint16)temp;
		}
	}
	
	for(int i=0;i<_filter_size;i++)
		if(GaussianFilter[i])
			free(GaussianFilter[i]);

	if(GaussianFilter)
		free(GaussianFilter);
	

	return result_img;
}


void MakeSobelMagnitudeImage(CSize _img_size, uint16* _src_image, uint16* _dist_mag_image, /*int16* _gx, int16* _gy,*/ int16* _dir)
{
	int hy[3][3] , hx[3][3];
	
	//sobel mask
	hx[0][0]=-1;			   //  -1  0  1
	hx[0][1]=0;				   //  -2  0  2
	hx[0][2]=1;				   //  -1  0  1 
	hx[1][0]=-2;
	hx[1][1]=0;
	hx[1][2]=2;
	hx[2][0]=-1;
	hx[2][1]=0;
	hx[2][2]=1;

	hy[0][0]=1;				   //  1  2	 1
	hy[0][1]=2;				   //  0  0	 0
	hy[0][2]=1;				   // -1 -2 -1 
	hy[1][0]=0;
	hy[1][1]=0;
	hy[1][2]=0;
	hy[2][0]=-1;
	hy[2][1]=-2;
	hy[2][2]=-1;

	
	
#pragma omp parallel for
	for(long int i=0;i<_img_size.height;i++)
	{
		int16 _gx, _gy;
		for(long int j=0;j<_img_size.width;j++)
		{
			if(j == 0 || j == _img_size.width-1 || i == 0 || i == _img_size.height-1)
			{
				_gx = 0;
				_gy = 0;
				_dist_mag_image[i*_img_size.width + j] = 0;
				_dir[i*_img_size.width + j] = 0;
			}
			else
			{
				double temp_ptr_X=((hx[0][0]*_src_image[(i-1)*_img_size.width + (j-1)]+hx[0][1]*_src_image[(i-1)*_img_size.width + j]+hx[0][2]*_src_image[(i-1)*_img_size.width + (j+1)]+
									hx[1][0]*_src_image[	i*_img_size.width + (j-1)]+hx[1][1]*_src_image[	   i*_img_size.width + j]+hx[1][2]*_src_image[	  i*_img_size.width + (j+1)]+
									hx[2][0]*_src_image[(i+1)*_img_size.width + (j-1)]+hx[2][1]*_src_image[(i+1)*_img_size.width + j]+hx[2][2]*_src_image[(i+1)*_img_size.width + (j+1)]));

				double temp_ptr_Y=((hy[0][0]*_src_image[(i-1)*_img_size.width + (j-1)]+hy[0][1]*_src_image[(i-1)*_img_size.width + j]+hy[0][2]*_src_image[(i-1)*_img_size.width + (j+1)]+
									hy[1][0]*_src_image[	i*_img_size.width + (j-1)]+hy[1][1]*_src_image[	   i*_img_size.width + j]+hy[1][2]*_src_image[	  i*_img_size.width + (j+1)]+
									hy[2][0]*_src_image[(i+1)*_img_size.width + (j-1)]+hy[2][1]*_src_image[(i+1)*_img_size.width + j]+hy[2][2]*_src_image[(i+1)*_img_size.width + (j+1)]));

				double temp_ptr=sqrt(temp_ptr_X*temp_ptr_X + temp_ptr_Y*temp_ptr_Y);

				_gx = (int16)(temp_ptr_X + 0.5);
				_gy = (int16)(temp_ptr_Y + 0.5);							
				_dir[i*_img_size.width + j]= (int16)(floor(atan2(temp_ptr_Y,temp_ptr_X)*RadToDeg));
				_dist_mag_image[i*_img_size.width + j] = (uint16)(temp_ptr+0.5);
				
				
			}
		}
	}
	
}


void Orientation(CSize imagesize, uint16* Gmag, int16* Gdir, uint8 Template_size, uint8* plhs)
{
	//Procedure
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int Half_template_size = Template_size / 2;
	
	long int main_col, main_row;
	main_col = imagesize.width;
	main_row = imagesize.height;
	
	//malloc gu_weight_pre_computed matrix
	int size = 2 * Half_template_size - 1;
	double **gu_weight_pre_computed = (double **) malloc(size * sizeof(double *));
	for (int row = -Half_template_size + 1; row <= Half_template_size - 1; row++) {
		gu_weight_pre_computed[Half_template_size - 1 + row] = (double *) malloc(size * sizeof(double));
	}
	// calculate gu_weight_pre_computed values
	for (int row = -Half_template_size + 1; row <= Half_template_size - 1; row++) {
		for (int col = -Half_template_size + 1; col <= Half_template_size - 1; col++) {
			double angle_sigma = 1.5;
			gu_weight_pre_computed[Half_template_size - 1 + row][Half_template_size - 1 + col] = exp(-(double) (row * row + col * col) / (2 * angle_sigma * angle_sigma));
		}
	}
#pragma omp parallel for
	for (long int mask_row = 0; mask_row < main_row; mask_row++) {
		for (long int mask_col = 0; mask_col < main_col; mask_col++) {
			const double sub_ratio = 18.0 / 4.0;
			const double bin_angle = 360.0 / 18.0;
			double bin[18] = {0.0};
			double sub_bin[4] = {0.0};
			
			//calulation of major direction, refering to SIFT descriptor based on gradient
			for (int row = -Half_template_size + 1; row <= Half_template_size - 1; row++) {
				for (int col = -Half_template_size + 1; col <= Half_template_size - 1; col++) {
					double gu_weight = gu_weight_pre_computed[Half_template_size - 1 + row][Half_template_size - 1 + col];
					long int radius2 = (row * row + col * col);
					long int pixel_row = mask_row + row;
					long int pixel_col = mask_col + col;
					if (radius2 <= (Half_template_size - 1) * (Half_template_size - 1) && pixel_row > 0 && pixel_row < main_row - 1 && pixel_col > 0 && pixel_col < main_col - 1) {
						double mag, theta;
						mag = Gmag[pixel_row * main_col + pixel_col];
						theta = (double) (Gdir[pixel_row * main_col + pixel_col]);
						
						if (theta < 0) {
							theta += 360;
						}
						long int index = (long int) (theta / bin_angle);
						if (index >= 18) {
							printf("ERROR: Orientation array out of bounds access\n");
						}
						bin[index] += mag * gu_weight;
						sub_bin[(int) (index / sub_ratio)] += mag * gu_weight;
					}
				}
			}
			
			double max_th = -100;
			int sub_th = 0, th = 0;
			for (int i = 0; i < 4; i++) {
				if (sub_bin[i] > max_th) {
					max_th = sub_bin[i];
					sub_th = i;
				}
			}
			
			max_th = -100;
			for (int i = (int) (sub_th * sub_ratio); i < (int) (sub_th * sub_ratio + sub_ratio); i++) {
				if (bin[i] > max_th) {
					max_th = bin[i];
					th = i;
				}
			}
			plhs[mask_row * main_col + mask_col] = (uint8) th;
		}
	}
	// free all of gu_weight_pre_computed
	for (int j = 0; j < size; j++) {
		if (gu_weight_pre_computed[j]) {
			free(gu_weight_pre_computed[j]);
		}
	}
	
	if (gu_weight_pre_computed) {
		free(gu_weight_pre_computed);
		
	}
}

bool VerticalLineLocus(NCCresult* nccresult, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs, CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size,
					   CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT* Grid_wgs, UGRID *GridPT3, NCCflag flag,
					   uint8 NumofIAparam, double* ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration, uint8* left_ori, uint8* right_ori,
					   double bin_angle, uint8 NumOfCompute, uint8 peak_level, FILE* fid, bool IsPar, bool Hemisphere, char* save_filepath, uint8 tile_row, uint8 tile_col, double* Boundary,
					   bool pre_DEMtif, char* tmpdir,double *meters_per_pixel, bool IsRA,double mag_avg,double mag_var)
{
    if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(im_resolution*pow(2,Pyramid_step)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
        
        printf("VerticalLineLocus : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }
    
	int Half_template_size = (int)(Template_size/2);
    
	double subBoundary[4];

	double h_divide;
	double height_step;
	
	D2DPOINT temp_p1, temp_p2;
	D3DPOINT temp_GrP;
	double temp_LIA[2] = {0,0};
	
	int numofpts;
	
	int pixel_buffer = 1000;
	int GNCC_level	= 3;
	bool check_ortho = false;
	int Mag_th = 0;
	double ncc_alpha, ncc_beta;
	int N_th;
	
	N_th = 110;
	if(Pyramid_step == 3)
	{
		Half_template_size = Half_template_size - 1;
		N_th = 75;

	}
	else if(Pyramid_step < 3)
	{
		Half_template_size = Half_template_size - 2;
		N_th = 45;
	}
	
	subBoundary[0]	  = Boundary[0];
	subBoundary[1]	  = Boundary[1];
	subBoundary[2]	  = Boundary[2];
	subBoundary[3]	  = Boundary[3];
	
	numofpts = Size_Grid2D.height*Size_Grid2D.width;

	temp_GrP.m_X = Grid_wgs[(int)(numofpts/2)].m_X;
	temp_GrP.m_Y = Grid_wgs[(int)(numofpts/2)].m_Y;
	temp_GrP.m_Z = minmaxHeight[0];
	temp_GrP.flag = 0;
	temp_p1		= GetObjectToImageRPC_single_mpp(LRPCs,NumofIAparam,temp_LIA,temp_GrP);
	
	temp_GrP.m_Z = minmaxHeight[1];
	temp_p2		= GetObjectToImageRPC_single_mpp(LRPCs,NumofIAparam,temp_LIA,temp_GrP);
	
	*meters_per_pixel = (minmaxHeight[1] - minmaxHeight[0]) / sqrt( pow(temp_p1.m_X - temp_p2.m_X,2.0) + pow(temp_p1.m_Y - temp_p2.m_Y,2.0));
	
	printf("mpp = %f\n",*meters_per_pixel);
	
    if(*meters_per_pixel > 5)
        *meters_per_pixel = 5;
	
    if(*meters_per_pixel < im_resolution)
        *meters_per_pixel = im_resolution;
    
	printf("h,x,y =%f\t %f\t %f\t\n%f\t %f\t %f\n",minmaxHeight[0], temp_p1.m_X, temp_p1.m_Y, minmaxHeight[1],temp_p2.m_X, temp_p2.m_Y);
	printf("mpp = %f\n",*meters_per_pixel);
	
	im_resolution = im_resolution*pow(2,Pyramid_step);
	
    
    if((Pyramid_step == 4 && iteration == 1) || IsRA == true)
        check_ortho = false;
    else
        check_ortho = true;
	   
    if(pre_DEMtif)
        check_ortho = true;
	
    //orthoimage pixel information save
    F2DPOINT *all_left_im_cd, *all_right_im_cd;
    int sub_imagesize_w, sub_imagesize_h;
    
    if(check_ortho)
    {
        sub_imagesize_w = (int)((subBoundary[2] - subBoundary[0])/im_resolution)+1;
        sub_imagesize_h = (int)((subBoundary[3] - subBoundary[1])/im_resolution)+1;
        
        printf("imsize %d %d %u %u\n",sub_imagesize_w,sub_imagesize_h,Size_Grid2D.width,Size_Grid2D.height);

        long int sub_imagesize_total = (long int)sub_imagesize_w * (long int)sub_imagesize_h;
        printf("sub_imagesize_total %ld\n",sub_imagesize_total);
        
        all_left_im_cd = (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
		if (all_left_im_cd == NULL)
		{
			printf("ERROR: Out of memory - all_left_im_cd is NULL\n");
			exit(1);
		}
        all_right_im_cd= (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
		if (all_right_im_cd == NULL)
		{
			printf("ERROR: Out of memory - all_right_im_cd is NULL\n");
			exit(1);
		}

#pragma omp parallel for schedule(guided)
		for(long int iter_count = 0 ; iter_count < sub_imagesize_total ; iter_count++)
		{
           	int pts_row = (int)(floor(iter_count/sub_imagesize_w));
           	int pts_col = iter_count % sub_imagesize_w;
           	int pt_index;
           	double t_X, t_Y;
           	int t_col, t_row;
           	long int pt_index_im;
           
           	t_X		= subBoundary[0] + pts_col*im_resolution;
           	t_Y		= subBoundary[1] + pts_row*im_resolution;
           
           	t_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
           	t_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
           
           	pt_index	= t_row*Size_Grid2D.width + t_col;
           	pt_index_im = pts_row*(long int)sub_imagesize_w + pts_col;
           
           	if(pt_index < Size_Grid2D.width * Size_Grid2D.height && t_col < Size_Grid2D.width && t_row < Size_Grid2D.height)
           	{
               	if(GridPT3[pt_index].Height != -1000)
               	{
                   	double temp_LIA[2];
                   	D3DPOINT temp_GP;
                   	D2DPOINT temp_GP_p;
                   	D2DPOINT Left_Imagecoord, Right_Imagecoord;
                   	D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
                   
                   	temp_GP_p.m_X = t_X;
                   	temp_GP_p.m_Y = t_Y;
                   
                   	temp_GP_p	  = ps2wgs_single(param,temp_GP_p);
                   	temp_GP.m_X	  = temp_GP_p.m_X;
                   	temp_GP.m_Y	  = temp_GP_p.m_Y;
                   	temp_GP.m_Z	  = GridPT3[pt_index].Height;
                   
                   	temp_LIA[0] = 0.0;
                   	temp_LIA[1] = 0.0;
                   	Left_Imagecoord		= GetObjectToImageRPC_single(LRPCs,NumofIAparam,temp_LIA,temp_GP);
                   	Right_Imagecoord	= GetObjectToImageRPC_single(RRPCs,NumofIAparam,ImageAdjust,temp_GP);
                   
                   	Left_Imagecoord_py	= OriginalToPyramid_single(Left_Imagecoord,Lstartpos,Pyramid_step);
                   	Right_Imagecoord_py = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,Pyramid_step);
                   
                   	all_left_im_cd[pt_index_im].m_X = Left_Imagecoord_py.m_X;
                   	all_left_im_cd[pt_index_im].m_Y = Left_Imagecoord_py.m_Y;
                   	all_right_im_cd[pt_index_im].m_X= Right_Imagecoord_py.m_X;
                   	all_right_im_cd[pt_index_im].m_Y= Right_Imagecoord_py.m_Y;
               	}
           	}
		}
    }
    
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
	
	if(pre_DEMtif)
	{
		if (Pyramid_step == 4) {
			ncc_alpha = 0.9 - ((4-Pyramid_step)*0.3 + (iteration-1)*0.05);
			if(ncc_alpha < 0.7)
				ncc_alpha = 0.7;
		}
		else if(Pyramid_step == 3)
		{
			ncc_alpha = 1.0 - ((4-Pyramid_step)*0.3 + (iteration-1)*0.05);
			if(ncc_alpha < 0.5)
				ncc_alpha = 0.5;
		}
		else if(Pyramid_step == 2)
		{
			ncc_alpha = 1.0 - (0.4 + (iteration-1)*0.05);
			if(ncc_alpha < 0.4)
				ncc_alpha = 0.4;
		}
		else
			ncc_alpha = 1.0 - ((4-Pyramid_step)*0.2 + (iteration-1)*0.05);
	}

	if(ncc_alpha < 0.1)
		ncc_alpha = 0.1;
	
	
	
	double kernal_w[3];
	kernal_w[1] = 0.3;
	kernal_w[0] = 0.3;
	kernal_w[2] = 0.4;
	if(Pyramid_step >= 2)
	{
		if(iteration < 5)
		{
			kernal_w[0] = 0.3 - (iteration-1)*0.02;
			kernal_w[2] = 0.4 + (iteration-1)*0.02;
		}
		else
		{
			kernal_w[0] = 0.2;
			kernal_w[2] = 0.5;
		}
	}
	
	if(IsRA == 1)
		ncc_alpha = 1.0;
	
	ncc_beta = 1.0 - ncc_alpha;
	
	printf("ncc_alpha ncc_beta %f %f\n",ncc_alpha,ncc_beta);
	
	if(Pyramid_step >= 3)
		h_divide = 3;
	else if(Pyramid_step == 2)
	{
		h_divide = 3;
	}
	else if(Pyramid_step == 1)
		h_divide = 3;
	else {
		h_divide = 6;
	}
	
	
	height_step = (double)(im_resolution/h_divide);
    
    if(height_step > 3)
        height_step = 3;
    
	printf("height step %f\n",height_step);

#pragma omp parallel for schedule(guided)
	for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
	{
		int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
		int pts_col = iter_count % Size_Grid2D.width;
		int pt_index = pts_row*Size_Grid2D.width + pts_col;
		
		{
			double pre_rho	= -1.0;
			double pre_INCC_roh = -1.0;
			double pre_GNCC_roh = -1.0;
			double max_1stroh = -1.0;
			double max_2ndroh = -1.0;
			double max_1stheight = - 1000;
			double max_2ndheight = -1000;
			
			
			double pre_height= 0.0;
			int direction	= 0;
			bool check_rho	= false;
			int count_height;
			int NumOfHeights;
			double start_H, end_H;
			int i,j;
			double ortho_th = 0.8 - (4 - Pyramid_step)*0.10;
			
			start_H		= GridPT3[pt_index].minHeight;
			end_H		= GridPT3[pt_index].maxHeight;

			
			if(check_image_boundary(LRPCs,RRPCs,NumofIAparam,ImageAdjust,Lstartpos,Rstartpos,Grid_wgs[pt_index],start_H,end_H,LImagesize,RImagesize,Half_template_size,Pyramid_step))
			{
				char t_temp_path[500];
				bool check_blunder_cell = false;
				
				if(GridPT3[pt_index].Matched_flag == 0)
					h_divide = 2;
				
				
				NumOfHeights = (int)((end_H -  start_H)/height_step);
				
				if ( Pyramid_step >= 2)
					check_blunder_cell = false;
				else if(GridPT3[pt_index].Matched_flag == 0)
					check_blunder_cell = true;
				else if( (GridPT3[pt_index].Matched_flag != 0)	&& (fabs(GridPT3[pt_index].maxHeight - GridPT3[pt_index].minHeight) < 1000) )
					check_blunder_cell = false;
				else
					check_blunder_cell = true;
				
				if(!IsRA)
				{
					if(Pyramid_step == 1)
					{
						if(iteration > 1)
						{
							if ( (GridPT3[pt_index].Matched_flag != 0)	&& (fabs(GridPT3[pt_index].maxHeight - GridPT3[pt_index].minHeight) < 1000) )
								check_blunder_cell = false;
							else
								check_blunder_cell = true;
						}
						else
							check_blunder_cell = false;
					}
				}
				
				if(Pyramid_step <= 1)
				{
					if(GridPT3[pt_index].maxHeight > minmaxHeight[1] || GridPT3[pt_index].minHeight < minmaxHeight[0])
					{
						check_blunder_cell = true;
					}
					if(GridPT3[pt_index].Matched_flag == 0)
						check_blunder_cell = true;
				}

				if((nccresult[pt_index].result2 != 0) && (nccresult[pt_index].result2 >= GridPT3[pt_index].minHeight && nccresult[pt_index].result2 <= GridPT3[pt_index].maxHeight) &&
				   (GridPT3[pt_index].minHeight >= GridPT3[pt_index].t_minHeight && GridPT3[pt_index].maxHeight <= GridPT3[pt_index].t_maxHeight) )
				{
					check_blunder_cell = true;
				}
				
				if(Pyramid_step == 0 && GridPT3[pt_index].Matched_flag == 2 && iteration > 100)
				{
					//no calculation
				}
				if(!check_blunder_cell)
				{
					for(count_height = 0 ; count_height < NumOfHeights ; count_height++)
					{
						double iter_height;
						bool check_false_h = false;
						
						iter_height		= start_H + count_height*height_step;

						if(count_height == 0)
						{
							nccresult[pt_index].result0 = -1.0;
							nccresult[pt_index].result2 = -1000;
							nccresult[pt_index].result1 = -1.0;
							nccresult[pt_index].result3 = -1000;
							nccresult[pt_index].result4 = 0;
						}
						
						double temp_LIA[2];
						

						D3DPOINT temp_GP[1];
						D2DPOINT Left_Imagecoord[1], Right_Imagecoord[1];
						D2DPOINT Left_Imagecoord_py[1], Right_Imagecoord_py[1];
		
						

						temp_GP[0].m_X = Grid_wgs[pt_index].m_X;
						temp_GP[0].m_Y = Grid_wgs[pt_index].m_Y;
						temp_GP[0].m_Z = iter_height;
						
						temp_LIA[0] = 0.0;
						temp_LIA[1] = 0.0;
						
						Left_Imagecoord[0]		= GetObjectToImageRPC_single(LRPCs,NumofIAparam,temp_LIA,temp_GP[0]);
						Right_Imagecoord[0]		= GetObjectToImageRPC_single(RRPCs,NumofIAparam,ImageAdjust,temp_GP[0]);
						
						Left_Imagecoord_py[0]	= OriginalToPyramid_single(Left_Imagecoord[0],Lstartpos,Pyramid_step);
						Right_Imagecoord_py[0]	= OriginalToPyramid_single(Right_Imagecoord[0],Rstartpos,Pyramid_step);

						bool check_weight = false;
						double weight = 1.0;
						
						int ori_diff;
						
						ori_diff = left_ori[(int)Left_Imagecoord_py[0].m_Y*LImagesize.width + (int)Left_Imagecoord_py[0].m_X] - 
							right_ori[(int)Right_Imagecoord_py[0].m_Y*RImagesize.width + (int)Right_Imagecoord_py[0].m_X];
						{
							double Weight, Left_CR, Left_CC, Right_CR, Right_CC, diff_theta;
							bool b_count;
							double left_bin[36] = {0.0};
							int left_max_th, left_th;
							bool check_orientation = false;	
							
							b_count = false;
							left_max_th = -100;
							
							Weight		= weight;
							Left_CR		= Left_Imagecoord_py[0].m_Y;
							Left_CC		= Left_Imagecoord_py[0].m_X;
							Right_CR	= Right_Imagecoord_py[0].m_Y;
							Right_CC	= Right_Imagecoord_py[0].m_X;
							
							diff_theta	= (double)(ori_diff);

							if(IsRA != 1)
							{
								if(Pyramid_step >= 3 && diff_theta*bin_angle < 90 && diff_theta*bin_angle > - 90)
									check_orientation = true;
								else if(Pyramid_step <= 2)
								{
									check_orientation = true;
									
								}
							}
							else
							{
								check_orientation = true;
							}
							
							if(check_orientation)
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
								
								double Sum_LR_ortho = 0;
								double Sum_L_ortho = 0;
								double Sum_R_ortho = 0;
								double Sum_L2_ortho = 0;
								double Sum_R2_ortho = 0;
								double Sum_LR_2_ortho = 0;
								double Sum_L_2_ortho = 0;
								double Sum_R_2_ortho = 0;
								double Sum_L2_2_ortho = 0;
								double Sum_R2_2_ortho = 0;
								double Sum_LR_3_ortho = 0;
								double Sum_L_3_ortho = 0;
								double Sum_R_3_ortho = 0;
								double Sum_L2_3_ortho = 0;
								double Sum_R2_3_ortho = 0;
								int Count_N_ortho[3] = {0};	   
								
								double Sum_LR_mag_ortho = 0;
								double Sum_L_mag_ortho = 0;
								double Sum_R_mag_ortho = 0;
								double Sum_L2_mag_ortho = 0;
								double Sum_R2_mag_ortho = 0;
								double Sum_LR_2_mag_ortho = 0;
								double Sum_L_2_mag_ortho = 0;
								double Sum_R_2_mag_ortho = 0;
								double Sum_L2_2_mag_ortho = 0;
								double Sum_R2_2_mag_ortho = 0;
								double Sum_LR_3_mag_ortho = 0;
								double Sum_L_3_mag_ortho = 0;
								double Sum_R_3_mag_ortho = 0;
								double Sum_L2_3_mag_ortho = 0;
								double Sum_R2_3_mag_ortho = 0;
								
								
								int row, col;
								int N;
								
								double val1, val2, de, de2, ncc_1, ncc_2, ncc_3;
								double ncc_1_mag, ncc_2_mag, ncc_3_mag;
								double ncc_1_ortho, ncc_2_ortho, ncc_3_ortho;
								double ncc_1_mag_ortho, ncc_2_mag_ortho, ncc_3_mag_ortho;
								
								double temp_rho;
								double temp_INCC_roh = 0;
								double temp_GNCC_roh = 0;
								bool flag_value;
								int grid_index;
								double diff_rho;
								int t_direction;

																	
								uint16 mag_center_l, mag_center_r;
								
								if(flag.rotate_flag)
									rot_theta = (double)(diff_theta*bin_angle*PI/180.0);
								
								double cos0 = cos(-rot_theta);
								double sin0 = sin(-rot_theta);
								
								for(row = -Half_template_size; row <= Half_template_size ; row++)
								{
									for(col = -Half_template_size; col <= Half_template_size ; col++)
									{
										int radius2  =  row*row + col*col;
										if(radius2 <= (Half_template_size-1)*(Half_template_size-1))
										{
											double pos_row_left		 = (Left_CR + row);
											double pos_col_left		 = (Left_CC + col);
											
											double temp_col		   = (cos0*col - sin0*row);
											double temp_row		   = (sin0*col + cos0*row);
											double pos_row_right	 = (Right_CR + temp_row);
											double pos_col_right	 = (Right_CC + temp_col);
											
											
											if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right	>= 0 && pos_col_right+1	< RImagesize.width &&
												pos_row_left >= 0 && pos_row_left+1	  < LImagesize.height && pos_col_left	>= 0 && pos_col_left+1	< LImagesize.width)
											{
												//interpolate left_patch
												double dx		   =  pos_col_left - (int)(pos_col_left);
												double dy		   =  pos_row_left - (int)(pos_row_left);
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
												dx			=  pos_col_right - (int)(pos_col_right);
												dy			=  pos_row_right - (int)(pos_row_right);
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
												
												Sum_LR			  = Sum_LR + LR;
												Sum_L			  = Sum_L  + left_patch;
												Sum_R			  = Sum_R  + right_patch;
												Sum_L2			  = Sum_L2 + L2;
												Sum_R2			  = Sum_R2 + R2;
												
												Sum_LR_mag			  = Sum_LR_mag + LR_mag;
												Sum_L_mag			  = Sum_L_mag  + left_mag_patch;
												Sum_R_mag			  = Sum_R_mag  + right_mag_patch;
												Sum_L2_mag			  = Sum_L2_mag + L2_mag;
												Sum_R2_mag			  = Sum_R2_mag + R2_mag;
												
												int size_1, size_2;
												size_1		  = (int)(Half_template_size/2);
												if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
												{
													if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
													{
														Sum_LR_2  = Sum_LR_2 + LR;
														Sum_L_2	  = Sum_L_2	 + left_patch;
														Sum_R_2	  = Sum_R_2	 + right_patch;
														Sum_L2_2  = Sum_L2_2 + L2;
														Sum_R2_2  = Sum_R2_2 + R2;
														Count_N[1]++;
														
														Sum_LR_2_mag  = Sum_LR_2_mag + LR_mag;
														Sum_L_2_mag	  = Sum_L_2_mag	 + left_mag_patch;
														Sum_R_2_mag	  = Sum_R_2_mag	 + right_mag_patch;
														Sum_L2_2_mag  = Sum_L2_2_mag + L2_mag;
														Sum_R2_2_mag  = Sum_R2_2_mag + R2_mag;
													}
												}
														
												size_2		  = size_1 + (int)((size_1/2.0) + 0.5);
												if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
												{
													if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
													{
														Sum_LR_3  = Sum_LR_3 + LR;
														Sum_L_3	  = Sum_L_3	 + left_patch;
														Sum_R_3	  = Sum_R_3	 + right_patch;
														Sum_L2_3  = Sum_L2_3 + L2;
														Sum_R2_3  = Sum_R2_3 + R2;
														Count_N[2]++;
														
														Sum_LR_3_mag  = Sum_LR_3_mag + LR_mag;
														Sum_L_3_mag	  = Sum_L_3_mag	 + left_mag_patch;
														Sum_R_3_mag	  = Sum_R_3_mag	 + right_mag_patch;
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
											
                                            //ortho_NCC
                                            if(check_ortho)
                                            {
                                                double pos_row_left_ortho = -100;
                                                double pos_col_left_ortho = -100;
                                                double pos_row_right_ortho= -100;
                                                double pos_col_right_ortho= -100;
                                                
                                                if(row != 0 && col != 0)
                                                {
                                                    long int pt_index_temp,pt_index_dem;
                                                    double t_X, t_Y;
                                                    int t_col, t_row, tt_col, tt_row;
                                                    
                                                    t_X		= GridPts[pt_index].m_X + col*im_resolution;
                                                    t_Y		= GridPts[pt_index].m_Y + row*im_resolution;
                                                    
                                                    t_col	= (int)((t_X - subBoundary[0])/im_resolution);
                                                    t_row	= (int)((t_Y - subBoundary[1])/im_resolution);
                                                    
                                                    tt_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
                                                    tt_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
                                                    
                                                    pt_index_temp = t_row*(long int)sub_imagesize_w + t_col;
                                                    pt_index_dem  = tt_row*Size_Grid2D.width + tt_col;
                                                    
                                                    if(pt_index_temp >= 0 && pt_index_temp < sub_imagesize_w * sub_imagesize_h &&
                                                       t_col >= 0 && t_col < sub_imagesize_w && t_row >=0 && t_row < sub_imagesize_h && 
                                                       pt_index_dem >= 0 && pt_index_dem < Size_Grid2D.width * Size_Grid2D.height &&
                                                       tt_col >= 0 && tt_col < Size_Grid2D.width && tt_row >=0 && tt_row < Size_Grid2D.height &&
                                                        all_left_im_cd != NULL && all_right_im_cd != NULL)
                                                    {
                                                        if(GridPT3[pt_index_dem].Height != -1000)
                                                        {
                                                            pos_row_left_ortho	= all_left_im_cd[pt_index_temp].m_Y;
                                                            pos_col_left_ortho	= all_left_im_cd[pt_index_temp].m_X;
                                                            pos_row_right_ortho	= all_right_im_cd[pt_index_temp].m_Y;
                                                            pos_col_right_ortho	= all_right_im_cd[pt_index_temp].m_X;
                                                            
                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    pos_row_left_ortho = pos_row_left;
                                                    pos_col_left_ortho = pos_col_left;
                                                    pos_row_right_ortho= pos_row_right;
                                                    pos_col_right_ortho= pos_col_right;
                                                }
                                                
                                                pos_row_left = pos_row_left_ortho;
                                                pos_col_left = pos_col_left_ortho;
                                                pos_row_right= pos_row_right_ortho;
                                                pos_col_right= pos_col_right_ortho;
                                                
                                                if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right	>= 0 && pos_col_right+1	< RImagesize.width &&
                                                    pos_row_left >= 0 && pos_row_left+1	  < LImagesize.height && pos_col_left	>= 0 && pos_col_left+1	< LImagesize.width)
                                                {
                                                    //interpolate left_patch
                                                    double dx		   =  pos_col_left - (int)(pos_col_left);
                                                    double dy		   =  pos_row_left - (int)(pos_row_left);
                                                    double dxdy = dx * dy;
                                                    double left_patch;
                                                    double right_patch;
                                                    double left_mag_patch;
                                                    double right_mag_patch;
                                                    
                                                    long int position = (long int) (pos_col_left) + (long int) (pos_row_left) * LImagesize.width;
                                                    
                                                    left_patch =
                                                        (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
                                                        (double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) +
                                                        (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
                                                    
                                                    left_mag_patch = (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) +
                                                        (double) (MagImages_L[position + 1]) * (dx - dxdy) +
                                                        (double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) +
                                                        (double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
                                                    
                                                    //interpolate right_patch
                                                    dx			=  pos_col_right - (int)(pos_col_right);
                                                    dy			=  pos_row_right - (int)(pos_row_right);
                                                    dxdy = dx * dy;
                                                    position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
                                                    right_patch =
                                                        (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
                                                        (double) (RightImage[position + RImagesize.width]) * (dy - dxdy) +
                                                        (double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
                                                    
                                                    right_mag_patch = (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) +
                                                        (double) (MagImages_R[position + 1]) * (dx - dxdy) +
                                                        (double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) +
                                                        (double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
                                                    
                                                    //end
                                                    Count_N_ortho[0]++;
                                                    
                                                    double LR = left_patch * right_patch;
                                                    double L2 = left_patch * left_patch;
                                                    double R2 = right_patch * right_patch;
                                                    double LR_mag = left_mag_patch * right_mag_patch;
                                                    double L2_mag = left_mag_patch * left_mag_patch;
                                                    double R2_mag = right_mag_patch * right_mag_patch;
                                                    
                                                    Sum_LR_ortho			= Sum_LR_ortho + LR;
                                                    Sum_L_ortho				= Sum_L_ortho  + left_patch;
                                                    Sum_R_ortho				= Sum_R_ortho  + right_patch;
                                                    Sum_L2_ortho			= Sum_L2_ortho + L2;
                                                    Sum_R2_ortho			= Sum_R2_ortho + R2;
                                                    
                                                    Sum_LR_mag_ortho			= Sum_LR_mag_ortho + LR_mag;
                                                    Sum_L_mag_ortho				= Sum_L_mag_ortho  + left_mag_patch;
                                                    Sum_R_mag_ortho				= Sum_R_mag_ortho  + right_mag_patch;
                                                    Sum_L2_mag_ortho			= Sum_L2_mag_ortho + L2_mag;
                                                    Sum_R2_mag_ortho			= Sum_R2_mag_ortho + R2_mag;
                                                    
                                                    int size_1, size_2;
                                                    size_1		  = (int)(Half_template_size/2);
                                                    if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
                                                    {
                                                        if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
                                                        {
                                                            Sum_LR_2_ortho	= Sum_LR_2_ortho + LR;
                                                            Sum_L_2_ortho	= Sum_L_2_ortho	 + left_patch;
                                                            Sum_R_2_ortho	= Sum_R_2_ortho	 + right_patch;
                                                            Sum_L2_2_ortho	= Sum_L2_2_ortho + L2;
                                                            Sum_R2_2_ortho	= Sum_R2_2_ortho + R2;
                                                            Count_N_ortho[1]++;
                                                            
                                                            Sum_LR_2_mag_ortho	= Sum_LR_2_mag_ortho + LR_mag;
                                                            Sum_L_2_mag_ortho	= Sum_L_2_mag_ortho	 + left_mag_patch;
                                                            Sum_R_2_mag_ortho	= Sum_R_2_mag_ortho	 + right_mag_patch;
                                                            Sum_L2_2_mag_ortho	= Sum_L2_2_mag_ortho + L2_mag;
                                                            Sum_R2_2_mag_ortho	= Sum_R2_2_mag_ortho + R2_mag;
							}              
						    }
                                                    
                                                    
						    size_2		  = size_1 + (int)((size_1/2.0) + 0.5);
						    if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
						    {
						      if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
						      {
							Sum_LR_3_ortho	= Sum_LR_3_ortho + LR;
							Sum_L_3_ortho	= Sum_L_3_ortho	 + left_patch;
							Sum_R_3_ortho	= Sum_R_3_ortho	 + right_patch;
							Sum_L2_3_ortho	= Sum_L2_3_ortho + L2;
							Sum_R2_3_ortho	= Sum_R2_3_ortho + R2;
							Count_N_ortho[2]++;
                                                            
							Sum_LR_3_mag_ortho	= Sum_LR_3_mag_ortho + LR_mag;
							Sum_L_3_mag_ortho	= Sum_L_3_mag_ortho	 + left_mag_patch;
							Sum_R_3_mag_ortho	= Sum_R_3_mag_ortho	 + right_mag_patch;
							Sum_L2_3_mag_ortho	= Sum_L2_3_mag_ortho + L2_mag;
							Sum_R2_3_mag_ortho	= Sum_R2_3_mag_ortho + R2_mag;
                                                            
						      }
						    }
							
						    
                                                }
                                            }
										}
									}
								}

								N				= Count_N[0];
								val1		  = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
								val2		  = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								
								if( val1*val2 > 0)
								{
								  de			  = sqrt(val1*val2);
								  de2			  = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
								  ncc_1			= de2/de;
								}
								else
									ncc_1			= -1.0;
								
								val1		  = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/N;
								val2		  = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								if( val1*val2 > 0)
								{
								  de			  = sqrt(val1*val2);
								  de2			  = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/N;
								  ncc_1_mag			= de2/de;
								}
								else
									ncc_1_mag			= -1.0;
								
								N					= Count_N[1];
								val1				= (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
								val2				= (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								if( val1*val2 > 0)
								{
								  de					= sqrt(val1*val2);
								  de2					= (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
								  ncc_2		  = de2/de;
								}
								else
									ncc_2			= -1.0;
								
								val1				= (double)(Sum_L2_2_mag) - (double)(Sum_L_2_mag*Sum_L_2_mag)/N;
								val2				= (double)(Sum_R2_2_mag) - (double)(Sum_R_2_mag*Sum_R_2_mag)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								if( val1*val2 > 0)
								{
								  de					= sqrt(val1*val2);
								  de2					= (double)(Sum_LR_2_mag) - (double)(Sum_L_2_mag*Sum_R_2_mag)/N;
								  ncc_2_mag		  = de2/de;
								}
								else
									ncc_2_mag			= -1.0;
								
								
								N					= Count_N[2];
								val1				= (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
								val2				= (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								if( val1*val2 > 0)
								{
								  de					= sqrt(val1*val2);
								  de2					= (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
								  ncc_3		  = de2/de;
								}
								else
									ncc_3			= -1.0;
								
								val1				= (double)(Sum_L2_3_mag) - (double)(Sum_L_3_mag*Sum_L_3_mag)/N;
								val2				= (double)(Sum_R2_3_mag) - (double)(Sum_R_3_mag*Sum_R_3_mag)/N;
								if(Pyramid_step <= 1)
								{
									if(val1 == 0)
										val1 = 0.00001;
									if(val2 == 0)
										val2 = 0.00001;
								}
								
								if( val1*val2 > 0)
								{
								  de					= sqrt(val1*val2);
								  de2					= (double)(Sum_LR_3_mag) - (double)(Sum_L_3_mag*Sum_R_3_mag)/N;
								  ncc_3_mag		  = de2/de;
								}
								else
									ncc_3_mag			= -1.0;
								
								
								//ortho_NCC
					            N				= Count_N_ortho[0];
                                val1		  = (double)(Sum_L2_ortho) - (double)(Sum_L_ortho*Sum_L_ortho)/N;
                                val2		  = (double)(Sum_R2_ortho) - (double)(Sum_R_ortho*Sum_R_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de			  = sqrt(val1*val2);
				  de2			  = (double)(Sum_LR_ortho) - (double)(Sum_L_ortho*Sum_R_ortho)/N;
				  ncc_1_ortho			  = de2/de;
				}
                                else
                                    ncc_1_ortho			  = -1.0;
                                
                                
                                val1		  = (double)(Sum_L2_mag_ortho) - (double)(Sum_L_mag_ortho*Sum_L_mag_ortho)/N;
                                val2		  = (double)(Sum_R2_mag_ortho) - (double)(Sum_R_mag_ortho*Sum_R_mag_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de			  = sqrt(val1*val2);
				  de2			  = (double)(Sum_LR_mag_ortho) - (double)(Sum_L_mag_ortho*Sum_R_mag_ortho)/N;
				  ncc_1_mag_ortho			  = de2/de;
				}
                                else
                                    ncc_1_mag_ortho			  = -1.0;
                                
                                N					= Count_N_ortho[1];
                                val1				= (double)(Sum_L2_2_ortho) - (double)(Sum_L_2_ortho*Sum_L_2_ortho)/N;
                                val2				= (double)(Sum_R2_2_ortho) - (double)(Sum_R_2_ortho*Sum_R_2_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de					= sqrt(val1*val2);
				  de2					= (double)(Sum_LR_2_ortho) - (double)(Sum_L_2_ortho*Sum_R_2_ortho)/N;
				  ncc_2_ortho			= de2/de;
				}
                                else
                                    ncc_2_ortho			  = -1.0;
                                
                                val1				= (double)(Sum_L2_2_mag_ortho) - (double)(Sum_L_2_mag_ortho*Sum_L_2_mag_ortho)/N;
                                val2				= (double)(Sum_R2_2_mag_ortho) - (double)(Sum_R_2_mag_ortho*Sum_R_2_mag_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de					= sqrt(val1*val2);
				  de2					= (double)(Sum_LR_2_mag_ortho) - (double)(Sum_L_2_mag_ortho*Sum_R_2_mag_ortho)/N;
				  ncc_2_mag_ortho			= de2/de;
				}
                                else
                                    ncc_2_mag_ortho			  = -1.0;
                                
                                
                                N					= Count_N_ortho[2];
                                val1				= (double)(Sum_L2_3_ortho) - (double)(Sum_L_3_ortho*Sum_L_3_ortho)/N;
                                val2				= (double)(Sum_R2_3_ortho) - (double)(Sum_R_3_ortho*Sum_R_3_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de					= sqrt(val1*val2);
				  de2					= (double)(Sum_LR_3_ortho) - (double)(Sum_L_3_ortho*Sum_R_3_ortho)/N;
				  ncc_3_ortho			= de2/de;
				}
                                else
                                    ncc_3_ortho			  = -1.0;
                                
                                val1				= (double)(Sum_L2_3_mag_ortho) - (double)(Sum_L_3_mag_ortho*Sum_L_3_mag_ortho)/N;
                                val2				= (double)(Sum_R2_3_mag_ortho) - (double)(Sum_R_3_mag_ortho*Sum_R_3_mag_ortho)/N;
                                if(Pyramid_step <= 1)
                                {
                                    if(val1 == 0)
                                        val1 = 0.00001;
                                    if(val2 == 0)
                                        val2 = 0.00001;
                                }
                                
                                if( val1*val2 > 0)
				{
				  de					= sqrt(val1*val2);
				  de2					= (double)(Sum_LR_3_mag_ortho) - (double)(Sum_L_3_mag_ortho*Sum_R_3_mag_ortho)/N;
				  ncc_3_mag_ortho			= de2/de;
				}
                                else
                                    ncc_3_mag_ortho			  = -1.0;
                                
                                flag_value		= true;
								
								temp_INCC_roh = (double)(ncc_1 + ncc_2 + ncc_3 + ncc_1_mag + ncc_2_mag + ncc_3_mag)/6.0;
								temp_GNCC_roh = (double)(ncc_1_ortho + ncc_2_ortho + ncc_3_ortho + ncc_1_mag_ortho + ncc_2_mag_ortho + ncc_3_mag_ortho)/6.0;
							
							
								if(check_ortho && GridPT3[pt_index].ortho_ncc > ortho_th)
								{
									temp_rho = temp_INCC_roh*ncc_alpha + temp_GNCC_roh*ncc_beta;
								}
								else
									temp_rho = temp_INCC_roh;
						
								grid_index			 = pt_index;

								
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
									check_rho	= true;
									t_direction = -1;
								}
								else if(pre_rho != -1)
								{
									bool check_peak = false;
									
									if (Pyramid_step == 0 && iteration == 3)
									{
										if(direction > 0 && t_direction < 0 )
											check_peak = true;
									}
									else if(Pyramid_step == 4 && iteration < 3)
									{
										if (direction > 0 && t_direction < 0)
											check_peak = true;
									}
									else
									{
										
										if (direction > 0 && t_direction < 0)
											check_peak = true;
									}
									
									if( check_peak )
									{
										nccresult[grid_index].result4 += 1;
										if(nccresult[grid_index].result0 < pre_rho)
										{
											double temp_1, temp_2;
											temp_1 = nccresult[grid_index].result0;
											nccresult[grid_index].result0 = pre_rho;

											temp_2 = nccresult[grid_index].result2;
											nccresult[grid_index].result2 = pre_height;
											
											nccresult[grid_index].INCC = pre_INCC_roh;
											nccresult[grid_index].GNCC = pre_GNCC_roh;
							   
											if(nccresult[grid_index].result1 < temp_1)
											{
												nccresult[grid_index].result1 = temp_1;
												nccresult[grid_index].result3 = temp_2;
											}
										}
										else
										{
											if(nccresult[grid_index].result1 < pre_rho)
											{
												nccresult[grid_index].result1 = pre_rho;
												nccresult[grid_index].result3 = pre_height;
											}
										}
									}
								}

								pre_rho				   = temp_rho;
								
								pre_height			   = iter_height;
								direction			   = t_direction;
								
								pre_INCC_roh		   = temp_INCC_roh;
								pre_GNCC_roh		   = temp_GNCC_roh;
								nccresult[grid_index].roh_count ++;
							
							}
						}
					}
				}
			}
		}
	}
    
    if(check_ortho)
    {
        if (all_left_im_cd) free(all_left_im_cd);
        if (all_right_im_cd) free(all_right_im_cd);
    }
    
	return true;
}


double VerticalLineLocus_seeddem(uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs,
								CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
								CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT *Grid_wgs, UGRID *GridPT3,
								uint8 NumofIAparam, double* ImageAdjust, uint8 Pyramid_step, D2DPOINT Lstartpos, D2DPOINT Rstartpos, 
								char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,double* Boundary, double* minmaxHeight, double seedDEM_sigma)
{
	printf("minmax %f %f\n",minmaxHeight[0],minmaxHeight[1]);
	
	//double* nccresult1;
	//double* orthoheight;
	
	int Half_template_size = (int)(Template_size/2.0);
	
	if(Pyramid_step == 4)
	{
		Half_template_size = Template_size - iteration*2;
		if(Half_template_size < (int)(Template_size/2.0))
			Half_template_size = (int)(Template_size/2.0);
	}
	else if(Pyramid_step == 3)
	{
		Half_template_size = Template_size - 2	- (iteration)*2;
		if(Half_template_size < (int)(Template_size/2.0))
			Half_template_size = (int)(Template_size/2.0);
	}
	else
	{
		Half_template_size = (int)(Template_size/2.0);
	}
	
	double subBoundary[4];
	
	double h_divide;
	double height_step;
	
	D2DPOINT temp_p1, temp_p2;
	D3DPOINT temp_GrP;
	double temp_LIA[2] = {0,0};
	
	int numofpts;
	F2DPOINT *all_left_im_cd, *all_right_im_cd;
	int sub_imagesize_w, sub_imagesize_h;
	int pixel_buffer = 1000;
	int GNCC_level	= 3;
	bool check_ortho = false;
	int Mag_th = 0;
	double ncc_alpha, ncc_beta;
	
	subBoundary[0]	  = Boundary[0];
	subBoundary[1]	  = Boundary[1];
	subBoundary[2]	  = Boundary[2];
	subBoundary[3]	  = Boundary[3];
	
	numofpts = Size_Grid2D.height*Size_Grid2D.width;
	
	//nccresult1 = (double*)calloc(sizeof(double),numofpts);
	//orthoheight	 = (double*)calloc(sizeof(double),numofpts);
	
	im_resolution = im_resolution*pow(2,Pyramid_step);
	
	sub_imagesize_w = (int)((subBoundary[2] - subBoundary[0])/im_resolution)+1;
	sub_imagesize_h = (int)((subBoundary[3] - subBoundary[1])/im_resolution)+1;
	long int sub_imagesize_total = (long int)sub_imagesize_w * (long int)sub_imagesize_h;
	printf("sub_imagesize_total %ld\n",sub_imagesize_total);
	
	all_left_im_cd = (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
    if (all_left_im_cd == NULL)
	{
		printf("ERROR: Out of memory - all_left_im_cd is NULL\n");
		exit(1);
	}
	all_right_im_cd= (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
    if (all_right_im_cd == NULL)
	{
		printf("ERROR: Out of memory - all_right_im_cd is NULL\n");
		exit(1);
	}

#pragma omp parallel for schedule(guided)
	for(long int iter_count = 0 ; iter_count < sub_imagesize_total ; iter_count++)
	{
		int pts_row = (int)(floor(iter_count/sub_imagesize_w));
		int pts_col = iter_count % sub_imagesize_w;
		int pt_index;
		double t_X, t_Y;
		int t_col, t_row;
		long int pt_index_im;
		
		t_X		= subBoundary[0] + pts_col*im_resolution;
		t_Y		= subBoundary[1] + pts_row*im_resolution;
		
		t_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
		t_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
		
		pt_index	= t_row*Size_Grid2D.width + t_col;
		pt_index_im = pts_row*(long int)sub_imagesize_w + pts_col;
		
		if(pt_index < Size_Grid2D.width * Size_Grid2D.height && pts_row < sub_imagesize_h && pts_col < sub_imagesize_w && pts_row >= 0 && pts_col >= 0 &&
		   t_col >= 0 && t_row >= 0 && t_col < Size_Grid2D.width && t_row < Size_Grid2D.height)
		{
			if(GridPT3[pt_index].Height != -1000)
			{
				double temp_LIA[2];
				D3DPOINT temp_GP;
				D2DPOINT temp_GP_p;
				D2DPOINT Left_Imagecoord, Right_Imagecoord;
				D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
				
				temp_GP_p.m_X = t_X;
				temp_GP_p.m_Y = t_Y;
				
				temp_GP_p	  = ps2wgs_single(param,temp_GP_p);
				temp_GP.m_X	  = temp_GP_p.m_X;
				temp_GP.m_Y	  = temp_GP_p.m_Y;
				temp_GP.m_Z	  = GridPT3[pt_index].Height;
				
				//orthoheight[pt_index] = temp_GP.m_Z;
				
				temp_LIA[0] = 0.0;
				temp_LIA[1] = 0.0;
				Left_Imagecoord		= GetObjectToImageRPC_single(LRPCs,NumofIAparam,temp_LIA,temp_GP);
				Right_Imagecoord	= GetObjectToImageRPC_single(RRPCs,NumofIAparam,ImageAdjust,temp_GP);
				
				Left_Imagecoord_py	= OriginalToPyramid_single(Left_Imagecoord,Lstartpos,Pyramid_step);
				Right_Imagecoord_py = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,Pyramid_step);
				
				all_left_im_cd[pt_index_im].m_X = Left_Imagecoord_py.m_X;
                all_left_im_cd[pt_index_im].m_Y = Left_Imagecoord_py.m_Y;
				all_right_im_cd[pt_index_im].m_X= Right_Imagecoord_py.m_X;
                all_right_im_cd[pt_index_im].m_Y= Right_Imagecoord_py.m_Y;
			}
		}
	}
	
	int count_total = 0;
	int count_low = 0;
#pragma omp parallel for reduction(+:count_low, count_total) schedule(guided)
	for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
	{
		int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
		int pts_col = iter_count % Size_Grid2D.width;
		int pt_index;
		pt_index = pts_row*Size_Grid2D.width + pts_col;
		
		if(pt_index < Size_Grid2D.width * Size_Grid2D.height && pts_row < Size_Grid2D.height && pts_col < Size_Grid2D.width && pts_row >= 0 && pts_col >= 0)
		{
			
			int i,j;
			int row,col;
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
			int N;
			
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
			
			double ncc_1, ncc_1_mag, ncc_2, ncc_2_mag, ncc_3, ncc_3_mag;
			double val1, val2, de, de2;
			
			double nccresult = -10.0;
			//nccresult1[pt_index] = -10.0;
			
			if(GridPT3[pt_index].Height != -1000)
				count_total+=1;
			
			for(row = -Half_template_size; row <= Half_template_size ; row++)
			{
				for(col = -Half_template_size; col <= Half_template_size ; col++)
				{
					double pos_row_left = -100;
					double pos_col_left = -100;
					double pos_row_right= -100;
					double pos_col_right= -100;
					
					long int pt_index_temp,pt_index_dem;
					double t_X, t_Y;
					int t_col, t_row, tt_col, tt_row;
					
					t_X		= GridPts[pt_index].m_X + col*im_resolution;
					t_Y		= GridPts[pt_index].m_Y + row*im_resolution;
					
					t_col	= (int)((t_X - subBoundary[0])/im_resolution);
					t_row	= (int)((t_Y - subBoundary[1])/im_resolution);
					
					tt_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
					tt_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
					
					pt_index_temp = t_row*(long int)sub_imagesize_w + t_col;
					pt_index_dem  = tt_row*Size_Grid2D.width + tt_col;
					
					if(pt_index_temp >= 0 && pt_index_temp < sub_imagesize_w * sub_imagesize_h &&
					   t_col >= 0 && t_col < sub_imagesize_w && t_row >=0 && t_row < sub_imagesize_h && 
					   pt_index_dem >= 0 && pt_index_dem < Size_Grid2D.width * Size_Grid2D.height &&
					   tt_col >= 0 && tt_col < Size_Grid2D.width && tt_row >=0 && tt_row < Size_Grid2D.height &&
						all_left_im_cd != NULL && all_right_im_cd != NULL)
					{
						if(GridPT3[pt_index_dem].Height != -1000)
						{
							pos_row_left	= all_left_im_cd[pt_index_temp].m_Y;
							pos_col_left	= all_left_im_cd[pt_index_temp].m_X;
							pos_row_right	= all_right_im_cd[pt_index_temp].m_Y;
							pos_col_right	= all_right_im_cd[pt_index_temp].m_X;
							
							if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right	>= 0 && pos_col_right+1	< RImagesize.width &&
								pos_row_left >= 0 && pos_row_left+1	  < LImagesize.height && pos_col_left	>= 0 && pos_col_left+1	< LImagesize.width)
							{
								//interpolate left_patch
								double dx = pos_col_left - (int) (pos_col_left);
								double dy = pos_row_left - (int) (pos_row_left);
								double dxdy = dx * dy;
								long int position = (long int) pos_col_left + ((long int) pos_row_left) * LImagesize.width;
								
								double left_patch;
								double right_patch;
								double left_mag_patch;
								double right_mag_patch;
								
								left_patch = (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
									(double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) + (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
								
								left_mag_patch = (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_L[position + 1]) * (dx - dxdy) +
									(double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) +
									(double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
								
								
								//interpolate right_patch
								dx = pos_col_right - (int) (pos_col_right);
								dy = pos_row_right - (int) (pos_row_right);
								dxdy = dx * dy;
								position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
								right_patch = (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
									(double) (RightImage[position + RImagesize.width]) * (dy - dxdy) + (double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
								
								right_mag_patch = (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_R[position + 1]) * (dx - dxdy) +
									(double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) +
									(double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
								
								//end
								Count_N[0]++;
								
								double LR = left_patch * right_patch;
								double L2 = left_patch * left_patch;
								double R2 = right_patch * right_patch;
								double LR_mag = left_mag_patch * right_mag_patch;
								double L2_mag = left_mag_patch * left_mag_patch;
								double R2_mag = right_mag_patch * right_mag_patch;
								
								Sum_LR			  = Sum_LR + LR;
								Sum_L			  = Sum_L  + left_patch;
								Sum_R			  = Sum_R  + right_patch;
								Sum_L2			  = Sum_L2 + L2;
								Sum_R2			  = Sum_R2 + R2;
								
								Sum_LR_mag			  = Sum_LR_mag + LR_mag;
								Sum_L_mag			  = Sum_L_mag  + left_mag_patch;
								Sum_R_mag			  = Sum_R_mag  + right_mag_patch;
								Sum_L2_mag			  = Sum_L2_mag + L2_mag;
								Sum_R2_mag			  = Sum_R2_mag + R2_mag;
								
								int size_1, size_2;
								size_1		  = (int)(Half_template_size/2);
								if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
								{
									if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
									{
										Sum_LR_2  = Sum_LR_2 + LR;
										Sum_L_2	  = Sum_L_2	 + left_patch;
										Sum_R_2	  = Sum_R_2	 + right_patch;
										Sum_L2_2  = Sum_L2_2 + L2;
										Sum_R2_2  = Sum_R2_2 + R2;
										Count_N[1]++;
										
										Sum_LR_2_mag  = Sum_LR_2_mag + LR_mag;
										Sum_L_2_mag	  = Sum_L_2_mag	 + left_mag_patch;
										Sum_R_2_mag	  = Sum_R_2_mag	 + right_mag_patch;
										Sum_L2_2_mag  = Sum_L2_2_mag + L2_mag;
										Sum_R2_2_mag  = Sum_R2_2_mag + R2_mag;
										
									}
								}
								
								size_2		  = size_1 + (int)((size_1/2.0) + 0.5);
								if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
								{
									if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
									{
										Sum_LR_3  = Sum_LR_3 + LR;
										Sum_L_3	  = Sum_L_3	 + left_patch;
										Sum_R_3	  = Sum_R_3	 + right_patch;
										Sum_L2_3  = Sum_L2_3 + L2;
										Sum_R2_3  = Sum_R2_3 + R2;
										Count_N[2]++;
										
										Sum_LR_3_mag  = Sum_LR_3_mag + LR_mag;
										Sum_L_3_mag	  = Sum_L_3_mag	 + left_mag_patch;
										Sum_R_3_mag	  = Sum_R_3_mag	 + right_mag_patch;
										Sum_L2_3_mag  = Sum_L2_3_mag + L2_mag;
										Sum_R2_3_mag  = Sum_R2_3_mag + R2_mag;
									}
								}
							}
						}
					}
				}
			}
			
			if(Count_N[2] > 0)
			{
				N				= Count_N[0];
				val1		  = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
				val2		  = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de			  = sqrt(val1*val2);
				de2			  = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
				if( val1*val2 > 0)
					ncc_1			= de2/de;
				else
					ncc_1			= -1.0;
				
				val1		  = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/N;
				val2		  = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de			  = sqrt(val1*val2);
				de2			  = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/N;
				if( val1*val2 > 0)
					ncc_1_mag			= de2/de;
				else
					ncc_1_mag			= -1.0;
				
				
				N					= Count_N[1];
				val1				= (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
				val2				= (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
				if( val1*val2 > 0)
					ncc_2		  = de2/de;
				else
					ncc_2			= -1.0;
				
				val1				= (double)(Sum_L2_2_mag) - (double)(Sum_L_2_mag*Sum_L_2_mag)/N;
				val2				= (double)(Sum_R2_2_mag) - (double)(Sum_R_2_mag*Sum_R_2_mag)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_2_mag) - (double)(Sum_L_2_mag*Sum_R_2_mag)/N;
				if( val1*val2 > 0)
					ncc_2_mag		  = de2/de;
				else
					ncc_2_mag			= -1.0;
				
				
				N					= Count_N[2];
				val1				= (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
				val2				= (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
				if( val1*val2 > 0)
					ncc_3		  = de2/de;
				else
					ncc_3			= -1.0;
				
				val1				= (double)(Sum_L2_3_mag) - (double)(Sum_L_3_mag*Sum_L_3_mag)/N;
				val2				= (double)(Sum_R2_3_mag) - (double)(Sum_R_3_mag*Sum_R_3_mag)/N;
				if(Pyramid_step <= 4)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}
				
				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_3_mag) - (double)(Sum_L_3_mag*Sum_R_3_mag)/N;
				if( val1*val2 > 0)
					ncc_3_mag		  = de2/de;
				else
					ncc_3_mag			= -1.0;
				
				
				nccresult = (ncc_1 + ncc_1_mag + ncc_2 + ncc_2_mag + ncc_3_mag + ncc_3_mag)/6.0;
			}
			else 
			{
				nccresult = -1;
			}
			//nccresult1[pt_index] = nccresult;
			
			if(nccresult < 0.3)
			{
				GridPT3[pt_index].ortho_ncc	= nccresult;
	
				GridPT3[pt_index].minHeight		-= 100;
				if(GridPT3[pt_index].minHeight < minmaxHeight[0])
					GridPT3[pt_index].minHeight		= minmaxHeight[0];
				GridPT3[pt_index].maxHeight		+= 100;
				if(GridPT3[pt_index].maxHeight > minmaxHeight[1])
					GridPT3[pt_index].maxHeight		= minmaxHeight[1];
			
				if(nccresult < 0.1)
				{
					GridPT3[pt_index].minHeight	 = -9999;
					GridPT3[pt_index].maxHeight	 = -9999;
				}
				if(nccresult > -1)
					count_low += 1;
			}
			else
			{
				GridPT3[pt_index].ortho_ncc	= nccresult;
			}
		}
		
	}
	if (all_left_im_cd) free(all_left_im_cd);
	if (all_right_im_cd) free(all_right_im_cd);
			   
	//free(orthoheight);
	//free(nccresult1);
	printf("%d %d\n",count_low,count_total);
	
	return (double)(count_low)/(double)(count_total)*100;
}

bool VerticalLineLocus_blunder(double* nccresult, double* INCC, uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs,
							   CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
							   CSize Size_Grid2D, TransParam param, D2DPOINT* GridPts, D2DPOINT *Grid_wgs, UGRID *GridPT3,
							   uint8 NumofIAparam, double* ImageAdjust, uint8 Pyramid_step, D2DPOINT Lstartpos, D2DPOINT Rstartpos, 
							   char* save_filepath, uint8 tile_row, uint8 tile_col, uint8 iteration,uint8 bl_count,double* Boundary,uint8* left_ori, uint8* right_ori,int blunder_selected_level, bool bblunder)
{
  //uint16* orthoimage_l, *orthoimage_r;
	
  //double* orthoheight;
	
    //if(Pyramid_step >= 1)
    {
        double template_area = 5.0;
        int t_Template_size = (int)((template_area/(im_resolution*pow(2,blunder_selected_level)))/2.0)*2+1;
        if(Template_size < t_Template_size)
            Template_size = t_Template_size;
        
        printf("VerticalLineLocus_blunder : t Template_size %d\t%d\n",t_Template_size,Template_size);
    }

	int Half_template_size = (int)(Template_size/2.0);
	
	if(bblunder)
	{
		if(Pyramid_step == 4)
		{
			Half_template_size = Template_size - iteration*2;
			if(Half_template_size < (int)(Template_size/2.0))
				Half_template_size = (int)(Template_size/2.0);
		}
		else if(Pyramid_step == 3)
		{
			Half_template_size = Template_size - 2	- (iteration)*2;
			if(Half_template_size < (int)(Template_size/2.0))
				Half_template_size = (int)(Template_size/2.0);
		}
		else
		{
			Half_template_size = (int)(Template_size/2.0);
		}
	}
	
	double subBoundary[4];
	
	double h_divide;
	double height_step;
	
	D2DPOINT temp_p1, temp_p2;
	D3DPOINT temp_GrP;
	double temp_LIA[2] = {0,0};
	
	int numofpts;
	F2DPOINT *all_left_im_cd, *all_right_im_cd;
	int sub_imagesize_w, sub_imagesize_h;
	int pixel_buffer = 1000;
	int GNCC_level	= 3;
	bool check_ortho = false;
	int Mag_th = 0;
	
	subBoundary[0]	  = Boundary[0];
	subBoundary[1]	  = Boundary[1];
	subBoundary[2]	  = Boundary[2];
	subBoundary[3]	  = Boundary[3];
	
	numofpts = Size_Grid2D.height*Size_Grid2D.width;
	
	//orthoimage_l = (uint16*)calloc(sizeof(uint16),numofpts);
	//orthoimage_r = (uint16*)calloc(sizeof(uint16),numofpts);
	
	
	//orthoheight	 = (double*)calloc(sizeof(double),numofpts);
	
	im_resolution = im_resolution*pow(2,blunder_selected_level);
	
	printf("im_resolution = %f\n",im_resolution);
	
	sub_imagesize_w = (int)((subBoundary[2] - subBoundary[0])/im_resolution)+1;
	sub_imagesize_h = (int)((subBoundary[3] - subBoundary[1])/im_resolution)+1;
	long int sub_imagesize_total = (long int)sub_imagesize_w * (long int)sub_imagesize_h;
	printf("sub_imagesize_total %ld\n",sub_imagesize_total);
	
	all_left_im_cd = (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
    if (all_left_im_cd == NULL)
	{
		printf("ERROR: Out of memory - all_left_im_cd is NULL\n");
		exit(1);
	}
	all_right_im_cd= (F2DPOINT*)calloc(sizeof(F2DPOINT),sub_imagesize_total);
    if (all_right_im_cd == NULL)
	{
		printf("ERROR: Out of memory - all_right_im_cd is NULL\n");
		exit(1);
	}
	
#pragma omp parallel for schedule(guided)
	for(long int iter_count = 0 ; iter_count < sub_imagesize_total ; iter_count++)
	{
		int pts_row = (int)(floor(iter_count/sub_imagesize_w));
		int pts_col = iter_count % sub_imagesize_w;
		int pt_index;
		double t_X, t_Y;
		int t_col, t_row;
		long int pt_index_im;
		
		t_X		= subBoundary[0] + pts_col*im_resolution;
		t_Y		= subBoundary[1] + pts_row*im_resolution;
		
		t_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
		t_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
		
		pt_index	= t_row*Size_Grid2D.width + t_col;
		pt_index_im = pts_row*(long int)sub_imagesize_w + pts_col;
		
		if(pt_index < Size_Grid2D.width * Size_Grid2D.height && pts_row < sub_imagesize_h && pts_col < sub_imagesize_w && pts_row >= 0 && pts_col >= 0 &&
		   t_col >= 0 && t_row >= 0 && t_col < Size_Grid2D.width && t_row < Size_Grid2D.height)
		{
			if(GridPT3[pt_index].Height != -1000)
			{
				double temp_LIA[2];
				D3DPOINT temp_GP;
				D2DPOINT temp_GP_p;
				D2DPOINT Left_Imagecoord, Right_Imagecoord;
				D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
				
				temp_GP_p.m_X = t_X;
				temp_GP_p.m_Y = t_Y;
				
				temp_GP_p	  = ps2wgs_single(param,temp_GP_p);
				temp_GP.m_X	  = temp_GP_p.m_X;
				temp_GP.m_Y	  = temp_GP_p.m_Y;
				temp_GP.m_Z	  = GridPT3[pt_index].Height;
				
				//orthoheight[pt_index] = temp_GP.m_Z;
				
				temp_LIA[0] = 0.0;
				temp_LIA[1] = 0.0;
				Left_Imagecoord		= GetObjectToImageRPC_single(LRPCs,NumofIAparam,temp_LIA,temp_GP);
				Right_Imagecoord	= GetObjectToImageRPC_single(RRPCs,NumofIAparam,ImageAdjust,temp_GP);
				
				Left_Imagecoord_py	= OriginalToPyramid_single(Left_Imagecoord,Lstartpos,blunder_selected_level);
				Right_Imagecoord_py = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,blunder_selected_level);
				
				all_left_im_cd[pt_index_im].m_X = Left_Imagecoord_py.m_X;
                all_left_im_cd[pt_index_im].m_Y = Left_Imagecoord_py.m_Y;
				all_right_im_cd[pt_index_im].m_X= Right_Imagecoord_py.m_X;
                all_right_im_cd[pt_index_im].m_Y= Right_Imagecoord_py.m_Y;
				
			}
			else {
			  //orthoheight[pt_index]	= -1000;
			}
		}
	}
	
	double ncc_alpha, ncc_beta;
	
	if(Pyramid_step >= 3)
		ncc_alpha = 1.0 - ((4-Pyramid_step)*0.3 + (iteration-1)*0.10);
	else if(Pyramid_step == 2)
		ncc_alpha = 1.0 - (0.5 + (iteration-1)*0.10);
	else
		ncc_alpha = 1.0 - ((4-Pyramid_step)*0.2 + (iteration-1)*0.10);
	
	if(ncc_alpha < 0.1)
		ncc_alpha = 0.1;
	
#pragma omp parallel for schedule(guided)
	for(int iter_count = 0 ; iter_count < Size_Grid2D.height*Size_Grid2D.width ; iter_count++)
	{
		int pts_row = (int)(floor(iter_count/Size_Grid2D.width));
		int pts_col = iter_count % Size_Grid2D.width;
		int pt_index;
		pt_index = pts_row*Size_Grid2D.width + pts_col;
		
		if(pt_index < Size_Grid2D.width * Size_Grid2D.height && pts_row < Size_Grid2D.height && pts_col < Size_Grid2D.width && pts_row >= 0 && pts_col >= 0)
		{
			int i,j;
			int row,col;
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
			int N;
			
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

			double ncc_1, ncc_1_mag, ncc_2, ncc_2_mag, ncc_3, ncc_3_mag;
			double ncc_1_I, ncc_1_mag_I, ncc_2_I, ncc_2_mag_I, ncc_3_I, ncc_3_mag_I;
			double val1, val2, de, de2;
			
			nccresult[pt_index] = -1.0;
			INCC[pt_index] = -1.0;

			bool check_INCC = false;
			
			for(row = -Half_template_size; row <= Half_template_size ; row++)
			{
				for(col = -Half_template_size; col <= Half_template_size ; col++)
				{
					double radius  = sqrt((double)(row*row + col*col));
					if(radius <= Half_template_size-1)
					{
						double pos_row_left = -100;
						double pos_col_left = -100;
						double pos_row_right= -100;
						double pos_col_right= -100;
						
						int pt_index_temp,pt_index_dem,pt_index_temp_c;
						double t_X, t_Y, t_X_c,t_Y_c;
						int t_col, t_row, tt_col, tt_row,t_col_c, t_row_c;
						
						t_X		= GridPts[pt_index].m_X + col*im_resolution;
						t_Y		= GridPts[pt_index].m_Y + row*im_resolution;
						
						t_col	= (int)((t_X - subBoundary[0])/im_resolution);
						t_row	= (int)((t_Y - subBoundary[1])/im_resolution);
						
						tt_col	= (int)((t_X - subBoundary[0])/DEM_resolution);
						tt_row	= (int)((t_Y - subBoundary[1])/DEM_resolution);
						
						pt_index_temp = t_row*sub_imagesize_w + t_col;
						
						pt_index_dem  = tt_row*Size_Grid2D.width + tt_col;
						
						if(pt_index_temp >= 0 && pt_index_temp < sub_imagesize_w * sub_imagesize_h &&
						   t_col >= 0 && t_col < sub_imagesize_w && t_row >=0 && t_row < sub_imagesize_h && 
						   pt_index_dem >= 0 && pt_index_dem < Size_Grid2D.width * Size_Grid2D.height &&
						   tt_col >= 0 && tt_col < Size_Grid2D.width && tt_row >=0 && tt_row < Size_Grid2D.height &&
							all_left_im_cd != NULL && all_right_im_cd != NULL)
						{
							if(GridPT3[pt_index_dem].Height != -1000)
							{
								pos_row_left	= all_left_im_cd[pt_index_temp].m_Y;
								pos_col_left	= all_left_im_cd[pt_index_temp].m_X;
								pos_row_right	= all_right_im_cd[pt_index_temp].m_Y;
								pos_col_right	= all_right_im_cd[pt_index_temp].m_X;
								
								if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right	>= 0 && pos_col_right+1	< RImagesize.width &&
									pos_row_left >= 0 && pos_row_left+1	  < LImagesize.height && pos_col_left	>= 0 && pos_col_left+1	< LImagesize.width)
								{
									//interpolate left_patch
									double dx		   =  pos_col_left - (int)(pos_col_left);
									double dy		   =  pos_row_left - (int)(pos_row_left);
									double left_patch;
									double right_patch;
									double left_mag_patch;
									double right_mag_patch;
									double dxdy = dx * dy;
									long int position = (long int) (pos_col_left) + (long int) (pos_row_left) * LImagesize.width;
									
									left_patch = (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
										(double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) + (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
									
									left_mag_patch = (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_L[position + 1]) * (dx - dxdy) +
										(double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) +
										(double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
									
									
									//interpolate right_patch
									dx			=  pos_col_right - (int)(pos_col_right);
									dy			=  pos_row_right - (int)(pos_row_right);
									dxdy = dx * dy;
									position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
									
									right_patch = (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
										(double) (RightImage[position + RImagesize.width]) * (dy - dxdy) +
										(double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
									
									right_mag_patch = (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_R[position + 1]) * (dx - dxdy) +
										(double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) +
										(double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
		
									//end
									Count_N[0]++;
									
									double LR = left_patch * right_patch;
									double L2 = left_patch * left_patch;
									double R2 = right_patch * right_patch;
									double LR_mag = left_mag_patch * right_mag_patch;
									double L2_mag = left_mag_patch * left_mag_patch;
									double R2_mag = right_mag_patch * right_mag_patch;
									
									Sum_LR			  = Sum_LR + LR;
									Sum_L			  = Sum_L  + left_patch;
									Sum_R			  = Sum_R  + right_patch;
									Sum_L2			  = Sum_L2 + L2;
									Sum_R2			  = Sum_R2 + R2;
									
									Sum_LR_mag			  = Sum_LR_mag + LR_mag;
									Sum_L_mag			  = Sum_L_mag  + left_mag_patch;
									Sum_R_mag			  = Sum_R_mag  + right_mag_patch;
									Sum_L2_mag			  = Sum_L2_mag + L2_mag;
									Sum_R2_mag			  = Sum_R2_mag + R2_mag;
									
									int size_1, size_2;
									size_1		  = (int)(Half_template_size/2);
									if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
									{
										if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
										{
											Sum_LR_2  = Sum_LR_2 + LR;
											Sum_L_2	  = Sum_L_2	 + left_patch;
											Sum_R_2	  = Sum_R_2	 + right_patch;
											Sum_L2_2  = Sum_L2_2 + L2;
											Sum_R2_2  = Sum_R2_2 + R2;
											Count_N[1]++;
											
											Sum_LR_2_mag  = Sum_LR_2_mag + LR_mag;
											Sum_L_2_mag	  = Sum_L_2_mag	 + left_mag_patch;
											Sum_R_2_mag	  = Sum_R_2_mag	 + right_mag_patch;
											Sum_L2_2_mag  = Sum_L2_2_mag + L2_mag;
											Sum_R2_2_mag  = Sum_R2_2_mag + R2_mag;
											
										}
									}
									
									size_2		  = size_1 + (int)((size_1/2.0) + 0.5);
									if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
									{
										if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
										{
											Sum_LR_3  = Sum_LR_3 + LR;
											Sum_L_3	  = Sum_L_3	 + left_patch;
											Sum_R_3	  = Sum_R_3	 + right_patch;
											Sum_L2_3  = Sum_L2_3 + L2;
											Sum_R2_3  = Sum_R2_3 + R2;
											Count_N[2]++;
											
											Sum_LR_3_mag  = Sum_LR_3_mag + LR_mag;
											Sum_L_3_mag	  = Sum_L_3_mag	 + left_mag_patch;
											Sum_R_3_mag	  = Sum_R_3_mag	 + right_mag_patch;
											Sum_L2_3_mag  = Sum_L2_3_mag + L2_mag;
											Sum_R2_3_mag  = Sum_R2_3_mag + R2_mag;
										}
									}
									
									if(row == 0 && col == 0)
									{
									  //orthoimage_l[pt_index] = (int)(left_patch);
									  //orthoimage_r[pt_index] = (int)(right_patch);

									}
								}
							}
						}
					}
				}
			}
			
			if(Count_N[0] > 0)
			{
				N				= Count_N[0];
				val1		  = fabs((double)(Sum_L2) - (double)(Sum_L*Sum_L)/N);
				val2		  = fabs((double)(Sum_R2) - (double)(Sum_R*Sum_R)/N);
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de			  = sqrt(val1*val2);
				de2			  = ((double)(Sum_LR) - (double)(Sum_L*Sum_R)/N);
				if( val1*val2 > 0)
					ncc_1			= de2/de;
				else
					ncc_1			= -1.0;
				
				val1		  = fabs((double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/N);
				val2		  = fabs((double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/N);
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de			  = sqrt(val1*val2);
				de2			  = ((double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/N);
				if( val1*val2 > 0)
					ncc_1_mag			= de2/de;
				else
					ncc_1_mag			= -1.0;
				
				
				N					= Count_N[1];
				val1				= (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
				val2				= (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
				if( val1*val2 > 0)
					ncc_2		  = de2/de;
				else
					ncc_2			= -1.0;
				
				val1				= (double)(Sum_L2_2_mag) - (double)(Sum_L_2_mag*Sum_L_2_mag)/N;
				val2				= (double)(Sum_R2_2_mag) - (double)(Sum_R_2_mag*Sum_R_2_mag)/N;
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_2_mag) - (double)(Sum_L_2_mag*Sum_R_2_mag)/N;
				if( val1*val2 > 0)
					ncc_2_mag		  = de2/de;
				else
					ncc_2_mag			= -1.0;
				
				
				N					= Count_N[2];
				val1				= (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
				val2				= (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
				if( val1*val2 > 0)
					ncc_3		  = de2/de;
				else
					ncc_3			= -1.0;
				
				val1				= (double)(Sum_L2_3_mag) - (double)(Sum_L_3_mag*Sum_L_3_mag)/N;
				val2				= (double)(Sum_R2_3_mag) - (double)(Sum_R_3_mag*Sum_R_3_mag)/N;
				if(Pyramid_step <= 1)
				{
					if(val1 == 0)
						val1 = 0.00001;
					if(val2 == 0)
						val2 = 0.00001;
				}

				de					= sqrt(val1*val2);
				de2					= (double)(Sum_LR_3_mag) - (double)(Sum_L_3_mag*Sum_R_3_mag)/N;
				if( val1*val2 > 0)
					ncc_3_mag		  = de2/de;
				else
					ncc_3_mag			= -1.0;
				
				nccresult[pt_index] = (ncc_1 + ncc_2 + ncc_3 + ncc_1_mag + ncc_2_mag + ncc_3_mag)/6.0;
			}
			else 
			{
				nccresult[pt_index] = -1;
			}
			
			GridPT3[pt_index].ortho_ncc = nccresult[pt_index];
		}
	}

	if (all_left_im_cd) free(all_left_im_cd);
	if (all_right_im_cd) free(all_right_im_cd);
	
	//free(orthoimage_l);
	//free(orthoimage_r);
	//free(orthoheight);
	
	return true;
}

bool VerticalLineLocus_Ortho(double *F_Height,D3DPOINT ref1_pt, D3DPOINT ref2_pt, D3DPOINT target_pt,
							 uint16 *MagImages_L,uint16 *MagImages_R, uint16* LeftImage, uint16* RightImage,
							 double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs, 
							 CSize LImagesize,	CSize RImagesize, CSize Size_Grid2D, TransParam param, uint8 NumofIAparam, 
							 double* ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, double meters_per_pixel,
							 D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration,  UGRID *GridPT3, int target_index, int ref1_index, int ref2_index,
							 double* boundary)
{
	double h_divide;
	double height_step;
	int NumOfHeights;
	double start_H, end_H;
	double F_NCC = -1.0;
	bool check_ncc = false;
	double gridspace = DEM_resolution;
	uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
	
	*F_Height = -9999;
	
	TIN_Grid_Size_X = Size_Grid2D.width;
	TIN_Grid_Size_Y = Size_Grid2D.height;
	
	if(Pyramid_step >= 3)
		h_divide = 1;
	else if(Pyramid_step == 2)
	{
		h_divide = 1;
	}
	else if(Pyramid_step ==1)
		h_divide = 1;
	else {
		h_divide = 1;
	}

	if(meters_per_pixel > 3)
		meters_per_pixel = 3;
	
	height_step = (double)(pow(2.0,(double)Pyramid_step)*meters_per_pixel/h_divide);
	
    //if(height_step > 10)
    //    height_step = 10;
    
	start_H		= GridPT3[target_index].minHeight;
	end_H		= GridPT3[target_index].maxHeight;
	
	NumOfHeights = (int)((end_H -  start_H)/height_step) + 1;
	
	int th_count = 20;

	for(int count_height = 0 ; count_height < NumOfHeights ; count_height++)
	{
		double TriP1[3];
		double TriP2[3];
		double TriP3[3];
		
		double TriMinXY[2], TriMaxXY[2];
		int PixelMinXY[2]={0};
		int PixelMaxXY[2]={0};
		int Col, Row;
		
		double val1, val2, de, de2, ncc_1,ncc_2,ncc;
		double Sum_LR = 0;
		double Sum_L = 0;
		double Sum_R = 0;
		double Sum_L2 = 0;
		double Sum_R2 = 0;
		
		double Sum_LR_mag = 0;
		double Sum_L_mag = 0;
		double Sum_R_mag = 0;
		double Sum_L2_mag = 0;
		double Sum_R2_mag = 0;
		
		int Count_N = 0;
		
		double iter_height		= start_H + count_height*height_step;
		
		TriP1[0]		= ref1_pt.m_X;
		TriP2[0]		= ref2_pt.m_X;
		TriP3[0]		= target_pt.m_X;
		
		TriP1[1]		= ref1_pt.m_Y;
		TriP2[1]		= ref2_pt.m_Y;
		TriP3[1]		= target_pt.m_Y;
		
		TriP1[2]		= ref1_pt.m_Z;
		TriP2[2]		= ref2_pt.m_Z;
		TriP3[2]		= iter_height;
		
		// calculation on BoundingBox(MinMax XY) of triangle
		TriMinXY[0]	= min(min(TriP1[0],TriP2[0]),TriP3[0]);
		TriMinXY[1]	= min(min(TriP1[1],TriP2[1]),TriP3[1]);
		TriMaxXY[0]	= max(max(TriP1[0],TriP2[0]),TriP3[0]);
		TriMaxXY[1]	= max(max(TriP1[1],TriP2[1]),TriP3[1]);
		
		PixelMinXY[0] = (int)((TriMinXY[0] - boundary[0])/gridspace + 0.5);
		PixelMinXY[1] = (int)((TriMinXY[1] - boundary[1])/gridspace + 0.5);
		PixelMaxXY[0] = (int)((TriMaxXY[0] - boundary[0])/gridspace + 0.5);
		PixelMaxXY[1] = (int)((TriMaxXY[1] - boundary[1])/gridspace + 0.5);
		
		PixelMinXY[0] -= 1;		PixelMinXY[1] -= 1;
		PixelMaxXY[0] += 1;		PixelMaxXY[1] += 1;
		if (PixelMaxXY[0] >= (int)(TIN_Grid_Size_X))	
			PixelMaxXY[0] =	 (int)(TIN_Grid_Size_X-1);
		if (PixelMaxXY[1] >= (int)(TIN_Grid_Size_Y))   
			PixelMaxXY[1] =	 (int)(TIN_Grid_Size_Y-1);
		if (PixelMinXY[0] < 0)	
			PixelMinXY[0] = 0;
		if (PixelMinXY[1] < 0)	
			PixelMinXY[1] = 0;
		
		for (Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
		  {
		    for (Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
		      {	
				double CurGPXY[2]={0.};
				double Z = 0.0;
				bool rtn = false;
				double _p1[3], _p2[3], _p3[3], v12[2], v1P[2];
				double v23[2], v2P[2], v31[2], v3P[2];
				int Sum;
				
				CurGPXY[0]	= (Col)*gridspace + boundary[0];
				CurGPXY[1]	= (Row)*gridspace + boundary[1];
				
				_p1[0]		= TriP1[0];
				_p1[1]		= TriP1[1];
				_p1[2]		= TriP1[2];
				
				_p2[0]		= TriP2[0];
				_p2[1]		= TriP2[1];
				_p2[2]		= TriP2[2];
				
				_p3[0]		= TriP3[0];
				_p3[1]		= TriP3[1];
				_p3[2]		= TriP3[2];
				
				v12[0]		= _p2[0]-_p1[0];
				v12[1]		= _p2[1]-_p1[1];
				
				v1P[0]		= CurGPXY[0]-_p1[0];
				v1P[1]		= CurGPXY[1]-_p1[1];
				
				v23[0]		= _p3[0]-_p2[0];
				v23[1]		= _p3[1]-_p2[1];
				
				v2P[0]		= CurGPXY[0]-_p2[0];
				v2P[1]		= CurGPXY[1]-_p2[1];
				
				v31[0]		= _p1[0]-_p3[0];
				v31[1]		= _p1[1]-_p3[1];
				
				v3P[0]		= CurGPXY[0]-_p3[0];
				v3P[1]		= CurGPXY[1]-_p3[1];
				
				Sum = 3;
				if (v12[0]*v1P[1]-v12[1]*v1P[0]<=0)
					Sum--;
				if (v23[0]*v2P[1]-v23[1]*v2P[0]<=0)
					Sum--;
				if (v31[0]*v3P[1]-v31[1]*v3P[0]<=0)
					Sum--;
				
				if (Sum==0 || Sum==3)
				{
					double v12[3], v13[3], Len, A, B, C, D;
					double Normal[3]={0.};
					
					v12[0]	= _p2[0]-_p1[0];
					v12[1]	= _p2[1]-_p1[1];
					v12[2]	= _p2[2]-_p1[2];
					
					v13[0]	= _p3[0]-_p1[0];
					v13[1]	= _p3[1]-_p1[1];
					v13[2]	= _p3[2]-_p1[2];
					
					Normal[0]=v12[1]*v13[2] - v12[2]*v13[1];
					Normal[1]=v12[2]*v13[0] - v12[0]*v13[2];
					Normal[2]=v12[0]*v13[1] - v12[1]*v13[0];	  
					
					Len=sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
					if(Len > 0)
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
					else
						rtn = false;
				}
				
				if (rtn)
				{
					D3DPOINT temp_GP;
					D2DPOINT Left_Imagecoord, Right_Imagecoord, temp_GP_p;
					D2DPOINT Left_Imagecoord_py, Right_Imagecoord_py;
					double temp_LIA[2];
					
					double pos_row_left;
					double pos_col_left;
					double pos_row_right;
					double pos_col_right;
					
					temp_GP_p.m_X = CurGPXY[0];
					temp_GP_p.m_Y = CurGPXY[1];
					
					temp_GP_p	  = ps2wgs_single(param,temp_GP_p);

					temp_GP.m_X = temp_GP_p.m_X;
					temp_GP.m_Y = temp_GP_p.m_Y;
					temp_GP.m_Z = Z;
					
					temp_LIA[0] = 0.0;
					temp_LIA[1] = 0.0;
					Left_Imagecoord		= GetObjectToImageRPC_single(LRPCs,NumofIAparam,temp_LIA,temp_GP);
					Right_Imagecoord	= GetObjectToImageRPC_single(RRPCs,NumofIAparam,ImageAdjust,temp_GP);
					
					Left_Imagecoord_py	= OriginalToPyramid_single(Left_Imagecoord,Lstartpos,Pyramid_step);
					Right_Imagecoord_py = OriginalToPyramid_single(Right_Imagecoord,Rstartpos,Pyramid_step);
					
					pos_row_left	  = Left_Imagecoord_py.m_Y;
					pos_col_left	  = Left_Imagecoord_py.m_X;
					
					pos_row_right	  = Right_Imagecoord_py.m_Y;
					pos_col_right	  = Right_Imagecoord_py.m_X;
					
					if( pos_row_right >= 0 && pos_row_right+1 < RImagesize.height && pos_col_right	>= 0 && pos_col_right+1	< RImagesize.width &&
						pos_row_left >= 0 && pos_row_left+1	  < LImagesize.height && pos_col_left	>= 0 && pos_col_left+1	< LImagesize.width)
					{
						//interpolate left_patch
						double dx = pos_col_left - (int) (pos_col_left);
						double dy = pos_row_left - (int) (pos_row_left);
						double left_patch;
						double right_patch;
						double left_mag_patch;
						double right_mag_patch;
						double dxdy = dx * dy;
						long int position = (long int) (pos_col_left) + (long int) (pos_row_left) * LImagesize.width;
						left_patch = (double) (LeftImage[position]) * (1 - dx - dy + dxdy) + (double) (LeftImage[position + 1]) * (dx - dxdy) +
							(double) (LeftImage[position + LImagesize.width]) * (dy - dxdy) + (double) (LeftImage[position + 1 + LImagesize.width]) * (dxdy);
						
						left_mag_patch = (double) (MagImages_L[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_L[position + 1]) * (dx - dxdy) +
							(double) (MagImages_L[position + LImagesize.width]) * (dy - dxdy) + (double) (MagImages_L[position + 1 + LImagesize.width]) * (dxdy);
						
						//interpolate right_patch
						dx = pos_col_right - (int) (pos_col_right);
						dy = pos_row_right - (int) (pos_row_right);
						dxdy = dx * dy;
						position = (long int) (pos_col_right) + (long int) (pos_row_right) * RImagesize.width;
						
						right_patch = (double) (RightImage[position]) * (1 - dx - dy + dxdy) + (double) (RightImage[position + 1]) * (dx - dxdy) +
							(double) (RightImage[position + RImagesize.width]) * (dy - dxdy) + (double) (RightImage[position + 1 + RImagesize.width]) * (dxdy);
						
						right_mag_patch = (double) (MagImages_R[position]) * (1 - dx - dy + dxdy) + (double) (MagImages_R[position + 1]) * (dx - dxdy) +
							(double) (MagImages_R[position + RImagesize.width]) * (dy - dxdy) + (double) (MagImages_R[position + 1 + RImagesize.width]) * (dxdy);
						
						//end
						Count_N++;
						
						Sum_LR			  = Sum_LR + left_patch*right_patch;
						Sum_L			  = Sum_L  + left_patch;
						Sum_R			  = Sum_R  + right_patch;
						Sum_L2			  = Sum_L2 + left_patch*left_patch;
						Sum_R2			  = Sum_R2 + right_patch*right_patch;
						
						Sum_LR_mag			  = Sum_LR_mag + left_mag_patch*right_mag_patch;
						Sum_L_mag			  = Sum_L_mag  + left_mag_patch;
						Sum_R_mag			  = Sum_R_mag  + right_mag_patch;
						Sum_L2_mag			  = Sum_L2_mag + left_mag_patch*left_mag_patch;
						Sum_R2_mag			  = Sum_R2_mag + right_mag_patch*right_mag_patch;
					}
				}
			}
		}
		
		if(Count_N >= th_count)
		{
			val1		  = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/Count_N;
			val2		  = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/Count_N;
			de			  = sqrt(val1*val2);
			de2			  = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/Count_N;
			if( val1*val2 > 0)
				ncc_1			= de2/de;
			else
				ncc_1			= -1.0;
			
			val1		  = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/Count_N;
			val2		  = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/Count_N;
			de			  = sqrt(val1*val2);
			de2			  = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/Count_N;
			if( val1*val2 > 0)
				ncc_2			= de2/de;
			else
				ncc_2			= -1.0;
			
			ncc = (ncc_1 + ncc_2)/2.0;
			if(ncc > 0.5)
			{
				check_ncc = true;
				if(F_NCC < ncc)
				{
					F_NCC = ncc;
					*F_Height = iter_height;
				}
			}
		}
		else if(Count_N > 0)
		{
			val1		  = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/Count_N;
			val2		  = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/Count_N;
			de			  = sqrt(val1*val2);
			de2			  = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/Count_N;
			if( val1*val2 > 0)
				ncc_1			= de2/de;
			else
				ncc_1			= -1.0;
			
			val1		  = (double)(Sum_L2_mag) - (double)(Sum_L_mag*Sum_L_mag)/Count_N;
			val2		  = (double)(Sum_R2_mag) - (double)(Sum_R_mag*Sum_R_mag)/Count_N;
			de			  = sqrt(val1*val2);
			de2			  = (double)(Sum_LR_mag) - (double)(Sum_L_mag*Sum_R_mag)/Count_N;
			if( val1*val2 > 0)
				ncc_2			= de2/de;
			else
				ncc_2			= -1.0;
			
			ncc = (ncc_1 + ncc_2)/2.0;
			if(ncc > 0.7)
			{ 
				check_ncc = true;
				*F_Height = (ref1_pt.m_Z + ref2_pt.m_Z)/2.0;
			}
		}
	}
	
	return true;
}


D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
	int i;
	D2DPOINT* out;

	out = (D2DPOINT*)malloc(sizeof(D2DPOINT)*numofpts);
	for(i=0;i<numofpts;i++)
	{
		out[i].m_X		= (InCoord[i].m_X/pow(2,Pyramid_step)) - Startpos.m_X;
		out[i].m_Y		= (InCoord[i].m_Y/pow(2,Pyramid_step)) - Startpos.m_Y;
	}

	return out;
	
}

D2DPOINT OriginalToPyramid_single(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
	D2DPOINT out;

	out.m_X		 = (InCoord.m_X/pow(2,Pyramid_step)) - Startpos.m_X;
	out.m_Y		 = (InCoord.m_Y/pow(2,Pyramid_step)) - Startpos.m_Y;

	return out;
	
}

D2DPOINT* PyramidToOriginal(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
	int i;
	D2DPOINT* out;
	out = (D2DPOINT*)malloc(sizeof(D2DPOINT)*numofpts);

	if(Pyramid_step > 0)
	{
		for(i=0;i<numofpts;i++)
		{
			out[i].m_X		= (InCoord[i].m_X +	 Startpos.m_X)*pow(2,Pyramid_step) + pow(2,Pyramid_step)/2.0;
			out[i].m_Y		= (InCoord[i].m_Y +	 Startpos.m_Y)*pow(2,Pyramid_step) + pow(2,Pyramid_step)/2.0;
		}
	}
	else
	{
		for(i=0;i<numofpts;i++)
		{
			out[i].m_X		= InCoord[i].m_X +	Startpos.m_X;
			out[i].m_Y		= InCoord[i].m_Y +	Startpos.m_Y;
		}
	}

	return out;
}


int SelectMPs(NCCresult* roh_height, CSize Size_Grid2D, D2DPOINT *GridPts_XY, UGRID *GridPT3,
			  double Th_roh, double Th_roh_min, double Th_roh_start, double Th_roh_next, uint8 Pyramid_step, uint8 total_pyramid,
			  uint8 iteration, uint8 peak_level, char *filename_mps, int pre_DEMtif, int IsRA, double MPP, double DEM_resolution, double im_resolution, int final_level_iteration)
{
	int count_MPs = 0;

	double PPM			= MPP; 
	printf("PPM of select MPs = %f\n",PPM);
	double roh_next;

	double minimum_Th = 0.1;
	
	if(Pyramid_step == 0)
		minimum_Th = 0.2;
	
	printf("minimum TH %f\n",minimum_Th);
	bool check_iter_end = false;
	FILE* temp_fid;
	double h_divide;
	
	if(Pyramid_step >= 3)
		h_divide = 2;
	else if(Pyramid_step == 2)
	{
		h_divide = 2;
	}
	else if(Pyramid_step ==1)
		h_divide = 2;
	else {
		h_divide = 6;
	}

	temp_fid			= fopen(filename_mps,"w");

	if(Pyramid_step >= 5)
		roh_next	= 0;
	else
		roh_next	= (double)(0.05);
	
	if(Pyramid_step == 0)
	{
		if(Th_roh - 0.10 < Th_roh_min)
			check_iter_end	= true;
	}
	else if(Pyramid_step <= 1)
	{
		if(Th_roh - 0.10 < Th_roh_min)
			check_iter_end	= true;
	}
	else if(Pyramid_step == 2)
	{
		if(Th_roh - 0.10 < Th_roh_min)
			check_iter_end	= true;
	}
	else if(Pyramid_step == 3) 
	{
		if(Th_roh - 0.10 < Th_roh_min)
			check_iter_end	= true;
	}
	else 
	{
		if(IsRA)
		{
			if(Th_roh - 0.10 < Th_roh_min)
				check_iter_end	= true;

		}
		else
		{
			if(Th_roh - 0.06 < Th_roh_min)
				check_iter_end	= true;
		}
	}

	im_resolution = im_resolution*pow(2,Pyramid_step);
	
	//#pragma omp parallel for shared(Size_Grid2D,GridPT3,GridPts_XY,Th_roh,roh_height,Pyramid_step,peak_level,PPM,temp_fid,check_iter_end,roh_next,Th_roh_start,Th_roh_next,total_pyramid,iteration) private(row,col) reduction(+:count_MPs)
	for(int iter_index = 0 ; iter_index < Size_Grid2D.height*Size_Grid2D.width ; iter_index++)
	{
		int row,col;
		row		= (int)(floor(iter_index/Size_Grid2D.width));
		col		= iter_index % Size_Grid2D.width;
		
		if(row >= 0 && row < Size_Grid2D.height && col >= 0 && col < Size_Grid2D.width)
		{
			bool index,index_1,index_2,index_3, roh_index;
			double temp, ROR, temp_h;
			double roh_th;
			D3DPOINT temp_mp;

			int grid_index = row*Size_Grid2D.width + col;

			if(Pyramid_step == 4 && iteration == 1)
				GridPT3[grid_index].Height			= -1000.0;
			
			index			= false;
			index_1			= false;
			index_2			= false;
			index_3			= false;
			roh_index		= false;
			

			roh_th		= roh_next;

			//ratio of 1st peak roh / 2nd peak roh
			ROR			= (roh_height[grid_index].result0 - roh_height[grid_index].result1)/roh_height[grid_index].result0;
			
			if(Pyramid_step <= 2)	
			{
				if(ROR >= (0.1 - 0.03*(3 - Pyramid_step)) && roh_height[grid_index].result0 > minimum_Th)
					index_2	= true;
			}
			else
			{
				if(ROR >= 0.1 && roh_height[grid_index].result0 > minimum_Th)
					index_2	= true;
			}
			
			// threshold of 1st peak roh
			if(roh_height[grid_index].result0 > GridPT3[grid_index].roh - roh_th)
				index_3		= true;

			if(roh_height[grid_index].result4 > 0)
				index		= true;

			//distance between min and max height
			temp			= abs(GridPT3[grid_index].maxHeight - GridPT3[grid_index].minHeight); 
			temp_h			= temp/(PPM/h_divide);
			
			if(Pyramid_step >= 1)
				roh_index	= index & index_2 & index_3;
			else if(Pyramid_step == 0 && final_level_iteration < 3)
			{
				/*if( index_3	 && roh_height[grid_index].result0 > minimum_Th)
					index_1	= true;
				
				roh_index	= ((index_3 & index_2) | index_1);*/
                roh_index	= index & index_2 & index_3;
			}
			else
			{
				if( index_3)
					index_1	= true;
				
				roh_index	= (index_2 | index_1);
			}
			
			if(GridPT3[grid_index].Matched_flag != 0)
			{
				index			= false;
				index_1			= false;
				index_2			= false;
				
				
				if(Pyramid_step <= 2)
				{
					if( (ROR < (0.1 - 0.03*(3 - Pyramid_step)) ) && (roh_height[grid_index].result1 > GridPT3[grid_index].roh - roh_th) ) 
						index	= true;
				}
				else
				{
					if( (ROR < 0.1) && (roh_height[grid_index].result1 > GridPT3[grid_index].roh - roh_th) ) 
						index	= true;
				}
				
				if(roh_height[grid_index].result2 > roh_height[grid_index].result3)
					index_1		= true;
				index_2			= !index_1;
				index_1			= index & index_1;
				index_2			= index & index_2;
				
				if(index_1)
				{
					GridPT3[grid_index].minHeight = (roh_height[grid_index].result3 - 0.5);
					GridPT3[grid_index].maxHeight = (roh_height[grid_index].result2 + 0.5);
					GridPT3[grid_index].Matched_flag = 4;
					
				}

				if(index_2)
				{
					GridPT3[grid_index].minHeight = (roh_height[grid_index].result2 - 0.5);
					GridPT3[grid_index].maxHeight = (roh_height[grid_index].result3 + 0.5);
					GridPT3[grid_index].Matched_flag = 4;
				}
			}
		
			if(Pyramid_step == 0 && GridPT3[grid_index].Matched_flag == 2 && iteration > 100)
			{
				temp_mp.m_Z	= GridPT3[grid_index].Height;
				temp_mp.m_X	= GridPts_XY[grid_index].m_X;
				temp_mp.m_Y	= GridPts_XY[grid_index].m_Y;
				
				count_MPs++;
				
				fprintf(temp_fid,"%lf %lf %lf 0\n",temp_mp.m_X,temp_mp.m_Y,temp_mp.m_Z);
	
				if(GridPT3[grid_index].roh < minimum_Th)
				{
					GridPT3[grid_index].roh = minimum_Th;
				}
				
			}
			else
			{
				//Set the matched pts and information
				if(roh_index)
				{
					double pre_H, pre_range;
					uint8 pre_Match;
					bool index_1, index_2, index_3, index_4, index_5, index_41, index_42, index_6, index_7, index;
					index			= false;
					index_1			= false;
					index_2			= false;
					index_3			= false;
					index_4			= false;
					index_41		= false;
					index_42		= false;
					index_5			= false;
					index_6			= false;
					index_7			= false;

					temp_mp.m_Z	= roh_height[grid_index].result2;
					temp_mp.m_X	= GridPts_XY[grid_index].m_X;
					temp_mp.m_Y	= GridPts_XY[grid_index].m_Y;
					
					//update MPs by previsous level results
					if(Pyramid_step < 2)
					{
						pre_H			= GridPT3[grid_index].Height;
						pre_range		= GridPT3[grid_index].maxHeight - GridPT3[grid_index].minHeight;
						pre_Match		= GridPT3[grid_index].Matched_flag;
						
						if(pre_Match == 0) //non matched grid
							index_1		= true;
						else if(pre_Match == 2) //matched grid
							index_2		= true;
						else if(pre_Match == 1 || pre_Match == 4) //inside triangle and min,max change status by matched grid
							index_3		= true;
						
						
						if(fabs(pre_H - temp_mp.m_Z) <= PPM*2*pow(2,Pyramid_step))
							index_4		= true;
						if(fabs(pre_H - temp_mp.m_Z) <= pre_range/2.0)
							index_5		= true;
						if(fabs(temp_mp.m_Z - GridPT3[grid_index].minHeight) <= PPM*3*pow(2,Pyramid_step))
							index_41	= true;
						if(fabs(temp_mp.m_Z - GridPT3[grid_index].maxHeight) <= PPM*3*pow(2,Pyramid_step))
							index_42	= true;
						
						index_6			= (index_4 | index_5 | index_41 | index_42) & index_3; // inside triangle by matched grid
						index_7			= index_2 & index_4; // matched grid
						
						index			= index_1 | index_7 | index_6;
						
						if(index)
						{
							count_MPs++;
							
							fprintf(temp_fid,"%lf %lf %lf 0\n",temp_mp.m_X,temp_mp.m_Y,temp_mp.m_Z);

							// update max_roh value
							GridPT3[grid_index].roh		= roh_height[grid_index].result0;
							if(GridPT3[grid_index].roh < minimum_Th)
							{
								GridPT3[grid_index].roh = minimum_Th;
							}
							
						}
						else
						{
							GridPT3[grid_index].roh = Th_roh;
							if(check_iter_end && Th_roh_start == GridPT3[grid_index].roh)
								GridPT3[grid_index].roh		= Th_roh_next;
						}
						
					}
					else
					{
						count_MPs++;
						
						fprintf(temp_fid,"%lf %lf %lf 0\n",temp_mp.m_X,temp_mp.m_Y,temp_mp.m_Z);

						// update max_roh value
						GridPT3[grid_index].roh		= roh_height[grid_index].result0;
						if(GridPT3[grid_index].roh < minimum_Th)
						{
							GridPT3[grid_index].roh = minimum_Th;
						}
					}
				}
				else
				{
					GridPT3[grid_index].roh = Th_roh;
					// update max_roh value
					if(check_iter_end && Th_roh_start == GridPT3[grid_index].roh)
						GridPT3[grid_index].roh		= Th_roh_next;
				}
			}
		}
	}
	fclose(temp_fid);

	return count_MPs;
}


UI3DPOINT *TINgeneration(bool last_flag, char *savepath, uint8 level, CSize Size_Grid2D, double img_resolution, double grid_resolution,
						 double min_max[],
						 double *subBoundary, int total_point_count, D3DPOINT *ptslists, int *iter_row, int *iter_col,
						 int *re_total_tri_counts)
{
	double interval_X, interval_Y;
	double interval_buffer;
	double X_distance, Y_distance;
	int interval_col, interval_row;
	
	
	int iter_count;
	int total_tri_counts = 0;
	int total_mps_counts = 0;
	int total_ptslist_count = 0;
	int count_tri;
	
	int i = 0;
	
	X_distance	= Size_Grid2D.width*grid_resolution;
	Y_distance	= Size_Grid2D.height*grid_resolution;
	
	interval_X	= X_distance;
	interval_Y	= Y_distance;
	
	interval_buffer = 0;
	
	interval_col = 1;
	interval_row = 1;
	
	if(level <= 3 )
	{
		if(grid_resolution <= pow(2,level)*img_resolution)
		{
			if (level == 3)
			{
				interval_X = 2000;
				interval_Y = 2000;
			}
			else if(level == 2)
			{
				interval_X = 1000;
				interval_Y = 1000;
			}
			else
			{
				interval_X = 500;
				interval_Y = 500;
			}
		}
		else
		{
			if (level == 3)
			{
			}
			else if(level == 2)
			{
				interval_X = 2000;
				interval_Y = 2000;
				
			}
			else
			{
				interval_X = 1000;
				interval_Y = 1000;
			}
		}
		
		interval_buffer = 4*(level+1)*grid_resolution;
		interval_col = (int)(X_distance/interval_X) + 1;
		interval_row = (int)(Y_distance/interval_Y) + 1;
	}
	
	if(grid_resolution > 10) //RA calculation
	{
		interval_X = 3000;
		interval_Y = 3000;
	}
	
	interval_buffer = 4*(level+1)*grid_resolution;
	interval_col = (int)(X_distance/interval_X) + 1;
	interval_row = (int)(Y_distance/interval_Y) + 1;
	
	printf("interval %f\t%f\t%d\t%d\n",interval_X,interval_Y,interval_col,interval_row);
	*iter_row = interval_row;
	*iter_col = interval_col;
	
	
	
	UI3DPOINT *trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*total_point_count*4);
	int tri_counts = 0;
	//#pragma omp parallel shared(interval_col,interval_row,savepath,subBoundary,interval_X,interval_Y,interval_buffer,grid_resolution,ptslists) private(iter_count) reduction(+:total_mps_counts,total_tri_counts)
	{
		//#pragma omp for ordered
		for(iter_count = 0; iter_count < interval_row*interval_col ; iter_count++)
		{
			int count_MPs_nums;
			int count_ptslist_num;
			
			int t_x, t_y;
			double pt_boundary[4];
			
			int count_MPs = 0;
			int count_re_tri = 0;
			CSize t_Size_Grid2D;
			
			
			char temp_pts_file[500];
			FILE *p_tempfile;
			
			int t_i;
			
			
			char bufstr[500];
			
			t_y		= (int)(floor(iter_count/interval_col));
			t_x		= iter_count % interval_col;
			sprintf(temp_pts_file,"%s/txt/temp_pts_%d_%d_.txt",savepath,t_x,t_y);
			
			pt_boundary[0] = subBoundary[0] + t_x*interval_X - interval_buffer;
			pt_boundary[1] = subBoundary[1] + t_y*interval_Y - interval_buffer;
			
			t_Size_Grid2D.width		= (int)((interval_X + 2*interval_buffer)/grid_resolution);
			t_Size_Grid2D.height	= (int)((interval_Y + 2*interval_buffer)/grid_resolution);
			
			pt_boundary[2] = pt_boundary[0] + t_Size_Grid2D.width*grid_resolution;
			pt_boundary[3] = pt_boundary[1] + t_Size_Grid2D.height*grid_resolution;
			
			D3DPOINT* temp_selected_ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*total_point_count);
			
			count_MPs_nums = 0;
			if(!last_flag)
			{
				for(t_i = 0 ; t_i < total_point_count ; t_i++)
				{
					if(ptslists[t_i].m_X >= pt_boundary[0] && ptslists[t_i].m_X <= pt_boundary[2] &&
					   ptslists[t_i].m_Y >= pt_boundary[1] && ptslists[t_i].m_Y <= pt_boundary[3])
					{
						if(ptslists[t_i].flag != 1 && ptslists[t_i].flag != 2)
						{
							temp_selected_ptslists[count_MPs_nums].m_X = ptslists[t_i].m_X;
							temp_selected_ptslists[count_MPs_nums].m_Y = ptslists[t_i].m_Y;
							temp_selected_ptslists[count_MPs_nums].m_Z = t_i;
							
							count_MPs_nums++;
						}
					}
				}
			}
			else {
				for(t_i = 0 ; t_i < total_point_count ; t_i++)
				{
					if(ptslists[t_i].m_X >= pt_boundary[0] && ptslists[t_i].m_X <= pt_boundary[2] &&
					   ptslists[t_i].m_Y >= pt_boundary[1] && ptslists[t_i].m_Y <= pt_boundary[3])
					{
						if(ptslists[t_i].flag != 1)
						{
							temp_selected_ptslists[count_MPs_nums].m_X = ptslists[t_i].m_X;
							temp_selected_ptslists[count_MPs_nums].m_Y = ptslists[t_i].m_Y;
							temp_selected_ptslists[count_MPs_nums].m_Z = t_i;
							
							count_MPs_nums++;
						}
					}
				}
			}
			
			//printf("count mps %d\n",count_MPs_nums);
 
			D3DPOINT *selected_ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs_nums);
			
			uint32 *temp_array = (uint32*)malloc(sizeof(uint32)*count_MPs_nums);
			
			for(t_i=0;t_i<count_MPs_nums;t_i++)
			{
				double temp_X, temp_Y, temp_Z,temp_f,temp_id;
				
				selected_ptslists[t_i].m_X = temp_selected_ptslists[t_i].m_X;
				selected_ptslists[t_i].m_Y = temp_selected_ptslists[t_i].m_Y;
				temp_array[t_i]			   = temp_selected_ptslists[t_i].m_Z;
			}
			
			UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_MPs_nums*4);
			
			sprintf(bufstr,"%s/txt/tri_%d_%d.txt",savepath,t_x,t_y);
			TINCreate(selected_ptslists,bufstr,count_MPs_nums,t_trilists,min_max,&count_tri);
			
			i = 0;
			for(i=0;i<count_tri;i++)
			{
				trilists[tri_counts].m_X = temp_array[t_trilists[i].m_X];
				trilists[tri_counts].m_Y = temp_array[t_trilists[i].m_Y];
				trilists[tri_counts].m_Z = temp_array[t_trilists[i].m_Z];
				tri_counts++;
			}
			
			free(t_trilists);
			
			total_tri_counts			+= count_tri;
			
			count_tri = 0;
			
			free(temp_array);
			
			free(selected_ptslists);
			
			free(temp_selected_ptslists);
		}
	}
	
	
	*re_total_tri_counts = total_tri_counts;
	
	
	UI3DPOINT *trilists_f	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*total_tri_counts);
	
	printf("total_tri_counts %d\n",total_tri_counts);
	
	for(i=0;i<total_tri_counts;i++)
	{
		trilists_f[i].m_X = trilists[i].m_X;
		trilists_f[i].m_Y = trilists[i].m_Y;
		trilists_f[i].m_Z = trilists[i].m_Z;
	}
	
	free(trilists);
	
	return trilists_f;
}

int DecisionMPs(bool flag_blunder, int count_MPs_input, double* Boundary, UGRID *GridPT3, uint8 Pyramid_step, double grid_resolution,
				uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_mps, char *filename_tri, double Hinterval, 
				bool *p_flag, double *pre_3sigma, double *pre_mean, int *count_Results, double *minz_mp, double *maxz_mp, double *minmaxHeight,
				uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs,
				CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
				TransParam param, D2DPOINT* Grid_wgs,D2DPOINT* GridPts,
				uint8 NumofIAparam, double* ImageAdjust, D2DPOINT Lstartpos, D2DPOINT Rstartpos, 
				char* save_filepath, uint8 tile_row, uint8 tile_col, int pre_DEMtif, int IsRA,uint8* left_ori, uint8* right_ori,int blunder_selected_level,double seedDEMsigma, double f_demsize)
{
	char bufstr[500];
	uint16 count			= 0;
	*p_flag				= true;
	int count_MPs		= count_MPs_input;
	int TIN_split_level = 0;
	int count_tri;
	
	if (grid_resolution <= 8)
	{
		if(grid_resolution <= 8 && grid_resolution > 4)
			TIN_split_level = 0;
		else if(grid_resolution <= 4)
			TIN_split_level = 2;
	}
	if(IsRA)
	{
		TIN_split_level = 4;
	}
	
	if(count_MPs > 3)
	{
		double gridspace			= grid_resolution;
		uint8 max_count			= 30;
		if(!flag_blunder)
		{
			if(Pyramid_step == 4)
				max_count = 10;
			else if(Pyramid_step == 3)
				max_count = 10;
			else {
				max_count = 10;
			}
		}
		
		bool check_dh			= false;
		bool height_check_flag	= false;
		bool flag				= true;
		

		FILE *survey;
		int Th_pts;

		printf("DEM resolution %f\n",f_demsize);
	
		if(flag_blunder)
		{
			if(Pyramid_step == 3)
				max_count  = 20;
			else if(Pyramid_step == 2)
			{
				if(f_demsize < 8 )
					max_count = 20;
				else
					max_count = 10;
			}
			else if(Pyramid_step == 1)
				max_count = 10;
			else if(Pyramid_step == 0)
			{
				if( iteration == 1)
					max_count = 10;
				else if(iteration == 2)
					max_count = 10;
				else if( iteration >= 3)
					max_count = 10;
			}
		}
		
		survey	= fopen(filename_mps_pre,"r");
		if(survey)
		{
			int i=0;
			BL blunder_param;
			double blunder_dh = 0;
			int blunder_th_counts;
			int pre_count_blunder = -100;
			int pre_count_MPs = -100;
			uint32 blunder_count[2];
			int iter_row,iter_col;
			
			

			blunder_param.Boundary	= Boundary;
			blunder_param.gridspace	= grid_resolution;
			blunder_param.height_check_flag = true;
			blunder_param.Hinterval = Hinterval;
			blunder_param.iteration = iteration;
			blunder_param.Pyramid_step = Pyramid_step;
			
			blunder_param.Size_Grid2D.width = Size_Grid2D.width;
			blunder_param.Size_Grid2D.height = Size_Grid2D.height;
			
			Th_pts		= 1;
			
			if(Pyramid_step >= 3)
				blunder_th_counts = 1;
			else {
				blunder_th_counts = 1;
			}
			
			//TIN generation start
			D3DPOINT *ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*count_MPs);
			
			
			double maxX_ptslists = -100000000;
			double maxY_ptslists = -100000000;
			double minX_ptslists =	100000000;
			double minY_ptslists =	100000000;
			
			while( i < count_MPs && (fscanf(survey,"%lf %lf %lf %hhd\n",&ptslists[i].m_X,&ptslists[i].m_Y,&ptslists[i].m_Z,&ptslists[i].flag)) != EOF )
			{
				if(maxX_ptslists < ptslists[i].m_X)
					maxX_ptslists = ptslists[i].m_X;
				if(maxY_ptslists < ptslists[i].m_Y)
					maxY_ptslists = ptslists[i].m_Y;
				if(minX_ptslists > ptslists[i].m_X)
					minX_ptslists = ptslists[i].m_X;
				if(minY_ptslists > ptslists[i].m_Y)
					minY_ptslists = ptslists[i].m_Y;
				
				i++;
			}
			double min_max[4] = {minX_ptslists, minY_ptslists, maxX_ptslists, maxY_ptslists};
			if( !(Pyramid_step == 0 && iteration == 3) )
			{
				UI3DPOINT *trilists;
				
				
				if(Pyramid_step >= TIN_split_level || count_MPs < 10000)
				{
					UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_MPs*4);
					
					sprintf(bufstr,"%s/txt/tri_%d_%d.txt",filename_tri,flag_blunder,count);
					TINCreate(ptslists,bufstr,count_MPs,t_trilists,min_max,&count_tri);
					
					trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
					i = 0;
					for(i=0;i<count_tri;i++)
					{
						trilists[i].m_X = t_trilists[i].m_X;
						trilists[i].m_Y = t_trilists[i].m_Y;
						trilists[i].m_Z = t_trilists[i].m_Z;
					}
					
					free(t_trilists);
				}
				else
				{
					trilists = TINgeneration(false,save_filepath, Pyramid_step, Size_Grid2D, im_resolution, DEM_resolution,
											 min_max,
											 Boundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
				}
				// TIN generation end
				
				
				while(flag == true && count < max_count && count_tri > 20)
				{
					int count_blunders = 0;
					double mt_minmaxheight[2];
					double* ortho_ncc = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
					double* INCC = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
					blunder_count[0] = 0;
					blunder_count[1] = 0;
					
					if(IsRA != 1)
					{
						SetHeightRange_blunder(minmaxHeight,ptslists, count_MPs, trilists,count_tri, GridPT3, blunder_param,mt_minmaxheight,false);
						
						printf("mt_minmax %f %f\n",mt_minmaxheight[0],mt_minmaxheight[1]);
						
						VerticalLineLocus_blunder(ortho_ncc, INCC, MagImages_L,MagImages_R,DEM_resolution, im_resolution, LRPCs, RRPCs, LImagesize_ori, LImagesize,
												  LeftImage, RImagesize_ori, RImagesize, RightImage, Template_size, 
												  Size_Grid2D, param, GridPts, Grid_wgs, GridPT3,
												  NumofIAparam, ImageAdjust, Pyramid_step, Lstartpos, Rstartpos, 
												  save_filepath, tile_row, tile_col, iteration,count,Boundary,left_ori, right_ori,blunder_selected_level,true);
					}
					blunder_detection_TIN(pre_DEMtif,ortho_ncc,INCC,flag_blunder,count,&blunder_dh,filename_mps,ptslists, count_MPs,
										  trilists,count_tri, GridPT3, blunder_param, blunder_count,minz_mp,maxz_mp,mt_minmaxheight,IsRA,seedDEMsigma);
					
					free(ortho_ncc);
					free(INCC);
					free(trilists);
					
					printf("blunder detection end\tblunder_dh = %f\n",blunder_dh);

					if(count > 0)
						count_blunders = abs(blunder_count[1] - pre_count_blunder);
					else
						count_blunders = blunder_count[1];
						
					if(count_blunders < blunder_th_counts)
						flag = false;
					
					if(blunder_count[0] < Th_pts)
					{
						flag = false;
					}
					
					if(pre_count_blunder == blunder_count[1] && pre_count_MPs == blunder_count[0])
						flag = false;
					else {
						pre_count_MPs = blunder_count[0];
						pre_count_blunder = blunder_count[1];
					}

					count++;
					
					if(count >= max_count || count_tri <= 20)
						flag = false;
					
					printf("start TIN\n");
					if(Pyramid_step >= TIN_split_level || count_MPs < 10000)
					{
						D3DPOINT *input_tri_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),blunder_count[0]);
						uint32 *check_id		= (uint32*)calloc(sizeof(uint32),blunder_count[0]);
						FILE *pTri;
						
						int t_tri_counts = 0;
						for(i=0;i<count_MPs;i++)
						{
							if(flag)
							{
								if(ptslists[i].flag != 1 && ptslists[i].flag != 2)
								{
									input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
									input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
									check_id[t_tri_counts]			= i;
									
									t_tri_counts++;
								}
							}
							else
							{
								if(ptslists[i].flag != 1)
								{
									input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
									input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
									check_id[t_tri_counts]			= i;
									
									t_tri_counts++;
								}
							}

						}
						
						UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*t_tri_counts*4);
						
						sprintf(bufstr,"%s/txt/tri_aft_%d_%d.txt",filename_tri,flag_blunder,count);
						TINCreate(input_tri_pts,bufstr,t_tri_counts,t_trilists,min_max,&count_tri);
						
						free(input_tri_pts);
						
						trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
						i = 0;
						for(i=0;i<count_tri;i++)
						{
							trilists[i].m_X = check_id[t_trilists[i].m_X];
							trilists[i].m_Y = check_id[t_trilists[i].m_Y];
							trilists[i].m_Z = check_id[t_trilists[i].m_Z];
						}
						
						free(t_trilists);
						free(check_id);
					}
					else
					{
						trilists = TINgeneration(!flag,save_filepath, Pyramid_step, Size_Grid2D, im_resolution, DEM_resolution,
												 min_max,
												 Boundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
					}
					printf("end TIN\n");
					
					count_Results[0]	= count_MPs;
					count_Results[1]	= count_tri;
					
					printf("iter = %d\tGridsize = %f\tMPs = %d\tBlunder = %d\tcount_tri = %d\n",count,grid_resolution,blunder_count[0],count_blunders,count_tri);
					
					//blunder remove from TIN minmax height
					int floor_ieration = 2;
					if((flag == false || count >= max_count || count_tri <= 20) && flag_blunder && (Pyramid_step > 10))
					{
						SetHeightRange_blunder(minmaxHeight,ptslists, count_MPs, trilists,count_tri, GridPT3, blunder_param,mt_minmaxheight,false);
						
						double* ortho_ncc = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
						double* INCC = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
						
						VerticalLineLocus_blunder(ortho_ncc, INCC, MagImages_L,MagImages_R,DEM_resolution, im_resolution, LRPCs, RRPCs, LImagesize_ori, LImagesize,
												  LeftImage, RImagesize_ori, RImagesize, RightImage, Template_size,
												  Size_Grid2D, param, GridPts, Grid_wgs, GridPT3,
												  NumofIAparam, ImageAdjust, Pyramid_step, Lstartpos, Rstartpos,
												  save_filepath, tile_row, tile_col, iteration,count,Boundary,left_ori, right_ori,blunder_selected_level,false);
						
						free(trilists);
						free(INCC);
						
						int non_blunder_count = 0;
						double th_ref_ncc = 0.1 + (iteration-1)*0.05;
						if(th_ref_ncc > 0.5)
							th_ref_ncc = 0.5;
						for(int t_index=0;t_index<count_MPs;t_index++)
						{
							int t_col		  = (int)((ptslists[t_index].m_X - Boundary[0])/grid_resolution + 0.5);
							int t_row		  = (int)((ptslists[t_index].m_Y - Boundary[1])/grid_resolution + 0.5);
							int ref_index	  = Size_Grid2D.width*t_row + t_col;
							if(ptslists[t_index].m_Z < GridPT3[ref_index].t_minHeight || ptslists[t_index].m_Z > GridPT3[ref_index].t_maxHeight)
								ptslists[t_index].flag = 1;
							else
							{
								ptslists[t_index].flag = 0;
								non_blunder_count++;
							}
						}
						free(ortho_ncc);
						
						
						printf("2 start TIN\n");
						if(Pyramid_step >= TIN_split_level || count_MPs < 10000)
						{
							D3DPOINT *input_tri_pts = (D3DPOINT*)calloc(sizeof(D3DPOINT),non_blunder_count);
							uint32 *check_id		= (uint32*)calloc(sizeof(uint32),non_blunder_count);
							FILE *pTri;
							
							int t_tri_counts = 0;
							for(i=0;i<count_MPs;i++)
							{
								if(flag)
								{
									if(ptslists[i].flag != 1 && ptslists[i].flag != 2)
									{
										input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
										input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
										check_id[t_tri_counts]			= i;
										
										t_tri_counts++;
									}
								}
								else
								{
									if(ptslists[i].flag != 1)
									{
										input_tri_pts[t_tri_counts].m_X = ptslists[i].m_X;
										input_tri_pts[t_tri_counts].m_Y = ptslists[i].m_Y;
										check_id[t_tri_counts]			= i;
										
										t_tri_counts++;
									}
								}
							}
							
							UI3DPOINT* t_trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*t_tri_counts*4);
							
							sprintf(bufstr,"%s/txt/tri_aft_%d_%d.txt",filename_tri,flag_blunder,count);
							TINCreate(input_tri_pts,bufstr,t_tri_counts,t_trilists,min_max,&count_tri);
							
							free(input_tri_pts);
 
							trilists	= (UI3DPOINT*)malloc(sizeof(UI3DPOINT)*count_tri);
							i = 0;
							for(i=0;i<count_tri;i++)
							{
								trilists[i].m_X = check_id[t_trilists[i].m_X];
								trilists[i].m_Y = check_id[t_trilists[i].m_Y];
								trilists[i].m_Z = check_id[t_trilists[i].m_Z];
							}
							
							free(t_trilists);
							free(check_id);
						}
						else
						{
							trilists = TINgeneration(!flag,save_filepath, Pyramid_step, Size_Grid2D, im_resolution, DEM_resolution,
													 min_max,
													 Boundary, count_MPs, ptslists, &iter_row, &iter_col, &count_tri);
						}
						printf("2 end TIN\n");
						
						count_Results[0]	= count_MPs;
						count_Results[1]	= count_tri;
						
						char tt_str[500];
						sprintf(tt_str,"%s/txt/blunder_mps_%d_%d_%d.txt",filename_tri,Pyramid_step,flag_blunder,count);
						FILE *pfile = fopen(tt_str,"w");
						for(int tcnt=0;tcnt<count_MPs;tcnt++)
						{
							if(ptslists[tcnt].flag != 1)
							{
								fprintf(pfile,"%lf %lf %lf\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z);
							}
						}
						fclose(pfile);
					}
				}
				
				FILE *pfile = fopen(filename_mps,"w");
				int tcnt = 0;
				count_Results[0] = 0;
				for(tcnt=0;tcnt<count_MPs;tcnt++)
				{
					if(ptslists[tcnt].flag != 1)
					{
						fprintf(pfile,"%lf %lf %lf %d\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z,ptslists[tcnt].flag);
						count_Results[0]++;
					}
				}
				fclose(pfile);
				
				double mt_minmaxheight[2];
				SetHeightRange_blunder(minmaxHeight,ptslists, count_Results[0], trilists,count_tri, GridPT3, blunder_param,mt_minmaxheight,false);
				
				free(trilists);
				
				if(flag_blunder == true)
					count = 50;
				else {
					count = 40;
				}

				double* ortho_ncc = (double*)calloc(sizeof(double),Size_Grid2D.height*Size_Grid2D.width);
				double* INCC = (double*)calloc(sizeof(double),Size_Grid2D.height*Size_Grid2D.width);
				
				VerticalLineLocus_blunder(ortho_ncc, INCC, MagImages_L,MagImages_R,DEM_resolution, im_resolution, LRPCs, RRPCs, LImagesize_ori, LImagesize,
										  LeftImage, RImagesize_ori, RImagesize, RightImage, Template_size, 
										  Size_Grid2D, param, GridPts, Grid_wgs, GridPT3,
										  NumofIAparam, ImageAdjust, Pyramid_step, Lstartpos, Rstartpos, 
										  save_filepath, tile_row, tile_col, iteration,count,Boundary,left_ori, right_ori,blunder_selected_level,false);
				free(ortho_ncc);
				free(INCC);
				
			}
			else
			{
				FILE *pfile = fopen(filename_mps,"w");
				int tcnt = 0;
				for(tcnt=0;tcnt<count_MPs;tcnt++)
				{
					if(ptslists[tcnt].flag != 1)
					{
						fprintf(pfile,"%lf %lf %lf %d\n",ptslists[tcnt].m_X,ptslists[tcnt].m_Y,ptslists[tcnt].m_Z,ptslists[tcnt].flag);
					}
					
				}
				count_Results[0] = count_MPs;
				fclose(pfile);
			}
			free(ptslists);

		}
		fclose(survey);
		
	}
	else
		*p_flag		= false;

	return count;
}

int DecisionMPs_setheight(bool flag_blunder, int count_MPs_input, double* Boundary, UGRID *GridPT3, uint8 Pyramid_step, double grid_resolution,
						  uint8 iteration, CSize Size_Grid2D, char *filename_mps_pre, char *filename_tri, double Hinterval, 
						  bool *p_flag, double *pre_3sigma, double *pre_mean, int *count_Results, double *minz_mp, double *maxz_mp, double *minmaxHeight,
						  uint16 *MagImages_L,uint16 *MagImages_R,double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs,
						  CSize LImagesize_ori, CSize LImagesize, uint16* LeftImage, CSize RImagesize_ori, CSize RImagesize, uint16* RightImage, uint8 Template_size, 
						  TransParam param, D2DPOINT* Grid_wgs,D2DPOINT* GridPts,
						  uint8 NumofIAparam, double* ImageAdjust, D2DPOINT Lstartpos, D2DPOINT Rstartpos, 
						  char* save_filepath, uint8 tile_row, uint8 tile_col,D3DPOINT *ptslists, UI3DPOINT *trilists,int numoftri,uint8* left_ori, uint8* right_ori,int blunder_selected_level)
{
	char bufstr[500];
	uint16 count			= 0;
	int count_MPs		= count_MPs_input;
	double gridspace			= grid_resolution;
	int Th_pts;
	int i=0;
	BL blunder_param;
	double blunder_dh = 0;
	double mt_minmaxheight[2];
	double* ortho_ncc = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
	double* INCC = (double*)calloc(Size_Grid2D.height*Size_Grid2D.width,sizeof(double));
	
	blunder_param.Boundary	= Boundary;
	blunder_param.gridspace	= grid_resolution;
	blunder_param.height_check_flag = true;
	blunder_param.Hinterval = Hinterval;
	blunder_param.iteration = iteration;
	blunder_param.Pyramid_step = Pyramid_step;
	
	blunder_param.Size_Grid2D.width = Size_Grid2D.width;
	blunder_param.Size_Grid2D.height = Size_Grid2D.height;
	
	SetHeightRange_blunder(minmaxHeight,ptslists, count_MPs, trilists,numoftri, GridPT3, blunder_param,mt_minmaxheight,false);
	
	VerticalLineLocus_blunder(ortho_ncc,INCC, MagImages_L,MagImages_R,DEM_resolution, im_resolution, LRPCs, RRPCs, LImagesize_ori, LImagesize,
							  LeftImage, RImagesize_ori, RImagesize, RightImage, Template_size, 
							  Size_Grid2D, param, GridPts, Grid_wgs, GridPT3,
							  NumofIAparam, ImageAdjust, Pyramid_step, Lstartpos, Rstartpos, 
							  save_filepath, tile_row, tile_col, iteration,100,Boundary,left_ori, right_ori,blunder_selected_level,false);
	free(ortho_ncc);
	free(INCC);
	
	return count;
}

void TINCreate(D3DPOINT *ptslists, char *filename_tri,int numofpts,UI3DPOINT* trilists,double min_max[],int *count_tri)
{
double minX_ptslists = min_max[0];
double minY_ptslists = min_max[1];
double maxX_ptslists = min_max[2];
double maxY_ptslists = min_max[3];
#ifdef VLIBRARY
	double distX_ptslists, distY_ptslists;
	double Scale_ptslists = 1000;
	Site *(*next)();
	
	*count_tri = 0;

	initializeVoronoi();	
	D3DPOINT *scaled_ptslists;

	scaled_ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*numofpts);

	distX_ptslists = maxX_ptslists - minX_ptslists;
	distY_ptslists = maxY_ptslists - minY_ptslists;
	for(int i=0;i<numofpts;i++)
	{
		scaled_ptslists[i].m_X = (ptslists[i].m_X - minX_ptslists)/distX_ptslists*Scale_ptslists;
		scaled_ptslists[i].m_Y = (ptslists[i].m_Y - minY_ptslists)/distY_ptslists*Scale_ptslists;
	}

	readsites(scaled_ptslists,numofpts);
	next = nextone;
	geominit();
	voronoi(next,trilists,&(*count_tri));
	free(scaled_ptslists);
	free_all();

#elif TRILIBRARY
	//clock_t start = clock(), mid1, mid2, mid3, end;

	struct triangulateio in, out;
    int x;
	double v12[2], v13[2], Normal;
    D3DPOINT pt1, pt2, pt3;
	
	if(numofpts < 3){
        *count_tri = 0;
        return;
    }
	D3DPOINT *shifted_ptslists;
	shifted_ptslists = (D3DPOINT*)malloc(sizeof(D3DPOINT)*numofpts);
	
	double midX_ptslists = (maxX_ptslists + minX_ptslists) / 2.0;
	double midY_ptslists = (maxY_ptslists + minY_ptslists) / 2.0;
	for(int i=0;i<numofpts;i++)
	{
		shifted_ptslists[i].m_X = ptslists[i].m_X - midX_ptslists;
		shifted_ptslists[i].m_Y = ptslists[i].m_Y - midY_ptslists;
	}

	/* Define input points. */

    in.numberofpoints = numofpts;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

	#pragma omp parallel for private(x) schedule(guided)
    for (x = 0; x < numofpts; x++){
        in.pointlist[2*x] = shifted_ptslists[x].m_X;
        in.pointlist[2*x+1] = shifted_ptslists[x].m_Y;
    }

	in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                          in.numberofpointattributes *
                                          sizeof(REAL));
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));

	/* Make necessary initializations so that Triangle can return a */
    /*   triangulation in `out'.  */

    out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
    //mid1 = clock();
	triangulate("zQ", &in, &out, (struct triangulateio *) NULL);
	//mid2 = clock();
	/* Transfer output to trilists */
	#pragma omp parallel for private(x) schedule(guided)
	for (x = 0; x < out.numberoftriangles; x++) {
        trilists[x].m_X = (uint32)out.trianglelist[x * 3 + 0];
        trilists[x].m_Z = (uint32)out.trianglelist[x * 3 + 1];
        trilists[x].m_Y = (uint32)out.trianglelist[x * 3 + 2];
/*
        pt1 = ptslists[trilists[x].m_X];
        pt2 = ptslists[trilists[x].m_Y];
        pt3 = ptslists[trilists[x].m_Z];

        v12[0]  = pt2.m_X-pt1.m_X;
        v12[1]  = pt2.m_Y-pt1.m_Y;

        v13[0]  = pt3.m_X-pt1.m_X;
        v13[1]  = pt3.m_Y-pt1.m_Y;

        Normal = v12[0]*v13[1] - v12[1]*v13[0];
        if (Normal > 0) {
        	trilists[x].m_Z = (uint32)out.trianglelist[x * 3 + 1];
        	trilists[x].m_Y = (uint32)out.trianglelist[x * 3 + 2];
        }*/
	}
	
	*count_tri = out.numberoftriangles;
	//mid3 = clock();
	/* Free all allocated arrays, including those allocated by Triangle. */

    free(in.pointlist);
    free(in.pointattributelist);
    free(in.pointmarkerlist);
    free(in.regionlist);
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.pointmarkerlist);
    free(out.trianglelist);
    free(out.triangleattributelist);
	free(out.neighborlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(shifted_ptslists);

	//end = clock();
/*
    int msec = (mid1 - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken to initalize %d seconds %d milliseconds", msec/1000, msec%1000);
	msec = (mid2 - mid1) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken to triangulate %d seconds %d milliseconds", msec/1000, msec%1000);
	msec = (mid3 - mid2) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken to transfer out %d seconds %d milliseconds", msec/1000, msec%1000);
	msec = (end - mid3) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken to free alloc %d seconds %d milliseconds", msec/1000, msec%1000);
	msec = (end - start) * 1000 / CLOCKS_PER_SEC;
    printf("Time taken total %d seconds %d milliseconds", msec/1000, msec%1000);
*/
#endif
}

bool blunder_detection_TIN(int pre_DEMtif,double* ortho_ncc, double* INCC, bool flag_blunder,uint16 count_bl,double* blunder_dh,char *file_pts,
						   D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *Gridpts, BL BL_param, uint32 *blunder_count,
						   double *minz_mp, double *maxz_mp, double *minmaxHeight, int IsRA,double seedDEMsigma)
{
	uint8 pyramid_step, iteration;
	pyramid_step	= BL_param.Pyramid_step;
	iteration		= BL_param.iteration;
	double gridspace;
	gridspace = BL_param.gridspace;
	
	printf("count bo = %d\n",count_bl);
	if(!(pyramid_step == 0 && iteration == 3))
	{
		int count_mp;
		int i, index;
		int num_triangles,num_points;
		int tcnt;
		double count_height_th, height_th;
		double lw_3sigma, up_3sigma;
		bool check_dh		= false;
		double sum_data2	= 0.0;
		double sum_data	= 0.0;
		double mean,std;
		double GSD, height_th_1, height_th_2;
		int dh_count = 0;
		double small_th = 0;
		int count_small = 0;
		uint8 max_nodes;
		double *boundary;
		CSize gridsize;
		//uint32 **savenode;
		//savenode will now be a pointer to a uint
		uint32* savenode;
		uint32 *nodecount;
		uint32 *hdiffbin;
		uint32 th80 = 1000;
		uint32 total_dh = 0;
		bool check_total_dh = false;
		uint32 hdiffcount = (uint32)BL_param.Hinterval;
		int th_count_bl = 3;
		
		max_nodes		= 30;
		num_triangles	= numOfTri;
		num_points	= numOfPts;
		//savenode		= (uint32**)calloc(num_points,sizeof(uint32*));
		//make savenode a giant array with row size max_nodes
		savenode = malloc((long)sizeof(uint32)*(long)max_nodes*(long)num_points);
		nodecount		= (uint32*)calloc(num_points,sizeof(uint32));
		hdiffbin	  = (uint32*)calloc(hdiffcount+1,sizeof(uint32));
		i = 0;
		/*for(i=0;i<num_points;i++)
		{
			savenode[i]		= (uint32*)calloc(max_nodes,sizeof(uint32));
			}*/
		if(pyramid_step <= 4 && pyramid_step >= 3)
		{
			//if(iteration <= 4)
				check_dh		= true;
		}
		else if(pyramid_step == 2)
		{
			if(BL_param.iteration <= 4)
			{
				check_dh		= true;
			}
		}
		else if(pyramid_step <= 1)
		{
		        if(iteration <= 2)
			{
		                check_dh		= true;
			}
		}
		
		boundary	= BL_param.Boundary;
		gridsize.width	= BL_param.Size_Grid2D.width;
		gridsize.height	= BL_param.Size_Grid2D.height;
			
		if(pyramid_step >= 1)
			GSD	= gridspace;
		else
			GSD	= 0.5*pow(2.0,pyramid_step);

		height_th_1= BL_param.Hinterval;
		height_th_2= GSD*10;
			
		if(IsRA == 1)
			height_th_2= GSD*30;
			
		if(height_th_1 < height_th_2)
			height_th		= height_th_1;
		else
			height_th		= height_th_2;
		
		if(pyramid_step == 2)
			small_th = 3.0;
		else
			small_th = 2.0;
			
		small_th	= gridspace;
		
		double sum_oncc = 0;
		double sum2_oncc = 0;
		int total_oncc_count = 0;
 
		//#pragma omp parallel for shared(count_bl,pts,tris,savenode,nodecount,num_triangles,height_th,num_points,hdiffbin) private(tcnt) reduction(+:sum_data2, sum_data,dh_count)
		for(tcnt=0;tcnt<num_triangles;tcnt++)
		{
			double dh1,dh2,dh3;
			int node1, node2, node3;
			D3DPOINT pt0,pt1,pt2;
			bool check_pt_index = true;
			int t_col,t_row;
			node1	= tris[tcnt].m_X;
			node2	= tris[tcnt].m_Y;
			node3	= tris[tcnt].m_Z;

			if(node1 < num_points)
			{
				if(nodecount[node1] < max_nodes)
				{
					//savenode[node1][nodecount[node1]] = tcnt;
					savenode[(long)node1*(long)max_nodes+(long)nodecount[node1]] = tcnt;
                    //#pragma omp critical
					nodecount[node1]++;
				}
				pt0		= pts[node1];
				
				t_col		  = (int)((pt0.m_X - boundary[0])/gridspace + 0.5);
				t_row		  = (int)((pt0.m_Y - boundary[1])/gridspace + 0.5);
				int pt0_index	  = gridsize.width*t_row + t_col;
				sum_oncc += ortho_ncc[pt0_index];
				sum2_oncc += ortho_ncc[pt0_index]* ortho_ncc[pt0_index];
				total_oncc_count++;
			}
			else
				check_pt_index = false;
			  
			if(node2 < num_points)
			{
				if(nodecount[node2] < max_nodes)
				{
					//savenode[node2][nodecount[node2]] = tcnt;
					savenode[(long)node2*(long)max_nodes+(long)nodecount[node2]] = tcnt;
                    //#pragma omp critical
					nodecount[node2]++;
				}
				pt1		= pts[node2];
				
				t_col		  = (int)((pt1.m_X - boundary[0])/gridspace + 0.5);
				t_row		  = (int)((pt1.m_Y - boundary[1])/gridspace + 0.5);
				int pt1_index	  = gridsize.width*t_row + t_col;
				sum_oncc += ortho_ncc[pt1_index];
				sum2_oncc += ortho_ncc[pt1_index]* ortho_ncc[pt1_index];
				total_oncc_count++;
			}
			else
				check_pt_index = false;
			  
			if(node3 < num_points)
			{
				if(nodecount[node3] < max_nodes)
				{
					//savenode[node3][nodecount[node3]] = tcnt;
					savenode[(long)node3*(long)max_nodes+(long)nodecount[node3]] = tcnt;
					//#pragma omp critical
					nodecount[node3]++;
				}
				pt2		= pts[node3];
				
				t_col		  = (int)((pt2.m_X - boundary[0])/gridspace + 0.5);
				t_row		  = (int)((pt2.m_Y - boundary[1])/gridspace + 0.5);
				int pt2_index	  = gridsize.width*t_row + t_col;
				sum_oncc += ortho_ncc[pt2_index];
				sum2_oncc += ortho_ncc[pt2_index]* ortho_ncc[pt2_index];
				total_oncc_count++;
			}
			else
				check_pt_index = false;

			if(check_pt_index)
			{
				dh1		= (double)fabs(pt0.m_Z - pt1.m_Z);
				sum_data += dh1;
				sum_data2 += dh1*dh1;
				dh_count++;
				
				dh2		= (double)fabs(pt0.m_Z - pt2.m_Z);
				sum_data += dh2;
				sum_data2 += dh2*dh2;
				dh_count++;
		   
				dh3		= (double)fabs(pt1.m_Z - pt2.m_Z);
				sum_data += dh3;
				sum_data2 += dh3*dh3;
				dh_count++;
				
				dh1	  = (uint32)dh1;
				if (dh1 < hdiffcount)
				{
					//#pragma omp critical
					hdiffbin[(uint32)dh1]++;
				}
				dh2	  = (uint32)dh2;
				if (dh2 < hdiffcount)
				{
					//#pragma omp critical
					hdiffbin[(uint32)dh2]++;
				}
				dh3	  = (uint32)dh3;
				if (dh3 < hdiffcount)
				{
					//#pragma omp critical
					hdiffbin[(uint32)dh3]++;
				}
			}
		}
		
		double oncc_mean	=  sum_oncc/total_oncc_count;
		double oncc_std	=	sqrt(fabs(sum2_oncc - (sum_oncc)*(sum_oncc)/total_oncc_count)/total_oncc_count);
		
		printf("oncc mean %f\tstd%f\n",oncc_mean,oncc_std);
		
		int blunder_pyramid_step;
		blunder_pyramid_step = 3;
			
		if(pyramid_step >= blunder_pyramid_step)
		{
			if(check_dh)
			{
				mean	=	sum_data/dh_count;
				std	=	sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
				lw_3sigma		= mean - 3*std;
				up_3sigma		= mean + 3*std;
			}
			else
			{
				lw_3sigma		= -100;
				up_3sigma		= height_th;	
			}
			
            printf("height_th = %f\t3sigma = %f\theight_th_2 = %f\n",height_th,up_3sigma,height_th_2);
            
			if(pyramid_step	 <= 4)
			{
				/*if(pyramid_step <= 4 && pyramid_step >= 3 && iteration <= 4)
				{
					if(up_3sigma < height_th)
						height_th = up_3sigma;
				}
				else if(pyramid_step == 2 && iteration <= 2)
				{
					if(up_3sigma < height_th)
						height_th = up_3sigma;
				}
				else if(iteration <= 2)
				{
					if(up_3sigma < height_th)
						height_th = up_3sigma;
				}
                */
                if(up_3sigma < height_th)
                    height_th = up_3sigma;
			}
		}
		else
		{
			total_dh = 0;
			tcnt = 0;
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
				
			mean	=	sum_data/dh_count;
			std	=	sqrt((sum_data2 - (sum_data)*(sum_data)/dh_count)/dh_count);
			lw_3sigma		= mean - 3*std;
			up_3sigma		= mean + 3*std;

			if(up_3sigma < th80)
				height_th = up_3sigma;
			else {
				height_th = th80;
			}

			if(height_th > 50)
				height_th = 50;
			
			printf("height_th = %f\t3sigma = %f\tth90 = %d\theight_th_2 = %f\n",height_th,up_3sigma,th80,height_th_2);
		}
		
		free(hdiffbin);
	  
		count_height_th = GSD*0.7;
		if(gridspace >= 2)
		{
			if(count_height_th < 3)
				count_height_th = 3;
		}
		printf("height_th %f\tcount_height_th %f\titeration %d\n",height_th,count_height_th,iteration);
		  
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
			ortho_ncc_th = 0.7;
		else if(pyramid_step >= 3)
			ortho_ncc_th = 0.6 - (iteration - 1)*0.02;
		else if(pyramid_step == 2)
			ortho_ncc_th = 0.5 - (iteration - 1)*0.02;
		else if(pyramid_step == 1)
			ortho_ncc_th = 0.4;
		else
			ortho_ncc_th = 0.3;

		
		if(pre_DEMtif && pyramid_step == 2 && seedDEMsigma <= 20)
			ortho_ncc_th = 0.3;
		
		double ortho_ancc_th = 100.;
		double th_ref_ncc = 0.1 + (iteration-1)*0.05;
		if(th_ref_ncc > ortho_ncc_th)
			th_ref_ncc = ortho_ncc_th;
		
		printf("ortho_ncc_th %f\t ortho_ancc_th %f\t th_ref_ncc %f\n",ortho_ncc_th,ortho_ancc_th,th_ref_ncc);
		
	   
#pragma omp parallel for private(index) reduction(+:count_small) schedule(guided)
		for(index=0;index<num_points;index++)
		{
			if(pts[index].flag != 1)
			{
				int count_th_positive	= 0;
				int count_th_negative	= 0;
				int count_small_height	= 0;
				int count = 0;
				int iter = 0;
				int max_iter;
				int t_col,t_row;
				  
				bool check_neigh = false;
				D3DPOINT ref_index_pt;
				int ref_index;
				  
				if(nodecount[index] < max_nodes)
					max_iter	= nodecount[index];
				else
					max_iter	= max_nodes;

				ref_index_pt = pts[index];
				t_col		  = (int)((ref_index_pt.m_X - boundary[0])/gridspace + 0.5);
				t_row		  = (int)((ref_index_pt.m_Y - boundary[1])/gridspace + 0.5);
				ref_index	  = gridsize.width*t_row + t_col;

				Gridpts[ref_index].anchor_flag = 0;
				  
				if(!IsRA)
				{
					if(ortho_ncc[ref_index] < th_ref_ncc && pyramid_step >= 2)
						pts[index].flag = 1;
				}
				if(pyramid_step >= 3 && iteration >= 4)
				{
					if(ortho_ncc[ref_index] >= 0.95 && flag_blunder)
					{
						Gridpts[ref_index].anchor_flag = 1;
					}
				}
				  
				for(iter = 0 ; iter < max_iter ; iter++)
				{
					uint32 reference_index = 0;
					int t_count = 0;
					double U[3], V[3];
					int target_index_0 = 0, target_index_1 = 0;
					  
					double dh1, dh2, dh3, ddh_1, ddh_2, ddh_3;
					double N[3];
					double norm, angle = 0;
					double dh_max;
					D3DPOINT pt0,pt1,pt2;
					int kk;
					bool check_index = false;
					uint32 temp_tri[3];
					  
					//if(savenode[index][iter] < num_triangles && savenode[index][iter] > 0)
					if(savenode[(long)index*(long)max_nodes+(long)iter] < num_triangles && savenode[(long)index*(long)max_nodes+(long)iter] > 0)
					{
						/*
						temp_tri[0]	= tris[savenode[index][iter]].m_X;
						temp_tri[1]	= tris[savenode[index][iter]].m_Y;
						temp_tri[2]	= tris[savenode[index][iter]].m_Z;*/
						temp_tri[0] = tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_X;
						temp_tri[1] = tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Y;
						temp_tri[2] = tris[savenode[(long)index*(long)max_nodes+(long)iter]].m_Z;

						for(kk=0;kk<3;kk++)
						{
							if(temp_tri[kk] == index)
								reference_index = index;
							else
							{
								if(!check_index)
								{
									target_index_0 = temp_tri[kk];
									check_index	   = true;
								}
								else
									target_index_1 = temp_tri[kk];
							}	
						}
						  
						if(reference_index < num_points && target_index_0 < num_points && target_index_1 < num_points &&
						   target_index_0 >= 0 && target_index_1 >= 0)
						{
							pt0	= pts[reference_index];
							pt1	= pts[target_index_0];
							pt2	= pts[target_index_1];
							  
							dh1 = pt0.m_Z - pt1.m_Z;
							dh2 = pt0.m_Z - pt2.m_Z;
							dh3 = pt1.m_Z - pt2.m_Z;

							ddh_1	  = fabs(dh1) - fabs(dh2);
							ddh_2	  = fabs(dh1) - fabs(dh3);
							ddh_3	  = fabs(dh2) - fabs(dh3);
					   
							if(pyramid_step <= 2 && fabs(dh1) <= small_th && fabs(dh2) <= small_th && fabs(dh3) <= small_th)
							{
								count_small_height++;
							}
							  
							U[0] = pt1.m_X - pt0.m_X;
							V[0] = pt2.m_X - pt0.m_X;
							U[1] = pt1.m_Y - pt0.m_Y;
							V[1] = pt2.m_Y - pt0.m_Y;
							U[2] = pt1.m_Z - pt0.m_Z;
							V[2] = pt2.m_Z - pt0.m_Z;
								
							N[0]		=	U[1]*V[2] - V[1]*U[2];
							N[1]		= -(U[0]*V[2] - V[0]*U[2]);
							N[2]		=	U[0]*V[1] - V[0]*U[1];

							norm = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
							if (norm != 0)
							{
								angle =acos(fabs(N[2])/norm)*RadToDeg;
							}
								  
							if(angle <= 0 && angle >= -90)
								angle = fabs(angle);
							else if(angle <= -270 && angle >= -360)
								angle = 360 + angle;
							else if(angle >= 270 && angle <= 360)
								angle = 360 - angle;
								  
							if(fabs(dh1) > fabs(dh2))
								dh_max = fabs(dh1);
							else
								dh_max = fabs(dh2);

							if(dh1 >= 0 && dh2 >= 0 && angle > 30)
								count_th_positive ++;
							if(dh1 <  0 && dh2 < 0	&& angle > 30)
								count_th_negative ++;
							  
							bool check_match;

							check_match	   = false;
							if(pyramid_step	 > 1)
							{
								check_match		= true;
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
							  
							if((fabs(ddh_1) > height_th || fabs(ddh_2) > height_th || fabs(ddh_3) > height_th) && check_match)
							{
								int order[3]	= {reference_index,target_index_0,target_index_1};
								double height[3] = {pt0.m_Z,pt1.m_Z,pt2.m_Z};
								double h1,h2,dh;
								int t_o_min,t_o_max,t_o_mid;

								// all node check with height difference.
								int Col			= (int)((pt0.m_X - boundary[0])/gridspace + 0.5);
								int Row			= (int)((pt0.m_Y- boundary[1])/gridspace + 0.5);
								int Index[3];
								uint8 m_bMatch_1, m_bMatch_2, m_bMatch_3;
								
								Index[0]	 = gridsize.width*Row + Col;
								if(Index[0] < (int)(gridsize.width*gridsize.height) && Index[0] >= 0)
									m_bMatch_1	= Gridpts[Index[0]].Matched_flag;
								else
									m_bMatch_1	= 0;

								Col				= (int)((pt1.m_X - boundary[0])/gridspace + 0.5);
								Row				= (int)((pt1.m_Y - boundary[1])/gridspace + 0.5);	 
								Index[1]	 = gridsize.width*Row + Col;
								if(Index[1] < (int)(gridsize.width*gridsize.height) && Index[1] >= 0)
									m_bMatch_2	= Gridpts[Index[1]].Matched_flag;
								else
									m_bMatch_2	= 0;

								Col				= (int)((pt2.m_X - boundary[0])/gridspace + 0.5);
								Row				= (int)((pt2.m_Y - boundary[1])/gridspace + 0.5);	 
								Index[2]	 = gridsize.width*Row + Col;
								if(Index[2] < (int)(gridsize.width*gridsize.height) && Index[2] >= 0)
									m_bMatch_3	= Gridpts[Index[2]].Matched_flag;
								else
									m_bMatch_3	= 0;

								if(check_match)
								{
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

									h1		  = height[t_o_mid] - height[t_o_min];
									h2		  = height[t_o_max] - height[t_o_mid];
									dh		  = h1 - h2;
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
						{
							pts[index].flag = 1;
						}
					}
					if((count_th_positive >= (int)(count*0.7 + 0.5) || count_th_negative >= (int)(count*0.7 + 0.5)) && count > 1)
					{
						pts[index].flag = 1;
					}
				}
				else
				{
					if(check_neigh == true)
					{
						if(!flag_blunder)
						{
							//#pragma omp critical
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
						//#pragma omp critical
						{
							if(flag_blunder)
							{
								if(pyramid_step >= 1)
								{
									if(ortho_ncc[ref_index] < 0.6)
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
				}
			}
		}
		 
		printf("count small %d\n",count_small);
		
		for(tcnt=0;tcnt<num_points;tcnt++)
		{
			if(pts[tcnt].flag == 1)
			{
				blunder_count[1]++;
			}
			else if(pts[tcnt].flag == 3)
			{
			        blunder_count[1]++;
			        pts[tcnt].flag = 1;
			}
			else
			{
				blunder_count[0]++;
				  
				if(*minz_mp > pts[tcnt].m_Z)
					*minz_mp		= pts[tcnt].m_Z;
				if(*maxz_mp < pts[tcnt].m_Z)
					*maxz_mp		= pts[tcnt].m_Z;
			}
			//free(savenode[tcnt]);
		}
		free(savenode);
		free(nodecount);
	}
  
	return true;
}

int SetttingFlagOfGrid(double *subBoundary,UGRID *GridPT3, uint8 Pyramid_step,double grid_resolution,uint8 iteration,
					   CSize Size_Grid2D,char *filename_mps_anchor,char *filename_mps_aft,int count_MPs_anchor,int count_MPs_blunder, char *filename_mps)
{
	int total_count = 0;
	double X,Y,Z;
	int t_flag;
	int i = 0;
	FILE *fid = fopen(filename_mps_anchor,"r");
	FILE *fid_all = fopen(filename_mps,"w");
	
	while( i < count_MPs_anchor && (fscanf(fid,"%lf %lf %lf %d\n",&X,&Y,&Z,&t_flag)) != EOF )
	{
		int grid_index;
		int t_col, t_row;
		
		t_col		  = (int)((X - subBoundary[0])/grid_resolution + 0.5);
		t_row		  = (int)((Y - subBoundary[1])/grid_resolution + 0.5);
		grid_index	   = Size_Grid2D.width*t_row + t_col;
		if(grid_index >= 0 && grid_index < Size_Grid2D.width*Size_Grid2D.height && t_col >=0 && t_col < Size_Grid2D.width && t_row >= 0 && t_row < Size_Grid2D.height)
		{
			GridPT3[grid_index].anchor_flag = 1;
		}
		
		i++;
		
	}
	fclose(fid);
	
	i=0;
	fid = fopen(filename_mps_aft,"r");
	while( i < count_MPs_blunder && (fscanf(fid,"%lf %lf %lf %d\n",&X,&Y,&Z,&t_flag)) != EOF )
	{
		int grid_index;
		int t_col, t_row;
		
		t_col		  = (int)((X - subBoundary[0])/grid_resolution + 0.5);
		t_row		  = (int)((Y - subBoundary[1])/grid_resolution + 0.5);
		grid_index	   = Size_Grid2D.width*t_row + t_col;
		if(grid_index >= 0 && grid_index < Size_Grid2D.width*Size_Grid2D.height && t_col >=0 && t_col < Size_Grid2D.width && t_row >= 0 && t_row < Size_Grid2D.height)
		{
			fprintf(fid_all,"%lf %lf %lf %d\n",X,Y,Z,t_flag);
			total_count++;
			if(GridPT3[grid_index].anchor_flag != 1)
			{
				GridPT3[grid_index].anchor_flag = 2;
				
			}
		}
		
		i++;
		
	}
	fclose(fid);
	fclose(fid_all);
	
	return total_count;
}

int Ortho_blunder(D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, bool update_flag,double *minH_grid, double *maxH_grid, BL BL_param,
				  uint16 *MagImages_L,uint16 *MagImages_R, uint16* LeftImage, uint16* RightImage,
				  double DEM_resolution, double im_resolution, double** LRPCs, double** RRPCs, 
				  CSize LImagesize,	 CSize RImagesize, CSize Size_Grid2D, TransParam param, uint8 NumofIAparam, 
				  double* ImageAdjust, double* minmaxHeight, uint8 Pyramid_step, double meters_per_pixel,
				  D2DPOINT Lstartpos, D2DPOINT Rstartpos, uint8 iteration,	UGRID *GridPT3, char *filename_mps, char *save_path)
{
	uint32 num_triangles;
	int i, tcnt;
	double gridspace;
	double *boundary;
	CSize gridsize;
	uint8 pyramid_step = Pyramid_step;
	uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;

	uint8* tris_check;
	
	bool check_stop_TIN;
	int max_count = 100;
	int while_count = 0;
	
	num_triangles=numOfTri;
	pyramid_step	= BL_param.Pyramid_step;
	iteration		= BL_param.iteration;
	gridspace = BL_param.gridspace;
	boundary	= BL_param.Boundary;
	gridsize.width	= BL_param.Size_Grid2D.width;
	gridsize.height	= BL_param.Size_Grid2D.height;
	TIN_Grid_Size_X= gridsize.width;
	TIN_Grid_Size_Y= gridsize.height;

	tris_check = (uint8*)calloc(num_triangles,sizeof(uint8));
	check_stop_TIN = false;
	
	while(!check_stop_TIN && while_count < max_count)
	{
		bool check_ortho_cal = false;
	  
		while_count++;
	  
		check_stop_TIN = true;
	  
//#pragma omp parallel for shared(tris_check,while_count,TIN_Grid_Size_X,TIN_Grid_Size_Y,gridspace,pts,tris,numOfPts,check_ortho_cal,check_stop_TIN,MagImages_L,MagImages_R, LeftImage, RightImage,DEM_resolution, im_resolution, LRPCs, RRPCs,LImagesize,	RImagesize, Size_Grid2D, param, NumofIAparam, ImageAdjust, minmaxHeight, Pyramid_step, meters_per_pixel, Lstartpos, Rstartpos, iteration, GridPT3,boundary) private(tcnt)
		for(tcnt=0;tcnt<(int)(num_triangles);tcnt++)
		{
//#pragma omp ordered
			{
				if(tris_check[tcnt] == 0)
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
			
					t_tri	= tris[tcnt];
					pdex0 = t_tri.m_X;
					pdex1 = t_tri.m_Y;
					pdex2 = t_tri.m_Z;
			  
					if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
					{
						int node1_index,node2_index,node3_index;
						uint8 node1_F,node2_F,node3_F;
						int node_f_count = 0;
				  
						pt0	= pts[pdex0];
						pt1	= pts[pdex1];
						pt2	= pts[pdex2];

						TriP1[0]		= pt0.m_X;
						TriP2[0]		= pt1.m_X;
						TriP3[0]		= pt2.m_X;

						TriP1[1]		= pt0.m_Y;
						TriP2[1]		= pt1.m_Y;
						TriP3[1]		= pt2.m_Y;

						TriP1[2]		= pt0.m_Z;
						TriP2[2]		= pt1.m_Z;
						TriP3[2]		= pt2.m_Z;

						node1_index = (int)((TriP1[1] - boundary[1])/gridspace + 0.5)*TIN_Grid_Size_X + (int)((TriP1[0] - boundary[0])/gridspace + 0.5);
						node2_index = (int)((TriP2[1] - boundary[1])/gridspace + 0.5)*TIN_Grid_Size_X + (int)((TriP2[0] - boundary[0])/gridspace + 0.5);
						node3_index = (int)((TriP3[1] - boundary[1])/gridspace + 0.5)*TIN_Grid_Size_X + (int)((TriP3[0] - boundary[0])/gridspace + 0.5);
				  
						node1_F		=  GridPT3[node1_index].anchor_flag;
						node2_F		=  GridPT3[node2_index].anchor_flag;
						node3_F		=  GridPT3[node3_index].anchor_flag;

						if(node1_F == 1 || node1_F == 3)
							node_f_count++;
						if(node2_F == 1 || node2_F == 3)
							node_f_count++;
						if(node3_F == 1 || node3_F == 3)
							node_f_count++;
				  
						if(node_f_count == 3)
							//	#pragma omp critical
						{
							tris_check[tcnt] = 1;
						}
				  
						if(while_count <= 0)
						{
							if(node_f_count == 2)
							{
								D3DPOINT ref1_pt, ref2_pt, target_pt;
								int target_index,ref1_index,ref2_index;
								int target_pt_index;
						  
								if(node1_F == 2)
								{
									ref1_pt	  = pt1;
									ref2_pt	= pt2;
									target_pt = pt0;
									target_index = node1_index;
									target_pt_index = pdex0;
								}
								else if(node2_F == 2)
								{
									ref1_pt	= pt0;
									ref2_pt	= pt2;
									target_pt = pt1;
									target_index = node2_index;
									target_pt_index = pdex1;
								}
								else
								{
									ref1_pt	= pt0;
									ref2_pt	= pt1;
									target_pt = pt2;
									target_index = node3_index;
									target_pt_index = pdex2;
								}
						  
								//#pragma omp critical
								{
									GridPT3[target_index].anchor_flag = 1;
								}
							}
					  
							//#pragma omp critical
							{
								check_stop_TIN = false;
								check_ortho_cal = true;
							}
						}
						else
						{
							if(node_f_count == 2)
							{
								//ortho matching process
								D3DPOINT ref1_pt, ref2_pt, target_pt;
								int target_index,ref1_index, ref2_index;
								int target_pt_index;
								double F_height;
					  
								if(node1_F == 2)
								{
									ref1_pt	  = pt1;
									ref2_pt	= pt2;
									target_pt = pt0;
									target_index = node1_index;
									ref1_index = node2_index;
									ref2_index = node3_index;
									target_pt_index = pdex0;
								}
								else if(node2_F == 2)
								{
									ref1_pt	= pt0;
									ref2_pt	= pt2;
									target_pt = pt1;
									target_index = node2_index;
									ref1_index = node1_index;
									ref2_index = node3_index;
									target_pt_index = pdex1;
								}
								else
								{
									ref1_pt	= pt0;
									ref2_pt	= pt1;
									target_pt = pt2;
									target_index = node3_index;
									ref1_index = node1_index;
									ref2_index = node2_index;
									target_pt_index = pdex2;
								}

								VerticalLineLocus_Ortho(&F_height,ref1_pt,ref2_pt,target_pt,
														MagImages_L,MagImages_R, LeftImage, RightImage,
														DEM_resolution, im_resolution, LRPCs, RRPCs, 
														LImagesize,	 RImagesize, Size_Grid2D, param, NumofIAparam, 
														ImageAdjust, minmaxHeight, Pyramid_step, meters_per_pixel,
														Lstartpos, Rstartpos, iteration,  GridPT3,target_index,ref1_index,ref2_index,boundary);
						  
						  
								if(F_height != -9999 )
								{
									//#pragma omp critical
									{
										pts[target_pt_index].m_Z = F_height;
								  
										GridPT3[target_index].anchor_flag = 3;
							  
										check_stop_TIN = false;
										check_ortho_cal = true;
										tris_check[tcnt] = 1;
									}
								}
							}
						}
					}
				}
			}
		}

		//#pragma omp critical
		{
			if(check_ortho_cal == false)
			{
				check_stop_TIN = true;
			}
		}
	}
	free(tris_check);
	
	return numOfPts;
}

UGRID* SetHeightRange(bool pre_DEMtif, double* minmaxHeight,int numOfPts, int numOfTri, UGRID *GridPT3, bool update_flag,
					  double *minH_grid, double *maxH_grid, BL BL_param,D3DPOINT *pts, UI3DPOINT *tris,int IsRA, double MPP, char* save_path, uint8 row, uint8 col,bool check_level_end,double seedDEMsigma)
{
	UGRID *result;
	
	uint32 num_triangles;
	int i, tcnt;
	double gridspace;
	double *boundary;
	CSize gridsize;
	uint8 pyramid_step, iteration;
	uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
	uint8* m_bHeight;
	double Total_Min_Z		=  100000;
	double Total_Max_Z		= -100000;
	double meters_per_pixel = MPP;
	
	printf("MPP of setheightrange = %f\n",MPP);
	
	int Col_C, Row_R;
	
	num_triangles			= numOfTri;
	pyramid_step	= BL_param.Pyramid_step;
	iteration		= BL_param.iteration;
	gridspace = BL_param.gridspace;
	boundary	= BL_param.Boundary;
	gridsize.width	= BL_param.Size_Grid2D.width;
	gridsize.height	= BL_param.Size_Grid2D.height;
	TIN_Grid_Size_X= gridsize.width;
	TIN_Grid_Size_Y= gridsize.height;
	
	double DEM_error		= meters_per_pixel*4;
	
	double BufferOfHeight	= DEM_error*pow(2.0,pyramid_step);

	if (pyramid_step == 1)
	{
		if(iteration >= 2)
			BufferOfHeight = meters_per_pixel*2;
		else {
			BufferOfHeight = meters_per_pixel*3;
		}
	}	
	else if(pyramid_step == 0)
	{
		if(iteration == 1)
			BufferOfHeight = meters_per_pixel*2;
		else {
			BufferOfHeight = meters_per_pixel;
		}
		
		if (BufferOfHeight < 0.5)
			BufferOfHeight = 0.5;
	}

	printf("BufferOfHeight = %f\n",BufferOfHeight);
	
    if(BufferOfHeight > 100)
        BufferOfHeight = 100;
    
    printf("BufferOfHeight = %f\n",BufferOfHeight);
    
	m_bHeight		= (uint8*)calloc(TIN_Grid_Size_Y*TIN_Grid_Size_X,sizeof(uint8));
	
//	#pragma omp parallel shared(GridPT3,pts,tris,num_triangles,m_bHeight,pyramid_step,iteration,gridspace,boundary,gridsize,TIN_Grid_Size_X,TIN_Grid_Size_Y,DEM_error) private(tcnt)
	{
//		#pragma omp for
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
			
			t_tri	= tris[tcnt];
			pdex0 = t_tri.m_X;
			pdex1 = t_tri.m_Y;
			pdex2 = t_tri.m_Z;
			
			if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
			{
				int node1_index, node2_index, node3_index;
				bool check_anchor_flag = false;
				pt0	= pts[pdex0];
				pt1	= pts[pdex1];
				pt2	= pts[pdex2];
				
				TriP1[0]		= pt0.m_X;
				TriP2[0]		= pt1.m_X;
				TriP3[0]		= pt2.m_X;
				
				TriP1[1]		= pt0.m_Y;
				TriP2[1]		= pt1.m_Y;
				TriP3[1]		= pt2.m_Y;
				
				TriP1[2]		= pt0.m_Z;
				TriP2[2]		= pt1.m_Z;
				TriP3[2]		= pt2.m_Z;
				
				double U[3], V[3], N[3];
				double norm, angle;
				
				U[0] = pt1.m_X - pt0.m_X;
				V[0] = pt2.m_X - pt0.m_X;
				U[1] = pt1.m_Y - pt0.m_Y;
				V[1] = pt2.m_Y - pt0.m_Y;
				U[2] = pt1.m_Z - pt0.m_Z;
				V[2] = pt2.m_Z - pt0.m_Z;
				
				N[0]		=	U[1]*V[2] - V[1]*U[2];
				N[1]		= -(U[0]*V[2] - V[0]*U[2]);
				N[2]		=	U[0]*V[1] - V[0]*U[1];
				
				norm  = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
				angle = acos(fabs(N[2])/norm)*180/3.141592;
				
                if(angle <= 0 && angle >= -90)
                    angle = fabs(angle);
                else if(angle <= -270 && angle >= -360)
                    angle = 360 + angle;
                else if(angle >= 270 && angle <= 360)
                    angle = 360 - angle;
                
                double angle_weight;
                /*if(angle <= 5)
                    angle_weight = 1.0 - 0.2*(5.0-angle)/5.0;
                else*/ if(angle <= 30)
                    angle_weight = 1.0 + 0.05*(angle)/5.0;
                else
                    angle_weight = 1.3;
                
				temp_MinZ = min(min(TriP1[2],TriP2[2]),TriP3[2]);
				temp_MaxZ = max(max(TriP1[2],TriP2[2]),TriP3[2]);
				
				if(temp_MinZ < Total_Min_Z)
					Total_Min_Z = temp_MinZ;
				if(temp_MaxZ > Total_Max_Z)
					Total_Max_Z = temp_MaxZ;

				
				double diff_H = fabs(Total_Max_Z - Total_Min_Z)/3.0;
				double BF;

                if(pyramid_step >= 3)
                    BF = BufferOfHeight*angle_weight;
				else
                    BF = BufferOfHeight;
                
				// calculation on BoundingBox(MinMax XY) of triangle
				TriMinXY[0]	= min(min(TriP1[0],TriP2[0]),TriP3[0]);
				TriMinXY[1]	= min(min(TriP1[1],TriP2[1]),TriP3[1]);
				TriMaxXY[0]	= max(max(TriP1[0],TriP2[0]),TriP3[0]);
				TriMaxXY[1]	= max(max(TriP1[1],TriP2[1]),TriP3[1]);
				
				PixelMinXY[0] = (int)((TriMinXY[0] - boundary[0])/gridspace + 0.5);
				PixelMinXY[1] = (int)((TriMinXY[1] - boundary[1])/gridspace + 0.5);
				PixelMaxXY[0] = (int)((TriMaxXY[0] - boundary[0])/gridspace + 0.5);
				PixelMaxXY[1] = (int)((TriMaxXY[1] - boundary[1])/gridspace + 0.5);
				
				PixelMinXY[0] -= 1;		PixelMinXY[1] -= 1;
				PixelMaxXY[0] += 1;		PixelMaxXY[1] += 1;
				if (PixelMaxXY[0] >= (int)(TIN_Grid_Size_X))	
					PixelMaxXY[0] =	 (int)(TIN_Grid_Size_X-1);
				if (PixelMaxXY[1] >= (int)(TIN_Grid_Size_Y))   
					PixelMaxXY[1] =	 (int)(TIN_Grid_Size_Y-1);
				if (PixelMinXY[0] < 0)	
					PixelMinXY[0] = 0;
				if (PixelMinXY[1] < 0)	
					PixelMinXY[1] = 0;
				
				for (Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
				{
					for (Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
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
							
							CurGPXY[0]	= (Col)*gridspace + boundary[0];
							CurGPXY[1]	= (Row)*gridspace + boundary[1];
							
							_p1[0]		= TriP1[0];
							_p1[1]		= TriP1[1];
							_p1[2]		= TriP1[2];
							
							_p2[0]		= TriP2[0];
							_p2[1]		= TriP2[1];
							_p2[2]		= TriP2[2];
							
							_p3[0]		= TriP3[0];
							_p3[1]		= TriP3[1];
							_p3[2]		= TriP3[2];
							
							v12[0]		= _p2[0]-_p1[0];
							v12[1]		= _p2[1]-_p1[1];
							
							v1P[0]		= CurGPXY[0]-_p1[0];
							v1P[1]		= CurGPXY[1]-_p1[1];
							
							v23[0]		= _p3[0]-_p2[0];
							v23[1]		= _p3[1]-_p2[1];
							
							v2P[0]		= CurGPXY[0]-_p2[0];
							v2P[1]		= CurGPXY[1]-_p2[1];
							
							v31[0]		= _p1[0]-_p3[0];
							v31[1]		= _p1[1]-_p3[1];
							
							v3P[0]		= CurGPXY[0]-_p3[0];
							v3P[1]		= CurGPXY[1]-_p3[1];
							
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
								
								v12[0]	= _p2[0]-_p1[0];
								v12[1]	= _p2[1]-_p1[1];
								v12[2]	= _p2[2]-_p1[2];
								
								v13[0]	= _p3[0]-_p1[0];
								v13[1]	= _p3[1]-_p1[1];
								v13[2]	= _p3[2]-_p1[2];
								
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
								GridPT3[Index].angle = (double)angle;
							}
							
							if (rtn)
							{
								double diff1, diff2, diff3, t1, t2;	   
								double pre_minHeight, pre_maxHeight;
								
								GridPT3[Index].t_minHeight	= GridPT3[Index].minHeight;
								GridPT3[Index].t_maxHeight	= GridPT3[Index].maxHeight;
								
								CurGPXY[0]	= (Col)*gridspace + boundary[0];
								CurGPXY[1]	= (Row)*gridspace + boundary[1];
								
								//// IDW
								diff1 = sqrt((CurGPXY[0] - TriP1[0])*(CurGPXY[0] - TriP1[0]) + (CurGPXY[1] - TriP1[1])*(CurGPXY[1] - TriP1[1]));
								diff2 = sqrt((CurGPXY[0] - TriP2[0])*(CurGPXY[0] - TriP2[0]) + (CurGPXY[1] - TriP2[1])*(CurGPXY[1] - TriP2[1]));
								diff3 = sqrt((CurGPXY[0] - TriP3[0])*(CurGPXY[0] - TriP3[0]) + (CurGPXY[1] - TriP3[1])*(CurGPXY[1] - TriP3[1]));
								
								if(diff1 == 0)
								{
									Z	= TriP1[2];
									GridPT3[Index].Matched_flag = 2;
								}
								else if(diff2 == 0)
								{
									Z	= TriP2[2];
									GridPT3[Index].Matched_flag = 2;
								}
								else if(diff3 == 0)
								{
									Z	= TriP3[2];
									GridPT3[Index].Matched_flag = 2;
								}
								
								m_bHeight[Index] = 1;
								GridPT3[Index].Height = Z; 
								
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
								t1		 = min(temp_MinZ, Z);
								if(GridPT3[Index].Matched_flag == 4) //extension minHeight
								{
									if(t1 - BF <= GridPT3[Index].minHeight)
										GridPT3[Index].minHeight   = t1 - BF;
								}
								else
								{
									//reduce minHeight
									GridPT3[Index].minHeight   = t1 - BF;
								}
								
								t2		 = max(temp_MaxZ, Z);
								if(GridPT3[Index].Matched_flag == 4)
								{
									if(t2 + BF >= GridPT3[Index].maxHeight)
										GridPT3[Index].maxHeight   = t2 + BF;
								}
								else
								{
									GridPT3[Index].maxHeight   = t2 + BF;
								}
								
								if(GridPT3[Index].minHeight > GridPT3[Index].maxHeight)
								{
									double temp_H = GridPT3[Index].minHeight;
									GridPT3[Index].minHeight = GridPT3[Index].maxHeight;
									GridPT3[Index].maxHeight = temp_H;
								}
								
								if(GridPT3[Index].Matched_flag == 2)
								{
									double mintemp, maxtemp;
									mintemp = Z - BF;
									maxtemp = Z + BF;
									
									GridPT3[Index].minHeight = mintemp;
									GridPT3[Index].maxHeight = maxtemp;
								}
								else
								{
									GridPT3[Index].Matched_flag = 1;
									if(GridPT3[Index].ortho_ncc > ortho_ncc_th)
									{
										
										GridPT3[Index].minHeight = Z - BF;
										GridPT3[Index].maxHeight = Z + BF;
									}
								}
                            }
						}
					}
				}
			}
		}
	}

	*minH_grid	= 100000.0;
	*maxH_grid	= -100000.0;
	
	double minH_temp = *minH_grid;
	double maxH_temp = *maxH_grid;

	UGRID *GridPT3_temp = (UGRID*)malloc(sizeof(UGRID)*TIN_Grid_Size_X*TIN_Grid_Size_Y);
	for (int counter = 0; counter < TIN_Grid_Size_X*TIN_Grid_Size_Y; counter++)
	{
	  GridPT3_temp[counter] = GridPT3[counter];
	}
	
	// minmaxheight setup for no matched grids
#pragma omp parallel for shared(TIN_Grid_Size_X,TIN_Grid_Size_Y,GridPT3,m_bHeight,Total_Min_Z,Total_Max_Z,minH_grid,maxH_grid) private(Col_C,Row_R) reduction(max:maxH_temp) reduction(min:minH_temp)
		for (Row_R=0; Row_R < (int)(TIN_Grid_Size_Y); Row_R++)
		{
			for (Col_C=0; Col_C < (int)(TIN_Grid_Size_X); Col_C++)
			{
				int Index = TIN_Grid_Size_X*Row_R + Col_C;
				
				if(GridPT3[Index].Matched_flag == 0)
				{
					double min_H = 10000.0;
					double max_H = -10000.0;
					bool  t_flag= false;
					int t_r, t_c;
					for(t_r = -1 ; t_r <= 1 ; t_r++)
					{
						for(t_c = -1 ; t_c <= 1 ; t_c++)
						{
							if(Col_C + t_c >= 0 && Col_C + t_c < (int)(TIN_Grid_Size_X) && Row_R + t_r >= 0 && Row_R + t_r < (int)(TIN_Grid_Size_Y))
							{
								int t_index = (Col_C + t_c) + (Row_R + t_r)*TIN_Grid_Size_X;
								if(GridPT3_temp[t_index].Matched_flag != 0)
								{
									t_flag = true;
									if(min_H > GridPT3_temp[t_index].minHeight)
										min_H	= GridPT3_temp[t_index].minHeight;
									if(max_H < GridPT3_temp[t_index].maxHeight)
										max_H	= GridPT3_temp[t_index].maxHeight;										 
								}
							}
						}
					}
					
					if(t_flag)
					{
						GridPT3[Index].minHeight = min_H;
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
							
							if(BF < pow(2.0,pyramid_step)*0.5)
								BF = pow(2.0,pyramid_step)*0.5;
							
							if(GridPT3[Index].minHeight < Total_Min_Z - BF)
								GridPT3[Index].minHeight   =  (Total_Min_Z - BF);
							if(GridPT3[Index].maxHeight > Total_Max_Z + BF)
								GridPT3[Index].maxHeight   =  (Total_Max_Z + BF);
						}
						
					}
				}
				else
				{
					if(m_bHeight[Index] == 0)
						GridPT3[Index].Matched_flag = 0;
				}
				
//#pragma omp critical
				{
					if(minH_temp > GridPT3[Index].minHeight)
						minH_temp	= GridPT3[Index].minHeight;
					if(maxH_temp < GridPT3[Index].maxHeight)
						maxH_temp	= GridPT3[Index].maxHeight;
				}
			}
		}

	*minH_grid = minH_temp;
	*maxH_grid = maxH_temp;

	free(GridPT3_temp);
	
	printf("end grid set height!!\n");
	
	int k, j;
	
	result					= (UGRID*)calloc(TIN_Grid_Size_Y*TIN_Grid_Size_X,sizeof(UGRID));
	
#pragma omp parallel for shared(TIN_Grid_Size_X,TIN_Grid_Size_Y,GridPT3,result,m_bHeight,Total_Min_Z,Total_Max_Z) private(k,j)
	for(k=0;k<(int)(TIN_Grid_Size_Y);k++)
	{
		for(j=0;j<(int)(TIN_Grid_Size_X);j++)
		{
			int matlab_index	= k*TIN_Grid_Size_X + j;
			
			if(GridPT3[matlab_index].minHeight > GridPT3[matlab_index].maxHeight)
			{
				double temp;
				temp = GridPT3[matlab_index].maxHeight;
				GridPT3[matlab_index].maxHeight = GridPT3[matlab_index].minHeight;
				GridPT3[matlab_index].minHeight = temp;
				
			}

			result[matlab_index].Height						= -1000;
			
			result[matlab_index].Matched_flag				= GridPT3[matlab_index].Matched_flag;
			if(m_bHeight[matlab_index] == 2 && GridPT3[matlab_index].Matched_flag == 0)
				result[matlab_index].Matched_flag			= 1;
			
			result[matlab_index].roh						= GridPT3[matlab_index].roh;
			
			result[matlab_index].Height						= GridPT3[matlab_index].Height;
			
			result[matlab_index].anchor_flag				= 0;

			result[matlab_index].ortho_ncc					= GridPT3[matlab_index].ortho_ncc;
			result[matlab_index].minHeight					= GridPT3[matlab_index].minHeight;
			
			if(result[matlab_index].minHeight < -100)
				result[matlab_index].minHeight			= -100;
			
			
			result[matlab_index].maxHeight					= GridPT3[matlab_index].maxHeight;
			
//			result[matlab_index].false_h_count	= 0;
			
			if(pyramid_step >= 2)
			{
				if(!check_level_end)
				{
					
					
/*					if(GridPT3[matlab_index].false_h_count > 0)
					{
						result[matlab_index].false_h_count				= GridPT3[matlab_index].false_h_count;
						
						result[matlab_index].false_h		  = (double*)calloc(sizeof(double),result[matlab_index].false_h_count);
						for(int kkk = 0 ; kkk<result[matlab_index].false_h_count ; kkk++)
							result[matlab_index].false_h[kkk] = GridPT3[matlab_index].false_h[kkk];
						free(GridPT3[matlab_index].false_h);
					}
*/				}
				else
				{
					
//					if(GridPT3[matlab_index].false_h_count > 0)
//						free(GridPT3[matlab_index].false_h);
				}
			}
			else
			{
//				if(GridPT3[matlab_index].false_h_count > 0)
//					free(GridPT3[matlab_index].false_h);
			}
			
			if(m_bHeight[matlab_index] == 0)
				result[matlab_index].Height = -1000.0;
			
		}
	}
	
	if(m_bHeight)
		free(m_bHeight);
	
	printf("end updating grid set height!!\n");
	
	free(GridPT3);
	
	printf("end memory release updating grid set height!!\n");

	return result;
}


UGRID* ResizeGirdPT3(CSize preSize, CSize resize_Size, double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight)
{
	
	UGRID *resize_GridPT3 = (UGRID *)calloc(sizeof(UGRID),resize_Size.height*resize_Size.width);
	
	for(int row=0;row<resize_Size.height;row++)
	{
		for(int col=0;col<resize_Size.width;col++)
		{
			long int index = row*resize_Size.width + col;
			double X = resize_Grid[index].m_X;
			double Y = resize_Grid[index].m_Y;
			int pos_c = (int)((X - Boundary[0])/pre_gridsize);
			int pos_r = (int)((Y - Boundary[1])/pre_gridsize);
			long int pre_index = pos_r*preSize.width + pos_c;
			if(pos_c >= 0 && pos_c < preSize.width && pos_r >= 0 && pos_r < preSize.height)
			{
				resize_GridPT3[index].minHeight		= preGridPT3[pre_index].minHeight;
				resize_GridPT3[index].maxHeight		= preGridPT3[pre_index].maxHeight;
				resize_GridPT3[index].Height		= preGridPT3[pre_index].Height;
				resize_GridPT3[index].Matched_flag	= preGridPT3[pre_index].Matched_flag;
				resize_GridPT3[index].roh			= preGridPT3[pre_index].roh;
				resize_GridPT3[index].anchor_flag	= preGridPT3[pre_index].anchor_flag;
				resize_GridPT3[index].ortho_ncc		= preGridPT3[pre_index].ortho_ncc;
				resize_GridPT3[index].angle			= preGridPT3[pre_index].angle;
//				resize_GridPT3[index].false_h_count = 0;
			}
			else
			{
				resize_GridPT3[index].minHeight		= (double)(minmaxheight[0] - 0.5);
				resize_GridPT3[index].maxHeight		= (double)(minmaxheight[1] + 0.5);
				resize_GridPT3[index].Height		= -1000;
				resize_GridPT3[index].Matched_flag	= 0;
				resize_GridPT3[index].roh			= 0.0;
				resize_GridPT3[index].anchor_flag	= 0;
				resize_GridPT3[index].ortho_ncc		= 0;
				resize_GridPT3[index].angle			= 0;
//				resize_GridPT3[index].false_h_count = 0;
			}
		}
	}
	
	printf("before release preGirdPT3\n");
	
	free(preGridPT3);
	
	printf("after release preGirdPT3\n");
	
	return resize_GridPT3;
}

UGRID* ResizeGirdPT3_RA(CSize preSize, CSize resize_Size, double* preBoundary,double* Boundary, D2DPOINT *resize_Grid, UGRID *preGridPT3, double pre_gridsize, double* minmaxheight)
{
    
    UGRID *resize_GridPT3 = (UGRID *)calloc(sizeof(UGRID),resize_Size.height*resize_Size.width);
    
    for(int row=0;row<resize_Size.height;row++)
    {
        for(int col=0;col<resize_Size.width;col++)
        {
            long int index = row*resize_Size.width + col;
            double X = resize_Grid[index].m_X;
            double Y = resize_Grid[index].m_Y;
            int pos_c = (int)((X - preBoundary[0])/pre_gridsize);
            int pos_r = (int)((Y - preBoundary[1])/pre_gridsize);
            long int pre_index = pos_r*preSize.width + pos_c;
            if(pos_c >= 0 && pos_c < preSize.width && pos_r >= 0 && pos_r < preSize.height)
            {
                resize_GridPT3[index].minHeight		= preGridPT3[pre_index].minHeight;
                resize_GridPT3[index].maxHeight		= preGridPT3[pre_index].maxHeight;
                resize_GridPT3[index].Height		= preGridPT3[pre_index].Height;
                resize_GridPT3[index].Matched_flag	= preGridPT3[pre_index].Matched_flag;
                resize_GridPT3[index].roh			= preGridPT3[pre_index].roh;
                resize_GridPT3[index].anchor_flag	= preGridPT3[pre_index].anchor_flag;
                resize_GridPT3[index].ortho_ncc		= preGridPT3[pre_index].ortho_ncc;
                resize_GridPT3[index].angle			= preGridPT3[pre_index].angle;
                //				resize_GridPT3[index].false_h_count = 0;
            }
            else
            {
                resize_GridPT3[index].minHeight		= (double)(minmaxheight[0] - 0.5);
                resize_GridPT3[index].maxHeight		= (double)(minmaxheight[1] + 0.5);
                resize_GridPT3[index].Height		= -1000;
                resize_GridPT3[index].Matched_flag	= 0;
                resize_GridPT3[index].roh			= 0.0;
                resize_GridPT3[index].anchor_flag	= 0;
                resize_GridPT3[index].ortho_ncc		= 0;
                resize_GridPT3[index].angle			= 0;
                //				resize_GridPT3[index].false_h_count = 0;
            }
        }
    }
    
    printf("before release preGirdPT3\n");
    
    free(preGridPT3);
    
    printf("after release preGirdPT3\n");
    
    return resize_GridPT3;
}

bool SetHeightRange_blunder(double* minmaxHeight,D3DPOINT *pts, int numOfPts, UI3DPOINT *tris,int numOfTri, UGRID *GridPT3, BL BL_param, double *mt_minmaxheight,bool blunder_update)
{
	uint32 num_triangles;
	int i, tcnt;
	double gridspace;
	double *boundary;
	CSize gridsize;
	uint8 pyramid_step, iteration;
	uint32 TIN_Grid_Size_X, TIN_Grid_Size_Y;
	uint8* m_bHeight;
	double Total_Min_Z		=  100000;
	double Total_Max_Z		= -100000;
	double meters_per_pixel= 2.0;
	double DEM_error		= meters_per_pixel*2;
	int Col_C, Row_R;
	
	num_triangles=numOfTri;
	pyramid_step	= BL_param.Pyramid_step;
	iteration		= BL_param.iteration;
	gridspace = BL_param.gridspace;
	boundary	= BL_param.Boundary;
	gridsize.width	= BL_param.Size_Grid2D.width;
	gridsize.height	= BL_param.Size_Grid2D.height;
	TIN_Grid_Size_X= gridsize.width;
	TIN_Grid_Size_Y= gridsize.height;
	
	m_bHeight		= (uint8*)calloc(TIN_Grid_Size_Y*TIN_Grid_Size_X,sizeof(uint8));
	
	//#pragma omp parallel shared(Total_Min_Z,Total_Max_Z,GridPT3,pts,tris,num_triangles,m_bHeight,pyramid_step,gridspace,boundary,TIN_Grid_Size_X,TIN_Grid_Size_Y,DEM_error) private(tcnt)
	{
		//#pragma omp for
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
			
			t_tri	= tris[tcnt];
			pdex0 = t_tri.m_X;
			pdex1 = t_tri.m_Y;
			pdex2 = t_tri.m_Z;
			
			if(pdex0 < numPts && pdex1 < numPts && pdex2 < numPts)
			{
				int node1_index, node2_index, node3_index;
				bool check_anchor_flag = false;
				pt0	= pts[pdex0];
				pt1	= pts[pdex1];
				pt2	= pts[pdex2];
				
				TriP1[0]		= pt0.m_X;
				TriP2[0]		= pt1.m_X;
				TriP3[0]		= pt2.m_X;
				
				TriP1[1]		= pt0.m_Y;
				TriP2[1]		= pt1.m_Y;
				TriP3[1]		= pt2.m_Y;
				
				TriP1[2]		= pt0.m_Z;
				TriP2[2]		= pt1.m_Z;
				TriP3[2]		= pt2.m_Z;
				
				temp_MinZ = min(min(TriP1[2],TriP2[2]),TriP3[2]);
				temp_MaxZ = max(max(TriP1[2],TriP2[2]),TriP3[2]);
				
				if(temp_MinZ < Total_Min_Z)
					Total_Min_Z = temp_MinZ;
				if(temp_MaxZ > Total_Max_Z)
					Total_Max_Z = temp_MaxZ;
				
				// calculation on BoundingBox(MinMax XY) of triangle
				TriMinXY[0]	= min(min(TriP1[0],TriP2[0]),TriP3[0]);
				TriMinXY[1]	= min(min(TriP1[1],TriP2[1]),TriP3[1]);
				TriMaxXY[0]	= max(max(TriP1[0],TriP2[0]),TriP3[0]);
				TriMaxXY[1]	= max(max(TriP1[1],TriP2[1]),TriP3[1]);
				
				PixelMinXY[0] = (int)((TriMinXY[0] - boundary[0])/gridspace + 0.5);
				PixelMinXY[1] = (int)((TriMinXY[1] - boundary[1])/gridspace + 0.5);
				PixelMaxXY[0] = (int)((TriMaxXY[0] - boundary[0])/gridspace + 0.5);
				PixelMaxXY[1] = (int)((TriMaxXY[1] - boundary[1])/gridspace + 0.5);
				
				PixelMinXY[0] -= 1;		PixelMinXY[1] -= 1;
				PixelMaxXY[0] += 1;		PixelMaxXY[1] += 1;
				if (PixelMaxXY[0] >= (int)(TIN_Grid_Size_X))	
					PixelMaxXY[0] =	 (int)(TIN_Grid_Size_X-1);
				if (PixelMaxXY[1] >= (int)(TIN_Grid_Size_Y))   
					PixelMaxXY[1] =	 (int)(TIN_Grid_Size_Y-1);
				if (PixelMinXY[0] < 0)	
					PixelMinXY[0] = 0;
				if (PixelMinXY[1] < 0)	
					PixelMinXY[1] = 0;
				
				for (Row=PixelMinXY[1]; Row <= PixelMaxXY[1]; Row++)
				{
					for (Col=PixelMinXY[0]; Col <= PixelMaxXY[0]; Col++)
					{
						int Index= TIN_Grid_Size_X*Row + Col;
						
						if (m_bHeight[Index] == 0)	
						{
							double CurGPXY[2]={0.};
							double Z = -1000.0;
							bool rtn = false;
							double _p1[3], _p2[3], _p3[3], v12[2], v1P[2];
							double v23[2], v2P[2], v31[2], v3P[2];
							int Sum;
							int t_col_count = 0;
							int t_row_count = 0;
							
							CurGPXY[0]	= (Col)*gridspace + boundary[0];
							CurGPXY[1]	= (Row)*gridspace + boundary[1];
							
							_p1[0]		= TriP1[0];
							_p1[1]		= TriP1[1];
							_p1[2]		= TriP1[2];
							
							_p2[0]		= TriP2[0];
							_p2[1]		= TriP2[1];
							_p2[2]		= TriP2[2];
							
							_p3[0]		= TriP3[0];
							_p3[1]		= TriP3[1];
							_p3[2]		= TriP3[2];
							
							v12[0]		= _p2[0]-_p1[0];
							v12[1]		= _p2[1]-_p1[1];
							
							v1P[0]		= CurGPXY[0]-_p1[0];
							v1P[1]		= CurGPXY[1]-_p1[1];
							
							v23[0]		= _p3[0]-_p2[0];
							v23[1]		= _p3[1]-_p2[1];
							
							v2P[0]		= CurGPXY[0]-_p2[0];
							v2P[1]		= CurGPXY[1]-_p2[1];
							
							v31[0]		= _p1[0]-_p3[0];
							v31[1]		= _p1[1]-_p3[1];
							
							v3P[0]		= CurGPXY[0]-_p3[0];
							v3P[1]		= CurGPXY[1]-_p3[1];
							
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
								
								v12[0]	= _p2[0]-_p1[0];
								v12[1]	= _p2[1]-_p1[1];
								v12[2]	= _p2[2]-_p1[2];
								
								v13[0]	= _p3[0]-_p1[0];
								v13[1]	= _p3[1]-_p1[1];
								v13[2]	= _p3[2]-_p1[2];
								
								Normal[0]=v12[1]*v13[2] - v12[2]*v13[1];
								Normal[1]=v12[2]*v13[0] - v12[0]*v13[2];
								Normal[2]=v12[0]*v13[1] - v12[1]*v13[0];	  
								
								Len=sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
								if(Len > 0)
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
								double t1,t2;
								m_bHeight[Index] = 1;
								GridPT3[Index].Height = Z;
								
								GridPT3[Index].t_minHeight = temp_MinZ;
								GridPT3[Index].t_maxHeight = temp_MaxZ;
							}
						}
					}
				}
			}
		}
	}
	if(m_bHeight)
		free(m_bHeight);
	
	mt_minmaxheight[0] = Total_Min_Z;
	mt_minmaxheight[1] = Total_Max_Z;
	
	
	return true;
}

void echoprint_Gridinfo(char *save_path,int row,int col,int level, int iteration, double update_flag, CSize *Size_Grid2D, UGRID *GridPT3, char *add_str)
{
	FILE *outfile_h,*outfile_min, *outfile_max,	 *outfile_roh, *outfile_flag;
	CSize temp_S;
	char t_str[500];
	int k,j;

	//sprintf(t_str,"%s/txt/tin_min_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	//outfile_min	= fopen(t_str,"w");
	//sprintf(t_str,"%s/txt/tin_max_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	//outfile_max	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/tin_h_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_h	= fopen(t_str,"w");
	//sprintf(t_str,"%s/txt/tin_roh_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	//outfile_roh	= fopen(t_str,"w");
	/*sprintf(t_str,"%s/txt/tin_flag_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	  outfile_flag	= fopen(t_str,"w");*/
	
//	if(outfile_min && outfile_max && outfile_h && outfile_roh && outfile_flag)
	{
		if(update_flag)
		{
			temp_S.height	= Size_Grid2D->height*2;
			temp_S.width	= Size_Grid2D->width*2;
		}
		else
		{
			temp_S.height	= Size_Grid2D->height;
			temp_S.width	= Size_Grid2D->width;
		}
			
		for(k=0;k<temp_S.height;k++)
		{
			for(j=0;j<temp_S.width;j++)
			{
				int matlab_index	= k*temp_S.width + j;

				//fprintf(outfile_min,"%f\t",GridPT3[matlab_index].minHeight);
				//fprintf(outfile_max,"%f\t",GridPT3[matlab_index].maxHeight);
				//if(GridPT3[matlab_index].Matched_flag != 0)
				fprintf(outfile_h,"%f\t",GridPT3[matlab_index].Height);
				//fprintf(outfile_roh,"%f\t",GridPT3[matlab_index].roh);
				/*fprintf(outfile_flag,"%d\t",GridPT3[matlab_index].Matched_flag);*/
			}
			//fprintf(outfile_min,"\n");
			//fprintf(outfile_max,"\n");
			fprintf(outfile_h,"\n");
			//fprintf(outfile_roh,"\n");
			/*fprintf(outfile_flag,"\n");*/
		}

		//fclose(outfile_min);
		//fclose(outfile_max);
		fclose(outfile_h);
		//fclose(outfile_roh);
		/*fclose(outfile_flag);*/
	}
}

void echo_print_nccresults(char *save_path,int row,int col,int level, int iteration, NCCresult *nccresult, CSize *Size_Grid2D, char *add_str)
{
	int k,j;
	FILE *outfile_min, *outfile_max, *outfile_h, *outfile_roh, *outfile_diff, *outfile_peak, *outINCC, *outGNCC,*outcount;
	CSize temp_S;
	char t_str[500];

	sprintf(t_str,"%s/txt/nccresult_roh1_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_min	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/nccresult_roh2_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_max	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/nccresult_h1_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_h	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/nccresult_h2_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_roh	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/nccresult_peaks_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_peak	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/nccresult_diff_level_%d_%d_%d_iter_%d_%s.txt",save_path,row,col,level,iteration,add_str);
	outfile_diff	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/INCC_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
	outINCC	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/GNCC_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
	outGNCC	= fopen(t_str,"w");
	sprintf(t_str,"%s/txt/rohcount_level_%d_%d_%d_iter_%d.txt",save_path,row,col,level,iteration);
	outcount	= fopen(t_str,"w");
	
	temp_S.height	= Size_Grid2D->height;
	temp_S.width	= Size_Grid2D->width;
		
	for(k=0;k<temp_S.height;k++)
	{
		for(j=0;j<temp_S.width;j++)
		{
			int matlab_index	= k*temp_S.width + j;

			fprintf(outfile_min,"%f\t",nccresult[matlab_index].result0);
			fprintf(outfile_max,"%f\t",nccresult[matlab_index].result1);
			fprintf(outfile_h,"%f\t",nccresult[matlab_index].result2);
			fprintf(outfile_roh,"%f\t",nccresult[matlab_index].result3);
			fprintf(outfile_peak,"%f\t",nccresult[matlab_index].result4);
			
			
			fprintf(outINCC,"%f\t",nccresult[matlab_index].INCC);
			fprintf(outGNCC,"%f\t",nccresult[matlab_index].GNCC);
			fprintf(outcount,"%d\t",nccresult[matlab_index].roh_count);
		}
		fprintf(outfile_min,"\n");
		fprintf(outfile_max,"\n");
		fprintf(outfile_h,"\n");
		fprintf(outfile_roh,"\n");
		fprintf(outfile_peak,"\n");
		//fprintf(outfile_diff,"\n");
		
		fprintf(outINCC,"\n");
		fprintf(outGNCC,"\n");
		fprintf(outcount,"\n");
	}

	fclose(outfile_min);
	fclose(outfile_max);
	fclose(outfile_h);
	fclose(outfile_roh);
	fclose(outfile_peak);
	//fclose(outfile_diff);
	
	fclose(outINCC);
	fclose(outGNCC);
	fclose(outcount);
}

int AdjustParam(uint8 Pyramid_step, int NumofPts, char * file_pts, D2DPOINT Lstartpos, D2DPOINT Rstartpos, double **LRPCs, double **RRPCs, double *ImageAdjust, NCCflag _flag,
				uint8 Template_size, uint16 *LeftImage, CSize LImagesize, uint16 *RightImage, CSize RImagesize, uint8 *left_ori, uint8 *right_ori, TransParam param,
				double bin_angle, uint8 total_pyramid, bool Hemisphere, char* save_filepath, char* tmpdir)
{
	int i,iter_count;
	bool check_stop = false;
	
	double  left_IA[2] = {0.0};
	double **subA;
	double **TsubA;
	double **InverseSubA;

	FILE *fid_pts;
	D3DPOINT* MPs;
	D3DPOINT *Coord;
	
	int ii,kk;

	subA	= (double**)malloc(9*sizeof(double*));
	TsubA	= (double**)malloc(6*sizeof(double*));
	InverseSubA	= (double**)malloc(6*sizeof(double*));

	for(ii=0;ii<9;ii++)
	{
		subA[ii]	= (double*)malloc(6*sizeof(double));
		if(ii < 6)
		{
			TsubA[ii]		= (double*)malloc(9*sizeof(double));
			InverseSubA[ii]	= (double*)malloc(6*sizeof(double));
		}
	}

	for(ii=0;ii<9;ii++)
		subA[ii][0]	  = 1.0;

	subA[0][1] = -1.0; subA[0][2] = -1.0; subA[0][3] =	1.0; subA[0][4] =  1.0; subA[0][5] =  1.0;
	subA[1][1] =  0.0; subA[1][2] = -1.0; subA[1][3] =	0.0; subA[1][4] =  0.0; subA[1][5] =  1.0;
	subA[2][1] =  1.0; subA[2][2] = -1.0; subA[2][3] =	1.0; subA[2][4] = -1.0; subA[2][5] =  1.0;
	subA[3][1] = -1.0; subA[3][2] =	 0.0; subA[3][3] =	1.0; subA[3][4] =  0.0; subA[3][5] =  0.0;
	subA[4][1] =  0.0; subA[4][2] =	 0.0; subA[4][3] =	0.0; subA[4][4] =  0.0; subA[4][5] =  0.0;
	subA[5][1] =  1.0; subA[5][2] =	 0.0; subA[5][3] =	1.0; subA[5][4] =  0.0; subA[5][5] =  0.0;
	subA[6][1] = -1.0; subA[6][2] =	 1.0; subA[6][3] =	1.0; subA[6][4] = -1.0; subA[6][5] =  1.0;
	subA[7][1] =  0.0; subA[7][2] =	 1.0; subA[7][3] =	0.0; subA[7][4] =  0.0; subA[7][5] =  1.0;
	subA[8][1] =  1.0; subA[8][2] =	 1.0; subA[8][3] =	1.0; subA[8][4] =  1.0; subA[8][5] =  1.0;

	for(ii=0;ii<6;ii++)
		for(kk=0;kk<9;kk++)
			TsubA[ii][kk]		= subA[kk][ii];

	InverseSubA[0][0] =	 0.555556; InverseSubA[0][1] =	0.000000; InverseSubA[0][2] =  0.000000; InverseSubA[0][3] = -0.333333; InverseSubA[0][4] =	 0.000000; InverseSubA[0][5] = -0.333333;
	InverseSubA[1][0] =	 0.000000; InverseSubA[1][1] =	0.166667; InverseSubA[1][2] =  0.000000; InverseSubA[1][3] =  0.000000; InverseSubA[1][4] =	 0.000000; InverseSubA[1][5] =	0.000000;
	InverseSubA[2][0] =	 0.000000; InverseSubA[2][1] =	0.000000; InverseSubA[2][2] =  0.166667; InverseSubA[2][3] =  0.000000; InverseSubA[2][4] =	 0.000000; InverseSubA[2][5] =	0.000000;
	InverseSubA[3][0] = -0.333333; InverseSubA[3][1] =	0.000000; InverseSubA[3][2] =  0.000000; InverseSubA[3][3] =  0.500000; InverseSubA[3][4] =	 0.000000; InverseSubA[3][5] =	0.000000;
	InverseSubA[4][0] =	 0.000000; InverseSubA[4][1] =	0.000000; InverseSubA[4][2] =  0.000000; InverseSubA[4][3] =  0.000000; InverseSubA[4][4] =	 0.250000; InverseSubA[4][5] =	0.000000;
	InverseSubA[5][0] = -0.333333; InverseSubA[5][1] =	0.000000; InverseSubA[5][2] =  0.000000; InverseSubA[5][3] =  0.000000; InverseSubA[5][4] =	 0.000000; InverseSubA[5][5] =	0.500000;

	fid_pts		= fopen(file_pts,"r");
	MPs = (D3DPOINT*)malloc(sizeof(D3DPOINT)*NumofPts);
    
	for(i=0;i<NumofPts;i++)
    {
		fscanf(fid_pts,"%lf %lf %lf %hhd\n",&MPs[i].m_X,&MPs[i].m_Y,&MPs[i].m_Z,&MPs[i].flag);
        //printf("t_coord %f\t%f\t%f\n",MPs[i].m_X,MPs[i].m_Y,MPs[i].m_Z);
    }
 
	fclose(fid_pts);

	iter_count = 1;

	Coord			= ps2wgs_3D(param,NumofPts,MPs);

	while(!check_stop && iter_count < 10)
	{
		uint8	Half_template_size;
		double b_factor;
		bool flag_boundary = false;
		int i;
		int count_pts = 0;
		double sum_weight_X		= 0;
		double sum_weight_Y		= 0;
		double sum_max_roh		= 0;
		double t_sum_weight_X		= 0;
		double t_sum_weight_Y		= 0;
		double t_sum_max_roh		= 0;
		double shift_X, shift_Y;

		//calculation image coord from object coord by RFM in left and right image
		b_factor			 = pow(2.0,(total_pyramid-Pyramid_step))*2;
		Half_template_size	 = (int)(Template_size/2.0);

		#pragma omp parallel for shared(Pyramid_step,Lstartpos,Rstartpos,ImageAdjust,NumofPts,LRPCs,RRPCs,left_IA,Coord,Half_template_size,b_factor,LImagesize,RImagesize,left_ori,right_ori,subA,TsubA,InverseSubA) private(i,t_sum_weight_X,t_sum_weight_Y,t_sum_max_roh) reduction(+:count_pts,sum_weight_X,sum_weight_Y,sum_max_roh)
		for(i = 0; i<NumofPts ; i++)
		{
			D2DPOINT Left_Imagecoord,Left_Imagecoord_p;
			D2DPOINT Right_Imagecoord, Right_Imagecoord_p;

			Left_Imagecoord_p	= GetObjectToImageRPC_single(LRPCs,2,left_IA,Coord[i]);
			Left_Imagecoord		= OriginalToPyramid_single(Left_Imagecoord_p,Lstartpos,Pyramid_step);

			Right_Imagecoord_p	= GetObjectToImageRPC_single(RRPCs,2,ImageAdjust,Coord[i]);
			Right_Imagecoord	= OriginalToPyramid_single(Right_Imagecoord_p,Rstartpos,Pyramid_step);

			if(	  Left_Imagecoord.m_Y  > Half_template_size*b_factor	+ 10					&& Left_Imagecoord.m_X	> Half_template_size*b_factor + 10
				  && Left_Imagecoord.m_Y  < LImagesize.height - Half_template_size*b_factor - 10	&& Left_Imagecoord.m_X	< LImagesize.width - Half_template_size*b_factor - 10
				  && Right_Imagecoord.m_Y > Half_template_size*b_factor	+ 10					&& Right_Imagecoord.m_X > Half_template_size*b_factor + 10
				  && Right_Imagecoord.m_Y < RImagesize.height - Half_template_size*b_factor - 10	&& Right_Imagecoord.m_X < RImagesize.width - Half_template_size*b_factor - 10)
			{
				double Left_X = Left_Imagecoord.m_X;
				double Left_Y = Left_Imagecoord.m_Y;
				double Right_X = Right_Imagecoord.m_X;
				double Right_Y = Right_Imagecoord.m_Y;
				int index_l	= ((int)Left_Y)*LImagesize.width + (int)Left_X;
				int index_r	= ((int)Right_Y)*RImagesize.width + (int)Right_X;
				double ori_diff;
				if( (index_l > 0 && index_l < LImagesize.height*LImagesize.width) && (index_r > 0 && index_r < RImagesize.height*RImagesize.width) )
				{
					ori_diff = left_ori[index_l] - right_ori[index_r];
					
					if(postNCC(Pyramid_step, ori_diff, Left_Y,	Left_X, Right_Y, Right_X,
							   subA,TsubA,InverseSubA,Template_size,_flag,bin_angle,LImagesize,RImagesize,LeftImage,RightImage,&t_sum_weight_X,&t_sum_weight_Y,&t_sum_max_roh))
					{
						//#pragma omp critical
						{
							sum_weight_X += t_sum_weight_X;
							sum_weight_Y += t_sum_weight_Y;
							sum_max_roh	 += t_sum_max_roh;
						}
						//#pragma omp atomic
						count_pts++;
					}
				}
			}
		}
		

		if(count_pts > 10)
		{
			shift_X				= sum_weight_X/sum_max_roh*pow(2.0,Pyramid_step);
			shift_Y				= sum_weight_Y/sum_max_roh*pow(2.0,Pyramid_step);
			if(fabs(shift_Y) < 0.1 && fabs(shift_X) < 0.1)
				check_stop = true;
 
			shift_X				+= ImageAdjust[1];
			shift_Y				+= ImageAdjust[0];

			ImageAdjust[1]		= shift_X;
			ImageAdjust[0]		= shift_Y;
		}
		else
		{
			check_stop = true;

		}

		iter_count++;
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
	
	free(Coord);
	free(MPs);

	return iter_count;
}


bool postNCC(uint8 Pyramid_step, double Ori_diff, double Left_CR,  double Left_CC, double Right_CR, double Right_CC, double **subA,double **TsubA,double **InverseSubA, uint8 Template_size, 
			 NCCflag _flag, double bin_angle, CSize leftsize, CSize rightsize, uint16* _leftimage, uint16* _rightimage, double *sum_weight_X, double *sum_weight_Y, double *sum_max_roh)
{

	// input order: NumOfPt, LeftImage, RightImage, Template_size, center_row_left, center_col_left, center_row_right, center_col_right, ncc_weight
	bool check_pt = false;
	int count_max_roh		= 0;

	uint8 rotate_flag,multi_flag,multi_flag_sum,inter_flag,weight_flag;

	//Image info input
	uint16 L_rowsize   = leftsize.height;
	uint16 L_colsize   = leftsize.width;

	uint16 R_rowsize   = rightsize.height;
	uint16 R_colsize   = rightsize.width;

	double t_weight_X	= 0;
	double t_weight_Y	= 0;
	double t_max_roh	= 0;
	double diff_theta;
	int mask_row, mask_col;

	int Half_template_size,half_mask_size;
	int count = 0;

	Half_template_size	= (int)(Template_size/2);
	half_mask_size		= 1;
	
	rotate_flag = _flag.rotate_flag;
	multi_flag	= _flag.multi_flag;
	multi_flag_sum	= _flag.multi_flag_sum;
	inter_flag	= _flag.inter_flag;
	weight_flag = _flag.weight_flag;

	
	diff_theta = Ori_diff;

	double *result_rho	= (double*)calloc(9,sizeof(double));
	double *XX			= (double*)calloc(6,sizeof(double));
	double *ATLT		= (double*)calloc(6,sizeof(double));
	int i, j, k;
	uint8 cell_count = 0;

	for(j=0;j<9;j++)
		result_rho[j]		= -1.00;

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
						double pos_row_left		 = (Left_CR + row);
						double pos_col_left		 = (Left_CC + col);

						double temp_col		   = (cos(-rot_theta)*col - sin(-rot_theta)*row);
						double temp_row		   = (sin(-rot_theta)*col + cos(-rot_theta)*row);
						double pos_row_right	 = (Right_CR + temp_row + mask_row);
						double pos_col_right	 = (Right_CC + temp_col + mask_col);

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
							
							//end
							Count_N[0]++;

							Sum_LR			  = Sum_LR + left_patch*right_patch;
							Sum_L			  = Sum_L  + left_patch;
							Sum_R			  = Sum_R  + right_patch;
							Sum_L2			  = Sum_L2 + left_patch*left_patch;
							Sum_R2			  = Sum_R2 + right_patch*right_patch;

							if(multi_flag == 1)
							{
								int size_1, size_2;
								size_1		  = (int)(Half_template_size/2);
								if( row >= -Half_template_size + size_1 && row <= Half_template_size - size_1)
								{
									if( col >= -Half_template_size + size_1 && col <= Half_template_size - size_1)
									{
										Sum_LR_2  = Sum_LR_2 + left_patch*right_patch;
										Sum_L_2	  = Sum_L_2	 + left_patch;
										Sum_R_2	  = Sum_R_2	 + right_patch;
										Sum_L2_2  = Sum_L2_2 + left_patch*left_patch;
										Sum_R2_2  = Sum_R2_2 + right_patch*right_patch;
										Count_N[1]++;
									}
								}

								size_2		  = size_1 + (int)((size_1/2.0) + 0.5);
								if( row >= -Half_template_size + size_2 && row <= Half_template_size - size_2)
								{
									if( col >= -Half_template_size + size_2 && col <= Half_template_size - size_2)
									{
										Sum_LR_3  = Sum_LR_3 + left_patch*right_patch;
										Sum_L_3	  = Sum_L_3	 + left_patch;
										Sum_R_3	  = Sum_R_3	 + right_patch;
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
				N				= Count_N[0];
				val1		  = (double)(Sum_L2) - (double)(Sum_L*Sum_L)/N;
				val2		  = (double)(Sum_R2) - (double)(Sum_R*Sum_R)/N;
				de			  = sqrt(val1*val2);
				de2			  = (double)(Sum_LR) - (double)(Sum_L*Sum_R)/N;
				if( val1*val2 > 0)
					ncc_1			= de2/de;
				else
					ncc_1			= -1.0;

				if(multi_flag == 1)
				{
					if(Count_N[1] > 0)
					{
						N					= Count_N[1];
						val1				= (double)(Sum_L2_2) - (double)(Sum_L_2*Sum_L_2)/N;
						val2				= (double)(Sum_R2_2) - (double)(Sum_R_2*Sum_R_2)/N;
						de					= sqrt(val1*val2);
						de2					= (double)(Sum_LR_2) - (double)(Sum_L_2*Sum_R_2)/N;
						if( val1*val2 > 0)
							ncc_2		  = de2/de;
						else
							ncc_2			= -1.0;
					}

					if(Count_N[2] > 0)
					{
						N					= Count_N[2];
						val1				= (double)(Sum_L2_3) - (double)(Sum_L_3*Sum_L_3)/N;
						val2				= (double)(Sum_R2_3) - (double)(Sum_R_3*Sum_R_3)/N;
						de					= sqrt(val1*val2);
						de2					= (double)(Sum_LR_3) - (double)(Sum_L_3*Sum_R_3)/N;
						if( val1*val2 > 0)
							ncc_3		  = de2/de;
						else
							ncc_3			= -1.0;
					}

				}

				if(multi_flag == 1)
				{
					if(Count_N[1] > 0 && Count_N[2] > 0)
						temp_rho	  = ((ncc_1 + ncc_2 + ncc_3)/3.0);
					else if(Count_N[1] > 0)
						temp_rho	  = ((ncc_1 + ncc_2)/2.0);
					else if(Count_N[2] > 0)
						temp_rho	  = ((ncc_1 + ncc_3)/2.0);
					else
						temp_rho		= ncc_1;
				}
				else
				{
					temp_rho	  = ncc_1;
				}
				

				grid_index			 = (mask_row+1)*3 + (mask_col+1);
				if(grid_index < 9)
					result_rho[grid_index] = temp_rho;
				cell_count++;
			}

			
		}
	}

	if(cell_count == 9)
	{
		double demnum;
		double max_X		= 100;
		double max_Y		= 100;
		double max_roh		= 0;
		bool find_index_1	= false;
		bool find_index_2	= false;
		bool find_index		= false;

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
				XX[i*1 + j]	= sum;
			}
		}

		demnum		= -pow(XX[4],2.0) + 4*XX[3]*XX[5];
		if(demnum > 0 && XX[3] < 0)
		{
			max_X = (- 2*XX[5]*XX[1] + XX[2]*XX[4])/demnum;
			max_Y = (- 2*XX[2]*XX[3] + XX[1]*XX[4])/demnum;
			max_roh =  XX[0]				+ XX[1]*max_X			+ XX[2]*max_Y
				+ XX[3]*max_X*max_X	+ XX[4]*max_X*max_Y	+ XX[5]*max_Y*max_Y;
			if(fabs(max_X) <= 1.0)
				find_index_1 = true;
			if(fabs(max_Y) <= 1.0)
				find_index_2 = true;
			if (Pyramid_step >= 2) 
				find_index	= find_index_1 & find_index_2 & (max_roh > 0.80);
			else
				find_index	= find_index_1 & find_index_2 & (max_roh > 0.90);

			if(find_index)
			{
				t_weight_X += max_X*max_roh;
				t_weight_Y += max_Y*max_roh;
				t_max_roh  += max_roh;

				check_pt = true;
			}
		}
	}
	free(result_rho);
	free(ATLT);
	free(XX);

	if(check_pt)
	{
		*sum_weight_X	= t_weight_X;
		*sum_weight_Y	= t_weight_Y;
		*sum_max_roh	= t_max_roh;
	}
	else
	{
		*sum_weight_X	= 0;
		*sum_weight_Y	= 0;
		*sum_max_roh	= 0;
	}
	
	return check_pt;
}

bool check_image_boundary(double **lrpc, double **rrpc, uint8 numofparam, double *rimageparam, D2DPOINT Lstartpos, D2DPOINT Rstartpos,
						  D2DPOINT pos_xy, double minH, double maxH, CSize Lsize, CSize Rsize, int H_template_size, int pyramid_step)
{
	bool check = true;

	bool bleft_s, bright_s, bleft_e, bright_e;
	D3DPOINT temp_gp;
	D2DPOINT temp;
	double limageparam[2] = {0.0};

	int buff_pixel = 1;
	//left start
	temp_gp.m_X = pos_xy.m_X;
	temp_gp.m_Y = pos_xy.m_Y;
	temp_gp.m_Z = minH;
	temp		= GetObjectToImageRPC_single(lrpc,numofparam,limageparam,temp_gp);
	temp		= OriginalToPyramid_single(temp,Lstartpos,pyramid_step);

	if(temp.m_X > H_template_size +buff_pixel && temp.m_X < Lsize.width - H_template_size -buff_pixel &&
	   temp.m_Y > H_template_size +buff_pixel && temp.m_Y < Lsize.height - H_template_size -buff_pixel)
		bleft_s		= true;
	else 
		bleft_s		= false;


	//left end
	temp_gp.m_Z = maxH;
	temp		= GetObjectToImageRPC_single(lrpc,numofparam,limageparam,temp_gp);
	temp		= OriginalToPyramid_single(temp,Lstartpos,pyramid_step);

	if(temp.m_X > H_template_size +buff_pixel && temp.m_X < Lsize.width - H_template_size -buff_pixel &&
	   temp.m_Y > H_template_size +buff_pixel && temp.m_Y < Lsize.height - H_template_size -buff_pixel)
		bleft_e		= true;
	else 
		bleft_e		= false;


	//right start
	temp_gp.m_Z = minH;
	temp		= GetObjectToImageRPC_single(rrpc,numofparam,rimageparam,temp_gp);
	temp		= OriginalToPyramid_single(temp,Rstartpos,pyramid_step);

	if(temp.m_X > H_template_size +buff_pixel && temp.m_X < Rsize.width - H_template_size -buff_pixel &&
	   temp.m_Y > H_template_size +buff_pixel && temp.m_Y < Rsize.height - H_template_size -buff_pixel)
		bright_s	 = true;
	else 
		bright_s	 = false;


	//left end
	temp_gp.m_Z = maxH;
	temp		= GetObjectToImageRPC_single(rrpc,numofparam,rimageparam,temp_gp);
	temp		= OriginalToPyramid_single(temp,Rstartpos,pyramid_step);

	if(temp.m_X > H_template_size +buff_pixel && temp.m_X < Rsize.width - H_template_size -buff_pixel &&
	   temp.m_Y > H_template_size +buff_pixel && temp.m_Y < Rsize.height - H_template_size -buff_pixel)
		bright_e	 = true;
	else 
		bright_e	 = false;

	if((bleft_s && bleft_e) && (bright_s && bright_e))
		check	= true;
	else
		check	= false;


	return check;

}
void RemoveFiles(char *save_path, char *lfilename, char *rfilename, int py_level, bool flag)
{
	int status;
	int count;

	char *filename_py;
	char t_str[500];
	int start_lv, end_lv;

	start_lv	= py_level;
	
	if(flag)
		end_lv		= py_level - 1;
	else
		end_lv		= -2;

	if(py_level == 0)
	{
		start_lv	= 0;
		end_lv		= -2;
	}

	for(count = 4 ; count >= 0 ; count--)
	{
		bool bfile = false;

		filename_py		= GetFileName(lfilename);
		filename_py		= remove_ext(filename_py);

		sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,count);

		status = remove(t_str);
		filename_py		= GetFileName(rfilename);
		filename_py		= remove_ext(filename_py);

		sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,count);
		status = remove(t_str);

		filename_py		= GetFileName(lfilename);
		filename_py		= remove_ext(filename_py);
		sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,count);
		status = remove(t_str);

		filename_py		= GetFileName(rfilename);
		filename_py		= remove_ext(filename_py);
		sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,count);
		status = remove(t_str);
		
		filename_py		= GetFileName(lfilename);
		filename_py		= remove_ext(filename_py);
		sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,count);
		status = remove(t_str);

		filename_py		= GetFileName(rfilename);
		filename_py		= remove_ext(filename_py);
		sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,count);
		status = remove(t_str);
 
		filename_py		= GetFileName(lfilename);
		filename_py		= remove_ext(filename_py);

		sprintf(t_str,"%s/%s.raw",save_path,filename_py);

		status = remove(t_str);
		
		filename_py		= GetFileName(rfilename);
		filename_py		= remove_ext(filename_py);

		sprintf(t_str,"%s/%s.raw",save_path,filename_py);
		status = remove(t_str);
	}
}

double MergeTiles(ProInfo info, int iter_row_start, int t_col_start, int iter_row_end,int t_col_end, int buffer,int final_iteration)
{
	FILE *poutDEM;
	FILE *poutMatchtag;
	FILE *poutheader;
	
    int header_line = 7;
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

	double *DEM, *DEMinter;
	//bool *Matchtag;
	char DEM_str[500];	  

	bool check_gs = false;

	//find boundary of DEM
	boundary[0] = 10000000.0;
	boundary[1] = 10000000.0;
	boundary[2] = -10000000.0;
	boundary[3] = -10000000.0;

#pragma omp parallel for private(index_file) schedule(guided)
	for(index_file = 0 ; index_file <= row_end*col_end ; index_file++)
	{
		int row,col;
		
		row = (int)(floor(index_file/col_end));
		col = index_file%col_end;
		
        if(row >= iter_row_start && row <= row_end && col >= t_col_start &&  col <= col_end)
        {
            FILE *pfile;
            char t_str[500];
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info.save_filepath,row,col,find_level,find_iter);
            pfile	= fopen(t_str,"r");
            if(pfile)
            {
                printf("matched tiles %s\n",t_str);
                fseek(pfile,0,SEEK_END);
                size = ftell(pfile);
                if(size > 0)
                {
                    char h_t_str[500];
                    FILE *p_hfile;
                    //sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col,find_level,find_iter);
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col);
                    p_hfile		= fopen(h_t_str,"r");
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
                                        boundary[0]		= t_boundary[0];
                                    if(boundary[1] > t_boundary[1])
                                        boundary[1]		= t_boundary[1];

                                    if(boundary[2] < t_boundary[2])
                                        boundary[2]		= t_boundary[2];
                                    if(boundary[3] < t_boundary[3])
                                        boundary[3]		= t_boundary[3];
                                }
                            }
                        }
                        fclose(p_hfile);
                    }
                }
                fclose(pfile);
            }
        }
	}
/*
	boundary[0] = (int)(boundary[0]/40.0)*40 - 40;
	boundary[1] = (int)(boundary[1]/40.0)*40 - 40;
	boundary[2] = (int)(boundary[2]/40.0)*40 + 40;
	boundary[3] = (int)(boundary[3]/40.0)*40 + 40;
*/
    boundary[0] -= 10*grid_size;//= (int)(boundary[0]/40.0)*40 - 40;
    boundary[1] -= 10*grid_size;//= (int)(boundary[1]/40.0)*40 - 40;
    boundary[2] += 10*grid_size;//= (int)(boundary[2]/40.0)*40 + 40;
    boundary[3] += 10*grid_size;//= (int)(boundary[3]/40.0)*40 + 40;
    
    printf("boundary %f\t%f\t%f\t%f\n",boundary[0],boundary[1],boundary[2],boundary[3]);
	//if(grid_size < info.DEM_resolution)
	check_gs = false;

	buffer	= (int)(floor(buffer/grid_size));

	DEM_size.width		= (int)(ceil( (double)(boundary[2] - boundary[0]) /grid_size ));
	DEM_size.height		= (int)(ceil( (double)(boundary[3] - boundary[1]) /grid_size ));
	DEM = (double*)malloc((long)DEM_size.height*(long)DEM_size.width*sizeof(double));
	
    printf("dem size %d\t%d\n",DEM_size.width,DEM_size.height);
    
#pragma omp parallel for schedule(guided)
	for(long index = 0 ; index < (long)DEM_size.height*(long)DEM_size.width ; index++)
	{
		DEM[index] = -9999;
	}

	if(check_gs)
	{
		DEMinter_size.width = (int)(ceil((boundary[2] - boundary[0])/info.DEM_resolution));
		DEMinter_size.height= (int)(ceil((boundary[3] - boundary[1])/info.DEM_resolution));
		DEMinter = (double*)malloc((long)DEMinter_size.height*(long)DEMinter_size.width*sizeof(double));

#pragma omp parallel for schedule(guided)
		for(long index = 0 ; index < (long)DEMinter_size.height*(long)DEMinter_size.width ; index++)
		{
			DEMinter[index] = -9999;
		}
	}

	//setting DEM value
#pragma omp parallel for private(index_file) schedule(guided)
	for(index_file = 0 ; index_file < row_end*col_end ; index_file++)
	{
		int row,col;
		FILE *pfile;
		char t_str[500];

		row = (int)(floor(index_file/col_end));
		col = index_file%col_end;
		
        if(row >= iter_row_start && row <= row_end && col >= t_col_start &&  col <= col_end)
        {
        
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_%d_%d.txt",info.save_filepath,row,col,find_level,find_iter);
            pfile	= fopen(t_str,"r");
            if(pfile)
            {
                fseek(pfile,0,SEEK_END);
                size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    char h_t_str[500];
                    FILE *p_hfile, *p_hvfile;
                    
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",info.save_filepath,row,col);
                    p_hfile		= fopen(h_t_str,"r");
                    if(p_hfile)
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
                        sprintf(hv_t_str,"%s/txt/tin_h_level_%d_%d_%d_iter_%d_final.txt",info.save_filepath,row,col,find_level,find_iter);
                        
                        p_hvfile	= fopen(hv_t_str,"r");
                        if(p_hvfile)
                        {
                            long index_total;
                            for(index_total = 0; index_total < row_size*col_size ; index_total++)
                            {
                                int iter_row,iter_col;
                                iter_row = floor(index_total/col_size);
                                iter_col = index_total%col_size;
                                {

                                    double t_col = ( (double)(t_boundary[0] + grid_size*iter_col - boundary[0])  /grid_size);
                                    double t_row = ( (double)(boundary[3] - (t_boundary[1] + grid_size*iter_row))/grid_size);
                                    long index = (long)(t_row*DEM_size.width + t_col + 0.01);
                                    
                                    double DEM_value;
                                    fscanf(p_hvfile,"%lf\t",&DEM_value);
                                    
                                    //if(t_row > 1400 && t_row < DEM_size.height - 1000 && t_col > 1200 && t_col < 1500)
                                    //    printf("buffer %d\t row %f\t col %f\t DEM %f\t %d\n",buffer,t_row,t_col,DEM_value,index);
                                    
                                    if(index >= 0 && index < (long)DEM_size.width*(long)DEM_size.height &&
                                       iter_row > buffer && iter_row < row_size - buffer &&
                                       iter_col > buffer && iter_col < col_size - buffer)
                                    {
                                        if(DEM_value > -1000 && DEM_value != 0)
                                            DEM[index] = DEM_value;
                                        //else if(t_row > 1400 && t_row < DEM_size.height - 1000 && t_col > 1200 && t_col < 1500)
                                        //    printf("buffer %d\t row %f\t col %f\t DEM %f\t %d\n",buffer,t_row,t_col,DEM[index],index);
                                        
                                        
                                    }
                                }
                                fscanf(p_hvfile,"\n");
                            }
                            
                            fclose(p_hvfile);
                        }
                        
                        fclose(p_hfile);
                    }
                }
                fclose(pfile);
            }
        }
	}
	
	if(check_gs)
	{
		long index;
		for(index = 0; index < (long)DEMinter_size.height*(long)DEMinter_size.width ; index++)
		{
			int row		= (int)(floor(index/DEMinter_size.width));
			int col		= index%DEMinter_size.width;

			double x		= boundary[0] + info.DEM_resolution*col;
			double y		= boundary[3] - info.DEM_resolution*row;

			double col_o = (x - boundary[0])/grid_size;
			double row_o = (boundary[3] - y)/grid_size;

			int	  row_o_i = (int)(row_o);
			int	  col_o_i = (int)(col_o);
			double dx	= row_o - (int)(row_o);
			double dy	= col_o - (int)(col_o);

			double s1	= dx*dy;
			double s2	= (1-dx)*dy;
			double s3	= dx*(1-dy);
			double s4	= (1-dx)*(1-dy);

			if(row_o >= 0 && row_o + 1 < DEM_size.height && col_o >= 0 && col_o +1 < DEM_size.width)
			{
				long index_1		=  (long)row_o_i	 *(long)DEM_size.width +	(long)col_o_i;
				long index_2		=  (long)row_o_i	 *(long)DEM_size.width + (long)(col_o_i + 1);
				long index_3		= (long)(row_o_i+1)*(long)DEM_size.width +	(long)col_o_i;
				long index_4		= (long)(row_o_i+1)*(long)DEM_size.width + (long)(col_o_i + 1);

				double value = s1*DEM[index_4] + s2*DEM[index_3] + s3*DEM[index_2] + s4*DEM[index_1];
				DEMinter[(long)row*(long)DEMinter_size.width + (long)col] = value;
			}
			else
			{
				double value = DEM[(long)(row_o)*(long)DEM_size.width + (long)(col_o)];
				DEMinter[(long)row*(long)DEMinter_size.width + (long)col] = value;
			}
		}

		sprintf(DEM_str, "%s/%s_dem_tin.txt", info.save_filepath, info.Outputpath_name);
		poutDEM	= fopen(DEM_str,"w");
		for (row = 0; row < DEMinter_size.height; row++) 
		{
			for (col = 0; col < DEMinter_size.width; col++) 
			{
				fprintf(poutDEM,"%f\t",DEMinter[(long)row*(long)DEMinter_size.width + (long)col]);
			}
			fprintf(poutDEM,"\r\n");
		}
		fclose(poutDEM);
		
		
		sprintf(DEM_str, "%s/%s_dem_header_tin.txt", info.save_filepath, info.Outputpath_name);
		
		printf("name %s\t%f\t%f\t%f\t%d\t%d\n",DEM_str,boundary[0],boundary[3],grid_size,DEM_size.width,DEM_size.height);
		
		poutheader = fopen(DEM_str,"w");
		fprintf(poutheader,"%f\t%f\t%f\t%d\t%d\n",boundary[0],boundary[3],info.DEM_resolution,DEMinter_size.width,DEMinter_size.height);
		fclose(poutheader);

		free(DEMinter);
		free(DEM);
		//free(Matchtag);
	}
	else
	{
		sprintf(DEM_str, "%s/%s_dem_tin.txt", info.save_filepath, info.Outputpath_name);
		poutDEM	= fopen(DEM_str,"w");
		for (row = 0; row < DEM_size.height; row++) 
		{
			for (col = 0; col < DEM_size.width; col++) 
			{
				fprintf(poutDEM,"%f\t",DEM[(long)row*(long)DEM_size.width + (long)col]);
			}
			fprintf(poutDEM,"\r\n");
		}
		fclose(poutDEM);
		
		
		sprintf(DEM_str, "%s/%s_dem_header_tin.txt", info.save_filepath, info.Outputpath_name);
		
		printf("name %s\t%f\t%f\t%f\t%d\t%d\n",DEM_str,boundary[0],boundary[3],grid_size,DEM_size.width,DEM_size.height);
		
		poutheader = fopen(DEM_str,"w");
		fprintf(poutheader,"%f\t%f\t%f\t%d\t%d\n",boundary[0],boundary[3],grid_size,DEM_size.width,DEM_size.height);
		fclose(poutheader);

		free(DEM);
		//free(Matchtag);
	}
	
	return grid_size;
}

void NNA_M(TransParam _param, char *save_path, char* Outputpath_name, char *iterfile, int row_start, int col_start,int row_end, int col_end, double grid_resolution, double mt_grid_resolution, int buffer_clip, int Hemisphere,int final_iteration)
{
	double dummy;
	int i0, cnthold,i,j,index;
	double minX,maxX, minY, maxY,DEM_minX,DEM_maxY,minHeight,maxHeight;
	double dis_X, dis_Y;
	double grid;
	double mt_grid;
	float *value;
	unsigned char	 *value_pt;
	int index_file;
	NNXY *pt_save;
	int row_count, col_count,row_count_mt, col_count_mt, DEM_rows,DEM_cols;
	FILE* fout;
	FILE* fheader;
	FILE	 *afile;
	FILE* t_file;
	
	char t_savefile[500];
	char outfile[500];
	char DEM_header[500];
	
	double X_center;
	double Y_center;
	
    int header_line = 7;
	int find_iter  = final_iteration;
	int interval_row = 20;
	int buffer_row	 = 5;
	int t_ndata;
	int ndim,ndim1,ndata;
	long total_search_count = 0;
	long total_mt_count = 0;
	long total_check_count = 0;
	long total_null_cell = 0;
	int threads_num;
	int guided_size;
	double *cal_gridpts_X = NULL;
	double *cal_gridpts_Y = NULL;
	struct calgrid
	{
		double X, Y;
		struct calgird *next;
	};
	
	struct calgrid *head = NULL;
	struct calgrid *cur = NULL;
	struct calgrid *tail = NULL;
	struct calgrid *newNode = NULL;
	
	ndim	 = 2;
	ndim1 = ndim + 1;
	
	minX = 10000000;
	minY = 10000000;
	maxX = -10000000;
	maxY = -10000000;
	
	minHeight = 100000;
	maxHeight = -100000;
	
	grid	 = grid_resolution;
	mt_grid	 = mt_grid_resolution;
	
	printf("DEM grid = %lf\t mt grid = %lf\n",grid,mt_grid);
	
	ndata = 0;
//#pragma omp parallel for shared(row_end,col_end,save_path) private(index_file) reduction(+:ndata)
	for(index_file = 0 ; index_file < row_end*col_end ; index_file++)
	{
		int row,col;
		FILE *pfile;
		char t_str[500];
		
		row = (int)(floor(index_file/col_end));
		col = index_file%col_end;
		
        if(row >= row_start && row <= row_end && col >= col_start &&  col <= col_end)
        {
            
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_0_%d.txt",save_path,row,col,find_iter);
            
            pfile	= fopen(t_str,"r");
            if(pfile)
            {
                //printf("row %d\tcol %d\n",row,col);
                long int size;
                fseek(pfile,0,SEEK_END);
                size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    double X,Y,Z;
                    while(fscanf(pfile,"%lf\t%lf\t%lf\n",&X,&Y,&Z) != EOF)
                    {
    //#pragma omp critical
                        {
                            if(minX > X)
                                minX	 = X;
                            if(minY > Y)
                                minY	 = Y;
                            
                            if(maxX < X)
                                maxX	 = X;
                            if(maxY < Y)
                                maxY	 = Y;
                            if(minHeight > Z)
                                minHeight = Z;
                            if(maxHeight < Z)
                                maxHeight = Z;
                        }
                        
                        ndata = ndata + 1;
                    }
                }
                fclose(pfile);
            }
        }
	}
	
    printf("end loading matched pts\n");
    
	/*minX = (int)(minX/40.0)*40 - 40;
	maxX = (int)(maxX/40.0)*40 + 40;
	minY = (int)(minY/40.0)*40 - 40;
	maxY = (int)(maxY/40.0)*40 + 40;
    */
    minX -= 10*grid;//= (int)(minX/40.0)*40 - 40;
    maxX += 10*grid;//= (int)(maxX/40.0)*40 + 40;
    minY -= 10*grid;//= (int)(minY/40.0)*40 - 40;
    maxY += 10*grid;//= (int)(maxY/40.0)*40 + 40;
	
	dis_X = maxX - minX;
	dis_Y = maxY - minY;
	
	col_count = (int)((maxX - minX)/grid) + 1;
	row_count = (int)((maxY - minY)/grid) + 1;
	
	col_count_mt = (int)((maxX - minX)/mt_grid) + 1;			
	row_count_mt = (int)((maxY - minY)/mt_grid) + 1;
	
	
	value = (float*)malloc(sizeof(float)*(long)row_count*(long)col_count);
	value_pt = (unsigned char*)malloc(sizeof(unsigned char)*(long)row_count*(long)col_count);
	pt_save	 = (NNXY*)malloc(sizeof(NNXY)*(long)row_count*(long)col_count);
	
	t_ndata = ndata;
	
	sprintf(DEM_header, "%s/%s_dem_header_tin.txt", save_path,Outputpath_name);
	fheader = fopen(DEM_header,"r");
	fscanf(fheader,"%lf\t%lf\t%lf\t%d\t%d",&DEM_minX,&DEM_maxY,&dummy,&DEM_cols,&DEM_rows);
	fclose(fheader);
	
	remove(DEM_header);
	
	
	for(i=0;i<row_count;i++)
	{
		for(j=0;j<col_count;j++)
			value[(long)i*(long)col_count + (long)j] = -9999;
	}
	
	if ((fheader = fopen(iterfile,"r")) != NULL)
	{
		sprintf(t_savefile,"%s/txt/tfile.txt",save_path);
		
		t_file = fopen(t_savefile,"w");
		
		for(long ix=0;ix<(long)DEM_cols*(long)DEM_rows;ix++)
		{
			int row,col;
			double t_z;
			double t_x, t_y;
			row = (int)(floor(ix/DEM_cols));
			col = ix%DEM_cols;
			
			t_x = DEM_minX + col*mt_grid;
			t_y = DEM_maxY - row*mt_grid;
			
			fscanf(fheader,"%lf",&t_z);
			
            double d_row,d_col;
			if(t_x >= minX && t_x < maxX && t_y >= minY && t_y < maxY)
			{
				d_row = ((maxY - t_y)/grid);
				d_col = ((t_x - minX)/grid);
                
                long t_index = (long)((long)d_row*(long)col_count + (long)d_col + 0.01);
                
                if(t_index < (long)col_count*(long)row_count && d_row < row_count && d_col < col_count)
                {
                    value[t_index] = t_z;
                    value_pt[t_index] = 0;
                    if(t_z > -1000)
                    {
                        fprintf(t_file,"%12.4lf\t%12.4lf\n",t_x,t_y);
                        total_search_count++;
                    }
                }
			}
			
		}
		fclose(t_file);
	}
	fclose(fheader);
	remove(iterfile);
	
	if(total_search_count > 0)
	{
		cal_gridpts_X	 = (double*)malloc(sizeof(double)*total_search_count);
		cal_gridpts_Y	 = (double*)malloc(sizeof(double)*total_search_count);
		t_file = fopen(t_savefile,"r");
		for(long ix=0;ix<total_search_count;ix++)
		{
			fscanf(t_file,"%lf\t%lf\n",&cal_gridpts_X[ix],&cal_gridpts_Y[ix]);
		}
		fclose(t_file);
		remove(t_savefile);
	}
	
	for(i=0;i<row_count;i++)
	{
		for(j=0;j<col_count;j++)
		{
			pt_save[(long)i*(long)col_count + (long)j].X = -999999999;
			pt_save[(long)i*(long)col_count + (long)j].Y = -999999999;
			pt_save[(long)i*(long)col_count + (long)j].Z = -9999;
		}
	}
	
	X_center = minX + dis_X/2.0;
	Y_center = minY + dis_Y/2.0;
	
	//printf("minXY %f\t%f\t%f\t%f\t %f\t%f\t %d\t%d\n",minX,maxX,minY,maxY,((maxX - minX)/grid) + 1,((maxY - minY)/grid) + 1,col_count_mt,row_count_mt);
    
	buffer_clip	 = (int)(floor(buffer_clip/mt_grid));
//#pragma omp parallel for shared(buffer_clip, row_end,col_end,save_path,pt_save,value,value_pt,X_center,Y_center,minX,maxY,grid,mt_grid,col_count,row_count,col_count_mt,row_count_mt,total_mt_count) private(index_file)
	for(index_file = 0 ; index_file < row_end*col_end ; index_file++)
	{
		int row,col;
		FILE *pfile;
		char t_str[500];
		
		row = (int)(floor(index_file/col_end));
		col = index_file%col_end;
		
        if(row >= row_start && row <= row_end && col >= col_start &&  col <= col_end)
        {
            
            sprintf(t_str,"%s/txt/matched_pts_%d_%d_0_%d.txt",save_path,row,col,find_iter);
            
            pfile	= fopen(t_str,"r");
            if(pfile)
            {
                //printf("matched row %d\tcol %d\n",row,col);
                
                long int size;
                fseek(pfile,0,SEEK_END);
                size = ftell(pfile);
                fseek(pfile,0L,SEEK_SET);
                if(size > 0)
                {
                    double pos_row, pos_col, clip_pos_row,clip_pos_col;
                    double t_x, t_y, t_z;
                    char h_t_str[500];
                    
                    FILE* p_hfile;
                    
                    sprintf(h_t_str,"%s/txt/headerinfo_row_%d_col_%d.txt",save_path,row,col);
                    p_hfile		= fopen(h_t_str,"r");
                    if(p_hfile)
                    {
                        //printf("header row %d\tcol %d\n",row,col);
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
                            
                            //printf("%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\n",t_row,t_col,t_level,t_boundary[0],t_boundary[1],t_grid_size,col_size,row_size);
                        }	
                        
                        while(fscanf(pfile,"%lf\t%lf\t%lf\n",&t_x,&t_y,&t_z) != EOF)
                        {
                            pos_col = (t_x - minX)/grid;
                            pos_row = (maxY - t_y)/grid;
                            
                            clip_pos_col = (t_x - t_boundary[0])/mt_grid;
                            clip_pos_row = (t_y - t_boundary[1])/mt_grid;
                            
                            
                            if(pos_row >= 0 && pos_row < row_count && pos_col >= 0 && pos_col < col_count && 
                               clip_pos_col > buffer_clip && clip_pos_col < col_size - buffer_clip && 
                               clip_pos_row > buffer_clip && clip_pos_row < row_size - buffer_clip)
                            {
                                long t_index = (long)((pos_row)*col_count + pos_col + 0.01);
                                value[t_index] = t_z;
                                value_pt[t_index] = 1;
                                pt_save[t_index].X = t_x;
                                pt_save[t_index].Y = t_y;
                                pt_save[t_index].Z = t_z;
                            }
                        }
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

    
	printf("start null\n");
	long count_null_cell = 0;
	int check_while = 0;
	while(check_while == 0)
	{
		count_null_cell = 0;
	
//#pragma omp parallel for shared(value,value_pt,col_count,row_count) private(index) reduction(+:count_null_cell)
		for (long index = 0; index < col_count*row_count; index++) 
		{
			int row, col;
			int check_size;
			long count_cell;
			int t_i, t_j;
			double sum_h;
			row = (int)(floor(index/col_count));
			col = index%col_count;
		
			if (value_pt[(long)row*(long)col_count + (long)col] == 0 && value[(long)row*(long)col_count + (long)col] > -9999)
			{
			
				check_size = 1;
				count_cell = 0;
				sum_h = 0;
				for (t_i = -check_size; t_i <= check_size;t_i++ ) 
				{
					for (t_j = -check_size; t_j <= check_size; t_j++) 
					{
						int index_row = row + t_i;
						int index_col = col + t_j;
						if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
						{
							if(value_pt[(long)index_row*(long)col_count + (long)index_col] == 1)
							{
								count_cell++;
								sum_h += value[(long)index_row*(long)col_count + (long)index_col];
							}
						}
					}
				}
				if (count_cell >= 6) 
				{
					double t_x, t_y;
					double pos_col,pos_row;
					
					value_pt[(long)row*(long)col_count + (long)col] = 1;
					value[(long)row*(long)col_count + (long)col]	  = (sum_h)/(count_cell);
					
					t_x = minX + col*grid;
					t_y = maxY - row*grid;
					
					if(pt_save[(long)row*(long)col_count + (long)col].Z == -9999)
					{
						pt_save[(long)row*(long)col_count + (long)col].X = t_x;
						pt_save[(long)row*(long)col_count + (long)col].Y = t_y;
						pt_save[(long)row*(long)col_count + (long)col].Z = value[(long)row*(long)col_count + (long)col]; 
					}
					count_null_cell ++;
				}
			}
		}
		
		if(count_null_cell == 0)
			check_while = 1;
		
		total_null_cell += count_null_cell;
	}
	total_mt_count = 0;
	
    printf("end null\n");
    
	if(total_search_count > 0)
	{
#pragma omp parallel for schedule(guided)
		for(long count = 0;count < total_search_count;count++)
		{
			double pos_row, pos_col, pos_row_mt, pos_col_mt;
			int check = 0;
			int row_interval, col_interval;
			double t_value;
			double query[2];
			double t_q[2];
			int iteration;
			char save_file[500];
			
			query[0]		= cal_gridpts_X[count];
			query[1]		= cal_gridpts_Y[count];
			
			t_q[0]			= query[0];
			t_q[1]			= query[1];
			
			row_interval	= 100;
			col_interval	= row_interval;
			
			pos_col = ((query[0] - minX)/grid);
			pos_row = ((maxY - query[1])/grid);
			
			pos_col_mt = ((query[0] - minX)/mt_grid);
			pos_row_mt = ((maxY - query[1])/mt_grid);
			
            int pos_index = (int)(pos_row*col_count + pos_col + 0.01);
            
			if (pos_col >= 0 && pos_col < col_count && pos_row >= 0 && pos_row < row_count)
			{
				t_value	 = value[pos_index];
				
				if(value_pt[pos_index] == 1)
				{
					check = 1;
				}
			}
			else {
				check = 1;
			}
			
			
			if(!check)
			{
				int check_while_stop = 0;
				iteration = row_interval;
				double height;
				int rendata;
				
				//IDW Interpolation
				height = FindNebPts_F_M_IDW(pt_save, row_count, col_count, grid, minX, minY, maxX, maxY, t_q[0], t_q[1],&rendata,iteration,iteration,ndim1,save_file);
				if(height == -9999)
					value[pos_index] = height;
				else if(height < minHeight-100 || height > maxHeight + 100)
				{
				}
				else 
				{
					value[pos_index] = height;
				}
			}
		}
	}
	
	if(cal_gridpts_X)
		free(cal_gridpts_X);
	if(cal_gridpts_Y)
		free(cal_gridpts_Y);
	
    printf("end interpolation\n");
	//smoothing
	float *value_sm = (float*)malloc(sizeof(float)*(long)col_count*(long)row_count);
				
//#pragma omp parallel for shared(value,value_sm,col_count,row_count) private(index)
	for (long index = 0; index < (long)col_count*(long)row_count; index++) 
	{
		int row, col;
		int check_size;
		long count_cell;
		long null_count_cell;
		int t_i, t_j;
		double sum_h, null_sum_h;
		
		row = (int)(floor(index/col_count));
		col = index%col_count;
		count_cell = 0;
		null_count_cell = 0;
		sum_h = 0;
		null_sum_h = 0;
		
		for (t_i = -1; t_i <= 1;t_i++ ) 
		{
			for (t_j = -1; t_j <= 1; t_j++) 
			{
				int index_row = row + t_i;
				int index_col = col + t_j;
				if(index_row >= 0 && index_row < row_count && index_col >= 0 && index_col < col_count)
				{
					if(value[(long)index_row*(long)col_count + (long)index_col] != -9999)
					{
						count_cell++;
						sum_h += value[(long)index_row*(long)col_count + (long)index_col];
					}
					
					if (value[(long)row*(long)col_count + (long)col] == -9999)
					{
						if(value[(long)index_row*(long)col_count + (long)index_col] != -9999 && t_i != 0 && t_j != 0)
						{
							null_count_cell++;
							null_sum_h += value[(long)index_row*(long)col_count + (long)index_col];
						}
					}

				}
			}
		}
		
		if(count_cell > 0)
			value_sm[(long)row*(long)col_count + (long)col] = sum_h / count_cell;
		else
			value_sm[(long)row*(long)col_count + (long)col] = value[(long)row*(long)col_count + (long)col];
		
		if(null_count_cell > 6)
			value_sm[(long)row*(long)col_count + (long)col] = null_sum_h / null_count_cell;
	}
	printf("end smoothing\n");
    
	free(value);
	value = (float*)malloc(sizeof(float)*(long)col_count*(long)row_count);
	memcpy(value,value_sm,sizeof(float)*(long)col_count*(long)row_count);
	free(value_sm);
	
	printf("end loading mpts, %ld\t%d\t%f\t%f\t%ld\t%ld!!\n",total_search_count,buffer_clip,minHeight,maxHeight,total_mt_count,total_check_count);
	
	sprintf(outfile, "%s/%s_dem.raw", save_path, Outputpath_name);
	fout	= fopen(outfile,"wb");
	fwrite(value,sizeof(float),(long)row_count*(long)col_count,fout);
	fclose(fout);
	
	sprintf(outfile, "%s/%s_matchtag.raw", save_path, Outputpath_name);
	fout	= fopen(outfile,"wb");
	fwrite(value_pt,sizeof(unsigned char),(long)row_count*(long)col_count,fout);
	fclose(fout);
	
	sprintf(DEM_header, "%s/%s_dem.hdr", save_path, Outputpath_name);
	Envihdr_writer(_param,DEM_header, col_count, row_count, grid, minX, maxY, Hemisphere,4);
	
	sprintf(DEM_header, "%s/%s_matchtag.hdr", save_path, Outputpath_name);
	Envihdr_writer(_param,DEM_header, col_count, row_count, grid, minX, maxY, Hemisphere,1);
	
	
	free(pt_save);
	free(value);
	free(value_pt);
}

void Envihdr_writer(TransParam _param, char *filename, int col_size, int row_size, double grid_size, double minX, double maxY, int NS_flag, int data_type)
{
	FILE* fid;
	fid = fopen(filename,"w");

	fprintf(fid,"ENVI\n");
	fprintf(fid,"description = {\n");
	fprintf(fid,"\tFile Imported into ENVI.}\n");
	fprintf(fid,"samples\t= %d\n",col_size);
	fprintf(fid,"lines\t= %d\n",row_size);
	fprintf(fid,"bands\t= 1\n");
	fprintf(fid,"header offset = 0\n");
	fprintf(fid,"file type = ENVI Standard\n");
	fprintf(fid,"data type = %d\n",data_type);
	fprintf(fid,"interleave = bsq\n");
	fprintf(fid,"sensor type = WorldView\n");
	fprintf(fid,"byte order = 0\n");

	if(_param.projection == 1)
	{
		fprintf(fid, "map info = {Polar Stereographic, 1, 1, %f, %f, %f, %f, WGS-84, units=Meters}\n", minX, maxY, grid_size, grid_size);
		
		
		if (NS_flag == true) {
			fprintf(fid, "projection info = {31, 637813.0, 635674.5, 70.000000, -45.000000, 0.0, 0.0, WGS-84, Polar Stereographic, units=Meters}\n");
		} else {
			fprintf(fid, "projection info = {31, 637813.0, 635674.5, -71.000000, 0.000000, 0.0, 0.0, WGS-84, Polar Stereographic, units=Meters}\n");
		}
	}
	else
	{
		if (NS_flag == true) {
			fprintf(fid, "map info = {UTM, 1, 1, %f, %f, %f, %f, %d, North, WGS-84, units=Meters}\n", minX, maxY, grid_size, grid_size,_param.zone);
		} else {
			fprintf(fid, "map info = {UTM, 1, 1, %f, %f, %f, %f, %d, South, WGS-84, units=Meters}\n", minX, maxY, grid_size, grid_size,_param.zone);
		}
	}

	fprintf(fid,"wavelength units = Unknown\n");
	fclose(fid);
}

CSize Envihdr_reader_seedDEM(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size)
{
	FILE* fid;
	CSize image_size;
	char bufstr[500];
	char t_str1[500],t_str2[500],t_str[500];
	int t_int1, t_int2;
	char* pos1;
	char* pos2;
	char* token = NULL;
	
	fid = fopen(filename,"r");
	
	while(!feof(fid))
	{
		fgets(bufstr,300,fid);
		if (strstr(bufstr,"samples")!=NULL)
			sscanf(bufstr,"samples = %d\n",&image_size.width);
		else if (strstr(bufstr,"lines")!=NULL)
			sscanf(bufstr,"lines = %d\n",&image_size.height);
		else if (strstr(bufstr,"map")!=NULL)
		{
			if(_param.projection == 1)
			{
				sscanf(bufstr, "map info = {%s %s %d%s %d%s %lf%s %lf%s %lf%s %lf%s %s %s\n", t_str1, t_str2, &t_int1, t_str, &t_int2, t_str, &(*minX), t_str, &(*maxY), t_str,
					   &(*grid_size), t_str, &(*grid_size), t_str, t_str, t_str);
			}
			else
			{
				sscanf(bufstr, "map info = {%s %d%s %d%s %lf%s %lf%s %lf%s %lf%s %d%s %s %s %s\n", t_str2, &t_int1, t_str, &t_int2, t_str, &(*minX), t_str, &(*maxY), t_str,
					   &(*grid_size), t_str, &(*grid_size), t_str, &t_int1, t_str, t_str, t_str, t_str);
			}
		}
	}
	fclose(fid);
	
	return image_size;
}

bool TFW_reader_seedDEM(char *filename, double *minX, double *maxY, double *grid_size)
{
	FILE* fid;
	double garbage;
	printf("%s\n",filename);
	fid = fopen(filename,"r");
	double temp;
	fscanf(fid,"%lf\n",&temp);
	printf("%f\n",temp);
	*grid_size = temp;
	fscanf(fid,"%lf\n",&garbage);
	fscanf(fid,"%lf\n",&garbage);
	fscanf(fid,"%lf\n",&garbage);
	
	fscanf(fid,"%lf\n",&temp);
	*minX = temp;
	fscanf(fid,"%lf\n",&temp);
	*maxY = temp;
	
	*minX = (*minX) - (*grid_size)/2.0;
	*maxY = (*maxY) - (*grid_size)/2.0;
	
	
	printf("%f\n",*minX);
	printf("%f\n",*maxY);
	printf("%f\n",*grid_size);
	
	
	fclose(fid);
	
	return true;
}

CSize Envihdr_reader(char *filename)
{
	FILE* fid;
	CSize image_size;
	char bufstr[500];
	
	fid = fopen(filename,"r");
	
	while(!feof(fid))
	{
		fgets(bufstr,300,fid);
		if (strstr(bufstr,"samples")!=NULL)
			sscanf(bufstr,"samples = %d\n",&image_size.width);
		else if (strstr(bufstr,"lines")!=NULL)
			sscanf(bufstr,"lines = %d\n",&image_size.height);
	}
	
	fclose(fid);
	
	return image_size;
}

double FindNebPts_F_M_IDW(NNXY *input, int row_size, int col_size, double grid, double minX, double minY, double maxX, double maxY, double X, double Y, int *numpts, int row_interval, int col_interval, int ndim1, char* path)
{
	double result;
    double row_pos, col_pos;
    int row,col;
	int size_pts;
	int check_stop = 0;
	int interval;
	int final_interval;
	double* diff;
	double* height;
	int count1,count2,count3,count4;
	
	col_pos = ((X - minX)/grid);
	row_pos = ((maxY - Y)/grid);
	
	interval = 3;
	
	
	while (check_stop == 0)
	{
		*numpts = 0;
		count1 = 0;
		count2 = 0;
		count3 = 0;
		count4 = 0;
		for(row = -interval;row<interval;row++)
		{
			for(col = -interval;col < interval ; col++)
			{
				int grid_pos = (int)((row_pos+row)*col_size + (col_pos+col));
				if(grid_pos >= 0 && grid_pos < row_size*col_size && 
				   row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size)// && radius > distance)
				{
					if(input[grid_pos].Z != -9999)
					{
						(*numpts)++;
						
						if (row >= 0 && row <interval && col >= 0 && col < interval) {
							count1++;
						}
						if (row >= 0 && row <interval && col < 0 && col >= -interval) {
							count2++;
						}
						if (row < 0 && row >= -interval && col < 0 && col >= -interval) {
							count3++;
						}
						if (row < 0 && row >= -interval && col >= 0 && col < interval) {
							count4++;
						}
					}
				}
			}
		}
		
		if (interval >= row_interval || ((*numpts) >= 20 && count1 > 2 && count2 > 2 && count3 > 2 && count4 > 2))
		{
			check_stop = 1;
			final_interval = interval;
		}
		else
			interval = interval + 2;
		
	}
	
	diff = (double*)malloc(sizeof(double)*(*numpts));
	height = (double*)malloc(sizeof(double)*(*numpts));
	*numpts = 0;
	for(row = -final_interval;row<final_interval;row++)
	{
		for(col = -final_interval;col < final_interval ; col++)
		{
			int grid_pos = (int)((row_pos+row)*col_size + (col_pos+col));
			int grid_pos_ref = (int)(row_pos*col_size + col_pos);
			if(grid_pos >= 0 && grid_pos < row_size*col_size && 
			   row_pos+row >= 0 && row_pos+row < row_size && col_pos+col >= 0 && col_pos+col < col_size && col != 0 && row != 0)// && radius > distance)
			{
				if(input[grid_pos].Z != -9999)
				{
					double x,y,x_ref,y_ref;
					x = input[grid_pos].X;
					y = input[grid_pos].Y;
					x_ref = X;
					y_ref = Y;
					diff[*numpts] = sqrt((x-x_ref)*(x-x_ref) + (y-y_ref)*(y-y_ref));
					height[*numpts] = input[grid_pos].Z;
					(*numpts)++;
				}
			}
		}
	}
	
	{
		double sum1, sum2;
		double p = 1.5;
	
		sum1 = 0;
		sum2 = 0;
		for(row=0;row<*numpts;row++)
		{
			sum1 += (height[row]/pow(diff[row],p));
			sum2 += (1.0/pow(diff[row],p));
		}
		
		if(sum2 > 0)
			result = sum1/sum2;
		else {
			result = -9999;
		}

	}
	
	free(diff);
	free(height);
	
	return result;
}

//orthogeneration
void orthogeneration(TransParam _param, ARGINFO args, char *ImageFilename, char *DEMFilename, char *Outputpath,int pair)
{
	if(args.RA_only) {
	    return;
	}
	char DEM_header[500];
	char RPCFilename[500];
	char OrthoFilename[500];
	char Ortho_header[500];
	
	char* tmp_chr; 
	double DEM_resolution, Image_resolution, Ortho_resolution;
	double OrthoGridFactor;
	double resolution[3];
	int impyramid_step;
	double** RPCs;
	double row_grid_size, col_grid_size;
	bool Hemisphere;
	double DEM_minX, DEM_maxY;
	double minLat,minLon;
	CSize DEM_size, Image_size;
	TransParam param;
	
	tmp_chr = remove_ext_ortho(ImageFilename);
	sprintf(RPCFilename,"%s.xml",tmp_chr);
	
	tmp_chr = remove_ext_ortho(DEMFilename);
	sprintf(DEM_header,"%s.hdr",tmp_chr);
	
	if(pair == 1)
	{
		sprintf(OrthoFilename, "%s/%s_ortho_image1.raw", Outputpath,args.Outputpath_name);
		sprintf(Ortho_header, "%s/%s_ortho_image1.hdr", Outputpath, args.Outputpath_name);
	}
	else
	{
		sprintf(OrthoFilename, "%s/%s_ortho_image2.raw", Outputpath,args.Outputpath_name);
		sprintf(Ortho_header, "%s/%s_ortho_image2.hdr", Outputpath, args.Outputpath_name);
	}
	
	printf("image = %s\n",ImageFilename);
	printf("rpc = %s\n",RPCFilename);
	printf("save = %s\n",Outputpath);
	printf("DEM = %s\n",DEMFilename);
	printf("DEM hdr= %s\n",DEM_header);
	printf("ortho = %s\n",OrthoFilename);
	printf("ortho hdr= %s\n",Ortho_header);
	
	
	Image_resolution = 0.5;
	
    // load RPCs info from xml file
    if(args.sensor_provider == 1)
    {
        RPCs			= OpenXMLFile_ortho(RPCFilename, &row_grid_size, &col_grid_size);
    }
    else
    {
        RPCs            = OpenXMLFile_Pleiades(RPCFilename);
    }
    
    
    if(!args.check_imageresolution)
    {
        if(args.sensor_provider == 1)
        {
            Image_resolution = (int)(((row_grid_size + col_grid_size)/2.0)*10 + 0.5)/10.0;
            
            if (Image_resolution < 0.75) {
                Image_resolution = 0.5;
            }
            else if(Image_resolution < 1.25)
                Image_resolution = 1.0;
            else
                Image_resolution = floor(Image_resolution);
        }
        else
            Image_resolution = 0.5;
    }
    else
    {
        Image_resolution  = args.image_resolution;
    }

    
	// load RPCs info from xml file
	/*RPCs			= OpenXMLFile_ortho(RPCFilename, &row_grid_size, &col_grid_size);
	
	  Image_resolution = (int)(((row_grid_size + col_grid_size)/2.0)*10 + 0.5)/10.0;
	
	  if (Image_resolution < 0.75) {
	  Image_resolution = 0.5;
	  }
	  else if(Image_resolution < 1.25)
	  Image_resolution = 1.0;
	*/
	printf("Image resolution %f\n",Image_resolution);
	
	minLat = RPCs[0][3];
	minLon = RPCs[0][2];
	
	//set projection conversion info.
    param.utm_zone   = args.utm_zone;
    _param.utm_zone  = args.utm_zone;
    
	SetTransParam(minLat, minLon, &Hemisphere, &_param);
	SetTransParam(minLat, minLon, &Hemisphere, &param);
	
	param.projection = _param.projection;
	printf("Hemis projection %d %d\n",Hemisphere, param.projection);

	// load DEM infor from hdr file.
	DEM_size = Envihdr_reader_DEM_ortho(param, DEM_header, &DEM_minX, &DEM_maxY, &DEM_resolution);
	
    if (!args.check_DEM_space)
        Ortho_resolution = DEM_resolution;
    else
        Ortho_resolution = args.DEM_space;
	
	resolution[0] = DEM_resolution;
	resolution[1] = Image_resolution;
	resolution[2] = Ortho_resolution;
    
	if(fabs(resolution[1] - resolution[2]) <= 1)
    {
        resolution[1] = Ortho_resolution;
    }
    
	OrthoGridFactor = resolution[2]/resolution[0];
	
	impyramid_step	= ceil(log(resolution[2]/resolution[1])/log(2));
	printf("impyramid_step %d\n",impyramid_step);
	
	CSize Orthoimagesize;
	double OrthoBoundary[4];
	
	// set generated orthoimage info by comparing DEM info
	SetOrthoBoundary_ortho(&Orthoimagesize, OrthoBoundary, RPCs, DEM_resolution, DEM_size, DEM_minX, DEM_maxY, param, Ortho_resolution);
	
	// set saving pointer for orthoimage
	uint16 *result_ortho	= (uint16*)calloc(Orthoimagesize.width*Orthoimagesize.height,sizeof(uint16));
	
	float *DEM_value	= NULL;
	
	double minmaxHeight[2];
	minmaxHeight[0]		= 99999;
	minmaxHeight[1]		= -99999;
	// load DEM value;
	DEM_value = LoadDEM_ortho(param, DEMFilename, DEM_header);
	
	printf("%d\n",DEM_size.width);
	printf("%d\n",DEM_size.height);
	printf("%f\n",DEM_minX);
	printf("%f\n",DEM_maxY);
	printf("%f\n",DEM_resolution);
	printf("%f\n",Ortho_resolution);		
	
	int i, j;
	for(i=0;i<DEM_size.height;i++)
	{
		for(j=0;j<DEM_size.width;j++)
		{
			if( (minmaxHeight[0] > DEM_value[i*DEM_size.width + j]) && (DEM_value[i*DEM_size.width + j] > -1000))
				minmaxHeight[0] = DEM_value[i*DEM_size.width + j];
			
			if( (minmaxHeight[1] < DEM_value[i*DEM_size.width + j]) && (DEM_value[i*DEM_size.width + j] > -1000))
				minmaxHeight[1] = DEM_value[i*DEM_size.width + j];
 		}
	}
	
	double subfactor					  = pow(6-impyramid_step,2.0);
	subfactor						= 1;
	int sub_height					= ceil(Orthoimagesize.height/subfactor);
	int sub_width					= ceil(Orthoimagesize.width/subfactor);
	int height_interval				= Ortho_resolution*sub_height;
	int width_interval				= Ortho_resolution*sub_width;
	
	int buffer_x					= Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
	int buffer_y					= Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
	
	int tile_count;
	for(i=0; i < subfactor;i++)
	{
		for(j=0; j<subfactor;j++)
		{
			char MEG[500];
			sprintf(MEG,"Tile %d %d, processing:%6.1f\n",j+1,i+1,((i)*subfactor+j+1)/subfactor/subfactor*100);
			printf("%s",MEG);

			double Y_size[2];
			Y_size[0]		 = OrthoBoundary[3] - (i+1)*height_interval;
			Y_size[1]		 = OrthoBoundary[3] - i*height_interval;
			double X_size[2];
			X_size[0]		 = OrthoBoundary[0] + j*width_interval;
			X_size[1]		 = OrthoBoundary[0] + (j+1)*width_interval;
			
			
			X_size[0]		= X_size[0] - buffer_x;
			Y_size[0]		= Y_size[0] - buffer_y;
			Y_size[1]		= Y_size[1] + buffer_y;
			X_size[1]		= X_size[1] + buffer_x;
			
			if( X_size[0] < OrthoBoundary[0])
				X_size[0]	= OrthoBoundary[0];
			
			if (X_size[1] > OrthoBoundary[2])
				X_size[1]	= OrthoBoundary[2];
			
			if (Y_size[0] < OrthoBoundary[1])
				Y_size[0]	= OrthoBoundary[1];
			
			if (Y_size[1] > OrthoBoundary[3])
				Y_size[1]	= OrthoBoundary[3];
			
            double subBoundary[4];
			subBoundary[0] = X_size[0];
			subBoundary[1] = Y_size[0];
			subBoundary[2] = X_size[1];
			subBoundary[3] = Y_size[1];
			
			D2DPOINT startpos_ori;
			char subsetImageFilename[500];
			char subsetImageFilename_hdr[500];
			CSize subsetsize;
			uint16 *subimage;
			bool check_subsetImage = false;
			
			subimage = subsetImage_ortho(param, RPCs, ImageFilename,
										 subBoundary,  minmaxHeight, &startpos_ori, subsetImageFilename, &subsetsize,&check_subsetImage);
			if(check_subsetImage)
			{
                int ori_impyramid_step = impyramid_step;
                if(impyramid_step > 0)
                    impyramid_step = 1;
                
				CSize data_size[impyramid_step+1];
				D2DPOINT startpos;
				char t_str[500];
				
				
				SetPySizes_ortho(data_size, subsetsize, impyramid_step);
				startpos.m_X	   = (double)(startpos_ori.m_X/pow(2,impyramid_step));		startpos.m_Y	   = (double)(startpos_ori.m_Y/pow(2,impyramid_step));
				
				uint16 *pyimg;
				
				if(impyramid_step > 0)
					pyimg = Preprocessing_ortho(ori_impyramid_step,data_size,subimage);
				
				Image_size	= data_size[impyramid_step];
				
				int col_size	= (int)((X_size[1] - X_size[0])/Ortho_resolution + 0.5);
				int row_size	= (int)((Y_size[1] - Y_size[0])/Ortho_resolution + 0.5);
				
//#pragma omp parallel for schedule(guided)
				for(int count = 0; count < col_size*row_size ; count++)
				{
					double row	= ((int)(floor(count/col_size)))*Ortho_resolution + Y_size[0];
					double col	= (count % col_size)*Ortho_resolution + X_size[0];

					int index1,index2,index3, index4;
					double t_col, t_row;
					int t_col_int, t_row_int;
					double dcol,drow;
					
					t_col		= (col - DEM_minX)/DEM_resolution;
					t_row		= (DEM_maxY - row)/DEM_resolution;
					
					t_col_int	= (int)(t_col + 0.01);
					t_row_int	= (int)(t_row + 0.01);
					
					dcol		= t_col - t_col_int;
					drow		= t_row - t_row_int;
					
                 	if(t_col_int >= 0 && t_col_int +1 < DEM_size.width && t_row_int >= 0 && t_row_int +1 < DEM_size.height)
					{
						double value1, value2, value3, value4, value;
						
						index1	= (t_col_int   ) + (t_row_int	)*DEM_size.width;
						index2	= (t_col_int +1) + (t_row_int	)*DEM_size.width;
						index3	= (t_col_int   ) + (t_row_int +1)*DEM_size.width;
						index4	= (t_col_int +1) + (t_row_int +1)*DEM_size.width;
						
						value1		= DEM_value[index1];
						value2		= DEM_value[index2];
						value3		= DEM_value[index3];
						value4		= DEM_value[index4];
						
						
						value		= value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
							+ value3*(1-dcol)*drow + value4*dcol*drow;
						
						D3DPOINT object;
						D2DPOINT objectXY;
						double imageparam[2] = {0.};
						object.m_X	= col;
						object.m_Y	= row;
						object.m_Z	= value;
						
						objectXY.m_X  = col;
						objectXY.m_Y  = row;
						
						if(value > -1000)
						{
							D2DPOINT wgsPt = ps2wgs_single(param, objectXY);
							D2DPOINT temp_pt;
							
							object.m_X	= wgsPt.m_X;
							object.m_Y	= wgsPt.m_Y;
							D2DPOINT image = GetObjectToImageRPC_single_ortho(RPCs, 2, imageparam, object);
							
							temp_pt		= OriginalToPyramid_single_ortho(image, startpos, impyramid_step);
							
							t_col		= temp_pt.m_X;
							t_row		= temp_pt.m_Y;
							
							t_col_int	= (int)(t_col + 0.01);
							t_row_int	= (int)(t_row + 0.01);
							
							dcol		= t_col - t_col_int;
							drow		= t_row - t_row_int;
							
							if(t_col_int >= 0 && t_col_int +1 < Image_size.width && t_row_int >= 0 && t_row_int +1 < Image_size.height
							   && (t_col_int +1) + (t_row_int +1)*Image_size.width < Image_size.width*Image_size.height)
							{
								double value1, value2, value3, value4, value;
								int index;
								index1	= (t_col_int   ) + (t_row_int	)*Image_size.width;
								index2	= (t_col_int +1) + (t_row_int	)*Image_size.width;
								index3	= (t_col_int   ) + (t_row_int +1)*Image_size.width;
								index4	= (t_col_int +1) + (t_row_int +1)*Image_size.width;
								
								if(impyramid_step > 0)
								{
									value1		= pyimg[index1];
									value2		= pyimg[index2];
									value3		= pyimg[index3];
									value4		= pyimg[index4];
								}
								else {
									value1		= subimage[index1];
									value2		= subimage[index2];
									value3		= subimage[index3];
									value4		= subimage[index4];
								}

								value		= value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
									+ value3*(1-dcol)*drow + value4*dcol*drow;
								
								t_col_int		= (int)((col - OrthoBoundary[0])/Ortho_resolution + 0.01);
								t_row_int		= (int)((OrthoBoundary[3] - row)/Ortho_resolution + 0.01);
								index			= t_col_int + t_row_int*Orthoimagesize.width;
								
                                if(t_col_int >= 0 && t_col_int < Orthoimagesize.width && t_row_int >= 0 && t_row_int < Orthoimagesize.height && index >= 0 && index < Orthoimagesize.width*Orthoimagesize.height)
								{
									result_ortho[index] = value1;
								}
							}
						}
					}
				}
				
				if(impyramid_step > 0)
					free(pyimg);
			}
			free(subimage);
		}
	}
	free(RPCs);
	free(DEM_value);
	
	FILE *fout	= fopen(OrthoFilename,"wb");
	fwrite(result_ortho,sizeof(uint16),Orthoimagesize.width*Orthoimagesize.height,fout);
	fclose(fout);
	
	Envihdr_writer(_param, Ortho_header, Orthoimagesize.width, Orthoimagesize.height, Ortho_resolution, OrthoBoundary[0], OrthoBoundary[3], Hemisphere, 12);
	free(result_ortho);
}

D2DPOINT OriginalToPyramid_single_ortho(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
	D2DPOINT out;
	
	out.m_X		 = (InCoord.m_X/pow(2,Pyramid_step)) - Startpos.m_X;
	out.m_Y		 = (InCoord.m_Y/pow(2,Pyramid_step)) - Startpos.m_Y;
	
	return out;
	
}

uint16 *Preprocessing_ortho(uint8 py_level, CSize *data_size, uint16 *subimg)
{
    uint16 *pyimg;
    int filter_size = pow(2,py_level)-1;
    if(filter_size < 3)
        filter_size = 3;
    
    printf("level %d\tfilter size %d\n",py_level,filter_size);
    
    pyimg = CreateImagePyramid_ortho(subimg,data_size[0],filter_size,(double)(filter_size/2.0));
    
    return pyimg;
    /*
	uint16 **pyimg;
	
	int i;
	
	FILE *pFile;
	
	pyimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
	
	for(i=0;i<py_level;i++)
	{
		if(i == 0)
			pyimg[i+1] = CreateImagePyramid_ortho(subimg,data_size[i],9,(double)(1.5));
		else
		{
			pyimg[i+1] = CreateImagePyramid_ortho(pyimg[i],data_size[i],9,(double)(1.5));
			free(pyimg[i]);
		}
	}
	
	return pyimg[py_level];
	*/
}

uint16* CreateImagePyramid_ortho(uint16* _input, CSize _img_size, int _filter_size, double _sigma)
{
	//_filter_size = 7, sigma = 1.6
	//long int i,j,r,c,l,k;
	double sigma = _sigma;
	double temp,scale;
	double sum = 0;
	double** GaussianFilter;
	CSize result_size;
	uint16* result_img;
	
	GaussianFilter = (double**)malloc(sizeof(double*)*_filter_size);
	for(int i=0;i<_filter_size;i++)
		GaussianFilter[i] = (double*)malloc(sizeof(double)*_filter_size);
	
	
	result_size.width = _img_size.width/2;
	result_size.height = _img_size.height/2;
	scale=sqrt(2*PI)*sigma;
	
	result_img = (uint16*)malloc(sizeof(uint16)*result_size.height*result_size.width);
	
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
		}
	}
	
#pragma omp parallel for private(temp) schedule(guided)
	for(long int r=0;r<result_size.height;r++)
	{
		for(long int c=0;c<result_size.width;c++)
		{
			temp = 0;
			
			for(int l=0;l<_filter_size;l++)
			{
				for(int k=0;k<_filter_size;k++)
				{
					//r'->2r+m, c'->2c+n
					if( (2*r + l-(int)(_filter_size/2)) >= 0 && (2*c + k-(int)(_filter_size/2)) >= 0 && 
						(2*r + l-(int)(_filter_size/2)) < _img_size.height && (2*c + k-(int)(_filter_size/2)) < _img_size.width)
					{
						temp += GaussianFilter[l][k]*_input[(2*r + l-(int)(_filter_size/2))*_img_size.width +(2*c + k-(int)(_filter_size/2))];
						
					}
				}
			}
			
			result_img[r*result_size.width + c] = (uint16)temp;
		}
	}
	
	for(int i=0;i<_filter_size;i++)
		if(GaussianFilter[i])
			free(GaussianFilter[i]);
	
	if(GaussianFilter)
		free(GaussianFilter);
	
	
	return result_img;
}


void SetPySizes_ortho(CSize *data_size, CSize subsetsize, int level)
{
	int i;
	
	data_size[0].height		= subsetsize.height;
	data_size[0].width		= subsetsize.width;
	for(i=0;i<level;i++)
	{
		data_size[i+1].width  = data_size[i].width/2;
		data_size[i+1].height = data_size[i].height/2;
	}
}


uint16 *subsetImage_ortho(TransParam transparam, double **RPCs, char *ImageFilename, 
						  double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, char *subsetImage, CSize* subsetsize, bool *ret)
{
	*ret = false;
	
	CSize Imagesize;
	uint16 *leftimage = NULL;
	if(GetImageSize_ortho(ImageFilename,&Imagesize))
	{
		int cols[2], rows[2];
		if(GetsubareaImage_ortho(transparam,RPCs,ImageFilename,&Imagesize,subBoundary,minmaxHeight,cols,rows) )
		{
			
			
			leftimage	= Readtiff_ortho(ImageFilename,Imagesize,cols,rows,subsetsize);
			
			
			startpos->m_X	= (double)(cols[0]);
			startpos->m_Y	= (double)(rows[0]);

			*ret		= true;
		}
	}
	
	return leftimage;
}

bool GetImageSize_ortho(char *filename, CSize *Imagesize)
{
	bool ret = false;
	
	char *ext;
	ext = strrchr(filename,'.');
	
	if(!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		
		TIFF *tif  = TIFFOpen(filename,"r");
		
		if(tif)
		{
			TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&Imagesize->width);
			TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&Imagesize->height);
			TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&Imagesize->width);
			TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&Imagesize->height);
			ret = true;
		}
		
		TIFFClose(tif);
	}
	else if(!strcmp("bin",ext+1))
	{
		char *tmp;
		tmp = remove_ext_ortho(filename);
		sprintf(tmp,"%s.hdr",tmp);
		*Imagesize = Envihdr_reader_ortho(tmp);
		
		ret = true;
		
	}
	return ret;
}

CSize Envihdr_reader_ortho(char *filename)
{
	FILE* fid;
	CSize image_size;
	char bufstr[500];
	
	fid = fopen(filename,"r");
	
	while(!feof(fid))
	{
		fgets(bufstr,300,fid);
		if (strstr(bufstr,"samples")!=NULL)
			sscanf(bufstr,"samples = %d\n",&image_size.width);
		else if (strstr(bufstr,"lines")!=NULL)
			sscanf(bufstr,"lines = %d\n",&image_size.height);
	}
	
	fclose(fid);
	
	return image_size;
}

bool GetsubareaImage_ortho(TransParam transparam, double **RPCs, char *ImageFilename, CSize *Imagesize,
						   double *subBoundary, double *minmaxHeight, int *cols, int *rows)
{
	bool ret = false;
	
	if(GetImageSize_ortho(ImageFilename,Imagesize))
	{
		int i;
		
		D3DPOINT t_pts[8];
		D3DPOINT *t_pts1;
		D2DPOINT *ImageCoord;
		int buffer, null_buffer;
		double ImageParam[2] = {0.};
		double minX =  1000000;
		double maxX = -1000000;
		double minY =  1000000;
		double maxY = -1000000;
		
		t_pts[0].m_X	= subBoundary[0];
		t_pts[1].m_X	= subBoundary[2];
		t_pts[2].m_X	= subBoundary[0];
		t_pts[3].m_X	= subBoundary[2];
		t_pts[4].m_X	= subBoundary[0];
		t_pts[5].m_X	= subBoundary[2];
		t_pts[6].m_X	= subBoundary[0];
		t_pts[7].m_X	= subBoundary[2];
		
		t_pts[0].m_Y	= subBoundary[1];
		t_pts[1].m_Y	= subBoundary[3];
		t_pts[2].m_Y	= subBoundary[1];
		t_pts[3].m_Y	= subBoundary[3];
		t_pts[4].m_Y	= subBoundary[3];
		t_pts[5].m_Y	= subBoundary[1];
		t_pts[6].m_Y	= subBoundary[3];
		t_pts[7].m_Y	= subBoundary[1];
		
		t_pts[0].m_Z	= minmaxHeight[0];
		t_pts[1].m_Z	= minmaxHeight[0];
		t_pts[2].m_Z	= minmaxHeight[1];
		t_pts[3].m_Z	= minmaxHeight[1];
		t_pts[4].m_Z	= minmaxHeight[0];
		t_pts[5].m_Z	= minmaxHeight[0];
		t_pts[6].m_Z	= minmaxHeight[1];
		t_pts[7].m_Z	= minmaxHeight[1];
		
		t_pts1			= ps2wgs_3D(transparam,8,t_pts);
		
		ImageCoord		= GetObjectToImageRPC_ortho(RPCs, 2, ImageParam, 8, t_pts1);
		
		for(i=0;i<8;i++)
		{
			if(minX > ImageCoord[i].m_X)
				minX	= ImageCoord[i].m_X;
			if(maxX < ImageCoord[i].m_X)
				maxX	= ImageCoord[i].m_X;
			if(minY > ImageCoord[i].m_Y)
				minY	= ImageCoord[i].m_Y;
			if(maxY < ImageCoord[i].m_Y)
				maxY	= ImageCoord[i].m_Y;
		}
		
		buffer				= 200;
		cols[0]				= (int)(ceil(minX)-buffer);
		cols[1]				= (int)(ceil(maxX)+buffer);
		rows[0]				= (int)(ceil(minY)-buffer);
		rows[1]				= (int)(ceil(maxY)+buffer);
		
		null_buffer			= 1;
		// Null pixel value remove
		if(cols[0]			<= null_buffer)
			cols[0]			= null_buffer;
		if(rows[0]			<= null_buffer)
			rows[0]			= null_buffer;
		if(cols[0]			> Imagesize->width - null_buffer)
			cols[0]			= Imagesize->width - null_buffer;
		if(rows[0]			> Imagesize->height - null_buffer)
			rows[0]			= Imagesize->height - null_buffer;
		
		if(cols[1]			<= null_buffer)
			cols[1]			= null_buffer;
		if(rows[1]			<= null_buffer)
			rows[1]			= null_buffer;
		if(cols[1]			> Imagesize->width - null_buffer)
			cols[1]			= Imagesize->width - null_buffer;
		if(rows[1]			> Imagesize->height - null_buffer)
			rows[1]			= Imagesize->height - null_buffer;
		
		free(t_pts1);
		free(ImageCoord);
		
		ret	= true;
	}
	
	return ret;
}

D2DPOINT* GetObjectToImageRPC_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP)
{
	D2DPOINT *IP;
	
	IP		= (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);

#pragma omp parallel for schedule(guided)
	for(int i=0;i<_numofpts;i++)
	{
		double L, P, H, Line, Samp;
		double deltaP = 0.0, deltaR = 0.0;
		double *Coeff;
		
		L		= (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
		P		= (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
		H		= (_GP[i].m_Z - _rpc[0][4])/_rpc[1][4];
        
        if(L < -10.0 || L > 10.0)
        {
            if(_GP[i].m_X > 0)
                _GP[i].m_X = _GP[i].m_X - 360;
            else
                _GP[i].m_X = _GP[i].m_X + 360;
            
            L		= (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
        }
        
        if(P < -10.0 || P > 10.0)
        {
            if(_GP[i].m_Y > 0)
                _GP[i].m_Y = _GP[i].m_Y - 360;
            else
                _GP[i].m_Y = _GP[i].m_Y + 360;
            
            P		= (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
        }
        
		Coeff	= (double*)malloc(sizeof(double)*4);

		for(int j=0;j<4;j++)
		{
			Coeff[j]	= _rpc[j+2][0]*1.0			+ _rpc[j+2][1]*L			+ _rpc[j+2][2]*P
				+ _rpc[j+2][3]*H			+ _rpc[j+2][4]*L*P			+ _rpc[j+2][5]*L*H
				+ _rpc[j+2][6]*P*H			+ _rpc[j+2][7]*L*L			+ _rpc[j+2][8]*P*P
				+ _rpc[j+2][9]*H*H			+ _rpc[j+2][10]*(P*L)*H		+ _rpc[j+2][11]*(L*L)*L
				+ _rpc[j+2][12]*(L*P)*P		+ _rpc[j+2][13]*(L*H)*H		+ _rpc[j+2][14]*(L*L)*P
				+ _rpc[j+2][15]*(P*P)*P		+ _rpc[j+2][16]*(P*H)*H		+ _rpc[j+2][17]*(L*L)*H
				+ _rpc[j+2][18]*(P*P)*H		+ _rpc[j+2][19]*(H*H)*H;
		}
		
		Line	 = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
		Samp	 = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
		
		switch(_numofparam)
		{
		case 2:
			deltaP		= _imageparam[0];
			deltaR		= _imageparam[1];
			break;
		case 6:
			deltaP		= _imageparam[0] + _imageparam[1]*Samp + _imageparam[2]*Line;
			deltaR		= _imageparam[3] + _imageparam[4]*Samp + _imageparam[5]*Line;
			break;
		}
		
		IP[i].m_Y		= deltaP + Line;
		IP[i].m_X		= deltaR + Samp;
		free(Coeff);
	}

	return IP;
}


uint16 *Readtiff_ortho(char *filename, CSize Imagesize, int *cols, int *rows, CSize *data_size)
{
	uint16 *out;
	FILE *bin;
	int check_ftype = 1; // 1 = tif, 2 = bin
	TIFF *tif = NULL;
	char *ext;
	ext = strrchr(filename,'.');
	
	if (!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
	{
		printf("tif open\n");
		tif	 = TIFFOpen(filename,"r");
		check_ftype = 1;
		printf("tif open end\n");
	}
	else if(!strcmp("bin",ext+1))
	{
		bin	 = fopen(filename,"rb");
		check_ftype = 2;
	}
	
	if(check_ftype == 1 && tif)
	{
		int i,j,row, col, tileW;
		
		tileW = -1;
		TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
		if(tileW < 0)
		{
			tsize_t scanline;
			tdata_t buf;
			uint16 s,nsamples;
			
			// scanline read
			data_size->width	= cols[1] - cols[0];
			data_size->height	= rows[1] - rows[0];
			
			long int data_length = data_size->height*data_size->width;
			out				= (uint16*)malloc(sizeof(uint16)*data_length);
			scanline		= TIFFScanlineSize(tif);
			
			buf				= _TIFFmalloc(scanline);
			
			TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&nsamples);
			
			for(s =0;s< nsamples;s++)
			{
				for (row=0;row<rows[0];row++)
					TIFFReadScanline(tif,buf,row,s);
				for (row=rows[0];row<rows[1];row++)
				{
					uint16* t_data;
					TIFFReadScanline(tif,buf,row,s);
					t_data = (uint16*)buf;
#pragma omp parallel for schedule(guided)
					for(int a = cols[0];a<cols[1];a++)
					{
						long int pos = (row-rows[0])*data_size->width + (a-cols[0]);
						out[pos] = t_data[a];
					}
				}
			}
			
			_TIFFfree(buf);
		}
		else
		{
			int tileL,count_W,count_L,starttileL,starttileW;
			int start_row,start_col,end_row,end_col;
			tdata_t buf;
			uint16* t_data;
			
			TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileW);
			TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileL);
			
			starttileL		= (int)(rows[0]/tileL);
			start_row		= starttileL*tileL;
			end_row			= ((int)(rows[1]/tileL)+1)*tileL;
			if(end_row > Imagesize.height)
				end_row = Imagesize.height;
			
			starttileW		= (int)(cols[0]/tileW);
			start_col		= starttileW*tileW;
			end_col			= ((int)(cols[1]/tileW)+1)*tileW;
			if(end_col > Imagesize.width)
				end_col = Imagesize.width;
			
			
			cols[0]			= start_col;
			cols[1]			= end_col;
			rows[0]			= start_row;
			rows[1]			= end_row;
			
			data_size->width = end_col - start_col;
			data_size->height= end_row - start_row;
			
			long int data_length = data_size->height*data_size->width;
			
			out				= (uint16*)malloc(sizeof(uint16)*data_length);
			
			buf				= _TIFFmalloc(TIFFTileSize(tif));
			
			count_L = (int)(data_size->height/tileL);
			count_W = (int)(data_size->width/tileW);
			
			for (row = 0; row < count_L; row ++)
			{
				for (col = 0; col < count_W; col ++)
				{
					TIFFReadTile(tif, buf, (col+starttileW)*tileW, (row+starttileL)*tileL, 0,0);
					t_data = (uint16*)buf;
#pragma omp parallel for private(i,j) schedule(guided)
					for (i=0;i<tileL;i++)
					{
						for (j=0;j<tileW;j++)
						{
							long int pos = ((row*tileL) + i)*data_size->width + ((col*tileL) + j);
							out[pos] = t_data[i*tileW + j];
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
		int r,c,a;
		data_size->width	= cols[1] - cols[0];
		data_size->height	= rows[1] - rows[0];
		
		long int data_length = data_size->height*data_size->width;
		
		out				= (uint16*)malloc(sizeof(uint16)*data_length);
		
		for(r = rows[0]; r < rows[1] ; r++)
		{
			fseek(bin,sizeof(uint16)*(r*Imagesize.width + cols[0]),SEEK_SET);
		
			uint16* t_data = (uint16*)malloc(sizeof(uint16)*data_size->width);
			
			fread(t_data,sizeof(uint16),data_size->width,bin);
			
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

D2DPOINT GetObjectToImageRPC_single_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, D3DPOINT _GP)
{
	D2DPOINT IP;
	
	int j;
	
	
	double L, P, H, Line, Samp, deltaP, deltaR;
	double *Coeff;
	deltaP = 0.0;
	deltaR = 0.0;
	
	L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
	P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
	H		= (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
    if(L < -10.0 || L > 10.0)
    {
        if(_GP.m_X > 0)
            _GP.m_X = _GP.m_X - 360;
        else
            _GP.m_X = _GP.m_X + 360;
        
        L		= (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    }
    
    if(P < -10.0 || P > 10.0)
    {
        if(_GP.m_Y > 0)
            _GP.m_Y = _GP.m_Y - 360;
        else
            _GP.m_Y = _GP.m_Y + 360;
        
        P		= (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    }
    
	Coeff	= (double*)malloc(sizeof(double)*4);

	for(j=0;j<4;j++)
	{
		Coeff[j]	= _rpc[j+2][0]*1.0			+ _rpc[j+2][1]*L			+ _rpc[j+2][2]*P
			+ _rpc[j+2][3]*H			+ _rpc[j+2][4]*L*P			+ _rpc[j+2][5]*L*H
			+ _rpc[j+2][6]*P*H			+ _rpc[j+2][7]*L*L			+ _rpc[j+2][8]*P*P
			+ _rpc[j+2][9]*H*H			+ _rpc[j+2][10]*(P*L)*H		+ _rpc[j+2][11]*(L*L)*L
			+ _rpc[j+2][12]*(L*P)*P		+ _rpc[j+2][13]*(L*H)*H		+ _rpc[j+2][14]*(L*L)*P
			+ _rpc[j+2][15]*(P*P)*P		+ _rpc[j+2][16]*(P*H)*H		+ _rpc[j+2][17]*(L*L)*H
			+ _rpc[j+2][18]*(P*P)*H		+ _rpc[j+2][19]*(H*H)*H;
	}
	
	Line	 = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
	Samp	 = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
	
	switch(_numofparam)
	{
	case 2:
		deltaP		= _imageparam[0];
		deltaR		= _imageparam[1];
		break;
	case 6:
		deltaP		= _imageparam[0] + _imageparam[1]*Samp + _imageparam[2]*Line;
		deltaR		= _imageparam[3] + _imageparam[4]*Samp + _imageparam[5]*Line;
		break;
	}
	
	IP.m_Y		= deltaP + Line;
	IP.m_X		= deltaR + Samp;
	free(Coeff);
	
	return IP;
}

bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary, 
							double **RPCs, double gridspace, CSize DEM_size, double minX, double maxY, TransParam param, double Ortho_resolution)
{
	double TopLeft[2];
	TopLeft[0]	= minX;
	TopLeft[1]	= maxY;
	double DEMboundary[4];
	DEMboundary[0]	= TopLeft[0];
	DEMboundary[1]	= TopLeft[1]-DEM_size.height*gridspace;
	DEMboundary[2]	= TopLeft[0]+DEM_size.width*gridspace;
	DEMboundary[3]	= TopLeft[1];
	double minLon, maxLon, minLat, maxLat;
	minLon			= -1.15*RPCs[1][2] + RPCs[0][2];
	maxLon			=  1.15*RPCs[1][2] + RPCs[0][2];
	minLat			= -1.15*RPCs[1][3] + RPCs[0][3];
	maxLat			=  1.15*RPCs[1][3] + RPCs[0][3];
	
	D2DPOINT *XY	= (D2DPOINT*)malloc(sizeof(D2DPOINT)*4);
	D2DPOINT *LonLat= (D2DPOINT*)malloc(sizeof(D2DPOINT)*4);
	double t_minX, t_maxX, t_minY, t_maxY;
	
	LonLat[0].m_X = minLon;
	LonLat[0].m_Y = minLat;
	LonLat[1].m_X = minLon;
	LonLat[1].m_Y = maxLat;
	LonLat[2].m_X = maxLon;
	LonLat[2].m_Y = maxLat;
	LonLat[3].m_X = maxLon;
	LonLat[3].m_Y = minLat;
	
	XY			= wgs2ps(param,4, LonLat);
	
	t_minX		= min(min(min(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
	t_maxX		= max(max(max(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
	t_minY		= min(min(min(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
	t_maxY		= max(max(max(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
	
	double ImageBoundary[4];
	ImageBoundary[0]	= floor(t_minX)-1;
	ImageBoundary[1]	= floor(t_minY)-1;
	ImageBoundary[2]	= ceil(t_maxX)+1;
	ImageBoundary[3]	= ceil(t_maxY)+1;
	
	Boundary[0]	= (max(DEMboundary[0],ImageBoundary[0]));
	Boundary[1]	= (max(DEMboundary[1],ImageBoundary[1]));
	Boundary[2]	= (min(DEMboundary[2],ImageBoundary[2]));
	Boundary[3]	= (min(DEMboundary[3],ImageBoundary[3]));
	
	Imagesize->height	 = ceil(fabs(Boundary[3] - Boundary[1])/Ortho_resolution);
	Imagesize->width	 = ceil(fabs(Boundary[2] - Boundary[0])/Ortho_resolution);
    
	
	printf("orthoimage height width %d \t%d\t %f\t%f\n",Imagesize->height,Imagesize->width,fabs(DEMboundary[3] - DEMboundary[1])/Ortho_resolution,fabs(DEMboundary[2] - DEMboundary[0])/Ortho_resolution);
	return true;
}


float *LoadDEM_ortho(TransParam param, char *DEM_path, char* hdr_path)
{
	double minX, maxX, minY,maxY,grid_size;
	CSize seeddem_size;
	float *seeddem = NULL;
	printf("hdr path %s\n",hdr_path);
	seeddem_size  = Envihdr_reader_DEM_ortho(param, hdr_path, &minX, &maxY, &grid_size);
	maxX	= minX + grid_size*seeddem_size.width;
	minY	= maxY - grid_size*seeddem_size.height;
	
	printf("seeddem cal\n");
	
	FILE *bin;
	
	bin = fopen(DEM_path,"rb");
	seeddem = (float*)malloc(sizeof(float)*seeddem_size.width*seeddem_size.height);
	fread(seeddem,sizeof(float),seeddem_size.width*seeddem_size.height,bin);
	
	printf("seeddem end\n");
	
	return seeddem;
}

double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c)
{
	double** out = NULL;
	
	FILE *pFile;
	char temp_str[1000];
	char linestr[1000];
	int i;
	char* pos1;
	char* pos2;
	char* token = NULL;
	
	double aa;
	
	pFile			= fopen(_filename,"r");
	if(pFile)
	{
		out = (double**)malloc(sizeof(double*)*7);
		while(!feof(pFile))
		{
			fscanf(pFile,"%s",temp_str);
			if(strcmp(temp_str,"<IMAGE>") == 0)
			{
				for(i=0;i<21;i++)
					fgets(temp_str,sizeof(temp_str),pFile);
				
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				*gsd_r			= atof(pos2);
				
				//for(i=0;i<2;i++)
				//	fgets(temp_str,sizeof(temp_str),pFile);
				
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				*gsd_c			= atof(pos2);
				
			}
			
			if(strcmp(temp_str,"<RPB>")==0)
			{
				for(i=0;i<5;i++)
					fgets(temp_str,sizeof(temp_str),pFile);
				
				out[6] = (double*)malloc(sizeof(double)*2);
				for(i=0;i<2;i++)
				{
					
					fgets(linestr,sizeof(linestr),pFile);
					pos1 = strstr(linestr,">")+1;
					pos2 = strtok(pos1,"<");
					out[6][i]			= atof(pos2);
					
				}
				
				out[0] = (double*)malloc(sizeof(double)*5);
				for(i=0;i<5;i++)
				{
					fgets(linestr,sizeof(linestr),pFile);
					pos1 = strstr(linestr,">")+1;
					pos2 = strtok(pos1,"<");
					out[0][i]			= atof(pos2);
				}
				aa						= out[0][2];
				out[0][2]				= out[0][3];
				out[0][3]				= aa;
				
				out[1] = (double*)malloc(sizeof(double)*5);
				for(i=0;i<5;i++)
				{
					fgets(linestr,sizeof(linestr),pFile);
					pos1 = strstr(linestr,">")+1;
					pos2 = strtok(pos1,"<");
					out[1][i]			= atof(pos2);
				}
				
				aa						= out[1][2];
				out[1][2]				= out[1][3];
				out[1][3]				= aa;
				
				fgets(temp_str,sizeof(temp_str),pFile);
				out[2] = (double*)malloc(sizeof(double)*20);
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				token = strtok(pos2," ");
				i=0;
				while(token != NULL)
				{
					out[2][i]			= atof(token);
					//printf("%e\n",out[2][i]);
					token = strtok(NULL," ");
					i++;
				}
				
				fgets(temp_str,sizeof(temp_str),pFile);
				fgets(temp_str,sizeof(temp_str),pFile);
				out[3] = (double*)malloc(sizeof(double)*20);
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				token = strtok(pos2," ");
				i=0;
				while(token != NULL)
				{
					out[3][i]			= atof(token);
					//printf("%e\n",out[3][i]);
					token = strtok(NULL," ");
					i++;
				}
				
				fgets(temp_str,sizeof(temp_str),pFile);
				fgets(temp_str,sizeof(temp_str),pFile);
				out[4] = (double*)malloc(sizeof(double)*20);
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				token = strtok(pos2," ");
				i=0;
				while(token != NULL)
				{
					out[4][i]			= atof(token);
					//printf("%e\n",out[4][i]);
					token = strtok(NULL," ");
					i++;
				}
				
				fgets(temp_str,sizeof(temp_str),pFile);
				fgets(temp_str,sizeof(temp_str),pFile);
				out[5] = (double*)malloc(sizeof(double)*20);
				fgets(linestr,sizeof(linestr),pFile);
				pos1 = strstr(linestr,">")+1;
				pos2 = strtok(pos1,"<");
				token = strtok(pos2," ");
				i=0;
				while(token != NULL)
				{
					out[5][i]			= atof(token);
					//printf("%e\n",out[5][i]);
					token = strtok(NULL," ");
					i++;
				}
			}
		}
		fclose(pFile);
		
		if(pos1)
			pos1 = NULL;
		if(pos2)
			pos2 = NULL;
		if(token)
			token = NULL;
	}
	return out;
}


CSize Envihdr_reader_DEM_ortho(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size)
{
	FILE* fid;
	CSize image_size;
	char bufstr[500];
	char t_str1[500],t_str2[500],t_str[500];
	int t_int1, t_int2;
	
	fid = fopen(filename,"r");
	
	while(!feof(fid))
	{
		fgets(bufstr,300,fid);
		if (strstr(bufstr,"samples")!=NULL)
			sscanf(bufstr,"samples = %d\n",&image_size.width);
		else if (strstr(bufstr,"lines")!=NULL)
			sscanf(bufstr,"lines = %d\n",&image_size.height);
		else if (strstr(bufstr,"map")!=NULL)
		{
			if(_param.projection == 1)
			{
				sscanf(bufstr, "map info = {%s %s %d%s %d%s %lf%s %lf%s %lf%s %lf%s %s %s\n", t_str1, t_str2, &t_int1, t_str, &t_int2, t_str, &(*minX), t_str, &(*maxY), t_str,
					   &(*grid_size), t_str, &(*grid_size), t_str, t_str, t_str);
			}
			else
			{
				sscanf(bufstr, "map info = {%s %d%s %d%s %lf%s %lf%s %lf%s %lf%s %d%s %s %s %s\n", t_str2, &t_int1, t_str, &t_int2, t_str, &(*minX), t_str, &(*maxY), t_str,
					   &(*grid_size), t_str, &(*grid_size), t_str, &t_int1, t_str, t_str, t_str, t_str);
			}
		}
	}
	fclose(fid);
	
	return image_size;
	
}

char* remove_ext_ortho(char* mystr)
{
	char *retstr;
	char *lastdot;
	if (mystr == NULL)
		return NULL;
	if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
		return NULL;
	strcpy (retstr, mystr);
	lastdot = strrchr (retstr, '.');
	if (lastdot != NULL)
		*lastdot = '\0';
	return retstr;
}
