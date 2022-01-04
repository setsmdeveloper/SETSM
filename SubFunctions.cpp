
//  SubFunctions.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//
#include "SubFunctions.hpp"
#include "readtiff.hpp"

int numcols[7] = {6,5,5,5,4,3,0};

char* remove_ext(const char* mystr)
{
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
        return NULL;
    if ((retstr = (char*)malloc(strlen (mystr) + 1)) == NULL)
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
        int lenth = end  - start+1;
        t_name = (char*)(malloc(sizeof(char)*(lenth+1)));
        int path_size = strlen(t_name);
        for (int i = 0; i < lenth; i++) {
            t_name[i] = lastSlash[i+start];
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

int Maketmpfolders(ProInfo *info)
{
    int check_folder = 1;
    
    if(!info->check_checktiff)
    {
        int status = mkdir(info->save_filepath,0777);
        if (opendir(info->save_filepath) == NULL)
        {
            if (status == -1)
            {
                check_folder = 0;
                printf("Outpath of '%s' cannot make, please check outpath!!\n",info->save_filepath);
                exit(1);
            }
        }
        
        char temp_filepath[500];
        sprintf(temp_filepath,"%s/txt",info->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(temp_filepath,"%s/tif",info->save_filepath);
        mkdir(temp_filepath,0777);
        sprintf(temp_filepath,"%s/tmp",info->save_filepath);
        mkdir(temp_filepath,0777);
    }
    return check_folder;
}

double getSystemMemory()
{
    struct rlimit lim;
    if( getrlimit(RLIMIT_RSS, &lim) == 0) {
      if(lim.rlim_cur == RLIM_INFINITY) {
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        return (pages * page_size) / (double)pwrtwo(30) ;
      }
      else {
        return (double)lim.rlim_cur/pwrtwo(30);
      }
    }
    else {
        return -1.0;
    }
}

/** Print maximum memory usage so far in human readable units */
void printMaxMemUsage() {
    // rusage.ru_maxrss in KB
    struct rusage ru;
    if(getrusage(RUSAGE_SELF, &ru)) {
        printf("WARNING: failed to get maximum memory usage\n");
        return;
    }

    auto max_rss = ru.ru_maxrss;
    double scaled = 0;
    const char *unit = NULL;

    if(max_rss < 1024) {
        scaled = max_rss / 1.0;
        unit = "KB";
    } else if(max_rss < 1024 * 1024) {
        scaled = max_rss / 1024.0;
        unit = "MB";
    } else {
        scaled = max_rss / (1024.0 * 1024.0);
        unit = "GB";
    }
    printf("Maximum memory usage: %.2f %s\n", scaled, unit);
}

bool GetRAinfo(ProInfo *proinfo, const char* RAfile, double **Imageparams, CPairInfo &pairinfo)
{
    bool check_load_RA = false;
    bool check_load_RA_each = false;
    
    FILE *pFile_echo  = fopen(RAfile,"r");
    if(pFile_echo)
    {
        int ti=0;
        while(!feof(pFile_echo))
        {
            int temp_ref,temp_tar;
            double bh_ratio;
            //for(int ti = 0; ti < proinfo->number_of_images ; ti++)
            {
                fscanf(pFile_echo,"%d\t%d\t%lf\t%lf\t%lf\n",&temp_ref,&temp_tar,&Imageparams[ti][0],&Imageparams[ti][1],&bh_ratio);
                /*if(ti == 0)
                    proinfo->check_selected_image[ti] = true;
                else*/
                {
                    if(Imageparams[ti][0] == 0 || Imageparams[ti][1] == 0)
                        check_load_RA_each = true;
                }
                printf("RAinfo %d\t%d\t%d\t%f\t%f\n",ti,temp_ref,temp_tar,Imageparams[ti][0],Imageparams[ti][1]);
                pairinfo.SetPairs(ti, temp_ref, temp_tar);
                pairinfo.SetBHratio(ti, bh_ratio);
            }
            ti++;
            pairinfo.SetSelectNumberOfPairs(ti);
        }
        fclose(pFile_echo);
        
        if(!check_load_RA_each && ti > 0)
            check_load_RA = true;
    }
    
    return check_load_RA;
}

bool GetRAinfoFromEcho(ProInfo *proinfo, const char* echofile, double **Imageparams, CPairInfo &pairinfo)
{
    bool check_load_RA = false;
    bool check_load_RA_each = false;
    FILE* pFile_echo  = fopen(echofile,"r");
    if(pFile_echo)
    {
        int ti=0;
        while(!feof(pFile_echo))
        {
            int temp_ref,temp_tar;
            double bh_ratio;
            //for(int ti = 0; ti < proinfo->number_of_images ; ti++)
            {
                char bufstr[500];
                fgets(bufstr,500,pFile_echo);
                if (strstr(bufstr,"RA param X")!=NULL)
                {
                    sscanf(bufstr,"RA param X = %lf\tY = %lf\t%d\t%d\t%lf\n",&Imageparams[ti][0],&Imageparams[ti][1],&temp_ref,&temp_tar,&bh_ratio);
                    /*if(ti == 0)
                        proinfo->check_selected_image[ti] = true;
                    else*/
                    {
                        if(Imageparams[ti][0] == 0 || Imageparams[ti][1] == 0)
                            check_load_RA_each = true;
                    }
                    printf("echo %d\t%d\t%d\t%f\t%f\n",ti,temp_ref,temp_tar,Imageparams[ti][0],Imageparams[ti][1]);
                    pairinfo.SetPairs(ti, temp_ref, temp_tar);
                    pairinfo.SetBHratio(ti, bh_ratio);
                    ti++;
                }
            }
            
            pairinfo.SetSelectNumberOfPairs(ti);
        }
        fclose(pFile_echo);
        
        if(!check_load_RA_each && ti > 0)
            check_load_RA = true;
    }
    
    return check_load_RA;
}

bool GetsubareaImage(const int sensor_type, const FrameInfo m_frameinfo, const int ti, const TransParam transparam, const double *ImageParam, const double * const *RPCs, char *ImageFilename, CSize Imagesize, const double *subBoundary, const double *minmaxHeight, long *cols, long *rows)
{
    bool ret = false;
    
    D3DPOINT t_pts[8];
    
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
    
    t_pts[0].m_Z    = minmaxHeight[0];
    t_pts[1].m_Z    = minmaxHeight[0];
    t_pts[2].m_Z    = minmaxHeight[1];
    t_pts[3].m_Z    = minmaxHeight[1];
    t_pts[4].m_Z    = minmaxHeight[0];
    t_pts[5].m_Z    = minmaxHeight[0];
    t_pts[6].m_Z    = minmaxHeight[1];
    t_pts[7].m_Z    = minmaxHeight[1];
    
    printf("subBoundary %f\t%f\t%f\t%f\n%d\t%f\t%f\t%f\t%f\n",subBoundary[0],subBoundary[1],subBoundary[2],subBoundary[3],transparam.pm,transparam.sa,transparam.sb,transparam.e2,transparam.c);
    
    D2DPOINT *ImageCoord = NULL;
    if(sensor_type == SB)
    {
        D3DPOINT *t_pts1          = ps2wgs_3D(transparam,8,t_pts);
        ImageCoord      = GetObjectToImageRPC(RPCs, 2, ImageParam, 8, t_pts1);
        
        free(t_pts1);
    }
    else
    {
        D2DPOINT *t_image     = GetPhotoCoordinate(t_pts, m_frameinfo.Photoinfo[ti], 8, m_frameinfo.m_Camera, m_frameinfo.Photoinfo[ti].m_Rm);
        ImageCoord  = PhotoToImage(t_image,8,m_frameinfo.m_Camera.m_CCDSize,m_frameinfo.m_Camera.m_ImageSize);
        
        free(t_image);
    }
    
    double minX =  1000000;
    double maxX = -1000000;
    double minY =  1000000;
    double maxY = -1000000;
    for(int i=0;i<8;i++)
    {
        if(minX > ImageCoord[i].m_X)
            minX    = ImageCoord[i].m_X;
        if(maxX < ImageCoord[i].m_X)
            maxX    = ImageCoord[i].m_X;
        if(minY > ImageCoord[i].m_Y)
            minY    = ImageCoord[i].m_Y;
        if(maxY < ImageCoord[i].m_Y)
            maxY    = ImageCoord[i].m_Y;
    }
    
    const int buffer              = 200;
    cols[0]             = (int)(ceil(minX)-buffer);
    cols[1]             = (int)(ceil(maxX)+buffer);
    rows[0]             = (int)(ceil(minY)-buffer);
    rows[1]             = (int)(ceil(maxY)+buffer);
    
    const int null_buffer         = 1;
    // Null pixel value remove
    if(cols[0]          <= null_buffer)
        cols[0]         = null_buffer;
    if(rows[0]          <= null_buffer)
        rows[0]         = null_buffer;
    if(cols[0]          > Imagesize.width - null_buffer)
        cols[0]         = Imagesize.width - null_buffer;
    if(rows[0]          > Imagesize.height - null_buffer)
        rows[0]         = Imagesize.height - null_buffer;
    
    if(cols[1]          <= null_buffer)
        cols[1]         = null_buffer;
    if(rows[1]          <= null_buffer)
        rows[1]         = null_buffer;
    if(cols[1]          > Imagesize.width - null_buffer)
        cols[1]         = Imagesize.width - null_buffer;
    if(rows[1]          > Imagesize.height - null_buffer)
        rows[1]         = Imagesize.height - null_buffer;
    
    printf("cols rows %ld\t%ld\t%ld\t%ld\n",cols[0],cols[1],rows[0],rows[1]);
    
    free(ImageCoord);
    
    ret = true;
    
    return ret;
}

bool GetImageSize(char *filename, CSize *Imagesize)
{
    bool ret = false;
    
    char *ext = strrchr(filename,'.');
    
    if(!strcmp("tif",ext+1) || !strcmp("TIF",ext+1))
    {
        *Imagesize = ReadGeotiff_info(filename, NULL, NULL, NULL);
        ret = true;
    }
    else if(!strcmp("bin",ext+1))
    {
        char *tmp = remove_ext(filename);
        char extension[] = ".hdr";
        char *tmp2;
        int length = strlen(tmp) + strlen(extension) + 1;
        tmp2 = (char *)malloc(length);
        strcpy(tmp2, tmp);
        strcat(tmp2, extension);
        free(tmp);
        tmp = tmp2;
        *Imagesize = Envihdr_reader(tmp);
        free(tmp);
        
        ret = true;
        
    }
    return ret;
}

bool GetsubareaImage_GeoTiff(ProInfo proinfo, char *ImageFilename, CSize Imagesize, double *subBoundary, long *cols,long *rows)
{
    bool ret = false;
    
    D2DPOINT t_pts[8];
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
    
    D2DPOINT *ImageCoord        = GetObjectToImage(8, t_pts,subBoundary,proinfo.resolution);
    
    double minX =  1000000;
    double maxX = -1000000;
    double minY =  1000000;
    double maxY = -1000000;
    
    for(int i=0;i<8;i++)
    {
        if(minX > ImageCoord[i].m_X)
            minX    = ImageCoord[i].m_X;
        if(maxX < ImageCoord[i].m_X)
            maxX    = ImageCoord[i].m_X;
        if(minY > ImageCoord[i].m_Y)
            minY    = ImageCoord[i].m_Y;
        if(maxY < ImageCoord[i].m_Y)
            maxY    = ImageCoord[i].m_Y;
    }
    
    const int buffer                = 0;
    cols[0]                = (int)(ceil(minX)-buffer);
    cols[1]                = (int)(ceil(maxX)+buffer);
    rows[0]                = (int)(ceil(minY)-buffer);
    rows[1]                = (int)(ceil(maxY)+buffer);
    
    const int null_buffer            = 0;
    // Null pixel value remove
    if(cols[0]            <= null_buffer)
        cols[0]            = null_buffer;
    if(rows[0]            <= null_buffer)
        rows[0]            = null_buffer;
    if(cols[0]            > Imagesize.width - null_buffer)
        cols[0]            = Imagesize.width - null_buffer;
    if(rows[0]            > Imagesize.height - null_buffer)
        rows[0]            = Imagesize.height - null_buffer;
    
    if(cols[1]            <= null_buffer)
        cols[1]            = null_buffer;
    if(rows[1]            <= null_buffer)
        rows[1]            = null_buffer;
    if(cols[1]            > Imagesize.width - null_buffer)
        cols[1]            = Imagesize.width - null_buffer;
    if(rows[1]            > Imagesize.height - null_buffer)
        rows[1]            = Imagesize.height - null_buffer;
    
    printf("cols rows %ld\t%ld\t%ld\t%ld\n",cols[0],cols[1],rows[0],rows[1]);
    
    free(ImageCoord);
    
    ret    = true;
    
    return ret;
}

bool GetsubareaImage_GeoTiff_tile(ProInfo proinfo, char *ImageFilename, CSize Imagesize, double *ImgBoundary, double *subBoundary, long *cols,long *rows)
{
    bool ret = false;
    
    D2DPOINT t_pts[8];
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
    
    D2DPOINT *ImageCoord        = GetObjectToImage(8, t_pts,ImgBoundary,proinfo.resolution);
    
    //printf("subBoundary = %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]);
    //printf("ImgBoundary = %f\t%f\t%f\t%f\n", ImgBoundary[0], ImgBoundary[1], ImgBoundary[2], ImgBoundary[3]);
    
    double minX =  1000000;
    double maxX = -1000000;
    double minY =  1000000;
    double maxY = -1000000;
    
    for(int i=0;i<8;i++)
    {
        if(minX > ImageCoord[i].m_X)
            minX    = ImageCoord[i].m_X;
        if(maxX < ImageCoord[i].m_X)
            maxX    = ImageCoord[i].m_X;
        if(minY > ImageCoord[i].m_Y)
            minY    = ImageCoord[i].m_Y;
        if(maxY < ImageCoord[i].m_Y)
            maxY    = ImageCoord[i].m_Y;
    }
    
    const int buffer                = 0;
    cols[0]                = (int)(ceil(minX)-buffer);
    cols[1]                = (int)(ceil(maxX)+buffer);
    rows[0]                = (int)(ceil(minY)-buffer);
    rows[1]                = (int)(ceil(maxY)+buffer);
    
    const int null_buffer            = 0;
    // Null pixel value remove
    if(cols[0]            <= null_buffer)
        cols[0]            = null_buffer;
    if(rows[0]            <= null_buffer)
        rows[0]            = null_buffer;
    if(cols[0]            > Imagesize.width - null_buffer)
        cols[0]            = Imagesize.width - null_buffer;
    if(rows[0]            > Imagesize.height - null_buffer)
        rows[0]            = Imagesize.height - null_buffer;
    
    if(cols[1]            <= null_buffer)
        cols[1]            = null_buffer;
    if(rows[1]            <= null_buffer)
        rows[1]            = null_buffer;
    if(cols[1]            > Imagesize.width - null_buffer)
        cols[1]            = Imagesize.width - null_buffer;
    if(rows[1]            > Imagesize.height - null_buffer)
        rows[1]            = Imagesize.height - null_buffer;
    
    //printf("cols rows %ld\t%ld\t%ld\t%ld\n",cols[0],cols[1],rows[0],rows[1]);
    
    free(ImageCoord);
    
    ret    = true;
    
    return ret;
}

bool GetsubareaImage_GeoTiff_ortho(double grid_size, CSize Imagesize, double *ImgBoundary, double *subBoundary, long *cols,long *rows)
{
    bool ret = false;
    
    D2DPOINT t_pts[8];
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
    
    D2DPOINT *ImageCoord        = GetObjectToImage(8, t_pts,ImgBoundary,grid_size);
    
    //printf("subBoundary = %f\t%f\t%f\t%f\n", subBoundary[0], subBoundary[1], subBoundary[2], subBoundary[3]);
    //printf("ImgBoundary = %f\t%f\t%f\t%f\n", ImgBoundary[0], ImgBoundary[1], ImgBoundary[2], ImgBoundary[3]);
    
    double minX =  1000000;
    double maxX = -1000000;
    double minY =  1000000;
    double maxY = -1000000;
    
    for(int i=0;i<8;i++)
    {
        if(minX > ImageCoord[i].m_X)
            minX    = ImageCoord[i].m_X;
        if(maxX < ImageCoord[i].m_X)
            maxX    = ImageCoord[i].m_X;
        if(minY > ImageCoord[i].m_Y)
            minY    = ImageCoord[i].m_Y;
        if(maxY < ImageCoord[i].m_Y)
            maxY    = ImageCoord[i].m_Y;
    }
    
    const int buffer                = 0;
    cols[0]                = (int)(ceil(minX)-buffer);
    cols[1]                = (int)(ceil(maxX)+buffer);
    rows[0]                = (int)(ceil(minY)-buffer);
    rows[1]                = (int)(ceil(maxY)+buffer);
    
    const int null_buffer            = 0;
    // Null pixel value remove
    if(cols[0]            <= null_buffer)
        cols[0]            = null_buffer;
    if(rows[0]            <= null_buffer)
        rows[0]            = null_buffer;
    if(cols[0]            > Imagesize.width - null_buffer)
        cols[0]            = Imagesize.width - null_buffer;
    if(rows[0]            > Imagesize.height - null_buffer)
        rows[0]            = Imagesize.height - null_buffer;
    
    if(cols[1]            <= null_buffer)
        cols[1]            = null_buffer;
    if(rows[1]            <= null_buffer)
        rows[1]            = null_buffer;
    if(cols[1]            > Imagesize.width - null_buffer)
        cols[1]            = Imagesize.width - null_buffer;
    if(rows[1]            > Imagesize.height - null_buffer)
        rows[1]            = Imagesize.height - null_buffer;
    
    //printf("cols rows %ld\t%ld\t%ld\t%ld\n",cols[0],cols[1],rows[0],rows[1]);
    
    free(ImageCoord);
    
    ret    = true;
    
    return ret;
}

float* GetDEMValue(char *GIMP_path,CSize seeddem_size)
{
    double minX, maxX, minY,maxY,a_minX,a_maxX,a_minY,a_maxY;
    
    float* seeddem = NULL;
    char* hdr_path;
    FILE *bin;
    TIFF *tif;
    char save_DEMfile[500];
    int i,j;
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
    
    //float *seeddem = NULL;
    printf("DEM file = %s\n",GIMP_path);
    if(check_ftype == 2)
    {
        long int data_length = (long)seeddem_size.width*(long)seeddem_size.height;
        seeddem = (float*)calloc(sizeof(float),data_length);
        bin = fopen(GIMP_path,"rb");
        fread(seeddem,sizeof(float),data_length,bin);
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
        seeddem = Readtiff_T(GIMP_path,LImagesize,cols,rows,&data_size,type);
        printf("tif open\n");

        free(LImagesize);
    }
    
    return seeddem;
}

bool OpenProject(char* _filename, ProInfo *info, ARGINFO args)
{
    bool bopened = false;
    FILE *pFile;
    char bufstr[500];

    info->IsRA          = true;
    info->IsSaveStep    = false;
    info->pre_DEMtif    = false;
    info->check_tiles_SR= false;
    info->check_tiles_ER= false;
    info->check_tiles_SC= false;
    info->check_tiles_EC= false;
    info->check_minH    = false;
    info->check_maxH    = false;
    info->check_checktiff   = false;
    info->tile_info[0]      = '\0';
    info->priori_DEM_tif[0] = '\0';
    info->metafilename[0]   = '\0';
    info->seedDEMsigma      = 0;
    info->check_Matchtag    = false;
    info->threads_num   = 0;
    
    FILE *limage;
    bopened = true;
    char* tmp_chr;
    
    for(int ti=0;ti < args.number_of_images ; ti++)
        info->check_selected_image[ti] = true;
    
    if(args.sensor_type == AB)
        info->IsRA = 0;
    if(args.sensor_provider == PT)
        info->IsRA = false;
    
    if (args.check_arg == 0) // project file open
    {
        if(args.sensor_type == SB) // RFM, pushbroom sensor
        {
            pFile       = fopen(_filename,"r");
            if( pFile == NULL)
            {
                printf("'default.txt' file doesn't exist in SETSM execution folder, please check it!!\n");
                exit(1);
            }
            else
            {
                int count_images = 0;
                        
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
                            for(int ti=0;ti < args.number_of_images ; ti++)
                            {
                                sscanf(bufstr,"Image_RA %hhd %lf %lf\n",&info->IsRA,&info->RA_param[ti][0],&info->RA_param[ti][1]);
                                if (info->IsRA == 1) {
                                    info->RA_param[ti][0] = 0.0;
                                    info->RA_param[ti][1] = 0.0;
                                }
                            }
                        }
                    }
                    else if (!args.check_Threads_num && strstr(bufstr,"Threads_num")!=NULL)
                        sscanf(bufstr,"Threads_num %d\n",&info->threads_num);
                    else if (args.check_arg == 0 && strstr(bufstr,"Image")!=NULL)
                    {
                        sscanf(bufstr,"Image1 %s\n",info->Imagefilename[count_images]);
                        
                        FILE *temp_ptif = fopen(info->Imagefilename[count_images],"r");
                        if(temp_ptif)
                            printf("image %d load completed!\n",count_images);
                        else
                        {
                            printf("image %d load faied. Please check filename!!\n",count_images);
                            exit(1);
                        }
                        
                        if(!args.check_sdm_ortho && (args.check_coreg != 1 && args.check_coreg != 2 && args.check_coreg != 3))
                        {
                            tmp_chr = remove_ext(info->Imagefilename[count_images]);
                            sprintf(info->RPCfilename[count_images],"%s.xml",tmp_chr);
                            
                            FILE *temp_pFile;
                            temp_pFile           = fopen(info->RPCfilename[count_images],"r");
                            //printf("xml file %s\n",info->LeftRPCfilename);
                            if(temp_pFile)
                                printf("image %d xml load completed!\n",count_images);
                            else
                            {
                                sprintf(info->RPCfilename[count_images],"%s.XML",tmp_chr);
                                //printf("xml file %s\n",info->LeftRPCfilename);
                                temp_pFile           = fopen(info->RPCfilename[count_images],"r");
                                if(temp_pFile)
                                    printf("image %d XML load completed!\n",count_images);
                                else
                                {
                                    printf("image %d xml/XML load failed!\n",count_images);
                                    exit(0);
                                }
                            }
                            free(tmp_chr);
                        }
                        count_images++;
                    }
                    else if (args.check_arg == 0 && strstr(bufstr,"outputpath")!=NULL)
                    {
                        sscanf(bufstr,"outputpath %s\n",info->save_filepath);
                        
                        
                        char *Outputpath_name  = SetOutpathName(info->save_filepath);
                        sprintf(info->Outputpath_name,"%s",Outputpath_name);
                        printf("after pathname %s\n",info->Outputpath_name);
                        
                        sprintf(info->tmpdir, "%s/tmp", info->save_filepath);

                        free(Outputpath_name);
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
                
                if(count_images > 1)
                    args.number_of_images = count_images;
                else
                {
                    printf("Please check input images!!\n");
                    exit(1);
                }
            }
            
            
        }
        else  // Collinear Equation, Frame sensor
        {
            printf("Load aerial info\n");
            bopened = OpenDMCproject(args.EO_Path, info, args);
        }
        
        
    }
    else
    {
        if(args.check_RA_line && args.check_RA_sample)
            info->IsRA = 0;
            
        sprintf(info->tmpdir,"%s/tmp",args.Outputpath);
        
        for(int ti= 0 ; ti < args.number_of_images ; ti++)
        {
            sprintf(info->Imagefilename[ti],"%s", args.Image[ti]);
        
            FILE *temp_ptif = fopen(info->Imagefilename[ti],"r");
            if(temp_ptif)
            {
                printf("image%d load completed!\n",ti);
                fclose(temp_ptif);
            }
            else
            {
                printf("image(%d/%d) load faied. Please check filename!! %s\n",ti,args.number_of_images,info->Imagefilename[ti]);
                exit(1);
            }
            
            if(!args.check_sdm_ortho && (args.check_coreg != 1 && args.check_coreg != 2 && args.check_coreg != 3))
            {
                tmp_chr = remove_ext(args.Image[ti]);
                sprintf(info->RPCfilename[ti],"%s.xml",tmp_chr);
                
                if(args.sensor_provider == PT)
                {
                    sprintf(info->RPCfilename[ti],"%s_RPC.TXT",tmp_chr);
                    
                    int chr_size = strlen(tmp_chr);
                    char metadata_file[500];
                    for(int i= 0 ; i < chr_size - 4 ; i++)
                        metadata_file[i] = tmp_chr[i];
                    metadata_file[chr_size - 4] = '\0';
                    
                    //sprintf(info->Imagemetafile[ti],"%s_metadata.xml",metadata_file);
                    sprintf(info->Imagemetafile[ti],"%s_metadata.xml",tmp_chr);
                    printf("imagemetafile %s\n",info->Imagemetafile[ti]);
                }
                                                                              
                FILE *temp_pFile;
                temp_pFile           = fopen(info->RPCfilename[ti],"r");
                printf("xml file %s\n",info->RPCfilename[ti]);
                if(temp_pFile)
                {
                    printf("image%d xml load completed!\n",ti);
                    fclose(temp_pFile);
                }
                else
                {
                    sprintf(info->RPCfilename[ti],"%s.XML",tmp_chr);
                    printf("xml file %s\n",info->RPCfilename[ti]);
                    temp_pFile           = fopen(info->RPCfilename[ti],"r");
                    if(temp_pFile)
                    {
                        printf("image%d XML load completed!\n",ti);
                        fclose(temp_pFile);
                    }
                    else
                    {
                        printf("image%d xml/XML load failed!\n",ti);
                        exit(0);
                    }
                }
                free(tmp_chr);
            }
        }
        
        sprintf(info->save_filepath,"%s",args.Outputpath);
        sprintf(info->Outputpath_name, "%s", args.Outputpath_name);
    }
    
    //default proinfo info loading
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
    
    for(int ti= 0 ; ti < args.number_of_images ; ti++)
    {
        printf("image%d = %s\n",ti,info->Imagefilename[ti]);
        printf("rpc%d = %s\n",ti,info->RPCfilename[ti]);
        
        //check image
        limage = fopen(info->Imagefilename[ti],"r");
        if(!limage)
        {
            bopened     = false;
            printf("Check input image%d filename!!\n",ti);
            exit(1);
        }
        if(limage)
            fclose(limage);
    }
    
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
 
    printf("%d\t%d\t%d\n",args.sensor_type,args.check_fl,args.check_ccd);
    if(args.sensor_type == AB) // Collinear Equation, Frame sensor
    {
        info->frameinfo.Photoinfo = (EO*)calloc(sizeof(EO),info->number_of_images);
        
        if(args.sensor_provider == PT)
        {
            if(args.check_ccd)
            {
                info->frameinfo.m_Camera.m_CCDSize = args.CCD_size;
                bopened = true;
            }
            else
            {
                bopened = false;
                printf("Please input CCD size!!\n");
                exit(1);
            }
        }
        else
        {
            if(args.check_fl && args.check_ccd)
            {
                info->frameinfo.m_Camera.m_focalLength = args.focal_length;
                info->frameinfo.m_Camera.m_CCDSize = args.CCD_size;
                info->frameinfo.m_Camera.m_ppx = 0.0;
                info->frameinfo.m_Camera.m_ppy = 0.0;
                
    //            printf("%f\t%d\t%d\t%f\n",info->frameinfo.m_Camera.m_focalLength,info->frameinfo.m_Camera.m_ImageSize.width,
    //                   info->frameinfo.m_Camera.m_ImageSize.height,info->frameinfo.m_Camera.m_CCDSize);
                
                //info->frameinfo.Photoinfo = (EO*)calloc(sizeof(EO),info->number_of_images);
                /*
                for(int ti = 0 ; ti < info->number_of_images ; ti++)
                {
                    FILE *p_xml;
                    p_xml = fopen(info->RPCfilename[ti],"r");
                    fscanf(p_xml,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&info->frameinfo.Photoinfo[ti].m_Xl,
                           &info->frameinfo.Photoinfo[ti].m_Yl,
                           &info->frameinfo.Photoinfo[ti].m_Zl,
                           &info->frameinfo.Photoinfo[ti].m_Wl,
                           &info->frameinfo.Photoinfo[ti].m_Pl,
                           &info->frameinfo.Photoinfo[ti].m_Kl);
                    fclose(p_xml);

                    sprintf(info->frameinfo.Photoinfo[ti].path,"%s",info->Imagefilename[ti]);
                    printf("%s\n",info->Imagefilename[ti]);
                    printf("%s\n",info->frameinfo.Photoinfo[ti].path);

                    printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                           info->frameinfo.Photoinfo[ti].path,
                           info->frameinfo.Photoinfo[ti].m_Xl,info->frameinfo.Photoinfo[ti].m_Yl,info->frameinfo.Photoinfo[ti].m_Zl,
                           info->frameinfo.Photoinfo[ti].m_Wl,info->frameinfo.Photoinfo[ti].m_Pl,info->frameinfo.Photoinfo[ti].m_Kl);


                }
    */
                bopened = true;
            }
            else
            {
                bopened = false;
                printf("Please input focal length and CCD size!!\n");
                exit(1);
            }
        }
    }
    //exit(1);
    return bopened;
}

bool OpenDMCproject(char* project_path, ProInfo *proinfo, ARGINFO args)
{
    bool ret = true;
    
    FILE *fp;
    char garbage[200];
    char temp_str[500];
    
    if(args.check_EO)
    {
        fp = fopen(project_path, "r");
        
        if(fp)
        {
            while(!feof(fp))
            {
                fscanf(fp, "%s\n", garbage);
                fscanf(fp, "%s\t%lf\n", garbage,&proinfo->frameinfo.m_Camera.m_focalLength);
                fscanf(fp, "%s\t%d\t%d\n", garbage,&proinfo->frameinfo.m_Camera.m_ImageSize.width,&proinfo->frameinfo.m_Camera.m_ImageSize.height);
                fscanf(fp, "%s\t%lf\n", garbage,&proinfo->frameinfo.m_Camera.m_CCDSize);
                proinfo->frameinfo.m_Camera.m_ppx = 0.0;
                proinfo->frameinfo.m_Camera.m_ppy = 0.0;
                
                printf("%f\t%d\t%d\t%f\n",proinfo->frameinfo.m_Camera.m_focalLength,proinfo->frameinfo.m_Camera.m_ImageSize.width,
                       proinfo->frameinfo.m_Camera.m_ImageSize.height,proinfo->frameinfo.m_Camera.m_CCDSize);
                
                fscanf(fp, "%s\n", garbage);
                fscanf(fp, "%s\t%d\t%d\t%d\t%d\n", garbage,&proinfo->frameinfo.NumberofStip,&proinfo->frameinfo.NumberofPhotos,&proinfo->frameinfo.start_stripID,&proinfo->frameinfo.end_stripID);
                printf("%d\t%d\t%d\t%d\n",proinfo->frameinfo.NumberofStip,proinfo->frameinfo.NumberofPhotos,proinfo->frameinfo.start_stripID,proinfo->frameinfo.end_stripID);
                
                fscanf(fp, "%s\n", garbage);
                
                proinfo->frameinfo.Photoinfo = (EO*)calloc(sizeof(EO),proinfo->frameinfo.NumberofPhotos);
                proinfo->number_of_images = proinfo->frameinfo.NumberofPhotos;
                
                int total_photos = 0;
                for(int i=0;i<proinfo->frameinfo.NumberofStip;i++)
                {
                    int image_number;
                    int strip_id;
                    fscanf(fp, "%s\t%d\t%d\n", garbage,&strip_id,&image_number);
                    
                    for(int j=0;j<image_number;j++)
                    {
                        proinfo->frameinfo.Photoinfo[total_photos].strip_ID = strip_id;
                        fscanf(fp, "%d", &proinfo->frameinfo.Photoinfo[total_photos].photo_ID);
                        
                        fscanf(fp, "%s", proinfo->frameinfo.Photoinfo[total_photos].path);
                        
                        
                        fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                               &proinfo->frameinfo.Photoinfo[total_photos].m_Xl,&proinfo->frameinfo.Photoinfo[total_photos].m_Yl,&proinfo->frameinfo.Photoinfo[total_photos].m_Zl,
                               &proinfo->frameinfo.Photoinfo[total_photos].m_Wl,&proinfo->frameinfo.Photoinfo[total_photos].m_Pl,&proinfo->frameinfo.Photoinfo[total_photos].m_Kl);
                        
                        printf("%d\t%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                               proinfo->frameinfo.Photoinfo[total_photos].strip_ID,proinfo->frameinfo.Photoinfo[total_photos].photo_ID,
                               proinfo->frameinfo.Photoinfo[total_photos].path,
                               proinfo->frameinfo.Photoinfo[total_photos].m_Xl,proinfo->frameinfo.Photoinfo[total_photos].m_Yl,proinfo->frameinfo.Photoinfo[total_photos].m_Zl,
                               proinfo->frameinfo.Photoinfo[total_photos].m_Wl,proinfo->frameinfo.Photoinfo[total_photos].m_Pl,proinfo->frameinfo.Photoinfo[total_photos].m_Kl);
                        double o = proinfo->frameinfo.Photoinfo[total_photos].m_Wl;
                        double p = proinfo->frameinfo.Photoinfo[total_photos].m_Pl;
                        double k = proinfo->frameinfo.Photoinfo[total_photos].m_Kl;
                        
                        proinfo->frameinfo.Photoinfo[total_photos].m_Rm = MakeRotationMatrix(o, p, k);
                        
                        total_photos++;
                    }
                }
            }
            fclose(fp);
        }
        else
        {
            ret = false;
        }
    }
    else
    {
        if(args.check_fl && args.check_ccd)
        {
            proinfo->frameinfo.m_Camera.m_focalLength = args.focal_length;
            proinfo->frameinfo.m_Camera.m_CCDSize = args.CCD_size;
            proinfo->frameinfo.m_Camera.m_ppx = 0.0;
            proinfo->frameinfo.m_Camera.m_ppy = 0.0;
            
            printf("%f\t%d\t%d\t%f\n",proinfo->frameinfo.m_Camera.m_focalLength,proinfo->frameinfo.m_Camera.m_ImageSize.width,
                   proinfo->frameinfo.m_Camera.m_ImageSize.height,proinfo->frameinfo.m_Camera.m_CCDSize);
            
            proinfo->frameinfo.Photoinfo = (EO*)calloc(sizeof(EO),proinfo->number_of_images);
            
            for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
            {
                sprintf(proinfo->frameinfo.Photoinfo[ti].path,"%s",proinfo->Imagefilename[ti]);
                printf("%s\n",proinfo->Imagefilename[ti]);
                printf("%s\n",proinfo->frameinfo.Photoinfo[ti].path);
            }
            
            ret = true;
        }
        else
        {
            ret = false;
            printf("Please input focal length and CCD size!!\n");
        }
    }

    
    return ret;
}

void SetPySizes(CSize *data_size_lr, const CSize Subsetsize, const int level)
{
    data_size_lr[0].height       = Subsetsize.height;
    data_size_lr[0].width        = Subsetsize.width;
    for(int i=0;i<level;i++)
    {
        data_size_lr[i+1].width  = data_size_lr[i].width/2;
        data_size_lr[i+1].height = data_size_lr[i].height/2;
    }
}

// Compute the correlation between two arrays using a numerically stable formulation
// Return the correlation coefficient or -99 if undefined
double Correlate(const double *L, const double *R, const int N)
{
    double Lmean = 0;
    double Rmean = 0;
    double rho;
    
    if(N > 0)
    {
        for (int i=0; i<N; i++)
        {
            Lmean += L[i];
            Rmean += R[i];
        }
        Lmean = Lmean / N;
        Rmean = Rmean / N;

        double SumLR = 0;
        double SumL2 = 0;
        double SumR2 = 0;

        for (int i=0; i<N; i++)
        {
            SumLR += (L[i]-Lmean)*(R[i]-Rmean);
            SumL2 += (L[i]-Lmean)*(L[i]-Lmean);
            SumR2 += (R[i]-Rmean)*(R[i]-Rmean);
        }

        if (SumL2 > 1e-8  &&  SumR2 > 1e-8)
        {
            rho = SumLR / (sqrt(SumL2*SumR2));
            //int rI = (int)round(rho*1000);
            //rho = (double)rI / 1000.0;
        }
        else
        {
            rho = (double) -99;
        }
    }
    else
    {
        rho = (double) -99;
    }
 
    return rho;
}

double Correlate(const vector<double> &L, const vector<double> &R, const int N)
{
    double Lmean = 0;
    double Rmean = 0;
    double rho;
    
    if(N > 0)
    {
        for (int i=0; i<N; i++)
        {
            Lmean += L[i];
            Rmean += R[i];
        }
        Lmean = Lmean / N;
        Rmean = Rmean / N;

        double SumLR = 0;
        double SumL2 = 0;
        double SumR2 = 0;

        for (int i=0; i<N; i++)
        {
            SumLR += (L[i]-Lmean)*(R[i]-Rmean);
            SumL2 += (L[i]-Lmean)*(L[i]-Lmean);
            SumR2 += (R[i]-Rmean)*(R[i]-Rmean);
        }

        if (SumL2 > 1e-8  &&  SumR2 > 1e-8)
        {
            rho = SumLR / (sqrt(SumL2*SumR2));
            //int rI = (int)round(rho*1000);
            //rho = (double)rI / 1000.0;
        }
        else
        {
            rho = (double) -99;
        }
    }
    else
    {
        rho = (double) -99;
    }
 
    return rho;
}

void SetTranParam_fromGeoTiff(TransParam *param, char* inputfile)
{
    TIFF *tif;
    GTIF *gtif;
    int cit_length;
    int  size;
    tagtype_t type;
    char *citation;
    geocode_t projCoordTransfCode;
    double projNatOriginLat;
    char ttt[100];
    char hem[100];
    
    printf("param projection %d\n",param->projection);
    
    tif  = XTIFFOpen(inputfile,"r");
    gtif = GTIFNew(tif);
    
    cit_length = GTIFKeyInfo( gtif, GTCitationGeoKey, &size, &type );
    //printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
    if (cit_length > 0)
    {
        //printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
        citation = (char*)malloc(sizeof(char)*cit_length );
        GTIFKeyGet(gtif, GTCitationGeoKey, citation, 0, cit_length);
        
        printf("1 Citation:%s\n",citation);
        sscanf(citation,"%s %s %d, %s %s",ttt,ttt,&param->utm_zone,hem,ttt);
        if(strcmp(hem,"Northern") && strcmp(hem,"Southern"))
        {
            char hemzone[500];
            printf("11 Citation:%s\n",citation);
            sscanf(citation,"%s %s / %s %s %d%s",ttt,ttt,ttt,ttt,&param->utm_zone,hem);
            printf("hemzone %s %d\n",hem,param->utm_zone);
        }
        
        free(citation);
    }
    
    //printf("111 citation %d\t%s\n",param.zone,hem);
    if(!strcmp(hem,"Northern") && cit_length > 0) //utm
    {
        param->bHemisphere = true;
        param->pm = 1;
        param->projection = 2;
    }
    else if(!strcmp(hem,"Southern") && cit_length > 0)
    {
        param->bHemisphere = false;
        param->pm = -1;
        param->projection = 2;
    }
    else if(!strcmp(hem,"N") && cit_length > 0) //utm
    {
        param->bHemisphere = true;
        param->pm = 1;
        param->projection = 2;
    }
    else if(!strcmp(hem,"S") && cit_length > 0)
    {
        param->bHemisphere = false;
        param->pm = -1;
        param->projection = 2;
    }
    else
    {
        cit_length = GTIFKeyInfo( gtif, ProjCoordTransGeoKey, &size, &type );
        //printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
        if (cit_length > 0)
        {
            //printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
            
            GTIFKeyGet(gtif, ProjCoordTransGeoKey, &projCoordTransfCode, 0, cit_length);
            
            //printf("Citation:%d\n",projCoordTransfCode);
            if(projCoordTransfCode == CT_PolarStereographic)
                param->projection = 1;
        }
        
        cit_length = GTIFKeyInfo( gtif, ProjNatOriginLatGeoKey, &size, &type );
        printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
        if (cit_length > 0)
        {
            //printf("length %d\tsize %d\t type %d\n",cit_length,size,type);
            citation = (char*)malloc(sizeof(char)*cit_length );
            
            GTIFKeyGet(gtif, ProjNatOriginLatGeoKey, &projNatOriginLat, 0, cit_length);
            
            printf("2 Citation:%lf\n",projNatOriginLat);
            
            if(projNatOriginLat > 0)
                param->bHemisphere = true;
            else
                param->bHemisphere = false;
            
            free(citation);
        }
    }
    GTIFFree(gtif);
    XTIFFClose(tif);
    
    printf("param %d\t%d\t%d\t%lf\n",param->projection,param->bHemisphere,projCoordTransfCode,projNatOriginLat);
    
    D2DPOINT minXmaxY;
    double t_minXY_X;
    double t_minXY_Y;
    
    double grid_size_c;
    CSize imagesize = ReadGeotiff_info(inputfile,&t_minXY_X,&t_minXY_Y,&grid_size_c);
    
    minXmaxY.m_X = (float)(t_minXY_X + grid_size_c*imagesize.width/2.0);
    minXmaxY.m_Y = (float)(t_minXY_Y - grid_size_c*imagesize.height/2.0);
    
    printf("coord %f\t%f\n",minXmaxY.m_X,minXmaxY.m_Y);
    SetTransParam_param(param,param->bHemisphere);
    D2DPOINT wgs_coord = ps2wgs_single(*param, minXmaxY);
    
    printf("coord %f\t%f\n",wgs_coord.m_X,wgs_coord.m_Y);
    SetTransParam((double)(wgs_coord.m_Y),(double)(wgs_coord.m_X),param);
    
    printf("param %d\t%d\t%d\t%lf\t%d\t%s\n",param->projection,param->bHemisphere,projCoordTransfCode,projNatOriginLat,param->utm_zone,param->direction);
    //exit(1);
}



void SetTransParam_param(TransParam *param, bool Hemisphere)
{
    if(Hemisphere)
    {
        param->a = (double)(6378137.0);
        param->e = (double)(0.08181919);
        param->phi_c = (double)(70.0*DegToRad);
        param->lambda_0 = (double)(-45.0*DegToRad);
    }
    else
    {
        param->a = (double)(6378137.0);
        param->e = (double)(0.08181919);
        param->phi_c = (double)(-71.0*DegToRad);
        param->lambda_0 = 0;
    }
    
    if(param->phi_c < 0)
    {
        param->pm       = -1;
        param->phi_c    = -param->phi_c;
        param->lambda_0 = -param->lambda_0;
    }
    else
        param->pm       = 1;
    
    param->t_c  = tan(PI/4.0 - param->phi_c/2.0)/(pow((1.0-param->e*sin(param->phi_c)) / (1.0+param->e*sin(param->phi_c)) ,(param->e/2.0)));
    param->m_c  = cos(param->phi_c)/sqrt(1.0-pow(param->e,2)*pow(sin(param->phi_c),2));
    
    //UTM param
    param->sa   = 6378137.000000;
    param->sb   = 6356752.314245;
    param->e2   = ( sqrt( ( param->sa*param->sa ) - ( param->sb*param->sb ) ) ) / param->sb;
    param->e2cuadrada = param->e2*param->e2;
    param->c    = (param->sa*param->sa)/param->sb;
}

void SetTransParam(double minLat, double minLon, TransParam *param)
{
    if(param->bHemisphere == 3)
    {
        if(minLat > 0)
            param->bHemisphere = true;
        else
            param->bHemisphere = false;
    }
    
    SetTransParam_param(param,param->bHemisphere);
 
    printf("%f %f %f %f %f %f %f\n",param->sa,param->sb,param->e2,param->e2cuadrada,param->c, minLat, minLon);
    double Lat = minLat;
    char direction[2];
    
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
        param->utm_zone = (int)( ( minLon / 6 ) + 31);
}

void RPCsFree(double **rpcs) {
    for(int i = 0; i < 7; i++) {
        free(rpcs[i]);
    }
    free(rpcs);
}

double** OpenXMLFile(ProInfo *proinfo, int ImageID, double* gsd_r, double* gsd_c, double* gsd, BandInfo* band)
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
    
    pFile           = fopen(proinfo->RPCfilename[ImageID],"r");
    if(pFile)
    {
        if(proinfo->sensor_type == SB) //RPCs info
        {
            out = (double**)calloc(7, sizeof(double*));
            out[0] = (double*)calloc(5, sizeof(double));
            out[1] = (double*)calloc(5, sizeof(double));
            out[6] = (double*)calloc(2, sizeof(double));
            
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
                    band->abscalfactor          = atof(pos2);
                    printf("abscalfactor %f\n",band->abscalfactor);
                }
                
                if(strcmp(token,"<EFFECTIVEBANDWIDTH") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    band->effbw         = atof(pos2);
                    printf("effbw %f\n",band->effbw);
                }
                
                if(strcmp(token,"<TDILEVEL") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    band->tdi           = atof(pos2);
                    printf("tdi %f\n",band->tdi);
                }
                
                if(strcmp(token,"<MEANCOLLECTEDROWGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd_r          = atof(pos2);
                    printf("collect row %f\n",*gsd_r);
                }
                if(strcmp(token,"<MEANCOLLECTEDCOLGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd_c          = atof(pos2);
                    printf("collect col %f\n",*gsd_c);
                }
                if(strcmp(token,"<MEANCOLLECTEDGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd            = atof(pos2);
                    printf("collect gsd %f\n",*gsd);
                }
                
                if(strcmp(token,"<MEANPRODUCTROWGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd_r          = atof(pos2);
                    printf("product row %f\n",*gsd_r);
                }
                if(strcmp(token,"<MEANPRODUCTCOLGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd_c          = atof(pos2);
                    printf("product col %f\n",*gsd_c);
                }
                if(strcmp(token,"<MEANPRODUCTGSD") == 0)
                {
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    *gsd            = atof(pos2);
                    printf("product gsd %f\n",*gsd);
                }
                
                
                if(strcmp(token,"<ERRBIAS") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[6][0]           = atof(pos2);
                    printf("ERRBIAS %f\n",out[6][0]);
                }
                if(strcmp(token,"<ERRRAND") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[6][1]           = atof(pos2);
                }
                if(strcmp(token,"<LINEOFFSET") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][0]           = atof(pos2);
                }
                if(strcmp(token,"<SAMPOFFSET") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][1]           = atof(pos2);
                }
                if(strcmp(token,"<LATOFFSET") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][2]           = atof(pos2);
                }
                if(strcmp(token,"<LONGOFFSET") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][3]           = atof(pos2);
                }
                if(strcmp(token,"<HEIGHTOFFSET") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[0][4]           = atof(pos2);
                }
                
                
                if(strcmp(token,"<LINESCALE") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][0]           = atof(pos2);
                }
                if(strcmp(token,"<SAMPSCALE") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][1]           = atof(pos2);
                }
                if(strcmp(token,"<LATSCALE") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][2]           = atof(pos2);
                }
                if(strcmp(token,"<LONGSCALE") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][3]           = atof(pos2);
                }
                if(strcmp(token,"<HEIGHTSCALE") == 0)
                {
                    //fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[1][4]           = atof(pos2);
                }
             
                
                if(strcmp(token,"<LINENUMCOEFList") == 0)
                {
                    out[2] = (double*)calloc(20, sizeof(double));
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    token2 = strtok(pos2," ");
                    
                    
                    i=0;
                    while(token2 != NULL && i<20)
                    {
                        out[2][i]           = atof(token2);
                        token2 = strtok(NULL," ");
                        
                        //printf("out[2][i] %f\n",out[2][i]);
                        
                        i++;
                        
                        
                    }
                }
                if(strcmp(token,"<LINEDENCOEFList") == 0)
                {
                    out[3] = (double*)calloc(20, sizeof(double));
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    token2 = strtok(pos2," ");
                    i=0;
                    while(token != NULL && i<20)
                    {
                        out[3][i]           = atof(token2);
                        token2 = strtok(NULL," ");
                        
                        //printf("out[3][i] %f\n",out[3][i]);
                        
                        i++;
                    }
                }
                if(strcmp(token,"<SAMPNUMCOEFList") == 0)
                {
                    out[4] = (double*)calloc(20, sizeof(double));
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    token2 = strtok(pos2," ");
                    
                    printf("token2 %s\n",token2);
                    
                    i=0;
                    while(token != NULL && i<20)
                    {
                        out[4][i]           = atof(token2);
                        token2 = strtok(NULL," ");
                        
                        //printf("out[4][i] %f\n",out[4][i]);
                        
                        i++;
                    }
                }
                if(strcmp(token,"<SAMPDENCOEFList") == 0)
                {
                    out[5] = (double*)calloc(20, sizeof(double));
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    token2 = strtok(pos2," ");
                    i=0;
                    while(token != NULL && i<20)
                    {
                        out[5][i]           = atof(token2);
                        token2 = strtok(NULL," ");
                        
                        //printf("out[5][i] %f\n",out[5][i]);
                        
                        i++;
                    }
                }
                
            }
            
            aa                      = out[0][2];
            out[0][2]               = out[0][3];
            out[0][3]               = aa;
            
            aa                      = out[1][2];
            out[1][2]               = out[1][3];
            out[1][3]               = aa;
            
            fclose(pFile);
            
            if(pos1)
                pos1 = NULL;
            if(pos2)
                pos2 = NULL;
            if(token)
                token = NULL;
        }
        else //Exterior orientation
        {
            if(proinfo->sensor_provider == PT)
            {
                printf("eofile\n");
                char eofile[500];
                char* t_str = remove_ext(proinfo->Imagefilename[ImageID]);
                
                sprintf(eofile,"%s_EO.TXT",t_str);
                printf("eofile %s\n",eofile);
                FILE *pfile = fopen(eofile,"r");
                if(pfile)
                {
                    //printf("%d/%d read existing EOs %s\n",ti,proinfo->number_of_images,eofile);
                    ReadEOs(eofile,proinfo->frameinfo.Photoinfo[ImageID],proinfo->frameinfo.m_Camera);
                    
                    printf("EO %f\t%f\t%f\t%f\t%f\t%f\n",
                           proinfo->frameinfo.Photoinfo[ImageID].m_Xl,proinfo->frameinfo.Photoinfo[ImageID].m_Yl,proinfo->frameinfo.Photoinfo[ImageID].m_Zl,
                           proinfo->frameinfo.Photoinfo[ImageID].m_Wl,proinfo->frameinfo.Photoinfo[ImageID].m_Pl,proinfo->frameinfo.Photoinfo[ImageID].m_Kl);
                    
                    printf("Camera %5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",
                           proinfo->frameinfo.m_Camera.m_focalLength,proinfo->frameinfo.m_Camera.m_ppx,proinfo->frameinfo.m_Camera.m_ppy,
                           proinfo->frameinfo.m_Camera.k1,proinfo->frameinfo.m_Camera.k2,proinfo->frameinfo.m_Camera.k3,
                           proinfo->frameinfo.m_Camera.p1,proinfo->frameinfo.m_Camera.p2,
                           proinfo->frameinfo.m_Camera.a1,proinfo->frameinfo.m_Camera.a2);
                    fclose(pfile);
                }
                
            }
            else
            {
                fscanf(pFile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                       &proinfo->frameinfo.Photoinfo[ImageID].m_Xl,&proinfo->frameinfo.Photoinfo[ImageID].m_Yl,&proinfo->frameinfo.Photoinfo[ImageID].m_Zl,
                       &proinfo->frameinfo.Photoinfo[ImageID].m_Wl,&proinfo->frameinfo.Photoinfo[ImageID].m_Pl,&proinfo->frameinfo.Photoinfo[ImageID].m_Kl);
                
                //printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                //       proinfo->frameinfo.Photoinfo[ImageID].path,
                //       proinfo->frameinfo.Photoinfo[ImageID].m_Xl,proinfo->frameinfo.Photoinfo[ImageID].m_Yl,proinfo->frameinfo.Photoinfo[ImageID].m_Zl,
                //       proinfo->frameinfo.Photoinfo[ImageID].m_Wl,proinfo->frameinfo.Photoinfo[ImageID].m_Pl,proinfo->frameinfo.Photoinfo[ImageID].m_Kl);
            }
        
            double o = proinfo->frameinfo.Photoinfo[ImageID].m_Wl;
            double p = proinfo->frameinfo.Photoinfo[ImageID].m_Pl;
            double k = proinfo->frameinfo.Photoinfo[ImageID].m_Kl;
        
            proinfo->frameinfo.Photoinfo[ImageID].m_Rm = MakeRotationMatrix(o, p, k);
            
            *gsd_r = proinfo->frameinfo.m_Camera.m_CCDSize*UMToMM/(proinfo->frameinfo.m_Camera.m_focalLength/proinfo->frameinfo.Photoinfo[ImageID].m_Zl);
            *gsd_c = *gsd_r;
            *gsd = *gsd_r;
                //printf("Rotation %f\t%f\t%f\n",proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m11,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m12,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m13);
                //printf("Rotation %f\t%f\t%f\n",proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m21,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m22,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m23);
                //printf("Rotation %f\t%f\t%f\n",proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m31,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m32,proinfo->frameinfo.Photoinfo[ImageID].m_Rm.m33);
        }
    }
    return out;
}

double** OpenXMLFile_Pleiades(char* _filename)
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
    
    pFile           = fopen(_filename,"r");
    if(pFile)
    {
        out = (double**)calloc(7, sizeof(double*));
        while(!feof(pFile))
        {
            fscanf(pFile,"%s",temp_str);
            if(strcmp(temp_str,"<Inverse_Model>") == 0)
            {
                fgets(linestr,sizeof(linestr),pFile);
                
                printf("sample num\n");
                out[4] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[4][i]           = atof(pos2);
                    
                    printf("%f\n",out[4][i]);
                }
                
                printf("sample den\n");
                out[5] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[5][i]           = atof(pos2);
                    
                    printf("%f\n",out[5][i]);
                }
                
                printf("line num\n");
                out[2] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[2][i]           = atof(pos2);
                    
                    printf("%f\n",out[2][i]);
                }
                
                printf("line den\n");
                out[3] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[3][i]           = atof(pos2);
                    
                    printf("%f\n",out[3][i]);
                }
                
                out[6] = (double*)calloc(2, sizeof(double));
                for(i=0;i<2;i++)
                {
                    
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    out[6][i]           = atof(pos2);
                    
                    printf("%f\n",out[6][i]);
                }
                
                for(i=0;i<14;i++)
                    fgets(linestr,sizeof(linestr),pFile);
                
                out[0] = (double*)calloc(5, sizeof(double)); //offset
                out[1] = (double*)calloc(5, sizeof(double)); //scale
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][2]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][2]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][3]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][3]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][4]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][4]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][1]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][1]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[1][0]           = atof(pos2);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                out[0][0]           = atof(pos2);
                
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

void WriteIRPCs_Planet(char* filename, double **IRPCs)
{
    FILE *pfile = fopen(filename,"w");
    fprintf(pfile,"LINE_OFF: %f\n",IRPCs[0][0]);
    fprintf(pfile,"SAMP_OFF: %f\n",IRPCs[0][1]);
    fprintf(pfile,"LAT_OFF: %f\n",IRPCs[0][3]);
    fprintf(pfile,"LONG_OFF: %f\n",IRPCs[0][2]);
    fprintf(pfile,"HEIGHT_OFF: %f\n",IRPCs[0][4]);
    
    fprintf(pfile,"LINE_SCALE: %f\n",IRPCs[1][0]);
    fprintf(pfile,"SAMP_SCALE: %f\n",IRPCs[1][1]);
    fprintf(pfile,"LAT_SCALE: %f\n",IRPCs[1][3]);
    fprintf(pfile,"LONG_SCALE: %f\n",IRPCs[1][2]);
    fprintf(pfile,"HEIGHT_SCALE: %f\n",IRPCs[1][4]);
    
    for(int i=0;i<20;i++)
        fprintf(pfile,"LINE_NUM_COEFF_%d: %16.15e\n",i+1,IRPCs[2][i]);
    
    for(int i=0;i<20;i++)
        fprintf(pfile,"LINE_DEN_COEFF_%d: %16.15e\n",i+1,IRPCs[3][i]);
    
    for(int i=0;i<20;i++)
        fprintf(pfile,"SAMP_NUM_COEFF_%d: %16.15e\n",i+1,IRPCs[4][i]);
    
    for(int i=0;i<20;i++)
        fprintf(pfile,"SAMP_DEN_COEFF_%d: %16.15e\n",i+1,IRPCs[5][i]);
    fclose(pfile);
}

double** OpenXMLFile_Planet(char* _filename)
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
    
    pFile           = fopen(_filename,"r");
    printf("RPC filename %s\n",_filename);
    if(pFile)
    {
        out = (double**)calloc(7, sizeof(double*));
        //while(!feof(pFile))
        {
            out[0] = (double*)calloc(5, sizeof(double)); //offset
            out[1] = (double*)calloc(5, sizeof(double)); //scale
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[0][0]);
                        
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[0][1]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[0][2]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[0][3]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[0][4]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[1][0]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[1][1]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[1][2]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[1][3]);
            
            fgets(linestr,sizeof(linestr),pFile);
            sscanf(linestr,"%s %lf\n",temp_str,&out[1][4]);
            /*
            for(i=0;i<2;i++)
            {
                for(int j=0;j<5;j++)
                    printf("%lf\n",out[i][j]);
            }
            */
            //exit(1);
            
            //fscanf(pFile,"%s",temp_str);
            //if(strcmp(temp_str,"<Inverse_Model>") == 0)
            {
                //printf("line num\n");
                out[2] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    sscanf(linestr,"%s %lf\n",temp_str,&out[2][i]);
                    //printf("%lf\n",out[2][i]);
                }
                
                //printf("line den\n");
                out[3] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    sscanf(linestr,"%s %lf\n",temp_str,&out[3][i]);
                    
                    //printf("%lf\n",out[3][i]);
                }
                
                //printf("sample num\n");
                out[4] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    sscanf(linestr,"%s %lf\n",temp_str,&out[4][i]);
                    
                    //printf("%lf\n",out[4][i]);
                }
                
                //printf("sample den\n");
                out[5] = (double*)calloc(20, sizeof(double));
                for(i=0;i<20;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    sscanf(linestr,"%s %lf\n",temp_str,&out[5][i]);
                    
                    //printf("%lf\n",out[5][i]);
                }
                
                
                
                out[6] = (double*)calloc(sizeof(double),2);
                //exit(1);
            }
            
            //printf("end read\n");
            
            aa                      = out[0][2];
            out[0][2]               = out[0][3];
            out[0][3]               = aa;
            
            aa                      = out[1][2];
            out[1][2]               = out[1][3];
            out[1][3]               = aa;
        }
        fclose(pFile);
        
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

    // Initializng variables that throw warnings
    Cloud = Nodata;
    MOffnadirangle = Nodata;
    MCrosstrackangle = Nodata;
    MIntrackangle = Nodata;
    
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
                    UL[i]           = atof(pos2);
                    Iinfo->UL[i]    = UL[i];
                    
                }
                //printf("UL %f %f \n",UL[0],UL[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    UR[i]           = atof(pos2);
                    Iinfo->UR[i]    = UR[i];
                }
                //printf("UR %f %f \n",UR[0],UR[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LR[i]           = atof(pos2);
                    Iinfo->LR[i]    = LR[i];
                }
                //printf("LR %f %f \n",LR[0],LR[1]);
                
                for(i=0;i<3;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LL[i]           = atof(pos2);
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
                        MSUNAz          = atof(pos2);
                    }
                    if(strcmp(token,"<MEANSUNEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSUNEl          = atof(pos2);
                    }
                    
                    if(strcmp(token,"<MEANSATAZ") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATAz          = atof(pos2);
                    }
                    if(strcmp(token,"<MEANSATEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATEl          = atof(pos2);
                    }
                    if(strcmp(token,"<MEANINTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MIntrackangle           = atof(pos2);
                    }
                    
                    if(strcmp(token,"<MEANCROSSTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MCrosstrackangle            = atof(pos2);
                    }
                    if(strcmp(token,"<MEANOFFNADIRVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MOffnadirangle          = atof(pos2);
                    }
                    if(strcmp(token,"<CLOUDCOVER") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        Cloud           = atof(pos2);
                        
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

void OpenXMLFile_orientation_planet(char* _filename, ImageInfo *Iinfo)
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
    double temp;
    
    // Initializng variables that throw warnings
    Cloud = Nodata;
    MOffnadirangle = Nodata;
    MCrosstrackangle = Nodata;
    MIntrackangle = Nodata;
    MSATEl = Nodata;
    MSATAz = Nodata;
    MSUNEl = Nodata;
    MSUNAz = Nodata;
    
    
    //printf("%s\n",_filename);
    
    pFile           = fopen(_filename,"r");
    if(pFile)
    {
        bool check_br = false;
        bool check_d = false;
        bool check_i = false;
        bool check_cen = false;
        bool check_cloud = false;
        while(!feof(pFile) && (!check_br || !check_d || !check_i || !check_cen || !check_cloud))
        {
            fscanf(pFile,"%s",temp_str);
            
            if(strcmp(temp_str,"<gml:centerOf>") == 0 && !check_cen)
            {
                check_cen = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                fgets(linestr,sizeof(linestr),pFile);
                pos1 = strstr(linestr,">")+1;
                pos2 = strtok(pos1,"<");
                //printf("pos2 %s\n",pos2);
                
                sscanf(pos2,"%lf %lf",&Iinfo->Center[0],&Iinfo->Center[1]);
                
                
                //printf("Center X Y %lf %lf \n",Iinfo->Center[0],Iinfo->Center[1]);
                //exit(1);
            }
            
            if(strcmp(temp_str,"<ps:geographicLocation>") == 0 && !check_br)
            {
                check_br = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                for(i=0;i<2;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    UL[i]           = atof(pos2);
                    Iinfo->UL[i]    = UL[i];
                    
                }
                temp = Iinfo->UL[0];
                Iinfo->UL[0] = Iinfo->UL[1];
                Iinfo->UL[1] = temp;
                //printf("UL %f %f \n",UL[0],UL[1]);
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                for(i=0;i<2;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    UR[i]           = atof(pos2);
                    Iinfo->UR[i]    = UR[i];
                }
                temp = Iinfo->UR[0];
                Iinfo->UR[0] = Iinfo->UR[1];
                Iinfo->UR[1] = temp;
                //printf("UR %f %f \n",UR[0],UR[1]);
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                for(i=0;i<2;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LR[i]           = atof(pos2);
                    Iinfo->LR[i]    = LR[i];
                }
                temp = Iinfo->LR[0];
                Iinfo->LR[0] = Iinfo->LR[1];
                Iinfo->LR[1] = temp;
                //printf("LR %f %f \n",LR[0],LR[1]);
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                for(i=0;i<2;i++)
                {
                    fgets(linestr,sizeof(linestr),pFile);
                    pos1 = strstr(linestr,">")+1;
                    pos2 = strtok(pos1,"<");
                    LL[i]           = atof(pos2);
                    Iinfo->LL[i]    = LL[i];
                }
                temp = Iinfo->LL[0];
                Iinfo->LL[0] = Iinfo->LL[1];
                Iinfo->LL[1] = temp;
                //printf("LL %f %f \n",LL[0],LL[1]);
            }
            
            if(strcmp(temp_str,"<eop:acquisitionParameters>") == 0 && !check_d)
            {
                check_d = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                bool check_end = false;
                int t_count = 0;
                while(!check_end && t_count < 10)
                {
                    t_count++;
                    fgets(linestr,sizeof(linestr),pFile);
                    strcpy(linestr1,linestr);
                    token1 = strstr(linestr,"<");
                    token = strtok(token1,"=");
                    
                    if(strcmp(token,"<opt:illuminationAzimuthAngle uom") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSUNAz          = atof(pos2);
                    }
                    if(strcmp(token,"<opt:illuminationElevationAngle uom") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSUNEl          = atof(pos2);
                    }
                    
                    if(strcmp(token,"<ps:azimuthAngle uom") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATAz          = atof(pos2);
                    }
                    if(strcmp(token,"<MEANSATEL") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MSATEl          = atof(pos2);
                    }
                    if(strcmp(token,"<MEANINTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MIntrackangle           = atof(pos2);
                    }
                    
                    if(strcmp(token,"<MEANCROSSTRACKVIEWANGLE") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MCrosstrackangle            = atof(pos2);
                    }
                    if(strcmp(token,"<ps:spaceCraftViewAngle uom") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        MOffnadirangle          = atof(pos2);
                    }
                    if(strcmp(token,"<CLOUDCOVER") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        Cloud           = atof(pos2);
                        
                        check_end = true;
                    }
                }
            }
            
            if(strcmp(temp_str,"<eop:product>") == 0 && !check_i)
            {
                check_i = true;
                
                fgets(temp_str,sizeof(temp_str),pFile);
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
                    
                    if(strcmp(token,"<ps:rowGsd") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        Iinfo->GSD.row_GSD    = atof(pos2);
                    }
                    if(strcmp(token,"<ps:columnGsd") == 0)
                    {
                        pos1 = strstr(linestr1,">")+1;
                        pos2 = strtok(pos1,"<");
                        Iinfo->GSD.col_GSD    = atof(pos2);
                        
                        check_end = true;
                    }
                }
            }
            
            if(strcmp(temp_str,"</eop:mask>") == 0 && !check_cloud)
            {
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                fgets(temp_str,sizeof(temp_str),pFile);
                
                fgets(linestr,sizeof(linestr),pFile);
                
                strcpy(linestr1,linestr);
                token1 = strstr(linestr,"<");
                token = strtok(token1,"=");
                
                if(strcmp(token,"<opt:cloudCoverPercentage uom") == 0)
                {
                    check_cloud = true;
                    pos1 = strstr(linestr1,">")+1;
                    pos2 = strtok(pos1,"<");
                    Cloud          = atof(pos2);
                }
            }
        }
        
        fclose(pFile);
        
        Iinfo->Mean_sun_azimuth_angle   = MSUNAz;
        Iinfo->Mean_sun_elevation       = MSUNEl;
        Iinfo->Mean_sat_azimuth_angle_xml   = MSATAz;
        Iinfo->Mean_sat_elevation       = MSATEl;
        Iinfo->Intrack_angle            = MIntrackangle;
        Iinfo->Crosstrack_angle         = MCrosstrackangle;
        Iinfo->Offnadir_angle_xml           = MOffnadirangle;
        Iinfo->cloud                    = Cloud;
        sprintf(Iinfo->imagetime,"%s",imagetime);
        sprintf(Iinfo->SatID,"%s",SatID);
        
        printf("MSUNAz MSUNEl MSATAz MSATEl MOffnadirangle %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",MSUNAz,MSUNEl,MSATAz,MSATEl,MOffnadirangle,Iinfo->GSD.row_GSD,Iinfo->GSD.col_GSD,Iinfo->cloud);
        
        if(pos1)
            pos1 = NULL;
        if(pos2)
            pos2 = NULL;
        if(token)
            token = NULL;
    }
}


void Open_planetmultiinfo(ProInfo *proinfo, char* _filename, ImageInfo *Iinfo)
{
    FILE *pFile;
    pFile           = fopen(_filename,"r");
    /*
    int total_line = 0;
    char t_str[500];
    while(!feof(pFile))
    {
        total_line++;
        fgets(t_str,500,pFile);
    }
    int total_images_end = total_line - 3;
    
    rewind(pFile);
    */
    int total_images;
    
    if(pFile)
    {
        char default_path[500];
        char full_image_path[500];
        fscanf(pFile,"%d",&total_images);
        fscanf(pFile,"%s",default_path);
        printf("image counts %d\tdefault path %s\n",total_images,default_path);
        int selected_count = 0;
        
        for(int i = 0; i<total_images ; i++)
        {
            int strip_ID;
            char temp_str[1000];
            char temp_meta[1000];
            ImageInfo temp_info;
            
            fscanf(pFile,"%s\t%d\t%f\t%f\t%f\n",temp_str,&strip_ID,&temp_info.Offnadir_angle,&temp_info.Mean_sat_azimuth_angle,&temp_info.GSD.pro_GSD);
            
            char temp_year[5];
            char temp_month[3];
            char temp_day[3];
            char temp_hour[3];
            char temp_min[3];
            
            for (int k = 0; k < 4; k++)
                temp_year[k] = temp_str[k];
            temp_year[4] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_month[k] = temp_str[k+4];
            temp_month[2] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_day[k] = temp_str[k+6];
            temp_day[2] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_hour[k] = temp_str[k+9];
            temp_hour[2] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_min[k] = temp_str[k+11];
            temp_min[2] = '\0';
            //printf("file name %s\t%d\t%f\t%f\t%f\n",temp_str,strip_ID,Iinfo[i].Offnadir_angle,Iinfo[i].Mean_sat_azimuth_angle,Iinfo[i].GSD.pro_GSD);
            printf("temp %s\t%s\t%s\t%s\t%s\n",temp_year,temp_month,temp_day,temp_hour,temp_min);
            sprintf(temp_meta,"%s/%s/%s/%s/%d/%s_metadata.xml",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
            printf("temp_meta %s\n",temp_meta);
            
            OpenXMLFile_orientation_planet(temp_meta,&temp_info);
            if(temp_info.cloud <= proinfo->Cloud_th)
            {
                /*
                sprintf(proinfo->Imagefilename[selected_count],"%s/%s/%s/%s/%d/%s_pan.tif",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                sprintf(proinfo->RPCfilename[selected_count],"%s/%s/%s/%s/%d/%s_pan_RPC.TXT",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                sprintf(proinfo->Imagemetafile[selected_count],"%s/%s/%s/%s/%d/%s_metadata.xml",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                */
                
                sprintf(proinfo->Imagefilename[selected_count],"%s/%s/%s/%s/%d/%s.tif",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                sprintf(proinfo->RPCfilename[selected_count],"%s/%s/%s/%s/%d/%s_RPC.TXT",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                sprintf(proinfo->Imagemetafile[selected_count],"%s/%s/%s/%s/%d/%s_metadata.xml",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
                
                Iinfo[selected_count] = temp_info;
                Iinfo[selected_count].Offnadir_angle = temp_info.Offnadir_angle;
                Iinfo[selected_count].Mean_sat_azimuth_angle = temp_info.Mean_sat_azimuth_angle;
                Iinfo[selected_count].GSD.pro_GSD = temp_info.GSD.pro_GSD;
                
                printf("selected_count %d\tpath %s\n%s\n%s\noffnadir %f\tazimuth %f\tGSD %f\t%f\n\n",selected_count,proinfo->Imagefilename[selected_count],proinfo->RPCfilename[selected_count],proinfo->Imagemetafile[selected_count],Iinfo[selected_count].Offnadir_angle,Iinfo[selected_count].Mean_sat_azimuth_angle,Iinfo[selected_count].GSD.pro_GSD,Iinfo[selected_count].cloud);
                
                Iinfo[selected_count].year = atoi(temp_year);
                Iinfo[selected_count].month = atoi(temp_month);
                Iinfo[selected_count].date = atoi(temp_day);
                Iinfo[selected_count].hour = atoi(temp_hour);
                Iinfo[selected_count].min = atoi(temp_min);
                
                sprintf(Iinfo[selected_count].filename,"%s",temp_str);
                sprintf(Iinfo[selected_count].fullpath,"%s/%s/%s/%s/%d",default_path,temp_year,temp_month,temp_day,strip_ID);
                Iinfo[selected_count].strip_ID = strip_ID;
                
                selected_count++;
            }
        }
        
        proinfo->number_of_images = selected_count;
        
        fclose(pFile);
    }
    //printf("total image %d\t%d\n",total_images_end,total_images);
    //exit(1);
}

void SetImageFileInfo(char* _filename, ImageInfo &Iinfo)
{
    char* temp_ext = SetOutpathName(_filename);
    char* temp_str = remove_ext(temp_ext);
    //printf("input str %s\n",temp_str);
    
    char temp_year[5];
    char temp_month[3];
    char temp_day[3];
    char temp_stripID[7];
    
    for (int k = 0; k < 4; k++)
        temp_year[k] = temp_str[k];
    temp_year[4] = '\0';
    
    for(int k=0 ; k < 2 ; k++)
        temp_month[k] = temp_str[k+4];
    temp_month[2] = '\0';
    
    for(int k=0 ; k < 2 ; k++)
        temp_day[k] = temp_str[k+6];
    temp_day[2] = '\0';
    
    for(int k=0 ; k < 6 ; k++)
        temp_stripID[k] = temp_str[k+9];
    temp_stripID[2] = '\0';
    
    //printf("fileinfo %s\t%s\t%s\t%s\n",temp_year,temp_month,temp_day,temp_stripID);
    
    Iinfo.year = atoi(temp_year);
    Iinfo.month = atoi(temp_month);
    Iinfo.date = atoi(temp_day);
    Iinfo.strip_ID = atoi(temp_stripID);
    sprintf(Iinfo.filename,"%s",temp_str);
    
    int filesize;
    char t_str[500];
    memset(t_str,'\0',sizeof(t_str));
    strcpy(t_str,_filename);
    char *fullpath = dirname(t_str);
    sprintf(Iinfo.fullpath,"%s",fullpath);
    
    printf("fileinfo %d\t%d\t%d\t%d\t%s\t%s\n",Iinfo.year,Iinfo.month,Iinfo.date,Iinfo.strip_ID,Iinfo.filename,Iinfo.fullpath);
    //printf("input str %s\n",_filename);
    //exit(1);
}

void Open_planetmultiinfo_args(ARGINFO *args)
{
    FILE *pFile;
    pFile           = fopen(args->Multi_input_file,"r");
    /*
    int total_line = 0;
    char t_str[500];
    while(!feof(pFile))
    {
        total_line++;
        fgets(t_str,500,pFile);
    }
    int total_images_end = total_line - 3;
    
    rewind(pFile);
    */
    int total_images;
    if(pFile)
    {
        char default_path[500];
        char full_image_path[500];
        fscanf(pFile,"%d",&total_images);
        fscanf(pFile,"%s",default_path);
        printf("image counts %d\tdefault path %s\n",total_images,default_path);
        args->number_of_images = total_images;
        for(int i = 0; i<total_images ; i++)
        {
            int strip_ID;
            double temp_value;
            char temp_str[1000];
            fscanf(pFile,"%s\t%d\t%lf\t%lf\t%lf\n",temp_str,&strip_ID,&temp_value,&temp_value,&temp_value);
            
            char temp_year[5];
            char temp_month[3];
            char temp_day[3];
            for (int k = 0; k < 4; k++)
                temp_year[k] = temp_str[k];
            temp_year[4] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_month[k] = temp_str[k+4];
            temp_month[2] = '\0';
            for(int k=0 ; k < 2 ; k++)
                temp_day[k] = temp_str[k+6];
            temp_day[2] = '\0';
            
            //sprintf(args->Image[i],"%s/%s/%s/%s/%d/%s_pan.tif",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
            sprintf(args->Image[i],"%s/%s/%s/%s/%d/%s.tif",default_path,temp_year,temp_month,temp_day,strip_ID,temp_str);
            printf("path %s\n\n",args->Image[i]);
        }
        fclose(pFile);
    }
    //printf("total image %d\t%d\n",total_images_end,total_images);
}

void WriteEOs(char* filename, EO eo, CAMERA_INFO ca)
{
    FILE *pfile = fopen(filename,"w");
    fprintf(pfile,"%f\t%f\t%f\t%f\t%f\t%f\n",eo.m_Xl,eo.m_Yl,eo.m_Zl,eo.m_Wl,eo.m_Pl,eo.m_Kl);
    fprintf(pfile,"%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\n",ca.m_focalLength,ca.m_ppx,ca.m_ppy,
            ca.k1,ca.k2,ca.k3,ca.p1,ca.p2,ca.a1,ca.a2);
    fclose(pfile);
}

void ReadEOs(char* filename, EO &eo, CAMERA_INFO &ca)
{
    FILE *pfile = fopen(filename,"r");
    //EO t_eo;
    fscanf(pfile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&eo.m_Xl,&eo.m_Yl,&eo.m_Zl,&eo.m_Wl,&eo.m_Pl,&eo.m_Kl);
    //eo = t_eo;
    fscanf(pfile,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&ca.m_focalLength,&ca.m_ppx,&ca.m_ppy,&ca.k1,&ca.k2,&ca.k3,&ca.p1,&ca.p2,&ca.a1,&ca.a2);
    fclose(pfile);
}

void CompareRPCs(double **RPCs, double **comRPCs)
{
    for(int i = 0 ; i< 5 ;i++)
    {
        printf("offset %5.4e\n",RPCs[0][i] - comRPCs[0][i]);
        printf("scale %5.4e\n",RPCs[1][i] - comRPCs[1][i]);
    }
    
    for(int i = 0 ; i < 20 ; i++)
    {
        printf("linenumcoef %5.4e\t%5.4e\t%5.4e\t",RPCs[2][i],comRPCs[2][i],RPCs[2][i] - comRPCs[2][i]);
        printf("linedencoef %5.4e\t%5.4e\t%5.4e\t",RPCs[3][i],comRPCs[3][i],RPCs[3][i] - comRPCs[3][i]);
        printf("samplenumcoef %5.4e\t%5.4e\t%5.4e\t",RPCs[4][i],comRPCs[4][i],RPCs[4][i] - comRPCs[4][i]);
        printf("sampledencoef %5.4e\t%5.4e\t%5.4e\n",RPCs[5][i],comRPCs[5][i],RPCs[5][i] - comRPCs[5][i]);
    }
}

EO simulatedEO(EO input_eo, CAMERA_INFO camera, D3DPOINT XYZ_center, EO rotate)
{
    EO out_eo;
    
    D3DPOINT shifted_XYZ;
    
    D3DPOINT dist;
    dist.m_X = input_eo.m_Xl - XYZ_center.m_X;
    dist.m_Y = input_eo.m_Yl - XYZ_center.m_Y;
    dist.m_Z = input_eo.m_Zl - XYZ_center.m_Z;
    
    double ori_distance = sqrt(dist.m_X*dist.m_X + dist.m_Y*dist.m_Y + dist.m_Z*dist.m_Z);
    double scale = ori_distance/camera.m_focalLength;
    //printf("scale ori %f\n",scale);
    D2DPOINT ori_img = GetPhotoCoordinate_single(XYZ_center, input_eo, camera, input_eo.m_Rm);
    //printf("ori_img %f\t%f\n",ori_img.m_X,ori_img.m_Y);
    
    //shifted_XYZ.m_X = tan(rotate.m_Pl*DegToRad)*ori_distance + rotate.m_Xl;
    //shifted_XYZ.m_Y = tan(rotate.m_Wl*DegToRad)*ori_distance + rotate.m_Yl;
    
    shifted_XYZ.m_X = tan(rotate.m_Pl*DegToRad)*ori_distance;
    shifted_XYZ.m_Y = tan(rotate.m_Wl*DegToRad)*ori_distance;
    
    printf("shifted_XYZ %f\t%f\n",shifted_XYZ.m_X,shifted_XYZ.m_Y);
    
    out_eo = input_eo;
    out_eo.m_Xl += shifted_XYZ.m_X;
    out_eo.m_Yl -= shifted_XYZ.m_Y;
    out_eo.m_Wl += rotate.m_Wl;
    out_eo.m_Pl += rotate.m_Pl;
    out_eo.m_Kl += rotate.m_Kl;
    
    RM M = MakeRotationMatrix(out_eo.m_Wl,out_eo.m_Pl,out_eo.m_Kl);
    
    D3DPOINT ray_vector;
    ray_vector.m_X = XYZ_center.m_X - out_eo.m_Xl;
    ray_vector.m_Y = XYZ_center.m_Y - out_eo.m_Yl;
    ray_vector.m_Z = XYZ_center.m_Z - out_eo.m_Zl;
    
    double mag = SQRT(ray_vector);
    ray_vector.m_X /= mag;
    ray_vector.m_Y /= mag;
    ray_vector.m_Z /= mag;
    printf("ray %f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z);
    
    bool check_stop = false;
    //if(rotate.m_Wl == 0 && rotate.m_Pl == 0)
    //    check_stop = true;
    
    int max_iteration = 50;
    int iter = 0;
    double max_correction = camera.m_CCDSize*UMToMM;
    while(!check_stop && iter < max_iteration)
    {
        iter++;
    
        //Xl, Yl adjustment
        D2DPOINT img_center = GetPhotoCoordinate_single(XYZ_center, out_eo, camera, M);
        
        D2DPOINT img_diff;
        img_diff.m_X = ori_img.m_X - img_center.m_X;
        img_diff.m_Y = ori_img.m_Y - img_center.m_Y;
        
        
        printf("iter %d\timg center %f\t%f\t eo %f\t%f\t%f\t%f\t%f\t%f\n",iter,img_center.m_X,img_center.m_Y,
               out_eo.m_Xl,out_eo.m_Yl,out_eo.m_Zl,out_eo.m_Wl,out_eo.m_Pl,out_eo.m_Kl);
        printf("img diff %f\t%f\t%f\n",img_diff.m_X,img_diff.m_Y,max_correction);
        
        if(fabs(img_diff.m_X) < max_correction && fabs(img_diff.m_Y) < max_correction)
        {
            check_stop = true;
   
        }
        else
        {
            shifted_XYZ.m_X = img_diff.m_X*scale;
            shifted_XYZ.m_Y = img_diff.m_Y*scale;
            out_eo.m_Xl -= shifted_XYZ.m_X;
            out_eo.m_Yl -= shifted_XYZ.m_Y;
        }
    }
    
    dist.m_X = out_eo.m_Xl - XYZ_center.m_X;
    dist.m_Y = out_eo.m_Yl - XYZ_center.m_Y;
    dist.m_Z = out_eo.m_Zl - XYZ_center.m_Z;
    
    double shifted_distance = sqrt(dist.m_X*dist.m_X + dist.m_Y*dist.m_Y + dist.m_Z*dist.m_Z);
    double diff_dist = ori_distance - shifted_distance;
    //printf("scale aft %f\n",shifted_distance/camera.m_focalLength);
    
    
    ray_vector.m_X = XYZ_center.m_X - out_eo.m_Xl;
    ray_vector.m_Y = XYZ_center.m_Y - out_eo.m_Yl;
    ray_vector.m_Z = XYZ_center.m_Z - out_eo.m_Zl;
    
    mag = SQRT(ray_vector);
    ray_vector.m_X /= mag;
    ray_vector.m_Y /= mag;
    ray_vector.m_Z /= mag;
    
    printf("ray %f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z);
    
    
    //Zl adjustment
    if(fabs(diff_dist) > 1 )
    {
        bool bcheck = false;
        int max_iter = 50;
        int t_iter = 0 ;
        double adjust_H;
        while(!bcheck && t_iter < max_iter)
        {
            t_iter++;
            
            shifted_distance = sqrt(dist.m_X*dist.m_X + dist.m_Y*dist.m_Y + dist.m_Z*dist.m_Z);
            diff_dist = ori_distance - shifted_distance;
            
            if(fabs(diff_dist) < 1)
                bcheck = true;
            else
            {
                adjust_H = cos(out_eo.m_Pl*DegToRad)*diff_dist;
                
                out_eo.m_Zl += adjust_H;
                double t = (out_eo.m_Zl - XYZ_center.m_Z)/ray_vector.m_Z;
                out_eo.m_Xl = XYZ_center.m_X + t*ray_vector.m_X;
                out_eo.m_Yl = XYZ_center.m_Y + t*ray_vector.m_Y;
                
                dist.m_X = out_eo.m_Xl - XYZ_center.m_X;
                dist.m_Y = out_eo.m_Yl - XYZ_center.m_Y;
                dist.m_Z = out_eo.m_Zl - XYZ_center.m_Z;
            }
            /*
            printf("ray %f\t%f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z,adjust_H);
            printf("t_iter %d\tdiff_dist %f\tXYZ %f\t%f\t%f\n",t_iter,diff_dist, out_eo.m_Xl,out_eo.m_Yl,out_eo.m_Zl);
            printf("scale aft %f\n",shifted_distance/camera.m_focalLength);
             */
        }
    }
    
    ray_vector.m_X = XYZ_center.m_X - out_eo.m_Xl;
    ray_vector.m_Y = XYZ_center.m_Y - out_eo.m_Yl;
    ray_vector.m_Z = XYZ_center.m_Z - out_eo.m_Zl;
    
    mag = SQRT(ray_vector);
    ray_vector.m_X /= mag;
    ray_vector.m_Y /= mag;
    ray_vector.m_Z /= mag;
    
    printf("ray %f\t%f\t%f\n",ray_vector.m_X,ray_vector.m_Y,ray_vector.m_Z);
    
    return out_eo;
}

CSize SetSimulatedImageSize(ImageInfo image_info, double min_H, double max_H, EO &eo, D3DPOINT ori_center, CAMERA_INFO camera, TransParam param, double *objectBR)
{
    CSize SImagesize;
    
    D3DPOINT minXY, maxXY;
    D2DPOINT temp_latlong[4];
    temp_latlong[0].m_X = image_info.UL[0];
    temp_latlong[0].m_Y = image_info.UL[1];
    
    temp_latlong[1].m_X = image_info.UR[0];
    temp_latlong[1].m_Y = image_info.UR[1];
    
    temp_latlong[2].m_X = image_info.LR[0];
    temp_latlong[2].m_Y = image_info.LR[1];
    
    temp_latlong[3].m_X = image_info.LL[0];
    temp_latlong[3].m_Y = image_info.LL[1];
    
    for(int i=0;i<4;i++)
    {
        D2DPOINT temp_XYZ = wgs2ps_single(param,temp_latlong[i]);
        if(i == 0)
        {
            minXY = temp_XYZ;
            maxXY = temp_XYZ;
        }
        else
        {
            if(minXY.m_X > temp_XYZ.m_X)
                minXY.m_X = temp_XYZ.m_X;
            if(minXY.m_Y > temp_XYZ.m_Y)
                minXY.m_Y = temp_XYZ.m_Y;
            
            if(maxXY.m_X < temp_XYZ.m_X)
                maxXY.m_X = temp_XYZ.m_X;
            if(maxXY.m_Y < temp_XYZ.m_Y)
                maxXY.m_Y = temp_XYZ.m_Y;
        }
    }
    
    objectBR[0] = minXY.m_X;
    objectBR[1] = minXY.m_Y;// + Limagesize[ti].height*image_info[ti].GSD.col_GSD;
    objectBR[2] = maxXY.m_X;
    objectBR[3] = maxXY.m_Y;
    
    D3DPOINT PTs[8];
    PTs[0].m_X = objectBR[0];
    PTs[0].m_Y = objectBR[1];
    PTs[0].m_Z = ori_center.m_Z;
    PTs[1].m_X = objectBR[0];
    PTs[1].m_Y = objectBR[3];
    PTs[1].m_Z = ori_center.m_Z;
    PTs[2].m_X = objectBR[2];
    PTs[2].m_Y = objectBR[1];
    PTs[2].m_Z = ori_center.m_Z;
    PTs[3].m_X = objectBR[2];
    PTs[3].m_Y = objectBR[3];
    PTs[3].m_Z = ori_center.m_Z;
    
    /*
    PTs[4].m_X = objectBR[0];
    PTs[4].m_Y = objectBR[1];
    PTs[4].m_Z = max_H;
    PTs[5].m_X = objectBR[0];
    PTs[5].m_Y = objectBR[3];
    PTs[5].m_Z = max_H;
    PTs[6].m_X = objectBR[2];
    PTs[6].m_Y = objectBR[1];
    PTs[6].m_Z = max_H;
    PTs[7].m_X = objectBR[2];
    PTs[7].m_Y = objectBR[3];
    PTs[7].m_Z = max_H;
    */
    RM M = MakeRotationMatrix(eo.m_Wl,eo.m_Pl,eo.m_Kl);
    
    D2DPOINT min_photo, max_photo;
    for(int i=0;i<4;i++)
    {
        D2DPOINT photo  = GetPhotoCoordinate_single(PTs[i],eo,camera,M);
        if(i == 0)
        {
            min_photo = photo;
            max_photo = photo;
        }
        else
        {
            if(min_photo.m_X > photo.m_X)
                min_photo.m_X = photo.m_X;
            if(min_photo.m_Y > photo.m_Y)
                min_photo.m_Y = photo.m_Y;
            
            if(max_photo.m_X < photo.m_X)
                max_photo.m_X = photo.m_X;
            if(max_photo.m_Y < photo.m_Y)
                max_photo.m_Y = photo.m_Y;
        }
    }
    
    /*
    D2DPOINT img_center = GetPhotoCoordinate_single(ori_center, eo, camera, M);
    printf("img center %f\t%f\t eo %f\t%f\t%f\t%f\t%f\t%f\n",img_center.m_X,img_center.m_Y,
           eo.m_Xl,eo.m_Yl,eo.m_Zl,eo.m_Wl,eo.m_Pl,eo.m_Kl);
    
    D2DPOINT XY_center( ori_center.m_X, ori_center.m_Y);
    printf("XY_center %f\t%f\t eo %f\t%f\t%f\t%f\t%f\t%f\n",XY_center.m_X,XY_center.m_Y,
           eo.m_Xl,eo.m_Yl,eo.m_Zl,eo.m_Wl,eo.m_Pl,eo.m_Kl);
    */
    double width = max_photo.m_X - min_photo.m_X;
    double height = max_photo.m_Y - min_photo.m_Y;
    
    SImagesize.width = ceil(width/(camera.m_CCDSize*UMToMM)/2.0)*2;
    SImagesize.height = ceil(height/(camera.m_CCDSize*UMToMM)/2.0)*2;
    /*
    D2DPOINT image(SImagesize.width/2.0, SImagesize.height/2.0);
    D2DPOINT photo = ImageToPhoto_single(image, camera.m_CCDSize, SImagesize);
    printf("photo %f\t%f\n",photo.m_X,photo.m_Y);
    D3DPOINT object = GetObjectCoordinate_single(img_center, ori_center.m_Z, eo, camera, M);
    printf("object %f\t%f\t%f\n",object.m_X,object.m_Y,object.m_Z);
    exit(1);
     */
    return SImagesize;
}

void SimulatedImageGeneration(float *seeddem, CSize seeddem_size, double minX, double maxY, double grid_size, double min_H, double max_H, uint16 *oriimage, CSize imagesize, CSize new_imagesize, char *imagefile, EO ori_eo, EO eo, CAMERA_INFO camera, TransParam param)
{
    uint16 *outimage = NULL;
    
    long total_pixel_count = (long)new_imagesize.width*(long)new_imagesize.height;
    outimage = (uint16*)calloc(sizeof(uint16),total_pixel_count);
    
    printf("minmaxXY %f\t%f\tgrid_size %f\tminmax H %f\t%f\timagesize %d\t%d\n",minX,maxY,grid_size,min_H,max_H,new_imagesize.width,new_imagesize.height);
    RM M = MakeRotationMatrix(eo.m_Wl,eo.m_Pl,eo.m_Kl);
    RM ori_M = MakeRotationMatrix(ori_eo.m_Wl,ori_eo.m_Pl,ori_eo.m_Kl);
    
#pragma omp parallel for schedule(guided)
    for(long pos = 0 ; pos <  total_pixel_count; pos++)
    {
        long row = (long)floor(pos/new_imagesize.width);
        long col = (long)(pos%new_imagesize.width);
        
        D2DPOINT image(col,row);
        D2DPOINT photo = ImageToPhoto_single(image, camera.m_CCDSize, new_imagesize);
        D3DPOINT object;
        //printf("pos col row photo %d\t%d\t%d\t%f\t%f\n",pos,col,row,photo.m_X,photo.m_Y);
        double height;
        double diff_h = 1000000;
        bool check_h = false;
        
        height = min_H;
        while(height >= min_H && height <= max_H && !check_h)
        {
            int max_iteration = 100;
            int iter = 0;
            bool check_stop = false;
            
            while(!check_stop && iter < max_iteration)
            {
                iter++;
                
                object = GetObjectCoordinate_single(photo, height, eo, camera, M);
                //printf("object start %f\t%f\t%f\n",object.m_X,object.m_Y,object.m_Z);
                D3DPOINT dem_pos;
                dem_pos.m_X = (object.m_X - minX)/grid_size;
                dem_pos.m_Y = (maxY - object.m_Y)/grid_size;
                //printf("dem pos %f\t%f\n",dem_pos.m_X,dem_pos.m_Y);
                long dem_index = (long)dem_pos.m_Y * (long)seeddem_size.width + (long)dem_pos.m_X;
                if(dem_pos.m_X >= 0 && dem_pos.m_X < seeddem_size.width && dem_pos.m_Y >= 0 && dem_pos.m_Y < seeddem_size.height)
                {
                    dem_pos.m_Z = seeddem[dem_index];
                    object.m_Z = dem_pos.m_Z;
                    diff_h = fabs(height - dem_pos.m_Z);
                    if(diff_h < 0.1)
                    {
                        check_stop = true;
                    }
                    else
                        height = dem_pos.m_Z;
                    
                    check_h = true;
                }
                else
                {
                    check_stop = true;
                }
            }
            height += 10;
        }
        
        //printf("detected object pos %d\t%d\t%f\t%f\t%f\t%f\n",col,row,object.m_X,object.m_Y,object.m_Z,diff_h);
        //find ori image points from the object intersection of ray
        //printf("object final %f\t%f\t%f\n",object.m_X,object.m_Y,object.m_Z);
        photo  = GetPhotoCoordinate_single(object,ori_eo,camera,ori_M);
        image = PhotoToImage_single(photo, camera.m_CCDSize, imagesize);
        //printf("coord %f\t%f\t%f\t%f\n",photo.m_X,photo.m_Y,image.m_X,image.m_Y);
        
        if(image.m_X >= 0 && image.m_X < imagesize.width && image.m_Y >= 0 && image.m_Y < imagesize.height)
        {
            double t_col       = image.m_X;
            double t_row       = image.m_Y;
            
            long t_col_int   = (long int)(t_col + 0.01);
            long t_row_int   = (long int)(t_row + 0.01);
            
            double dcol        = t_col - t_col_int;
            double drow        = t_row - t_row_int;
            
            if(t_col_int >= 0 && t_col_int +1 < imagesize.width && t_row_int >= 0 && t_row_int +1 < imagesize.height
               && (t_col_int +1) + (t_row_int +1)*(long)imagesize.width < total_pixel_count)
            {
                double value1, value2, value3, value4;
                long index1  = (t_col_int   ) + (t_row_int   )*(long)imagesize.width;
                long index2  = (t_col_int +1) + (t_row_int   )*(long)imagesize.width;
                long index3  = (t_col_int   ) + (t_row_int +1)*(long)imagesize.width;
                long index4  = (t_col_int +1) + (t_row_int +1)*(long)imagesize.width;
                
                value1      = oriimage[index1];
                value2      = oriimage[index2];
                value3      = oriimage[index3];
                value4      = oriimage[index4];
                //printf("values %f\t%f\t%f\t%f\n",value1,value2,value3,value4);
                double value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow) + value3*(1-dcol)*drow + value4*dcol*drow;
                //printf("value %ld\t%f\n",(uint16)value,value);
                outimage[pos] = (uint16)value;
                //exit(1);
            }
        }
    }
    
    WriteGeotiff(imagefile,outimage,new_imagesize.width,new_imagesize.height,1,0,0,param.projection,param.utm_zone,param.bHemisphere,12);
    
    free(outimage);
}

float median(int n, float* x,float min, float max)
{
    float temp;
    int i, j;
    
    long hist_leng = ceil((max - min)/0.01);
    
    long *hist = (long*)calloc(sizeof(long),hist_leng);
    
    for(i=0;i<n;i++)
    {
        long pos = (long)((x[i] - min)*100);
        if(pos >= 0 && pos < hist_leng)
            hist[pos]++;
    }
    
    bool check = false;
    int count = 0;
    long start_pos, current_count;
    long total_count = 0;
    while(!check && count < hist_leng)
    {
        total_count += hist[count];
        
        if(total_count > (int)(n/2.0))
        {
            start_pos = count;
            current_count = total_count;
            check = true;
        }
        count++;
    }
    
    long start_count = current_count - hist[start_pos];
    long remain_count = (int)(n/2.0) - start_count;
    
    float *search_x = (float*)calloc(sizeof(float),hist[start_pos]);
    count = 0;
    for(i=0;i<n;i++)
    {
        long pos = (long)((x[i] - min)*100);
        if(pos == start_pos)
        {
            search_x[count] = x[i];
            count++;
        }
    }
    
    n = hist[start_pos];
    //for(i = 0 ; i< n ; i++)
    //    printf("id %d\t%f\t%d\n",i+start_count,search_x[i],hist[start_pos]);
    
    
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(search_x[j] < search_x[i]) {
                // swap elements
                temp = search_x[i];
                search_x[i] = search_x[j];
                search_x[j] = temp;
                
                //printf("temp xi xj %f\t%f\t%f\n",temp,x[i],x[j]);
            }
        }
    }
    
    float MED = search_x[remain_count];
    
    free(hist);
    free(search_x);
    
    //for(i = 0 ; i< n ; i++)
    //    printf("id %d\t%f\n",i+start_count,search_x[i]);
    
    return MED;
    /*
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        printf("median %f\n",(x[n/2] + x[n/2 - 1]) / 2.0);
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        printf("median %f\n",x[n/2]);
        return x[n/2];
    }
     */
}

float binmedian(int n, float *x)
{
    // Compute the mean and standard deviation
    float sum = 0;
    int i;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    float mu = sum/n;
    
    sum = 0;
    for (i = 0; i < n; i++) {
        sum += (x[i]-mu)*(x[i]-mu);
    }
    float sigma = sqrt(sum/n);
    
    // Bin x across the interval [mu-sigma, mu+sigma]
    int bottomcount = 0;
    int bincounts[1001];
    for (i = 0; i < 1001; i++) {
        bincounts[i] = 0;
    }
    float scalefactor = 1000/(2*sigma);
    float leftend =  mu-sigma;
    float rightend = mu+sigma;
    int bin;
    
    for (i = 0; i < n; i++) {
        if (x[i] < leftend) {
            bottomcount++;
        }
        else if (x[i] < rightend) {
            bin = (int)((x[i]-leftend) * scalefactor);
            bincounts[bin]++;
        }
    }
    
    // If n is odd
    if (n & 1) {
        // Recursive step
        int k, r, count, medbin;
        float oldscalefactor, oldleftend;
        int oldbin;
        float temp;
        
        k = (n+1)/2;
        r = 0;
        
        for (;;) {
            // Find the bin that contains the median, and the order
            // of the median within that bin
            count = bottomcount;
            for (i = 0; i < 1001; i++) {
                count += bincounts[i];
                
                if (count >= k) {
                    medbin = i;
                    k = k - (count-bincounts[i]);
                    break;
                }
            }
            
            bottomcount = 0;
            for (i = 0; i < 1001; i++) {
                bincounts[i] = 0;
            }
            oldscalefactor = scalefactor;
            oldleftend = leftend;
            scalefactor = 1000*oldscalefactor;
            leftend = medbin/oldscalefactor + oldleftend;
            rightend = (medbin+1)/oldscalefactor + oldleftend;
            
            // Determine which points map to medbin, and put
            // them in spots r,...n-1
            i = r; r = n;
            while (i < r) {
                oldbin = (int)((x[i]-oldleftend) * oldscalefactor);
                
                if (oldbin == medbin) {
                    r--;
                    SWAP(x[i],x[r]);
                    
                    // Re-bin on a finer scale
                    if (x[i] < leftend) {
                        bottomcount++;
                    }
                    else if (x[i] < rightend) {
                        bin = (int)((x[i]-leftend) * scalefactor);
                        bincounts[bin]++;
                    }
                }
                else {
                    i++;
                }
            }
            
            // Stop if all points in medbin are the same
            int samepoints = 1;
            for (i = r+1; i < n; i++) {
                if (x[i] != x[r]) {
                    samepoints = 0;
                    break;
                }
            }
            if (samepoints) {
                return x[r];
            }
            
            // Stop if there's <= 20 points left
            if (n-r <= 20) {
                break;
            }
        }
        
        // Perform insertion sort on the remaining points,
        // and then pick the kth smallest
        float a;
        int j;
        for (i = r+1; i < n; i++) {
            a = x[i];
            for (j = i-1; j >= r; j--) {
                if (x[j] > a) {
                    break;
                }
                x[j+1] = x[j];
            }
            x[j+1] = a;
        }
        
        return x[r-1+k];
    }
    
    // If n is even (not implemented yet)
    else {
        return 0;
    }
}

double quickselect(double *arr, int n, int k)
{
    unsigned long i,ir,j,l,mid;
    double a,temp;
    
    l=0;
    ir=n-1;
    for(;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir] < arr[l]) {
                SWAP(arr[l],arr[ir]);
            }
            return arr[k];
        }
        else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l],arr[l+1]);
            }
            i=l+1;
            j=ir;
            a=arr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
            }
            arr[l+1]=arr[j];
            arr[j]=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}

double quickselect(vector<double> &arr, int n, int k)
{
    unsigned long i,ir,j,l,mid;
    double a,temp;
    
    l=0;
    ir=n-1;
    for(;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir] < arr[l]) {
                SWAP(arr[l],arr[ir]);
            }
            return arr[k];
        }
        else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l],arr[l+1]);
            }
            i=l+1;
            j=ir;
            a=arr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i],arr[j]);
            }
            arr[l+1]=arr[j];
            arr[j]=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}

float quickselect(vector<PairCA> &arr, int n, int k)
{
    int i, j;
    for (i = 0; i < n-1; i++)
    {
        // Last i elements are already in place
        for (j = 0; j < n-i-1; j++)
        {
           if (arr[j].CA < arr[j+1].CA)
           {
              //swap(&arr[j], &arr[j+1]);
               
               PairCA temp = arr[j];
               arr[j] = arr[j+1];
               arr[j+1] = temp;
           }
        }
    }
    return 0;
    /*
    unsigned long i,ir,j,l,mid;
    double a;
    PairCA temp(0,0);
    l=0;
    ir=n-1;
    for(;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir].CA < arr[l].CA) {
                temp = arr[l];
                arr[l] = arr[ir];
                arr[ir] = temp;
            }
            return arr[k].CA;
        }
        else {
            mid=(l+ir) >> 1;
            temp = arr[mid];
            arr[mid] = arr[l+1];
            arr[l+1] = temp;
            
            if (arr[l].CA > arr[ir].CA) {
                temp = arr[l];
                arr[l] = arr[ir];
                arr[ir] = temp;
            }
            if (arr[l+1].CA > arr[ir].CA) {
                temp = arr[l+1];
                arr[l+1] = arr[ir];
                arr[ir] = temp;
            }
            if (arr[l].CA > arr[l+1].CA) {
                temp = arr[l];
                arr[l] = arr[l+1];
                arr[l+1] = temp;
            }
            i=l+1;
            j=ir;
            a=arr[l+1].CA;
            for (;;) {
                do i++; while (arr[i].CA < a);
                do j--; while (arr[j].CA > a);
                if (j < i) break;
                temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
            arr[l+1].CA=arr[j].CA;
            arr[j].CA=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
     */
}

double cal_var(vector<double> &arr)
{
    int arr_count = arr.size();
    double sum = 0;
    for(int i = 0 ; i < arr_count ; i++)
    {
        sum += arr[i];
        //printf("arr[i] %f\n",arr[i]);
    }
    double mean = sum/arr_count;
    //printf("mean %f\n",mean);
    double sum_var = 0;
    for(int i = 0 ; i < arr_count ; i++)
    {
        sum_var += (arr[i] - mean)*(arr[i] - mean);
    }
    return sqrt(sum_var/arr_count);
}

void GetOverlappedImages_LatLong(const double *tile_boundary, const int total_images_count, const ImageInfo *imageinfo, vector<int> &selected_ID)
{
    for(int count = 0 ; count < total_images_count ; count++)
    {
        D2DPOINT tile_lt(tile_boundary[0], tile_boundary[3]);
        D2DPOINT tile_rb(tile_boundary[2], tile_boundary[1]);
        D2DPOINT Image_lt(imageinfo[count].UL[0],imageinfo[count].UL[1]);
        D2DPOINT Image_rb(imageinfo[count].LR[0],imageinfo[count].LR[1]);
        bool check_overlap = CheckOverlap(tile_lt, tile_rb, Image_lt, Image_rb);
        
        if(check_overlap)
            selected_ID.push_back(count);
    }
}

void GetOverlappedImages_XY(const double *tile_boundary_XY, const TransParam param, const int total_images_count, const ImageInfo *imageinfo, vector<int> &selected_ID)
{
    
    for(int count = 0 ; count < total_images_count ; count++)
    {
        D2DPOINT tile_lt(tile_boundary_XY[0], tile_boundary_XY[3]);
        D2DPOINT tile_rb(tile_boundary_XY[2], tile_boundary_XY[1]);
        D2DPOINT Image_lt(imageinfo[count].UL[0],imageinfo[count].UL[1]);
        D2DPOINT Image_rb(imageinfo[count].LR[0],imageinfo[count].LR[1]);
        
        D2DPOINT Image_lt_xy = wgs2ps_single(param,Image_lt);
        D2DPOINT Image_rb_xy = wgs2ps_single(param,Image_rb);
        bool check_overlap = CheckOverlap(tile_lt, tile_rb, Image_lt_xy, Image_rb_xy);
        
        if(check_overlap)
            selected_ID.push_back(count);
    }
}

bool CheckOverlap(const D2DPOINT br1_lt, const D2DPOINT br1_rb, const D2DPOINT br2_lt, const D2DPOINT br2_rb)
{
    if (br1_lt.m_X > br2_rb.m_X || br2_lt.m_X > br1_rb.m_X)
        return false;
    
    if (br1_lt.m_Y < br2_rb.m_Y || br2_lt.m_Y < br1_rb.m_Y)
        return false;
    
    return true;
}


void MakeSobelMagnitudeImage(const CSize _img_size, const uint16* _src_image, uint16* _dist_mag_image, int16* _dir)
{
    int hy[3][3] , hx[3][3];
    
    //sobel mask
    hx[0][0]=-1;               //  -1  0  1
    hx[0][1]=0;                //  -2  0  2
    hx[0][2]=1;                //  -1  0  1
    hx[1][0]=-2;
    hx[1][1]=0;
    hx[1][2]=2;
    hx[2][0]=-1;
    hx[2][1]=0;
    hx[2][2]=1;

    hy[0][0]=1;                //  1  2  1
    hy[0][1]=2;                //  0  0  0
    hy[0][2]=1;                // -1 -2 -1
    hy[1][0]=0;
    hy[1][1]=0;
    hy[1][2]=0;
    hy[2][0]=-1;
    hy[2][1]=-2;
    hy[2][2]=-1;

#pragma omp parallel for
    for(long int i=0;i<_img_size.height;i++)
    {
        for(long int j=0;j<_img_size.width;j++)
        {
            if(j == 0 || j == _img_size.width-1 || i == 0 || i == _img_size.height-1)
            {
                _dist_mag_image[i*_img_size.width + j] = 0;
                _dir[i*_img_size.width + j] = 0;
            }
            else
            {
                double temp_ptr_X=((hx[0][0]*_src_image[(i-1)*_img_size.width + (j-1)]+hx[0][1]*_src_image[(i-1)*_img_size.width + j]+hx[0][2]*_src_image[(i-1)*_img_size.width + (j+1)]+
                                    hx[1][0]*_src_image[    i*_img_size.width + (j-1)]+hx[1][1]*_src_image[    i*_img_size.width + j]+hx[1][2]*_src_image[    i*_img_size.width + (j+1)]+
                                    hx[2][0]*_src_image[(i+1)*_img_size.width + (j-1)]+hx[2][1]*_src_image[(i+1)*_img_size.width + j]+hx[2][2]*_src_image[(i+1)*_img_size.width + (j+1)]));

                double temp_ptr_Y=((hy[0][0]*_src_image[(i-1)*_img_size.width + (j-1)]+hy[0][1]*_src_image[(i-1)*_img_size.width + j]+hy[0][2]*_src_image[(i-1)*_img_size.width + (j+1)]+
                                    hy[1][0]*_src_image[    i*_img_size.width + (j-1)]+hy[1][1]*_src_image[    i*_img_size.width + j]+hy[1][2]*_src_image[    i*_img_size.width + (j+1)]+
                                    hy[2][0]*_src_image[(i+1)*_img_size.width + (j-1)]+hy[2][1]*_src_image[(i+1)*_img_size.width + j]+hy[2][2]*_src_image[(i+1)*_img_size.width + (j+1)]));

                double temp_ptr=sqrt(temp_ptr_X*temp_ptr_X + temp_ptr_Y*temp_ptr_Y);

                _dir[i*_img_size.width + j]= (int16)(floor(atan2(temp_ptr_Y,temp_ptr_X)*RadToDeg));
                _dist_mag_image[i*_img_size.width + j] = (uint16)(temp_ptr+0.5);
            }
        }
    }
}


void Orientation(const CSize imagesize, const uint16* Gmag, const int16* Gdir, const uint8 Template_size, uint8* plhs)
{
    //Procedure
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int Half_template_size = Template_size / 2;

    long int main_col = imagesize.width;
    long int main_row = imagesize.height;

    //malloc gu_weight_pre_computed matrix
    int size = 2 * Half_template_size - 1;
    double **gu_weight_pre_computed = (double **) malloc(size * sizeof(double *));
    for (int row = -Half_template_size + 1; row <= Half_template_size - 1; row++)
        gu_weight_pre_computed[Half_template_size - 1 + row] = (double *) malloc(size * sizeof(double));

    // calculate gu_weight_pre_computed values
    for (int row = -Half_template_size + 1; row <= Half_template_size - 1; row++)
    {
        for (int col = -Half_template_size + 1; col <= Half_template_size - 1; col++)
        {
            //double angle_sigma = 1.5;
            gu_weight_pre_computed[Half_template_size - 1 + row][Half_template_size - 1 + col] = exp(-(double) (row * row + col * col) / 4.5);
        }
    }

    //const int bin_count = 18;
    //const double sub_ratio = 4.6; //(double)bin_count / 4.0;
    //const double bin_angle = 20.0; //360.0 / (double)bin_count;
#pragma omp parallel for collapse(2)
    for (long int mask_row = 0; mask_row < main_row; mask_row++)
    {
        for (long int mask_col = 0; mask_col < main_col; mask_col++)
        {
            double bin[18] = {0.0};
            double sub_bin[4] = {0.0};

            //calulation of major direction, refering to SIFT descriptor based on gradient
            for (int row = max<long>(-mask_row, -Half_template_size + 1); row <= min<long>(main_row-mask_row-1, Half_template_size - 1); row++)
            {
                int dcol = numcols[abs(row)];
                for (int col = max<long>(-mask_col, -dcol); col <= min<long>(dcol, main_col-mask_col-1); col++)
                {
                    double gu_weight = gu_weight_pre_computed[Half_template_size - 1 + row][Half_template_size - 1 + col];
                    //long int radius2 = (row * row + col * col);
                    long int pixel_row = mask_row + row;
                    long int pixel_col = mask_col + col;
                    if (/*radius2 <= (Half_template_size - 1) * (Half_template_size - 1) &&*/ pixel_row > 0 && pixel_row < main_row - 1 && pixel_col > 0 && pixel_col < main_col - 1)
                    // Possible OOB read in here
                    {
                        double mag = Gmag[pixel_row * main_col + pixel_col];
                        // overflow(?)
                        double theta = (double) (Gdir[pixel_row * main_col + pixel_col]);

                        if (theta < 0)
                            theta += 360;
                        else if(theta > 360)
                            theta -= 360;

                        long int index = (long int) (theta / 20.0);
                        if (index >= 0 && index < 18)
                        {
                            bin[index] += mag * gu_weight;
                            sub_bin[(int) (index / SUB_RATIO)] += mag * gu_weight;
                        }
                    }
                }
            }

            double max_th = -100;
            int sub_th = 0, th = 0;
            for (int i = 0; i < 4; i++)
            {
                if (sub_bin[i] > max_th)
                {
                    max_th = sub_bin[i];
                    sub_th = i;
                }
            }

            max_th = -100;
            for (int i = (int) (sub_th * SUB_RATIO); i < (int) (sub_th * SUB_RATIO + SUB_RATIO); i++)
            {
                if (bin[i] > max_th)
                {
                    max_th = bin[i];
                    th = i;
                }
            }
            plhs[mask_row * main_col + mask_col] = (uint8) th;
        }
    }
    // free all of gu_weight_pre_computed
    for (int j = 0; j < size; j++)
    {
        if (gu_weight_pre_computed[j])
            free(gu_weight_pre_computed[j]);
    }

    if (gu_weight_pre_computed)
        free(gu_weight_pre_computed);
}



GMA_double* GMA_double_create(uint32 size_row, uint32 size_col)
{
    long int cnt;
    GMA_double *out;
    out=(GMA_double*)malloc(sizeof(GMA_double));
    out->nrows=size_row;
    out->ncols=size_col;
    out->val=(long double**)calloc(sizeof(long double*),size_row);
    //out->val[0]=(float*)malloc(sizeof(float)*size_row*size_col);
    out->data=(long double*)calloc(sizeof(long double),size_row*size_col);
    for(cnt=0;cnt<size_row;cnt++)
    {
        //out->val[cnt]=out->val[0]+sizeof(float*)*size_col*cnt;
        out->val[cnt]=&(out->data[cnt*size_col]);
    }
    return out;
}

void GMA_double_destroy(GMA_double* in)
{
    /*
     int32_t cnt;
     for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
     free(in);
     */
    free(in->val);
    free(in->data);
    free(in);
}


void GMA_double_inv(GMA_double *a, GMA_double *I)
{
    long int cnt1,cnt2,cnt3;
    long double pivot,coeff;
    
    GMA_double *b=GMA_double_create(a->nrows,a->ncols);
    //printf("duplicate the matrix a\n");
    
    for(cnt1=0;cnt1<a->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            b->val[cnt1][cnt2]=a->val[cnt1][cnt2];
        }
    }
    
    //printf("initialize I\n");
    for(cnt1=0;cnt1<b->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<b->nrows;cnt2++)
        {
            I->val[cnt1][cnt2]=(cnt1==cnt2)?1:0;
        }
    }
    
    //printf("Gaussian elimation - forward\n");
    for(cnt1=0;cnt1<b->nrows-1;cnt1++)
    {
        pivot=b->val[cnt1][cnt1];
        for(cnt2=cnt1+1;cnt2<b->ncols;cnt2++)
        {
            coeff=b->val[cnt2][cnt1]/pivot;
            //printf("%d %d\tcoeff %f\n",cnt1, cnt2, coeff);
            for(cnt3=0;cnt3<b->nrows;cnt3++)
            {
                b->val[cnt2][cnt3]-=b->val[cnt1][cnt3]*coeff;
                //printf("%f\n",b->val[cnt2][cnt3]);
                I->val[cnt2][cnt3]-=I->val[cnt1][cnt3]*coeff;
                //printf("%f\t%f\n",b->val[cnt2][cnt3], I->val[cnt2][cnt3]);
            }
        }
    }
    /*
     printf("b matrix\n");
     for(cnt1=0;cnt1<4;cnt1++)
     {
     for(cnt2=0;cnt2<4;cnt2++)
     {
     printf("%f\t",I->val[cnt1][cnt2]);
     }
     printf("\n");
     }
     */
    
    //printf("Backward Elimination\n");
    for(cnt1=b->nrows-1;cnt1>=0;cnt1--)
    {
        pivot=b->val[cnt1][cnt1];
        for(cnt2=cnt1-1;cnt2>=0;cnt2--)
        {
            //printf("%d %d\n",cnt1, cnt2);
            coeff=b->val[cnt2][cnt1]/pivot;
            
            //printf("%d %d\tcoeff %f\n",cnt1, cnt2, coeff);
            for(cnt3=b->nrows-1;cnt3>=0;cnt3--)
            {
                //printf("%d\t%d\t%d\n",cnt1,cnt2,cnt3);
                b->val[cnt2][cnt3]-=b->val[cnt1][cnt3]*coeff;
                //printf("%f\n",b->val[cnt2][cnt3]);
                I->val[cnt2][cnt3]-=I->val[cnt1][cnt3]*coeff;
                
                //printf("%f\t%f\n",b->val[cnt2][cnt3], I->val[cnt2][cnt3]);
            }
            //printf("enf for loop\n");
        }
    }
    
    //printf("scaling\n");
    for(cnt1=0;cnt1<b->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<b->nrows;cnt2++)
        {
            I->val[cnt1][cnt2]/=b->val[cnt1][cnt1];
        }
        
    }
    
    //printf("release memory\n");
    GMA_double_destroy(b);
    
}

void GMA_double_mul(GMA_double *a, GMA_double *b, GMA_double *out)
{
    long int cnt1,cnt2,cnt3;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    {
        for(cnt2=0;cnt2<b->ncols;cnt2++)
        {
            out->val[cnt1][cnt2]=0;
            for(cnt3=0;cnt3<a->ncols;cnt3++)
            {
                out->val[cnt1][cnt2]+=a->val[cnt1][cnt3]*b->val[cnt3][cnt2];
            }
        }
    }
}

void GMA_double_Tran(GMA_double *a, GMA_double *out)
{
    long int cnt1,cnt2,cnt3;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            out->val[cnt2][cnt1] = 0;
            out->val[cnt2][cnt1] = a->val[cnt1][cnt2];
        }
    }
}

void GMA_double_sub(GMA_double *a, GMA_double *b, GMA_double *out)
{
    long int cnt1,cnt2;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            out->val[cnt1][cnt2]=0;
            out->val[cnt1][cnt2] = a->val[cnt1][cnt2] - b->val[cnt1][cnt2];
        }
    }
}

void GMA_double_sum(GMA_double *a, GMA_double *b, GMA_double *out)
{
    long int cnt1,cnt2;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            out->val[cnt1][cnt2]=0;
            out->val[cnt1][cnt2] = a->val[cnt1][cnt2] + b->val[cnt1][cnt2];
        }
    }
}


void GMA_double_printf(GMA_double *a)
{
    long int cnt1,cnt2;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            printf("[%ld,%ld] %e\t",cnt1,cnt2,a->val[cnt1][cnt2]);
        }
        printf("\n");
    }
}

//Returns created triangulation pointer
FullTriangulation *TINCreate_list(D3DPOINT *ptslists, long int numofpts, vector<UI3DPOINT> *trilists, double min_max[], long int *count_tri, double resolution)
{
    if (numofpts <= 2) {
        *count_tri = 0;
        return NULL;
    }
    
    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];
    
    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //Check to ensure grid width/height fit in 16 bit integers
    if(width > 32767 || height > 32767)
    {
        printf("ERROR: Grid is too large. width: %d height: %d\n", width, height);
        exit(1);
    }
    //printf("\tTINCreate: PTS = %d, width = %d, height = %d, resolution = %f\n", numofpts, width, height, resolution);
    
    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    GridPoint *grid_points = new GridPoint[numofpts];
    
    //printf("done s1\n");
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        grid_points[t].col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        grid_points[t].row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[grid_points[t].row * width + grid_points[t].col] = t;
    }
    
    //printf("done s2\n");
    GridPoint **points_ptrs = new GridPoint*[numofpts];
#pragma omp parallel for
    for (std::size_t t = 0; t < numofpts; ++t) points_ptrs[t] = grid_points + t;
    
    FullTriangulation *triangulation = new FullTriangulation(width, height);
    //double begin = omp_get_wtime();
    //printf("done s3\n");
    triangulation->Triangulate(points_ptrs, numofpts);
    //printf("done triangulation\n");
    //double end = omp_get_wtime();
    //printf("Triangulate took %lf with %d points\n", end - begin, numofpts);
    
    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    //printf("s3-1\n");
    *count_tri = (int)(triangulation->GetAllTris(&tris));
    //printf("count tri = %d\n",*count_tri);
    for(long t = 0 ; t < tris.size() ; t++)
    {
        int row, col;
        UI3DPOINT temp_pt;
        
        row = tris[t].pts[0].row;
        col = tris[t].pts[0].col;
        temp_pt.m_X = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[1].row;
        col = tris[t].pts[1].col;
        temp_pt.m_Y = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[2].row;
        col = tris[t].pts[2].col;
        temp_pt.m_Z = index_in_ptslists[row * width + col];
        
        trilists->push_back(temp_pt);
    }
    //printf("s4\n");
    //delete [] tris;
    delete [] points_ptrs;
    delete [] grid_points;
    
    //printf("s5\n");
    return triangulation;
}

FullTriangulation *TINCreate_list_vector(vector<D3DPOINT> &ptslists, long int numofpts, vector<UI3DPOINT> &trilists, double min_max[], long int *count_tri, double resolution)
{
    if (numofpts <= 2) {
        *count_tri = 0;
        return NULL;
    }
    
    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];
    
    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //Check to ensure grid width/height fit in 16 bit integers
    if(width > 32767 || height > 32767)
    {
        printf("ERROR: Grid is too large. width: %d height: %d\n", width, height);
        exit(1);
    }
    //printf("\tTINCreate: PTS = %d, width = %d, height = %d, resolution = %f\n", numofpts, width, height, resolution);
    
    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    GridPoint *grid_points = new GridPoint[numofpts];
    
    //printf("done s1\n");
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        grid_points[t].col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        grid_points[t].row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[grid_points[t].row * width + grid_points[t].col] = t;
    }
    
    //printf("done s2\n");
    GridPoint **points_ptrs = new GridPoint*[numofpts];
#pragma omp parallel for
    for (std::size_t t = 0; t < numofpts; ++t) points_ptrs[t] = grid_points + t;
    
    FullTriangulation *triangulation = new FullTriangulation(width, height);
    //double begin = omp_get_wtime();
    //printf("done s3\n");
    triangulation->Triangulate(points_ptrs, numofpts);
    //printf("done triangulation\n");
    //double end = omp_get_wtime();
    //printf("Triangulate took %lf with %d points\n", end - begin, numofpts);
    
    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    //printf("s3-1\n");
    *count_tri = (int)(triangulation->GetAllTris(&tris));
    //printf("count tri = %d\n",*count_tri);
    for(long t = 0 ; t < tris.size() ; t++)
    {
        int row, col;
        UI3DPOINT temp_pt;
        
        row = tris[t].pts[0].row;
        col = tris[t].pts[0].col;
        temp_pt.m_X = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[1].row;
        col = tris[t].pts[1].col;
        temp_pt.m_Y = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[2].row;
        col = tris[t].pts[2].col;
        temp_pt.m_Z = index_in_ptslists[row * width + col];
        
        trilists.push_back(temp_pt);
    }
    //printf("s4\n");
    //delete [] tris;
    delete [] points_ptrs;
    delete [] grid_points;
    
    //printf("s5\n");
    return triangulation;
}

void TINUpdate_list(D3DPOINT *ptslists, long int numofpts, vector<UI3DPOINT> *trilists, double min_max[], long int *count_tri, double resolution, FullTriangulation *oldTri, D3DPOINT *blunderlist, long int numblunders)
{
    
    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];
    
    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //printf("\tTINUpdate: PTS = %d, blunders = %d, width = %d, height = %d, resolution = %f\n", numofpts, numblunders, width, height, resolution);
    
    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        //index_in_ptslists[((ptslists[t].m_Y - minY_ptslists)/resolution)*width+((ptslists[t].m_X - minX_ptslists) / resolution)] = t;
        INDEX col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        INDEX row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[row * width + col] = t;
    }
    
    GridPoint *grid_blunders = new GridPoint[numblunders];
    for (std::size_t t = 0; t < numblunders; ++t)
    {
        grid_blunders[t].col = 0.5 + (blunderlist[t].m_X - minX_ptslists) / resolution;
        grid_blunders[t].row = 0.5 + (blunderlist[t].m_Y - minY_ptslists) / resolution;
    }
    
    GridPoint **blunder_ptrs = new GridPoint*[numblunders];
#pragma omp parallel for
    for (std::size_t t = 0; t < numblunders; ++t) blunder_ptrs[t] = grid_blunders + t;
    
    //double begin = omp_get_wtime();
    
    oldTri->Retriangulate(blunder_ptrs, numblunders);
    
    //double end = omp_get_wtime();
    //printf("Retriangulate took %lf with %d points, %d blunders\n", end - begin, numofpts, numblunders);
    
    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    *count_tri = (int)(oldTri->GetAllTris(&tris));
    for(long t = 0 ; t < tris.size() ; t++)
    {
        int row, col;
        UI3DPOINT temp_pt;
        
        row = tris[t].pts[0].row;
        col = tris[t].pts[0].col;
        
        temp_pt.m_X = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[1].row;
        col = tris[t].pts[1].col;
        temp_pt.m_Y = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[2].row;
        col = tris[t].pts[2].col;
        temp_pt.m_Z = index_in_ptslists[row * width + col];
        
        trilists->push_back(temp_pt);
    }
    delete [] blunder_ptrs;
    delete [] grid_blunders;
}

void TINUpdate_list_vector(vector<D3DPOINT> &ptslists, long int numofpts, vector<UI3DPOINT> &trilists, double min_max[], long int *count_tri, double resolution, FullTriangulation *oldTri, vector<D3DPOINT> &blunderlist, long int numblunders)
{
    
    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];
    
    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //printf("\tTINUpdate: PTS = %d, blunders = %d, width = %d, height = %d, resolution = %f\n", numofpts, numblunders, width, height, resolution);
    
    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        //index_in_ptslists[((ptslists[t].m_Y - minY_ptslists)/resolution)*width+((ptslists[t].m_X - minX_ptslists) / resolution)] = t;
        INDEX col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        INDEX row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[row * width + col] = t;
    }
    
    GridPoint *grid_blunders = new GridPoint[numblunders];
    for (std::size_t t = 0; t < numblunders; ++t)
    {
        grid_blunders[t].col = 0.5 + (blunderlist[t].m_X - minX_ptslists) / resolution;
        grid_blunders[t].row = 0.5 + (blunderlist[t].m_Y - minY_ptslists) / resolution;
    }
    
    GridPoint **blunder_ptrs = new GridPoint*[numblunders];
#pragma omp parallel for
    for (std::size_t t = 0; t < numblunders; ++t) blunder_ptrs[t] = grid_blunders + t;
    
    //double begin = omp_get_wtime();
    
    oldTri->Retriangulate(blunder_ptrs, numblunders);
    
    //double end = omp_get_wtime();
    //printf("Retriangulate took %lf with %d points, %d blunders\n", end - begin, numofpts, numblunders);
    
    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    *count_tri = (int)(oldTri->GetAllTris(&tris));
    for(long t = 0 ; t < tris.size() ; t++)
    {
        int row, col;
        UI3DPOINT temp_pt;
        
        row = tris[t].pts[0].row;
        col = tris[t].pts[0].col;
        
        temp_pt.m_X = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[1].row;
        col = tris[t].pts[1].col;
        temp_pt.m_Y = index_in_ptslists[row * width + col];
        
        row = tris[t].pts[2].row;
        col = tris[t].pts[2].col;
        temp_pt.m_Z = index_in_ptslists[row * width + col];
        
        trilists.push_back(temp_pt);
    }
    delete [] blunder_ptrs;
    delete [] grid_blunders;
}

void SetTinBoundary(LevelInfo &rlevelinfo, const D3DPOINT &TriP1, const D3DPOINT &TriP2, const D3DPOINT &TriP3, int *PixelMinXY, int *PixelMaxXY, double &Total_Min_Z, double &Total_Max_Z, double &temp_MinZ, double &temp_MaxZ)
{
    temp_MinZ = min(min(TriP1.m_Z,TriP2.m_Z),TriP3.m_Z);
    temp_MaxZ = max(max(TriP1.m_Z,TriP2.m_Z),TriP3.m_Z);
    
    if(temp_MinZ < Total_Min_Z)
        Total_Min_Z = temp_MinZ;
    if(temp_MaxZ > Total_Max_Z)
        Total_Max_Z = temp_MaxZ;
    
    // calculation on BoundingBox(MinMax XY) of triangle
    const double TriMinXY[2] = { min(min(TriP1.m_X,TriP2.m_X),TriP3.m_X), min(min(TriP1.m_Y,TriP2.m_Y),TriP3.m_Y)};
    const double TriMaxXY[2] = { max(max(TriP1.m_X,TriP2.m_X),TriP3.m_X), max(max(TriP1.m_Y,TriP2.m_Y),TriP3.m_Y)};
    PixelMinXY[0] = (int)((TriMinXY[0] - rlevelinfo.Boundary[0])/(*rlevelinfo.grid_resolution) + 0.5);
    PixelMinXY[1] = (int)((TriMinXY[1] - rlevelinfo.Boundary[1])/(*rlevelinfo.grid_resolution) + 0.5);
    PixelMaxXY[0] = (int)((TriMaxXY[0] - rlevelinfo.Boundary[0])/(*rlevelinfo.grid_resolution) + 0.5);
    PixelMaxXY[1] = (int)((TriMaxXY[1] - rlevelinfo.Boundary[1])/(*rlevelinfo.grid_resolution) + 0.5);
    
    PixelMinXY[0] -= 1;     PixelMinXY[1] -= 1;
    PixelMaxXY[0] += 1;     PixelMaxXY[1] += 1;
    
    if (PixelMaxXY[0] >= (int)(rlevelinfo.Size_Grid2D->width))
        PixelMaxXY[0] =  (int)(rlevelinfo.Size_Grid2D->width-1);
    if (PixelMaxXY[1] >= (int)(rlevelinfo.Size_Grid2D->height))
        PixelMaxXY[1] =  (int)(rlevelinfo.Size_Grid2D->height-1);
    if (PixelMinXY[0] < 0)
        PixelMinXY[0] = 0;
    if (PixelMinXY[1] < 0)
        PixelMinXY[1] = 0;
}

bool IsTinInside(const D3DPOINT &CurGPXY, const D3DPOINT &TriP1, const D3DPOINT &TriP2, const D3DPOINT &TriP3, float &Z)
{
    bool rtn = false;
    //float temp_Z;
    D3DPOINT v12(TriP2 - TriP1);
    D3DPOINT v1P(CurGPXY - TriP1);
    D3DPOINT v23(TriP3 - TriP2);
    D3DPOINT v2P(CurGPXY - TriP2);
    D3DPOINT v31(TriP1 - TriP3);
    D3DPOINT v3P(CurGPXY - TriP3);
    
    int Sum = 3;
    if (v12.m_X*v1P.m_Y-v12.m_Y*v1P.m_X <= 0)
        Sum--;
    if (v23.m_X*v2P.m_Y-v23.m_Y*v2P.m_X <= 0)
        Sum--;
    if (v31.m_X*v3P.m_Y-v31.m_Y*v3P.m_X <= 0)
        Sum--;
    
    if (Sum==0 || Sum==3)
    {
        D3DPOINT v12(TriP2 - TriP1);
        D3DPOINT v13(TriP3 - TriP1);
        D3DPOINT Normal(v12.m_Y*v13.m_Z - v12.m_Z*v13.m_Y, v12.m_Z*v13.m_X - v12.m_X*v13.m_Z, v12.m_X*v13.m_Y - v12.m_Y*v13.m_X);
        
        const double Len = SQRT(Normal);
        if(Len > 0 )
        {
            Normal.m_X/=Len;
            Normal.m_Y/=Len;
            Normal.m_Z/=Len;
            
            const double A = Normal.m_X;
            const double B = Normal.m_Y;
            const double C = Normal.m_Z;
            const double D = -(A*TriP1.m_X+B*TriP1.m_Y+C*TriP1.m_Z);
            
            if(C != 0)
            {
                Z = -1.0 * ((A * CurGPXY.m_X) + (B * CurGPXY.m_Y) + D) / C;
                //Z = temp_Z;
                rtn = true;
            }
            else
                rtn = false;
        }
    }
    
    return rtn;
}

double SetNormalAngle(const D3DPOINT &pts0, const D3DPOINT &pts1, const D3DPOINT &pts2)
{
    double angle;
    
    const D3DPOINT U(pts1 - pts0);
    const D3DPOINT V(pts2 - pts0);
    const D3DPOINT N(U.m_Y*V.m_Z - V.m_Y*U.m_Z, -(U.m_X*V.m_Z - V.m_X*U.m_Z), U.m_X*V.m_Y - V.m_X*U.m_Y,0);
    
    double norm  = SQRT(N);
    if (norm != 0)
        angle = acos(fabs(N.m_Z)/norm)*RadToDeg;
    
    SetAngle(angle);
    
    return angle;
}

void SetAngle(double &angle)
{
    if(angle <= 0 && angle >= -90)
        angle = fabs(angle);
    else if(angle <= -270 && angle >= -360)
        angle = 360 + angle;
    else if(angle >= 270 && angle <= 360)
        angle = 360 - angle;
}

//removable functions not used in SETSM DEM main
/*void TINUpdate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution, FullTriangulation *oldTri, D3DPOINT *blunderlist, int numblunders)
{

    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];

    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    printf("\tTINUpdate: PTS = %d, blunders = %d, width = %d, height = %d, resolution = %f\n", numofpts, numblunders, width, height, resolution);

    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        //index_in_ptslists[((ptslists[t].m_Y - minY_ptslists)/resolution)*width+((ptslists[t].m_X - minX_ptslists) / resolution)] = t;
        INDEX col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
    INDEX row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
    index_in_ptslists[row * width + col] = t;
    }

    GridPoint *grid_blunders = new GridPoint[numblunders];
    for (std::size_t t = 0; t < numblunders; ++t)
    {
        grid_blunders[t].col = 0.5 + (blunderlist[t].m_X - minX_ptslists) / resolution;
        grid_blunders[t].row = 0.5 + (blunderlist[t].m_Y - minY_ptslists) / resolution;
    }

    GridPoint **blunder_ptrs = new GridPoint*[numblunders];
    #pragma omp parallel for
    for (std::size_t t = 0; t < numblunders; ++t) blunder_ptrs[t] = grid_blunders + t;
    
    //double begin = omp_get_wtime();

    oldTri->Retriangulate(blunder_ptrs, numblunders);

    //double end = omp_get_wtime();
    //printf("Retriangulate took %lf with %d points, %d blunders\n", end - begin, numofpts, numblunders);

    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    *count_tri = (int)(oldTri->GetAllTris(&tris));
    std::size_t t = 0;
    for(it_tr = tris.begin(); it_tr != tris.end(); ++it_tr)
    {
        int row, col;

        row = it_tr->pts[0].row;
        col = it_tr->pts[0].col;

        trilists[t].m_X = index_in_ptslists[row * width + col];

        row = it_tr->pts[1].row;
        col = it_tr->pts[1].col;
        trilists[t].m_Y = index_in_ptslists[row * width + col];
        
        row = it_tr->pts[2].row;
        col = it_tr->pts[2].col;
        trilists[t].m_Z = index_in_ptslists[row * width + col];
        t++;
    }
    delete [] blunder_ptrs;
    delete [] grid_blunders;
}

FullTriangulation *TINCreate_float(F3DPOINT *ptslists, long numofpts, UI3DPOINT* trilists, double min_max[], long *count_tri, double resolution)
{
    if (numofpts <= 2) {
        *count_tri = 0;
        return NULL;
    }

    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];

    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //Check to ensure grid width/height fit in 16 bit integers
    if(width > 65535 || height > 65535)
    {
        printf("ERROR: Grid is too large. width: %d height: %d\n", width, height);
        exit(1);
    }
    printf("\tTINCreate: PTS = %d, width = %d, height = %d, resolution = %f\n", numofpts, width, height, resolution);

    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    GridPoint *grid_points = new GridPoint[numofpts];
    
    printf("done s1\n");
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        grid_points[t].col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        grid_points[t].row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[grid_points[t].row * width + grid_points[t].col] = t;
        //printf("ID %d\t coord %d\t%d\t%f\n",t,grid_points[t].col,grid_points[t].row,ptslists[t].m_Z);
    }
    printf("done s2\n");
    GridPoint **points_ptrs = new GridPoint*[numofpts];
    #pragma omp parallel for
    for (std::size_t t = 0; t < numofpts; ++t) points_ptrs[t] = grid_points + t;

    FullTriangulation *triangulation = new FullTriangulation(width, height);
    //double begin = omp_get_wtime();
    printf("done s3\n");
    triangulation->Triangulate(points_ptrs, numofpts);
    printf("done triangulation\n");
    //double end = omp_get_wtime();
    //printf("Triangulate took %lf with %d points\n", end - begin, numofpts);

    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    printf("s3-1\n");
    *count_tri = (int)(triangulation->GetAllTris(&tris));
    printf("count tri = %d\n",*count_tri);
    std::size_t t = 0;
    for(it_tr = tris.begin(); it_tr != tris.end(); ++it_tr)
    {
        int row, col;

        row = it_tr->pts[0].row;
        col = it_tr->pts[0].col;
        trilists[t].m_X = index_in_ptslists[row * width + col];

        row = it_tr->pts[1].row;
        col = it_tr->pts[1].col;
        trilists[t].m_Y = index_in_ptslists[row * width + col];
        
        row = it_tr->pts[2].row;
        col = it_tr->pts[2].col;
        trilists[t].m_Z = index_in_ptslists[row * width + col];
        t++;
    }
    printf("s4\n");
    delete [] points_ptrs;
    delete [] grid_points;
    
    printf("s5\n");
    return triangulation;
}

FullTriangulation *TINCreate(D3DPOINT *ptslists, int numofpts, UI3DPOINT* trilists, double min_max[], int *count_tri, double resolution)
{
    if (numofpts <= 2) {
        *count_tri = 0;
        return NULL;
    }

    double minX_ptslists = min_max[0];
    double minY_ptslists = min_max[1];
    double maxX_ptslists = min_max[2];
    double maxY_ptslists = min_max[3];

    INDEX width     = 1 + (maxX_ptslists - minX_ptslists) / resolution;
    INDEX height    = 1 + (maxY_ptslists - minY_ptslists) / resolution;
    //Check to ensure grid width/height fit in 16 bit integers
    if(width > 32767 || height > 32767)
    {
        printf("ERROR: Grid is too large. width: %d height: %d\n", width, height);
        exit(1);
    }
    printf("\tTINCreate: PTS = %d, width = %d, height = %d, resolution = %f\n", numofpts, width, height, resolution);

    std::unordered_map<std::size_t, std::size_t> index_in_ptslists;
    index_in_ptslists.reserve(numofpts);
    GridPoint *grid_points = new GridPoint[numofpts];

    printf("done s1\n");
    for (std::size_t t = 0; t < numofpts; ++t)
    {
        grid_points[t].col = 0.5 + (ptslists[t].m_X - minX_ptslists) / resolution;
        grid_points[t].row = 0.5 + (ptslists[t].m_Y - minY_ptslists) / resolution;
        index_in_ptslists[grid_points[t].row * width + grid_points[t].col] = t;
    }

    printf("done s2\n");
    GridPoint **points_ptrs = new GridPoint*[numofpts];
    #pragma omp parallel for
    for (std::size_t t = 0; t < numofpts; ++t) points_ptrs[t] = grid_points + t;

    FullTriangulation *triangulation = new FullTriangulation(width, height);
    //double begin = omp_get_wtime();
    printf("done s3\n");
    triangulation->Triangulate(points_ptrs, numofpts);
    printf("done triangulation\n");
    //double end = omp_get_wtime();
    //printf("Triangulate took %lf with %d points\n", end - begin, numofpts);

    vector<Tri> tris;
    vector<Tri>::iterator it_tr;
    printf("s3-1\n");
    *count_tri = (int)(triangulation->GetAllTris(&tris));
    printf("count tri = %d\n",*count_tri);
    size_t t=0;
    for(it_tr = tris.begin(); it_tr != tris.end(); ++it_tr)
    {
        int row, col;

        row = it_tr->pts[0].row;
        col = it_tr->pts[0].col;
        trilists[t].m_X = index_in_ptslists[row * width + col];

        row = it_tr->pts[1].row;
        col = it_tr->pts[1].col;
        trilists[t].m_Y = index_in_ptslists[row * width + col];
        
        row = it_tr->pts[2].row;
        col = it_tr->pts[2].col;
        trilists[t].m_Z = index_in_ptslists[row * width + col];
        t++;
    }
    printf("s4\n");
    delete [] points_ptrs;
    delete [] grid_points;
    
    printf("s5\n");
    return triangulation;
}
*/


void SetVecKernelValue(const KernelPatchArg &patch, const int row, const int col, const D2DPOINT &pos_left, const D2DPOINT &pos_right, const int radius2, int *Count_N)
{

    if( pos_right.m_Y >= 0 && pos_right.m_Y + 1 < patch.RImagesize.height && pos_right.m_X  >= 0 && pos_right.m_X + 1 < patch.RImagesize.width && pos_left.m_Y >= 0 && pos_left.m_Y + 1 < patch.LImagesize.height && pos_left.m_X >= 0 && pos_left.m_X + 1  < patch.LImagesize.width)
    {
        const long int position = (long int) pos_left.m_X + (long int) pos_left.m_Y *(long int)patch.LImagesize.width;
        const long int position_right = (long int) pos_right.m_X + (long int) pos_right.m_Y *(long int)patch.RImagesize.width;
        
        double left_patch, left_mag_patch, right_patch, right_mag_patch;
        
        double dx          =   pos_left.m_X - floor(pos_left.m_X);
        double dy          =   pos_left.m_Y - floor(pos_left.m_Y);
        double dx_r        =  pos_right.m_X - floor(pos_right.m_X);
        double dy_r        =  pos_right.m_Y - floor(pos_right.m_Y);

        //interpolate left_patch
        left_patch = InterpolatePatch(patch.left_image, position, patch.LImagesize, dx, dy);
        left_mag_patch = InterpolatePatch(patch.left_mag_image, position,patch.LImagesize, dx, dy);
        
        //interpolate right_patch
        right_patch = InterpolatePatch(patch.right_image, position_right,patch.RImagesize, dx_r, dy_r);
        right_mag_patch = InterpolatePatch(patch.right_mag_image, position_right,patch.RImagesize, dx_r, dy_r);

        auto &rkernel = patch.rkernel;
            
        if(left_patch > 0 && right_patch > 0)
        {
            rkernel.left_patch_vecs(0, Count_N[0]) = left_patch;
            rkernel.left_mag_patch_vecs(0, Count_N[0]) = left_mag_patch;
            rkernel.right_patch_vecs(0, Count_N[0]) = right_patch;
            rkernel.right_mag_patch_vecs(0, Count_N[0]) = right_mag_patch;
            Count_N[0]++;
            
            const int size_1        = (int)(rkernel.Half_template_size/2);
            if(radius2 <= (rkernel.Half_template_size - size_1 + 1)*(rkernel.Half_template_size - size_1 + 1))
            {
                if(row >= -rkernel.Half_template_size + size_1 && row <= rkernel.Half_template_size - size_1)
                {
                    if(col >= -rkernel.Half_template_size + size_1 && col <= rkernel.Half_template_size - size_1)
                    {
                        rkernel.left_patch_vecs(1, Count_N[1]) = left_patch;
                        rkernel.left_mag_patch_vecs(1, Count_N[1]) = left_mag_patch;
                        rkernel.right_patch_vecs(1, Count_N[1]) = right_patch;
                        rkernel.right_mag_patch_vecs(1, Count_N[1]) = right_mag_patch;
                        Count_N[1]++;
                    }
                }
            }
            
            const int size_2        = size_1 + (int)((size_1/2.0) + 0.5);
            if(radius2 <= (rkernel.Half_template_size - size_2 + 1)*(rkernel.Half_template_size - size_2 + 1))
            {
                if(row >= -rkernel.Half_template_size + size_2 && row <= rkernel.Half_template_size - size_2)
                {
                    if(col >= -rkernel.Half_template_size + size_2 && col <= rkernel.Half_template_size - size_2)
                    {
                        rkernel.left_patch_vecs(2, Count_N[2]) = left_patch;
                        rkernel.left_mag_patch_vecs(2, Count_N[2]) = left_mag_patch;
                        rkernel.right_patch_vecs(2, Count_N[2]) = right_patch;
                        rkernel.right_mag_patch_vecs(2, Count_N[2]) = right_mag_patch;
                        Count_N[2]++;
                    }
                }
            }
        }
    }
}

double InterpolatePatch(const uint16 *Image, const long int position, const CSize Imagesize, const double dx, const double dy)
{
    double patch =
        (double) (Image[position]) * (1 - dx) * (1 - dy) +
        (double) (Image[position + 1]) * dx * (1 - dy) +
        (double) (Image[position + Imagesize.width]) * (1 - dx) * dy +
        (double) (Image[position + 1 + Imagesize.width]) * dx * dy;

    return patch;
}

void ComputeMultiNCC(SetKernel &rsetkernel, const int Th_rho, const int *Count_N, double &count_NCC, double &sum_NCC_multi)
{
    if(Count_N[0] > Th_rho && Count_N[1] > Th_rho && Count_N[2] > Th_rho)
    {
        double temp_roh = 0;
        double count_roh = 0;
        double temp_NCC_roh;
        
        for (int k=0; k<3; k++)
        {
            const double ncc = Correlate(rsetkernel.left_patch_vecs.row(k), rsetkernel.right_patch_vecs.row(k), Count_N[k]);
            if (ncc != -99)
            {
                count_roh++;
                temp_roh += ncc;
            }
            
            const double ncc_mag = Correlate(rsetkernel.left_mag_patch_vecs.row(k), rsetkernel.right_mag_patch_vecs.row(k), Count_N[k]);
            if (ncc_mag != -99)
            {
                count_roh++;
                temp_roh += ncc_mag;
            }
        }
        if (count_roh > 0)
        {
            temp_NCC_roh = temp_roh/count_roh;
            sum_NCC_multi += temp_NCC_roh;
            count_NCC ++;
        }
    }
}



D2DPOINT *SetDEMGrid(const double *Boundary, const double Grid_x, const double Grid_y, CSize *Size_2D)
{
    D2DPOINT *GridPT = NULL;

    Size_2D->width  = (int)(ceil((double)(Boundary[2] - Boundary[0])/Grid_x));
    Size_2D->height = (int)(ceil((double)(Boundary[3] - Boundary[1])/Grid_y));

    GridPT  = (D2DPOINT*)malloc(sizeof(D2DPOINT)*(long)Size_2D->height*(long)Size_2D->width);
#pragma omp parallel for schedule(guided)
    for(long row = 0 ; row < Size_2D->height ; row++)
        for(long col = 0; col < Size_2D->width ; col++)
        {
            long index = row*(long)Size_2D->width + col;
            GridPT[index].m_X = Boundary[0] + col*Grid_x;
            GridPT[index].m_Y = Boundary[1] + row*Grid_y;
        }

    return GridPT;
}

void SetPyramidImages(const ProInfo *proinfo, const int py_level_set, const CSize * const *data_size_lr, uint16 ***SubImages, uint16 ***SubMagImages, uint8 ***SubOriImages)
{
    for(int iter_level = 0 ; iter_level < py_level_set; iter_level++)
    {
        for(int image_index = 0 ; image_index < proinfo->number_of_images ; image_index++)
        {
            if(proinfo->check_selected_image[image_index])
            {
                long int data_length = (long int)data_size_lr[image_index][iter_level].height*(long int)data_size_lr[image_index][iter_level].width;
   
                if(iter_level > 0)
                {
                    SubImages[iter_level][image_index] = CreateImagePyramid(SubImages[iter_level-1][image_index],data_size_lr[image_index][iter_level-1],9,(double)(1.5));
                }
                
                int16 *dirimg = (int16*)malloc(sizeof(int16)*data_length);
                MakeSobelMagnitudeImage(data_size_lr[image_index][iter_level],SubImages[iter_level][image_index],SubMagImages[iter_level][image_index],dirimg);
                Orientation(data_size_lr[image_index][iter_level],SubMagImages[iter_level][image_index],dirimg,15,SubOriImages[iter_level][image_index]);
                
                free(dirimg);
            }
        }
    }
}

uint16 *SubsetImageFrombitsToUint16(const int image_bits, char *imagefile, long *cols, long *rows, CSize *subsize)
{
    uint16 *out = NULL;
    
    CSize temp_size;
    GetImageSize(imagefile,&temp_size);
    switch(image_bits)
    {
        case 8:
        {
            uint8 type_t(0);
            uint8 *t_data8    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            if(t_data8)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint16*)malloc(sizeof(uint16)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                    out[index] = t_data8[index];
                free(t_data8);
            }
            break;
        }
        case 12:
        {
            uint16 type_t(0);
            uint16 *t_data16    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            if(t_data16)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint16*)malloc(sizeof(uint16)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                    out[index] = t_data16[index];
                free(t_data16);
            }
            break;
        }
    }
    
    return out;
}

uint8 *SubsetImageFrombitsToUint8(const int image_bits, char *imagefile, long *cols, long *rows, CSize *subsize)
{
    uint8 *out = NULL;
    
    CSize temp_size;
    GetImageSize(imagefile,&temp_size);
    switch(image_bits)
    {
        case 8:
        {
            uint8 type_t(0);
            uint8 *t_data8    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            if(t_data8)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint8*)malloc(sizeof(uint8)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                    out[index] = t_data8[index];
                free(t_data8);
            }
            break;
        }
        case 12:
        {
            uint16 type_t(0);
            uint16 *t_data16    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            if(t_data16)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint8*)malloc(sizeof(uint8)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                    out[index] = t_data16[index];
                free(t_data16);
            }
            break;
        }
        case 32:
        {
            float type_t(0);
            float *t_datafloat    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            int count = 0;
            if(t_datafloat)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint8*)malloc(sizeof(uint8)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                {
                    out[index] = t_datafloat[index];
                    if(out[index] == 1)
                        count++;
                }
                printf("count %d\ttotal count %ld\n",count,data_size);
                free(t_datafloat);
            }
            break;
        }
        case 64:
        {
            double type_t(0);
            double *t_datafloat    = Readtiff_T(imagefile,&temp_size,cols,rows,subsize, type_t);
            if(t_datafloat)
            {
                long data_size = (long)subsize->width*(long)subsize->height;
                out = (uint8*)malloc(sizeof(uint8)*data_size);
                #pragma omp parallel for schedule(guided)
                for(long index = 0 ; index < data_size ; index++)
                    out[index] = t_datafloat[index];
                free(t_datafloat);
            }
            break;
        }
    }
    
    return out;
}

uint16* LoadPyramidImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level)
{
    uint16 *out = NULL;

    char t_str[500];
    char *filename_py     = GetFileName(subsetfile);
    filename_py     = remove_ext(filename_py);
    sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,py_level);
    free(filename_py);

    FILE *pFile           = fopen(t_str,"rb");
    if(pFile)
    {
        out         = (uint16*)malloc(sizeof(uint16)*data_size.height*data_size.width);
        fread(out,sizeof(uint16),data_size.height*data_size.width,pFile);
        fclose(pFile);
    }

    return out;
}

uint8* LoadPyramidOriImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level)
{
    uint8 *out = NULL;
    char t_str[500];
    char *filename_py     = GetFileName(subsetfile);
    filename_py     = remove_ext(filename_py);
    sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,py_level);
    free(filename_py);
    
    FILE *pFile           = fopen(t_str,"rb");
    if(pFile)
    {
        out     = (uint8*)malloc(sizeof(uint8)*data_size.height*data_size.width);
        fread(out,sizeof(uint8),data_size.height*data_size.width,pFile);
        fclose(pFile);
    }

    return out;
}


uint16* LoadPyramidMagImages(const char *save_path,char *subsetfile,const CSize data_size,const uint8 py_level)
{
    uint16 *out = NULL;
    char t_str[500];
    char *filename_py     = GetFileName(subsetfile);
    filename_py     = remove_ext(filename_py);
    sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,py_level);
    free(filename_py);
    
    FILE *pFile           = fopen(t_str,"rb");
    if(pFile)
    {
        out     = (uint16*)malloc(sizeof(uint16)*data_size.height*data_size.width);
        fread(out,sizeof(uint16),data_size.height*data_size.width,pFile);
        fclose(pFile);
    }

    return out;
}
void RemoveFiles(const ProInfo *proinfo, const char *save_path, char **filename, const int py_level, const bool flag)
{
    int start_lv, end_lv;

    start_lv    = py_level;
    if(flag)
        end_lv      = py_level - 1;
    else
        end_lv      = -2;

    if(py_level == 0)
    {
        start_lv    = 0;
        end_lv      = -2;
    }

    for(int ti = 0 ; ti < proinfo->number_of_images ; ti++)
    {
        if(proinfo->check_selected_image[ti])
        {
            for(int count = py_level ; count >= 0 ; count--)
            {
                char *filename_py     = GetFileName(filename[ti]);
                char t_str[500];
                filename_py     = remove_ext(filename_py);
                sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,count);
                int status = remove(t_str);
                free(filename_py);

                filename_py     = GetFileName(filename[ti]);
                filename_py     = remove_ext(filename_py);
                sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,count);
                status = remove(t_str);
                free(filename_py);

                filename_py     = GetFileName(filename[ti]);
                filename_py     = remove_ext(filename_py);
                sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,count);
                status = remove(t_str);
                free(filename_py);

                filename_py     = GetFileName(filename[ti]);
                filename_py     = remove_ext(filename_py);
                sprintf(t_str,"%s/%s.raw",save_path,filename_py);
                status = remove(t_str);
                free(filename_py);
            }
        }
    }
}

void Preprocessing(const ProInfo *proinfo, const char *save_path,char **Subsetfile, const uint8 py_level, const CSize * const *data_size_lr, FILE *fid)
{
    for(uint8 count = 0; count<proinfo->number_of_images ; count++)
    {

        char t_str[500];
        FILE *pFile_raw   = fopen(Subsetfile[count],"rb");
        char *filename_py     = GetFileName(Subsetfile[count]);
        filename_py     = remove_ext(filename_py);

        sprintf(t_str,"%s/%s_py_5.raw",save_path,filename_py);
        FILE *pFile_check_file    = fopen(t_str,"rb");

        if(pFile_raw && !pFile_check_file && proinfo->check_selected_image[count])
        {
            uint16 **pyimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
            uint16 **magimg = (uint16**)malloc(sizeof(uint16*)*(py_level+1));
            int16 **dirimg = (int16**)malloc(sizeof(int16*)*(py_level+1));
            uint8 **oriimg = (uint8**)malloc(sizeof(uint8*)*(py_level+1));
            CSize *data_size = (CSize*)malloc(sizeof(CSize)*(py_level+1));

            fprintf(fid,"image %d\tlevel = 0\t width = %d\theight = %d\n",count,data_size_lr[count][0].width,data_size_lr[count][0].height);
            for(int i=0;i<py_level+1;i++)
                data_size[i]        = data_size_lr[count][i];
            
            long int data_length = (long int)data_size[0].height*(long int)data_size[0].width;
            pyimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
            magimg[0] = (uint16*)malloc(sizeof(uint16)*data_length);
            dirimg[0] = (int16*)malloc(sizeof(int16)*data_length);
            oriimg[0] = (uint8*)malloc(sizeof(uint8)*data_length);

            fread(pyimg[0],sizeof(uint16),data_length,pFile_raw);

            MakeSobelMagnitudeImage(data_size[0],pyimg[0],magimg[0],dirimg[0]);
            Orientation(data_size[0],magimg[0],dirimg[0],15,oriimg[0]);
            
            sprintf(t_str,"%s/%s_py_0.raw",save_path,filename_py);
            FILE *pFile   = fopen(t_str,"wb");
            fwrite(pyimg[0],sizeof(uint16),data_length,pFile);
            fclose(pFile);
            
            sprintf(t_str,"%s/%s_py_0_mag.raw",save_path,filename_py);
            pFile   = fopen(t_str,"wb");
            fwrite(magimg[0],sizeof(uint16),data_length,pFile);
            fclose(pFile);
            free(magimg[0]);

            free(dirimg[0]);

            sprintf(t_str,"%s/%s_py_0_ori.raw",save_path,filename_py);
            pFile   = fopen(t_str,"wb");
            fwrite(oriimg[0],sizeof(uint8),data_length,pFile);
            fclose(pFile);
            free(oriimg[0]);

            for(int i=0;i<py_level;i++)
            {
                pyimg[i+1] = CreateImagePyramid(pyimg[i],data_size[i],9,(double)(1.5));
                free(pyimg[i]);

                fprintf(fid,"image %d\tlevel = %d\t width = %d\theight = %d\n",count,i+1,data_size_lr[count][i+1].width,data_size_lr[count][i+1].height);
        
                long int data_length_array[6];
                data_length_array[i] = (long int)data_size[i+1].height*(long int)data_size[i+1].width;
                magimg[i+1] = (uint16*)malloc(sizeof(uint16)*data_length_array[i]);
                dirimg[i+1] = (int16*)malloc(sizeof(int16)*data_length_array[i]);
                oriimg[i+1] = (uint8*)malloc(sizeof(uint8)*data_length_array[i]);
                MakeSobelMagnitudeImage(data_size[i+1],pyimg[i+1],magimg[i+1],dirimg[i+1]);
                Orientation(data_size[i+1],magimg[i+1],dirimg[i+1],15,oriimg[i+1]);

                sprintf(t_str,"%s/%s_py_%d.raw",save_path,filename_py,i+1);
                pFile   = fopen(t_str,"wb");
                fwrite(pyimg[i+1],sizeof(uint16),data_length_array[i],pFile);
                fclose(pFile);
                if(i == py_level-1)
                    free(pyimg[i+1]);

                sprintf(t_str,"%s/%s_py_%d_mag.raw",save_path,filename_py,i+1);
                pFile   = fopen(t_str,"wb");
                fwrite(magimg[i+1],sizeof(uint16),data_length_array[i],pFile);
                fclose(pFile);
                free(magimg[i+1]);

                free(dirimg[i+1]);

                sprintf(t_str,"%s/%s_py_%d_ori.raw",save_path,filename_py,i+1);
                pFile   = fopen(t_str,"wb");
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

        free(filename_py);
    }
}

//removable by Geotiff reader
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

CSize Envihdr_reader_seedDEM(TransParam _param, char *filename, double *minX, double *maxY, double *grid_size)
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

bool TFW_reader_LSFDEM(char *filename, double *minX, double *maxY, double *grid_size, int *zone, char *dir)
{
    FILE* fid;
    double garbage;
    printf("%s\n",filename);
    fid = fopen(filename,"r");
    int count_line = 0;
    char bufstr[500];
    while(!feof(fid))
    {
        double temp;
        fgets(bufstr,300,fid);
        if(count_line == 0)
        {
            sscanf(bufstr,"%lf",&temp);
            *grid_size = temp;
        }
        else if(count_line == 4)
        {
            sscanf(bufstr,"%lf",&temp);
            *minX = temp;
        }
        else if(count_line == 5)
        {
            sscanf(bufstr,"%lf",&temp);
            *maxY = temp;
        }
        else if(count_line == 7)
        {
            int ttt;
            sscanf(bufstr,"%d",&ttt);
            *zone = ttt;
        }
        else if(count_line == 8)
        {
            char ttt[100];
            sscanf(bufstr,"%s",ttt);
            sprintf(dir,"%s",ttt);
            
        }
        count_line++;
    }
    
    *minX = (*minX) - (*grid_size)/2.0;
    *maxY = (*maxY) + (*grid_size)/2.0;
    
    
    printf("%f\n",*minX);
    printf("%f\n",*maxY);
    printf("%f\n",*grid_size);
    if(count_line > 5)
    {
        printf("%d\n",*zone);
        printf("%s\n",dir);
    }
    
    fclose(fid);
    
    return true;
}


