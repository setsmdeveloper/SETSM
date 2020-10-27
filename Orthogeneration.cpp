//
//  Orthogeneration.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "Orthogeneration.hpp"


//orthogeneration
void orthogeneration(const TransParam _param, const ARGINFO args, char *ImageFilename, char *DEMFilename, const char *Outputpath, const int pair, const int DEM_divide, const double * const *Imageparams)
{
    if(args.RA_only) {
        return;
    }
    
    ProInfo *proinfo = new ProInfo;
    proinfo->sensor_type = args.sensor_type;
    proinfo->number_of_images = 1;
    
    char DEM_header[500];
    char RPCFilename[500];
    char OrthoFilename[500];
    char OrthoGEOTIFFFilename[500];
    char Ortho_header[500];
    
    time_t ST = 0, ET = 0;
    double gap;
    
    ST = time(0);
    FrameInfo m_frameinfo;
    m_frameinfo.m_Camera.m_focalLength  = 0;
    m_frameinfo.m_Camera.m_CCDSize      = 0;
    
    if(args.sensor_type == AB)
    {
        m_frameinfo.m_Camera.m_focalLength  = args.focal_length;
        m_frameinfo.m_Camera.m_CCDSize      = args.CCD_size;
        m_frameinfo.m_Camera.m_ppx = 0;
        m_frameinfo.m_Camera.m_ppy = 0;
        
        m_frameinfo.Photoinfo = (EO*)malloc(sizeof(EO));
        sprintf(m_frameinfo.Photoinfo[0].path,"%s",ImageFilename);
    }
    
    printf("sensor %f\t%f\n",m_frameinfo.m_Camera.m_focalLength,m_frameinfo.m_Camera.m_CCDSize);
    
    char *tmp_chr = remove_ext(ImageFilename);
    sprintf(RPCFilename,"%s.xml",tmp_chr);
    sprintf(proinfo->RPCfilename[0],"%s",RPCFilename);
    FILE *fid_xml = fopen(RPCFilename,"r");
    if(!fid_xml)
    {
        sprintf(RPCFilename,"%s.XML",tmp_chr);
        sprintf(proinfo->RPCfilename[0],"%s",RPCFilename);
        FILE *fid_XML = fopen(RPCFilename,"r");
        if(!fid_XML)
        {
            printf("Please check xml file!! SETSM supports a format of 'xml' or 'XML'");
            exit(1);
        }
        else
            fclose(fid_XML);
    }
    else
        fclose(fid_xml);
    
    free(tmp_chr);

    tmp_chr = remove_ext(DEMFilename);
    sprintf(DEM_header,"%s.hdr",tmp_chr);
    char *Ifilename  = SetOutpathName(ImageFilename);
    char *tmp_no_ext = remove_ext(Ifilename);
    
    if(DEM_divide == 0)
    {
        sprintf(OrthoFilename, "%s/%s_ortho_%3.1f.raw",Outputpath,tmp_no_ext,args.DEM_space);
        sprintf(Ortho_header, "%s/%s_ortho_%3.1f.hdr",Outputpath, tmp_no_ext,args.DEM_space);
        sprintf(OrthoGEOTIFFFilename, "%s/%s_ortho_%3.1f.tif",Outputpath, tmp_no_ext,args.DEM_space);
    }
    else
    {
        sprintf(OrthoFilename, "%s/%s_%d_ortho_%3.1f.raw",Outputpath, tmp_no_ext,DEM_divide,args.DEM_space);
        sprintf(Ortho_header, "%s/%s_%d_ortho_%3.1f.hdr",Outputpath, tmp_no_ext,DEM_divide,args.DEM_space);
        sprintf(OrthoGEOTIFFFilename, "%s/%s_%d_ortho_%3.1f.tif",Outputpath, tmp_no_ext,DEM_divide,args.DEM_space);
    }

    free(tmp_no_ext);
    free(tmp_chr);
    free(Ifilename);

    printf("image = %s\n",ImageFilename);
    printf("rpc = %s\n",RPCFilename);
    printf("save = %s\n",Outputpath);
    printf("DEM = %s\n",DEMFilename);
    printf("DEM hdr= %s\n",DEM_header);
    printf("ortho = %s\n",OrthoFilename);
    printf("ortho hdr= %s\n",Ortho_header);
    printf("ortho geotiff = %s\n", OrthoGEOTIFFFilename);
    
    double Image_resolution = 0.5;
    
    // load RPCs info from xml file
    double row_grid_size, col_grid_size, product_grid_size;
    double** RPCs;
    BandInfo band;
    if(args.sensor_type == SB)
    {
        if(args.sensor_provider == DG)
            RPCs            = OpenXMLFile(proinfo, 0, &row_grid_size, &col_grid_size,&product_grid_size, &band);
        else if(args.sensor_provider == PL)
            RPCs            = OpenXMLFile_Pleiades(RPCFilename);
        else if(args.sensor_provider == PT)
            RPCs            = OpenXMLFile_Planet(RPCFilename);
    }
    else
    {
        FILE *pFile           = fopen(RPCFilename,"r");
        
        fscanf(pFile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &m_frameinfo.Photoinfo[0].m_Xl,&m_frameinfo.Photoinfo[0].m_Yl,&m_frameinfo.Photoinfo[0].m_Zl,
               &m_frameinfo.Photoinfo[0].m_Wl,&m_frameinfo.Photoinfo[0].m_Pl,&m_frameinfo.Photoinfo[0].m_Kl);
        
        printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               m_frameinfo.Photoinfo[0].path,
               m_frameinfo.Photoinfo[0].m_Xl,m_frameinfo.Photoinfo[0].m_Yl,m_frameinfo.Photoinfo[0].m_Zl,
               m_frameinfo.Photoinfo[0].m_Wl,m_frameinfo.Photoinfo[0].m_Pl,m_frameinfo.Photoinfo[0].m_Kl);
        double o = m_frameinfo.Photoinfo[0].m_Wl;
        double p = m_frameinfo.Photoinfo[0].m_Pl;
        double k = m_frameinfo.Photoinfo[0].m_Kl;
        
        m_frameinfo.Photoinfo[0].m_Rm = MakeRotationMatrix(o, p, k);
        
        printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",m_frameinfo.Photoinfo[0].m_Rm.m11,m_frameinfo.Photoinfo[0].m_Rm.m12,
               m_frameinfo.Photoinfo[0].m_Rm.m13,m_frameinfo.Photoinfo[0].m_Rm.m21,m_frameinfo.Photoinfo[0].m_Rm.m22,m_frameinfo.Photoinfo[0].m_Rm.m23,
               m_frameinfo.Photoinfo[0].m_Rm.m31,m_frameinfo.Photoinfo[0].m_Rm.m32,m_frameinfo.Photoinfo[0].m_Rm.m33);
        Image_resolution = m_frameinfo.m_Camera.m_CCDSize*UMToMM/(m_frameinfo.m_Camera.m_focalLength/m_frameinfo.Photoinfo[0].m_Zl);
        
        Image_resolution = (int)((Image_resolution)*10 + 0.5)/10.0;
    }
    
    if(args.sensor_type == SB)
    {
        if(!args.check_imageresolution)
        {
            if(args.sensor_provider == DG)
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
    }

    printf("Image resolution %f\n",Image_resolution);
    printf("Hemis projection %d %d\n",_param.bHemisphere, _param.projection);
    printf("param %s %d %d\n", _param.direction,_param.utm_zone,_param.projection);
    
    // load DEM infor from geotiff file.
    double Ortho_resolution, DEM_resolution, DEM_minX, DEM_maxY;
    CSize DEM_size = ReadGeotiff_info(DEMFilename, &DEM_minX, &DEM_maxY, &DEM_resolution);
    
    if (!args.check_DEM_space)
        Ortho_resolution = DEM_resolution;
    else
        Ortho_resolution = args.DEM_space;
    
    double resolution[3] = {DEM_resolution, Image_resolution, Ortho_resolution};
    
    if(fabs(resolution[1] - resolution[2]) <= 1)
        resolution[1] = Ortho_resolution;
    
    const double OrthoGridFactor = resolution[2]/resolution[0];
    
    int impyramid_step  = ceil(log(resolution[2]/resolution[1])/log(2));
    printf("impyramid_step %d\n ortho_Resolution %f\n",impyramid_step,Ortho_resolution);
    
    bool check_overlap = false;
    CSize Orthoimagesize_temp;
    double OrthoBoundary[4];
    // set generated orthoimage info by comparing DEM info
    if(args.sensor_type == SB)
        check_overlap = SetOrthoBoundary_ortho(&Orthoimagesize_temp, OrthoBoundary, RPCs, DEM_resolution, DEM_size, DEM_minX, DEM_maxY, _param, Ortho_resolution);
    else
    {
        GetImageSize(ImageFilename,&m_frameinfo.m_Camera.m_ImageSize);
        check_overlap = SetDEMBoundary_ortho_photo(&Orthoimagesize_temp, OrthoBoundary,DEM_resolution, DEM_size, DEM_minX, DEM_maxY, Ortho_resolution, m_frameinfo.Photoinfo[0], m_frameinfo.m_Camera, m_frameinfo.Photoinfo[0].m_Rm);
    }
    
    if(check_overlap)
    {
        // set saving pointer for orthoimage
        CSize Orthoimagesize(DEM_size.width, DEM_size.height);
        OrthoBoundary[0]  = DEM_minX;
        OrthoBoundary[1]  = DEM_maxY-DEM_size.height*DEM_resolution;
        OrthoBoundary[2]  = DEM_minX+DEM_size.width*DEM_resolution;
        OrthoBoundary[3]  = DEM_maxY;
        
        uint16 *result_ortho    = (uint16*)calloc((long)Orthoimagesize.width*(long)Orthoimagesize.height,sizeof(uint16));
        
        // load DEM value;
        float *DEM_value = GetDEMValue(DEMFilename, DEM_size);
        
        printf("%d\n",DEM_size.width);
        printf("%d\n",DEM_size.height);
        printf("%f\n",DEM_minX);
        printf("%f\n",DEM_maxY);
        printf("%f\n",DEM_resolution);
        printf("%f\n",Ortho_resolution);
        
        double minmaxHeight[2] = {99999.0, -99999.0};
        for(long i=0;i<DEM_size.height;i++)
        {
            for(long j=0;j<DEM_size.width;j++)
            {
                long index = i*(long)DEM_size.width + j;
                if( (minmaxHeight[0] > DEM_value[index]) && (DEM_value[index] > -1000))
                    minmaxHeight[0] = DEM_value[index];
                
                if( (minmaxHeight[1] < DEM_value[index]) && (DEM_value[index] > -1000))
                    minmaxHeight[1] = DEM_value[index];
            }
        }
        
        double subfactor                      = pow(4-impyramid_step,2.0);
        if(subfactor <= 1)
            subfactor                       = 1;
        const int sub_height                  = ceil(Orthoimagesize.height/subfactor);
        const int sub_width                   = ceil(Orthoimagesize.width/subfactor);
        const int height_interval             = Ortho_resolution*sub_height;
        const int width_interval              = Ortho_resolution*sub_width;
        
        const int buffer_x                    = Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
        const int buffer_y                    = Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
        
        double imageparam[2];
        if(pair == 1)
        {
            imageparam[0] = 0.0;
            imageparam[1] = 0.0;
        }
        else
        {
            imageparam[0] = Imageparams[1][0];
            imageparam[1] = Imageparams[1][1];
        }
        printf("Image ID %d\tRPCs bias %f\t%f\n",pair,imageparam[0],imageparam[1]);
        
        printf("Orthoimage info size %d\t%d\tBR %f\t%f\t%f\t%f\n",Orthoimagesize.width,Orthoimagesize.height,OrthoBoundary[0],OrthoBoundary[1],OrthoBoundary[2],OrthoBoundary[3]);
        
        long data_length_ortho = (long)Orthoimagesize.width*(long)Orthoimagesize.height;
        for(long i=0; i < subfactor;i++)
        {
            for(long j=0; j<subfactor;j++)
            {
                char MEG[500];
                sprintf(MEG,"Tile %d %d, processing:%6.1f\n",j+1,i+1,((i)*subfactor+j+1)/subfactor/subfactor*100);
                printf("%s",MEG);

                double Y_size[2];
                Y_size[0]        = OrthoBoundary[3] - (i+1)*height_interval;
                Y_size[1]        = OrthoBoundary[3] - i*height_interval;
                double X_size[2];
                X_size[0]        = OrthoBoundary[0] + j*width_interval;
                X_size[1]        = OrthoBoundary[0] + (j+1)*width_interval;
                
                X_size[0]       = X_size[0] - buffer_x;
                Y_size[0]       = Y_size[0] - buffer_y;
                Y_size[1]       = Y_size[1] + buffer_y;
                X_size[1]       = X_size[1] + buffer_x;
                
                if( X_size[0] < OrthoBoundary[0])
                    X_size[0]   = OrthoBoundary[0];
                
                if (X_size[1] > OrthoBoundary[2])
                    X_size[1]   = OrthoBoundary[2];
                
                if (Y_size[0] < OrthoBoundary[1])
                    Y_size[0]   = OrthoBoundary[1];
                
                if (Y_size[1] > OrthoBoundary[3])
                    Y_size[1]   = OrthoBoundary[3];
                
                double subBoundary[4];
                subBoundary[0] = X_size[0];
                subBoundary[1] = Y_size[0];
                subBoundary[2] = X_size[1];
                subBoundary[3] = Y_size[1];
                
                D2DPOINT startpos_ori;
                CSize subsetsize;
                bool check_subsetImage = false;
                
                uint16 *subimage = subsetImage_ortho(args.sensor_type, m_frameinfo, _param, imageparam, RPCs, ImageFilename,
                                             subBoundary,  minmaxHeight, &startpos_ori, &subsetsize, &check_subsetImage);
                if(check_subsetImage)
                {
                    const int ori_impyramid_step = impyramid_step;
                    if(impyramid_step > 0)
                        impyramid_step = 1;
                    
                    CSize* data_size = new CSize[impyramid_step+1];
                    D2DPOINT startpos;
                    
                    SetPySizes(data_size, subsetsize, impyramid_step);
                    startpos.m_X       = (double)(startpos_ori.m_X/pwrtwo(impyramid_step));      startpos.m_Y       = (double)(startpos_ori.m_Y/pwrtwo(impyramid_step));
                    
                    uint16 *pyimg = NULL;
                    if(impyramid_step > 0)
                        pyimg = Preprocessing_ortho(ori_impyramid_step,data_size,subimage);
                    
                    CSize Image_size  = data_size[impyramid_step];
                    long data_length_image = (long)Image_size.width*(long)Image_size.height;
                    const int col_size    = (int)((X_size[1] - X_size[0])/Ortho_resolution + 0.5);
                    const int row_size    = (int)((Y_size[1] - Y_size[0])/Ortho_resolution + 0.5);
                    
    #pragma omp parallel for schedule(guided)
                    for(long count = 0; count < (long)col_size*(long)row_size ; count++)
                    {
                        const double row  = (floor(count/col_size))*Ortho_resolution + Y_size[0];
                        const double col  = (count % col_size)*Ortho_resolution + X_size[0];

                        double t_col       = (col - DEM_minX)/DEM_resolution;
                        double t_row       = (DEM_maxY - row)/DEM_resolution;
                        
                        long t_col_int   = (long)(t_col + 0.01);
                        long t_row_int   = (long)(t_row + 0.01);
                        
                        if(t_col_int >= 0 && t_col_int +1 < DEM_size.width && t_row_int >= 0 && t_row_int +1 < DEM_size.height)
                        {
                            long index  = (t_col_int   ) + (t_row_int   )*(long)DEM_size.width;
                            double value = DEM_value[index];
                            D3DPOINT object;
                            D2DPOINT objectXY;
                            object.m_X  = col;
                            object.m_Y  = row;
                            object.m_Z  = value;
                            
                            objectXY.m_X  = col;
                            objectXY.m_Y  = row;
                            
                            if(value > -1000)
                            {
                                D2DPOINT image;
                                D2DPOINT temp_pt;
                                if(args.sensor_type == SB)
                                {
                                    D2DPOINT wgsPt = ps2wgs_single(_param, objectXY);
                                   
                                    object.m_X  = wgsPt.m_X;
                                    object.m_Y  = wgsPt.m_Y;
                                    image = GetObjectToImageRPC_single(RPCs, 2, imageparam, object);
                                }
                                else
                                {
                                    D2DPOINT photo  = GetPhotoCoordinate_single(object,m_frameinfo.Photoinfo[0],m_frameinfo.m_Camera,m_frameinfo.Photoinfo[0].m_Rm);
                                    image = PhotoToImage_single(photo, m_frameinfo.m_Camera.m_CCDSize, m_frameinfo.m_Camera.m_ImageSize);
                                }
                                
                                temp_pt     = OriginalToPyramid_single(image, startpos, impyramid_step);
                                
                                t_col       = temp_pt.m_X;
                                t_row       = temp_pt.m_Y;
                                
                                t_col_int   = (long int)(t_col + 0.01);
                                t_row_int   = (long int)(t_row + 0.01);
                                
                                double dcol        = t_col - t_col_int;
                                double drow        = t_row - t_row_int;
                                
                                if(t_col_int >= 0 && t_col_int +1 < Image_size.width && t_row_int >= 0 && t_row_int +1 < Image_size.height
                                   && (t_col_int +1) + (t_row_int +1)*(long)Image_size.width < data_length_image)
                                {
                                    double value1, value2, value3, value4;
                                    long index1  = (t_col_int   ) + (t_row_int   )*(long)Image_size.width;
                                    long index2  = (t_col_int +1) + (t_row_int   )*(long)Image_size.width;
                                    long index3  = (t_col_int   ) + (t_row_int +1)*(long)Image_size.width;
                                    long index4  = (t_col_int +1) + (t_row_int +1)*(long)Image_size.width;
                                    
                                    if(impyramid_step > 0)
                                    {
                                        value1      = pyimg[index1];
                                        value2      = pyimg[index2];
                                        value3      = pyimg[index3];
                                        value4      = pyimg[index4];
                                    }
                                    else
                                    {
                                        value1      = subimage[index1];
                                        value2      = subimage[index2];
                                        value3      = subimage[index3];
                                        value4      = subimage[index4];
                                    }

                                    value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow) + value3*(1-dcol)*drow + value4*dcol*drow;
                                    
                                    t_col_int       = (long )((col - OrthoBoundary[0])/Ortho_resolution + 0.01);
                                    t_row_int       = (long )((OrthoBoundary[3] - row)/Ortho_resolution + 0.01);
                                    index           = t_col_int + t_row_int*(long)Orthoimagesize.width;
                                    
                                    if(t_col_int >= 0 && t_col_int < Orthoimagesize.width && t_row_int >= 0 && t_row_int < Orthoimagesize.height && index >= 0 && index < data_length_ortho)
                                        result_ortho[index] = value1;
                                }
                            }
                        }
                    }
                    
                    if(impyramid_step > 0)
                        free(pyimg);
                    
                    free(subimage);
                    
                    delete[] data_size;
                }
                
            }
        }
        
        if(args.sensor_type == SB)
            RPCsFree(RPCs);
        free(DEM_value);
        
        WriteGeotiff(OrthoGEOTIFFFilename, result_ortho, Orthoimagesize.width, Orthoimagesize.height, Ortho_resolution, OrthoBoundary[0], OrthoBoundary[3], _param.projection, _param.utm_zone, _param.bHemisphere, 12);
        free(result_ortho);
        
        ET = time(0);
        
        gap = difftime(ET,ST);
        printf("ortho finish(time[m] = %5.2f)!!\n",gap/60.0);
    }
    else
    {
        printf("check overlap area between DEM and image, or match a projection type of input image based on DEM projection by adding '-projection' option\n");
        RPCsFree(RPCs);
    }

    delete proinfo;
}

uint16 *Preprocessing_ortho(const uint8 py_level, CSize *data_size, uint16 *subimg)
{
    int filter_size = pwrtwo(py_level)-1;
    if(filter_size < 3)
        filter_size = 3;
    
    uint16 *pyimg = CreateImagePyramid(subimg,data_size[0],filter_size,(double)1.6);
    
    return pyimg;
}

uint16 *subsetImage_ortho(const int sensor_type, const FrameInfo m_frameinfo, const TransParam transparam, const double *Imageparam, double **RPCs, char *ImageFilename, double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, CSize* subsetsize, bool *ret)
{
    *ret = false;
    
    CSize Imagesize;
    uint16 *leftimage = NULL;
    if(GetImageSize(ImageFilename,&Imagesize))
    {
        long cols[2], rows[2];
        if(GetsubareaImage(sensor_type, m_frameinfo, 0, transparam, Imageparam, RPCs,ImageFilename, Imagesize,subBoundary,minmaxHeight,cols,rows) )
        {
            uint16 type(0);
            leftimage   = Readtiff_T(ImageFilename,&Imagesize,cols,rows,subsetsize,type);
            startpos->m_X   = (double)(cols[0]);
            startpos->m_Y   = (double)(rows[0]);

            *ret        = true;
        }
    }
    
    return leftimage;
}

bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary, const double * const *RPCs, const double gridspace, const CSize DEM_size, const double minX, const double maxY, const TransParam param, const double Ortho_resolution)
{
    const double TopLeft[2] = {minX, maxY};
    double DEMboundary[4];
    DEMboundary[0]  = TopLeft[0];
    DEMboundary[1]  = TopLeft[1]-DEM_size.height*gridspace;
    DEMboundary[2]  = TopLeft[0]+DEM_size.width*gridspace;
    DEMboundary[3]  = TopLeft[1];
    
    printf("DEMBoundary %f\t%f\t%f\t%f\n",DEMboundary[0],DEMboundary[1],DEMboundary[2],DEMboundary[3]);
    
    const double minLon          = -1.15*RPCs[1][2] + RPCs[0][2];
    const double maxLon          =  1.15*RPCs[1][2] + RPCs[0][2];
    const double minLat          = -1.15*RPCs[1][3] + RPCs[0][3];
    const double maxLat          =  1.15*RPCs[1][3] + RPCs[0][3];
    
    printf("lon lat %f\t%f\t%f\t%f\n",minLon,maxLon,minLat,maxLat);
    
    D2DPOINT *XY = NULL;
    D2DPOINT LonLat[4];
    
    LonLat[0].m_X = minLon;
    LonLat[0].m_Y = minLat;
    LonLat[1].m_X = minLon;
    LonLat[1].m_Y = maxLat;
    LonLat[2].m_X = maxLon;
    LonLat[2].m_Y = maxLat;
    LonLat[3].m_X = maxLon;
    LonLat[3].m_Y = minLat;
    
    printf("param %s %d %d\n", param.direction,param.utm_zone,param.projection);
    XY          = wgs2ps(param,4, LonLat);
    
    const double t_minX      = min(min(min(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
    const double t_maxX      = max(max(max(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
    const double t_minY      = min(min(min(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
    const double t_maxY      = max(max(max(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
    
    double ImageBoundary[4];
    ImageBoundary[0]    = floor(t_minX)-1;
    ImageBoundary[1]    = floor(t_minY)-1;
    ImageBoundary[2]    = ceil(t_maxX)+1;
    ImageBoundary[3]    = ceil(t_maxY)+1;
    
    if(param.projection != 1)
    {
        ImageBoundary[0]    = DEMboundary[0];
        ImageBoundary[1]    = DEMboundary[1];
        ImageBoundary[2]    = DEMboundary[2];
        ImageBoundary[3]    = DEMboundary[3];
    }
    
    printf("ImageBoundary %f\t%f\t%f\t%f\n",ImageBoundary[0],ImageBoundary[1],ImageBoundary[2],ImageBoundary[3]);
    
    D2DPOINT DEM_lt,DEM_rb,Image_lt,Image_rb;
    DEM_lt.m_X = DEMboundary[0];
    DEM_lt.m_Y = DEMboundary[3];
    DEM_rb.m_X = DEMboundary[2];
    DEM_rb.m_Y = DEMboundary[1];
    
    Image_lt.m_X = ImageBoundary[0];
    Image_lt.m_Y = ImageBoundary[3];
    Image_rb.m_X = ImageBoundary[2];
    Image_rb.m_Y = ImageBoundary[1];
    
    const bool check_overlap = CheckOverlap(DEM_lt, DEM_rb, Image_lt, Image_rb);
    
    Boundary[0] = (max(DEMboundary[0],ImageBoundary[0]));
    Boundary[1] = (max(DEMboundary[1],ImageBoundary[1]));
    Boundary[2] = (min(DEMboundary[2],ImageBoundary[2]));
    Boundary[3] = (min(DEMboundary[3],ImageBoundary[3]));
    
    Imagesize->height    = ceil(fabs(Boundary[3] - Boundary[1])/Ortho_resolution);
    Imagesize->width     = ceil(fabs(Boundary[2] - Boundary[0])/Ortho_resolution);

    free(XY);

    printf("orthoimage height width %d \t%d\t %f\t%f\n",Imagesize->height,Imagesize->width,fabs(DEMboundary[3] - DEMboundary[1])/Ortho_resolution,fabs(DEMboundary[2] - DEMboundary[0])/Ortho_resolution);
    return check_overlap;
}

bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, const double gridspace, const CSize DEM_size, const double minX, const double maxY, const double Ortho_resolution, const EO Photo, const CAMERA_INFO m_Camera, const RM M)
{
    double TopLeft[2] = {minX, maxY};
    double DEMboundary[4];
    DEMboundary[0]  = TopLeft[0];
    DEMboundary[1]  = TopLeft[1]-DEM_size.height*gridspace;
    DEMboundary[2]  = TopLeft[0]+DEM_size.width*gridspace;
    DEMboundary[3]  = TopLeft[1];

    
    double MSL = 0;
    
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
    
    //printf("photo coord %f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n",top_left.m_X,top_left.m_Y,top_right.m_X,top_right.m_Y,bottom_right.m_X,bottom_right.m_Y,bottom_left.m_X,bottom_left.m_Y);
    
    top_left = ImageToPhoto_single(top_left,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    top_right = ImageToPhoto_single(top_right,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    bottom_right = ImageToPhoto_single(bottom_right,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    bottom_left = ImageToPhoto_single(bottom_left,m_Camera.m_CCDSize,m_Camera.m_ImageSize);
    
    top_left_3D     = GetObjectCoordinate_single(top_left,MSL,Photo, m_Camera, M);
    top_right_3D    = GetObjectCoordinate_single(top_right,MSL,Photo, m_Camera, M);
    bottom_right_3D = GetObjectCoordinate_single(bottom_right,MSL,Photo, m_Camera, M);
    bottom_left_3D  = GetObjectCoordinate_single(bottom_left,MSL,Photo, m_Camera, M);
    
    const double IminX = (top_left_3D.m_X < bottom_left_3D.m_X) ? top_left_3D.m_X : bottom_left_3D.m_X;
    const double IminY = (bottom_left_3D.m_Y < bottom_right_3D.m_Y) ? bottom_left_3D.m_Y : bottom_right_3D.m_Y;
    const double ImaxX = (top_right_3D.m_X > bottom_right_3D.m_X) ? top_right_3D.m_X : bottom_right_3D.m_X;
    const double ImaxY = (top_left_3D.m_Y > top_right_3D.m_Y) ? top_left_3D.m_Y : top_right_3D.m_Y;
    
    double ImageBoundary[4];
    
    ImageBoundary[0] =  floor(IminX);
    ImageBoundary[1] =  floor(IminY);
    ImageBoundary[2] =  ceil(ImaxX);
    ImageBoundary[3] =  ceil(ImaxY);
    //printf("photo coord %f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n",top_left.m_X,top_left.m_Y,top_right.m_X,top_right.m_Y,bottom_right.m_X,bottom_right.m_Y,bottom_left.m_X,bottom_left.m_Y);
    //printf("ImageBoundary %f\t%f\t%f\t%f\n",ImageBoundary[0],ImageBoundary[1],ImageBoundary[2],ImageBoundary[3]);
    //printf("DEMboundary %f\t%f\t%f\t%f\n",DEMboundary[0],DEMboundary[1],DEMboundary[2],DEMboundary[3]);
    
    D2DPOINT DEM_lt,DEM_rb,Image_lt,Image_rb;
    DEM_lt.m_X = DEMboundary[0];
    DEM_lt.m_Y = DEMboundary[3];
    DEM_rb.m_X = DEMboundary[2];
    DEM_rb.m_Y = DEMboundary[1];
    
    Image_lt.m_X = ImageBoundary[0];
    Image_lt.m_Y = ImageBoundary[3];
    Image_rb.m_X = ImageBoundary[2];
    Image_rb.m_Y = ImageBoundary[1];
    
    const bool check_overlap = CheckOverlap(DEM_lt, DEM_rb, Image_lt, Image_rb);
    
    Boundary[0] = (max(DEMboundary[0],ImageBoundary[0]));
    Boundary[1] = (max(DEMboundary[1],ImageBoundary[1]));
    Boundary[2] = (min(DEMboundary[2],ImageBoundary[2]));
    Boundary[3] = (min(DEMboundary[3],ImageBoundary[3]));
    
    Imagesize->height    = ceil(fabs(Boundary[3] - Boundary[1])/Ortho_resolution);
    Imagesize->width     = ceil(fabs(Boundary[2] - Boundary[0])/Ortho_resolution);
    
    
    printf("orthoimage height width %d \t%d\t %f\t%f\n",Imagesize->height,Imagesize->width,fabs(DEMboundary[3] - DEMboundary[1])/Ortho_resolution,fabs(DEMboundary[2] - DEMboundary[0])/Ortho_resolution);
    
    return true;
}

