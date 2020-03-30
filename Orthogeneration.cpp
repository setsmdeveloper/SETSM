//
//  Orthogeneration.cpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#include "Orthogeneration.hpp"


//orthogeneration
void orthogeneration(TransParam _param, ARGINFO args, char *ImageFilename, char *DEMFilename, char *Outputpath,int pair,int DEM_divide, double **Imageparams)
{
    if(args.RA_only) {
        return;
    }
    char DEM_header[500];
    char RPCFilename[500];
    char OrthoFilename[500];
    char OrthoGEOTIFFFilename[500];
    char Ortho_header[500];
    
    time_t ST = 0, ET = 0;
    double gap;
    
    char* tmp_chr;
    char* tmp_no_ext;
    double DEM_resolution, Image_resolution, Ortho_resolution;
    double OrthoGridFactor;
    double resolution[3];
    int impyramid_step;
    double** RPCs;
    double row_grid_size, col_grid_size, product_grid_size;
    bool Hemisphere = false;
    double DEM_minX, DEM_maxY;
    double minLat,minLon;
    CSize DEM_size, Image_size;
    TransParam param = _param;
    FrameInfo m_frameinfo;
    m_frameinfo.m_Camera.m_focalLength  = 0;
    m_frameinfo.m_Camera.m_CCDSize      = 0;
    
    
    ST = time(0);
    
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
    
    tmp_chr = remove_ext(ImageFilename);
    
    char *Ifilename  = SetOutpathName(ImageFilename);
    tmp_no_ext = remove_ext(Ifilename);
    
    sprintf(RPCFilename,"%s.xml",tmp_chr);
    
    FILE* fid_xml;
    fid_xml = fopen(RPCFilename,"r");
    if(!fid_xml)
    {
        sprintf(RPCFilename,"%s.XML",tmp_chr);
        FILE* fid_XML;
        fid_XML = fopen(RPCFilename,"r");
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
    
    tmp_chr = remove_ext(DEMFilename);
    sprintf(DEM_header,"%s.hdr",tmp_chr);
    
    //if(pair == 1)
    {
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
    }
    /*else
    {
        if(DEM_divide == 0)
        {
            sprintf(OrthoFilename, "%s/%s_ortho_image2.raw", Outputpath,args.Outputpath_name);
            sprintf(Ortho_header, "%s/%s_ortho_image2.hdr", Outputpath, args.Outputpath_name);
            sprintf(OrthoGEOTIFFFilename, "%s/%s_ortho_image2.tif", Outputpath, args.Outputpath_name);
        }
        else
        {
            sprintf(OrthoFilename, "%s/%s_%d_ortho_image2.raw", Outputpath,args.Outputpath_name,DEM_divide);
            sprintf(Ortho_header, "%s/%s_%d_ortho_image2.hdr", Outputpath, args.Outputpath_name,DEM_divide);
            sprintf(OrthoGEOTIFFFilename, "%s/%s_%d_ortho_image2.tif", Outputpath, args.Outputpath_name,DEM_divide);
        }
    }
    */
    printf("image = %s\n",ImageFilename);
    printf("rpc = %s\n",RPCFilename);
    printf("save = %s\n",Outputpath);
    printf("DEM = %s\n",DEMFilename);
    printf("DEM hdr= %s\n",DEM_header);
    printf("ortho = %s\n",OrthoFilename);
    printf("ortho hdr= %s\n",Ortho_header);
    printf("ortho geotiff = %s\n", OrthoGEOTIFFFilename);
    
    
    Image_resolution = 0.5;
    
    // load RPCs info from xml file
    if(args.sensor_type == SB)
    {
        if(args.sensor_provider == DG)
        {
            RPCs            = OpenXMLFile_ortho(RPCFilename, &row_grid_size, &col_grid_size,&product_grid_size);
        }
        else if(args.sensor_provider == PL)
        {
            RPCs            = OpenXMLFile_Pleiades(RPCFilename);
        }
        else if(args.sensor_provider == PT)
        {
            RPCs            = OpenXMLFile_Planet(RPCFilename);
        }
    }
    else
    {
        FILE *pFile;
        pFile           = fopen(RPCFilename,"r");
        
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

    
    // load RPCs info from xml file
    /*RPCs          = OpenXMLFile_ortho(RPCFilename, &row_grid_size, &col_grid_size);
    
      Image_resolution = (int)(((row_grid_size + col_grid_size)/2.0)*10 + 0.5)/10.0;
    
      if (Image_resolution < 0.75) {
      Image_resolution = 0.5;
      }
      else if(Image_resolution < 1.25)
      Image_resolution = 1.0;
    */
    printf("Image resolution %f\n",Image_resolution);
    /*
    if(args.sensor_type == SB)
    {
        minLat = RPCs[0][3];
        minLon = RPCs[0][2];
        
        //set projection conversion info.
        param.utm_zone   = args.utm_zone;
        _param.utm_zone  = args.utm_zone;
        
        SetTransParam(minLat, minLon, &Hemisphere, &_param);
        SetTransParam(minLat, minLon, &Hemisphere, &param);
        
        param.projection = _param.projection;
    }
     */
    printf("Hemis projection %d %d\n",Hemisphere, param.projection);
    printf("param %s %d %d\n", param.direction,param.zone,param.projection);
    
    // load DEM infor from geotiff file.
    DEM_size = ReadGeotiff_info(DEMFilename, &DEM_minX, &DEM_maxY, &DEM_resolution);
    
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
    
    impyramid_step  = ceil(log(resolution[2]/resolution[1])/log(2));
    printf("impyramid_step %d\n ortho_Resolution %f\n",impyramid_step,Ortho_resolution);
    
    CSize Orthoimagesize;
    double OrthoBoundary[4];
    bool check_overlap;
    // set generated orthoimage info by comparing DEM info
    if(args.sensor_type == SB)
        check_overlap = SetOrthoBoundary_ortho(&Orthoimagesize, OrthoBoundary, RPCs, DEM_resolution, DEM_size, DEM_minX, DEM_maxY, param, Ortho_resolution);
    else
    {
        GetImageSize_ortho(ImageFilename,&m_frameinfo.m_Camera.m_ImageSize);
        check_overlap = SetDEMBoundary_ortho_photo(&Orthoimagesize, OrthoBoundary,DEM_resolution, DEM_size, DEM_minX, DEM_maxY, Ortho_resolution, m_frameinfo.Photoinfo[0], m_frameinfo.m_Camera, m_frameinfo.Photoinfo[0].m_Rm);
    }
    
    if(check_overlap)
    {
        // set saving pointer for orthoimage
        uint16 *result_ortho    = (uint16*)calloc((long)Orthoimagesize.width*(long)Orthoimagesize.height,sizeof(uint16));
        
        float *DEM_value    = NULL;
        
        double minmaxHeight[2];
        minmaxHeight[0]     = 99999;
        minmaxHeight[1]     = -99999;
        // load DEM value;
        DEM_value = GetDEMValue(DEMFilename, DEM_size);
        
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
        
        double subfactor                      = pow(4-impyramid_step,2.0);
        if(subfactor <= 1)
            subfactor                       = 1;
        int sub_height                  = ceil(Orthoimagesize.height/subfactor);
        int sub_width                   = ceil(Orthoimagesize.width/subfactor);
        int height_interval             = Ortho_resolution*sub_height;
        int width_interval              = Ortho_resolution*sub_width;
        
        int buffer_x                    = Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
        int buffer_y                    = Ortho_resolution*5*ceil(1.0/OrthoGridFactor);
        
        int tile_count;
        
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
        
        for(i=0; i < subfactor;i++)
        {
            for(j=0; j<subfactor;j++)
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
                char subsetImageFilename[500];
                char subsetImageFilename_hdr[500];
                CSize subsetsize;
                uint16 *subimage;
                bool check_subsetImage = false;
                
                subimage = subsetImage_ortho(args.sensor_type, m_frameinfo, param, RPCs, ImageFilename,
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
                    startpos.m_X       = (double)(startpos_ori.m_X/pwrtwo(impyramid_step));      startpos.m_Y       = (double)(startpos_ori.m_Y/pwrtwo(impyramid_step));
                    
                    uint16 *pyimg;
                    
                    if(impyramid_step > 0)
                        pyimg = Preprocessing_ortho(ori_impyramid_step,data_size,subimage);
                    
                    Image_size  = data_size[impyramid_step];
                    
                    //printf("Image_size %d\t%d\n",Image_size.width,Image_size.height);
                    
                    int col_size    = (int)((X_size[1] - X_size[0])/Ortho_resolution + 0.5);
                    int row_size    = (int)((Y_size[1] - Y_size[0])/Ortho_resolution + 0.5);
                    
    #pragma omp parallel for schedule(guided)
                    for(int count = 0; count < col_size*row_size ; count++)
                    {
                        double row  = ((int)(floor(count/col_size)))*Ortho_resolution + Y_size[0];
                        double col  = (count % col_size)*Ortho_resolution + X_size[0];

                        long int index1,index2,index3, index4;
                        double t_col, t_row;
                        long int t_col_int, t_row_int;
                        double dcol,drow;
                        
                        t_col       = (col - DEM_minX)/DEM_resolution;
                        t_row       = (DEM_maxY - row)/DEM_resolution;
                        
                        t_col_int   = (long int)(t_col + 0.01);
                        t_row_int   = (long int)(t_row + 0.01);
                        
                        dcol        = t_col - t_col_int;
                        drow        = t_row - t_row_int;
                        
                        if(t_col_int >= 0 && t_col_int +1 < DEM_size.width && t_row_int >= 0 && t_row_int +1 < DEM_size.height)
                        {
                            double value1, value2, value3, value4, value;
                            
                            index1  = (t_col_int   ) + (t_row_int   )*(long)DEM_size.width;
                            index2  = (t_col_int +1) + (t_row_int   )*(long)DEM_size.width;
                            index3  = (t_col_int   ) + (t_row_int +1)*(long)DEM_size.width;
                            index4  = (t_col_int +1) + (t_row_int +1)*(long)DEM_size.width;
                            
                            value1      = DEM_value[index1];
                            value2      = DEM_value[index2];
                            value3      = DEM_value[index3];
                            value4      = DEM_value[index4];
                            
                            
                            value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                + value3*(1-dcol)*drow + value4*dcol*drow;
                            
                            D3DPOINT object;
                            D2DPOINT objectXY;
                            //double imageparam[2] = {0.};
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
                                    D2DPOINT wgsPt = ps2wgs_single(param, objectXY);
                                   
                                    object.m_X  = wgsPt.m_X;
                                    object.m_Y  = wgsPt.m_Y;
                                    image = GetObjectToImageRPC_single_ortho(RPCs, 2, imageparam, object);
                                }
                                else
                                {
                                    D2DPOINT photo  = GetPhotoCoordinate_single(object,m_frameinfo.Photoinfo[0],m_frameinfo.m_Camera,m_frameinfo.Photoinfo[0].m_Rm);
                                    image = PhotoToImage_single(photo, m_frameinfo.m_Camera.m_CCDSize, m_frameinfo.m_Camera.m_ImageSize);
                                    
                                    //printf("done photoToImage %f\t%f\n",image.m_X,image.m_Y);
                                }
                                
                                
                                temp_pt     = OriginalToPyramid_single_ortho(image, startpos, impyramid_step);
                                
                                t_col       = temp_pt.m_X;
                                t_row       = temp_pt.m_Y;
                                
                                t_col_int   = (long int)(t_col + 0.01);
                                t_row_int   = (long int)(t_row + 0.01);
                                
                                dcol        = t_col - t_col_int;
                                drow        = t_row - t_row_int;
                                
                                //printf("temp_pt %f\t%f\n",temp_pt.m_X,temp_pt.m_Y);
                                if(t_col_int >= 0 && t_col_int +1 < Image_size.width && t_row_int >= 0 && t_row_int +1 < Image_size.height
                                   && (t_col_int +1) + (t_row_int +1)*(long)Image_size.width < (long)Image_size.width*(long)Image_size.height)
                                {
                                    //printf("inside of image value\n");
                                    double value1, value2, value3, value4, value;
                                    long int index;
                                    index1  = (t_col_int   ) + (t_row_int   )*(long)Image_size.width;
                                    index2  = (t_col_int +1) + (t_row_int   )*(long)Image_size.width;
                                    index3  = (t_col_int   ) + (t_row_int +1)*(long)Image_size.width;
                                    index4  = (t_col_int +1) + (t_row_int +1)*(long)Image_size.width;
                                    
                                    if(impyramid_step > 0)
                                    {
                                        value1      = pyimg[index1];
                                        value2      = pyimg[index2];
                                        value3      = pyimg[index3];
                                        value4      = pyimg[index4];
                                    }
                                    else {
                                        value1      = subimage[index1];
                                        value2      = subimage[index2];
                                        value3      = subimage[index3];
                                        value4      = subimage[index4];
                                    }

                                    value       = value1*(1-dcol)*(1-drow) + value2*dcol*(1-drow)
                                        + value3*(1-dcol)*drow + value4*dcol*drow;
                                    
                                    t_col_int       = (long int)((col - OrthoBoundary[0])/Ortho_resolution + 0.01);
                                    t_row_int       = (long int)((OrthoBoundary[3] - row)/Ortho_resolution + 0.01);
                                    index           = t_col_int + t_row_int*(long)Orthoimagesize.width;
                                    
                                    if(t_col_int >= 0 && t_col_int < Orthoimagesize.width && t_row_int >= 0 && t_row_int < Orthoimagesize.height && index >= 0 && index < (long)Orthoimagesize.width*(long)Orthoimagesize.height)
                                    {
                                        //printf("inside of ortho\n");
                                        result_ortho[index] = value1;
                                    }
                                }
                            }
                        }
                    }
                    
                    if(impyramid_step > 0)
                        free(pyimg);
                    
                    free(subimage);
                }
                
            }
        }
        free(RPCs);
        free(DEM_value);
        
        WriteGeotiff(OrthoGEOTIFFFilename, result_ortho, Orthoimagesize.width, Orthoimagesize.height, Ortho_resolution, OrthoBoundary[0], OrthoBoundary[3], _param.projection, _param.zone, _param.bHemisphere, 12);
        free(result_ortho);
        
        ET = time(0);
        
        gap = difftime(ET,ST);
        printf("ortho finish(time[m] = %5.2f)!!\n",gap/60.0);
    }
    else
    {
        printf("check overlap area between DEM and image, or match a projection type of input image based on DEM projection by adding '-projection' option\n");
        free(RPCs);
    }
}

D2DPOINT OriginalToPyramid_single_ortho(D2DPOINT InCoord, D2DPOINT Startpos, uint8 Pyramid_step)
{
    D2DPOINT out;
    
    out.m_X      = (InCoord.m_X/pwrtwo(Pyramid_step)) - Startpos.m_X;
    out.m_Y      = (InCoord.m_Y/pwrtwo(Pyramid_step)) - Startpos.m_Y;
    
    return out;
    
}

uint16 *Preprocessing_ortho(uint8 py_level, CSize *data_size, uint16 *subimg)
{
    uint16 *pyimg;
    int filter_size = pwrtwo(py_level)-1;
    if(filter_size < 3)
        filter_size = 3;
    
    //printf("level %d\tfilter size %d\n",py_level,filter_size);
    
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
    
    data_size[0].height     = subsetsize.height;
    data_size[0].width      = subsetsize.width;
    for(i=0;i<level;i++)
    {
        data_size[i+1].width  = data_size[i].width/2;
        data_size[i+1].height = data_size[i].height/2;
    }
}


uint16 *subsetImage_ortho(int sensor_type, FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename,
                          double *subBoundary, double *minmaxHeight, D2DPOINT *startpos, char *subsetImage, CSize* subsetsize, bool *ret)
{
    *ret = false;
    
    CSize Imagesize;
    uint16 *leftimage = NULL;
    if(GetImageSize_ortho(ImageFilename,&Imagesize))
    {
        int cols[2], rows[2];
        if(GetsubareaImage_ortho(sensor_type,m_frameinfo,transparam,RPCs,ImageFilename,&Imagesize,subBoundary,minmaxHeight,cols,rows) )
        {
            
            
            leftimage   = Readtiff_ortho(ImageFilename,Imagesize,cols,rows,subsetsize);
            
            
            startpos->m_X   = (double)(cols[0]);
            startpos->m_Y   = (double)(rows[0]);

            *ret        = true;
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
        *Imagesize = ReadGeotiff_info(filename, NULL, NULL, NULL);
        ret = true;
    }
    else if(!strcmp("bin",ext+1))
    {
        char *tmp;
        tmp = remove_ext(filename);
        sprintf(tmp,"%s.hdr",tmp);
        *Imagesize = Envihdr_reader_ortho(tmp);
        free(tmp);
        
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

bool GetsubareaImage_ortho(int sensor_type, FrameInfo m_frameinfo, TransParam transparam, double **RPCs, char *ImageFilename, CSize *Imagesize,
                           double *subBoundary, double *minmaxHeight, int *cols, int *rows)
{
    bool ret = false;
    
    if(GetImageSize_ortho(ImageFilename,Imagesize))
    {
        int i;
        
        D3DPOINT t_pts[8];
        
        D2DPOINT *ImageCoord;
        int buffer, null_buffer;
        double ImageParam[2] = {0.};
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
        
        t_pts[0].m_Z    = minmaxHeight[0];
        t_pts[1].m_Z    = minmaxHeight[0];
        t_pts[2].m_Z    = minmaxHeight[1];
        t_pts[3].m_Z    = minmaxHeight[1];
        t_pts[4].m_Z    = minmaxHeight[0];
        t_pts[5].m_Z    = minmaxHeight[0];
        t_pts[6].m_Z    = minmaxHeight[1];
        t_pts[7].m_Z    = minmaxHeight[1];
        
        if(sensor_type == SB)
        {
            D3DPOINT *t_pts1;
            t_pts1          = ps2wgs_3D(transparam,8,t_pts);
            ImageCoord      = GetObjectToImageRPC(RPCs, 2, ImageParam, 8, t_pts1);
            
            free(t_pts1);
        }
        else
        {
            D2DPOINT *t_image;
            t_image     = GetPhotoCoordinate(t_pts, m_frameinfo.Photoinfo[0], 8, m_frameinfo.m_Camera, m_frameinfo.Photoinfo[0].m_Rm);
            ImageCoord  = PhotoToImage(t_image,8,m_frameinfo.m_Camera.m_CCDSize,m_frameinfo.m_Camera.m_ImageSize);
            
            free(t_image);
        }
        
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
        }
        
        buffer              = 200;
        cols[0]             = (int)(ceil(minX)-buffer);
        cols[1]             = (int)(ceil(maxX)+buffer);
        rows[0]             = (int)(ceil(minY)-buffer);
        rows[1]             = (int)(ceil(maxY)+buffer);
        
        null_buffer         = 1;
        // Null pixel value remove
        if(cols[0]          <= null_buffer)
            cols[0]         = null_buffer;
        if(rows[0]          <= null_buffer)
            rows[0]         = null_buffer;
        if(cols[0]          > Imagesize->width - null_buffer)
            cols[0]         = Imagesize->width - null_buffer;
        if(rows[0]          > Imagesize->height - null_buffer)
            rows[0]         = Imagesize->height - null_buffer;
        
        if(cols[1]          <= null_buffer)
            cols[1]         = null_buffer;
        if(rows[1]          <= null_buffer)
            rows[1]         = null_buffer;
        if(cols[1]          > Imagesize->width - null_buffer)
            cols[1]         = Imagesize->width - null_buffer;
        if(rows[1]          > Imagesize->height - null_buffer)
            rows[1]         = Imagesize->height - null_buffer;
        
        
        free(ImageCoord);
        
        ret = true;
    }
    
    return ret;
}

D2DPOINT* GetObjectToImageRPC_ortho(double **_rpc, uint8 _numofparam, double *_imageparam, uint16 _numofpts, D3DPOINT *_GP)
{
    D2DPOINT *IP;
    
    IP      = (D2DPOINT*)malloc(sizeof(D2DPOINT)*_numofpts);

#pragma omp parallel for schedule(guided)
    for(int i=0;i<_numofpts;i++)
    {
        double L, P, H, Line, Samp;
        double deltaP = 0.0, deltaR = 0.0;
        double Coeff[4];
        
        L       = (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
        P       = (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
        H       = (_GP[i].m_Z - _rpc[0][4])/_rpc[1][4];
        
        if(L < -10.0 || L > 10.0)
        {
            if(_GP[i].m_X > 0)
                _GP[i].m_X = _GP[i].m_X - 360;
            else
                _GP[i].m_X = _GP[i].m_X + 360;
            
            L       = (_GP[i].m_X - _rpc[0][2])/_rpc[1][2];
        }
        
        if(P < -10.0 || P > 10.0)
        {
            if(_GP[i].m_Y > 0)
                _GP[i].m_Y = _GP[i].m_Y - 360;
            else
                _GP[i].m_Y = _GP[i].m_Y + 360;
            
            P       = (_GP[i].m_Y - _rpc[0][3])/_rpc[1][3];
        }

        for(int j=0;j<4;j++)
        {
            Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
                + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
                + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
                + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
                + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
                + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
                + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
        }
        
        Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
        Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
        
        switch(_numofparam)
        {
        case 2:
            deltaP      = _imageparam[0];
            deltaR      = _imageparam[1];
            break;
        case 6:
            deltaP      = _imageparam[0] + _imageparam[1]*Samp + _imageparam[2]*Line;
            deltaR      = _imageparam[3] + _imageparam[4]*Samp + _imageparam[5]*Line;
            break;
        }
        
        IP[i].m_Y       = deltaP + Line;
        IP[i].m_X       = deltaR + Samp;
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
        //printf("tif open\n");
        tif  = TIFFOpen(filename,"r");
        check_ftype = 1;
        //printf("tif open end\n");
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
            tsize_t scanline;
            tdata_t buf;
            uint16 s,nsamples;
            
            // scanline read
            data_size->width    = cols[1] - cols[0];
            data_size->height   = rows[1] - rows[0];
            
            long int data_length = data_size->height*data_size->width;
            out             = (uint16*)malloc(sizeof(uint16)*data_length);
            scanline        = TIFFScanlineSize(tif);
            
            buf             = _TIFFmalloc(scanline);
            
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
            
            starttileL      = (int)(rows[0]/tileL);
            start_row       = starttileL*tileL;
            end_row         = ((int)(rows[1]/tileL)+1)*tileL;
            if(end_row > Imagesize.height)
                end_row = Imagesize.height;
            
            starttileW      = (int)(cols[0]/tileW);
            start_col       = starttileW*tileW;
            end_col         = ((int)(cols[1]/tileW)+1)*tileW;
            if(end_col > Imagesize.width)
                end_col = Imagesize.width;
            
            
            cols[0]         = start_col;
            cols[1]         = end_col;
            rows[0]         = start_row;
            rows[1]         = end_row;
            
            data_size->width = end_col - start_col;
            data_size->height= end_row - start_row;
            
            long int data_length = data_size->height*data_size->width;
            
            out             = (uint16*)malloc(sizeof(uint16)*data_length);
            
            buf             = _TIFFmalloc(TIFFTileSize(tif));
            
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
        data_size->width    = cols[1] - cols[0];
        data_size->height   = rows[1] - rows[0];
        
        long int data_length = data_size->height*data_size->width;
        
        out             = (uint16*)malloc(sizeof(uint16)*data_length);
        
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
    double Coeff[4];
    deltaP = 0.0;
    deltaR = 0.0;
    
    L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    H       = (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
    if(L < -10.0 || L > 10.0)
    {
        if(_GP.m_X > 0)
            _GP.m_X = _GP.m_X - 360;
        else
            _GP.m_X = _GP.m_X + 360;
        
        L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    }
    
    if(P < -10.0 || P > 10.0)
    {
        if(_GP.m_Y > 0)
            _GP.m_Y = _GP.m_Y - 360;
        else
            _GP.m_Y = _GP.m_Y + 360;
        
        P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    }

    for(j=0;j<4;j++)
    {
        Coeff[j]    = _rpc[j+2][0]*1.0          + _rpc[j+2][1]*L            + _rpc[j+2][2]*P
            + _rpc[j+2][3]*H            + _rpc[j+2][4]*L*P          + _rpc[j+2][5]*L*H
            + _rpc[j+2][6]*P*H          + _rpc[j+2][7]*L*L          + _rpc[j+2][8]*P*P
            + _rpc[j+2][9]*H*H          + _rpc[j+2][10]*(P*L)*H     + _rpc[j+2][11]*(L*L)*L
            + _rpc[j+2][12]*(L*P)*P     + _rpc[j+2][13]*(L*H)*H     + _rpc[j+2][14]*(L*L)*P
            + _rpc[j+2][15]*(P*P)*P     + _rpc[j+2][16]*(P*H)*H     + _rpc[j+2][17]*(L*L)*H
            + _rpc[j+2][18]*(P*P)*H     + _rpc[j+2][19]*(H*H)*H;
    }
    
    Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
    Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample
    
    switch(_numofparam)
    {
    case 2:
        deltaP      = _imageparam[0];
        deltaR      = _imageparam[1];
        break;
    case 6:
        deltaP      = _imageparam[0] + _imageparam[1]*Samp + _imageparam[2]*Line;
        deltaR      = _imageparam[3] + _imageparam[4]*Samp + _imageparam[5]*Line;
        break;
    }

    IP.m_Y      = deltaP + Line;
    IP.m_X      = deltaR + Samp;

    return IP;
}

bool SetOrthoBoundary_ortho(CSize *Imagesize, double *Boundary,
                            double **RPCs, double gridspace, CSize DEM_size, double minX, double maxY, TransParam param, double Ortho_resolution)
{
    double TopLeft[2];
    TopLeft[0]  = minX;
    TopLeft[1]  = maxY;
    double DEMboundary[4];
    DEMboundary[0]  = TopLeft[0];
    DEMboundary[1]  = TopLeft[1]-DEM_size.height*gridspace;
    DEMboundary[2]  = TopLeft[0]+DEM_size.width*gridspace;
    DEMboundary[3]  = TopLeft[1];
    
    printf("DEMBoundary %f\t%f\t%f\t%f\n",DEMboundary[0],DEMboundary[1],DEMboundary[2],DEMboundary[3]);
    
    double minLon, maxLon, minLat, maxLat;
    minLon          = -1.15*RPCs[1][2] + RPCs[0][2];
    maxLon          =  1.15*RPCs[1][2] + RPCs[0][2];
    minLat          = -1.15*RPCs[1][3] + RPCs[0][3];
    maxLat          =  1.15*RPCs[1][3] + RPCs[0][3];
    
    printf("lon lat %f\t%f\t%f\t%f\n",minLon,maxLon,minLat,maxLat);
    
    D2DPOINT *XY;
    D2DPOINT LonLat[4];
    double t_minX, t_maxX, t_minY, t_maxY;
    
    LonLat[0].m_X = minLon;
    LonLat[0].m_Y = minLat;
    LonLat[1].m_X = minLon;
    LonLat[1].m_Y = maxLat;
    LonLat[2].m_X = maxLon;
    LonLat[2].m_Y = maxLat;
    LonLat[3].m_X = maxLon;
    LonLat[3].m_Y = minLat;
    
    printf("param %s %d %d\n", param.direction,param.zone,param.projection);
    XY          = wgs2ps(param,4, LonLat);
    
    t_minX      = min(min(min(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
    t_maxX      = max(max(max(XY[0].m_X,XY[1].m_X),XY[2].m_X),XY[3].m_X);
    t_minY      = min(min(min(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
    t_maxY      = max(max(max(XY[0].m_Y,XY[1].m_Y),XY[2].m_Y),XY[3].m_Y);
    
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
    
    bool check_overlap = CheckOverlap(DEM_lt, DEM_rb, Image_lt, Image_rb);
    
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




bool SetDEMBoundary_ortho_photo(CSize *Imagesize, double *Boundary, double gridspace, CSize DEM_size, double minX, double maxY, double Ortho_resolution, EO Photo, CAMERA_INFO m_Camera, RM M)
{
    double TopLeft[2];
    TopLeft[0]  = minX;
    TopLeft[1]  = maxY;
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
    
    double IminX = (top_left_3D.m_X < bottom_left_3D.m_X) ? top_left_3D.m_X : bottom_left_3D.m_X;
    double IminY = (bottom_left_3D.m_Y < bottom_right_3D.m_Y) ? bottom_left_3D.m_Y : bottom_right_3D.m_Y;
    double ImaxX = (top_right_3D.m_X > bottom_right_3D.m_X) ? top_right_3D.m_X : bottom_right_3D.m_X;
    double ImaxY = (top_left_3D.m_Y > top_right_3D.m_Y) ? top_left_3D.m_Y : top_right_3D.m_Y;
    
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
    
    bool check_overlap = CheckOverlap(DEM_lt, DEM_rb, Image_lt, Image_rb);
    
    Boundary[0] = (max(DEMboundary[0],ImageBoundary[0]));
    Boundary[1] = (max(DEMboundary[1],ImageBoundary[1]));
    Boundary[2] = (min(DEMboundary[2],ImageBoundary[2]));
    Boundary[3] = (min(DEMboundary[3],ImageBoundary[3]));
    
    Imagesize->height    = ceil(fabs(Boundary[3] - Boundary[1])/Ortho_resolution);
    Imagesize->width     = ceil(fabs(Boundary[2] - Boundary[0])/Ortho_resolution);
    
    
    printf("orthoimage height width %d \t%d\t %f\t%f\n",Imagesize->height,Imagesize->width,fabs(DEMboundary[3] - DEMboundary[1])/Ortho_resolution,fabs(DEMboundary[2] - DEMboundary[0])/Ortho_resolution);
    
    return true;
}

double** OpenXMLFile_ortho(char* _filename, double* gsd_r, double* gsd_c, double* gsd)
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
    
    pFile           = fopen(_filename,"r");
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
                out[2] = (double*)malloc(sizeof(double)*20);
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
                out[3] = (double*)malloc(sizeof(double)*20);
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
                out[4] = (double*)malloc(sizeof(double)*20);
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
                out[5] = (double*)malloc(sizeof(double)*20);
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
