//
//  CoordConversion.hpp
//  
//
//  Created by Myoung-Jong Noh on 3/30/20.
//

#ifndef CoordConversion_hpp
#define CoordConversion_hpp

#include <cmath>
#include <stdlib.h>

#include "Typedefine.hpp"

//WGS to PS, WGS to UTM
D2DPOINT *wgs2ps(TransParam _param, int _numofpts, D2DPOINT *_wgs);
D2DPOINT wgs2ps_single(TransParam _param, D2DPOINT _wgs);
D3DPOINT *wgs2ps_3D(TransParam _param, int _numofpts, D3DPOINT *_wgs);
D2DPOINT *ps2wgs(const TransParam _param,const long int _numofpts,const D2DPOINT *_ps);
D2DPOINT ps2wgs_single(const TransParam _param, const D2DPOINT _ps);
D3DPOINT *ps2wgs_3D(TransParam _param, int _numofpts, D3DPOINT *_ps);
D3DPOINT *ps2wgs_3D_vector(TransParam _param, int _numofpts, vector<D3DPOINT> &_ps);

//RPC conversion : image to object
D2DPOINT* GetObjectToImageRPC(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const uint16 _numofpts, D3DPOINT *_GP);
static D2DPOINT GetObjectToImageRPC_single(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP)
{
    double L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    double P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    double H       = (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
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
    
    double Coeff[4];
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

    double Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
    double Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample

    double deltaP      = _imageparam[0];
    double deltaR      = _imageparam[1];

    D2DPOINT IP;
    IP.m_Y      = deltaP + Line;
    IP.m_X      = deltaR + Samp;

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

static D2DPOINT GetObjectToImageRPC_single_noBias(const double * const *_rpc, const uint8 _numofparam, D3DPOINT _GP)
{
    double L       = (_GP.m_X - _rpc[0][2])/_rpc[1][2];
    double P       = (_GP.m_Y - _rpc[0][3])/_rpc[1][3];
    double H       = (_GP.m_Z - _rpc[0][4])/_rpc[1][4];
    
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
    
    double Coeff[4];
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

    double Line     = ((Coeff[0]/Coeff[1])*_rpc[1][0] + _rpc[0][0]); //Line
    double Samp     = ((Coeff[2]/Coeff[3])*_rpc[1][1] + _rpc[0][1]); //Sample

    D2DPOINT IP;
    IP.m_Y      = Line;
    IP.m_X      = Samp;

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

D2DPOINT GetObjectToImageRPC_single_mpp(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP);

//Geotiff conversion : image to geo-coord
D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres);
D2DPOINT GetObjectToImage_single(uint16 _numofpts, D2DPOINT _GP, double *boundary, double imageres);

//Pyramid coord conversion
D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
static D2DPOINT OriginalToPyramid_single(const D2DPOINT InCoord, const D2DPOINT Startpos, const uint8 Pyramid_step)
{
    D2DPOINT out;

    out.m_X      = (InCoord.m_X/pwrtwo(Pyramid_step)) - Startpos.m_X;
    out.m_Y      = (InCoord.m_Y/pwrtwo(Pyramid_step)) - Startpos.m_Y;

    return out;
    
}
D2DPOINT* PyramidToOriginal(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);

//Air-borne frame sensor coord conversion : photo to object
RM MakeRotationMatrix(double o, double p, double k);
D2DPOINT *GetPhotoCoordinate(D3DPOINT *A, EO Photo, int _numofpts, CAMERA_INFO Camera, RM M);
D3DPOINT *GetObjectCoordinate(D2DPOINT *a, double z,EO Photo, int _numofpts, CAMERA_INFO Camera, RM M);
D2DPOINT *PhotoToImage(D2DPOINT *_photo, int _numofpts, float _CCDSize, CSize _imgsize);
D2DPOINT *ImageToPhoto(D2DPOINT *_image, int _numofpts, float _CCDSize, CSize _imgsize);

D2DPOINT GetPhotoCoordinate_single(const D3DPOINT A, const EO Photo, const CAMERA_INFO Camera, const RM M);
D3DPOINT GetObjectCoordinate_single(D2DPOINT a, double z,EO Photo, CAMERA_INFO Camera, RM M);
D2DPOINT PhotoToImage_single(const D2DPOINT _photo, const double _CCDSize, const CSize _imgsize);
D2DPOINT ImageToPhoto_single(D2DPOINT _image, float _CCDSize, CSize _imgsize);

#endif /* CoordConversion_hpp */
