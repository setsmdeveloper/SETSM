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

//RPC conversion : image to object
D2DPOINT* GetObjectToImageRPC(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, const uint16 _numofpts, D3DPOINT *_GP);
D2DPOINT GetObjectToImageRPC_single(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP);
D2DPOINT GetObjectToImageRPC_single_mpp(const double * const *_rpc, const uint8 _numofparam, const double *_imageparam, D3DPOINT _GP);

//Geotiff conversion : image to geo-coord
D2DPOINT* GetObjectToImage(uint16 _numofpts, D2DPOINT *_GP, double *boundary, double imageres);
D2DPOINT GetObjectToImage_single(uint16 _numofpts, D2DPOINT _GP, double *boundary, double imageres);

//Pyramid coord conversion
D2DPOINT* OriginalToPyramid(uint16 numofpts,D2DPOINT* InCoord, D2DPOINT Startpos, uint8 Pyramid_step);
D2DPOINT OriginalToPyramid_single(const D2DPOINT InCoord, const D2DPOINT Startpos, const uint8 Pyramid_step);
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
