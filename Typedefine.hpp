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
#ifndef _Typedefine_H_
#define _Typedefine_H_

//uint definition
#include "tiff.h"
#include "tiffio.h"
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define PI 3.141592653589793
#define DegToRad PI/180
#define RadToDeg 180/PI
#define UMToMM 0.001
#define MMToUM 1000
#define MaxImages 100
#define MaxNCC 700
#define Nodata -9999
#define Roh_min -30.0
#define Roh_max 30.0
#ifndef bool
//#define bool unsigned char
#define true 0x1
#define false 0x0
#define SQ(x)         (x) * (x)
#define SWAP(a,b) temp=a;a=b;b=temp;
#define MAXRAND     0x7fffffff
#define BIGNUM         1e37
#define WEEBIT         0.000000000001
#define MAXDIM         10
#define MAXSTR         48
#define pwrtwo(x) (1 << (x))
#define CLD_COV 10
#endif

enum SensorType {SB , AB};
enum SensorProvider {DG, PL, PT};
enum PyImageSelect {OR, BD, NX};

using std::vector;

typedef struct tagUI2DPoint
{
	uint32 m_X;
	uint32 m_Y;
    
    tagUI2DPoint()
    {
        m_X = 0;
        m_Y = 0;
    }
    tagUI2DPoint(uint32 m_X, uint32 m_Y):m_X(m_X),m_Y(m_Y)
    {
        
    }
    tagUI2DPoint(const tagUI2DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
    }
    
} UI2DPOINT;

typedef struct tagUI3DPoint
{
	uint32 m_X;
	uint32 m_Y;
	uint32 m_Z;
    
    tagUI3DPoint()
    {
        m_X = 0;
        m_Y = 0;
        m_Z = 0;
    }
    tagUI3DPoint(uint32 m_X, uint32 m_Y, uint32 m_Z):m_X(m_X),m_Y(m_Y),m_Z(m_Z)
    {
        
    }
    tagUI3DPoint(const tagUI3DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
        m_Z = p.m_Z;
    }
    
} UI3DPOINT;

typedef struct tagD2DPoint
{
    double m_X;
    double m_Y;
    
    tagD2DPoint()
    {
        m_X = 0;
        m_Y = 0;
    }
    tagD2DPoint(double m_X, double m_Y):m_X(m_X),m_Y(m_Y)
    {
        
    }
    tagD2DPoint(const tagD2DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
    }
    
    tagD2DPoint& operator=(const tagD2DPoint &p)
    {
        this->m_X = p.m_X;
        this->m_Y = p.m_Y;
        return *this;
    }
    
    tagD2DPoint operator+(const tagD2DPoint &p)
    {
        tagD2DPoint temp(this->m_X + p.m_X, this->m_Y + p.m_Y);
        return temp;
    }
    
    tagD2DPoint operator-(const tagD2DPoint &p)
    {
        tagD2DPoint temp(this->m_X - p.m_X, this->m_Y - p.m_Y);
        return temp;
    }
    
} D2DPOINT;

typedef struct tagD3DPoint
{
    double m_X;
	double m_Y;
	double m_Z;
    short m_roh;
    uint8 flag;
    
    tagD3DPoint()
    {
        m_X = 0;
        m_Y = 0;
        m_Z = 0;
        m_roh = 0;
        flag = false;
    }
    tagD3DPoint(double m_X, double m_Y, double m_Z, short roh = 0, uint8 flag = 0):m_X(m_X),m_Y(m_Y),m_Z(m_Z),m_roh(roh),flag(flag)
    {
        
    }
    tagD3DPoint(const tagD3DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
        m_Z = p.m_Z;
        m_roh = p.m_roh;
        flag = p.flag;
    }
    
    tagD3DPoint& operator=(const tagD2DPoint &p)
    {
        this->m_X = p.m_X;
        this->m_Y = p.m_Y;
        return *this;
    }
    
    tagD3DPoint& operator=(const tagD3DPoint &p)
    {
        this->m_X = p.m_X;
        this->m_Y = p.m_Y;
        this->m_Z = p.m_Z;
        this->m_roh = p.m_roh;
        this->flag = p.flag;
        return *this;
    }
    
    tagD3DPoint operator+(const tagD3DPoint &p) const
    {
        tagD3DPoint temp(this->m_X + p.m_X, this->m_Y + p.m_Y, this->m_Z + p.m_Z, this->m_roh + p.m_roh, 0);
        return temp;
    }
    
    tagD3DPoint operator-(const tagD3DPoint &p) const
    {
        tagD3DPoint temp(this->m_X - p.m_X, this->m_Y - p.m_Y, this->m_Z - p.m_Z, this->m_roh - p.m_roh, 0);
        return temp;
    }
    
} D3DPOINT;

typedef struct tagD3DPointSave
{
    float m_X;
    float m_Y;
    float m_Z;
} D3DPOINTSAVE;

typedef struct tagTransParam
{
	double t_c, m_c;
	//UTM param
	double sa, sb, e2, e2cuadrada, c;
	// a and b are radius and eccentricity of WGS84 ellipsoid, repectively.
	double a, e;
	// phi_c and lambda_0 are the latitude of true scale of standard parallel and meridian along positive Y axis, respectively.
	double phi_c, lambda_0;
	
	int pm;
	int projection;
    int utm_zone;
	
	char direction[2];
	bool bHemisphere;
} TransParam;

typedef struct tagCSize
{
	unsigned int width;
	unsigned int height;
    
    tagCSize(){}
    
    tagCSize(const int width, const int height):width(width), height(height)
    {
    }
    
    tagCSize(const tagCSize &size)
    {
        width = size.width;
        height = size.height;
    }
} CSize;

typedef struct tagNCCflag
{
	uint8 rotate_flag;
	uint8 multi_flag; 
	uint8 multi_flag_sum;
	uint8 inter_flag;
	uint8 weight_flag;
} NCCflag;

typedef struct tagNCCresult
{
    float result2; //first peak height
    float result3; //second peak height
    short minHeight;
    short maxHeight;
    
	short result0; //first peak roh
	short result1; //second peak roh
	short max_WNCC;
    
    //short GNCC;
    unsigned short NumOfHeight;
    //float *GNCC_multi;
	unsigned char result4; //peak count
    
    //int max_WNCC_pos;
	
    
    bool check_height_change;
	//int roh_count;
    
    tagNCCresult()
    {
        result2 = Nodata;
        result3 = Nodata;
        minHeight = -1000;
        maxHeight = -1000;
        result0 = -1.0;
        result1 = -1.0;
        max_WNCC = -1.0;
        NumOfHeight = 0;
        result4 = 0;
        check_height_change = false;
    }
} NCCresult;

typedef struct tagMultiMPs
{
    short peak_roh; //first peak roh
    float peak_height; //first peak height
    bool check_matched;
} MultiMPs;

typedef struct UpdateGrid{
    UpdateGrid(long len, long num_pairs) : len(len), num_pairs(num_pairs),
        _Height(len, 0), _minHeight(len, 0), _maxHeight(len, 0), _roh(len, 0),
        /* HACK for bug TODO fixme (remove +1) */ _ortho_ncc(len * (num_pairs + 1), 0),
        _Mean_ortho_ncc(len, 0), _Matched_flag(len, 0), _anchor_flag(len, 0),
        _selected_pair(len, 0), _total_images(len, 0), _ncc_seleceted_pair(len, 0)
        {}
    UpdateGrid() : UpdateGrid(0, 0) {}

    float &Height(long i) { return _Height[i]; }
    const float &Height(long i) const { return _Height[i]; }
    short &minHeight(long i) { return _minHeight[i]; }
    const short &minHeight(long i) const { return _minHeight[i]; }
    short &maxHeight(long i) { return _maxHeight[i]; }
    const short &maxHeight(long i) const { return _maxHeight[i]; }

    short &roh(long i) { return _roh[i]; }
    const short &roh(long i) const { return _roh[i]; }

    /* HACK for bug TODO fixme (remove +1s) */
    short &ortho_ncc(long i, long n) { return _ortho_ncc[i * (num_pairs + 1) + n + 1]; }
    const short &ortho_ncc(long i, long n) const { return _ortho_ncc[i]; }

    short &Mean_ortho_ncc(long i) { return _Mean_ortho_ncc[i]; }
    const short &Mean_ortho_ncc(long i) const { return _Mean_ortho_ncc[i]; }

    unsigned char &Matched_flag(long i) { return _Matched_flag[i]; }
    const unsigned char &Matched_flag(long i) const { return _Matched_flag[i]; }
    unsigned char &anchor_flag(long i) { return _anchor_flag[i]; }
    const unsigned char &anchor_flag(long i) const { return _anchor_flag[i]; }

    unsigned char &selected_pair(long i) { return _selected_pair[i]; }
    const unsigned char &selected_pair(long i) const { return _selected_pair[i]; }
    unsigned char &total_images(long i) { return _total_images[i]; }
    const unsigned char &total_images(long i) const { return _total_images[i]; }
    signed char &ncc_seleceted_pair(long i) { return _ncc_seleceted_pair[i]; }
    const signed char &ncc_seleceted_pair(long i) const { return _ncc_seleceted_pair[i]; }


private:
    long len;
    long num_pairs;
    vector<float> _Height; //after blunder detection
    vector<short> _minHeight;
    vector<short> _maxHeight;

    vector<short> _roh;
    vector<short> _ortho_ncc;
    vector<short> _Mean_ortho_ncc; // selected peak pair ortho ncc

    vector<unsigned char> _Matched_flag;
    vector<unsigned char> _anchor_flag;

    vector<unsigned char> _selected_pair; //reference image
    vector<unsigned char> _total_images;
    vector<signed char> _ncc_seleceted_pair; //selected peak pair
//    float height_counts;
}UGRID;

typedef struct tagVoxelinfo
{
    //float ANCC;
    //float height;
    bool flag_cal[MaxNCC];
    short INCC[MaxNCC];
    //float *INCC_multi;
}VOXEL;

typedef struct LSFinfo{
    //unsigned char lsf_std;
    unsigned char lsf_kernel;
}LSFINFO;

class CPairInfo
{
private:
    
    int m_NumberOfPairs;
    int m_SelectNumberOfPairs;
    int m_MinOffImage;
    float m_HeightStep;
    
    std::vector<UI2DPOINT> m_pairs;
    std::vector<float> m_BHratio;
    std::vector<float> m_ConvergenceAngle;
    std::vector<float> m_CenterDist;
    std::vector<float> m_SigmaZ;

    void allocate(int numberofpairs)
    {
        m_pairs = std::vector<UI2DPOINT>(numberofpairs);
        m_BHratio = std::vector<float>(numberofpairs);
        m_ConvergenceAngle = std::vector<float>(numberofpairs);
        m_CenterDist = std::vector<float>(numberofpairs);
        m_SigmaZ = std::vector<float>(numberofpairs);
        
        printf("allocate pairinfo array\n");
    }
    
public:
    CPairInfo() : m_NumberOfPairs(0), m_MinOffImage(-1), m_HeightStep(0), m_SelectNumberOfPairs(0)
    {
    }
    
    CPairInfo(int numberofpairs) : m_NumberOfPairs(numberofpairs)
    {
        m_SelectNumberOfPairs = m_NumberOfPairs;
        allocate(m_NumberOfPairs);
    }
    
    void SetNumberOfPairs(int numberofpairs)
    {
        m_NumberOfPairs = numberofpairs;
        m_SelectNumberOfPairs = m_NumberOfPairs;
        
        allocate(m_NumberOfPairs);
    }
    
    void SetSelectNumberOfPairs(int numberofpairs)
    {
        m_SelectNumberOfPairs = numberofpairs;
    }
    
    void SetMinOffImage(int minoffimage)
    {
        m_MinOffImage = minoffimage;
    }
    
    void SetHeightStep(int heightStep)
    {
        m_HeightStep = heightStep;
    }
    
    void SetPairs(int pos, UI2DPOINT value)
    {
        m_pairs[pos] = value;
    }
    
    void SetPairs(int pos, int X, int Y)
    {
        m_pairs[pos].m_X = X;
        m_pairs[pos].m_Y = Y;
    }
    
    void SetBHratio(int pos, float value)
    {
        m_BHratio[pos] = value;
    }
    
    void SetConvergenceAngle(int pos, float value)
    {
        m_ConvergenceAngle[pos] = value;
    }
    
    void SetCenterDist(int pos, float value)
    {
        m_CenterDist[pos] = value;
    }
    
    void SetSigmaZ(int pos, float value)
    {
        m_SigmaZ[pos] = value;
    }
    
    int& NumberOfPairs()
    {
        return m_NumberOfPairs;
    }
    
    int& SelectNumberOfPairs()
    {
        return m_SelectNumberOfPairs;
    }
    
    int& MinOffImageID()
    {
        return m_MinOffImage;
    }
    
    float& HeightStep()
    {
        return m_HeightStep;
    }
    
    UI2DPOINT& pairs(int pos)
    {
        return m_pairs[pos];
    }
    
    float& BHratio(int pos)
    {
        return m_BHratio[pos];
    }
    
    float& ConvergenceAngle(int pos)
    {
        return m_ConvergenceAngle[pos];
    }
    
    float& CenterDist(int pos)
    {
        return m_CenterDist[pos];
    }
    
    float& SigmaZ(int pos)
    {
        return m_SigmaZ[pos];
    }
};

typedef struct BlunderIP{
	CSize Size_Grid2D;
	double gridspace;
	double Hinterval; 
	const double* Boundary;
	uint8 Pyramid_step;
	uint8 iteration;
	bool height_check_flag;
}BL;

//Aerial frame Camera
typedef struct tagCameraInfo
{
    double m_focalLength;
    CSize m_ImageSize;
    double m_CCDSize;
    double m_ppx;
    double m_ppy;
} CAMERA_INFO;

typedef struct tagRotationMatrix
{
    double m11, m12, m13, m21, m22, m23, m31, m32, m33;
} RM;

typedef struct tagEO
{
    int strip_ID;
    int photo_ID;
    char path[500];
    double m_Xl;
    double m_Yl;
    double m_Zl;
    double m_Wl;
    double m_Pl;
    double m_Kl;
    RM m_Rm;
} EO;


typedef struct tagFrameInfo
{
    CAMERA_INFO m_Camera;
    int NumberofStip;
    int NumberofPhotos;
    int start_stripID;
    int end_stripID;
    EO* Photoinfo;
} FrameInfo;


typedef struct ProjectInfo{
	double resolution;
	double DEM_resolution;
	double preDEM_space;
	double cal_boundary[4];
	double RA_param[MaxNCC][2];
	double seedDEMsigma;
    double LBoundary[4];
    double RBoundary[4];
    double GCP_spacing;
    double CA_th;
    double Cloud_th;
    
	double minHeight;
	double maxHeight;
	double System_memory;
    double required_memory;
    
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	int threads_num;
    int number_of_images;
    enum SensorType sensor_type; // 1 is for RFM (default), 2 is for Collinear Equation (Frame)
    enum SensorProvider sensor_provider; //DG = DG, Pleiades = PL if sensor_type = 1
    uint8 pyramid_level;
    uint8 end_level;
    uint8 SDM_SS;
    double SDM_AS;
    double SDM_days;
    uint8 image_bits;
    
	char Imagefilename[MaxImages][500];
	char RPCfilename[MaxImages][500];
    char Imagemetafile[MaxImages][500];
    
	char save_filepath[500];
	char Outputpath_name[500];
	char tmpdir[500];
	char tile_info[500];
	char priori_DEM_tif[500];
	char metafilename[500];
	
	bool check_minH;
	bool check_maxH;
	bool check_tiles_SR;
	bool check_tiles_ER;
	bool check_tiles_SC;
	bool check_tiles_EC;
	bool check_gridonly;
	bool check_boundary;
	bool check_checktiff;
    bool check_ortho;
	bool IsRA, IsSP, IsRR, IsSaveStep, Overall_DEM, Affine_RA, pre_DEMtif, check_tile_array;
    bool check_Matchtag;
    bool check_selected_image[MaxImages];
    bool check_full_cal;
    
    //SGM test flag
    bool check_SNCC;
    bool check_updateheight;
    bool check_blunderdetection;
    bool check_NCCpeak;
    bool check_minTH;
    bool check_adaptive_P2;
    bool check_orthoblunder;
    bool check_8d;
    int SGM_py;
    
    FrameInfo frameinfo;
    
	uint8 SPnumber[2],NumOfTile_row, NumOfTile_col;	
} ProInfo;

typedef struct ArgumentInfo{
	double DEM_space;
	double seedDEMsigma;
	double minHeight;
	double maxHeight;
	double ra_line[MaxNCC];
	double ra_sample[MaxNCC];
	double Min_X, Max_X, Min_Y, Max_Y;
	double image_resolution;
	double overlap_length;
    double focal_length;
    double CCD_size;
    double System_memory;
    double SDM_days;
    double SDM_AS;
    double DS_sigma;
    double DS_tx;
    double DS_ty;
    double GCP_spacing;
    double CA_th;
    double Cloud_th;
    
	int check_arg; // 0 : no input, 1: 3 input
	int Threads_num;
	int start_row;
	int end_row;
	int start_col;
	int end_col;
	int RA_row;
	int RA_col;
	int tilesize;
	int projection; //PS = 1, UTM = 2
    TransParam param;
    int utm_zone;
    enum SensorType sensor_type; // SB is for RFM (default), AB is for Collinear Equation (Frame)
    enum SensorProvider sensor_provider; //DG = DG, Pleiades = PL if sensor_type = 1
    int ortho_count;
    int RA_only;
    int number_of_images; // 2 is for stereo (default), n is for multi more than 3
    uint8 pyramid_level;
    uint8 SDM_SS;
    int DS_kernel;
    int RA_line_count;
    int RA_sample_count;
    
    char Image[MaxImages][500];
	char Outputpath[500];
	char Outputpath_name[500];
	char seedDEMfilename[500];
	char metafilename[500];
    char EO_Path[500];
    char DEM_input_file[500];
    char Multi_input_file[500];
    
    bool check_DS_txy;
    bool check_downsample;
	bool check_DEM_space;
	bool check_Threads_num;
	bool check_seeddem;
	bool check_minH;
	bool check_maxH;
	bool check_tiles_SR;
	bool check_tiles_ER;
	bool check_tiles_SC;
	bool check_tiles_EC;
	bool check_RA_line;
	bool check_RA_sample;
	bool check_gridonly; 
	bool check_RA_tileR;
	bool check_RA_tileC;
	bool check_tilesize;
	bool check_boundary;
	bool check_checktiff;
	bool check_RA_only;
    bool check_ortho;
    bool check_imageresolution;
    bool check_LSF;
    bool check_LSF_DEM;
    bool check_LSFDEMpath;
    int check_LSF2;
    bool check_Matchtag;
    bool check_EO;
    bool check_fl;
    bool check_ccd;
    bool check_full_cal;
    int check_txt_input;
    int check_coreg;
    int check_sdm_ortho;
    int check_DEM_coreg_output;
    
    //SGM test flag
    bool check_SNCC;
    bool check_updateheight;
    bool check_blunderdetection;
    bool check_NCCpeak;
    bool check_minTH;
    bool check_adaptive_P2;
    bool check_orthoblunder;
    bool check_8d;
    int SGM_py;
    
} ARGINFO;

typedef struct tagImageGSD
{
    float row_GSD;
    float col_GSD;
    float pro_GSD;
} ImageGSD;

typedef struct tagBandInfo
{
    float abscalfactor;
    float effbw;
    float tdi;
} BandInfo;

typedef struct tagImageInfo
{
    float Mean_sun_azimuth_angle;
    float Mean_sun_elevation;
    float Mean_sat_azimuth_angle;
    float Mean_sat_elevation;
    float Intrack_angle;
    float Crosstrack_angle;
    float Offnadir_angle;
    float Offnadir_angle_xml;
    float Mean_sat_azimuth_angle_xml;
    float cloud;
    float Image_ori;
    float Image_ori_azi;
    float dx,dy,f;
    float UL[3], UR[3],LR[3],LL[3];
    double Center[2];
    float convergence_angle;
	
    ImageGSD GSD;
    
	int month;
    int date;
    int year;
    int scandirection;
    
    char filename[500];
    char Ofilename[500];
    char imagetime[500];
    char SatID[500];
} ImageInfo;

typedef struct
{
    long int ncols;
    long int nrows;
    long double **val;
    long double *data;
} GMA_double;


typedef struct UpdateGridSDM{
    float roh;
    float ortho_ncc;
    float col_shift;
    float row_shift;
}UGRIDSDM;

typedef struct tagNCCresultSDM
{
    float result0;
    D2DPOINT result2;
    D2DPOINT result3;
} NCCresultSDM;

typedef struct tagTINinfo
{
    D3DPOINT *normal;
    uint16 *slope;
    uint16 *aspect;
    float *ncc;
    float *dem;
} TINinfo;

typedef struct tagConformalparam
{
    float scale;
    float omega;
    float phi;
    float kappa;
    float Tx;
    float Ty;
    float Tz;
} Conformalparam;

typedef struct taglevelinfo
{
    const uint16 * const *py_Images;
    const uint8 * const *py_OriImages;
    const uint16 * const *py_MagImages;
    
    const uint16 * const *py_BImages;
    const uint16 * const *py_BMagImages;

    const uint16 * const *py_Images_next;
    const uint8 * const *py_OriImages_next;
    const uint16 * const *py_MagImages_next;
    
    const vector<D2DPOINT> *py_Startpos;
    const vector<D2DPOINT> *py_BStartpos;
    const vector<D2DPOINT> *py_Startpos_next;
    
    const ImageInfo *imageinfo;
    
    const double * const * const *RPCs;
    const double *Boundary;
    const int *Pyramid_step;
    const CSize * const *py_Sizes;
    const unsigned char *Template_size;
    const double * const *ImageAdjust;
    const TransParam *param;
    const unsigned char *NumOfIAparam;
    const double *bin_angle;
    const int *blunder_selected_level; 
    const CSize *Size_Grid2D;
    const long int *Grid_length;
    const D2DPOINT* GridPts;
    const D2DPOINT* Grid_wgs;
    int reference_id;
    const double *height_step;
    const double *grid_resolution;
    const double *Hinterval;
    //double *minmaxHeight; //iteratively changed
    const int *Py_combined_level;
    const unsigned char *iteration;
    bool check_ortho;
    bool *check_matching_rate;
    CPairInfo *pairinfo;
    const CSize *Imagesize_ori;
    bool check_SGM;
    double MPP;
} LevelInfo;

template <typename T>
class Matrix {
public:
    Matrix(size_t rows, size_t cols)
      : rows_ {rows},
        cols_ {cols},
        data_(rows * cols)
    {
    }

    Matrix(size_t rows, size_t cols, const T& initial_value)
      : rows_ {rows},
        cols_ {cols},
        data_(rows * cols, initial_value)
    {
    }

    T& operator() (unsigned row, unsigned col)
    {
        return data_[cols_*row + col];
    }

    const T& operator() (unsigned row, unsigned col) const
    {
        return data_[cols_*row + col];
    }

    inline
    const T * row(size_t i) const
    {
        return data_.data()  + cols_*i;
    }

private:
    size_t rows_ = 0;
    size_t cols_ = 0;
    vector<T> data_;
};

typedef struct tagSetKernel
{
    int reference_id;
    int ti;
    const int Half_template_size;
    unsigned patch_size;

    Matrix<double> left_patch_vecs;
    Matrix<double> left_mag_patch_vecs;
    Matrix<double> right_patch_vecs;
    Matrix<double> right_mag_patch_vecs;


    tagSetKernel(const int reference_id,const int ti,const int Half_template_size):
        reference_id(reference_id),
        ti(ti),
        Half_template_size(Half_template_size),
        patch_size((2*Half_template_size+1) * (2*Half_template_size+1)),
        left_patch_vecs{3, patch_size},
        left_mag_patch_vecs{3, patch_size},
        right_patch_vecs{3, patch_size},
        right_mag_patch_vecs{3, patch_size}
    {
    }

    ~tagSetKernel() {
    }
} SetKernel;

template <typename T>
class MixedHeight3DGrid{
public:
    MixedHeight3DGrid() = default; // needed for use in containers
    MixedHeight3DGrid(long length) : _indices(length, 0) {}

    void reserve(long index, long length) {
        _indices[index] = length;
    }

    void allocate(const T &iv) {
        // get the total length needed
        long len = 0;
        for(auto i : _indices)
        {
            len += i;
        }

        // allocate the vector
        _data = vector<T>(len, iv);

        // update the indices to hold offsets
        long curr_offset = 0;
        for(auto &i : _indices)
        {
            auto curr_len = i;
            i = curr_offset;
            curr_offset += curr_len;
        }
    }

    void allocate() {
        T iv{};
        allocate(iv);
    }

    T& value(long i, long j) {
        return _data[_indices[i] + j];
    }

    const T& value(long i, long j) const {
        return _data[_indices[i] + j];
    }

private:

    vector<long> _indices;
    vector<T> _data;
};

typedef MixedHeight3DGrid<short> SumCostContainer;

class GridPairs {
public:
    GridPairs(long length) : pair_ids(length, -1), next_id(0) {}

    void add_pairs(long index, const vector<short> &pairs) {
        short id = find_pair_id(pairs);
        if(id != -1) {
            pair_ids[index] = id;
        } else {
            // not in our list, so add it
            pair_lists[next_id] = pairs;
            pair_ids[index] = next_id;
            next_id++;
       }
    }

    void remap_pairs(const std::map<short, short> &pair_map) {
        for(auto &it : pair_lists) {
            for(short &p_num : it.second) {
                if(pair_map.find(p_num) != pair_map.end()) {
                    p_num = pair_map.at(p_num);
                }
            }
        }
    }

    const vector<short>& get_pairs(long index) const {
        int id = pair_ids[index];
        if(id < 0) {
            return empty;
        } else {
            return pair_lists.at(id);
        }
    }

private:

    short find_pair_id(const vector<short> &pairs) {
        for(auto &it : pair_lists) {
            if(it.second == pairs) {
                return it.first;
            }
        }
        return -1;
    }

    std::map<int, vector<short>> pair_lists;
    vector<int> pair_ids;
    vector<short> empty;
    int next_id;
};

#endif

