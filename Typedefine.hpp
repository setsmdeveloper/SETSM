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
#include <algorithm>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "log.hpp"
#include <stdexcept>

#define PI 3.141592653589793
#define DegToRad PI/180
#define RadToDeg 180/PI
#define UMToMM 0.001
#define MMToUM 1000
#define MaxImages 1000
#define MaxNCC 20000
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
enum SensorProvider {DG, PL, PT, BS}; //Digital Globe, Pleiades, Planet, BlackSky
enum PyImageSelect {OR, BD, NX};

using std::vector;

typedef struct tagCLB
{
    double b11,b12,b13,b14,b15,b16;
    double b21,b22,b23,b24,b25,b26;
    double J, K;
    /*void printf()
    {
        printf("CLB\n%f\t%f\t%f\t%f\t%f\t%f\n",b11,b12,b13,b14,b15,b16);
        printf("%f\t%f\t%f\t%f\t%f\t%f\n",b21,b22,b23,b24,b25,b26);
        printf("%f\t%f\n",J, K);
    }*/
} CLB;

typedef struct tagUParams
{
    double f;//focal length
    double xp, yp;//principal point
    double k1, k2, k3;//ratdial distortion
    double p1, p2;//tangential distortion
    double a1, a2;//affinity
    
    double X, Y, Z;//perspective center
    double omega, phi, kappa;//attitude
    
    double XA, YA, ZA;//Ground Coordinate
    
    double xa, ya;//image point
} UPARAMS;

typedef struct tagPDCs
{
    double a[2][2];//image point term
    double b[2][6];//EO parameters and Ground Coordinates term
    double c[2][10];//Calibration data term
    double Fo, Go;//initial approximation of F, G
} PDCS;

typedef struct tagDLTparam
{
    double A1,A2,A3,B1,B2,B3,C1,C2,C3,D1,D2;
} DLTPARAM;

typedef struct tagUI8_2DPoint
{
    uint8 m_X;
    uint8 m_Y;
    
    tagUI8_2DPoint()
    {
        m_X = 0;
        m_Y = 0;
    }
    tagUI8_2DPoint(uint8 m_X, uint8 m_Y):m_X(m_X),m_Y(m_Y)
    {
        
    }
    tagUI8_2DPoint(const tagUI8_2DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
    }
    
} UI8_2DPoint;

typedef struct tagUI16_2DPoint
{
    uint16 m_X;
    uint16 m_Y;
    
    tagUI16_2DPoint()
    {
        m_X = 0;
        m_Y = 0;
    }
    tagUI16_2DPoint(uint16 m_X, uint16 m_Y):m_X(m_X),m_Y(m_Y)
    {
        
    }
    tagUI16_2DPoint(const tagUI16_2DPoint &p)
    {
        m_X = p.m_X;
        m_Y = p.m_Y;
    }
    
} UI16_2DPoint;

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

typedef struct tagF3DPoint
{
    float m_X;
    float m_Y;
    float m_Z;
} F3DPOINT;

typedef struct tagF2DPoint
{
    float m_X;
    float m_Y;
} F2DPOINT;

typedef struct tagD3DPointSave
{
    float m_X;
    float m_Y;
    float m_Z;
} D3DPOINTSAVE;

typedef struct tagDate
{
    int year;
    int month;
    int day;
    int period;
} DATE;

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
	int bHemisphere;
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
    short ortho_roh;
    //float peak_height; //first peak height
    bool check_matched;
} MultiMPs;

typedef struct UpdateGrid{
    UpdateGrid(long len, long num_pairs) : len(len), num_pairs(num_pairs),
        _angle(len,0), _Height(len, 0), _minHeight(len, 0), _maxHeight(len, 0), _roh(len, 0),
        /* HACK for bug TODO fixme (remove +1) */ _ortho_ncc(len * (num_pairs + 1), 0),
        _Mean_ortho_ncc(len, 0), _Matched_flag(len, 0), _anchor_flag(len, 0),
        _selected_pair(len, 0), _total_images(len, 0), _ncc_seleceted_pair(len, 0), PairCheck(len), _ndvi_pixel(len,0)//, PairArray(len)
        {}
    UpdateGrid() : UpdateGrid(0, 0) {}

    unsigned char &NDVI_pixel(long i) {return _ndvi_pixel[i]; };
    const unsigned char &NDVI_pixel(long i) const {return _ndvi_pixel[i]; };
    
    float &Angle(long i) { return _angle[i]; }
    const float &Angle(long i) const { return _angle[i]; }
    
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

    short &selected_pair(long i) { return _selected_pair[i]; }
    const short &selected_pair(long i) const { return _selected_pair[i]; }
    unsigned char &total_images(long i) { return _total_images[i]; }
    const unsigned char &total_images(long i) const { return _total_images[i]; }
    short &ncc_seleceted_pair(long i) { return _ncc_seleceted_pair[i]; }
    const short &ncc_seleceted_pair(long i) const { return _ncc_seleceted_pair[i]; }
    
    //vector<vector<short>> PairArray;
    vector<vector<short>> PairCheck;

private:
    long len;
    long num_pairs;
    vector<float> _angle;
    vector<float> _Height; //after blunder detection
    vector<short> _minHeight;
    vector<short> _maxHeight;

    vector<short> _roh;
    vector<short> _ortho_ncc;
    vector<short> _Mean_ortho_ncc; // selected peak pair ortho ncc

    vector<unsigned char> _Matched_flag;
    vector<unsigned char> _anchor_flag;

    vector<short> _selected_pair; //reference image
    vector<unsigned char> _total_images;
    vector<short> _ncc_seleceted_pair; //selected peak pair
    
    vector<unsigned char> _ndvi_pixel;
    
//    float height_counts;
}UGRID;
/*
typedef struct tagVoxelinfo
{
    //float ANCC;
    //float height;
    bool flag_cal[MaxNCC];
    short INCC[MaxNCC];
    //float *INCC_multi;
}VOXEL;
*/
typedef struct LSFinfo{
    //unsigned char lsf_std;
    unsigned char lsf_kernel;
}LSFINFO;

class GridPairs {
public:
    vector<float> grid_sigmaZ; //NR noise variance
    vector<float> grid_max_sigmaZ; // height step for Verticallinelocus and AWNCC_MPs,and MPP
    vector<unsigned short> grid_height_step;
    vector<float> grid_mean_sigmaZ; //NR noise variance weight
    
    GridPairs(long length) : pair_ids(length, -1), next_id(0), grid_sigmaZ(length,0), grid_max_sigmaZ(length, 0), grid_mean_sigmaZ(length, 0), grid_height_step(length, 0) {}

    GridPairs()
    {
    }
    
    void Setlength(long length)
    {
        for(size_t i = 0 ; i < length ; i++)
            pair_ids.push_back(-1);
        next_id = 0;
        
        //grid_sigmaZ = std::vector<float>(length);
        //grid_max_sigmaZ = std::vector<float>(length);
        //grid_height_step = std::vector<unsigned short>(length);
        //grid_mean_sigmaZ = std::vector<float>(length);
    }
    
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

    void remove_pairs(long index, std::vector<short> pairs_to_remove) {
        std::vector<short> current_pairs = get_pairs(index);
        std::vector<short> new_pairs;

        // build list of new pairs for index
        for(auto pnum : current_pairs) {
            if(std::find(pairs_to_remove.begin(), pairs_to_remove.end(), pnum) == pairs_to_remove.end()) {
                // pair is not in removal list, so add to new list
                new_pairs.push_back(pnum);
            }
        }
        add_pairs(index, new_pairs);
    }

    void clear_pairs(long length)
    {
        
        pair_ids.clear();
        empty.clear();
        
        grid_sigmaZ.clear();
        grid_max_sigmaZ.clear();
        grid_height_step.clear();
        grid_mean_sigmaZ.clear();
        
        for(size_t i = 0 ; i < length ; i++)
        {
            pair_lists[i].clear();
            //pair_ids.push_back(-1);
        }
        next_id = 0;
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

typedef struct tagGNCCGrid
{
    // Constructor
    tagGNCCGrid()
            : gncc_grid(0), Gridpairs(0) {}

    tagGNCCGrid(size_t length, const GridPairs &Gridpairs)
            : gncc_grid(length), Gridpairs(Gridpairs)
    {
        for(size_t index = 0 ; index < length ; index++)
        {
            allocate(index);
        }
    }
    
    // Const and nonconst accessors for GNCC
    short& GNCC(size_t grid_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        //size_t num_pairs = pairs.size();
        return gncc_grid[grid_index][pair_index];
    }

    // const versions of the above
    const short& GNCC(size_t grid_index, size_t pair_id) const {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        //size_t num_pairs = pairs.size();
        return gncc_grid[grid_index][pair_index];
    }
    
    
    bool has_pair(size_t grid_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        return pair_id_to_index(pairs, pair_id) >= 0;
    }
    
    /** Set size to zero and clear memory */
    void clear(size_t grid_index) {
        std::vector<short>().swap(gncc_grid[grid_index]);
    }
    
    void allocate(size_t grid_index) {
        auto num_pairs = Gridpairs.get_pairs(grid_index).size();
        gncc_grid[grid_index] = std::vector<short>(num_pairs, INCC_UNSET);
    }
    
private:
    std::vector<std::vector<short>> gncc_grid;
    GridPairs Gridpairs;
    
    static constexpr short INCC_UNSET = -1000;

    int pair_id_to_index(const std::vector<short> &pairs, size_t pair_id) const {
        auto it = std::find(pairs.begin(), pairs.end(), pair_id);
        if(it == pairs.end()) {
            return -1;
        }
        return std::distance(pairs.begin(), it);
    }

    void dbg_crash(size_t grid_index, const std::vector<short> &pairs, size_t pair_id) const {
        LOG("GNCC Failed to find pair id %lu for grid index %lu\n", pair_id, grid_index);
        LOG("number of pairs: %d\n", pairs.size());

        LOG("pairs:\n");
        for(int i = 0; i < pairs.size(); i++) {
            LOG("    pair[%d] = %d\n", i, pairs[i]);
        }
        abort();
    }
}GNCCGrid;

typedef struct tagGridroh
{
    // Constructor
    tagGridroh()
            : grid_rohth(0), Gridpairs(0), length(0) {}

    /*
    tagGridroh(size_t length, const GridPairs &Gridpairs)
            : grid_rohth(length), Gridpairs(Gridpairs)
    {
        for(size_t index = 0 ; index < length ; index++)
        {
            allocate(index);
        }
    }
    */
    // Const and nonconst accessors for GNCC
    short& Grid_roh(size_t grid_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        //size_t num_pairs = pairs.size();
        return grid_rohth[grid_index][pair_index];
    }

    // const versions of the above
    const short& Grid_roh(size_t grid_index, size_t pair_id) const {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        //size_t num_pairs = pairs.size();
        return grid_rohth[grid_index][pair_index];
    }
    
    
    bool has_pair(size_t grid_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        return pair_id_to_index(pairs, pair_id) >= 0;
    }
    
    /** Set size to zero and clear memory */
    void clear(size_t grid_index) {
        std::vector<short>().swap(grid_rohth[grid_index]);
    }
    
    void clear_all() {
        if(length > 0)
        {
            for(int i = 0 ; i < length ; i++)
                std::vector<short>().swap(grid_rohth[i]);
            grid_rohth.clear();
            
            Gridpairs.clear_pairs(length);
        }
        
        length = 0;
    }
    
    
    void allocate(long total_grid, short val, GridPairs m_Gridpairs) {
        if(length > 0)
        {
            for(int i = 0 ; i < length ; i++)
                std::vector<short>().swap(grid_rohth[i]);
            grid_rohth.clear();
            
            Gridpairs.clear_pairs(length);
        }
        
        length = total_grid;
        
        grid_rohth = std::vector<std::vector<short>>(length);
        Gridpairs.Setlength(length);
        
        for(long i = 0 ; i < length ; i++)
        {
            auto &pairs = m_Gridpairs.get_pairs(i);
            auto num_pairs = pairs.size();
            if(pairs.size() > 0)
            {
                //printf("pairsize %d\n",pairs.size());
                Gridpairs.add_pairs(i,pairs);
                grid_rohth[i] = std::vector<short>(num_pairs, val);
            }
        }
    }
    
private:
    std::vector<std::vector<short>> grid_rohth;
    GridPairs Gridpairs;
    long length;
    static constexpr short INCC_UNSET = -1000;

    int pair_id_to_index(const std::vector<short> &pairs, size_t pair_id) const {
        auto it = std::find(pairs.begin(), pairs.end(), pair_id);
        if(it == pairs.end()) {
            return -1;
        }
        return std::distance(pairs.begin(), it);
    }

    void dbg_crash(size_t grid_index, const std::vector<short> &pairs, size_t pair_id) const {
        LOG("gridroh Failed to find pair id %lu for grid index %lu\n", pair_id, grid_index);
        LOG("number of pairs: %d\n", pairs.size());

        LOG("pairs:\n");
        for(int i = 0; i < pairs.size(); i++) {
            LOG("    pair[%d] = %d\n", i, pairs[i]);
        }
        abort();
    }
}Gridroh;

class CPairInfo
{
private:
    
    int m_NumberOfPairs;
    int m_SelectNumberOfPairs;
    int m_MinOffImage;
    float m_HeightStep;
    int m_max_countMPs_pair;
    int m_ref_imageID;
    
    std::vector<short> m_cal;
    std::vector<D2DPOINT> m_RBias;
    std::vector<unsigned char> m_cloud;
    std::vector<float> m_Tz;
    std::vector<float> m_Tz_var;
    std::vector<float> m_Tz_res;
    std::vector<vector<unsigned char>> m_linked_pair;
    
    std::vector<D2DPOINT> m_referencepairs;
    std::vector<uint16> m_oripairnumber;
    
    std::vector<UI2DPOINT> m_pairs;
    std::vector<float> m_BHratio;
    std::vector<float> m_ConvergenceAngle;
    std::vector<float> m_ConvergenceAngle_IRPC;
    std::vector<float> m_ConvergenceAngle_given;
    std::vector<float> m_CenterDist;
    std::vector<float> m_SigmaZ;
    std::vector<float> m_AzimuthDiff;
    std::vector<unsigned char> m_diff_day;
    std::vector<unsigned char> m_diff_time;
    std::vector<float> m_MPs_P;
    std::vector<float> m_resolution;
    
    std::vector<float> m_AE;
    std::vector<float> m_BIE;
    std::vector<D3DPOINT> m_BaseRay;
    
    std::vector<float> m_AE_IRPC;
    std::vector<float> m_BIE_IRPC;
    std::vector<D3DPOINT> m_BaseRay_IRPC;
    
    std::vector<float> m_ConvergenceAngle_EQ;
    std::vector<float> m_ConvergenceAngle_EQ_IRPC;
    
    std::vector<unsigned char> m_Coverage;
    std::vector<unsigned short> m_pairHeightStep;
    
    void allocate(int numberofpairs)
    {
        m_max_countMPs_pair = -1;
        m_ref_imageID = -1;
        
        m_cal = std::vector<short>(numberofpairs,1);
        m_RBias = std::vector<D2DPOINT>(numberofpairs);
        m_cloud = std::vector<unsigned char>(numberofpairs);
        m_Tz = std::vector<float>(numberofpairs,0);
        m_Tz_var = std::vector<float>(numberofpairs,1.0);
        m_Tz_res = std::vector<float>(numberofpairs,1.0);
        D2DPOINT temp(0,0);
        m_referencepairs = std::vector<D2DPOINT>(numberofpairs,temp);
        m_oripairnumber = std::vector<uint16>(numberofpairs,0);
        
        m_linked_pair = std::vector<vector<unsigned char>>(numberofpairs);
        for(int index = 0 ; index < numberofpairs ; index++)
            m_linked_pair[index] = std::vector<unsigned char>(numberofpairs,0);
        
        m_pairs = std::vector<UI2DPOINT>(numberofpairs);
        m_BHratio = std::vector<float>(numberofpairs);
        m_ConvergenceAngle = std::vector<float>(numberofpairs);
        m_ConvergenceAngle_given = std::vector<float>(numberofpairs);
        m_CenterDist = std::vector<float>(numberofpairs);
        m_SigmaZ = std::vector<float>(numberofpairs);
        m_AzimuthDiff = std::vector<float>(numberofpairs);
        m_diff_day = std::vector<unsigned char>(numberofpairs);
        m_diff_time = std::vector<unsigned char>(numberofpairs);
        m_MPs_P = std::vector<float>(numberofpairs);
        m_resolution = std::vector<float>(numberofpairs);
        
        m_AE = std::vector<float>(numberofpairs);
        m_BIE = std::vector<float>(numberofpairs);
        m_BaseRay = std::vector<D3DPOINT>(numberofpairs);
        
        m_ConvergenceAngle_IRPC = std::vector<float>(numberofpairs);
        m_AE_IRPC = std::vector<float>(numberofpairs);
        m_BIE_IRPC = std::vector<float>(numberofpairs);
        m_BaseRay_IRPC = std::vector<D3DPOINT>(numberofpairs);
        
        m_ConvergenceAngle_EQ = std::vector<float>(numberofpairs);
        m_ConvergenceAngle_EQ_IRPC = std::vector<float>(numberofpairs);
        
        m_Coverage = std::vector<unsigned char>(numberofpairs);
        m_pairHeightStep = std::vector<unsigned short>(numberofpairs);
        
        printf("allocate pairinfo array\n");
    }
    
    
    
public:
    CPairInfo() : m_NumberOfPairs(0), m_MinOffImage(-1), m_HeightStep(0), m_SelectNumberOfPairs(0), m_max_countMPs_pair(-1), m_ref_imageID(-1)
    {
    }
    
    CPairInfo(int numberofpairs) : m_NumberOfPairs(numberofpairs)
    {
        m_SelectNumberOfPairs = m_NumberOfPairs;
        allocate(m_NumberOfPairs);
    }
    
    void initialize()
    {
        m_NumberOfPairs = 0;
        m_SelectNumberOfPairs = 0;
        m_MinOffImage = -1;
        m_HeightStep = 0;
        m_max_countMPs_pair = -1;
        m_ref_imageID = -1;
        
        m_cal.clear();
        vector<short>().swap(m_cal);
        
        m_RBias.clear();
        vector<D2DPOINT>().swap(m_RBias);
        
        m_cloud.clear();
        vector<unsigned char>().swap(m_cloud);
        
        m_Tz.clear();
        vector<float>().swap(m_Tz);
        
        m_Tz_var.clear();
        vector<float>().swap(m_Tz_var);
        
        m_Tz_res.clear();
        vector<float>().swap(m_Tz_res);
        
        m_referencepairs.clear();
        vector<D2DPOINT>().swap(m_referencepairs);
        
        m_oripairnumber.clear();
        vector<uint16>().swap(m_oripairnumber);
        
        m_pairs.clear();
        vector<UI2DPOINT>().swap(m_pairs);
        
        m_BHratio.clear();
        vector<float>().swap(m_BHratio);
        
        m_ConvergenceAngle.clear();
        vector<float>().swap(m_ConvergenceAngle);
        
        m_ConvergenceAngle_given.clear();
        vector<float>().swap(m_ConvergenceAngle_given);
        
        m_CenterDist.clear();
        vector<float>().swap(m_CenterDist);
        
        m_SigmaZ.clear();
        vector<float>().swap(m_SigmaZ);
        
        m_AzimuthDiff.clear();
        vector<float>().swap(m_AzimuthDiff);
        
        m_diff_day.clear();
        vector<unsigned char>().swap(m_diff_day);
        
        m_diff_time.clear();
        vector<unsigned char>().swap(m_diff_time);
        
        m_MPs_P.clear();
        vector<float>().swap(m_MPs_P);
        
        m_resolution.clear();
        vector<float>().swap(m_resolution);
        
        m_AE.clear();
        vector<float>().swap(m_AE);
        
        m_BIE.clear();
        vector<float>().swap(m_BIE);
        
        m_BaseRay.clear();
        vector<D3DPOINT>().swap(m_BaseRay);
        
        m_ConvergenceAngle_IRPC.clear();
        vector<float>().swap(m_ConvergenceAngle_IRPC);
        
        m_AE_IRPC.clear();
        vector<float>().swap(m_AE_IRPC);
        
        m_BIE_IRPC.clear();
        vector<float>().swap(m_BIE_IRPC);
        
        m_BaseRay_IRPC.clear();
        vector<D3DPOINT>().swap(m_BaseRay_IRPC);
        
        m_ConvergenceAngle_EQ.clear();
        vector<float>().swap(m_ConvergenceAngle_EQ);
        
        m_ConvergenceAngle_EQ_IRPC.clear();
        vector<float>().swap(m_ConvergenceAngle_EQ_IRPC);
        
        m_Coverage.clear();
        vector<unsigned char>().swap(m_Coverage);
        
        m_pairHeightStep.clear();
        vector<unsigned short>().swap(m_pairHeightStep);
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
    
    void SetMaxCountMPs_pair(int pairID)
    {
        m_max_countMPs_pair = pairID;
    }
    
    void SetRefImageID(int pairID)
    {
        m_ref_imageID = pairID;
    }
    
    void SetCal(int pos, short value)
    {
        m_cal[pos] = value;
    }
    
    void SetRBias(int pos, D2DPOINT value)
    {
        m_RBias[pos] = value;
    }
    
    void SetCloud(int pos, unsigned char value)
    {
        m_cloud[pos] = value;
    }
    
    void SetTz(int pos, float value)
    {
        m_Tz[pos] = value;
    }
    
    void SetTz_var(int pos, float value)
    {
        m_Tz_var[pos] = value;
    }
    
    void SetTz_res(int pos, float value)
    {
        m_Tz_res[pos] = value;
    }
    
    void Setreferencepairs(int pos, D2DPOINT value)
    {
        m_referencepairs[pos] = value;
    }
    
    void Setoripairnumber(int pos, uint16 value)
    {
        m_oripairnumber[pos] = value;
    }
    
    void SetLinked_pair(int pos, short linked_pair, unsigned char value)
    {
        m_linked_pair[pos][linked_pair] = value;
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
    
    void SetAzimuth(int pos, float value)
    {
        m_AzimuthDiff[pos] = value;
    }
    
    void SetDiffDay(int pos, unsigned char value)
    {
        m_diff_day[pos] = value;
    }
    
    void SetDiffTime(int pos, unsigned char value)
    {
        m_diff_time[pos] = value;
    }
    
    void SetMatchingP(int pos, float value)
    {
        m_MPs_P[pos] = value;
    }
    
    void SetResolution(int pos, float value)
    {
        float convert = (value + 0.05)*10;
        int integer = (int)convert;
        convert = integer/10.0;
        
        m_resolution[pos] = convert;
    }
    
    void SetAE(int pos, float value)
    {
        m_AE[pos] = value;
    }
    
    void SetBIE(int pos, float value)
    {
        m_BIE[pos] = value;
    }
    
    void SetBaseRay(int pos, D3DPOINT value)
    {
        m_BaseRay[pos] = value;
    }
    
    void SetConvergenceAngle_IRPC(int pos, float value)
    {
        m_ConvergenceAngle_IRPC[pos] = value;
    }
    
    void SetAE_IRPC(int pos, float value)
    {
        m_AE_IRPC[pos] = value;
    }
    
    void SetBIE_IRPC(int pos, float value)
    {
        m_BIE_IRPC[pos] = value;
    }
    
    void SetBaseRay_IRPC(int pos, D3DPOINT value)
    {
        m_BaseRay_IRPC[pos] = value;
    }
    
    void SetConvergenceAngle_given(int pos, float value)
    {
        m_ConvergenceAngle_given[pos] = value;
    }
    
    void SetConvergenceAngle_EQ(int pos, float value)
    {
        m_ConvergenceAngle_EQ[pos] = value;
    }
    
    void SetConvergenceAngle_EQ_IRPC(int pos, float value)
    {
        m_ConvergenceAngle_EQ_IRPC[pos] = value;
    }
    
    void SetCoverage(int pos, unsigned char value)
    {
        m_Coverage[pos] = value;
    }
    
    void SetpairHeightStep(int pos, unsigned short value)
    {
        m_pairHeightStep[pos] = value;
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
    
    int& MaxCountMPs_pair()
    {
        return m_max_countMPs_pair;
    }
    
    int& RefImageID()
    {
        return m_ref_imageID;
    }
    
    short& cal(int pos)
    {
        return m_cal[pos];
    }
    
    unsigned char& Cloud(int pos)
    {
        return m_cloud[pos];
    }
    
    D2DPOINT& RBias(int pos)
    {
        return m_RBias[pos];
    }
    
    float& Tz(int pos)
    {
        return m_Tz[pos];
    }
    
    float& Tz_var(int pos)
    {
        return m_Tz_var[pos];
    }
    
    float& Tz_res(int pos)
    {
        return m_Tz_res[pos];
    }
    
    D2DPOINT& referencepairs(int pos)
    {
        return m_referencepairs[pos];
    }
    
    uint16& oripairnumber(int pos)
    {
        return m_oripairnumber[pos];
    }
    
    unsigned char& linked_pair(int pos, short query_pair)
    {
        return m_linked_pair[pos][query_pair];
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
    
    float& Azimuth(int pos)
    {
        return m_AzimuthDiff[pos];
    }
    
    unsigned char DiffDay(int pos)
    {
        return m_diff_day[pos];
    }
    
    unsigned char DiffTime(int pos)
    {
        return m_diff_time[pos];
    }
    
    float MatchingP(int pos)
    {
        return m_MPs_P[pos];
    }
    
    float Resolution(int pos)
    {
        return m_resolution[pos];
    }
    
    float& AE(int pos)
    {
        return m_AE[pos];
    }
    
    float& BIE(int pos)
    {
        return m_BIE[pos];
    }
    
    D3DPOINT& BaseRay(int pos)
    {
        return m_BaseRay[pos];
    }
    
    
    float& AE_IRPC(int pos)
    {
        return m_AE_IRPC[pos];
    }
    
    float& BIE_IRPC(int pos)
    {
        return m_BIE_IRPC[pos];
    }
    
    float& ConvergenceAngle_IRPC(int pos)
    {
        return m_ConvergenceAngle_IRPC[pos];
    }
    
    D3DPOINT& BaseRay_IRPC(int pos)
    {
        return m_BaseRay_IRPC[pos];
    }
    
    float& ConvergenceAngle_given(int pos)
    {
        return m_ConvergenceAngle_given[pos];
    }
    
    float& ConvergenceAngle_EQ(int pos)
    {
        return m_ConvergenceAngle_EQ[pos];
    }
    
    float& ConvergenceAngle_EQ_IRPC(int pos)
    {
        return m_ConvergenceAngle_EQ_IRPC[pos];
    }
    
    unsigned char& Coverage(int pos)
    {
        return m_Coverage[pos];
    }
    
    unsigned short& pairHeightStep(int pos)
    {
        return m_pairHeightStep[pos];
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
    
    double k1, k2, k3;//ratdial distortion
    double p1, p2;//tangential distortion
    double a1, a2;//affinity
    
    tagCameraInfo()
    {
        m_focalLength = 0;
        m_ImageSize.width = 0;
        m_ImageSize.height = 0;
        m_CCDSize = 0;
        m_ppx = 0;
        m_ppy = 0;
        
        k1 = k2 = k3 = 0;//ratdial distortion
        p1 = p2 = 0;//tangential distortion
        a1 = a2 = 0;//affinity
    };
    
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
    double CA_max_th;
    double Max_daygap;
    int pair_max_th;
    int pair_Azimuth_th;
    int pair_time_dif;
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
    int pair_options;
    int awnccmp;
    int merging_option;
    int NR_level;
    int NR_kernel_size;
    bool check_awncc;;
    int check_Planet_RA;
    int Planet_VC_level;
    int NDVI_th;
    int NDVI_pair_number;
    double NR_noise_var;
    
    
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
    char pre_shiftX[500];
    char pre_shiftY[500];
    char pre_roh[500];
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
	bool IsRA, IsSP, IsRR, IsSaveStep, Overall_DEM, Affine_RA, pre_DEMtif, check_tile_array, pre_SDM;
    bool check_Matchtag;
    bool check_selected_image[MaxImages];
    bool check_full_cal;
    bool check_pairinfo_only;
    
    bool check_NDVIPan;
    bool check_NDVI_pair_nochange;
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
    double CA_max_th;
    double Max_daygap;
    int pair_max_th;
    int pair_Azimuth_th;
    int pair_time_dif;
    double Cloud_th;
    double sim_shiftX;
    double sim_shiftY;
    double gamma;
    int min_S;
    int max_S;
    int NDVI_th;
    
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
    uint8 Coreg_PL; // 2 is default.
    uint8 SDM_SS;
    int DS_kernel;
    int RA_line_count;
    int RA_sample_count;
    double NR_noise_var;
    int NR_kernel_size;
    
    int Phi_start;
    int Phi_end;
    double Phi_interval;
    int Kappa_start;
    int Kappa_end;
    double Kappa_interval;
    int omega_start;
    int omega_end;
    double omega_interval;
    double sim_scale_start;
    double sim_scale_end;
    double sim_scale_interval;
    
    char Image[MaxImages][500];
    char Outputpath[500];
	char Outputpath_name[500];
	char seedDEMfilename[500];
	char metafilename[500];
    char EO_Path[500];
    char DEM_input_file[500];
    char Multi_input_file[500];
    char TIF_lists[500];
    
    bool check_minmaxStretch;
    bool check_gamma;
    bool check_DS_txy;
    bool check_downsample;
    int check_NDVIPAN;
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
    bool check_simulate;
    bool check_sdm_file;
    bool check_pairinfo_only;
    
    int check_txt_input;
    int check_coreg;
    int check_sdm_ortho;
    bool check_sdm_seed;
    int check_DEM_coreg_output;
    int pair_options;
    int awnccmp;
    int merging_option;
    int NR_level;
    bool check_awncc;
    int check_Planet_RA;
    int Planet_VC_level;
    
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
    
    DATE start_date;
} ARGINFO;

typedef struct tagImageGSD
{
    float row_GSD;
    float col_GSD;
    float pro_GSD;
} ImageGSD;

typedef struct tagBandInfo
{
    double abscalfactor;
    double effbw;
    double tdi;
    
    //multiband info
    double abscalfactor_multi[8];
    double effbw_multi[8];
    double tdi_multi[8];

    double calibrated_abscal_multi[8];
    double calibrated_effbw_multi[8];
} BandInfo;

typedef struct tagEsunDictWV02
{
    // Spectral Irradiance in W/m2/um
	// (from Thuillier 2003 - used by DG calibration team as of 2016)
    float BAND_P = 1571.36;
	float BAND_C = 1773.81;
    float BAND_B = 2007.27;
	float BAND_G = 1829.62;
	float BAND_Y = 1701.85;
	float BAND_R = 1538.85;
	float BAND_RE = 1346.09;
	float BAND_N = 1053.21;
	float BAND_N2 = 856.599;
} EsunDictWV02;

typedef struct tagGainDictWV02
{
    // Spectral Irradiance in W/m2/um
    float BAND_P = 0.942;
	float BAND_C = 1.151;
	float BAND_B = 0.988;
	float BAND_G = 0.936;
	float BAND_Y = 0.949;
	float BAND_R = 0.952;
	float BAND_RE = 0.974;
	float BAND_N = 0.961;
	float BAND_N2 = 1.002;
} GainDictWV02;

typedef struct tagBiasDictWV02
{
    // Spectral Irradiance in W/m2/um
    float BAND_P = -2.704;
	float BAND_C = -7.478;
	float BAND_B = -5.736;
	float BAND_G = -3.546;
	float BAND_ = -3.564;
	float BAND_R = -2.512;
	float BAND_RE = -4.120;
	float BAND_N = -3.300;
	float BAND_N2 = -2.891;
} BiasDictWV02;


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
    D2DPOINT min_XY;
    D2DPOINT max_XY;
    double Center[2];
    float convergence_angle;
    float convergence_angle_json;
    float AZ_ray[2]; //1 for IRPC vector, 2 for EO vector
    float EL_ray[2];
    
    ImageGSD GSD;
    
	int month;
    int date;
    int year;
    int hour;
    int min;
    int scandirection;
    int strip_ID;
    
    char filename[500];
    char fullpath[500];
    char Ofilename[500];
    char imagetime[500];
    char SatID[500];
    
    bool check_EO;
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
    uint8 mask;
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
    double **ImageAdjust;
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
    int max_covergae_pair;
    bool check_Height_update;
    double anchor_nccth;
    double blunder_nccth;
    double GNCC_th;
    
    GNCCGrid *GNCCGrid;
} LevelInfo;

typedef struct PairCA
{
    short pair_ID;
    float CA;
    float cloud;
public:
    PairCA()
    {
        
    }
    
    PairCA(int pairid, float ca, float cl)
    {
        pair_ID = pairid;
        CA = ca;
        cloud  = cl;
    }
} PairCA;

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

typedef MixedHeight3DGrid<float> SumCostContainer;



// Used in vectors when thread safety is needed
// just replace std::vector<bool> with std::vector<ConcurrentBool>
struct ConcurrentBool {
    bool val;

    ConcurrentBool(bool val) : val(val) {}
    operator bool() const { return val; }
};

/** Compute a running average of floating point values
 *      Example usage:
 *          RunningAverage avg;
 *          avg.update(2.0);
 *          avg.update(6.0);
 *          avg.update(120.0);
 *          // will print "average is 42.6667"
 *          std::cout << "average is " << avg.average() << std::endl;
 * */
class RunningAverage {
    public:
        RunningAverage() : n(0), v(0) {}

        /* Returns true if no values added, false otherwise */
        bool is_empty()
        {
            return n == 0;
        }

        /* Returns average of values if not empty,
         * throws domain_error otherwise */
        float average()
        {
            if(!n)
            {
                throw std::domain_error("average of empty list");
            }
            return v / n;
        }

        /* Add value to the running average */
        void update(float x) {
            ++n;
            v += x;
        };

    private:
        float v;
        int n;

};


#endif

