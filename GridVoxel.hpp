#ifndef GRID_VOXEL_H
#define GRID_VOXEL_H

#include <vector>

inline short DoubleToSignedChar_voxel(double val)
{
    return (short)(val*1000.0);
}

inline double SignedCharToDouble_voxel(short val)
{
    return (double)(val)/1000.0;
}



/*
 * Allocated like: 
        grid_voxel = (VOXEL**)calloc(sizeof(VOXEL*),Size_Grid2D.width*Size_Grid2D.height);
 *
 * That is, each point in the voxel grid is a VOXEL * (VOXEL pointer)
 * Passed to InitializeVoxel. First thing this function does is free older
 * version of the voxel arrays:
 *
        free(grid_voxel[t_i]);
  
 * Then, allocates them:
        grid_voxel[t_i] = (VOXEL*)calloc(sizeof(VOXEL),NumberofHeightVoxel);
 *
 * And, for each height index, updates the INCC and flag_cal:
        for(int h = 0 ; h < NumberofHeightVoxel ; h++)
        {
            for(int pair = 0 ; pair < plevelinfo.pairinfo->NumberOfPairs; pair++)
            {
                grid_voxel[t_i][h].flag_cal[pair] = false;
                grid_voxel[t_i][h].INCC[pair] = DoubleToSignedChar_voxel(-1);
            }
        }
 *
 *
 * Used in AWNCC_[AWNCC, SGM, single, multi), InitializeVoxel,
 * SGM_[con, start]_pos, VerticalLineLocus, SignedCharToDouble_voxel
 */

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



class VoxelTower {
public:
    VoxelTower(int max_ncc) : ncc_len(max_ncc) {}

    bool& flag_cal(size_t h_index, size_t index) {
        return _flag_cal[h_index * ncc_len + index].val;
    }
    short& INCC(size_t h_index, size_t index) {
        return _INCC[h_index * ncc_len + index];
    }

    // const versions of the above
    const bool& flag_cal(size_t h_index, size_t index) const {
        return _flag_cal[h_index * ncc_len + index].val;
    }
    const short& INCC(size_t h_index, size_t index) const {
        return _INCC[h_index * ncc_len + index];
    }

    /** Set size to zero and clear memory */
    void clear() {
        std::vector<BoolWrapper>().swap(_flag_cal);
        std::vector<short>().swap(_INCC);
    }

    /** Allocate n elements.
     *
     * Initialize flag_cal values to false.
     * Initialize INCC values to -1
     */
    void allocate(size_t n) {
        _flag_cal = std::vector<BoolWrapper>(n * ncc_len, BoolWrapper(true));
        _INCC = std::vector<short>(n * ncc_len, DoubleToSignedChar_voxel(-1));
    }


private:
    struct BoolWrapper {
        BoolWrapper(bool b) { val = b; }
        bool val;
    };

    // Cannot use bool here, it's not thread safe. Need to wrap it instead.
    std::vector<BoolWrapper> _flag_cal;
    std::vector<short> _INCC;
    size_t ncc_len;

};

class GridVoxel {
public:
    GridVoxel(size_t length, int max_ncc) : towers(length, VoxelTower(max_ncc)) {}
    VoxelTower & operator[](size_t n) {
        return towers[n];
    }
    const VoxelTower & operator[](size_t n) const {
        return towers[n];
    }
    void clearall() {
        std::vector<VoxelTower>().swap(towers);
    }

private:
    std::vector<VoxelTower> towers;
};

#endif

