#ifndef GRID_VOXEL_H
#define GRID_VOXEL_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>

inline short DoubleToSignedChar_voxel(double val)
{
    if(val > 30 || val < -30)
    {
        printf("DoubleToSignedChar_voxel overflow %f\n",val);
        exit(1);
    }
    return (short)(val*1000.0);
}

inline double SignedCharToDouble_voxel(short val)
{
    if(val > 30000 || val < -30000)
    {
        printf("SignedCharToDouble_voxel overflow %f\n",val);
        exit(1);
    }
    return (double)(val)/1000.0;
}

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

// Thin wrapper around vector, passes through the max_ncc value
class GridVoxel {
public:
    typedef std::vector<VoxelTower>::size_type size_type;

    GridVoxel(size_t length, int max_ncc) : towers(length, VoxelTower(max_ncc)) {}

    VoxelTower & operator[](size_type n) {
        return towers[n];
    }
    const VoxelTower & operator[](size_type n) const {
        return towers[n];
    }

    void clearall() {
        std::vector<VoxelTower>().swap(towers);
    }

private:
    std::vector<VoxelTower> towers;
};

#endif
