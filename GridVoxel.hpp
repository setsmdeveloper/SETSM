#ifndef GRID_VOXEL_H
#define GRID_VOXEL_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include "log.hpp"
#include <cstdlib>
#include <stdio.h>

constexpr inline short _DoubleToSignedChar_voxel(double val)
{
    return (short)(val*1000.0);
}

constexpr inline double _SignedCharToDouble_voxel(short val)
{
    return (double)(val)/1000.0;
}

inline short DoubleToSignedChar_voxel(double val)
{
    if(val > 30 || val < -30)
    {
        printf("DoubleToSignedChar_voxel overflow %f\n",val);
        exit(1);
    }
    return _DoubleToSignedChar_voxel(val);
}

inline double SignedCharToDouble_voxel(short val)
{
    if(val > 30000 || val < -30000)
    {
        printf("SignedCharToDouble_voxel overflow %f\n",val);
        exit(1);
    }
    return _SignedCharToDouble_voxel(val);
}

class VoxelTower {
public:
    VoxelTower() : _num_pairs(0) {}

    short& INCC(size_t h_index, int pair_id) {
        int pair_index = get_index(pair_id);
        if(pair_index < 0)
            dbg_crash(pair_id);
        return _INCC[h_index * _num_pairs + pair_index];
    }

    // const versions of the above
    const short& INCC(size_t h_index, int pair_id) const {
        int pair_index = get_index(pair_id);
        if(pair_index < 0)
            dbg_crash(pair_id);
        return _INCC[h_index * _num_pairs + pair_index];
    }

    bool is_cal(size_t h_index, int pair_id) {
        int pair_index = get_index(pair_id);
        if(pair_index < 0)
            return false;
        return _INCC[h_index * _num_pairs + pair_index] != INCC_UNSET;
    }

    bool has_pair(int pair_id) {
        return get_index(pair_id) >= 0;
    }

    /** Set size to zero and clear memory */
    void clear() {
        std::vector<short>().swap(_INCC);
    }

    /** Allocate n elements.
     *
     * Initialize INCC values to -1
     */
    void allocate(size_t n, const std::vector<short> &pairs) {
        _pairs = pairs;
        _num_pairs = _pairs.size();
        _INCC = std::vector<short>(n * _num_pairs, INCC_UNSET);
    }


private:
    int get_index(int pair_id) const {
        auto it = std::find(_pairs.begin(), _pairs.end(), pair_id);
        if(it == _pairs.end()) {
            return -1;
        }
        return std::distance(_pairs.begin(), it);
    }

    void dbg_crash(int pair_id) const {
        LOG("Failed to find pair id %d\n", pair_id);
        LOG("number of pairs: %d\n", _pairs.size());
        LOG("pairs:\n");
        for(int i = 0; i < _pairs.size(); i++) {
            LOG("    pair[%d] = %d\n", i, _pairs[i]);
        }
        abort();
    }

    struct BoolWrapper {
        BoolWrapper(bool b) { val = b; }
        bool val;
    };

    static constexpr short INCC_UNSET = _DoubleToSignedChar_voxel(-1);
    // Cannot use bool here, it's not thread safe. Need to wrap it instead.
    std::vector<short> _INCC;
    std::vector<short> _pairs;
    size_t _num_pairs;

};

// Thin wrapper around vector, passes through the max_ncc value
class GridVoxel {
public:
    typedef std::vector<VoxelTower>::size_type size_type;

    GridVoxel(size_t length) : towers(length) {}

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
