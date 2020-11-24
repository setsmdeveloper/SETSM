#ifndef GRID_VOXEL_H
#define GRID_VOXEL_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cstdlib>

inline short DoubleToSignedChar_voxel(double val)
{
    return (short)(val*1000.0);
}

inline double SignedCharToDouble_voxel(short val)
{
    return (double)(val)/1000.0;
}

class VoxelTower {
public:
    VoxelTower() : _num_pairs(0) {}

    bool& flag_cal(size_t h_index, int pair_id) {
        return _flag_cal[h_index * _num_pairs + get_index(pair_id)].val;
    }
    short& INCC(size_t h_index, int pair_id) {
        return _INCC[h_index * _num_pairs + get_index(pair_id)];
    }

    // const versions of the above
    const bool& flag_cal(size_t h_index, int pair_id) const {
        return _flag_cal[h_index * _num_pairs + get_index(pair_id)].val;
    }
    const short& INCC(size_t h_index, int pair_id) const {
        return _INCC[h_index * _num_pairs + get_index(pair_id)];
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
    void allocate(size_t n, const std::vector<uint8_t> &pairs) {
        _pairs = pairs;
        _num_pairs = _pairs.size();
        _flag_cal = std::vector<BoolWrapper>(n * _num_pairs, BoolWrapper(true));
        _INCC = std::vector<short>(n * _num_pairs, DoubleToSignedChar_voxel(-1));
    }


private:
    int get_index(int pair_id) const {
        auto it = std::find(_pairs.begin(), _pairs.end(), pair_id);
        if(it == _pairs.end()) {
            abort();
        }
        return std::distance(_pairs.begin(), it);
    }

    struct BoolWrapper {
        BoolWrapper(bool b) { val = b; }
        bool val;
    };

    // Cannot use bool here, it's not thread safe. Need to wrap it instead.
    std::vector<BoolWrapper> _flag_cal;
    std::vector<short> _INCC;
    std::vector<uint8_t> _pairs;
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
