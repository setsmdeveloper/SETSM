#ifndef GRID_VOXEL_H
#define GRID_VOXEL_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include "log.hpp"
#include <cstdlib>
#include <stdio.h>
#include "Typedefine.hpp"

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
        printf("SignedCharToDouble_voxel overflow %d\n",val);
        exit(1);
    }
    return _SignedCharToDouble_voxel(val);
}

class GridVoxel {
public:

    // Constructor
    GridVoxel(size_t length, const GridPairs &Gridpairs) 
            : incc_grid(length), Gridpairs(Gridpairs) {}

    // Const and nonconst accessors for INCC
    short& INCC(size_t grid_index, size_t h_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(grid_index, pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        size_t num_pairs = pairs.size();
        return incc_grid[grid_index][h_index * num_pairs + pair_index];
    }

    // const versions of the above
    const short& INCC(size_t grid_index, size_t h_index, size_t pair_id) const {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(grid_index, pairs, pair_id);
        if(pair_index < 0)
            dbg_crash(grid_index, pairs, pair_id);
        size_t num_pairs = pairs.size();
        return incc_grid[grid_index][h_index * num_pairs + pair_index];
    }

    bool is_cal(size_t grid_index, size_t h_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        int pair_index = pair_id_to_index(grid_index, pairs, pair_id);
        if(pair_index < 0)
            return false;
        size_t num_pairs = pairs.size();
        return incc_grid[grid_index][h_index * num_pairs + pair_index] != INCC_UNSET;
    }

    bool has_pair(size_t grid_index, size_t pair_id) {
        auto &pairs = Gridpairs.get_pairs(grid_index);
        return pair_id_to_index(grid_index, pairs, pair_id) >= 0;
    }

    /** Set size to zero and clear memory */
    void clear(size_t grid_index) {
        std::vector<short>().swap(incc_grid[grid_index]);
    }

    /** Allocate n elements.
     *
     * Initialize INCC values to -1
     */
    void allocate(size_t grid_index, size_t n) {
        auto num_pairs = Gridpairs.get_pairs(grid_index).size();
        incc_grid[grid_index] = std::vector<short>(n * num_pairs, INCC_UNSET);
    }

private:
    std::vector<vector<short>> incc_grid;
    GridPairs Gridpairs;

    static constexpr short INCC_UNSET = _DoubleToSignedChar_voxel(-1);

    int pair_id_to_index(size_t grid_index, const vector<short> &pairs, size_t pair_id) const {
        auto it = std::find(pairs.begin(), pairs.end(), pair_id);
        if(it == pairs.end()) {
            return -1;
        }
        return std::distance(pairs.begin(), it);
    }

    void dbg_crash(size_t grid_index, const vector<short> &pairs, size_t pair_id) const {
        LOG("Failed to find pair id %lu for grid index %lu\n", pair_id, grid_index);
        LOG("number of pairs: %d\n", pairs.size());

        LOG("pairs:\n");
        for(int i = 0; i < pairs.size(); i++) {
            LOG("    pair[%d] = %d\n", i, pairs[i]);
        }
        abort();
    }
};

#endif
