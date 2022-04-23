/**
 * @file    parallel_heat_solver.h
 * @author  xforto00 <xforto00@stud.fit.vutbr.cz>
 *
 * @brief   Course: PPP 2021/2022 - Project 1
 *
 * @date    2022-04-11
 */

#ifndef PARALLEL_HEAT_SOLVER_H
#define PARALLEL_HEAT_SOLVER_H

#include "base_heat_solver.h"
#include <cmath>
#include <string>
#include <algorithm>
#include <assert.h>
#include <bits/stdc++.h>
#include <list>
#include <vector>

using namespace std;

/**
 * @brief The ParallelHeatSolver class implements parallel MPI based heat
 *        equation solver in 2D using 1D and 2D block grid decomposition.
 */
class ParallelHeatSolver : public BaseHeatSolver
{
    //============================================================================//
    //                            *** BEGIN: NOTE ***
    //
    // Modify this class declaration as needed.
    // This class needs to provide at least:
    // - Constructor which passes SimulationProperties and MaterialProperties
    //   to the base class. (see below)
    // - Implementation of RunSolver method. (see below)
    //
    // It is strongly encouraged to define methods and member variables to improve
    // readability of your code!
    //
    //                             *** END: NOTE ***
    //============================================================================//

public:
    /**
     * @brief Constructor - Initializes the solver. This should include things like:
     *        - Construct 1D or 2D grid of tiles.
     *        - Create MPI datatypes used in the simulation.
     *        - Open SEQUENTIAL or PARALLEL HDF5 file.
     *        - Allocate data for local tile.
     *        - Initialize persistent communications?
     *
     * @param simulationProps Parameters of simulation - passed into base class.
     * @param materialProps   Parameters of material - passed into base class.
     */
    ParallelHeatSolver(SimulationProperties &simulationProps, MaterialProperties &materialProps);
    virtual ~ParallelHeatSolver();

    /**
     * @brief Run main simulation loop.
     * @param outResult Output array which is to be filled with computed temperature values.
     *                  The vector is pre-allocated and its size is given by dimensions
     *                  of the input file (edgeSize*edgeSize).
     *                  NOTE: The vector is allocated (and should be used) *ONLY*
     *                        by master process (rank 0 in MPI_COMM_WORLD)!
     */
    virtual void RunSolver(std::vector<float, AlignedAllocator<float> > &outResult);

    //vector<float> ParallelHeatSolver::SplitVector(vector<float> input_vector, int n, int size)
    void print_array(int* arr, int width, int height);
    void print_array(float* arr, int width, int height);
    list<vector<int>> SplitRows(int *input_arr, int local_tile_size, int local_tile_size_cols);
    vector<int> EnlargeTile(list<vector<int>> input_list, int local_tile_size_cols);
    list<vector<float>> SplitRows(float *input_arr, int local_tile_size, int local_tile_size_cols);
    vector<float> EnlargeTile(list<vector<float>> input_list, int local_tile_size_cols);
    vector<int> Enlarge1DTile(int *input_arr, int local_tile_size);
    vector<float> Enlarge1DTile(float *input_arr, int local_tile_size);
    float ComputeMiddleColAvgTemp(const float *data, int enlarged_tile_size_rows, int enlarged_tile_size_cols, int middle_item_tile_col_id) const;
    int count_1D_index(int row, int length_of_row, int column);
    float* TrimTileWithoutBorders(float* arr, int enlarged_tile_size_rows, int enlarged_tile_size_cols, float* result);
    void UpdateTileHybrid(const float *oldTemp, float *newTemp, const float *params, const int *map, size_t offset_rows_begin, size_t offset_rows_end, size_t offset_cols_begin, size_t offset_cols_end, size_t enlarged_tile_size_rows, size_t enlarged_tile_size_cols, float airFlowRate, float coolerTemp) const;



protected:
    int m_rank;     ///< Process rank in global (MPI_COMM_WORLD) communicator.
    int m_size;     ///< Total number of processes in MPI_COMM_WORLD.

    AutoHandle<hid_t> m_fileHandle;                             ///< Output HDF5 file handle.
    std::vector<float, AlignedAllocator<float> > m_tempArray;
};

#endif // PARALLEL_HEAT_SOLVER_H
