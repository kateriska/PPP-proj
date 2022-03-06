/**
 * @file    parallel_heat_solver.cpp
 * @author  xlogin00 <xlogin00@stud.fit.vutbr.cz>
 *
 * @brief   Course: PPP 2020/2021 - Project 1
 *
 * @date    2021-MM-DD
 */

#include "parallel_heat_solver.h"

using namespace std;

//============================================================================//
//                            *** BEGIN: NOTE ***
//
// Implement methods of your ParallelHeatSolver class here.
// Freely modify any existing code in ***THIS FILE*** as only stubs are provided
// to allow code to compile.
//
//                             *** END: NOTE ***
//============================================================================//

ParallelHeatSolver::ParallelHeatSolver(SimulationProperties &simulationProps,
                                       MaterialProperties &materialProps)
    : BaseHeatSolver (simulationProps, materialProps),
     m_fileHandle(H5I_INVALID_HID, static_cast<void (*)(hid_t )>(nullptr))
{
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    // Creating EMPTY HDF5 handle using RAII "AutoHandle" type
    //
    // AutoHandle<hid_t> myHandle(H5I_INVALID_HID, static_cast<void (*)(hid_t )>(nullptr))
    //
    // This can be particularly useful when creating handle as class member!
    // Handle CANNOT be assigned using "=" or copy-constructor, yet it can be set
    // using ".Set(/* handle */, /* close/free function */)" as in:
    // myHandle.Set(H5Fopen(...), H5Fclose);

    // Requested domain decomposition can be queried by
    // m_simulationProperties.GetDecompGrid(/* TILES IN X */, /* TILES IN Y */)
    // 1. Open output file if its name was specified.

    if (m_rank == 0 && simulationProps.IsUseParallelIO() == false)
    {
    if (!m_simulationProperties.GetOutputFileName().empty())
    {
        m_fileHandle.Set(H5Fcreate(simulationProps.GetOutputFileName("par").c_str(),H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    }
    }

}

ParallelHeatSolver::~ParallelHeatSolver()
{

}

void ParallelHeatSolver::RunSolver(std::vector<float, AlignedAllocator<float> > &outResult)
{
    // UpdateTile(...) method can be used to evaluate heat equation over 2D tile
    //                 in parallel (using OpenMP).
    // NOTE: This method might be inefficient when used for small tiles such as
    //       2xN or Nx2 (these might arise at edges of the tile)
    //       In this case ComputePoint may be called directly in loop.

    // ShouldPrintProgress(N) returns true if average temperature should be reported
    // by 0th process at Nth time step (using "PrintProgressReport(...)").

    // Finally "PrintFinalReport(...)" should be used to print final elapsed time and
    // average temperature in column.
    //std::vector<float> local_a(local_size);
      //local_a.reserve(local_size);
    //std::vector<float> loca_b(local_size);

    int outSizeX;
    int outSizeY;
    int local_tile_size_x;
    int local_tile_size_y;
    int local_tile_size;

    int domain_length;
    //vector<float> init_temp(m_materialProperties.GetInitTemp().size());

    if (m_rank == 0)
    {
      m_simulationProperties.GetDecompGrid(outSizeX, outSizeY);
      local_tile_size_x = sqrt(m_materialProperties.GetInitTemp().size()) / outSizeX;
      local_tile_size_y = sqrt(m_materialProperties.GetInitTemp().size()) / outSizeY;
      cout << local_tile_size_x << endl;
      cout << local_tile_size_y << endl;
      local_tile_size = local_tile_size_x * local_tile_size_y;
      domain_length = m_materialProperties.GetInitTemp().size();

      // compute position of each process in grid
      // 0 - left upper corner

    }
    MPI_Bcast(&outSizeX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&outSizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&domain_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << local_tile_size << endl;

    vector<float> init_temp;
    vector<float> tmp_vector;
    vector<float> out_result;
    vector<int> domain_map;
    vector<float> domain_params;

    if (m_rank == 0)
    {
    for (size_t i = 0; i < m_materialProperties.GetInitTemp().size(); ++i) {
        init_temp.push_back(m_materialProperties.GetInitTemp().at(i));
    }
    for (size_t i = 0; i < m_materialProperties.GetDomainParams().size(); ++i) {
        domain_params.push_back(m_materialProperties.GetDomainParams().at(i));
    }

    for (size_t i = 0; i < m_materialProperties.GetDomainMap().size(); ++i) {
        domain_map.push_back(m_materialProperties.GetDomainMap().at(i));
    }
    }


    string vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < domain_map.size(); ++i)
    {
      int item = domain_map.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }


    cout << vector_res << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    vector<float> init_temp_local(local_tile_size);
    vector<float> tmp_vector_local(local_tile_size);
    vector<float> out_result_local(local_tile_size);
    vector<int> domain_map_local(local_tile_size);
    vector<float> domain_params_local(local_tile_size);

    cout << init_temp_local.data() <<endl;
    MPI_Scatter(init_temp.data(), local_tile_size, MPI_FLOAT,  init_temp_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(domain_params.data(), local_tile_size, MPI_FLOAT,  domain_params_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(domain_map.data(), local_tile_size, MPI_INT,  domain_map_local.data(), local_tile_size, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(tmp_vector.data(), local_tile_size, MPI_FLOAT,  tmp_vector_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(out_result.data(), local_tile_size, MPI_FLOAT,  out_result_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);


    vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < domain_map_local.size(); ++i)
    {
      int item = domain_map_local.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }
    cout << vector_res << endl;

    if (m_rank == 0)
    {
      double startTime = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // 2. Copy initial temperature into both working arrays
    /*
    std::copy(init_temp_local.begin(),
              init_temp_local.end(), m_tempArray.begin());
    std::copy(init_temp_local.begin(),
              init_temp_local.end(), outResult.begin());

    float *workTempArrays[] = { m_tempArray.data(), outResult.data() };
    */
    
    /*
    for(size_t iter = 0; iter < m_simulationProperties.GetNumIterations(); ++iter)
    {
        // 4. Compute new temperature for each point in the domain (except borders)
        // border temperatures should remain constant (plus our stencil is +/-2 points).
        for(size_t i = 2; i < m_materialProperties.GetEdgeSize() - 2; ++i)
        {
            for(size_t j = 2; j < m_materialProperties.GetEdgeSize() - 2; ++j)
            {
                ComputePoint(workTempArrays[1], workTempArrays[0],
                        m_materialProperties.GetDomainParams().data(),
                        m_materialProperties.GetDomainMap().data(),
                        i, j,
                        m_materialProperties.GetEdgeSize(),
                        m_simulationProperties.GetAirFlowRate(),
                        m_materialProperties.GetCoolerTemp());
            }
        }
        */


    /*
    int local_tile_size_x;
    int local_tile_size_y;
    int local_tile_size;

    //int local_tile_size_x = m_materialProperties.GetInitTemp().size() / outSizeX;
    //int local_tile_size_y =  m_materialProperties.GetInitTemp().size() / outSizeY;
    //int local_tile_size = local_tile_size_x * local_tile_size_y;

    cout << m_materialProperties.GetInitTemp().size() << endl;
    local_tile_size_x = sqrt(m_materialProperties.GetInitTemp().size()) / outSizeX;
    local_tile_size_y = sqrt(m_materialProperties.GetInitTemp().size()) / outSizeY;
    local_tile_size = local_tile_size_x * local_tile_size_y;
    cout << local_tile_size << endl;




    //cout << local_tile_size << endl;
    //vector<float> local_tmp(outSizeX * outSizeY);
    //vector<float> local_out(outSizeX * outSizeY);

    vector<float> init_temp;
    vector<float> tmp_vector(m_materialProperties.GetInitTemp().size());
    vector<float> out_result(m_materialProperties.GetInitTemp().size());
    vector<int> domain_map(m_materialProperties.GetDomainMap().size());
    vector<float> domain_params(m_materialProperties.GetDomainParams().size());

    vector<float> init_temp_local(local_tile_size);
    vector<float> tmp_vector_local(local_tile_size);
    vector<float> out_result_local(local_tile_size);
    vector<int> domain_map_local(local_tile_size);
    vector<float> domain_params_local(local_tile_size);

    //vector<float> init_temp_gathered(m_materialProperties.GetInitTemp().size());


    //cout << init_temp_local.size() <<endl;



    if (m_rank == 0)
    {
    for (size_t i = 0; i < m_materialProperties.GetInitTemp().size(); ++i) {
        init_temp.push_back(m_materialProperties.GetInitTemp().at(i));
    }

    cout << init_temp.size() << endl;

    for (size_t i = 0; i < init_temp.size(); ++i) {
        std::cout << init_temp.at(i) << "; ";
    }

    for (size_t i = 0; i < m_materialProperties.GetDomainParams().size(); ++i) {
        domain_params.push_back(m_materialProperties.GetDomainParams().at(i));
    }

    for (size_t i = 0; i < m_materialProperties.GetDomainMap().size(); ++i) {
        domain_map.push_back(m_materialProperties.GetDomainMap().at(i));
    }
    cout << "End of process 0" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    /*
    cout << m_materialProperties.GetInitTemp().size() << endl;
    cout << m_materialProperties.GetDomainMap().size() << endl;
    cout << m_materialProperties.GetDomainParams().size() << endl;


    //int local_size = m_materialProperties.GetInitTemp().size() / m_size;
    cout << "Call Scatter" << endl;
    cout << init_temp_local.data() << endl;
    int res = MPI_Scatter(init_temp.data(), local_tile_size, MPI_FLOAT,  init_temp_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //cout << res << endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    /*
    MPI_Scatter(domain_map.data(), local_tile_size, MPI_INT,  domain_map_local.data(), local_tile_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(domain_params.data(), local_tile_size, MPI_FLOAT,  domain_params.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(tmp_vector.data(), local_tile_size, MPI_FLOAT,  tmp_vector_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(out_result.data(), local_tile_size, MPI_FLOAT,  out_result_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //cout << init_temp_local.size() << endl;

    //MPI_Barrier(MPI_COMM_WORLD);
    string vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < init_temp_local.size(); ++i)
    {
      int item = init_temp.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }
    cout << vector_res << endl;

    float *init_temp_gathered = NULL;
    if (m_rank == 0)
    {
      init_temp_gathered = malloc(sizeof(float) * m_size * m_materialProperties.GetInitTemp().size());
      //vector<float> init_temp_gathered(m_materialProperties.GetInitTemp().size());

      MPI_Gather(init_temp_local.data(), local_tile_size, MPI_FLOAT, init_temp_gathered, local_tile_size, MPI_FLOAT, 0,
           MPI_COMM_WORLD);

    }
    */



    /*
    std::copy(m_materialProperties.GetInitTemp().begin(),
                m_materialProperties.GetInitTemp().end(), m_tempArray.begin());
    std::copy(m_materialProperties.GetInitTemp().begin(),
                m_materialProperties.GetInitTemp().end(), outResult.begin());
    cout << m_materialProperties.GetInitTemp().size();





    if (m_simulationProperties.GetDecompMode() == DECOMP_MODE_1D)
    {

    }
    else
    {

    }
    */





}
