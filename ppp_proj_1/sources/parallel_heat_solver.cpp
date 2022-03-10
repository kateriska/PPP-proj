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
/*
vector<float> ParallelHeatSolver::SplitVector(vector<float> input_vector, int n, int size)
{




    // create an array of vectors to store the sub-vectors
    std::vector<float> vec[size];

    // each iteration of this loop process the next set of `n` elements
    // and store it in a vector at k'th index in `vec`
    for (int k = 0; k < size; ++k)
    {
        // get range for the next set of `n` elements
        auto start_itr = std::next(v.cbegin(), k*n);
        auto end_itr = std::next(v.cbegin(), k*n + n);

        // allocate memory for the sub-vector
        vec[k].resize(n);

        // copy elements from the input range to the sub-vector
        std::copy(start_itr, end_itr, vec[k].begin());
    }

    return vec;

}
*/

void ParallelHeatSolver::print_array(int* arr, int width, int height)
{
    int i;
    for (i = 0; i< width * height;i++)
    {
        if((i != 0) && (i % width == 0))
            printf("\n");
        printf("%4d ", arr[i]);
    }
    putchar('\n');
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
    }

    MPI_Bcast(&outSizeX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&outSizeY, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
        cout << m_materialProperties.GetInitTemp().at(i) << endl;
        float value = m_materialProperties.GetInitTemp().at(i);
        //float edited_value = round( value * 10.0 ) / 10.0;
        //cout << edited_value << endl;
        //string value_str = to_string(value);
        //float edited_value = stof(value_str);
        //cout << value_str << endl;
        //cout << value;
        init_temp.push_back(value);


      }

      //list<vector<float> >& init_temp_original(outSizeX);

      int n = local_tile_size * outSizeX;
      int size = outSizeY;

      // create an array of vectors to store the sub-vectors
    std::vector<float> init_temp_original[size];


    // each iteration of this loop process the next set of `n` elements
    // and store it in a vector at k'th index in `vec`
    // https://www.techiedelight.com/split-vector-into-subvectors-cpp/
    /*
    for (int k = 0; k < size; ++k)
    {
        // get range for the next set of `n` elements
        auto start_itr = std::next(init_temp.cbegin(), k*n);
        auto end_itr = std::next(init_temp.cbegin(), k*n + n);

        // allocate memory for the sub-vector
        init_temp_original[k].resize(n);

        // code to handle the last sub-vector as it might
        // contain fewer elements

        // copy elements from the input range to the sub-vector
        std::copy(start_itr, end_itr, init_temp_original[k].begin());

    }
    */

    for (size_t i = 0; i < m_materialProperties.GetDomainParams().size(); ++i) {
        domain_params.push_back(m_materialProperties.GetDomainParams().at(i));
    }

    for (size_t i = 0; i < m_materialProperties.GetDomainMap().size(); ++i) {
      int value = m_materialProperties.GetDomainMap().at(i);
        domain_map.push_back(value);
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

    int m_num;

    MPI_Comm_size(MPI_COMM_WORLD, &m_num);
    int displacements[m_num];
    int counts[m_num];
    fill_n(counts, m_num, 1);

    int res = 0;
    bool skip_to_next_row = false;

    if (outSizeX != 1 && outSizeY == 1)
    {
      for(int x = 0; x < outSizeX; x++)
      {
        if (x == 0)
        {
          res = 0;
        }
        else
        {
          cout << "Local tile size x " << local_tile_size_x << endl;
          res = res + local_tile_size_y;
        }
        displacements[x] = res;
        cout << "Result is " << res << endl;
      }
    }
    else {
    for(int x = 0; x < outSizeX; x++){
        for(int y = 0; y < outSizeY; y++){
            cout << x << ", " << y<< endl;
            cout << "Process on index" << x*outSizeX+y <<  endl;
            //int res =  x * sqrt(domain_length) * local_tile_size_x + y * local_tile_size_y;
            if (x == 0 & y == 0)
            {
              res = 0;
            }
            else if (skip_to_next_row == true)
            {
              res = local_tile_size_x * sqrt(domain_length);
              skip_to_next_row = false;
            }
            else
            {
              res =  res + local_tile_size_x;
            }

            displacements[x*outSizeX+y] = res;
            cout << "Result is " << res << endl;


        }
        skip_to_next_row = true;
        //res = sqrt(domain_length) * local_tile_size_y;
    }
  }



    for (auto d : counts)
{
    std::cout << d << std::endl;
}



  int enlarged_tile_size = (local_tile_size_x + 4) * (local_tile_size_y + 4);
  int enlarged_tile_size_x = local_tile_size_x + 4;
  int enlarged_tile_size_y = local_tile_size_y + 4;

  float *init_temp_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  int *domain_map_local = (int*)malloc(enlarged_tile_size * sizeof(int));
  float *domain_params_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  assert(domain_map_local != NULL);

  MPI_Datatype tile_t, resized_tile_t;
  MPI_Type_vector(local_tile_size_x, local_tile_size_y, sqrt(domain_length), MPI_FLOAT, &tile_t);
  MPI_Type_create_resized(tile_t, 0, sizeof(float), &resized_tile_t);
  MPI_Type_commit(&tile_t);
  MPI_Type_commit(&resized_tile_t);


    /*
    vector<float> init_temp_local(local_tile_size);
    vector<float> tmp_vector_local(local_tile_size);
    vector<float> out_result_local(local_tile_size);
    vector<int> domain_map_local(local_tile_size);
    vector<float> domain_params_local(local_tile_size);
    */

    //cout << init_temp_local.data() <<endl;
    //MPI_Scatter(init_temp.data(), local_tile_size, MPI_FLOAT,  init_temp_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(domain_params.data(), local_tile_size, MPI_FLOAT,  domain_params_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(domain_map.data(), local_tile_size, MPI_INT,  domain_map_local.data(), local_tile_size, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(tmp_vector.data(), local_tile_size, MPI_FLOAT,  tmp_vector_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(out_result.data(), local_tile_size, MPI_FLOAT,  out_result_local.data(), local_tile_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //float *init_temp_local;



    //int displacements[m_num] = {0,1,16,17};
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(&(init_temp[0]), counts, displacements, resized_tile_t, &(init_temp_local[0]), local_tile_size_x * local_tile_size_y, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_map[0]), counts, displacements, resized_tile_t, &(domain_map_local[0]), local_tile_size_x * local_tile_size_y, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_params[0]), counts, displacements, resized_tile_t, &(domain_params_local[0]), local_tile_size_x * local_tile_size_y, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < m_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (i == m_rank) {
            printf("\nRank: %d\n", m_rank);
            print_array(domain_map_local, local_tile_size_x, local_tile_size_y);
        }
    }

    vector<int> middle_ranks;
    int middle_column = outSizeY / 2;

    for (int i = 0; i < m_num; i++) {
      int column_rank = i % outSizeY;
      if (middle_column == column_rank || i == 0)
      {
        middle_ranks.push_back(i);
      }
    }



    vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < middle_ranks.size(); ++i)
    {
      int item = middle_ranks.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }


    cout << vector_res << endl;
    MPI_Group WORLD_GROUP;
    MPI_Group MIDDLE_COLUMN_GROUP;
    MPI_Comm MPI_COMM_MIDDLE_COLUMN;


    MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP);
    MPI_Group_incl(WORLD_GROUP, middle_ranks.size(), middle_ranks.data(), &MIDDLE_COLUMN_GROUP);
    MPI_Comm_create(MPI_COMM_WORLD, MIDDLE_COLUMN_GROUP, &MPI_COMM_MIDDLE_COLUMN);

    cout << "Hhhh" << endl;
    cout << enlarged_tile_size_x << "; " << enlarged_tile_size_y << endl;;
    vector<int> domain_map_local_with_borders(enlarged_tile_size);
    vector<float> domain_params_with_borders(enlarged_tile_size);
    cout << enlarged_tile_size << endl;


    for (int i = 0; i < local_tile_size_y; i++) {
        domain_map_local_with_borders.insert(domain_map_local_with_borders.begin() + ((i + 2) * enlarged_tile_size_x) + 2, &domain_map_local[i * local_tile_size_x],&domain_map_local[i * local_tile_size_x + local_tile_size_x]);
        //block_b_map.insert(block_b_map.begin() + ((i + 2) * b_block_w) + 2, &block_map[i * _block_w],&block_map[i * _block_w + _block_w]);
    }

    list<vector<int>> domain_map_local_list;

    bool start_new_vector = false;
    vector<int> domain_map_local_row(local_tile_size_x);
    for (size_t i = 0; i < local_tile_size; ++i)
    {
      cout << local_tile_size << endl;
      int item = domain_map_local[i];
      if (start_new_vector == true)
      {
        cout << "Starting new vector " << endl;
        domain_map_local_row.clear();
        start_new_vector = false;
      }
      else
      {
        domain_map_local_row.push_back(item);
      }
      if (i % local_tile_size_x == 1 && i != 1)
      {
        start_new_vector = true;
        domain_map_local_list.push_back(domain_map_local_row);
      }

    }

    for (auto vect : domain_map_local_list) {
        // Each element of the list is
        // a vector itself
        vector<int> currentVector = vect;

        cout << "[ ";

        // Printing vector contents
        for (auto element : currentVector)
            cout << element << ' ';

        cout << ']';
        cout << '\n';
    }






    vector_res = "";
    cout << m_rank << endl;
    cout << domain_map_local_with_borders.size() << endl;
    for (size_t i = 0; i < domain_map_local_with_borders.size(); ++i)
    {
      int item = domain_map_local_with_borders.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }

    cout << vector_res << endl;








    vector_res = "";
    //cout << m_rank << endl;
    /*
    for (size_t i = 0; i < domain_map_local.size(); ++i)
    {
      int item = domain_map_local.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }
    //cout << vector_res << endl;
    */

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
