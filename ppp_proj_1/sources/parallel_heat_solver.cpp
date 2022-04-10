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

#define FROM_DOWN_RANK_TAG_ROW_1 1
#define FROM_DOWN_RANK_TAG_ROW_2 2
#define FROM_UPPER_RANK_TAG_ROW_1 3
#define FROM_UPPER_RANK_TAG_ROW_2 4

#define FROM_LEFT_RANK_TAG_COL_1 5
#define FROM_LEFT_RANK_TAG_COL_2 6

#define FROM_RIGHT_RANK_TAG_COL_1 7
#define FROM_RIGHT_RANK_TAG_COL_2 8

#define AVERAGE_TEMP_TAG 9

#define FROM_DOWN_LEFT_RANK_TAG 10
#define FROM_DOWN_RIGHT_RANK_TAG 11
#define FROM_UPPER_LEFT_RANK_TAG 12
#define FROM_UPPER_RIGHT_RANK_TAG 13
//#define FROM_LEFT_RANK_TAG 3
//#define FROM_RIGHT_RANK_TAG 4

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

    else if (simulationProps.IsUseParallelIO() == true)
    {
      hid_t access_property_list = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(access_property_list, MPI_COMM_WORLD, MPI_INFO_NULL);
      m_fileHandle.Set(H5Fcreate(simulationProps.GetOutputFileName("par").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_property_list), H5Fclose);
      H5Pclose(access_property_list);
    }

}

ParallelHeatSolver::~ParallelHeatSolver()
{

}

// split vector via rows to list of vectors - each vector is one row
list<vector<int>> ParallelHeatSolver::SplitRows(int *input_arr, int local_tile_size, int local_tile_size_cols)
{

  list<vector<int>> output_list;

  bool start_new_vector = false;
  vector<int> row_items;
  int row_items_count = 0;
  for (size_t i = 0; i < local_tile_size; ++i)
  {
    int item = input_arr[i];
    //cout << "ITEM " << item << endl;

    if (start_new_vector == true)
    {
      row_items_count = 0;
      //cout << "Starting new vector " << endl;
      row_items.clear();
      start_new_vector = false;
    }

    row_items.push_back(item);
    row_items_count++;

    if (row_items_count == local_tile_size_cols)
    {
      start_new_vector = true;
      output_list.push_back(row_items);

    }
  }
    return output_list;

}

// split vector via rows to list of vectors - each vector is one row
list<vector<float>> ParallelHeatSolver::SplitRows(float *input_arr, int local_tile_size, int local_tile_size_cols)
{

  list<vector<float>> output_list;

  bool start_new_vector = false;
  vector<float> row_items;
  int row_items_count = 0;
  for (size_t i = 0; i < local_tile_size; ++i)
  {
    float item = input_arr[i];
    //cout << "ITEM " << item << endl;

    if (start_new_vector == true)
    {
      row_items_count = 0;
      //cout << "Starting new vector " << endl;
      row_items.clear();
      start_new_vector = false;
    }

    row_items.push_back(item);
    row_items_count++;

    if (row_items_count == local_tile_size_cols)
    {
      start_new_vector = true;
      output_list.push_back(row_items);

    }
  }
    return output_list;

}

vector<int> ParallelHeatSolver::EnlargeTile(list<vector<int>> input_list, int local_tile_size_cols)
{
  vector<int> output_vector_with_borders;
  list<vector<int>> enlarged_list;

  //cout << m_rank << endl;
  for (auto vect : input_list) {
      // Each element of the list is
      // a vector itself
      vector<int> currentVector = vect;

      cout << "[ ";// Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';
      currentVector.insert(currentVector.begin(), 0);
      currentVector.insert(currentVector.begin(), 0);
      currentVector.push_back(0);
      currentVector.push_back(0);
      enlarged_list.push_back(currentVector);
      cout << ']';
      cout << '\n';
  }

  vector<int> zeros_vector((local_tile_size_cols + 4), 0);
  enlarged_list.push_front(zeros_vector);
  enlarged_list.push_front(zeros_vector);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  for (auto vect : enlarged_list) {
      // Each element of the list is
      // a vector itself
      vector<int> currentVector = vect;
      for (size_t i = 0; i < currentVector.size(); ++i)
      {
        output_vector_with_borders.push_back(currentVector.at(i));
      }

      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';


      cout << ']';
      cout << '\n';
  }

  string vector_res = "";
  for (size_t i = 0; i < output_vector_with_borders.size(); ++i)
  {
    int item = output_vector_with_borders.at(i);
    vector_res.append(to_string(item));
    vector_res.append(", ");
  }
  //cout << "Print result map" << endl;
  //cout << vector_res << endl;
  return output_vector_with_borders;
}

vector<float> ParallelHeatSolver::EnlargeTile(list<vector<float>> input_list, int local_tile_size_cols)
{
  vector<float> output_vector_with_borders;
  list<vector<float>> enlarged_list;
  for (auto vect : input_list) {
      // Each element of the list is
      // a vector itself
      vector<float> currentVector = vect;

      cout << "[ ";// Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';
      currentVector.insert(currentVector.begin(), 0);
      currentVector.insert(currentVector.begin(), 0);
      currentVector.push_back(0);
      currentVector.push_back(0);
      enlarged_list.push_back(currentVector);
      cout << ']';
      cout << '\n';
  }

  vector<float> zeros_vector((local_tile_size_cols + 4), 0);
  enlarged_list.push_front(zeros_vector);
  enlarged_list.push_front(zeros_vector);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  for (auto vect : enlarged_list) {
      // Each element of the list is
      // a vector itself
      vector<float> currentVector = vect;
      for (size_t i = 0; i < currentVector.size(); ++i)
      {
        output_vector_with_borders.push_back(currentVector.at(i));
      }

      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';


      cout << ']';
      cout << '\n';
  }

  string vector_res = "";
  for (size_t i = 0; i < output_vector_with_borders.size(); ++i)
  {
    float item = output_vector_with_borders.at(i);
    vector_res.append(to_string(item));
    vector_res.append(", ");
  }
  //cout << "Print result map" << endl;
  //cout << vector_res << endl;
  return output_vector_with_borders;
}

vector<int> ParallelHeatSolver::Enlarge1DTile(int *input_arr, int local_tile_size)
{
  list<vector<int>> enlarged_list;
  vector<int> output_vector_with_borders;

  vector<int> zeros_vector((local_tile_size + 4), 0);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  vector<int> all_values_tile;
  all_values_tile.push_back(0);
  all_values_tile.push_back(0);

  for (size_t i = 0; i < local_tile_size; ++i)
  {
    int item = input_arr[i];
    all_values_tile.push_back(item);
  }

  all_values_tile.push_back(0);
  all_values_tile.push_back(0);

  enlarged_list.push_back(all_values_tile);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  for (auto vect : enlarged_list) {
      // Each element of the list is
      // a vector itself
      vector<int> currentVector = vect;
      for (size_t i = 0; i < currentVector.size(); ++i)
      {
        output_vector_with_borders.push_back(currentVector.at(i));
      }

      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';


      cout << ']';
      cout << '\n';
  }

  return output_vector_with_borders;


}

vector<float> ParallelHeatSolver::Enlarge1DTile(float *input_arr, int local_tile_size)
{
  list<vector<float>> enlarged_list;
  vector<float> output_vector_with_borders;

  vector<float> zeros_vector((local_tile_size + 4), 0.0);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  vector<float> all_values_tile;
  all_values_tile.push_back(0);
  all_values_tile.push_back(0);

  for (size_t i = 0; i < local_tile_size; ++i)
  {
    float item = input_arr[i];
    all_values_tile.push_back(item);
  }

  all_values_tile.push_back(0);
  all_values_tile.push_back(0);

  enlarged_list.push_back(all_values_tile);
  enlarged_list.push_back(zeros_vector);
  enlarged_list.push_back(zeros_vector);

  for (auto vect : enlarged_list) {
      // Each element of the list is
      // a vector itself
      vector<float> currentVector = vect;
      for (size_t i = 0; i < currentVector.size(); ++i)
      {
        output_vector_with_borders.push_back(currentVector.at(i));
      }

      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';


      cout << ']';
      cout << '\n';
  }

  return output_vector_with_borders;


}




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

void ParallelHeatSolver::print_array(float* arr, int width, int height)
{
    int i;
    for (i = 0; i< width * height;i++)
    {
        if((i != 0) && (i % width == 0))
            printf("\n");
        printf("%f ", arr[i]);
    }
    putchar('\n');
}

int ParallelHeatSolver::count_1D_index(int row, int length_of_row, int column)
{
  int index_1D = (row * length_of_row) + column;
  return index_1D;
}

float* ParallelHeatSolver::TrimTileWithoutBorders(float* arr, int enlarged_tile_size_rows, int enlarged_tile_size_cols, float* result)
{
  //float result[(enlarged_tile_size_rows - 4) * (enlarged_tile_size_cols - 4)];
  int result_index = 0;
  for(size_t i = 2; i < enlarged_tile_size_rows - 2; ++i)
  {
        for(size_t j = 2; j < enlarged_tile_size_cols - 2; ++j)
        {
            int index_1D = count_1D_index(i, enlarged_tile_size_cols, j);
            result[result_index] = arr[index_1D];
            result_index++;
        }
  }

  print_array(result, enlarged_tile_size_rows - 4, enlarged_tile_size_cols - 4);

  return result;
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

    int out_size_cols;
    int out_size_rows;
    int local_tile_size_cols; // size from x direction (cols for matrix meaning)
    int local_tile_size_rows; // size from y direction (rows in matrix meaning)
    int local_tile_size;

    int domain_length;
    //vector<float> init_temp(m_materialProperties.GetInitTemp().size());

    if (m_rank == 0)
    {
      m_simulationProperties.GetDecompGrid(out_size_cols, out_size_rows);
      //cout << "COLS " << out_size_cols << endl;
      //cout << "ROWS " << out_size_rows << endl;
      local_tile_size_cols = sqrt(m_materialProperties.GetInitTemp().size()) / out_size_cols;
      local_tile_size_rows = sqrt(m_materialProperties.GetInitTemp().size()) / out_size_rows;
      //cout << local_tile_size_rows << endl;
      local_tile_size = local_tile_size_cols * local_tile_size_rows;
      domain_length = m_materialProperties.GetInitTemp().size();

      /* recalculating params in case of 1D decomposition and getting correct composition of ranks
      e.g. 1D decomposition for 32 ranks and 16 x 16 grid, each rank has half of row of grid (16^2 / 32)

      +----+----+
      | 0  | 1  |
      +----+----+
      | 2  | 3  |
      +----+----+
      | 4  | 5  |
      +----+----+
      | 6  | 7  |
      +----+----+
      | 8  | 9  |
      +----+----+
      | 10 | 11 |
      +----+----+
      | 12 | 13 |
      +----+----+
      | 14 | 15 |
      +----+----+
      | 16 | 17 |
      +----+----+
      | 18 | 19 |
      +----+----+
      | 20 | 21 |
      +----+----+
      | 22 | 23 |
      +----+----+
      | 24 | 25 |
      +----+----+
      | 26 | 27 |
      +----+----+
      | 28 | 29 |
      +----+----+
      | 30 | 31 |
      +----+----+
*/
      if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D)
      {
        local_tile_size_cols = domain_length / out_size_cols;
        local_tile_size_rows = 1;
        out_size_cols = sqrt(m_materialProperties.GetInitTemp().size()) / local_tile_size_cols;
        out_size_rows = m_size / out_size_cols;

        cout << "Out size" << endl;
        cout << out_size_cols << endl;
        cout << out_size_rows << endl;
      }
    }

    MPI_Bcast(&out_size_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&out_size_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&domain_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //cout << local_tile_size << endl;

    vector<float> init_temp;
    vector<float> tmp_vector;
    vector<float> out_result;
    vector<int> domain_map;
    vector<float> domain_params;

    if (m_rank == 0)
    {

      for (size_t i = 0; i < m_materialProperties.GetInitTemp().size(); i++) {
        //cout << m_materialProperties.GetInitTemp().at(i) << endl;
        float value = m_materialProperties.GetInitTemp().at(i);
        init_temp.push_back(value);
      }

      //std::copy(m_materialProperties.GetInitTemp().begin(),
      //          m_materialProperties.GetInitTemp().end(), init_temp.begin());

    for (size_t i = 0; i < m_materialProperties.GetDomainParams().size(); i++) {
      float value = m_materialProperties.GetDomainParams().at(i);
      domain_params.push_back(value);
    }

    for (size_t i = 0; i < m_materialProperties.GetDomainMap().size(); i++) {
      int value = m_materialProperties.GetDomainMap().at(i);
      domain_map.push_back(value);
    }

    string vector_res = "";
    //cout << m_rank << endl;
    for (size_t i = 0; i < domain_map.size(); i++)
    {
      int item = domain_map.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }

    cout << "DOMAIN PARAMS ================" << endl;
    cout << vector_res << endl;

    }



    MPI_Barrier(MPI_COMM_WORLD);

    int m_num;

    MPI_Comm_size(MPI_COMM_WORLD, &m_num);
    int displacements[m_num];
    int counts[m_num];
    fill_n(counts, m_num, 1);

    int res = 0;
    bool skip_to_next_row = false;

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D)
    {
      for(int x = 0; x < m_size; x++)
      {
        if (x == 0)
        {
          res = 0;
        }
        else
        {
          //cout << "Local tile size x " << local_tile_size_cols << endl;
          res = res + local_tile_size_cols;
        }
        displacements[x] = res;
        //cout << "Result is " << res << endl;
      }
    }
    else {
    for(int x = 0; x < out_size_rows; x++){
        for(int y = 0; y < out_size_cols; y++){
            cout << x << ", " << y<< endl;
            cout << "Process on index" << x*out_size_cols+y <<  endl;
            //int res =  x * sqrt(domain_length) * local_tile_size_cols + y * local_tile_size_rows;
            if (x == 0 & y == 0)
            {
              res = 0;
            }
            else if (skip_to_next_row == true)
            {
              res = local_tile_size_rows * sqrt(domain_length) * x;
              skip_to_next_row = false;
            }
            else if (skip_to_next_row == false)
            {
              res =  res + local_tile_size_cols;
            }

            displacements[x*out_size_cols+y] = res;
            cout << "Result is " << res << endl;


        }
        skip_to_next_row = true;
        //res = sqrt(domain_length) * local_tile_size_rows;
    }
  }



  int enlarged_tile_size = (local_tile_size_cols + 4) * (local_tile_size_rows + 4);
  int enlarged_tile_size_cols = local_tile_size_cols + 4;
  int enlarged_tile_size_rows = local_tile_size_rows + 4;

  cout << "Enlarged sizes cols " << enlarged_tile_size_cols << endl;
  cout << "Enlarged sizes rows " << enlarged_tile_size_rows << endl;

  float *init_temp_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  int *domain_map_local = (int*)malloc(enlarged_tile_size * sizeof(int));
  float *domain_params_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  assert(domain_map_local != NULL);

  // float *init_temp_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  MPI_Win win;
  float* win_memory;
  //MPI_Win_create(init_temp_local, enlarged_tile_size * sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  MPI_Win_allocate(2 * enlarged_tile_size * sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_memory, &win);

  MPI_Datatype tile_t, resized_tile_t;
  MPI_Type_vector(local_tile_size_rows, local_tile_size_cols, sqrt(domain_length), MPI_FLOAT, &tile_t);
  MPI_Type_create_resized(tile_t, 0, sizeof(float), &resized_tile_t);
  MPI_Type_commit(&tile_t);
  MPI_Type_commit(&resized_tile_t);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(&(init_temp[0]), counts, displacements, resized_tile_t, &(init_temp_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_map[0]), counts, displacements, resized_tile_t, &(domain_map_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_params[0]), counts, displacements, resized_tile_t, &(domain_params_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < m_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (i == m_rank) {
            printf("\nRank: %d\n", m_rank);
            print_array(domain_params_local, local_tile_size_cols, local_tile_size_rows);
        }
    }

    vector<int> middle_ranks;

    // FIND RANKS THAT HOLDS MIDDLE COLUMN OF BOARD AND CREATE NEW COMMUNICATOR FOR THEM
    int middle_item_col_id = (sqrt(domain_length) / 2); // middle column of values
    int middle_item_tile_col_id = middle_item_col_id % local_tile_size_cols;
    int middle_col_id = middle_item_col_id / local_tile_size_cols;

    for (int i = 0; i < m_num; i++) {
      int column_rank = i % out_size_cols;
      if (middle_col_id == column_rank || i == 0)
      {
        cout << "Middle process rank is " << i << endl;
        middle_ranks.push_back(i);
      }
    }


    string vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < middle_ranks.size(); i++)
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

    vector<float> init_temp_local_with_borders;
    vector<float> domain_params_local_with_borders;
    vector<int> domain_map_local_with_borders;

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D)
    {
      domain_map_local_with_borders = Enlarge1DTile(domain_map_local, local_tile_size_cols);
      domain_params_local_with_borders = Enlarge1DTile(domain_params_local, local_tile_size_cols);
      init_temp_local_with_borders = Enlarge1DTile(init_temp_local, local_tile_size_cols);
    }
    else {
    list<vector<int>> domain_map_local_list = SplitRows(domain_map_local, local_tile_size, local_tile_size_rows);
    domain_map_local_with_borders = EnlargeTile(domain_map_local_list, local_tile_size_cols);
    list<vector<float>> domain_params_local_list = SplitRows(domain_params_local, local_tile_size, local_tile_size_cols);
    domain_params_local_with_borders = EnlargeTile(domain_params_local_list, local_tile_size_cols);
    list<vector<float>> init_temp_local_list = SplitRows(init_temp_local, local_tile_size, local_tile_size_cols);
    init_temp_local_with_borders = EnlargeTile(init_temp_local_list, local_tile_size_cols);
    }

    MPI_Datatype tile_row_t, tile_col_t, resized_tile_col_t;
    MPI_Type_contiguous(enlarged_tile_size_cols, MPI_FLOAT, &tile_row_t);
    MPI_Type_commit(&tile_row_t);

    MPI_Type_vector(enlarged_tile_size_rows, 2, enlarged_tile_size_cols, MPI_FLOAT, &tile_col_t);
    //MPI_Type_create_resized(tile_col_t, 1, sizeof(float), &resized_tile_col_t);
    MPI_Type_commit(&tile_col_t);
    //MPI_Type_commit(&resized_tile_col_t);

    int row_id = m_rank / out_size_cols;
    int col_id = m_rank % out_size_cols;

    int down_rank = m_rank + out_size_cols;
    int upper_rank = m_rank - out_size_cols;
    int left_rank = m_rank - 1;
    int right_rank = m_rank + 1;

    int down_left_rank = m_rank + out_size_cols - 1;
    int down_right_rank = m_rank + out_size_cols + 1;
    int upper_left_rank = m_rank - out_size_cols - 1;
    int upper_right_rank = m_rank - out_size_cols + 1;

    vector<float> init_temp_local_with_borders_recieved;
    vector<float> init_temp_recieved(enlarged_tile_size);
    for (size_t i = 0; i < init_temp_local_with_borders.size(); i++) {
      float value = init_temp_local_with_borders.at(i);
      //cout << "Value ===" << value << endl;
      init_temp_local_with_borders_recieved.push_back(value);
    }

    // each enlarge tile can be divided into 2 x 2 tiles
    MPI_Datatype enlarge_tile_part_t, resized_enlarge_tile_part_t;
    MPI_Type_vector(2, 2, enlarged_tile_size_cols, MPI_FLOAT, &enlarge_tile_part_t);
    MPI_Type_create_resized(enlarge_tile_part_t, 0, sizeof(float), &resized_enlarge_tile_part_t);
    MPI_Type_commit(&enlarge_tile_part_t);
    MPI_Type_commit(&resized_enlarge_tile_part_t);

    MPI_Barrier(MPI_COMM_WORLD);

    int total_request_count = 0;


    if (row_id == 0 || row_id == out_size_rows - 1)
    {

      total_request_count += 2;
    }
    else
    {
      total_request_count += 4;
    }

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_2D || (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && out_size_cols != 1)) {
    if (col_id == 0 || col_id == out_size_cols - 1)
    {
      total_request_count += 2;
    }
    else
    {
      total_request_count += 4;
    }
  }


  int edges_request_count = 0;
  if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_2D || (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && out_size_cols != 1)) {
  if ((row_id == 0 || row_id == out_size_rows - 1) && (col_id == 0 || col_id == out_size_cols - 1))
  {
    edges_request_count += 2;
  }
  else if (col_id == 0 || col_id == out_size_cols - 1 || row_id == 0 || row_id == out_size_rows - 1)
  {
    edges_request_count += 4;
  }
  else
  {
    edges_request_count += 8;
  }
}

    cout << "Upcoming requests for rank " << m_rank << ":   " << total_request_count << endl;


    MPI_Request requests[total_request_count];
    int num_requests = 0;
    MPI_Barrier(MPI_COMM_WORLD);

    if (row_id != out_size_rows - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&init_temp_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
    }


    MPI_Waitall(total_request_count, requests, NULL);
    cout << "HHHH fff " << endl;

    MPI_Request requests_init_temp_recieved[edges_request_count];
    num_requests = 0;

    vector<float> init_temp_local_with_borders_recieved_with_edges;
    for (size_t i = 0; i < init_temp_local_with_borders_recieved.size(); i++) {
      float value = init_temp_local_with_borders_recieved.at(i);
      init_temp_local_with_borders_recieved_with_edges.push_back(value);
    }

    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, 0);
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }

    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, 2);
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, 2);
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);

    }


    MPI_Waitall(edges_request_count, requests_init_temp_recieved, NULL);

    vector_res = "";
    cout << "Process with received values " << m_rank << endl;

  /*  for (size_t i = 0; i < init_temp_local_with_borders_recieved.size(); ++i)
    {
      float item = init_temp_local_with_borders_recieved.at(i);
      vector_res.append(to_string(item));
      vector_res.append(";; ");
    }

    cout << vector_res << endl;*/

    vector<int> domain_map_local_with_borders_recieved;
    for (size_t i = 0; i < domain_map_local_with_borders.size(); i++) {
      int value = domain_map_local_with_borders.at(i);
      domain_map_local_with_borders_recieved.push_back(value);
    }

    MPI_Request requests_domain_map[total_request_count];
    num_requests = 0;
    MPI_Barrier(MPI_COMM_WORLD);

    if (row_id != out_size_rows - 1)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }


    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_map_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    MPI_Waitall(total_request_count, requests_domain_map, NULL);

    vector<int> domain_map_local_with_borders_recieved_with_edges;

    for (size_t i = 0; i < domain_map_local_with_borders_recieved.size(); i++) {
      int value = domain_map_local_with_borders_recieved.at(i);
      domain_map_local_with_borders_recieved_with_edges.push_back(value);
    }


    MPI_Request requests_domain_map_recieved[edges_request_count];
    num_requests = 0;

    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }

    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&domain_map_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, 2);
      MPI_Isend(&domain_map_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);
    }

    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&domain_map_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, 2);
      MPI_Isend(&domain_map_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);
    }


    MPI_Waitall(edges_request_count, requests_domain_map_recieved, NULL);


    vector<float> domain_params_local_with_borders_recieved;

    for (size_t i = 0; i < domain_params_local_with_borders.size(); i++) {
      float value = domain_params_local_with_borders.at(i);
      domain_params_local_with_borders_recieved.push_back(value);
    }


    MPI_Request requests_domain_params[total_request_count];
    num_requests = 0;
    MPI_Barrier(MPI_COMM_WORLD);

    if (row_id != out_size_rows - 1)
    {
        // store to last two rows
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != 0)
    {
        // store to last two rows
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }


    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_params_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
       /// -------->
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    MPI_Waitall(total_request_count, requests_domain_params, NULL);

    vector<float> domain_params_local_with_borders_recieved_with_edges;

    for (size_t i = 0; i < domain_params_local_with_borders_recieved.size(); i++) {
      float value = domain_params_local_with_borders_recieved.at(i);
      domain_params_local_with_borders_recieved_with_edges.push_back(value);
    }




    MPI_Request requests_domain_params_recieved[edges_request_count];
    num_requests = 0;

    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }

    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&domain_params_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, 2);
      MPI_Isend(&domain_params_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);
    }

    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      MPI_Isend(&domain_params_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, 2);
      MPI_Isend(&domain_params_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);
    }


    MPI_Waitall(edges_request_count, requests_domain_params, NULL);

    for (int i = 0; i < m_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (i == m_rank) {
            printf("\nRANK: %d\n", m_rank);
            print_array(&domain_params_local_with_borders_recieved_with_edges[0], enlarged_tile_size_cols, enlarged_tile_size_rows);
        }
    }

    if (m_rank == 6) {
    //cout << "Process with received values params" << m_rank << endl;

    for (size_t i = 0; i < domain_params_local_with_borders_recieved.size(); i++)
    {
      float item = domain_params_local_with_borders_recieved.at(i);
      vector_res.append(to_string(item));
      vector_res.append(";; ");
    }

    cout << vector_res << endl;}

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

    MPI_Win_fence(0, win);


    if (m_rank == 0)
    {
      double startTime = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    float middleColAvgTemp = 0.0f;
    float *workTempArrays[] = {init_temp_local_with_borders_recieved_with_edges.data(), init_temp_local_with_borders_recieved_with_edges.data()};

    int offset_cols_begin = 0;
    int offset_rows_begin = 0;
    int offset_cols_end = 0;
    int offset_rows_end = 0;

    if (row_id == 0)
    {
      offset_rows_begin += 2;
    }
    if (row_id == out_size_rows - 1)
    {
      offset_rows_end += 2;
    }
    if (col_id == 0)
    {
      offset_cols_begin += 2;
    }
    if (col_id == out_size_cols - 1)
    {
      offset_cols_end += 2;
    }

    if (m_rank == 2)
    {
      cout << "Offsets ==" << endl;
      cout << offset_cols_begin << endl;
      cout << "Offsets ==" << endl;
      cout << offset_rows_begin << endl;
      cout << "Offsets ==" << endl;
      cout << offset_cols_end << endl;
      cout << "Offsets ==" << endl;
      cout << offset_rows_end << endl;
    }


    /*
    MPI_Type_vector(local_tile_size_rows, local_tile_size_cols, sqrt(domain_length), MPI_FLOAT, &tile_t);
    MPI_Type_create_resized(tile_t, 0, sizeof(float), &resized_tile_t);
    MPI_Type_commit(&tile_t);
    MPI_Type_commit(&resized_tile_t);
    */
    MPI_Datatype worker_tile_t;
    int tile[2];
    int worker_size[2];
    int worker_start[2];
    worker_size[0] = enlarged_tile_size_rows;
    worker_size[1] = enlarged_tile_size_cols;
    tile[0] = local_tile_size_rows;
    tile[1] = local_tile_size_cols;
    worker_start[0] = 2;
    worker_start[1] = 2;
    float result_domain[domain_length];

    MPI_Type_create_subarray(2, worker_size, tile, worker_start, MPI_ORDER_C, MPI_FLOAT, &worker_tile_t);
    MPI_Type_commit(&worker_tile_t);

    MPI_Datatype farmer_matrix_t;
    MPI_Datatype resized_farmer_matrix_t;
    int farmer_size[2];
    int farmer_start[2];
    farmer_size[0] = sqrt(domain_length);
    farmer_size[1] = sqrt(domain_length);
    farmer_start[0] = 0;
    farmer_start[1] = 0;

    MPI_Type_create_subarray(2, farmer_size, tile, farmer_start, MPI_ORDER_C, MPI_FLOAT, &farmer_matrix_t);
    MPI_Type_create_resized(farmer_matrix_t, 0, sizeof(float), &resized_farmer_matrix_t);
    MPI_Type_commit(&farmer_matrix_t);
    MPI_Type_commit(&resized_farmer_matrix_t);

    // MPI_Scatterv(&(init_temp[0]), counts, displacements, resized_tile_t, &(init_temp_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &result_domain[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);

    const char* fileName    = "File1.h5";
      const char* datasetName = "Dataset-1";

      if (m_rank == 0)
      {
        // 1. Declare an HDF5 file.
        hid_t file = H5I_INVALID_HID;

        // 2. Create a file with write permission. Use such a flag that overrides existing file.
        //    The list of flags is in the header file called H5Fpublic.h
        cout << "Creating file" << endl;
        file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


        // 3. Create file and memory spaces. We will only write a single value.
        const hsize_t rank = 1;
        const hsize_t size = 1;

        hid_t filespace = H5Screate_simple(rank, &size, nullptr);
        hid_t memspace  = H5Screate_simple(rank, &size, nullptr);

        // 4. Create a dataset of a size [1] and int datatype.
        //    The list of predefined datatypes can be found in H5Tpublic.h
        hid_t dataset   = H5Dcreate(file, datasetName, H5T_NATIVE_INT, filespace,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        //5. Write value into the dataset.
        cout << "Write scalar value" << endl;
        const int value  = 128;
        H5Dwrite(dataset, H5T_NATIVE_INT, filespace, filespace, H5P_DEFAULT, &value);

        // 6. Close dataset.
        H5Dclose(dataset);

        // 7. Close file

        H5Fclose(file);
      }

      hid_t filespace;
      hid_t memspace;
      hid_t dataset;
      hid_t xferPList;

    if (!m_simulationProperties.GetOutputFileName().empty() && m_simulationProperties.IsUseParallelIO())
    {
        hsize_t datasetSize[] = {hsize_t(sqrt(domain_length)), hsize_t(sqrt(domain_length))};
        hsize_t memSize[]     = {hsize_t(local_tile_size_rows), hsize_t(local_tile_size_cols)};

        filespace = H5Screate_simple(2, datasetSize, nullptr);
        memspace  = H5Screate_simple(2, memSize,     nullptr);

        xferPList = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xferPList, H5FD_MPIO_COLLECTIVE);

        dataset = H5Dcreate(m_fileHandle, "snapshot_dataset", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        hsize_t slabStart[] = {hsize_t(row_id * local_tile_size_rows), hsize_t(col_id * local_tile_size_cols)};
        hsize_t slabSize[]  = {hsize_t(local_tile_size_rows), hsize_t(local_tile_size_cols)};

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, slabStart, nullptr, slabSize, nullptr);
    }

    if (m_rank == 0)
    {
      cout << "Printing whole domain" << endl;
      print_array(result_domain, sqrt(domain_length), sqrt(domain_length));

      cout << "Printing my tile without borders" << endl;
      float *result_tile;
      float result[local_tile_size];
      result_tile = TrimTileWithoutBorders(&domain_params_local_with_borders_recieved_with_edges[0], enlarged_tile_size_rows, enlarged_tile_size_cols, result);
      cout << "HHHH " << endl;
      cout << result_tile << endl;
      print_array(result_tile, local_tile_size_rows, local_tile_size_cols);
    }

    float whole_domain_temp_snapshot[domain_length]; // snapshot of whole domain result for outputing into file by rank 0 sequentially
    MPI_Win_fence(0, win);

    for (size_t iter = 0; iter < m_simulationProperties.GetNumIterations(); ++iter)
    {

      for(size_t i = 2 + offset_rows_begin; i < enlarged_tile_size_rows - offset_rows_end - 2; ++i)
      {
            for(size_t j = 2 + offset_cols_begin; j < enlarged_tile_size_cols - offset_cols_end - 2; ++j)
            {
                ComputePoint(workTempArrays[1], workTempArrays[0],
                        domain_params_local_with_borders_recieved_with_edges.data(),
                        domain_map_local_with_borders_recieved_with_edges.data(),
                        i, j,
                        enlarged_tile_size_cols,
                        m_simulationProperties.GetAirFlowRate(),
                        m_materialProperties.GetCoolerTemp());
            }
      }

      MPI_Barrier(MPI_COMM_WORLD);

      if (m_simulationProperties.IsRunParallelP2P()) {
      MPI_Request requests_simulation[total_request_count];
      int num_requests_simulation = 0;


      if (row_id != out_size_rows - 1)
      {
          // store to last two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != 0)
      {
          // store to last two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != out_size_cols - 1)
      {
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != 0)
      {
          MPI_Irecv(&workTempArrays[1][0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != 0)
      {
          // send two rows up from down rank
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != out_size_rows - 1)
      {
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }


      if (col_id != 0)
      {
           // <--------
          MPI_Isend(&workTempArrays[0][2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != out_size_cols - 1)
      {
        /// -------->
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      MPI_Waitall(total_request_count, requests_simulation, NULL);
    }

      else
      {
        cout << "Using RMA " << endl;
        int offset = enlarged_tile_size * (( iter + 1) % 2);
        cout << "OFFSET " << offset << endl;
        MPI_Win_fence(0, win);

        // The execution of a put operation is similar to the execution of a send by the origin process and a matching receive by the target process. The obvious difference is that all arguments are provided by one call --- the call executed by the origin process.

        if (row_id != 0)
        {
            // send two rows up from down rank
            MPI_Put(&workTempArrays[0][enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, (enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)) + offset, 2, tile_row_t, win);
        }

        if (row_id != out_size_rows - 1)
        {
            MPI_Put(&workTempArrays[0][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, (enlarged_tile_size_cols * 0) + offset, 2, tile_row_t, win);
        }
        if (col_id != 0)
        {
             // <--------
            MPI_Put(&workTempArrays[0][2], 1, tile_col_t, left_rank, enlarged_tile_size_cols - 1 - 1 + offset, 1, tile_col_t, win);
        }
        if (col_id != out_size_cols - 1)
        {
          /// -------->
            MPI_Put(&workTempArrays[0][enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, 0 + offset, 1, tile_col_t, win);
        }

        MPI_Win_fence(0, win);

      }
      MPI_Barrier(MPI_COMM_WORLD);


      float final_iteration_temp = 0.0f;
      if (count(middle_ranks.begin(), middle_ranks.end(), m_rank))
      {
        cout << "Middle " << m_rank << endl;

        if (m_rank == 0)
        {
          middleColAvgTemp = 0.0f;
        }
        else
        {
          middleColAvgTemp = ComputeMiddleColAvgTemp(workTempArrays[0], enlarged_tile_size_rows, enlarged_tile_size_cols, out_size_rows, out_size_cols, middle_item_tile_col_id, m_rank);
          //cout << "Temperature " << middleColAvgTemp << endl;
        }

        vector<float> all_average_temp(out_size_rows + 1);
        MPI_Gather(&middleColAvgTemp, 1, MPI_FLOAT, &all_average_temp[0], 1, MPI_FLOAT, 0, MPI_COMM_MIDDLE_COLUMN);

        if (m_rank == 0)
        {
          float all_average_temp_sum = 0.0f;

          for (size_t i = 0; i < all_average_temp.size(); ++i)
          {
            if (i == 0)
            {
              continue;
            }
            float item = all_average_temp.at(i);
            //cout << "VAAAL " << item << endl;
            all_average_temp_sum += item;
          }

          //cout << "ALLLL SUM " << all_average_temp_sum << endl;
          final_iteration_temp = all_average_temp_sum / (all_average_temp.size() - 1);
        }

      }

      MPI_Barrier(MPI_COMM_WORLD);



      if (!m_simulationProperties.GetOutputFileName().empty()) {
            if (!m_simulationProperties.IsUseParallelIO() && ((iter % m_simulationProperties.GetDiskWriteIntensity()) == 0))
            {
              MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &whole_domain_temp_snapshot[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);
              if (m_rank == 0)
              {
                StoreDataIntoFile(m_fileHandle, iter, &whole_domain_temp_snapshot[0]);
              }
            }
            else if (m_simulationProperties.IsUseParallelIO())
            {
              cout << "Iter " << iter << endl;
              if (iter ==  m_simulationProperties.GetNumIterations() - 1) {
              cout << "JJJJ" << endl;
              float *tile_snapshot;
              float tile_result[local_tile_size];
              tile_snapshot = TrimTileWithoutBorders(&workTempArrays[0][0], enlarged_tile_size_rows, enlarged_tile_size_cols, tile_result);
              print_array(tile_snapshot, local_tile_size_rows, local_tile_size_cols);
              H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, xferPList, &tile_snapshot[0]);
            }
            }
      }

      swap(workTempArrays[0], workTempArrays[1]);

      if (m_rank == 0)
      {
        PrintProgressReport(iter, final_iteration_temp);
      }





    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (!m_simulationProperties.GetOutputFileName().empty() && m_simulationProperties.IsUseParallelIO())
    {
      H5Dclose(dataset);
    }


}

float ParallelHeatSolver::ComputeMiddleColAvgTemp(const float *data, int enlarged_tile_size_rows, int enlarged_tile_size_cols, int out_size_rows, int out_size_cols, int middle_item_tile_col_id, int m_rank) const
{
    float middleColAvgTemp = 0.0f;
    float result = 0.0f;


      // index of column which will have influence on middle temperature
      int middle_column_tile_id = 2 + middle_item_tile_col_id;
      cout << "Middle column tile id " << middle_column_tile_id <<  endl;
      cout << "Enlarge rows " << enlarged_tile_size_rows << endl;
      cout << "Enlarge cols " << enlarged_tile_size_cols << endl;

      for(size_t i = 2; i < enlarged_tile_size_rows - 2; ++i)
      {
            for(size_t j = 2; j < enlarged_tile_size_cols - 2; ++j)
            {
                int index_1D = i * enlarged_tile_size_cols + j;

                int row_id = index_1D / enlarged_tile_size_cols;
                int col_id = index_1D % enlarged_tile_size_cols;

                if (m_rank == 2)
                {
                  float value = data[index_1D];
                  cout << "TEMP " << value << endl;
                }

                if ((middle_column_tile_id == col_id && row_id != 0) && (middle_column_tile_id == col_id && row_id != 1) && (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 2) && (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 1))
                {
                  cout << "IMDEX" << index_1D << endl;
                  cout << data[index_1D] << endl;
                  float value = data[index_1D];
                  middleColAvgTemp += value;
                }
                /*
                cout << "IMDEX" << index_1D << endl;
                cout << "vaaaal" << data[index_1D] << endl;
                middleColAvgTemp += data[index_1D];
                */
            }
      }


    result = middleColAvgTemp / (enlarged_tile_size_rows - 4);
    return result;
}
