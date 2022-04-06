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
  vector<int> output_vector_with_borders;

  output_vector_with_borders.push_back(0);
  output_vector_with_borders.push_back(0);

  for (size_t i = 0; i < local_tile_size; ++i)
  {
    int item = input_arr[i];
    output_vector_with_borders.push_back(item);
  }

  output_vector_with_borders.push_back(0);
  output_vector_with_borders.push_back(0);

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

vector<float> ParallelHeatSolver::Enlarge1DTile(float *input_arr, int local_tile_size)
{
  vector<float> output_vector_with_borders;

  output_vector_with_borders.push_back(0);
  output_vector_with_borders.push_back(0);

  for (size_t i = 0; i < local_tile_size; ++i)
  {
    float item = input_arr[i];
    output_vector_with_borders.push_back(item);
  }

  output_vector_with_borders.push_back(0);
  output_vector_with_borders.push_back(0);

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
      int init_temp_size = sqrt(m_materialProperties.GetInitTemp().size());
      //cout << "MAT " << init_temp_size << endl;
      local_tile_size_cols = init_temp_size / out_size_cols;
      local_tile_size_rows = sqrt(m_materialProperties.GetInitTemp().size()) / out_size_rows;
      //cout << "L" << local_tile_size_cols << endl;
      //cout << local_tile_size_rows << endl;
      local_tile_size = local_tile_size_cols * local_tile_size_rows;
      domain_length = m_materialProperties.GetInitTemp().size();

      if (out_size_cols != 1 && out_size_rows == 1)
      {
        local_tile_size_cols = domain_length / out_size_cols;
        local_tile_size_rows = 1;
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



    }


    string vector_res = "";
    //cout << m_rank << endl;
    for (size_t i = 0; i < domain_params.size(); i++)
    {
      float item = domain_params.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }

    cout << "DOMAIN PARAMS ================" << endl;
    cout << vector_res << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    int m_num;

    MPI_Comm_size(MPI_COMM_WORLD, &m_num);
    int displacements[m_num];
    int counts[m_num];
    fill_n(counts, m_num, 1);

    int res = 0;
    bool skip_to_next_row = false;
    //cout << "LOCAL TILE SIZE X " << local_tile_size_cols << endl;
    if (out_size_cols != 1 && out_size_rows == 1)
    {
      cout << "HHHHHH" << endl;
      for(int x = 0; x < out_size_cols; x++)
      {
        if (x == 0)
        {
          res = 0;
        }
        else
        {
          cout << "Local tile size x " << local_tile_size_cols << endl;
          res = res + local_tile_size_cols;
        }
        displacements[x] = res;
        cout << "Result is " << res << endl;
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

  float *init_temp_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  int *domain_map_local = (int*)malloc(enlarged_tile_size * sizeof(int));
  float *domain_params_local = (float*)malloc(enlarged_tile_size * sizeof(float));
  assert(domain_map_local != NULL);

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

    if (out_size_cols != 1 && out_size_rows == 1)
    {
      cout << "TODO" << endl;
    }
    else {
      //int middle_column_id = (sqrt(domain_length) / 2) - 1;


      //int middle_col_id = middle_column_id / local_tile_size_cols;



    }

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


    vector_res = "";
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
    if (out_size_cols != 1 && out_size_rows == 1)
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

    MPI_Type_vector(enlarged_tile_size_rows, 1, enlarged_tile_size_cols, MPI_FLOAT, &tile_col_t);
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

    if (out_size_cols != 1 && out_size_rows == 1)
    {
      if (col_id == 0 || col_id == out_size_cols - 1)
      {
        total_request_count += 4;
      }
      else
      {
        total_request_count += 8;
      }
    }
    else {
    if (row_id == 0 || row_id == out_size_rows - 1)
    {

      total_request_count += 4;
    }
    else
    {
      total_request_count += 8;
    }


    if (col_id == 0 || col_id == out_size_cols - 1)
    {
      total_request_count += 4;
    }
    else
    {
      total_request_count += 8;
    }
  }

  int edges_request_count = 0;

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

    cout << "Upcoming requests for rank " << m_rank << ":   " << total_request_count << endl;


    MPI_Request requests[total_request_count];
    int num_requests = 0;
    MPI_Barrier(MPI_COMM_WORLD);

    if (row_id != out_size_rows - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 1], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }


    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Irecv(&init_temp_local_with_borders_recieved[1], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * 3], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 2)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }


    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&init_temp_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Isend(&init_temp_local_with_borders[3], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests[num_requests++]);
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols - 1 - 2], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests[num_requests++]);
    }


    MPI_Waitall(total_request_count, requests, NULL);

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
      //cout << "1D index is " << index_1D << endl;
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      //cout << "1D index is " << index_1D << endl;
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      //cout << "1D index is " << index_1D << endl;
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      //cout << "1D index is " << index_1D << endl;
      MPI_Irecv(&init_temp_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }

    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      //cout << "1D index is " << index_1D << endl;
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      int index_1D = count_1D_index(2, enlarged_tile_size_cols, 2);
      //cout << "1D index is " << index_1D << endl;
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, enlarged_tile_size_cols - 4);
      //cout << "1D index is " << index_1D << endl;
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);
    }
    if ((row_id != out_size_rows - 1 && col_id != 0) && (row_id != out_size_rows - 1) && (col_id != 0))
    {
      int index_1D = count_1D_index(enlarged_tile_size_rows - 4, enlarged_tile_size_cols, 2);
      cout << "1D index is " << index_1D << endl;
      MPI_Isend(&init_temp_local_with_borders_recieved[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_init_temp_recieved[num_requests++]);

    }


    MPI_Waitall(edges_request_count, requests_init_temp_recieved, NULL);






    //MPI_Barrier(MPI_COMM_WORLD);

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
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 1], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Irecv(&domain_map_local_with_borders_recieved[1], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * 3], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 2)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }


    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_map_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Isend(&domain_map_local_with_borders[3], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols - 1 - 2], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
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
      //cout << "Rank " << m_rank << endl;
      //cout << down_left_rank << endl;
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << down_right_rank << endl;
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << upper_left_rank << endl;
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << upper_right_rank << endl;
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_map_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_right_rank, FROM_UPPER_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_map_recieved[num_requests++]);

    }

    //if (row_id != 0 || col_id != out_size_cols - 1) {
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
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != 0)
    {
        // store to last two rows
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 1], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }


    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (col_id != 0)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Irecv(&domain_params_local_with_borders_recieved[1], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * 3], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 2)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }


    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_params_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Isend(&domain_params_local_with_borders[3], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (col_id != out_size_cols - 1)
    {
       /// -------->
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols - 1 - 2], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
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
      //cout << "Rank " << m_rank << endl;
      //cout << down_left_rank << endl;
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_left_rank, FROM_DOWN_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != out_size_rows - 1 && col_id != out_size_cols - 1) && (row_id != out_size_rows - 1) && (col_id != out_size_cols - 1))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << down_right_rank << endl;
      int index_1D = count_1D_index(enlarged_tile_size_rows - 2, enlarged_tile_size_cols, enlarged_tile_size_cols - 2);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, down_right_rank, FROM_DOWN_RIGHT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != 0) && (row_id != 0) && (col_id != 0))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << upper_left_rank << endl;
      int index_1D = count_1D_index(0, enlarged_tile_size_cols, 0);
      MPI_Irecv(&domain_params_local_with_borders_recieved_with_edges[index_1D], 1, enlarge_tile_part_t, upper_left_rank, FROM_UPPER_LEFT_RANK_TAG, MPI_COMM_WORLD, &requests_domain_params_recieved[num_requests++]);

    }
    if ((row_id != 0 && col_id != out_size_cols - 1) && (row_id != 0) && (col_id != out_size_cols - 1))
    {
      //cout << "Rank " << m_rank << endl;
      //cout << upper_right_rank << endl;
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
            printf("\nRank: %d\n", m_rank);
            print_array(&domain_params_local_with_borders_recieved_with_edges[0], enlarged_tile_size_cols, enlarged_tile_size_rows);
        }
    }

    //cout << "LLLLLLLLLLLLLL " << endl;

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

    if (m_rank == 0)
    {
      double startTime = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    float middleColAvgTemp = 0.0f;
    float *workTempArrays[] = {init_temp_local_with_borders_recieved_with_edges.data(), init_temp_local_with_borders_recieved_with_edges.data()};
    //print_array(workTempArrays[0],8,8);


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

    if (m_rank == 4)
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

      MPI_Request requests_simulation[total_request_count];
      int num_requests_simulation = 0;


      if (row_id != out_size_rows - 1)
      {
          // store to last two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != 0)
      {
          // store to last two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * 1], 1, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != out_size_cols - 1)
      {
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != 0)
      {
          MPI_Irecv(&workTempArrays[1][0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Irecv(&workTempArrays[1][1], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != 0)
      {
          // send two rows up from down rank
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 3], 1, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (row_id != out_size_rows - 1)
      {
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 2)], 1, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }


      if (col_id != 0)
      {
           // <--------
          MPI_Isend(&workTempArrays[0][2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Isend(&workTempArrays[0][3], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      if (col_id != out_size_cols - 1)
      {
        /// -------->
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_1, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols - 1 - 2], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL_2, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
      }

      MPI_Waitall(total_request_count, requests_simulation, NULL);
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
          middleColAvgTemp = ComputeMiddleColAvgTemp(workTempArrays[0], enlarged_tile_size_rows, enlarged_tile_size_cols, out_size_rows, out_size_cols, middle_item_tile_col_id);
          //cout << "Temperature " << middleColAvgTemp << endl;
        }

        vector<float> all_average_temp(out_size_rows + 1);
        MPI_Gather(&middleColAvgTemp, 1, MPI_FLOAT, &all_average_temp[0], 1, MPI_FLOAT, 0, MPI_COMM_MIDDLE_COLUMN);

        if (m_rank == 0)
        {
          float all_average_temp_sum = 0;

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

      if (m_rank == 0)
      {
        PrintProgressReport(iter, final_iteration_temp);
      }

      swap(workTempArrays[0], workTempArrays[1]);






    }
    MPI_Barrier(MPI_COMM_WORLD);







}

float ParallelHeatSolver::ComputeMiddleColAvgTemp(const float *data, int enlarged_tile_size_rows, int enlarged_tile_size_cols, int out_size_rows, int out_size_cols, int middle_item_tile_col_id) const
{
    float middleColAvgTemp = 0.0f;
    float result = 0.0f;

    if (out_size_cols != 1 && out_size_rows == 1)
    {

    }
    else
    {
      // index of column which will have influence on middle temperature
      int middle_column_tile_id = 2 + middle_item_tile_col_id;
      cout << "Middle column tile id " << middle_column_tile_id <<  endl;
      cout << "Enlarge rows " << enlarged_tile_size_rows << endl;
      cout << "Enlarge cols " << enlarged_tile_size_cols << endl;
      //int row_id = m_rank / out_size_cols;
      //int col_id = m_rank % out_size_cols;


      for(size_t i = 2; i < enlarged_tile_size_rows - 2; ++i)
      {
            for(size_t j = 2; j < enlarged_tile_size_cols - 2; ++j)
            {
                int index_1D = i * enlarged_tile_size_cols + j;

                int row_id = index_1D / enlarged_tile_size_cols;
                int col_id = index_1D % enlarged_tile_size_cols;

                if ((middle_column_tile_id == col_id && row_id != 0) || (middle_column_tile_id == col_id && row_id != 1) || (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 2) || (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 1))
                {
                  cout << "IMDEX" << index_1D << endl;
                  middleColAvgTemp += data[index_1D];
                }
                /*
                cout << "IMDEX" << index_1D << endl;
                cout << "vaaaal" << data[index_1D] << endl;
                middleColAvgTemp += data[index_1D];
                */
            }
      }
      /*
      for (int i = 0; i < m_num; i++)
      {
        int column_rank = i % out_size_rows;
        if (middle_column == column_rank || i == 0)
        {
          middle_ranks.push_back(i);
        }
      }
      */
    }
    /*
    for(size_t i = 2; i < enlarged_tile_size_rows - 2; ++i)
    {
          for(size_t j = 2; j < enlarged_tile_size_cols - 2; ++j)
          {
              int index_1D = i * enlarged_tile_size_rows + j;
              cout << "IMDEX" << index_1D << endl;
              cout << "vaaaal" << data[index_1D] << endl;
              middleColAvgTemp += data[index_1D];
          }
    }*/

    result = middleColAvgTemp / (enlarged_tile_size_rows - 4);
    return result;
}
