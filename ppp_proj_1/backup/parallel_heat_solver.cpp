/**
 * @file    parallel_heat_solver.cpp
 * @author  xforto00 <xforto00@stud.fit.vutbr.cz>
 *
 * @brief   Course: PPP 2021/2022 - Project 1
 *
 * @date    2022-04-11
 */

#include "parallel_heat_solver.h"

using namespace std;

#define FROM_DOWN_RANK_TAG_ROW 1
#define FROM_UPPER_RANK_TAG_ROW 2

#define FROM_LEFT_RANK_TAG_COL 3
#define FROM_RIGHT_RANK_TAG_COL 4

#define AVERAGE_TEMP_TAG 5

#define FROM_DOWN_LEFT_RANK_TAG 6
#define FROM_DOWN_RIGHT_RANK_TAG 7
#define FROM_UPPER_LEFT_RANK_TAG 8
#define FROM_UPPER_RIGHT_RANK_TAG 9

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
     m_fileHandle(H5I_INVALID_HID, static_cast<void (*)(hid_t )>(nullptr)),
     m_fileHandle_iteration(H5I_INVALID_HID, static_cast<void (*)(hid_t )>(nullptr))
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

      /*
      tuning of writing to file parallel
      Metadata Read Storm Problem (II)
      • Metadata read operations are treated by the library as independent read operations.
      • Consider a very large MPI job size where all processes want to open a dataset that already exists in the file.
      • All processes
        – Call H5Dopen(“/G1/G2/D1”);
        – Read the same metadata to get to the dataset and the metadata of the dataset itself
        • IF metadata not in cache, THEN read it from disk.
        – Might issue read requests to the file system for the same small metadata.

      Avoiding a Read Storm
      • Hint that metadata access is done collectively
      – H5Pset_coll_metadata_write, H5Pset_all_coll_metadata_ops
      • A property on an access property list
      • If set on the file access property list, then all metadata read operations will be
      required to be collective
      • Can be set on individual object property list
      • If set, MPI rank 0 will issue the read for a metadata entry to the file system and
      broadcast to all other ranks

      source: https://www.hdfgroup.org/wp-content/uploads/2020/02/20200206_ECPTutorial-final.pdf
      */
      H5Pset_coll_metadata_write(access_property_list, true);
      H5Pset_all_coll_metadata_ops(access_property_list, true);
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

    if (start_new_vector == true)
    {
      row_items_count = 0;
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

    if (start_new_vector == true)
    {
      row_items_count = 0;
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
/*
// 2D tile is enlarged with borders set to 0
// eg.

original tile:
0.002136 0.002136 0.000051 0.000051
0.000051 0.000051 0.000051 0.000051
0.002136 0.002136 0.000051 0.000051
0.000051 0.000051 0.000051 0.000051

enlarged tile with their borders:
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
0.000000 0.000000 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
*/
vector<int> ParallelHeatSolver::EnlargeTile(list<vector<int>> input_list, int local_tile_size_cols)
{
  vector<int> output_vector_with_borders;
  list<vector<int>> enlarged_list;

  for (auto vect : input_list) {
      // Each element of the list is
      // a vector itself
      vector<int> currentVector = vect;

      currentVector.insert(currentVector.begin(), 0);
      currentVector.insert(currentVector.begin(), 0);
      currentVector.push_back(0);
      currentVector.push_back(0);
      enlarged_list.push_back(currentVector);
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
      /*
      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';

      cout << ']';
      cout << '\n';
      */
  }

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

      currentVector.insert(currentVector.begin(), 0);
      currentVector.insert(currentVector.begin(), 0);
      currentVector.push_back(0);
      currentVector.push_back(0);
      enlarged_list.push_back(currentVector);
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
      /*
      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';

      cout << ']';
      cout << '\n';
      */
  }

  return output_vector_with_borders;
}

/*
// 1D tile is enlarged with borders set to 0
// eg.

original tile:
0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051

enlarged tile with their borders:
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
*/
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

      /*
      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';

      cout << ']';
      cout << '\n';
      */
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
      /*
      cout << "[ ";

      // Printing vector contents
      for (auto element : currentVector)
          cout << element << ' ';

      cout << ']';
      cout << '\n';
      */
  }

  return output_vector_with_borders;
}


// debug function for printing 1D array in visual form of 2D array
void ParallelHeatSolver::print_array(int* arr, int width, int height)
{
    for (int i = 0; i < width * height; i++)
    {
        if ((i != 0) && (i % width == 0))
        {
          printf("\n");
        }
        printf("%d ", arr[i]);
    }
    putchar('\n');
}

void ParallelHeatSolver::print_array(float* arr, int width, int height)
{
    for (int i = 0; i < width * height; i++)
    {
        if ((i != 0) && (i % width == 0))
        {
          printf("\n");
        }
        printf("%f ", arr[i]);
    }
    putchar('\n');
}

// function for counting 1D index from 2D index
int ParallelHeatSolver::count_1D_index(int row, int length_of_row, int column)
{
  int index_1D = (row * length_of_row) + column;
  return index_1D;
}

// function for converting tile with borders to tile without borders, used before paralell writing into output file
float* ParallelHeatSolver::TrimTileWithoutBorders(float* arr, int enlarged_tile_size_rows, int enlarged_tile_size_cols, float* result)
{
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

  //print_array(result, enlarged_tile_size_rows - 4, enlarged_tile_size_cols - 4);
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
    int local_tile_size; // local_tile_size_cols * local_tile_size_rows
    int domain_length; // length of whole domain (e.g. if domain is 16 x 16, domain_length is 256)

    if (m_rank == 0)
    {
      m_simulationProperties.GetDecompGrid(out_size_cols, out_size_rows);
      //cout << "COLS " << out_size_cols << endl;
      //cout << "ROWS " << out_size_rows << endl;
      local_tile_size_cols = sqrt(m_materialProperties.GetInitTemp().size()) / out_size_cols;
      local_tile_size_rows = sqrt(m_materialProperties.GetInitTemp().size()) / out_size_rows;
      local_tile_size = local_tile_size_cols * local_tile_size_rows;
      domain_length = m_materialProperties.GetInitTemp().size();

      /* recalculating params in case of 1D decomposition and getting correct composition of ranks
      e.g. 1D decomposition for 32 ranks and 16 x 16 grid, each rank has half of row of grid (16^2 / 32)

0       +----+----+
        | 0  | 1  |
1       +----+----+
        | 2  | 3  |
2       +----+----+
        | 4  | 5  |
3       +----+----+
        | 6  | 7  |
4       +----+----+
        | 8  | 9  |
5       +----+----+
        | 10 | 11 |
6       +----+----+
        | 12 | 13 |
7       +----+----+
        | 14 | 15 |
8       +----+----+
        | 16 | 17 |
9       +----+----+
        | 18 | 19 |
10      +----+----+
        | 20 | 21 |
11      +----+----+
        | 22 | 23 |
12      +----+----+
        | 24 | 25 |
13      +----+----+
        | 26 | 27 |
14      +----+----+
        | 28 | 29 |
15      +----+----+
        | 30 | 31 |
16      +----+----+
*/
      if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
      {
        local_tile_size_cols = domain_length / out_size_cols;
        local_tile_size_rows = 1;
        out_size_cols = sqrt(domain_length) / local_tile_size_cols;
        out_size_rows = m_size / out_size_cols;

        /*
        cout << "Out sizeh" << endl;
        cout << out_size_cols << endl;
        cout << out_size_rows << endl;
        cout << local_tile_size_cols << endl;
        cout << local_tile_size_rows << endl;
        */
      }

      // 1 rank has more rows of domain
      // e.g. 8 ranks for 16 x 16 grid, each rank has 2 rows of 16 x 16 grid, tile has dimensions 16 x 2, ranks are ordered into 1 x 8 decomposition of domain
      /*
0       +---+
1       | 1 |
2       +---+
3       | 2 |
4       +---+
5       | 3 |
6       +---+
7       | 4 |
8       +---+
9       | 5 |
10      +---+
11      | 6 |
12      +---+
13      | 7 |
14      +---+
15      | 8 |
16      +---+
*/
      else if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size < sqrt(domain_length))
      {
        local_tile_size_cols = sqrt(domain_length);
        local_tile_size_rows = sqrt(domain_length) / m_size;
        out_size_cols = 1;
        out_size_rows = sqrt(domain_length) / local_tile_size_rows;

        /*
        cout << "Out size" << endl;
        cout << out_size_cols << endl;
        cout << out_size_rows << endl;
        cout << local_tile_size_cols << endl;
        cout << local_tile_size_rows << endl;
        */
      }
    }

    MPI_Bcast(&out_size_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&out_size_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_tile_size_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&domain_length, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<float> init_temp;
    vector<float> tmp_vector;
    vector<float> out_result;
    vector<int> domain_map;
    vector<float> domain_params;

    if (m_rank == 0)
    {
      for (size_t i = 0; i < m_materialProperties.GetInitTemp().size(); i++)
      {
        float value = m_materialProperties.GetInitTemp().at(i);
        init_temp.push_back(value);
      }

      for (size_t i = 0; i < m_materialProperties.GetDomainParams().size(); i++)
      {
        float value = m_materialProperties.GetDomainParams().at(i);
        domain_params.push_back(value);
      }

      for (size_t i = 0; i < m_materialProperties.GetDomainMap().size(); i++)
      {
        int value = m_materialProperties.GetDomainMap().at(i);
        domain_map.push_back(value);
      }

    string vector_res = "";
    for (size_t i = 0; i < domain_map.size(); i++)
    {
      int item = domain_map.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }

    //cout << "DOMAIN PARAMS ================" << endl;
    //cout << vector_res << endl;

    }

    int m_num;

    MPI_Comm_size(MPI_COMM_WORLD, &m_num);
    int displacements[m_num];
    int counts[m_num];
    fill_n(counts, m_num, 1);

    int res = 0;
    bool skip_to_next_row = false;

    // counting of displacements for scatterv, displacements are different for 1D and 2D decomposition
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
    else
    {
      for(int x = 0; x < out_size_rows; x++)
      {
        for(int y = 0; y < out_size_cols; y++)
        {
            //cout << x << ", " << y<< endl;
            //cout << "Process on index" << x*out_size_cols+y <<  endl;
            //int res =  x * sqrt(domain_length) * local_tile_size_cols + y * local_tile_size_rows;
            if (x == 0 && y == 0)
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
            //cout << "Result is " << res << endl;
        }
        skip_to_next_row = true;
        //res = sqrt(domain_length) * local_tile_size_rows;
    }
  }

    // compute size of tile with their borders
    int enlarged_tile_size = (local_tile_size_cols + 4) * (local_tile_size_rows + 4);
    int enlarged_tile_size_cols = local_tile_size_cols + 4;
    int enlarged_tile_size_rows = local_tile_size_rows + 4;

    //cout << "Enlarged sizes cols " << enlarged_tile_size_cols << endl;
    //cout << "Enlarged sizes rows " << enlarged_tile_size_rows << endl;

    // allocation for tile with borders for rank
    float *init_temp_local = (float*)malloc(enlarged_tile_size * sizeof(float));
    int *domain_map_local = (int*)malloc(enlarged_tile_size * sizeof(int));
    float *domain_params_local = (float*)malloc(enlarged_tile_size * sizeof(float));

    // datatype for tile without borders
    MPI_Datatype tile_t, resized_tile_t;
    MPI_Type_vector(local_tile_size_rows, local_tile_size_cols, sqrt(domain_length), MPI_FLOAT, &tile_t);
    MPI_Type_create_resized(tile_t, 0, sizeof(float), &resized_tile_t);
    MPI_Type_commit(&tile_t);
    MPI_Type_commit(&resized_tile_t);

    // rank 0 send each rank their tile values
    MPI_Scatterv(&(init_temp[0]), counts, displacements, resized_tile_t, &(init_temp_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_map[0]), counts, displacements, resized_tile_t, &(domain_map_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(domain_params[0]), counts, displacements, resized_tile_t, &(domain_params_local[0]), local_tile_size_cols * local_tile_size_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /*
    for (int i = 0; i < m_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (i == m_rank) {
            printf("\nRank: %d\n", m_rank);
            print_array(domain_params_local, local_tile_size_cols, local_tile_size_rows);
        }
    }
    */

    vector<int> middle_ranks;

    // FIND RANKS THAT HOLDS MIDDLE COLUMN OF BOARD AND CREATE NEW COMMUNICATOR FOR THEM
    int middle_item_col_id = (sqrt(domain_length) / 2); // middle column of values
    int middle_item_tile_col_id = middle_item_col_id % local_tile_size_cols; // index of middle column of values in their tile
    int middle_col_id = middle_item_col_id / local_tile_size_cols; // column of ranks that will be in special communicator

    for (int i = 0; i < m_num; i++) {
      int column_rank = i % out_size_cols;
      if (middle_col_id == column_rank || i == 0)
      {
        //cout << "Middle process rank is " << i << endl;
        middle_ranks.push_back(i);
      }
    }

    /*
    string vector_res = "";
    cout << m_rank << endl;
    for (size_t i = 0; i < middle_ranks.size(); i++)
    {
      int item = middle_ranks.at(i);
      vector_res.append(to_string(item));
      vector_res.append(", ");
    }

    cout << vector_res << endl;
    */

    MPI_Group WORLD_GROUP;
    MPI_Group MIDDLE_COLUMN_GROUP;
    MPI_Comm MPI_COMM_MIDDLE_COLUMN;

    MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP);
    MPI_Group_incl(WORLD_GROUP, middle_ranks.size(), middle_ranks.data(), &MIDDLE_COLUMN_GROUP);
    MPI_Comm_create(MPI_COMM_WORLD, MIDDLE_COLUMN_GROUP, &MPI_COMM_MIDDLE_COLUMN);

    vector<float> init_temp_local_with_borders;
    vector<float> domain_params_local_with_borders;
    vector<int> domain_map_local_with_borders;

    // tile is enlarged with constant borders
    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
    {
      domain_map_local_with_borders = Enlarge1DTile(domain_map_local, local_tile_size_cols);
      domain_params_local_with_borders = Enlarge1DTile(domain_params_local, local_tile_size_cols);
      init_temp_local_with_borders = Enlarge1DTile(init_temp_local, local_tile_size_cols);
    }
    else
    {
      list<vector<int>> domain_map_local_list = SplitRows(domain_map_local, local_tile_size, local_tile_size_rows);
      domain_map_local_with_borders = EnlargeTile(domain_map_local_list, local_tile_size_cols);
      list<vector<float>> domain_params_local_list = SplitRows(domain_params_local, local_tile_size, local_tile_size_cols);
      domain_params_local_with_borders = EnlargeTile(domain_params_local_list, local_tile_size_cols);
      list<vector<float>> init_temp_local_list = SplitRows(init_temp_local, local_tile_size, local_tile_size_cols);
      init_temp_local_with_borders = EnlargeTile(init_temp_local_list, local_tile_size_cols);
    }

    // datatypes for row of matrix and block of 2 columns of matrix
    MPI_Datatype tile_row_t, tile_col_t;
    MPI_Type_contiguous(enlarged_tile_size_cols, MPI_FLOAT, &tile_row_t);
    MPI_Type_commit(&tile_row_t);

    MPI_Type_vector(enlarged_tile_size_rows, 2, enlarged_tile_size_cols, MPI_FLOAT, &tile_col_t);
    MPI_Type_commit(&tile_col_t);

    // compute 2D position of each rank
    int row_id = m_rank / out_size_cols;
    int col_id = m_rank % out_size_cols;

    // compute 1D position of rank's neighbors
    int down_rank = m_rank + out_size_cols;
    int upper_rank = m_rank - out_size_cols;
    int left_rank = m_rank - 1;
    int right_rank = m_rank + 1;

    int down_left_rank = m_rank + out_size_cols - 1;
    int down_right_rank = m_rank + out_size_cols + 1;
    int upper_left_rank = m_rank - out_size_cols - 1;
    int upper_right_rank = m_rank - out_size_cols + 1;

    /*
      SENDING NEIGHBOUR VALUES TO BORDERS OF TILE FOR EACH INIT TEMP, DOMAIN PARAMS AND DOMAIN MAP (process for sending values to tile borders is same for all 3 arrays: init temp, domain map and domain params)
      FIRSTLY ROWS AND COLUMNS FROM NEIGHBOURS ARE SENT AND INSERTED INTO BORDERS
      THEN UPPER LEFT, UPPER RIGHT, DOWN LEFT AND DOWN RIGHT CORNERS ARE SEND

      e.g. domain params for rank 5, 16 x 16 domain and 16 ranks in grid:

      enlarged tile with borders:
      0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
      0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
      0.000000 0.000000 0.002136 0.002136 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.000051 0.000051 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.002136 0.002136 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.000051 0.000051 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
      0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000

      after sending rows and columns:
      0.000000 0.000000 0.002916 0.002916 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.002136 0.002136 0.002916 0.002916 0.000000 0.000000
      0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916
      0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916
      0.000000 0.000000 0.002136 0.002136 0.002916 0.002916 0.000000 0.000000
      0.000000 0.000000 0.000051 0.000051 0.002916 0.002916 0.000000 0.000000

      after sending upper left, upper right, down left and down right corners from correct neighbors:
      0.000051 0.000051 0.002916 0.002916 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916
      0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916
      0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916
      0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916

    */

    vector<float> init_temp_local_with_borders_recieved;
    for (size_t i = 0; i < init_temp_local_with_borders.size(); i++)
    {
      float value = init_temp_local_with_borders.at(i);
      init_temp_local_with_borders_recieved.push_back(value);
    }

    // each enlarge tile can be divided into 2 x 2 tiles
    MPI_Datatype enlarge_tile_part_t;
    MPI_Type_vector(2, 2, enlarged_tile_size_cols, MPI_FLOAT, &enlarge_tile_part_t);
    MPI_Type_commit(&enlarge_tile_part_t);

    // counting of rank's requests during Isend and Irecv for sending columns and rows
    int total_request_count = 0;

    if (row_id == 0 || row_id == out_size_rows - 1)
    {

      total_request_count += 2;
    }
    else
    {
      total_request_count += 4;
    }

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_2D || (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && out_size_cols != 1))
    {
      if (col_id == 0 || col_id == out_size_cols - 1)
      {
        total_request_count += 2;
      }
      else
      {
        total_request_count += 4;
      }
    }

    // counting of rank's requests during Isend and Irecv upper left, upper right, down left and down right corners
    int edges_request_count = 0;
    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_2D || (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && out_size_cols != 1))
    {
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

    //cout << "Upcoming requests for rank " << m_rank << ":   " << total_request_count << endl;

    // SENDING NEIGHBOUR VALUES TO BORDERS OF TILE FOR INIT TEMP

    // sending rows and columns from upper, down, left and right neighbors
    MPI_Request requests[total_request_count];
    int num_requests = 0;

    if (row_id != out_size_rows - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (row_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (col_id != 0)
    {
        MPI_Irecv(&init_temp_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    if (row_id != 0)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&init_temp_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&init_temp_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests[num_requests++]);
    }

    MPI_Waitall(total_request_count, requests, MPI_STATUSES_IGNORE);

    /*
    specially for 1D decomposition when m_size >= sqrt(domain_length) (tile is one row (y dimension is 1)):
    enlarged tile with their borders (domain params example):
    // this is row overlapping my neighbor's upper neighbor rank (the rank is not same of these two upper borders)  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    // this is row overlapping my upper neighbor rank (as in 2D)                                                    0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    // this is my original tile with my values                                                                      0.000000 0.000000 0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
    // this is row overlapping my upper neighbor rank (as in 2D)                                                    0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    // this is row overlapping my neighbor's down neighbor rank (the rank is not same of these two down borders)    0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000

    during this process ranks recieve correct row from neighbor's upper neighbor rank or neighbor's down neighbor rank and save them to first or last row of their enlarged tile,
    two borders from upper and down neighbor ranks were already send before and are simmilar to 2D

    e.g. for domain 16 x 16, 1D decomposition for 16 ranks (each rank has one row of domain), domain params, situation for rank 4:
    already sended rows from neighbors ranks 3 and 5:
    RANK: 4
    0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916 0.000051 0.000051 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000

    now we send rows from neighbor's neighbors rank 2 and 6:
    RANK: 4
    0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916 0.002916 0.002916 0.002916 0.002916 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.000051 0.000051 0.000051 0.000051 0.002916 0.002916 0.002916 0.002916 0.000051 0.000051 0.000051 0.000051 0.000051 0.000051 0.000000 0.000000
    0.000000 0.000000 0.000051 0.000051 0.002136 0.002136 0.002136 0.002136 0.002916 0.002916 0.002916 0.002916 0.002136 0.002136 0.002136 0.002136 0.000051 0.000051 0.000000 0.000000
    */
    int rows_second_neighbor_request_count = 0;
    // requests for Irecv
    if (row_id == 0 || row_id == out_size_rows - 1 || row_id - 1 == 0 || row_id + 1 == out_size_rows - 1)
    {
      rows_second_neighbor_request_count += 1;
    }
    else
    {
      rows_second_neighbor_request_count += 2;
    }
    // requests for Isend
    if (row_id == 0 || row_id == out_size_rows - 1 || row_id - 1 == 0 || row_id + 1 == out_size_rows - 1)
    {
      rows_second_neighbor_request_count += 1;
    }
    else
    {
      rows_second_neighbor_request_count += 2;
    }

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
    {
      num_requests = 0;
      MPI_Request requests_init_temp_recieved_1D[rows_second_neighbor_request_count];

      vector<float> init_temp_local_with_borders_recieved_1D;

      for (size_t i = 0; i < init_temp_local_with_borders_recieved.size(); i++)
      {
        float value = init_temp_local_with_borders_recieved.at(i);
        init_temp_local_with_borders_recieved_1D.push_back(value);
      }

      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Irecv(&init_temp_local_with_borders_recieved_1D[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank + out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_init_temp_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Irecv(&init_temp_local_with_borders_recieved_1D[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank - out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_init_temp_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Isend(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank - out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_init_temp_recieved_1D[num_requests++]);
      }
      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Isend(&init_temp_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, down_rank + out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_init_temp_recieved_1D[num_requests++]);
      }

      MPI_Waitall(rows_second_neighbor_request_count, requests_init_temp_recieved_1D, MPI_STATUSES_IGNORE);

      init_temp_local_with_borders_recieved.clear();
      for (size_t i = 0; i < init_temp_local_with_borders_recieved_1D.size(); i++)
      {
        float value = init_temp_local_with_borders_recieved_1D.at(i);
        init_temp_local_with_borders_recieved.push_back(value);
      }
    }

    // sending upper left, upper right, down left and down right corners from correct neighbors:
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

    MPI_Waitall(edges_request_count, requests_init_temp_recieved, MPI_STATUSES_IGNORE);


    // SENDING NEIGHBOUR VALUES TO BORDERS OF TILE FOR DOMAIN MAP
    vector<int> domain_map_local_with_borders_recieved;
    for (size_t i = 0; i < domain_map_local_with_borders.size(); i++)
    {
      int value = domain_map_local_with_borders.at(i);
      domain_map_local_with_borders_recieved.push_back(value);
    }

    MPI_Request requests_domain_map[total_request_count];
    num_requests = 0;

    if (row_id != out_size_rows - 1)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (row_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (col_id != 0)
    {
        MPI_Irecv(&domain_map_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_map_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
      /// -------->
        MPI_Isend(&domain_map_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_map[num_requests++]);
    }

    MPI_Waitall(total_request_count, requests_domain_map, MPI_STATUSES_IGNORE);

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
    {
      num_requests = 0;
      MPI_Request requests_domain_map_recieved_1D[rows_second_neighbor_request_count];

      vector<int> domain_map_local_with_borders_recieved_1D;

      for (size_t i = 0; i < domain_map_local_with_borders_recieved.size(); i++)
      {
        int value = domain_map_local_with_borders_recieved.at(i);
        domain_map_local_with_borders_recieved_1D.push_back(value);
      }

      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Irecv(&domain_map_local_with_borders_recieved_1D[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank + out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Irecv(&domain_map_local_with_borders_recieved_1D[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank - out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Isend(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank - out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map_recieved_1D[num_requests++]);
      }
      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Isend(&domain_map_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, down_rank + out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_map_recieved_1D[num_requests++]);
      }

      MPI_Waitall(rows_second_neighbor_request_count, requests_domain_map_recieved_1D, MPI_STATUSES_IGNORE);

      domain_map_local_with_borders_recieved.clear();
      for (size_t i = 0; i < domain_map_local_with_borders_recieved_1D.size(); i++)
      {
        int value = domain_map_local_with_borders_recieved_1D.at(i);
        domain_map_local_with_borders_recieved.push_back(value);
      }
    }

    vector<int> domain_map_local_with_borders_recieved_with_edges;

    for (size_t i = 0; i < domain_map_local_with_borders_recieved.size(); i++)
    {
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

    MPI_Waitall(edges_request_count, requests_domain_map_recieved, MPI_STATUSES_IGNORE);

    // SENDING NEIGHBOUR VALUES TO BORDERS OF TILE FOR DOMAIN PARAMS
    vector<float> domain_params_local_with_borders_recieved;

    for (size_t i = 0; i < domain_params_local_with_borders.size(); i++) {
      float value = domain_params_local_with_borders.at(i);
      domain_params_local_with_borders_recieved.push_back(value);
    }

    MPI_Request requests_domain_params[total_request_count];
    num_requests = 0;

    if (row_id != out_size_rows - 1)
    {
        // store to last two rows
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (row_id != 0)
    {
        // store to last two rows
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (col_id != 0)
    {
        MPI_Irecv(&domain_params_local_with_borders_recieved[0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    if (row_id != 0)
    {
        // send two rows up from down rank
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (row_id != out_size_rows - 1)
    {
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (col_id != 0)
    {
         // <--------
        MPI_Isend(&domain_params_local_with_borders[2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }
    if (col_id != out_size_cols - 1)
    {
       /// -------->
        MPI_Isend(&domain_params_local_with_borders[enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_domain_params[num_requests++]);
    }

    MPI_Waitall(total_request_count, requests_domain_params, MPI_STATUSES_IGNORE);

    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
    {
      num_requests = 0;
      MPI_Request requests_domain_params_recieved_1D[rows_second_neighbor_request_count];

      vector<float> domain_params_local_with_borders_recieved_1D;

      for (size_t i = 0; i < domain_params_local_with_borders_recieved.size(); i++)
      {
        float value = domain_params_local_with_borders_recieved.at(i);
        domain_params_local_with_borders_recieved_1D.push_back(value);
      }

      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Irecv(&domain_params_local_with_borders_recieved_1D[enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank + out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Irecv(&domain_params_local_with_borders_recieved_1D[enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank - out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params_recieved_1D[num_requests++]);
      }
      if (row_id != 0 && row_id - 1 != 0)
      {
          MPI_Isend(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank - out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params_recieved_1D[num_requests++]);
      }
      if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
      {
          MPI_Isend(&domain_params_local_with_borders_recieved[enlarged_tile_size_cols * 2], 1, tile_row_t, down_rank + out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_domain_params_recieved_1D[num_requests++]);
      }

      MPI_Waitall(rows_second_neighbor_request_count, requests_domain_params_recieved_1D, MPI_STATUSES_IGNORE);

      domain_params_local_with_borders_recieved.clear();
      for (size_t i = 0; i < domain_params_local_with_borders_recieved_1D.size(); i++)
      {
        float value = domain_params_local_with_borders_recieved_1D.at(i);
        domain_params_local_with_borders_recieved.push_back(value);
      }
    }

    vector<float> domain_params_local_with_borders_recieved_with_edges;

    for (size_t i = 0; i < domain_params_local_with_borders_recieved.size(); i++)
    {
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

    MPI_Waitall(edges_request_count, requests_domain_params, MPI_STATUSES_IGNORE);

    /*
    // debug part of code for printing tile for each rank and see whether each rank has correct tile values
    for (int i = 0; i < m_num; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (i == m_rank) {
            printf("\nRANK: %d\n", m_rank);
            print_array(&domain_params_local_with_borders_recieved_with_edges[0], enlarged_tile_size_cols, enlarged_tile_size_rows);
        }
    }
    */


    /*
    if (m_rank == 6) {
    //cout << "Process with received values params" << m_rank << endl;

    for (size_t i = 0; i < domain_params_local_with_borders_recieved.size(); i++)
    {
      float item = domain_params_local_with_borders_recieved.at(i);
      vector_res.append(to_string(item));
      vector_res.append(";; ");
    }

    cout << vector_res << endl;}
    */
    // allocation for window for RMA
    MPI_Win win;
    float* win_memory;
    MPI_Win_allocate(2 * enlarged_tile_size * sizeof(float), sizeof(float), MPI_INFO_NULL, MPI_COMM_WORLD, &win_memory, &win);

    float *workTempArrays[2];
    // init arrays for new and old values
    if (m_simulationProperties.IsRunParallelP2P())
    {
      workTempArrays[0] = init_temp_local_with_borders_recieved_with_edges.data();
      workTempArrays[1] = init_temp_local_with_borders_recieved_with_edges.data();
    }
    else if (m_simulationProperties.IsRunParallelRMA())
    {
      float* win_memory_data1 = win_memory;
      float* win_memory_data2 = win_memory + enlarged_tile_size;

      for (int i = 0; i < enlarged_tile_size; i++)
      {
        win_memory_data1[i] = init_temp_local_with_borders_recieved_with_edges[i];
      }
      copy(win_memory_data1, win_memory_data1 + enlarged_tile_size, win_memory_data2);
      workTempArrays[0] = win_memory_data1;
      workTempArrays[1] = win_memory_data2;
    }

    // compute offset in case of corner ranks because corners of size of 2 of whole domain are static
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

    // specially for 1D decomposition if m_size >= sqrt(domain_length) cause not only first and last row ranks but also second and last but one (predposledni) ranks contains static border corners
    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && row_id - 1 == 0 && m_size >= sqrt(domain_length))
    {
      offset_rows_begin += 1;
    }
    if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && row_id + 1 == out_size_rows - 1 && m_size >= sqrt(domain_length))
    {
      offset_rows_end += 1;
    }

    /*
    if (m_rank == 2)
    {
      cout << "offset beg row " << offset_rows_begin << endl;
      cout << "offset end row " << offset_rows_end << endl;
      cout << "offset beg col " << offset_cols_begin << endl;
      cout << "offset end col " << offset_cols_end << endl;
    }
    */


    // datatypes for whole domain and tile without their borders - used for correct gatherv when rank 0 writes whole result into output file
    MPI_Datatype worker_tile_t;
    int tile[2];
    int worker_size[2];
    int worker_start[2];
    worker_size[0] = enlarged_tile_size_rows;
    worker_size[1] = enlarged_tile_size_cols;
    tile[0] = local_tile_size_rows;
    tile[1] = local_tile_size_cols;
    worker_start[0] = 2; // cut borders of size of 2
    worker_start[1] = 2;
    float result_domain[domain_length]; // result of whole domain temperature to write by rank 0 sequentially to file

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

    //MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &result_domain[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);

    // setting of paralel writing into output file
    hid_t filespace;
    hid_t memspace;
    hid_t dataset;
    hid_t xferPList;

    hsize_t datasetSize[] = {hsize_t(sqrt(domain_length)), hsize_t(sqrt(domain_length))};
    hsize_t memSize[]     = {hsize_t(local_tile_size_rows), hsize_t(local_tile_size_cols)};
    hsize_t slabStart[] = {hsize_t(row_id * local_tile_size_rows), hsize_t(col_id * local_tile_size_cols)};
    hsize_t slabSize[]  = {hsize_t(local_tile_size_rows), hsize_t(local_tile_size_cols)};

    if (!m_simulationProperties.GetOutputFileName().empty() && m_simulationProperties.IsUseParallelIO())
    {
        filespace = H5Screate_simple(2, datasetSize, nullptr);
        memspace  = H5Screate_simple(2, memSize,     nullptr);

        xferPList = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(xferPList, H5FD_MPIO_COLLECTIVE);

        dataset = H5Dcreate(m_fileHandle, "result_domain_temperature", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, slabStart, nullptr, slabSize, nullptr);
    }

    /*
    if (m_rank == 0)
    {
      cout << "Printing whole domain" << endl;
      print_array(result_domain, sqrt(domain_length), sqrt(domain_length));

      cout << "Printing my tile without borders" << endl;
      float *result_tile;
      float result[local_tile_size];
      result_tile = TrimTileWithoutBorders(&domain_params_local_with_borders_recieved_with_edges[0], enlarged_tile_size_rows, enlarged_tile_size_cols, result);
      print_array(result_tile, local_tile_size_rows, local_tile_size_cols);
    }
    */

    float whole_domain_temp_snapshot[domain_length]; // snapshot of whole domain result for outputing into file by rank 0 sequentially

    int iteration_offset = 0; // offset of iteration when using RMA

    /*
    ACTUAL SIMULATION
    */
    double simulation_start_time;
    if (m_rank == 0)
    {
      simulation_start_time = MPI_Wtime();
    }

    float middleColAvgTemp = 0.0f;
    float final_iteration_temp = 0.0f;

    for (size_t iter = 0; iter < m_simulationProperties.GetNumIterations(); ++iter)
    {
      // compute tile
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

      // change tile borders in p2p mode
      if (m_simulationProperties.IsRunParallelP2P())
      {
        MPI_Request requests_simulation[total_request_count];
        int num_requests_simulation = 0;

        if (row_id != out_size_rows - 1)
        {
          // store to last two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 1)], 2, tile_row_t, down_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (row_id != 0)
        {
          // store to first two rows
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * 0], 2, tile_row_t, upper_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (col_id != out_size_cols - 1)
        {
          MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols - 1 - 1], 1, tile_col_t, right_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (col_id != 0)
        {
          MPI_Irecv(&workTempArrays[1][0], 1, tile_col_t, left_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }

        if (row_id != 0)
        {
          // send two rows up from down rank
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 2], 2, tile_row_t, upper_rank, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (row_id != out_size_rows - 1)
        {
          // send two rows down
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1 - 3)], 2, tile_row_t, down_rank, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (col_id != 0)
        {
           // <--------
          MPI_Isend(&workTempArrays[0][2], 1, tile_col_t, left_rank, FROM_RIGHT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }
        if (col_id != out_size_cols - 1)
        {
        /// -------->
          MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols - 1 - 3], 1, tile_col_t, right_rank, FROM_LEFT_RANK_TAG_COL, MPI_COMM_WORLD, &requests_simulation[num_requests_simulation++]);
        }

        MPI_Waitall(total_request_count, requests_simulation, MPI_STATUSES_IGNORE);

        // for 1D if m_size >= sqrt(domain_length): tile is not only send to borders to my upper or down rank but also to borders of neighbor's upper or down rank
        // borders are basically overlapping not two ranks as in 2D but three ranks, because size of border is 2 and size of actual tile is 1
        if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
        {
          num_requests = 0;
          MPI_Request requests_simulation_1D[rows_second_neighbor_request_count];

          if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
          {
              MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * (enlarged_tile_size_rows - 1)], 1, tile_row_t, down_rank + out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation_1D[num_requests++]);
          }
          if (row_id != 0 && row_id - 1 != 0)
          {
              MPI_Irecv(&workTempArrays[1][enlarged_tile_size_cols * 0], 1, tile_row_t, upper_rank - out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation_1D[num_requests++]);
          }
          if (row_id != 0 && row_id - 1 != 0)
          {
              MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank - out_size_cols, FROM_DOWN_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation_1D[num_requests++]);
          }
          if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
          {
              MPI_Isend(&workTempArrays[0][enlarged_tile_size_cols * 2], 1, tile_row_t, down_rank + out_size_cols, FROM_UPPER_RANK_TAG_ROW, MPI_COMM_WORLD, &requests_simulation_1D[num_requests++]);
          }

          MPI_Waitall(rows_second_neighbor_request_count, requests_simulation_1D, MPI_STATUSES_IGNORE);
        }

      }

      else // change tile borders in rma mode
      {
        int offset = enlarged_tile_size * iteration_offset;
        MPI_Win_fence(0, win);

        // The execution of a put operation is similar to the execution of a send by the origin process and a matching receive by the target process. The obvious difference is that all arguments are provided by one call --- the call executed by the origin process.
        if (row_id != 0)
        {
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

        if (m_simulationProperties.GetDecompMode() == SimulationProperties::DECOMP_MODE_1D && m_size >= sqrt(domain_length))
        {
          MPI_Win_fence(0, win);

          if (row_id != 0 && row_id - 1 != 0)
          {
              MPI_Put(&workTempArrays[0][enlarged_tile_size_cols * 2], 1, tile_row_t, upper_rank - out_size_cols, enlarged_tile_size_cols * (enlarged_tile_size_rows - 1) + offset, 1, tile_row_t, win);
          }
          if (row_id != out_size_rows - 1 && row_id + 1 != out_size_rows - 1)
          {
              MPI_Put(&workTempArrays[0][enlarged_tile_size_cols * 2], 1, tile_row_t, down_rank + out_size_cols, enlarged_tile_size_cols * 0 + offset, 1, tile_row_t, win);
          }

          MPI_Win_fence(0, win);
        }

        // switch offset for next iteration
        if (iteration_offset == 1)
        {
          iteration_offset = 0;
        }
        else if (iteration_offset == 0)
        {
          iteration_offset = 1;
        }
      }

      /*
      if (m_rank == 5)
      {
        cout << "Rank 5" << endl;
        print_array(workTempArrays[0], enlarged_tile_size_rows, enlarged_tile_size_cols);
      }
      if (m_rank == 6)
      {
        cout << "Rank 6" << endl;
        print_array(workTempArrays[1], enlarged_tile_size_rows, enlarged_tile_size_cols);
      }
      */

      final_iteration_temp = 0.0f;
      // compute middle column temperature if rank is in communicator and holds middle column of domain
      if (count(middle_ranks.begin(), middle_ranks.end(), m_rank))
      {
        //cout << "Middle " << m_rank << endl;

        if (m_rank == 0 && out_size_cols != 1) // when rank 0 doesnt hold some values from middle column
        {
          //cout << "Rank 0 is not a part of middle column computing" << endl;
          middleColAvgTemp = 0.0f;
        }
        else
        {
          middleColAvgTemp = ComputeMiddleColAvgTemp(workTempArrays[0], enlarged_tile_size_rows, enlarged_tile_size_cols, middle_item_tile_col_id);
        }

        vector<float> all_average_temp(out_size_rows + 1);
        // rank 0 gathers all average middle column temperatures from ranks in communicator
        MPI_Gather(&middleColAvgTemp, 1, MPI_FLOAT, &all_average_temp[0], 1, MPI_FLOAT, 0, MPI_COMM_MIDDLE_COLUMN);

        if (m_rank == 0)
        {
          float all_average_temp_sum = 0.0f;

          for (size_t i = 0; i < all_average_temp.size(); ++i)
          {
            if (i == 0 && out_size_cols != 1) // when rank 0 doesnt hold some values from middle column
            {
              continue;
            }
            float item = all_average_temp.at(i);
            all_average_temp_sum += item;
          }

          final_iteration_temp = all_average_temp_sum / (all_average_temp.size() - 1);
        }
      }

      // writing to output file during simulation
      hid_t filespace_iteration;
      hid_t memspace_iteration;
      hid_t dataset_iteration;
      hid_t xferPList_iteration;
      if (!m_simulationProperties.GetOutputFileName().empty())
      {
        // rank 0 gathers whole domain temp from all ranks (their tiles without borders) for sequentially saving into file
        //MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &whole_domain_temp_snapshot[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);
        /*
        if (m_rank == 0) {
          //cout << "Printing snapshot" << endl;
          //print_array(whole_domain_temp_snapshot, sqrt(domain_length), sqrt(domain_length));
        }*/
            if (!m_simulationProperties.IsUseParallelIO() && ((iter % m_simulationProperties.GetDiskWriteIntensity()) == 0))
            {
              // sequential writing by rank 0 - gathers all ranks tiles without borders to one array which is inserted into output file
              MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &whole_domain_temp_snapshot[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);
              if (m_rank == 0)
              {
                StoreDataIntoFile(m_fileHandle, iter, &whole_domain_temp_snapshot[0]);
              }
            }
            // when parallel writing is enabled
            else if (m_simulationProperties.IsUseParallelIO() && ((iter % m_simulationProperties.GetDiskWriteIntensity()) == 0))
            {
              // create new separate file for result temperature of this iteration
              string input_file_name = m_simulationProperties.GetOutputFileName("par");
              string iteration_file_name = "";
              iteration_file_name.append(input_file_name.substr(0,input_file_name.size() - 1 - 2));
              iteration_file_name.append("_iteration_");
              iteration_file_name.append(to_string(iter));
              iteration_file_name.append(".h5");

              hid_t access_property_list = H5Pcreate(H5P_FILE_ACCESS);
              H5Pset_coll_metadata_write(access_property_list, true);
              H5Pset_all_coll_metadata_ops(access_property_list, true);

              H5Pset_fapl_mpio(access_property_list, MPI_COMM_WORLD, MPI_INFO_NULL);
              m_fileHandle_iteration.Set(H5Fcreate(iteration_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_property_list), H5Fclose);
              H5Pclose(access_property_list);

              // get current temperatures in tile without their borders which are shared
              float *tile_snapshot;
              float tile_result[local_tile_size];
              tile_snapshot = TrimTileWithoutBorders(&workTempArrays[0][0], enlarged_tile_size_rows, enlarged_tile_size_cols, tile_result);
              //print_array(tile_snapshot, local_tile_size_rows, local_tile_size_cols);

              string dataset_name = "";
              dataset_name.append("iteration_");
              dataset_name.append(to_string(iter));

              // seting file for this iteration
              filespace_iteration = H5Screate_simple(2, datasetSize, nullptr);
              memspace_iteration  = H5Screate_simple(2, memSize,     nullptr);

              xferPList_iteration = H5Pcreate(H5P_DATASET_XFER);
              H5Pset_dxpl_mpio(xferPList_iteration, H5FD_MPIO_COLLECTIVE);

              dataset_iteration = H5Dcreate(m_fileHandle_iteration, dataset_name.c_str(), H5T_NATIVE_FLOAT, filespace_iteration, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

              H5Sselect_hyperslab(filespace_iteration, H5S_SELECT_SET, slabStart, nullptr, slabSize, nullptr);
              // write to separate file this iteration
              H5Dwrite(dataset_iteration, H5T_NATIVE_FLOAT, memspace_iteration, filespace_iteration, xferPList_iteration, &tile_snapshot[0]);
              // write current temparature to file which name was input via command line
              H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, filespace, xferPList, &tile_snapshot[0]);
            }

      }

      swap(workTempArrays[0], workTempArrays[1]);

      if (m_rank == 0 && ShouldPrintProgress(iter))
      {
        PrintProgressReport(iter, final_iteration_temp);
      }

      if (!m_simulationProperties.GetOutputFileName().empty() && m_simulationProperties.IsUseParallelIO() && ((iter % m_simulationProperties.GetDiskWriteIntensity()) == 0))
      {
        H5Dclose(dataset_iteration);
      }

    } // end of simulation

    double simulation_end_time;
    if (m_rank == 0)
    {
      simulation_end_time = MPI_Wtime();
      double simulation_duration = simulation_end_time - simulation_start_time;
      //cout << final_iteration_temp << endl;
      PrintFinalReport(simulation_duration, final_iteration_temp, "par");

    }

    if (!m_simulationProperties.GetOutputFileName().empty() && m_simulationProperties.IsUseParallelIO())
    {
      H5Dclose(dataset);
    }

    // temperatures from all tiles are gathered in rank 0 and returned as outResult
    MPI_Gatherv(&workTempArrays[0][0], 1, worker_tile_t, &outResult[0], counts, displacements, resized_farmer_matrix_t, 0, MPI_COMM_WORLD);

    // cleaning
    MPI_Win_free(&win);

    MPI_Type_free(&tile_t);
    MPI_Type_free(&resized_tile_t);
    MPI_Type_free(&tile_row_t);
    MPI_Type_free(&tile_col_t);
    MPI_Type_free(&enlarge_tile_part_t);
    MPI_Type_free(&worker_tile_t);
    MPI_Type_free(&farmer_matrix_t);
    MPI_Type_free(&resized_farmer_matrix_t);

    if (count(middle_ranks.begin(), middle_ranks.end(), m_rank))
    {
      MPI_Comm_free(&MPI_COMM_MIDDLE_COLUMN);
    }

    MPI_Group_free(&WORLD_GROUP);
    MPI_Group_free(&MIDDLE_COLUMN_GROUP);
}

// computing of average temparature of middle column in tile
float ParallelHeatSolver::ComputeMiddleColAvgTemp(const float *data, int enlarged_tile_size_rows, int enlarged_tile_size_cols, int middle_item_tile_col_id) const
{
    float middleColAvgTemp = 0.0f;
    float result = 0.0f;

    // index of column which will have influence on middle temperature (index computed for tile + 2 because of enlarged borders)
    int middle_column_tile_id = 2 + middle_item_tile_col_id;
    //cout << "Middle column tile id " << middle_column_tile_id <<  endl;
    //cout << "Enlarge rows " << enlarged_tile_size_rows << endl;
    //cout << "Enlarge cols " << enlarged_tile_size_cols << endl;

    for(size_t i = 2; i < enlarged_tile_size_rows - 2; ++i)
    {
        for(size_t j = 2; j < enlarged_tile_size_cols - 2; ++j)
        {
            int index_1D = i * enlarged_tile_size_cols + j;

            int row_id = index_1D / enlarged_tile_size_cols;
            int col_id = index_1D % enlarged_tile_size_cols;

            /*
            if (m_rank == 2)
            {
                float value = data[index_1D];
                //cout << "TEMP " << value << endl;
            }
            */

            // ignore border values of enlarged tile
            if ((middle_column_tile_id == col_id && row_id != 0) && (middle_column_tile_id == col_id && row_id != 1) && (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 2) && (middle_column_tile_id == col_id && row_id != enlarged_tile_size_rows - 1))
            {
                  //cout << "IMDEX" << index_1D << endl;
                  //cout << data[index_1D] << endl;
                  float value = data[index_1D];
                  //cout << "TEMP OF MIDDLE COL OF "  << m_rank << ": " <<  value << endl;
                  middleColAvgTemp += value;
            }
        }
    }

    result = middleColAvgTemp / (enlarged_tile_size_rows - 4);
    return result;
}
