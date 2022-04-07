/**
 * @file    sequential_heat_solver.cpp
 * @authors Filip Vaverka <ivaverka@fit.vutbr.cz>
 *          Jiri Jaros <jarosjir@fit.vutbr.cz>
 *          Kristian Kadlubiak <ikadlubiak@fit.vutbr.cz>
 *
 * @brief   Course: PPP 2021/2022 - Project 1
 *          This file contains implementation of sequential heat equation solver.
 *
 * @date    2022-02-03
 */

#include "sequential_heat_solver.h"

int SequentialHeatSolver::count_1D_index(int row, int length_of_row, int column)
{
  int index_1D = (row * length_of_row) + column;
  return index_1D;
}

SequentialHeatSolver::SequentialHeatSolver(SimulationProperties &simulationProps,
                                           MaterialProperties &materialProps)
    : BaseHeatSolver(simulationProps, materialProps),
      m_tempArray(materialProps.GetGridPoints()),
      m_fileHandle(H5I_INVALID_HID, static_cast<void (*)(hid_t )>(nullptr))
{
    // 1. Open output file if its name was specified.
    if(!m_simulationProperties.GetOutputFileName().empty())
        m_fileHandle.Set(H5Fcreate(simulationProps.GetOutputFileName("seq").c_str(),
                                   H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
}

void SequentialHeatSolver::RunSolver(std::vector<float, AlignedAllocator<float> > &outResult)
{
    // 2. Copy initial temperature into both working arrays
    std::copy(m_materialProperties.GetInitTemp().begin(),
              m_materialProperties.GetInitTemp().end(), m_tempArray.begin());

    string vector_res = "";

    for (size_t i = 0; i < m_tempArray.size(); i++)
    {
        float item = m_tempArray.at(i);
        vector_res.append(to_string(item));
        vector_res.append(", ");
    }

    cout << "TEMP PARAMS ================" << endl;
    cout << vector_res << endl;
    std::copy(m_materialProperties.GetInitTemp().begin(),
              m_materialProperties.GetInitTemp().end(), outResult.begin());

    float *workTempArrays[] = { m_tempArray.data(), outResult.data() };
    float middleColAvgTemp = 0.0f;
    double startTime = MPI_Wtime();

    // 3. Begin iterative simulation main loop
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

        vector_res = "";
        for (size_t i = 0; i < 256; i++)
        {
            float item = workTempArrays[0][i];
            vector_res.append(to_string(item));
            vector_res.append(", ");
        }

        cout << "TEMP PARAMS ================" << endl;
        cout << vector_res << endl;
        /*
        for(size_t i = 0; i < m_materialProperties.GetEdgeSize(); ++i)
        {
            for(size_t j = 0; j < m_materialProperties.GetEdgeSize(); ++j)
            {
                int index_1D = count_1D_index(i, m_materialProperties.GetEdgeSize(), j);
                if (i >= 0 && )
            }
        }
        */

        // 5. Compute average temperature in the middle column of the domain.
        middleColAvgTemp = ComputeMiddleColAvgTemp(workTempArrays[0]);

        // 6. Store the simulation state if appropriate (ie. every N-th iteration)
        if(m_fileHandle != H5I_INVALID_HID && ((iter % m_simulationProperties.GetDiskWriteIntensity()) == 0))
            StoreDataIntoFile(m_fileHandle, iter, workTempArrays[0]);

        // 7. Swap source and destination buffers
        std::swap(workTempArrays[0], workTempArrays[1]);

        // 8. Print current progress (prints progress only every 10% of the simulation).
        PrintProgressReport(iter, middleColAvgTemp);
    }

    // 9. Measure total execution time and report
    double elapsedTime = MPI_Wtime() - startTime;
    PrintFinalReport(elapsedTime, middleColAvgTemp, "seq");

    // 10. Copy result over if necessary (even/odd number of buffer swaps).
    if(m_simulationProperties.GetNumIterations() & 1)
        std::copy(m_tempArray.begin(), m_tempArray.end(), outResult.begin());
}

float SequentialHeatSolver::ComputeMiddleColAvgTemp(const float *data) const
{
    float middleColAvgTemp = 0.0f;
    for(size_t i = 0; i < m_materialProperties.GetEdgeSize(); ++i) {
        middleColAvgTemp += data[i*m_materialProperties.GetEdgeSize() + m_materialProperties.GetEdgeSize() / 2];
        int index = i*m_materialProperties.GetEdgeSize() + m_materialProperties.GetEdgeSize() / 2;
        cout << index << endl;
        cout << data[i*m_materialProperties.GetEdgeSize() + m_materialProperties.GetEdgeSize() / 2] << endl;
      }

    return middleColAvgTemp / float(m_materialProperties.GetEdgeSize());
}
