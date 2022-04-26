#! /usr/bin/env python3
import csv
import matplotlib.pyplot as plt

def weak_scaling(data_dict, out_prefix, suffixes, proc_counts):
    """
    Params:
    data_dict: Dictionary s puvodnimi nactenymi hodnotami ze souboru (promenna data_dict)
    out_prefix: prefix vystupniho souboru
    suffixes: typy mereni (promenna suffixes)
    proc_counts: pocty procesoru u mereni (promenna x_axis_points)
    """
    split_dict = {}
    for i in suffixes:
        split_dict[i] = {}

    for key, data_series in data_dict.items():
        key_split = key.split('_', maxsplit=1)
        prefix = int(key_split[0])
        suffix = "_" + key_split[1]
        split_dict[suffix][prefix] = data_series

    weak_scale_dict = {}
    for i in suffixes:
        weak_scale_dict[i] = {}

    for key, data in split_dict.items():
        for num_items, elems in data.items():
            for j_idx, j in enumerate(proc_counts):
                if (int(num_items) // j) not in weak_scale_dict[key].keys():
                    weak_scale_dict[key][int(num_items) // j] = [None for i in range(len(x_axis_points))]

                weak_scale_dict[key][int(num_items) // j][j_idx] = elems[j_idx]

        weak_scale_dict[key] = dict(sorted(weak_scale_dict[key].items()))

    markers = ["o", "v", "^", "<", "8", "s", "p", "P", "*", "h", "X", "D", "d", "1", "H", "4", "x", "+"]
    for suffix, data in weak_scale_dict.items():
        plt.figure()
        marker_idx = 0
        for key, datapoint in data.items():
            plt.plot(proc_counts[:len(datapoint)], datapoint, marker=markers[marker_idx], label=key)
            marker_idx += 1
        plt.xlabel('Number of cores')
        plt.ylabel('Iteration time [ms]')
        plt.yscale('log')
        plt.xscale('log', basex = 2)
        plt.legend( prop={'size': 6})
        plt.grid(True)
        plt.savefig("{}{}.svg".format(out_prefix, suffix))


suffixes=["_p2p", "_p2p_seq_IO", "_p2p_par_IO" ,"_rma", "_rma_seq_IO", "_rma_par_IO"]

ploth_width  = 4
ploth_height = 3
legend_font_size = 6

x_axis_points = [1, 16, 32, 64, 128]
#x_axis_points = [1,2,4]
data = []
with open('run_full_mpi_2d_out.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=';')
    data = list(csvreader)

# print(data)
# exit()

data_dict = {}
headerlen = 1

for i in range(0,len(data)-headerlen):
    suffix = suffixes[i % 6]
    row = data[i + headerlen]
    datakey = str(row[4]) + suffix
    iterationtime = float(row[13])
    if not datakey in data_dict.keys():
        data_dict[datakey] = []
    data_dict[datakey].append(iterationtime)

weak_scaling(data_dict, "ppp_weak_scaling_mpi_", suffixes, x_axis_points)
for comm_type in ['p2p', 'rma']:


    # scaling
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"

        if key.find(comm_type) != -1:
            scaling_values      = data_series
            plt.plot(x_axis_points[0:len(scaling_values)], scaling_values, marker=plot_marker, label=key)

    plt.xlabel('Number of cores')
    plt.ylabel('Iteration time [ms]')
    plt.yscale('log')
    plt.xscale('log', basex = 2)
    plt.legend( prop={'size': legend_font_size})
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_scaling_mpi_" + comm_type + ".svg", bbox_inches='tight')

    # speedup
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]

            plt.plot(x_axis_points[0:len(speedup_values)], speedup_values, marker=plot_marker, label=key)

    plt.legend(prop={'size': legend_font_size})
    plt.xlabel('Number of cores')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_speedup_mpi_" + comm_type + ".svg", bbox_inches='tight')

    # efficiency
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]
            efficiency_values = []
            for i in range(0,len(speedup_values)):
                efficiency_values.append(speedup_values[i]/x_axis_points[i]*100)

            plt.plot(x_axis_points[0:len(efficiency_values)], efficiency_values, marker=plot_marker, label=key)


    plt.xlabel('Number of cores')
    plt.ylabel('Efficiency [%]')
    plt.legend(prop={'size': legend_font_size})
    plt.grid(True)
#plt.show()
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_efficiency_mpi_" + comm_type + ".svg", bbox_inches='tight')



x_axis_points = [1, 2*9, 4*9, 8*9, 16*9, 32*9]
#x_axis_points = [1,2,4]
data = []
with open('run_full_hybrid_2d_out.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=';')
    data = list(csvreader)


data_dict = {}


for i in range(0,len(data)-headerlen):
    suffix = suffixes[i % 6]
    row = data[i + headerlen]
    datakey = str(row[4]) + suffix
    iterationtime = float(row[13])
    if not datakey in data_dict.keys():
        data_dict[datakey] = []
        data_dict[datakey].append(iterationtime)
    else:
        data_dict[datakey].append(iterationtime)

weak_scaling(data_dict, "ppp_weak_scaling_hybrid_2d_", suffixes, x_axis_points)
for comm_type in ['p2p', 'rma']:

    # scaling
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"

        if key.find(comm_type) != -1:
            scaling_values      = data_series
            plt.plot(x_axis_points[0:len(scaling_values)], scaling_values, marker=plot_marker, label=key)

    plt.xlabel('Number of cores')
    plt.ylabel('Iteration time [ms]')
    plt.yscale('log')
    plt.xscale('log', basex = 2)
    plt.legend(prop={'size': legend_font_size})
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_scaling_hybrid_" + comm_type + ".svg", bbox_inches='tight')

    # speedup
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]

            plt.plot(x_axis_points[0:len(speedup_values)], speedup_values, marker=plot_marker, label=key)

    plt.legend(prop={'size': legend_font_size})
    plt.xlabel('Number of cores')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_speedup_hybrid_" + comm_type + ".svg", bbox_inches='tight')

    # efficiency
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]
            efficiency_values = []
            for i in range(0,len(speedup_values)):
                efficiency_values.append(speedup_values[i]/x_axis_points[i]*100)

            plt.plot(x_axis_points[0:len(efficiency_values)], efficiency_values, marker=plot_marker, label=key)


    plt.xlabel('Number of cores')
    plt.ylabel('Efficiency [%]')
    plt.legend(prop={'size': legend_font_size})
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_efficiency_hybrid_" + comm_type + ".svg", bbox_inches='tight')


x_axis_points = [1, 2*9, 4*9, 8*9, 16*9, 32*9]
#x_axis_points = [1,2,4]
data = []
with open('run_full_hybrid_1d_out.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=';')
    data = list(csvreader)


data_dict = {}

for i in range(0,len(data)-headerlen):
    suffix = suffixes[i % 6]
    row = data[i + headerlen]
    datakey = str(row[4]) + suffix
    iterationtime = float(row[13])
    if not datakey in data_dict.keys():
        data_dict[datakey] = []
        data_dict[datakey].append(iterationtime)
    else:
        data_dict[datakey].append(iterationtime)

weak_scaling(data_dict, "ppp_weak_scaling_hybrid_1d_", suffixes, x_axis_points)
for comm_type in ['p2p', 'rma']:

    # scaling
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"

        if key.find(comm_type) != -1:
            scaling_values      = data_series
            plt.plot(x_axis_points[0:len(scaling_values)], scaling_values, marker=plot_marker, label=key)

    plt.xlabel('Number of cores')
    plt.ylabel('Iteration time [ms]')
    plt.yscale('log')
    plt.xscale('log', basex = 2)
    plt.legend(prop={'size': legend_font_size})
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_scaling_hybrid_1D_" + comm_type + ".svg", bbox_inches='tight')

    # speedup
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]

            plt.plot(x_axis_points[0:len(speedup_values)], speedup_values, marker=plot_marker, label=key)

    plt.legend(prop={'size': legend_font_size})
    plt.xlabel('Number of cores')
    plt.ylabel('Speedup')
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_speedup_hybrid_1D_" + comm_type + ".svg", bbox_inches='tight')

    # efficiency
    plt.figure(figsize=[ploth_width, ploth_height])
    for key, data_series in data_dict.items():
        plot_marker = ""
        if key.find('256')  != -1:
            plot_marker = "o"
        if key.find('512')  != -1:
            plot_marker = "v"
        if key.find('1024') != -1:
            plot_marker = "s"
        if key.find('2048') != -1:
            plot_marker = "X"
        if key.find('4096') != -1:
            plot_marker = "d"
        if key.find(comm_type) != -1:
            scaling_values      = data_series
            speedup_values      = [data_series[0]/data_point for data_point in data_series]
            efficiency_values = []
            for i in range(0,len(speedup_values)):
                efficiency_values.append(speedup_values[i]/x_axis_points[i]*100)

            plt.plot(x_axis_points[0:len(efficiency_values)], efficiency_values, marker=plot_marker, label=key)


    plt.xlabel('Number of cores')
    plt.ylabel('Efficiency [%]')
    plt.legend(prop={'size': legend_font_size})
    plt.grid(True)
    plt.gca().set_position([0, 0, 1, 1])
    plt.savefig("ppp_efficiency_hybrid_1D_" + comm_type + ".svg", bbox_inches='tight')
