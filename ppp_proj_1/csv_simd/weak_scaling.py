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

    markers = ["o", "v", "s", "X", "d", "+", "x", "H", "*", "p", "1", ".", "^", "<", ">", "2", "3", "4"]
    for suffix, data in weak_scale_dict.items():
        plt.figure()
        marker_idx = 0
        for key, datapoint in data.items():
            plt.plot(proc_counts[:len(datapoint)], datapoint, marker=markers[marker_idx], label=key)
            marker_idx += 1
        plt.xlabel('Number of cores')
        plt.ylabel('Iteration time [ms]')
        plt.yscale('log')
        plt.xscale('log', base = 2)
        plt.legend( prop={'size': 6})
        plt.grid(True)
        plt.savefig("{}{}.svg".format(out_prefix, suffix))
