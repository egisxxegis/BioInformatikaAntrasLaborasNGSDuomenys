import pylab
from matplotlib import pyplot as plot
from tabulate import tabulate


def prepare_plot(feed_count, title):
    plot.close()
    plot.cla()
    plot.clf()
    coord_min_max_x = [0, 100]
    coord_min_max_y = [1, feed_count]

    # axis limits
    pylab.xlim(coord_min_max_x)
    pylab.ylim(coord_min_max_y)

    plot.title(title)
    plot.xlabel('C/G nukleotidų dalis read\'e procentais (%)')
    plot.ylabel('Read\'ų skaičius')

    # show values on contour
    plot.rcParams.update({'font.size': 12})


def plot_peaks(peaks, the_ys):
    color = "-m"
    index = 1
    dec = (max(the_ys) - min(the_ys)) / 10
    index_rev = 10 - index
    for peak in peaks:
        plot.plot([peak[0], peak[0]], [the_ys[0], the_ys[-1]], color)
        plot.annotate(f"{index}", [peak[0], dec * index_rev])
        plot.annotate(f"peak", [peak[0], dec * index_rev - dec / 2])
        index += 1
        index_rev -= 1
    return


def print_data_in_table(data, top_down_headers=None, left_right_headers=None, dimensions=1):
    if left_right_headers is None:
        left_right_headers = ["data"]
    if dimensions < 1 or dimensions > 2:
        print(f'{dimensions} dimensions are not supported in print_data_in_table.')
        return
    if dimensions == 1:
        if top_down_headers is None:
            print(tabulate(data, left_right_headers))
        else:
            new_data = []
            for index in range(len(data)):
                new_data.append([top_down_headers[index], data[index]])
            print(tabulate(new_data, ["name"] + left_right_headers))
        return
    elif dimensions == 2:
        # we substract one because top_down_headers will be added
        if len(data) % (len(left_right_headers) - 1) > 0:
            print(f'Given data cannot be converted to {len(left_right_headers)} column table.')
            return
        elif len(left_right_headers) == 1:
            print(f'1 column header can not form 2 dimensional table.')
            return
        if top_down_headers is None:
            print(f'top down headers are must for 2 dimensional table. Else use 1 dimensional')
            return
        new_data = []
        the_index = 0
        the_index_of_header = 0
        the_length = len(left_right_headers) - 1
        while the_index + the_length < len(data)+1:
            new_data.append(
                [top_down_headers[the_index_of_header]]
                + data[the_index: the_index + the_length])
            the_index += the_length
            the_index_of_header += 1
        print(tabulate(new_data, left_right_headers, tablefmt="psql"))
        return
