from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plot
import pylab
from tabulate import tabulate


ENCODINGS = [
    ("Sanger Phred+33", "!", "I"),
    ("Solexa Solexa+64", ";", "h"),
    ("Illumina 1.3+ Phred+64", "@", "h"),
    ("Illumina 1.5+ Phred+64", "B", "i"),
    ("Illumina 1.8+ Phred+33", "!", "J")
]


class HandRecord:

    def __init__(self):
        self.id = ""
        self.description = ""
        self.seq = ""
        self.rating = ""
        return


def hand_read_fastq(file_name):
    the_return = []
    state_id = 0
    state_seq = 1
    state_plus = 2
    state_rating = 3
    state = 0
    temp_record = HandRecord()
    with open(file_name, "r") as file:
        for line in file:
            # print(f"line: ''{line}'' len:{len(line)} true? {True if line else False}")
            line = line.replace("\r", "").replace("\n", "")

            if state == state_id:
                if line[0] == '@':
                    temp_record.id, temp_record.description = line[1:].split(" ", 1)
                    state = state_seq
                    continue
            elif state == state_seq:
                temp_record.seq = line
                state = state_plus
                continue
            elif state == state_plus:
                if line[0] == '+':
                    state = state_rating
                    continue
            elif state == state_rating:
                temp_record.rating = line
                the_return.append(temp_record)
                temp_record = HandRecord()
                state = state_id
    return the_return


def identify_encoding(rating: str, options: [(str, str, str)]):
    # example options = [("Sanger Phred+33", "!", "I")]
    low = min(rating)
    high = max(rating)
    to_return = []
    for option in options:
        if low < option[1]:
            continue
        if high > option[2]:
            continue
        to_return.append(option)
    return to_return


def calc_distribution(source: str, what: [str]):
    total = len(source)
    the_sum = 0
    for i in range(0, len(what)):
        the_sum += source.count(what[i])
    return the_sum / total


def prepare_plot(feed_count, title):
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


def find_pikas(data: [float], min_len_in_pikas=5, max_amount_of_pikas=5):
    data = [x for x in data]  # copy
    data.sort()
    to_return = []
    pikas = []
    target_len = 9999999
    max_diff_in_pikas = max(data) / 2 / max_amount_of_pikas

    append = False
    carry = False
    for i_rev in range(0, len(data))[-1::-1]:
        pikas.append(data[i_rev])
        if len(pikas) < min_len_in_pikas:
            continue

        if pikas[0] - pikas[-1] > max_diff_in_pikas and len(pikas) - 1 >= min_len_in_pikas:
            pikas = pikas[0:-1]
            append = True
            carry = True
        elif len(pikas) >= target_len and data[i_rev] != pikas[-2]:
            append = True
            i_rev -= 1

        if append:
            to_return.append(pikas)
            target_len = len(to_return[-1]) * 5
            append = False
            if i_rev > -1 and len(to_return) < max_amount_of_pikas:
                pikas = [data[i_rev]] if carry else []
                carry = False
            else:
                return to_return

    to_return.append(pikas)
    return to_return


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


def peaks_to_seqs(peaks: [float], the_xs: [float], seqs: [HandRecord]):
    to_return = []
    the_xs = [x for x in the_xs]  # copy
    for peak in peaks:
        the_seq = []
        for dist in peak:
            the_i = the_xs.index(dist)
            the_xs[the_i] = -1999.99991
            the_seq.append(seqs[the_i])
        to_return.append(the_seq)
    return to_return


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


def main():
    global ENCODINGS
    file_name = "data\\reads_for_analysis.fastq"
    seqs = hand_read_fastq(file_name)

    # find final encoding
    counter = [0] * len(ENCODINGS)
    for seq in seqs:
        encoding = identify_encoding(seq.rating, ENCODINGS)
        for enco in encoding:
            counter[ENCODINGS.index(enco)] += 1
    encoding_index = counter.index(max(counter))  # if there are multi answers, we pick the top one (leftist)
    encoding = ENCODINGS[encoding_index]
    print("Quality rating encoding: " + encoding[0])

    # calculate cg - s and vizualize
    the_xs = [100 * calc_distribution(sq.seq, ["C", "G", "S"]) for sq in seqs]
    the_ys = [read for read in range(1, len(seqs) + 1)]
    prepare_plot(len(the_xs), "C/G nukleotidų pasiskirstymas read'uose")
    plot.scatter(the_xs, the_ys, 2, "r")
    plot.show()

    # find pika'chu
    pikas = find_pikas(the_xs, 5, 4)
    prepare_plot(len(the_xs), "C/G nukleotidų pasiskirstymas read'uose ir rasti pikai")
    plot.scatter(the_xs, the_ys, 2, "r")
    plot_peaks(pikas, the_ys)
    plot.show()

    # perform blast'a search for each 5 seq of peaks
    print("Preparing data for BLAST search")
    the_seqs = peaks_to_seqs(pikas, the_xs, seqs)
    matches = []
    c = 0
    for peak in the_seqs:
        for i in range(0, min(5, len(peak))):
            c += 1
            print("performing BLAST search " + f"({c}/{len(pikas*5)})")
            print(f"seq: {peak[i].seq}")
            result_handle = NCBIWWW.qblast(program="blastn",
                                           database="nt",
                                           sequence=f";the test sample\r\n{peak[i].seq}",
                                           gapcosts="5 2",
                                           hitlist_size=1,
                                           nucl_penalty=-3,
                                           nucl_reward=2,
                                           entrez_query='"bacteria"[organism]')
            result_xml = result_handle.read()
            print(result_xml)
            result_title = ""
            try:
                result_title = ET.fromstring(result_xml)\
                    .find("BlastOutput_iterations")\
                    .find("Iteration")\
                    .find("Iteration_hits")\
                    .find("Hit")\
                    .find("Hit_def")\
                    .text
            except AttributeError:
                result_title = "???"
            print(f"best match: {result_title}")
            matches.append((peak[i].id, result_title))
            print_data_in_table(matches, left_right_headers=["Read'o id", "Rastas organizmas"])
    print("BLAST search done")

    # print search results
    print_data_in_table(matches, left_right_headers=["Read'o id", "Rastas organizmas"])

    # for seq_record in SeqIO.parse(file_name, "fastq"):
    #     # print(seq_record.letter_annotations["phred_quality"])
    #     print(seq_record.id)
    #     return


if __name__ == '__main__':
    main()
else:
    print("not main. I am quitting")
