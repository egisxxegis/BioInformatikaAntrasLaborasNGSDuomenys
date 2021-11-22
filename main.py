from Bio import SeqIO
import matplotlib.pyplot as plot
import pylab


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


def prepare_plot(feed_count):
    coord_min_max_x = [0, 100]
    coord_min_max_y = [1, feed_count]

    # axis limits
    pylab.xlim(coord_min_max_x)
    pylab.ylim(coord_min_max_y)

    plot.title('C/G nukleotidų pasiskirstymas read\'uose')
    plot.xlabel('C/G nukleotidų dalis read\'e')
    plot.ylabel('Read\'ų skaičius')

    # show values on contour
    plot.rcParams.update({'font.size': 12})


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
    prepare_plot(len(the_xs))
    plot.plot(the_xs, the_ys)
    plot.show()

    # find pika

    # for seq_record in SeqIO.parse(file_name, "fastq"):
    #     # print(seq_record.letter_annotations["phred_quality"])
    #     print(seq_record.id)
    #     return


if __name__ == '__main__':
    main()
else:
    print("not main. I am quitting")
