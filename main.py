from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as elTree
import matplotlib.pyplot as plot

from fileHandler import hand_read_fastq
from representationTools import prepare_plot, plot_peaks, print_data_in_table
from sequenceTools import identify_encoding, calc_distribution, find_pikas, peaks_to_seqs

ENCODINGS = [
    ("Sanger Phred+33", "!", "I"),
    ("Solexa Solexa+64", ";", "h"),
    ("Illumina 1.3+ Phred+64", "@", "h"),
    ("Illumina 1.5+ Phred+64", "B", "i"),
    ("Illumina 1.8+ Phred+33", "!", "J")
]


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
                result_title = elTree.fromstring(result_xml)\
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
