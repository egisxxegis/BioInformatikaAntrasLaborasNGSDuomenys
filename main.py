from Bio import SeqIO
import copy


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
                temp_record = HandRecord
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

def main():
    file_name = "data\\reads_for_analysis.fastq"
    seqs = hand_read_fastq(file_name)

    for seq_record in SeqIO.parse(file_name, "fastq"):
        # print(seq_record.letter_annotations["phred_quality"])
        print(seq_record.id)
        return


if __name__ == '__main__':
    main()
else:
    print("not main. I am quitting")

