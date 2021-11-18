from main import *

test_counter = 0


def print_if_false(title, boolean):
    global test_counter
    test_counter += 1
    if not boolean:
        print("----- TEST. " + title + " failed.")
        return


if __name__ == "__main__":
    data = hand_read_fastq("data\\reads_for_analysis.fastq")
    print_if_false("number of seqs read", len(data) == 25749)
    data = ''

    print(f"Number of successful tests: {test_counter}")
