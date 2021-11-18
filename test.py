from main import *

test_counter = 0


def print_if_false(title, boolean):
    global test_counter
    test_counter += 1
    if not boolean:
        print("----- TEST. " + title + " failed.")
        test_counter -= 1
        return


if __name__ == "__main__":
    data = hand_read_fastq("data\\reads_for_analysis.fastq")
    out = len(data)
    answer = 25749
    print_if_false("number of seqs read", out == answer)
    data = ''

    data = "!$$H"
    param = [("Sanger Phred+33", "!", "I"), ("Solexa Solexa+64", ";", "h")]
    out = identify_encoding(data, param)
    answer = [("Sanger Phred+33", "!", "I")]
    print_if_false("solexa +33 and +64", out == answer)

    data = "!$$H"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J")]
    out = identify_encoding(data, param)
    answer = param
    print_if_false("solexa +33 illumina 1.8+ +33", out == answer)

    data = "=>?"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J"), ("Illumina 1.5+ Phred+64", "B", "i")]
    out = identify_encoding(data, param)
    answer = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J")]
    print_if_false("solexa +33 and illumina 1.8+ +33 out of three", out == answer)

    data = "g"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J"), ("Illumina 1.5+ Phred+64", "B", "i")]
    out = identify_encoding(data, param)
    answer = [("Illumina 1.5+ Phred+64", "B", "i")]
    print_if_false("solexa +33 and illumina 1.8+ +33 out of three", out == answer)

    print(f"Number of successful tests: {test_counter}")
