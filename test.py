from main import *

test_counter = 0


def print_if_not_equal(title, left, right):
    global test_counter
    test_counter += 1
    if left != right:
        print("----- TEST. " + title + " failed: ")
        print(f"----------- {left}")
        print(f"-----------  !=")
        print(f"----------- {right}")
        print("-----------------------------------------")
        test_counter -= 1
        return


if __name__ == "__main__":
    data = hand_read_fastq("data\\reads_for_analysis.fastq")
    out = len(data)
    answer = 25749
    print_if_not_equal("number of seqs read", out, answer)

    data = "!$$H"
    param = [("Sanger Phred+33", "!", "I"), ("Solexa Solexa+64", ";", "h")]
    out = identify_encoding(data, param)
    answer = [("Sanger Phred+33", "!", "I")]
    print_if_not_equal("solexa +33 and +64", out, answer)

    data = "!$$H"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J")]
    out = identify_encoding(data, param)
    answer = param
    print_if_not_equal("solexa +33 illumina 1.8+ +33", out, answer)

    data = "=>?"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J"), ("Illumina 1.5+ Phred+64", "B", "i")]
    out = identify_encoding(data, param)
    answer = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J")]
    print_if_not_equal("solexa +33 and illumina 1.8+ +33 out of three", out, answer)

    data = "g"
    param = [("Sanger Phred+33", "!", "I"), ("Illumina 1.8+ Phred+33", "!", "J"), ("Illumina 1.5+ Phred+64", "B", "i")]
    out = identify_encoding(data, param)
    answer = [("Illumina 1.5+ Phred+64", "B", "i")]
    print_if_not_equal("solexa +33 and illumina 1.8+ +33 out of three", out, answer)

    data = "GACTGACT"
    param = ["G", "C"]
    out = calc_distribution(data, param)
    answer = 0.5
    print_if_not_equal("50 dist", out, answer)

    data = "GSCTGSCT"
    param = ["G", "C", "S"]
    out = calc_distribution(data, param)
    answer = 3/4
    print_if_not_equal("3/4 dist", out, answer)

    data = "GAAAAA"
    param = ["C"]
    out = calc_distribution(data, param)
    answer = 0
    print_if_not_equal("0 dist", out, answer)

    data = [90, 80, 80], [60]
    param = [[60, 70, 80, 80, 90],
             [HandRecord(), HandRecord(), HandRecord(), HandRecord(), HandRecord()]
             ]
    out = peaks_to_seqs(data, *param)
    answer = [[param[1][4], param[1][2], param[1][3]], [param[1][0]]]
    print_if_not_equal("peaks to seqs", out, answer)

    print(f"Number of successful tests: {test_counter}")
