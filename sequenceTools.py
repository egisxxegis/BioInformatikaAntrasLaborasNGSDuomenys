from HandRecord import HandRecord


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
