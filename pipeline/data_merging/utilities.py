from math import floor, log10

def round_sigfigs(num, sig_figs):
    if num != 0:
        return round(num, -int(floor(log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0


def isEmpty(value):
    return value == '-' or value is None or value == [] or value == ['-'] or value == ''

