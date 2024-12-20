import numpy as np
from ..sequence.position_weight_matrix import PWM


def parse_jaspar_line(line):
    letter, rest = line.split(maxsplit=1)
    rest = rest.strip()[1:-1].split()
    counts = [float(n) for n in rest]
    return letter.strip(), counts


def read_jaspar_matrix(filename):
    f = open(filename)
    _ = f.readline()
    pwm = dict((parse_jaspar_line(line) for line in f))
    return PWM.from_dict(pwm)
    # alphabet, matrix = zip(*(parse_jaspar_line(line) for line in f))
    # 
    # alphabet = "".join(alphabet)
    # matrix = np.array(matrix, dtype="float")
    # return PWM.from_dict({char: row for
    # return alphabet, matrix


def read_csv_motif(filename: str) -> PWM:
    '''
    Read a PWM from a CSV file. The first line should be the alphabet, and the rest should be the matrix with probabilities.

    Parameters
    ----------
    filename

    Returns
    -------
    PWM

    '''
    f = open(filename)
    alphabet = f.readline().strip().split(",")
    pwm = {letter: [] for letter in alphabet}
    for line in f:
        parts = line.strip().split(",")
        for i, letter in enumerate(alphabet):
            pwm[letter].append(float(parts[i]))
    return PWM.from_dict(pwm)
