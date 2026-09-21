from Bio import Align
from Bio.Align import substitution_matrices
# matrix = substitution_matrices.load("BLOSUM62")

GAP_PENALTY = 1

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    # raise NotImplementedError()

    rows = len(seq2) + 1
    cols = len(seq1) + 1

    # initialise matrix
    alignment_matrix = [["" for _ in range(cols)] for _ in range(rows)]
    alignment_matrix[0][0] = 0

    # set matrix initial row
    for i in range(rows):
        alignment_matrix[i][0] = -i*GAP_PENALTY


    # set matrix initial col
    for j in range(cols):
        alignment_matrix[0][j] = -j*GAP_PENALTY

    # recurrence
    for i in range(1, rows):
        for j in range(1, cols):
            recurr_one = alignment_matrix[i-1][j-1] + scoring_function(seq1[j-1], seq2[i-1])
            recurr_two = alignment_matrix[i-1][j] - GAP_PENALTY
            recurr_three = alignment_matrix[i][j-1] - GAP_PENALTY

            alignment_matrix[i][j] = max(recurr_one, recurr_two, recurr_three)

    # TO DO - add backtracking



def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

def main():
    global_alignment('abracadabra', 'dabarakadara',lambda x, y: [-1, 1][x == y])

if __name__ == '__main__':
    main()
