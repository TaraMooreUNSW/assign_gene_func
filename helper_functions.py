from Bio import Align
from Bio.Align import substitution_matrices

GAP_PENALTY = 8

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
    aligner = Align.PairwiseAligner()
    blosum_matrix = substitution_matrices.load("BLOSUM62")
    aligner.substitution_matrix = blosum_matrix

    rows = len(seq2) + 1
    cols = len(seq1) + 1

    # initialise matrix
    alignment_matrix = [[{} for _ in range(cols)] for _ in range(rows)]
    alignment_matrix[0][0] = 0

    # set matrix initial row
    for i in range(rows):
        alignment_matrix[i][0] = {"score": -i*GAP_PENALTY,
                                  "prev_cell_row": "",
                                  "prev_cell_col": ""}

    # set matrix initial col
    for j in range(cols):
        alignment_matrix[0][j] = {"score": -j*GAP_PENALTY,
                                  "prev_cell_row": "",
                                  "prev_cell_col": ""}

    # recurrence
    for i in range(1, rows):
        for j in range(1, cols):
            # recurr_one = alignment_matrix[i-1][j-1]["score"] + scoring_function(seq1[j-1], seq2[i-1])
            recurr_one = alignment_matrix[i-1][j-1]["score"] + aligner.score(seq1[j-1], seq2[i-1])
            recurr_two = alignment_matrix[i-1][j]["score"] - GAP_PENALTY
            recurr_three = alignment_matrix[i][j-1]["score"] - GAP_PENALTY

            scores = [recurr_one, recurr_two, recurr_three]
            max_score = max(scores)
            max_score_index = scores.index(max_score)

            alignment_matrix[i][j]["score"] = max_score

            # traceback pointer
            if max_score_index == 0:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j - 1
            elif max_score_index == 1:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j
            elif max_score_index == 2:
                alignment_matrix[i][j]["prev_cell_row"] = i
                alignment_matrix[i][j]["prev_cell_col"] = j - 1


    # TO DO - add backtracking from bottom right
    print(alignment_matrix[rows-1][cols-1])

    # print(alignment_matrix)
    # traceback_i = len(seq2) + 1
    # traceback_j = len(seq1) + 1
    # traceback_k = 0
    # traceback_id = 0
    # while


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
    global_alignment('HEAGAWGHEE', 'PAWHEAE',lambda x, y: [-1, 1][x == y])

if __name__ == '__main__':
    main()
