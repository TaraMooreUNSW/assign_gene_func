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
    rows = len(seq2) + 1
    cols = len(seq1) + 1

    # initialise matrix
    alignment_matrix = [[{} for _ in range(cols)] for _ in range(rows)]
    alignment_matrix[0][0] = 0

    # set matrix initial row
    for i in range(rows):
        alignment_matrix[i][0] = {"score": -i*GAP_PENALTY,
                                  "prev_cell_row": i - 1,
                                  "prev_cell_col": 0}

    # set matrix initial col
    for j in range(cols):
        alignment_matrix[0][j] = {"score": -j*GAP_PENALTY,
                                  "prev_cell_row": 0,
                                  "prev_cell_col": j - 1}

    # recurrence
    for i in range(1, rows):
        for j in range(1, cols):
            recurr_one = (alignment_matrix[i-1][j-1]["score"]
                        + scoring_function(seq1[j-1], seq2[i-1]))
            recurr_two = alignment_matrix[i-1][j]["score"] - GAP_PENALTY
            recurr_three = alignment_matrix[i][j-1]["score"] - GAP_PENALTY

            scores = [recurr_one, recurr_two, recurr_three]
            max_score = max(scores)
            max_score_index = scores.index(max_score)

            alignment_matrix[i][j]["score"] = max_score

            # traceback recording
            if max_score_index == 0:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j - 1
            elif max_score_index == 1:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j
            elif max_score_index == 2:
                alignment_matrix[i][j]["prev_cell_row"] = i
                alignment_matrix[i][j]["prev_cell_col"] = j - 1

    # Traceback from bottom right
    seq1_aligned = list(seq1)
    seq2_aligned = list(seq2)

    t_i = len(seq2)
    t_j = len(seq1)
    while t_i > 0 or t_j > 0:
        i_backpoint = alignment_matrix[t_i][t_j]["prev_cell_row"]
        j_backpoint = alignment_matrix[t_i][t_j]["prev_cell_col"]

        if i_backpoint == t_i:
            t_j = j_backpoint           # gap in seq2
            seq2_aligned.insert(t_i, '-')
        elif j_backpoint == t_j:
            t_i = i_backpoint           # gap in seq1
            seq1_aligned.insert(t_j, '-')
        else:
            t_i = i_backpoint
            t_j = j_backpoint

    seq1_final = "".join(seq1_aligned)
    seq2_final = "".join(seq2_aligned)

    return(seq1_final, seq2_final, f"{alignment_matrix[rows - 1][cols - 1]["score"]:.1f}")


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
    rows = len(seq2) + 1
    cols = len(seq1) + 1

    # initialise matrix
    alignment_matrix = [[{} for _ in range(cols)] for _ in range(rows)]
    alignment_matrix[0][0] = 0

    # set matrix initial row
    for i in range(rows):
        alignment_matrix[i][0] = {"score": 0}

    # set matrix initial col
    for j in range(cols):
        alignment_matrix[0][j] = {"score": 0}

    max_cell_score = 0
    max_cell_row = 0
    max_cell_col = 0

    for i in range(1, rows):
        for j in range(1, cols):
            recurr_one = (alignment_matrix[i-1][j-1]["score"]
                        + scoring_function(seq1[j-1], seq2[i-1]))
            recurr_two = alignment_matrix[i-1][j]["score"] - GAP_PENALTY
            recurr_three = alignment_matrix[i][j-1]["score"] - GAP_PENALTY

            scores = [recurr_one, recurr_two, recurr_three, 0]
            max_score = max(scores)
            max_score_index = scores.index(max_score)

            alignment_matrix[i][j]["score"] = max_score
            if max_score >= max_cell_score:
                max_cell_score = max_score
                max_cell_row = i
                max_cell_col = j

            # traceback recording
            if max_score_index == 0 and max_score > 0:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j - 1
            elif max_score_index == 1 and max_score > 0:
                alignment_matrix[i][j]["prev_cell_row"] = i - 1
                alignment_matrix[i][j]["prev_cell_col"] = j
            elif max_score_index == 2 and max_score > 0:
                alignment_matrix[i][j]["prev_cell_row"] = i
                alignment_matrix[i][j]["prev_cell_col"] = j - 1

    # Traceback from cell with highest score
    seq1_aligned = []
    seq2_aligned = []

    t_i = max_cell_row
    t_j = max_cell_col
    while alignment_matrix[t_i][t_j]["score"] > 0:
        i_backpoint = alignment_matrix[t_i][t_j]["prev_cell_row"]
        j_backpoint = alignment_matrix[t_i][t_j]["prev_cell_col"]

        if i_backpoint == t_i:
            t_j = j_backpoint           # gap in seq2
            seq1_aligned.append(seq1[t_j])
            seq2_aligned.append('-')
        elif j_backpoint == t_j:
            t_i = i_backpoint           # gap in seq1
            seq2_aligned.append(seq2[t_i])
            seq1_aligned.append('-')
        else:
            t_i = i_backpoint
            t_j = j_backpoint
            seq1_aligned.append(seq1[t_j])
            seq2_aligned.append(seq2[t_i])

    seq1_aligned.reverse()
    seq2_aligned.reverse()
    seq1_final = "".join(seq1_aligned)
    seq2_final = "".join(seq2_aligned)

    return(seq1_final, seq2_final, f"{max_cell_score:.1f}")


## This is an example scoring function, you should implement a version which uses a scoring matrix
def scoring_function_simple(aa_i,aa_j):
    """
    Scoring function using simple substitution matrix

    Parameters
        ----------
        aa_i: str
            First amino acid.
        aa_j: str
            Second amino acid.

        Returns
        -------
        int
            Score of match between amino acids.

        Examples
        --------
        >>> scoring_function_simple('W', 'W')
        (1)
    """
    score = [-1, 1][aa_i == aa_j]
    return (score)

def scoring_BLOSUM62(aa_i, aa_j):
    """
    Scoring function using the BLOSUM62 substitution matrix

    Parameters
        ----------
        aa_i: str
            First amino acid.
        aa_j: str
            Second amino acid.

        Returns
        -------
        float
            Score of match between amino acids.

        Examples
        --------
        >>> scoring_BLOSUM('W', 'W')
        (11.0)
    """
    aligner = Align.PairwiseAligner()
    blosum_matrix = substitution_matrices.load("BLOSUM50")
    aligner.substitution_matrix = blosum_matrix

    score = aligner.score(aa_i, aa_j)
    return score

# def main():
    # global_alignment('HEAGAWGHEE', 'PAWHEAE', scoring_BLOSUM62)
    # global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    # local_alignment("HEAGAWGHEE", "PAWHEAE", scoring_BLOSUM62)
    # local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])

# if __name__ == '__main__':
#     main()
