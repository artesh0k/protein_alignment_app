from itertools import islice

from django.http import JsonResponse
from django.shortcuts import render
import numpy as np
from Bio.Align import substitution_matrices
import os
import time

from django.template.loader import render_to_string


# Create your views here.

def substitution_matrices_info(request):
    pam250 = substitution_matrices.load('PAM250')
    pam70 = substitution_matrices.load('PAM70')
    pam30 = substitution_matrices.load('PAM30')
    blosum45 = substitution_matrices.load('BLOSUM45')
    blosum50 = substitution_matrices.load('BLOSUM50')
    blosum62 = substitution_matrices.load('BLOSUM62')
    blosum80 = substitution_matrices.load('BLOSUM80')
    blosum90 = substitution_matrices.load('BLOSUM90')

    sub_matrices = {
        "pam250": pam250,
        "pam70": pam70,
        "pam30": pam30,
        "blosum45": blosum45,
        "blosum50": blosum50,
        "blosum62": blosum62,
        "blosum80": blosum80,
        "blosum90": blosum90,
    }
    return render(request, 'sub_mat.html', {'sub_matrices': sub_matrices})


def get_alignment_results(request):
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sub_mat_name = request.POST.get('sub_mat_name', '')

        results = {}
        base_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(base_dir, '../insulin_proteins-40.fasta')

        start_time = time.time()

        for key, value in fasta_database_reader(file_path).items():
            output_1, output_2 = smith_waterman(sequence_1=sequence_1, sequence_2=value, sub_mat_name=sub_mat_name)
            similarity = calculate_similarity(output_1, output_2, sub_mat_name)
            score = calculate_score(output_1, output_2, sub_mat_name)
            identity = calculate_identity(output_1, output_2)
            gaps = calculate_gaps(output_1, output_2)
            highlighted_1, highlighted_2 = highlight_matches(output_1, output_2, sub_mat_name)
            results[output_2] = (similarity, identity, highlighted_1, highlighted_2, key, gaps, len(output_1), score)

        end_time = time.time()
        time_results = round(end_time - start_time, 5)

        results = dict(sorted(results.items(), key=lambda item: item[1][0], reverse=True))

        results = dict(islice(results.items(), 100))  # first 100 sequences results

        html = render_to_string("alignment_results.html",
                                {'results': results, 'sequence_1': sequence_1, 'time': time_results})
        return JsonResponse({'html': html})

    return JsonResponse({'error': 'Invalid request method'}, status=400)


def home(request):
    return render(request, "home.html")


def highlight_matches(seq1, seq2, sub_mat_name):
    sub_mat = substitution_matrices.load(sub_mat_name)
    highlighted_seq1 = ""
    highlighted_seq2 = ""

    for a, b in zip(seq1, seq2):
        if a == b:
            highlighted_seq1 += f"<span style='color: green; '>{a}</span>"
            highlighted_seq2 += f"<span style='color: green; '>{b}</span>"
        elif a != '-' and b != '-' and sub_mat[(a, b)] >= 0:
            highlighted_seq1 += f"<span style='color: blue; '>{a}</span>"
            highlighted_seq2 += f"<span style='color: blue; ;'>{b}</span>"
        else:
            highlighted_seq1 += f"<span style='color: red; font-weight: bold;'>{a}</span>"
            highlighted_seq2 += f"<span style='color: red; font-weight: bold;'>{b}</span>"

    return highlighted_seq1, highlighted_seq2


def fasta_line_reader(sequence_file):
    lines = open(sequence_file).readlines()
    sequence_name_row = lines[0][1:]
    sequence = ''.join(line.strip() for line in lines[1:])
    return sequence_name_row.strip(), sequence


def fasta_database_reader(sequence_file):
    sequences = {}
    with open(sequence_file, 'r') as file:
        current_name = None
        current_sequence = []

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = ''.join(current_sequence)
                current_name = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_name:
            sequences[current_name] = ''.join(current_sequence)

    return sequences


def calculate_identity(sequence_1, sequence_2):
    aligned_positions = 0
    identical_positions = 0
    for res1, res2 in zip(sequence_1, sequence_2):
        aligned_positions += 1
        if res1 == res2:
            identical_positions += 1

    identity = (identical_positions / aligned_positions) * 100 if aligned_positions > 0 else 0
    return round(identity, 2)


def calculate_similarity(sequence_1, sequence_2, sub_mat_name):
    sub_mat = substitution_matrices.load(sub_mat_name)
    aligned_positions = 0
    similar_positions = 0

    for res1, res2 in zip(sequence_1, sequence_2):
        aligned_positions += 1

        if res1 == '-' or res2 == '-':
            continue

        if sub_mat[(res1, res2)] >= 0:
            similar_positions += 1

    similarity = (similar_positions / aligned_positions) * 100 if aligned_positions > 0 else 0
    return round(similarity, 2)


def calculate_gaps(sequence_1, sequence_2):
    gaps = 0
    for res1, res2 in zip(sequence_1, sequence_2):
        if res1 == '-':
            gaps = gaps + 1
        if res2 == '-':
            gaps = gaps + 1

    return gaps


def calculate_score(sequence_1, sequence_2, sub_mat_name):
    sub_mat = substitution_matrices.load(sub_mat_name)
    gap_penalty = -2
    score = 0
    for i in range(len(sequence_1) - 1):
        if (sequence_1[i] == '-' and sequence_2[i] != '-') or (sequence_1[i] != '-' and sequence_2[i] == '-'):
            score = score + gap_penalty
        else:
            score = score + sub_mat[(sequence_1[i], sequence_2[i])]

    return score


def smith_waterman(sequence_1, sequence_2, sub_mat_name):
    # Creat Matrices
    main_matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))
    match_checker_matrix = np.zeros((len(sequence_1), len(sequence_2)))

    # Providing the scores for match ,mismatch and gap
    sub_mat = substitution_matrices.load(sub_mat_name)
    gap_penalty = -2

    # Fill the match checker matrix according to match or mismatch
    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):
            match_checker_matrix[i][j] = sub_mat[(sequence_1[i], sequence_2[j])]

    # Filling up the matrix using Smith_Waterman algorithm
    # STEP 1 : Initialisation
    # fill the first row and column with 0

    # STEP 2 : Matrix Filling (using a linear gap penalty)
    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            left_up_value = main_matrix[i - 1][j - 1] + match_checker_matrix[i - 1][j - 1]
            up_value = main_matrix[i - 1][j] + gap_penalty  # a linear gap penalty
            left_value = main_matrix[i][j - 1] + gap_penalty  # a linear gap penalty

            left_up_value = 0 if left_up_value < 0 else left_up_value
            left_value = 0 if left_value < 0 else left_value
            up_value = 0 if up_value < 0 else up_value

            main_matrix[i][j] = max(left_up_value, left_value, up_value)

    # STEP 3 : Traceback
    max_value = max(max(row) for row in main_matrix)
    positions = [(i, j) for i, row in enumerate(main_matrix) for j, value in enumerate(row) if value == max_value]

    alignments = []
    for ti, tj in positions:
        actual_value = max_value
        temp_aligned_1 = ""
        temp_aligned_2 = ""

        while actual_value != 0:
            if main_matrix[ti][tj] == main_matrix[ti - 1][tj - 1] + match_checker_matrix[ti - 1][tj - 1]:
                temp_aligned_1 = sequence_1[ti - 1] + temp_aligned_1
                temp_aligned_2 = sequence_2[tj - 1] + temp_aligned_2

                ti = ti - 1
                tj = tj - 1
                actual_value = main_matrix[ti][tj]

            elif ti > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj] + gap_penalty:
                temp_aligned_1 = sequence_1[ti - 1] + temp_aligned_1
                temp_aligned_2 = "-" + temp_aligned_2

                ti = ti - 1
                actual_value = main_matrix[ti][tj]

            else:
                temp_aligned_1 = "-" + temp_aligned_1
                temp_aligned_2 = sequence_2[tj - 1] + temp_aligned_2

                tj = tj - 1
                actual_value = main_matrix[ti][tj]

        score = calculate_score(temp_aligned_1, temp_aligned_2, sub_mat_name)
        alignments.append((temp_aligned_1, temp_aligned_2, score))

    # choose the best score in different alignment
    best_alignment = max(alignments, key=lambda x: x[2])

    best_aligned_1 = best_alignment[0]
    best_aligned_2 = best_alignment[1]

    return best_aligned_1, best_aligned_2
