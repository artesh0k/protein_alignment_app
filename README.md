### README.md

# Protein Alignment App

This project implements a local alignment algorithm for protein sequences using the Smith-Waterman algorithm. The alignment results are influenced by various substitution matrices (e.g., **PAM30**, **BLOSUM62**, **PAM250**). The application allows users to input a protein sequence, select a substitution matrix, and perform alignments against a dataset of protein sequences.

The project is built with Python and Django and utilizes **NumPy** and **Biopython** libraries for computations and bioinformatics functionalities.

---

## Table of Contents

1. [Features](#features)  
2. [Technologies Used](#technologies-used)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [API Endpoints](#api-endpoints)  
6. [Example Results](#example-results)  

---

## Features

- **Local Alignment:** Implements the Smith-Waterman algorithm for local sequence alignment.
- **Substitution Matrices:** Supports a variety of substitution matrices, including:
  - PAM250, PAM70, PAM30
  - BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90
- **Alignment Metrics:** Provides the following metrics:
  - **Score:** Alignment score based on the selected matrix.
  - **Identity:** Percentage of identical positions.
  - **Similarity:** Percentage of similar positions (based on positive scores in the matrix).
  - **Gaps:** Number of gap characters in the alignment.
- **Visualization:** Highlights matches, mismatches, and gaps in the alignment results.

---

## Technologies Used

- **Programming Language:** Python
- **Web Framework:** Django
- **Libraries:**
  - **NumPy:** For handling numerical operations and matrices.
  - **Biopython:** For bioinformatics functions and substitution matrices.

---

## Installation

### Prerequisites

1. **Python Installation:**  
   Download and install Python from the [official website](https://www.python.org/downloads/).

2. **Dependencies:**  
   Install Django, NumPy, and Biopython using `pip`:

   ```bash
   pip install django
   pip install numpy
   pip install biopython
   ```

### Project Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/artesh0k/protein_alignment_app.git
   cd protein_alignment_app
   ```

2. **Run the Server:**

   ```bash
   python manage.py runserver
   ```

3. **Access the Application:**  
   Open a web browser and go to `http://127.0.0.1:8000/`.

---

## Usage

1. **Home Page:**  
   Enter a protein sequence and select a substitution matrix.

2. **Submit Alignment:**  
   Click the **"Submit"** button to perform the local alignment against the predefined dataset.

3. **View Results:**  
   The results page will display:
   - Aligned sequences
   - Alignment metrics (score, identity, similarity, gaps)
   - Highlighted regions (matches in green, similar residues in blue, mismatches in red)

---

![image](https://github.com/user-attachments/assets/b61825d2-d637-40ff-8c6e-47b709765a3d)
![image](https://github.com/user-attachments/assets/0ade8abe-6cbb-4e24-b8ed-8f0e3dd7ab3c)

## API Endpoints

1. **Home Page (`/`):**  
   Displays the input form for entering sequences and selecting substitution matrices.

2. **Substitution Matrices Info (`/substitution_matrices_info`):**  
   Returns available substitution matrices.

3. **Get Alignment Results (`/get_alignment_results`):**  
   Accepts a POST request with the following parameters:
   - **`sequence_1`**: The input protein sequence.
   - **`sub_mat_name`**: The selected substitution matrix.

   Returns a JSON response containing the alignment results and metrics.

---

## Example Results

After submitting an input sequence, the results might look like:

- **Sequence 1:**  
  `MKTIIALSYIFCLVFA`

- **Aligned Sequence 1:**  
  `MK-TIIALSYIFCLVFA`

- **Aligned Sequence 2:**  
  `MKGSIFTL-FLFSVLFA`

- **Metrics:**  
  - **Identity:** 41.18%
  - **Similarity:** 82.35%
  - **Score:** 41.0  
  - **Gaps:** 2
  - **Length:** 17  

---
