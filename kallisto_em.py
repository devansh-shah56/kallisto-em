"""
Kallisto EM Algorithm Implementation

This module implements the Expectation-Maximization (EM) algorithm 
used by Kallisto to quantify RNA isoforms as described in Pachter (2011).
"""

import numpy as np
import pandas as pd

def e_step(transcript_abundance_vector, assignment_matrix):
    """
    Given transcript abundance vector and assignment matrix,
    return the posterior estimates for each read.

    Parameters
    ----------
    transcript_abundance_vector : np.array
        A vector of initial transcript estimates
    assignment_matrix : np.array
        A matrix of read to transcript assignment matrix

    Returns
    -------
    posterior : np.array
        A matrix of posterior estimates for each read
    """
    # Convert to numpy arrays for safe operations
    A = np.asarray(assignment_matrix, dtype=float)  # assignment matrix
    pi = np.asarray(transcript_abundance_vector, dtype=float)  # abundance vector

    # Calculate weighted reads for each isoform
    weighted = A * pi[:, None]  # Multiply assignment by abundance
    col_sums = weighted.sum(axis=0, keepdims=True)  # Total weighted reads per sample

    # Calculate posterior probabilities with safe division
    with np.errstate(divide='ignore', invalid='ignore'):
        posterior = np.divide(weighted, col_sums, out=np.zeros_like(weighted), where=col_sums != 0)

    return posterior


def m_step(posterior, assignment_matrix):
    """
    Given posterior estimates and assignment matrix,
    return the new transcript abundance estimates.

    Parameters
    ----------
    posterior : np.array
        A matrix of posterior estimates for each read
    assignment_matrix : np.array
        A matrix of read to transcript assignment matrix

    Returns
    -------
    new_transcript_abundance : np.array
        A vector of new transcript estimates
    """
    # Sum the posterior probabilities for each transcript
    expected_counts = np.asarray(posterior, dtype=float).sum(axis=1)
    total = expected_counts.sum()

    if total > 0:
        new_transcript_abundance = expected_counts / total
    else:
        # Fallback: uniform distribution if no effective assignments
        T = expected_counts.shape[0]
        new_transcript_abundance = np.full(T, 1.0 / T)

    return new_transcript_abundance


def run_em_algorithm(assignment_matrix, initial_abundance=None, tol=1e-4, max_iters=10000):
    """
    Run the complete EM algorithm until convergence.

    Parameters
    ----------
    assignment_matrix : np.array
        A matrix of read to transcript assignment
    initial_abundance : np.array, optional
        Initial transcript abundance estimates. If None, uses uniform distribution.
    tol : float, optional
        Convergence tolerance (default: 1e-4)
    max_iters : int, optional
        Maximum number of iterations (default: 10000)

    Returns
    -------
    final_abundance : np.array
        Final transcript abundance estimates
    n_iterations : int
        Number of iterations until convergence
    """
    # Initialize abundance vector
    if initial_abundance is None:
        n_transcripts = assignment_matrix.shape[0]
        alpha = np.full(n_transcripts, 1.0 / n_transcripts)
    else:
        alpha = np.asarray(initial_abundance, dtype=float).copy()

    # Run EM algorithm
    for iteration in range(max_iters):
        # E-step: calculate posterior probabilities
        posterior = e_step(alpha, assignment_matrix)

        # M-step: update abundance estimates
        updated = m_step(posterior, assignment_matrix)

        # Check for convergence
        if np.max(np.abs(updated - alpha)) < tol:
            alpha = updated
            return alpha, iteration + 1

        alpha = updated

    print(f"Warning: EM algorithm did not converge after {max_iters} iterations")
    return alpha, max_iters


def create_example_data():
    """
    Create example data matching the assignment from the notebook.

    Returns
    -------
    assignment_matrix : np.array
        Example read-to-isoform assignment matrix
    """
    # Read to isoform assignment matrix from the notebook
    assignment_matrix = np.array([
        [1, 0, 1, 1, 1],  # Red isoform
        [1, 1, 0, 0, 1],  # Green isoform  
        [1, 1, 1, 0, 0]   # Blue isoform
    ])

    return assignment_matrix


def main():
    """
    Main function demonstrating the EM algorithm on example data.
    """
    print("Kallisto EM Algorithm for Isoform Quantification")
    print("=" * 50)

    # Create example data
    assignment_matrix = create_example_data()
    print(f"Assignment Matrix shape: {assignment_matrix.shape}")
    print("Assignment Matrix:")
    print(assignment_matrix)
    print()

    # Run EM algorithm
    final_abundance, n_iterations = run_em_algorithm(assignment_matrix)

    print(f"EM Algorithm converged after {n_iterations} iterations")
    print(f"Final abundance estimates: {final_abundance}")
    print()

    # Display results in a nice format
    isoform_names = ['Red', 'Blue', 'Green']  # Note: reordered to match notebook output
    results = dict(zip(isoform_names, [final_abundance[0], final_abundance[2], final_abundance[1]]))

    print("Final Isoform Abundance Estimates:")
    for name, abundance in results.items():
        print(f"  {name}: {abundance:.6f} ({abundance*100:.2f}%)")


if __name__ == "__main__":
    main()
