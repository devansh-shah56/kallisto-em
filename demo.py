#!/usr/bin/env python3
"""
Demo script showing how to use the Kallisto EM algorithm implementation.
"""

import numpy as np
from kallisto_em import run_em_algorithm, e_step, m_step

def demo_basic_usage():
    """Demonstrate basic usage of the EM algorithm."""
    print("🧬 Kallisto EM Algorithm Demo")
    print("=" * 40)

    # Create a simple example
    print("Creating a simple 2-isoform example...")
    assignment_matrix = np.array([
        [1, 1, 0, 1],  # Isoform A
        [0, 1, 1, 1]   # Isoform B  
    ])

    print(f"Assignment matrix shape: {assignment_matrix.shape}")
    print("Assignment matrix:")
    print(assignment_matrix)
    print()

    # Run EM algorithm
    print("Running EM algorithm...")
    abundance, iterations = run_em_algorithm(assignment_matrix, tol=1e-6)

    print(f"✅ Converged after {iterations} iterations")
    print(f"Final abundances: {abundance}")
    print(f"Isoform A: {abundance[0]:.4f} ({abundance[0]*100:.1f}%)")
    print(f"Isoform B: {abundance[1]:.4f} ({abundance[1]*100:.1f}%)")
    print()

def demo_step_by_step():
    """Demonstrate step-by-step execution."""
    print("🔍 Step-by-step EM demonstration")
    print("=" * 40)

    # Simple 2x2 example
    assignment_matrix = np.array([[1, 0], [1, 1]])
    alpha = np.array([0.6, 0.4])  # Initial guess

    print(f"Assignment matrix:\n{assignment_matrix}")
    print(f"Initial abundance: {alpha}")
    print()

    # One iteration
    print("E-step:")
    posterior = e_step(alpha, assignment_matrix)
    print(f"Posterior probabilities:\n{posterior}")
    print()

    print("M-step:")
    new_alpha = m_step(posterior, assignment_matrix)
    print(f"Updated abundances: {new_alpha}")
    print(f"Change: {np.abs(new_alpha - alpha)}")
    print()

if __name__ == "__main__":
    demo_basic_usage()
    demo_step_by_step()
