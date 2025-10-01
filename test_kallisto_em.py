"""
Unit tests for the Kallisto EM algorithm implementation.
"""

import numpy as np
import pytest
from kallisto_em import e_step, m_step, run_em_algorithm, create_example_data


class TestKallistoEM:

    def test_e_step_basic(self):
        """Test E-step with simple data."""
        assignment_matrix = np.array([[1, 0], [1, 1]])
        abundance = np.array([0.5, 0.5])

        posterior = e_step(abundance, assignment_matrix)

        # Check shape
        assert posterior.shape == (2, 2)

        # Check that columns sum to 1 (probabilities)
        col_sums = posterior.sum(axis=0)
        np.testing.assert_array_almost_equal(col_sums, [1.0, 1.0])

    def test_m_step_basic(self):
        """Test M-step with simple data."""
        posterior = np.array([[0.5, 0.0], [0.5, 1.0]])
        assignment_matrix = np.array([[1, 0], [1, 1]])

        new_abundance = m_step(posterior, assignment_matrix)

        # Check that abundances sum to 1
        assert abs(new_abundance.sum() - 1.0) < 1e-10

        # Check shape
        assert len(new_abundance) == 2

    def test_run_em_algorithm_convergence(self):
        """Test that EM algorithm converges."""
        assignment_matrix = create_example_data()

        final_abundance, n_iterations = run_em_algorithm(
            assignment_matrix, tol=1e-6, max_iters=1000
        )

        # Check convergence
        assert n_iterations < 1000

        # Check that abundances sum to 1
        np.testing.assert_almost_equal(final_abundance.sum(), 1.0)

        # Check that all abundances are non-negative
        assert np.all(final_abundance >= 0)

    def test_create_example_data(self):
        """Test example data creation."""
        assignment_matrix = create_example_data()

        # Check shape
        assert assignment_matrix.shape == (3, 5)

        # Check that matrix contains only 0s and 1s
        assert np.all(np.isin(assignment_matrix, [0, 1]))

    def test_em_reproducibility(self):
        """Test that EM algorithm gives consistent results."""
        assignment_matrix = create_example_data()

        # Run multiple times
        results = []
        for _ in range(3):
            abundance, _ = run_em_algorithm(assignment_matrix, tol=1e-8)
            results.append(abundance)

        # Check that results are consistent
        for i in range(1, len(results)):
            np.testing.assert_array_almost_equal(results[0], results[i], decimal=6)

    def test_uniform_initialization(self):
        """Test uniform initialization."""
        assignment_matrix = create_example_data()

        abundance, _ = run_em_algorithm(assignment_matrix)

        # Should still converge to same result regardless of initialization
        assert abundance is not None
        np.testing.assert_almost_equal(abundance.sum(), 1.0)

    def test_edge_cases(self):
        """Test edge cases."""
        # Single isoform
        single_isoform = np.array([[1, 1, 1]])
        abundance, _ = run_em_algorithm(single_isoform)
        np.testing.assert_array_almost_equal(abundance, [1.0])

        # Equal assignment matrix  
        equal_matrix = np.array([[1, 1], [1, 1]])
        abundance, _ = run_em_algorithm(equal_matrix)
        np.testing.assert_array_almost_equal(abundance, [0.5, 0.5])


if __name__ == "__main__":
    pytest.main([__file__])
