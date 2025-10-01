# Kallisto EM Algorithm Implementation

[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python implementation of the Expectation-Maximization (EM) algorithm used by Kallisto for RNA isoform quantification, as described in Pachter (2011).

## 🧬 Overview

This repository contains a clean, well-documented implementation of the EM algorithm that Kallisto uses to assign RNA sequencing reads to different isoforms of a gene. The algorithm iteratively estimates transcript abundances by:

1. **E-step**: Computing posterior probabilities of read assignments
2. **M-step**: Updating transcript abundance estimates based on these probabilities

## 🚀 Features

- **Pure Python implementation** with NumPy for numerical computations
- **Educational focus** with clear documentation and examples  
- **Modular design** allowing easy integration into other projects
- **Convergence detection** with configurable tolerance
- **Example data** to demonstrate the algorithm

## 📁 Repository Structure

```
kallisto-em-algorithm/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── kallisto_em.py           # Main implementation
├── example_notebook.ipynb   # Original Jupyter notebook
├── test_kallisto_em.ppy     # Unit tests
└── algorithm_explanation.md  # Detailed algorithm explanation
```

## 🛠️ Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/devansh-shah56/kallisto-em
   cd kallisto-em
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

## 📊 Usage

### Basic Usage

```python
import numpy as np
from kallisto_em import run_em_algorithm, create_example_data

# Create example assignment matrix
assignment_matrix = create_example_data()

# Run EM algorithm
final_abundance, n_iterations = run_em_algorithm(assignment_matrix)

print(f"Final abundance estimates: {final_abundance}")
print(f"Converged in {n_iterations} iterations")
```

### Command Line

```bash
python kallisto_em.py
```

This will run the algorithm on example data and display the results.

### Custom Data

```python
import numpy as np
from kallisto_em import run_em_algorithm

# Your custom read-to-isoform assignment matrix
# Rows = isoforms, Columns = reads
assignment_matrix = np.array([
    [1, 0, 1, 1],  # Isoform 1 assignments
    [1, 1, 0, 1],  # Isoform 2 assignments  
    [0, 1, 1, 1]   # Isoform 3 assignments
])

# Run with custom parameters
abundance, iterations = run_em_algorithm(
    assignment_matrix, 
    tol=1e-5,        # Convergence tolerance
    max_iters=1000   # Maximum iterations
)
```

## 📈 Algorithm Details

The EM algorithm alternates between two steps:

### E-step (Expectation)
Calculate posterior probabilities that read *n* belongs to isoform *k*:

```
p(Z_n = k | Y, α^(t)) = (y_{k,n} * α_k^(t)) / (Σ_l y_{l,n} * α_l^(t))
```

### M-step (Maximization)  
Update abundance estimates:

```
α_k^(t+1) = (1/N) * Σ_n p(Z_n = k | Y, α^(t))
```

Where:
- `Y` is the assignment matrix
- `α` represents transcript abundances
- `Z_n` indicates which isoform read *n* comes from

## 🧪 Example Results

Running on the provided example data with 3 isoforms and 5 reads:

```
Assignment Matrix:
[[1 0 1 1 1]  # Red isoform
 [1 1 0 0 1]  # Green isoform
 [1 1 1 0 0]] # Blue isoform

Final Isoform Abundance Estimates:
  Red: 0.640296 (64.03%)
  Blue: 0.179852 (17.99%)
  Green: 0.179852 (17.99%)
```

## 🔬 Background

This implementation is based on the seminal paper:
> Pachter, L. (2011). Models for transcript quantification from RNA-Seq. *arXiv preprint arXiv:1104.3889*.

The EM algorithm is fundamental to many RNA-seq quantification tools including:
- Kallisto
- Sailfish  
- RSEM
- eXpress

## 🧑‍💻 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Original assignment implementation for DH607 course
- Pachter (2011) for the theoretical foundation
- The Kallisto development team for the practical implementation

## 📚 Further Reading

- [Kallisto paper](https://www.nature.com/articles/nbt.3519)
- [EM Algorithm Tutorial](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm)
- [RNA-seq quantification review](https://www.annualreviews.org/doi/abs/10.1146/annurev-biodata-073014-020548)

---

⭐ **If you find this implementation helpful, please consider giving it a star!**
