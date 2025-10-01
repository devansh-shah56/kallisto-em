# Kallisto EM Algorithm Explanation

## Overview

The Expectation-Maximization (EM) algorithm is at the heart of modern RNA-seq quantification tools like Kallisto. This document provides a detailed explanation of how the algorithm works to estimate transcript abundances from ambiguous read mappings.

## Problem Setup

### The Challenge
When RNA sequencing reads map to multiple isoforms of a gene, we need to probabilistically assign each read to its most likely source isoform. This is complicated because:

1. **Reads are ambiguous** - A single read may map equally well to multiple isoforms
2. **Isoform abundances are unknown** - We want to estimate these
3. **Read assignments depend on abundances** - More abundant isoforms should receive more reads

### Mathematical Formulation

Let:
- `Y` = assignment matrix (which isoforms each read can map to)
- `α` = transcript abundance vector (what we want to estimate)  
- `Z` = latent variable indicating true isoform assignment for each read

## The EM Algorithm

### Initialization
Start with uniform abundance estimates:
```
α₀ = [1/K, 1/K, ..., 1/K]
```
where K is the number of isoforms.

### Iteration t

#### E-step: Estimate Read Assignments
For each read n and isoform k, calculate the posterior probability:

```
p(Zₙ = k | Y, α⁽ᵗ⁾) = (y_{k,n} × α_k⁽ᵗ⁾) / (∑ᵢ y_{l,n} × α_l⁽ᵗ⁾)
```

**Intuition**: Reads are assigned proportionally to:
1. Whether they CAN map to the isoform (y_{k,n})
2. How abundant that isoform currently is (α_k⁽ᵗ⁾)

#### M-step: Update Abundance Estimates
```
α_k⁽ᵗ⁺¹⁾ = (1/N) × ∑ₙ p(Zₙ = k | Y, α⁽ᵗ⁾)
```

**Intuition**: New abundance = average probability that reads come from this isoform

### Convergence
Continue until abundance estimates stabilize:
```
max |α_k⁽ᵗ⁺¹⁾ - α_k⁽ᵗ⁾| < tolerance
```

## Worked Example

Consider 3 isoforms (Red, Green, Blue) and 5 reads with this assignment matrix:

```
        Read1 Read2 Read3 Read4 Read5
Red   [  1,    0,    1,    1,    1 ]
Green [  1,    1,    0,    0,    1 ]  
Blue  [  1,    1,    1,    0,    0 ]
```

### Iteration 0
Start with uniform abundances: α = [1/3, 1/3, 1/3]

### Iteration 1

**E-step**: Calculate posterior probabilities
For Read1 (can map to all isoforms):
- P(Z₁ = Red) = (1 × 1/3) / (1×1/3 + 1×1/3 + 1×1/3) = 1/3
- P(Z₁ = Green) = 1/3  
- P(Z₁ = Blue) = 1/3

For Read2 (maps to Green and Blue only):
- P(Z₂ = Red) = 0
- P(Z₂ = Green) = (1 × 1/3) / (0×1/3 + 1×1/3 + 1×1/3) = 1/2
- P(Z₂ = Blue) = 1/2

... and so on for all reads.

**M-step**: Update abundances
```
α_Red = (1/5) × (1/3 + 0 + 1 + 1 + 1/2) = 14/30 ≈ 0.467
α_Green = (1/5) × (1/3 + 1/2 + 0 + 0 + 1/2) = 8/30 ≈ 0.267  
α_Blue = (1/5) × (1/3 + 1/2 + 1 + 0 + 0) = 8/30 ≈ 0.267
```

The algorithm continues until convergence...

### Final Result
After convergence:
- Red: 64.03%
- Green: 17.99%  
- Blue: 17.99%

## Key Insights

1. **Red isoform dominates** because it can uniquely claim Read3 and Read4
2. **Green and Blue tie** because they have symmetric assignment patterns for the ambiguous reads
3. **The algorithm balances** between explaining the data and maintaining parsimony

## Implementation Details

### Numerical Stability
- Use log-space computations for very small probabilities
- Handle division by zero when no reads map to an isoform
- Ensure probabilities sum to 1 within numerical precision

### Convergence Criteria
- Typical tolerance: 1e-4 to 1e-6
- Maximum iterations: 1000-10000
- Can also check likelihood convergence

### Computational Complexity
- Time: O(iterations × reads × isoforms)  
- Space: O(reads × isoforms)

## Extensions

### Effective Length Correction
Real implementations account for isoform lengths:
```
p(Zₙ = k) ∝ α_k × effective_length_k
```

### Bootstrapping
Multiple runs with resampled reads provide uncertainty estimates.

### Fragment-level Models
Paired-end reads require more sophisticated assignment models.

## References

1. Pachter, L. (2011). Models for transcript quantification from RNA-Seq
2. Bray et al. (2016). Near-optimal probabilistic RNA-seq quantification (Kallisto paper)
3. Dempster et al. (1977). Maximum likelihood from incomplete data via the EM algorithm
