# mathematics-for-computer-vision

A structured study repository working through the mathematical foundations required for computer vision and autonomous navigation engineering.

---

## Purpose

This repository covers the core mathematics that underpin computer vision, state estimation, and robotic systems: linear algebra, convex optimization, probability, and Lie theory. It implements exercises and derivations drawn from the primary source material listed below.

All meaningful changes are tracked through a structured issue and pull request workflow, producing an auditable record of progress over time.

---

## Source Material

### High Priority

- Strang, G., *Linear Algebra and Its Applications*, 5th ed. — with MIT OCW 18.06
- Boyd, S. and Vandenberghe, L., *Convex Optimization*
- Bertsekas, D. and Tsitsiklis, J., *Introduction to Probability*, 2nd ed.
- Solà, J., Deray, J., Atchuthan, D., *A Micro Lie Theory for State Estimation in Robotics*, arXiv:1812.01537

### Supplemental

- Trefethen, L. and Bau, D., *Numerical Linear Algebra* — SVD, QR, conditioning, stability
- Hubbard, J. and Hubbard, B., *Vector Calculus, Linear Algebra, and Differential Forms*, 5th ed.
- Nocedal, J. and Wright, S., *Numerical Optimization*, 2nd ed. — Gauss–Newton, Levenberg–Marquardt
- MIT OCW 18.065 — Matrix Methods in Data Analysis, Signal Processing, and ML

### Additional

- Axler, S., *Linear Algebra Done Right*, 4th ed.
- Wasserman, L., *All of Statistics*
- Hall, B., *Lie Groups, Lie Algebras, and Representations*, 2nd ed.
- MIT OCW 18.02 — Multivariable Calculus

---

## Repository Structure

```
mathematics-for-computer-vision/
├── .github/
│   ├── ISSUE_TEMPLATE/
│   └── PULL_REQUEST_TEMPLATE.md
├── CHANGELOG.md
├── CONTRIBUTING.md
└── README.md
```

Subdirectories will be added as work on each topic begins.

---

## Workflow

See [`CONTRIBUTING.md`](./CONTRIBUTING.md) for the full issue, branch, and PR protocol.

---

## Build Standard

Code in this repository targets **C++20** and is compiled with:

```bash
g++ -std=c++20 -Wall -Wextra -Werror -Wpedantic -o <output> <source>.cpp
```

Multi-file projects use CMake. No source file in this repository should produce compiler warnings.
