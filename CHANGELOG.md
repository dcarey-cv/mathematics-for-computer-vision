# Changelog

All notable changes to this repository are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased]

---

## [2026-06-14]

### strang-linear-algebra/chapter_01.cpp

#### Changed
- 'gaussianElimination()' function
  - added abilty handle non-pivot-ordered matrices
  - added singular matrix error handling

#### Added
- bool partialPivot()
- Eigen::MatrixXd reconstructL()
- Eigen::MatrixXd reconstructD()
- Eigen::MatrixXd reconstructU()
- bool gaussElimComparison()


---

## [2026-06-13]

### strang-linear-algebra/chapter_01.cpp

#### Added
- 'gaussianElimination()' function which automates the Gaussian Elimination algorithm over dynamically-sized, pivot-ordered matrices
  - this algorithm is compared to the `Eigen` inbuilt solve function