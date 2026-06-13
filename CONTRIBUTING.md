# Contributing

This document describes the self-review workflow used throughout this repository. As this repo is targetted towards self-improvement, external contributions are unlikely to be included. That said -- discussion, recommendations, questions, etcetera are all welcome so long as they flow through the proper channels.

---

## Guiding Philosophy

Each domain in the curriculum produces implementations, derivations, or algorithm prototypes. Most contributions in this repo are algorithms which are first created manually and then compared to common library implementations. This strategy aims to provide me, dcarey-cv, with a deeper understanding of the mathematical foundations of computer vision, as well as their associated computational overhead under the hood.

---

## Issue Types

| Template | When to Use |
|---|---|
| `Code Quality / Refactor` | An implementation that runs but is incorrect, inefficient, or non-idiomatic |
| `Concept Documentation` | A mathematical concept, language feature, or algorithm pattern worth recording directly |
| `Exercise Extension` | A textbook exercise or paper derivation worth expanding to deepen understanding |

Open issues from the **Issues** tab → **New Issue** → select the appropriate template.

---

## Branch Naming

```
refactor/<short-description>
docs/<short-description>
feat/<short-description>
```

Examples:
- `feat/svd-image-compression`
- `refactor/kalman-filter-eigen-types`
- `docs/lie-group-so3-notes`

---

## Pull Request Protocol

Every change goes through a PR, which are then reviewed and merged by the author, dcarey-cv, after confirming the checklist in the PR template is complete.

---

## Changelog

All merged PRs must produce an entry in `CHANGELOG.md` under the `[Unreleased]` section before merging. Follow the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

---

## Build Standard

All C++ code in this repository targets **C++20**. Single-file programs compile with:

```bash
g++ -std=c++20 -Wall -Wextra -Werror -Wpedantic -o <output> <source>.cpp
```

Multi-file projects use CMake. No PR should introduce new compiler warnings.
