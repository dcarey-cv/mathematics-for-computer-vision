# Contributing

This document describes the self-review workflow used throughout this repository.

---

## Guiding Philosophy

Each domain in the curriculum produces implementations, derivations, or algorithm prototypes. This workflow exists to ensure every meaningful piece of work is preceded by a written intent, reviewed against a checklist, and logged — not to simulate a team environment, but to produce a legible, honest record of growth.

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

Every change goes through a PR. PRs are self-merged after confirming the checklist in the PR template is complete.

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
