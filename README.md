# 🚁 Helicopter Relief Routing Solver

A solver for the **Helicopter Relief Problem**, modeling disaster relief logistics where helicopters deliver aid to stranded villages while obeying all logistical and helicopter-specific constraints.  
Written in **C++17** and utilizing advanced local search with restarts, it aims to maximize the value of delivered aid within feasible, efficient routes.

---

## ✨ Features

- **Constraint-aware Search**
  - Enforces village supply caps: max 9 meals/person and 1 other supply/person
  - Respects helicopter trip weight and distance limits
  - Per-helicopter maximum distance (`DMax`)
- **Robust Search Strategy**
  - Multi-restart hill climbing, including up to 8 sideways moves
  - Diverse, trimmed neighborhood operators (add, remove, tweak, shuffle)
  - Time-bounded search, halts at 90% of allowed time
- **Smart Post-processing**
  - Oversupply trimming to avoid excess per village
  - Trip merging for route simplification
  - Strict, final feasibility verification
- **Easy Build & Portability**
  - Provided `Makefile` for out-of-the-box compilation
  - Tested on Ubuntu 20.04 with g++

---

## 📂 Project Structure

.
├── solver.cpp      # Solver implementation
├── solver.h        # Data structures & declarations
├── main.cpp        # Entry point, input/output
├── Makefile        # Build instructions
├── writeup.txt     # Course writeup
└── README.md       # This file


---

## ⚙️ Algorithm Overview

### Initialization
- Each helicopter starts with one trip to a village.
- Package allocation is demand-aware, prioritizing perishable food.

### Local Search
- Neighborhood moves: add, remove, tweak, shuffle
- Hill climbing, allowing up to 8 sideways steps per restart
- Multiple restarts until 90% of time expires

### Post-processing
- Trim oversupply beyond village demand
- Merge trips when feasible
- Final scoring and strict feasibility checks

---

## 📝 Full Problem Statement

The goal is to model and solve a new AI search problem: disaster relief logistics with helicopters. Given villages, cities, helicopters, and package types, allocate relief to maximize aid delivered—subject to all practical constraints and minimizing logistical strain cost.

### Constraints
- **Villages:** max 9 meals/person & 1 other supply/person
- **Per helicopter:** trip weight/distance limits, total distance ≤ `DMax`

### Objective

\[
\text{Objective} = (\text{Total delivered value}) - (\text{Total cost})
\]

Trip cost: \( F + \alpha \times \text{distance} \)

---

## 🔧 Future Improvements

- Smarter initialization heuristics
- Richer, more diverse local moves (e.g., cross-helicopter package swaps)
- Parallel restarts for faster, more robust search

---

## 🤝 Author

Developed by **Himanshu Singh**  
💻 Focus: Algorithms, Optimization, AI Search

---
