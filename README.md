# 🚁 Helicopter Relief Routing Solver

This project implements a solver for the **Helicopter Relief Problem**, a disaster relief logistics scenario where helicopters deliver supplies to stranded villages while respecting helicopter and logistical constraints.  

The solver is written in **C++17**, uses **local search with restarts**, and maximizes the total value of delivered aid while keeping all constraints valid.

---

## ✨ Features

- **Constraint-aware search**  
  - Per-village supply caps (9 meals/person, 1 other/person).  
  - Per-trip helicopter limits (weight and distance).  
  - Per-helicopter maximum distance (`DMax`).  

- **Search strategy**  
  - Multi-restart hill climbing with sideways moves (up to 8).  
  - Trimmed neighborhood operators (add/remove/tweak/shuffle).  
  - Time-bounded search that stops at **90% of the allowed time**.  

- **Post-processing**  
  - Oversupply trimming.  
  - Trip merging to simplify routes.  
  - Strict feasibility checks.  

- **Portable build system**  
  - Provided `Makefile` compiles into an executable `main`.  
  - Tested on Ubuntu 20.04 with g++.  

---

## 📂 Project Structure

.
├── solver.cpp # Solver implementation
├── solver.h # Data structures & declarations
├── main.cpp # Entry point, input/output
├── Makefile # Build instructions
├── writeup.txt # Course writeup
└── README.md # This file

yaml
Copy code

---

## 🚀 Build & Run

### Build
```bash
make
Run
bash
Copy code
./main input.txt output.txt
Clean
bash
Copy code
make clean
📄 Input Format (Sample)
yaml
Copy code
1
100
0.01 1 0.1 2 0.005 0.1
2 0 0 10 10
2 0 5 1000 0 10 1000
2 1 100 25 10 1 2 100 50 10 1
📤 Output Format (Sample)
yaml
Copy code
1 2
9000 0 2000 2 1 9000 0 1000 2 0 0 1000
8889 111 0 1 2 8889 111 0
-1
2 0
-1
⚙️ Algorithm Overview
Initialization

Each helicopter starts with one trip to a village.

Package allocation is demand-aware, with preference for perishable food.

Local Search

Neighborhood moves (add/remove/tweak/shuffle).

Hill climbing with sideways moves (max 8).

Multiple restarts until 90% of the time limit.

Post-Processing

Trim oversupply beyond village demand.

Merge trips when possible.

Final scoring and feasibility checks.

📜 Full Problem Statement
Goal: The goal of this assignment is to take a complex new problem and formulate and solve it as search. Formulation as search is an integral skill of AI that will come in handy whenever you are faced with a new problem. Heuristic search will allow you to find optimal solutions. Local search may not find the optimal solution, but is usually able to find good solutions for really large problems.

Scenario
There are floods in a region. There is an urgent need to carry out relief operations from neighboring cities. There are a fixed set of relief packages to deliver to affected villages using a fleet of helicopters. Each village has a number of people stranded. The goal is to allocate relief deliveries to helicopters such that maximum aid is delivered, while respecting capacity and range limits of each helicopter, and while simultaneously minimizing a “logistical strain” cost caused by routing aid inefficiently.

We are provided with:

Villages V (coordinates and stranded population).

Cities C (coordinates).

Helicopters H (home city, weight and distance capacities, fixed and per-distance costs).

Packages T = {d, p, o} with weights and values.

Global maximum distance DMax.

Constraints:

Villages: max 9 meals/person and 1 other supply/person.

Helicopter trips: respect weight capacity and distance capacity.

Each helicopter’s total travel ≤ DMax.

Cost model:
Trip cost = F + α * distance.
Objective = (Total value delivered) – (Total cost).

Input Specification
Processing time in minutes.

DMax (max distance).

Six numbers: w(d) v(d) w(p) v(p) w(o) v(o).

Cities: C followed by 2C coordinates.

Villages: V followed by 3V numbers (x, y, stranded population).

Helicopters: H followed by 5H numbers (home id, wcap, dcap, F, α).

Output Specification
For each helicopter (1..H):

Helicopter number and number of trips.

For each trip:

Packages picked up (d, p, o).

Number of villages visited.

For each village: village id + dropped packages.

End with -1.

Example
Sample Input

yaml
Copy code
1
100
0.01 1 0.1 2 0.005 0.1
2 0 0 10 10
2 0 5 1000 0 10 1000
2 1 100 25 10 1 2 100 50 10 1
Sample Output

yaml
Copy code
1 2
9000 0 2000 2 1 9000 0 1000 2 0 0 1000
8889 111 0 1 2 8889 111 0
-1
2 0
-1
🔧 Future Improvements
Smarter initialization heuristics.

More diverse neighborhood moves (cross-helicopter package swaps).

Parallel restarts for faster performance.

🙋 Author
Developed by [Your Name]

💻 Focus: Algorithms, Optimization, AI Search

🌐 GitHub: yourusername
