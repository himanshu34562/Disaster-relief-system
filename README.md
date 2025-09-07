# Disaster-relief-system
Implementing a disaster relief system using local search with random restarts.
# 🚁 Helicopter Relief Routing Solver

This repository contains my solution for the **Helicopter Relief Problem**.  
The problem models a disaster relief scenario where helicopters deliver supplies to stranded villages while respecting helicopter and logistical constraints.  

The solver is written in **C++17** and uses **local search with restarts** to maximize the total value of delivered aid.

---

## ✨ Features

- Constraint-aware search:
  - Per-village supply caps (9 meals/person, 1 other/person).
  - Per-trip helicopter limits (weight and distance).
  - Per-helicopter maximum total distance (`DMax`).
- Multi-restart hill climbing with sideways moves (up to 8).
- Trimmed neighborhood operators (add, remove, tweak, shuffle).
- Time-bounded search:
  - Algorithm stops at **90% of the allowed time limit**.
- Post-processing:
  - Trip merging and oversupply trimming.
- Tested on Ubuntu 20.04 with g++.

---

## 📂 Project Structure

.
├── solver.cpp # Solver implementation
├── solver.h # Data structures & declarations
├── main.cpp # Entry point, handles input/output
├── Makefile # Build instructions (creates executable main)
├── writeup.txt # Course writeup (for submission)
└── README.md # This file


---

## 🚀 Build & Run

### Build
```bash
make

./main input.txt output.txt

make clean

📄 Input Format (Sample)
1
100
0.01 1 0.1 2 0.005 0.1
2 0 0 10 10
2 0 5 1000 0 10 1000
2 1 100 25 10 1 2 100 50 10 1

📤 Output Format (Sample)
1 2
9000 0 2000 2 1 9000 0 1000 2 0 0 1000
8889 111 0 1 2 8889 111 0
-1
2 0
-1

⚙️ Algorithm Overview

Initialization

Each helicopter starts with one feasible trip to a single village.

Supply allocation favors perishable food when possible.

Local Search

Neighborhood moves adjust trips (add/remove/tweak/shuffle).

Hill climbing with sideways moves up to 8.

Multiple restarts until 90% of the time limit is reached.

Post-Processing

Oversupply is trimmed to respect village caps.

Trips are merged when feasible.

Final evaluation ensures validity.

📜 Full Problem Statement

Goal: The goal of this assignment is to take a complex new problem and formulate and solve it as search. Formulation as search is an integral skill of AI that will come in handy whenever you are faced with a new problem. Heuristic search will allow you to find optimal solutions. Local search may not find the optimal solution, but is usually able to find good solutions for really large problems.

Scenario
There are floods in a region. There is an urgent need to carry out relief operations from neighboring cities. There are a fixed set of relief packages to deliver to affected villages using a fleet of helicopters. Each village has a number of people stranded. The goal is to allocate relief deliveries to helicopters such that maximum aid is delivered, while respecting capacity and range limits of each helicopter, and while simultaneously minimizing a “logistical strain” cost caused by routing aid inefficiently.

For this we assume, we are provided a list of villages V that need relief. Each village v is associated with its 2D coordinates (xv, yv). Similarly, we have a few cities C, along with 2D coordinates for each city. The distance between any two locations is given by the straight line distance between their coordinates (in kms). For each village, we are also provided nv: the number of people expected to be stranded in that village.

There are three types of packets T = {d, p, o}, which stand for dry food, wet (perishable) food, and other supplies. The goal is to send about 9 meals per stranded person, and about 1 unit of other supplies per stranded person. Each packet of type t weighs a fixed amount w(t). Each packet of type t has a fixed value v(t). The meals may be dry or perishable food, but wet food is preferable, i.e., v(p) > v(d).

There are H helicopters. Each helicopter h has a home city home(h). It also has a weight capacity wcap(h), and a distance capacity dcap(h) per trip. Each trip starts at the home city and ends at the home city, and package weight on the trip cannot exceed wcap(h) and total distance traveled cannot exceed dcap(h). Moreover, total distance traveled by any helicopter across trips cannot exceed a given DMax.

Each trip costs F + alpha*distance. Here F is the fixed cost per trip (for takeoff and landing). And alpha represents some notion of fuel efficiency.
Total value of solution = Total value achieved – total trip cost.

The goal of the assignment is to produce a plan for each helicopter, such that the total value is maximized. How many trips they do. How many total packages of each type they start with per trip. Which villages do they visit and in what order. How many packages of each type do they drop per village.

Input Format

Total processing time available in minutes.

DMax: max distance in kilometres.

Six numbers: w(d) v(d) w(p) v(p) w(o) v(o).

Number of cities C, followed by 2C coordinates.

Number of villages V, followed by 3V numbers: x, y, stranded population.

Number of helicopters H, followed by 5H numbers: home city id, weight cap, distance cap, F, alpha.

Output Format

For each helicopter (in order 1..H):

Print helicopter number and number of trips.

For each trip:

Packages picked up (d, p, o).

Number of villages visited.

For each village: village id and dropped packages.

End with -1.

Sample Input

1
100
0.01 1 0.1 2 0.005 0.1
2 0 0 10 10
2 0 5 1000 0 10 1000
2 1 100 25 10 1 2 100 50 10 1


Sample Output

1 2
9000 0 2000 2 1 9000 0 1000 2 0 0 1000
8889 111 0 1 2 8889 111 0
-1
2 0
-1

🙋 Author

Developed by [Your Name]

💻 Focus: Algorithms, Optimization, AI Search

🌐 GitHub: yourusername

