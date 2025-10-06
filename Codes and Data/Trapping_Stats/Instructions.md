# 📊 Trapping Time Statistics for Glassy and Conventional Molecular Clutch Models

This folder contains all simulation codes, analysis scripts, and raw figures used to compute **trapping time statistics** for both the **Glassy** and **Conventional** molecular clutch models.

---

## 🧠 Overview

The *molecular clutch model* describes how stochastic binding and unbinding of clutches mediate the transmission of cytoskeletal forces to the substrate.  
A **trapping event** corresponds to a stable period where clutches remain bound before a complete collapse (detachment) occurs.

---

## ⚙️ Trapping Time Definition

For both the **Glassy** and **Conventional** models:

1. The number of bound clutches (`Ncbound`) is monitored independently for the **left** and **right** sides of the system.  
2. A **collapse event** (or **failure**) occurs whenever the `Ncbound` value of **either side** becomes zero.  
3. The **time interval between two consecutive collapse events** is defined as a **trapping time**.  
4. Collecting and analyzing these time intervals provides a statistical description of the trapping dynamics.

---

## 🧩 Folder Contents

This directory includes the following key files:
  
- `Glassy_Trapping_Statistics.m` — Runs 50 independent glassy model simulations and aggregates trapping times into one histogram.   
- `Conventional_Trapping_Statistics.m` — Runs 50 independent conventional model simulations and aggregates trapping times into one histogram.  
- `*_One/50_Trajectory.fig` — Figure files containing raw output figures from both single and multi-trajectory runs.

---

## 📈 Output Description

Each multi-trajectory script produces:

- A **histogram of trapping times**, representing the probability distribution of durations between consecutive clutch collapses.  


---

## 🧪 Reproducibility Notes

- All simulations include stochastic clutch dynamics; therefore, individual trajectories vary.  
- Running multiple (e.g., 50) trajectories ensures reliable ensemble statistics and smoother distributions.  
- Random seeds can be fixed for exact reproducibility across runs.

---

## 🗂️ Suggested Workflow

1. Run the single-trajectory script in the main folder to visualize clutch dynamics in one run.  
2. Run the 50-trajectory script (`*_Trapping_Statistics.m`) to gather statistical distributions.  
3. Analyze and compare the trapping time histograms of the **Glassy** and **Conventional** models.

---

