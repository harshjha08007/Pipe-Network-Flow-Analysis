# Pipe-Network-Flow-Analysis
Includes Python source code, EPANET input files, and result comparison plots for multi-loop pipe flow simulation and validation.

Implemented a Newton–Raphson based solver to compute steady-state pipe flows in a four-loop water distribution network, given inflow (5 units) and outflows (2 and 3 units). Pipe resistances were defined as 
R=k+p, with p=10x+y (based on roll number). Results were validated by simulating the same network in EPANET.
