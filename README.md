# ODE/PDE Solver for MATLAB using GPUs

Suppose you really want to run a (2D) turbulence simulation but don't have access to C, Fortran, Julia, Python etc. It turns out that MATLAB is actually surprisingly fast when running on a GPU and generally outperforms it's multithreaded CPU performance by an order of magnitude.

This repository contains MATLAB functions for building ODEs and PDEs (on periodic grids) and solving them numerically using various timesteppers. Calculations can be performed on a GPU if available. Some examples are also included.