# Aerodynamics of Biplanes and Tandem Wings at Low Reynolds Number

This project reproduces selected results from  
Jones et al. (2015) – *Aerodynamics of Biplane and Tandem Wings at Low Reynolds Numbers*  
using analytical aerodynamic theories.

The goal was to evaluate how far inviscid models can replicate lift trends at Re = 100,000.

## Objectives

- Replicate lift curve of a single flat plate (Figure 3)
- Replicate lift ratio of tandem/biplane configuration (Figure 5a)
- Compare theoretical predictions with experimental trends
- Assess limitations of classical aerodynamic models at low Reynolds number

## Models Used

### 1) Thin Airfoil Theory (2D)

- Computed Fourier coefficients of extracted camber line
- Lift slope ≈ 2π (ideal inviscid result)
- Correctly predicts linear trend
- Overestimates lift slope at low Re
- Cannot predict stall

### 2) Prandtl Lifting Line Theory (3D)

- Included induced effects and finite aspect ratio
- Improved lift slope prediction
- Extended to tandem configuration by computing mutual induced velocities

## Main Results

Single wing:
- Linear lift region well captured
- Real lift slope ≈ 0.5–0.6 of ideal 2π
- Stall cannot be predicted (no viscous modeling)

Tandem configuration:
- Qualitative lift-ratio trend reproduced
- ~0.55 lift ratio near zero x-gap
- Divergence at larger gaps
- Symmetry discrepancy (~0.25 gap instead of 0)

## Key Takeaways

- Thin Airfoil Theory captures 2D trends but is too optimistic
- Lifting Line Theory improves realism via 3D induced effects
- Neither model predicts stall or post-stall behavior
- Strong insight into limits of inviscid aerodynamics at low Reynolds numbers

## Tools

- MATLAB
- Fourier series implementation
- Numerical integration
- Iterative solvers with relaxation


