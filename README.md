# Optimization

Gradient-based optimisation in MATLAB, with contour visualisation of the search.

| File | What it does |
| --- | --- |
| `Conjugate_Gradient_Branin_hoo_function.m` | Conjugate gradient method with a line search over the step length. Takes an optional starting point (default `[2 1]`), runs to a gradient-norm tolerance of 1e-6 or 10,000 iterations, and overlays the objective's contours for visualisation. The Branin-Hoo function is included; a quadratic test objective is active by default, so swap the commented `f` to switch between them. Returns `xopt`, `fopt`, the iteration count, the gradient norm, and the final step size. |
| `contour_test_shade.m` | Contour plot of a quadratic objective subject to three inequality constraints, shading the infeasible side of each so the feasible region and the constrained optimum are visible. Exports to PNG or PDF. |

## Requirements

MATLAB. Nothing beyond base plotting is needed.
