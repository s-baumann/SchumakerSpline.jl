```@meta
CurrentModule = SchumakerSpline
```

# Internal Functions

```@index
Pages = ["api.md"]
```

### Main Struct

```@docs
    Schumaker
    evaluate
    evaluate_integral
```

### Extrapolation Schemes

```@docs
    Schumaker_ExtrapolationSchemes
```

### Working with Splines

```@docs
find_derivative_spline
find_roots
find_optima
get_crossover_in_interval
get_intersection_points
reshape_values
splice_splines
```



### Two dimensional splines

```@docs
Schumaker2d
```

### Plotting of splines

```@docs
plot
```

### Internal Functions

```@docs
schumakerIndInterval
imputeGradients
quadratic_formula_roots
test_if_intercept_in_interval
getCoefficientMatrix
extrapolate
```

