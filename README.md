# StapletonMathPackage

[![NuGet](https://img.shields.io/nuget/v/StapletonMathPackage.svg)](https://www.nuget.org/packages/StapletonMathPackage/)

**StapletonMathPackage** is a C# numerical methods library providing tools for linear algebra, optimization, root-finding, interpolation, and ODE solving.  It is pretty old, and I mostly started it as a wrapper around other math packages because I wanted to stick with double arrays in my program rather than specialized matrix classes. 
This version wraps around Math.Numerics for many of the linear algebra operations.   

Licensed under the [MIT License](LICENSE).

---

##  Features

### Linear Algebra
- Matrix and vector operations (`MatrixMath`, `VectorMath`)
- Determinant, trace, inverse, transpose
- Eigenvalues & eigenvectors (general & 3×3 symmetric closed form)
- Gaussian elimination, LU, banded solvers
- Utilities for stacking, extracting, reshaping, swapping rows/columns

### Calculus
- Numerical integration: trapezoid, Romberg
- Double integration
- Vector and matrix-valued integrands

### Root Finding
- **Newton-Raphson solvers**:
  - With Jacobian supplied
  - With numerical tangent
  - With local secant updates
  - With initial slope
- **Regula Falsi (False Position)** method

### Interpolation
- Cubic spline interpolation with derivatives
- Tridiagonal solver for spline coefficients

### Optimization
- Genetic Algorithm framework:
  - Chromosome representation
  - Mutation & perturbation
  - Tournament-based crossover & population evolution
  - Multicore evaluation via `ThreadLocal<Random>`

### ODE System Solvers
- Matrix exponential (Taylor, Padé, scaling & squaring)
- Constant step matrix exponential integration
- Variable coefficient systems
- Non-homogeneous systems with forcing functions

---

##  Installation

Install via NuGet:

```bash
dotnet add package StapletonMathPackage
```

Or in Visual Studio:

```
PM> Install-Package StapletonMathPackage
```

---

##  Usage Examples

### Matrix Operations
```csharp
using RandomMath;

double[,] A = { { 1, 2 }, { 3, 4 } };
double[,] B = { { 5, 6 }, { 7, 8 } };

double[,] C = MatrixMath.Multiply(A, B);
double detA = MatrixMath.Determinant(A);
```

### Cubic Spline Interpolation
```csharp
using RandomMath.CubicSplineFit;

var xs = new double[] { 0, 1, 2, 3 };
var ys = new double[] { 0, 1, 0, 1 };

var spline = new CubicSplineFit(xs, ys);

double dy = 0, d2y = 0;
double yInterp = spline.Interpolate(1.5, ref dy, ref d2y);
```

### Genetic Algorithm
```csharp
using RandomMath.GeneticAlgorithm;

public class SphereFunction : IOptimize {
    public int nX => 3;
    public double Eval(double[] X) => X.Sum(x => x*x);
    public IOptimize DeepCopy() => new SphereFunction();
}

var ga = new GeneticAlgorithm(new SphereFunction());
ga.Run();
Console.WriteLine($"Best fitness: {ga.Best.Fitness}");
```

### Newton-Raphson (Jacobian supplied)
```csharp
using RandomMath.NewtonRaphson;

public class SimpleFunc : IMatrixFunctionAndDerivative {
    public double[] Eval(double[] X) => new[] { X[0]*X[0] - 2 };
    public double[,] DEval(double[] X) => new double[,] { { 2*X[0] } };
}

var solver = new NewtonRaphsonJacobian(
    new double[] { 0 },    // target y
    new double[] { 1.0 },  // initial guess
    new SimpleFunc(),
    1e-10, 50
);
solver.Solve();
Console.WriteLine($"Root ≈ {solver.X[0]}");
```

### Matrix Exponential ODE Solve
```csharp
using RandomMath.ODESystemSolver;

double[,] A = { { 0, 1 }, { -1, 0 } };
var mexp = new MatrixExponential(A, 1.0);
var expAx = mexp.Solve(1.0); // exp(A*1.0)
```

---

##  Documentation

- All public APIs include XML comments visible in IntelliSense.
- Example usage included in this README.
- For advanced methods (ODE solvers, GA tuning), see source files in `RandomMath.*` namespaces.

---

##  Contributing

Contributions are welcome!  
Please open an issue or submit a PR with improvements, bug fixes, or new numerical routines.

Guidelines:
- Follow existing naming conventions.
- Provide XML documentation for new public APIs.
- Add unit tests where possible.

---

##  License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
