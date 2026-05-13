# 💻 Numerical Analysis Course Projects

**Course:** Fundamentals of Numerical Analysis

**Status:** Complete

This repository contains algorithmic implementations focused on Polynomial and Spline Interpolation. The project demonstrates the ability to approximate complex functions using various mathematical methods and analyze their accuracy through error visualization.

---

## 📂 Project Overview

| Project | Concept | Description |
| --- | --- | --- |
| **Interpolation Suite** | Numerical Approximation | A MATLAB tool to calculate Newton Polynomials and Cubic Splines with error analysis. |

---

## 1️⃣ Function Interpolator & Error Analyzer

**Location:** `/Numerical-Interpolation`

This project is a comprehensive MATLAB application that takes a symbolic function and a set of distinct knot points to generate several types of interpolating functions. It provides a visual comparison between the original function and its approximations.

### Key Features

* **Newton Divided Differences**: Calculates the unique polynomial of degree $n$ that passes through $n+1$ points using the divided difference table method.
* **Natural Cubic Spline**: Constructs a piecewise-polynomial approximation where the second derivatives at the endpoints are set to zero ($S''(x_0) = S''(x_n) = 0$).
* **Clamped (Restricted) Cubic Spline**: Generates a spline where the first derivatives at the boundaries match the original function's derivatives ($S'(x_0) = f'(x_0)$ and $S'(x_n) = f'(x_n)$).
* **Comparative Visualization**: Generates high-resolution plots comparing:
* The original function $f(x)$.
* The Newton interpolating polynomial.
* Both Natural and Clamped Cubic Splines.


* **Error Analysis**: Specifically plots $|f(x) - P(x)|$ for each method, allowing for a direct visual assessment of which approximation is most stable and accurate for the given knot points.

### Technical Implementation

* **`AmirNoshadiNew.m`**: The primary script handling symbolic math, matrix construction for spline coefficients, and the plotting engine.
* **Symbolic Processing**: Uses MATLAB's `syms` and `str2func` to allow users to input any mathematical rule (e.g., `sin(x) + exp(x)`) directly from the console.
* **Divided Difference Table**: Efficiently computes coefficients for the Newton form of the interpolating polynomial.
* **Tri-diagonal Matrix Solver**: Solves the system of equations required to find the $c_j$ coefficients for cubic splines.

### 💻 How to Run

1. **Open MATLAB**.
2. **Run the script**:
```matlab
AmirNoshadiNew

```


3. **Provide Inputs**:
* **Function**: e.g., `exp(x)`
* **Knot Points**: e.g., `[0, 0.5, 1, 1.5, 2]`



---

## 🛠 Tech Stack

* **Language:** MATLAB
* **Toolboxes:** Symbolic Math Toolbox
* **Approach:** The project emphasizes the **Runge's Phenomenon** observation—demonstrating how higher-order polynomials might oscillate and why piecewise cubic splines often provide a superior fit for numerical data.


**Date:** January 2024
