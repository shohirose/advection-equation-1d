# 1-D Advection Equation
C++ source code to solve one-dimensional scalar advection equation. Problems are taken from [1].

1-D scalar advection equation with constant velocity

$$
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0
$$

is solved by using the finie difference method. Currently, the following spacial reconstruction schemes are implemented:

- Lax-Wendroff scheme
- 1st-order upwind scheme
- Beam-Warming scheme
- Fromm scheme
- TVD scheme with minmod, Superbee, van Leer, and van Albada slope limiters.

In addition, periodic boundaries, the Roe-Riemann solver, and the explicit Euler scheme for time integration are used.

Please refer to [1] for the details of each scheme.

# How to compile

Run the following commands under the root directory of the project:

```
$ cmake -S . -B build
$ cmake --build build
```

Please note that the project depends on the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library, which is automatically downloaded and built by CMake using the `FetchContent` module.

Then, run all simulators and get results. For Linux,

```
$ ./run_all.sh
```

For Windows,

```
$ .\run_all.ps1
```

To visualize results, open [`plot.ipynb`](./plot.ipynb) with Jupyter Lab, and run all cells.

# References
1. 肖鋒・長﨑孝夫　2020　数値流体解析の基礎 －Visual C++とgnuplotによる圧縮性・非圧縮性流体解析－　コロナ社