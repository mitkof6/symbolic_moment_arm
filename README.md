Symbolic Derivation of the Muscle Moment Arm Matrix 
---

[OpenSim](https://github.com/opensim-org/opensim-core) is a framework for
modeling and simulation of musculoskeletal systems. The muscle moment arm is an
important variable for evaluating the effectiveness of a muscle to actuate a
particular joint. Calculating the muscle moment arm requires knowledge of the
muscle path and wrapping surfaces. OpenSim is restricted to evaluate the muscle
moment arm at an arbitrary configuration (e.g., <a
href="https://www.codecogs.com/eqnedit.php?latex=R(q)&space;\in&space;\mathcal{R}^{n&space;x&space;m},&space;q&space;\in&space;\mathcal{R}^{n}"
target="_blank"><img
src="https://latex.codecogs.com/gif.latex?R(q)&space;\in&space;\mathcal{R}^{n&space;x&space;m},&space;q&space;\in&space;\mathcal{R}^{n}"
title="R(q) \in \mathcal{R}^{n x m}, q \in \mathcal{R}^{n}" /></a> with *n*
degrees of freedom and *m* muscles), lacking the information for calculating
higher order derivatives (e.g., <a
href="https://www.codecogs.com/eqnedit.php?latex=\partial&space;R(q)&space;/&space;\partial&space;q"
target="_blank"><img
src="https://latex.codecogs.com/gif.latex?\partial&space;R(q)&space;/&space;\partial&space;q"
title="\partial R(q) / \partial q" /></a>). This project evaluates the moment
arm at different configurations and approximates its terms using multivariate
polynomial fitting, thus a symbolic expression is derived. Examples are provided
for the gait2392 model (23 DoFs and 92 muscles). Results are calculated and
stored in the *.dat* files which can be loaded using python's pickle
utility. The *R.dat* file is the symbolic expression (sympy Matrix) of the
muscle moment arm matrix. Visual inspection of the polynomial fitting is
provided below.


![Moment arm of vas_int_r at knee joint](data/gait2392/fig/vas_int_r_knee_angle_r.png)

![Moment arm of tib_ant_r at knee joint](data/gait2392/fig/tib_ant_r_ankle_angle_r.png)

Dependencies
---

There are two versions that work either with OpenSim v3.3 or OpenSim v4.0.

- OpenSim v3.3: Python 2.7 bindings 
- OpenSim v4.0: Python 3.7 bindings
- sympy: `pip install sympy`
- numpy: `pip install numpy`
- matplotlib: `pip install matplotlib` (use custom implementation that
  fixes some bugs with Python 3.7)
- multipolyfit: `pip install multipolyfit` for multivariate polynomial fitting


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img
alt="Creative Commons License" style="border-width:0"
src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is
licensed under a <a rel="license"
href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution
4.0 International License</a>.
