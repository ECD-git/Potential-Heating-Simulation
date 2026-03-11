# Potential-Well-Heating-Simulation
A classical simulation to detail the heating of an anti-hydrogen atom in a potential well driven by a potential, in one dimension.

Utilises symplectic leapfrog algorithm for energy conserving numerical integration of equations of motion;

$$a_i(x_i) = A(x_i) = \frac{F(x_i)}{m} = \frac{1}{m}\bigg(-\frac{dU(x_i)}{dx} \bigg).$$
$$v_{i+1/2} = v_{i-1/2} + a_i\Delta t,$$
$$x_{i+1} = x_i + v_{i+1/2}\Delta t.$$

where $x,v,a$ are the position, velocity and acceleration of the particle in one dimension, $U$ is the potential the particle is under and $\Delta t$ is the time step of the simulation.

Simulation done in C++ with results plotted in Python (No technical reason just good for practicing both)

---

### Spring Potential (Static) Test:

$$F = -kx, \hspace{1cm}\phi = \frac{1}{2}kx^2.$$

Resulting kinetic energy graph and $\Delta t = 0.0001\text{s}$ and $T = 1\text{s}$, for a particle with no initial velocity at $x_i = -0.25\text{m}$.
![[SpringPotentialStatic-1second-0.0001.png]]
![image](https://github.com/ECD-git/Potential-Heating-Simulation/blob/main/Figures%2011%3A3%3A26/F%3D-kx.png)
Shows energy conservation, oscillating back and forth over a potential proportional to $x^2$.

---