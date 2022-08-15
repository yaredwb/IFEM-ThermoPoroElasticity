# IFEM-ThermoPoroElasticity
Thermo-poro-elasticity solver based on the Isogeometric Analysis (IGA) toolbox [IFEM](https://github.com/OPM/IFEM)

The governing momentum, mass and energy balance equations that are solved by the IFEM-ThermoPoroElasticity solver are:

The governing equations may be summarized in terms of the field variables $ \bm u $, $ p^\rmw $ and $ T $ as

$$
\begin{aligned}
\bm L^\intercal \left[ \tilde{\bm D} \bm L \bm u - \alpha p^\rmw \tilde{\bm I} - \tilde{\bm D} \alpha_\rms T \tilde{\bm I} \right]  + \rho \bm b &= \bm 0 \\
\alpha \tilde{\bm I}^\intercal \bm L \frac{\partial {\bm u}}{\partial t} + c\frac{\partial p^\rmw}{\partial t} + \nabla \cdot \left[ -\frac{1}{\gamma_\rmw}\bm k(\nabla p^\rmw - \rho_\rmw \bm b) \right] &= 0 \\
(\rho c)_\rmef \frac{\partial T}{\partial t} + \rho_\rmw c_\rmw \left[ -\frac{1}{\gamma_\rmw}\bm k(\nabla p^\rmw - \rho_\rmw \bm b) \right]  \cdot \nabla T + \nabla \cdot (-\lambda \nabla T) &= Q.
\end{aligned}
$$

**Note:** This repository is not upto date with the latest version of the IFEM toolbox.
