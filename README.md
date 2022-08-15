# IFEM-ThermoPoroElasticity
Thermo-poro-elasticity solver based on the Isogeometric Analysis (IGA) toolbox [IFEM](https://github.com/OPM/IFEM)

The governing momentum, mass and energy balance equations that are solved by the IFEM-ThermoPoroElasticity solver are:

$$
\begin{aligned}
\pmb L^\intercal \left[ \tilde{\pmb D} \pmb L \pmb u - \alpha p^w \tilde{\pmb I} - \tilde{\pmb D} \alpha_s T \tilde{\pmb I} \right]  + \rho \pmb b &= \pmb 0 \\
\alpha \tilde{\pmb I}^\intercal \pmb L \frac{\partial {\pmb u}}{\partial t} + c\frac{\partial p^w}{\partial t} + \nabla \cdot \left[ -\frac{1}{\gamma_w}\pmb k(\nabla p^w - \rho_w \pmb b) \right] &= 0 \\
(\rho c)_{eff} \frac{\partial T}{\partial t} + \rho_w c_w \left[ -\frac{1}{\gamma_w}\pmb k(\nabla p^w - \rho_w \pmb b) \right]  \cdot \nabla T + \nabla \cdot (-\lambda \nabla T) &= Q.
\end{aligned}
$$

**Note:** This repository is not upto date with the latest version of the IFEM toolbox.
