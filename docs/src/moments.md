# Moments

Or, rotational motion.

::: info

###  An aside on rotations

[There are many ways to describe rotation in three dimensions](https://en.wikipedia.org/wiki/Rotation_formulations_in_three_dimensions). I prefer to use the direction cosine matrix (DCM for short). It is a linear operator on vectors

A DCM ``\boldsymbol{C}_{BI}`` is a ``3\times 3`` matrix of real numbers with the following properties:

```math
\begin{align}
  \boldsymbol{C}_{BI}^\mathsf{T}\boldsymbol{C}_{BI} &= \boldsymbol{1}  \label{eqn:dcm-ident} \\
  \det \boldsymbol{C}_{BI} &= 1
\end{align}
```

Note that equation $\eqref{eqn:dcm-ident}$ implies that the transpose of a DCM is equal to its inverse.

We often know a vector in some reference frame and want to express that vector in a different reference frame. For example denote the spacecraft body frame ``B`` and the Earth-centered, inertial frame ``I``. If some vector ``\boldsymbol{u}_I`` is expressed in ``I``, it can be transformed to ``B`` by left-multiplication by the DCM:

The set of all valid DCMs is known as the special orthogonal group of order 3, or ``\mathrm{SO}(3)``. This is the same group that rotation matrices belong to. A rotation matrix is the transpose of a DCM. This is relevant to the order in which a DCM vs. rotation matrix operates on a vector to transform it.


:::

## Dynamics

Rotational motion is simulated by propagating the state

```math
(\boldsymbol{C}_{BI}, \boldsymbol{\omega}_B^{BI})
```

where ``\boldsymbol{C}_{BI}`` is the direction cosine matrix, and ``\boldsymbol{\omega}_B^{BI}`` is the angular velocity. The angular velocity is propagated under Euler's equation for rigid bodies:

```math
\dot{\boldsymbol{\omega}}_B^{BI} = \boldsymbol{I}_B^{-1}\big(\sum_k \boldsymbol{m}^k_B - \boldsymbol{\omega}_B^{BI}\times\boldsymbol{I}_B\boldsymbol{\omega}_B^{BI} \big)
```

and the DCM evolves kinematically following the angular velocity:

```math
\dot{\boldsymbol{C}}_{BI} = -\boldsymbol{\omega}_B^{BI^\times}\boldsymbol{C}_{BI}
```

## Control

Attitude control is implemented using a geometric PD law. Let the reference attitude state be 

```math
(\boldsymbol{C}_{BI}^r, \boldsymbol{\omega}_B^{BI^r})
```

Then the control law has the form

```math
\boldsymbol{u}_B = -k_p\big(\log(\boldsymbol{C}_{BI}^r\boldsymbol{C}_{BI})\big)^\vee - k_d\big(\boldsymbol{\omega}_B^{BI} - \boldsymbol{\omega}_B^{BI^r}\big)
```