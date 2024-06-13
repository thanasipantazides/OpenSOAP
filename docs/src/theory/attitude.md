# Attitude
Attitude means orientation of the spacecraft. There are many ways to define this in three dimensions: quaternions, axis-angle parameters, Euler angles, Direction Cosine Matrix/Rotation Matrix, and Rodrigues parameters. I choose Direction Cosine Matrix here.

## Direction Cosine Matrices
Say you have two right-handed, 3-dimensional reference frames, ``\mathcal{A}`` and ``\mathcal{B}``. Recall that a *vector* only specifies a direction, it is note located at a specific point in a coordinate system. Say you have a vector ``\boldsymbol{x}``. This direction exists regardless of the reference frame you *resolve* it in. You won't be able to write down specific numbers for each vector component without resolving it, but you can see that it points in a certain direction. Denote the resolution of ``\boldsymbol{x}`` in frame ``\mathcal{A}`` as ``\boldsymbol{x}_\mathcal{A}``. Then, the direction cosine matrix (DCM) transforming ``\boldsymbol{x}_\mathcal{A}`` into frame ``\mathcal{B}`` is ``\boldsymbol{C}_{\mathcal{B}\mathcal{A}}``:

```math
\boldsymbol{x}_\mathcal{B} = \boldsymbol{C}_{\mathcal{B}\mathcal{A}}\boldsymbol{x}_\mathcal{A}
```

Note how the subscripts for each frame are aligned left/right. 

All valid direction cosine matrices have the following properties:
```math
\begin{aligned}
    \boldsymbol{C}_{\mathcal{B}\mathcal{A}}^\mathsf{T}\boldsymbol{C}_{\mathcal{B}\mathcal{A}} &= \boldsymbol{1} \\
    \mathrm{det}(\boldsymbol{C}_{\mathcal{B}\mathcal{A}}) &= 1
\end{aligned}
```
The first condition effectively specifies that vector length is preserved under rotations, and the second condition restricts us to right-handed coordinate systems. Note that the first condition implies that the inverse of the DCM is equal to its transpose. Because this inverse operation transforms vectors in the opposite direction, the following are equivalent:
```math
\boldsymbol{C}_{\mathcal{B}\mathcal{A}}^\mathsf{T} = \boldsymbol{C}_{\mathcal{B}\mathcal{A}}^{-1} = \boldsymbol{C}_{\mathcal{A}\mathcal{B}}
```

The set of all ``3\times 3`` matrices that satisfies both these conditions is referred to as the special orthogonal group of dimension 3, or ``\mathrm{SO}(3)``:

```math
\mathrm{SO}(3) \equiv \{\boldsymbol{C}_{\mathcal{B}\mathcal{A}} \ | \ \boldsymbol{C}_{\mathcal{B}\mathcal{A}}^\mathsf{T}\boldsymbol{C}_{\mathcal{B}\mathcal{A}} = \boldsymbol{1}, \
    \mathrm{det}(\boldsymbol{C}_{\mathcal{B}\mathcal{A}}) = 1 \}
```
