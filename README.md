# SimpleRotations.jl

A minimal package for converting between rotation representations, generating random rotations, projecting to SO(3). To use in Julia, run:

```Julia
pkg> add https://github.com/lopenguin/SimpleRotations.jl
```

## Representations
- **Matrix**: Any 3x3 matrix (Julia's built in type) may be treated as a rotation matrix. The matrix $R\in\mathrm{SO}(3)$ if it is orthogonal ($R'R=RR'=I$) and has determinant $+1$.
- **Axis-Angle**: any axis `ω` (3-vector, will be normalized) and angle `θ` (float). 
- **Quaternion**: any 4-vector `[qw, qx, qy, qz]` with `qw` (the first element) treated as the scalar part.

## Conversion Table
Currently only conversions to and from rotation matrices are implemented.

|            | Quaternion     | Axis-Angle          | Matrix          |
|------------|----------------|---------------------|-----------------|
| Quaternion |                |                     | `rotm2quat(R)`  |
| Axis-Angle |                |                     | `rotm2axang(R)` |
| Matrix     | `quat2rotm(q)` | `axang2rotm (ω, θ)` |                 |

## Other functions
- `randrotation()` generates a uniformly random rotation matrix.
- `project2SO3(M)` projects a 3x3 matrix `M` to $\mathrm{SO}(3)$ using singular value decomposition.
- `roterror(R₁, R₂)` computes the angular error in degrees between two rotation matrices