module SimpleRotations

using LinearAlgebra

export axang2rotm, rotm2axang, rotm2quat, quat2rotm, randrotation, unifrandrotation, project2SO3, project2SOd, roterror


"""
    axang2rotm(ω, θ)
    
Convert an axis `ω` and angle `θ` to a 3x3 rotation matrix.

Normalizes `ω` and returns identity when `θ = 0`.
"""
function axang2rotm(ω, θ)
    normalize!(ω)
    if θ == 0
        return diagm(ones(3))
    end
    K = [0. -ω[3] ω[2]; ω[3] 0. -ω[1]; -ω[2] ω[1] 0.]
    R = I + sin(θ)*K + (1 - cos(θ))*(K*K)
end

# this is to prevent errors when close to 1
@noinline my_slow_acos(x) = x ≈ 1 ? zero(x) : x ≈ -1 ? one(x)*π : acos(x)
my_acos(x) = abs(x) <= one(x) ? acos(x) : my_slow_acos(x)

"""
    rotm2axang(R)

Convert a rotation matrix `R` to axis-angle pair (`ω`, `θ`)
"""
function rotm2axang(R)
    θ = my_acos((tr(R) - 1) / 2)
    ω = [R[3,2] - R[2,3]; R[1,3] - R[3,1]; R[2,1] - R[1,2]] ./ (2*sin(θ))
    return ω, θ
end

"""
    rotm2quat(R)

Convert rotation matrix `R` to quaternion `[qw,qx,qy,qz]`.
"""
function rotm2quat(R)
    ω, θ = rotm2axang(R)
    q = [cos(θ/2); ω*sin(θ/2)]
end

"""
    quat2rotm(q)

Convert quaternion `[qw,qx,qy,qz]` to rotation matrix `R`.
"""
function quat2rotm(q)
    R = [q[1]^2+q[2]^2-q[3]^2-q[4]^2  2*(q[2]*q[3]-q[1]*q[4])  2*(q[2]*q[4]+q[1]*q[3]);
         2*(q[2]*q[3]+q[1]*q[4])  q[1]^2-q[2]^2+q[3]^2-q[4]^2  2*(q[3]*q[4]-q[1]*q[2]);
         2*(q[2]*q[4]-q[1]*q[3])  2*(q[3]*q[4]+q[1]*q[2])  q[1]^2-q[2]^2-q[3]^2+q[4]^2]
end

"""
    randrotation()
    
Generate uniformly random rotation matrix via quaternion sampling.
"""
function randrotation()
    q = randn(4)
    normalize!(q)
    # convert to rotation matrix
    R = quat2rotm(q)
end

"""
    project2SO3(M)
    
Project the 3x3 matrix `M` to SO(3) via SVD.
"""
function project2SO3(M)
    F = svd(M);
    R = F.U*F.V';
    if (det(R)<0)
        R = F.U * Diagonal([1.; 1; -1]) * F.V';  
    end
    return R
end

"""
    roterror(R₁, R₂)

Calculate angular difference (degrees) between `R₁` and `R₂`.

Essentially just a copy of `rotm2axang` without axis computation.
"""
function roterror(R₁, R₂)
    R = R₁'*R₂
    θ = my_acos((tr(R) - 1) / 2)
    # ω = [R[3,2] - R[2,3]; R[1,3] - R[3,1]; R[2,1] - R[1,2]] ./ (2*sin(θ))
    return θ*180/π
end


end # module SimpleRotations
