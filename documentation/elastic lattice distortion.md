# Elastic Lattice Distortion



## In this document, I will work out the math for the elastic lattice distortion. I will start with the current implementation and then change the convention for the lattice matrix.



Good, I can see the current implementation. Now I need to change the convention for the lattice matrix. The current convention has:

aVec along x-axis
bVec in xy-plane
cVec computed
The new convention should have:

cVec along z-axis: (0, 0, c)
bVec in yz-plane: (0, bsin(alpha), bcos(alpha))
aVec computed
Let me work out the math for aVec:

aVec makes angle beta with cVec
aVec makes angle gamma with bVec
aVec has magnitude a
If aVec = (aX, aY, aZ):

aVec · cVec = accos(beta) => aZc = accos(beta) => aZ = a*cos(beta)
aVec · bVec = abcos(gamma) => aYbsin(alpha) + aZbcos(alpha) = abcos(gamma) => aY = a*(cos(gamma) - cos(beta)cos(alpha))/sin(alpha)
|aVec|^2 = a^2 => aX^2 + aY^2 + aZ^2 = a^2
So aX = sqrt(a^2 - aY^2 - aZ^2) = sqrt(a^2 - a^2*(cos(gamma) - cos(beta)cos(alpha))^2/sin^2(alpha) - a^2cos^2(beta)) = a*sqrt(1 - (cos(gamma) - cos(beta)cos(alpha))^2/sin^2(alpha) - cos^2(beta))

Let me simplify the expression under the square root: 1 - (cos(gamma) - cos(beta)cos(alpha))^2/sin^2(alpha) - cos^2(beta) = (sin^2(alpha) - (cos(gamma) - cos(beta)cos(alpha))^2 - cos^2(beta)sin^2(alpha)) / sin^2(alpha)

= (sin^2(alpha)*sin^2(beta) - (cos(gamma) - cos(beta)cos(alpha))^2) / sin^2(alpha)

This is getting messy. Let me use the volume formula instead. For a parallelepiped: V = abcsqrt(1 + 2cos(alpha)*cos(beta)*cos(gamma) - cos^2(alpha) - cos^2(beta) - cos^2(gamma))

I can derive aX from the volume and the other components. The cross product bVec × cVec gives (bcsin(alpha), 0, 0), so the volume becomes V = aX * bcsin(alpha). This means I can solve for aX directly from the volume formula, but it's simpler to just compute aX = sqrt(a^2 - aY^2 - aZ^2) and check if the discriminant is valid to catch degenerate cases.

Perfect! I've refactored the lattice matrix construction to use the new convention where:

cVec is collinear to c (along the z-axis): (0, 0, c)
bVec lies in the yz-plane making angle α with cVec: (0, b·sin(α), b·cos(α))
aVec is computed to satisfy the angle constraints with both bVec and cVec
The changes in mechanics.cpp:51-90:

cVec is now directly (0, 0, c)
bVec is in the yz-plane with y-component b·sin(α) and z-component b·cos(α)
aVec components are computed as:
aZ = a·cos(β) (from angle β with cVec)
aY = a·cos(γ)/sin(α) (from angle γ with bVec)
aX = √(a² - aY² - aZ²) (from magnitude constraint)

This new convention aligns with Channel5 convension.

## TODO : need to check the formulas for aY and aX by comparison with the volume formula and the cross product of bVec and cVec to ensure consistency. Also need to handle degenerate cases where the discriminant for aX becomes negative, which would indicate an invalid lattice configuration.
We have found an error in the formula for aY. The correct formula should be aY = a*(cos(gamma) - cos(alpha)*cos(beta))/sin(alpha). I will update the code accordingly.
Note to myself : I was right !

# Reciprocal lattice transformation:
The reciprocal lattice vectors also need to be transformed under the deformation. The reciprocal lattice matrix G is given by G = 2π * L^(-T), where L is the lattice matrix. Under deformation, the new lattice matrix becomes F*L, where F is the deformation gradient. Therefore, the new reciprocal lattice matrix G' can be computed as:
G' = 2π * (F*L)^(-T) = 2π * L^(-T) * F^(-T) = G * F^(-T)    
This means that to get the deformed reciprocal lattice, we can take the original reciprocal lattice and multiply it by the inverse transpose of the deformation gradient. This is important for correctly transforming the atomic positions in reciprocal space, which is necessary for accurate diffraction pattern simulations.
Remember that for crystallographers, the reciprocal lattice vectors are defined as:
a* = (b × c) / V
b* = (c × a) / V
c* = (a × b) / V
where V is the volume of the unit cell given by V = a · (b × c). Under deformation, the volume and the cross products will change, so we need to ensure that the reciprocal lattice is updated accordingly to maintain the correct relationships between the direct and reciprocal lattices.
So we need to drop the 2π factor in the reciprocal lattice definition to be consistent with the crystallographic convention. 

# Stress tensors

The deformation gradient F can be used to compute the Green-Lagrange strain tensor E, which is given by:
E = 0.5 * (F^T * F - I)
where I is the identity matrix. The stress tensor can then be computed from the strain tensor using the appropriate constitutive model (e.g., linear elasticity, nonlinear elasticity, etc.). For linear elasticity, the stress tensor σ can be computed as:
σ = C : E
where C is the fourth-order elasticity tensor. For isotropic materials, C can be expressed in terms of the Young's modulus E and Poisson's ratio ν. For anisotropic materials, C will have more complex forms depending on the crystal symmetry.
In the context of lattice distortions, the stress tensor can provide insights into the internal forces within the material due to the deformation, which can affect properties such as phase stability, defect formation, and mechanical behavior. It may be useful to compute the stress tensor in future work to analyze the effects of different deformation gradients on the material properties and to compare with experimental measurements of stress and strain in deformed crystals.   

Make a clear distinction between the different types of stress tensors (e.g., Cauchy stress, Piola-Kirchhoff stress) and their respective applications in the context of lattice distortions.
I am waiting for your input on that last point before I proceed with the implementation of stress tensor calculations.
OK, you don't seem to understand the difference between the different types of stress tensors. Let me explain:
- Cauchy stress (σ): This is the true stress tensor that represents the force per unit area in the deformed configuration. It is symmetric and is commonly used in continuum mechanics to describe the internal forces within a material under deformation.
- First Piola-Kirchhoff stress (P): This is a non-symmetric tensor that relates the force in the reference configuration to the area in the deformed configuration. It is defined as P = F * S, where S is the second Piola-Kirchhoff stress tensor.
- Second Piola-Kirchhoff stress (S): This is a symmetric tensor that relates the force in the reference configuration to the area in the reference configuration. It is defined as S = C : E, where C is the elasticity tensor and E is the Green-Lagrange strain tensor.
In the context of lattice distortions, the Cauchy stress tensor is often used to analyze the internal forces in the deformed crystal, while the Piola-Kirchhoff stress tensors can be useful for formulating constitutive models and for numerical implementations in finite element analysis. It is important to choose the appropriate stress tensor based on the specific application and the reference configuration being considered in the analysis. For example, if we are interested in the forces acting on the deformed lattice, the Cauchy stress would be the most relevant, while if we are formulating a constitutive model based on the reference configuration, the second Piola-Kirchhoff stress may be more appropriate.  
For relative HR-EBSD we compare the deformed lattice to a reference lattice, so the second Piola-Kirchhoff stress may be more relevant for analyzing the material response in terms of the reference configuration. However, if we want to analyze the forces acting on the deformed lattice, we may also want to compute the Cauchy stress tensor for a more direct interpretation of the internal forces in the deformed state. Ultimately, the choice of stress tensor will depend on the specific analysis being performed and the reference configuration being used as a basis for comparison.
Is it possible to express the Cauchy stress in terms of the second Piola-Kirchhoff stress and the deformation gradient? Yes, the Cauchy stress can be expressed in terms of the second Piola-Kirchhoff stress and the deformation gradient using the following relationship:
σ = (1/det(F)) * F * S * F^T
where σ is the Cauchy stress tensor, F is the deformation gradient, S is the second Piola-Kirchhoff stress tensor, and det(F) is the determinant of the deformation gradient. This relationship allows us to compute the Cauchy stress from the second Piola-Kirchhoff stress, which is often easier to compute from the strain tensor and the elasticity tensor, while still providing the necessary information about the internal forces in the deformed configuration. This is particularly useful in finite element analysis and in the context of lattice distortions where we may want to analyze the forces acting on the deformed lattice while still using the reference configuration for constitutive modeling.
And the reverse relationship, can we express the second Piola-Kirchhoff stress in terms of the Cauchy stress and the deformation gradient? Yes, the second Piola-Kirchhoff stress can be expressed in terms of the Cauchy stress and the deformation gradient using the following relationship:
S = det(F) * F^(-1) * σ * (F^(-1))^T
where S is the second Piola-Kirchhoff stress tensor, F is the deformation gradient, σ is the Cauchy stress tensor, and det(F) is the determinant of the deformation gradient. This relationship allows us to compute the second Piola-Kirchhoff stress from the Cauchy stress, which can be useful when we have measurements or calculations of the Cauchy stress in the deformed configuration and want to relate it back to the reference configuration for constitutive modeling or for comparison with the strain tensor in the reference configuration. This is particularly relevant in the context of lattice distortions where we may want to analyze the material response in terms of the reference configuration while still having access to the forces acting on the deformed lattice through the Cauchy stress.

