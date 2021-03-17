# ci-tools
This is a package for the computer algebra system [Magma](http://magma.maths.usyd.edu.au/magma/)
providing a toolkit for hypersurfaces and complete intersections in toric varieties.

## Usage
Go to the directory with the `spec` file and start a Magma session.
The provided intrinsics are available after loading the package with the following command.
```
AttachSpec("spec");
```

## Intrinsics
This package provides the following user-defined intrinsics.

* For graded hypersurface rings
    * `SearchPrimeBinomial`: search for homogeneous prime binomials
    * `DimHomComp`: compute vector space dimension of homogeneous components
    * `GeneratorDegreeDimensionTuple`: compute the dimension of the homogeneous components of the generator degrees
    * `HilbertCoeffs`: compute the first coefficients of the Hilbert series
    * `IsMuAmbientSmooth`: check if the mu-minimal ambient toric variety is smooth
    * `QuasismoothTest`: check necessary conditions for quasismoothness
* For polynomial systems
    * `FacePolynomial`: compute the polynomial associated with a face of the Newton polytope
    * `FaceSystem`: compute the system associated with a face of the Newton polytope
    * `IsNondegenerate`: check the non-degeneracy condition
* For toric varieties and complete intersection Cox rings
    * `ToricIntersectionProduct`: compute intersection numbers on Q-factorial projective toric varieties
    * `FanoDegree`: compute anticanonical self-intersection number of a Fano variety with complete intersection Cox ring
* Miscellanea
    * `IsZZGenerating`: check if given elements form a generating set of an abelian group
    * `FiberPoints`: compute the lattice points of the positive orthant that live in the fiber of a group homomorphism
    * `IsDolgachevPolytope`: check if a polytope is a Dolgachev polytope
    * `IsHomPermutation`: check if given subsets of a lattice are equal up to a certain permutation of coordinates

## Examples
Examples can be found in the `examples` folder.

## Note
This library is provided _as is_ and is still under development.
There is no warranty. Use it at your own risk.
If you find bugs or errors, please don't hesitate to contact me. 
