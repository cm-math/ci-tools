S<[T]> := PolynomialAlgebra(Rationals(), 5);
f := (T[1] - T[2])*(T[3] - T[4]) + T[5]^2;
ver, B0 := IsNondegenerate([f]);
ver;
// output: false

B0;
/* output:
2-dimensional polytope B0 with 4 vertices:
    (1, 0, 1, 0, 0),
    (1, 0, 0, 1, 0),
    (0, 1, 1, 0, 0),
    (0, 1, 0, 1, 0)
*/

FacePolynomial(f, B0);
/* output:
T[1]*T[3] - T[1]*T[4] - T[2]*T[3] + T[2]*T[4]
*/
