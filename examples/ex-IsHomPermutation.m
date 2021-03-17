S := PolynomialAlgebra(Rationals(), 6);
Q := [[1,1,1,0,0,0], [0,0,0,1,1,1]];
g1 := S.1*S.4 + S.2*S.5 + S.3*S.6;
g2 := S.1*S.6 + S.2*S.4 + S.3*S.5;
E := [Exponents(f) : f in Monomials(g1)];
F := [Exponents(f) : f in Monomials(g2)];
IsHomPermutation(Q, E, F);

/* output
true [ 1, 2, 3, 6, 4, 5 ]
*/
