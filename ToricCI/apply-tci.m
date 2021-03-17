load "tci.m";

FWPS := function(Q : T := [])
	x := Q[1];

	r := #x;

	Qrat := [];	

	error if #Q - #T gt 1, "wrong number of torsion";

	if #T gt 0 then
		// express torsion in Q by rational numbers
		for i := 1 to #T do
			Append(~Qrat, [q/T[i] : q in Q[i+1]]);
		end for;


		Z := FakeProjectiveSpace(Rationals(), x, Qrat);
	else
		Z := WeightedProjectiveSpace(Rationals(), x);
	end if;

	return Z;
end function;



RandomSystemFWPS := function(Q, M : T := [], N := 19, sign := true)
	x := Q[1];
	u := [mu[1] : mu in M];

	r := #x;

	Qrat := [];	

	error if #Q - #T gt 1, "wrong number of torsion";

	if #T gt 0 then
		// express torsion in Q by rational numbers
		for i := 1 to #T do
			Append(~Qrat, [q/T[i] : q in Q[i+1]]);
		end for;


		Z := FakeProjectiveSpace(Rationals(), x, Qrat);
	else
		Z := WeightedProjectiveSpace(Rationals(), x);
	end if;

	WDiv := DivisorGroup(Z);

	DL := [ u[i]*Divisor(Z, 1) : i in [1..#u] ];
	BL := [ MoveToOrthant(RiemannRochPolytope(D)) : D in DL ];

	R<[S]> := PolynomialAlgebra(Rationals(), Dimension(Ambient(Rep(BL))));

	// flag if coefficientes have sign
	Sign := sign select [-1, 1] else [1];

	F := [\
		&+[ Random(Sign)*Random(1,N)*Monomial(R, Eltseq(v)) : v in Vertices(B)]\
		: B in BL ];

	return F;
end function;


BuildNiceSystem := function(Q, M : T := [])
	counter := 0;
	repeat                         
		counter +:= 1;
		F := RandomSystemFWPS(Q, M: T := T,\
			N := 1 + (counter div 10), sign := (counter mod 10 ge 5));
	until IsNondegenerate(F);

	return F, counter;
end function;
//IsNondegenerate(F);

/* 
Habs geschafft, etwas Entartetes zu wuerfeln:

x := [1,2,2,2,3,3];
u := [6, 6];

Loading "apply-tci.m"
[
    Weil divisor with coefficients:
        6, 0, 0, 0, 0, 0,
    Weil divisor with coefficients:
        6, 0, 0, 0, 0, 0
]
[
    7*T[1]^3*T[2]^3*T[4]^3 - 6*T[1]^3*T[3]^3*T[4]^3 - 7*T[1]^3*T[4]^3 +
        8*T[1]^2*T[2]*T[3]*T[4]^4*T[5]^2 - T[1]^2*T[2]*T[3]*T[4]^4 -
        4*T[2]^3*T[3]^3,
    4*T[1]^3*T[2]^3*T[4]^3 + 8*T[1]^3*T[3]^3*T[4]^3 - 4*T[1]^3*T[4]^3 +
        8*T[1]^2*T[2]*T[3]*T[4]^4*T[5]^2 + 7*T[1]^2*T[2]*T[3]*T[4]^4 -
        T[2]^3*T[3]^3
]
false 1-dimensional polytope with 2 vertices:
    (6, 6, 0, 6, 0),
    (6, 0, 0, 6, 0)
>
*/
