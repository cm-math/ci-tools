freeze;
import "basics.m": ZZIndex, IsCollinear;



//============================================================================= 
/* ToricIntersectionProduct
 * ----------------------------------------------------------------------------
 * Computes the intersection number of n classes D = [v_1, ..., v_n]
 * on an n-dimensional complete toric variety given by its
 * Cox ring generator degrees Q = [w_1, ..., w_r] and an ample class u.
 *
 * Input.
 * Q: Seq[Seq], sequences are rows of Q seen as matrix
 * u: ample class
 * D: Seq[Seq], sequences are intersecting classes
 *
 * Optional Input.
 * T: Seq[RngIntElt] order of cyclic factors
 * 
 * Ouput.
 * Intersection number of classes given
*/

intrinsic ToricIntersectionProduct(Q::[], u::[RngIntElt], D::[] : T := [])\
-> RngIntElt
{ Compute intersection number of divisor classes given by D }
	// help function
	coordinates := func<V, u |\
		Solution(Matrix(Rationals(), V), Vector(Rationals(), u))>;

	Q := Matrix(Q);
	r := Ncols(Q);
	rho := Nrows(Q) - #T;
	Qcols := Rows(Transpose(Q));

	require #D eq r - rho:\
		"Number of divisor classes must match dimension of toric variety";

	// discard torsion parts of D and u (if existing)
	D := [mu[1..rho] : mu in D];
	u := u[1..rho];

	// find primitive ray generators
	V := [];
	N := [];
	W := [];
	W0 := [];

	for wtup in Qcols do
		w := ElementToSequence(wtup);
		w0 := w[1..rho];
		wprim := [w0[i]/GCD(w0) : i in [1..#w0]];

		i := Position(V, wprim);
		if i eq 0  then
			Append(~V, wprim);
			Append(~N, 1);
			Append(~W, w);
		else
			N[i] +:= 1;
			Insert(~W, &+[N[j] : j in [1..i]], w);
		end if;
	end for;

	W0 := [w[1..rho] : w in W];

	// number of primitive ray generators of GIT-fan
	t := #V;

	// compute a first naive presentation of D_i as linear combination
	S := PolynomialAlgebra(Rationals(), t);
	f := Identity(S);
	for Y in D do
		c := coordinates(V, Y);
		f *:= &+[c[i]*S.i : i in [1..t]];
	end for;

	K := ToricLattice(rho);
	u := K!u;

	repeat
		for M in Exclude(Monomials(f), Identity(S)) do
			E := Exponents(M);	


			// can be computed
			if &and[ E[i] le N[i] : i in [1..t] ] then
				// expand monomials to find complementary weights
				I := {};
				for i in [1..t] do
					k := (i gt 1) select &+[N[j] : j in [1..i-1]] else 0;
					I join:= {k + l : l in [1..E[i]]};
				end for;

				J := {1..r} diff I;

				tau := Cone([K!W0[j] : j in J]);

				if IsInInterior(u, tau) then
					p := 1/(ZZIndex(W[Setseq(J)] : T := T) * \
						&*[GCD((W0[i])) : i in I]);
				else
					p := 0;
				end if;

				f +:= MonomialCoefficient(f, M)*(p - M);
			else
				I := [i : i in [1..t] | E[i] lt N[i]];
				J := [j : j in [1..t] | E[j] gt N[j]];
				j := Rep(J);
				d := E[j] - N[j];


				Mpr := 0;

				sg := Cone( [K!V[i] : i in I]);

				if IsInInterior(u, sg) then
					// present V[j] as linear combination over
					// generators of non over used rays
					c := coordinates(V[I], V[j]);

					E[j] -:= 1;
					Mpr := Monomial(S, E)*(&+[c[i]*S.I[i] : i in [1..#I]]);
				end if;

				f +:= MonomialCoefficient(f, M)*(Mpr - M);	
			end if;

		end for;
	until LeadingTotalDegree(f) le 0;	

	return LeadingCoefficient(f);
end intrinsic;
//============================================================================= 



intrinsic FanoDegree(Q::[], Mu::[] : T := []) -> SeqEnum[RngIntElt]
{Anticanonical self-intersection number provided\
 specifying data defines a complete intersection.}
	r := #Rep(Q);
	s := #Mu;
	rho := #Q - #T;

	require rho le 2: "currently works only for Picard number at most two :(";

	wantican := AnticanonicalClass(Q, Mu : T := T);

	// H := useFundamentalSystem select\
		[a/d : a in wantican | true where d is GCD(wantican)] else wantican;
	u := wantican;


	// -----------------------------------------------------------------------
	// choose ample class u1 of a suitable Q-factorial small quasimod
	// -----------------------------------------------------------------------
	if rho eq 1 then
		u1 := [1] cat [0 : i in [1..#T]];
	else
		W0 := Rows(Transpose(RowSubmatrixRange(Matrix(Q), 1, rho)));
		u0 := u[1..rho];

		// reorder columns of Q s.t. they are ordered counter-clockwise	
		orientation := func<u,v | Sign(Determinant(Matrix([v,u])))>;
		Sort(~W0, orientation);	

		J := [j : j in [1..#W0] | IsCollinear(W0[j], u0)];
		u10 := (IsEmpty(J)) select Vector(u0) else Vector(u0) + W0[Max(J) + 1];

		u1 := Eltseq(u10) cat [0 : i in [1..#T]];
	end if;
	// -----------------------------------------------------------------------

	return ToricIntersectionProduct(Q, u1, [u : k in [1..r-s-rho]] cat Mu : T := T);
end intrinsic;
