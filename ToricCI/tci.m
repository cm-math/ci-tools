freeze;

//============================================================================= 
/* MoveToOrthant
 * ----------------------------------------------------------------------------
 * Translates polytope B into positive orthant
 * 
 * Input.
 * B: TorPol
 * 
 * Ouput.
 * TorPol
*/

MoveToOrthant := function(B)
	n := Dimension(Ambient(B));
	V := Vertices(B);

	v := [ -Min([0] cat [v.i : v in V]) : i in [1..n] ];
	v := Ambient(B)!v;

	return B + v, v;
end function;
//============================================================================= 


//============================================================================= 
/* MinkowskiDecomposition
 * ----------------------------------------------------------------------------
 * For a Minkowski sum B = B_1 + ... B_s and
 * and face B' of B computes the unique decomposition
 * B' = B_1' + ... + B_s' with faces B_i' of B_i.
 * 
 * Input.
 * BL: Seq[TorPol], sequence with entries B_1, ..., B_s
 * F: TorPol, face of B
 * 
 * Ouput.
 * Seq[TorPol], sequence with entries B_1', ..., B_s'
*/

MinkowskiDecomposition := function(BL, F)
	B := &+BL;

	// IsFace(B, F) is not used due to mysterious side effects	
	error if not F in Faces(B, Dimension(F)), "F ist not a face of B";

	// interior point of normal cone
	N := Dual(Cone([u-v : u in Vertices(B), v in Vertices(F)]));
	h := &+RGenerators(N);

	F0L := [];

	for B0 in BL do
		V := Vertices(B0);
		a := Min([h*v : v in V]);
		F0 := Polytope([v : v in V | h*v eq a]);
		Append(~F0L, F0);
	end for;

	return F0L;
end function;
//============================================================================= 



//============================================================================= 
/* IsSmoothTorus
 * ----------------------------------------------------------------------------
 * Returns true if and only if the differential DF(z) is of full
 * rank at every point of V(F) inside the torus.
 * 
 * Input.
 * F: Seq[RngMPolElt]
 * 
 * Ouput.
 * BoolElt
*/

IsSmoothTorus := function(F)
	S := Parent(Rep(F));
	I := Ideal(F cat Minors(JacobianMatrix(F), #F));
	return IsInRadical(&*[S.i : i in [1..Rank(S)]], I);
end function;
//============================================================================= 


intrinsic NewtonPolytope(F::[RngMPolElt]) -> TorPol
{ Returns the Newton polytope of F i.e. the Minkowski sum of all
  Newton polytopes of members of F. }
	S := Parent(Rep(F));
	Q := RationalFunctionField(BaseRing(S), Rank(S));

	BL := [NewtonPolytope(Q!f) : f in F];

	return &+BL;
end intrinsic;


//============================================================================= 
/* FacePolynomial
 * ----------------------------------------------------------------------------
 * Computes the polynomial associated with a face
 * of the Newton polytope of f.
 * 
 * Input.
 * f: RngMPolElt
 * B0: TorPol, face of Newton polytope of f
 * 
 * Ouput.
 * RngMPolElt
*/

intrinsic FacePolynomial(f::RngMPolElt, B0::TorPol) -> RngMPolElt
{ Returns the polynomial associated with a face of the Newton polytope of f. }
	S := Parent(f);
	return &+[ MonomialCoefficient(f, Eltseq(v)) * Monomial(S, Eltseq(v)) : v in Points(B0) ];
end intrinsic;
//============================================================================= 



//============================================================================= 
/* FacePolynomial
 * ----------------------------------------------------------------------------
 * Computes the system of polynomials associated with a face
 * of the Minkowski sum of the Newton polytops of f_1, ..., f_s
 * 
 * Input.
 * F: Seq[RngMPolElt], sequence with entries f_1, ..., f_s
 * B0: TorPol, face of Minkowski sum of Newton polytops
 * 
 * Ouput.
 * Seq[RngMPolElt]
*/

intrinsic FaceSystem(F::[RngMPolElt], B0::TorPol) -> SeqEnum[RngMPol]
{ Returns the system associated with a face of the Newton polytope of f. }
	S := Parent(Rep(F));
	Q := RationalFunctionField(BaseRing(S), Rank(S));

	// BL = [B_1, ..., B_s]
	BL := [NewtonPolytope(Q!f) : f in F];

	B0L := MinkowskiDecomposition(BL, B0);

	return [FacePolynomial(F[i], B0L[i]) : i in [1..#F]];
end intrinsic;
//============================================================================= 



//============================================================================= 
/* IsNondegenerate
 * ----------------------------------------------------------------------------
 * Checks if a system of polynomials is non-degenerate
 * 
 * Input.
 * F: Seq[RngMPolElt], sequence with entries f_1, ..., f_s
 * 
 * Ouput.
 * BoolElt, TorPol
*/

intrinsic IsNondegenerate(F::[RngMPolElt]) -> BoolElt, TorPol
{ True iff F is non-degenerate.  If F is not non-degenerate,
  also returns a face of the Minkowski sum of the Newton polytopes
  for which the differential of the associated face system of F
  is not of rank s at every point point of the torus. }
	S := Parent(Rep(F));
	Q := RationalFunctionField(BaseRing(S), Rank(S));

	// BL = [B_1, ..., B_s]
	BL := [NewtonPolytope(Q!f) : f in F];

	B := &+BL;
	
	for B0 in &cat(Faces(B)) do
		if not IsSmoothTorus(FaceSystem(F, B0)) then
			return false, B0;
		end if;
	end for;

	return true, _;
end intrinsic;
//============================================================================= 
