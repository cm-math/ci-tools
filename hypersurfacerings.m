freeze;

import "utilities.m": SeqSupp, TransposeInSeq;
import "basics.m": MakePrimitive;


//-----------------------------------------------------------------------------
// SearchPrimeBinomial
//-----------------------------------------------------------------------------
// Tests if T_i defines a prime element in S/<g> for general g of degree mu
// by looking for prime binomials.
// If true, also returns exponents of the found prime binomial.
//
// Input.
// Q: Seq[Seq], sequences are rows of Q seen as matrix
// mu: relation degree
// i: index of variable T_i 
//
// Parameters.
// T: Seq[RngIntElt], torsion sequence
// 
// Ouput.
// BoolElt, Seq[Seq]
//
intrinsic SearchPrimeBinomial(Q::SeqEnum, mu::SeqEnum[RngIntElt], i : T := [])\
-> BoolElt, SeqEnum
{ Returns true if and only if there is a mu-homogeneous prime
  binomial not depending on i-th variable. }
	r := #Rep(Q);
	Q := Matrix(Q);

	// monomials of degree mu not depending on T_i
	Qsub := [Eltseq(w) : w in Rows(RemoveColumn(Q, i))];
	B := FiberPoints(Qsub, mu : T := T);

	// run through all mu-homogeneous binomials
	for Iset in Subsets({1..#B}, 2) do
		I := Setseq(Iset);
		nu1 := B[I[1]];
		nu2 := B[I[2]]; 

		// Apply primeness criterion for binomials
		ExpDiff := [ nu2[j] - nu1[j] : j in [1..r-1] ];
		if IsDisjoint(Seqset(SeqSupp(nu1)), Seqset(SeqSupp(nu2)))\
		   and GCD(ExpDiff) eq 1\
		then
			// adjust nu1, nu2 to original variables
			Insert(~nu1, i, 0);
			Insert(~nu2, i, 0);

			return true, [nu2, nu1];
		end if;
	end for;

	return false, _;
end intrinsic;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// SearchPrimeBinomials
// ----------------------------------------------------------------------------
// Tests if variables define prime elements in S/<g> for general g of degree mu
// by looking for suitable prime binomials.
//
// Input.
// Q: Seq[Seq], sequences are rows of Q seen as matrix
// mu: relation degree
// 
// Ouput.
// BoolElt, Seq[Seq[Seq]]]
//
intrinsic SearchPrimeBinomials(Q::SeqEnum, mu::SeqEnum[RngIntElt] : T := [])
-> BoolElt, SeqEnum
{ Returns true if and only if for each i there is a mu-homogeneous
  prime binomial not depending on i-th variable.
  If true, returns also a list of suitable prime binomials. }
	r := #Rep(Q);
	W := [Eltseq(w) : w in Rows(Transpose(Matrix(Q)))];
	Wset := Seqset(W);

	// prepare Seq for verifying binomials
	verEx := [ [ [-1 : j in [1..r]], [-1 : j in [1..r]] ] : i in [1..r] ];

	for w in Wset do
		i := Index(W, w);
		isPrime, primeBinomial := SearchPrimeBinomial(Q, mu, i : T := T);

		if not isPrime then
			return false, _;
		end if;

		// generate verifying binomial for each T_k
		// with deg(T_k) = deg(T_i).
		// nu1, nu2 might depend on T_k when i ne k.
		// In this case: swap roles of T_i and T_k.
		for k in [k : k in [1..r] | W[k] eq w] do
			verEx[k] := [TransposeInSeq(primeBinomial[1], i, k), TransposeInSeq(primeBinomial[2], i, k)]; 
		end for;
	end for;

	return true, verEx;
end intrinsic;
//-----------------------------------------------------------------------------



intrinsic DimHomComp(Q::SeqEnum, mu::SeqEnum, w::SeqEnum : T := [])
-> RngIntElt
{ Dimension of homogeneous component R_w of R. }
	v := [w[i] - mu[i] : i in [1..#w]];

	n := #FiberPoints(Q, w : T := T);
	m := #FiberPoints(Q, v : T := T);

	return n-m;
end intrinsic;



intrinsic GeneratorDegreeDimensionTuple(Q::SeqEnum, mu::SeqEnum : T := [])
-> SeqEnum
{ Returns the generator degree dimension tuple of R. }
	return Sort([DimHomComp(Q, mu, Eltseq(w) : T := [])\
	       : w in Seqset(Rows(Transpose(Matrix(Q))))]);
end intrinsic;



intrinsic AnticanonicalClass(Q::SeqEnum, Mu::SeqEnum : T := [])
-> SeqEnum[RngIntElt]
{Anticanonical class provided SD defines a complete intersection.}
	w := Eltseq(&+Rows(Transpose(Matrix(Q))) - &+Rows(Matrix(Mu)));

	for i->t in T do
		w[#w - #T + i] mod:= t;
	end for;

	return w;
end intrinsic;



intrinsic HilbertCoeffs(Q::SeqEnum, mu::SeqEnum, n::RngIntElt\
	  : T := [], useFundamentalSystem := false)
-> SeqEnum
{ Computes the first n coefficients of the Hilbert series of R. }
	require n gt 1: "n must be greater than one.";

	wak := AnticanonicalClass(Q, [mu] : T := T);	
	
	// primitive element might be not unique
	w0 := (useFundamentalSystem) select MakePrimitive(wak) else wak;

	coeffs := [DimHomComp(Q, mu, w : T := T) where w is [k*a : a in w0]\
	           : k in [1..n-1]];

	return [1] cat coeffs;
 end intrinsic;



//============================================================================= 
// Geometric tools for Picard number two 
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/* ToricRelevantFaces
 *-----------------------------------------------------------------------------
 * Compute faces gamma_0 of gamma where u lies in the relative interior
 * of Q(gamma_0) in case of Picard number one or two.
 * 
 * Input.
 * Q: Seq[Seq], Matrix
 * u: Seq[RngIntElt], ample class
 * 
 * Ouput.
 * Seq[Seq], faces sorted by dimension
*/
ToricRelevantFaces := function(Q, u)
	error if #Q gt 2, "Picard number #Q must be one or two";

	// help function
	NegativeOriented := func<u,v | Sign(Determinant(Matrix([v,u])))>;

	// sort columns by orientation
	Q := Matrix(Q);
	r := Ncols(Q);
	k := Nrows(Q);
	Qcols := Rows(Transpose(Q));

	// Picard number one
	if k eq 1 then
		return [ Setseq(I) : I in Subsets({1..r}) | I ne {} ];
	// Picar number two
	else
		Sort(~Qcols, NegativeOriented, ~p);	
		u := Vector(u);

		p := Eltseq(p);
		i1 := Max([ i : i in {1..r} | NegativeOriented(u, Qcols[i]) ge 0]);
		i2 := Min([ i : i in {i1..r} | NegativeOriented(Qcols[i], u) ge 0 ]);

		return Sort([ p[Setseq(I join J)]\
               : I in Subsets({1..i1}), J in Subsets({i2..r})\
               | not (IsEmpty(I) or IsEmpty(J)) ], func<A, B | Sign(#A - #B)>);
	end if;
end function;



//-----------------------------------------------------------------------------
/* RelevantFacesLS
 *-----------------------------------------------------------------------------
 * Compute releveant faces i.e.
 * toric relevant faces which are
 * F-faces.
 * The relation g is a general member
 * of the linear system given by M,
 * i.e. the entries for M are precisely
 * the exponents of the monomials of g
 * with non-zero coefficient.
 *
 * Restricted to Picard number 2
 * since ToricRelevantFaces is so.
 * 
 * Input.
 * Q: Seq[Seq], Matrix
 * mu: Seq[RngIntElt], relation degree
 * u: Seq[RngIntElt], ample class
 * 
 * Ouput.
 * Seq[Seq], faces sorted by dimension
*/
RelevantFacesLS := function(Q, M, u : T := [])
	Q0 := Q[1..#Q - #T];
	u0 := u[1..#Q - #T];
	torrlv := ToricRelevantFaces(Q0, u0);
	return [ F : F in torrlv | #I ne 1\
           where I is {w : w in M | Seqset(SeqSupp(w)) subset Seqset(F)} ];
end function;



//-----------------------------------------------------------------------------
// RelevantFaces
//-----------------------------------------------------------------------------
// Compute releveant faces i.e. toric relevant faces which are
// barX-faces. Restricted to Picard number 2 since ToricRelevantFaces is so.
// 
// Input.
// Q: Seq[Seq], Matrix
// mu: Seq[RngIntElt], relation degree
// u: Seq[RngIntElt], ample class
// 
// Ouput.
// Seq[Seq], faces sorted by dimension
//
RelevantFaces := function(Q, mu, u : T := [])
	M := FiberPoints(Q, mu : T := T);
	return RelevantFacesLS(Q, M, u : T := T);
end function;



//-----------------------------------------------------------------------------
/* IsMinimalAmbSmooth
 *-----------------------------------------------------------------------------
 * Check if the minimal ambient toric variety of any hypersurface defined by
 * a polynomial with precisely the monomials M is smooth.
 * 
 * Input.
 * Q: Seq[Seq], degree matrix as row sequence
 * M: Seq[Seq], exponent vectors of monomials
 * u: Seq[RngIntElt], ample class
 * 
 * Ouput.
 * BoolElt
*/
intrinsic IsMinimalAmbSmooth(Q::SeqEnum[SeqEnum[RngIntElt]],\
                               M::SeqEnum[SeqEnum[RngIntElt]],\
                               u::SeqEnum[RngIntElt]\
                               :\
                               T := []) -> BoolElt
{ Returns true iff the minimal ambient toric variety of any polynomial
  with monomials M is smooth.
  Only available for Picard number <= 2. }
	QM := Matrix(Q);
	k := Nrows(QM);
	r := Ncols(QM);	

	require k - #T le 2: "Picard number must be less or equal than 2";

	for F in RelevantFacesLS(Q, M, u : T := T) do 
		A := Submatrix(QM, [1..k], F);
		
		if not IsZZGenerating([Eltseq(w) : w in Rows(Transpose(A))] : T := T) then
			return false, F;
		end if;
	end for;

	return true;
end intrinsic;



intrinsic IsMuAmbientSmooth(Q, mu, u : T := []) -> BoolElt
{ Returns true iff mu-minimal ambient toric variety\
  specified by (Q, mu, u) is smooth.\
  Only available for Picard number <= 2. }
	M := FiberPoints(Q, mu : T := T);
	return IsMinimalAmbSmooth(Q, M, u : T := T);
end intrinsic;



intrinsic QuasismoothTest(Q, mu, u : T := []) -> BoolElt
{ Returns true iff relation degree of a variety X with hypersurface
  Cox ring and given specifying data satisfies necessary\
  conditions for X being quasismooth.
  Only available for Picard number <= 2. }
	B := FiberPoints(Q, mu : T := T);
	rlv := RelevantFaces(Q, mu, u : T := T);

	for F in rlv do 
		// check if mu lies in Q(gamma_0 \cap ZZ^r)
		found := false;

		for v in B do
			if Seqset(SeqSupp(v)) subset Seqset(F) then
				found := true;
				break;
			end if;
		end for;

		if found then
			continue;
		end if;


		// check if mu lies in w + Q(gamma_0 \cap ZZ^r)
		// for some Cox ring generator degree w
		for w in [Eltseq(w) : w in Rows(Transpose(Matrix(Q)))] do
			diffmuw := [mu[i] - w[i] : i in [1..#mu]];
			B2 := FiberPoints(Q, diffmuw : T := T);

			for v in B2 do
				if Seqset(SeqSupp(v)) subset Seqset(F) then
					found := true;
					break;
				end if;
			end for;

			if found then
				break;
			end if;
		end for;

		if found then
			continue;
		else
			return false, F;
		end if;
	end for;

	return true, [];
end intrinsic;
