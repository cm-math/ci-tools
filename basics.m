freeze;

//============================================================================= 
/* ZZIndex
 * ----------------------------------------------------------------------------
 * Compute the index of the subgroup spanned by a sequence
 * of elements in ZZ^k x (ZZ/t_1*ZZ) x ... (ZZ/t_q*ZZ) 
 * 
 * Input.
 * W: Seq[Seq]
 * T: Seq[RngIntElt]
 * 
 * Ouput.
 * RngIntElt - index of subgroup if finite, zero otherwise
*/
ZZIndex := function(W : T := [])
	r := #Rep(W);
	q := #T;

	error if q gt r, "#T must be less or equal than k";

	// Append (i-th entry of T) * e_(i-th torsion coordinate) to W
	for i := r - q + 1 to r do
		w := [ (i eq j) select T[i - r + q] else 0 : j in [1..r] ];
		Append(~W, w);
	end for;

	D := ElementaryDivisors(Matrix(W));
	
	// index of <W> is finite
	if #D eq r then
		return &*D;
	else
		return 0;
	end if;
end function;
//=============================================================================



//============================================================================= 
/* IsZZGenerating
 * ----------------------------------------------------------------------------
 * Check if a given sequence of elements
 * is a generating set of ZZ^k x (ZZ/m_1*ZZ) x ... (ZZ/m_t*ZZ) 
 * 
 * Input.
 * W: Seq[Seq]
 *
 * Optional input.
 * M: Seq[RngIntElt]
 * 
 * Ouput.
 * BoolElt
*/
intrinsic IsZZGenerating(W : T := []) -> BoolElt
{ Returns true iff W forms a generating set. }
	k := #Rep(W);
	q := #T;

	require q le k: "#T must be less or equal than k!";

	// Append (i-th entry of T) * e_(i-th torsion coordinate) to W
	for i := k - q + 1 to k do
		w := [ (i eq j) select T[i - k + q] else 0 : j in [1..k] ];
		Append(~W, w);
	end for;

	A := Transpose(Matrix(W));
	
	return SmithForm(A) eq HorizontalJoin(\
		IdentityMatrix(Integers(), k),\
		ZeroMatrix(Integers(), k, #W-k));
end intrinsic;
//============================================================================= 



//============================================================================= 
/* IsNwiseGenerating
 * ----------------------------------------------------------------------------
 * Check if each n columns of Q
 * form a generating set of ZZ^k x (ZZ/m_1*ZZ) x ... (ZZ/m_t*ZZ) 
 * 
 * Input.
 * Q: Seq[Seq] - matrix as a row sequence
 * n: RngIntElt
 * 
 * Optional input.
 * M: Seq[RngIntElt]
 * 
 * Ouput.
 * BoolElt
*/
IsNwiseGenerating := function(Q, n : M := [])
	W := [Eltseq(w) : w in Rows(Transpose(Matrix(Q)))]; 

	for I in Subsets({1..#W}, n) do                      
		if not IsZZGenerating(W[Setseq(I)] : T := M) then
			return false;
		end if;
	end for;

	return true;
end function;
//============================================================================= 



//============================================================================= 
/* IsAlmostFree
 * ----------------------------------------------------------------------------
 * Check if each n-1 columns of Q
 * form a generating set of ZZ^k x (ZZ/t_1*ZZ) x ... (ZZ/t_q*ZZ) 
 * 
 * Input.
 * Q: Seq[Seq] - matrix as a row sequence
 * 
 * Optional input.
 * T: Seq[RngIntElt]
 * 
 * Ouput.
 * BoolElt
*/
intrinsic IsAlmostFree(Q : T := []) -> BoolElt
{ Returns true iff and only if each n-1 columns of Q form a generating set. }
    r := #Rep(Q);
	return IsNwiseGenerating(Q, r-1 : T := T);
end intrinsic;
//============================================================================= 



//============================================================================= 
/* IsNwiseCoprime
 * ----------------------------------------------------------------------------
 * Returns true if and only if each n
 * elements of A are coprime.
 * 
 * Input.
 * A: Seq[RngIntElt]
 * n: RngIntElt
 * 
 * Ouput.
 * BoolElt
*/
IsNwiseCoprime := function(A, n)
	for I in Subsets({1..#A}, n) do
		if GCD(A[Setseq(I)]) ne 1 then
			return false;
		end if;
	end for;

	return true;
end function;
//============================================================================= 



//============================================================================= 
/* IsCollinear
 * ----------------------------------------------------------------------------
 * Check if two vectors v, w
 * are collinear.
 * 
 * Input.
 * v: Seq[RngIntElt]
 * w: Seq[RngIntElt]
 * 
 * Ouput.
 * BoolElt
*/
IsCollinear := function(v,w)
	v := Vector(Rationals(), v);
	w := Vector(Rationals(), w);

	n := #ElementToSequence(v);

	I := {i : i in [1..n] | v[i] ne 0};	
	if #I eq 0 then
		return true;
	end if;

	i := Rep(I);
	
	u := [w[j]-w[i]/v[i]*v[j] : j in [1..n]];
	return u eq ZeroSequence(Rationals(), n);
end function;
//============================================================================= 



//============================================================================= 
// IsDivisibleByMod
//-----------------------------------------------------------------------------
// Test if x is divisible by k modulo n.
// 
// Input.
// x, k, n: RngIntElt
// 
// Ouput.
// BoolElt
//
IsDivisibleByMod := function(x, k, n)
	// k | x in Z/nZ means n | k*y - x for some 0 <= a < n.
	for y := 0 to n - 1 do
		if IsDivisibleBy(k*y - x, n) then
			return true, y;
		end if;
	end for;

	return false, _ ;
end function;
//============================================================================= 



//============================================================================= 
// MakePrimitive
//-----------------------------------------------------------------------------
// Given an element w of a f.g. abelian group K = Z^k x G,
// returns the largest integer k such that w = k*v holds for some v in K.
// 
// Input.
// w: SeqEnum[RngIntElt]
//
// Paramters.
// T: SeqEnum[RngIntElt]
// 
// Ouput.
// SeqEnum, v
// BoolElt, k
//
MakePrimitive := function(w : T := [])
	w0 := [w[i] : i in [1..#w-#T]];
	wtor := w[#w-#T+1 .. #w];

	error if #w0 eq 0 or w0 eq [0 : i in [1..#w0]],\
		"free part of w must not be empty nor zero";

	for k in Reverse(Divisors(GCD(w0))) do
		v := [x div k : x in w0];

		for i -> t in T do
			isDiv, y := IsDivisibleByMod(wtor[i], k, t);
			if isDiv then
				Append(~v, y);
			else
				break;
			end if;
		end for;

		if #v eq #w then
			return v, k;
		end if;
	end for;
end function;


//============================================================================= 
/* FiberPoly
 * ----------------------------------------------------------------------------
 * Input.
 * Q: Seq[Seq], Matrix
 * w: Seq[RngIntElt] 
 * 
 * Ouput.
 * Polyhedron, intersection of fiber Q^{-1}(w) with positive orthant
*/
intrinsic FiberPoly(Q::SeqEnum[SeqEnum[RngIntElt]],
                   w::SeqEnum[RngIntElt]
                   : T := []) -> TorPol
{Intersection of fiber Q^{-1\}(w) with positive orthant.}
	Q := Matrix(Q);

	// only consider free parts
	if #T gt 0 then
		Q := Submatrix(Q, [1..Nrows(Q) - #T], [1..Ncols(Q)]);
		w := w[1..Nrows(Q)];
	end if;

	E := ToricLattice(Ncols(Q));
	K := ToricLattice(Nrows(Q));

	Qmap := LatticeMap(E, K, Transpose(Q));	
	w := K!w;

	// Kernel of Q as cone
	M := Cone(KernelBasis(Qmap));
	M := M + (-M);

	// If fiber of w is empty, then Preimage throws an error.
	// We avoid this by using exception handling
	// and let the function return an empty set.
	try
		fib := Preimage(Qmap, w) + ConeToPolyhedron(M);
	catch e
		return [];
	end try;

	return fib meet ConeToPolyhedron(PositiveQuadrant(E));
end intrinsic;
//============================================================================= 



//============================================================================= 
/* FiberPoints
 * ----------------------------------------------------------------------------
 * Computes the set of integer points of
 * the fiber polytope of w.
 * Note that these points correspond to
 * the monomials of degree w.
 * 
 * Input.
 * A: Seq[Seq], Matrix
 * w: Seq[RngIntElt] 
 * 
 * Ouput.
 * Seq
*/
intrinsic FiberPoints(Q::SeqEnum[SeqEnum[RngIntElt]],
                   w::SeqEnum[RngIntElt]
                   : T := []) -> SeqEnum
{Returns sequence of all points a with Q(a) = w.}
	// ensure matrix Q is represented as row sequence
	Q := [Eltseq(w) : w in Rows(Matrix(Q))];

	B0 := [ Eltseq(v) : v in Points( FiberPoly(Q, w : T := T) )];
	
	if IsEmpty(T) then
		return B0;
	else
		r := #Rep(Q);
		B := [ v : v in B0 | &and[\
			(&+[Q[i][j]*v[j] : j in [1..r]] - w[i]) mod T[i - #Q + #T] eq 0\
			: i in [#Q - #T + 1 .. #Q] ] ];
		return B;
	end if;
end intrinsic;
//============================================================================= 



//============================================================================= 
/* IsHomPermutation
 * ----------------------------------------------------------------------------
 * Checks if sequences of exponent vector (representing a set of monomials)
 * are equal up to a homogeneous permuation of variables w.r.t Q.

W = [w_1, ..., w_r]
 * 
 * Input.
 * Q: degree matrix as row sequence
 * E, F: sequences of exponent vectors to compare
 * 
 * Ouput.
 * BoolElt
 * If returns true, also returns according
 * permutation sequence Seq[RngIntElt]
*/

intrinsic IsHomPermutation(Q, E, F) -> BoolElt
{ check things }
	Q := Matrix(Q);
	r := NumberOfColumns(Q);

	W :=  [Eltseq(w) : w in Rows(Transpose(Q))];
	Wset := Seqset(W);

	// Sequence of index blocks of same degree
	I := [ [ i : i in [1..r] | W[i] eq Eltseq(w) ] : w in Wset];
	
	// Mapping sequence to I, i.e., i in I[N[i]]
	N := [];
	for i := 1 to r do
		for j := 1 to #I do
			if i in I[j] then
				Append(~N, j);
				break;
			end if;
		end for;
	end for;

	// go trough possible homogeneous permutations
	for P in CartesianProduct([Sym(#I[i]) : i in [1..#I]]) do
		// realize permutation as a sequence
		Pseq := [ I[k][Eltseq(P[k])[l]] : i in [1..r] | true \
			where l is Index(I[k], i) where k is N[i] ];

		// apply permutation
		Fpermuted := [v[Pseq] : v in F];

		if Seqset(E) eq Seqset(Fpermuted) then
			return true, Pseq;
		end if;
	end for;

	return false, _;
end intrinsic;
//=============================================================================
