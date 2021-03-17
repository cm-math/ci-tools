//============================================================================= 
/* TupleSeq
 * ----------------------------------------------------------------------------
 * Given a tuple T with all entries of same type,
 * construct a sequence whose entries are the elements of T in the same order.
 * 
 * Input.
 * w: Tup
 * 
 * Ouput.
 * Seq
*/
TupleSeq := function(T)
	return [ T[i] : i in [1..#T] ];
end function;
//============================================================================= 



//============================================================================= 
// IsWellFormed
//-----------------------------------------------------------------------------
// Returns true if and only if each #x-1 elements of x are coprime.
// 
// Input.
// x: Seq[RngIntElt]
// 
// Ouput.
// BoolElt
//============================================================================= 
IsWellFormed := function(x)
	return &and[ Gcd(Remove(x,i)) eq 1 : i in [1..#x] ];
end function;
//============================================================================= 



//============================================================================= 
/* IsNwiseCoprime
 * ----------------------------------------------------------------------------
 * Returns true if and only if each n elements of x are coprime.
 * 
 * Input.
 * x: Seq[RngIntElt]
 * n: RngIntElt
 * 
 * Ouput.
 * BoolElt
*/
IsNwiseCoprime := function(x, n)
	for I in Subsets({1..#x}, n) do
		if GCD(x[Setseq(I)]) ne 1 then
			return false;
		end if;
	end for;

	return true;
end function;
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
/* TerminalConstraint
 * ----------------------------------------------------------------------------
 * Check if the condition from Lemma 1.7.3 (ii) is fullfilled for
 * each n-dimensional cone.
 * 
 * Input.
 * x: Seq[RngIntElt]
 * n: RngIntElt
 * 
 * Ouput.
 * BoolElt
*/
TerminalConstraint := function(x, n)
	for I in Subsets({1..#x}, n) do
		J := {1..#x} diff I;
		if GCD(x[Setseq(J)]) ge &+x[Setseq(I)]  then
			return false;
		end if;
	end for;

	return true;
end function;
//============================================================================= 



//============================================================================= 
/* TorsionBound
 * ----------------------------------------------------------------------------
 * Compute bound on order of cyclic factors of Cl(Z) as in Lemma 1.7.2.
 * 
 * Input.
 * x: Seq[RngIntElt]
 * u: RngIntElt
 *
 * Output.
 * RngIntElt
*/
TorsionBound := function(x, u)
	L := [u div n : n in x];
	return GCD([LCM(L[Setseq(I)]) : I in Subsets({2..#L}, #x-2)]);
end function;
//============================================================================= 



//============================================================================= 
/* AreGradingsIsomorphic
 * ----------------------------------------------------------------------------
 * Check if two rank one gradings are equivalent
 * up to certain changing coordinates and reordering variables
 * 
 * Input.
 * Q1: matrix as row sequence
 * Q2: matrix as row sequence
 *
 * Optional input.
 * T: Seq[RngIntElt] Orders of cyclic factors
 *
 * Output.
 * RngIntElt
*/
AreGradingsIsomorphic := function(Q1, Q2 : T := [])
	w1 := Q1[1];
	Q1tor := Q1[2..#Q1];
	w2 := Q2[1];
	Q2tor := Q2[2..#Q2];

	// free parts are not compatible
	if w1 ne w2 and w1 ne [-x : x in w2] then
		return false;
	end if;

	// here it is not checked if a row is multiplied by an unit
	// e.g. [ [3, 3, 3], [1,1,1] ] corresponds to
	// [ [3, 3, 3], [2, 2, 2] ] w.r.t K = Z x Z/3Z

	// run through admissible changes of torsion part
	for K in CartesianProduct([{0..T[i] - 1} : i in [1..#T]]) do
		Q1torpr := [];
		for i := 1 to #T do
			wtor := Q1tor[i];
			wtorpr := [(wtor[j] + K[i]*w1[j]) mod T[i] : j in [1..#w1]];
			Append(~Q1torpr, wtorpr);
		end for;

		if IsHomPermutation([w1], Q1torpr, Q2tor) then
			return true;
		end if;
	end for;

	return false;
end function;
//============================================================================= 



//============================================================================= 
// FilterGradings
// ----------------------------------------------------------------------------
// Remove (some) redundant entries from a sequence of ranke one gradings
// 
// Input.
// Qs: Seq[Seq[RngIntElt]],
//  sequence of degree matrices each given as row sequence
// T: Seq[RngIntElt], orders of finite cyclic factors
//
// Output.
// Seq[Seq[RngIntElt]]
//
FilterGradings := function(Qs, T)
	for i := 1 to #Qs do
		if not IsDefined(Qs, i) then
			continue;
		end if;

		for j := i+1 to #Qs do
			if not IsDefined(Qs, j) then
				continue;
			end if;

			if AreGradingsIsomorphic(Qs[i], Qs[j] : T := T) then
				Undefine(~Qs, j);
			end if;
		end for;
	end for;

	return [Qs[i] : i in [1..#Qs] | IsDefined(Qs, i)];
end function;
//============================================================================= 



//============================================================================= 
// IsACCterminal
//-----------------------------------------------------------------------------
// Returns true and only if the FWPS determined by the degree matrix Q
// has only terminal three-dimensional cones
// 
// Input.
// Q: Seq[RngIntElt]
//
// Optional input.
// T: Seq[RngIntElt] Orders of cyclic factors
// 
// Output.
// BoolElt
//
IsACCterminal := function(Q : T := [])
	Sigma := FanOfFakeWPS(Q, T);
	hasIntPoints := false;

	// go through 3-dimensional cones of Sigma
	for sigma in Cones(Sigma, 3) do
		Asigma := Polytope(Append(Rays(sigma), Zero(Ambient(sigma))));

		if NumberOfPoints(Asigma) gt #Vertices(Asigma) then
			hasIntPoints := true;
			break;
		end if;
	end for;

	return not hasIntPoints;
end function;
//============================================================================= 



//============================================================================= 
// ExtendsToBasepointfreeClass
// ----------------------------------------------------------------------------
// Returns true and only if the FWPS associated with Q (and T) admits a
// basepoint free class mu having u as its ZZ-part.
// If true, also returns mu.
// 
// Input.
// Q: Seq[RngIntElt]
// u: RngIntElt
// T: Seq[RngIntElt] Orders of cyclic factors
// 
// Output.
// BoolElt, [RngIntElt]
//
ExtendsToBasepointfreeClass := function(Q, u, T)
	r := #Rep(Q);
	q := #Q;

	// numbers l_i with l_i * x_i = u
	L := [u div Q[1][j] : j in [1..r]];

	// set of all l_i * w_i
	M := { [L[j]*Q[1][j]] cat [L[j]*Q[i][j] mod T[i-1]\
	     : i in [2..q]] : j in [1..r] };

	// if always l_i w_i = l_j w_j then mu = l_i w_i 
	// otherwise there is no mu as desired
	if #M eq 1 then
		return true, Rep(M);
	else
		return false, _;
	end if;
end function;
//============================================================================= 



//============================================================================= 
// CandidateRows
// ----------------------------------------------------------------------------
// Returns all rows representing elements in Z/tZ that might show up in
// the torsion part of degree matrix Q with Z-part x when mu has Z-part u
// Assumes x_0 = 1.
// 
// Input.
// x: Seq[RngIntElt]
// u: Seq[RngIntElt]
// t: RngIntElt Orders of cyclic factors
// 
// Output.
// SeqEnum
//
CandidateRows := function(x, u, t)
	r := #x;
	Qtors := [];

	// use x_0 = 1 to assume first column of Q is always 0
	for Y in CartesianPower({0..t-1}, r-1) do
		Qtor := [0] cat TupleSeq(Y);
		Q := [x, Qtor];

		if IsNwiseGenerating(Q, r-2 : M := [t]) and\
		   &and[ExtendsToBasepointfreeClass(Q, u[j], [t]) : j in [1..#u]]\
		then
			Append(~Qtors, Qtor);
		end if;
	end for;

	return Qtors;
end function;
//============================================================================= 



//============================================================================= 
// ProduceTorsionGradings
// ----------------------------------------------------------------------------
// Produces specifying data for Fano terminal hypersurfaces in fake weighted
// projective spaces Z arising from a non-degenerate Laurent polynomial
// with Cl(Z) having q finite cyclic factors.
// The Z-parts x and u = (u_1, ..., u_s) of the degree matrix resp.
// the relation degree matrix are given.
// 
// Input.
// x: Seq[RngIntElt]
// u: Seq[RngIntElt]
// q: RngIntElt
//
// Paramters.
// verbose: print counter and intermediate results
//
// Output.
// [*SpecData*] (list of specifying data)
//
ProduceTorsionGradings := function(x, u, q : verbose := false)
	r := #x;
	s := #u;
	SD := [* *];

	// upper bound on order of finite cyclic factors	
	tupper := GCD([TorsionBound(x, u[j]) : j in [1..s]]);
		
	// possible values for the t_i
	D := Exclude(Divisors(tupper), 1);

	if verbose then
		printf "D = %o\nCompute candidate rows...\n", D;
	end if;

	// possible rows in the torsion part of Q
	Qtors := [ CandidateRows(x, u, t) : t in D];

	for Ttup in CartesianPower({d : d in Divisors(tupper) | d gt 1}, q) do

		// discard tuple T if it is not in increasing order ...
		T := TupleSeq(Ttup);
		if T ne Sort(T)  then
			continue;
		end if;

		// .. or the factors of the corresponding finite group can be 
		// reduced by the Chinese Remainder Theorem
		if q gt 1 and \
		Min([GCD(T[Setseq(I)]) : I in Subsets({1..q}, 2)])\
		eq 1 then
			continue;
		end if;


		if verbose then
			printf "T = %o\n", T;
		end if;


		// bring the first torsion row of Q in a more specific constellation
		// by applying a suitably sorting the columns of Q
		Qtors1 := [Q[2] : Q in FilterGradings([ [x] cat [Qtor]\
		          : Qtor in Qtors[Index(D, T[1])] ], [T[1]])];

		// tuples in G are all the candidates for the whole torsion part of Q
		G := CartesianProduct([Qtors1] cat [ Qtors[Index(D, T[k])] : k in [2..q]]);

		goodQs := [];

		dbgCounter := 1;

		// run over all candidates for degree matrices and save those
		// that stem from a terminal Fano threefold in goodQs
		for Qtor in G do
			if verbose then
				printf "%o/%o\n", dbgCounter, #G;
				dbgCounter +:= 1;
			end if;

			Q := [x] cat [Qtor[k] : k in [1..q]];

			// essential tests if specifying data lead
			// to terminal Fano threefold
			if &and[ExtendsToBasepointfreeClass(Q, u0, T) : u0 in Seqset(u)]\
				 and IsAlmostFree(Q : T := T)\
				 and IsACCterminal(Q : T := T)\
			then
				Append(~goodQs, Q);

				if verbose then
					printf "Found degree matrix: %o\n\n",  Q;
				end if;
			end if;
		end for;

		if verbose then
			print "filter gradings...";
		end if;

		goodQs := FilterGradings(goodQs, T);
		for Q in goodQs do
			newSD := [* Q,\
			         [mu where _, mu is ExtendsToBasepointfreeClass(Q, u[j],T)\
					 : j in [1..s]],\
			         T *];
			Append(~SD, newSD);
		end for;
	end for;

	return SD;
end function;
//============================================================================= 



//============================================================================= 
// FanoBpfDegs
//-----------------------------------------------------------------------------
// Compute all u = (u_1, ..., u_s) such that 
// (1) u_1 + ... u_s > x_0 + ... + x_n
// (2) for any i, j we have u_j = l_ij * x_i for some l_ij >= 2
// 
// Input.
// x: Seq[RngIntElt]
// s: RngIntElt
// 
// Output. 
// Seq[Seq[RngIntElt]]
//
FanoBpfDegs := function(x, s)
	ratio := &+x/LCM(x);
	kupper := (IsIntegral(ratio)) select ratio - 1 else Floor(ratio);

	// each u_i is of the form k_i * lcm(x)
	// and k_1 + ... + k_s = k where k * lcm(x) < x_0 + ... + x_n
	if Max(x) eq LCM(x) then
		return &cat[ [[L[i]*LCM(x) : i in [1..s]]\
		       : L in Partitions(k, s) | 1 notin L]\
		       : k in [1..kupper] ];
	else
		return &cat[ [ [L[i]*LCM(x) : i in [1..s]]\
		       : L in Partitions(k, s)]\
		       : k in [1..kupper]];
	end if;
end function;
//============================================================================= 



//============================================================================= 
// FirstTuples
// ----------------------------------------------------------------------------
// Compute all candidates for Z-parts x of specifying data of
// a s-codimensional terminal Fano toric CI in a FWPS within a given bound M
// 
// Input.
// dim: RngIntElt
// s: RngIntElt
// M: RngIntElt
//
// Output.
// [* [* [RngIntElt], [RngIntElt] *] *]
// (a list of lists of all pairs (x, u) satisfying conditions stated
// in the proof of Theorem 1.1.3)
//
FirstTuples := function(dim, s, M : verbose := false)
	SpecifyingData := [* *];
	r := dim+s+1;

	// the while loops allows us go through all ordered r-tuples x without
	// keeping all of them in memory
	x := [1 : j in [1..r]];
	i := 1;	

	while x[r] lt M+1 do
		if i eq 1 then

	//------------------------------------------------------------------------
	// what to do for each tuple
	//------------------------------------------------------------------------
			Degs := FanoBpfDegs(x, s);

			if #Degs ge 1 and IsNwiseCoprime(x, r-2)\
			   and TerminalConstraint(x, 3)\
			then
				for u in Degs do
					newData := [* x, [u[j] : j in [1..s]], [] *];
					Append(~SpecifyingData, newData);

					if verbose then
						// PrintSpecifyingData(newData);
					end if;
				end for;
			end if;
	//------------------------------------------------------------------------

		end if;

		if i eq r then
			x[r] +:= 1;
			i := 1;	
		elif x[i] lt x[i+1] then
			x[i] +:= 1;
			i := 1;	
		else
			x[i] := 1;
			i +:= 1;
		end if;
	end while;

	return SpecifyingData;
end function;
//============================================================================= 



//============================================================================= 
// Classification of terminal Fano threefolds
//-----------------------------------------------------------------------------
// Perform the computions for the proof of Theorem 1.1.3.
// 
// Output.
// A list of lists [* Q, mu, T *] where
// - Q is the (generator) degree matrix
// - mu = [mu_1, ..., mu_s] is the relation degree matrix
// - T is a sequence 
// (a list of specifying data for all threedimensional terminal Fano
// toric complete intersections in a fake weighted projective space)
ClassificationTerminalFano3foldsFWPS := function()
	print "Start classification ...";
	print "x;u;nVerified;time";

	// bounds for s = 1, 2, 3
	M := [41, 21, 1];
	Data := [**];

	for s := 1 to 3 do
		// compute the candidates for the ZZ-parts x and u = (u_1, ..., u_s)
		ZZdata := FirstTuples(3, s, M[s]); 

		// for each par (x,u) go through the finite number of Q, mu
		// having x resp. u as ZZ-part
		for D in ZZdata do
			x := D[1];
			u := D[2];	

			printf "%o;%o;", x, u;

			// affirmatively tested candidates for this configuration of x, u
			newData := [* *];

			T0 := Time();

			// consider torsion-free constellation
			if IsACCterminal([x]) then
 				Append(~newData, [* [x], [[u[j]] : j in [1..s]], [] *]);
			end if;

			// look for candidates with torsion
			for q := 1 to s+1 do
				newData cat:= ProduceTorsionGradings(x, u, q);
			end for;

			printf "%o;%o\n", #newData, Time(T0);
			Data cat:= newData;
		end for;
	end for;

	return Data;
end function;
//============================================================================= 
