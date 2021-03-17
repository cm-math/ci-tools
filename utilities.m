freeze; 

//============================================================================= 
/* SeqSupp
 * ----------------------------------------------------------------------------
 * Input.
 * S: SeqEnum
 * 
 * Ouput.
 * Seq[RngIntElt], entries are inidices i where S[i] ne 0
*/
SeqSupp := function(S)
	return [i : i in [1..#S] | IsDefined(S, i) and S[i] ne 0];
end function;
//============================================================================= 



//============================================================================= 
/* TupleSeq
 * ----------------------------------------------------------------------------
 * Given a tuple T with all elements of same type,
 * construct a sequence where the elements
 * are the entries of T taken in the same order.
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
// TransposeInSeq
//-----------------------------------------------------------------------------
// Swaps i-th and j-th element of a sequence.
// 
// Input.
// S: SeqEnum
// i: RngIntElt, index 
// j: RngIntElt, index
// 
// Ouput.
// Sequence S' with S'[i] = S[j], S'[j] = S[i] and S'[k] = S[k] otherwise.
//
TransposeInSeq := function(S, i, j);
	x := S[i];
	S[i] := S[j];
	S[j] := x;

	return S;
end function;
//============================================================================= 
