//=============================================================================
/* FanOfFakeWPS
 * ----------------------------------------------------------------------------
 * This is a wrapper for FanOfFakeProjectiveSpace and FanOfWPS
 * from the Toric Varieties package.
*/
intrinsic FanOfFakeWPS(Q::[SeqEnum[RngIntElt]], T::[RngIntElt]) -> TorFan
{ The fan of the fake weighted projective space determined by Q. } 
	require #Q - #T ge 1: "T has too many cyclic factors";
	require IsAlmostFree(Q : T := T):\
		"This data does not define a fake weighted projective space";

	x := Q[1];
	r := #x;

	if #T gt 0 then
		// express torsion rows of Q by rational numbers
		Qrat := [ [q/T[i] : q in Q[i+1]] : i in [1..#T] ];	

		Sigma := FanOfFakeProjectiveSpace(x, Qrat);
	else
		Sigma := FanOfWPS(x);
	end if;

	return Sigma;
end intrinsic;
//============================================================================= 
