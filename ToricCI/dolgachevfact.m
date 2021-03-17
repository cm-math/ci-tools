freeze;

// Returns true iff P is nonsingular in every one-dimensional face
IsNonSingularInOneDimFaces := function(P)
	for F in Faces(P, 1) do
		if IsSingular(DualFaceInDualFan(P,F)) then
			return false, F;
		end if;
	end for;

	return true;
end function;



intrinsic IsDolgachevPolytope(P::TorPol) -> BoolElt
{ Returns true if and only if P is a Dolgachev poyltope. }
	if Dimension(P) lt 3 then
		//printf "polytope is of dimension %o!\n", Dimension(P);
		return false;
	end if;

	// check that P intersects every coordinate hyperplane
	V := [Eltseq(v) : v in Vertices(P)];
	for i := 1 to Dimension(Ambient(P)) do
		meetsHyperplane := false;
		for v in V do
			if v[i] eq 0 then
				meetsHyperplane := true;
				break;
			end if;
		end for;

		if not meetsHyperplane then
			//printf "polytope does not intersect [x_%o = 0]!\n", i;
			return false;
		end if;
	end for;

	if not IsNonSingularInOneDimFaces(P) then
		// printf "polytpe has singular one-dimensional face!\n";
		return false;
	end if;

	return true;
end intrinsic;



intrinsic HasFactorialQuotient(f::RngMPolElt) -> BoolElt
{ Returns true if Dolgachev criterion applies to f and false otherwise. }
	S := Parent(f);
	P := NewtonPolytope(f);

	return IsDolgachevPolytope(P) and IsNondegenerate(f);
end intrinsic;
