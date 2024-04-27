intrinsic ToricModel(gammalist::SeqEnum[RngIntElt]) -> Any
  {Returns a toric model for the hypergeometric motive with given gamma list.
   Only a little attempt is made to minimize the equation.}

  gammaseq := [gamma : gamma in gammalist];
  H := HypergeometricData(gammalist);
  _, i0 := Min([Abs(gamma) : gamma in gammaseq]);
  gammami0 := [gammalist[i] : i in [1..#gammalist] | i ne i0];
  mker := Kernel(Transpose(Matrix([gammami0])));
  mmat0 := Matrix(Basis(mker));
  mmat := HorizontalJoin(
          HorizontalJoin(ColumnSubmatrix(mmat0,i0-1),Matrix([[0] : i in [1..Nrows(mmat0)]])),
          ColumnSubmatrix(mmat0,i0,Ncols(mmat0)-i0+1));
  mmat := LLL(mmat);
  g, k := Xgcd(gammaseq);
  assert g eq 1;  // I only understand what's happening in the primitive case

  _<u> := FunctionField(Rationals());
  R<[x]> := PolynomialRing(Parent(u), Nrows(mmat));
  mmatT := Transpose(mmat);
  f := 0;
  for i := 1 to Nrows(mmatT) do
    ei := Eltseq(mmatT[i]);
    f +:= u^k[i] * &*[(1/(1/x[j]))^ei[j] : j in [1..#ei]];
  end for;
  t := (-1)^(&+[gamma : gamma in gammaseq | gamma gt 0])*MValue(H)*u;
  return f, t;
end intrinsic;