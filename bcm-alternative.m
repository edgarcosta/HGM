SublistsSummingToZeroHelper := function(S);
  n := #S;
  if n eq 0 then
    return [];
  end if;
  K := Kernel(Transpose(Matrix(Vector(S))));
  vecs := [x[1] : x in ShortVectors(Lattice(K), n) | &and[xi in [0,1] : xi in Eltseq(x[1])]];
  sublists := {Sort([[S[i] : i in [1..n] | v[i] eq 1], [S[i] : i in [1..n] | v[i] eq 0]]) : v in vecs};
  return SetToSequence(sublists);
end function;

SublistsSummingToZero := function(Slist);
  S := [x : x in Slist];
  sublistFrontier := {[S]};
  sublists := {[S]};
  repeat
    // try to split 
    Sd := Random(sublistFrontier);  // wish I could just get the first one
    splits := CartesianProduct([SublistsSummingToZeroHelper(Sd[i]) : i in [1..#Sd] | Sd[i] ne []]);
    if #splits gt 1 then
      sublistFrontier join:= {&cat[c : c in spli] : spli in splits | [] notin spli[1]};
    end if;
    sublists join:= {[Sdi : Sdi in Sd | Sdi ne []]};
    sublistFrontier diff:= {Sd};
  until sublistFrontier eq {};
  return SetToSequence(sublists);
end function;

// Hs := PossibleHypergeometricData(4 : Weight := 1);
// [SublistsSummingToZero(GammaList(HypergeometricData(H))) : H in Hs];   

// munged from Geometry/HypGeomMot/identify.m
function BCMAlternativeSchemes(H,t)
  Gfull := H`Glist; 
  
  schemes := [* *];
  for G in SublistsSummingToZero(Gfull) do
    // if not &and[Gcd(Gi) eq 1 : Gi in G] then continue; end if;  // BCM condition
    Gtot := &cat(G);  // need to ensure everything in the right order
    assert 0 notin Gtot;
    A := AffineSpace(Parent(t), #Gtot);
    eqs := [&*[A.i^(-Gtot[i]) : i in [1..#Gtot] | Gtot[i] lt 0]-
            &*[A.i^(Gtot[i]) : i in [1..#Gtot] | Gtot[i] gt 0]/MValue(H)/t];
    cntr := 0;
    for Gi in G do
      Append(~eqs, &+[A.(cntr+j) : j in [1..#Gi] | Gi[j] lt 0]-1);
      Append(~eqs, &+[A.(cntr+j) : j in [1..#Gi] | Gi[j] gt 0]-1);
      cntr +:= #Gi;
    end for;
    Append(~schemes, <G, Scheme(A, eqs)>);
  end for;
  return schemes;
end function;

/*
_<t> := FunctionField(Rationals());
for H in Hs do
  Hdat := HypergeometricData(H);
  for X in BCMAlternativeSchemes(Hdat,t) do
    if Dimension(X[2]) eq 1 then
      Xc := Curve(X[2]);
      printf "A=%o,B=%o: g = %o, #defeq = %o\n", Hdat`cycA, Hdat`cycB, Genus(Curve(Xc)), #DefiningEquations(Xc);
      if Genus(Xc) eq 3 then
        print Image(CanonicalMap(Xc));
      end if;
    end if;
  end for;
end for;
*/
