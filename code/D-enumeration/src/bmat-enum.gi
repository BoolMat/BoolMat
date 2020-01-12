IsRowReducedBooleanMat := function(mat)
  local emptyblist, list, union, row, x;
  emptyblist := BlistList([1 .. DimensionOfMatrixOverSemiring(mat)], []);
  list := AsList(mat);
  list := Filtered(list, x -> x <> emptyblist);
  if not Size(Set(list)) = Size(list) then
    return false;
  fi;
  for row in list do 
    union := emptyblist;
    for x in list do
      if IsSubsetBlist(row, x) and not row = x then
        union := UnionBlist(union, x);
      fi;
    od;
    if union = row then
      return false;
    fi;
  od;
  return true;
end;

IsColReducedBooleanMat := function(mat)
  return IsRowReducedBooleanMat(TransposedMat(mat));
end;

bmat_enum := function(n) 
  local max, reps, mat, ht, nums, row_in_orb, row_orb_by_row, first_row, out,
  empty_blist, stabs, inv_stabs, forbidden, forbidden_at_stage, orbit_reps, bt,
  nr_ones;
 
  max := 2 ^ n;
  reps := [];

  mat := List([1 .. n], x -> BlistList([1 .. n], []));
  ht := HTCreate(0);
  nums := List([1 .. n], x -> 0);
  row_in_orb := List([1 .. max], x -> false);
  row_orb_by_row := List([1 .. n], x -> []);
  first_row := 0;
  out := [];
  empty_blist := BlistNumber(1, n);

  stabs := [];
  inv_stabs := [];
  forbidden := BlistList([1 .. max], []);
  forbidden_at_stage := List([1 .. n], x -> BlistList([1 .. max], []));

  bt := function(k)
    local bl, seen, new_row, x, orb, y, bmat, i, j, l, rown, m, row;
    if k = 0 then
      for i in [1 .. n - 1] do
        first_row := i;
        for j in [2 .. n] do
          bl      := BlistList([1 .. n], [j .. n]); 
          stabs := List([1 .. i], x -> Stabilizer(SymmetricGroup(n), [j .. n],
                                                  OnSets));
          mat[i]  := bl;
          nums[i] := NumberBlist(bl);
          nr_ones := SizeBlist(bl);

          row_orb_by_row[i] := [empty_blist, bl];
          row_in_orb[1] := true;
          row_in_orb[nums[i]] := true;

          forbidden := BlistList([1 .. max], []);

          seen := ListWithIdenticalEntries(max, false);
          orbit_reps := List([1 .. n], x -> []);
          for l in [nums[i] + 1 .. max] do
            if not seen[l] then
              Add(orbit_reps[i + 1], l);
              for x in Orbit(stabs[i], ListBlist([1 .. n], BlistNumber(l, n)), OnSets) do
                y := NumberBlist(BlistList([1 .. n], x));
                seen[y] := true;
              od;
            fi;
          od;

          bt(i + 1);
          row_in_orb[nums[i]] := false;
        od;
        mat[i] := empty_blist;
      od; 
    fi;

    if k in [2 .. n - 1] then
      for rown in orbit_reps[k] do
        bl := BlistNumber(rown, n);
        if not forbidden[rown] and not row_in_orb[rown] and SizeBlist(bl) >= nr_ones then
          mat[k] := bl;
          nums[k] := rown;
 
          #UPDATE ROW ORB
          for m in [first_row .. k - 1] do 
            for row in row_orb_by_row[m] do
              new_row := UnionBlist(row, bl);
              x := NumberBlist(new_row);
              if not row_in_orb[x] then
                Add(row_orb_by_row[k], new_row);
                row_in_orb[x] := true;
              fi;
            od;
          od;

          # UPDATE STABS
          stabs[k] := Stabilizer(stabs[k - 1], ListBlist([1 .. n], bl), OnSets);

          # GET NEXT ORBIT REPS
          seen := ListWithIdenticalEntries(max, false);
          for i in [1 .. rown] do
            if SizeBlist(BlistNumber(i, n)) < nr_ones then
              continue;
            fi;
            if not seen[i] then
              orb := Orbit(stabs[k], ListBlist([1 .. n], BlistNumber(i, n)), OnSets);
              for x in orb do
                y := NumberBlist(BlistList([1 .. n], x));
                seen[y] := true;
              od;
            fi;
          od;

          for i in [rown + 1 .. max] do
            if SizeBlist(BlistNumber(i, n)) < nr_ones then
              continue;
            fi;
            if not seen[i] then
              Add(orbit_reps[k + 1], i);
            fi;
            for x in Orbit(stabs[k], ListBlist([1 .. n], BlistNumber(i, n)), OnSets) do
              seen[NumberBlist(BlistList([1 .. n], x))] := true;
            od;
          od;
          
          #UPDATE FORBIDDEN
          seen := ListWithIdenticalEntries(max, false);
          for i in [1 .. rown - 1] do
            if SizeBlist(BlistNumber(i, n)) < nr_ones then
              continue;
            fi;
            if not seen[i] then
              orb := Orbit(stabs[k - 1], ListBlist([1 .. n], BlistNumber(i, n)), OnSets);
              for x in orb do
                y := NumberBlist(BlistList([1 .. n], x));
                if not forbidden[y] then
                  forbidden_at_stage[k][y] := true;
                fi;
                forbidden[y] := true;
                seen[y] := true;
              od;
            fi;
          od;

          #DIVE
          bt(k + 1);
          
          #RESET ROW ORB
          for row in row_orb_by_row[k] do
            row_in_orb[NumberBlist(row)] := false;
          od;
          row_orb_by_row[k] := [];
          
          # RESET ORBIT REPS
          orbit_reps[k + 1] := [];

          # RESET FORBIDDEN
          SubtractBlist(forbidden, forbidden_at_stage[k]);
          forbidden_at_stage[k] := BlistList([1 .. max], []);
        fi;
      od;
      mat[k] := BlistNumber(1, n);
      nums[k] := 1;
    fi;

    if k = n then
      for rown in orbit_reps[k] do
        bl := BlistNumber(rown, n);
        if not row_in_orb[rown] and SizeBlist(bl) >= nr_ones then
          mat[k] := bl;
          bmat := MatrixNC(BooleanMatType, ShallowCopy(mat));
          if IsColReducedBooleanMat(bmat) then
            x := NumberBooleanMat(CanonicalBooleanMat(bmat));
            if HTValue(ht, x) = fail then
              HTAdd(ht, x, true);
              Add(out, x);
            fi;
          fi;
        fi;
      od;
      mat[k] := BlistNumber(0, n);
      nums[k] := 0;
    fi;
  end;

  bt(0);
  
  #A special case
  Append(out, [2, 1]);
  return out; 
end;
