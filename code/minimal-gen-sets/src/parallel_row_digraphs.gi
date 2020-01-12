filter := function(filename, md, rem)
  local D, quo, range, n, valid, count, partial_map, j, result, i, k;
  D := ReadDigraphs(filename);;
  quo := QuoInt(Size(D), md);
  range := [rem, rem + md .. (quo - 1) * md + rem];
  Sort(D);
  n := Length(D);
  valid := BlistList([1 .. n], range);
  count := 0;
  for i in range do
    if count mod 100 = 0 then
      Print("Process ", rem, " on ", count, " out of ", Size(range), "\n");
    fi;
    partial_map := [1 .. DigraphNrVertices(D[i]) - 9];
    Append(partial_map,  [1 .. 9]);
    for k in  [1 .. DigraphNrVertices(D[i]) - 9] do
      Unbind(partial_map[k]);
    od;
    for k in [1 .. n] do
      j := n - k + 1;
      if i <> j and Size(HomomorphismDigraphsFinder(D[i],
                                            D[j],
                                            fail,
                                            [],
                                            1,
                                            fail,
                                            2,
                                            DigraphVertices(D[j]),
                                            partial_map + DigraphNrVertices(D[j]) - 9,
                                            fail,
                                            fail)) > 0 then
        valid[i] := false;
        break;
      fi;
    od;
    count := count + 1;
  od;
  result := D{ListBlist([1 .. Length(valid)], valid)};
  WriteDigraphs(Concatenation(filename, "_filt_", String(rem), ".d6.gz"), result);
end;
