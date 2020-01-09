augmented_digraph := function(blists)
  local n, out, i, j;
  n := Length(blists);
  out := List([1 .. n + 9], x -> []);
  for i in [1 .. n] do
    for j in [1 .. n] do
      if IsSubsetBlist(blists[i], blists[j]) then
        Add(out[i], j);
      fi;
    od;
    Add(out[i], n + 1 + SizeBlist(blists[i]));
  od;
  return Digraph(out);
end;

augmented_digraph_2 := function(blists)
  local n, out, i, j;
  n := Length(blists);
  dim := Length(blists[1]);
  out := List([1 .. n + 9 + dim], x -> []);
  for i in [1 .. n] do
    for j in [1 .. n] do
      if IsSubsetBlist(blists[i], blists[j]) then
        Add(out[i], j);
      fi;
    od;
    Add(out[i], n + 1 + SizeBlist(blists[i]));
    for j in ListBlist([1 .. n], blists[i]) do
      Add(out[i], n + 1 + Length(blists[i]);
  od;
  return Digraph(out);
end;
