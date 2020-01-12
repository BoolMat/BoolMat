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

# Converts digraphs!
# From: one list of numbers per line, each number representing a binary row
# To:   augmented digraphs 
write_augmented_digraphs := function(in, out, dim) 
  fin := IO_File(in);
  x := IO_ReadLine(fin);
  fout := DigraphFile(out, "r");
  while Length(x) > 0 do
    x := EvalString(x);
    Perform(x, y -> BlistNumber(y + 1, dim);
    WriteDigraphs(fout, augmented_digraph(x));
    x := IO_ReadLine(fin);
  od;
  IO_Close(fin);
  IO_Close(fout);
end;
