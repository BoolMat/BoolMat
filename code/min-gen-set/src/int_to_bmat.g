IntToBMat := function(n, b)
  b := BlistNumber(b + 1, 64);
  return List([0 .. n - 1], i -> b{[8 * i + 1 .. 8 * i + n]});
end;

FileOfIntsToListOfBMats := function(n)
  local x;
  x := StringFile(Concatenation("build/output/bmat_filtered_", String(n), ".txt"));
  x := EvalString(Concatenation("[", ReplacedString(x, "\n", ","), "]"));
  x := List(x, y -> BooleanMat(IntToBMat(n, y)));
  Print("writing to ", Concatenation("build/output/bmat-gens-", String(n), ".gz"), "\n");
  return WriteGenerators(Concatenation("build/output/bmat-gens-", String(n), ".gz"), x);
end;
