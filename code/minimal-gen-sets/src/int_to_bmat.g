IntToBMat := function(n, b)
  b := BlistNumber(b + 1, 64);
  return List([0 .. n - 1], i -> b{[8 * i + 1 .. 8 * i + n]});
end;

x := StringFile("build/output/bmat_filtered_6.txt");
x := EvalString(Concatenation("[", ReplacedString(x, "\n", ","), "]"));
x := List(x, y -> BooleanMat(IntToBMat(6, y)));
WriteGenerators("~/Desktop/bmat-gens-6.gz", x);
