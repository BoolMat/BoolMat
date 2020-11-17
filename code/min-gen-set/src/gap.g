#############################################################################
##
##  gap.g
##  Copyright (C) 2020                                   James D. Mitchell
##                                                         Finlay L. Smith
##
##  Licensing information can be found in the README file of this package.
##
#############################################################################

# This file contains various GAP functions required to find the minimum size
# generating sets for the monoid of all n x n Boolean matrices.

# Convert the integer to a blist (the integer is obtained from HPCombi::BMat8
# using the to_int mem fn.
IntToBMat := function(n, b)

  b := BlistNumber(b + 1, 64);
  return List([0 .. n - 1], i -> b{[8 * i + 1 .. 8 * i + n]});
end;

# Convert a file of ints into a file of GAP boolean matrices.
FileOfIntsToListOfBMats := function(n, process_id, sublist, p)
  local x, cycle, id;

  x := StringFile(Concatenation("build/output/bmat_full_",
                                String(n),
                                "_prefiltered.txt"));
  x := EvalString(Concatenation("[", ReplacedString(x, "\n", ","), "]"));
  x := x{OnTuples(sublist, p ^ -1)};
  x := List(x, y -> BooleanMat(IntToBMat(n, y)));
  if process_id = 1 then
    Add(x, AsBooleanMat((1, 2), n));
    cycle := [2 .. n];
    Add(cycle, 1);
    Add(x, AsBooleanMat(PermList(cycle)));
    id := List([1 .. n - 1], x -> BlistList([1 .. n], [x]));
    Add(x, BooleanMat(Concatenation(id, [BlistList([1 .. n], [])])));
    Add(x, BooleanMat(Concatenation(id, [BlistList([1 .. n], [n-1, n])])));
  
  fi;

  if not IsEmpty(x) then
    WriteGenerators(Concatenation("build/output/bmat-gens-",
                                  String(n),
                                  "-",
                                  String(process_id),
                                  ".gz"), x);
  fi;
end;

FileOfBoolMatsToFileOfInts := function(n, fin, fout)
  local x, f, padded, blist, mat, i;
  x := ReadGenerators(fin);  
  f := IO_CompressedFile(fout, "w");
  for mat in x do
    padded := [];
    for i in [1 .. 8] do
      blist := BlistList([1 .. 8], []);
      if i <= n then
        blist{[1 .. n]} := mat[i];
      fi;
      Add(padded, blist);
    od;
    IO_WriteLine(f, String(NumberBooleanMat(BooleanMat(padded)) - 1));
  od;
  IO_Close(f);
end;

# TODO(FLS) what is this?
AugmentedDigraph := function(dim, blists)
  local n, out, i, j;
  n   := Length(blists);
  out := List([1 .. n + dim + 1], x -> []);
  for i in [1 .. n] do
    for j in [1 .. n] do
      if IsSubsetBlist(blists[i], blists[j]) then
        Add(out[i], j);
      fi;
    od;
    # Add(out[i], n + 1 + SizeBlist(blists[i]));
  od;
  # TODO use DigraphTransitiveReduction??
  return Digraph(out);
end;

# Converts digraphs!
# From: one list of numbers per line, each number representing a binary row
# To:   augmented digraphs
WriteAugmentedDigraphs := function(in_filename, out_filename, dim)
  local in_file, line, out_file;

  in_file  := IO_File(in_filename);
  line     := IO_ReadLine(in_file);
  out_file := DigraphFile(out_filename, "w");
  while Length(line) > 0 do
    line := EvalString(line);
    Apply(line, x -> BlistNumber(x + 1, 8));
    WriteDigraphs(out_file, AugmentedDigraph(dim, line));
    line := IO_ReadLine(in_file);
  od;
  IO_Close(in_file);
  IO_Close(out_file);
end;

FilterByHomomorphisms := function(dim, in_filename, out_filename, md, rem)
  local D, quo, range, p, n, valid, count, partial_map, result, i, k;
  D     := ReadDigraphs(in_filename);
  quo   := QuoInt(Size(D), md);
  range := Intersection([rem, rem + md .. quo * md + rem], [1 .. Size(D)]);
  p     := Sortex(D);
  n     := Length(D);
  valid := BlistList([1 .. n], range);
  count := 0;
  for i in range do
    if dim >= 7 and count mod 100 = 0 then
      Print("Process ", rem, " on ", count, " out of ", Size(range), "\n");
    fi;
    partial_map := [1 .. DigraphNrVertices(D[i]) - (dim + 1)];
    Append(partial_map,  [1 .. dim + 1]);
    for k in  [1 .. DigraphNrVertices(D[i]) - (dim + 1)] do
      Unbind(partial_map[k]);
    od;
    valid[n] := false;
    for k in [n - 1, n - 2 .. 1] do
      if i <> k and Size(HomomorphismDigraphsFinder(D[i],
                                            D[k],
                                            fail,
                                            [],
                                            1,
                                            fail,
                                            2,
                                            DigraphVertices(D[k]),
                                            partial_map
                                            + DigraphNrVertices(D[k]) - dim - 1,
                                            fail,
                                            fail)) > 0 then
        valid[i] := false;
        break;
      fi;
    od;
    count := count + 1;
  od;
  result := D{ListBlist([1 .. Length(valid)], valid)};
  # Print(Length(result), "\n");
  WriteDigraphs(out_filename, result, "w");
  FileOfIntsToListOfBMats(dim, rem, ListBlist([1 .. Length(valid)], valid), p);
end;
