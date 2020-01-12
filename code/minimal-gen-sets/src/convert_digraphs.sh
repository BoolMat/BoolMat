#!/bin/bash
~/gap/bin/gap.sh -A -T -q -m  <<< "LoadPackage(\"digraphs\");
Read(\"digraph_maker.gi\"); write_augmented_digraphs(\"row_space_numbers_8.txt\", \"augmented_digraphs_8.d6.gz\", 60, $rem); quit;"
exit 0
