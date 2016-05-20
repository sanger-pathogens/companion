Name "aragorn: aragorn_to_gff3 empty data"
Keywords "aragorn aragorn_to_gff3"
Test do
  run_test("echo "" | #{$bin}gt #{$bindir}/aragorn_to_gff3.lua")
  run("diff #{last_stdout} #{$testdata}/aragorn.empty.gff3")
end

Name "aragorn: aragorn_to_gff3 regular file"
Keywords "aragorn aragorn_to_gff3"
Test do
  run("#{$bin}gt gff3 -sort -tidy -retainids #{$testdata}/aragorn.gff3 > ref")
  run_test("#{$bin}gt #{$bindir}/aragorn_to_gff3.lua < #{$testdata}/aragorn.out")
  run("#{$bin}gt gff3 -sort -tidy -retainids #{last_stdout} | tee gff3")
  run("diff #{last_stdout} ref")
  run_test("#{$bin}gt gff3validator gff3")
end

Name "aragorn: aragorn_to_gff3 special case ('???' amino acid, '.XX' anticodon)"
Keywords "aragorn aragorn_to_gff3"
Test do
  run_test("#{$bin}gt #{$bindir}/aragorn_to_gff3.lua < #{$testdata}/aragorn.unreliable_codons.out")
end

Name "aragorn: aragorn_to_gff3 malformed file (incomplete lines)"
Keywords "aragorn aragorn_to_gff3"
Test do
  run_test("#{$bin}gt #{$bindir}/aragorn_to_gff3.lua < #{$testdata}/aragorn.malformed.out", :retval => 1)
end

Name "aragorn: aragorn_to_gff3 malformed file (no seqs)"
Keywords "aragorn aragorn_to_gff3"
Test do
  run_test("#{$bin}gt #{$bindir}/aragorn_to_gff3.lua < #{$testdata}/aragorn.malformed2.out", :retval => 1)
end