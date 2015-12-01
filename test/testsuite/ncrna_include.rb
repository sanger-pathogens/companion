Name "ncrna: empty data"
Keywords "ncrna infernal infernal_to_gff3"
Test do
  run_test("echo '' | #{$bin}gt #{$bindir}/infernal_to_gff3.lua")
  run("diff #{last_stdout} #{$testdata}/infernal.empty.gff3")
end

Name "ncrna: infernal_to_gff3 regular file"
Keywords "ncrna infernal infernal_to_gff3"
Test do
  run("#{$bin}gt gff3 -sort -tidy -retainids  #{$testdata}/infernal.gff3 > ref")
  run_test("#{$bin}gt #{$bindir}/infernal_to_gff3.lua < #{$testdata}/infernal.out")
  run("#{$bin}gt gff3 -sort -tidy -retainids #{last_stdout} | tee gff3")
  run("diff #{last_stdout} ref")
  run_test("#{$bin}gt gff3validator gff3")
end

# TODO add tests for malformed files