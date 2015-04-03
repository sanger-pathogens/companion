Name "gtscripts are runnable"
Keywords "misc"
Test do
  run_test("#{$bin}gt #{$bindir}/test.lua")
  grep last_stdout, /travis test runs/
end
