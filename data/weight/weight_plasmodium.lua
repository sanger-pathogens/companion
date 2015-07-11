function get_weight(gene)
  local fac = 1
  -- apply reward for being annotated by RATT
  if gene:get_source() == "RATT" then
    fac = 5
  end
  return gene:get_range():length() * fac
end