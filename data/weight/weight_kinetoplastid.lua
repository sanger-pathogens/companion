function get_weight(gene, regionmapping)
  local fac = 1
  local nof_cds = 0
  -- count the number of CDS/exons
  for c in gene:children() do
    if c:get_type() == "CDS" then
      nof_cds = nof_cds + 1
    end
  end
  -- apply reward for being annotated by RATT
  -- unless RATT would produce a multi-exon gene, then only use that
  -- if there's no other gene
  if gene:get_source() == "RATT" then
    if nof_cds == 1 then
      fac = 3
    else
      fac = .1
    end
    -- disregard pseudogenes in reference
    if gene:get_attribute("is_pseudo_in_ref") == 'true' then
      fac = 0
    end
  end
  return gene:get_range():length() * fac
end