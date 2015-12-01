function count_stopcodons(s)
 local pat = "[+*#]"
 local n = 0
 local percent = function(s) n = n + 1 end
 s:gsub(pat, percent)
 return n
end

-- XXX make this more generic
valid_splice_sites_fwd = { five_prime = { GT = true, gt = true },
                           three_prime = { AG = true, ag = true }}

valid_splice_sites_rev = { three_prime = { AC = true, ac = true },
                           five_prime = { CT = true, ct = true }}

function gene_has_valid_splice_site(gene, regionmapping)
  for c in gene:children() do
    if c:get_type() == 'mRNA' then
      local status, rval = pcall(GenomeTools_genome_node.extract_sequence, c, "mRNA", false, regionmapping)
      if not status then
        -- transcripts for which no sequence can be extracted
        io.stderr:write("warning: " .. rval .. "\n")
        return false
      else
        local i = 0
        local n = 0
        local mrna_range = c:get_range()
        if gene:get_strand() == '-' then
          rval = revcomp(rval)
        end
        -- count CDS features in this transcript
        for cds in c:children() do
          if cds:get_type() == 'CDS' then
            n = n + 1
          end
        end
        -- test splice sites
        for cds in c:children() do
          if n > 1 and cds:get_type() == 'CDS' then
            local five_prime_splice_site = '  '
            local left_ok = true
            local three_prime_splice_site = '  '
            local right_ok = true
            if i > 0 then
              three_prime_splice_site = rval:sub(cds:get_range():get_start() - mrna_range:get_start() - 1,
                             cds:get_range():get_start() - mrna_range:get_start())
              if three_prime_splice_site:len() > 0 then
                if gene:get_strand() == '-' then
                  if not valid_splice_sites_rev.three_prime[three_prime_splice_site] then
                    right_ok = false
                  end
                else
                  if not valid_splice_sites_fwd.three_prime[three_prime_splice_site] then
                    right_ok = false
                  end
                end
              end
            end
            if i < n then
              five_prime_splice_site = rval:sub(cds:get_range():get_end() - mrna_range:get_start() + 2,
                             cds:get_range():get_end() - mrna_range:get_start() + 3)
              if five_prime_splice_site:len() > 0 then
                if gene:get_strand() == '-' then
                  if not valid_splice_sites_rev.five_prime[five_prime_splice_site] then
                    left_ok = false
                  end
                else
                  if not valid_splice_sites_fwd.five_prime[five_prime_splice_site] then
                    left_ok = false
                  end
                end
              end
            end
            if not (left_ok and right_ok) then
              --print(gene:get_strand() .. ":    " .. tostring(left_ok) .. "(" .. five_prime_splice_site .. ")" .. " .. " .. tostring(right_ok) .. "(" .. three_prime_splice_site .. ")")
              return false
            end
            i = i + 1
          end
        end
      end
    end
  end
  return true
end

function get_weight(gene, regionmapping)
  local fac = 1
  -- apply reward for being annotated by RATT
  if gene:get_source() == "RATT" then
    fac = 1.5
    local valid_splice_site = gene_has_valid_splice_site(gene, regionmapping)
    if not valid_splice_site then
       -- ... unless a RATT gene has invalid splice site boundaries,
       -- in which case disregard it
      fac = 0
    end
    for c in gene:children() do
      if c:get_type() == 'mRNA' then
        local status, rval = pcall(GenomeTools_genome_node.extract_and_translate_sequence, c, "CDS", true, regionmapping)
        if status then
          prot = rval
          -- ... unless a RATT gene is broken or exceeding boundaries,
          -- in which case disregard it
          if count_stopcodons(prot) > 1 then
            fac = 0
          end
        else
          fac = 0
        end
      end
    end
    -- disregard pseudogenes in reference
    if gene:get_attribute("is_pseudo_in_ref") == 'true' then
      fac = 0
    end
  end
  return gene:get_range():length() * fac
end