--[[  TEMPLATE CHECKS
      ===============  ]]

check_parent = function (n, parent_type)
  it("appears as part of a " .. parent_type, function()
    expect(n:appears_as_child_of_type(parent_type)).should_be(true)
  end)
end
is_a_lone_feature = function (n)
  it("appears as a root node", function()
    expect(n:appears_as_root_node()).should_be(true)
  end)

  it("should not have children", function()
    expect(count(n:direct_children())).should_be(0)
  end)
end
does_not_cross_a_contig_boundary = function(n)
  it("does not span more than one contig", function()
    local ovl = feature_index:get_features_for_range(n:get_seqid(),
                                                     n:get_range())
    local contigs = {}
    for _,v in ipairs(ovl) do
      if v:get_type() == "contig" or v:get_type() == "region" then
        table.insert(contigs, v)
      end
    end
    expect(#contigs).should_be_smaller_than(2)
  end)
end

--[[  GENE
      ====  ]]

describe.feature("gene", function(gene)

  does_not_cross_a_contig_boundary(gene)

  it("appears as a root node", function()
    expect(gene:appears_as_root_node()).should_be(true)
  end)

  it("contains a transcript", function()
    expect(gene:has_child_of_supertype("transcript")).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in gene:children() do
      expect(gene:get_range():contains(child:get_range())).should_be(true)
    end
  end)

  it("has consistent strands across all children", function()
    for child in gene:children() do
      expect(gene:get_strand()).should_be(child:get_strand())
    end
  end)

  it("does not overlap another protein coding gene", function()
    local nof_genes_at_site = 0
    local overlapping = feature_index:get_features_for_range(gene:get_seqid(),
                                                             gene:get_range())
    if gene:has_child_of_type("mRNA") then
      for _,f in ipairs(overlapping) do
        if f:get_type() == "gene" and f:has_child_of_type("mRNA") then
          nof_genes_at_site = nof_genes_at_site + 1
        end
      end
      expect(nof_genes_at_site).should_equal(1)
    end
  end)
end)

--[[  PSEUDOGENE
      ==========  ]]

describe.feature("pseudogene", function(pseudogene)

  does_not_cross_a_contig_boundary(pseudogene)

  it("contains a pseudogenic_transcript", function()
    expect(pseudogene:has_child_of_type("pseudogenic_transcript")).should_be(true)
  end)

  it("appears as a root node", function()
    expect(pseudogene:appears_as_root_node()).should_be(true)
  end)

  it("contains all child features within its coordinates", function()
    for child in pseudogene:children() do
      expect(pseudogene:get_range():overlap(child:get_range())).should_be(true)
    end
  end)

  it("has consistent strands across all children", function()
    for child in pseudogene:children() do
      expect(pseudogene:get_strand()).should_be(child:get_strand())
    end
  end)
end)

--[[  PSEUDOGENIC_TRANSCRIPT
      ======================  ]]

describe.feature("pseudogenic_transcript", function(ptranscript)
  check_parent(ptranscript, "pseudogene")

  it("contains at least one pseudogenic_exon", function()
    expect(ptranscript:has_child_of_type("pseudogenic_exon")).should_be(true)
  end)

end)

--[[  PSEUDOGENIC_EXON
      ================  ]]

describe.feature("pseudogenic_exon", function(pexon)
  check_parent(pexon, "pseudogenic_transcript")

  it("should not have children", function()
    expect(count(pexon:direct_children())).should_be(0)
  end)
end)

--[[  MRNA
      ====  ]]

describe.feature("mRNA", function(mrna)
  local dnaseq = mrna:extract_sequence("CDS", true, region_mapping):lower()
  local protseq = mrna:extract_and_translate_sequence("CDS", true,
                                                      region_mapping)

  check_parent(mrna, "gene")

  it("has at least one CDS child", function()
    expect(mrna:has_child_of_type("CDS")).should_be(true)
  end)

  it("has only allowed child types", function()
    for f in mrna:direct_children() do
      expect({"CDS", "exon", "five_prime_UTR",
              "three_prime_UTR"}).should_contain(f:get_type())
    end
  end)

  it("has a coding sequence >= 3bp", function()
    expect(dnaseq:len()).should_be_larger_than(2)
  end)

  it("has a protein product >= 10aa", function()
    expect(protseq:len()).should_be_larger_than(10)
  end)

  it("has non-partial CDS ending on a stop codon", function()
    local in_partial_gene = false
    for _,g in ipairs(mrna:get_path()) do
      if g:get_type() == "gene" then
        if ((g:get_strand() == "+" and g:get_attribute("End_range")) or
            (g:get_strand() == "-" and g:get_attribute("Start_range")) or
             g:get_attribute("threeEndPartial")) then
          in_partial_gene = true
          break
        end
      end
    end
    if not in_partial_gene then
      expect(protseq:sub(-1)).should_match("[*+#]")
    end
  end)

  it("has non-partial CDS beginning with a start codon", function()
    local in_partial_gene = false
    for _,g in ipairs(mrna:get_path()) do
      if g:get_type() == "gene" then
        if ((g:get_strand() == "+" and g:get_attribute("Start_range")) or
            (g:get_strand() == "-" and g:get_attribute("End_range")) or
             g:get_attribute("fiveEndPartial")) then
          in_partial_gene = true
          break
        end
      end
    end
    if not in_partial_gene then
      expect(protseq:sub(1,1)).should_match("[Mm]")
    end
  end)

  it("agrees exactly with CDS/UTR coordinates of its children", function()
    local rng = nil
    -- collect and join CDS ranges
    for c in mrna:children() do
      if c:get_type() == "CDS" or string.match(c:get_type(), "UTR") then
        if not rng then
          rng = c:get_range()
        else
          rng = rng:join(c:get_range())
        end
      end
    end
    -- should overlap with at least one feature
    expect(rng).should_be_truthy()
    -- check if coordinates match
    if rng then
      expect(rng:get_start() == mrna:get_range():get_start() and
             rng:get_end() == mrna:get_range():get_end()).should_be_truthy()
    end
  end)
end)

--[[  FIVE_PRIME_UTR
      ==============  ]]

describe.feature("five_prime_UTR", function(futr)
  check_parent(futr, "gene")

  it("should not have children", function()
    expect(#(collect(futr:direct_children()))).should_be(0)
  end)
end)

--[[  THREE_PRIME_UTR
      ===============  ]]

describe.feature("three_prime_UTR", function(tutr)
  check_parent(tutr, "gene")

  it("should not have children", function()
    expect(#(collect(tutr:direct_children()))).should_be(0)
  end)
end)

--[[  CDS
      ===  ]]

describe.feature("CDS", function(cds)
  it("appears as child of an mRNA", function()
    expect(cds:appears_as_child_of_type("mRNA")).should_be(true)
  end)

  it("should not have children", function()
    expect(#(collect(cds:direct_children()))).should_be(0)
  end)
end)

--[[  POLYPEPTIDE
      ===========  ]]

derives_from = {}
describe.feature("polypeptide", function(pp)
  does_not_cross_a_contig_boundary(pp)

  it("must be derived_from a unique mRNA", function()
    local dfrom = pp:get_attribute("Derives_from")
    expect(dfrom).should_not_be(nil)
    expect(derives_from).should_not_have_key(dfrom)
    derives_from[dfrom] = true
  end)

  it("appears as a root node", function()
    expect(pp:appears_as_root_node()).should_be(true)
  end)

  it("has a required product attribute", function()
    expect(pp:get_attribute("product")).should_not_be(nil)
  end)

  local overlapping = feature_index:get_features_for_range(pp:get_seqid(),
                                                           pp:get_range())

  it("overlaps the transcript it derives_from", function()
    local num_transcripts = 0
    expect(#overlapping).should_be_larger_than(0)
    if #overlapping > 0 then
      for _,ovl_feat in ipairs(overlapping) do
        for c in ovl_feat:children() do
          if c:get_attribute("ID") == pp:get_attribute("Derives_from") then
            num_transcripts = num_transcripts + 1
          end
        end
      end
      expect(num_transcripts).should_be_larger_than(0)
    end
  end)

  it("agrees exactly with CDS of >=1 overlapping coding transcripts", function()
    local nof_possible = 0
    local nof_correct = 0
    -- check every feature in the range
    for _,ovl_feat in ipairs(overlapping) do
      for n in ovl_feat:children() do
        if n:get_type() == "mRNA"
            or n:get_type() =="pseudogenic_transcript" then
          -- locate transcript (no pseudogene etc)
          local rng = nil
          -- collect and join CDS ranges
          for c in n:children() do
            if c:get_type() == "CDS"
                or c:get_type() == "pseudogenic_exon" then
              if not rng then
                rng = c:get_range()
              else
                rng = rng:join(c:get_range())
              end
            end
          end
          -- check if coordinates match
          if rng then
            nof_possible =  nof_possible + 1
            if rng:get_start() == pp:get_range():get_start() and
               rng:get_end() == pp:get_range():get_end() then
              nof_correct = nof_correct + 1
            end
          end
        end
      end
    end
    if nof_possible > 0 then
      expect(nof_correct).should_be_larger_than(0)
    end
  end)
end)

--[[  NCRNA
      =====  ]]

describe.feature("ncRNA", function(node)
  check_parent(node, "gene")
end)

--[[  TRNA
      ====  ]]

describe.feature("tRNA", function(node)
  check_parent(node, "gene")
end)

--[[  RRNA
      ====  ]]

describe.feature("rRNA", function(node)
  check_parent(node, "gene")
end)
