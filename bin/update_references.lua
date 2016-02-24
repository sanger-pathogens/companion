#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2015 Genome Research Ltd

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
]]

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("table_save")
local json = require ("dkjson")

ncrna_visitor = gt.custom_visitor_new()
function ncrna_visitor:visit_feature(fn)
  for fn2 in fn:get_children() do
    if fn2:get_type() == 'ncRNA' or
         fn2:get_type() == 'snRNA' or
         fn2:get_type() == 'snoRNA' or
         fn2:get_type() == 'rRNA' or
         fn2:get_type() == 'lncRNA' or
         fn2:get_type() == 'snRNA' or
         fn2:get_type() == 'scRNA' or
         fn2:get_type() == 'tRNA' then
      for n in fn2:get_children() do
        if n:get_type() == "CDS" or n:get_type():match("UTR") then
          if n:get_attribute("Dbxref") and not fn2:get_attribute("Dbxref") then
            -- copy over Dbxref from CDS to ncRNA
            fn2:add_attribute("Dbxref", n:get_attribute("Dbxref"))
          end
          if n:get_attribute("product") and not fn2:get_attribute("product") then
            -- copy over product from CDS to ncRNA
            fn2:add_attribute("product", n:get_attribute("product"))
          end
          -- remove CDS from ncRNA
          fn:remove_leaf(n)
        end
      end
      if fn2:get_attribute("score") then
        fn2:remove_attribute("score")
      end
      if fn2:get_attribute("model_range") then
        fn2:remove_attribute("model_range")
      end
      if fn2:get_attribute("model_name") then
        fn2:remove_attribute("model_name")
      end
      if fn2:get_attribute("evalue") then
        fn2:remove_attribute("evalue")
      end
      if fn2:get_attribute("anticodon") then
        fn2:remove_attribute("anticodon")
      end
      if fn2:get_attribute("gc_content") then
        fn2:remove_attribute("gc_content")
      end
      if fn2:get_attribute("gc") then
        fn2:remove_attribute("gc")
      end
      if fn2:get_attribute("aa") then
        fn2:remove_attribute("aa")
      end
    end
  end
  return 0
end

quotes_visitor = gt.custom_visitor_new()
function quotes_visitor:visit_feature(fn)
  for fn2 in fn:get_children() do
    for k,v in fn2:attribute_pairs() do
      if v:match('"') then
        v = v:gsub('"','')
        fn2:set_attribute(k, v)
      end
    end
  end
  return 0
end

attrs_to_remove = {"translation",
                   "controlled_curation",
                   "gPI_anchor_cleavage_site",
                   "orthologous_to",
                   "ortholog_cluster",
                   "comment",
                   "membrane_structure",
                   "polypeptide_domain",
                   "cytoplasmic_polypeptide_region",
                   "non_cytoplasmic_polypeptide_region",
                   "transmembrane_polypeptide_region",
                   "signal_peptide",
                   "paralogous_to",
                   "timelastmodified",
                   "full_GO",
                   "colour",
                   "feature_id",
                   "isObsolete"}
remove_attrs_visitor = gt.custom_visitor_new()
function remove_attrs_visitor:visit_feature(fn)
  for n in fn:get_children() do
    for _,v in ipairs(attrs_to_remove) do
      if n:get_attribute(v) then
        n:remove_attribute(v)
      end
    end
  end
  return 0
end

rank_fixer_visitor = gt.custom_visitor_new()
function rank_fixer_visitor:visit_feature(fn)
  for n in fn:get_children() do
    local p = n:get_attribute("product")
    if p then
      local products = gff3_extract_structure(p)
      local ranked = {}
      local seen_preferred = false
      -- gather product names
      for _,v in ipairs(products) do
        table.insert(ranked, v)
      end
      local s = n:get_attribute("product_synonym")
      if s then
        local syns = split(s, ",")
        for _,v in ipairs(syns) do
          table.insert(ranked, {term=v, rank=1})
        end
        n:remove_attribute("product_synonym")
      end
      s = n:get_attribute("synonym")
      if s then
        local syns = split(s, ",")
        for _,v in ipairs(syns) do
          table.insert(ranked, {term=v, rank=1})
        end
        n:remove_attribute("synonym")
      end
      -- assign "is_preferred" flag
      if #ranked > 1 then
        for _,v in ipairs(ranked) do
          if not rank and not seen_preferred then
            v.rank = nil
            v.is_preferred = "true"
            seen_preferred = true
          else
            v.rank = 1
            v.is_preferred = nil
          end
        end
      end
      if string.len(gff3_explode(ranked)) > 0 then
        n:set_attribute("product", gff3_explode(ranked))
      end
    end
  end
  return 0
end

synonym_mover_visitor = gt.custom_visitor_new()
function synonym_mover_visitor:visit_feature(fn)
  if fn:get_type() == "gene" or fn:get_type() == "pseudogene" then
    for fn2 in fn:get_children() do
      if fn2:get_type():match("RNA")
          or fn2:get_type():match("transcript")
          or fn2:get_type():match("CDS") then
        local syn = fn2:get_attribute("synonym")
        if syn then
          if not fn:get_attribute("synonym") then
            fn:set_attribute("synonym", syn)
          end
          fn2:remove_attribute("synonym")
        end
      end
    end
  end
  return 0
end

polypeptide_child_trimmer_visitor = gt.custom_visitor_new()
function polypeptide_child_trimmer_visitor:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    for fn2 in fn:direct_children() do
      if fn2:get_type() == "membrane_structure" then
        for fn3 in fn2:direct_children() do
          fn:remove_leaf(fn3)
        end
      end
      fn:remove_leaf(fn2)
    end
  end
  return 0
end

exon_remover_visitor = gt.custom_visitor_new()
function exon_remover_visitor:visit_feature(fn)
  if fn:get_type() == "gene" then
    for fn2 in fn:children() do
      if fn2:get_type() == "exon" then
        fn:remove_leaf(fn2)
      end
    end
  end
  return 0
end

stat_visitor = gt.custom_visitor_new()
stat_visitor.stats = {}
stat_visitor.stats.nof_genes = 0
stat_visitor.stats.nof_coding_genes = 0
stat_visitor.stats.nof_regions = 0
stat_visitor.stats.nof_chromosomes = 0
function stat_visitor:visit_feature(fn)
  local seqid = fn:get_seqid()
  if fn:get_type() == 'gene' then
    local coding = false
    self.stats.nof_genes = self.stats.nof_genes + 1
    for n in fn:get_children() do
      if n:get_type() == 'mRNA' then
        if not coding then
          coding = true
        end
      end
    end
    -- is coding?
    if coding then
      self.stats.nof_coding_genes = self.stats.nof_coding_genes + 1
    end
  end
  return 0
end
function stat_visitor:visit_region(rn)
  local seqid = rn:get_seqid()
  self.stats.nof_regions = self.stats.nof_regions + 1
  -- how many sequences are full chromosomes?
  if self.chromosome_pattern and string.match(seqid, self.chromosome_pattern) then
    self.stats.nof_chromosomes = self.stats.nof_chromosomes + 1
  elseif self.chromosome_mapping and self.chromosome_mapping[seqid] then
    self.stats.nof_chromosomes = self.stats.nof_chromosomes + 1
  end
  return 0
end


-- =========================================

-- load local 'references.json' file
local reffile = io.open("references-in.json", "rb")
if not reffile then
  error("could not open references-in.json in current directory")
end
local refcontent = reffile:read("*all")
reffile:close()
refs = json.decode(refcontent)

-- import all defined references
for name, values in pairs(refs.species) do
  print("Importing ".. name .. "...")
  local has_chr_mapping = false
  if values.chr_mapping then
    has_chr_mapping = true
  end

  -- make clean organism directory
  os.execute("rm -rf " .. name)
  os.execute("mkdir " .. name)

  -- tidy annotation
  if file_exists(values.gff) then
    os.execute("gt gff3 -tidy -retainids -sort -fixregionboundaries "
                  .. values.gff .. " > " .. name .. "/annotation_preclean.gff3")
  end

  -- fix up annotation
  stat_visitor.stats = {}
  stat_visitor.stats.nof_genes = 0
  stat_visitor.stats.nof_coding_genes = 0
  stat_visitor.stats.nof_regions = 0
  stat_visitor.stats.nof_chromosomes = 0
  stat_visitor.chromosome_pattern = values.chromosome_pattern
  stat_visitor.chromosome_mapping = values.chr_mapping
  stat_visitor.name = values.name
  fixup_stream = gt.custom_stream_new_unsorted()
  fixup_stream.instream = gt.gff3_in_stream_new_sorted(name .. "/annotation_preclean.gff3")
  function fixup_stream:next_tree()
    local node = self.instream:next_tree()
    if node then
      node:accept(ncrna_visitor)
      node:accept(quotes_visitor)
      node:accept(remove_attrs_visitor)
      node:accept(rank_fixer_visitor)
      node:accept(synonym_mover_visitor)
      node:accept(polypeptide_child_trimmer_visitor)
      node:accept(exon_remover_visitor)
      node:accept(stat_visitor)
    end
    return node
  end
  out_stream = gt.gff3_out_stream_new_retainids(fixup_stream, name .. "/annotation.gff3")
  local gn = out_stream:next_tree()
  while (gn) do
    gn = out_stream:next_tree()
  end
  if file_exists(name .. "/annotation_preclean.gff3") then
    os.remove(name .. "/annotation_preclean.gff3")
  end
  values.gff = nil
  out_stream = nil
  collectgarbage() -- hack, forces close of streamed file
  io.flush()

  -- prepare genome FASTA
  if file_exists(values.genome) then
    if values.genome:match("%.gz$") then
      os.execute("zcat " .. values.genome .. " > " .. name .. "/genome.fasta.in")
    else
      os.execute("cp " .. values.genome .. " " .. name .. "/genome.fasta.in")
    end
  end

  -- filter out chromosomes
  if file_exists(name .. "/genome.fasta.in") then
    -- realize mapping as hash if not given before
    if not has_chr_mapping then
      values.chr_mapping = {}
    end
    local keys, seqs = get_fasta_nosep(name .. "/genome.fasta.in")
    if values.chromosome_pattern or has_chr_mapping then
      local outfile = io.open(name .. "/chromosomes.fasta", "w+")
      for hdr, seq in pairs(seqs) do
        local trans_id = hdr:split(' ')[1]
        local m = nil
        if values.chr_mapping[trans_id] then
          m = trans_id
        elseif values.chromosome_pattern then
          m = hdr:match(values.chromosome_pattern)
        end
        if m then
          outfile:write(">" .. trans_id .. "\n")
          print_max_width(seq, outfile, 60)
          if not has_chr_mapping then
            values.chr_mapping[trans_id] = m
          end
        end
      end
    end
    local genomeoutfile = io.open(name .. "/genome.fasta", "w+")
    for hdr, seq in pairs(seqs) do
      local trans_id = hdr:split(' ')[1]
      genomeoutfile:write(">" .. trans_id .. "\n")
      print_max_width(seq, genomeoutfile, 60)
    end
  end
  values.genome = nil
  values.chromosomes = nil
  os.remove(name .. "/genome.fasta.in")
  collectgarbage() -- hack, forces close of streamed file
  io.flush()

  -- prepare GAF
  if values.gaf and file_exists(values.gaf) then
    os.execute("cp " .. values.gaf .. " " .. name .. "/go.gaf")
  else
    io.stderr:write("warning: GAF for " .. name .. " does not exist in " .. tostring(values.gaf) .. "\n")
  end
  values.gaf = nil

  -- prepare models
  -- SNAP
  if values.snap_model and file_exists(values.snap_model) then
    os.execute("cp " .. values.snap_model .. " " .. name .. "/snap.hmm")
    values.snap_model = true
  else
    values.snap_model = nil
  end

  -- AUGUSTUS
  AUG_MODEL_FILES = {'_parameters.cfg', '_parameters.cfg.orig1'}
  if values.augustus_model and file_exists(values.augustus_model) then
    os.execute("mkdir -p " .. lfs.currentdir() .. "/" .. name .. "/species")
    os.execute("cp -pr " .. values.augustus_model
               .. " " .. lfs.currentdir() .. "/"
               ..name .. "/species/augustus_species")
    os.execute("cp -pr " .. gt.script_dir .. "/../data/augustus/model "
               .. lfs.currentdir() .. "/" ..name .. "/model")
    os.execute("cp -pr " .. gt.script_dir .. "/../data/augustus/profile "
               .. lfs.currentdir() .. "/" ..name .. "/profile")
    values.augustus_model = true
    local species_path = lfs.currentdir() .. "/" .. name .. "/species/augustus_species"
    for file in lfs.dir(species_path) do
      --file is the current file or directory name
      for _,pat in ipairs(AUG_MODEL_FILES) do
        if file:match(pat .. "$") then
          os.rename(species_path .."/" .. file,
                    species_path .. "/augustus_species" .. pat)
        end
      end
    end
  else
    values.augustus_model = nil
  end

  -- extract proteins
  -- XXX TODO: check for applicability of mapping
  os.execute("gt extractfeat -type CDS -join -retainids -translate -seqfile "
    .. name .. "/genome.fasta -matchdescstart -o ".. name .. "/proteins_preclean.fasta "
    .. name .. "/annotation.gff3")

  -- open ggfile
  local ggfile = io.open(name .. "/ggline.gg", "w+")
  ggfile:write(name .. ": ")

  -- truncate proteins, make gg line
  if file_exists(name .. "/proteins_preclean.fasta") then
    local keys, seqs = get_fasta_nosep(name .. "/proteins_preclean.fasta")
    local outfile = io.open(name .. "/proteins.fasta", "w+")
    for hdr, seq in pairs(seqs) do
      local trans_id = hdr:split(' ')[1]
      outfile:write(">" .. trans_id .. "\n")
      print_max_width(seq, outfile, 60)
      ggfile:write(trans_id .. " ")
    end
  end
  ggfile:write("\n")
  if file_exists(name .. "/proteins_preclean.fasta") then
    os.remove(name .. "/proteins_preclean.fasta")
  end
  values.pep = nil -- lfs.currentdir() .. "/" .. name .. "/proteins.fasta"

  values.nof_genes = stat_visitor.stats.nof_genes
  values.nof_coding_genes = stat_visitor.stats.nof_coding_genes
  values.nof_regions = stat_visitor.stats.nof_regions
  values.nof_chromosomes = stat_visitor.stats.nof_chromosomes

  -- write out table with metadata (number of genes, etc.)
  metadata_json_out = io.open(name .. "/metadata.json", "w+")
  metadata_json_out:write(json.encode(values,{ indent = true }))
  metadata_json_out:write("\n")
  assert(table.save(stat_visitor.stats, name .. "/metadata.lua" ) == nil )
end

full_json_out = io.open("references.json", "w+")
full_json_out:write(json.encode(refs,{ indent = true }))
full_json_out:write("\n")