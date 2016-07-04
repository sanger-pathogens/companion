#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2015-2016 Genome Research Ltd

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

function usage()
  io.stderr:write(string.format("Usage: %s <GFF annotation> <seqfile>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

seqkeys, seqs = get_fasta_nosep(arg[2])
rm = gt.region_mapping_new_seqfile_matchdescstart(arg[2])

ftvis = gt.custom_visitor_new()
ftvis.mrnas = {}
ftvis.pseudogenic_transcripts = {}
ftvis.other = {}
function ftvis:visit_feature(fn)
  for n in fn:get_children() do
    if n:get_type() == 'mRNA' then
      local id = n:get_attribute('ID')
      if id then
        ftvis.mrnas[id] = true
      end
      break
    elseif n:get_type() == 'pseudogenic_transcript' then
      local id = n:get_attribute('ID')
      if id then
        ftvis.pseudogenic_transcripts[id] = true
      end
      break
    elseif n:get_type():match('.+RNA$') then
      local id = n:get_attribute('ID')
      if id then
        ftvis.other[id] = true
      end
      break
    end
  end
end

cv = gt.custom_visitor_new()
cv.nof_genes = 0
cv.nof_pseudogenes = 0
cv.nof_trnas = 0
cv.nof_coding_genes = 0
cv.nof_genes_with_introns = 0
cv.nof_genes_with_function = 0
cv.nof_pseudogenes_with_function = 0
cv.nof_genes_with_transferred = 0
cv.nof_singleton_genes = 0
cv.nof_singleton_genes_with_function = 0
cv.nof_regions = 0
cv.nof_chromosomes = 0
cv.overall_length = 0
cv.coding_length = 0
cv.gc_overall = 0
cv.gc_coding = 0
function cv:visit_feature(fn)
  local seqid = fn:get_seqid()
  if fn:get_type() == 'gene' then
    local nof_exons = 0
    local coding = false
    cv.nof_genes = cv.nof_genes + 1
    for n in fn:get_children() do
      if n:get_type() == 'mRNA' then
        if not coding then
          coding = true
        end
        -- calculate coding GC content
        local seq = string.lower(n:extract_sequence('CDS', true, rm))
        for c in seq:gmatch("[gc]") do
          cv.gc_coding = cv.gc_coding + 1
        end
        for c in seq:gmatch("[agct]") do
          cv.coding_length = cv.coding_length + 1
        end
      elseif n:get_type() == 'tRNA' then
        cv.nof_trnas = cv.nof_trnas + 1
      elseif n:get_type() == 'CDS' then
        nof_exons = nof_exons + 1
      end
    end
    -- is coding?
    if coding then
      cv.nof_coding_genes = cv.nof_coding_genes + 1
    end
    -- has more than one exon?
    if nof_exons > 1 then
      cv.nof_genes_with_introns = cv.nof_genes_with_introns + 1
    end
  elseif fn:get_type() == 'pseudogene' then
    cv.nof_pseudogenes = cv.nof_pseudogenes + 1
  elseif fn:get_type() == 'polypeptide' then
    local orths = fn:get_attribute("orthologous_to")
    local product = fn:get_attribute("product")
    -- is species-specific?
    if not orths then
      cv.nof_singleton_genes = cv.nof_singleton_genes + 1
    end
    -- has defined function?
    if product and not string.match(product, "hypothetical") then
      local df = fn:get_attribute("Derives_from")
      if not df then
        warning("polypeptide feature on line " .. fn:get_line_number()
                  .. " does not have a Derives_from attribute, skipping")
      else
        if ftvis.pseudogenic_transcripts[df] then
          cv.nof_pseudogenes_with_function = cv.nof_pseudogenes_with_function + 1
        else
          cv.nof_genes_with_function = cv.nof_genes_with_function + 1
        end
        if not orths then
          cv.nof_singleton_genes_with_function = cv.nof_singleton_genes_with_function + 1
        end
        -- function transferred from reference?
        if string.match(product, 'with.3DGeneDB:') then
          cv.nof_genes_with_transferred = cv.nof_genes_with_transferred + 1
        end
      end
    end
  end
  return 0
end
function cv:visit_region(rn)
  local seqid = rn:get_seqid()
  cv.nof_regions = cv.nof_regions + 1
  -- how many sequences are full chromosomes?
  if string.match(seqid, '^[^.]+_[0-9]+') then
    cv.nof_chromosomes = cv.nof_chromosomes + 1
  end
  return 0
end
function cv:calc_gc_overall()
  return (cv.gc_overall/cv.overall_length)
end
function cv:calc_gc_coding()
  return (cv.gc_coding/cv.coding_length)
end

-- tool code starts here

for h,s in pairs(seqs) do
  local seq = string.lower(s)
  -- calculate overall GC content
  for c in seq:gmatch("[GCgc]") do
    cv.gc_overall = cv.gc_overall + 1
  end
  for c in seq:gmatch("[ATGCagtc]") do
    cv.overall_length = cv.overall_length + 1
  end
end

visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.visitor)
  end
  return node
end

visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.visitor = ftvis
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.visitor = cv
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end


print("nof_regions: " .. cv.nof_regions)
--print("nof_chromosomes: " .. cv.nof_chromosomes)
print("overall_length: " .. cv.overall_length)
print("gc_overall: " .. string.format("%.2f", cv:calc_gc_overall()*100))
print("nof_genes: " .. cv.nof_genes)
print("nof_pseudogenes: " .. cv.nof_pseudogenes)
print("gene_density: " .. string.format("%.2f", cv.nof_coding_genes/(cv.overall_length/1000000)))
print("avg_coding_length: " .. string.format("%d", cv.coding_length/cv.nof_genes))
print("nof_coding_genes: " .. cv.nof_coding_genes)
print("nof_genes_with_mult_cds: " .. cv.nof_genes_with_introns)
print("nof_genes_with_function: " .. cv.nof_genes_with_function)
print("nof_pseudogenes_with_function: " .. cv.nof_pseudogenes_with_function)
--print("nof_genes_with_transferred\: " .. cv.nof_genes_with_transferred)
print("nof_trnas: " .. cv.nof_trnas)
print("gc_coding: " .. string.format("%.2f", cv:calc_gc_coding()*100))
