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

function usage()
  io.stderr:write(string.format("Usage: %s <ref GFF annotation> " ..
                                "<target GFF annotation> <blast output> " ..
                                "<outpath> <chr_pattern> <bin_chr> <ref_dir> " ..
                                "<ref_name>\n" , arg[0]))
  os.exit(1)
end

if #arg < 6 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
local json = require ("dkjson")

genes_out = io.open(arg[4] .. "/genes.txt", "w+")
gaps_out = io.open(arg[4] .. "/gaps.txt", "w+")
links_out = io.open(arg[4] .. "/links.txt", "w+")
chr_out = io.open(arg[4] .. "/chromosomes.txt", "w+")
karyotype_out = io.open(arg[4] .. "/karyotype.txt", "w+")
bin_out = io.open(arg[4] .. "/bin.txt", "w+")

chr_pattern = arg[5]
bin_chr = arg[6]
ref_dir = arg[7]
ref_name = arg[8]
local reffile = io.open(ref_dir .. '/' .. "references.json", "rb")
if not reffile then
  error("could not open references.json in " .. ref_dir)
end
local refcontent = reffile:read("*all")
reffile:close()
refs = json.decode(refcontent)
if not refs.species[ref_name] then
  error("reference " .. ref_name .. " not found in JSON definition")
end

has_bin = false

-- gather real lengths of reference chromsomes
refsizes = {}
refhdr, refseqs = get_fasta_nosep(ref_dir .. '/' .. ref_name .. '/chromosomes.fasta')
for k, v in pairs(refseqs) do
  refsizes[k] = v:len()
end

-- holds reference chromosomes
ref_chr = {}

visitor = gt.custom_visitor_new()
function visitor:visit_feature(fn)
  if fn:get_type() == 'gene' or fn:get_type() == 'gap' then
    local color = "red"
    if fn:get_strand() == "-" then
      color = "blue"
    end
    if fn:get_type() == "gap" then
      color = "yellow"
      gaps_out:write(fn:get_seqid() .. "  " .. fn:get_range():get_start()
          .. "  " .. fn:get_range():get_end() .. "  color=" .. color .. "\n")
    else
      genes_out:write(fn:get_seqid() .. "  " .. fn:get_range():get_start()
          .. "  " .. fn:get_range():get_end() .. "  color=" .. color .. "\n")
    end
  end
  return 0
end
function visitor:visit_region(rn)
  local endrng = rn:get_range():get_end()
  if self.is_ref and refsizes[rn:get_seqid()] then
    endrng = refsizes[rn:get_seqid()]
  end
  if tonumber(endrng) == nil then
   error("non-numeric end range encountered for chromosome " .. rn:get_seqid() .. ": " .. endrng)
  end
  karyotype_out:write("chr - " .. rn:get_seqid() .. " " .. rn:get_seqid()
                        .. " 1 " .. endrng .. " ")
  if rn:get_seqid() == bin_chr then
    karyotype_out:write("black")
  else
    karyotype_out:write("grey")
  end
  karyotype_out:write("\n")
  if rn:get_seqid() == bin_chr then
    has_bin = true
  end
  if self.print_chr then
    m = rn:get_seqid():match(chr_pattern)
    if m and rn:get_seqid() ~= bin_chr then
      chr_out:write(m .. "\n")
    end
  end
  if self.is_ref then
    local chr = nil
    if refs.species[ref_name].chromosome_pattern then
      chr = rn:get_seqid():match(refs.species[ref_name].chromosome_pattern)
    end
    if refs.species[ref_name].chr_mapping then
      chr = refs.species[ref_name].chr_mapping[rn:get_seqid()]
    end
    if chr then
      table.insert(ref_chr, rn:get_seqid())
    end
  end
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = visitor
visitor.print_chr = true
visitor.is_ref = true
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[2])
visitor_stream.vis = visitor
visitor.print_chr = false
visitor.is_ref = false
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

for l in io.lines(arg[3]) do
  targetid, refid, _, _, _, _, tfrom, tto, rfrom, rto, eval = unpack(split(l, "%s+"))
  if math.abs(tonumber(tto) - tonumber(tfrom)) >= 1500 then
    links_out:write(targetid .. " " .. tfrom .. " " .. tto .. " " .. refid .. " " .. rfrom .. " " .. rto .."\n")
  end
end

if has_bin then
  bin_out:write(bin_chr .. "\t" .. table.concat(ref_chr, ";") .. "\n")
end
