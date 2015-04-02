#!/usr/bin/env gt

--[[
  Copyright (c) 2015 Sascha Steinbiss <ss34@sanger.ac.uk>
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
                                "<outpath> <chr_prefix> <bin_chr>\n" , arg[0]))
  os.exit(1)
end

if #arg < 6 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

genes_out = io.open(arg[4] .. "/genes.txt", "w+")
gaps_out = io.open(arg[4] .. "/gaps.txt", "w+")
links_out = io.open(arg[4] .. "/links.txt", "w+")
chr_out = io.open(arg[4] .. "/chromosomes.txt", "w+")
karyotype_out = io.open(arg[4] .. "/karyotype.txt", "w+")
bin_out = io.open(arg[4] .. "/bin.txt", "w+")

chr_prefix = arg[5]
bin_chr = arg[6]

has_bin = false

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
  karyotype_out:write("chr - " .. rn:get_seqid() .. " " .. rn:get_seqid()
                        .. " " .. rn:get_range():get_start() .. " "
                        .. rn:get_range():get_end() .. " ")
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
    m = rn:get_seqid():match(chr_prefix)
    if m and rn:get_seqid() ~= bin_chr then
      chr_out:write(m .. "\n")
    end
  end
  if self.is_ref then
    table.insert(ref_chr, rn:get_seqid())
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
  links_out:write(targetid .. " " .. tfrom .. " " .. tto .. " " .. refid .. " " .. rfrom .. " " .. rto .."\n")
end

if has_bin then
  bin_out:write(bin_chr .. "\t" .. table.concat(ref_chr, ";") .. "\n")
end
