#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd

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
  io.stderr:write(string.format("Usage: %s <GFF annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

snap_exons = {}

for l in io.lines(arg[1]) do
  if l:sub(1,1) ~= '#' then
    seqid, src, type, start, stop, score, strand, _,  id = unpack(split(l, "\t"))
    if id then
      if not snap_exons[id] then
        snap_exons[id] = {}
      end
      table.insert(snap_exons[id], {seqid, src, type, start, stop, score, strand})
    end
  end
end

queue = {}
geneno = 1

for k,v in pairs(snap_exons) do
  local totalrange = nil
  exnr = 1
  cds = {}
  exons = {}
  for _,e in ipairs(v) do
    seqid, src, type, start, stop, score, strand = unpack(e)
    if start <= stop then
      local thisrng = gt.range_new(tonumber(start), tonumber(stop))
      if not totalrange then
        totalrange = thisrng
      else
        totalrange = totalrange:join(thisrng)
      end
      node = gt.feature_node_new(seqid, "CDS", tonumber(start), tonumber(stop), strand)
      node:set_attribute("ID",  "gene" .. geneno .. ".CDS"..exnr)
      node:set_source(src)
      table.insert(cds, node)
      node = gt.feature_node_new(seqid, "exon", tonumber(start), tonumber(stop), strand)
      node:set_attribute("ID", "gene" .. geneno .. ".exon"..exnr)
      node:set_source(src)
      table.insert(exons, node)
      exnr = exnr + 1
    end
  end
  if totalrange then
    gene = gt.feature_node_new(seqid, "gene",
                               totalrange:get_start(),
                               totalrange:get_end(),
                               strand)
    gene:set_attribute("ID", "gene"..geneno)
    gene:set_source(src)
    mrna = gt.feature_node_new(seqid, "mRNA",
                               totalrange:get_start(),
                               totalrange:get_end(),
                               strand)
    mrna:set_attribute("ID", "mrna"..geneno)
    mrna:set_source(src)
    for _,n in pairs(cds) do
      mrna:add_child(n)
    end
    for _,n in pairs(exons) do
      mrna:add_child(n)
    end
    gene:add_child(mrna)
    table.insert(queue, gene)
    geneno = geneno + 1
  end
end

vis_stream = gt.custom_stream_new_unsorted()
vis_stream.queue = queue
function vis_stream:next_tree()
  if table.getn(self.queue) > 0 then
    return table.remove(self.queue, 1)
  else
    return nil
  end
end

out_stream = gt.gff3_out_stream_new(vis_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end