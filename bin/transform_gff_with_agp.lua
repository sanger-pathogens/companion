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
  io.stderr:write(string.format("Usage: %s <base GFF> <AGP/GFF> "
                                .. "<base FASTA> <target FASTA> "
                                .. "[gaps] [check]\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

gff = arg[1]
agp = arg[2]
gaps = arg[5]
if arg[3] then
  bkeys, bseqs = get_fasta(arg[3])
end
if arg[4] then
  tkeys, tseqs = get_fasta(arg[4])
end

function get_seq(seqs, id)
  for k,v in pairs(seqs) do
    if k:match("^"..id.."$") then
      return v
    end
  end
  io.stderr:write("could not find sequence matching object " .. id
                    .. " in input\n")
  os.exit(1)
end

gapqueue = {}

-- load and index AGP
mappings = {}

for l in io.lines(arg[2]) do
  if string.sub(l, 1, 1) ~= '#' then
    obj, obj_s, obj_e, part_n, type, c6, c7, c8, c9 = unpack(split(l, "%s+"))
    if type == 'N' then
      fn = gt.feature_node_new(obj, "gap", obj_s, obj_e, ".")
      fn:add_attribute("estimated_length", obj_e-obj_s+1)
      if gaps then
        fn:add_attribute("gap_type", gaps)
        table.insert(gapqueue, fn)
      end
    elseif type == 'U' then
      fn = gt.feature_node_new(obj, "gap", obj_s, obj_e, ".")
      if gaps then
        fn:add_attribute("gap_type", gaps)
        table.insert(gapqueue, fn)
      end
    else
      mappings[c6] = {obj=obj, obj_s=tonumber(obj_s),
                                        obj_e=tonumber(obj_e),
                                        s_s=tonumber(c7), s_e=tonumber(c8),
                                        strand=c9}
    end
  end
end

-- visitor for feature transformation
visitor = gt.custom_visitor_new()
function visitor:visit_feature(fn)
  seqid = fn:get_seqid() --:gsub("%.%d+$","")
  -- is this sequence/fragment part of a layout?
  if mappings[seqid] then
    -- yes, transform all children
    for node in fn:children() do
      -- assign new seqid
      node:change_seqid(mappings[seqid].obj)
      rng = node:get_range()
      -- was the fragment used in reverse orientation?
      if mappings[seqid].strand == "-" then
        -- transform coordinates for reverse
        new_rng = gt.range_new((mappings[seqid].s_e - rng:get_end() + 1)
                                 + mappings[seqid].obj_s - 1,
                               (mappings[seqid].s_e - rng:get_start() + 1)
                                 + mappings[seqid].obj_s - 1)
        node:set_range(new_rng)
        -- flip strands
        if node:get_strand() == "+" then
          node:set_strand("-")
        elseif node:get_strand() == "-" then
          node:set_strand("+")
        end
      else
        -- otherwise just transform coordinates by offsetting
        new_rng = gt.range_new(rng:get_start() + mappings[seqid].obj_s - 1,
                               rng:get_end() + mappings[seqid].obj_s - 1)
        node:set_range(new_rng)
      end
      if arg[6] then   -- optionally, check seqs -- must stay the same!
        seq1 = get_seq(bseqs, seqid):sub(rng:get_start(), rng:get_end())
        seq2 = get_seq(tseqs, mappings[seqid].obj):sub(new_rng:get_start(),
                       new_rng:get_end())
        if  mappings[seqid].strand == "-" then
          seq2 = revcomp(seq2)
        end
        if seq1 ~= seq2 then
          print(seq1)
          print(seq2)
        end
        assert(seq1 == seq2)
      end
    end
  else
    -- this sequence is unassembled, do not change its coordinates
    io.stderr:write("no mapping for seqid " .. seqid .. "\n")
    --for node in fn:children() do
    --  node:change_seqid(seqid)
    --end
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
visitor_stream.instream = gt.gff3_in_stream_new_sorted(gff)
visitor_stream.vis = visitor
visitor_stream.queue = gapqueue
visitor_stream.gaps = false
function visitor_stream:next_tree()
  local node = nil
  if not self.gaps then
    node = self.instream:next_tree()
    if node then
      node:accept(self.vis)
    else
      self.gaps = true
    end
  end
  if self.gaps then
    if table.getn(self.queue) > 0 then
      return table.remove(self.queue)
    else
      return nil
    end
  end
  return node
end

-- attach to output stream and pull
out_stream = gt.gff3_out_stream_new(visitor_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
