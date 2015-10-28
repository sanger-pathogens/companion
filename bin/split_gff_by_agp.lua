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
                                .. "[check]\n" , arg[0]))
  os.exit(1)
end

if #arg < 4 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

gff = arg[1]
agp = arg[2]
bkeys, bseqs = get_fasta(arg[3])
tkeys, tseqs = get_fasta(arg[4])

function get_seq(seqs, id)
  for k,v in pairs(seqs) do
    if k:match(id) then
      return v
    end
  end
  io.stderr:write("could not find sequence matching object " .. id
                    .. " in input\n")
  os.exit(1)
end

-- load and index AGP
-- TODO: make this a proper 'class'
components = {}
-- collect components, AGP must be sorted!
for l in io.lines(arg[2]) do
  if string.sub(l, 1, 1) ~= '#' then
    obj, obj_s, obj_e, part_n, type, c6, c7, c8, c9 = unpack(split(l, "%s+"))
    if type == 'U' or type == 'N' then
      if not components[obj] then
        components[obj] = {}
      end
      table.insert(components[obj], {type="gap", seqid=obj,
                               start=tonumber(obj_s),
                               stop=tonumber(obj_e)})
    else
      if not components[obj] then
        components[obj] = {}
      end
      table.insert(components[obj], {type="seq", seqid=obj, start=tonumber(obj_s),
                                   stop=tonumber(obj_e), subject=c6,
                                   s_s=tonumber(c7), s_e=tonumber(c8),
                                   strand=c9})
    end
  end
end
-- connect gaps to flanking seqs
for k,v in pairs(components) do
  for i = 1,#v do
    if v[i].type == "gap" then
      assert(v[i+1].type == "seq")
      v[i].next = v[i+1]
      assert(v[i-1].type == "seq")
      v[i].prev = v[i-1]
    end
  end
end

hdrcache = {}
function get_real_seqhdr(hdr)
  if not hdrcache[hdr] then
    for _,v in ipairs(tkeys) do
      v = v:split("%s")[1]
      if v == hdr then
        hdrcache[hdr] = v
        return v
      end
    end
    return hdr
  else
    return hdrcache[hdr]
  end
end

gff3vis = gt.gff3_visitor_new()

-- visitor for feature transformation
visitor = gt.custom_visitor_new()
function visitor:visit_feature(fn)
  local seqid = fn:get_seqid()
  -- is this sequence/fragment part of a layout?
  if components[seqid] then
    --print("have components for seqid " .. seqid)
    for _,m in ipairs(components[seqid]) do
      local next_gene = nil
      local rng = gt.range_new(m.start, m.stop)
      if m.type == "seq" and rng:contains(fn:get_range()) then
        -- fully contained, no splitting necessary
        for c in fn:children() do
          old_rng = c:get_range()
          if m.strand == "-" then
            new_rng = gt.range_new(m.s_e - (c:get_range():get_end() - m.start),
                                   m.s_e - (c:get_range():get_start() - m.start))
            c:set_range(new_rng)
            c:change_seqid(get_real_seqhdr(m.subject))
            if c:get_strand() == "+" then
              c:set_strand("-")
            elseif c:get_strand() == "-" then
              c:set_strand("+")
            end
          else
            -- otherwise just transform coordinates by offsetting
            new_rng = gt.range_new(c:get_range():get_start() - m.start + 1,
                                   c:get_range():get_end() - m.start + 1)
            c:set_range(new_rng)
            c:change_seqid(get_real_seqhdr(m.subject))
          end
          if arg[5] then   -- optionally, check seqs -- must stay the same!
            seq1 = get_seq(bseqs, seqid):sub(old_rng:get_start(),
                                             old_rng:get_end())
            seq2 = get_seq(tseqs, m.subject):sub(c:get_range():get_start(),
                                                 c:get_range():get_end())
            if m.strand == "-" then
              seq2 = revcomp(seq2)
            end
            assert(seq1 == seq2)
          end
        end
        fn:accept(gff3vis)
      end
    end
  else
    -- this sequence is unassembled, do not change its coordinates
    newseqid = get_real_seqhdr(fn:get_seqid())
    for c in fn:children() do
      c:change_seqid(newseqid)
    end
    io.stderr:write("no mapping for seqid " .. seqid .. "\n")
    fn:accept(gff3vis)
  end

  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
visitor_stream.instream = gt.gff3_in_stream_new_sorted(gff)
visitor_stream.vis = visitor
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

out_stream = visitor_stream --  gt.gff3_out_stream_new(visitor_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end


  --         if (seq1 ~= seq2) then
  --            print("1 (" .. seqid .. ", " .. old_rng:get_start()
  --                     .. "-" ..  old_rng:get_end() .. "):" .. seq1 ..
  --                  "\n2:" ..  m.subject .. ", " .. c:get_range():get_start()
  --                     .. "-" ..  c:get_range():get_end() .. "):" .. seq2
  --                     .. "\n\n")
  --            os.exit(1)
  --          end
