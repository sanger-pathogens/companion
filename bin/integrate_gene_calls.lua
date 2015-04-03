#!/usr/bin/env gt

--[[
  Copyright (c) 2014-2015 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2015 Genome Research Ltd

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
  io.stderr:write("Selects the 'best' gene models from a pooled set of " ..
                  "GFF3 annotations.\n")
  io.stderr:write(string.format("Usage: %s <GFF with overlapping annotation>\n",
                                arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("SimpleChainer")

function get_weight(gene)
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
  end
  return gene:get_range():length() * fac
end

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
function stream:process_current_cluster()
  local bestset = nil
  local max = 0

-- XXX debug
--  print("before: " .. #self.curr_gene_set)
--  if #self.curr_gene_set > 50 then
--    for _,v in ipairs(self.curr_gene_set) do
--      v:accept(vis)
--    end
--  end

  -- keep only non-overlapping chain with highest weight
  bestset = SimpleChainer.new(self.curr_gene_set):chain()

-- XXX debug
--  print("after: " .. #bestset .. "\n")
--  for _,v in ipairs(bestset) do
--    v:accept(vis)
--  end


  for _,v in ipairs(bestset) do
    table.insert(self.outqueue, v)
  end
end

function stream:next_tree()
  local complete_cluster = false
  local mygn = nil

  if #self.outqueue > 0  then
    return table.remove(self.outqueue, 1)
  else
    complete_cluster = false
  end

  while not complete_cluster do
    mygn = self.instream:next_tree()
    if mygn then
      rval, err = pcall(GenomeTools_genome_node.get_type, mygn)
      if rval then
        --print(mygn)
        local fn = mygn
        local new_rng = mygn:get_range()
        if fn:get_type() == "gene" then
          if #self.curr_gene_set == 0 then
            table.insert(self.curr_gene_set, fn)
            self.curr_rng = new_rng
          else
            if self.last_seqid == fn:get_seqid()
                and self.curr_rng:overlap(new_rng) then
              table.insert(self.curr_gene_set, fn)
              self.curr_rng = self.curr_rng:join(new_rng)
            else
              -- no more overlap
              self:process_current_cluster()
              self.curr_gene_set = {}
              table.insert(self.curr_gene_set, fn)
              self.curr_rng = new_rng
              if #self.outqueue > 0  then
                outgn = table.remove(self.outqueue, 1)
                complete_cluster = true
              end
            end
          end
          self.last_seqid = mygn:get_seqid()
        end
      else
        -- no feature node
        self:process_current_cluster()
        self.curr_gene_set = {}
        table.insert(self.outqueue, fn)
        if #self.outqueue > 0  then
          outgn = table.remove(self.outqueue, 1)
          complete_cluster = true
        end
      end
    else
      -- end of annotation
      outgn = mygn
      break
    end
  end
  return outgn
end

stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
stream.idx = feature_index

out_stream = gt.gff3_out_stream_new(stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

