#!/usr/bin/env gt

--[[
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
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
  io.stderr:write(string.format("Usage: %s <GFF overlapping annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

function powerset(s)
   local t = {{}}
   for i = 1, #s do
      for j = 1, #t do
         t[#t+1] = {s[i],unpack(t[j])}
      end
   end
   return t
end

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
-- find non-overlapping subset with maximal total length
function stream:process_current_cluster()
  local bestset = nil
  local max = 0
  -- enumerate all subsets
  for _,s in ipairs(powerset(self.curr_gene_set)) do
    local overlaps = false
    local rng = nil
    -- check for overlaps
    for _,v in ipairs(s) do
      if not rng then
        rng = v:get_range()
      else
        if rng:overlap(v:get_range()) then
          overlaps = true
          break
        else
          rng = rng:join(v:get_range())
        end
      end
    end
    -- score this cluster
    if not overlaps then
      local sum = 0
      if bestset == nil then
        bestset = s
      end
      for _,v in ipairs(s) do
        -- apply reward for being annotated by RATT
        local fac = 1
        if v:get_source() == "RATT" then
          fac = 3
        end
        --if self.idx:has_seqid(v:get_seqid()) then
          --fac = 1
         -- ptu = self.idx:get_features_for_range(v:get_seqid(), v:get_range())
         -- if ptu and #ptu == 1 and ptu[1]:get_strand() ~= "?" then
         --   if ptu[1]:get_strand() == v:get_strand() then
         --     fac = 2.0
         --   else
         --     fac = -1
         --   end
         -- end
        --end
        sum = sum + (v:get_range():length() * fac)
      end
      if sum > max then
        max = sum
        bestset = s
      end
    end
  end
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
            if self.last_seqid == fn:get_seqid() and self.curr_rng:overlap(new_rng) then
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

-- create feature index
--feature_index = gt.feature_index_memory_new()
-- add GFF3 file to index
--feature_index:add_gff3file(arg[2])

stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
stream.idx = feature_index

out_stream = gt.gff3_out_stream_new(stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

