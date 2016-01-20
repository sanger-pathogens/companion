#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
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
  io.stderr:write(string.format("Usage: %s <GFF overlapping annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

function ReverseTable(t)
  local reversedTable = {}
  local itemCount = #t
  for k, v in ipairs(t) do
    reversedTable[itemCount + 1 - k] = v
  end
  return reversedTable
end

function clear_partial_attribs(cur_gene)
  if cur_gene:get_attribute("threeEndPartial") then
    cur_gene:remove_attribute("threeEndPartial")
  end
  if cur_gene:get_attribute("fiveEndPartial") then
    cur_gene:remove_attribute("fiveEndPartial")
  end
  if cur_gene:get_attribute("Start_range") then
    cur_gene:remove_attribute("Start_range")
  end
  if cur_gene:get_attribute("End_range") then
    cur_gene:remove_attribute("End_range")
  end
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
function stream:process_current_cluster()
  local gaps = {}

  -- count and collect gaps between scaffolds -- these always need to be broken
  for _,n in ipairs(self.curr_gene_set) do
    if n:get_type() == "gap" then
      table.insert(stream.outqueue, n)
      if n:get_attribute("gap_type")  == "between scaffolds" then
        table.insert(gaps, n)
      end
    end
  end

  -- no relevant gaps in cluster, so just pass along genes
  if #gaps == 0 then
    for _,n in ipairs(self.curr_gene_set) do
      if n:get_type() ~= "gap" then
        table.insert(stream.outqueue, n)
      end
    end
    return 0
  end

  -- iteratively split (sorted!) gene models at (sorted!) gaps
  for _,n in ipairs(self.curr_gene_set) do
    local cur_gene = n
    local outbuf = {}
    local skip_gene = false
    if n:get_type() ~= "gap" then
      for _,g in ipairs(gaps) do
        local rest = clone_cc(cur_gene)
        for c in cur_gene:children() do
          if c:get_range():get_start() > g:get_range():get_end() then
            cur_gene:remove_leaf(c)
          elseif c:get_range():overlap(g:get_range()) then
            if c:get_range():get_start() > g:get_range():get_start() - 1 then
              io.stderr:write("splitting would create invalid gene model, " ..
                              "skipping " .. cur_gene:get_attribute("ID") ..
                              "\n")
              skip_gene = true
            else
              local new_rng = gt.range_new(c:get_range():get_start(),
                                           g:get_range():get_start() - 1)
              c:set_range(new_rng)
              cur_gene:set_attribute("End_range",".,.")
              if c:get_strand() == "-" then
                cur_gene:set_attribute("fiveEndPartial", "true")
              else
                cur_gene:set_attribute("threeEndPartial", "true")
              end
            end
          end
        end
        for c in rest:children() do
          if c:get_range():get_end() < g:get_range():get_start() then
            rest:remove_leaf(c)
          elseif c:get_range():overlap(g:get_range()) then
            -- XXX make sure we don't create invalid genes
            if g:get_range():get_end() + 1 > c:get_range():get_end() then
              io.stderr:write("splitting would create invalid gene model, " ..
                              "skipping " .. cur_gene:get_attribute("ID") ..
                              "\n")
              skip_gene = true
            else
              local rest_rng = gt.range_new(g:get_range():get_end() + 1,
                                            c:get_range():get_end())
              c:set_range(rest_rng)
              rest:set_attribute("Start_range",".,.")
              if c:get_strand() == "-" then
                rest:set_attribute("threeEndPartial", "true")
              else
                rest:set_attribute("fiveEndPartial", "true")
              end
            end
          end
        end
        if not skip_gene then
          table.insert(outbuf, cur_gene)
        end
        cur_gene = rest
      end
      if not skip_gene then
        table.insert(outbuf, cur_gene)
      end

     if not skip_gene then
        -- repair CDS entries in gene models to fix frameshifts
        -- obtain and order CDS and gap features
        local cdslist = {}
        for _,v in ipairs(outbuf) do
          for q in v:children() do
            if q:get_type() == "CDS" or q:get_type() == "gap" then
              table.insert(cdslist, q)
            end
          end
        end
        -- reverse order for reverse strand features
        assert(cdslist[1]:get_type() == "CDS")
        local startphase = cdslist[1]:get_phase()
        if n:get_strand() == "-" then
          cdslist = ReverseTable(cdslist)
          cdslist[1]:set_phase(startphase)
        end

        -- assign correct strands
        assert(cdslist[1]:get_type() == "CDS")
        local phase = cdslist[1]:get_phase()
        local gaplen = 0
        -- determine starting phase (mostly 0, but you never know...)
        phase = (3 - (cdslist[1]:get_range():length() - phase) % 3) % 3
        for i = 2,#cdslist do
          if cdslist[i]:get_type() == "gap" then
            gaplen = cdslist[i]:get_range():length()
          elseif cdslist[i]:get_type() == "CDS" then
            -- adjust phase for shift introduced by gap length
            phase = (phase - (gaplen % 3)) % 3
            cdslist[i]:set_phase(phase)
            -- update running phase
            phase = (3 - (cdslist[i]:get_range():length() - phase) % 3) % 3
          end
        end

        -- split finished, pass along split gene models and gaps
        for _,v in ipairs(outbuf) do
          table.insert(stream.outqueue, v)
        end
      end
    end
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
        local fn = mygn
        local new_rng = mygn:get_range()
        if fn:has_child_of_type("CDS") or fn:get_type() == "gap" then
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
        else
          table.insert(self.outqueue, fn)
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
stream.vis = visitor
stream.idx = feature_index

out_stream = gt.gff3_out_stream_new(stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

