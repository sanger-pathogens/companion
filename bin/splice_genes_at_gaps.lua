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
  io.stderr:write(string.format("Usage: %s <GFF w/ genes + gaps>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function ReverseTable(t)
  local reversedTable = {}
  local itemCount = #t
  for k, v in ipairs(t) do
    reversedTable[itemCount + 1 - k] = v
  end
  return reversedTable
end

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
function stream:process_current_cluster()
  local gaps = {}

  -- count and collect intra-scaffold gaps
  for _,n in ipairs(self.curr_gene_set) do
    if n:get_type() == "gap"
        and n:get_attribute("gap_type")  == "within scaffold" then
      table.insert(gaps, n)
    end
  end

  -- no inter-contig gaps in cluster, so just pass along genes
  if #gaps == 0 then
    for _,n in ipairs(self.curr_gene_set) do
      table.insert(stream.outqueue, n)
    end
    return 0
  end

  -- split CDS of each gene
  for _,n in ipairs(self.curr_gene_set) do
    local cur_gene = n
    local outbuf = {}
    local cur_cds = nil
    if n:get_type() ~= "gap" then
      for _,g in ipairs(gaps) do
        local mrna = nil
        for c in cur_gene:children() do
          if c:get_type() == "mRNA" then
            mrna = c
          end
          if c:get_type() == "CDS" then
            local new_cds = nil
            if c:get_range():overlap(g:get_range()) then
              if not cur_cds then
                cur_cds = c
              end
              local new_rng = gt.range_new(c:get_range():get_start(),
                                           g:get_range():get_start() - 1)
              local rest_rng = gt.range_new(g:get_range():get_end() + 1,
                                            c:get_range():get_end())
              c:set_range(new_rng)
              new_cds = gt.feature_node_new(c:get_seqid(), "CDS",
                                                  rest_rng:get_start(),
                                                  rest_rng:get_end(),
                                                  c:get_strand())
              mrna:add_child(new_cds)
              new_cds:set_source(c:get_source())
              table.insert(outbuf, cur_cds)
              table.insert(outbuf, g)
            end
            cur_cds = new_cds
          end
        end
        table.insert(outbuf, cur_cds)
      end



      -- repair CDS entries in gene models to fix frameshifts
      local cdslist = {}
      for _,v in ipairs(outbuf) do
        for q in v:children() do
          if q:get_type() == "CDS" or q:get_type() == "gap" then
            table.insert(cdslist, q)
          end
        end
      end

      if #cdslist > 0 then
        n:set_attribute("internalGap", "true")

        -- reverse order for reverse strand features
        assert(cdslist[1]:get_type() == "CDS")
        local startphase = cdslist[1]:get_phase()
        if n:get_strand() == "-" then
          cdslist = ReverseTable(cdslist)
          cdslist[1]:set_phase(startphase)
        end

        -- update phases
        assert(cdslist[1]:get_type() == "CDS")
        local gaplen = 0
        local nofgaps = 0
        for i = 2,#cdslist do
          if cdslist[i]:get_type() == "gap" then
            gaplen = cdslist[i]:get_range():length()
            nofgaps = nofgaps + 1
          elseif cdslist[i]:get_type() == "CDS" then
            -- adjust coords for shift introduced by gap length, if necessary
            local shift = 0
            if (gaplen % 3) ~= 0 then
              shift = 3 - (gaplen % 3)
            end
            local newrng = nil
            if cdslist[i]:get_strand() == "-" then
              newrng = gt.range_new(cdslist[i]:get_range():get_start(),
                                    cdslist[i]:get_range():get_end() - shift)
            else
              newrng = gt.range_new(cdslist[i]:get_range():get_start() + shift,
                                    cdslist[i]:get_range():get_end())
            end
            cdslist[i]:set_range(newrng)
          end
        end

        -- take care of boundary case where after frameshift adjustment, the
        -- terminal CDS ends up too short to still contain a stop codon on
        -- its own
        if #cdslist > 2 and cdslist[#cdslist]:get_range():length() < 3 then
          local last_cds_range = cdslist[#cdslist-2]:get_range()
          -- in this case, drop it
          n:remove_leaf(cdslist[#cdslist])
          -- if there was only one gap, remove internal gap flag
          if nofgaps == 1 then
            n:remove_attribute("internalGap")
          end
          -- add appropriate partial flag
          if n:get_strand() == "-" then
            n:set_attribute("Start_range",".,.")
          else
            n:set_attribute("End_range",".,.")
          end
          -- readjust gene coordinates across this CC
          for this_child in n:children() do
            if this_child:get_range():overlap(last_cds_range) then
              if last_cds_range:get_start() > this_child:get_range():get_start() then
                local rng = gt.range_new(last_cds_range:get_start(),
                                         this_child:get_range():get_end())
                this_child:set_range(rng)
              end
              if last_cds_range:get_end() < this_child:get_range():get_end() then
                local rng = gt.range_new(this_child:get_range():get_start(),
                                         last_cds_range:get_end())
                this_child:set_range(rng)
              end
            end
          end
        end
      end
    end
    table.insert(self.outqueue, n)
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
        if fn:get_type() == "gene" or fn:get_type() == "gap" then
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

