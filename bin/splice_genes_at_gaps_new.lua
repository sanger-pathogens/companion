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
  io.stderr:write(string.format("Usage: %s <GFF w/ genes + gaps>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
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

  -- count and collect intra-scaffold gaps
  for _,n in ipairs(self.curr_gene_set) do
    if n:get_type() == "gap"
        and n:get_attribute("gap_type") == "within scaffold" then
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

  -- sort gaps by start position just to be safe
  table.sort(gaps, function (a,b)
    local a_s = a:get_range():get_start()
    local b_s = b:get_range():get_start()
    if a_s < b_s then
      return true
    else
      return false
    end
  end)

  -- adjust CDS according to gaps
  for _,n in ipairs(self.curr_gene_set) do
    local cur_gene = n
    local nof_cds = 0
    local skip_gene = false

    if n:get_type() ~= "gap" then
      for child in cur_gene:children() do
        if skip_gene then
          break
        end
        -- we are at a transcript
        if child:get_type() == "mRNA" then
          local mrna = child

          -- merge gaps if necessary to address small scattered gap clusters
          local relevant_gaps = {}
          local last_rng = nil
          local relevant_s = nil
          local relevant_e = nil
          for _,g in ipairs(gaps) do
            this_rng = g:get_range()
            if not relevant_s then
              relevant_s = this_rng:get_start()
            end
            if last_rng and (this_rng:get_start() - last_rng:get_end()) > 5 then
              relevant_e = last_rng:get_end()
              local new_gap = gt.feature_node_new(g:get_seqid(), "gap",
                                                  relevant_s, relevant_e,
                                                  g:get_strand())
              table.insert(relevant_gaps, new_gap)
              relevant_s = this_rng:get_start()
            end
            last_rng = this_rng
          end
          if relevant_s then
            relevant_e = last_rng:get_end()
            local new_gap = gt.feature_node_new(n:get_seqid(), "gap",
                                                    relevant_s, relevant_e,
                                                    n:get_strand())
            table.insert(relevant_gaps, new_gap)
          end

          -- handle gaps for this transcript
          for _,g in ipairs(relevant_gaps) do
            if skip_gene then
              break
            end
            local gap_rng = g:get_range()
            -- traverse all CDS
            for cds in mrna:children() do
              if cds:get_type() == "CDS" then
                local cds_rng = cds:get_range()
                nof_cds = nof_cds + 1
                -- is this CDS affected by the gap?
                if cds_rng:overlap(gap_rng) and not gap_rng:contains(cds_rng) then
--                  io.stderr:write(">>>>> " .. tostring(n:get_attribute("ID")) .. "\n")
--                  io.stderr:write(">>>>> overlap " .. tostring(g)
--                                  .. "  " .. tostring(cds) .. "\n")
                  -- handle overlap situation:
                  if cds_rng:contains(gap_rng) then
                    -- ================  cds
                    --       ====        gap
--                    io.stderr:write(">>>>> containment\n")
                    local shift = 0
                    if (gap_rng:length() % 3) ~= 0 then
                      shift = 3 - (gap_rng:length() % 3)
                    end
--                    io.stderr:write(">>> shift " .. shift .. "\n")
                    local new_rng = nil
                    local rest_rng = nil
                    if n:get_strand() == "-" then
                      if (cds_rng:get_start() > gap_rng:get_start() - 1 - shift) or
                         (gap_rng:get_end() + 1 > cds_rng:get_end()) then
                        io.stderr:write("splitting would create invalid gene model, " ..
                              "skipping " .. cur_gene:get_attribute("ID") ..
                              "\n")
                        skip_gene = true
                        break
                      else
                        new_rng = gt.range_new(cds_rng:get_start(),
                                                   gap_rng:get_start() - 1 - shift)
                        rest_rng = gt.range_new(gap_rng:get_end() + 1,
                                                    cds_rng:get_end())
                      end
                    else
                      if (cds_rng:get_start() > gap_rng:get_start() - 1) or
                         (gap_rng:get_end() + 1 + shift > cds_rng:get_end()) then
                        io.stderr:write("splitting would create invalid gene model, " ..
                              "skipping " .. cur_gene:get_attribute("ID") ..
                              "\n")
                        skip_gene = true
                        break
                      else
                        new_rng = gt.range_new(cds_rng:get_start(),
                                                   gap_rng:get_start() - 1)
                        rest_rng = gt.range_new(gap_rng:get_end() + 1 + shift,
                                                    cds_rng:get_end())
                      end
                    end
                    assert(new_rng)
                    assert(rest_rng)
                    cds:set_range(new_rng)
                    local new_cds = gt.feature_node_new(cds:get_seqid(), "CDS",
                                                        rest_rng:get_start(),
                                                        rest_rng:get_end(),
                                                        cds:get_strand())
                    mrna:add_child(new_cds)
                    new_cds:set_source(cds:get_source())
                  elseif gap_rng:get_start() <= cds_rng:get_start() and
                         gap_rng:get_end() >= cds_rng:get_start() then
                    --    =============  cds
                    --  ======           gap
--                    io.stderr:write(">>> case 1\n")
                    local new_rng = gt.range_new(gap_rng:get_end() + 1,
                                                 cds_rng:get_end())
                    cds:set_range(new_rng)
                    -- TODO: fixup phases
                  elseif gap_rng:get_end() >= cds_rng:get_end() and
                         gap_rng:get_start() <= cds_rng:get_end() then
                    --  =============    cds
                    --            =====  gap
--                    io.stderr:write(">>> case 2\n")
                    local new_rng = gt.range_new(cds_rng:get_start(),
                                                 gap_rng:get_start() - 1)
                    cds:set_range(new_rng)
                    -- TODO: fixup phases
                  else
                    -- any other case
                    -- e.g. containment of CDS in gap (should not happen)
                    if gap_rng:contains(cds_rng) then
                      io.stderr:write("gap (" .. tostring(g)
                                      .. ") contains CDS ( " .. tostring(cds)
                                      .. "), removing this CDS")
                      mrna:remove_leaf(cds)
                    end
                  end
                end
              end -- if CDS
            end -- for mRNA children
          end  -- for gaps
        end -- if mRNA

      end -- for gene children
    end
    if not skip_gene then
      table.insert(self.outqueue, n)
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

