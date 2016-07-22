#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2015-2016 Genome Research Ltd

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

to_gene = {pseudogene = 'gene', pseudogenic_transcript = 'mRNA',
           pseudogenic_exon = 'CDS' }

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("optparse")

op = OptionParser:new({usage="%prog <options> merged_gff_file.gff3 sequence.fas",
                       oneliner="Adds pseudogene annotations to existing set of gene models.",
                       version="0.1"})
op:option{"-t", action='store', dest='threshold',
                help="minimum gene coverage ratio threshold for pseudogene "
                  .. "required to replace a gene (e.g. 0.6 for 60%)"}
op:option{"-m", action='store', dest='minratio',
                help="minimal 'stitching' gene coverage ratio "
                  .. "required to replace multiple genes (e.g. 0.8 for 80%)"}
op:option{"-x", action='store', dest='maxratio',
                help="maximal 'stitching' gene coverage ratio "
                  .. "required to replace multiple genes (e.g. 1.2 for 120%)"}
op:option{"-i", action='store_true', dest='ignore_frame_agreement',
                help="do not require frame agreement for pseudogene "
                  .. "calling (default: false)"}
op:option{"-r", action='store', dest='rejectfile',
                help="file to store 'rejected' pseudogene candidates in GFF3 format"}
op:option{"-l", action='store', dest='longest_until',
                help="pick candidate with median length if at least "
                  .. "<LONGEST_UNTIL> candidates per locus"}
op:option{"-d", action='store_true', dest='debug',
                help="output debug information on stderr"}
options,args = op:parse({threshold=0.6, minratio=0.8, maxratio=1.2,
                         debug=false, rejectfile=nil,
                         ignore_frame_agreement=false, longest_until=4})

function usage()
  op:help()
  os.exit(1)
end

if #args < 2 then
  usage()
end

local rejected_pseudogenes = {}
local rm = gt.region_mapping_new_seqfile_matchdescstart(args[2])

-- returns table of transcripts and their coding lengths
function get_transcript_lengths(g, transcript_type, cds_type)
  local transcripts = {}
  local has_transcripts = false
  for c in g:children() do
    if c:get_type() == tostring(transcript_type) then
      transcripts[c] = 0
      has_transcripts = true
      for c2 in c:children() do
        if c2:get_type() == tostring(cds_type) then
          transcripts[c] = transcripts[c] + c2:get_range():length()
        end
      end
    end
  end
  local max = 0
  local best = nil
  for k,v in pairs(transcripts) do
    if v > max then
      max = v
      best = k
    end
  end
  return has_transcripts, transcripts, best
end

-- compares two genes by the coding lengths of their longest transcript
function gene_cmp_by_length(a, b)
  local has_a, transcripts_a, longest_a =
                     get_transcript_lengths(a, "pseudogene", "pseudogenic_exon")
  local has_b, transcripts_b, longest_b =
                     get_transcript_lengths(b, "pseudogene", "pseudogenic_exon")
  if not has_a and has_b then
    return true
  elseif not has_b and has_a then
    return false
  else
    assert(has_a and has_b)
    assert(longest_a and longest_b)
    if transcripts_a[longest_a]  < transcripts_b[longest_b] then
      return true
    else
      return false
    end
  end
end

function frame_agreement(f1, f1type, f2, f2type)
  local okay = true
  if f1:get_strand() ~= f2:get_strand() then
    if options.debug then
      io.stderr:write("found strands in disagreement in frame check\n")
    end
    return false
  else
    local c1 = nil
    local c2 = nil
    if f1:get_strand() == '-' then
      -- minus strand: take stop position of last coding feature
      for c in f1:children() do
        if c:get_type() == f1type then
          c1 = c:get_range():get_end()
        end
      end
      for c in f2:children() do
        if c:get_type() == f2type then
          c2 = c:get_range():get_end()
        end
      end
    else
      -- plus/unknown strand: take start position of first coding feature
      for c in f1:children() do
        if c:get_type() == f1type then
          c1 = c:get_range():get_start()
          break
        end
      end
      for c in f2:children() do
        if c:get_type() == f2type then
          c2 = c:get_range():get_start()
          break
        end
      end
    end
    if not (c1 and c2) then
      if options.debug then
        io.stderr:write("could not find feature type to determine frame agreement\n")
      end
      return false
    end
    if (c1 % 3) ~= (c2 % 3) then
      if options.debug then
        io.stderr:write("different frames: " ..
                        tostring(f2:get_attribute("ID")) ..
                        " " .. c1 % 3 .. " vs " .. c2 % 3 ..
                        " (" .. tostring(f2:get_strand()) ..")\n")
      end
      return false
    else
      if options.debug then
        io.stderr:write("same frame " .. tostring(f2:get_attribute("ID")) ..
                        "\n")
      end
      return true
    end
  end
end

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
function stream:find_best_by_length(set)
  if not set or #set == 0 then
    return nil
  else
    local pseudogenes = {}
    for _,v in ipairs(set) do
      if v:get_type() == "pseudogene" then
        table.insert(pseudogenes, v)
      end
    end
    table.sort(pseudogenes, gene_cmp_by_length)
    -- select pseudogene with median length if more than a set number of hits,
    -- otherwise pick longest one
    local out = nil
    if #pseudogenes < options.longest_until then
      out = pseudogenes[#pseudogenes]
    else
      out = pseudogenes[math.max(1,math.floor(#pseudogenes/2))]
    end
    return out
  end
end
function stream:find_best_by_score(set)
  if not set or #set == 0 then
    return nil
  else
    local best = nil
    local maxscore = 0
    for _,v in ipairs(set) do
      if v:get_type() == "pseudogene" and v:get_score() then
          if tonumber(v:get_score()) > maxscore then
            best = v
            maxscore = v:get_score()
          end
      end
    end
    return best
  end
end
function stream:process_current_cluster()
  -- pull out genes
  local genes = {}
  for _,v in ipairs(self.curr_gene_set) do
    for c in v:children() do
      if c:get_type() == "CDS" then
        table.insert(genes, v)
        break
      end
    end
  end

  -- find 'best' pseudogene match
  best = self:find_best_by_length(self.curr_gene_set)

  if not best then
     -- no pseudogenes found, just output gene
    for _,g in ipairs(genes) do
      table.insert(self.outqueue, g)
    end
  else
    -- we have an aligned region
    local frameshift = best:get_attribute("has_frameshift")
    local internal_stop = best:get_attribute("has_internal_stop")
    local has_start = best:get_attribute("has_start")
    local has_stop = best:get_attribute("has_stop")
    if #genes == 0 then
      -- no overlaps, just create a new feature
      if frameshift == 'true' or internal_stop == 'true' then
        -- this is a new pseudogene in a yet undetected spot
        if options.debug then
          io.stderr:write("new pseudogene " .. best:get_attribute("ID") .. "\n")
        end
        table.insert(self.outqueue, best)
      else
        if has_start  == 'true' and has_stop == 'true' then
          -- this is a gene with protein homology missed by the other predictors
          if options.debug then
            io.stderr:write("new gene " .. best:get_attribute("ID") .. "\n")
          end
          table.insert(self.outqueue, deep_copy(best, nil, to_gene))
        else
          -- this is an aligned part with no mutations but no start/stop
          -- codons either
          table.insert(self.outqueue, best)
        end
      end
    elseif #genes == 1 then
      -- we have exactly one overlapping gene
      if genes[1]:get_range():get_start() == best:get_range():get_start()
             and genes[1]:get_range():get_end() == best:get_range():get_end() then
             -- best alignment matches gene perfectly, only output gene
        --io.stderr:write("HSP and gene identical " .. best:get_attribute("ID") .. "\n")
        table.insert(self.outqueue, genes[1])
      elseif not genes[1]:get_range():contains(best:get_range()) then
        if options.ignore_frame_agreement or frame_agreement(genes[1], "CDS", best, "pseudogenic_exon") then
          -- only handle in-frame cases
          local ratio = best:get_range():length()/genes[1]:get_range():length()
          if options.debug then
            io.stderr:write(ratio .. "\n")
          end
          if ratio > tonumber(options.threshold) then
            -- if in frame and length ratio is below threshold
            if options.debug then
              io.stderr:write("replaced gene " .. genes[1]:get_attribute("ID") .. "\n")
            end
            table.insert(self.outqueue, best)
          else
            table.insert(self.outqueue, genes[1])
            best:add_attribute("reject_reason", "ratio_" .. tostring(ratio))
            table.insert(rejected_pseudogenes, best)
          end
        else
          table.insert(self.outqueue, genes[1])
          best:add_attribute("reject_reason", "frame_agreement")
          table.insert(rejected_pseudogenes, best)
        end
      else
        -- otherwise ignore the pseudogene
        table.insert(self.outqueue, genes[1])
      end
    else
      -- we have more than one overlapping gene
      -- check for gaps
      local has_gap = false
      if best:extract_sequence("pseudogene", false, rm):match("[nN]+") then
        has_gap = true
      end
      if not has_gap then
        -- count frame switches
        local nof_pexons = 0
        for c in best:children() do
          if c:get_type() == 'pseudogenic_exon' then
            nof_pexons = nof_pexons + 1
          end
        end
        -- determine size of gene
        local total_gene_length = 0
        for _,g in ipairs(genes) do
          total_gene_length = total_gene_length + g:get_range():length()
          --  table.insert(self.outqueue, g)
        end
        local ratio = (total_gene_length)/(best:get_range():length())
        if ratio <= tonumber(options.maxratio)
                                    and ratio >= tonumber(options.minratio) then
          if options.debug then
            io.stderr:write("replaced genes (")
            for _,g in ipairs(genes) do
              io.stderr:write(g:get_attribute("ID") .. " " )
            end
            io.stderr:write(") with " .. best:get_attribute("ID") .. "\n")
          end
          if frameshift == 'true' or internal_stop == 'true' then
            table.insert(self.outqueue, best)
          else
            if has_start  == 'true' and has_stop == 'true' then
              -- this is a gene with protein homology missed by the other predictors
              table.insert(self.outqueue, deep_copy(best, nil, to_gene))
            else
              table.insert(self.outqueue, best)
            end
          end
        else
          -- keep genes
          for _,g in ipairs(genes) do
            table.insert(self.outqueue, g)
          end
          best:add_attribute("reject_reason","ratio_" .. tostring(ratio))
          table.insert(rejected_pseudogenes, best)
        end
      else
        -- keep genes
        for _,g in ipairs(genes) do
          table.insert(self.outqueue, g)
        end
        if options.debug then
          io.stderr:write("has gaps " .. best:get_attribute("ID") .. "\n")
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
        if fn:get_type() == "gene" or fn:get_type() == "pseudogene" then
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

stream.instream = gt.gff3_in_stream_new_sorted(args[1])
stream.idx = feature_index

out_stream = gt.gff3_out_stream_new(stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

if options.rejectfile then
  local arr_in_stream = gt.custom_stream_new_unsorted()
  arr_in_stream.arr = rejected_pseudogenes
  function arr_in_stream:next_tree()
    if #self.arr > 0  then
      return table.remove(self.arr, 1)
    else
      return nil
    end
  end
  local rejected_out_stream = gt.gff3_out_stream_new(arr_in_stream,
                                                     options.rejectfile)
  local gn = rejected_out_stream:next_tree()
  while (gn) do
    gn = rejected_out_stream:next_tree()
  end
end
