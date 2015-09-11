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
  io.stderr:write(string.format("Usage: %s <GFF with CDS annotations> <prefix> <chromosome pattern>\n" , arg[0]))
  os.exit(1)
end

if #arg < 3 then
  usage()
  os.exit(1)
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

increment = 100

cv = gt.custom_visitor_new()
cv.queue = {}
cv.numbers = {}
cv.suffixes = {}
cv.suffixes.gene = ''
function cv:visit_feature(fn)
  local transcript_counter = 0
  -- leave gaps/contigs untouched
  if fn:get_type() == "gap" or fn:get_type() == "contig" then
    return 0
  end
  local seqid = fn:get_seqid()
  -- default for 'bin' seqs
  local chr = "00"
  m = seqid:match(arg[3])
  if m then
    if tonumber(m) then
      chr = string.format("%02d", tonumber(m))
    else
      chr = m
    end
  end
  -- initialize counter for feature numbers
  if not self.numbers[chr] then
    self.numbers[chr] = 5000
  end
  -- determine chromosome number
  local base_id = string.format("%s_%s%07d", arg[2], chr, self.numbers[chr])
  for n in fn:get_children() do
    local this_id = base_id
    --delete "Name"s
    if n:get_type() ~= "contig" and n:get_attribute("Name") then
      n:remove_attribute("Name")
    end
    -- for not yet seen features...
    if n:get_attribute("_done") then
      n:remove_attribute("_done")
    else
      n:set_attribute("ID", this_id)
      -- process transcript, TODO: make this SO aware
      if n:get_type() == "mRNA" or n:get_type() == "pseudogenic_transcript" then
        transcript_counter = transcript_counter + 1
        n:set_attribute("ID", base_id .. "." .. transcript_counter)
        transcript_region_ctr = {}
        local mrna_id = n:get_attribute("ID")
        local rng = nil
        -- collect range covered by CDSs etc.
        for c in n:get_children() do
          local t = c:get_type()
          if t == "CDS" or t == "pseudogenic_exon" or t == "exon" or t:match("_codon") then
            -- these are `coding' regions
            if t == "CDS" or t == "pseudogenic_exon" then
              if not rng then
                rng = c:get_range()
              else
                rng = rng:join(c:get_range())
              end
            end
            -- count number of subfeatures per transcript
            if not transcript_region_ctr[t] then
              transcript_region_ctr[t] = 1
            else
              transcript_region_ctr[t] = transcript_region_ctr[t] + 1
            end
            -- construct subfeature ID
            this_id = mrna_id .. ":" .. t ..":" .. transcript_region_ctr[t]
            local typestring = t
            if self.suffixes[typestring] then
                typestring = self.suffixes[typestring]
            else
              typestring = ":" .. typestring
            end
            -- set it
            c:set_attribute("ID", this_id)
            -- mark as already handled
            c:set_attribute("_done", "true")
          end
          n:set_attribute("_done", "true")
        end
        if rng then
          -- enqueue new polypeptide feature
          pp = gt.feature_node_new(n:get_seqid(), "polypeptide",
                                   rng:get_start(), rng:get_end(),
                                   n:get_strand())
          pp:add_attribute("Derives_from", mrna_id)
          pp:add_attribute("ID", string.format("%s:pep", mrna_id))
          table.insert(cv.queue, 1, pp)
        end
      else
        -- any other feature, just number and tag it
        local typestring = n:get_type()
        if self.suffixes[typestring] then
            typestring = self.suffixes[typestring]
        else
          typestring = ":" .. typestring
        end
        n:set_attribute("ID", this_id .. typestring)
      end
    end
  end
  -- clean up
  for c in fn:get_children() do
    if c:get_attribute("_done") then
      c:remove_attribute("_done")
    end
  end
  self.numbers[chr] = self.numbers[chr] + increment
  return 0
end

vis_stream = gt.custom_stream_new_unsorted()
vis_stream.queue = cv.queue
vis_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vis_stream.visitor = cv
function vis_stream:next_tree()
  if table.getn(self.queue) > 0 then
    return table.remove(self.queue)
  else
    local node = self.instream:next_tree()
    if node then
      node:accept(self.visitor)
    end
    return node
  end
end

out_stream = gt.gff3_out_stream_new(vis_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end