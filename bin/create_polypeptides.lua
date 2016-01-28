#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2016 Genome Research Ltd

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

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("optparse")

default_length = 4
default_increment = 100
default_startval = 5000
default_seed = os.time()

op = OptionParser:new({usage="%prog [options] <GFF with CDS annotations> <prefix> <chromosome pattern>",
                       oneliner="Adds polypeptide features and systematic IDs to raw gene features.",
                       version="0.1"})
op:option{"-r", action='store_true', dest='random_alphanumeric',
                help="use random alphanumeric locus_tags/geneIDs"}
op:option{"-l", action='store', dest='random_length', default=default_length,
                help="length of random alphanumeric IDs"}
op:option{"-i", action='store', dest='increment', default=default_increment,
                help="increment between consecutive gene IDs"}
op:option{"-s", action='store', dest='startval', default=default_startval,
                help="start value per chromosome for gene IDs"}
op:option{"-e", action='store', dest='seed', default=default_seed,
                help="seed value for PRNG"}
options,args = op:parse({random_alphanumeric = false,
                         random_length = tonumber(default_length),
                         increment = tonumber(default_increment),
                         startval = tonumber(default_startval),
                         seed = tonumber(default_seed)})

function usage()
  op:help()
  os.exit(1)
end

if #args ~= 3 then
  usage()
end

-- implementation for numeric gene IDs, with chromosome numbers
IDProviderNumeric = {}
IDProviderNumeric.__index = IDProviderNumeric
function IDProviderNumeric.new(startval, increment, prefix, chr_pattern)
  idp = {}
  setmetatable(idp, IDProviderNumeric)
  idp.startval = startval
  idp.increment = increment
  idp.prefix = prefix
  idp.chr_pattern = chr_pattern
  idp.numbers = {}
  return idp
end
function IDProviderNumeric:next_for_feature(feature)
  local seqid = feature:get_seqid()
  local m = seqid:match(self.chr_pattern)
  local chr = '00'
  if m then
    if tonumber(m) then
      chr = string.format("%02d", tonumber(m))
    else
      chr = m
    end
  end
  assert(chr)
  if not self.numbers[chr] then
    self.numbers[chr] = startval
  end
  local rval = string.format("%s_%s%07d", self.prefix, chr, self.numbers[chr])
  self.numbers[chr] = self.numbers[chr] + increment
  return rval
end

-- implementation for random alphanumeric IDs
IDProviderRandomAlpha = {}
IDProviderRandomAlpha.__index = IDProviderRandomAlpha
function IDProviderRandomAlpha.new(prefix, length, seed)
  idp = {}
  setmetatable(idp, IDProviderRandomAlpha)
  idp.length = length
  idp.alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
  idp.maxcount = math.pow(idp.alpha:len(), length)
  idp.ids = {}
  idp.count = 0
  idp.prefix = prefix
  math.randomseed(seed)
  return idp
end
function IDProviderRandomAlpha:next_for_feature(feature)
  local id = ""
  if self.count + 1 > self.maxcount then
    error("can't handle more than " .. self.maxcount
            .. " distinct IDs with length " .. self.length)
  end
  -- This is just sample+test, but should still be fast enough for typical
  -- feature numbers as |idp.alpha|^4 (default) is 1679616 already. If we get
  -- annotations with more features, we have more serious problems anyway.
  while id == "" or self.ids[id] do
    id = ""
    for i = 1, self.length do
      local pos = math.random(1, self.alpha:len())
      id = id .. self.alpha:sub(pos, pos)
    end
  end
  local rval = string.format("%s_%s", self.prefix, id)
  self.ids[id] = true
  self.count = self.count + 1
  return rval
end

increment = tonumber(options.increment)
if increment < 1 then
  error("increment value must be >0")
end
startval = tonumber(options.startval)
if startval < 0 then
  error("startval value must be >0")
end
if options.random_length and options.random_alphanumeric then
  random_length = tonumber(options.random_length)
  if random_length < 4 then
    error("startval value must be >3")
  end
end

cv = gt.custom_visitor_new()
if options.random_alphanumeric then
  cv.id_provider = IDProviderRandomAlpha.new(args[2],
                                             options.random_length,
                                             options.seed)
else
  cv.id_provider = IDProviderNumeric.new(startval, increment, args[2], args[3])
end
cv.queue = {}
cv.suffixes = {}
cv.suffixes.gene = ''
function cv:visit_feature(fn)
  local transcript_counter = 0
  -- leave gaps/contigs untouched
  if fn:get_type() == "gap" or fn:get_type() == "contig" then
    return 0
  end
  local base_id = self.id_provider:next_for_feature(fn)
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
  return 0
end

vis_stream = gt.custom_stream_new_unsorted()
vis_stream.queue = cv.queue
vis_stream.instream = gt.gff3_in_stream_new_sorted(args[1])
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