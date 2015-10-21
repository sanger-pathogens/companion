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

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("SimpleChainer")
require("optparse")

op = OptionParser:new({usage="%prog <options> < merged.gff3",
                       oneliner="Selects the 'best' gene models "
                         .. "from a pooled set of GFF3 annotations.",
                       version="0.1"})
op:option{"-w", action='store', dest='weight_func',
                help="Lua script defining the weight function 'get_weight()'"}
op:option{"-s", action='store', dest='sequence',
                help="sequence file for the given annotation"}
options,args = op:parse({weight_func=nil, sequence=nil})

function usage()
  op:help()
  os.exit(1)
end

if not options.sequence then
  usage()
end
regmap = gt.region_mapping_new_seqfile_matchdescstart(options.sequence)

-- default weight function: gene length
-- this should most of the time be overridden by a more specific,
-- evidence aware function
function _get_weight(gene)
  return gene:get_range():length()
end

-- load weight function, if given
if not options.weight_func then
  get_weight = _get_weight
else
  dofile(options.weight_func)
end

if not get_weight then
  error("no weight function ('get_weight()') found in file "
           .. options.weight_func)
  os.exit(1)
end

stream = gt.custom_stream_new_unsorted()
stream.outqueue = {}
stream.curr_gene_set = {}
stream.curr_rng = nil
stream.last_seqid = nil
function stream:process_current_cluster()
  local bestset = nil
  local max = 0

  -- keep only non-overlapping chain with highest weight
  bestset = SimpleChainer.new(self.curr_gene_set, get_weight, regmap):chain()

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

stream.instream = gt.gff3_in_stream_new_sorted()
stream.idx = feature_index

out_stream = gt.gff3_out_stream_new(stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end

