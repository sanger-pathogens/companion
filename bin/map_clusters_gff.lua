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
  io.stderr:write("Annotates GFF3 files with OrthoMCL results.\n")
  io.stderr:write(string.format("Usage: %s <GFF3_file> <orthomcl output> "
                                .."<map_file>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function make_gene_name(str)
  if string.match(str, "%.%.") then
    return split(str, "%.%.")[1]
  elseif string.match(str, ":") then
    return split(str, ":")[1]
  else
    return str
  end
end

orth_vis = gt.custom_visitor_new()
function orth_vis:visit_feature(fn)
  -- TODO: also consider pseudogenes, also make use of Derives_from to link
  -- mRNAs and polypeptides
  if fn:get_type() == "polypeptide" then
    id = fn:get_attribute("ID")
    if not id then
      error("node " .. tostring(fn) .. " encountered with no ID in "
             .. fn:get_filename() .. ":" .. fn:get_line_number())
    end
    id = make_gene_name(id)
    if self.clindex[id] then
      orthos = {}
      for _, member in ipairs(self.clindex[id].members) do
        if member ~= id then
          table.insert(orthos, member)
        end
      end
      if #orthos > 0 then
        fn:set_attribute("orthologous_to", table.concat(orthos, ","))
        fn:set_attribute("ortholog_cluster", self.clindex[id].name)
      end
    end
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

-- parse mapfile and build mappings to original genes, if given
if arg[3] then
  map_old_new = {}
  map_new_old = {}
  for l in io.lines(arg[3]) do
    k, v = unpack(split(l, "\t"))
    if not k or not v then
      error("could not parse mapfile line: " .. l)
    end
    map_new_old[k] = v
    map_old_new[v] = k
  end
end

-- build cluster index
clindex = {}
for l in io.lines(arg[2]) do
  local name, nofgenes, members = string.match(l,
                                  "^(ORTHOMCL%d+)%((%d+) genes[^)]+%):%s+(.+)")
  if not name or not members then
    error("could not parse cluster or members from line " .. l)
  end
  local n = 0
  thisclust = {name = name, members = {}}
  for mname, mspecies in members:gmatch("([^(]+)%(([^)]+)%)%s?") do
    if map_new_old then
      if not map_new_old[mname] then
        error("missing mapping for member " .. mname .. " from cluster " .. name)
      end
      oldname = make_gene_name(map_new_old[mname])
      table.insert(thisclust.members, oldname)
      clindex[oldname] = thisclust
    else
      genename = make_gene_name(mname)
      table.insert(thisclust.members, genename)
      clindex[genename] = thisclust
    end
    n = n + 1
  end
  if n ~= tonumber(nofgenes) then
    error("parsed " .. n .. " members from cluster " .. name
            .. ", but expected " .. nofgenes)
  end
end

-- annotate nodes with clusters
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = orth_vis
orth_vis.clindex = clindex

out_stream = gt.gff3_out_stream_new_retainids(visitor_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
