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

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
local json = require ("dkjson")

function usage()
  io.stderr:write(string.format("Usage: %s <in.gff3> <in.protein.fasta> "
                                .. "<refdir> <speciesprefix> "
                                .. "<reference prefix>\n" , arg[0]))
  os.exit(1)
end

if #arg < 5 then
  usage()
end

ingff = arg[1]
inprots = arg[2]
refdir = arg[3]
--refgroup = arg[4]
speciesprefix = arg[4]
refprefix = arg[5]

-- load reference info
local reffile = io.open(refdir .. "/references.json", "rb")
if not reffile then
  error("invalid reference directory '" .. refdir .. "' -- missing references.json file")
end
local refcontent = reffile:read("*all")
reffile:close()
refs = json.decode(refcontent)

-- check for reference prefix
if not refs.species[refprefix] then
  error("invalid reference prefix '" .. refprefix .. "' not found in reference directory " .. refdir)
end

-- parse clusters
clindex, clusters = get_clusters(refdir .. "/_all/all_orthomcl.out")

-- global core clusters, for tree drawing
global_core_clusters = {}
-- group core clusters, for gene loss/gain comparison
group_core_clusters = {}
-- determine core
for clustername,cluster in pairs(clusters) do
  -- collect core clusters from ALL references
  local is_global_core = true
  for sp,_ in pairs(refs.species) do
    if is_global_core and not cluster.specidx[sp] then
      is_global_core = false
    end
  end
  if is_global_core then
    table.insert(global_core_clusters, cluster)
  end
  for refgroup,_ in pairs(refs.groups) do
    if not group_core_clusters[refgroup] then
      group_core_clusters[refgroup] = {}
    end
    -- collect group core clusters (those with all species in group)
    local is_group_core = true
    for _,sp in ipairs(refs.groups[refgroup]) do
      if is_group_core and not cluster.specidx[sp] then
        is_group_core = false
      end
    end
    if is_group_core then
      table.insert(group_core_clusters[refgroup], cluster)
    end
  end
end

cv = gt.custom_visitor_new()
singletons = {}
function cv:visit_feature(fn)
  for n in fn:get_children() do
    if n:get_type() == 'polypeptide' then
      local orths = n:get_attribute("orthologous_to")
      local dfrom = n:get_attribute("Derives_from")
      local found = false
      local groups = {}
      -- collect orthologs for this gene product
      if orths then
        for _,orth in ipairs(split(orths,",")) do
          if clindex[orth] then
            table.insert(groups, clindex[orth])
          end
        end
      end
      -- collect singletons
      if not n:get_attribute("ortholog_cluster") then
        singletons[dfrom] = true
      end
      if #groups > 0 then
        groups = table_unique(groups)
        -- only handle genes with all orthologs in the same cluster
        if #groups == 1 and dfrom then
          -- add this gene to cluster
          local cluster = groups[1]
          table.insert(cluster.members, {dfrom, speciesprefix})
          table.insert(cluster.species, speciesprefix)
          cluster.specidx[speciesprefix] = true
        end
      end
    end
  end
  return 0
end

-- visitor to output CC roots with members in its 'nididx' field to a text file
-- vor circos visualization
circos_visitor = gt.custom_visitor_new()
function circos_visitor:visit_feature(fn)
  for n in fn:get_children() do
    local nid = n:get_attribute("ID")
    if nid and self.nididx[nid] then
      self.io:write(fn:get_seqid() .. "  " .. fn:get_range():get_start()
                .. "  " .. fn:get_range():get_end() .. "  color=" .. self.color .. "\n")
    end
  end
end

-- perform core comparison
local vstream = visitor_stream_new(gt.gff3_in_stream_new_sorted(ingff), cv)
local gn = vstream:next_tree()
while (gn) do
  gn = vstream:next_tree()
end

-- extract core set for tree drawing
prots = {}
specseqs = {}
for s,v in pairs(refs.species) do
  local hdr, seq = get_fasta_nosep(refdir .. "/" .. s .. "/proteins.fasta")
  prots[s] = seq
  specseqs[s] = ""
end
_, nseq = get_fasta_nosep(inprots)
prots[speciesprefix] = nseq
specseqs[speciesprefix] = ""
i = 0
treegenes = io.open("tree_selection.genes", "w+")
treeclusters = io.open("tree_selection.clusters", "w+")
treeout = io.open("tree_selection.fasta", "w+")
for _,cl in ipairs(global_core_clusters) do
  local seen_species = {}
  if cl.specidx[speciesprefix] then
    for _,m in ipairs(cl.members) do
      local seq = prots[m[2]][m[1]]
      if not seen_species[m[2]] and seq then
        seen_species[m[2]] = true
        treegenes:write(m[1] .. "\t")
        specseqs[m[2]] = specseqs[m[2]] .. seq
      end
    end
    treegenes:write("\n")
    treeclusters:write(cl.name .. "\n")
    i = i + 1
    if i == 50 then
      break
    end
  end
end
if i > 0 then
  for s,seq in pairs(specseqs) do
    treeout:write(">" .. s .. "\n")
    print_max_width(seq, treeout, 80)
  end
end

-- searching for global core clusters with missing members in new species
missing = {}
local global_outfile = io.open("core_comparison.txt", "w+")
for _,cluster in ipairs(global_core_clusters) do
  if not cluster.specidx[speciesprefix] then
    for _, v in ipairs(cluster.members) do
      global_outfile:write("global\t" .. cluster.name .. "\t" .. v[1] ..
                           "\t" .. v[2] .. "\n")
      if refprefix == v[2] then
        missing[v[1]] = true
      end
    end
  end
end

-- searching for group core clusters with missing members in new species
for refgroup,members in pairs(refs.groups) do
  for _,cluster in ipairs(group_core_clusters[refgroup]) do
    if not cluster.specidx[speciesprefix] then
      for _, v in ipairs(cluster.members) do
        global_outfile:write(refgroup .. "\t" .. cluster.name .. "\t" .. v[1] ..
                             "\t" .. v[2] .. "\n")
      end
    end
  end
end

-- circos output
local stream = visitor_stream_new(gt.gff3_in_stream_new_sorted(refdir .. "/" .. refprefix .. "/annotation.gff3"), circos_visitor)
circos_visitor.nididx = missing
circos_visitor.color = 'green'
circos_visitor.io = io.open("core_comp_circos.txt", "w+")
local gn = stream:next_tree()
while (gn) do
  gn = stream:next_tree()
end


local stream = visitor_stream_new(gt.gff3_in_stream_new_sorted(ingff), circos_visitor)
circos_visitor.nididx = singletons
circos_visitor.color = 'black'
local gn = stream:next_tree()
while (gn) do
  gn = stream:next_tree()
end
