#!/usr/bin/env gt

--[[
  Copyright (c) 2015 Sascha Steinbiss <ss34@sanger.ac.uk>
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
                                .. "<refdir> <speciesprefix>\n" , arg[0]))
  os.exit(1)
end

if #arg < 4 then
  usage()
end

ingff = arg[1]
inprots = arg[2]
refdir = arg[3]
--refgroup = arg[4]
speciesprefix = arg[4]

-- load reference info
local reffile = io.open(refdir .. "/references.json", "rb")
local refcontent = reffile:read("*all")
reffile:close()
refs = json.decode(refcontent)

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
for _,cl in ipairs(global_core_clusters) do
  local seen_species = {}
  for _,m in ipairs(cl.members) do
    if not seen_species[m[2]] then
      seen_species[m[2]] = true
      local seq = prots[m[2]][m[1]]
      treegenes:write(m[1] .. "\t")
      specseqs[m[2]] = specseqs[m[2]] .. seq
    end
  end
  treegenes:write("\n")
  i = i + 1
  if i == 50 then
    break
  end
end
treeout = io.open("tree_selection.fasta", "w+")
for s,seq in pairs(specseqs) do
  treeout:write(">" .. s .. "\n")
  print_max_width(seq, treeout, 80)
end

-- searching for global core clusters with missing members in new species
missing = {}
for _,cluster in ipairs(global_core_clusters) do
  if not cluster.specidx[speciesprefix] then
    for _, v in ipairs(cluster.members) do
      print("global\t" .. cluster.name .. "\t" .. v[1] .. "\t" .. v[2])
    end
  end
end

-- searching for group core clusters with missing members in new species
for refgroup,members in pairs(refs.groups) do
  print("missing ".. refgroup .. " core clusters (" .. #group_core_clusters[refgroup] .. " total)")
  for _,cluster in ipairs(group_core_clusters[refgroup]) do
    if not cluster.specidx[speciesprefix] then
      for _, v in ipairs(cluster.members) do
        print(refgroup .. "\t" .. cluster.name .. "\t" .. v[1] .. "\t" .. v[2])
      end
    end
  end
end
