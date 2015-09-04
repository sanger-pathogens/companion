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
  io.stderr:write("Transfer functional annotations from orthologs in "
                  .. "reference genomes.\n")
  io.stderr:write(string.format("Usage: %s <target GFF3_file> "
                                .. "<source GFF file>\n", arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
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

-- this visitor collects the parent gene nodes of coding transcripts in a table
-- indexed by the mRNA ID so it can be queried using the Derives_from attribute
-- to find the gene a polypeptide derives from
get_parent_gene_vis = gt.custom_visitor_new()
function get_parent_gene_vis:visit_feature(fn)
  if fn:get_type() == "gene" then
    for fn2 in fn:children() do
      if fn2:get_type() == "mRNA" then
        if fn2:get_attribute("ID") then
          self.genes[fn2:get_attribute("ID")] = fn
        end
      end
    end
  end
end

annotate_vis = gt.custom_visitor_new()
function annotate_vis:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    local orths = fn:get_attribute("orthologous_to")
    local nof_orths = 0
    if orths then
      local gos = {}
      local names = {}
      local dbxrefs = {}
      local genenames = {}
      local synonyms = {}
      for _,orth in ipairs(split(orths, ",")) do
        nof_orths = nof_orths + 1
        orth = orth:gsub("GeneDB:","")     -- to be safe
        if self.store[orth] then
          local onode = self.store[orth]
          local orig_orth = orth
          -- product
          mem_name = onode:get_attribute("product")
          if mem_name then
            prod_a = gff3_extract_structure(mem_name)
            if prod_a[1] and prod_a[1].term then
              mem_name = prod_a[1].term
              -- remove pseudogene/fragment tags from product
              mem_name = mem_name:gsub(" ?%(pseudogene%)","")
              mem_name = mem_name:gsub(",? ?pseudogene","")
              mem_name = mem_name:gsub(" ?%(fragment%)","")
              mem_name = mem_name:gsub(",? ?fragment","")
              if not string.match(mem_name, "hypothetic") then
                table.insert(names, {mem_name, orig_orth})
              end
            end
          end
          -- Dbxref
          mem_dbxref = onode:get_attribute("Dbxref")
          if mem_dbxref then
            for _,dbxref in ipairs(split(mem_dbxref, ",")) do
              table.insert(dbxrefs, gff3_decode(dbxref))
            end
          end
        end
        if self.refgenestore[orth] then
          mem_genename = self.refgenestore[orth]:get_attribute("Name")
          if mem_genename then
            for _,genename in ipairs(split(mem_genename, ",")) do
              genename = split(genename, "%%3B")[1]
              table.insert(genenames, gff3_decode(genename))
            end
          end
          mem_synonym = self.refgenestore[orth]:get_attribute("synonym")
          if mem_synonym then
            for _,synonym in ipairs(split(mem_synonym, ",")) do
              synonym = split(synonym, "%%3B")[1]
              table.insert(synonyms, gff3_decode(synonym))
            end
          end
        end
      end
      -- transfer products
      if #names > 0 then
        local sorted_names = {}
        local dfrom = fn:get_attribute("Derives_from")
        -- group orthologs by name
        for _, name in ipairs(names) do
          if not sorted_names[name[1]] then
            sorted_names[name[1]] = {}
          end
          table.insert(sorted_names[name[1]], "GeneDB:"..name[2])
        end
        local full_names = {}
        local rank = 0
        local nof_names = 0
        -- count names
        for k,v in pairs(sorted_names) do
          nof_names = nof_names + 1
        end
        -- handle each name/orthologset pair
        for k,v in pairs(sorted_names) do
          local ref_gene = self.refgenestore[make_gene_name(dfrom)]
          if ref_gene and (ref_gene:get_attribute("Start_range")
              or ref_gene:get_attribute("End_range")) then
            str = "term=" .. k .. " (fragment);evidence=IEA;with=" .. table.concat(v, "|")
          else
            str = "term=" .. k .. ";evidence=IEA;with=" .. table.concat(v, "|")
          end
          if nof_names > 1 then
            if rank == 0 then
              str = str .. ";is_preferred=true"
              rank = rank + 1
            else
              str = str .. ";rank=" .. rank
            end
          end
          table.insert(full_names, gff3_encode(str))
        end
        fn:set_attribute("product", table.concat(full_names, ","))
      else
        str = "term=hypothetical protein"
        if nof_orths > 0 then
          str = str .. ", conserved"
        end
        fn:set_attribute("product", gff3_encode(str))
      end
      if #dbxrefs > 0 then
        -- XXX add evidence/WITH?
        local newdbxrefs = {}
        for _,v in ipairs(dbxrefs) do
          table.insert(newdbxrefs, gff3_encode(v))
        end
        fn:set_attribute("Dbxref", table.concat(table_unique(newdbxrefs), ","))
      end
      genenames = table_unique(genenames)
      -- for now, only handle n:1 orthology
      if #genenames == 1 then
        local dfrom = fn:get_attribute("Derives_from")
        if dfrom then
          local target_gene = self.destgenestore[dfrom]
          if target_gene then
            target_gene:set_attribute("Name", gff3_encode(genenames[1]))
            for i = 2,#genenames do
              table.insert(synonyms, genenames[i])
            end
          end
        end
        synonyms = table_unique(synonyms)
        if #synonyms > 0 then
          local dfrom = fn:get_attribute("Derives_from")
          if dfrom then
            local target_gene = self.destgenestore[dfrom]
            if target_gene then
              target_gene:set_attribute("synonym", table.concat(synonyms, ","))
            end
          end
        end
      end
    else
      fn:set_attribute("product", gff3_encode("term=hypothetical protein"))
    end
  end
  return 0
end

get_pp_vis = gt.custom_visitor_new()
function get_pp_vis:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    id = fn:get_attribute("ID")
    if not id then
      error("node " .. tostring(fn) .. " encountered with no ID in "
             .. fn:get_filename() .. ":" .. fn:get_line_number())
    end
    id = make_gene_name(id)
    self.store[make_gene_name(id)] = fn
  elseif fn:get_type() == "gene" then
    id = fn:get_attribute("ID")
    if not id then
      error("node " .. tostring(fn) .. " encountered with no ID in "
             .. fn:get_filename() .. ":" .. fn:get_line_number())
    end
    id = make_gene_name(id)
    self.genestore[make_gene_name(id)] = fn
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    for _,v in ipairs(self.visitors) do
      node:accept(v)
    end
  end
  return node
end

-- set up required indices
store = {}
refgenestore = {}
destgenestore = {}

-- get annotations from reference(s)
visitor_stream.visitors = {get_parent_gene_vis, get_pp_vis}
get_pp_vis.store = store
get_pp_vis.genestore = refgenestore
get_parent_gene_vis.genes = refgenestore
for i = 2,#arg do
  visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[i])
  local gn = visitor_stream:next_tree()
  while (gn) do
    gn = visitor_stream:next_tree()
  end
end

outnodes = {}
-- annotate target genome
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.visitors = {get_parent_gene_vis, annotate_vis}
get_parent_gene_vis.genes = destgenestore
annotate_vis.store = store
annotate_vis.refgenestore = refgenestore
annotate_vis.destgenestore = destgenestore
local gn = visitor_stream:next_tree()
while (gn) do
  table.insert(outnodes, gn)
  gn = visitor_stream:next_tree()
end

-- XXX do the retainids here
gff3vis = gt.gff3_visitor_new()
for _,n in ipairs(outnodes) do
  n:accept(gff3vis)
end
