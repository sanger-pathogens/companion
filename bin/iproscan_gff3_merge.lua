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

function usage()
  io.stderr:write("Infer and annotate GFFs with functional data from " ..
                  "annotate_eukaryote GFF3 output.\n")
  io.stderr:write(string.format("Usage: %s <GFF3 with polypeptide annotations> "
                                .. "<InterproScan output GFF3>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
end

FILTERED_ATTRIBS = {date = true, Target = true, status = true}
FILTERED_SOURCES = {TMHMM = true, Phobius = true}

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

-- adjusts product string to make more sense in a protein context
function make_product_name(str)
  words = split(str, "%s+")
  lastword = words[#words]
  if lastword == 'domain' or lastword == 'motif' or lastword == 'region'
      or string.match(lastword, "finger.?$") then
    str = str .. " containing protein, putative"
  else
    str = str .. ", putative"
  end
  return str
end

-- returns array of hit descriptions and hit IDs as evidence
function extract_hit_strings(hits, dbtag)
  local with = {}
  local descs = {}
  local dbxrefs = {}
  for _,node in ipairs(hits) do
    local name = gff3_decode(node:get_attribute("Name"))
    local desc = gff3_decode(node:get_attribute("signature_desc"))
    local dbxref = gff3_decode(node:get_attribute("Dbxref"))
    if desc then
      if name then
        table.insert(with, dbtag .. ":" .. name:gsub("\"",""))
      end
      if string.match(desc, ": ") then
        table.insert(descs, split(desc, ": ")[2])
      else
        table.insert(descs, desc)
      end
    end
    if dbxref then
      for _, dbx in ipairs(split(dbxref:gsub("\"",""), ",")) do
        table.insert(dbxrefs, dbx)
      end
    end
  end
  with = table_unique(with)
  descs = table_unique(descs)
  dbxrefs = table_unique(dbxrefs)
  if #descs > 0 then
    prod = "term=" .. make_product_name(table.concat(descs,"/"))
             .. ";evidence=IEA;with=" .. table.concat(with, "|")
  end
  return prod, descs, with, dbxrefs
end

-- extracts GO IDs with evidence
function extract_gos(hits, dbtag, gos)
  if hits then
    for _,hit in ipairs(hits) do
      local ot = hit:get_attribute("Ontology_term")
      local mygos = nil
      if ot then
        mygos = table_unique(split(ot:gsub("\"", ""), ","))
      end
      if mygos then
        local with = hit:get_attribute("Name")
        if with then
          with = dbtag .. ":" .. with:gsub("\"","")
        end
        for _,mygo in ipairs(mygos) do
          if not gos[mygo] then
            gos[mygo] = {}
          end
          table.insert(gos[mygo], with)
          gos[mygo] = table_unique(gos[mygo])
        end
      end
    end
  end
end

-- converts amino acid coords relative to <parent> to DNA coordinates
function aminoloc_to_dnaloc(parent, range, strand)
  par_rng = parent:get_range()
  ret_range = nil
  if strand ~= '+' then -- iproscan GFF feats should always be on the + strand
    error("Strand is #{strand} but should be +")
  end
  if parent:get_strand() == '-' then
    ret_range = {par_rng:get_end() - (range:get_end()-1)*3+1,
                 par_rng:get_end() - (range:get_start()-1)*3}
  else
    ret_range = {(par_rng:get_start() + (range:get_start()-1)*3),
                 (par_rng:get_start() + (range:get_end()-1)*3)-1}
  end
  return ret_range
end

annotate_vis = gt.custom_visitor_new()
function annotate_vis:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    local ppid = fn:get_attribute("Derives_from")
    if not ppid then
      error("polypeptide node (" .. fn.get_filename() .. ":"
              .. fn:get_line_number() .. ") has no ID attribute")
    end
    local hits = self.hits[ppid]
    local gos = {}
    if hits then
      -- copy over Interpro features to polypeptide
      for k,v in pairs(hits) do
        if not FILTERED_SOURCES[k] then
          for _,n in ipairs(v) do
            local rng = aminoloc_to_dnaloc(fn, n:get_range(), n:get_strand())
            if rng[1] <= rng[2] then
              if fn:get_range():contains(gt.range_new(rng[1],rng[2])) then
                local new_node = gt.feature_node_new(fn:get_seqid(),
                                                     "protein_match",
                                                     rng[1], rng[2],
                                                     fn:get_strand())
                new_node:set_source(k)
                for attr, attrv in n:attribute_pairs() do
                  if not FILTERED_ATTRIBS[attr] then
                    new_node:set_attribute(attr, string.gsub(attrv, "\"",""))
                  end
                end
                fn:add_child(new_node)
              else
                io.stderr:write("coordinates for feature outside of parent: "
                                 .. tostring(fn:get_range()) .. " vs. "
                                 .. tostring(rng) .. " -- not attaching to "
                                 .. "polypeptide parent")
              end
            else
              io.stderr:write("warning: invalid range " .. rng[1] .. " > "
                                 .. rng[2] .. "\n")
            end
          end
        end
      end
      local had_product = (fn:get_attribute("product"):match("hypothetical") == nil)
      local had_dbxref = fn:get_attribute("Dbxref")
      -- get SMART annots first
      if hits["SMART"] then
        prod, descs, with, dbxrefs = extract_hit_strings(hits["SMART"], "SMART")
        if #descs > 0 and not had_product then
          -- they sometimes have a trailing period
          fn:set_attribute("product", gff3_encode(prod:gsub("%.,",",")))
        end
        if #dbxrefs > 0 and not had_dbxref then
          local escaped_dbxrefs = {}
          for _,v in ipairs(dbxrefs) do
            table.insert(escaped_dbxrefs, gff3_encode(v))
          end
          fn:set_attribute("Dbxref", table.concat(escaped_dbxrefs, ","))
        end
        extract_gos(hits["SMART"], "SMART", gos)
      end
      -- if we have Pfam hits, they are preferred...
      if hits["Pfam"] then
        prod, descs, with, dbxrefs = extract_hit_strings(hits["Pfam"], "Pfam")
        if #descs > 0 and not had_product then
          fn:set_attribute("product", gff3_encode(prod))
        end
        if #dbxrefs > 0 and not had_dbxref then
          local escaped_dbxrefs = {}
          for _,v in ipairs(dbxrefs) do
            table.insert(escaped_dbxrefs, gff3_encode(v))
          end
          fn:set_attribute("Dbxref", table.concat(escaped_dbxrefs, ","))
        end
        extract_gos(hits["Pfam"], "Pfam", gos)
      end
      -- and TIGRFAM will override everything
      if hits["TIGRFAM"] then
        prod, descs, with, dbxrefs = extract_hit_strings(hits["TIGRFAM"],
                                                         "JCVI_TIGRFAMS")
        if #descs > 0 and not had_product then
          fn:set_attribute("product", gff3_encode(prod))
        end
        if #dbxrefs > 0 and not had_dbxref then
          local escaped_dbxrefs = {}
          for _,v in ipairs(dbxrefs) do
            table.insert(escaped_dbxrefs, gff3_encode(v))
          end
          fn:set_attribute("Dbxref", table.concat(escaped_dbxrefs, ","))
        end
        extract_gos(hits["TIGRFAM"], "JCVI_TIGRFAMS", gos)
      end
      oterms = {}
      full_gos = {}
      for k,v in pairs(gos) do
        table.insert(oterms, gff3_encode(k))
       -- table.insert(full_gos, gff3_encode("GOid=" .. k
       --                     .. ";evidence=IEA;with=" .. table.concat(v, "|")))
      end
      if #oterms > 0 and not fn:get_attribute("Ontology_term") then
        fn:set_attribute("Ontology_term", table.concat(oterms, ","))
      end
      -- if #full_gos > 0 then
      --   fn:set_attribute("full_GO", table.concat(full_gos, ","))
      -- end
    end
  end
  return 0
end

-- visitor to collect interproscan results into indexed storage structure
get_ipro_vis = gt.custom_visitor_new()
get_ipro_vis.hits = {}
get_ipro_vis.counts = {}
function get_ipro_vis:visit_feature(fn)
  for node in fn:get_children() do
    if node:get_type() == "protein_match" then
      local nid = node:get_seqid()
      local src = node:get_source()
      if not self.hits[nid] then
        self.hits[nid] = {}
      end
      if not self.hits[nid][src] then
        self.hits[nid][src] = {}
      end
      if not self.counts[nid] then
        self.counts[nid] = {}
      end
      if not self.counts[nid][src] then
        self.counts[nid][src] = 0
      end
      table.insert(self.hits[nid][src], node)
      self.counts[nid][src] = self.counts[nid][src] + 1
    end
  end
  return 0
end

visitor_stream = visitor_stream_new(gt.gff3_in_stream_new_sorted(arg[2]),
                                    get_ipro_vis)
-- get and index iproscan results
store = {}
get_ipro_vis.store = store
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

visitor_stream = visitor_stream_new(gt.gff3_in_stream_new_sorted(arg[1]),
                                    annotate_vis)
-- use the results to annotate an input GFF file
annotate_vis.store = store
annotate_vis.hits = get_ipro_vis.hits
out_stream = gt.gff3_out_stream_new_retainids(visitor_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
