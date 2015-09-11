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
  io.stderr:write("Output GO annotation in GAF format from annotate_eukaryote "
                  .. "GFF3 output.\n")
  io.stderr:write(string.format("Usage: %s <GFF3 with polypeptide annotations> "
                                .. "<InterproScan output GFF3> <db> <taxonid> "
                                .. "<GO OBO>\n" , arg[0]))
  os.exit(1)
end

if #arg < 5 then
  usage()
end

FILTERED_ATTRIBS = {date = true, Target = true, status = true}
FILTERED_SOURCES = {TMHMM = true, Phobius = true}

db = arg[3]
taxonid = arg[4]
go_obo = arg[5]

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function parse_obo(filename)
  local gos = {}
  local this_id = {}
  local this_name = nil
  for l in io.lines(filename) do
    if string.match(l, "%[Term%]") then
      if #this_id and this_name then
        for _, n in ipairs(this_id) do
          gos[n] = {name=this_name, aspect=this_aspect}
        end
      end
      this_id = {}
      this_term = nil
      this_aspect = nil
    end
    m = string.match(l, "id: (.+)")
    if m then
      table.insert(this_id, m)
    end
    m = string.match(l, "^name: (.+)")
    if m then
      this_name = m
    end
    m = string.match(l, "^namespace: (.+)")
    if m then
      if m == "biological_process" then
        this_aspect = "P"
      elseif m == "molecular_function" then
        this_aspect = "F"
      elseif m == "cellular_component" then
        this_aspect = "C"
      end
    end
  end
  return gos
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
          if mygo ~= "GO:0005515" then  -- only wih IPI!
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
      -- get SMART annots first
      if hits["SMART"] then
        extract_gos(hits["SMART"], "SMART", gos)
      end
      -- if we have Pfam hits, they are preferred...
      if hits["Pfam"] then
        extract_gos(hits["Pfam"], "Pfam", gos)
      end
      -- and TIGRFAM will override everything
      if hits["TIGRFAM"] then
        extract_gos(hits["TIGRFAM"], "JCVI_TIGRFAMS", gos)
      end
      full_gos = {}
      prod = gff3_extract_structure(fn:get_attribute("product"))
      for k,v in pairs(gos) do
        if self.obos[k] then
          aspect = self.obos[k].aspect
          print(table.concat({tostring(db),
                              fn:get_attribute("ID"):split("%.")[1],
                              fn:get_attribute("ID"):split("%.")[1],
                              "",
                              k,
                              "GO_REF:0000002", -- transferred from Interpro
                              "IEA",
                              table.concat(v, "|"),
                              aspect,
                              tostring(prod[1].term),
                              "",
                              "gene",
                              "taxon:" .. tostring(taxonid),
                              os.date("%Y%m%d"),
                              tostring(db)
                             }, "\t"))
        end
      end
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
    if node:get_type() ~= "polypeptide" then
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

obos = parse_obo(go_obo)

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
annotate_vis.obos = obos
annotate_vis.hits = get_ipro_vis.hits
out_stream =  visitor_stream -- gt.gff3_out_stream_new_retainids(visitor_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
