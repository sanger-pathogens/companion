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
  io.stderr:write("Transforms a GFF3 file for loading into TriTrypDB.\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

function gff3_encode(s)
  s = string.gsub(s, "[\t\n\r;=%&,]", function (c)
        return string.format("%%%02X", string.byte(c))
      end)
  return s
end

function gff3_decode(s)
  s = string.gsub(s, "%%([0-9a-fA-F][1-9a-fA-F])", function (n)
        return string.char(tonumber("0x" .. n))
      end)
  return s
end

function gff3_extract_structure(str)
  local ret = {}
  for _,v in ipairs(split(str, ", ?")) do
    local res = {}
    local v = gff3_decode(v)
    for _,pair in ipairs(split(v, ";")) do
      if string.len(pair) > 0 then
        local key, value = unpack(split(pair, "="))
        res[key] = value
      end
    end
    table.insert(ret, res)
  end
  return ret
end

function gff3_explode(tab)
  local ret = {}
  for _,item in ipairs(tab) do
    local _tmp = {}
    for k,v in pairs(item) do
      table.insert(_tmp, k .. "=" .. v)
    end
    table.insert(ret, gff3_encode(table.concat(_tmp, ";")))
  end
  return table.concat(ret, ",")
end

function table_count(tt, value)
  local count = 0
  for ii,xx in pairs(tt) do
    if xx == value then count = count + 1 end
  end
  return count
end

function table_unique(tt)
  local newtable = {}
  for ii,xx in ipairs(tt) do
    if table_count(newtable, xx) == 0 then
      newtable[#newtable+1] = xx
    end
  end
  return newtable
end

function table_add(t1,t2)
    for i=1,#t2 do
        t1[#t1+1] = t2[i]
    end
    return t1
end

function split(str, pat)
   local t = {}
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
        table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end

partial_visitor = gt.custom_visitor_new()
function partial_visitor:visit_feature(fn)
  if fn:get_type() == 'gene' then
    if fn:get_attribute("Start_range") then
      if fn:get_strand() == "-" then
        fn:set_attribute("threeEndPartial", "1")
      else
        fn:set_attribute("fiveEndPartial", "1")
      end
    end
    if fn:get_attribute("End_range") then
      if fn:get_strand() == "+" then
        fn:set_attribute("threeEndPartial", "1")
      else
        fn:set_attribute("fiveEndPartial", "1")
      end
    end
  end
  return 0
end

seleno_tag_visitor = gt.custom_visitor_new()
function seleno_tag_visitor:visit_feature(fn)
  if fn:get_type() == 'polypeptide' then
    if fn:get_attribute("product")
      and string.match(fn:get_attribute("product"), "selenoprotein") then
      -- we have a selenoprotein
      if not fn:get_attribute("stop_codon_redefined_as_selenocysteine") then
        fn:add_attribute("stop_codon_redefined_as_selenocysteine", "1")
      end
    end
  end
  return 0
end

ncrna_visitor = gt.custom_visitor_new()
function ncrna_visitor:visit_feature(fn)
  for fn2 in fn:get_children() do
    if fn2:get_type() == 'ncRNA' or
         fn2:get_type() == 'snRNA' or
         fn2:get_type() == 'snoRNA' or
         fn2:get_type() == 'rRNA' or
         fn2:get_type() == 'lncRNA' or
         fn2:get_type() == 'snRNA' or
         fn2:get_type() == 'scRNA' or
         fn2:get_type() == 'tRNA' then
      for n in fn2:get_children() do
        if n:get_type() == "CDS" or n:get_type():match("UTR") then
          if n:get_attribute("Dbxref") and not fn2:get_attribute("Dbxref") then
            -- copy over Dbxref from CDS to ncRNA
            fn2:add_attribute("Dbxref", n:get_attribute("Dbxref"))
          end
          if n:get_attribute("product") and not fn2:get_attribute("product") then
            -- copy over product from CDS to ncRNA
            fn2:add_attribute("product", n:get_attribute("product"))
          end
          -- remove CDS from ncRNA
          fn:remove_leaf(n)
        end
      end
      if fn2:get_attribute("score") then
        fn2:remove_attribute("score")
      end
      if fn2:get_attribute("model_range") then
        fn2:remove_attribute("model_range")
      end
      if fn2:get_attribute("model_name") then
        fn2:remove_attribute("model_name")
      end
      if fn2:get_attribute("evalue") then
        fn2:remove_attribute("evalue")
      end
      if fn2:get_attribute("anticodon") then
        fn2:remove_attribute("anticodon")
      end
      if fn2:get_attribute("gc_content") then
        fn2:remove_attribute("gc_content")
      end
      if fn2:get_attribute("gc") then
        fn2:remove_attribute("gc")
      end
      if fn2:get_attribute("aa") then
        fn2:remove_attribute("aa")
      end
    end
  end
  return 0
end

quotes_visitor = gt.custom_visitor_new()
function quotes_visitor:visit_feature(fn)
  for fn2 in fn:get_children() do
    for k,v in fn2:attribute_pairs() do
      if v:match('"') then
        v = v:gsub('"','')
        fn2:set_attribute(k, v)
      end
    end
  end
  return 0
end

uc_visitor = gt.custom_visitor_new()
function uc_visitor:visit_feature(fn)
  ucs = {}
  if not fn:get_attribute("tritryp_uc") then
    for n in fn:get_children() do
      hist = n:get_attribute("history")
      if hist then
        for m in hist:gmatch('[Tt][Rr][Ii][Tt][Rr][Yy][Pp]_[Uu][Cc]:"?([^;%%"]+)') do
          table.insert(ucs, m)
        end
        for m in hist:gmatch('[Tt][Rr][Ii][Tt][Rr][Yy][Pp]_[Uu][Cc]%%3."?([^;%%"]+)') do
          table.insert(ucs, m)
        end
        fn:remove_attribute("history")
      end
    end
    ucs = table_unique(ucs)
    if #ucs > 0 then
      fn:add_attribute("tritryp_uc", table.concat(ucs, ","))
    end
  end
  return 0
end


attrs_to_remove = {"translation",
                   "controlled_curation",
                   "gPI_anchor_cleavage_site",
                   "orthologous_to",
                   "ortholog_cluster",
                   "comment",
                   "membrane_structure",
                   "polypeptide_domain",
                   "cytoplasmic_polypeptide_region",
                   "non_cytoplasmic_polypeptide_region",
                   "transmembrane_polypeptide_region",
                   "signal_peptide",
                   "paralogous_to",
                   "timelastmodified",
                   "full_GO",
                   "colour",
                   "ratt_ortholog",
                   "feature_id",
                   "isObsolete"}
remove_attrs_visitor = gt.custom_visitor_new()
function remove_attrs_visitor:visit_feature(fn)
  for n in fn:get_children() do
    for _,v in ipairs(attrs_to_remove) do
      if n:get_attribute(v) then
        n:remove_attribute(v)
      end
    end
  end
  return 0
end

polypeptide_child_trimmer_visitor = gt.custom_visitor_new()
function polypeptide_child_trimmer_visitor:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    for fn2 in fn:direct_children() do
      if fn2:get_type() == "membrane_structure" then
        for fn3 in fn2:direct_children() do
          fn:remove_leaf(fn3)
        end
      end
      fn:remove_leaf(fn2)
    end
  end
  return 0
end

exon_remover_visitor = gt.custom_visitor_new()
function exon_remover_visitor:visit_feature(fn)
  if fn:get_type() == "gene" then
    for fn2 in fn:children() do
      if fn2:get_type() == "exon" then
        fn:remove_leaf(fn2)
      end
    end
  end
  return 0
end

-- set up streams

fixup_stream = gt.custom_stream_new_unsorted()
fixup_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
function fixup_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(partial_visitor)
    node:accept(seleno_tag_visitor)
    node:accept(ncrna_visitor)
    node:accept(quotes_visitor)
    node:accept(remove_attrs_visitor)
    node:accept(polypeptide_child_trimmer_visitor)
    node:accept(exon_remover_visitor)
  end
  return node
end

-- run streams

out_stream = gt.gff3_out_stream_new(fixup_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end
