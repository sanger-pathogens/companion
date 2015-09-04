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

function usage()
  io.stderr:write(string.format("Usage: %s <pfam2go file>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function get_pfam2go(file)
  ret = {}
  for line in io.lines(file) do
    if string.sub(line, 1, 1) ~= "!" then
      local la = line:split_ws()
      local lb = line:split_sep(';')
      pfamid = (la[1]:split_sep(':'))[2]
      if not ret[pfamid] then
        ret[pfamid] = {}
      end
      local goid = (lb[2]:gsub("%s+",""))
      table.insert(ret[pfamid], goid)
    end
  end
  return ret
end

function string:split_sep(sep)
  local sep, fields = sep or ":", {}
  local pattern = string.format("([^%s]+)", sep)
  self:gsub(pattern, function(c) fields[#fields+1] = c end)
  return fields
end

function string:split_ws()
  local fields = {}
  self:gsub("([^%s]+)", function(c) fields[#fields+1] = c end)
  return fields
end

pfam2go = get_pfam2go(arg[1])
gff3v = gt.gff3_visitor_new()

while true do
  local line = io.read()
  if line == nil then break end

  if string.sub(line, 1, 1) ~= "#" then
    local la = line:split_ws()
    local targetname = la[1]
    local targetacc = la[2]
    local qryname = la[4]
    local qryacc = la[5]
    local evalue = tonumber(la[7])
    local score = la[8]
    local hmmfrom = tonumber(la[16])
    local hmmto = tonumber(la[17])
    local alifrom = tonumber(la[18])
    local alito = tonumber(la[19])
    -- stitch the description back together again
    local desc = {}
    for i = 23,#la do
      table.insert(desc, la[i])
    end
    local targetdesc = table.concat(desc, " ")
    local gostr = ""
    local pfid = targetacc:split_sep('.')[1]
    if pfam2go[pfid] then
      gostr = table.concat(pfam2go[pfid], ",")
    end

    fn = gt.feature_node_new(qryname, "protein_match", alifrom, alito, "+")
    if targetname and string.len(targetname) > 0 then
      fn:add_attribute("Name", pfid)
    end
    if targetdesc and string.len(targetdesc) > 0 then
      fn:add_attribute("signature_desc", gff3_encode(targetdesc))
    end
    if gostr and string.len(gostr) > 0 then
      fn:add_attribute("Ontology_term", gostr)
    end
    fn:set_source("Pfam")
    fn:add_attribute("Target", targetname .. " " .. hmmfrom .. " " .. hmmto)
    fn:set_score(evalue)
    fn:accept(gff3v)
  end
end