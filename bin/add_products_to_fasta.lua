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
  io.stderr:write("Add preferred product terms to protein sequence headers.\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation> <protein seqs>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

function gff3_extract_structure(str)
  ret = {}
  for _,v in ipairs(split(str, ", ?")) do
    res = {}
    v = gff3_decode(v)
    for _,pair in ipairs(split(v, ";")) do
      if string.len(pair) > 0 then
        key, value = unpack(split(pair, "="))
        res[key] = value
      end
    end
    table.insert(ret, res)
  end
  return ret
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

function gff3_decode(s)
  if not s then
    return s
  end
  s = string.gsub(s, "%%([0-9a-fA-F][1-9a-fA-F])", function (n)
        return string.char(tonumber("0x" .. n))
      end)
  return s
end

function get_fasta(filename, sep)
  local keys = {}
  local seqs = {}
  local cur_hdr = nil
  local cur_seqs = {}
  for l in io.lines(filename) do
    hdr = l:match(">(.*)")
    if hdr then
      table.insert(keys, hdr)
      if #cur_seqs > 0 and cur_hdr then
        if not seqs[cur_hdr] then
          seqs[cur_hdr] = table.concat(cur_seqs, sep)
        end
      end
      cur_hdr = hdr
      cur_seqs = {}
    else
      table.insert(cur_seqs, l)
    end
  end
  if cur_hdr and not seqs[cur_hdr] then
    seqs[cur_hdr] = table.concat(cur_seqs, sep)
  end
  return keys, seqs
end

function get_fasta_nosep(filename)
  return get_fasta(filename, "")
end

function print_max_width(str, ioo, width)
  local i = 1
  while (i + width - 1) < str:len() do
    ioo:write(str:sub(i, i + width - 1))
    ioo:write("\n")
    i = i + width
  end
  ioo:write(str:sub(i))
  ioo:write("\n")
end

products = {}

visitor = gt.custom_visitor_new()
visitor.last_seqid = nil
function visitor:visit_feature(fn)
  for node in fn:children() do
    if node:get_type() == "polypeptide" then
      local dfrom = node:get_attribute("Derives_from")
      local product = node:get_attribute("product")
      if (dfrom and product) then
        local p_data = gff3_extract_structure(product)
        local mterm = p_data[1].term
        for _,v in ipairs(p_data) do
          if v.is_preferred then
            mterm = v.term
          end
        end
        products[dfrom] = mterm
      end
    end
  end
  return 0
end

visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = visitor
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

keys, seqs = get_fasta_nosep(arg[2])
for i,k in ipairs(keys) do
  local rid = split(k, "%s+")[1]
  if products[rid] then
    print(">".. rid .. " " .. products[rid])
  else
    print(">".. rid)
  end
  print_max_width(seqs[k], io.stdout, 60)
end
