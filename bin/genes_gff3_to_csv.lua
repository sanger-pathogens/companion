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
  io.stderr:write(string.format("Usage: %s <GFF>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

peps = {}
gather_v = gt.custom_visitor_new()
function gather_v:visit_feature(fn)
  for node in fn:children() do
    if node:get_type() == "polypeptide" then
      local dfrom = node:get_attribute("Derives_from")
      if dfrom then
        peps[dfrom] = node
      end
    end
  end
  return 0
end

output_v = gt.custom_visitor_new()
function output_v:visit_feature(fn)
  for node in fn:children() do
    local gprod
    local gtype
    if node:get_type() == "mRNA" then
      local id = node:get_attribute("ID")
      if id and peps[id] then
        pr = peps[id]:get_attribute("product")
        if pr then
          pr_a = gff3_extract_structure(peps[id]:get_attribute("product"))
          gtype = "coding"
          gprod = pr_a[1].term
        end
      end
    elseif node:get_type() == "pseudogenic_transcript" then
      local id = node:get_attribute("ID")
      if id and peps[id] then
        pr = peps[id]:get_attribute("product")
        if pr then
          pr_a = gff3_extract_structure(peps[id]:get_attribute("product"))
          gtype = "pseudogene"
          gprod = pr_a[1].term
        end
      end
    elseif node:get_type():match("RNA$") then
      gtype = "non_coding"
      gprod = node:get_type()
    end
    if gtype and gprod then
      print(fn:get_attribute("ID") .. "\t" .. gtype .. "\t" .. gprod .. "\t"
               .. fn:get_seqid() .. "\t" .. fn:get_range():get_start()
                .. "\t" .. fn:get_range():get_end() .. "\t" .. fn:get_strand())
    end
  end
  return 0
end

-- setup generic visitor stream
vis_stream = gt.custom_stream_new_unsorted()
function vis_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.v)
  end
  return node
end

vis_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vis_stream.v = gather_v
local gn = vis_stream:next_tree()
while (gn) do
  gn = vis_stream:next_tree()
end

vis_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vis_stream.v = output_v
local gn = vis_stream:next_tree()
while (gn) do
  gn = vis_stream:next_tree()
end