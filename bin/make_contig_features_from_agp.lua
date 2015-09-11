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
  io.stderr:write(string.format("Usage: %s <AGP> [gap_type]\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

agp = arg[1]
gap_type = arg[2]

gff3v = gt.gff3_visitor_new()

for l in io.lines(agp) do
  if string.sub(l, 1, 1) ~= '#' then
    obj, obj_s, obj_e, part_n, type, c6, c7, c8, c9 = unpack(split(l, "%s+"))
    -- obj = obj:gsub("%.%d+$","") -- no longer required!?
    if type == 'N' then
      fn = gt.feature_node_new(obj, "gap", obj_s, obj_e, ".")
      --fn:add_attribute("gap_type", type)
      fn:add_attribute("estimated_length", obj_e-obj_s+1)
      if gap_type then
        fn:add_attribute("gap_type", gap_type)
      end
      fn:accept(gff3v)
    elseif type == 'U' then
      fn = gt.feature_node_new(obj, "gap", obj_s, obj_e, ".")
      --fn:add_attribute("gap_type", type)
      if gap_type then
        fn:add_attribute("gap_type", gap_type)
      end
      fn:accept(gff3v)
    else
      fn = gt.feature_node_new(obj, "contig", obj_s, obj_e, c9)
      fn:add_attribute("Name", c6)
      fn:accept(gff3v)
    end
  end
end