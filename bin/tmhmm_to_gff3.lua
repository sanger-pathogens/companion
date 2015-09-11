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
  io.stderr:write(string.format("Usage: %s <TMHMM2 output> "
                                .." <GFF with polypeptide annotations>\n" ,
                                arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

cv = gt.custom_visitor_new()
cv.tmhs = {}
function cv:visit_feature(fn)
  local pp_id = nil
  if fn:get_type() == "polypeptide" then
    pp_id = fn:get_attribute("ID")
    local tm = nil
    -- search for a TMH in the TMHMM results
    for seq,v in pairs(cv.tmhs) do
      if string.match(pp_id, seq) then
        tm = v
        break
      end
    end
    -- if we have a TMH for this feature then attach it
    if tm and tm.number > 0 then
      node_rng = fn:get_range()
      whole_region_rng = {tonumber(tm.parts[1][2])-1, tonumber(tm.parts[#tm.parts][3])-1}
      -- print(node_rng:length() - 3)  -- minus stop codon
      -- print((whole_region_rng[2]-whole_region_rng[1]+1) * 3)
      if fn:get_strand() == '-' then
        whole_region_rng = {node_rng:get_end() - (whole_region_rng[2])*3+1,
                            node_rng:get_end() - (whole_region_rng[1])*3}
      else
        whole_region_rng = {node_rng:get_start() + (whole_region_rng[1])*3,
                            node_rng:get_start() + (whole_region_rng[2])*3}
      end
      mnode = gt.feature_node_new(fn:get_seqid(), "membrane_structure",
                                  whole_region_rng[1], whole_region_rng[2],
                                  fn:get_strand())
      mnode:set_attribute("ID", pp_id..".membrane")
      mnode:set_source("TMHMM2.0")
      fn:add_child(mnode)
      for i,part in ipairs(tm.parts) do
        ptype = nil
        if part[1] == "inside" then
          ptype = "cytoplasmic_polypeptide_region"
        elseif part[1] == "outside" then
          ptype = "non_cytoplasmic_polypeptide_region"
        else
          ptype = "transmembrane_helix"
        end
        if fn:get_strand() == '-' then
          whole_region_rng = {node_rng:get_end() - (tonumber(part[3])-1)*3,
                              node_rng:get_end() - (tonumber(part[2])-1)*3}
        else
          whole_region_rng = {node_rng:get_start() + (tonumber(part[2])-1)*3,
                              node_rng:get_start() + (tonumber(part[3])-1)*3}
        end
        pt = gt.feature_node_new(fn:get_seqid(), ptype,
                                 whole_region_rng[1], whole_region_rng[2],
                                 fn:get_strand())
        mnode:add_child(pt)
        pt:set_source("TMHMM2.0")
      end
    end
  end
  return 0
end

-- extract TMH data from TMHMM output and store per protein
for l in io.lines(arg[1]) do
  a,b = string.match(l, "^# (.+)%s+Number of predicted TMHs:%s+(%d+)")
  if a and b then
    v = split(a, ':')[1]
    if not cv.tmhs[v] then
      cv.tmhs[v] = {}
    end
    cv.tmhs[v].number = tonumber(b)
  else
    if string.sub(l, 1, 2) ~= "# " then
      local seq, _, state, start, stop = unpack(split(l, "%s+"))
      seq = split(seq, ':')[1]
      if cv.tmhs[seq] and cv.tmhs[seq].number > 0 then
        if cv.tmhs[seq].parts == nil then
          cv.tmhs[seq].parts = {}
        end
        table.insert(cv.tmhs[seq].parts, {state, start, stop})
      end
    end
  end
end

vis_stream = gt.custom_stream_new_unsorted()
vis_stream.instream = gt.gff3_in_stream_new_sorted(arg[2])
vis_stream.visitor = cv
function vis_stream:next_tree()
   local node = self.instream:next_tree()
  if node then
    node:accept(self.visitor)
  end
  return node
end

out_stream = gt.gff3_out_stream_new(vis_stream)

local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end