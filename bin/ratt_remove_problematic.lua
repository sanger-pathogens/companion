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
  io.stderr:write("\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation> <*.Report.txt> " ..
                                "[<*.Report.txt> ...]\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

problematic_genes = {}

-- parse problematic genes from report file
for i = 2,#arg do
  for l in io.lines(arg[i]) do
    geneid,_,_,_,_,_,_,_,errorStill,StartStillBad,StopStillBad,
      frameshiftsStill,JoinExons,PossiblePseudo,CorrectionLog = unpack(split(l, "\t"))
    if geneid ~= 'Gene_ID' and PossiblePseudo then
      if tonumber(errorStill) + tonumber(StartStillBad) + tonumber(StopStillBad)
        + tonumber(frameshiftsStill) + tonumber(JoinExons) + tonumber(PossiblePseudo) > 0 then
        problematic_genes[geneid] = true
      end
      if CorrectionLog:match("MANUAL") then
        problematic_genes[geneid] = true
      end
    end
  end
end

-- find them in stream of nodes
vis = gt.custom_visitor_new()
function vis:visit_feature(fn)
  self.ok = true
  if fn:get_attribute("ratt_ortholog")
    and problematic_genes[fn:get_attribute("ratt_ortholog")] then
    self.ok = false
  end
  -- check for overlapping exons per transcripts -> do not allow those
  if self.ok then
    for n in fn:children() do
      if not self.ok then
        break
      end
      if n:get_type() == "mRNA" then
        for c1 in n:children() do
          if not self.ok then
            break
          end
          if c1:get_type() == 'CDS' then
            for c2 in n:children() do
              if c2:get_type() == 'CDS' then
                if c1 ~= c2 and c1:get_range():overlap(c2:get_range()) then
                  self.ok = false
                  break
                end
              end
            end
          end
        end
      end
    end
  end
end
function vis:visit_sequence(gn)
  self.ok = true
end
function vis:visit_region(gn)
  self.ok = true
end
function vis:visit_meta(gn)
  self.ok = true
end
function vis:visit_comment(gn)
  self.ok = true
end

-- skip these nodes in output stream
flt_stream = gt.custom_stream_new_unsorted()
flt_stream.queue = queue
flt_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
flt_stream.vis = vis
function flt_stream:next_tree()
  local gn = self.instream:next_tree()
  if gn then
    gn:accept(self.vis)
  end
  while gn and not self.vis.ok do
    gn = self.instream:next_tree()
    if gn then
      gn:accept(self.vis)
    end
  end
  return gn
end

out_stream = gt.gff3_out_stream_new(flt_stream)
local gn = out_stream:next_tree()
while gn do
  gn = out_stream:next_tree()
end