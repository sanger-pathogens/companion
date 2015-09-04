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
  io.stderr:write("Given an annotation in GFF3 format, outputs a GFF2 file " ..
                  "accepted by AUGUSTUS's gff2gbSmallDNA.pl tool.\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation>\n", arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

cv = gt.custom_visitor_new()
function cv:visit_feature(fn)
  mrna_context = false
  transcript_id = ""
  for n in fn:get_children() do
    if n:get_type() == 'mRNA' then
      mrna_context = true
      transcript_id = n:get_attribute("ID")
    end
    if n:get_type() == 'CDS' then
      if mrna_context then
        score = 0
        if n:get_score() then
          score = n:get_score()
        end
        print(n:get_seqid() .. "\t.\t"
              .. n:get_type() .. "\t"
              .. n:get_range():get_start() .. "\t"
              .. n:get_range():get_end() .. "\t"
              .. score .. "\t"
              .. n:get_strand() .. "\t"
              .. n:get_phase() .. "\t"
              .. "transcript_id \"" .. transcript_id .. "\"")
      end
    end
  end
  return 0
end

local cds_out_stream = visitor_stream_new(gt.gff3_in_stream_new_sorted(arg[1]), cv)
local gn = cds_out_stream:next_tree()
while (gn) do
  gn = cds_out_stream:next_tree()
end
