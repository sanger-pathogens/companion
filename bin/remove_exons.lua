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
  io.stderr:write("Removes exon features from gene features in a GFF3 file.\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
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
