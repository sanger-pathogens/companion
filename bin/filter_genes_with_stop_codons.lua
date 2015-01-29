#!/usr/bin/env gt

--[[
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
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
  io.stderr:write(string.format("Usage: %s <GFF annotation> <sequence>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

region_mapping = gt.region_mapping_new_seqfile_matchdescstart(arg[2])

outvis = gt.gff3_visitor_new()

visitor = gt.custom_visitor_new()
visitor.last_seqid = nil
function visitor:visit_feature(fn)
  local ignore = false
  for node in fn:children() do
    if node:get_type() == "mRNA" then
      local protseq = node:extract_and_translate_sequence("CDS", true,
                                                          region_mapping)
      if protseq:sub(1, -2):match("[*+#]") then
        ignore = true
        break
      end
    end
  end
  if not ignore then
    fn:accept(outvis)
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
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

