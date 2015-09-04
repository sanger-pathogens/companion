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

math.randomseed(os.time())

function usage()
  io.stderr:write("Randomly samples a number of single transcript protein coding gene CCs.\n")
  io.stderr:write(string.format("Usage: %s <GFF with gene annotations> " ..
                                "<number of genes to sample>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
end

cv = gt.custom_visitor_new()
cv.out = nil
function cv:visit_feature(fn)
  local gene = false
  local mrna = false
  local nof_transcripts = 0
  local cds = false
  for n in fn:get_children() do
    if n:get_type() == "gene" then
      gene = true
    elseif n:get_type() == "mRNA" then
      mrna = true
      nof_transcripts = nof_transcripts + 1
    elseif n:get_type() == "CDS" then
      cds = true
    end
  end
  if gene and mrna and nof_transcripts == 1 and cds then
    self.out = fn
  else
    self.out = nil
  end
  return 0
end

sample_stream = gt.custom_stream_new_unsorted()
sample_stream.queue = cv.queue
sample_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
sample_stream.visitor = cv
sample_stream.k = tonumber(arg[2])
sample_stream.rnd = {}
sample_stream.i = 1
sample_stream.filled = false
function sample_stream:next_tree()
  if not self.filled then
    -- do reservoir sampling from stream
    local gn = self.instream:next_tree()
    while (gn and #self.rnd < self.k) do
      gn:accept(self.visitor)
      if self.visitor.out then
        table.insert(self.rnd, gn)
        self.i = self.i + 1
      end
      gn = self.instream:next_tree()
    end
    while (gn) do
      gn:accept(self.visitor)
      if self.visitor.out then
        local j = math.random(self.i)
        if j <= self.k then
          self.rnd[j] = gn
        end
        self.i = self.i + 1
      end
      gn = self.instream:next_tree()
    end
    self.filled = true
  end
  if self.filled then
    if table.getn(self.rnd) > 0 then
      return table.remove(self.rnd)
    end
  end
end

out_stream = gt.gff3_out_stream_new(sample_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end