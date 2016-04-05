#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2016 Genome Research Ltd

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

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

-- overlapping cluster handler, picks feature with best e-value
function reconcile_overlaps(cluster)
  local best = nil
  local min_eval = nil
  for _,n in ipairs(cluster) do
    for c in n:children() do
      this_eval = tonumber(c:get_attribute('evalue'))
      if this_eval then
        if not min_eval or this_eval < min_eval then
          best = n
          min_eval = this_eval
        end
      end
    end
  end
  assert(best)
  return {best}
end

instream = infernal_in_stream_new(io.stdin)
ovlstream = overlap_stream_new(instream, nil, reconcile_overlaps)
outstream = gt.gff3_out_stream_new_retainids(ovlstream)

local n = outstream:next_tree()
while n do
  n = outstream:next_tree()
end