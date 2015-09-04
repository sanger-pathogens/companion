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
  io.stderr:write(string.format("Usage: %s <GFF annotation> [wrong]\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

if arg[2] then
  only_wrong = true
else
  only_wrong = false
end

-- simple linear score for strand membership
function strandscore(strand)
  if strand == '+' then
    return 1
  elseif strand == '-' then
    return -1
  else
    return 0
  end
end

-- sets implemented as tables, "show" contains the accepted genes, "out" rejected genes
show = {}
out = {}
-- XXX TODO: handle polypeptide features + children

function handle_wrong(fn)
  show[fn] = nil -- delete it from the set of accepted genes
  out[fn] = true -- put it into the set of rejected genes
end

gff3vis = gt.gff3_visitor_new()

visitor = gt.custom_visitor_new()
visitor.last_seqid = nil
visitor.last_feat = nil
visitor.last_score = 0
-- window size, must be evenly divisible by 2
visitor.r = 6
visitor.i = 0
visitor.polycis_start = nil
visitor.polycis_end = nil
visitor.queue = {}
function visitor:visit_feature(fn)
  local coding = false
  show[fn] = true

  -- re-initialize for each input sequence
  if fn:get_seqid() ~= visitor.last_seqid then
    visitor.strandscore = 0
    visitor.i = 1
    visitor.last_seqid = fn:get_seqid()
    visitor.queue = {}
    visitor.last_feat = nil
  end
  -- only act on coding genes
  for node in fn:children() do
    if node:get_type() == "mRNA" then
      coding = true
    end
  end
  if coding then
    local wrong = false
    local tp = false
    -- is our window full?
    visitor.strandscore = visitor.strandscore + strandscore(fn:get_strand())
    if visitor.i >= visitor.r then
      r = table.remove(visitor.queue, 1)
      visitor.strandscore = visitor.strandscore - strandscore(r:get_strand())
      -- determine if middle node agrees with strand context
      if (visitor.queue[visitor.r/2]:get_strand() == '-'
            and visitor.strandscore > 0) or
         (visitor.queue[visitor.r/2]:get_strand() == '+'
            and visitor.strandscore < 0) then
         wrong = true
      end
      if wrong then
        handle_wrong(visitor.queue[visitor.r/2])
      end
    end
    table.insert(visitor.queue, fn)

    visitor.i = visitor.i + 1
    visitor.last_feat = this
    visitor.last_score = visitor.strandscore
  end
  return 0
end
function visitor:visit_comment(n)
  n:accept(gff3vis)
end
function visitor:visit_region(n)
  n:accept(gff3vis)
end
function visitor:visit_sequence(n)
  n:accept(gff3vis)
end
function visitor:visit_meta(n)
  n:accept(gff3vis)
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

-- decide what to show: removed or fixed annotations
if only_wrong then
  coll = out
else
  coll = show
end

for k,_ in pairs(coll) do
  k:accept(gff3vis)
end

