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
  io.stderr:write(string.format("Usage: %s <GFF>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

remove_codons_visitor = gt.custom_visitor_new()
function remove_codons_visitor:visit_feature(fn)
  for node in fn:children() do
    if node:get_type() == "exon"
         or node:get_type() == "intron"
         or node:get_type() == "start_codon"
         or node:get_type() == "stop_codon" then
      fn:remove_leaf(node)
    end
  end
  return 0
end

function join_first(fn)
-- join intron with first CDS/exon, drop first intron
  if fn:get_strand() == '+' then
    local rng = nil
    for i in fn:children() do
      if i:get_type() == 'intron' then
        rng = i:get_range()
        fn:remove_leaf(i)
        break
      end
    end
    if rng then
      for c in fn:children() do
        if c:get_type() == 'CDS' then
          c:set_range(c:get_range():join(rng))
          break
        end
      end
      for c in fn:children() do
        if c:get_type() == 'exon' then
          c:set_range(c:get_range():join(rng))
          break
        end
      end
    end
  end
  if fn:get_strand() == '-' then
    local lintron, lcds, lexon = nil, nil, nil
    for i in fn:children() do
      if i:get_type() == 'intron' then
        lintron = i
      end
    end
    if lintron then
      for c in fn:children() do
        if c:get_type() == 'CDS' then
          lcds = c
        end
      end
      for c in fn:children() do
        if c:get_type() == 'exon' then
          lexon = c
        end
      end
      if lcds and lexon then
        lcds:set_range(lexon:get_range():join(lintron:get_range()))
        lexon:set_range(lexon:get_range():join(lintron:get_range()))
        fn:remove_leaf(lintron)
      end
    end
  end
end

function join_last(fn)
  -- join first intron with first CDS/exon, drop it afterwards
  if fn:get_strand() == '-' then
    rng = nil
    for i in fn:children() do
      if i:get_type() == 'intron' then
        rng = i:get_range()
        fn:remove_leaf(i)
        break
      end
    end
    if rng then
      for c in fn:children() do
        if c:get_type() == 'CDS' then
          c:set_range(c:get_range():join(rng))
          break
        end
      end
      for c in fn:children() do
        if c:get_type() == 'exon' then
          c:set_range(c:get_range():join(rng))
          break
        end
      end
    end
  end
  if fn:get_strand() == '+' then
    local lintron, lcds, lexon = nil, nil, nil
    for i in fn:children() do
      if i:get_type() == 'intron' then
        lintron = i
      end
    end
    if lintron then
      for c in fn:children() do
        if c:get_type() == 'CDS' then
          lcds = c
        end
      end
      for c in fn:children() do
        if c:get_type() == 'exon' then
          lexon = c
        end
      end
      if lcds and lexon then
        lcds:set_range(lexon:get_range():join(lintron:get_range()))
        lexon:set_range(lexon:get_range():join(lintron:get_range()))
        fn:remove_leaf(lintron)
      end
    end
  end
end

mark_partial_visitor = gt.custom_visitor_new()
function mark_partial_visitor:visit_feature(fn)
  has_start, has_stop, has_intron = false, false, false
  if fn:get_type() == 'gene' then
    -- determine if gene is partial
    for fn2 in fn:children() do
      if fn2:get_type() == 'mRNA' then
        local nof_children = 0
        for n in fn2:children() do
          nof_children = nof_children + 1
        end
        local i = 1
        for n in fn2:children() do
          if n:get_type() == 'start_codon' then
            has_start = true
          elseif n:get_type() == 'stop_codon' then
            has_stop = true
          elseif n:get_type() == 'intron' then
            has_intron = true
          end
        end
        -- set corresponding attributes
        if not has_start then
          fn:add_attribute("fiveEndPartial","true")
          if fn:get_strand() == "+" then
            fn:add_attribute("Start_range",".,.")
          elseif fn:get_strand() == "-" then
            fn:add_attribute("End_range",".,.")
          end
        end
        if not has_stop then
          fn:add_attribute("threeEndPartial","true")
          if fn:get_strand() == "-" then
            fn:add_attribute("Start_range",".,.")
          elseif fn:get_strand() == "+" then
            fn:add_attribute("End_range",".,.")
          end
        end
      end
      -- disable the leading/trailing intron joining for now
      if not has_start and has_intron then
        -- join_first(fn2)
      end
      if not has_stop and has_intron then
        -- join_last(fn2)
      end
    end
  end
  return 0
end

-- setup generic visitor stream
vis_stream = gt.custom_stream_new_unsorted()
vis_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
vis_stream.mark_partial_visitor = mark_partial_visitor
vis_stream.remove_codons_visitor = remove_codons_visitor
function vis_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.mark_partial_visitor)
    node:accept(self.remove_codons_visitor)
  end
  return node
end

out_stream = gt.gff3_out_stream_new(vis_stream)
local gn = out_stream:next_tree()
while (gn) do
  gn = out_stream:next_tree()
end