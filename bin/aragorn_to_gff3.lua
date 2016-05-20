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

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
local src = "Aragorn_1.2.36"

print("##gff-version\t 3")

line1 = nil
seqid = nil
trnano = 1

while true do
  local line2 = io.read()
  if line2 == nil then break end

  m = string.match(line2, "GC = ([0-9.]+)")
  if m then
    gc = m
  end
  c, spos, epos = string.match(line2, "Sequence (%a?)%[(%d+),(%d+)%]")
  m1, m2 = string.match(line2, "tRNA%-([a-zA-Z|().?]+)%(([a-z.]+)%)")
  if m1 then
    aa = m1
  end
  if m2 then
    anticodon = m2
  end
  if string.len(line2) > 0 and string.sub(line2, 1, 2) ~= '--' then
    if not line1 then
      line1 = line2
    elseif string.match(line2, "nucleotides in") then
      seqid = line1
      trnano = 1
    elseif spos or epos then
      local geneid = seqid..".tRNA."..trnano
      trnano = trnano + 1
      if c == 'c' then
        strand = '-'
      else
        strand = '+'
      end
      print(seqid .. "\t" .. src .. "\tgene\t" .. spos .. "\t" .. epos .. "\t.\t" .. strand .. "\t.\tID=" .. geneid)
      print(seqid .. "\t" .. src .. "\ttRNA\t" .. spos .. "\t" .. epos .. "\t.\t" .. strand .. "\t.\tID=" .. geneid ..".trna;Parent=" .. geneid .. ";aa=" .. gff3_encode(aa) .. ";anticodon=" .. anticodon .. ";gc_content=" .. gc)
      aa = nil
      anticodon = nil
    end
  end
  line1 = line2
end
