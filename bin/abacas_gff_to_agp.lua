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
  io.stderr:write(string.format("Usage: %s <ABACAS gff> [<ABACAS gff> ...]\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

print("##agp-version\t2.0")
for i = 1,#arg do
--  print("# ABACAS source file: " .. arg[i])
  local seqid = nil
  seqid = arg[i]:match("Res.(.+).contigs.gff")
  j = 1
  for l in io.lines(arg[i]) do
    if string.sub(l, 1, 1) ~= '#' then
      _, type1, type2, start, stop, frame, strand, _, attr = unpack(split(l, "%s+"))
      if not attr then
        io.stderr:write("non-comment line with less than 9 columns: " .. l .. "\n")
        os.exit(1)
      end
      if type1 == "Contig" then
        contig_id = attr:match("contig=([^\"]+)\"")
        size = attr:match("+%((%d+)%)")
        if attr:match("REVERSED") then
          strand = "-"
        else
          strand = "+"
        end
        print(seqid .. "\t" .. start .. "\t" .. stop .. "\t" .. j .. "\t"
              .. "P" .. "\t" .. contig_id .. "\t" .. 1 .. "\t" .. size .. "\t"
              .. strand)
--        print(tonumber(stop) - tonumber(start) + 1)
 --       assert(tonumber(size) == tonumber(stop) - tonumber(start) )
        j = j + 1
      elseif type1 == "GAP" then
        if l:match("label=\"GAP2\"") then
          gaptype = "U"
        else
          gaptype = "N"
        end
        print(seqid .. "\t" .. start .. "\t" .. stop .. "\t" .. j .. "\t"
              .. gaptype .. "\t" .. tonumber(stop) - tonumber(start) + 1 .. "\t"
              .. "contig" .. "\t" .. "no" .. "\t" .. "na")
      end
      j = j + 1
    end
  end
end