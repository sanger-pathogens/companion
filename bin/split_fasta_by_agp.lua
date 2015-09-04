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
  io.stderr:write(string.format("Usage: %s <FASTA> <AGP>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
local lfs = require("lfs")

-- targetdir = arg[3]
keys, seqs = get_fasta(arg[1])

--[[
if not lfs.attributes(targetdir, "mode") then
  ret, err = lfs.mkdir(targetdir)
  if not ret then
    io.stderr:write(err .. "\n")
    os.exit(1)
  end
end]]

-- can be made MUCH faster, but will do for now
function get_seq(seqs, id)
  for k,v in pairs(seqs) do
    if k:match(id) then
      return v
    end
  end
  io.stderr:write("could not find sequence matching object " .. id
                    .. " in input\n")
  os.exit(1)
end

for l in io.lines(arg[2]) do
  if string.sub(l, 1, 1) ~= '#' then
    obj, obj_s, obj_e, part_n, comp_type, c6, c7, c8, c9 = unpack(split(l, "%s+"))
    if not c9 then
      io.stderr:write("non-comment line with less than 9 columns: " .. l .. "\n")
      os.exit(1)
    end
    if comp_type == 'U' or comp_type == 'N' then
      -- gap
    else
      print(">" .. c6)
      -- pad output file if not the whole contig is part of the scaffold
      for i = 1,(tonumber(c7)-1) do
        io.write("N")
      end
      seq = string.sub(get_seq(seqs, obj), obj_s, obj_e)
      if c9 == '-' then
        seq = revcomp(seq)
      end
      io.write(seq)
      io.write("\n")
    end
  end
end
