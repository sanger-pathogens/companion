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

function chomp(s)
  return string.gsub(s, "\n$", "")
end

function usage()
  io.stderr:write(string.format("Usage: %s [file] [n] [outdirprefix]\n", arg[0]))
  io.stderr:write("Splits FASTA files in to [n] smaller ones.\n")
  os.exit(1)
end

if #arg ~= 3 then
  usage()
end

filename = arg[1]
n = tonumber(arg[2])
outdirprefix = arg[3]

local i, j = 0, 0
hdr = ""
seq = {}

f = io.open(outdirprefix .. "." .. j, "w")
io.output(f)
for l in io.lines(filename) do
  if string.sub(l,1,1) == '>' then
    if i > 0 then
      if i % n == 0 then
        f:close()
        j = j + 1
        i = 0
        f = io.open(outdirprefix .. "." .. j, "w")
        io.output(f)
      end
      io.write(">"..hdr.."\n")
      io.write(table.concat(seq,"\n") .. "\n")
    end
    hdr = string.sub(l, 2)
    seq = {}
    i = i + 1
  else
    table.insert(seq, tostring(chomp(l)))
  end
end
if #seq > 0 then
  io.write(">"..hdr.."\n")
  io.write(table.concat(seq,"\n") .. "\n")
end
f:close()
