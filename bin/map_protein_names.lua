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
  io.stderr:write("Assigns simple numbered headers to FASTA files, saves originals in map file.\n")
  io.stderr:write(string.format("Usage: %s <speciescode> <protein FASTA file> <mapfile>\n" , arg[0]))
  os.exit(1)
end

if #arg < 3 then
  usage()
  os.exit(1)
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function string:chomp()
  return self:gsub("\n$", "")
end

speciescode = arg[1]
mapfile = arg[3]

prots = {}
i = 1

for l in io.lines(arg[2]) do
  l = l:chomp()
  if string.sub(l, 1, 1) == '>' then
    if string.len(l) == 1 then
      io.stderr:write("empty header")
      os.exit(1)
    end
    desc = string.sub(l, 2)
    if prots[desc] then
      io.stderr:write("duplicate product: " .. desc)
      os.exit(1)
    end
    prots[desc] = speciescode .. "_" .. i
    print(">" .. prots[desc])
    i = i + 1
  else
    print(l)
  end
end

local f,err = io.open(mapfile,"w")
if not f then
  return print(err)
end
for k, v in pairs(prots) do
  f:write(v .. "\t" .. k .. "\n")
end
f:close()