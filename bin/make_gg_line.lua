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
  io.stderr:write("Creates a line in a .gg file suitable as input for OrthoMCL.\n")
  io.stderr:write(string.format("Usage: #{$0} <speciescode> <protein FASTA file>\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
  os.exit(1)
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function string:chomp()
  return self:gsub("\n$", "")
end

speciescode = arg[1]
prots = {}

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
    if string.match(l, " ") then
      io.stderr:write("space in header: " .. desc)
      os.exit(1)
    end
    table.insert(prots,desc)
  end
end

table.sort(prots)
print(speciescode .. ": " .. table.concat(prots," "))
