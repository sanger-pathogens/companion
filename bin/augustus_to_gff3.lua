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

geneid = nil
transid = nil
exno = 1

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function string:split_sep(sep)
  local sep, fields = sep or ":", {}
  local pattern = string.format("([^%s]+)", sep)
  self:gsub(pattern, function(c) fields[#fields+1] = c end)
  return fields
end

function string:split_ws()
  local fields = {}
  self:gsub("([^%s]+)", function(c) fields[#fields+1] = c end)
  return fields
end

function get_id(attrib)
  id = string.match(attrib, "ID=([^;]+);?")
  return id
end

while true do
  local line = io.read()
  if line == nil then break end

  if string.sub(line, 1, 2) == "##" then
    print(line)
  else
    la = line:split_sep('\t')
    seqid = la[1]
    src = la[2]
    ftype = la[3]
    from = la[4]
    to = la[5]
    score = la[6]
    strand = la[7]
    phase = la[8]
    attrib = la[9]

    if seqid and ftype then
      if ftype == "gene" then
        geneid = get_id(attrib)
        print(line)
      elseif ftype == "transcript" or ftype == "mRNA" then
        print(seqid .. "\t" .. src .. "\tmRNA\t" .. from .. "\t" .. to
                .. "\t" .. score .. "\t" .. strand .. "\t"
                .. phase .. "\t" .. attrib)
        transid = get_id(attrib)
        exno = 1
      elseif ftype == "CDS" then
        print(seqid .. "\t" .. src .. "\texon\t" .. from .. "\t" .. to
                .. "\t" .. score .. "\t" .. strand .. "\t" .. phase
                .. "\tID=" .. transid .. ".exon." .. exno .. ";Parent="
                .. transid)
        exno = exno + 1
        print(line)
      else
        print(line)
      end
    end
  end
end