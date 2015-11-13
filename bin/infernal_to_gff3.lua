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

-- TODO: make this dynamic
models = {}
models.rRNA = {"SSU_rRNA_eukarya","5S_rRNA", "5_8S_rRNA"}
models.snRNA = {"U1", "U2", "U4", "U5" "U6"}
models.snoRNA = {"snoTBR5", "snoTBR17", "snoTBR7"}

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

function qry2type(qry)
  for m in pairs(models) do
    for _,s in ipairs(models[m]) do
      if qry == s then
        return m
      end
    end
  end
  return qry
end

print("##gff-version\t 3")
local i = 1
while true do
  local line = io.read()
  if line == nil then break end

  if string.len(line) > 0 and string.sub(line, 1, 1) ~= "#" then
    la = split(line, '%s+')
    seqid = la[1]
    seqacc = la[2]
    qry = la[3]
    qryacc = la[4]
    mfrom = la[6]
    mto = la[7]
    sfrom = la[8]
    sto = la[9]
    strand = la[10]
    trunc = la[11]
    gc = la[13]
    score = la[15]
    evalue = la[16]
    inc = la[17]

    if strand == '-' then
      sfrom, sto = sto, sfrom
    end
    print(seqid .. "\tGenomeTools\tgene\t" .. sfrom .. "\t" .. sto .. "\t"
            .. score .. "\t" .. strand .. "\t.\tID=ncRNA" .. i)
    print(seqid .. "\tINFERNAL\t" .. qry2type(qry) .. "\t" .. sfrom .. "\t"
            .. sto .. "\t" .. score .. "\t" .. strand .. "\t.\tID=ncRNA" .. i
            .. ":" .. qry2type(qry) .. ";Parent=ncRNA"
            .. i ..";Name=" .. qry .. ";gc=" .. gc .. ";evalue=" .. evalue
            ..";score=" .. score ..";model_name=" .. qryacc .. ";model_range="
            .. mfrom .. "-" .. mto)
    i = i + 1
  end
end