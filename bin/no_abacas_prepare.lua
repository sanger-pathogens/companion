#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2015 Genome Research Ltd

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
  io.stderr:write("Creates a set of input files without ABACAS input.\n")
  io.stderr:write(string.format("Usage: %s <file> <outfileprefix>\n" , arg[0]))
  os.exit(1)
end

if #arg ~= 2 then
  usage()
end

seqfile = arg[1]
outfileprefix = arg[2]

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

-- scaffold datastructures
scafs = {}

start = 1
keys, seqs = get_fasta_nosep(arg[1])
for hdr, seq in pairs(seqs) do
  if not scafs[hdr] then
    scafs[hdr] = {}
  end
  local seqlen = string.len(seq)
  this_scaf = {start = start,
               stop = start + seqlen - 1,
               seq = string.sub(seq, start, stop),
               contigs = {}}
--print(hdr .. " " .. this_scaf.start .. " " .. this_scaf.stop )
  local i_end = 0
  -- gather all 'contigs'
  while true do
    i_start, i_end = string.find(this_scaf.seq, "[^Nn]+", i_end + 1)
    if i_start == nil then break end
    this_contig = {start = i_start,
                   stop = i_end,
                   seq = string.sub(this_scaf.seq, i_start, i_end)}
--print("    " .. this_contig.start .. " " .. this_contig.stop )
    -- sanity check
    if this_contig.seq:match("[Nn]") then
      io.stderr:write("contig with N character encountered")
      os.exit(1)
    end
    table.insert(this_scaf.contigs, this_contig)
  end
  table.insert(scafs[hdr], this_scaf)
end

-- open files
pseudochr_fasta_out = io.open(outfileprefix .. ".pseudochr.fasta", "w+")
scaf_fasta_out = io.open(outfileprefix .. ".scafs.fasta", "w+")
scaf_agp_out = io.open(outfileprefix .. ".pseudochr.agp", "w+")
scaf_agp_out:write("##agp-version\t2.0\n")
ctg_fasta_out = io.open(outfileprefix .. ".contigs.fasta", "w+")
ctg_agp_out = io.open(outfileprefix .. ".scafs.agp", "w+")
ctg_agp_out:write("##agp-version\t2.0\n")

-- do the output
scaf_i = 1
contig_i = 1
table.sort(keys)
-- for all toplevel seqs...
for _,seqid in ipairs(keys) do
  pseudochr_fasta_out:write(">" .. seqid .. "\n")
  print_max_width(seqs[seqid], pseudochr_fasta_out, 60)
  print(seqid .. ":  " .. #scafs[seqid])
  local i = 1
  local s_last_stop = 0
  for _, s in ipairs(scafs[seqid]) do
    local scafname = seqid
    scaf_fasta_out:write(">" .. scafname .. "\n")
    print_max_width(s.seq, scaf_fasta_out, 60)
    if s_last_stop > 0 then
      scaf_agp_out:write(seqid .. "\t" .. tonumber(s_last_stop)+1 .. "\t"
                        .. tonumber(s.start)-1 .. "\t" .. i .. "\tU\t"
                        .. (tonumber(s.start)-1)-(tonumber(s_last_stop)+1) + 1
                        .. "\tcontig\tno\talign_xgenus\n")
      i = i + 1
    end
    scaf_agp_out:write(seqid .. "\t" .. s.start .. "\t"
                       .. s.stop .. "\t" .. i .. "\tW\t" .. scafname
                       .. "\t1\t" .. string.len(s.seq) .. "\t+\n")
    scaf_i = scaf_i + 1
    local j = 1
    local c_last_stop = 0
    for _, c in ipairs(s.contigs) do
      local ctgname = seqid .. "_CTG" .. string.format("%06d", contig_i)
      ctg_fasta_out:write(">" .. ctgname .. "\n")
      print_max_width(c.seq, ctg_fasta_out, 60)
      if c_last_stop > 0 then
        ctg_agp_out:write(scafname .. "\t" .. tonumber(c_last_stop)+1 .. "\t"
                          .. tonumber(c.start)-1 .. "\t" .. j .. "\tN\t"
                          .. (tonumber(c.start)-1)-(tonumber(c_last_stop)+1) + 1
                          .. "\tscaffold\tyes\tunspecified\n")
        j = j + 1
      end
      ctg_agp_out:write(scafname .. "\t" .. c.start .. "\t"
                        .. c.stop .. "\t" .. j .. "\tF\t" .. ctgname
                        .. "\t1\t" .. string.len(c.seq) .. "\t+\n")
      contig_i = contig_i + 1
      j = j + 1
      c_last_stop = c.stop
    end
    i = i + 1
    s_last_stop = tonumber(s.stop)
  end
end
