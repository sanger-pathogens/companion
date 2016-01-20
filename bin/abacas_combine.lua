#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2015 Genome Research Ltd

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
require("lfs")
require("optparse")
local json = require ("dkjson")

op = OptionParser:new({usage="%prog <ABACAS directory> <outfileprefix> "
                                .. "<reference directory> <reference name> <chr prefix> "
                                .. "<bin seqid> <seq prefix>",
                       oneliner="Converts ABACAS2 output file into Companion "
                                  .. "specific formats.",
                       version="0.1"})
options,args = op:parse()

function usage()
  op:help()
  os.exit(1)
end

if #args ~= 7 then
  usage()
end

abacas_dir = args[1]
outfileprefix = args[2]
refdir = args[3]
refname = args[4]
local reffile = io.open(refdir .. '/' .. "references.json", "rb")
if not reffile then
  error("could not open references.json in " .. refdir)
end
local refcontent = reffile:read("*all")
reffile:close()
refs = json.decode(refcontent)
if not refs.species[refname] then
  error("reference " .. refname .. " not found in JSON definition")
end
chrprefix = args[5]
binseqid = args[6]
seqprefix = args[7]

-- scaffold datastructures
scafs = {}
-- pseudochromosomes
pseudochr_seq = {}

-- parse contigs and gaps
for file in lfs.dir(abacas_dir) do
  if file:match("^Res\.(.+)\.gff") then
    local seqid = file:match("^Res\.(.+).contigs\.gff")
    local keys, seqs = get_fasta_nosep(abacas_dir .. "/Res." .. seqid .. ".fna")
    pseudochr_seq[seqid] = seqs[seqid]
    scafs[seqid] = {}
    for l in io.lines(file) do
      if string.sub(l, 1, 1) ~= '#' then
        _, type1, type2, start, stop, frame, strand, _, attr =
                                                         unpack(split(l, "%s+"))
        if not attr then
          io.stderr:write("non-comment line with <9 columns: " .. l .. "\n")
          os.exit(1)
        end
        if type1 == "Contig" then
          local contig_name = attr:match('contig=([^"]+)')
          assert(contig_name)
          -- this is a 'scaffold'
          this_scaf = {start = start,
                       stop = stop,
                       name = contig_name,
                       seq = string.sub(seqs[seqid], start, stop),
                       contigs = {}}
          local i_end = 0
          -- gather all 'contigs'
          while true do
            i_start, i_end = string.find(this_scaf.seq, "[^Nn]+", i_end + 1)
            if i_start == nil then break end
            this_contig = {start = i_start,
                           stop = i_end,
                           seq = string.sub(this_scaf.seq, i_start, i_end)}
            -- sanity check
            if this_contig.seq:match("[Nn]") then
              io.stderr:write("contig with N character encountered\n")
              os.exit(1)
            end
            table.insert(this_scaf.contigs, this_contig)
          end
          table.insert(scafs[seqid], this_scaf)
        elseif type1 == "GAP" then
          -- just a sanity check
          gapseq = string.sub(seqs[seqid], start, stop)
          if gapseq:match("[^Nn]") then
            io.stderr:write("Gap with non-N character encountered\n")
            os.exit(1)
          end
        end
      end
    end
  end
end

ref_target_chromosome_map = {}

-- ensure correct naming for all sequences
newscafs = {}
newkeys = {}
newpseudochr_seq = {}
for k,v in pairs(scafs) do
  local chr = nil
  if refs.species[refname].chr_mapping then
    chr = refs.species[refname].chr_mapping[k]
  end
  if refs.species[refname].chromosome_pattern then
    chr = k:match(refs.species[refname].chromosome_pattern)
  end
  if not chr then
    error("chromosome without mapping encountered: " .. k)
  end
  local newid = chrprefix .. "_" .. chr
  if newscafs[newid] then
    error("new ID " .. newid .. " assigned more than once, "
                     .. "check your refpattern")
  end
  newscafs[newid] = v
  newpseudochr_seq[newid] = pseudochr_seq[k]
  table.insert(newkeys, newid)
  ref_target_chromosome_map[chr] = {k, newid}
end
scafs = newscafs
keys = newkeys
pseudochr_seq = newpseudochr_seq

-- handle 'bin' seqs
scafs[binseqid] = {}
start = 1
stop = 1
binkeys, binseqs = get_fasta_nosep(abacas_dir .. "/Res.abacasBin.fna")
tmp = {}
if #binkeys > 0 then
  table.insert(newkeys, binseqid)
  for k,v in pairs(binseqs) do
    stop = start + string.len(v) - 1
    this_scaf ={start = tonumber(start),
                stop = tonumber(stop),
                name = k,
                seq = v,
                contigs = {}}
    local i_end = 0
    while true do
      i_start, i_end = string.find(this_scaf.seq, "[^Nn]+", i_end + 1)
      if i_start == nil then break end
      this_contig = {start = i_start,
                     stop = i_end,
                     seq = string.sub(this_scaf.seq, i_start, i_end)}
      table.insert(this_scaf.contigs, this_contig)
    end
    table.insert(scafs[binseqid], this_scaf)
    if start > 1 then
      table.insert(tmp, string.rep("N",100))
    end
    start = stop + 1 + 100
    table.insert(tmp, v)
  end
  pseudochr_seq[binseqid] = table.concat(tmp,"")
end

-- open files
pseudochr_fasta_out = io.open(outfileprefix .. ".pseudochr.fasta", "w+")
scaf_fasta_out = io.open(outfileprefix .. ".scafs.fasta", "w+")
scaf_agp_out = io.open(outfileprefix .. ".pseudochr.agp", "w+")
scaf_agp_out:write("##agp-version\t2.0\n")
ctg_fasta_out = io.open(outfileprefix .. ".contigs.fasta", "w+")
ctg_agp_out = io.open(outfileprefix .. ".scafs.agp", "w+")
ctg_agp_out:write("##agp-version\t2.0\n")
ref_target_mapping_out = io.open("ref_target_mapping.txt", "w+")

-- do the output
scaf_i = 1
contig_i = 1
table.sort(newkeys)
-- for all toplevel seqs...
for _,seqid in ipairs(newkeys) do
  pseudochr_fasta_out:write(">" .. seqid .. "\n")
  print_max_width(trim_ns(pseudochr_seq[seqid]), pseudochr_fasta_out, 60)
  print(seqid .. ":  " .. #scafs[seqid])
  local i = 1
  local s_last_stop = 0
  for _, s in ipairs(scafs[seqid]) do
    local scafname = tostring(s.name) -- seqprefix .. "_SCAF" .. string.format("%06d", scaf_i)
    scaf_fasta_out:write(">" .. scafname .. "\n")
    print_max_width(s.seq, scaf_fasta_out, 60)
    if s_last_stop > 0 then
      -- XXX: use no:na for bin chromosomes, otherwise yes:align_xgenus
      if seqid == binseqid then
        scaf_agp_out:write(seqid .. "\t" .. tonumber(s_last_stop)+1 .. "\t"
                        .. tonumber(s.start)-1 .. "\t" .. i .. "\tU\t"
                        .. (tonumber(s.start)-1)-(tonumber(s_last_stop)+1) + 1
                        .. "\tcontig\tno\tna\n")
      else
        scaf_agp_out:write(seqid .. "\t" .. tonumber(s_last_stop)+1 .. "\t"
                          .. tonumber(s.start)-1 .. "\t" .. i .. "\tU\t"
                          .. (tonumber(s.start)-1)-(tonumber(s_last_stop)+1) + 1
                          .. "\tscaffold\tyes\talign_xgenus\n")
      end
      i = i + 1
    end
    scaf_agp_out:write(seqid .. "\t" .. s.start .. "\t"
                       .. s.stop .. "\t" .. i .. "\tW\t" .. scafname
                       .. "\t1\t" .. string.len(s.seq) .. "\t+\n")
    scaf_i = scaf_i + 1
    local j = 1
    local c_last_stop = 0
    for _, c in ipairs(s.contigs) do
      local ctgname = seqprefix .. "_CTG" .. string.format("%06d", contig_i)
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

for k,v in pairs(ref_target_chromosome_map) do
  ref_target_mapping_out:write(k .. "\t" .. v[1] .. "\t" .. v[2] .. "\n")
end
