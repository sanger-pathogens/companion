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
  io.stderr:write(string.format("Usage: %s <base GTF> <AGP> "
                                .. "<base FASTA> <target FASTA> "
                                .. "[check]\n" , arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

gff = arg[1]
agp = arg[2]
if arg[3] then
  bkeys, bseqs = get_fasta(arg[3])
end
if arg[4] then
  tkeys, tseqs = get_fasta(arg[4])
end

function get_seq(seqs, id)
  for k,v in pairs(seqs) do
    if k:match("^"..id.."$") then
      return v
    end
  end
  io.stderr:write("could not find sequence matching object " .. id
                    .. " in input\n")
  os.exit(1)
end

-- load and index AGP
local mappings = {}
local lineno = 1
for l in io.lines(arg[2]) do
  if string.sub(l, 1, 1) ~= '#' then
    obj, obj_s, obj_e, part_n, type, c6, c7, c8, c9 = unpack(split(l, "%s+"))
    if not c6 then
      error("could not parse line " .. lineno .. "of AGP file " .. arg[2])
    end
    if type == 'N' or type == 'U' then
      -- ignore gaps
    else
      mappings[c6] = {obj=obj, obj_s=tonumber(obj_s),
                                        obj_e=tonumber(obj_e),
                                        s_s=tonumber(c7), s_e=tonumber(c8),
                                        strand=c9}
    end
  end
  lineno = lineno + 1
end

-- process GTF
lineno = 1
for l in io.lines(arg[1]) do
  if string.sub(l, 1, 1) ~= '#' then
    seqid, src, type, start, stop, score, strand, phase, attrs = unpack(split(l, "\t"))
    if not attrs then
      error("could not parse line " .. lineno .. "of GTF file " .. arg[1])
    end
    -- is this sequence/fragment part of a layout?
    if mappings[seqid] then
      local outstrand = strand
      local new_rng = nil
      -- was the fragment used in reverse orientation?
      if mappings[seqid].strand == "-" then
        -- transform coordinates for reverse
        new_rng = gt.range_new((mappings[seqid].s_e - stop + 1)
                                 + mappings[seqid].obj_s - 1,
                               (mappings[seqid].s_e - start + 1)
                                 + mappings[seqid].obj_s - 1)
        -- flip strands
        if strand == "+" then
          outstrand = "-"
        elseif strand == "-" then
          outstrand = "+"
        end
      else
        -- otherwise just transform coordinates by offsetting
        new_rng = gt.range_new(start + mappings[seqid].obj_s - 1,
                               stop + mappings[seqid].obj_s - 1)
      end
      assert(new_rng)
      assert(outstrand)
      print(mappings[seqid].obj .. "\t" .. src .. "\t" .. type .. "\t" .. new_rng:get_start() .. "\t" .. new_rng:get_end() .. "\t"
             .. score .. "\t" .. outstrand .. "\t" .. phase .. "\t" .. attrs)
      if arg[5] then   -- optionally, check seqs -- must stay the same!
        seq1 = get_seq(bseqs, seqid):sub(rng:get_start(), rng:get_end())
        seq2 = get_seq(tseqs, mappings[seqid].obj):sub(new_rng:get_start(),
                       new_rng:get_end())
        if  mappings[seqid].strand == "-" then
          seq2 = revcomp(seq2)
        end
        if seq1 ~= seq2 then
          print(seq1)
          print(seq2)
        end
        assert(seq1 == seq2)
      end
    else
      -- this sequence is unassembled, do not change its coordinates
      io.stderr:write("no mapping for seqid " .. seqid .. "\n")
      --for node in fn:children() do
      --  node:change_seqid(seqid)
      --end
    end
  end
  lineno = lineno + 1
end
