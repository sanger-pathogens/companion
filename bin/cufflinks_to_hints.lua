#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2015-2016 Genome Research Ltd

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

-- This code is inspired by https://gist.github.com/epaule/1216188

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")
require("optparse")

op = OptionParser:new({usage="%prog <options> < merged.gtf",
                       oneliner="Creates AUGUSTUS hint files from "
                         .. "cufflinks/cuffmerge transcripts.",
                       version="0.1"})
op:option{"-t", action='store', dest='hint_type',
                help="hints type to use (default: E)"}
options,args = op:parse({hint_type='E'})

function usage()
  op:help()
  os.exit(1)
end

function print_gff(seqid, type, start, stop, strand, name)
  print(tostring(seqid) .. "\tCufflinks\t" .. tostring(type)
           .. "\t" .. tostring(start)
           .. "\t" .. tostring(stop) .. "\t.\t" .. tostring(strand)
           .. "\t.\tgrp=" .. name .. ";src=" .. options.hint_type)
end


feats = {}
-- parse from GTF input
for f in gtf_lines(io.stdin) do
  if f.type == 'exon' and f.attribs.transcript_id then
    if not feats[f.attribs.transcript_id] then
      feats[f.attribs.transcript_id] = {}
    end
    table.insert(feats[f.attribs.transcript_id], {seq = f.seqid,
                                                  f.attribs.gene_id,
                                                 tid = f.attribs.transcript_id,
                                                 start = f.from,
                                                 stop = f.to,
                                                 strand=f.strand})
  end
end

-- create hints file
for n,fs in pairs(feats) do
  -- sort by seqid and start position
  table.sort(fs, function(a,b)
    if a.seq == b.seq then
      return (a.start < b.start)
    else
      if a.seq < b.seq then
        return true
      else
        return false
      end
    end
  end)

  for i,f in ipairs(fs) do
    if i > 1 then
      local intron_start = fs[i-1].stop + 1
      local intron_end = fs[i].start - 1
      print_gff(f.seq, "intron", intron_start, intron_end, f.strand, f.tid)
    end
    if i == 0 or i == #fs then
      print_gff(f.seq, "exonpart", f.start, f.stop, f.strand, f.tid)
    else
      print_gff(f.seq, "exonpart", f.start, f.stop, f.strand, f.tid)
    end
  end
end

