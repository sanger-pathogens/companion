#!/usr/bin/env gt

--[[
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
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
  io.stderr:write("For a GFF3 file produced by annotation pipeline, output EMBL file suitable for Chado loader.\n")
  io.stderr:write("Requires an OBO file with a GO definition to add GO terms (pipeline only has IDs).\n")
  io.stderr:write(string.format("Usage: %s <GFF annotation> <GO OBO file> <organism name> <sequence>\n" , arg[0]))
  os.exit(1)
end

if #arg < 3 then
  usage()
end

package.path = gt.script_dir .. "/?.lua;" .. package.path
require("lib")

region_mapping = gt.region_mapping_new_seqfile_matchdescstart(arg[4])

function parse_obo(filename)
  local gos = {}
  local this_id = {}
  local this_name = nil
  for l in io.lines(filename) do
    if string.match(l, "%[Term%]") then
      if #this_id and this_name then
        for _, n in ipairs(this_id) do
          gos[n] = this_name
        end
      end
      this_id = {}
      this_term = nil
    end
    m = string.match(l, "id: (.+)")
    if m then
      table.insert(this_id, m)
    end
    m = string.match(l, "^name: (.+)")
    if m then
      this_name = m
    end
  end
  return gos
end

collect_vis = gt.custom_visitor_new()
collect_vis.lengths = {}
collect_vis.seqs = {}
collect_vis.pps = {}
function collect_vis:visit_feature(fn)
  if fn:get_type() == "polypeptide" then
    local df = fn:get_attribute("Derives_from")
    if df then
      self.pps[df] = fn
    end
  end
  return 0
end
function collect_vis:visit_sequence(sn)
  if sn:get_sequence_length() > 0 then
    self.lengths[sn:get_seqid()] = sn:get_sequence_length()
    self.seqs[sn:get_seqid()] = sn:get_sequence()
  end
  return 0
end

function format_embl_attrib(node, attrib, qualifier, fct)
  if node and node:get_attribute(attrib) then
    for _,v in ipairs(split(node:get_attribute(attrib), ",")) do
      if fct then
        s = fct(v)
      else
        s = gff3_decode(v)
      end
      if s then
        io.write("FT                   /" .. qualifier .. "=\"" .. s .."\"\n")
      end
    end
  end
end

function format_embl_sequence(sequence)
  local a = 0
  local c = 0
  local g = 0
  local t = 0
  local other = 0
  local l = 0
  -- count char distribution
  for ch in sequence:gmatch("%a") do
    ch = string.lower(ch)
    l = l + 1
    if ch == "a" then
      a = a + 1
    elseif ch == "c" then
      c = c + 1
    elseif ch == "g" then
      g = g + 1
    elseif ch == "t" then
      t = t + 1
    else
      other = other + 1
    end
  end
  -- show statistics
  io.write("SQ   Sequence " .. l .. " BP; " .. a .. " A; ".. c .. " C; "..
                                            g .. " G; ".. t .. " T; "..
                                            other .. " other;\n")
  local i = 1
  local pos = 0
  -- format and output sequence
  io.write("     ")
  for c in sequence:gmatch("%a", 10) do
    io.write(c)
    if i % 10 == 0 then
      io.write(" ")
    end
    if i % 60 == 0 then
      io.write(string.format("%9s\n     ", i))
    end
    i = i + 1
  end
  io.write(string.format(string.rep(' ',(80-i%60-(i%60)/10-13)) .. "%10d\n", i-1))
end

embl_vis = gt.custom_visitor_new()
embl_vis.pps = collect_vis.pps
embl_vis.gos = parse_obo(arg[2])
embl_vis.last_seqid = nil
function embl_vis:visit_feature(fn)
  if embl_vis.last_seqid ~= fn:get_seqid() then
    if embl_vis.last_seqid then
      -- sequence output disabled for now
      format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
      io.write("//\n")
      io.output():close()
    end
    embl_vis.last_seqid = fn:get_seqid()
    io.output(fn:get_seqid()..".embl", "w+")
    io.write("ID   XXX; XXX; linear; XXX; XXX; XXX; XXX.\n")
    io.write("XX   \n")
    io.write("DE   " .. arg[3] .. ", " .. fn:get_seqid() .. "\n")
    io.write("XX   \n")
    io.write("AC   XXX.\n")
    io.write("XX   \n")
    io.write("PR   Project:00000000;\n")
    io.write("XX   \n")
    io.write("KW   \n")
    io.write("XX   \n")
    io.write("RN   [1]\n")
    io.write("RA   Authors;\n")
    io.write("RT   Title;\n")
    io.write("RL   Unpublished.\n")
    io.write("XX   \n")
    io.write("OS   " .. arg[3] .. "\n")
    io.write("XX   \n")
    io.write("FH   Key             Location/Qualifiers\n")
    io.write("FH   \n")
    io.write("FT   source          1.." .. collect_vis.lengths[fn:get_seqid()] .. "\n")
    io.write("FT                   /organism=\"" .. arg[3] .. "\"\n")
    io.write("FT                   /mol_type=\"genomic DNA\"\n")
  end

  for node in fn:get_children() do
    if node:get_type() == "mRNA" or node:get_type() == "pseudogenic_transcript" then
      local cnt = 0
      for cds in node:get_children() do
        if cds:get_type() == "CDS" or cds:get_type() == "pseudogenic_exon" then
          cnt = cnt + 1
        end
      end
      io.write("FT   CDS             ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      if cnt > 1 then
        io.write("join(")
      end
      local i = 1
      for cds in node:get_children() do
        if cds:get_type() == "CDS" then
          if i == 1 and fn:get_attribute("Start_range") then
            io.write("<")
          end
          io.write(cds:get_range():get_start())
          io.write("..")
          if i == cnt and fn:get_attribute("End_range") then
            io.write(">")
          end
          io.write(cds:get_range():get_end())
          if i ~= cnt then
            io.write(",")
          end
          i = i + 1
        end
      end
      if cnt > 1 then
        io.write(")")
      end
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      local pp = self.pps[node:get_attribute("ID")]
      format_embl_attrib(pp, "product", "product",
          function (s)
            local pr_a = gff3_extract_structure(s)
            local gprod = pr_a[1].term
            if gprod then
              return gprod
            else
              return nil
            end
          end)
      format_embl_attrib(fn , "ID", "locus_tag", nil)
      if fn:get_type() == "pseudogene" then
        io.write("FT                   /pseudo\n")
        io.write("FT                   /pseudogene=\"unknown\"\n")
      end
      format_embl_attrib(pp, "Dbxref", "db_xref", nil)
      format_embl_attrib(pp, "Dbxref", "EC_number",
          function (s)
            m = string.match(s, "EC:([0-9.-]+)")
            if m then
              return m
            else
              return nil
            end
          end)
      local protseq = node:extract_and_translate_sequence("CDS", true,
                                                          region_mapping)
      io.write("FT                   /translation=\"" .. protseq:sub(1,-2) .."\"\n")
      io.write("FT                   /transl_table=1\n")
      -- not really used in EMBL
      -- format_embl_attrib(pp, "Ontology_term", "Ontology_term")
      -- disabled ortholog output for now (ask maa if its even implemented)
      if pp and false and pp:get_attribute("orthologous_to") then
        for _,v in ipairs(split(pp:get_attribute("orthologous_to"), ",")) do
          io.write("FT                   /ortholog=\"" ..
                          v .. " " .. v .. ";program=OrthoMCL;rank=0\"\n")
        end
      end
      if pp and pp:get_attribute("full_GO") then
        for _,v in ipairs(split(pp:get_attribute("full_GO"), ",")) do
          id = string.match(v, "(GO:[0-9]+)")
          if self.gos[id] then
            io.write("FT                   /GO=\"" ..
                            gff3_decode(v) .. ";term=" .. self.gos[id].."\"\n")
          else
            print("ID not found: " .. id)
          end
        end
      end
    elseif node:get_type() == "tRNA" then
      io.write("FT   tRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      if node:get_attribute("aa") then
        io.write("FT                   /product=\"" .. node:get_attribute("aa") .. " transfer RNA")
        if node:get_attribute("anticodon") then
          io.write(" (" .. node:get_attribute("anticodon") .. ")")
        end
        io.write("\"\n")
      end
      format_embl_attrib(fn , "ID", "locus_tag", nil)
    elseif string.match(node:get_type(), "snRNA") or string.match(node:get_type(), "snoRNA") then
      io.write("FT   ncRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      io.write("FT                   /ncRNA_class=\"" .. node:get_type() .. "\"\n")
      format_embl_attrib(fn , "ID", "locus_tag", nil)
    elseif string.match(node:get_type(), "rRNA") then
      io.write("FT   rRNA            ")
      if node:get_strand() == "-" then
        io.write("complement(")
      end
      io.write(node:get_range():get_start() .. ".." .. node:get_range():get_end())
      if node:get_strand() == "-" then
        io.write(")")
      end
      io.write("\n")
      io.write("FT                   /product=\"" .. node:get_type() .. "\"\n")
      format_embl_attrib(fn , "ID", "locus_tag", nil)
    end
  end
  return 0
end

-- make simple visitor stream, just applies given visitor to every node
visitor_stream = gt.custom_stream_new_unsorted()
function visitor_stream:next_tree()
  local node = self.instream:next_tree()
  if node then
    node:accept(self.vis)
  end
  return node
end

-- get sequences, sequence lengths and derives_from relationships
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = collect_vis
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

-- output EMBL code as we traverse the GFF
visitor_stream.instream = gt.gff3_in_stream_new_sorted(arg[1])
visitor_stream.vis = embl_vis
embl_vis.obo = gos
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end
-- output last seq
format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
io.write("//\n")
io.output():close()

