#!/usr/bin/env gt

--[[
  Author: Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014-2016 Genome Research Ltd

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
require("optparse")

op = OptionParser:new({usage="%prog <options> <GFF annotation> <GO OBO file> <organism name> <sequence> [GAF]",
                       oneliner="For a Sanger Pathogens style GFF3 file, output EMBL format.",
                       version="0.2"})
op:option{"-e", action='store_true', dest='embl_compliant',
                help="output reduced 'ENA compliant' format"}
op:option{"-o", action='store_true', dest='only_on_seq',
                help="ignore features for which no sequence is found"}
op:option{"-s", action='store_true', dest='contig_as_source',
                help="emit 'source' features for 'contig' input features"}
op:option{"-p", action='store', dest='projectid',
                help="ENA project ID (e.g. 'PRJEB1234')"}
op:option{"-f", action='store', dest='gaf_level',
                help="value of GAF 'DB Object Type' column (transcript/gene)"}
op:option{"-1", action='store_true', dest='add_onebase_gap',
                help="add gap features for single-base Ns"}
op:option{"-d", action='store_true', dest='remove_duplicate_gene',
                help="remove duplicate gene qualifiers"}
op:option{"-n", action='store', dest='note',
                help="add custom note to features"}
op:option{"-t", action='store_true', dest='trim_flanking_n',
                help="trim leading/trailing Ns in sequences, " ..
                     "adjusting coordinates"}
options,args = op:parse({embl_compliant=false, only_on_seq=false,
                         contig_as_source=false, projectid='00000000',
                         add_onebase_gap=false, remove_duplicate_gene=false,
                         trim_flanking_n=false, gaf_level='gene', note=nil})

function usage()
  op:help()
  os.exit(1)
end

if #args < 4 then
  usage()
end

embl_compliant = options.embl_compliant ~= false
contig_as_source = options.contig_as_source ~= false
seqfile = args[4]
gaffile = args[5]
obofile = args[2]
gff3file = args[1]
organismname = args[3]

region_mapping = gt.region_mapping_new_seqfile_matchdescstart(seqfile)

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
  local last_line = true
  -- format and output sequence
  io.write("     ")
  for c in sequence:gmatch("%a", 10) do
    io.write(c)
    if i % 10 == 0 then
      io.write(" ")
    end
    if i % 60 == 0 then
      io.write(string.format("%9s\n", i))
      -- this is to address the situation where a sequence fills the last line
      -- completely
      if i == l then
        last_line = false
      else
        io.write("     ")
      end
    end
    i = i + 1
  end
  if last_line then
    io.write(string.format(string.rep(' ',(80-i%60-(i%60)/10-13)) .. "%10d\n", i-1))
  end
end

embl_vis = gt.custom_visitor_new()
embl_vis.pps = collect_vis.pps
embl_vis.gos = parse_obo(obofile)
embl_vis.last_seqid = nil

function print_generic(fn, mtype, etype, adjust)
  local cnt = 0
  for c in fn:children() do
    if c:get_type() == mtype then
      cnt = cnt + 1
    end
  end
  if cnt == 0 then
    return false
  end
  -- TODO: make spacing automatic
  io.write("FT   " .. etype .. "           ")
  if fn:get_strand() == "-" then
    io.write("complement(")
  end
  if cnt > 1 then
    io.write("join(")
  end
  local i = 1
  local coding_length = 0
  local start_phase = 0
  local end_phase = 0
  local offset = 0
  if adjust then
    offset = adjust[c:get_seqid()]
  end
  for c in fn:children() do
    if c:get_type() == mtype then
      if i == 1 and fn:get_attribute("Start_range") then
        io.write("<")
      end
      io.write(c:get_range():get_start() - offset)
      io.write("..")
      if i == cnt and fn:get_attribute("End_range") then
        io.write(">")
      end
      io.write(c:get_range():get_end() - offset)
      if i ~= cnt then
        io.write(",")
      end
      coding_length = coding_length + c:get_range():length()
      i = i + 1
    end
  end
  if cnt > 1 then
    io.write(")")
  end
  if fn:get_strand() == "-" then
    io.write(")")
  end
  io.write("\n")
  return true
end

function output_note()
  if options.note then
    io.write("FT                   /note=\"".. options.note .. "\"\n")
  end
end

function is_multi_transcript_gene(fn)
  if fn:get_type() == 'gene' then
    local nof_mrnas = 0
    for n in fn:children() do
      if n:get_type() == 'mRNA' then
        nof_mrnas = nof_mrnas + 1
      end
    end
    return (nof_mrnas > 1)
  else
    return false
  end
end

function output_multi_transcript_qualifiers(gene, mrna)
  local other_mrnas = {}
  local gene_id = gene:get_attribute("ID")
  local transcript_id = mrna:get_attribute("ID")
  for n in gene:children() do
    if n:get_type() == 'mRNA' then
      if mrna:get_attribute("ID") ~= n:get_attribute("ID") then
        table.insert(other_mrnas, n:get_attribute("ID"))
      end
    end
  end
  io.write("FT                   /shared_id=\"" .. gene_id .. "\"\n")
  io.write("FT                   /systematic_id=\"" .. transcript_id .. "\"\n")
  for _, v in ipairs(other_mrnas) do
    io.write("FT                   /other_transcript=\"" .. v .. "\"\n")
  end
end

function output_single_base_gaps(seq)
  local i_end = 0
  while true do
    i_start, i_end = string.find(seq, "[Nn]+", i_end + 1)
    if i_start == nil then break end
    if i_start == i_end then
      io.write("FT   assembly_gap    ")
      io.write(i_start .. ".." .. i_end)
      io.write("\n")
      io.write("FT                   /estimated_length=" .. (i_end - i_start + 1) .. "\n")
      io.write("FT                   /gap_type=\"unknown\"\n")
    end
  end
end

function output_embl_header(fn)
  if embl_vis.last_seqid ~= fn:get_seqid() then
    if embl_vis.last_seqid then
      if not collect_vis.seqs[embl_vis.last_seqid] then
        warning("sequence " .. embl_vis.last_seqid .. " not in collected set")
      end
      if options.add_onebase_gap then
        output_single_base_gaps(collect_vis.seqs[embl_vis.last_seqid])
      end
      format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
      io.write("//\n")
      io.output():close()
    end
    embl_vis.last_seqid = fn:get_seqid()
    io.output(fn:get_seqid()..".embl", "w+")
    if embl_compliant then
      io.write("ID   XXX; XXX; linear; "
                 .. "genomic DNA; STD; UNC; "
                 .. tostring(collect_vis.lengths[fn:get_seqid()]) .." BP.\n")
      io.write("XX   \n")
      io.write("AC   XXX;\n")
      io.write("XX   \n")
      io.write("AC * _" .. fn:get_seqid() .. ";\n")
      io.write("XX   \n")
    else
      io.write("ID   " .. fn:get_seqid() .. "; SV 1; linear; "
                 .. "genomic DNA; STD; UNC; "
                 .. tostring(collect_vis.lengths[fn:get_seqid()]) .." BP.\n")
      io.write("XX   \n")
      io.write("AC   " .. fn:get_seqid() .. ";\n")
      io.write("XX   \n")
    end
    io.write("PR   Project:" .. options.projectid .. ";\n")
    io.write("XX   \n")
    io.write("DE   " .. organismname .. ", " .. fn:get_seqid() .. ".\n")
    io.write("XX   \n")
    io.write("KW   \n")
    io.write("XX   \n")
    io.write("RN   [1]\n")
    io.write("RA   Authors;\n")
    io.write("RT   Title;\n")
    io.write("RL   Unpublished.\n")
    io.write("XX   \n")
    io.write("OS   " .. organismname .. "\n")
    io.write("XX   \n")
    io.write("FH   Key             Location/Qualifiers\n")
    io.write("FH   \n")
    if not contig_as_source then
      io.write("FT   source          1.." .. collect_vis.lengths[fn:get_seqid()] .. "\n")
      io.write("FT                   /organism=\"" .. organismname .. "\"\n")
      io.write("FT                   /mol_type=\"genomic DNA\"\n")
    end
  end
end

function embl_vis:output_gene_qualifier(genesym, geneid)
  if not genesym then
    genesym = split(geneid,":")[1]
  end
  if not options.remove_duplicate_gene or not self.seen_gene_names[genesym] then
    io.write("FT                   /gene=\"".. genesym .. "\"\n")
  end
end

function embl_vis:visit_feature(fn)
  if collect_vis.seqs[fn:get_seqid()] then
    output_embl_header(fn)
    local offset = 0
    if self.adjust then
      offset = self.adjust[fn:get_seqid()]
    end
    for node in fn:get_children() do
      if node:get_type() == "mRNA" or node:get_type() == "pseudogenic_transcript" then
        local cnt = 0
        for cds in node:get_children() do
          if cds:get_type() == "CDS" or cds:get_type() == "pseudogenic_exon" then
            cnt = cnt + 1
          end
        end
        if cnt == 0 then
          break
        end
        io.write("FT   CDS             ")
        if node:get_strand() == "-" then
          io.write("complement(")
        end
        if cnt > 1 then
          io.write("join(")
        end
        local i = 1
        local coding_length = 0
        local start_phase = 0
        local end_phase = 0
        for cds in node:get_children() do
          if cds:get_type() == "CDS" or cds:get_type() == "pseudogenic_exon" then
            if i == 1 and fn:get_attribute("Start_range") then
              if fn:get_strand() == '-' then
                end_phase = tonumber(cds:get_phase())
              else
                start_phase = tonumber(cds:get_phase())
              end
              io.write("<")
            end
            io.write(cds:get_range():get_start() - offset)
            io.write("..")
            if i == cnt and fn:get_attribute("End_range") then
              if fn:get_strand() == '-' then
                start_phase = tonumber(cds:get_phase())
              else
                end_phase = tonumber(cds:get_phase())
              end
              io.write(">")
            end
            io.write(cds:get_range():get_end() - offset)
            if i ~= cnt then
              io.write(",")
            end
            coding_length = coding_length + cds:get_range():length()
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
        format_embl_attrib(fn , "ID", "locus_tag",
          function (s)
            return split(s,':')[1]
          end)
        if fn:get_type() == "pseudogene" then
          io.write("FT                   /pseudo\n")
          --io.write("FT                   /pseudogene=\"unknown\"\n")
        end
        output_note()
        local target_product = "product"
        if fn:get_type() == "pseudogene" and embl_compliant then
          target_product = "note"
        end
        format_embl_attrib(pp, "product", target_product,
          function (s)
            local pr_a = gff3_extract_structure(s)
            local gprod = nil
            if #pr_a > 0 and pr_a[1].term then
              gprod = pr_a[1].term
            end
            -- use preferred term if defined
            for _,pr in ipairs(pr_a) do
              if pr.is_preferred and pr.term then
                gprod = pr.term
                break
              end
            end
            if gprod then
              return gprod
            else
              return nil
            end
        end)
        format_embl_attrib(pp, "Dbxref", "EC_number",
            function (s)
              m = string.match(s, "EC:([0-9.-]+)")
              if m then
                return m
              else
                return nil
              end
            end)
        -- add gene to 'unroll' multiple spliceforms
        local geneid = fn:get_attribute("ID")
        local genesym = fn:get_attribute("Name")
        self:output_gene_qualifier(geneid, genesym)

        -- output multi-transcript qualifiers, as outlined in
        -- http://mediawiki.internal.sanger.ac.uk/index.php/Annotation_of_alternately_spliced_genes
        if not embl_compliant then
          if is_multi_transcript_gene(fn) then
            output_multi_transcript_qualifiers(fn, node)
          else
            if node:get_attribute("ID") then
              io.write("FT                   /systematic_id=\"" .. node:get_attribute("ID") .. "\"\n")
            end
          end
        end

        -- translation
        local protseq = nil
        if node:get_type() == "mRNA" then
          cdnaseq = node:extract_sequence("CDS", true, region_mapping)
          protseq = node:extract_and_translate_sequence("CDS", true,
                                                        region_mapping)
          --local startcodon = cdnaseq:sub(start_phase + 1, start_phase + 3):lower()
          --if embl_compliant and startcodon:match('^[tc]tg$') then
            -- The ENA expects non-M start codons (e.g. Ls) to be represented as
            -- Ms in the literal translation (according to their validator
            -- output). Emulating this behaviour here, diverging from the actual
            -- translation table.
          --  if protseq:sub(1,1):upper() == 'L' then
          --    protseq = 'M' .. protseq:sub(2)
          --  end
          --end
          if embl_compliant and cdnaseq:len() % 3 > 0 then
            -- The ENA expects translations for DNA sequences with a length
            -- which is not a multiple of three in a different way from the one
            -- produced by GenomeTools. While GenomeTools cuts off the
            -- remainder and does not translate the partial codon, the ENA
            -- validator fills it up with Ns and translates it. This behaviour
            -- is emulated here to make the \translation value validate
            -- properly.
            cdnaseq = cdnaseq .. string.rep('n', 3 - (cdnaseq:len() % 3))
            if start_phase > 0 then
              cdnaseq = cdnaseq:sub(start_phase + 1)
            end
            protseq = gt.translate_dna(cdnaseq)
          end

          -- clip off stop codon
          if protseq:sub(-1,-1) == "*" then
            protseq = protseq:sub(1,-2)
          end
          io.write("FT                   /translation=\"" .. protseq .."\"\n")
          io.write("FT                   /codon_start=" .. start_phase + 1 .. "\n")
        end
        io.write("FT                   /transl_table=1\n")
        -- orthologs
        local nof_orths = 0
        if pp and pp:get_attribute("orthologous_to") and not embl_compliant then
          for _,v in ipairs(split(pp:get_attribute("orthologous_to"), ",")) do
            io.write("FT                   /ortholog=\"" ..
                            v .. " " .. v .. ";program=OrthoMCL;rank=0\"\n")
            nof_orths = nof_orths + 1
          end
        end
        -- assign colours
        if node:get_type() == "mRNA" and not embl_compliant then
          if pp and pp:get_attribute("product") then
            local prod = pp:get_attribute("product")
            if  prod ~= "term%3Dhypothetical protein" then
              if nof_orths > 0 or prod:match("conserved") then
                io.write("FT                   /colour=10\n")   -- orange: conserved
              else
                io.write("FT                   /colour=7\n")    -- yellow: assigned Pfam
              end
            else
              if coding_length < 500 then
                io.write("FT                   /colour=6\n")    -- dark pink: short unlikely
              else
                io.write("FT                   /colour=8\n")    -- light green: hypothetical
              end
            end
          end
        elseif node:get_type() == "pseudogenic_transcript" and not embl_compliant then
          io.write("FT                   /colour=13\n")     -- pseudogene
        end
        -- add name
        local name = fn:get_attribute("Name")
        if name and not embl_compliant then
          io.write("FT                   /primary_name=\"".. name .. "\"\n")
        end
        -- GO terms
        local gaf_relevant_id = fn:get_attribute("ID")
        if options.gaf_level == 'transcript' then
          gaf_relevant_id = node:get_attribute("ID")
        end
        if gaf_relevant_id then
          if self.gaf and self.gaf[gaf_relevant_id] and not embl_compliant then
            for _,v in ipairs(self.gaf[gaf_relevant_id]) do
              if self.gos[v.goid] then
                io.write("FT                   /GO=\"aspect=" .. v.aspect ..
                                                    ";GOid=" .. v.goid ..
                                                    ";term=" .. self.gos[v.goid])
                if v.withfrom and string.len(v.withfrom) > 0 then
                  io.write(";with=" .. v.withfrom)
                end
                if v.evidence then
                  io.write(";evidence=" .. v.evidence)
                end
                io.write("\"\n")
              else
                io.stderr:write('no GO entry found for ' .. v.goid ..
                                ' -- possibly outdated GO OBO file\n')
              end
            end
          end
        end
        -- add UTRs
        if print_generic(node, 'five_prime_UTR', "5'UTR", self.adjust) then
          self:output_gene_qualifier(geneid, genesym)
          format_embl_attrib(node , "ID", "locus_tag",
            function (s)
              return split(s,':')[1]
            end)
          if is_multi_transcript_gene(fn) and not embl_compliant then
            io.write("FT                   /shared_id=\"" .. geneid .. "\"\n")
          end
        end
        if print_generic(node, 'three_prime_UTR', "3'UTR", self.adjust) then
          self:output_gene_qualifier(geneid, genesym)
          format_embl_attrib(node , "ID", "locus_tag",
            function (s)
              return split(s,':')[1]
            end)
          if is_multi_transcript_gene(fn) and not embl_compliant then
            io.write("FT                   /shared_id=\"" .. geneid .. "\"\n")
          end
        end
      elseif node:get_type() == "tRNA" then
        io.write("FT   tRNA            ")
        if node:get_strand() == "-" then
          io.write("complement(")
        end
        io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
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
          io.write("FT                   /gene=\"" .. fn:get_attribute("ID") .. "\"\n")
        end
        format_embl_attrib(node , "ID", "locus_tag", nil)
        output_note()
      elseif string.match(node:get_type(), "snRNA") or string.match(node:get_type(), "snoRNA") then
        io.write("FT   ncRNA            ")
        if node:get_strand() == "-" then
          io.write("complement(")
        end
        io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
        if node:get_strand() == "-" then
          io.write(")")
        end
        io.write("\n")
        io.write("FT                   /ncRNA_class=\"" .. node:get_type() .. "\"\n")
        io.write("FT                   /gene=\"" .. fn:get_attribute("ID") .. "\"\n")
        format_embl_attrib(node , "ID", "locus_tag", nil)
        output_note()
      elseif string.match(node:get_type(), "ncRNA") then
        io.write("FT   ncRNA            ")
        if node:get_strand() == "-" then
          io.write("complement(")
        end
        io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
        if node:get_strand() == "-" then
          io.write(")")
        end
        io.write("\n")
        io.write("FT                   /ncRNA_class=\"other\"\n")
        io.write("FT                   /gene=\"" .. fn:get_attribute("ID") .. "\"\n")
        format_embl_attrib(node , "ID", "locus_tag", nil)
        output_note()
      elseif string.match(node:get_type(), "rRNA") then
        io.write("FT   rRNA            ")
        if node:get_strand() == "-" then
          io.write("complement(")
        end
        io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
        if node:get_strand() == "-" then
          io.write(")")
        end
        io.write("\n")
        io.write("FT                   /product=\"" .. node:get_type() .. "\"\n")
        io.write("FT                   /gene=\"" .. fn:get_attribute("ID") .. "\"\n")
        format_embl_attrib(node , "ID", "locus_tag", nil)
        output_note()
      elseif string.match(node:get_type(), "gap") then
        if node:get_range():get_start() >= offset
            and node:get_range():get_end() <= collect_vis.seqs[node:get_seqid()]:len() then
          io.write("FT   assembly_gap    ")
          io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
          io.write("\n")
          io.write("FT                   /estimated_length=" .. node:get_range():length() .. "\n")
          if node:get_attribute('gap_type') then
            io.write("FT                   /gap_type=\"" .. node:get_attribute('gap_type') .. "\"\n")
            if node:get_attribute('gap_type') == 'within scaffold' then
              io.write("FT                   /linkage_evidence=\"align genus\"\n")
            end
          else
            io.write("FT                   /gap_type=\"unknown\"\n")
          end
        end
      elseif string.match(node:get_type(), "contig") and contig_as_source then
        io.write("FT   source          ")
        io.write(node:get_range():get_start() - offset .. ".." .. node:get_range():get_end() - offset)
        io.write("\n")
        io.write("FT                   /organism=\"" .. organismname .. "\"\n")
        io.write("FT                   /mol_type=\"unassigned DNA\"\n")
      end
    end
  else
    if not options.only_on_seq then
      error("sequence " .. tostring(fn:get_seqid()) .. " for " .. tostring(fn) .. " cannot be found in input")
      os.exit(1)
    end
  end
  return 0
end
function embl_vis:visit_region(rn)
  if collect_vis.seqs[rn:get_seqid()] then
    if self.last_seqid ~= rn:get_seqid() then
      output_embl_header(rn)
      self.last_seqid = rn:get_seqid()
    end
  end
end

function clip_flanking_ns(seqs)
  local adjust = {}
  for k,v in pairs(seqs) do
    local leading = 0
    local trailing = 0
    local m = v:match("^[Nn]+")
    if m then
      io.stderr:write("sequence " .. k .. " has " .. m:len() .. " leading Ns, trimming...\n")
      leading = m:len()
    end
    if leading == v:len() then
      io.stderr:write("sequence " .. k .. " is empty, removing it...\n")
      seqs[k] = nil
    else
      m = v:match("[Nn]+$")
      if m then
        io.stderr:write("sequence " .. k .. " has " .. m:len() .. " trailing Ns, trimming...\n")
        trailing = m:len()
      end
      if leading > 0 or trailing > 0 then
        seqs[k] = v:sub(1+leading, v:len()-trailing)
      end
    end
    adjust[k] = leading
  end
  return adjust, seqs
end

-- load GAF
gaf_store = {}
if gaffile then
  for l in io.lines(gaffile) do
    if l:sub(1,1) ~= '!' then
      local db,dbid,dbobj,qual,goid,dbref,evidence,withfrom,aspect,dbobjname,
        dbobjsyn,dbobjtype,taxon,data,assignedby = unpack(split(l, '\t'))
      if not gaf_store[dbid] then
        gaf_store[dbid] = {}
      end
      table.insert(gaf_store[dbid], {db=db, dbid=dbid, qual=qual, goid=goid,
                                 dbref=dbref, evidence=evidence, withfrom=withfrom,
                                 aspect=aspect, dbobjname=dbobjname,
                                 dbobjsyn=dbobjsyn, dbobjtype=dbobjtype,
                                 taxon=taxon, data=data, assignedby=assignedby})
    end
  end
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
visitor_stream.instream = gt.gff3_in_stream_new_sorted(gff3file)
visitor_stream.vis = collect_vis
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end

_, seqs = get_fasta_nosep(seqfile)
local adjust = nil
if options.trim_flanking_n then
  adjust, seqs = clip_flanking_ns(seqs)
end
collect_vis.seqs = {}
collect_vis.lengths = {}
for k,v in pairs(seqs) do
  collect_vis.seqs[k] = v
  collect_vis.lengths[k] = v:len()
end

-- output EMBL code as we traverse the GFF
visitor_stream.instream = gt.gff3_in_stream_new_sorted(gff3file)
visitor_stream.vis = embl_vis
embl_vis.obo = gos
embl_vis.gaf = gaf_store
embl_vis.seen_gene_names = {}
embl_vis.adjust = adjust
local gn = visitor_stream:next_tree()
while (gn) do
  gn = visitor_stream:next_tree()
end
-- output last seq
if options.add_onebase_gap then
  output_single_base_gaps(collect_vis.seqs[embl_vis.last_seqid])
end
format_embl_sequence(collect_vis.seqs[embl_vis.last_seqid])
io.write("//\n")
io.output():close()

