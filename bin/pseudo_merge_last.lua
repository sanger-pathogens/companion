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
  io.stderr:write("Extracts potential pseudogenes from LAST output as GFF3 "
                    .. "features, annotated with frameshifts and internal "
                    .. "stop codons.\n")
  io.stderr:write(string.format("Usage: %s <lastout> <sequence>\n", arg[0]))
  os.exit(1)
end

if #arg < 2 then
  usage()
end

P_MAXGAP = 100
C_MAXGAP = 500

rm = gt.region_mapping_new_seqfile_matchdescstart(arg[2])

package.path = gt.script_dir .. "/./?.lua;" .. package.path
require("lib")

function precedes(l1, l2, cgap, pgap)
  return ((l2.crange:get_start() >= l1.crange:get_end()
            and l2.crange:get_start() - l1.crange:get_end() + 1 < cgap) and
         (l2.prange:get_start() >= l1.prange:get_end()
            and l2.prange:get_start() - l1.prange:get_end() + 1 < pgap))
end

gff3v = gt.gff3_visitor_new()
pno = 0

function round(num, idp)
  local mult = 10^(idp or 0)
  return math.floor(num * mult + 0.5) / mult
end

function process_set_for_query(qry, set)
  local locs = {}
  -- sort hits by chromosome
  for _,i in ipairs(set) do
    if not locs[i.chr] then
      locs[i.chr] = {}
    end
    table.insert(locs[i.chr], i)
  end
  -- join into superhits (like PseudoPipe)
  for k,v in pairs(locs) do
    local last_hit = nil
    local cur_superhit = {}
    local superhits = {cur_superhit}
    table.sort(v, function(a,b)
      return a.crange:get_start() < b.crange:get_start()
    end)
    for _,l in ipairs(v) do
      if #cur_superhit == 0 then
        table.insert(cur_superhit, l)
      else
        if last_hit and precedes(last_hit, l, C_MAXGAP, P_MAXGAP) then
          table.insert(cur_superhit, l)
        else
          cur_superhit = {l}
          table.insert(superhits, cur_superhit)
        end
        last_hit = l
      end
    end
    -- output superhits as potential pseudogenes in annotated GFF3
    for _,h in ipairs(superhits) do
      local prot_len = h[1].prange:get_end() - h[#h].prange:get_start() + 1
      if prot_len/h[1].plength >= 0.8 then
        local frameshift = false
        local internal_stop = false
        -- parent feature
        local pp = gt.feature_node_new(k, "pseudogene", h[1].crange:get_start(),
                                       h[#h].crange:get_end(), h[1].cstrand)
        local pt = gt.feature_node_new(k, "pseudogenic_transcript",
                                       h[1].crange:get_start(),
                                       h[#h].crange:get_end(), h[1].cstrand)
        pp:add_attribute("ID", "pseudo_" .. qry .. tostring(pno))
        pt:add_attribute("ID", "pseudo_trans_" .. qry .. tostring(pno))
        pno = pno + 1
        pp:add_child(pt)
        -- gather frameshifts from LAST blocks and create subfeatures
        for _,hit in ipairs(h) do
          local lengths = {}
          local start = hit.crange:get_start()
          local cur_length = 0
          for _,b in ipairs(split(hit.blocks,',')) do
            -- if there is a colon in the block, we have an indel with base
            -- level gap lengths -- these are frameshifted if not divisible
            -- by 3
            if b:match(':') then
              local pdis
              local ddis
              pdis, ddis = unpack(split(b, ':'))
              if tonumber(ddis) % 3 == 0 then
                cur_length = cur_length + tonumber(ddis)
              else
                -- frameshift: start a new pseudogenic exon
                if cur_length > 3 then
                  local nnode = gt.feature_node_new(k, "pseudogenic_exon",
                                                    start,
                                                    start + cur_length - 1,
                                                    hit.cstrand)
                  pt:add_child(nnode)
                end
                start = start + cur_length + tonumber(ddis)
                cur_length = 0
                frameshift = true
              end
            else
              -- no indel, matches only -- these lengths are in amino acids,
              -- not bases
              cur_length = cur_length + 3*tonumber(b)
            end
          end
          -- finish last started pseudogenic exon
          local nnode = gt.feature_node_new(k, "pseudogenic_exon",
                                            start, start + cur_length - 1,
                                            hit.cstrand)
          pt:add_child(nnode)
          -- reverse the exons if this is a minus strand feature
          if pp:get_strand() == '-' then
            local pt_end = pt:get_range():get_end()
            local pt_start = pt:get_range():get_start()
            for c in pt:children() do
              local new_end = pt_end - (c:get_range():get_start() - pt_start)
              local new_start = pt_end - (c:get_range():get_end() - pt_start)
              c:set_range(gt.range_new(new_start, new_end))
            end
          end
        end

        -- adjust overlapping pseudogenic exons
        local last_pe = nil
        local i = 0
        for pe in pt:children() do
          if pe:get_type() == 'pseudogenic_exon' then
            if last_pe then
              if last_pe:get_range():overlap(pe:get_range()) then
                if pt:get_strand() == '-' then
                  local newend = pe:get_range():get_end() - 3
                  if pe:get_range():get_start() <= newend then
                    local new_range = gt.range_new(pe:get_range():get_start(), newend)
                    pe:set_range(new_range)
                  end
                else
                  local newstart = pe:get_range():get_start() + 3
                  if newstart <= pe:get_range():get_end() then
                    local new_range = gt.range_new(newstart, pe:get_range():get_end())
                    pe:set_range(new_range)
                  end
                end
              end
            end
            last_pe = pe
            i = i + 1
          end
        end

        -- determine internal stop
        local aaseq = pt:extract_and_translate_sequence("pseudogenic_exon",
                                                        true, rm)
        if aaseq:sub(1, -2):match("[*+#]") then
          internal_stop = true
        end
        pp:add_attribute("has_internal_stop", tostring(internal_stop))
        pp:add_attribute("has_frameshift", tostring(frameshift))
        -- check presence of start/stop codons
        if aaseq:sub(1,1):match("[Mm]") then
          pp:add_attribute("has_start", "true")
        end
        if aaseq:sub(-1):match("[*+#]") then
          pp:add_attribute("has_stop", "true")
        end
        -- add misc details
        pp:add_attribute("Target", qry .. " " .. h[1].prange:get_start()
                           .. " " .. h[#h].prange:get_end())
        pp:add_attribute("original_prot_length", h[1].plength)
        pp:set_score(h[1].score)
        -- output
        pp:accept(gff3v)
      end
    end
  end
end

function parse_last_out(lastoutfile)
  local lastchr = nil
  local matches_for_prot = {}
  local score = nil
  local prot = nil
  local chr = nil
  local prange = nil
  local crange = nil
  for l in io.lines(lastoutfile) do
    if l:sub(1,1) ~= '#' then
      score, prot, pstart, palnsize, pstrand, pseqsize, chr, cstart, calnsize,
        cstrand, cseqsize, blocks  = unpack(split(l, "%s+"))
        if not blocks then
          io.stderr:write("invalid/incomplete LAST line: '" .. l .."'\n")
          os.exit(1)
        end
      if cstrand == '-' then
        cstart = cseqsize - cstart - calnsize
      end
      local prange = gt.range_new(tonumber(pstart)+1,
                            tonumber(pstart) + tonumber(palnsize))
      local crange = gt.range_new(tonumber(cstart)+1,
                            tonumber(cstart) + tonumber(calnsize))
      if not matches_for_prot[prot] then
        matches_for_prot[prot] = {}
      end
      table.insert(matches_for_prot[prot], {chr = chr, crange = crange,
                                            prange = prange,
                                            prot = prot, pstrand = pstrand,
                                            cstrand = cstrand,
                                            plength = pseqsize,
                                            clength = cseqsize,
                                            score = score,
                                            blocks = blocks})
      if lastchr ~= chr then
        io.stderr:write(chr.."\n")
        lastchr = chr
      end
    end
  end
  return matches_for_prot
end

mfp = parse_last_out(arg[1])
for k,v in pairs(mfp) do
  process_set_for_query(k, v)
end