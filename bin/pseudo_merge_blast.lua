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
  io.stderr:write("Extracts best hits from tblastn output as GFF3 "
                    .. "'protein_match' features.\n")
  io.stderr:write(string.format("Usage: %s <blastout>\n" , arg[0]))
  os.exit(1)
end

if #arg < 1 then
  usage()
end

package.path = gt.script_dir .. "/./?.lua;" .. package.path
require("lib")

function make_ranges_for_hit(hit)
  return({qryrng = gt.range_new(tonumber(hit.qfrom), tonumber(hit.qto)),
          subjrng = gt.range_new(tonumber(hit.sfrom), tonumber(hit.sto)),
          subjstrand = hit.strand})
end

function process_set_for_query(qry, set)
  local locs = {}
  for _,i in ipairs(set) do
    -- location (chromosome)
    if not locs[i.subj] then
      locs[i.subj] = {}
    end
    -- pile up + join overlapping chromosomal ranges
    if #locs[i.subj] == 0 then
      table.insert(locs[i.subj], make_ranges_for_hit(i))
    else
      local joined = false
      local r1 = make_ranges_for_hit(i)
      for _,r2 in ipairs(locs[i.subj]) do
        if r1.subjrng:overlap(r2.subjrng) and not joined then
          r2.subjrng:join(r1.subjrng)
          joined = true
        end
      end
      if not joined then
        table.insert(locs[i.subj], make_ranges_for_hit(i))
      end
    end
  end
  -- final merge
  for k,v in pairs(locs) do
    for _,l in ipairs(v) do
      -- TODO
      print(k .. "\t" .. ".\tregion\t" .. tostring(l.subjrng:get_start())
             .. "\t" .. tostring(l.subjrng:get_end()) .. "\t.\t".. l.subjstrand .. "\t.\tName=pseudogene_"..qry)
    end
  end
end

function process_blast_out(blastoutfile)
  local lastqry = nil
  local curr_set = {}

  for l in io.lines(blastoutfile) do
    local strand = nil
    qry, subj, id, len, mism, gaps, qfrom,
                           qto, sfrom, sto, eval, score = unpack(split(l, "\t"))
    if qry ~= lastqry then
      if lastqry and #curr_set > 0 then
        process_set_for_query(lastqry, curr_set)
      end
      curr_set = {}
    else
      if tonumber(sfrom) > tonumber(sto) then
        from = tonumber(sto)
        to = tonumber(sfrom)
        strand = "-"
      else
        from = tonumber(sfrom)
        to = tonumber(sto)
        strand = "+"
      end
      if tonumber(id) > 30 and tonumber(eval) < 1e-4 then
        table.insert(curr_set, {subj = subj, id = id, len = len, mis = mism,
                                gaps = gaps, qfrom = qfrom, qto = qto,
                                sfrom = from , sto = to,
                                eval = eval, score = score,
                                strand = strand })
      end
    end
    lastqry = qry
  end
end

process_blast_out(arg[1])

--print("##gff-version\t3")