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

function gff3_encode(s)
  s = string.gsub(s, "[\t\n\r;=%&,]", function (c)
        return string.format("%%%02X", string.byte(c))
      end)
  return s
end

function gff3_decode(s)
  if not s then
    return s
  end
  s = string.gsub(s, "%%([0-9a-fA-F][1-9a-fA-F])", function (n)
        return string.char(tonumber("0x" .. n))
      end)
  return s
end

function gff3_explode(tab)
  local ret = {}
  for _,item in ipairs(tab) do
    local _tmp = {}
    for k,v in pairs(item) do
      table.insert(_tmp, k .. "=" .. v)
    end
    table.insert(ret, gff3_encode(table.concat(_tmp, ";")))
  end
  return table.concat(ret, ",")
end

function gff3_extract_structure(str)
  ret = {}
  for _,v in ipairs(split(str, ", ?")) do
    res = {}
    v = gff3_decode(v)
    for _,pair in ipairs(split(v, ";")) do
      if string.len(pair) > 0 then
        key, value = unpack(split(pair, "="))
        res[key] = value
      end
    end
    table.insert(ret, res)
  end
  return ret
end

function file_exists(name)
   local f=io.open(name,"r")
   if f~=nil then io.close(f) return true else return false end
end

function split(str, pat)
   local t = {}
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
        table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end

function string:split(sep)
  local sep, fields = sep or ":", {}
  local pattern = string.format("([^%s]+)", sep)
  self:gsub(pattern, function(c) fields[#fields+1] = c end)
  return fields
end

function is_experimental(evidence)
   exp_strings={"Experimen", "Direct Assa", "Physical Interactio",
                "Mutant Phenotyp", "Genetic Interactio", "Expression Patter",
                "EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}
   for _,str in ipairs(exp_strings) do
     if string.match(tostring(evidence), str) then
      return true
     end
   end
   return false
end

function is_curated(evidence)
   cur_strings={"Structural Similarit", "Sequence Ortholog",
                "Sequence Alignmen",
                "Sequence Mode", "Genomic Contex", "aspect of Ancesto",
                "aspect of Descendan", "Key Residue", "Rapid Divergenc",
                "Reviewed Computational Analysi", "Author Statemen", "Curato"}
   for _,str in ipairs(cur_strings) do
     if string.match(tostring(evidence), str) then
      return true
     end
   end
   return false
end

function map(func, array)
  local new_array = {}
  for i,v in ipairs(array) do
    new_array[i] = func(v)
  end
  return new_array
end

function table_count(tt, value)
  local count = 0
  for ii,xx in pairs(tt) do
    if xx == value then count = count + 1 end
  end
  return count
end

function table_unique(tt)
  local newtable = {}
  for ii,xx in ipairs(tt) do
    if table_count(newtable, xx) == 0 then
      newtable[#newtable+1] = xx
    end
  end
  return newtable
end

function table_concat(t1,t2)
    for i=1,#t2 do
        t1[#t1+1] = t2[i]
    end
    return t1
end

function revcomp(str)
  local comp_table = {}
  comp_table[string.byte('A')] = 'T';
  comp_table[string.byte('T')] = 'A';
  comp_table[string.byte('C')] = 'G';
  comp_table[string.byte('G')] = 'C';
  comp_table[string.byte('a')] = 't';
  comp_table[string.byte('t')] = 'a';
  comp_table[string.byte('c')] = 'g';
  comp_table[string.byte('g')] = 'c';
  comp_table[string.byte('n')] = 'n';
  comp_table[string.byte('N')] = 'N';
  local x, s = {}, str:reverse()
  for i = 1, #s do x[i] = comp_table[s:byte(i)] end
  return table.concat(x)
end

-- remove trailing and leading whitespace from string.
function trim(s)
  -- from PiL2 20.4
  return (s:gsub("^%s*(.-)%s*$", "%1"))
end

-- remove leading whitespace from string.
function ltrim(s)
  return (s:gsub("^%s*", ""))
end

-- remove trailing whitespace from string.
function rtrim(s)
  local n = #s
  while n > 0 and s:find("^%s", n) do n = n - 1 end
  return s:sub(1, n)
end

function get_fasta(filename, sep)
  local keys = {}
  local seqs = {}
  local cur_hdr = nil
  local cur_seqs = {}
  if not file_exists(filename) then
    error("file " .. filename .. " does not exist")
  end
  for l in io.lines(filename) do
    hdr = l:match(">(.*)")
    if hdr then
      hdr = trim(split(hdr,"%s+")[1])
      table.insert(keys, hdr)
      if #cur_seqs > 0 and cur_hdr then
        if not seqs[cur_hdr] then
          seqs[cur_hdr] = table.concat(cur_seqs, sep)
        end
      end
      cur_hdr = hdr
      cur_seqs = {}
    else
      table.insert(cur_seqs, l)
    end
  end
  if cur_hdr and not seqs[cur_hdr] then
    seqs[cur_hdr] = table.concat(cur_seqs, sep)
  end
  return keys, seqs
end

function get_fasta_nosep(filename)
  return get_fasta(filename, "")
end

function gtf_lines(in_io)
  -- returns, for each GTF line, a table with the 8 first column values as
  -- strings/numbers; attributes are given in a sub-table called 'attribs'
  assert(in_io)
  return function()
    local line = in_io:read()
    while (line and line:sub(1,1) == '#') do
      line = in_io:read()
    end
    if line then
      local larr = split(line, "\t")
      if #larr ~= 9 then
        error('could not parse 9 columns from input line "' .. line .. '"')
      end
      local seqid, source, type, from, to, score, strand, phase,
            attr = unpack(larr)
      local feat = {seqid = seqid, source = source, type = type,
                    from = tonumber(from), to = tonumber(to),
                    score = tonumber(score), strand = strand, phase = phase}
      local attr_table = {}
      attr:gsub('([%a_]+) "([^"]+)";', function(k,v)
        attr_table[k] = v
      end)
      feat.attribs = attr_table
      return feat
    else
      return nil
    end
  end
end

function visitor_stream_new(instream, vis)
  -- make simple visitor stream, just applies given visitor to every node
  local visitor_stream = gt.custom_stream_new_unsorted()
  visitor_stream.instream = instream
  visitor_stream.vis = vis
  function visitor_stream:next_tree()
    local node = self.instream:next_tree()
    if node then
      node:accept(self.vis)
    end
    return node
  end
  return visitor_stream
end

-- completely clones a CC
local _clone_id = 0
function deep_copy(root, copy, table)
  if not copy then
    local ntype = root:get_type()
    if table and table[ntype] then
      ntype = table[ntype]
    end
    copy = gt.feature_node_new(root:get_seqid(), ntype,
                                   root:get_range():get_start(),
                                   root:get_range():get_end(),
                                   root:get_strand())
    for k,v in root:attribute_pairs() do
      if k == "ID" then
        copy:set_attribute(k, v .. ":" .. _clone_id)
        _clone_id = _clone_id + 1
      elseif k ~= "Parent" then
        copy:set_attribute(k, v)
      end
    end
  end
  for c in root:direct_children() do
    local ntype = c:get_type()
    if table and table[ntype] then
      ntype = table[ntype]
    end
    new_node = gt.feature_node_new(c:get_seqid(), ntype,
                                   c:get_range():get_start(),
                                   c:get_range():get_end(),
                                   c:get_strand())
    for k,v in c:attribute_pairs() do
      if k == "ID" then
        new_node:set_attribute(k, v .. ":" .. _clone_id)
        _clone_id= _clone_id + 1
      elseif k ~= "Parent" then
        new_node:set_attribute(k, v)
      end
    end
    copy:add_child(new_node)
    deep_copy(c, new_node, table)
  end
  return copy
end

function trim_ns(inseq)
  inseq = inseq:gsub("^[Nn]+", ""):gsub("[Nn]+$", "")
  return inseq
end

function print_max_width(str, ioo, width)
  local i = 1
  while (i + width - 1) < str:len() do
    ioo:write(str:sub(i, i + width - 1))
    ioo:write("\n")
    i = i + width
  end
  ioo:write(str:sub(i))
  ioo:write("\n")
end

function clone_cc(incc)
  return deep_copy(incc, nil)
end

function make_gene_name(str)
  if string.match(str, "%.%.") then
    return split(str, "%.%.")[1]
  elseif string.match(str, ":") then
    return split(str, ":")[1]
  else
    return str
  end
end

function get_clusters(filename)
  local clindex = {}
  local clusters = {}
  if not file_exists(filename) then
    error("file " .. filename .. " can not be loaded")
  end
  for l in io.lines(filename) do
    local name, nofgenes, noftaxa, members = string.match(l,
                             "^([^(]+)%((%d+) genes,(%d+) taxa%):%s+(.+)")
    if not name or not members then
      error("could not parse cluster or members from line '" .. l .. "'")
    end
    clusters[name] = {}
    local n = 0
    local thisclust = {name = name, members = {}, specidx = {}, species = {}}
    for mname, mspecies in members:gmatch("([^(]+)%(([^)]+)%)%s?") do
      table.insert(thisclust.members, {mname, mspecies})
      thisclust.specidx[mspecies] = true
      table.insert(thisclust.species, mspecies)
      clindex[mname] = thisclust
      clindex[make_gene_name(mname)] = thisclust
      n = n + 1
    end
    thisclust.species = table_unique(thisclust.species)
    if n ~= tonumber(nofgenes) then
      error("parsed " .. n .. " members from cluster " .. name
              .. ", but expected " .. nofgenes)
    end
    clusters[name] = thisclust
  end
  return clindex, clusters
end

function print_r ( t )
    local print_r_cache={}
    local function sub_print_r(t,indent)
        if (print_r_cache[tostring(t)]) then
            print(indent.."*"..tostring(t))
        else
            print_r_cache[tostring(t)]=true
            if (type(t)=="table") then
                for pos,val in pairs(t) do
                    if (type(val)=="table") then
                        print(indent.."["..pos.."] => "..tostring(t).." {")
                        sub_print_r(val,indent..string.rep(" ",string.len(pos)+8))
                        print(indent..string.rep(" ",string.len(pos)+6).."}")
                    elseif (type(val)=="string") then
                        print(indent.."["..pos..'] => "'..val..'"')
                    else
                        print(indent.."["..pos.."] => "..tostring(val))
                    end
                end
            else
                print(indent..tostring(t))
            end
        end
    end
    if (type(t)=="table") then
        print(tostring(t).." {")
        sub_print_r(t,"  ")
        print("}")
    else
        sub_print_r(t,"  ")
    end
    print()
end

function overlap_stream_new(instream, types, func)
  local overlap_stream = gt.custom_stream_new_unsorted()
  overlap_stream.instream = instream
  overlap_stream.func = func
  overlap_stream.outqueue = {}
  overlap_stream.curr_gene_set = {}
  overlap_stream.curr_rng = nil
  overlap_stream.last_seqid = nil
  overlap_stream.inkeys = nil
  overlap_stream.end_reached = false

  -- register relevant types as keys for set queries later
  if types then
    for _,v in ipairs(types) do
      if not self.inkeys then
        self.inkeys = {}
      end
      overlap_stream.inkeys[v] = true
    end
  end

  -- handler function for each cluster
  function overlap_stream:process_current_cluster()
    bestset = self.func(self.curr_gene_set)
    for _,v in ipairs(bestset) do
      table.insert(self.outqueue, v)
    end
    self.curr_gene_set = {}
  end

  function overlap_stream:next_tree()
    local complete_cluster = false
    local mygn = nil
    if #self.outqueue > 0  then
      return table.remove(self.outqueue, 1)
    else
      if self.end_reached then
        return nil
      else
        complete_cluster = false
      end
    end
    while not complete_cluster do
      mygn = self.instream:next_tree()
      if mygn then
        rval, err = pcall(GenomeTools_genome_node.get_type, mygn)
        if rval then
          local fn = mygn
          local new_rng = mygn:get_range()
          if not self.inkeys or self.inkeys[fn:get_type()] then
            if #self.curr_gene_set == 0 then
              table.insert(self.curr_gene_set, fn)
              self.curr_rng = new_rng
            else
              if self.last_seqid == fn:get_seqid()
                  and self.curr_rng:overlap(new_rng) then
                table.insert(self.curr_gene_set, fn)
                self.curr_rng = self.curr_rng:join(new_rng)
              else
                -- no more overlap
                self:process_current_cluster()
                table.insert(self.curr_gene_set, fn)
                self.curr_rng = new_rng
                if #self.outqueue > 0  then
                  outgn = table.remove(self.outqueue, 1)
                  complete_cluster = true
                end
              end
            end
            self.last_seqid = mygn:get_seqid()
          end
        else
          -- no feature node
          self:process_current_cluster()
          table.insert(self.outqueue, fn)
          if #self.outqueue > 0  then
            outgn = table.remove(self.outqueue, 1)
            complete_cluster = true
          end
        end
      else
        -- end of annotation
        outgn = mygn
        if #self.curr_gene_set > 0 then
          self:process_current_cluster()
          table.insert(self.outqueue, fn)
        end
        outgn = table.remove(self.outqueue, 1)
        complete_cluster = true
        self.end_reached = true
      end
    end
    return outgn
  end
  return overlap_stream
end

function infernal_in_stream_new(myio)
  local assign_type = function (name)
    if name:match("rRNA") then
      return 'rRNA'
    elseif name:match("[Ss][Nn][Oo][RT]") then
      return 'snoRNA'
    elseif name:match("U%d") then
      return 'snRNA'
    end
    -- catch all
    return 'ncRNA'
  end

  local feats = {}
  local i = 1
  while true do
    local line = myio:read()
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
      if qry ~= 'tRNA' then
        local newgene = gt.feature_node_new(seqid, "gene", sfrom, sto, strand)
        newgene:set_score(tonumber(score))
        newgene:set_attribute("ID", "ncRNA" .. tostring(i))
        newgene:set_source("GenomeTools")

        local toptype = assign_type(qry)
        local newrna = gt.feature_node_new(seqid, toptype, sfrom, sto, strand)
        newrna:set_score(tonumber(score))
        newrna:set_attribute("ID", "ncRNA" .. tostring(i) .. ":" .. gff3_encode(toptype))
        newrna:set_attribute("Name", gff3_encode(qry))
        newrna:set_attribute("gc", gff3_encode(gc))
        newrna:set_attribute("evalue", gff3_encode(evalue))
        newrna:set_attribute("score", gff3_encode(score))
        newrna:set_attribute("model_name", gff3_encode(qry))
        newrna:set_attribute("model_acc", gff3_encode(qryacc))
        newrna:set_attribute("model_range", tostring(mfrom) .. "-" .. tostring(mto))
        newrna:set_source("INFERNAL")
        newgene:add_child(newrna)
        table.insert(feats, newgene)
        i = i + 1
      end
    end
  end

  local s = gt.custom_stream_new_unsorted()
  s.feats = feats
  table.sort(s.feats, function(a,b)
    -- XXX: this should best be done in gt by wrapping genome_node_compare()
    if a:get_seqid() < b:get_seqid() then
      return true
    elseif  a:get_seqid() > b:get_seqid() then
      return false
    end
    local arng = a:get_range()
    local brng = b:get_range()
    if arng:get_start() < brng:get_start() then
      return true
    elseif arng:get_start() > brng:get_start() then
      return false
    else
      assert(arng:get_start() == brng:get_start())
      if arng:get_end() < brng:get_end() then
        return true
      else
        return false
      end
    end
  end)
  function s:next_tree()
    if #self.feats > 0 then
      local t = table.remove(self.feats, 1)
      return t
    else
      return nil
    end
  end
  return s
end