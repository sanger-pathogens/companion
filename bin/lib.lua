--[[
  Copyright (c) 2014-2015 Sascha Steinbiss <ss34@sanger.ac.uk>
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
    error("file " .. filename .. " can not be loaded")
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
