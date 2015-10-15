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

SimpleChainer = {}
SimpleChainer.__index = SimpleChainer

function precedes(a,b)
  return a:get_range():get_end() < b:get_range():get_start()
end

function SimpleChainer.new(arr, weight_func, data)
	local sc = {}
	local mc = {}
	setmetatable(sc, SimpleChainer)
	for _,v in ipairs(arr) do
		local nval = {item = v, weight = weight_func(v, data)}
		table.insert(mc, nval)
	end
	sc.mc = mc
	return sc
end

function SimpleChainer:chain()
	local overallmaxscore = 0
	local bestmatch = nil
	local score = {}
	local prec = {}
	-- dynamic programming
	for j = 1,#self.mc do
	  local maxscore = self.mc[j].weight
	  for i = 1,j do
	    if precedes(self.mc[i].item, self.mc[j].item)
	      and maxscore < self.mc[j].weight + score[self.mc[i]] then
	      maxscore = self.mc[j].weight + score[self.mc[i]]
	      prec[self.mc[j]] = self.mc[i]
	    end
	  end
	  score[self.mc[j]] = maxscore
	  if overallmaxscore < maxscore then
	    overallmaxscore = maxscore
	    bestmatch = self.mc[j]
	  end
	end
	-- traceback
	local nextbest = bestmatch
    local bestchain = {}
	while nextbest do
		table.insert(bestchain, 1, nextbest.item)
	    nextbest = prec[nextbest]
	end
	return bestchain
end
