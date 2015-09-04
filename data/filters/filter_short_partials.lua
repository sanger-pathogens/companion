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

name        = "Small partial gene filter"
author      = "Sascha Steinbiss"
version     = "1.0"
email       = "ss34@sanger.ac.uk"
short_descr = "Remove partial genes shorter than a given amount of bases."
description = "Selects only features which are not partial genes "
              .. "(Start_range/End_range set or {three,five}EndPartial set) "
              .. "shorter than a defined length."

FILTER_LENGTH = 300

function filter(fn)
  if fn:get_type() == "gene" then
    if fn:get_attribute("Start_range") or fn:get_attribute("End_range") or
      fn:get_attribute("fiveEndPartial") or fn:get_attribute("threeEndPartial") then
      if fn:get_range():length() < FILTER_LENGTH then
        return true
      end
    end
  end
  return false
end
