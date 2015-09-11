#!/usr/bin/env ruby
#
# Author: Sascha Steinbiss <ss34@sanger.ac.uk>
# Copyright (c) 2014 Genome Research Ltd
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

if $0 == __FILE__
  $:<< "."            # favor the local stest version
  require 'stest'
  at_exit do
    OnError do exit 1 end
  end
end

# set some global variables
if $arguments["path"] then
  $path=File.join($arguments["path"], "")
else
  $path=""
end

if $arguments["bin"] then
  $bin=File.join($arguments["bin"], "")
else
  $bin=""
end

if $arguments["cur"] then
  $cur=$arguments["cur"]
else
  $cur=File.join(Dir.pwd, "..", "")
end

$systemname=`uname -s`
$systemname.chomp!

if $arguments["seed"] then
  $SEED = $arguments["seed"]
else
  $SEED = rand(2**31)
end

# define helper function
def run_test(str, opts = {})
  run("env G_DEBUG=gc-friendly G_SLICE=always-malloc #{$path}#{str}", opts)
end

def with_environment(variables={})
  if block_given?
    old_values = variables.map{ |k,v| [k,ENV[k]] }
    begin
       variables.each{ |k,v| ENV[k] = v }
       result = yield
    ensure
      old_values.each{ |k,v| ENV[k] = v }
    end
    result
  else
    variables.each{ |k,v| ENV[k] = v }
  end
end

$rootdir=File.join(Dir.pwd , "..", "..")
$bindir=File.join(Dir.pwd , "..", "..", "bin")
$testdata=File.join(Dir.pwd , "..", "testdata", "")

require 'aragorn_include'
require 'misc_include'
require 'ncrna_include'

#we now have all tests in $testsuite.

if $arguments["threads"] then
  $testsuite.nof_threads = $arguments["threads"].to_i
end

#start tests
$testsuite.run
