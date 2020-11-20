-- -*- mode: Lua; mode: fold -*- 
-- init/grapham_init.lua
-- 
-- Copyright (c) Matti Vihola 2009-2013
-- 
-- This file is part of Grapham.
-- 
-- Grapham is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
-- 
-- Grapham is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
-- 
-- You should have received a copy of the GNU General Public License
-- along with Grapham.  If not, see <http://www.gnu.org/licenses/>.

local _io = io
local _string = string

-- Make math library functions visible also as globals
for name,value in pairs(math) do
   _G[name] = value
end

const = {}
model = {}
para = {}

-- REPEAT_BLOCK  Helper for creating repeated blocks in the model {{{
--
-- Usage: repeat_block(blocks, N)
--        repeat_block(blocks, data_1, ...)
-- 
-- In: 
--   blocks -- A table of strings of variable names
--   N      -- How many times the block is repeated
--   data_i -- Alternatively, one can give directly the data
--             of the i:th variable in "blocks".
--
-- Out: none
-- 
-- The function repeats a given block N times. Then, if called by 
-- data, the corresponding block elements are instantiated with
-- the data.
function repeat_block(block, ...)
   if arg.n < 1 then
      error("Error: repeat_block: not enough arugments.")
   end
   local k, j, var, vari, varj, inst
   local N

   local function length(arg)
      local Nargn
      if type(arg) == "table" then
         Nargn = #arg
      elseif type(arg) == "userdata" then
         Nargn = arg:size()
      end
      return Nargn
   end
   
   if type(arg[1]) == "number" then
      N = arg[1]; inst = false
   else 
      N = length(arg[1]); inst = true
   end
   local function add_suffix(tb, strs, n)
      local k, j, m, str, found
      local tb2 = {}
      if tb ~= nil then
         for k=1,#tb do
            str = tb[k]
            tb2[k] = str
            m = 0
            found = true
            -- Count the occurrences of trailing "-"'s
            while found do
               found = false
               str = str:gsub("-$", function() 
                                       m = m+1; found = true; 
                                       return ""
                                    end)
            end
            for j=1,#strs do
               if str == strs[j] then tb2[k] = str .. n-m end
            end
         end
      end
      return tb2
   end
   
   local function copy_table(tb)
      local tb2 = {};
      for i,v in pairs(tb) do
         tb2[i] = v;
      end
      return tb2
   end
   
   for k=2,arg.n do
      if length(arg[k]) ~= N then
         error("Error: repeat_block: invalid arguments.")
      end
   end
   
   if data == nil then
      data = {}
   end
   
   for k=1,#block do
      var = block[k];
      for j=1,N do
         varj = var .. j
         if model[var] == nil then
            error("Error: repeat_block: cannot find variable '" ..
                  var .. "' from model.")
         end
         if model[varj] ~= nil then
            error("Error: repeat_block: '"..varj.."' already in model.")
         end
         if inst and k <= arg.n then
            data[varj] = arg[k][j]
         end
         model[varj] = copy_table(model[var])
         model[varj].parents = add_suffix(model[varj].parents, block, j)
      end
      -- Remove element
      model[var] = nil
   end
   for var in pairs(model) do   
      if model[var].parents ~= nil then
         for k=1,#model[var].parents do
            for j=1,#block do
               if model[var].parents[k] == block[j] then
                  error("Error: repeat_block: cannot repeat a block"
                        .." with children outside the block.")
               end
            end
         end
      end
   end
end
-- }}}

-- READ_CSV  Helper for reading a CSV (comma separated values) data file {{{
-- 
-- Usage: vars, data = read_csv(fname)
-- 
-- In:
--   fname -- The name of the CSV file. Note: the path is given relative
--            to where grapham is executed.
-- 
-- Out:
--   vars -- Names of the variables read from the CSV header
--           (Empty table, if no variable names have been read)
--   data -- A table with each element having the data for the corresponding
--           variable (i.e. column in the CSV data)
-- 
-- Reads a CSV (comma separated values) file in. Note that the function
-- is not perfect, i.e. it does *NOT* handle well escaped elements.
function read_csv(fname)
   local var, line, data, vars, k, head
   -- Note: this pattern _does not_ handle escaping at all!
   local csv_pat = '([^,]+)(,*)'

   local fh = assert(_io.open(fname, "r"))
   vars = {}; data = {}
   line =  fh:read("*line")
   head = true
   for var in _string.gfind(line, csv_pat) do
      var = _string.gsub(var, '"', '')
      if tonumber(var) ~= nil then 
         head = false
      else
         vars[#vars+1] = var
      end
   end
   if head and #vars ~= 0 then
      for k=1,#vars do data[k] = {} end
   else
      k=1
      for var in _string.gfind(line, csv_pat) do
         data[k] = {tonumber(var)}
         k = k+1
      end
   end
   line = fh:read("*line")
   
   while line ~= nil do
      k=1
      for var in _string.gfind(line, csv_pat) do
         data[k][#data[k]+1] = tonumber(var)
         k = k+1
      end
      line = fh:read("*line")
   end
   fh:close()
   return vars, data
end
-- }}}

-- REQUIRES_NUMLUA  Helper for checking for Numeric Lua {{{
function requires_numlua()
   if matrix == nil then
      error("This model requires Numeric Lua!")
   end
end
-- }}}
