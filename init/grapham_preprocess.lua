-- -*- mode: Lua; mode: fold -*- 
-- init/grapham_preprocess.lua
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

-- This script "preprocesses" the user-supplied Grapham configuration
-- files. This includes checking the integrity, supplying the default 
-- values for several fields, and "normalising" the file for easier
-- reading into C.
if _PREPROCESS ~= nil then
   error("\nGlobal variable _PREPROCESS already defined.\n");
end

-- Do everything in a function to avoid messing things up.
function _PREPROCESS() 
   local MODEL_VNAME = "model"
   local PARA_VNAME = "para"
   local CONST_VNAME = "const"
   local INST_VNAME = "data"
   local FUNC_VNAME = "functional"
   
   local PARA_FIELDS = { --{{{
      NITER = "niter",
      NBURN = "nburn",
      NTHIN = "nthin",
      ACC1  = "acc_opt",
      ACC2  = "acc_opt2",
      SEED = "seed",
      ALG = "algorithm",
      BLOCKING = "blocking",
      BLOCKS = "blocks",
      BLOCKS_CHOL = "blocks_chol",
      BLOCKS_CHOL_SC = "blocks_chol_sc",
      BLOCKS_SC = "blocks_scaling",
      OUTF = "outfile",
      OUTADAPT = "adapt_outfile",
      OUTADAPTFMT = "adapt_outfmt",
      OUTFMT = "outfmt",
      OUTVARS = "outvars",
      CLOSE = "close_hook",
      CUSTOM_SCALING_FUN = "scaling_adapt",
      OUTCFG = "outcfg",
      DR_SCALING = "dr",
      CLIB = "clib",
      ADAPT_WEIGHT = "adapt_weight",
      ADAPT_WEIGHT_SC = "adapt_weight_sc",
      INIT = "init",
      MIX_WEIGHT = "p_mix",
      PROPOSAL = "proposal",
      RANDOM_SCAN = "random_scan",
   } --}}}
   local PARA_DEFAULTS = { --{{{
      NITER = 10e3,
      NBURN = 0,
      NTHIN = 1,
      ACC1  = 0.44,
      ACC2  = 0.234,
      SEED = nil,
      ALG = "am",
      BLOCKING = "node",
      BLOCKS = nil,
      BLOCKS_CHOL = nil,
      BLOCKS_CHOL_SC = nil,
      BLOCKS_SC = nil,
      OUTF = nil,
      OUTFMT = "bin",
      OUTADAPT = nil,
      OUTADAPTFMT = "bin",
      OUTVARS = nil,
      CLOSE = nil,
      CUSTOM_SCALING_FUN = nil,
      OUTCFG = nil,
      DR_SCALING = 0.0,
      CLIB = nil,
      ADAPT_WEIGHT = 1.0,
      ADAPT_WEIGHT_SC = 0.66,
      INIT = "greedy",
      MIX_WEIGHT = nil,
      PROPOSAL = "norm",
      RANDOM_SCAN = 0,
   } --}}}
   local PARA_IN_BLOCK = "in_block"
   local function check_integer(x) --{{{
      if type(x) == "number" and math.floor(x) == x then return true
      else return false end
   end --}}}
   local function check_number_inclusive(lower, upper, int) --{{{
      local function fun(x)
        if type(x) ~= "number" then return false
        elseif lower ~= nil and x<lower then return false
        elseif upper ~= nil and x>upper then return false
        else return int == nil or check_integer(x) end  
      end 
      return fun
   end --}}}
   local function check_number_exclusive(lower, upper, int) --{{{
      local function fun(x)
         if type(x) ~= "number" then return false
         elseif lower ~= nil and x<=lower then return false
         elseif upper ~= nil and x>=upper then return false
         else return int == nil or check_integer(x) end  
      end 
      return fun
   end --}}}
   local function check_enumerated(vals) --{{{
      local function fun(x)
         local val
         local found = false
         for _,val in ipairs(vals) do
            if val == x then found = true end
         end
         return found
      end
      return fun
   end --}}}
   local function check_type(xtype) --{{{
     local function fun(x) return type(x) == xtype end
     return fun
   end --}}}
   local function check_table_of_strings(x) --{{{
      local function fun(x)
         local k,val,n
         if type(x) ~= "table" then 
            return false
         else
            n = #x
            for k,val in pairs(x) do
               if not check_number_inclusive(1,n,1)(k) then
                  return false
               elseif type(val) ~= "string" then 
                  return false 
               end
            end
         end
         return true
      end
      return fun
   end --}}}
   -- Simple check functions for various parameters
   local PARA_CHECK = { --{{{
      NITER = check_number_inclusive(0,nil,1),
      NBURN = check_number_inclusive(0,nil,1),
      NTHIN = check_number_inclusive(0,nil,1),
      ACC1 = check_number_exclusive(0,1),
      ACC2 = check_number_exclusive(0,1),
      SEED = check_number_inclusive(0,nil,1),
      ALG = check_enumerated{"asm","am","aswam","rbam","rbaswam","ram","metropolis"},
      BLOCKING = check_enumerated{"sc","node","full"},
      OUTF = check_type("string"),
      OUTADAPT = check_type("string"),
      OUTFMT = check_enumerated{"ascii","bin"},
      OUTADAPTFMT = check_enumerated{"ascii","bin"},
      CLOSE = check_type("function"),
      CUSTOM_SCALING_FUN = check_type("function"),
      OUTCFG = check_type("string"),
      DR_SCALING = check_number_inclusive(0),
      CLIB = check_type("string"),
      ADAPT_WEIGHT = function(x)
         return type(x) == "function" or check_number_inclusive(0)(x)
      end,
      INIT = check_enumerated{"greedy","freeze","trad"},
      MIX_WEIGHT = function(x)
         return type(x) == "function" or check_number_inclusive(0)(x)
      end,
      PROPOSAL = check_enumerated{"norm","laplace","unif","student"},
   } --}}}
   
   local FUNC_FIELDS = { --{{{
      DIM = "dim",
      NAME = "name",
      ARGS = "args",
   } --}}}
   local FUNC_CHECK = { --{{{
      DIM = check_type("number"),
      NAME = check_type("string"),
      ARGS = check_table_of_strings,
   } --}}}
   
   local MODEL_FIELDS = { --{{{
      TYPE = "type",
      KIND = "kind",
      DIM = "dim",
      DENSITY = "density",
      PARENTS = "parents",
      INIT = "init_val",
      INSTANTIATED = "instantiated",
      LIMITS = "limits"
   } --}}}
   local MODEL_CHECK = { --{{{
      TYPE = check_enumerated{"number","vector","custom","matrix"},
      KIND = check_enumerated{"real","integer","symmetric"},
      DIM = check_type("table"),
      DENSITY = function(x) 
         return type(x) == "function" or type(x) == "string"
      end,
      PARENTS = check_table_of_strings,
   } --}}}
   
   local MODEL_DEFAULTS = { --{{{
      KIND = "real",
      INSTANTIATED = 0,
   } --}}}
   
   local TYPES = { --{{{
      NUMBER = "number",
      VECTOR = "vector",
      CUSTOM = "custom",
      MATRIX = "matrix",
   } --}}}

   local field, var, vtype, dim, dim2, init_val, value, k, block
   
   local function gerror(...) --{{{
      error("\nError: " .. unpack(arg))
   end --}}}
   local function gwarning(...) --{{{
      print("\nWarning: " .. unpack(arg))
   end --}}}
   
   -- Read INST_VNAME into model
   local data = _G[INST_VNAME] --{{{
   local model = _G[MODEL_VNAME]
   if data == nil then 
      data = {}
   else
      if type(data) ~= "table" then
         gerror("'" ..INST_VNAME.."' must be a table.")
      end
      for var,value in pairs(data) do
         if model[var] == nil then
            model[var] = {}
         end
         if model[var][MODEL_FIELDS.INIT] ~= nil then
            gerror("the variable '" ..var.."' in '"
                   ..INST_VNAME.. "' has already an initial value.")
         else
            model[var][MODEL_FIELDS.INIT] = value
            model[var][MODEL_FIELDS.INSTANTIATED] = 1
         end
      end
   end 
   _G[MODEL_VNAME] = model --}}}
  
   -- Read CONST_VNAME into model
   local const = _G[CONST_VNAME] --{{{
   if const ~= nil then   
      if type(const) ~= "table" then
         gerror("'" ..CONST_VNAME.."' must be a table.")
      end
      for var, value in pairs(const) do
         if _G[MODEL_VNAME][var] == nil then
            _G[MODEL_VNAME][var] = {}
            _G[MODEL_VNAME][var][MODEL_FIELDS.INIT] = value
            _G[MODEL_VNAME][var][MODEL_FIELDS.DENSITY] = "duniform"
            _G[MODEL_VNAME][var][MODEL_FIELDS.INSTANTIATED] = 1
         else
            gerror("constant variable '" ..var.."' already in " 
                   .. MODEL_VNAME ..".")
         end
         if data[var] ~= nil then
            gerror("constant variable '" ..var.. "' already in " 
                   ..INST_VNAME..".");
         end
      end
   end --}}}
   
   -- Check whether there is a given field
   local function isfield(names, field) --{{{
      if field == nil then return true end
      local v
      local found = false
      for _, v in pairs(names) do
         if v == field then 
            found = true
            break
         end
      end
      return found
   end --}}}
   
   -- Infer the type and dimension from a value
   local function determine_type_dim(value) --{{{
         local vtype,dim,dim2
         if type(value) == "number" then
            vtype = TYPES.NUMBER
            dim = 1
         elseif type(value) == "table" then
            vtype = TYPES.VECTOR
            dim = #value
            if type(value[1]) == "table" then
               dim2 = 0
               for k = 1,dim do
                  if type(value[k]) ~= "table" then
                     return nil
                  end
                  dim2 = math.max(#value[k], dim2)
               end
               dim = {dim, dim2}
            end
         elseif matrix ~= nil and getmetatable(value) == matrix._MT then
            vtype = TYPES.MATRIX
            dim = value:size()
         else
            return nil
         end
         if type(dim) == "number" then dim = {dim} end
         return vtype,dim
   end --}}}

   -- Check MODEL_VNAME 
   local model = _G[MODEL_VNAME] --{{{
   if type(model) ~= "table" then
      gerror("'" .. MODEL_VNAME .. "' must be a table.")
   else
      for var in pairs(model) do
         -- Check the variable name
         if not (type(var) == "string" and string.match(var,"^%w+[%w_]*$"))  then
            gerror("invalid variable '" .. MODEL_VNAME .. "." .. var .. "'.")
         end
         
         -- Set defaults
         for field, value in pairs(MODEL_DEFAULTS) do
            if model[var][MODEL_FIELDS[field]] == nil then
               model[var][MODEL_FIELDS[field]] = value
            end
         end
         
         if type(model[var][MODEL_FIELDS.DIM]) == "number" then
            model[var][MODEL_FIELDS.DIM] = {model[var][MODEL_FIELDS.DIM]}
         end      
         -- Simple checks of the fields
         for field, fun in pairs(MODEL_CHECK) do
            if model[var][MODEL_FIELDS[field]] ~= nil then
               if not fun(model[var][MODEL_FIELDS[field]]) then
                  gerror("invalid value of '" ..MODEL_VNAME.."."..var.."."
                         ..MODEL_FIELDS[field].."'.")
               end
            end
         end
         
         -- Check (and guess if necessary) the type/dimension
         vtype = model[var][MODEL_FIELDS.TYPE]
         dim = model[var][MODEL_FIELDS.DIM]
         init_val = model[var][MODEL_FIELDS.INIT]
         if init_val == nil and dim == nil and vtype == nil then
            dim = {1}; vtype = TYPES.NUMBER; init_val = 0.0
         else
            if init_val ~= nil then
               if vtype == nil then
                  vtype,dim2 = determine_type_dim(init_val)
                  if dim == nil then dim = dim2 end
               end
            elseif dim ~= nil and vtype == nil then
               if #dim == 1 and dim[1] == 1 then
                  vtype = TYPES.NUMBER
               else
                  vtype = TYPES.VECTOR
               end
            end
         end
         if vtype == nil then
            gerror("could not infer the type of variable '" ..var.. "'.");
         end
         if dim == nil then
            gerror("could not infer the dimension of variable '" ..var.. "'.");
         end
         if vtype == TYPES.CUSTOM and init_val == nil then
            gerror("the type '" .. vtype .. 
                   "' of variable '" .. var .. 
                   "' requires an initial value.")
         end
         if not isfield(TYPES, vtype) then
            gerror("variable '" ..var.. "' has unknown type '" ..vtype.. "'.");
         end
         
         model[var][MODEL_FIELDS.TYPE] = vtype
         model[var][MODEL_FIELDS.DIM] = dim 
         _G[MODEL_VNAME] = model
         
         -- Check the rest of the fields
         for field in pairs(model[var]) do
            if not isfield(MODEL_FIELDS, field) then
               gwarning("unrecognised field '" .. MODEL_VNAME .. "." 
                     .. var .. "." .. field .. "'");
            end
         end
      end
   end --}}}
   
   -- Find the uninstantiated variables.
   local function uninstantiated(model) --{{{
      local uvars = {}
      local var,val
      for var,val in pairs(model) do
         if model[var][MODEL_FIELDS.INSTANTIATED] == 0 then
            table.insert(uvars,var)
         end
      end
      return uvars
   end --}}}
   
   -- Parse the string with variable name and possible indices
   local function var_index(varstr,model) --{{{
      local var, i, j, dim, dim1, dim2
      var, i, j = string.match(varstr, "^%s*([_%w]+)%s*(%b[])%s*(%b[])%s*$")
      if j == nil then
         var, i  = string.match(varstr, "^%s*([_%w]+)%s*(%b[])%s*$")
         if i == nil then
            var  = string.match(varstr, "^%s*([_%w]+)%s*$")
         end
      else
         j = tonumber(string.sub(j,2,-2))
         if j == nil then return nil end
      end
      if i ~= nil then 
         i = tonumber(string.sub(i, 2, -2)) 
         if i == nil then return nil end
      end
      if model[var] == nil then return nil
      else
         dim1 = model[var][MODEL_FIELDS.DIM][1]
         dim2 = model[var][MODEL_FIELDS.DIM][2]
         if j ~= nil then
            if type(dim1) == "number" and type(dim2) == "number" then
               if i>dim1 or j>dim2 then return end
               i = i + (j-1)*dim
            else return
            end
         elseif i ~= nil then
            if type(dim1) ~= "number" then return
            elseif i>dim1 then return
            end
         end
      end
      return var, i
   end --}}}
   
   -- Check PARA_VNAME
   local para = _G[PARA_VNAME] --{{{
   local model = _G[MODEL_VNAME] 
   local fun
   if para == nil then
      para = {}
   elseif type(para) ~= "table" then
      gerror("'" .. PARA_VNAME .. "' must be a table.")
   else      
      -- Set the default values.
      for field, value in pairs(PARA_DEFAULTS) do
         if para[PARA_FIELDS[field]] == nil then
            para[PARA_FIELDS[field]] = value
         end
      end
      
      -- Check the fields
      for field, fun in pairs(PARA_CHECK) do
         if para[PARA_FIELDS[field]] ~= nil then
            if not fun(para[PARA_FIELDS[field]]) then
               gerror("invalid value of '" ..PARA_VNAME.."."..PARA_FIELDS[field]..
                      "'.")
            end
         end
      end
      
      -- Set the default for saved variables.
      if para[PARA_FIELDS.OUTFILE] ~= nil then
         if para[PARA_FIELDS.OUTVARS] == nil then
            para[PARA_FIELDS.OUTVARS] = uninstantiated(model)
         elseif type(para[PARA_FIELDS.OUTVARS]) ~= "table" then
            gerror("'" .. PARA_VNAME .. "." .. PARA_FIELDS.OUTVARS .. 
                   "' must be a table.")
         end
      end
      
      -- Do the blocking
      if para[PARA_FIELDS.BLOCKS] == nil then
         para[PARA_FIELDS.BLOCKS] = {}
      else
         if type(para[PARA_FIELDS.BLOCKS]) ~= "table" then
            gerror("'" ..PARA_VNAME .."." ..PARA_FIELDS.BLOCKS ..
                   "' must be a table.")
         else
            local i, j, k, block, varstr   
            local Blocks = {}
            local N = 0
            for k, block in ipairs(para[PARA_FIELDS.BLOCKS]) do
               if type(block) ~= "table" then
                  gerror("'" ..PARA_VNAME .."." ..PARA_FIELDS.BLOCKS ..
                   "[" .. k .. "]' must be a table .")
               else
                  local n = 0
                  if #block ~= 0 then
                     N = N+1
                     Blocks[N] = {}
                  end
                  for _, varstr in ipairs(block) do
                     var,i = var_index(varstr, model)
                     if var == nil then
                        gerror("failed to find variable '" ..varstr.. 
                               "' requested in '"
                               ..PARA_VNAME.. "." ..PARA_FIELDS.BLOCKS.."'.")
                     end
                     n = n+1
                     Blocks[N][n] = {var, i}
                  end
               end
            end
            --para[PARA_FIELDS.BLOCKS] = Blocks
         end
      end
      
      -- Finally, check for any unknown fields
      for field in pairs(para) do
         if not isfield(PARA_FIELDS, field) then
            gwarning("unrecognised field '" .. PARA_VNAME .. "." 
                     .. field .. "'");
         end
      end
   end 
   _G[PARA_VNAME] = para
   --}}}
   
   -- Check FUNC_VNAME
   local functional = _G[FUNC_VNAME] --{{{
   if functional == nil or type(functional) == "function" then
      -- do nothing
   elseif type(functional) ~= "table" then
      gerror("'" .. FUNC_VNAME .. "' must be a function or a table.")
   else
      -- Simple check
      for field, fun in pairs(FUNC_CHECK) do
         if functional[FUNC_FIELDS[field]] ~= nil then
            if not fun(functional[FUNC_FIELDS[field]]) then
               gerror("invalid value of '" ..FUNC_VNAME.."."..FUNC_FIELDS[field]..
                      "'.")
            end
         end
      end
      
      -- Check for unknown fields
      for field, _ in pairs(functional) do
         if not isfield(FUNC_FIELDS, field) then
            gwarning("Warning: unrecognised field '" .. FUNC_VNAME .. "." 
                     .. field .. "'");
         end
      end
   end --}}}

end

-- Call the preprocessing function
_PREPROCESS()

