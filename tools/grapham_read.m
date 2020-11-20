function [data,names] = grapham_read(fname, varargin)
% GRAPHAM_READ  Read binary data stored in Grapham to Matlab
%  
% Usage: [data,names] = grapham_read(fname, ...)
%               datas = grapham_read(fname, ...)
%
% In:
%   fname -- String containing the name of the file to be read.
%   ...   -- The optional arguments given as a list "NAME",value, where
%            "NAME" is one of the following:
%            NTHIN -- Thinning, i.e. only the value at every NTHIN'th
%                     iteration is read (default: 1).
%            NMAX  -- The maximum number of elements that are read.
%                     (default: 1e6).
%            START -- The number of records omitted from the start
%                     (default: 0).
%            BUF   -- The buffer size, i.e. the size of the block
%                     that is read at a time. (default: 1e7).
%
% Out:
%   data  -- An M-by-N matrix containing the data from "fname".
%            Each column corresponds to a variable.
%   names -- Cell array of variable names.
%   datas -- Struct array, with field corresponding to the variables
%            read from the file "fname". The variable names may be
%            altered to obtain valid field names

% tools/grapham_read.m
% 
% Copyright (c) Matti Vihola 2009-2013
% 
% This file is part of Grapham.
% 
% Grapham is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Grapham is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Grapham.  If not, see <http://www.gnu.org/licenses/>.
p = optional_arguments(...
  struct('BUF', 1e7,... 
         'NMAX', 1e6,... 
         'NTHIN', 1,... 
         'START', 0), ...
  varargin);

[fid,msg] = fopen(fname, 'r');

if fid == -1
  error(msg)
end

hdr = fgetl_nocr(fid);
if length(hdr) == 0
  error('Cannot parse the header!');
end
vars = parse_vnames(hdr);
N = length(vars);
nblock = max(N*p.NTHIN, p.BUF - rem(p.BUF, N*p.NTHIN));
dd = zeros(0,N);
fread(fid, p.START, 'double');
while true
  [dd_,Ndd]= fread(fid, nblock, 'double');
  if (Ndd == 0)
    break
  elseif rem(Ndd, N) ~= 0  
    warning('the last record appears to be broken -- data omitted!');
    M = floor(Ndd/N);
    dd_=dd_(1:M*N);
  else
    M = Ndd/N;
  end
  dd_ = reshape(dd_,N,M)';
  dd = [dd; dd_(1:p.NTHIN:end,:)];
  if size(dd,2)>p.NMAX
    dd = dd(1:p.NMAX,:);
    break;
  end
end
fclose(fid);

if nargout > 1 
  data = dd; 
  names = vars; 
else
  if N>=1
    data = struct(normalise_name(vars{1}), dd(:,1));
  else
    data = [];
  end
  for k=2:N
    data = setfield(data, normalise_name(vars{k}), dd(:,k));
  end
end

function vars = parse_vnames(hdr)
vars = {};
hdrlen = length(hdr);
if hdrlen == 0 
  return
end
s = 1;
for k=[findstr(hdr, ',') length(hdr)+1]    
  str = hdr(s+1:(k-2));
  s = k+1;
  vars{end+1} = str;
end    

function name = normalise_name(name)
for k = 1:length(name)
  c = name(k);
  v = double(c);
  if ~( (48 <= v & v <= 57) ... % 0-9
      | (65 <= v & v <= 90) ... % A-Z
      | v == 95 ...             % _
      | (97 <= v & v <= 122))   % a-z
    c = '_';
  end
  name(k) = c;
end

function p = optional_arguments(p, args)
N = length(args);
if rem(N,2) ~= 0
  error('You must specify optional arguments in pairs!');
end
for k=1:2:N
  if ~isfield(p, args{k})
    error(['Invalid parameter: ''' args{k} '''']);
  end
  p = setfield(p, args{k}, args{k+1});
end

function [str,err] = fgetl_nocr(fid)
str = fgets(fid);
if double(str(end)) == 10
  str = str(1:end-1);
elseif length(str)>=2 & double(str(end)) == 13 & double(str(end-1)) == 10
  str = str(1:end-2);
  fseek(fid, -1, 'cof');
else
  str = '';
end
