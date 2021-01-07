function[Out]=Uncellfun(f,varargin)
%% Same as cellfun but automatically puts in the 'uniformoutput' part
Out=cellfun(f,varargin{:},'UniformOutput',0);
end