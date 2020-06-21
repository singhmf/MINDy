function[Out]=Uncellfun(f,X)
%% Same as cellfun but automatically puts in the 'uniformoutput' part
Out=cellfun(f,X,'UniformOutput',0);
end