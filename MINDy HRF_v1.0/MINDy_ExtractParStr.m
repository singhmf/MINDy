ParStrFields=fields(ParStr);
for iPSF=1:numel(ParStrFields)
    assignin('caller',ParStrFields{iPSF},ParStr.(ParStrFields{iPSF}));
end