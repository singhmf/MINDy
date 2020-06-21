%% Does the recording of variables in a separate module to clean code
%% Doesn't record error--need to do that separate
    if strcmpi(ParStr.RecCorr(1),'y')
        xPt=(W0*(PolyZ)-((Decay).*X)+C);
        yPt=dY(:,bInd);
     Rec.RecCorr{1}(:,ceil(iBatch/RecWrate))=DiagCorr(yPt',xPt');
     Rec.RecTotCorr{1}(:,ceil(iBatch/RecWrate))=corr(xPt(:),yPt(:));
    end
     if strcmpi(RecW(1),'y')
         %% Records W0 instead of W1!
        Rec.RecW{iRep}(:,:,ceil(iBatch/RecWrate))=W0;
     end
    if iBatch==1&&(size(A,2)~=1)&&(ismatrix(ndims(Rec.RecA{iRep})))
        ss=size(Rec.RecA{iRep});
        Rec.RecA{1}=zeros(ss(1),size(A,2),ss(2));
    end
    if strcmpi(RecA(1),'y')
        if size(A,2)==1
    Rec.RecA{1}(:,ceil(iBatch/RecWrate))=A;
        else            
    Rec.RecA{1}(:,:,ceil(iBatch/RecWrate))=A;
        end
            
    end                                             
        
    %% I changed to do B1 instead!!!
    if strcmpi(RecB2(1),'y')
    Rec.RecB2{1}(:,ceil(iBatch/RecWrate))=B1(:,1);
    end
    if strcmpi(RecDdiff(1),'y')
        Rec.RecD{1}(:,ceil(iBatch/RecWrate))=D;
    end