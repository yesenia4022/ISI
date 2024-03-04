function [tcOut tcourseHi tcourseLo tcourseHi_var tcourseLo_var tcourseBlank primDom blank legStr] = getpixeldata2(pos,W)

global Tens Tens_var f0m f0m_var Analyzer symbolInfo G_handles

varflag = get(G_handles.EbarFlag,'Value');

if isempty(Tens_var{1})
    varflag = 0;
    %disp('To make things run a little faster, variance tensor not generated')
end


xran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));
yran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));
nopix = length(yran)*length(xran);

nc = getnoconditions;

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
Nloop = nc;
blank = [];
if bflag
    Nloop = nc-1;    
    dum = f0m{end}(yran,xran);
    blank = mean(dum(:));
end

Nsym = length(Analyzer.L.param);  %number of looping parameters
for i = 1:Nsym
    allDom{i} = eval(Analyzer.L.param{i}{2});
    Pdim(i) = length(allDom{i});
end


%Make tuning curve
for i = 1:Nloop %Will not include blank
    
    dum = f0m{i}(yran,xran);
    tc(i) = mean(dum(:)); 
    
end

if length(Pdim) == 1
    Pdim = [1 Pdim];
end

tcmat = reshape(tc,Pdim);


%Don't use getdomain.m. We want the vector in the order that it was input in
%the looper.
primDom = eval(Analyzer.L.param{symbolInfo.ID(1)}{2});


tcmatC = tcmat;
for i = 2:Nsym
    collapseDim = symbolInfo.ID(i); %dimension to collapse
    switch symbolInfo.Collapse(i-1)         
        case 1  %Max           
            tcmatC = max(tcmatC,[],collapseDim);
        case 2  %Mean  
            tcmatC = mean(tcmatC,collapseDim);
    end
end

tcmatC = squeeze(tcmatC);
[dum idBest] = max(tcmatC);
[dum idWorst] = min(tcmatC);

tcOut = tcmatC;

%%

%Insert values at the correct location
Nt = size(Tens{1},3);
tcourseArray = zeros(Nt,length(tc));
for i = 1:length(tc)
    tcoursedum = squeeze(sum(sum(Tens{i}(yran,xran,:),1),2))/nopix;
    tcourseArray(:,i) = tcoursedum(:);
end


for i = 1:size(tcourseArray(:,1))   %loop each time point
    tcourseMat{i} = reshape(tcourseArray(i,:),Pdim);
end


for j = 1:length(tcourseMat)
    
    Matdum = tcourseMat{j};
    for i = 2:Nsym
        collapseDim = symbolInfo.ID(i); %dimension to collapse
        switch symbolInfo.Collapse(i-1)
            case 1  %Max
                Matdum = max(Matdum,[],collapseDim);
            case 2  %Mean
                Matdum = mean(Matdum,collapseDim);
                            case 3  %Mean
                Matdum = mean(Matdum,collapseDim);
        end
    end
    tcourseAll(j,:) = squeeze(Matdum);
    
end

tcourseHi = tcourseAll(:,idBest(1));
tcourseLo = tcourseAll(:,idWorst(1));

tcourseBlank = [];
if bflag
    tcourseBlank = squeeze(sum(sum(Tens{end}(yran,xran,:),1),2))/nopix;
end


tcourseHi_var = [];
tcourseLo_var = [];
legStr = '';
