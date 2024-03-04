function [tcmat tcourseHi tcourseLo tcourseHi_var tcourseLo_var primDom blank legStr tcourseBlank] = getpixeldata(pos,W)

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

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters

%Make tuning curve
for i = 1:Nloop
    
    dum = f0m{i}(yran,xran);
    tc(i) = mean(dum(:)); 
    if varflag
        dum = f0m_var{i+1}(yran,xran);
        tc_var(i) = mean(dum(:));
    end
    
end


for i = 1:Nsym
    allDom{i} = getdomain(symbolInfo.str{i});
    dim(i) = length(allDom{i});
end
primDom = allDom{1};

%Create Ndim kernel template

switch Nsym
    case 1
        tcmat = zeros(dim(1),1);
        tcourseArray = cell(dim(1),1);
        tcourseArray_var = cell(dim(1),1);
    case 2
        tcmat = zeros(dim(1),dim(2));
        tcourseArray = cell(dim(1),dim(2));
        tcourseArray_var = cell(dim(1),dim(2));
    case 3
        tcmat = zeros(dim(1),dim(2),dim(3));
        tcourseArray = cell(dim(1),dim(2),dim(3));
        tcourseArray_var = cell(dim(1),dim(2),dim(3));
end


%Insert values at the correct location
for i = 1:length(tc)
    tcoursedum = squeeze(sum(sum(Tens{i}(yran,xran,:),1),2))/nopix;
    if varflag
        tcoursedum_var = squeeze(sum(sum(Tens_var{i}(yran,xran,:),1),2))/nopix;
    else 
        tcoursedum_var = zeros(size(tcoursedum));
    end
    vals = Analyzer.loops.conds{i}.val;
    clear loc
    for j = 1:Nsym
        loc(j) = find(allDom{j} == vals{symbolInfo.ID(j)});
    end
    
    switch Nsym
        case 1
            tcmat(loc(1)) = tc(i);
            tcourseArray{loc(1)} = tcoursedum;
            tcourseArray_var{loc(1)} = tcoursedum_var;
        case 2
            tcmat(loc(1),loc(2)) = tc(i);
            tcourseArray{loc(1),loc(2)} = tcoursedum;
            tcourseArray_var{loc(1),loc(2)} = tcoursedum_var;
        case 3
            tcmat(loc(1),loc(2),loc(3)) = tc(i);
            tcourseArray{loc(1),loc(2),loc(3)} = tcoursedum;
            tcourseArray_var{loc(1),loc(2),loc(3)} = tcoursedum_var;
    end
end

tcourseBlank = [];
if bflag
    tcourseBlank = squeeze(sum(sum(Tens{end}(yran,xran,:),1),2))/nopix;
end

legStr{1} = 'blank';
if Nsym == 3
    
    oppCollapse = symbolInfo.Collapse(2);
    
    switch oppCollapse
        
        case 1  %Take slice at maximum
            
            [v id] = find(tcmat(:) == max(tcmat(:)));
            zloc = ceil(id/(dim(1)*dim(2)));
            
            tcmat = squeeze(tcmat(:,:,zloc));            
            tcourseArray = tcourseArray(:,:,zloc);   %Cell arrays can be indexed like this apparently
            tcourseArray_var = tcourseArray_var(:,:,zloc);   %Cell arrays can be indexed like this apparently
            
        case 2  %Take mean over opposing parameters
            
            tcmat = squeeze(mean(tcmat,3));  %Take mean across last dimension       
            
            for i = 1:dim(1)
                for j = 1:dim(2)
                    tcourseNew{i,j} = 0;
                    tcourseNew_var{i,j} = 0;
                    for k = 1:dim(3)
                        tcourseNew{i,j} = tcourseNew{i,j} + tcourseArray{i,j,k}/dim(3);
                        tcourseNew_var{i,j} = tcourseNew_var{i,j} + tcourseArray_var{i,j,k}/dim(3);
                    end
                end
            end
            tcourseArray = tcourseNew;
            tcourseArray_var = tcourseNew_var;
            clear tcourseNew tcourseNew_var
            
    end
       
end

%'tcmat' should have at most 2 dimensions at this point
if Nsym > 1
    oppCollapse = symbolInfo.Collapse(1);
    
    %switch oppCollapse
    
    %         case 1  %Take slice at maximum
    %
    %             [idy idx] = find(tcmat == max(tcmat(:)));
    %             tcmat = tcmat(:,idx);
    %             tcourseArray = tcourseArray(:,idx);
    %             tcourseArray_var = tcourseArray_var(:,idx);
    %
    %             legStr{2} = num2str(round(allDom{2}(idx)*100)/100);
    
    if oppCollapse == 1  %Take mean over opposing parameters
        
        tcmat = squeeze(mean(tcmat,2));  %Take mean across last dimension
        
        for i = 1:dim(1)
            tcourseNew{i} = 0;
            tcourseNew_var{i} = 0;
            for j = 1:dim(2)
                tcourseNew{i} = tcourseNew{i} + tcourseArray{i,j}/dim(2);
                tcourseNew_var{i} = tcourseNew_var{i} + tcourseArray_var{i,j}/dim(2);
            end
        end
        tcourseArray = tcourseNew;
        tcourseArray_var = tcourseNew_var;
        clear tcourseNew tcourseNew_var
        
        legStr{2} = 'All';
        
    else 2  %Display all tuning curves, and take mean time course
        for k = 1:length(allDom{2})
            legStr{k+1} = num2str(round(allDom{2}(k)*100)/100);
        end
        %Don't do anything
        
    end
    
end

 
d = size(tcmat);
if d(1) == 1 || d(2) == 1
    [xma] = find(tcmat == max(tcmat(:)));
    [xmi] = find(tcmat == min(tcmat(:)));
    
    tcourseHi = tcourseArray{xma};
    tcourseLo = tcourseArray{xmi};
    
    tcourseHi_var = tcourseArray_var{xma};
    tcourseLo_var = tcourseArray_var{xmi};
else
    [yma xma] = find(tcmat == max(tcmat(:)));
    [ymi xmi] = find(tcmat == min(tcmat(:)));
    
    tcourseHi = tcourseArray{yma,xma};
    tcourseLo = tcourseArray{ymi,xmi};
    
    tcourseHi_var = tcourseArray_var{yma,xma};
    tcourseLo_var = tcourseArray_var{ymi,xmi};
end
