function [kernPop_neuron kernPop_glia popBlank_neuron popBlank_glia CoM_neuron CoM_glia pardom] = GetCellKernels(f0,cellmask1,cellmask2)

%Returns the kernel response matrices and corresponding locations for cells
%in cellmask

%'nflag' if 1 then it takes the glia.  If 0, then the neurons

global Analyzer
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition


Np = length(Analyzer.loops.conds{1}.symbol);  %number of looping params

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    f0blank = f0{end};
    f0(end) = [];
end

for z = 1:Np  %loop through each loop parameter
    for i = 1:length(f0)
        paramlist{z}(i) = Analyzer.loops.conds{i}.val{z};
    end
end

for z = 1:Np    
    paramlist{z} = round(paramlist{z}*1000)/1000;    
    pardom{z} = unique(paramlist{z});    
end

cmdum = cellmask1;  %Cell mask within this ROI
masklabel = bwlabel(cmdum);
celldom = unique(masklabel);
celldom = celldom(2:end); %first element is the neuropil
Ncell = length(celldom);

dim = size(cellmask1);
CoM = zeros(Ncell,2);
Nneuron = 0;
Nglia = 0;
for p = 1:Ncell
    [idcelly idcellx] = find(masklabel == celldom(p));
    
    idvec = dim(1)*(idcellx-1) + idcelly;
    
        if sum(cellmask2(idvec)) == 0  %if there is no overlap, call it a neuron
            Nneuron = Nneuron+1;
            CoM_neuron(Nneuron,:) = [mean(idcelly) mean(idcellx)];  %center of mass
            id_neuron{Nneuron} = find(masklabel(:) == celldom(p));
        else
            Nglia = Nglia+1;
            CoM_glia(Nglia,:) = [mean(idcelly) mean(idcellx)];  %center of mass
            id_glia{Nglia} = find(masklabel(:) == celldom(p));
        end

end


for i = 1:length(f0)
   
        for z = 1:Np
            kernloc(z) = find(pardom{z} == paramlist{z}(i));  %location in kernel
        end
        
        %Get kernel value for each cell
        
        %Neurons...
        for p = 1:Nneuron

            R = mean(f0{i}(id_neuron{p}));
            
            switch Np
                
                case 1
                    kernPop_neuron(kernloc(1),p) = R;                
                case 2                    
                    kernPop_neuron(kernloc(1),kernloc(2),p) = R;
            end

        end
        
        %Glia...
        for p = 1:Nglia

            R = mean(f0{i}(id_glia{p}));
            
            switch Np
                
                case 1
                    kernPop_glia(kernloc(1),p) = R;                
                case 2                    
                    kernPop_glia(kernloc(1),kernloc(2),p) = R;
            end

        end
        
end

if bflag
    for p = 1:Nneuron
        popBlank_neuron(p) = mean(f0blank(id_neuron{p}));
    end
    
    for p = 1:Nglia
        popBlank_glia(p) = mean(f0blank(id_glia{p}));
    end
else
    popBlank = [];
end


