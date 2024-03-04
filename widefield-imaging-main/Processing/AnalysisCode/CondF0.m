function mu = CondF0(Tens,Tlim)

Flim = getframeidx(Tlim,1)
mu = cell(1,length(Tens));
for i = 1:length(Tens)
   mu{i} = zeros(size(Tens{1}(:,:,1))); 
end

for i = 1:length(Tens)  %loop through each condition
        i
    mu{i} = mean(Tens{i}(:,:,Flim(1):Flim(2)),3);
    
end
