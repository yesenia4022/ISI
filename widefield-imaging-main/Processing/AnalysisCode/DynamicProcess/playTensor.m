function playTensor(Tens,fp,varargin)


mu = mean(Tens,3);
sig = std(Tens,[],3);
for i = 1:length(Tens(1,1,:))
    Tens(:,:,i) = (Tens(:,:,i)-mu);
end

if ~isempty(varargin)
    mi = max(Tens(:));
    ma = min(Tens(:));
     
    bw = varargin{1};

    for i = 1:length(Tens(1,1,:))
        
        kdum = Tens(:,:,i);
        madum = max(kdum(find(bw)));
        midum = min(kdum(find(bw)));

        ma = max([madum(1) ma]);
        mi = min([midum(1) mi]);

    end
    
else
    ma = max(Tens(:));
    mi = min(Tens(:));
end

figure(1)
for i = 1:length(Tens(1,1,:))
    
    imagesc(Tens(:,:,i),[mi ma]), colormap gray
    pause(fp)
    
end

