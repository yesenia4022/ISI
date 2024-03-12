%For interleaved adaptation experiment
global Tens Analyzer f0m Flim
nFrames=size(Tens{1},3);
adaptTens=zeros(size(Tens{1},1),size(Tens{1},2), nFrames);
for i=1:nFrames
    fr= cat(3, Tens{2}(:,:,i), Tens{4}(:,:,i), Tens{6}(:,:,i), Tens{8}(:,:,i), Tens{10}(:,:,i), Tens{12}(:,:,i), Tens{14}(:,:,i), Tens{16}(:,:,i));
    mimg= mean(fr,3);
    adaptTens(:,:,i)=mimg;

end

figure()
for i=1:nFrames
    subplot(4,4,i)
    imagesc(adaptTens(:,:,i), [0, .1])
    colormap jet
    colorbar
    
end

Tens_update=cell(1,9); %actual number of unique conditions
Tens_update{1}=Tens{1};
Tens_update{2}=adaptTens;
z=3;
for i=3:9
   Tens_update{i}= Tens{z};
   z=z+2 
    
end
Analyzer.L.param{1}{2}= '[0 1 45 90 135 180 225 270 315]';
Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); 
f0m = CondF0(Tens_update,Flim);

Analyzer.loops.conds=cell(1,9);
Analyzer.loops.conds{1}.val={0};
Analyzer.loops.conds{1}.symbol={'ori'}; 
Analyzer.loops.conds{2}.val={1}; 
Analyzer.loops.conds{2}.symbol={'ori'}; 
z=45;
for i=3:9
    Analyzer.loops.conds{i}.val={z};
    Analyzer.loops.conds{i}.symbol={'ori'}; 
    z=z+45
end
%plot time course for each condition
for i =3:9
   figure()
   for j=1:nFrames
    subplot(4,4,j)
    imagesc(Tens_update{i}(:,:,j), [0, .1])
    colormap jet
    colorbar
    
   end
    
end

% %plot 14th frame (approx. peak reponse) for each condition 
% figure()
% for i=1:9
%    subplot(2,5,i)
%    imagesc(Tens_update{i}(350:450,300:600,14), [0, .1])
%    colormap jet
%    colorbar
% %    title(Analyzer.loops.conds{i}.val)
% end
f_mu=zeros(1,9);
for i=1:9
%     ori_vec(i)=Analyzer.loops.conds{i}.val{1};
    f_mu(i)=mean2(mean(Tens_update{i}(400:500,400:500,6:end),3));
end
ori_vec=[0 1 45 90 135 180 225 270 315 ];
figure(), plot(ori_vec, f_mu), ylim([0 .05]), xlabel('ori'), ylabel('delta F'), title('xl2 V1 Ori Adaptation biased at 1 degree, contrast 100%')
hold on

%unbiased 

ori_vec=zeros(1,8);
f_mu=zeros(1,8);
for i=1:8
    ori_vec(i)=Analyzer.loops.conds{i}.val{1};
    f_mu(i)=mean2(mean(Tens{i}(400:500,400:500,6:end),3));
end
plot(ori_vec, f_mu), ylim([0 .05]), xlabel('ori'), ylabel('delta F'), %title('xl2 V1 Ori unbiased, contrast 100%'), 
hold off




















