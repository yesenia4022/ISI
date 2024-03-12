%ori=zero  ws5 u007 imrate=20hz 
global repDom Tens G_handles

b = str2double(get(G_handles.bstart,'String')); %in msec as well
b(2) = str2double(get(G_handles.bstop,'String'));
slowMo = get(G_handles.slowMotionFlag,'Value');



nt = getnotrials;
nc= getnoconditions;
nreps=nt/nc;
Tens_0= cell(1,nreps);

for i=1:nreps
    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
    Tens_0{i}=Tens{1};
    
    
end

figure()
for i=1:21
    subplot(2, 11, i)
    imagesc(Tens_0{2}(:,:,i), [0 800])
    colormap jet
    colorbar
end

% figure(), imagesc(Tens_0{40}(:,:,15), [0 .1]), colormap jet

f_mu0=zeros(1,nreps);
for i=1:nreps
    f_mu0(i)=mean2(mean(Tens_0{i}(500:550,500:550,15:20),3));
   
end


%% 


figure(), scatter(1:nreps, f_mu0), ylim([-.1 .1]), xlabel('trial number'), ylabel('delta F'), title('ws5 u007 V1 Ori 0 only ')

%long exposure
nfr= size(Tens_0{1},3);
f_mu=zeros(1,nfr);
for i=1:nfr
    f_mu(i)=mean2(Tens_0{1}(500:550,500:550,i));
    
end
% time= nfr/10; %10 Hz imrate   
figure(), plot(0:.1:121, f_mu),   ylim([-.1 .1]), xlabel('time (sec)'), ylabel('delta F'), title('xl2 V1 ori 0 60 sec continuous') 
    
    

%%0/90 adaptation

figure()
for i=1:21
    subplot(4,6,i)
    imagesc(Tens{3}(:,:,i), [0 .1])
    colormap jet
    colorbar
end

f_mu=zeros(1,length(Tens));
for i=1:length(Tens)
    f_mu(i)=mean2(mean(Tens{i}(500:550,500:550,15:20),3));
end

%looking at responses as trials (and thus adaptation) progress
f_mu=zeros(10,7);
for i=1:10

    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
%     Tens_0{i}=cat(3, Tens{1},Tens{2}, Tens{3},Tens{5},Tens{6}, Tens{7}  ;
    for j=1:length(Tens)
        f_mu(i,j)=mean2(mean(Tens{j}(500:550,500:550,15:20),3));
    end
    
end



figure()
subplot(2,1,1)
plot(1:10, f_mu(:,4))
xlabel('trial #')
title('ori=90')
subplot(2,1,2)
plot(1:10, f_mu(:,5))
xlabel('trial #')
title('ori=0 adapter preceding 90 deg presentation')

figure()
for i=1:7
    subplot(3,3,i)
    plot(1:10, f_mu(:,i))



end


%long Adapt 40 sec

nt = getnotrials;
nc= getnoconditions;
nreps=nt/nc;
Tens_longAdapt= cell(1,nreps);

for i=1:nreps
    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
    Tens_longAdapt{i}=Tens;
    
    
end

for j=1:20
figure()
for i=1:16
   subplot(4,5,i)
   imagesc(Tens_longAdapt{j}{3}(:,:,i), [0, .1])
   colorbar
end
end
%%
f_muLM=zeros(20,2);
for i=1:20
    
    for j=1:2
        f_muLM(i,j)=mean2(mean(Tens_longAdapt{i}{j+1}(350:480, 270:330,16:21),3));
    end
    
end
ori_std= std(f_muLM,1);
ori_mu=mean(f_muLM,1);
figure(), errorbar([90 0], ori_mu, ori_std), ylim([-.02 .04]), xlim([-10 100])
xlabel('orientation(deg)'), ylabel('dF/F'),title('xb7 u006 longAdapt 90 deg LM oreintation tuning')

f_muRL=zeros(20,2);
for i=1:20
    
    for j=1:2
        f_muRL(i,j)=mean2(mean(Tens_longAdapt{i}{j+1}(600:630,360:400,16:21),3));
    end
    
end
ori_std= std(f_muRL,1);
ori_mu=mean(f_muRL,1);
figure(), errorbar([90 0], ori_mu, ori_std), ylim([-.02 .04]), xlim([-10 100])
xlabel('orientation(deg)'), ylabel('dF/F'),title('xb7 u006 longAdapt 90 deg RL oreintation tuning')

f_muV1=zeros(20,2);
for i=1:20
    
    for j=1:2
        f_muV1(i,j)=mean2(mean(Tens_longAdapt{i}{j+1}(400:450,450:500,16:21),3));
    end
    
end
ori_std= std(f_muV1,1);
ori_mu=mean(f_muV1,1);
figure(), errorbar([90 0], ori_mu, ori_std), ylim([-.02 .05]), xlim([-10 100])
xlabel('orientation(deg)'), ylabel('dF/F'),title('xb7 u006 longAdapt 90 deg V1 oreintation tuning')



figure()
for i=1:20
    subplot(10,2,i)
    plot(f_mu(i,:))
    ylim([0 .1])
end

ori_std= std(f_mu,1);
ori_mu=mean(f_mu,1);
figure(),
errorbar(0:45:315, ori_mu, ori_std)
ylim([-.02 .12])
xlim([-10 350])
xlabel('orientation(deg)')
ylabel('dF/F')
title('xb7 u005 longAdapt 0 deg oreintation tuning')

for i=2:2:10
    figure()
    for j=1:31
        subplot(3,11,j)
        imagesc(Tens_longAdapt{i}{4}(:,:,j), [0, .2])
        colorbar
    end
    title(['90deg trial ' num2str(i)])
end

ori_vec=[0 90];
for i=1:2
    figure() 
    
    for j=1:20
        subplot(5,4,j)
        mu_frame=zeros(1,16);
        for z=1:16
            mu_frame(z)=mean2(Tens_longAdapt{j}{i+1}(400:450,450:500,z));
        end
        plot(0:.1:1.54, mu_frame)
        ylim([-.02 .08])
        xlabel('t(sec)')
        title([num2str(ori_vec(i)) 'deg trial' num2str(j)])
    end
end


ori_vec=[90 0];
mu_frame_cond=cell(1,2);
for i=1:2
    figure() 
    mu_frame_trials=zeros(20,26);
    for j=1:20
        subplot(5,4,j)
        mu_frame=zeros(1,16);
        for z=1:26
            mu_frame(z)=mean2(Tens_longAdapt{j}{i+1}(350:480, 270:330,z));
        end
        plot(0:.1:2.54, mu_frame)
        ylim([-.02 .06])
        xlabel('t(sec)')
        title([num2str(ori_vec(i)) 'deg trial' num2str(j)])
        mu_frame_trials(j,:)=mu_frame;
    end
    mu_frame_cond{i}=mu_frame_trials;
end

figure(), plot(0:.1:2.54,mean(mu_frame_cond{1}, 1)), hold on, plot(0:.1:2.54,mean(mu_frame_cond{2}, 1)), hold off





%interleaved long adapt
Tens_update_trials=cell(1,10);
for i=1:10
    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
    nFrames=size(Tens{1},3);
    adaptTens=zeros(size(Tens{1},1),size(Tens{1},2), nFrames);
    for j=1:nFrames
        fr= cat(3, Tens{2}(:,:,j), Tens{4}(:,:,j), Tens{6}(:,:,j), Tens{8}(:,:,j), Tens{10}(:,:,j), Tens{12}(:,:,j), Tens{14}(:,:,j), Tens{16}(:,:,j));
        mimg= mean(fr,3);
        adaptTens(:,:,j)=mimg;
        
    end
    Tens_update=cell(1,9); %actual number of unique conditions
    Tens_update{1}=adaptTens;
    for j=2:9
        Tens_update{j}= Tens{j};
    end
    Tens_update_trials{i}=Tens_update;
end


f_mu=zeros(10,8);
for i=1:10
    
    for j=1:8
        f_mu(i,j)=mean2(mean(Tens_update_trials{i}{j+1}(400:450,400:450  ,12:16),3));
    end
    
end

ori_std= std(f_mu,1);
ori_mu=mean(f_mu,1);
figure(),
errorbar(0:45:315, ori_mu, ori_std)
ylim([-.02 .12])
xlim([-10 350])
xlabel('orientation(deg)')
ylabel('dF/F')
title('xb7 u005 interleave Adapt 400:450,400:450 ')






ori_vec=[0:45:315];
for i=1:8
    figure() 
    
    for j=1:10 %reps
        subplot(2,5,j)
        mu_frame=zeros(1,31);
        for z=1:31
            mu_frame(z)=mean2(Tens_update_trials{j}{i+1}(500:550,500:550,z));
        end
        plot(0:1/20:1.54, mu_frame)
        ylim([-.02 .2])
        xlabel('t(sec)')
        title([num2str(ori_vec(i)) 'deg trial' num2str(j)])
    end
end









nFrames=size(Tens{1},3);
adaptTens=zeros(size(Tens{1},1),size(Tens{1},2), nFrames);
for i=1:nFrames
    fr= cat(3, Tens{2}(:,:,i), Tens{4}(:,:,i), Tens{6}(:,:,i), Tens{8}(:,:,i), Tens{10}(:,:,i), Tens{12}(:,:,i), Tens{14}(:,:,i), Tens{16}(:,:,i));
    mimg= mean(fr,3);
    adaptTens(:,:,i)=mimg;

end

Tens_update=cell(1,9); %actual number of unique conditions
Tens_update{1}=adaptTens;

for i=2:9
   Tens_update{i}= Tens{i};    
end

figure()
ori_vec=0:45:315;
for i=1:8
    subplot(2,4,i)
    mu_frame=zeros(1,10);
    for z=1:10
        mu_frame(z)=mean2(Tens_update{i+1}(500:550,500:550,z));
    end
    plot(1:10, mu_frame)
    ylim([-.02 .2])
    xlabel('t(sec)')
    title([num2str(ori_vec(i)) 'deg'])
end

figure()
for i=1:31
    subplot(3,11,i)
    imagesc(Tens_update{4}(:,:,i), [0 .1])
    colorbar
end


f_mu=zeros(1,8);
for i=1:8
%     ori_vec(i)=Analyzer.loops.conds{i}.val{1};
    f_mu(i)=mean2(mean(Tens_update{i+1}(500:550,500:550,8:10),3));
end
ori_vec=[0 45 90 135 180 225 270 315 ];
figure(), scatter(ori_vec, f_mu),  xlabel('ori'), ylabel('delta F'), title('xl2 V1 Ori Adaptation biased at 1 degree, contrast 100%'), ylim([0 .07])




