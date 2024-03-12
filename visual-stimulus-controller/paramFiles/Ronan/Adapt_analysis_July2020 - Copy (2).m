%0 vs 90 30 reps
global repDom Tens G_handles

b = str2double(get(G_handles.bstart,'String')); %in msec as well
b(2) = str2double(get(G_handles.bstop,'String'));
slowMo = get(G_handles.slowMotionFlag,'Value');



nt = getnotrials;
nc= getnoconditions;
nreps=nt/nc;
Tens_0= cell(1,nreps);
Tens_90=cell(1,nreps);
for i=1:nreps
    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
    Tens_0{i}=Tens{1};
    Tens_90{i}=Tens{2};
    
end

figure(), imagesc(Tens_0{1}(400:500,400:500,16),[0 .1]), colormap jet

figure()
for i =1:33
   subplot(3,11,i) 
   imagesc(Tens_0{1}(:,:,i),[0 .1])
   colormap jet
   colorbar
end


f_mu0=zeros(1,nreps);
f_mu90=zeros(1,nreps);

for i=1:nreps
    f_mu0(i)=mean2(mean(Tens_0{i}(400:500,400:500,8:33),3));
    f_mu90(i)=mean2(mean(Tens_90{i}(400:500,400:500,8:33),3));
end

figure(), scatter(1:nreps, f_mu0), hold on, scatter(1:nreps, f_mu90), ylim([-.1 .1]), xlabel('trial number'), ylabel('delta F'), title('xl2 V1 Ori 0/90 '), hold off

%zero deg only
nt = getnotrials;
nc= getnoconditions;
nreps=nt/nc;
Tens_0= cell(1,nreps);
for i=1:nreps
    repDom=i;
    [Tens Tens_var] = CondTensor3(b,slowMo);
    Tens_0{i}=Tens{1};
    
    
end

f_mu0=zeros(1,nreps);
for i=1:nreps
    f_mu0(i)=mean2(mean(Tens_0{i}(400:500,400:500,8:33),3));
   
end

figure(), scatter(1:nreps, f_mu0), ylim([-.1 .1]), xlabel('trial number'), ylabel('delta F'), title('xl2 V1 Ori 0 only ')

%zero deg for 1 min continuous
nfr= size(Tens{1},3);
f_mu=zeros(1,nfr);
for i=1:nfr
    f_mu(i)=mean2(Tens{1}(500:550,500:550,i));
    
end
% time= nfr/10; %10 Hz imrate   
figure(), plot(0:.1:120.5,f_mu),   ylim([-.02 .12]), xlabel('time (sec)'), ylabel('delta F'), title('xl2 V1 ori 0 120 sec continuous') 
    
    
