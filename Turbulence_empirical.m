clear all;

load schaefercog.mat;

NSUB=1000;
NPARCELLS=1000;
Tmax=1200;
lambda=0.18;

% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
    end
end

enstrophy=zeros(NPARCELLS,Tmax);
enstrophy_su=zeros(NPARCELLS,Tmax);
signal_filt=zeros(NPARCELLS,Tmax);
Phases=zeros(NPARCELLS,Tmax);
Phases_su=zeros(NPARCELLS,Tmax);
Rspatime=zeros(1,NSUB);
Rspa=zeros(NSUB,NPARCELLS);
Rtime=zeros(NSUB,Tmax);
acfspa=zeros(NSUB,101);
acftime=zeros(NSUB,101);
Rspatime_su=zeros(1,NSUB);
Rspa_su=zeros(NSUB,NPARCELLS);
Rtime_su=zeros(NSUB,Tmax);
acfspa_su=zeros(NSUB,101);
acftime_su=zeros(NSUB,101);

Cexp=zeros(NPARCELLS,NPARCELLS);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        Cexp(i,j)=exp(-lambda*rr(i,j));
    end
    Cexp(i,i)=1;
end


for sub=1:NSUB
    sub
    if sub==1
        load hcp1003_REST1_LR_schaefer_1-100.mat;
    end
    if sub==101
        clear subject;
        load hcp1003_REST1_LR_schaefer_101-200.mat;
    end
    if sub==201
        clear subject;
        load hcp1003_REST1_LR_schaefer_201-300.mat;
    end
    if sub==301
        clear subject;
        load hcp1003_REST1_LR_schaefer_301-400.mat;
    end
    if sub==401
        clear subject;
        load hcp1003_REST1_LR_schaefer_401-500.mat;
    end
    if sub==501
        clear subject;
        load hcp1003_REST1_LR_schaefer_501-600.mat;
    end
    if sub==601
        clear subject;
        load hcp1003_REST1_LR_schaefer_601-700.mat;
    end
    if sub==701
        clear subject;
        load hcp1003_REST1_LR_schaefer_701-800.mat;
    end
    if sub==801
        clear subject;
        load hcp1003_REST1_LR_schaefer_801-900.mat;
    end
    if sub==901
        clear subject;
        load hcp1003_REST1_LR_schaefer_901-1003.mat;
    end
    
    ts=subject{sub}.schaeferts;
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
        Phases_su(seed,:)=Phases(seed,randperm(Tmax));
    end
    
    Phases_su=Phases_su(randperm(NPARCELLS),:);
    for i=1:NPARCELLS
        sumphases=nansum(repmat(Cexp(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)))/nansum(Cexp(i,:));
        enstrophy(i,:)=abs(sumphases);
        sumphases=nansum(repmat(Cexp(i,:)',1,Tmax).*complex(cos(Phases_su),sin(Phases_su)))/nansum(Cexp(i,:));
        enstrophy_su(i,:)=abs(sumphases);
    end
    
    Rspatime(sub)=nanstd(enstrophy(:));
    Rspatime_su(sub)=nanstd(enstrophy_su(:));
    
    Rspa(sub,:)=nanstd(enstrophy,[],2)';
    Rtime(sub,:)=nanstd(enstrophy,[],1);
    
    acfspa(sub,:)=autocorr(Rspa(sub,:),100);
    acftime(sub,:)=autocorr(Rtime(sub,:),100);
   
    Rspa_su(sub,:)=nanstd(enstrophy_su,[],2)';
    Rtime_su(sub,:)=nanstd(enstrophy_su,[],1);
    
    acfspa_su(sub,:)=autocorr(Rspa_su(sub,:),100);
    acftime_su(sub,:)=autocorr(Rtime_su(sub,:),100);
end


save turbu_emp.mat Rspatime Rspatime_su Rspa Rtime Rspa_su Rtime_su acfspa acftime acfspa_su acftime_su;


save turbu_emp_Phases.mat Phases Phases_su;

load turbu_emp_Phases.mat;
load turbu_emp.mat;
NPARCELLS=1000;
Turbu=nanmean(Rspatime)
Turbu_su=nanmean(Rspatime_su)

p=ranksum(Rspatime,Rspatime_su)

figure(1)
boxplot([Rspatime' Rspatime_su']);

%% std across tome vs space
figure(2)
shadedErrorBar(1:NPARCELLS,nanmean(Rspa),nanstd(Rspa),'-r',0.7)
hold on;
shadedErrorBar(1:NPARCELLS,nanmean(Rspa_su),nanstd(Rspa_su),'-k',0.7)

%% std across space as a function of time (I cut border effects...thus 100:1100 (time)
figure(3)
shadedErrorBar(100:1100,nanmean(Rtime(:,100:1100)),nanstd(Rtime(:,100:1100)),'-r',0.7)
hold on;
shadedErrorBar(100:1100,nanmean(Rtime_su(:,100:1100)),nanstd(Rtime_su(:,100:1100)),'-k',0.7)

%% autocorr space
figure(4)
shadedErrorBar(1:101,nanmean(acfspa),nanstd(acfspa),'-r',0.7)
hold on;
shadedErrorBar(1:101,nanmean(acfspa_su),nanstd(acfspa_su),'-k',0.7)

%% autocorr time 
figure(5)
shadedErrorBar(1:101,nanmean(acftime),nanstd(acftime),'-r',0.7)
hold on;
shadedErrorBar(1:101,nanmean(acftime_su),nanstd(acftime_su),'-k',0.7)

%% Phases
figure(6)
plot(Phases(:,500),'r*');
hold on;
plot(Phases_su(:,500),'k*');

%% check:
for i=400:2:600
    i
    plot(Phases(:,i),'r*');
    pause;
end
