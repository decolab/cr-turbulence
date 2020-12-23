clear all;
load empirical_spacorr_rest.mat;
load schaefercog.mat;

NSUB=1003;
NPARCELLS=1000;
NR=400;
NRini=20;
NRfin=80;
NSUBSIM=200;

empcorrfcn=corrfcn;
empgrandcorrfcn=grandcorrfcn;

for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end

Markov_SC=1;
if Markov_SC==1
    lambda=0.18;
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C(i,j)=exp(-lambda*rr(i,j));
        end
        C(i,i)=0;
    end
else
    load sc_schaefer.mat
    C=sc_schaefer;
    C=C/max(max(C));
end

C1=C;
for i=1:NPARCELLS
    C1(i,i)=1;
end

neighbours=cell(1,NPARCELLS);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        r=rr(i,j);
        index=floor(r/delta)+1;
        if index==NR+1
            index=NR;
        end
        if index>1 && index<=35
            neighbours{i}=[neighbours{i} j];
        end
    end
end
%%%
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

%%%%
fce1=zeros(100,NPARCELLS,NPARCELLS);
tss=zeros(NPARCELLS,1200);
PowSpect=zeros(600,NPARCELLS,100);

% skip preprocessing of 1003 participants
load hpcdata1003_f_diff_fce.mat

% Parameters for Hopf model
Tmax=1200;
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

linfunc = @(A, x)(A(1)*x+A(2));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');

fcsimul=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
corrfcn=zeros(NSUBSIM,NPARCELLS,NR);
ensspasub=zeros(NSUBSIM,NPARCELLS);
ensspasub1=zeros(NSUBSIM,NPARCELLS);
DTspatime=zeros(NPARCELLS,Tmax);
DTspatime1=zeros(NPARCELLS,Tmax);
Rsub=zeros(1,NSUBSIM);
DTsub=zeros(1,NSUBSIM);

G_range=0.5:0.05:2;
enstrophy_r=zeros(length(G_range),NR);

ii=1;
for G=G_range
    
    G
    wC = G*C;
    sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
    
    fcsimul=zeros(NSUBSIM,NPARCELLS,NPARCELLS);
    ensspasub=zeros(NSUBSIM, NPARCELLS);
    corrfcn=zeros(NSUBSIM,NPARCELLS,NR);
    %% Hopf Simulation
    for sub=1:NSUBSIM
        a=-0.02*ones(NPARCELLS,2);
        xs=zeros(Tmax,NPARCELLS);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
        z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 2000 time steps
        for t=0:dt:2000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        end
        % actual modelling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts=xs';
        signal_filt=zeros(NPARCELLS,Tmax);
        Xanalytic=zeros(NPARCELLS,Tmax);
        Phases=zeros(NPARCELLS,Tmax);
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,1:1200)=filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
        
        %%% Perturbation
        a=(-0.02+0.02*rand(NPARCELLS,1)).*ones(NPARCELLS,2);
        nn=0;
        
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        ts=xs';
        Rspatime1=zeros(NPARCELLS,Tmax);
        Tspatime1=zeros(NPARCELLS,Tmax);
        
        signal_filt2=zeros(NPARCELLS,Tmax);
        Xanalytic2=zeros(NPARCELLS,Tmax);
        Phases2=zeros(NPARCELLS,Tmax);
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt2(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic2 = hilbert(demean(signal_filt2(seed,:)));
            Phases2(seed,:) = angle(Xanalytic2);
        end
        
        parfor i=1:NPARCELLS
            % compute enstrophy
            enstrophy=nansum(C1(i,:)'.*complex(cos(Phases2),sin(Phases2)))/sum(C1(i,:));
            Rspatime1(i,:)=abs(enstrophy);
            Tspatime1(i,:)=angle(enstrophy);
        end
        
        DTspatime1=zeros(NPARCELLS,Tmax);
        for i=1:NPARCELLS
            DTspatime1(i,:)=Tspatime1(i,:)-nanmean(Tspatime1(neighbours{i},:));
        end
        ensspasub1(sub,:)=(nanmean(Rspatime1,2))';
        
        %%%  end perturbation
        
        fcsimulc=zeros(NPARCELLS,NPARCELLS);
        fcsimulc(:,:)=corrcoef(signal_filt'); %compute current corr of fc
        fcsimul(sub,:,:)=zeros(NPARCELLS,NPARCELLS);
        fcsimul(sub,:,:)=corrcoef(signal_filt'); % save this fc
        
        
        % Calculate enstrophy
        Rspatime=zeros(NPARCELLS,Tmax);
        Tspatime=zeros(NPARCELLS,Tmax);
        corrfcnt=zeros(NPARCELLS,NR);
        parfor i=1:NPARCELLS
            numind=zeros(1,NR);
            corrfcn_1=zeros(1,NR);
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                %                mcc=fcsimul(sub,i,j);
                mcc=fcsimulc(i,j); % get the current value
                if ~isnan(mcc)
                    corrfcn_1(index)=corrfcn_1(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
            
            corrfcnt(i,:)=corrfcn_1./numind;
            
            %           corrfcn(sub,i,:)=corrfcn_1./numind; % save current corrfcn
            
            %%% enstrophy
            enstrophy=nansum(C1(i,:)'.*complex(cos(Phases),sin(Phases)))/sum(C1(i,:));
            Rspatime(i,:)=abs(enstrophy);
            Tspatime(i,:)=angle(enstrophy);
        end
        corrfcn(sub,:,:)=zeros(NPARCELLS,NR);
        corrfcn(sub,:,:)=corrfcnt(:,:);
        
        DTspatime=zeros(NPARCELLS,Tmax);
        for i=1:NPARCELLS
            DTspatime(i,:)=Tspatime(i,:)-nanmean(Tspatime(neighbours{i},:));
        end
        
        Rsub(sub)=nanstd(Rspatime(:));
        DTsub(sub)=nanstd(DTspatime(:));
        ensspasub(sub,:)=(nanmean(Rspatime,2))';
        save(['Progress_' ndate '_G_' num2str(G) '_sub_' num2str(sub)],'G','sub');
    end % subsim

    %% compute measures 

    fcsim=squeeze(nanmean(fcsimul,1));

    % compute segint, ie mean segregation/integration
    [MM QQ]=community_louvain(fcsim,[],[],'negative_sym');
    modularity=QQ;
    integration=nanmean(nanmean(abs(fcsim)-eye(NPARCELLS)));
    segint(ii)=modularity*integration

    % compute mean and std for each simulated matrix
    asegint=zeros(1,NSUBSIM);
    for t=1:NSUBSIM
        currfcsim=squeeze(fcsimul(t,:,:));
        [MM QQ]=community_louvain(currfcsim,[],[],'negative_sym');
        modularity=QQ;
        integration=nanmean(nanmean(abs(fcsim)-eye(NPARCELLS)));
        asegint(t)=modularity*integration;
    end;
    asegintmean(ii)=nanmean(asegint);
    asegintstd(ii)=nanstd(asegint);

    % compute fcfitt
    ccaux=corrcoef(fce(Isubdiag),fcsim(Isubdiag),'rows','pairwise');
    fcfitt(ii)=ccaux(2)
    % compute mean and std for each simulated matrix
    afcfitt=zeros(1,NSUBSIM);
    for t=1:NSUBSIM
        currfcsim=squeeze(fcsimul(t,:,:));
        ccaux=corrcoef(fce(Isubdiag),currfcsim(Isubdiag),'rows','pairwise');
        afcfitt(t)=ccaux(2);
    end;
    afcfittmean(ii)=nanmean(afcfitt)
    afcfittstd(ii)=nanstd(afcfitt)
    

    % compute fcfitt_hete
    fce2=fce-eye(NPARCELLS);
    GBCemp=nanmean(fce2,2);
    fcsim2=fcsim-eye(NPARCELLS);
    GBCsim=nanmean(fcsim2,2);
    for i=1:NPARCELLS
        errgbc1(i)=sqrt((GBCsim(i)-GBCemp(i))^2);
    end
    fcfitt_hete(ii)=nanmean(errgbc1)
 
    % compute Rmeta, mean and std
    Rmeta(ii)=nanmean(Rsub)
    Rmetastd(ii)=nanstd(Rsub)

    % compute DTmeta, mean and std
    DTmeta(ii)=nanmean(DTsub)
    DTmetastd(ii)=nanstd(DTsub)
    
    % compute infocapacity, mean and std
    infocapacity(ii)=nanmean(nanstd(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)))
    infocapacitystd(ii)=nanstd(nanstd(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)))

    % compute susceptibility, mean and std
    susceptibility(ii)=nanmean(nanmean(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)))
    susceptibilitystd(ii)=nanstd(nanmean(ensspasub1-ones(NSUBSIM,1)*nanmean(ensspasub)))
    
    corrfcn_he=squeeze(nanmean(corrfcn));
    grandcorrfcn=squeeze(nanmean(corrfcn_he));
    %% Comparison
    for i=1:NPARCELLS
        for k=NRini:NR
            err11(k)=sqrt((corrfcn_he(i,k)-empcorrfcn(i,k))^2);
        end
        err1(i)=nanmean(err11);
    end
    err_hete(ii)=nanmean(err1)
    
    for k=NRini:NR
        errg1(k)=sqrt((grandcorrfcn(k)-empgrandcorrfcn(k))^2);
    end
    err_grand(ii)=nanmean(errg1)
    
    %%% Powerlaw a
   
    clear xcoor;
    clear ycoor;
    nn=1;
    for k=NRini:NRfin
        if grandcorrfcn(k)>0
            xcoor(nn)=log(xrange(k));
            ycoor(nn)=log(grandcorrfcn(k)/grandcorrfcn(NRini));
            nn=nn+1;
        end
    end
    A0=[-1 1];
    [Afit Residual]= lsqcurvefit(linfunc,A0,xcoor,ycoor,[-4 -10],[4 10],options);
    GoFpow(ii)=Residual
    SlopeCorrPow(ii)=abs(Afit(1));
    
    %%%
    
    ii=ii+1;
end

%% plot figures

figure(1)
plot(G_range,err_grand);
figure(2)
plot(G_range,err_hete);
figure(3)
plot(G_range,fcfitt);
figure(4)
plot(G_range,fcfitt_hete);
figure(5)
plot(G_range,segint);
figure(6)
plot(G_range,GoFpow);
figure(7)
plot(G_range,SlopeCorrPow);
figure(8)
plot(G_range,infocapacity*5);
figure(9)
plot(G_range,susceptibility);
figure(10)
plot(G_range,Rmeta*10);
figure(11)
plot(G_range,DTmeta);

figure(12)
shadedErrorBar(G_range, asegintmean, asegintstd./sqrt(NSUBSIM))
figure(13)
shadedErrorBar(G_range, afcfittmean, afcfittstd./sqrt(NSUBSIM))

figure(14)
shadedErrorBar(G_range, infocapacity*5, (infocapacitystd./sqrt(NSUBSIM))*5)
figure(15)
shadedErrorBar(G_range, susceptibility, (susceptibilitystd./sqrt(NSUBSIM)))

figure(16)
shadedErrorBar(G_range, Rmeta, Rmeta./sqrt(NSUBSIM))
figure(17)
shadedErrorBar(G_range, DTmeta, DTmeta./sqrt(NSUBSIM))
