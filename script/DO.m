
addpath ~steve/NETWORKS/FSLNets
load workspace8.mat    % precomputed UK Biobank IDPs (etc.) data. See previous papers from FMRIB for more information

X = [ ALL_IDPs(:,IDP_modality_types(1:size(ALL_IDPs,2))>0) NODEamps25 NODEamps100 NET25 NET100 ];
X=X-repmat(nanmedian(X),size(X,1),1); X=X./repmat(nanmedian(abs(X)),size(X,1),1); X(abs(X)>6)=NaN;
X=nets_normalise(X);  Xgood=sum(isnan(X),2)<50;  X(isnan(X))=randn(size(X(isnan(X))))*0.01;
Y=age;
confX=conf(:,(size(confA1,2)+size(confA2,2)+1):end);
Kp=find((sum(isnan([Y confX]),2)==0)&(Xgood==1));  KpN=length(Kp);  Kp=Kp(randperm(KpN));  100*KpN/N
Y=Y(Kp);  X=nets_inormal(X(Kp,:)); confX=nets_normalise(confX(Kp,:));

%%% real data analyses for Tables 1 and 2
po1=fopen(sprintf('~steve/BrainAge1.txt'),'w');
po2=fopen(sprintf('~steve/BrainAge2.txt'),'w');
for J = [ 0 1 2 5 10 20 50 100 1000 2635 ]
  [delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X,Y,J,confX);
  page1=Y+delta1; page2=Y+delta2; page3=Y+delta3;
  [grotR,grotP]=corr([delta1 delta2 delta3 delta2q delta3q],vars_i_deconf(Kp,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99);
  %%% linear table
  fprintf(po1,...
    '%4d &%.2f &%.2f &%.2f &%.2f &%.1f &%.1f &%.1f &%.2f &%.2f &%.2f \\\\',...
    [J corr(Y,[delta1 page1 page2 page3]) mean(abs([delta1 delta2 delta3])) Q(1:3)]);
  %%% nonlinear table
  fprintf(po2,...
    '%4d &%.1f &%.1f &%.1f &%.1f &%.1f &%.2f &%.2f &%.2f &%.2f &%.2f \\\\',...
    [J mean(abs([delta1 delta2 delta3 delta2q delta3q])) Q]);
end
fclose(po1); fclose(po2);

%%% real data example figure, from running J=50.  For Fig 3
J=50; [delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X,Y,J,confX);
[grotR,grotP]=corr([delta1 delta2 delta3 delta2q delta3q],vars_i_deconf(Kp,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99);
[grotR0,grotP0]=corr([delta1 delta2 delta3 delta2q delta3q],vars_i_deconfnoage(Kp,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99);
grot=-log10([grotP0' grotP']);  grot(isinf(grot))=max(grot(~isinf(grot)));
subplot(3,3,1); dscatter(Y,delta1); hold on; plot([0 100],[0 0],'k');
set(gca,'Xlim',[45 80],'Ylim',[-20 20]); xlabel('Age'); ylabel('\delta_1');
subplot(3,3,2); dscatter(Y,delta2); hold on; plot([0 100],[0 0],'k');
set(gca,'Xlim',[45 80],'Ylim',[-20 20]); xlabel('Age'); ylabel('\delta_2');
subplot(3,3,3); hold off; dscatter(grot(:,6),grot(:,7)); hold on; plot([0 40],[0 40]);
set(gca,'Xlim',[0 40],'Ylim',[0 40]); xlabel('\delta_1  -log_{10}P'); ylabel('\delta_2  -log_{10}P');
subplot(3,3,4:6); hold off; dscatter(grot(:,1),grot(:,6)); hold on; plot([0 40],[0 40]);
set(gca,'Xlim',[0 315],'Ylim',[0 22]); xlabel('\delta_1  -log_{10}P  (no age deconfounding on non-imaging variables)'); ylabel('\delta_1  -log_{10}P');
subplot(3,3,7:9); hold off; dscatter(grot(:,1),grot(:,7)); hold on; plot([0 40],[0 40]);
set(gca,'Xlim',[0 315],'Ylim',[0 35]); xlabel('\delta_1  -log_{10}P  (no age deconfounding on non-imaging variables)'); ylabel('\delta_2  -log_{10}P');
set(gcf,'Position',[0 0 950 900],'PaperPositionMode','auto');
OutFile='UKB_P_scatter'; print('-dpdf',sprintf('/home/fs0/steve/%s',OutFile));
system(sprintf('cd /home/fs0/steve; pdfcrop %s.pdf %s.pdf ; tar cvfz ~/grot %s.pdf',OutFile,OutFile,OutFile));


%%% delta scaling (Section 2.8)
Ymin=min(Y); Ymax=max(Y); y0=(Y-Ymin)/(Ymax-Ymin);
lambda=exp( pinv(nets_demean(y0)) * nets_demean(log(abs(delta2q))) )-1;
lambda*mean(abs(delta2q))
[grot,grotstats]=robustfit(nets_demean(y0),nets_demean(log(abs(delta2q))),'ols',[])
grotstats.p

%%% correlations against nonimaging measures (Tables 3 and 4)

SEX=sex(Kp);  J=50;
[delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X,Y,J,confX);
[delta1F,delta2F,delta3F,delta2qF,delta3qF] = DeltaCV(X(SEX==0,:),Y(SEX==0,:),J,confX(SEX==0,:));
[delta1M,delta2M,delta3M,delta2qM,delta3qM] = DeltaCV(X(SEX==1,:),Y(SEX==1,:),J,confX(SEX==1,:));

grot0=delta2q; grot0(SEX==1)=NaN; grot1=delta2q; grot1(SEX==0)=NaN; nanmean(grot0)-nanmean(grot1)
hist([grot0 grot1],100);
legend(sprintf('female (mean=%.2fy std=%.1fy)',nanmean(grot0),nanstd(grot0)),sprintf('male (mean=%.2fy std=%.1fy)',nanmean(grot1),nanstd(grot1)));

[ corr(delta2q(SEX==0),delta2qF) corr(delta2q(SEX==1),delta2qM) ]

VarsDeconf=vars_i_deconf(Kp,:);
[grotR,grotP]=corr(delta2q,VarsDeconf,'rows','pairwise');
[grotRF1,grotPF1]=corr(delta2q(SEX==0),VarsDeconf(SEX==0,:),'rows','pairwise');
[grotRM1,grotPM1]=corr(delta2q(SEX==1),VarsDeconf(SEX==1,:),'rows','pairwise');
[grotRF2,grotPF2]=corr(delta2qF,VarsDeconf(SEX==0,:),'rows','pairwise');
[grotRM2,grotPM2]=corr(delta2qM,VarsDeconf(SEX==1,:),'rows','pairwise');
grotLOGP=-log10([grotP' grotPF1' grotPM1' grotPF2' grotPM2']);
corr(grotLOGP,'rows','pairwise')

grotPP=grotP'; grotRR=grotR';
%grotPP=grotPF2'; grotRR=grotRF2';
%grotPP=grotPM2'; grotRR=grotRM2';
[grotY,grotI]=sort(-log10(grotPP),'descend'); 
for i= grotI(grotY>8)'
  redstring=''; if grotRR(i)>0, redstring='\RED '; end;
  disp(sprintf('%2.1f &%.2f &%s &%s%s',-log10(grotPP(i)),grotRR(i),varsVARS{varskeep(i)},redstring,varsHeader{varskeep(i)}));
end

%%% now for IDPs (Table 5)
XX=ALL_IDPs_i_deconf(Kp,:); XN=IDPnames(IDP_modality_types>0);
[grotR,grotP]=corr(delta2q,XX,'rows','pairwise'); Q=prctile(-log10(grotP'),99)
[grotRF1,grotPF1]=corr(delta2q(SEX==0),XX(SEX==0,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99)
[grotRM1,grotPM1]=corr(delta2q(SEX==1),XX(SEX==1,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99)
[grotRF2,grotPF2]=corr(delta2qF,XX(SEX==0,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99)
[grotRM2,grotPM2]=corr(delta2qM,XX(SEX==1,:),'rows','pairwise'); Q=prctile(-log10(grotP'),99)
grotLOGP=-log10([grotP' grotPF1' grotPM1' grotPF2' grotPM2']);  grotLOGP(isinf(grotLOGP)) = max(grotLOGP(~isinf(grotLOGP)));
corr(grotLOGP,'rows','pairwise')

[grotY,grotI]=sort(abs(grotR'),'descend'); 
for i= grotI'
  if max(max(abs(grotR(i)),abs(grotRF2(i))),abs(grotRM2(i))) > 0.3
    redstring=''; if abs(grotRF2(i)-grotRM2(i))>0.05, redstring='\RED '; end;
    grot=strrep(XN{i},'_',' ');
    disp(sprintf('%.2f &%.2f &%.2f  &%s%s \\\\',grotR(i),grotRF2(i),grotRM2(i),redstring,grot));
  end
end


%%%%%%%%%%%%% non-additive pre-simulation (Section 2.8)
poop=[];
for lambda=[-1:0.01:1]
  N=20000; Y=rand(N,1)+0.1; % Y=randn(N,1)*10;  % Y=(rand(N,1)-0.5)+0.02*randn(N,1);  Y=60+Y*25;
  minY=min(Y); maxY=max(Y); Y=(Y-minY)/(maxY-minY); % rescale Y from 0:1
  d=randn(N,1);  dd=d.*(1+lambda*Y);
  grot3= pinv(nets_demean(Y)) * nets_demean(log(dd.^2)/2);
  grot4= exp(grot3)-1;
  grot6 = grot3+ (lambda-log(1+lambda)); % cheating option
  poop=[poop;lambda  grot3  grot4 grot6]; poop(isnan(poop))=0;
end
X=[poop(:,3) poop(:,3).^2 ones(length(poop(:,1)),1)]; Xbeta=pinv(X)*poop(:,1); grot0=X*Xbeta; poop=[poop grot0]; %cheating
plot(poop(:,1),poop); corr(poop)



%%% linear simulations (Table 1)
clear all;
Nsimulations=20; for sim=1:Nsimulations, j=1, sim
  N=10000;  Y=(rand(N,1)-0.5)+0.02*randn(N,1);  Y=60+Y*25;  deltaTRUE=2*randn(N,1);
  minY=min(Y); maxY=max(Y); Y0=(Y-minY)/(maxY-minY); 
  NOISE=0.5; ND=100;  %%% SIM 1
  %NOISE=10.0; ND=1;   %%% SIM 2
  %NOISE=0.5; ND=100;  deltaTRUE=deltaTRUE.*(1 + 0.5*Y0);         %%% SIM 3 non-additive delta
  Yb=Y+deltaTRUE;
  X0=nets_normalise([Yb randn(N,ND-1)]);  Xmix=randn(ND,3000).^5;  X1=X0*Xmix; X2=nets_normalise(X1);
  X3=nets_demean(X2+NOISE*randn(size(X2)));
  for J = [ 0 1 10 50 100 1000 2990 ]    % permuting pcaU gives same results as J=1 (ie null model ~ bad model)
    [delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X3,Y,J,[]);
    page1=Y+delta1; page2=Y+delta2; page3=Y+delta3;
    % acosd(dot(Y,-delta1)/(norm(Y)*norm(delta1)))
    %disp(sprintf('%4d  %.2f %.2f %.2f %.2f  %.2f %.2f %.2f  %.2f %.2f %.2f',J,corr(Y,[delta1 page1 page2 page3]), ...
    %             mean(abs([delta1 delta2 delta3])),corr(deltaTRUE,[delta1 delta2 delta3])));
    AllResults(j,:,sim)=[J corr(Y,[delta1 page1 page2 page3]) mean(abs([delta1 delta2 delta3])) corr(deltaTRUE,[delta1 delta2 delta3])]; j=j+1;

    %%%% non-additive delta results Sim3
    %poop=exp( pinv(nets_demean(Y0)) * nets_demean(log(abs([deltaTRUE delta1 delta2 delta3]))) )-1;
    %disp(sprintf('%.2f %.2f %.2f %.2f',poop(1,:)));
  end
end
grot1=mean(AllResults,3); grot2=std(AllResults,[],3); 
for j = 1:size(AllResults,1)
  grot=[grot1(j,2:end); grot2(j,2:end)]; grot=grot(:)';
  disp(sprintf(...
   '%4d &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} \\\\',...
   grot1(j,1),grot));
end


%%% quadratic simulations (Table 2)
clear all;
Nsimulations=20; for sim=1:Nsimulations, jj=1, sim
for Yquad=[0 0.01 0.025]
  j=1;
  N=10000;  Y=(rand(N,1)-0.5)+0.02*randn(N,1);  Y=60+Y*25;  deltaTRUE=2*randn(N,1);
  Y2=nets_demean(Yquad*nets_demean(Y).^2); Y2=Y2-nets_demean(Y)*(pinv(nets_demean(Y))*Y2);
  Yb=Y+deltaTRUE+Y2;
  NOISE=0.5; ND=100;
  X0=nets_normalise([Yb randn(N,ND-1)]);  Xmix=randn(ND,3000).^5;  X1=X0*Xmix; X2=nets_normalise(X1);
  X3=nets_demean(X2+NOISE*randn(size(X2)));
  for J = [ 0 10 50 100 1000 2990 ]
    [delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X3,Y,J,[]);
    page1=Y+delta1; page2=Y+delta2; page3=Y+delta3;
    AllResults(j,:,sim,jj)=[J mean(abs([delta1 delta2 delta3   delta2q delta3q])) corr(deltaTRUE,[delta1  delta2 delta3  delta2q delta3q])]; j=j+1;
  end
  jj=jj+1;
end
end
for jj=1:size(AllResults,4)
  grot1=mean(AllResults(:,:,:,jj),3); grot2=std(AllResults(:,:,:,jj),[],3); 
  for j = 1:size(AllResults,1)
    grot=[grot1(j,2:end); grot2(j,2:end)]; grot=grot(:)';
    disp(sprintf(...
     '%4d &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.1f\\PM{%.1f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f} &%.2f\\PM{%.2f}\\\\',...
     grot1(j,1),grot));
  end
end

