FigSize=[0 0 600 350]; MyStart

Yquad=0.05; N=20000;  Y=(rand(N,1)-0.5)+0.02*randn(N,1);  Y=60+Y*25;
deltaTRUE=2*randn(N,1);  
Y2=nets_demean(Yquad*nets_demean(Y).^2); Y2=Y2-nets_demean(Y)*(pinv(nets_demean(Y))*Y2);
Yb=Y+deltaTRUE+Y2;  Y=nets_demean(Y);
NOISE=0.2; ND=5;
X0=nets_normalise([Yb randn(N,ND-1)]);  Xmix=randn(ND,100).^5;  X1=X0*Xmix; X2=nets_normalise(X1);
X=nets_demean(X2+NOISE*randn(size(X2))); X=nets_svds(X,10);

beta1=pinv(X)*Y;  page=X*beta1;  delta1=page-Y;   beta2=pinv(Y)*delta1;  delta2=delta1-Y*beta2;  page2=Y+delta2;
YY=nets_demean([Y Y.^2]); YY=nets_demean([Y YY(:,2)-YY(:,1)*(pinv(YY(:,1))*YY(:,2))]);
beta1=pinv(X)*Y;  page=X*beta1;  delta1q=page-Y; beta2=pinv(YY)*delta1q;  delta2q=delta1q-YY*beta2; page2q=Y+delta2q;

Xplot=Y; Yplot=page; MyTitle='A.\rm  Predicted age \itY_{B1}\rm vs. age \itY'; 
L=20; LL=1; subplot(2,3,1);    MyPlot2;

Xplot=Y; Yplot=delta1; MyTitle='B.\rm  Initial delta \it\delta_1\rm vs. age \itY';
L=15; LL=2; subplot(2,3,2);  MyPlot2;

Xplot=deltaTRUE; Yplot=delta1; MyTitle='C.\rm  Initial delta \it\delta_1\rm vs. true delta';
L=10; LL=1; subplot(2,3,3);  MyPlot2;

Xplot=Y; Yplot=page2q; MyTitle='D.\rm  Predicted age \itY_{B2q}\rm vs. age \itY'; 
L=20; LL=1; subplot(2,3,4);    MyPlot2;

Xplot=Y; Yplot=delta2q; MyTitle='E.\rm  Corrected delta \it\delta_{2q}\rm vs. age \itY';
L=15; LL=2; subplot(2,3,5);  MyPlot2;

Xplot=deltaTRUE; Yplot=delta2q; MyTitle='F.\rm  Corrected delta \it\delta_{2q}\rm vs. true delta';
L=10; LL=1; subplot(2,3,6);  MyPlot2;

OutFile='QuadraticSimulations'; MyStop;
