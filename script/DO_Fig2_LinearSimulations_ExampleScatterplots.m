FigSize=[0 0 600 350]; MyStart

N=20000;  Y=(rand(N,1)-0.5)+0.01*randn(N,1);  Y=60+Y*25;  deltaTRUE=2*randn(N,1);  Yb=Y+deltaTRUE;
ND=5; X=nets_normalise([Yb randn(N,ND-1)]);
X=nets_normalise(X*randn(ND,100).^5);
X=nets_demean(X+0.2*randn(size(X)));
Y=nets_demean(Y);

beta1=pinv(X)*Y;  page=X*beta1;  delta1=page-Y;   beta2=pinv(Y)*delta1;  delta2=delta1-Y*beta2;  page2=Y+delta2;

Xplot=Y; Yplot=page; MyTitle='A.\rm  Predicted age \itY_{B1}\rm vs. age \itY'; 
L=20; LL=1; subplot(2,3,1);    MyPlot2;

Xplot=Y; Yplot=delta1; MyTitle='B.\rm  Initial delta \it\delta_1\rm vs. age \itY';
L=15; LL=2; subplot(2,3,2);  MyPlot2;

Xplot=deltaTRUE; Yplot=delta1; MyTitle='C.\rm  Initial delta \it\delta_1\rm vs. true delta';
L=10; LL=1; subplot(2,3,3);  MyPlot2;

Xplot=Y; Yplot=page2; MyTitle='D.\rm  Predicted age \itY_{B2}\rm vs. age \itY'; 
L=20; LL=1; subplot(2,3,4);    MyPlot2;

Xplot=Y; Yplot=delta2; MyTitle='E.\rm  Corrected delta \it\delta_2\rm vs. age \itY';
L=15; LL=2; subplot(2,3,5);  MyPlot2;

Xplot=deltaTRUE; Yplot=delta2; MyTitle='F.\rm  Corrected delta \it\delta_2\rm vs. true delta';
L=10; LL=1; subplot(2,3,6);  MyPlot2;

OutFile='LinearSimulations'; MyStop;

