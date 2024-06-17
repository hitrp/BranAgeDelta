FigSize=[0 0 700 500]; MyStart

N=10000;               % number of samples
X0=randn(N,1);         % true something
Y0=X0+randn(N,1)*0.3;  % true something plus noise

X=X0; Y=Y0; beta=pinv(demean(X))*demean(Y);
MyTitle='A.  Gaussian A, B and noise'; L=4; subplot(2,3,1); MyPlot;

X=X0; Y=Y0+randn(N,1)*1; beta=pinv(demean(X))*demean(Y);
MyTitle='B.  Extra noise on B'; L=4; subplot(2,3,2); MyPlot;

Y=Y0(abs(X0)<1); X=X0(abs(X0)<1); beta=pinv(demean(X))*demean(Y);
MyTitle='C.  Truncation of A'; L=2; subplot(2,3,3); MyPlot;

Y=demean(Y0); X=demean(X0); beta=inv(X'*X + 5000)*X'*Y;
MyTitle='D.  Regularised fit'; L=4; subplot(2,3,4); MyPlot;

X=X0+randn(N,1)*1; Y=Y0; beta=pinv(demean(X))*demean(Y);
MyTitle='E.  Noise on A (regression dilution)'; L=4; subplot(2,3,5); MyPlot;

X=X0(abs(Y0)<0.75); Y=Y0(abs(Y0)<0.75); beta=pinv(demean(X))*demean(Y);
MyTitle='F.  Truncation of B'; L=2; subplot(2,3,6); MyPlot;

OutFile='BiasExamples'; MyStop;
