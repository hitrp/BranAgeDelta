
function [delta1,delta2,delta3,delta2q,delta3q] = DeltaCV(X,Y,J,CONF);

delta1=0*Y/0; delta2=delta1; delta3=delta1; delta2q=delta1; delta3q=delta1;
Ncv=10;                    % number of cross-validation folds
I=randi(Ncv,length(Y),1);  % create random cross-validation folds

for i=1:Ncv
  OUT=(I==i);  IN=(I~=i);

  x=X(IN,:);  y=Y(IN);  ym=mean(y);  y=y-ym;
  yy=nets_demean([y y.^2]);  yy=nets_demean([y yy(:,2)-yy(:,1)*(pinv(yy(:,1))*yy(:,2))]);
  if size(CONF,1)>0,  CONFbeta=pinv(CONF(IN,:))*x; x=x-CONF(IN,:)*CONFbeta;  end;
  if J>0,  [x,pcaS,pcaV]=nets_svds(x,J);  end;
  beta1=pinv(x)*y; deltaIN=x*beta1-y;  beta2=pinv(y)*deltaIN;  beta2q=pinv(yy)*deltaIN;
  gamma1=pinv(y)*x;  gamma2=pinv(yy)*x;

  x=X(OUT,:);  y=Y(OUT)-ym;
  yy=nets_demean([y y.^2]);  yy=nets_demean([y yy(:,2)-yy(:,1)*(pinv(yy(:,1))*yy(:,2))]);
  if size(CONF,1)>0,  x=x-CONF(OUT,:)*CONFbeta;  end;
  if J>0, x=x*pcaV*inv(pcaS);  end;
  delta1(OUT)=x*beta1-y;  delta2(OUT)=delta1(OUT)-y*beta2;  delta2q(OUT)=delta1(OUT)-yy*beta2q;
  d=x-y*gamma1; delta3(OUT)=d*pinv(gamma1);
  d=x-yy*gamma2; delta3q(OUT)=d*pinv(gamma2(1,:));
end

% for slightly improved results, run this many times and average the outputs

