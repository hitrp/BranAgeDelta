sigma=sqrt(mean((Y-X*beta).^2));
hold off;  p1=plot([-L L],[-L L],'k','LineWidth',1);
hold on;  dscatter(X,Y);
p2=plot(0.8*[-L L],0.8*[-L L]*beta,'r','LineWidth',2);
xlabel('A'); ylabel('B');
axis equal; set(gca,'XLim',[-L L],'YLim',[-L L]);
title({MyTitle,sprintf('\\rm\\beta=%.2f    \\sigma=%.2f',beta,sigma)});
sub_pos = get(gca,'position'); % get subplot axis position
xf=1.02; yf=1.07;
set(gca,'position',sub_pos.*[1/xf 1/yf xf yf]) % stretch subplot width and height
