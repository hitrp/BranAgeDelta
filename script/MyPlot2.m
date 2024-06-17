hold off; dscatter(Xplot,Yplot); hold on;   
if LL==1
  p1=plot([-L L],[-L L],'k','LineWidth',1);
else
  p1=plot([-L L],[0 0],'k','LineWidth',1);
end
axis equal; set(gca,'XLim',[-L L],'YLim',[-L L]);
title(MyTitle);
sub_pos = get(gca,'position'); % get subplot axis position
xf=0.95; yf=1.0;
set(gca,'position',sub_pos.*[1/xf 1/yf xf yf]) % stretch subplot width and height
