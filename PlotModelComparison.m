function PlotModelComparison(sumprobability1,sumprobability2)
% This function takes the output of RunModelComparison and produces a
% contour plot visualization. The axes are the probabilities of enabling and
% disabling change and the contours are the log ratio of probabilities.
perange=10.^([-6:.1:-2]');
pdrange=perange;
pemat=perange*ones(1,length(pdrange));
pdmat=ones(length(perange),1)*(pdrange');
figure
contourf(pemat,pdmat,log2(sumprobability1./sumprobability2),400,'LineStyle','None')
set(gca,'xScale','log','yScale','log','TickLength',[.025 .025],'LineWidth',3,'FontSize',18,'xTick',10.^[-7:1:-2],'yTick',10.^[-7:1:-2]);
xlabel('Probability of enabling change','FontSize',24);
ylabel('Probability of disabling change','FontSize',24);
c=colorbar('FontSize',18);
c.Label.String='log_2 ratio probability pathway1/pathway2';
colormat=[.6 0 0; 
    .7 0 0;
    .8 0 0;
    .9 0 0;
    1 0 0;
    1 .1 .1;
    1 .2 .2;
    1 .3 .3;
    1 .4 .4;
    1 .5 .5;
    1 .6 .6;
    1 .7 .7;
    1 .8 .8;
    1 .9 .9;
    1 1 1;
    .9 .9 1;
    .8 .8 1;
    .7 .7 1;
    .6 .6 1;
    .5 .5 1;
    .4 .4 1;
    .3 .3 1;
    .2 .2 1;
    .1 .1 1;
    0 0 1;
    0 0 .9;
    0 0 .8;
    0 0 .7;
    0 0 .6];
colormap(colormat);
caxis([-2 2]);
axis square
end
    