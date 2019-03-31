function [h1,h2]=CoolHistogram(x1,x2,bins,varname)
% CoolHistogram plots a nice looking red and grey histogram for two group
% comparisons

edges=linspace(min([x1;x2]),max([x1;x2]),bins);
h1=histogram(x1,edges,'Normalization','probability');hold on
h2=histogram(x2,edges,'Normalization','probability');
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel(varname)
ylabel('Probability')
end