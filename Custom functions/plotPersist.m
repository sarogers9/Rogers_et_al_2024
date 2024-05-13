function [] = plotPersist(dat,groups,groupNames,conditions,labels)
%%
colors = {'r','b','k'};
xs = {'Hab','Acq','Ext1','Ext2','Ext3'};
figure
for g = 1:length(groups)
for r=1:3
    subplot(3,length(groups),g+length(groups)*(r-1))
    for dir = 1:length(conditions)
        data = dat(groups{g},:,r,dir);
        hold on
        errorbar(1:length(mean(data)),mean(data),std(data)./sqrt(length(data)-1),colors{dir},'Linewidth',3)
    end
    if r<2
        title(groupNames{g})
    elseif r==3
        xticks([1 2 3 4 5])
    end
    if isempty(labels)
        ylabel(conditions{r})
    else
        ylabel(labels{r})
    end
    xlabel('Number of Sessions')
    ylim([0 .8])

end
end
end