function [] = plot2metrics(data1,data2,titles,groupNames,groupMembers,ylab1,ylab2,ylim1,ylim2)

nGroups = length(groupNames);
figure
for n=1:length(groupNames)
    subplot(2,nGroups,n)
    barWithError(data1(:,groupMembers{n})')
    ylabel(ylab1)
    title(groupNames{n})
    ylim(ylim1)
    xticklabels(titles)
end
for n=1:length(groupNames)
    subplot(2,nGroups,n+nGroups)
    barWithError(data2(:,groupMembers{n})')
    ylabel(ylab2)
    title(groupNames{n})
    ylim(ylim2)
    xticklabels(titles)
end
end