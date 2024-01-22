function [mu,sem] = barWithError(data)
%this function makes bar graphs with SEM error bars
%INPUTS =   data - your data matrix, with variables to plot in columns and
%                  observations/replicates in rows
%OUTPUTS =  mu - mean
%           sem - standard error
%           bar plot
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%
    mu = mean(data);
    coords = length(data(1,:));
    sem = std(data)/sqrt(length(data(:,1))-1);
    
    bar(mu)
    hold on
    errorbar(1:coords,mu,sem)