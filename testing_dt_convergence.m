clear
close all
clc

%% Data

dt = [0.1, 0.05, 0.01, 0.005,0.001];

%for-loop implementation
h1 = [9.99465275695814626,...
    9.99455322465442109,...
    9.99446128138663070,...
    9.99444962527489444,...
    9.99443814714951451];

h2 = [9.96644518845135607,...
    9.96433435544512669,...
    9.96244104144809484,...
    9.96218944533607953,...
    9.96193403608311634];

%% Evaluation and plotting

thing = @(h) abs(h-h(end));
h1 = thing(h1);
h2 = thing(h2);

figure()
loglog(dt,h1,'-x')
dt = log(dt);
h1 = log(h1)
order = fit(dt(1:end-1)',h1(1:end-1)','poly1')