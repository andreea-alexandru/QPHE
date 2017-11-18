filename = 'Data/512_64_results_30.txt';
C = strsplit(filename,{'/','_','.'});
key = C(2);
msg = C(3);
K = C(5);
samples = load(filename);
n = samples(:,1);
m = samples(:,2);
t_on = samples(:,3);
t_off= samples(:,4);

h = figure('rend','painters','pos',[10 10 900 600]);
scatter(n,t_on,70,m,'filled');
xlim([0,90]);ylim([0,115]);
xlabel('Number of variables')
ylabel('Online time [sec]')
title(strcat(K,' iterations, ',msg,' bits, ', key,' key-size'))
grid on
colormap winter
c = colorbar;
c.Label.String = 'Number of constraints';
name = strcat('results_',msg,'_',K,'_',key,'_on');
saveas(h,name{1},'epsc');

h = figure('rend','painters','pos',[10 10 900 600]);
scatter(n,t_off,70,m,'filled');
xlim([0,90]);ylim([0,2.5]);
xlabel('Number of variables')
ylabel('Offline time [sec]')
title(strcat(K,' iterations, ',msg,' bits, ', key,' key-size'))
grid on
colormap winter
c = colorbar;
c.Label.String = 'Number of constraints';
name = strcat('results_',msg,'_',K,'_',key,'_off');
saveas(h,name{1},'epsc');