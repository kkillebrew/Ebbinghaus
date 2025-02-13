clear
close all
% filename_list={'S1','S2','KEN','HXH'};
 filename_list={'S2','KEN','HXH'};
color = {'g','m','r','b'};

xaxis= [-.25,0,.15,.25,.35,.5,1];

xx_axis = -.25:.01:1;

for i=1:length(filename_list)
    load(sprintf('%s%s',filename_list{i},'_1_ThreeCon'));
    
    results(:,:,i) = done_data/20;
    figure(i)
    for j = 1:3
        datafit(:,:,j) = [results(j,:,i)',ones(7,1)];
        b(i,:,j) = glmfit(xaxis',datafit(:,:,j),'binomial','logit');
        PSE(i,j) = -b(i,1,j)/b(i,2,j);
        fitdata(j,:,i) = 100* exp(b(i,1,j)+b(i,2,j)*xx_axis')./(1+exp(b(i,1,j)+b(i,2,j)*xx_axis'));
        plot(xaxis,100*results(j,:,i)',color{j},'LineWidth',2);
        hold on
        plot(xx_axis,fitdata(j,:,i),color{j});
        plot(xaxis,50*ones(length(xaxis),1),'r--','LineWidth',2);
        set(gca,'ylim',[0,100]);
        set(gca,'xlim',[-.25,1]);
        
    end
    ylabel('% Perceived Shrinking');
    xlabel('Growth Factor');
    
end
mean_PSE = mean(PSE);
stderr_PSE = std(PSE)/sqrt(length(filename_list)-1);
mean_fit =mean(fitdata,3);
mean_results = mean(results,3);
figure(i+1)
for j = 1:3
    plot(xaxis,100*mean_results(j,:)',color{j},'LineWidth',2);
    hold on
    plot(xx_axis,mean_fit(j,:),color{j});
    plot(xaxis,50*ones(length(xaxis),1),'r--','LineWidth',2);
    set(gca,'ylim',[0,100]);
    set(gca,'xlim',[-.25,1]);
    
end
figure(i+2)
for j = 1:3
    bar(j,mean_PSE(j),color{j});
    hold on;
    errorbar(j,mean_PSE(j),stderr_PSE(j),'k.');
    plot(1:4,120*ones(4,1),'r--','LineWidth',2);
    set(gca,'ylim',[0,.5]);
end