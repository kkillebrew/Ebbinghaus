clear
close all

%  file_list={'BKN_1_revamp','DNS_1_revamp'};
% file_list={'DEC_1_qtangles1'};
 file_list={'LXL_1_revampHORRM','TMD_1_revampHORRM','UMD_1_revampHORRM','KEN_1_revampHORRM'};
%file_list={'S3_1_revampHORRM'};

growth_list = linspace(-1.2,1.2,18);
final_list=zeros(5,18,length(file_list));
xaxis=growth_list;
xx_axis=-1.2:.01:1.2;
color = {'g','m','r','b','c'};

for a=1:length(file_list)
    
    load(['../data/' file_list{a}]);
    
    for i=1:length(rawdata)        
        for j=1:18
            if rawdata(i,4)==growth_list(j)
                rawdata(i,4)=j;
            end
        end
        
    end
    
    pres_list=zeros(5,18);
    shrink_list=zeros(5,18);
    
    
    for i=1:length(rawdata)
        for j=1:18
            if rawdata(i,4)==j
                
                pres_list(rawdata(i,1),j)=pres_list(rawdata(i,1),j)+1;
                
                if rawdata(i,5)==1
                    shrink_list(rawdata(i,1),j)=shrink_list(rawdata(i,1),j)+1;
                end
                
            end
        end
    end
    
    final_list(:,:,a)=shrink_list./pres_list;
    
    figure(a)
    for j=1:5
        datafit(:,:,j) = [shrink_list(j,:)',pres_list(j,:)']; %        datafit(:,:,j) = [final_list(j,:,a)',ones(18,1)];
        b(a,:,j) = glmfit(xaxis',datafit(:,:,j),'binomial','logit');
        PSE(a,j) = -b(a,1,j)/b(a,2,j);
        fitdata(j,:,a) = 100* exp(b(a,1,j)+b(a,2,j)*xx_axis')./(1+exp(b(a,1,j)+b(a,2,j)*xx_axis'));
        plot(xaxis,100*final_list(j,:,a)',color{j},'LineWidth',2);
        hold on
        plot(xx_axis,fitdata(j,:,a),color{j});
        plot(xaxis,50*ones(length(xaxis),1),'r--','LineWidth',2);
        set(gca,'ylim',[0,100]);
        set(gca,'xlim',[-1.2,1.2]);
    end
    
    figure(length(file_list)+a)
    for j = 1:5
        bar(j,PSE(a,j),color{j});
        hold on;
        plot(1:5,120*ones(5,1),'r--','LineWidth',2);
        set(gca,'ylim',[0,.5]);
    end
    
end


mean_PSE = mean(PSE);
stderr_PSE = std(PSE)/sqrt(length(file_list)-1);
mean_fit =mean(fitdata,3);
mean_results = mean(final_list,3);
figure(a+a+1)
for j = 1:5
    plot(xaxis,100*mean_results(j,:)',color{j},'LineWidth',2);
    hold on
    plot(xx_axis,mean_fit(j,:),color{j});
    plot(xaxis,50*ones(length(xaxis),1),'r--','LineWidth',2);
    set(gca,'ylim',[0,100]);
    set(gca,'xlim',[-1.2,1.2]);
    
end
figure(a+a+2)
for j = 1:5
    bar(j,mean_PSE(j),color{j});
    hold on;
    errorbar(j,mean_PSE(j),stderr_PSE(j),'k.');
    plot(1:5,120*ones(5,1),'r--','LineWidth',2);
    set(gca,'ylim',[0,.5]);
end