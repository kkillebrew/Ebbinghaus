
clear
close all
rev_num=5;
stair_num=3;
% file_list={'BKN_1_revamp','DNS_1_revamp'};
% file_list={'DEC_1_qtangles1'};
  file_list={'S2_1_revampHORRM','LXL_1_revampHORRM','DBD_1_revampHOR','TMD_1_revampHORRM','UMD_1_revampHORRM','KEN_1_revampHORRM'};
%  file_list={'KEN_1_revampHORRM'};
revmeans=zeros(length(file_list),5);
for a=1:length(file_list)
    load(file_list{a});
    for i=1:length(rawdata)
        switch rawdata(i,1)
            case -9
                rawdata(i,1)=1;
            case 10
                rawdata(i,1)=2;
            case 11
                rawdata(i,1)=3;
            case 12
                rawdata(i,1)=4;
            case 13
                rawdata(i,1)=5;
                case 0
                rawdata(i,1)=1;
            case 22.5
                rawdata(i,1)=2;
            case 45
                rawdata(i,1)=3;
            case 67.5
                rawdata(i,1)=4;
            case 90
                rawdata(i,1)=5;
        end
    end
    
    thereversal_list=[];
    for i=1:5
        for j=1:2
            for k=1:stair_num
                z=1;
                previous=0;
                for l=1:length(rawdata)
                    if rawdata(l,1)==i && rawdata(l,2)==j && rawdata(l,3)==k
                        if previous~=0;
                            if rawdata(l,5)~=previous
                                thereversal_list(i,j,k,z)=rawdata(l,4);
                                z=z+1;
                            end
                        end
                        previous=rawdata(l,5);
                    end
                    
                end
            end
        end
    end
    
    end_rev_list=zeros(5,2,stair_num,rev_num);
    for i=1:5
        for j=1:2
            for k=1:2
                z=1;
                for l=length(thereversal_list(i,j,k,:)):-1:1
                    if thereversal_list(i,j,k,l)~=0
                    end_rev_list(i,j,k,z)=thereversal_list(i,j,k,l);
                    z=z+1;
                    end
                    if z>rev_num
                        break
                    end
                end
            end
        end
    end


 counter=zeros(5,2,stair_num);
    
    for i=1:5
        for j=1:2
            for k=1:stair_num
                for l=1:length(end_rev_list(i,j,k,:))
                    if end_rev_list(i,j,k,l)~=0
                        counter(i,j,k)=counter(i,j,k)+1;
                    end
                end
            end
        end
    end
    
    
    if stair_num==2;
        for i=1:5
            revmeans(a,i)=(sum(end_rev_list(i,1,1,:))+sum(end_rev_list(i,1,2,:))+...
                sum(end_rev_list(i,2,1,:))+sum(end_rev_list(i,2,2,:)))...
                /(counter(i,1,1)+counter(i,1,2)+counter(i,2,1)+counter(i,2,2));
        end
    else
        for i=1:5
            revmeans(a,i)=(sum(end_rev_list(i,1,1,:))+sum(end_rev_list(i,1,2,:))+sum(end_rev_list(i,1,3,:))+...
                sum(end_rev_list(i,2,1,:))+sum(end_rev_list(i,2,2,:))+sum(end_rev_list(i,2,3,:)))...
                /(counter(i,1,1)+counter(i,1,2)+counter(i,1,3)+counter(i,2,1)+counter(i,2,2)+counter(i,2,3));
        end
    end

    revmeans(a,:)=revmeans(a,:)*(1/max(revmeans(a,:)));
    figure(a)
    bar(revmeans(a,:))
end
ste_array=std(revmeans)/sqrt(length(file_list)-1);
figure(a+1)
bar(mean(revmeans));
hold on
errorbar(mean(revmeans),ste_array,'k.');


