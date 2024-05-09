
clear

% load('DEC_1_5contets');
% load('BKN_1_revamp');

file_list={'BKN_1_revamp','DNS_1_revamp' 'RM_revamp5'};
%file_list={'RM_qtangles'};

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
        end
    end
    
    reversal_list=[];
    for i=1:5
        for j=1:2
            for k=1:2
                z=1;
                previous=0;
                for l=1:length(rawdata)
                    if rawdata(l,1)==i && rawdata(l,2)==j && rawdata(l,3)==k
                        if previous~=0;
                            if rawdata(l,5)~=previous
                                reversal_list(i,j,k,z)=rawdata(l,4);
                                z=z+1;
                            end
                        end
                        previous=rawdata(l,5);
                    end
                    
                end
            end
        end
    end
    
    counter=zeros(5,2,2);
    
    for i=1:5
        for j=1:2
            for k=1:2
                for l=1:length(reversal_list(i,j,k,:))
                    if reversal_list(i,j,k,l)~=0
                        counter(i,j,k)=counter(i,j,k)+1;
                    end
                end
            end
        end
    end
    
    
    for i=1:5
        revmeans(a,i)=(sum(reversal_list(i,1,1,:))+sum(reversal_list(i,1,2,:))+sum(reversal_list(i,2,1,:))+sum(reversal_list(i,2,2,:)))...
            /(counter(i,1,1)+counter(i,1,2)+counter(i,2,1)+counter(i,2,2));
    end
    
    figure(a)
    bar(revmeans(a,:))
end
% ste_array=std(revmeans)/sqrt(length(file_list)-1);
figure(a+1)
bar(mean(revmeans));


