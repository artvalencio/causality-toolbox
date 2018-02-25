function [results] = camir(x,y,tau)
% CAMIR Calculator of Causal Mutual Information Rate
%--------------------------------------------------------
% Limitations: 
%                (i) time-series in excess of 10^8 points may experience
%                issues
%--------------------------------------------------------
% Inputs:
%           -x,y: time-series, provide as single columns of the data
%           values, without timestamp
%           -tau: delay
% Assumption:
%	    - initial partition division is binary, 
%       - location of the division is allowed to change until find the optimum         
%--------------------------------------------------------
% Outputs:
%           -camirval_xy: Maximal Causal Mutual Information Rate X->Y
%           -camirval_xy: Maximal Causal Mutual Information Rate Y->X
%           -optimal_partition_xy: Partition maximizing Causal Mutual
%            Information Rate X->Y
%           -optimal_partition_yx: Partition maximizing Causal Mutual
%            Information Rate Y->X
%           -camir_xy: Causal Mutual Information Rate X->Y for each L (symbolic
%            sequence length) and position of partition division lines
%           -camir_yx: Causal Mutual Information Rate Y->X for each L (symbolic
%            sequence length) and position of partition division lines
%           -cami_xy: Causal Mutual Information X->Y for each L (symbolic
%            sequence length) and position of partition division lines
%           -cami_yx: Causal Mutual Information Y->X for each L and
%            partition divison location
%           -mi: Mutual Information of the of X and Y symbolic sequences
%            (I_L(X;Y)), for each sequence length L and partition division
%            location
%           -diridx: Directionality Index (positive for X->Y) for each L 
%            and partition divisons
%           -te_xy: Transfer Entropy X->Y for each L and partition divisons
%           -te_yx: Transfer Entropy Y->X for each L and partition divisons
%--------------------------------------------------------
% (C) Arthur Valencio* and Dr Murilo S Baptista 25 Feb 2018
%     ICSMB, University of Aberdeen
%     * Support: CNPq (Brazil) 
       
    %Calculate CaMI for varying partition location and symbolic length
    for i=0.1:0.1:0.9
        for j=0.1:0.1:0.9
            for L=1:1:4
                try
                    lx=L;
                    ly=2*L;
                    [tcami_xy,tcami_yx,tmi_past,tdiridx,tte_xy,tte_yx]=cami(x,y,lx,ly,i,j,tau,0);
                    a=floor(10*i);
                    b=floor(10*j);
                    cami_xy(a,b,L)=tcami_xy;
                    cami_yx(a,b,L)=tcami_yx;
                    mi(a,b,L)=tmi_past;
                    diridx(a,b,L)=tdiridx;
                    te_xy(a,b,L)=tte_xy;
                    te_yx(a,b,L)=tte_yx;
                catch
                    break;
                end
            end
        end
    end
    save('cami_xy.mat','cami_xy')
    save('mi.mat','mi')
    
    %Plot MI
    figure(1)
    for i=0.1:0.1:0.9
        for j=0.1:0.1:0.9
            data=[];
            a=floor(10*i);
            b=floor(10*j);
            maxlen=length(mi(a,b,:));
            data(1:maxlen)=mi(a,b,:);
            plot(1:maxlen,data);
            hold on;
            [posy,posx]=max(mi(a,b,:));
            text(posx,posy,strcat('[',num2str(i),',',num2str(j),']'));
        end
    end
    hold off;
    xlabel('L');
    ylabel('Mutual Information');
    

    %Plot CaMI
    figure(2)
    for i=0.1:0.1:0.9
        for j=0.1:0.1:0.9
            data=[];
            a=floor(10*i);
            b=floor(10*j);
            maxlen=length(cami_xy(a,b,:));
            data(1:maxlen)=cami_xy(a,b,:);
            plot(1:maxlen,data);
            hold on;
            [posy,posx]=max(cami_xy(a,b,:));
            text(posx,posy,strcat('[',num2str(i),',',num2str(j),']'));
        end
    end
    hold off;
    xlabel('L');
    ylabel('CaMI_{X\rightarrow Y}');
    
    %Get CaMIR plot from derivative

    for i=0.1:0.1:0.9
        for j=0.1:0.1:0.9
            a=floor(10*i);
            b=floor(10*j);
            tempxy(1:length(cami_xy(a,b,:)))=cami_xy(a,b,:);
            tempyx(1:length(cami_yx(a,b,:)))=cami_yx(a,b,:);
            coefsxy=polyfit(1:length(cami_xy(a,b,:)),tempxy,1);
            coefsyx=polyfit(1:length(cami_yx(a,b,:)),tempyx,1);
            camir_xy(floor(10*i),floor(10*j),:)=coefsxy(1);
            camir_yx(floor(10*i),floor(10*j),:)=coefsyx(1);
        end
    end
    
    [camirval_xy]=max(camir_xy(:));
    [row,column]=find(camir_xy == camirval_xy);
    optimal_partition_xy=[row,column]./10;
    [camirval_yx]=max(camir_yx(:));
    [row2,column2]=find(camir_yx == camirval_yx);
    optimal_partition_yx=[row2,column2]./10;
    
    disp('Maximal CaMIR X->Y:')
    disp(camirval_xy)
    
    results.camirval_xy=camirval_xy;
    results.camirval_yx=camirval_yx;
    results.optimal_partition_xy=optimal_partition_xy;
    results.optimal_partition_yx=optimal_partition_yx;
    results.camir_xy=camir_xy;
    results.camir_yx=camir_yx;
    results.cami_xy=cami_xy;
    results.cami_yx=cami_yx;
    results.mi=mi;
    results.diridx=diridx;
    results.te_xy=te_xy;
    results.te_yx=te_yx;
    
end
