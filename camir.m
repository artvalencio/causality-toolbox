function [rates,fullresults] = camir(x,y,tau)
% CAMIR Calculator of Causal Mutual Information Rate
%--------------------------------------------------------
% Function inputs:
%           -x,y: time-series, provide as single columns of the data
%           values, without timestamp
%           -tau: delay
% Aditional inputs (dialog box):
%           - Partition division lines (minimuam, maximum, step):
%             The method will find the optimum partition based on data, but
%             the user needs to provide the minimum and maximum thresholds 
%             and the step to which the partition will be changed
%           - Max L: Maximum symbolic sequence length for initial
%           computation, dependent on the available memory/time.
%           Sugested between 5 and 10.
%           - L_min and L_max of linear part: simply select the most
%           adequate interval of L for your data.
%--------------------------------------------------------
% Assumption:
%	        - initial partition division is binary      
%--------------------------------------------------------
% Outputs:
%           Rates:
%           -camir_xy: Maximal Causal Mutual Information Rate X->Y
%           -camir_yx: Maximal Causal Mutual Information Rate Y->X
%           -mir: Mutual Information Rate in the optimal partition (from CaMIR X->Y)
%           -ter_xy: Transfer Entropy Rate X->Y in the optimal partition
%           -ter_yx: Transfer Entropy Rate Y->X in the optimal partition
%           -optimal_partition_xy: Optimal Partition in the direction X->Y;
%           -optimal_partition_yx: Optimal Partition in the direction Y->X;
%
%           Full results:
%           -camirval_xy: Maximal Causal Mutual Information Rate X->Y
%           -camirval_xy: Maximal Causal Mutual Information Rate Y->X
%           -camirval_xy: Mutual Information Rate in the optimal partition
%           -terval_xy: Transfer Entropy Rate X->Y in the optimal partition
%           -terval_xy: Transfer Entropy Rate Y->X in the optimal partition
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
% (C) Arthur Valencio* and Dr Murilo S Baptista, version update 17 May 2018
%     ICSMB, University of Aberdeen
%     * Support: CNPq (Brazil) [206246/2014-5]
       
    %prompts user for initial configuration
    prompt = {'Min X div line','Max X div line','X div line step',...
              'Min Y div line','Max Y div line','Y div line step',...
              'Max L'};
    title = 'Partition divisions';
    answer = inputdlg(prompt,title);
    xpartmin=str2num(answer{1});
    xpartmax=str2num(answer{2});
    xpartstep=str2num(answer{3});
    ypartmin=str2num(answer{4});
    ypartmax=str2num(answer{5});
    ypartstep=str2num(answer{6});
    Lmax=str2num(answer{7});
    
    for i=1:length(x)
        xerror(i)=rand;
        yerror(i)=rand;
    end
    
    %Calculation of CaMI, MI and TE for each L and each partition configuration
    for i=xpartmin:xpartstep:xpartmax
        for j=ypartmin:ypartstep:ypartmax
            for L=1:Lmax
                try
                    %Calculation of CaMI, MI and TE for each L
                    [tcami_xy,tcami_yx,tmi,tdiridx,tte_xy,tte_yx]=cami(x,y,L,L,i,j,tau,'bits');
                    a=floor(10*i);
                    b=floor(10*j);
                    cami_xy(a,b,L)=tcami_xy;
                    cami_yx(a,b,L)=tcami_yx;
                    mi(a,b,L)=tmi;
                    diridx(a,b,L)=tdiridx;
                    te_xy(a,b,L)=tte_xy;
                    te_yx(a,b,L)=tte_yx;
                    %Calculation of error levels of CaMI, MI and TE
                    [etcami_xy,etcami_yx,etmi,etdiridx,ette_xy,ette_yx]=cami(xerror,yerror,L,L,i,j,tau,'bits');
                    ea=floor(10*i);
                    eb=floor(10*j);
                    ecami_xy(a,b,L)=etcami_xy;
                    ecami_yx(a,b,L)=etcami_yx;
                    emi(a,b,L)=etmi;
                    ediridx(a,b,L)=etdiridx;
                    ete_xy(a,b,L)=ette_xy;
                    ete_yx(a,b,L)=ette_yx;
                catch
                    break;
                end
            end
        end
    end
        
    %Plot MI
    colorsize=length(xpartmin:xpartstep:xpartmax)*length(ypartmin:ypartstep:ypartmax);
    figure(1)
    colidx=0;
    for i=xpartmin:xpartstep:xpartmax
        for j=ypartmin:ypartstep:ypartmax         
            %colorcode for plot lines
            col=hsv(colorsize);
            colidx=colidx+1;
            %plot MI for each partition
            data=[];
            a=floor(10*i);
            b=floor(10*j);
            maxlen=length(mi(a,b,:));
            data(1:maxlen)=mi(a,b,:);
            plot(1:maxlen,data,'color',col(colidx,:));
            hold on;
            [posy,posx]=max(mi(a,b,:));
            text(posx,posy,strcat('[',num2str(i),',',num2str(j),']'));
            %plot MI errors for each partition
            edata=[];
            ea=floor(10*i);
            eb=floor(10*j);
            emaxlen=length(emi(a,b,:));
            edata(1:emaxlen)=emi(a,b,:);
            plot(1:emaxlen,edata,'--','color',col(colidx,:));
            [eposy,eposx]=max(emi(a,b,:));
            text(eposx,eposy,strcat('error: [',num2str(i),',',num2str(j),']'));            
        end
    end
    hold off;
    xlabel('L');
    ylabel('Mutual Information');    
    %prompts user for linear interval
    prompt = {'Linear interval starts at L=',' And ends at L='};
    title = 'Linear interval';
    answ = inputdlg(prompt,title);
    lmin=str2num(answ{1});
    lmax=str2num(answ{2}); 

    %Plot CaMI
    figure(2)
    colidx=0;
    for i=xpartmin:xpartstep:xpartmax
        for j=ypartmin:ypartstep:ypartmax
            %colorcode for plot lines
            col=hsv(colorsize);
            colidx=colidx+1;
            %Plot CaMI for each partition           
            data=[];
            a=floor(10*i);
            b=floor(10*j);
            maxlen=length(cami_xy(a,b,:));
            data(1:maxlen)=cami_xy(a,b,:);
            plot(1:maxlen,data,'color',col(colidx,:));
            hold on;
            [posy,posx]=max(cami_xy(a,b,:));
            text(posx,posy,strcat('[',num2str(i),',',num2str(j),']'));
            %Plot errors
            edata=[];
            ea=floor(10*i);
            eb=floor(10*j);
            emaxlen=length(ecami_xy(a,b,:));
            edata(1:emaxlen)=ecami_xy(a,b,:);
            plot(1:emaxlen,edata,'--','color',col(colidx,:));
            hold on;
            [eposy,eposx]=max(ecami_xy(a,b,:));
            text(eposx,eposy,strcat('error: [',num2str(i),',',num2str(j),']'));
        end
    end
    hold off;
    xlabel('L');
    ylabel('CaMI_{X\rightarrow Y}');
    %prompts user for linear interval
    prompt = {'Linear interval starts at L=',' And ends at L='};
    title = 'Linear interval';
    answ2 = inputdlg(prompt,title);
    lmin2=str2num(answ2{1});
    lmax2=str2num(answ2{2}); 
    
    %Plot TE
    figure(3)
    colidx=0;
    for i=xpartmin:xpartstep:xpartmax
        for j=ypartmin:ypartstep:ypartmax
            %colorcode for plot lines
            col=hsv(colorsize);
            colidx=colidx+1;
            %Plot TE
            data=[];
            a=floor(10*i);
            b=floor(10*j);
            maxlen=length(te_xy(a,b,:));
            data(1:maxlen)=te_xy(a,b,:);
            plot(1:maxlen,data,'color',col(colidx,:));
            hold on;
            [posy,posx]=max(te_xy(a,b,:));
            text(posx,posy,strcat('[',num2str(i),',',num2str(j),']'));
            %Plot errors
            edata=[];
            ea=floor(10*i);
            eb=floor(10*j);
            emaxlen=length(ete_xy(a,b,:));
            edata(1:emaxlen)=ete_xy(a,b,:);
            plot(1:emaxlen,edata,'--','color',col(colidx,:));
            hold on;
            [eposy,eposx]=max(ete_xy(a,b,:));
            text(eposx,eposy,strcat('error: [',num2str(i),',',num2str(j),']'));
        end
    end
    hold off;
    xlabel('L');
    ylabel('TE_{X\rightarrow Y}');
    %prompts user for linear interval
    prompt = {'Linear interval starts at L=',' And ends at L='};
    title = 'Linear interval';
    answ3 = inputdlg(prompt,title);
    lmin3=str2num(answ3{1});
    lmax3=str2num(answ3{2}); 
    
    %Get CaMIR
    for i=xpartmin:xpartstep:xpartmax
        for j=ypartmin:ypartstep:ypartmax
            a=floor(10*i);
            b=floor(10*j);
            tempxy(1:length(cami_xy(a,b,lmin2:lmax2)))=cami_xy(a,b,lmin2:lmax2);
            tempyx(1:length(cami_yx(a,b,lmin2:lmax2)))=cami_yx(a,b,lmin2:lmax2);
            coefsxy=polyfit(lmin2:lmax2,tempxy,1);
            coefsyx=polyfit(lmin2:lmax2,tempyx,1);
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
    %Get TER
    tempxy(1:length(cami_xy(row,column,lmin3:lmax3)))=te_xy(row,column,lmin3:lmax3);
    tempyx(1:length(cami_yx(row,column,lmin3:lmax3)))=te_yx(row,column,lmin3:lmax3);
    coefsxy=polyfit(lmin3:lmax3,tempxy,1);
    coefsyx=polyfit(lmin3:lmax3,tempyx,1);
    terval_xy=coefsxy(1);
    terval_yx=coefsyx(1);
    %Get MIR
    tempmi(1:length(mi(row,column,lmin:lmax)))=mi(a,b,lmin:lmax);
    coefs=polyfit(lmin:lmax,tempmi,1);
    mirval=coefs(1);   
    
    %Build rates output
    rates.mir=mirval;
    rates.camir_xy=camirval_xy;
    rates.camir_yx=camirval_yx;
    rates.ter_xy=terval_xy;
    rates.ter_yx=terval_yx;
    rates.optimal_partition_xy=optimal_partition_xy;
    rates.optimal_partition_yx=optimal_partition_yx;
    
    %Build full results output    
    fullresults.mirval=mirval;
    fullresults.camirval_xy=camirval_xy;
    fullresults.camirval_yx=camirval_yx;
    fullresults.optimal_partition_xy=optimal_partition_xy;
    fullresults.optimal_partition_yx=optimal_partition_yx;
    fullresults.camir_xy=camir_xy;
    fullresults.camir_yx=camir_yx;
    fullresults.cami_xy=cami_xy;
    fullresults.cami_yx=cami_yx;
    fullresults.mi=mi;
    fullresults.diridx=diridx;
    fullresults.te_xy=te_xy;
    fullresults.te_yx=te_yx;
    fullresults.terval_xy=terval_xy;
    fullresults.terval_yx=terval_yx;
    
end
