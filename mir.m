function [mirval]=mir(x,y,xpartition,ypartition,tau)
%MIR Calculates the Mutual Information Rate
%---------------------------------------------
%Input: x,y,: time-series
%       xpartition, ypartition: partition locations in x and y
%       tau: embedding delay (tau=1 for maps)
%----------------------------------------------
%Output:
%       mirval: value of MIR in bits/iteration
%----------------------------------------------
%Example:
%       mirval=mir(x,y,0.5,0.5,1)
%       For a binary initial partition of a map with time-series given by
%       x and y
%----------------------------------------------
%(C) Arthur Valencio* and Dr Murilo S. Baptista, 18 Mar 2018
%    ICSMB University of Aberdeen
%    * Support: CNPq, Brazil
%----------------------------------------------
%   If useful, please cite:
%   A. Valencio and M.S. Baptista (2018). Causality Toolbox: functions for 
%       calculating information theory measures from time-series. Available
%       at: https://github.com/artvalencio/causality-toolbox


    for L=1:10
        mival(L)=mi(x,y,xpartition,ypartition,L,tau);
    end
    figure
    plot(1:7,mival)
    xlabel('L')
    ylabel('Mutual Information')
    title('Find the linear part:')
    list2={'L=1','L=2','L=3','L=4','L=5','L=6','L=7','L=8','L=9','L=10'};
    lstart = listdlg('PromptString',{'Select starting L','(linear part begins):'},...
        'SelectionMode','single','ListString',list2);
    lend = listdlg('PromptString',{'Select finishing L', '(linear part ends):'},...
        'SelectionMode','single','ListString',list2);
    
    coefs=polyfit(lstart:1:lend,mival(lstart:lend),1);
    mirval=coefs(1);
    msgbox(strcat('MIR=',num2str(mirval)));
    
end

function [mival]=mi(x,y,xpartition,ypartition,L,tau)
%calculates the mutual information

    if length(x)~=length(y)
        error('x and y must be of same size!')
    end
    if length(xpartition)~=length(ypartition)
        error('initial partitions in x and in y must be of same size!')
    end
    ns=length(xpartition)+1; %partition resolution
    
%calculating symbols
    for n=1:length(x) %assign data points to partition symbols in x
        Sx(n)=-1;
        for i=1:length(xpartition)
            if x(n)<xpartition(i)
                Sx(n)=i-1;
                break;
            end
        end
        if Sx(n)==-1
            Sx(n)=ns-1;
        end
    end
    for n=1:length(y) %assign data points to partition symbols in y
        Sy(n)=-1;
        for i=1:length(ypartition)
            if y(n)<ypartition(i)
                Sy(n)=i-1;
                break;
            end
        end
        if Sy(n)==-1
            Sy(n)=ns-1;
        end
    end
    
    [p_x,p_y,p_xy]=getprobabilities(Sx,Sy,L,ns,tau,length(x));
    
    %Calculating mutual information
    mival=0;
    for i=1:ns^L
        for j=1:ns^L
            if (p_x(i)*p_y(j)>1e-14)&&(p_xy(i,j)>1e-14)
                mival=mival+p_xy(i,j)*(log(p_xy(i,j)/(p_x(i)*p_y(j))))/log(2);
            end
        end
    end

end

function [p_xp,p_yp,p_xyp]=getprobabilities(Sx,Sy,L,ns,tau,len)
% calculates the symbolic sequences and probabilities

    %initializing phi: removing points out-of-reach (start)
    phi_x(1:tau*L)=NaN;
    phi_yp(1:tau*L)=NaN;
    
    %initializing probabilities of boxes
    p_xp(1:ns^L+1)=0;
    p_yp(1:ns^L+1)=0;
    p_xyp(1:ns^L+1,1:ns^L+1)=0;
    %calculating phi_x, about the past of x
    for n=tau*L+1:len
        phi_x(n)=0;
        k=n-L;%running index for sum over tau-spaced elements
        for i=n-tau*L:tau:n-tau
            phi_x(n)=phi_x(n)+Sx(k)*ns^((n-1)-k);
            k=k+1;
        end
        p_xp(phi_x(n)+1)=p_xp(phi_x(n)+1)+1;
    end
    p_xp=p_xp/sum(p_xp);
    %calculating phi_yp, about the past of y
    for n=tau*L+1:len
        phi_yp(n)=0;
        k=n-L;
        for i=n-tau*L:tau:n-tau
            phi_yp(n)=phi_yp(n)+Sy(k)*ns^((n-1)-k);
            k=k+1;
        end
        p_yp(phi_yp(n)+1)=p_yp(phi_yp(n)+1)+1;
    end
    p_yp=p_yp/sum(p_yp);
    %calculating joint probabilities
    for n=tau*L+1:len
        p_xyp(phi_x(n)+1,phi_yp(n)+1)=p_xyp(phi_x(n)+1,phi_yp(n)+1)+1;
    end
    p_xyp=p_xyp/sum(sum(p_xyp));
end

