function [totalcorrel,totalentropy]=totals(x,l,ns,linepos,tau)
%TOTALS Calculates the total correlation and total joint entropy of a multivariate system
%---------------------------------------
%Inputs:
%       A: time-series of multiple nodes (each variable as a column)
%       l: length of the symbolic sequence
%       ns: number of initial symbols (symbolic base)
%       linepos: position of partition division lines (each variable as a column)
%       tau: lag to be applied if flow, 1 if map
%---------------------------------------
%Outputs: 
%       correl: Total Correlation
%       entropy: Joint Entropy
%---------------------------------------
%Example:
%       
%       [totalcorrel,totalentropy]=totals(randn(5e6,5),4,3,(ones(5,1)*[-0.4 0.4])',1);
%       
%       This example calculates the total correlation for
%       a system made of 5 gaussian distributions.
%       The user selected a partition with 3 symbols (0,1,2) defined by the
%       splitting lines traced at -0.4 and 0.4. The selected symbolic sequence 
%       length is of 4 symbols. 
%       It means that each 'box' is given by a sequence such as '0-2-1-1', 
%       where the digit '0','1' or '2' represents if the data point is 
%       smaller than -0.4, between -0.4 and 0.4 or greater than 0.4.
%       
%---------------------------------------
%LaTeX Expression of the definitions:
%       Total Correlation:
%       $C(x_1,x_2,...,x_n)=\sum_{x_1,x_2,...,x_n} p(x_1,x_2,...,x_n) \log{\frac{p(x_1,x_2,...,x_n)}{\prod_i p(x_i)}}$
%
%       Joint Entropy:
%       $H(x_1,x_2,...,x_n)=\sum_{x_1,x_2,...,x_n} p(x_1,x_2,...,x_n) \log{p(x_1,x_2,...,x_n)}$
%---------------------------------------
%(C) Arthur Valencio(1)* and Murilo Baptista(1), 11 December 2017
%  incorporating discussions with Nicolas Rubido(2)
%  (1)ICSMB, University of Aberdeen,UK
%  (2)Universidad de la Republica, Uruguay
%   *Support: CNPq, Brazil

    %initial defs
    len=length(x(:,1));
    numelements=length(x(1,:));
    S(1:len,1:numelements)=-1;    
    
    %calculating symbols
    for j=1:numelements
        for n=1:len %assign data points to partition symbols in x
            S(n,j)=-1;
            for i=1:length(linepos(:,j))
                if x(n,j)<linepos(i,j)
                    S(n,j)=i-1;
                    break;
                end
            end
            if S(n,j)==-1
                S(n,j)=ns-1;
            end
        end
    end
    
    %get probs
    [pjoint,pindiv]=getprobabilities2(S,l,ns,tau,len,numelements);
    
    %calculate total correlation
    totalcorrel=0;
    for i=1:length(pjoint)
        %calc product term
        product=0;
        for k=1:numelements
            product=product+pindiv(k);
        end
        %calc correl
        if product>0 && pjoint(i)>0
            totalcorrel=totalcorrel+pjoint(i)*log(pjoint(i)/product);
        end
    end
    
    %calculate total entropy
    totalentropy=0;
    for i=1:length(pjoint)
        if pjoint(i)>0
            totalentropy=totalentropy+pjoint(i)*log(pjoint(i));
        end
    end
    
end

function [pjoint,pindiv]=getprobabilities2(S,l,ns,tau,len,nodes)

    %initializing symbolic box-counter
    phi(1:len,1:nodes)=NaN;
    
    %calculating individual probabilities
    pindiv(1:ns^l+1,1:nodes)=0;
    for j=1:nodes
        for n=tau*l+1:len
            phi(n,j)=0;
            k=n-l;%running index for sum over tau-spaced elements
            for i=n-tau*l:tau:n-tau
                phi(n,j)=phi(n,j)+S(k,j)*ns^((n-1)-k);
                k=k+1;
            end
            pindiv(phi(n,j)+1,j)=pindiv(phi(n,j)+1,j)+1;
        end
        pindiv(:,j)=pindiv(:,j)/sum(pindiv(:,j));
    end
    
    %calculating joint probabilities
    maxidx=max(phi(:))+1+nodes*len;
    pjoint(1:maxidx)=0;
    for n=tau*l+1:len
        for j=1:nodes
            pjoint(phi(:,j)+1)=pjoint(phi(:,j)+1)+1;
        end
    end
    pjoint=pjoint/sum(pjoint);
end