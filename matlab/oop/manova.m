classdef manova < explanatory
    %MANOVA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimensions;
        pval;
        statistics;
        
        
    end
    
    methods
        function this = manova(group,options)
            %MANOVA Construct an instance of this class
            %   Detailed explanation goes here
%             addpath lib;
            this = this@explanatory(group,options);
        end
        function output=HottelingT(this)
            addpath HotellingT2;
            
            inp2 = this.diagnosis;
            male = find(inp2>0);
            female = find(inp2==0);
            n1= size(male,1);
            n2=size(female,1);
            this.pval = zeros(this.num_edge,length(unique(inp2))-1);
            
            for j_edge = 1 : this.num_edge
                a=squeeze(this.all_edges(j_edge, :,:));
                inp1=a;
                X = [inp1(male,:);inp1(female,:)];
                [p] = T2Hot2ihe(X,n1,n2,0.05);
                this.pval(j_edge,:) = p(1);
            end
            [output]=this.myNBSfdr(this.pval(:,1),0.05,100);
        end
        function output=run(this)
            inp2 = this.diagnosis;
            this.pval = zeros(this.num_edge,length(unique(inp2))-1);
            for j_edge = 1 : this.num_edge
                inp1 = squeeze(this.all_edges(j_edge, :,6))';
                [d,p,~] = manova1(inp1,inp2);
                this.dimensions(j_edge) = d;
                this.pval(j_edge,:) = p(1);
            end
            [pID,pN] = this.my_fdr(this.pval(:,1),0.01);
        end
        function [mask]=myNBSfdr(this,network,eta,iterations)
            [matrix,biggestArea,lcc] = this.permTest(this.pval(:,1),eta,iterations);
            [~,L,bigval] = this.my_largestCC(network<eta);
            stats = regionprops(L,'area');
            avec=[stats.Area];
            [in val] = find(avec>biggestArea);
            mask = zeros(this.num_node,this.num_node);
            for each=val
                mask =mask+L==each;
            end
            
        end
        
        function [matrix,biggestArea,lcc]=permTest(this,network,eta,iterations)
            matrix = zeros(iterations+1,this.num_edge);
            matrix(1,:)=network;
            lcc=zeros(iterations,1);
            for i=1:iterations
                disp(i);
                inp2 = this.diagnosis;
                inp2 = inp2(randperm(length(inp2)));
                this.pval = zeros(this.num_edge,length(unique(inp2))-1);
                for j_edge = 1 : this.num_edge
                    inp1 = squeeze(this.all_edges(j_edge, :,:));
                    [d,p,~] = manova1(inp1,inp2);
                    this.dimensions(j_edge) = d;
                    this.pval(j_edge,:) = p(1);
                end
                matrix(i,:)=this.pval(:,1);
                a = this.pval(:,1)<eta;
                [~,~,bigval] = this.my_largestCC(a);
                lcc(i)=bigval;
            end
            [f,x] = ecdf(lcc);
            pos = find((1-eta)<=f,1);
            biggestArea = x(pos);
        end
        function clusterPlot(this,mask,X,Y)
            [cx,cy] = find(mask);
            lindices = tril(ones(268));
            lindices(find(eye(268)))=0;
            lindices = find(lindices);
            P = zeros(length(unique(Y)),length(cx),length(X(1,1,:))); % 4*25*6
            for g=1:length(unique(Y))
                indices = Y==g;
                Xg=squeeze(mean(X(:,indices,:),2));
                
                for i=1:length(cx)
                    for j=1:length( Xg(i,:))
                        network = zeros(268,268);
                        network(lindices) = Xg(:,j);
                        network = network'+network;
                        P(g,i,j)=network(cx(i),cy(i));
                    end
                end
            end
            %             subplot(1,3,3);
            d1=length(unique(Y));
            d2=5;%length(cx);
            for g=1:d1
                for j=1:d2
                    disp((g-1)*d2+j);
                    subplot(d1,d2,(g-1)*d2+j); %length(unique(Y))
                    bar(squeeze(P(g,j,:)))
                end
            end
            disp(P);
        end
    end
end

