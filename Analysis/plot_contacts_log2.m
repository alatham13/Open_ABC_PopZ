clear all;
close all;

% Plot FL_salt300 -----------------------------------------------------------

% import data
PopZ1=importdata('inter_map.txt');
PopZ2=importdata('intra_map.txt');

% N is the number of chains
N=max(size(PopZ1));

% contact matrix over all chains
contanct_FL=zeros(N,N);

% loop over all chains
for i=1:N
    % set inter-molecular along diagonal to prevent loss of information
    contanct_FL(i,i)=2*PopZ1(i,i);
    for j=i+1:N
        % top diagonal is intra molecular. Normalize based on the symetery
        % of the chain
        contanct_FL(i,j)=(PopZ2(i,j)+PopZ2(j,i))/2;
        % bottom diagonal is inter molecular. Sum together to account for
        % interactions affecting both chains, which is intentionally not
        % included in the orignal script for speed
        contanct_FL(j,i)=(PopZ1(i,j)+PopZ1(j,i));
    end
end

% plot
figure;
heatmap(contanct_FL,'CellLabelColor','none','GridVisible','off');
colormap(flipud(hot));
set(gca,'FontSize',52,'FontName','Helvetica');
%set(gca,'FontSize',52,'FontName','Helvetica','Colorscaling','log');
ax=gca;
caxis([0 0.01]);