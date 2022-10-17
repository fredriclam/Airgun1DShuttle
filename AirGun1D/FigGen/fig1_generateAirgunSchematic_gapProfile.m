% Generate `metadata` from performing one run
L =  metadata.discretization.chambers.total_travel_length;
x = linspace(0,L,1000);

dr = nan(size(x));
for i = 1:length(x)
    dr(i) = metadata.discretization.chambers.gapProfile(x(i));
end

figure(1); clf;
set(gcf,'position',[2079, 136, 601, 280]);