function [dummy] = draw_cs(coord_x,irr,af,alg_num,nx,ny)
%
dummy = 1;
%
% figure; %
figure(alg_num) %
clf reset;
hold on;
axis equal
axis off;
%
nm = size(irr,1);
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    mem0(1,:,i) = coord_x(j1,:);
    mem0(2,:,i) = coord_x(j2,:);
end
%
% af = af / (max(af) / 10); % original
% af = af / (max(af) / 4); %

for i=1:nm
    if af(i) > 0%1.5e-8 % Threshold value (1.0*10^(-4) originally) 1.5e-8
        plot(mem0(:,1,i), mem0(:,2,i), 'k-', 'LineWidth',af(i)/(max(af)/15));
%     elseif af(i) > 0 %
%         plot(mem0(:,1,i), mem0(:,2,i), 'r-', 'LineWidth',af(i)/(max(af)/15)); %
    end
    hold on;
end

for j=1:size(coord_x,1)
    plot(coord_x(j,1), coord_x(j,2), 'ok',...
        'MarkerFaceColor', 'w', 'MarkerSize',6);
end

plot(coord_x(1,1), coord_x(1,2), 'ok',...
    'MarkerFaceColor', 'k', 'MarkerSize',16);
% h(1) = plot(coord_x(nx+1,1), coord_x(nx+1,2), 'ok',...
%     'MarkerFaceColor', 'k', 'MarkerSize',16);
h(1) = plot(coord_x(nx+1,1), coord_x(1,2), 'ok',...
    'MarkerFaceColor', 'k', 'MarkerSize',16);
% h(2) = plot(coord_x((nx+1)*(ny+1),1), coord_x((nx+1)*(ny+1),2), 'ok',...
%     'MarkerFaceColor', 'r', 'MarkerSize',16); % Non-struct top-right
% h(2) = plot(coord_x(nx/2+1,1), coord_x(nx/2+1,2), 'ok',...
%     'MarkerFaceColor', 'r', 'MarkerSize',16); % Non-struct bottom-center
h(2) = plot(coord_x(nx/2+1,1), coord_x(1,2), 'ok',...
    'MarkerFaceColor', 'r', 'MarkerSize',16); % Non-struct mass middle-right

% legend(h([1,2]),{'Fixed','Mass'},'interpreter','latex','Location','northeast')
% set(gca,'FontSize',28)
% ax = gca;
% ax.FontSize = 28;
% ax.TickLabelInterpreter = 'latex';
% ax.YTickLabel = cell(size(ax.YTickLabel));
% ax.XTickLabel = cell(size(ax.XTickLabel));

axis equal