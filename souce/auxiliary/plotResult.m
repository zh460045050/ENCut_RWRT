function [  ] = plotResult( X, final_lab, line_list )
%PLOTRESULT 此处显示有关此函数的摘要
%   此处显示详细说明

num_lab = max(max(final_lab));
num_line = length(line_list);
v_lab = cell(1, num_lab);
for i = 1:num_lab
    v_lab{i} = X(final_lab ==i, :);
end
%v_lab3 = v(lab ==3, :);
figure();
hold on;
color = ['r','b','k','c','m','g','y'];
for i = 1:num_lab
    scatter(v_lab{i}(:,1), v_lab{i}(:,2), 'o', color(mod(i, 7) + 1));
end


% p_x = [min(X(:,1)):0.1:max(X(:,1))];
% for i = 1:num_line
%     normV = line_list{i}.V / sum(line_list{i}.V);
%     %normV = line_list{i}.V;
%     b = line_list{i}.b;
%     p_y = (-normV(1) * p_x + b) / normV(2);
%     plot(p_x, p_y, '-');
% end



xlim([min(X(:,1)), max(X(:,1))])
ylim([min(X(:,2)), max(X(:,2))])
hold off;

end

