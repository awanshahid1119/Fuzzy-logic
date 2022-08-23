for i = 1:size(train_data,1)
    for j = 1:size(train_data,2)
        [tt2] = find(x_new==(round(train_data(i,j)*1000))); % input at x-axis
        disp(tt2);
        fprintf('Rule done for i = %d and j = %d.\n',i,j);
        for k = 1:size(tt2,2)
            tt3 = 
        end
        [Degree(i,j), Rule(i,j)] = max(x_R(tt2,:));
    end
end