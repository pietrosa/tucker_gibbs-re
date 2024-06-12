function ret = stack_images(X_list, n, m)
ret = cat(2,X_list{1:m});
for i = 2:n
    ret_ = cat(2,X_list{(1:m)+(i-1)*m});
    ret = [ret; ret_];
end
