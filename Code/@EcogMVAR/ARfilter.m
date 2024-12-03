function [x_hat, total_er, REV] = ARfilter(x, A, Q)
[N_elec, N_trial, Nt] = size(x);
p = size(A,3);
AA = reshape(A,N_elec,N_elec*p);
x_hat = nan*ones(N_elec, N_trial, Nt-p);
for i = 1:Nt-p
	xx = reshape(permute(x(:,:,i:i+p-1),[1,3,2]), N_elec*p,N_trial);
	er = transpose(mvnrnd(zeros(N_elec,1),Q,N_trial));
	x_hat(:,:,i) = AA*xx + er;
end
total_er = x(:,:,p+1:end) - x_hat;
total_er = sum(total_er.^2,'all')/numel(total_er);
var_x = sum((x(:,:,p+1:end)-mean(x(:,:,p+1:end),'all')).^2,'all')/numel(total_er);
REV = total_er/var_x;
x_hat = cat(3, x(:,:,1:p), x_hat);
end
