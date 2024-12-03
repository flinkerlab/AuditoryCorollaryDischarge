function y=Af(A,f,dt)
[nd,~,np]=size(A);
y=eye(nd);
for tp=1:np
    y=y+A(:,:,tp)*exp(-2*pi*1i*tp*f*dt);
end
end