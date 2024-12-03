function [y,Hout]=pdc(A,Q,f,dt, mode)
% computes partial directed coherence
if ~exist('mode', 'var'); mode = 1; end

F = EcogMVAR.Af(A,f,dt);
if mode==2
	F=F'*F;
end
if mode==3
	Hf = inv(F);
	F = Hf*Hf';
end
if mode==4
	Hf = inv(F);
	F = F'*(Q\F);
end
%F=inv(Af(A,f,dt)); %DTF
[n,~]=size(F);
y=zeros(n);
for a=1:n
    for b=1:n
		if mode==1
			y(a,b)=abs(F(a,b))/sqrt(F(:,b)'*F(:,b));
		elseif mode==2
			y(a,b) = abs(F(a,b))/sqrt(abs(F(a,a).*F(b,b)));
		elseif mode==3
			y(a,b) = abs(F(a,b))/sqrt(abs(F(a,a).*F(b,b)));
		elseif mode==4
			y(a,b) = abs(F(a,b))/sqrt(abs(F(a,a).*F(b,b)));
		elseif mode==5
			y(a,b)=abs((1/Q(a,a)).*F(a,b))/...
				   sqrt((diag(Q).*F(:,b))'*(diag(Q).*F(:,b)));
		else
			error('mode not implemented');
		end
        %y(a,b)=abs(F(a,b))^2/(F(a,:)*F(a,:)'); % DTF
    end
end
if (mode==4)&&(nargout>1)
	Hout = Hf;
end
end
