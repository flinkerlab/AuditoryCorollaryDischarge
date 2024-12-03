function [xr] = WTDN(x, wname, thresh)
	level = floor(log2(length(x)));
	cA = cell(level+1,1);
	cD = cell(level+1,1);
	cA{1} = x;
	cD{1} = zeros(size(x));
	for i=1:level
		[cA{i+1,1},cD{i+1,1}] = dwt(cA{i},wname);
		cD{i+1,1} = wthresh(cD{i+1,1},'h',thresh);
	end
	cAr = cell(level+1,1);
	cAr{level+1} = cA{level+1};
	for i = level+1:-1:2
		cAr{i-1,1} = idwt(cAr{i,1}(1:length(cA{i,1})),cD{i,1},wname);
	end
	xr = cAr{1,1};
end
