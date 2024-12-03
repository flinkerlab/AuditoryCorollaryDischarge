function Visual_view2d(x, w)
        if nargin<2
			N = size(x,2)*0.05;
            w = (1/N)*ones(N,1);
        end
        [z,m] = size(x);
        slice = 1;

        fig1 = figure;
		if size(w,1)>1
			w = w';
		end
		x_smooth = conv2(x,w);
        mini = min(min(x(:)), min(x_smooth(:)));
        maxi = max(max(x(:)), max(x_smooth(:)));
		sig = squeeze(x(slice,:));
		sig_mean = squeeze(mean(sig,'all'));
		sig_smooth = x_smooth(slice,:);

        plt= plot(sig, 'b-','LineWidth', 1);
		hold on;
        plt_smooth = plot(sig_smooth, 'r-','LineWidth', 1.5);
        plt_mean = plot([1:length(sig)]*sig_mean, 'g--','LineWidth', 1);
        ylim([mini-0.1,maxi+0.1]);
		grid on;
        drawnow;

        txt = uicontrol(fig1,...
                'Style','text',...
                'String','slice 1',...
                'Units', 'Normalized',...
                'Position', [0.2 0.05 0.6 0.05]);
        sld = uicontrol(fig1,...
                'Style', 'slider',...
                'Min',1,'Max',z,'Value',1,...
                'Units', 'Normalized',...
                'Position', [0.2 0.0 0.6 0.05],...
                'Callback', {@fun1, x, x_smooth, txt, plt, plt_smooth, plt_mean },...
                'SliderStep', [1/(z-1) 1]);

end

function fun1(hObject, eventdata, x, x_smooth, txt, plt, plt_smooth, plt_mean)
slice = floor(get(hObject,'Value'));
plt.YData = x(slice,:);
plt_mean.YData = squeeze(mean(x(slice,:),'all'))*ones(1,size(x,2));
plt_smooth.YData = x_smooth(slice,:);
txt.String = ['slice ',num2str(slice)];

end

