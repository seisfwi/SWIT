% pick_line : pick a line from a figure;
%
%
% Based on doc(ginput);

function [xy,xys]=pick_line(name,LineWidth,LineColor);


if nargin<1, name='line';end
if nargin<2, LineWidth=3;end
if nargin<3, LineColor=[1 1 0];end

hold on
% Initially, the list of points is empty.
xy = [];
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);
    n = n+1;
    xy(:,n) = [xi;yi];
    if n ==1;
        p=plot(xy(1),xy(2),'-','LineWidth',LineWidth,'Color',LineColor);
        p_point=plot(xy(1),xy(2),'*','LineWidth',LineWidth,'Color',LineColor);
    else
        set(p_point,'Xdata',xy(1,:));
        set(p_point,'Ydata',xy(2,:));

        try
            %SPLINE INTERPOLATION
            t = 1:n;
            ts = 1: 0.1: n;
            xys = spline(t,xy,ts);
            set(p,'Xdata',xys(1,:));
            set(p,'Ydata',xys(2,:));
            %        p=plot(xy(1,:),xy(2,:),'w-');
        end
    end
    drawnow;
    %plot(xi,yi,'ro')
end

hold off

eval(sprintf('save %s.mat xy xys',name));



