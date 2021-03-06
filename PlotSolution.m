function  PlotSolution(sol,model )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    for uav=1:model.UAV

    XS =sol(uav).XS;
    YS =sol(uav).YS;
    ZS =sol(uav).ZS;
    xx=sol(uav).xx;
    yy=sol(uav).yy;
    zz=sol(uav).zz;
    global Scene;
    %alpha =0:pi/50:2*pi;
    %[x1 y1 z1]=sphere;       %将球体数据写入三矩阵中
    %a=[8 -2 2  4];                  %设置球体参数
    %s1=surf(x1*a(1,4)+a(1,1),y1*a(1,4)+a(1,2),z1*a(1,4)+a(1,3),'FaceColor',[0,0,1]);
    figure(Scene);
    view(0,90);
    hold on;
    if(model.std_ga)
    plot3(xx,yy,zz,'k','LineWidth',2);
    else
     plot3(xx,yy,zz,'r','LineWidth',2);
    end
    plot3(XS(2:model.dim+1),YS(2:model.dim+1),ZS(2:model.dim+1),'ro');
    end
    hold off;
    title('GA')
    grid on;
end

