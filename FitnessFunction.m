function [ cost,sol ] = FitnessFunction( chromosome,model )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    for uav=1:model.UAV
    x= zeros(1,model.dim);
    y= zeros(1,model.dim);
    z = zeros(1,model.dim);
    %取第uav个航路的坐标
    for i=1:model.dim
    x(i) = chromosome.pos(i,1,uav);
    y(i) = chromosome.pos(i,2,uav);
    z(i) = chromosome.pos(i,3,uav);
    end
    sx = model.sx(uav);
    sy = model.sy(uav);
    sz = model.sz(uav);
    ex = model.endp(1);
    ey =model.endp(2);
    ez=model.endp(3);
        
    
    xobs = model.xobs;
    yobs = model.yobs;
    zobs = model.zobs;
    robs = model.robs;
    
    XS=[sx x ex];
    YS=[sy y ey];
    ZS=[sz z ez];
    k =numel(XS);
    TP =linspace(0,1,k);
    tt =linspace(0,1,50);
    xx =[];
    yy =[];
    zz=[];
    for i=1:k-1
    %每一段向量分成10个点
    x_r = linspace(XS(i),XS(i+1),10);
    y_r= linspace(YS(i),YS(i+1),10);
    z_r =linspace(ZS(i),ZS(i+1),10);
    xx = [xx,x_r];
    yy = [yy,y_r];
    zz =[zz ,z_r];
    end
    
    %calc L
    dx =diff(xx);
    dy =diff(yy);
    dz = diff(zz);
    Length = sum(sqrt(dx.^2+dy.^2+dz.^2));
    nobs = numel(xobs);
     violation=0;
    for i=1:nobs
       d = sqrt( (xx-xobs(i)).^2+(yy-yobs(i)).^2 );
       v = max(1-d/robs(i),0);
       violation = violation + mean(v);
    end
    sol(uav).TP=TP;
    sol(uav).XS =XS;
    sol(uav).YS=YS;
    sol(uav).ZS=ZS;
    sol(uav).tt=tt;
    sol(uav).xx=xx;
    sol(uav).yy=yy;
    sol(uav).zz=zz;
    sol(uav).dx=dx;
    sol(uav).dy=dy;
    sol(uav).dz=dz;
    sol(uav).Length=Length;
    sol(uav).violation=violation;
    sol(uav).IsFeasible=(violation==0);
    
    %计算协调适应值
    % 3、飞行高度限制
     high=0;
     for k=1:numel(XS)
        x=XS(k);
        y=YS(k);
        h=terrain(x,y);        
        if ZS(k)<=(h+10)  %限制飞行最低高度
            high=high+10000;          
        elseif ZS(k)>375   %限制飞行最高高度              
            high=high+10000;           
        else  
            high=high+abs(ZS(k)-287); %计算与理想高度差距和      
        end        
    end
    
    %z
 
     w1 =0.03;
     w2=0.3;
     w3=0.1;
     w4=0.1;
     %markov evaluatea
     %获取所有维度的坐标
     r_xx=[];r_yy=[];r_zz=[];
    for i=2:numel(XS)-2
    %每一段向量分成10个点
    r_x = linspace(XS(i),XS(i+1),3);
    r_y= linspace(YS(i),YS(i+1),3);
    r_z =linspace(ZS(i),ZS(i+1),3);
    r_xx = [r_xx,r_x];
    r_yy = [r_yy,r_y];
    r_zz =[r_zz ,r_z];
    end
     
   Allpos = [r_xx',r_yy',r_zz'];
   [stateProbabilityProcess, expectedCostProcess]=MarkovEvaluate(Allpos,model);
   sol(uav).MarkovState = stateProbabilityProcess;
   sol(uav).MarkovCost = expectedCostProcess;
   sol(uav).costs=[w1*Length,w3*high,w4*mean(expectedCostProcess)*150,w2*abs(Length/model.vel -chromosome.ETA)*50];
   single_cost(uav)= w1*Length +w3*high+w4*mean(expectedCostProcess)*150+w2*abs( Length/model.vel -chromosome.ETA)*50;
   end
   cost = mean(single_cost);
end

