function [ flag,AttackAlpha,AttackBeta] = IsReasonble( chromosome,model,uav )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%检查航路是否合理
%位置坐标越界则该航路不满足要求，重新生成
sum_alpha =zeros(1,model.UAV);
sum_beta =zeros(1,model.UAV);
flag=0;
   for i=1:model.dim
      if  chromosome.pos(i,1,uav) <model.Xmin || chromosome.pos(i,1,uav)<0 || ...
          chromosome.pos(i,2,uav)<model.Ymin || chromosome.pos(i,2,uav)>model.Ymax||...
          chromosome.pos(i,3,uav)<model.Zmin || chromosome.pos(i,3,uav) > model.Zmax
          AttackAlpha=0;
          AttackBeta =0;
          flag =0;
          return
      end
   end

   
  %检查最后的偏角能否符合要求
     %航路最后一个点
   lastpoint=[chromosome.pos(model.dim,1,uav),chromosome.pos(model.dim,2,uav),chromosome.pos(model.dim,3,uav)];
   endpoint =[model.ex,model.ey,model.ez];
   last2end = endpoint -lastpoint;
   
  %计算最终偏角方向和最后一个点到终点的方向的夹角
   %分别计算航偏角和俯仰角 
   for i=1:model.dim
      sum_alpha(uav) =sum_alpha(uav) + chromosome.alpha(i,uav);
      sum_beta(uav)  = sum_beta(uav) + chromosome.beta(i,uav);
   end
   
    %计算起始到目标的向量
   st = [model.ex-model.sx(uav),model.ey-model.sy(uav),model.ez-model.sz(uav)];
   %水平向量
   vhorizontal=[1,0];
   %计算起始到目标的航偏角
    st_alpha = rad2deg( acos(dot(st(1:2),vhorizontal)/norm(st(1:2))/norm(vhorizontal) )  );
    %如果正弦值小于0
    if st(2)/norm(st(1:2)) <0
        st_alpha =360 - st_alpha;        
    end
    %计算起始到目标的俯仰角
    st_beta = rad2deg(asin(st(3)/norm(st)));
   %角度转换弧度
    sum_alpha(uav) = sum_alpha(uav) + st_alpha;
    sum_beta(uav) = sum_beta(uav) +st_beta;
    sum_alpha(uav) = deg2rad(sum_alpha(uav));
    sum_beta(uav) = deg2rad(sum_beta(uav));
    %总的航偏角的方向向量
    lastdeg =[cos(sum_alpha(uav)),sin(sum_alpha(uav))];
    %投影到XOY计算航偏角的最后变化值
    theta = rad2deg(acos(dot(last2end(1:2),lastdeg)/norm(last2end(1:2))/norm(lastdeg)));
    %计算last2end的俯仰角
    ag1 = rad2deg(asin(last2end(3)/norm(last2end)));
    %用last2end的俯仰角 - 总的俯仰角 = 从最后一个点到终点的俯仰角变化 
    ag2 =abs( ag1 - sum_beta(uav));
    %计算最后的攻击角
    AttackAlpha(uav) = rad2deg(acos(last2end(1)/norm(last2end(1:2))));
    AttackBeta(uav) = ag1;
    %根据指定攻击角计算每个航偏角平均增加的角度值
%     average_value(uav) = (model.attack_alpha(uav) -  AttackAlpha(uav))/(model.dim+1);
  
    if theta >0 && theta < model.alpha_max &&...
       ag2 >0 && ag2 <model.beta_max
        flag=1;
    else
        flag= 0;
    end

%    %检查能否达到时间上的协同 
%   [flag_time ,ETA_r] =EstimateTime( chromosome,model ); 
%   %两者都满足说明该解符合要求
%   flag_r = flag_time ;
%   ETA =ETA_r;
%   index =flag;
end

