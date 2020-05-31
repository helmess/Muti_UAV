function [ alpha,T,beta ] = InitialChromosome( model,num,uav)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


   %计算起始到目标的向量 
   st = [model.ex-model.sx(uav),model.ey-model.sy(uav),model.ez-model.sz(uav)];
   dist = norm(st);
   %计算从起始到目标的平均时间
   deltaT = dist/model.vel/(model.dim+1);
   %计算起始到目标的航偏角
    if model.sx - model.ex~=0
    st_alpha = rad2deg(atan((model.ey-model.sy(uav))/(model.ex-model.sx(uav))));
    else
    st_alpha = rad2deg(asin(1));
    end
    %计算起始到目标的俯仰角
    st_beta = rad2deg(asin(st(3)/norm(st)));
%初始化染色体基因    
alpha = zeros(model.dim,1);
beta = zeros(model.dim,1);
T=zeros(model.dim,1);
%%前num个染色体初始化
if num < model.num
sum_alpha =st_alpha;
sum_beta =st_beta;
for i =1:model.dim 
   if sum_alpha < st_alpha
       %alpha(i)在(0,alpha_max)
      alpha(i) = rand*(model.alpha_max - (model.alpha_max+model.alpha_min)/2)+(model.alpha_max+model.alpha_min)/2;
   else
       %alpha(i)在(alpha_min,0)
      alpha(i) = rand*( (model.alpha_max+model.alpha_min)/2-model.alpha_min )+model.alpha_min;
   end
   %同上
   if sum_beta < st_beta
      beta(i) = rand*(model.beta_max - (model.beta_max+model.beta_min)/2)+(model.beta_max+model.beta_min)/2;
   else
      beta(i) = rand*((model.beta_max+model.beta_min)/2 -model.beta_min )+model.beta_min;
   end
      sum_beta = sum_beta + beta(i);
      sum_alpha = sum_alpha+ alpha(i);
      T(i)= (rand*(2*deltaT)-deltaT)*0.1 + deltaT;
end
%%num后染色体随机初始化
else
   for i=1:model.dim
        alpha(i)=  (rand*2 -1)*model.alpha_max;
        beta(i)=  (rand*2 -1)*model.beta_max;
        T(i)=(rand*(2*deltaT)-deltaT)*0.1+ deltaT;
    end
 
end


end

