function [ vel,alpha,beta,T ] = Update_vel_pos( next_chromosome,model,uav )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    w=0.98;
    c1=1;
    c2=1;
    vel=zeros(3,model.dim);
    %更新航偏角速度
    vel(1,:) =  w*next_chromosome.vel(1,:,uav)+...
        c1*rand(1,model.dim).*(next_chromosome.best.alpha(:,uav)' - next_chromosome.alpha(:,uav)')+...
        c2*rand(1,model.dim).*(model.p_global.alpha(:,uav)' - next_chromosome.alpha(:,uav)');
    %更新俯仰角速度
     vel(2,:) =  w*next_chromosome.vel(2,:,uav)+...
        c1*rand(1,model.dim).*(next_chromosome.best.beta(:,uav)' - next_chromosome.beta(:,uav)')+...
        c2*rand(1,model.dim).*(model.p_global.beta(:,uav)' - next_chromosome.beta(:,uav)');
    %跟新时间速度
    vel(3,:) =  w*next_chromosome.vel(3,:,uav)+...
        c1*rand(1,model.dim).*(next_chromosome.best.T(:,uav)' - next_chromosome.T(:,uav)')+...
        c2*rand(1,model.dim).*(model.p_global.T(:,uav)' - next_chromosome.T(:,uav)');
    %速度约束
    vel_alpha_max =0.1*(model.alpha_max-model.alpha_min);
    vel_alpha_min =-vel_alpha_max;
    vel_beta_max =0.1*(model.beta_max-model.beta_min);
    vel_beta_min =-vel_beta_max;
    vel_T_max =0.1*(model.Tmax-model.Tmin);
    vel_T_min =-vel_T_max;
    %约束
    vel(1,:) =max(vel(1,:),vel_alpha_min);
    vel(1,:) =min(vel(1,:),vel_alpha_max);
    
    vel(2,:) =max(vel(2,:),vel_beta_min);
    vel(2,:) =min(vel(2,:),vel_beta_max);

    vel(3,:) =max(vel(3,:),vel_T_min);
    vel(3,:) =min(vel(3,:),vel_T_max);    
    %跟新alpha,beta,T
    alpha = next_chromosome.alpha(:,uav) + vel(1,:)';
    beta = next_chromosome.beta(:,uav) + vel(2,:)';
    T = next_chromosome.T(:,uav) + vel(3,:)';
    alpha =max(alpha,model.alpha_min);
    alpha =min(alpha,model.alpha_max);
    
    beta =max(beta,model.beta_min);
    beta =min(beta,model.beta_max);
    
%     T =max(T,model.Tmin);
%     T =min(T,model.Tmax);
    
end

