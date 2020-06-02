function [ globel ] = Muti_Uav_Ga(model )
%定义染色体
%123
model.endp =[model.ex,model.ey,model.ez];
my_chromosome.pos=[];
my_chromosome.alpha=[];
my_chromosome.beta=[];
my_chromosome.atkalpha=[];
my_chromosome.atkbeta=[];
my_chromosome.T=[];
my_chromosome.sol=[];
my_chromosome.cost=[];
my_chromosome.costs=[];
my_chromosome.ETA=[];
my_chromosome.IsFeasible=[];
my_chromosome.initialized_uav=[];
my_chromosome.Paths=[];
%初始染色体个数
chromosome = repmat(my_chromosome,model.NP,1);
%子代染色体
next_chromosome = repmat(my_chromosome,model.NP,1);
%两代染色体
AllChromosome = repmat(my_chromosome,model.NP*2,1);
%种群的适应度值
seeds_fitness=zeros(1,model.NP);
%局部最优
local =repmat(my_chromosome,model.UAV,1);
for uav=1:model.UAV
local(uav).cost =inf;
end
globel.cost =inf;
%种群初始化
h= waitbar(0,'initial chromosome');
for i=1:model.NP
  chromosome(i).initialized_uav=zeros(model.UAV,1);
  while sum(chromosome(i).initialized_uav) ~=model.UAV
  %需要初始化的无人机序号
  index = find(chromosome(i).initialized_uav==0);
  %初始化角度和时间
  for j =1:numel(index)
  uav =index(j);
  [chromosome(i).alpha(:,uav),chromosome(i).T(:,uav),chromosome(i).beta(:,uav)] = InitialChromosome(model,i,uav);
  %根据角度获得对应坐标
  [chromosome(i).pos(:,:,uav)] = Angel2Pos(chromosome(i),model,uav);
  %形成可执行路径后,由于实际的路径可能比起始到目标的直线距离远,调整运行时间T
   [chromosome(i).T(:,uav),chromosome(i).Paths(uav)] =Modify_Chromosom_T(chromosome(i),model,uav);
   %重新计算新的pos
  [chromosome(i).pos(:,:,uav)] = Angel2Pos(chromosome(i),model,uav);
  %检查各无人机航路是否合理
  [flag,chromosome(i).atkalpha(uav),chromosome(i).atkbeta(uav)] = IsReasonble(chromosome(i),model,uav);
  chromosome(i).initialized_uav(uav)=flag;
  end
  %计算协同满足要求
  max_length = max(chromosome(i).Paths);
  %以最长距离的航路作为协同时间
  chromosome(i).ETA = max_length / model.vel;
  end
  chromosome(i).IsFeasible =1;
  %计算每个符合协调函数解的适应度值和每个解的具体解决方案
  [chromosome(i).cost,chromosome(i).sol,chromosome(i).costs] = FitnessFunction(chromosome(i),model);
  %记录所有解的适应度值，作为轮盘赌的集合
  seeds_fitness(i) = chromosome(i).cost;
  h=waitbar(i/model.NP,h,[num2str(i),':chromosomes finished']);
  
end
close(h)
% for i=1:model.NP
% PlotSolution(chromosome(i).sol,model);
% end
%开始迭代进化
for it=1:model.MaxIt
    %得到最大和平均适应度值
    model.f_max =max(seeds_fitness);
    model.f_avg =mean(seeds_fitness);
    %由于适应度值越小越好
    seeds_fitness = 1./seeds_fitness;
    total_fitness = sum(seeds_fitness);
    seeds_probability = seeds_fitness/ total_fitness;
    %计算累计概率
    seeds_accumulate_probability = cumsum(seeds_probability, 2);
    %根据轮盘赌选择父母,总共选择出NP个子代
    
    for seed=1:2:model.NP
    flag =0;
    %保证父母和子代都符合要求
    while flag~=1
    [parents,flag] = SelectChromosome(seeds_accumulate_probability,model,chromosome);
    %在父母染色体进行基因重组和变异操作，
    %并获得保证每个子代都符合约束条件
    end
    
    papa=randi(model.UAV,1,1);
    if local(papa).cost~=inf && it <5
        parents(2) = local(papa);
    end
    [ sons] = CrossoverAndMutation( parents,model );
    
    %符合要求以后计算子代的适应度值
    [sons(1).cost,sons(1).sol] = FitnessFunction(sons(1),model);
    [sons(2).cost,sons(2).sol] = FitnessFunction(sons(2),model);
    next_chromosome(seed) = (sons(1));
    next_chromosome(seed+1) = (sons(2));
    end
   %把新旧合并同一种群
    AllChromosome(1:model.NP) = chromosome(1:model.NP);
    AllChromosome(model.NP+1:model.NP*2) = next_chromosome(1:model.NP);
    %精英保留,新旧种群一起比较
    
    for i=1:model.NP*2
    eval_array(i,:) = [i,AllChromosome(i).cost];
    end
    %以cost从小到大进行排序
    eval_array =sortrows(eval_array,2);
    last_cost=eval_array(1,2);
    cnt =1;
    chromosome(cnt) = AllChromosome(eval_array(1,1));
    %下次迭代的染色体为不重复cost的最优染色体
    for i=2:model.NP*2
        current_cost = eval_array(i,2);
        if current_cost ~= last_cost
        cnt = cnt+1;
        chromosome(cnt) = AllChromosome(eval_array(i,1));
        last_cost = current_cost;
        end
    end
    %如果下次迭代的染色体数目不够，就根据轮盘赌补染色体。
    cnt_r =cnt;
    while cnt <model.NP
        cnt= cnt+1;
        chromosome(cnt) = AllChromosome(eval_array(cnt - cnt_r,1));
    end
    %选出迭代的染色体和全局最优染色体
    for index =1:model.NP
        seeds_fitness(index) =chromosome(index).cost; 
        for uav=1:model.UAV
            if local(uav).cost >chromosome(index).costs(uav)
                local(uav) = chromosome(index);
                local(uav).cost =chromosome(index).costs(uav);
            end
        end
        %全局最优
        if globel.cost >chromosome(index).cost
            globel =chromosome(index);
        end
    end

    best(it) = globel.cost;

     disp(['it: ',num2str(it),'   best value:',num2str(globel.cost)]);
    
    
    
end
PlotSolution(globel.sol,model )
figure;
plot(best);
end


