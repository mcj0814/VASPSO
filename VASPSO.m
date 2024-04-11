%%VASPSO
function [gbestValue]=VASPSO(N,T_Max,LB,UB,dim,fun)
% function [cg_curve]=VASPSO(N,T_Max,LB,UB,dim,fun)

Vmax0=0.5*(UB-LB);
w_max=0.9;
w_min=0.4;
x=rand;

%% Learning Factor

c1=2;
c2=c1;

%% Initialize Population

HN = round(N/2);%Divide
Xmax=repmat(UB,N,1);
Xmin=repmat(LB,N,1);
HXmax=repmat(UB,HN,1);
HXmin=repmat(LB,HN,1);
X=Xmin+(Xmax-Xmin).*rand(N,dim);
HX=HXmin+(HXmax-HXmin).*rand(HN,dim);
%Vmax=repmat(Vmax0,N,1);
HVmax=repmat(Vmax0,HN,1);
HV = -HVmax+2*HVmax.*rand(HN,dim);

%% Evaluation of particles.
fX(N)=inf;
for i=1:N
    fX(i)=cec17_func(X(i,:)', fun);
end
FEs=N;
%% Initialize Pbest and Gbest

Pbest=X;   
fPbest=fX;
[gbestValue, gbestIndex]=min(fPbest); 
Gbest=Pbest(gbestIndex,:); 
%Fbest=[Fbest gbestValue];   


%% Iteration
for t=2:T_Max%开始迭代，从第2次迭代开始。
    %% Update Position X and Velocity V and check the responding bounds更新位置X和速度V，并检查边界值。

     w(t)=exp(-(2.5 * t / T_Max)^2.5);
   
    %Sort 
     for i = 1:N
       fX(i) = cec17_func(X(i,:)', fun);
    end    
    fX = sort(fX);
    
    u=randperm(N,2);
  
     for i=1:HN 
        [U,index]=min(fPbest(u));
        if min(fPbest(u))<fPbest(i)
            Ubest(i,:)=Pbest(index,:);
        else
            Ubest(i,:)=Pbest(i,:);
        end
     end
    %Velocity pausing
    if rand< 0.3
    HV=HV.^(w(t)*rand)+c1*rand(HN,dim).*(Ubest-HX)+c2*rand(HN,dim).*(repmat(mean(Pbest),HN,1)-HX);
   HV=max(-HVmax,min(HVmax,HV));
    end
    
    % First swarm
    for i=1:HN
        if exp(fX(i))/exp(mean(fX))>rand
            HX(i,:)=w*HX(i,:)+(1-w)*HV(i,:)+Gbest;
        else
            HX(i,:)=HX(i,:)+HV(i,:);
        end
    end
    HX=max(HXmin,min(HXmax,HX));
    
    %% Second swarm
      for i = HN:N
      for j = 1:dim
           BackorFrong = w(t) * rand * abs(Gbest(j))^w(t); 
           if rand < 0.49
              X(i, j) = Gbest(j) + BackorFrong;
           else
               X(i, j) = Gbest(j) - BackorFrong;
           end
       end
    
    end
    
    %% Terminal replacement mechanism
    for i=1:N
        fX(i)=cec17_func(X(i,:)', fun);
    end
    [worst,index]=max(fX);
    z=[1:index-1,index+1:N];
    d=randperm(length(z),2);
    NewX=Gbest+rand*(Pbest(z(d(2)),:)-Pbest(z(d(1)),:));
    NewX=max(LB,min(UB,NewX));
    fNewX=cec17_func(NewX',fun);
    if rand< 0.99   
    if fNewX<worst
        X(index,:)=NewX;
        fX(index)=fNewX;
    end
    end
        
    %% Evaluate Fitness and Change Pbest
    for i=1:N
        if fX(i) < fPbest(i)
            Pbest(i,:)=X(i,:);
            fPbest(i)=fX(i);
        end
        if fPbest(i)<gbestValue
            Gbest=Pbest(i,:);
            gbestValue=fPbest(i);
        end
    end
    FEs = FEs + N;
    %% Plot
    
       cg_curve(t) = gbestValue;

end
end  