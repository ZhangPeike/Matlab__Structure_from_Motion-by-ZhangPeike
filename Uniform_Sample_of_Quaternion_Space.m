%A New Approach for Solving the 5-Point Relative Pose Problem
%Part A producing Q sample
%%
clear;
close all;
disp('Start!');
time1=datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(time1);
frames.images{1}='templeR0013.png';
frames.images{2}='templeR0014.png';
N=5000;
%
f=zeros(4,N);
r=zeros(4*N,N);
%step1

Q_Set=(rand(4,N)-0.5)*2;

%step2
for i=1:N
    Q_Set(:,i)=Q_Set(:,i)/norm(Q_Set(:,i));
end
disp('Quaternion is created initially...');
disp(Q_Set);
%%
%step3
Q_Set_ori=Q_Set;
for i=1:N
   for j=1:N
       if j~=i
           r(4*j-3:4*j,i)=Q_Set(:,i)-Q_Set(:,j);
       end
   end
end
for i=1:N
    count=0;    
    for j=1:N        
        count=count+1;        
        if j~=i
            if count<2
            f(:,i)=r(4*j-3:4*j,i)/(norm(r(4*j-3:4*j,i)))^3;
            else
            f(:,i)=f(:,i)+r(4*j-3:4*j,i)/(norm(r(4*j-3:4*j,i)))^3;
            end
        end
    end
end
%step4
Q_Set=Q_Set+f;
disp('Move q...');
disp('Without normalizing');
disp(Q_Set);
for i=1:N
    Q_Set(:,i)=Q_Set(:,i)/norm(Q_Set(:,i));
end
disp('Normalizing');
Q_Set_update=Q_Set;   
disp(Q_Set_update);

Bmax=zeros(1,N);
for i=1:N
   Bmax(i)=norm(Q_Set_ori(:,i)-Q_Set_update(:,i));
end
Max=max(Bmax);
disp('residual');
disp(Bmax);
disp(Max);
%%
while Max>0.005,
    disp('Iterating...');
    Q_Set_ori=Q_Set_update;
    Q_Set=Q_Set_update;
    for i=1:N
       for j=1:N
           if j~=i
               r(4*j-3:4*j,i)=Q_Set(:,i)-Q_Set(:,j);
           end
       end
    end
    for i=1:N
        count=0;    
        for j=1:N        
            count=count+1;        
            if j~=i
                if count<2
                f(:,i)=r(4*j-3:4*j,i)/(norm(r(4*j-3:4*j,i)))^3;
                else
                f(:,i)=f(:,i)+r(4*j-3:4*j,i)/(norm(r(4*j-3:4*j,i)))^3;
                end
            end
        end
    end
    %step4
    Q_Set=Q_Set+f;
    for i=1:N
        Q_Set(:,i)=Q_Set(:,i)/norm(Q_Set(:,i));
    end
    Q_Set_update=Q_Set;   
    for i=1:N
       Bmax(i)=norm(Q_Set_ori(:,i)-Q_Set_update(:,i));
    end
    Max=max(Bmax);
    disp('Iterating residual:');
    disp(Max);
end
%%
%Show the sampling result
%{
for i=1:N
    Rt(:,:,i)=[Qua2R(Q_Set_update(:,i))';zeros(1,3)]';
    drawCamera(Rt(:,:,i),0.01,0.01,1,0.2,0.1);
end
%}
for i=1:N
    ang(i)=2*acos(Q_Set_update(1,i));
    Vec(:,i)=Q_Set_update(2:4,i)/sin(0.5*ang(i));
    plot3(Vec(1,i),Vec(2,i),Vec(3,i),'g+');
    hold on
end
time2=datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(time2);
disp('Program OK');
