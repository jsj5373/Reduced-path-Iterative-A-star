close all;
clear all;
clc;


% for i = 1 : q_size+1

%             q_V(:,i) = [-5*sin(pi/25*(i-1))*pi/25,0,0];
%             q_A(:,i) = [-5*cos(pi/25*(i-1))*pi/25*pi/25,0,0];
% end

%         q = zeros(3,q_size+1);
%         q_V = zeros(3,q_size+1);
%         q_A = zeros(3,q_size+1);
%         q(:,1) = [0,0,0];  %%straight trajectory
%         for i = 2 : q_size/5+1
%             q(:,i) = [0,0,0.1*(i-1)*(i-1)];
%             q_V(:,i) = [0,0,0.2*(i-1)];
%             q_A(:,i) = [0,0,0.2];
%         end
%         for i = q_size/5+2 : 2*q_size/5+1
%             q(:,i) = [0.1*(i-11)*(i-11),0,10];
%             q_V(:,i) = [0.2*(i-11),0,0];
%             q_A(:,i) = [0.2,0,0];
%         end
%         for i = 2*q_size/5+2 : q_size+1
%             q(:,i) = [10+2*(i-21),0,10];
%             q_V(:,i) = [2,0,0];
%             q_A(:,i) = [0,0,0];
%         end



t1 = linspace(1,100,51);

[P_r,R_r,V_r,V0_r,W_r,A_r,A0_r,Q,Q_V,Q_A,j_t] = generate_trajectory(0.1,5);
t = linspace(1,100,501);
tt = linspace(1,100,51);
figure(1)
subplot(1,3,1)
plot(t,P_r(1,:));
title('Position');
hold on
plot(tt, Q(1,:),'r.')
legend('trajectory','checkpoint')

subplot(1,3,2)
plot(t,V0_r(1,:));
title('Velocity');
hold on
plot(tt,Q_V(1,:),'r.')
legend('trajectory','checkpoint')

subplot(1,3,3)
plot(t,A_r(1,:));
title('Acceleration');
hold on
plot(tt,Q_A(1,:),'r.')
legend('trajectory','checkpoint')

figure(2)
plot(t,j_t);


%%%%% generate_trajectory %%%%
% 인수로 시간간격 dt와 traj_no 만을 받음
% 기존 인수로 받았던 N,speedx,scale,speed 등은 구간을 어떻게 설정하느냐에 따라 달라지므로 변수로 받지 않음
% 각 traj_no에서 check point의 수와 좌표, check point 사이 시간간격(m*dt)을 설정해 구간사이 속도를 설정할 수 있도록 함


function [P_D,R_D,V_D,V0_D,W_D,A_D,A0_D,q,Q_V,Q_A,j_t] = generate_trajectory(dt,traj_no)  
    if traj_no == 1
        q_size = 50;
        m = 10;
        dT = m * dt;
        q = zeros(3,q_size+1);
        q_V = zeros(3,q_size+1);
        q_A = zeros(3,q_size+1);
        q(:,1) = [0,0,0];  %%straight trajectory
        for i = 2 : q_size/5+1
            q(:,i) = [0,0,0.1*(i-1)*(i-1)];
           
        end
        for i = q_size/5+2 : 2*q_size/5+1
            q(:,i) = [0.1*(i-11)*(i-11),0,10];
           
        end
        for i = 2*q_size/5+2 : q_size+1
            q(:,i) = [10+2*(i-21),0,10];
           
        end
       for i = 1 : q_size+1
            if i ==1
                q_V(:,i) =  (q(:,i+1)-q(:,i))/dT;                                                
            elseif i == q_size+1
                q_V(:,i) =  (q(:,i)-q(:,i-1))/dT;
            else
                q_V(:,i) =  (q(:,i+1)-q(:,i-1))/2/dT;
            end
        end
        for i = 1 : q_size+1                  
            if i ==1
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i))/dT;                                                
            elseif i == q_size+1
                q_A(:,i) =  (q_V(:,i)-q_V(:,i-1))/dT;
            else
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i-1))/2/dT;
            end                            
        end

    elseif traj_no ==2
        q_size = 50;
        m = 10;
        dT = m * dt;
        q = zeros(3,q_size+1);
        q_V = zeros(3,q_size+1);
        q_A = zeros(3,q_size+1);
        for i = 1 : q_size/2 +1
            q(:,i) = [2*(i-1),0,0];  
        end
        for i = q_size/2+2 : q_size+1 
            q(:,i) = [50,0,0];  
        end

        for i = 1 : q_size+1
            if i ==1
                q_V(:,i) =  (q(:,i+1)-q(:,i))/dT;                                                
            elseif i == q_size+1
                q_V(:,i) =  (q(:,i)-q(:,i-1))/dT;
            else
                q_V(:,i) =  (q(:,i+1)-q(:,i-1))/2/dT;
            end
        end
        for i = 1 : q_size+1                  
            if i ==1
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i))/dT;                                                
            elseif i == q_size+1
                q_A(:,i) =  (q_V(:,i)-q_V(:,i-1))/dT;
            else
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i-1))/2/dT;
            end                            
        end
   
    elseif traj_no ==3
        q_size = 50;
        m = 20;
        q = zeros(3,q_size+1);
        q(:,1) = [0,0,0];  %%abrupt start/stop
        for i = 2 : q_size/5+1
            q(:,i) = [0,0,0];  
        end
        for i = q_size/5+2 : 2*q_size/5+1
            q(:,i) = [10,0,10];
        end
        for i = 2*q_size/5+2 : q_size+1
            q(:,i) = [0,0,0];
        end
    elseif traj_no == 4
        q_size = 100;
        m = 10;
        q = zeros(3,q_size+1);
        q(:,1) = [0,0,0];  %%straight trajectory
        for i = 2 : q_size/5+1
            q(:,i) = [0,0,0.1*(i-1)*(i-1)];  
        end
        for i = q_size/5+2 : q_size+1
            q(:,i) = [0.1*(i-11)*(i-11),0,10];
        end
     
    elseif traj_no == 5
        q_size = 50;
        m = 10;
        dT = m * dt;
        q = zeros(3,q_size+1);
        q_V = zeros(3,q_size+1);
        q_A = zeros(3,q_size+1);

        for i = 1 : q_size+1
           
            q(:,i) = [5*cos(pi/25*(i-1)),0,0];
        end
       
        for i = 1 : q_size+1
            if i ==1
                q_V(:,i) =  (q(:,i+1)-q(:,i))/dT;                                                
            elseif i == q_size+1
                q_V(:,i) =  (q(:,i)-q(:,i-1))/dT;
            else
                q_V(:,i) =  (q(:,i+1)-q(:,i-1))/2/dT;
            end
        end
        for i = 1 : q_size+1                  
            if i ==1
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i))/dT;                                                
            elseif i == q_size+1
                q_A(:,i) =  (q_V(:,i)-q_V(:,i-1))/dT;
            else
                q_A(:,i) =  (q_V(:,i+1)-q_V(:,i-1))/2/dT;
            end                            
        end
       
           


    end
        N = q_size*m;

        delta_x = zeros(3,q_size);
        delta_y = zeros(3,q_size);
        delta_z = zeros(3,q_size);
        coeff_x = zeros(3,q_size+1);
        coeff_y = zeros(3,q_size+1);
        coeff_z = zeros(3,q_size+1);

        s_t = zeros(3,N+1);
        v_t = zeros(3,N+1);
        a_t = zeros(3,N+1);
        j_t = zeros(3,N+1);
        s_t(:,1) =  q(:,1);

       
        for i = 1 : q_size
               
                delta_x(:,i) = [q(1,i+1)-q(1,i)-q_V(1,i)*dT-(1/2)*q_A(1,i)*(dT^2); q_V(1,i+1)-q_V(1,i)-q_A(1,i)*dT; q_A(1,i+1)-q_A(1,i)];
                delta_y(:,i) = [q(2,i+1)-q(2,i)-q_V(2,i)*dT-(1/2)*q_A(2,i)*(dT^2); q_V(2,i+1)-q_V(2,i)-q_A(2,i)*dT; q_A(2,i+1)-q_A(2,i)];
                delta_z(:,i) = [q(3,i+1)-q(3,i)-q_V(3,i)*dT-(1/2)*q_A(3,i)*(dT^2); q_V(3,i+1)-q_V(3,i)-q_A(3,i)*dT; q_A(3,i+1)-q_A(3,i)];
                coeff_x(:,i) = (1/(dT^5))*[720 -360*dT 60*(dT^2);-360*dT 168*(dT^2) -24*(dT^3);60*(dT^2) -24*(dT^3) 3*(dT^4)]*delta_x(:,i);
                coeff_y(:,i) = (1/(dT^5))*[720 -360*dT 60*(dT^2);-360*dT 168*(dT^2) -24*(dT^3);60*(dT^2) -24*(dT^3) 3*(dT^4)]*delta_y(:,i);
                coeff_z(:,i) = (1/(dT^5))*[720 -360*dT 60*(dT^2);-360*dT 168*(dT^2) -24*(dT^3);60*(dT^2) -24*(dT^3) 3*(dT^4)]*delta_z(:,i);
%                
%                 q_V(1,i+1) = (dT^4)*coeff_x(1,i)/24+(dT^3)*coeff_x(2,i)/6+(dT^2)*coeff_x(3,i)/2+q_A(1,i)*dT+q_V(1,i);
%                 q_V(2,i+1) = (dT^4)*coeff_y(1,i)/24+(dT^3)*coeff_y(2,i)/6+(dT^2)*coeff_y(3,i)/2+q_A(2,i)*dT+q_V(2,i);
%                 q_V(3,i+1) = (dT^4)*coeff_z(1,i)/24+(dT^3)*coeff_z(2,i)/6+(dT^2)*coeff_z(3,i)/2+q_A(3,i)*dT+q_V(3,i);
%
%                 q_A(1,i+1) = (dT^3)*coeff_x(1,i)/6+(dT^2)*coeff_x(2,i)/2+dT*coeff_x(3,i)+q_A(1,i);
%                 q_A(2,i+1) = (dT^3)*coeff_y(1,i)/6+(dT^2)*coeff_y(2,i)/2+dT*coeff_y(3,i)+q_A(2,i);
%                 q_A(3,i+1) = (dT^3)*coeff_z(1,i)/6+(dT^2)*coeff_z(2,i)/2+dT*coeff_z(3,i)+q_A(3,i);
               

        end
       
        for i = 1 : q_size+1
            Q_V(:,i) = q_V(:,i);
            Q_A(:,i) = q_A(:,i);
        end

        for i = 1 : N+1
            t= mod((i-1),m)*dt;
            j = 1+floor((i-1)/m);
           
            s_t(1,i) = (t^5)*coeff_x(1,j)/120+(t^4)*coeff_x(2,j)/24+(t^3)*coeff_x(3,j)/6+(1/2)*q_A(1,j)*(t^2)+q_V(1,j)*t+q(1,j);
            s_t(2,i) = (t^5)*coeff_y(1,j)/120+(t^4)*coeff_y(2,j)/24+(t^3)*coeff_y(3,j)/6+(1/2)*q_A(2,j)*(t^2)+q_V(2,j)*t+q(2,j);
            s_t(3,i) = (t^5)*coeff_z(1,j)/120+(t^4)*coeff_z(2,j)/24+(t^3)*coeff_z(3,j)/6+(1/2)*q_A(3,j)*(t^2)+q_V(3,j)*t+q(3,j);            
            v_t(1,i) = (t^4)*coeff_x(1,j)/24+(t^3)*coeff_x(2,j)/6+(t^2)*coeff_x(3,j)/2+q_A(1,j)*t+q_V(1,j);            
            v_t(2,i) = (t^4)*coeff_y(1,j)/24+(t^3)*coeff_y(2,j)/6+(t^2)*coeff_y(3,j)/2+q_A(2,j)*t+q_V(2,j);
            v_t(3,i) = (t^4)*coeff_z(1,j)/24+(t^3)*coeff_z(2,j)/6+(t^2)*coeff_z(3,j)/2+q_A(3,j)*t+q_V(3,j);
            a_t(1,i) = (t^3)*coeff_x(1,j)/6+(t^2)*coeff_x(2,j)/2+t*coeff_x(3,j)+q_A(1,j);
            a_t(2,i) = (t^3)*coeff_y(1,j)/6+(t^2)*coeff_y(2,j)/2+t*coeff_y(3,j)+q_A(2,j);
            a_t(3,i) = (t^3)*coeff_z(1,j)/6+(t^2)*coeff_z(2,j)/2+t*coeff_z(3,j)+q_A(3,j);
            j_t(1,i) = (t^2)*coeff_x(1,j)/2+coeff_x(2,j)*t+coeff_x(3,j);
        end
%           v_t(1,:) = lowpass(v_t(1,:),0.9);
         a_t(1,:) = lowpass(a_t(1,:),0.01);
        
        
        
%         for i = 1 : N+1
%             t= mod((i-1),m)*dt;
%             j = 1+floor((i-1)/m);
%             a_t(1,i) = (t^3)*coeff_x(1,j)/6+(t^2)*coeff_x(2,j)/2+t*coeff_x(3,j)+q_A(1,j);
%             a_t(2,i) = (t^3)*coeff_y(1,j)/6+(t^2)*coeff_y(2,j)/2+t*coeff_y(3,j)+q_A(2,j);
%             a_t(3,i) = (t^3)*coeff_z(1,j)/6+(t^2)*coeff_z(2,j)/2+t*coeff_z(3,j)+q_A(3,j);
%             j_t(1,i) = (t^2)*coeff_x(1,j)/2+coeff_x(2,j)*t+coeff_x(3,j);
%         end
        
%         for i = 1 : N
%             v_t(:,i) = (s_t(:,i+1)-s_t(:,i))/dt;    
%         end    
%         v_t(:,N+1) = (s_t(:,N+1)-s_t(:,N))/dt;
% 
% %         v_t(1,:) = lowpass(v_t(1,:),0.0001);
% 
%         for i = 1 : N        
%             a_t(:,i) = (v_t(:,i+1)-v_t(:,i))/dt;
%         end
%         a_t(:,N+1) = (v_t(:,N+1)-v_t(:,N))/dt;
% %         a_t(1,:) = lowpass(a_t(1,:),0.0001);
% 
%         for i = 1 : N        
%             j_t(:,i) = (a_t(:,i+1)-a_t(:,i))/dt;
%         end
%         j_t(:,N+1) = (a_t(:,N+1)-a_t(:,N))/dt;

        for i = 1 : N+1
            P_D(:,i) = [s_t(1,i); s_t(2,i); s_t(3,i)];
            R_D(:,:,i) = [1 0 0;
                          0 1 0;
                          0 0 1];
            V0_D(:,i) = [v_t(1,i); v_t(2,i); v_t(3,i)];
            V_D(:,i) = R_D(:,:,i)'*V0_D(:,i);
            A0_D(:,i) = [a_t(1,i); a_t(2,i); a_t(3,i)];
            A_D(:,i) = R_D(:,:,i)*A0_D(:,i);
            W_D(:,i) = [0 0 0]';

        end


end
