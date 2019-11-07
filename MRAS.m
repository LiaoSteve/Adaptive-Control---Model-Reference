% MARS 
% J(theta1,theta2)=0.5*error^2
clc,clear;
a=1;b=0.5;
am=2;bm=2;
R=[0.2 1 5];%gamma
uc(1:100)=1;
uc(101:200)=-1;
uc=[uc uc uc uc uc];
uc=[uc uc uc uc uc];  
uc=[uc uc uc uc uc]; 
ts=0.1;
time=500;
for s=1:3
    r=R(s);
    theta1(1)=0;
    theta2(1)=0;

    theta1_prime=0;
    theta2_prime=0;

    y(1)=0;
    ym(1)=0;    
    i=0;
    for k=0:ts:(time-ts)
        i=i+1;
        u(i)=theta1(i)*uc(i)-theta2(i)*y(i);

        y_prime=-a*y(i)+b*u(i);
        ym_prime=-am*ym(i)+bm*uc(i);

        y(i+1)=y(i)+y_prime*ts;
        ym(i+1)=ym(i)+ym_prime*ts;

        e=y(i+1)-ym(i+1);    

        theta1_double_prime=( -am*theta1_prime-r*am*uc(i+1)*e );
        theta2_double_prime=( -am*theta2_prime+r*am*y(i+1)*e  );

        theta1_prime=theta1_prime+theta1_double_prime*ts;
        theta2_prime=theta2_prime+theta2_double_prime*ts;

        theta1(i+1)=theta1(i)+theta1_prime*ts;
        theta2(i+1)=theta2(i)+theta2_prime*ts;
    end
    figure,
    subplot(411)
    plot(0:ts:time,y); 
    hold on
    plot(0:ts:time,ym);
    suptitle(['MRAS : \gamma=',num2str(r)])
    ylabel('y,ym')
    text(1.5,0.2145,' y')
    text(1.7, 1.4,'ym')
    %axis([-inf, inf, -1.5, 1.5])
    
    subplot(412)
    plot(0:ts:time-ts,u);
    ylabel('u')
    text(65, 1,' u')
    axis([-inf, inf, -7, 7])
    
    subplot(413)
    plot(0:ts:time,theta1);
    ylabel('\theta_1'),title(['\theta_1(100)=',num2str(theta1(1001))])
    
    subplot(414)
    plot(0:ts:time,theta2);
    xlabel('time'),ylabel('\theta_2'),title(['\theta_2(100)=',num2str(theta2(1001))])
    if r==1 && time>=500
        figure
        plot(theta1,theta2),xlabel('\theta_1'),ylabel('\theta_2')
        title('Relation between \theta_1 and \theta_2')
        hold on,plot(theta1,theta1-(a/b))
        axis([-inf, inf, -1, inf])
    end    
end
