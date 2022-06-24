%% plot of II_2
j=0;
figure(1+j);
    a1=plot(0.5*t1, sin(0.5*t1+pi/2), 'k-','LineWidth',1);
    hold on;
    a2=plot(X1(:,1), X1(:,2), '-','LineWidth', 2);
    hold on;
    a3=plot(X2(:,1), X2(:,2), '-','LineWidth', 2);
    hold on;
    a4=plot(X3(:,1), X3(:,2), '-','LineWidth', 2);
    hold on;
    a5=plot(Xa1(:,1), Xa1(:,2), '--','LineWidth', 2);
    hold on;
   % a5=plot(hat_Xa1(:,1), hat_Xa1(:,2), '--','LineWidth', 2);
    hold on;
    a6=plot(Xa2(:,1), Xa2(:,2), '--','LineWidth', 2);
    hold on;
    a7=plot(Xa3(:,1), Xa3(:,2), '--','LineWidth', 2);   
    xlabel('x[m]','FontSize',14);
    ylabel('y[m]','FontSize',14);
    
%marker 
    plot(X1(1,1), X1(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X2(1,1), X2(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X3(1,1), X3(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    
    plot(X1(80,1), X1(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X2(80,1), X2(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X3(80,1), X3(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X1(150,1), X1(150,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X2(150,1), X2(150,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(X3(150,1), X3(150,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    
% marker
    plot(hat_Xa1(1,1), hat_Xa1(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa2(1,1), hat_Xa2(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa3(1,1), hat_Xa3(1,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa1(80,1), hat_Xa1(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa2(80,1), hat_Xa2(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa3(80,1), hat_Xa3(80,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa1(151,1), hat_Xa1(151,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa2(151,1), hat_Xa2(151,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
    plot(hat_Xa3(151,1), hat_Xa3(151,2),'ks', 'MarkerSize',3,'MarkerFaceColor',[1 1 1]);
% lianxian LDMPC
    plot([X1(1,1), X2(1,1)],[X1(1,2), X2(1,2)],'r-.');
    plot([X2(1,1), X3(1,1)],[X2(1,2), X3(1,2)],'r-.');
    plot([X3(1,1), X1(1,1)],[X3(1,2), X1(1,2)],'r-.');

    plot([X1(80,1), X2(80,1)],[X1(80,2), X2(80,2)],'r-.');
    plot([X2(80,1), X3(80,1)],[X2(80,2), X3(80,2)],'r-.');
    plot([X3(80,1), X1(80,1)],[X3(80,2), X1(80,2)],'r-.');

    plot([X1(end,1), X2(end,1)],[X1(end,2), X2(end,2)],'r-.');
    plot([X2(end,1), X3(end,1)],[X2(end,2), X3(end,2)],'r-.');
    plot([X3(end,1), X1(end,1)],[X3(end,2), X1(end,2)],'r-.');
% lianxian 
    plot([hat_Xa1(1,1), hat_Xa2(1,1)],[hat_Xa1(1,2), hat_Xa2(1,2)],'b--');
    plot([hat_Xa2(1,1), hat_Xa3(1,1)],[hat_Xa2(1,2), hat_Xa3(1,2)],'b--');
    plot([hat_Xa3(1,1), hat_Xa1(1,1)],[hat_Xa3(1,2), hat_Xa1(1,2)],'b--');

    plot([hat_Xa1(80,1), hat_Xa2(80,1)],[hat_Xa1(80,2), hat_Xa2(80,2)],'b--');
    plot([hat_Xa2(80,1), hat_Xa3(80,1)],[hat_Xa2(80,2), hat_Xa3(80,2)],'b--');
    plot([hat_Xa3(80,1), hat_Xa1(80,1)],[hat_Xa3(80,2), hat_Xa1(80,2)],'b--');

    plot([hat_Xa1(151,1), hat_Xa2(151,1)],[hat_Xa1(151,2), hat_Xa2(151,2)],'b--');
    plot([hat_Xa2(151,1), hat_Xa3(151,1)],[hat_Xa2(151,2), hat_Xa3(151,2)],'b--');
    plot([hat_Xa3(151,1), hat_Xa1(151,1)],[hat_Xa3(151,2), hat_Xa1(151,2)],'b--');
    grid on;
% 放大
%     axes('Position',[0.68 0.17 0.2 0.18]);
%     p1=plot(x1(:,1), x1(:,2),'-');
%     set(p1,'LineWidth',2);
%     hold on
%     p2=plot(x2(:,1), x2(:,2),'-');
%     set(p2,'LineWidth',2);
%     p3=plot(x3(:,1), x3(:,2),'-');
%     set(p3,'LineWidth',2);
%     p5=plot(0.5*t1, sin(0.5*t1+pi/2),'k-');
%     set(p5,'LineWidth',1);
% 
%     plot(x1(6,1), x1(6,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[1 1 1]);
%     plot(x2(6,1), x2(6,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[1 1 1]);
%     plot(x3(6,1), x3(6,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[1 1 1]);
%     plot(x1(7,1), x1(7,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[0 0 0]);
%     plot(x2(7,1), x2(7,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[0 0 0]);
%     plot(x3(7,1), x3(7,2),'ks', 'MarkerSize',6,'MarkerFaceColor',[0 0 0]);
%     axis([-0.45 -0.2 -0.3 0.6])
    hold off;
    legend([a1,a2,a3,a4,a5,a6,a7],{'reference','AUV1', 'AUV2', 'AUV3','AUVa1', 'AUVa2', 'AUVa3'},'FontSize',12); 

     figure(5+j);
     subplot(3,1,1)
     plot(T*(1:150),(D1(:,1)-hat_Xa1(1:end-1,7)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D2(:,1)-hat_Xa2(1:end-1,7)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D3(:,1)-hat_Xa3(1:end-1,7)).^2,'LineWidth', 2);
     grid on;
     ylabel('||e_{1i}||^2','FontSize',14);
    axes('Position',[1.68 1.17 1.2 2.18]);
    p1=plot(T*(1:150), (D1(:,1)-hat_Xa1(1:end-1,7)).^2,'-');
    set(p1,'LineWidth',2);
    hold on
    p2=plot(T*(1:150),(D2(:,1)-hat_Xa2(1:end-1,7)).^2,'-');
    set(p2,'LineWidth',2);
    p3=plot(T*(1:150),(D3(:,1)-hat_Xa3(1:end-1,7)).^2,'-');
    set(p3,'LineWidth',2);
     

     subplot(3,1,2)
     plot(T*(1:150),(D1(:,2)-hat_Xa1(1:end-1,8)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D2(:,2)-hat_Xa2(1:end-1,8)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D3(:,2)-hat_Xa3(1:end-1,8)).^2,'LineWidth', 2);
     grid on;
     ylabel('||e_{2i}||^2','FontSize',14);
     subplot(3,1,3)
     plot(T*(1:150),(D1(:,3)-hat_Xa1(1:end-1,9)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D2(:,3)-hat_Xa2(1:end-1,9)).^2,'LineWidth', 2);
     hold on;
     plot(T*(1:150),(D3(:,3)-hat_Xa3(1:end-1,9)).^2,'LineWidth', 2);
     grid on;
     ylabel('||e_{3i}||^2','FontSize',14);
     xlabel('Time[s]','FontSize',14);
     legend({'AUV1','AUV2','AUV3'},'FontSize',12);

    figure(2+j);
    subplot(3,1,1)
    plot(t1, 0.5*t1, 'k-','LineWidth',1);
    hold on;
    plot(t1, X1(:,1), '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,1), '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,1), '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,1), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,1), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,1), '--','LineWidth', 2);
    hold on;
    grid on;
%     xlabel('Time[s]','interpreter','latex');
    ylabel('x[m]','FontSize',14);
    
    figure(2+j);
    subplot(3,1,2)
    plot(t1, sin(0.5*t1+pi/2), 'k-','LineWidth',1);
    hold on;
    plot(t1, X1(:,2), '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,2), '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,2), '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,2), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,2), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,2), '--','LineWidth', 2);
    hold on;
    grid on;
%     xlabel('Time[s]','interpreter','latex');
    ylabel('y[m]','FontSize',14);
    figure(2); 
    subplot(3,1,3)
    plot(t1, atan2(0.5*cos(0.5*t1+pi/2),0.5), 'k-','LineWidth',1);
    hold on;
    plot(t1, X1(:,3), '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,3), '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,3), '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,3), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,3), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,3), '--','LineWidth', 2);
    hold on;
    grid on;
%     leg2=legend('$AUV_1$', '$AUV_2$', '$AUV_3$','$AUV_{a1}$', '$AUV_{a2}$',...
%         '$AUV_{a3}$'); 
%     set(leg2,'interpreter','latex')
     legend({'AUV1','AUV2','AUV3','AUVa1','AUVa2','AUVa3'},'FontSize',12,'Orientation','horizontal');
    xlabel('Time[s]','FontSize',14);
    ylabel('\psi[rad]','FontSize',14);

    figure(3+j); % control input
     subplot(3,1,1)
     plot(t1(1:length(t1)-1),u1(:,1),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u2(:,1),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u3(:,1),'LineWidth', 2);
     grid on;
     hold on;
     plot(t1(1:length(t1)-1),Ua1(:,1), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua2(:,1), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua3(:,1), '--','LineWidth', 2);

%      xlabel('Time[s]','interpreter','latex');
     ylabel('u_1[N]','FontSize',14);
     
     subplot(3,1,2)
     plot(t1(1:length(t1)-1),u1(:,2),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u2(:,2),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u3(:,2),'LineWidth', 2);
     grid on;
     hold on;
     plot(t1(1:length(t1)-1),Ua1(:,2), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua2(:,2), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua3(:,2), '--','LineWidth', 2);
     ylabel('u_2[N]','FontSize',14);
     
     subplot(3,1,3)
     plot(t1(1:length(t1)-1),u1(:,3),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u2(:,3),'LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),u3(:,3),'LineWidth', 2);
     grid on;
     hold on;
     plot(t1(1:length(t1)-1),Ua1(:,3), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua2(:,3), '--','LineWidth', 2);
     hold on;
     plot(t1(1:length(t1)-1),Ua3(:,3), '--','LineWidth', 2);
%      leg3=legend('AUV_1', 'AUV_2', 'AUV_3','AUV_{a1}', 'AUV_{a2}',...
%         'AUV_{a3}'); 
%      set(leg3,'interpreter','latex')
     legend({'AUV1','AUV2','AUV3','AUVa1','AUVa2','AUVa3'},'FontSize',12,'Orientation','horizontal');
     xlabel('Time[s]','FontSize',14);
     ylabel('u_3[N\cdotm]','FontSize',14);
     
    figure(4+j);% error-x
    subplot(3,1,1)
    plot(t1, X1(:,1)-0.5*t1, '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,1)-0.5*t1, '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,1)-0.5*t1+1, '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,1)-0.5*t1, '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,1)-0.5*t1, '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,1)-0.5*t1+1, '--','LineWidth', 2);
    hold on;
    grid on;
    %xlabel('Time[s]','interpreter','latex');
    ylabel('x[m]','FontSize',14);
    
    figure(4);
    subplot(3,1,2)
    plot(t1, X1(:,2)-sin(0.5*t1+pi/2)-1, '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,2)-sin(0.5*t1+pi/2)+1, '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,2)-sin(0.5*t1+pi/2), '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,2)-sin(0.5*t1+pi/2)-1, '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,2)-sin(0.5*t1+pi/2)+1, '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,2)-sin(0.5*t1+pi/2), '--','LineWidth', 2);
    grid on;
    %legend('AUV1', 'AUV2', 'AUV3','AUVa1', 'AUVa2', 'AUVa3');
%     xlabel('Time[s]','interpreter','latex');
    ylabel('y[m]','FontSize',14);
    
    figure(4+j);
    subplot(3,1,3)
    plot(t1, X1(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '-','LineWidth', 2);
    hold on;
    plot(t1, X2(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '-','LineWidth', 2);
    hold on;
    plot(t1, X3(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '-','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa1(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa2(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '--','LineWidth', 2);
    hold on;
    plot(t1, hat_Xa3(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5), '--','LineWidth', 2);
    grid on;
    legend({'AUV1','AUV2','AUV3','AUVa1','AUVa2','AUVa3'},'FontSize',12,'Orientation','horizontal');
    xlabel('Time[s]','FontSize',14);
    ylabel('\psi[rad]','FontSize',14);
    
    figure(6+j); % the lyapunov function value
    plot(t1(1:end-1), mv1, '-','LineWidth', 2);
    hold on;
    plot(t1(1:end-1), mv2, '-','LineWidth', 2);
    hold on;
    plot(t1(1:end-1), mv3, '-','LineWidth', 2);
    grid on;
    hold on;
    plot(t1(1:end-1), lv1, '--','LineWidth', 2);
    hold on;
    plot(t1(1:end-1), lv2, '--','LineWidth', 2);
    hold on;
    plot(t1(1:end-1), lv3, '--','LineWidth', 2);
    hold on;
    plot(t1(1:end-1), lv2+lv1+lv3, 'k-.','LineWidth', 3);
    hold on;
    plot(t1(1:end-1), mv1+mv2+mv3, 'r-.','LineWidth', 3);
    legend({'AUV1','AUV2','AUV3','AUVa1','AUVa2','AUVa3','AUVa','AUV'},'FontSize',12);
    xlabel('Time[s]','FontSize',14);
    ylabel('V','FontSize',14);
    
    %formaiton tracking error of auxiliary contorller
    ea11 = hat_Xa1(:,1)-0.5*t1;
    ea21 = hat_Xa2(:,1)-0.5*t1;
    ea31 = hat_Xa3(:,1)-0.5*t1+1;
    ea12 = hat_Xa1(:,2)-sin(0.5*t1+pi/2)-1;
    ea22 = hat_Xa2(:,2)-sin(0.5*t1+pi/2)+1;
    ea32 = hat_Xa3(:,2)-sin(0.5*t1+pi/2);
    ea13 = hat_Xa1(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    ea23 = hat_Xa2(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    ea33 = hat_Xa3(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    ea1 = mean(ea11.^2+ea21.^2+ea31.^2);
    ea2 = mean(ea12.^2+ea22.^2+ea32.^2);
    ea3 = mean(ea13.^2+ea23.^2+ea33.^2);
    

    % formation tracking error of the LDMPC
    e11 = X1(:,1)-0.5*t1;
    e21 = X2(:,1)-0.5*t1;
    e31 = X3(:,1)-0.5*t1+1;
    e12 = X1(:,2)-sin(0.5*t1+pi/2)-1;
    e22 = X2(:,2)-sin(0.5*t1+pi/2)+1;
    e32 = X3(:,2)-sin(0.5*t1+pi/2);
    e13 = X1(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    e23 = X2(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    e33 = X3(:,3)-atan2(0.5*cos(0.5*t1+pi/2),0.5);
    e1 = mean(e11.^2+e21.^2+e31.^2);
    e2 = mean(e12.^2+e22.^2+e32.^2);
    e3 = mean(e13.^2+e23.^2+e33.^2);
% function [] = plot1( x,y,r )
% theta=0:0.1:2*pi;
% Circle1=x+r*cos(theta);
% Circle2=y+r*sin(theta);
% c=[123,14,52];
% plot(Circle1,Circle2,'c','linewidth',1);
% %axis equal
% end