%% plot of II_2
% Define the color
cor0 = [0, 0.4470, 0.7410];
cor1 = [0.8500, 0.3250, 0.0980];
cor2 = [0.9290, 0.6940, 0.1250];
cor3 = [0.4940, 0.1840, 0.5560];
cor4 = [0.4660, 0.6740, 0.1880];
cor5 = [0.3010, 0.7450, 0.9330];
cor6 = [0.6350, 0.0780, 0.1840];
[T0, T1, T2, T3, T4, T5, T6] = deal(TIme(1:Tstep)');
% Specify the shading area along the optimal nominal state trajectory
up0 = X0(1,1:Tstep) + 0.5;
up1 = X1(1,1:Tstep) + 0.5;
up2 = X2(1,1:Tstep) + 0.5;
up3 = X3(1,1:Tstep) + 0.5;
up4 = X4(1,1:Tstep) + 0.5;
up5 = X5(1,1:Tstep) + 0.5;
up6 = X6(1,1:Tstep)+ 0.5;
down0 = X0(1,1:Tstep) - 0.5;
down1 = X1(1,1:Tstep) - 0.5;
down2 = X2(1,1:Tstep)- 0.5;
down3 = X3(1,1:Tstep) - 0.5;
down4 = X4(1,1:Tstep) - 0.5;
down5 = X5(1,1:Tstep) - 0.5;
down6 = X6(1,1:Tstep) - 0.5;


% fig1: state trajectory of the actual and nominal
figure(1);
subplot(6,1,1:4)
c0= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down0 fliplr(up0)],cor0,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c1= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down1 fliplr(up1)],cor1,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c2= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down2 fliplr(up2)],cor2,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c3= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down3 fliplr(up3)],cor3,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c4= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down4 fliplr(up4)],cor4,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c5= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down5 fliplr(up5)],cor5,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
c6= fill([T0(1:Tstep) fliplr(T0(1:Tstep))],[down6 fliplr(up6)],cor6,'FaceAlpha',0.2,'EdgeColor','none');
hold all;
plot(T1, Xreal0(1,1:Tstep),'-', 'color',cor0,'LineWidth', 1);
hold on;
a1=plot(T1, Xreal1(1,1:Tstep), '-', 'color',cor1,'LineWidth', 1);
hold on;
a2=plot(T1, Xreal2(1,1:Tstep), '-', 'color',cor2,'LineWidth', 1);
hold on;
a3=plot(T1, Xreal3(1,1:Tstep), '-', 'color',cor3,'LineWidth', 1);
hold on;
a4=plot(T1, Xreal4(1,1:Tstep), '-', 'color',cor4,'LineWidth', 1);
hold on;
a5=plot(T1, Xreal5(1,1:Tstep), '-', 'color',cor5,'LineWidth', 1);
hold on;
a6=plot(T1, Xreal6(1,1:Tstep), '-', 'color',cor6,'LineWidth', 1);
hold on;
b0=plot(T1, X0(1,1:Tstep), ':', 'color',cor0,'LineWidth', 1.5);
hold on;
b1=plot(T1, X1(1,1:Tstep), ':', 'color',cor1,'LineWidth', 1.5);
hold on;
b2=plot(T1, X2(1,1:Tstep), ':', 'color',cor2,'LineWidth', 1.5);
hold on;
b3=plot(T1, X3(1,1:Tstep), ':', 'color',cor3,'LineWidth', 1.5);
hold on;
b4=plot(T1, X4(1,1:Tstep), ':', 'color',cor4,'LineWidth', 1.5);
hold on;
b5=plot(T1, X5(1,1:Tstep), ':', 'color',cor5,'LineWidth', 1.5);
hold on;
b6=plot(T1, X6(1,1:Tstep), ':', 'color',cor6,'LineWidth', 1.5);
grid on;
set(gca,'xticklabel',[]);
ylabel('${\boldmath x_i}$[m]','FontSize',16,'interpreter','latex');

subplot(6,1,5:6)
a0=plot(T1, Xreal0(2,1:Tstep),'-', 'color',cor0,'LineWidth', 1);
hold on;
a1=plot(T1, Xreal1(2,1:Tstep), '-', 'color',cor1,'LineWidth', 1);
hold on;
a2=plot(T1, Xreal2(2,1:Tstep), '-', 'color',cor2,'LineWidth', 1);
hold on;
a3=plot(T1, Xreal3(2,1:Tstep), '-', 'color',cor3,'LineWidth', 1);
hold on;
a4=plot(T1, Xreal4(2,1:Tstep), '-', 'color',cor4,'LineWidth', 1);
hold on;
a5=plot(T1, Xreal5(2,1:Tstep), '-', 'color',cor5,'LineWidth', 1);
hold on;
a6=plot(T1, Xreal6(2,1:Tstep), '-', 'color',cor6,'LineWidth', 1);
hold on;
b0=plot(T1, X0(2,1:Tstep), ':', 'color',cor0,'LineWidth', 1.5);
hold on;
b1=plot(T1, X1(2,1:Tstep), ':', 'color',cor1,'LineWidth', 1.5);
hold on;
b2=plot(T1, X2(2,1:Tstep), ':', 'color',cor2,'LineWidth', 1.5);
hold on;
b3=plot(T1, X3(2,1:Tstep), ':', 'color',cor3,'LineWidth', 1.5);
hold on;
b4=plot(T1, X4(2,1:Tstep), ':', 'color',cor4,'LineWidth', 1.5);
hold on;
b5=plot(T1, X5(2,1:Tstep), ':', 'color',cor5,'LineWidth', 1.5);
hold on;
b6=plot(T1, X6(2,1:Tstep), ':', 'color',cor6,'LineWidth', 1.5);
patch([1.6 2 2 1.6],[0 0 2 2],cor1,'FaceAlpha',0.25,'EdgeColor','none');
text(2,0.5,'\leftarrow Acceleration');
patch([3.9 4.3 4.3 3.9],[0 0 2 2],cor4,'FaceAlpha',0.25,'EdgeColor','none');
text(4.3,0.4,'\leftarrow Deceleration');
grid on;
legend([a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6],...
    {'$q_0$','$q_1$','$q_2$','$q_3$','$q_4$','$q_5$','$q_6$','$\tilde{q}_0$',...
    '$\bar{q}_1$','$\bar{q}_2$','$\bar{q}_3$','$\bar{q}_4$','$\bar{q}_5$','$\bar{q}_6$'},...
    'FontSize',14,'interpreter','latex');
xlabel('Time[s]','FontSize',16,'interpreter','latex');
ylabel('$\mu_i$[m/s]','FontSize',16,'interpreter','latex');
% Fig.2: control input
figure(2);
plot(T1,U0(:),'color',cor0,'LineWidth', 1.5);
hold on;
plot(T1,U1(:),'color',cor1,'LineWidth', 1.5);
hold on;
plot(T1,U2(:),'color',cor2,'LineWidth', 1.5);
hold on;
plot(T1,U3(:),'color',cor3,'LineWidth', 1.5);
hold on;
plot(T1,U4(:),'color',cor4,'LineWidth', 1.5);
hold on;
plot(T1,U5(:),'color',cor5,'LineWidth', 1.5);
hold on;
plot(T1,U6(:),'color',cor6,'LineWidth', 1.5);
grid on;
ylabel('$u_i$[N]','FontSize',14,'interpreter','latex');
legend({'ASV0','ASV1','ASV2','ASV3','ASV4','ASV5','ASV6'},...
    'FontSize',12,'interpreter','latex');
xlabel('Time[s]','FontSize',14,'interpreter','latex');



% Fig.3: error of the nominal states error-x
figure(3);% error-x
subplot(8,1,1:2);
plot(T1, (X0(1,1:Tstep)-Ref(1,1:Tstep)), '-','color',cor0,'LineWidth', 1.5);
hold on;
plot(T1, (X1(1,1:Tstep)-Ref(1,1:Tstep))+2.5, '-','color',cor1,'LineWidth', 1.5);
hold on;
plot(T1, (X2(1,1:Tstep)-Ref(1,1:Tstep))+5, '-','color',cor2,'LineWidth', 1.5);
hold on;
plot(T1, (X3(1,1:Tstep)-Ref(1,1:Tstep))+7.5, '-','color',cor3,'LineWidth', 1.5);
hold on;
plot(T1, (X4(1,1:Tstep)-Ref(1,1:Tstep))+10, '-','color',cor4,'LineWidth', 1.5);
hold on;
plot(T1, (X5(1,1:Tstep)-Ref(1,1:Tstep))+12.5, '-','color',cor5,'LineWidth', 1.5);
hold on;
plot(T1, (X6(1,1:Tstep)-Ref(1,1:Tstep))+15, '-','color',cor6,'LineWidth', 1.5);
grid on;
set(gca,'xticklabel',[]);
ylabel('$\bar{x}_i^e$[m]','FontSize',14,'interpreter','latex');

% Fig.: error of the nominal velocity error-x
figure(3)
subplot(8,1,3:4);
plot(T1, (X0(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X1(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X2(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X3(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X4(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X5(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
hold on;
plot(T1, (X6(2,1:Tstep)-Ref(2,1:Tstep)), '-','LineWidth', 1.5);
grid on;
set(gca,'xticklabel',[]);
ylabel('$\bar{\mu}_i^e$[m/s]','FontSize',14,'interpreter','latex');
figure(3);
subplot(8,1,5:6);
plot(T1, (Xreal0(1,1:Tstep)-X0(1,1:Tstep)), '-','color',cor0,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal1(1,1:Tstep)-X1(1,1:Tstep)), '-','color',cor1,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal2(1,1:Tstep)-X2(1,1:Tstep)), '-','color',cor2,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal3(1,1:Tstep)-X3(1,1:Tstep)), '-','color',cor3,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal4(1,1:Tstep)-X4(1,1:Tstep)), '-','color',cor4,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal5(1,1:Tstep)-X5(1,1:Tstep)), '-','color',cor5,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal6(1,1:Tstep)-X6(1,1:Tstep)), '-','color',cor6,'LineWidth', 1.5);
grid on;
set(gca,'xticklabel',[]);
ylabel('$e_i^x$[m]','FontSize',14,'interpreter','latex');

figure(3);
subplot(8,1,7:8);
plot(T1, (Xreal0(2,1:Tstep)-X0(2,1:Tstep)), '-','color',cor0,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal1(2,1:Tstep)-X1(2,1:Tstep)), '-','color',cor1,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal2(2,1:Tstep)-X2(2,1:Tstep)), '-','color',cor2,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal3(2,1:Tstep)-X3(2,1:Tstep)), '-','color',cor3,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal4(2,1:Tstep)-X4(2,1:Tstep)), '-','color',cor4,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal5(2,1:Tstep)-X5(2,1:Tstep)), '-','color',cor5,'LineWidth', 1.5);
hold on;
plot(T1, (Xreal6(2,1:Tstep)-X6(2,1:Tstep)), '-','color',cor6,'LineWidth', 1.5);
grid on;
legend({'ASV0','ASV1','ASV2','ASV3','ASV4','ASV5','ASV6'},'FontSize',12,'interpreter','latex');
xlabel('{Time[s]}','FontSize',14,'interpreter','latex');
ylabel('${e}_i^{\mu}$[m/s]','FontSize',14,'interpreter','latex');

