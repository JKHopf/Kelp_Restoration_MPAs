function [N, Nbio, Ybio, Rs] = Func_OverTimePlot(ND, tx, lincol, SpParas)

subplot(5,1,1)
        hold on
%         plot(tx, sum(ND(:,:),1), lincol,...
%             'DisplayName' , "Yr = "+yrDist+", Lgth = "+rDist)
        plot(tx, sum(ND(:,:),1)./sum(ND(:,1),1),lincol)
        title('Total metapopulation')
        xlim([0,100])
%         ylim([0.9,1.3])
    subplot(5,1,2)
        hold on
%         plot(tx, sum(ND(1:SpParas.A_max,:),1),lincol)
        plot(tx, sum(ND(1:SpParas.A_max,:),1)./sum(ND(1:SpParas.A_max,1),1),lincol)
        title('Reserve population (restore)')
        xlim([0,100])
%         ylim([0.9,1.3])
    subplot(5,1,3)
        hold on
%         plot(tx, sum(ND(SpParas.A_max+1:SpParas.A_max*2,:),1),lincol)
        plot(tx, sum(ND(SpParas.A_max+1:SpParas.A_max*2,:),1)./sum(ND(SpParas.A_max+1:SpParas.A_max*2,1),1),lincol)
        title('Reserve population')
        xlim([0,100])
%         ylim([0.9,1.3])
    subplot(5,1,4)
        hold on
%         plot(tx, sum(ND(2*SpParas.A_max+1:SpParas.A_max*3,:),1),lincol)
        plot(tx, sum(ND(2*SpParas.A_max+1:SpParas.A_max*3,:),1)./sum(ND(2*SpParas.A_max+1:SpParas.A_max*3,1),1),lincol)
        title('Fished population (restore)')
        xlim([0,100])
%         ylim([0.9,1.3])
    subplot(5,1,5)
        hold on
%         plot(tx, sum(ND(3*SpParas.A_max+1:SpParas.A_max*4,:),1),lincol)
        plot(tx, sum(ND(3*SpParas.A_max+1:SpParas.A_max*4,:),1)./sum(ND(3*SpParas.A_max+1:SpParas.A_max*4,1),1),lincol)
        title('Fished population')
        xlim([0,100])
%         ylim([0.9,1.3])