
%% LOADING DATA 
% =========================================================================

% Loading SMF, PRO-seq and clustering data for drosophila 
Dd = readtable('data/Modelling/Drosophila_state_frequency.txt');
Td = readtable('data/Modelling/Drosophila_S2_TRP_promoter_counts_SpikeIn_norm_top5.txt');
Id = readtable('data/Modelling/Drosophila_K_means_clustering_TSS_cluster_mapping.txt');
Td = innerjoin(Td,Id,'Keys','gene');
Dd = innerjoin(Dd,Id,'Keys','gene');

% Loading SMF, PRO-seq and clustering data for mouse 
Dm = readtable('data/Modelling/Mouse_state_frequency.txt');
Tm = readtable('data/Modelling/mouse_TKO_TRP_promoter_counts_SpikeIn_norm_top5.txt');
Im = readtable('data/Modelling/mouse_K_means_clustering_TSS_cluster_mapping.txt');
Tm = innerjoin(Tm,Im,'Keys','gene');
Dm = innerjoin(Dm,Im,'Keys','gene');

% Calculating total counts for PRO-seq averageing over replicates (2.5min has only one replicate)
Ad = [(Td.S2_0min_R1+Td.S2_0min_R2)/2,   Td.S2_2_5min_R2, (Td.S2_5min_R1+Td.S2_5min_R2)/2, (Td.S2_10min_R1+Td.S2_10min_R2)/2, (Td.S2_20min_R1+Td.S2_20min_R2)/2];
Am = [(Tm.TKO_0min_R1+Tm.TKO_0min_R2)/2, Tm.TKO_2_5min_R1, (Tm.TKO_5min_R1+Tm.TKO_5min_R2)/2, (Tm.TKO_10min_R1+Tm.TKO_10min_R2)/2, (Tm.TKO_20min_R1+Tm.TKO_20min_R2)/2];

% Extracting cluster IDs into an array 
Cd = Td.cluster;
Cm = Tm.cluster;

% PRO-seq experimental times 
times = [0 2.5 5 10 20];
tmodel = 0:0.1:20;
numclus = 4;

% Species colors to match the style in other figures
color_m = [107,37,36]/256;
color_d = [41,42,97]/256;


%% Figure 4 (q vs kt/ki)
% =========================================================================

% Calculating probability of promoter free of nucleosome (Pf)
% probability of promoter bound by Pol II (Pp) and q = Pp/Pf
Pf_d = (100-Dd.nucleosome)./(100-Dd.unassigned);
Pp_d = (Dd.PIC_polII+Dd.polII)./(100-Dd.unassigned);
qd = Pp_d./Pf_d;

Pf_m = (100-Dm.nucleosome)./(100-Dm.unassigned);
Pp_m = (Dm.PIC_polII+Dm.polII)./(100-Dm.unassigned);
qm = Pp_m./Pf_m;

% Calculating the ratios r = kt/ki from q's 
rm = (1-qm)./qm;
rd = (1-qd)./qd;
q = 0.001:0.01:1;
r = (1-q)./q;
medqm = median(qm(isfinite(qm)));
medqd = median(qd(isfinite(qd)));
medrm = (1-medqm)/medqm;
medrd = (1-medqd)/medqd;

% Calculating distributions of q's and r's 
[yqm,xqm]=ksdensity(qm,'Support','nonnegative','BoundaryCorrection','reflection');
[yqd,xqd]=ksdensity(qd,'Support','nonnegative','BoundaryCorrection','reflection');
[yrm,xrm]=ksdensity(log10(rm(isfinite(rm))));
[yrd,xrd]=ksdensity(log10(rd(isfinite(rd))));

% Figure 4
fig4 = figure('Position',[0 0 400 300]);
subplot(3,3,[4 5 7 8])
semilogx(r,q,'Color',[0.5,0.5,0.5],'LineWidth',2)
hold on;
plot([medrd,medrd],[medqd,1],'--','Color',[0.5,0.5,0.5])
plot([medrm,medrm],[medqm,1],'--','Color',[0.5,0.5,0.5])
plot([medrd,medrd],[medqd,1],'--','Color',[0.5,0.5,0.5])
plot([medrd,10^3],[medqd,medqd],'--','Color',[0.5,0.5,0.5])
plot([medrm,10^3],[medqm,medqm],'--','Color',[0.5,0.5,0.5])
scatter(medrm,medqm,100,color_m,'filled','MarkerFaceAlpha',0.75)
scatter(medrd,medqd,100,color_d,'filled','MarkerFaceAlpha',0.75)
hold on
axis([10^-2 10^3 0 1])
ylabel('q = P_{PolII}/P_{open}')
xlabel('r = k_t / k_i')

subplot(3,3,[1 2])
hold on
plot(xrd,yrd,'LineWidth',3,'Color',color_d)
plot(xrm,yrm,'LineWidth',3,'Color',color_m)
xticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'})

subplot(3,3,[6 9])
hold on
plot(yqd,xqd,'LineWidth',3,'Color',color_d)
plot(yqm,xqm,'LineWidth',3,'Color',color_m)
axis([0 15 0 1])
saveas(fig4,'FiguresKasit2024/Kasit2024_figure4.eps','epsc');
saveas(fig4,'FiguresKasit2024/Kasit2024_figure4.pdf','pdf');


%% Figure 5 (Fitting turnover rates and initation rates)
% =========================================================================

% Fitting turnover rates in mouse 
fprintf('Mouse:\n')
[to_clus_m,hl_clus_m,sig_clus_m] = fit_turnover_rate(Am,Cm,times,numclus,'Mouse','FiguresKasit2024/Kasit2024_figsup_turnover_fits_mouse');

% Fitting turnover rates in drosophila 
fprintf('Drosophila:\n')
[to_clus_d,hl_clus_d,sig_clus_d] = fit_turnover_rate(Ad,Cd,times,numclus,'Drosophila','FiguresKasit2024/Kasit2024_figsup_turnover_fits_droso');

% Turnover time for each gene 
to_m = to_clus_m(Dm.cluster);
to_d = to_clus_d(Dd.cluster);
hl_m = hl_clus_m(Dm.cluster);
hl_d = hl_clus_d(Dd.cluster);

% Calulating initiation rates using q = 1/(1 + ti/to)
ti_m = to_m.*(1-qm)./qm;
ti_d = to_d.*(1-qd)./qd; 

% Kernel density estiamtion of the distribution of initation times
[yd,xd] = ksdensity(log10(ti_d(isfinite(ti_d))),'Bandwidth',0.15);
[ym,xm] = ksdensity(log10(ti_m(isfinite(ti_m))),'Bandwidth',0.15);

% Figure 5 (distributions of initiation times)
fig5 = figure('Position',[0 0 300 200]);
hold on
plot(xd,yd,'LineWidth',3,'Color',color_d);
plot(xm,ym,'LineWidth',3,'Color',color_m);
xticks([-1,0,1,2,3])
xticklabels({'0.1','1','10','100','1000'})
xlabel(' 1/k_i (min)')
ylabel('Probability density')
legend({'Drosophila','Mouse'})
legend('Location', 'northwest');
saveas(fig5,'FiguresKasit2024/Kasit2024_figure5.eps','epsc');
saveas(fig5,'FiguresKasit2024/Kasit2024_figure5.pdf','pdf');


%% Writing a table
% =========================================================================
table_mouse_q_to_ti_r = table(Dm.gene,qm,to_m,hl_m,ti_m,rm,'VariableNames',{'gene','q','turnover_half_life_[min]','turnover_time_[min]','initation_time_[min]','ratio'});
table_droso_q_to_ti_r = table(Dd.gene,qd,to_d,hl_d,ti_d,rd,'VariableNames',{'gene','q','turnover_half_life_[min]','turnover_time_[min]','initation_time_[min]','ratio'});
writetable(table_mouse_q_to_ti_r,'data/Modelling/output/Kasit2024_table_parameters_mouse.csv');
writetable(table_droso_q_to_ti_r,'data/Modelling/output/Kasit2024_table_parameters_drosophila.csv');

