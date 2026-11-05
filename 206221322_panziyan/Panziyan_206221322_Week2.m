%% 潘子岩 206221322：储能选址定容与交流校核
% 独立实现。依据课程第二周模板：无光伏、LinDistFlow MILP、Big-M线性化。
% 时段t对应(t-1):00至t:00；功率MW，能量MWh，费用元/年。
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'data'));
prepare_tools(root);
yalmip('clear');
Stu_ID = 206221322; student_name = '潘子岩';
alpha = 1 + 0.002*mod(Stu_ID,100);
sites = [4 9 14];
network = loadcase('case_nongcun');
profile = load(fullfile(root,'data','profiles_24h.mat'));
cfg = struct('etaC',0.90,'etaD',0.95,'socMin',0.20,'socMax',0.80,...
    'soc0',0.50,'pOptions',1:5,'eOptions',1:5,'vMin',0.960,'vMax',1.050,...
    'powerCost',80000,'energyCost',550000,'omCost',200000,'days',365,...
    'discount',0.05,'life',20,'dt',1);
fprintf('%s | %d | alpha=%.3f | 候选母线 %s\n',student_name,Stu_ID,alpha,mat2str(sites));

%% 1. 基准方案与规划优化
baseline = ac_series(network,profile,alpha,sites,zeros(3,profile.T));
plan = solve_planning(network,profile,alpha,sites,cfg);
assert(plan.exitCode==0,'规划未达到求解器最优终止条件');
checked = ac_series(network,profile,alpha,sites,plan.discharge-plan.charge);
audit = verify_plan(network,profile,alpha,sites,cfg,plan,checked);
assert(audit.passed,'数值校验失败，请检查 audit 输出');
planning = table(sites',round(plan.installed),plan.power,plan.energy,...
    'VariableNames',{'Bus','Installed','Power_MW','Energy_MWh'});
disp(planning); disp(audit);

%% 2. 费用采用两个明确口径：规划目标不含网损，交流核算含网损
price = profile.price_MWh(:); hours=(1:profile.T)';
baseGrid = sum(network.bus(:,3))*alpha*profile.load_curve_pu(:);
baseCost = cfg.days*sum(price.*baseGrid)*cfg.dt;
lossBenefit = cfg.days*sum(price.*(baseline.loss-checked.loss))*cfg.dt;
trueBaseCost = cfg.days*sum(price.*baseline.grid)*cfg.dt;
truePlanCost = sum(plan.cost(1:2))+cfg.days*sum(price.*checked.grid)*cfg.dt;
costs = table([0;plan.cost(1)],[0;plan.cost(2)],[baseCost;plan.cost(3)],...
    [trueBaseCost;truePlanCost],'RowNames',{'Base','ESS'},...
    'VariableNames',{'AnnualInvestment','AnnualOM','PlanningGridCost','AC_TotalCost'});
disp(costs);

%% 3. 保留课程规定的标准结果字段
Result = struct('StuID',Stu_ID,'Obj_Value',sum(plan.cost),...
    'Opt_Sites',find(round(plan.installed)>0),...
    'Opt_P_rated',plan.power,'Opt_E_rated',plan.energy,...
    'Max_V_Err',audit.maxVoltageError);
save(fullfile(root,'result_summary.mat'),'Result');
save(fullfile(root,'calculation_details.mat'),'plan','baseline','checked',...
    'audit','cfg','planning','costs','lossBenefit','alpha','sites','student_name');
writetable(planning,fullfile(root,'Planning_206221322.csv'));
hourly = table(hours,baseGrid,plan.grid,checked.grid,sum(plan.charge,1)',...
    sum(plan.discharge,1)',min(baseline.voltage,[],1)',min(checked.voltage,[],1)',...
    'VariableNames',{'Hour','BaseLoad_MW','LinGrid_MW','ACGrid_MW',...
    'Charge_MW','Discharge_MW','BaseMinVoltage_pu','ESSMinVoltage_pu'});
writetable(hourly,fullfile(root,'Hourly_206221322.csv'));

%% 4. 独立的四组结果图
draw_results(root,Stu_ID,student_name,sites,cfg,plan,baseline,checked,price,baseCost);
fprintf('\n完成：选址 %s，总功率 %.0f MW，总容量 %.0f MWh。\n',...
    mat2str(sites(round(plan.installed)>0)),sum(plan.power),sum(plan.energy));
fprintf('交流最低电压 %.6f p.u.，最大线性误差 %.6f p.u.\n',...
    min(checked.voltage,[],'all'),audit.maxVoltageError);
fprintf('规划目标 %.2f 元/年；交流核算总成本 %.2f 元/年。\n',sum(plan.cost),truePlanCost);
fprintf('数值自检通过；图、程序和结果位于 %s\n',root);

function prepare_tools(root)
% 优先使用已配置的工具箱；支持把依赖放在本文件夹 dependencies 下。
    candidates={fullfile(root,'dependencies'),...
        '/Users/apple/Documents/Codex/MATLAB_Toolboxes',...
        '/Users/apple/Documents/Codex/2026-10-30/meta-lab/outputs/ESS_Day3_Mac/dependencies'};
    for k=1:numel(candidates)
        if isfolder(candidates{k}), addpath(genpath(candidates{k})); end
    end
    addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');
    assert(exist('yalmip','file')==2,'需要安装并添加 YALMIP 路径');
    assert(exist('runpf','file')==2,'需要安装并添加 MATPOWER 路径');
end

function out=solve_planning(m,prof,alpha,sites,c)
% 按档位展开的SoC乘积，四条Big-M约束保证二进制乘积精确等价。
    T=prof.T; n=size(m.bus,1); b=size(m.branch,1); ns=numel(sites); nk=numel(c.pOptions);
    fb=m.branch(:,1); tb=m.branch(:,2); sb=m.baseMVA;
    loadP=m.bus(:,3)*alpha/sb*prof.load_curve_pu(:)';
    loadQ=m.bus(:,4)*alpha/sb*prof.load_curve_pu(:)';
    choose=binvar(ns,nk,'full'); built=binvar(ns,1);
    state=sdpvar(ns,T,'full'); product=sdpvar(ns,nk,T,'full');
    charge=sdpvar(ns,T,'full'); discharge=sdpvar(ns,T,'full'); mode=binvar(ns,T,'full');
    squaredV=sdpvar(n,T,'full'); flowP=sdpvar(b,T,'full'); flowQ=sdpvar(b,T,'full');
    gridP=sdpvar(1,T,'full'); gridQ=sdpvar(1,T,'full');
    power=choose*c.pOptions'; energy=choose*c.eOptions'; maxPower=max(c.pOptions);
    con=[sum(choose,2)==built,c.socMin<=state<=c.socMax,state(:,T)==c.soc0,...
        0<=charge,0<=discharge,charge<=repmat(power,1,T),...
        discharge<=repmat(power,1,T),charge<=maxPower*mode,...
        discharge<=maxPower*(1-mode)];
    for i=1:ns
        for k=1:nk
            d=reshape(product(i,k,:),1,T);
            con=[con,c.socMin*choose(i,k)<=d<=c.socMax*choose(i,k),...
                state(i,:)-c.socMax*(1-choose(i,k))<=d,...
                d<=state(i,:)-c.socMin*(1-choose(i,k))]; %#ok<AGROW>
        end
        stored=c.eOptions*reshape(product(i,:,:),nk,T);
        previous=[c.soc0*energy(i),stored(1:T-1)];
        con=[con,stored-previous==c.dt*(c.etaC*charge(i,:)-discharge(i,:)/c.etaD)]; %#ok<AGROW>
    end
    incidence=sparse(tb,1:b,1,n,b)-sparse(fb,1:b,1,n,b);
    injection=sparse(sites,1:ns,1,n,ns); slack=sparse(1,1,1,n,1);
    con=[con,incidence*flowP+slack*gridP+injection*(discharge-charge)/sb==loadP,...
        incidence*flowQ+slack*gridQ==loadQ,squaredV(1,:)==1,...
        c.vMin^2<=squaredV(2:end,:)<=c.vMax^2];
    for j=1:b
        con=[con,squaredV(tb(j),:)==squaredV(fb(j),:)-...
            2*(m.branch(j,3)*flowP(j,:)+m.branch(j,4)*flowQ(j,:))]; %#ok<AGROW>
    end
    % 关口出力同时满足课程算例中的发电机上下限。
    con=[con,m.gen(1,10)/sb<=gridP<=m.gen(1,9)/sb,...
        m.gen(1,5)/sb<=gridQ<=m.gen(1,4)/sb];
    crf=c.discount*(1+c.discount)^c.life/((1+c.discount)^c.life-1);
    investment=crf*sum(c.powerCost*power+c.energyCost*energy);
    maintenance=c.omCost*sum(power);
    gridCost=c.days*c.dt*sum(prof.price_MWh(:)'.*(sb*gridP));
    settings=sdpsettings('solver','intlinprog','verbose',0);
    settings.intlinprog.RelativeGapTolerance=1e-7;
    settings.intlinprog.AbsoluteGapTolerance=1e-6;
    timer=tic; result=optimize(con,investment+maintenance+gridCost,settings);
    assert(result.problem==0,'MILP求解失败：%s',result.info);
    out=struct('installed',value(built),'power',value(power),'energy',value(energy),...
        'soc',value(state),'charge',value(charge),'discharge',value(discharge),...
        'voltage',sqrt(value(squaredV)),'grid',sb*value(gridP)',...
        'flowP',value(flowP),'flowQ',value(flowQ),'gridQ',sb*value(gridQ)',...
        'cost',[value(investment),value(maintenance),value(gridCost)],...
        'exitCode',result.problem,'solveSeconds',toc(timer));
end

function r=ac_series(m,prof,alpha,sites,net)
    T=prof.T; r.voltage=zeros(size(m.bus,1),T); r.loss=zeros(T,1); r.grid=zeros(T,1);
    opt=mpoption('verbose',0,'out.all',0);
    for t=1:T
        sample=m; sample.bus(:,3:4)=m.bus(:,3:4)*alpha*prof.load_curve_pu(t);
        sample.bus(sites,3)=sample.bus(sites,3)-net(:,t);
        pf=runpf(sample,opt); assert(pf.success,'第%d小时交流潮流未收敛',t);
        r.voltage(:,t)=pf.bus(:,8); r.grid(t)=sum(pf.gen(:,2));
        r.loss(t)=sum(real(get_losses(pf)));
    end
end

function a=verify_plan(m,prof,alpha,sites,c,p,ac)
    energy=p.energy.*p.soc; previous=[p.energy*c.soc0,energy(:,1:end-1)];
    a.energyResidual=max(abs(energy-previous-c.dt*(c.etaC*p.charge-p.discharge/c.etaD)),[],'all');
    a.cycleResidual=max(abs(p.soc(:,end)-c.soc0));
    a.simultaneousPower=max(min(p.charge,p.discharge),[],'all');
    a.maxVoltageError=max(abs(p.voltage-ac.voltage),[],'all');
    n=size(m.bus,1); b=size(m.branch,1);
    A=sparse(m.branch(:,2),1:b,1,n,b)-sparse(m.branch(:,1),1:b,1,n,b);
    J=sparse(sites,1:3,1,n,3); S=sparse(1,1,1,n,1);
    residual=A*p.flowP+S*p.grid'/m.baseMVA+J*(p.discharge-p.charge)/m.baseMVA-...
        m.bus(:,3)*alpha/m.baseMVA*prof.load_curve_pu(:)';
    a.powerBalanceResidual=max(abs(residual),[],'all');
    a.passed=a.energyResidual<1e-5 && a.cycleResidual<1e-5 && a.simultaneousPower<1e-5 &&...
        a.powerBalanceResidual<1e-7 && all(p.soc>=c.socMin-1e-5 & p.soc<=c.socMax+1e-5,'all') &&...
        all(p.charge>=-1e-5 & p.discharge>=-1e-5,'all') &&...
        all(p.charge<=p.power+1e-5 & p.discharge<=p.power+1e-5,'all') &&...
        min(ac.voltage,[],'all')>=0.95-1e-5 && max(ac.voltage,[],'all')<=1.05+1e-5;
end

function draw_results(root,id,name,sites,c,p,base,ac,price,baseCost)
    set(groot,'defaultAxesFontName','Songti SC','defaultTextFontName','Songti SC',...
        'defaultAxesFontSize',11,'defaultLineLineWidth',1.7);
    blue=[0.12 0.36 0.55]; orange=[0.86 0.38 0.14]; green=[0.20 0.55 0.42];
    T=numel(price); t=1:T; active=find(round(p.installed)>0);
    f=new_figure([120 80 1200 760]); layout=tiledlayout(f,2,2,'TileSpacing','compact');
    nexttile; bar(sites,[p.power,p.energy],'grouped'); xlabel('候选母线'); ylabel('MW / MWh');
    xticks(sites); legend('额定功率 / MW','额定能量 / MWh','Location','northwest'); title('选址与额定配置'); grid on;
    nexttile; bar(t,[sum(p.discharge,1)',-sum(p.charge,1)'],'stacked');
    xlabel('时段 / h'); ylabel('功率 / MW'); legend('放电','充电'); title('储能总充放电功率'); grid on;
    nexttile; hold on;
    for i=active', plot(0:T,[c.soc0,p.soc(i,:)]*100,'DisplayName',sprintf('母线%d',sites(i))); end
    yline(c.socMin*100,':','HandleVisibility','off'); yline(c.socMax*100,':','HandleVisibility','off');
    ylim([0 100]); xlim([0 T]); xlabel('时刻 / h'); ylabel('SOC / %'); legend('Location','best'); title('SOC：含初始时刻'); grid on;
    nexttile; yyaxis left; plot(t,p.grid,'Color',blue); ylabel('关口购电 / MW');
    yyaxis right; stairs(t,price,'Color',orange); ylabel('电价 / (元/MWh)'); xlabel('时段 / h');
    title('购电计划与分时电价'); grid on;
    title(layout,sprintf('%s · %d | 储能配置与运行',name,id)); emit(f,root,sprintf('Fig_Week2_%d',id));

    f=new_figure([120 80 1200 720]); layout=tiledlayout(f,2,2,'TileSpacing','compact');
    [~,weak]=min(min(base.voltage,[],2));
    nexttile; plot(t,base.voltage(weak,:),'--','Color',[.45 .45 .45]); hold on;
    plot(t,p.voltage(weak,:),'-','Color',blue); plot(t,ac.voltage(weak,:),'-o','Color',orange,'MarkerSize',3);
    yline(.95,':r'); ylabel('电压 / p.u.'); xlabel('时段 / h'); grid on;
    title(sprintf('薄弱母线%d：电压校核',weak)); legend('无储能AC','规划Lin','配置后AC','下限','Location','best');
    nexttile; plot(t,min(base.voltage,[],1),'--','Color',[.45 .45 .45]); hold on;
    plot(t,min(ac.voltage,[],1),'-','Color',green); yline(.95,':r'); grid on;
    title('全网逐时最低电压'); ylabel('电压 / p.u.'); xlabel('时段 / h'); legend('无储能','配置后','下限');
    nexttile; imagesc(t,1:size(ac.voltage,1),ac.voltage); axis xy; colorbar;
    xlabel('时段 / h'); ylabel('母线编号'); title('配置后交流电压分布');
    nexttile; imagesc(t,1:size(ac.voltage,1),1000*(p.voltage-ac.voltage)); axis xy; colorbar;
    xlabel('时段 / h'); ylabel('母线编号'); title('线性模型减交流电压 / 10^{-3} p.u.');
    title(layout,sprintf('%s · %d | 交流潮流验证',name,id)); emit(f,root,sprintf('Fig_Week2_Voltage_%d',id));

    f=new_figure([120 80 1100 500]); layout=tiledlayout(f,1,2,'TileSpacing','compact');
    nexttile; bar([0 0 baseCost;p.cost]/1e4,'stacked'); xticklabels({'无储能','储能方案'});
    ylabel('万元/年'); title('规划目标口径：不含线路损耗'); legend('投资年值','年运维','购电','Location','northwest'); grid on;
    nexttile; bar(p.cost(1:2)/1e4); xticklabels({'投资年值','年运维'}); ylabel('万元/年');
    title('储能固定费用明细'); grid on;
    title(layout,sprintf('%s · %d | 年费用构成',name,id)); emit(f,root,sprintf('Fig_Week2_Cost_%d',id));

    f=new_figure([120 80 1100 500]); layout=tiledlayout(f,1,2,'TileSpacing','compact');
    nexttile; plot(t,base.loss,'--','Color',[.45 .45 .45]); hold on; plot(t,ac.loss,'Color',green);
    xlabel('时段 / h'); ylabel('网损 / MW'); title('逐时有功损耗'); legend('无储能','配置后'); grid on;
    nexttile; baseAC=c.days*sum(price.*base.grid); essAC=sum(p.cost(1:2))+c.days*sum(price.*ac.grid);
    bar([baseAC essAC]/1e4); xticklabels({'无储能','储能方案'}); ylabel('万元/年');
    title('交流核算总费用：含实际线路损耗'); grid on;
    title(layout,sprintf('%s · %d | 网损与实际费用',name,id)); emit(f,root,sprintf('Fig_Week2_ACCost_%d',id));
end

function f=new_figure(position)
    f=figure('Color','white','Position',position,'Visible','off');
    if isprop(f,'Theme'), theme(f,'light'); end
end
function emit(f,root,stem)
    exportgraphics(f,fullfile(root,[stem '.png']),'Resolution',240,'BackgroundColor','white');
    exportgraphics(f,fullfile(root,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
    savefig(f,fullfile(root,[stem '.fig'])); close(f);
    fprintf('已输出 %s\n',stem);
end
