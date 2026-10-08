% This code reproduces results and graphs in the paper. Curve 1 refers to fitting data 2010-17, Curve 2 to fitting data 2018-2026. 
% Cubic Spline Fits are based on F10.7 being the independent variable on x-axis and the image mean values (mn...) being the dependent 
% variable on y-axis. For display purporses in the paper the x and y axis are switched as one can not fit functions that
% begin to have infinite or close to it slopes.

%% Load data and define global variables
clear; close all;
load('F10_mnImages_B10Paper.mat'); % Load Daily F10.7 at 20:00UT and mean Image intensity at 20:00UT in 94A, 171A, 193A, 304A
evlpts  = 64:1:345;      % min to max F10.7 values 2010 - 2026
options = optimset('Display','iter','PlotFcns',@optimplotfval);  % options set for optimixation routine fmincon
f10c1   = f10(1:2754);   % F10.7 2010-2017
f10c2   = f10(2755:end); % F10.7 2018-2026

%% Step 1: Scale Curve 2 to Curve 1 for 193A by spline fitting each and multiplying curve 2 data by scale factor

% ----------------- Fit 193 Data 2010 - 2017: Curve 1 ---------------------
datac1 = mn193(1:2754); % mean 193

knots_0 = [66,73,92,195,210]';  % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,datac1,f10c1);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots1,L1errc1] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);  % call function minization routine
spline1  = csape(xknots1,datac1.'/fnval(csape(xknots1,eye(length(xknots1)),'var'),f10c1(:).'),'var');

% ----------------- Fit 193 Data 2017 - 2026: Curve 2 ---------------------
datac2 = mn193(2755:end);

knots_0 = [70,75,174,266,314]'; % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,datac2,f10c2);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots2,L1errc2] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minization routine
spline_c2  = csape(xknots2,datac2.'/fnval(csape(xknots2,eye(length(xknots2)),'var'),f10c2(:).'),'var');

% ----------------- Scale data of Curve 1 to Curve 2 193 and get mapping factor for data ---------------------
ScaleFactor_c1toc2      = (ppval(spline_c2,f10c1)./ppval(spline1,f10c1));
mapped_datac1_to_datac2 = ScaleFactor_c1toc2.*datac1;

% ----------------- Final Fit to All 193 Data ---------------------
mn193all = [mapped_datac1_to_datac2;datac2]; % different from mn193 in that it is mapped

knots_0 = [69,110,180,260,310]';  % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,mn193all,f10);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknotsall,L1err193] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minization routine
spline193 = csape(xknotsall,mn193all.'/fnval(csape(xknotsall,eye(length(xknotsall)),'var'),f10(:).'),'var');

%% Step 2: Scale Curve 2 to Curve 1 for 304A by spline fitting each and multiplying Curve 2 data by scale factor

%  ----------------- Fit 304 Data 2010 - 2017: Curve 1 ---------------------
datac1 = mn304(1:2754); % mean 193

knots_0 = [80,100,180,210]';  % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,datac1,f10c1);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots1,L1errc1] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minimization routine
spline1  = csape(xknots1,datac1.'/fnval(csape(xknots1,eye(length(xknots1)),'var'),f10c1(:).'),'var');

%  ----------------- Fit 304 Data 2017 - 2026: Curve 2 ---------------------
datac2 = mn304(2755:end);

knots_0 = [80,100,250,308]';   % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,datac2,f10c2);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots2,L1errc2] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minization routine
spline_c2  = csape(xknots2,datac2.'/fnval(csape(xknots2,eye(length(xknots2)),'var'),f10c2(:).'),'var');

%  ----------------- Scale data of Curve 1 to Curve 2 304 and get mapping factor for data  ---------------------
ScaleFactor_c1toc2       = (ppval(spline_c2,f10c1)./ppval(spline1,f10c1));
mapped_datac1_to_datac2  = ScaleFactor_c1toc2.*datac1;

%   ----------------- Fit Final 304 Data ---------------------
mn304all = [mapped_datac1_to_datac2;datac2];

knots_0 = [80,110,190,308]';   % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,mn304all,f10);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknotsall,L1err304] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minization routine
spline304 = csape(xknotsall,mn304all.'/fnval(csape(xknotsall,eye(length(xknotsall)),'var'),f10(:).'),'var');

%% Step 3: ----------------- Fit 94 and 171 Data ---------------------

% ----------------- Fit 171 ----------------------------------
datac1 = mn171(1:2754); %severe data spread of 171 for F10>230 causes large errors in fit

knots_0 = [68,90,175,205]';    % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,datac1,f10c1);        % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots1,L1err171] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);  % call function minization routine
spline171  = csape(xknots1,datac1.'/fnval(csape(xknots1,eye(length(xknots1)),'var'),f10c1(:).'),'var'); % create spline with knots from minimization algorithm

% ----------------- Fit 94 ----------------------------------
data = mn94;

knots_0 =  [72,90,140,300]';               % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,data,f10);            % L1 Error function, fmincon will call it
lb = [knots_0(1) -inf -inf  knots_0(end)]; % Fix lower bound lb 
ub = [knots_0(1)  inf  inf  knots_0(end)]; % Fix upper bound ub 

[xknots94,L1err94] = fmincon(fn,knots_0,[],[],[],[],lb,ub,[],options);   % call function minization routine
spline94  = csape(xknots94,data.'/fnval(csape(xknots94,eye(length(xknots94)),'var'),f10(:).'),'var');

%% Third Step: Scale all DATA (171,193,304) to 94

map171_to_94 = ppval(spline94,f10)./ppval(spline171,f10).* mn171;
map193_to_94 = ppval(spline94,f10)./ppval(spline193,f10).* mn193all;
map304_to_94 = ppval(spline94,f10)./ppval(spline304,f10).* mn304all;

%% Fit Spline to all data scaled to 94
alldata = [mn94; map171_to_94; map193_to_94; map304_to_94];
f10rep  = repmat(f10,4,1); % repeat F10.7 vector 4x to match length of 'data' in line above

knots_0_allData = [64,100,155,345];         % intial guess for knots of cubic spline; determined by minimization algorithm
fn = @(knots)f(knots,alldata,f10rep);                      % L1 Error function, fmincon will call it
lb = [knots_0_allData(1) -inf -inf  knots_0_allData(end)]; % Fix lower bound lb 
ub = [knots_0_allData(1)  inf  inf  knots_0_allData(end)]; % Fix upper bound ub 

[xknotsall,L1errall] = fmincon(fn,knots_0_allData,[],[],[],[],lb,ub,[],options);   % call function minization routine
splineall  = csape(xknotsall,alldata.'/fnval(csape(xknotsall,eye(length(xknotsall)),'var'),f10rep(:).'),'var'); % Final spline fit to all data

%% B10.7 values

y = linspace(min(f10),75); % ymax value of 75 satisfies criteria in paper
x = ppval(splineall,y);    % Sampled values along approx. linear part to <1% RMSE and R2 = 1
mdlmap = fitlm(x,y)        % Linear regression fit to check it meets criteria in paper
imean  = ppval(splineall,evlpts); % image brightness means from spline approx. given evlpts.
coef   = polyfit(x, y, 1);        % check polynomial linear fit same as a coefficients of regression 
b10    = polyval(coef, imean);  % Produce new B10.7
max(b10)

figure  % Plot Fig. 8 in paper
plot(alldata,f10rep,'.','MarkerSize',8','Color',[0.07,0.62,1.00]); hold on; grid;
plot(ppval(splineall,evlpts),evlpts,'LineWidth',2,'Color','k'); 
xplot = min(ppval(splineall,evlpts)):0.01:max(ppval(splineall,evlpts)); 
yplot = polyval(coef , xplot);
plot(xplot,yplot,'m-','LineWidth',2); %will plot on all pts plot

figure   % Plot Fig. 9 in paper
plot(evlpts,round(b10),'.') ;  

%% Function called by minization routine
function L1 = f(knots,data,f10)
    sp  = csape(knots,data.'/fnval (csape(knots,eye(length(knots)),'var'),f10(:).'),'var');
    for k = 1:162        % number of unique f10 values in least degraded part of SDO/AIA
        ind = find(f10 == min(f10)+k-1); % find the indices for the current value of f10
        id{k} = ind;                     % store those indices in a cell array
        for j = 1:numel(id{k})
            res{k}(j) = abs(data(id{k}(j)) - ppval(sp,64+k-1));
            d{k}(j)   = data(id{k}(j));
        end
    end

rescol = [res{:}]'; truth = [d{:}]';
L1 = abs(sum((rescol).^1))/abs(sum((truth).^1));   %calculate L1 Error of cubic spline fit
end