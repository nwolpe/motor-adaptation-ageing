%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Script for (robust) fitting of trajectory errors in Cam-CAN 
%%%% visuomotor learning task with two exponentials 
%%%%
%%%% Written by James Ingram & Noham Wolpe 2018-2019
%%%% correspondence to: nw305@medschl.cam.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
FilePath = '~/Dropbox/MotorLearning/';
FileName = 'CamCAN-VML.dat';

load(sprintf('%s%s',FilePath,FileName), '-mat');

SubjectCount = length(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotFlag = true;
% PlotFlag = false;
PDF_Flag = false;

% for (sub)plotting:
fn = 0;
pm = 4;
pn = 4;
pc = pm*pn;
pp = 0;

opt = statset;
opt.RobustWgtFun = 'bisquare'; % robust weighting function = bisquare
beta0 = [0.5 10 10]; % initial values for fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fitting

for k=1:SubjectCount
    % Trial values:
    s1=data(k).PhaseTrialTrajectoryError(2,1:120)'; % Exposure phase
    s2=data(k).PhaseTrialTrajectoryError(3,1:48)';  % Post-Exposure phase
    
    % Cycle values (mean of 4 trials):
    y1=mean(reshape(s1',4,length(s1)/4))';          % Exposure phase
    y2=mean(reshape(s2',4,length(s2)/4))';          % Post-Exposure phase
    
    y = [y1; y2];         % Trajectory Error - cycle-wise
    x = (0:length(y)-1)';   % Cycle Number
    
    [beta, r] = nlinfit(x,y,@exponentials_fit_constrained,beta0,opt);
    
    beta(1) = beta(1);
    beta(2) = abs(beta(2));         % Adaptation time-constant
    beta(3) = abs(beta(3));         % De-adaptation time-constant
    
    yhat = exponentials_fit_constrained(beta,x);
    
    % Calculate final de-adaptation
    xi = linspace(0,41,1000);
    yp = exponentials_fit_constrained(beta,xi);
    beta(5) = yp(end);              % Final de-adaptation.
    
    % Time to half adaptation
    j = find(xi <= 29);
    mid_a = (beta(1)+yp(j(1)))/2;
    i = find(yp(j) <= mid_a,1);
    beta(4) = 1+xi(j(i));           % Time-to-half adaptation.
    
    % Time to half de-adaptation
    j = find(xi >= 30);
    mid_d = (beta(5)+yp(j(1)))/2;
    i = find(yp(j) >= mid_d,1);
    beta(6) = xi(j(i))-29;          % Time-to-half de-adaptation.
    
    params(k,:) = beta;
    rsq(k) = 1 - sum(r.^2) / sum((y - mean(y)).^2); % goodness of fit R2
    
    if PlotFlag
        if( (pp == pc) || (pp == 0) )
            fn = fn + 1;
            figure(fn);
            clf;
            
            pp = 0;
        end
        
        pp = pp + 1;
        subplot(pm,pn,pp);
        plot(x+1,y,'k.')
        hold on
        plot(x+1,yhat,'r-')
        axis([ 0 42 -30 60 ])
        xy = axis;
        X = [ xy(1) xy(2) ];
        Y = [ 0 0 ];
        plot(X,Y,'k:');
        plot(beta(4),mid_a,'r.');
        plot(30+beta(6),mid_d,'r.');
        
        str = sprintf('%i Aend=%2.1f t(0.5)=%2.1f Dend=%.2f',k,beta(1),beta(4),beta(5));
        title(str);
        
%         if( PDF_Flag && ((rem(k,16) == 0) || (k == SubjectCount)) )
%             shg
%             pause
%             %export_fig -append -painters results.pdf
%             clf
%         end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot parameter histograms

close all;

figure(1);
clf;

pm = 2;
pn = 3;
pp = 0;

pp = pp + 1;
subplot(pm,pn,pp);
x = params(:,1);
[y,bx] = hist(x,30);
bar(bx,y,'k');
axis([ (min(bx)-2) (max(bx)+2) 0 50 ]);
title(sprintf('mean=%.1f std=%.1f (deg)',mean(x),std(x)));
axis([ (min(x)-2) (max(x)+2) 0 50 ]);
xlabel('Exposure Final Offset (deg)');
ylabel('Count');

pp = pp + 1;
subplot(pm,pn,pp);
x = params(:,4);
[y,bx] = hist(x,30);
bar(bx,y,'k');
axis([ (min(bx)-0.5) (max(bx)+0.5) 0 50 ]);
title(sprintf('mean=%.1f std=%.1f (cycles)',mean(x),std(x)));
xlabel('Exposure Time-to-half (cycles)');
ylabel('Count');

pp = pp + 1;
subplot(pm,pn,pp);
x = params(:,5);
[y,bx] = hist(x,30);
bar(bx,y,'k');
axis([ (min(bx)-2) (max(bx)+2) 0 50 ]);
title(sprintf('mean=%.1f std=%.1f (deg)',mean(x),std(x)));
axis([ (min(x)-2) (max(x)+2) 0 50 ]);
xlabel('Post-Exposure Final Offset (deg)');
ylabel('Count');

pp = pp + 1;
subplot(pm,pn,pp);
x = params(:,6);
[y,bx] = hist(x,30);
bar(bx,y,'k');
axis([ (min(bx)-0.5) (max(bx)+0.5) 0 50 ]);
title(sprintf('mean=%.1f std=%.1f (cycles)',mean(x),std(x)));
xlabel('Post-Exposure Time-to-half (cycles)');
ylabel('Count');

pp = pp + 1;
subplot(pm,pn,pp);
x = rsq;
[y,bx] = hist(x,30);
bar(bx,y,'k');
axis([ (min(bx)-0.1) (max(bx)+0.1) 0 50 ]);
title(sprintf('mean=%.1f std=%.1f (deg)',mean(x),std(x)));
xlabel('Rsquared of fit');
ylabel('Count');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot time to half adaptation against exponential time constant

figure(2);
clf;

pm = 2;
pn = 2;
pp = 0;

pp = pp + 1;
subplot(pm,pn,pp);
hold on;
x = params(:,2); % Exposure time-to-half.
y = params(:,4); % Post-Exposure time-to-half.
i = find((x < 100) & (y < 100));
x = x(i);
y = y(i);
plot(x,y,'k.');

xlabel('Exponential Time Constant (cycles)');
ylabel('Time-to-half Adpatation (cycles)');
title('Exposure');

pp = pp + 1;
subplot(pm,pn,pp);
hold on;
x = params(:,3); % Exposure time-to-half.
y = params(:,6); % Post-Exposure time-to-half.
i = find((x < 100) & (y < 100));
x = x(i);
y = y(i);
plot(x,y,'k.');

xlabel('Exponential Time Constant (cycles)');
ylabel('Time-to-half Adpatation (cycles)');
title('Post-Exposure');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% store parameters
D = zeros(SubjectCount,7);
% D(:,1) = CCid
% D(:,2) = Final adaptation (deg).
% D(:,3) = Adaptation time-constant (cycles).
% D(:,4) = De-Adaptation time-constant (cycles).
% D(:,5) = Time-to-half final adaptation (cycles).
% D(:,6) = Final de-apdatation (deg).
% D(:,7) = Time-to-half final de-adaptation (cycles).
% D(:,8) = R^2 of fit.

for i=1:SubjectCount
    D(i,1) = str2num(data(i).CCid);
end

D(:,2:7) = params;
D(:,8) = rsq;

% FileName = 'CamCAN-VML-FitParameters.mat';
% save(sprintf('%s%s',FilePath,FileName),'D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%