function  y = exponentials_fit_constrained(b,x) 

final = b(1);
tau1=abs(b(2));
tau2=abs(b(3));
ang=30; % trajectory error on first cycle constrained to be 30 = perturbation angle


if b(1)<=0 || b(1)>=30 % constraining final adaptation offset to be 0<b<30
    y = zeros(length(x),1);
else
    F=(final-ang*exp(-30/tau1))/(1-exp(-30/tau1));
    
    % Exposure
    i = find((x >= 0) & (x <= 29)); % exposure phase 1<cycle<30
    y(i)=F - (F-ang).*exp(-x(i)./tau1);
    yend = y(i(end));
    
    % Post-Exposure
    i = find(x > 29); % post exposure phase cycle>30
    y(i)=(yend-ang)*exp(-(x(i)-30)/tau2);
    y=y';
end