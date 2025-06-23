clear
%close all

ct = 5000000;
snr=10000;
Mvec=[1:5];
D = 5;
height = 3; %d 
snrvec = [100:5 : 120];
M=1;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =4*D;
leaky = 1;
blockage = 0.1;
Rtar = 5;
eps = 2^Rtar-1; 
Dlvec = [1:1:5];

for mi = 1 : length(Dlvec) 
    D_leng =Dlvec(mi)*D;
    snr = 10^(snrvec(1)/10); 
    tau1 = sqrt(eta*snr/eps);
    sum1=0;sum2=0;sum3=0;sum4=0;
    for i = 1 : ct    
        loc = zeros(M,2);
        loc(:,1) = D_leng*rand(M,1)-D_leng/2; %length        
        loc(:,2) = D*rand(M,1)-D/2; %width,      

        %conventional antenna
        dist = max(1,sqrt(loc(:,1).^2+ loc(:,2).^2 +height^2));
        blockage_probability = exp(-blockage*dist); % the larger the distance is, the smaller this probability is
        blockage_phi = rand(M,1) <= blockage_probability; % smaller probability, i.e., 0, give a number close to 0, always blocked
        r_all(:,i) = log2(1+eta.*snr./dist.^2)/M;      
        r_all_blockage(:,i) = log2(1+eta*blockage_phi.*snr./dist.^2)/M;  
        if r_all(:,i)< Rtar
            sum1 = sum1 + 1;
        end
        if r_all_blockage(:,i)< Rtar
            sum2 = sum2 + 1;
        end
        
        %pinching antenna
        %case III
        wd_loss = 10.^(0.08*(abs(loc(:,1)+D_leng/2))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %case II
        %wd_loss = 10.^(0.08*(abs(loc(:,1)))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %wd_loss = 1; %no loss for now
        dist2 = max(1,sqrt(loc(:,2).^2 +height^2));
        blockage_probability2 = exp(-blockage*dist2); % the smaller the distance is, the larger this probability is
        blockage_phi2 = rand(M,1) <= blockage_probability2; % larger probability, i.e., 1, give a number close to 1, not blocked
  
        r_all_pin(:,i) = log2(1+eta*snr./dist2.^2./wd_loss)/M;   
        r_all_pin_blockage(:,i) = log2(1+eta*blockage_phi2.*snr./dist2.^2./wd_loss)/M;   
        if r_all_pin(:,i)< Rtar
            sum3 = sum3 + 1;
        end
        if r_all_pin_blockage(:,i)< Rtar
            sum4 = sum4 + 1;
        end        
    end

    %conventional antenna
    r_ave_sum(mi) = mean(sum(r_all,1));%  sum rate  
    r_ave_sum_blockage(mi) = mean(sum(r_all_blockage,1));%  sum rate  
    %pinching antenna 
    r_ave_pin_sum(mi) = mean(sum(r_all_pin,1)); %sum rate  
    r_ave_pin_sum_blockage(mi) = mean(sum(r_all_pin_blockage,1)); %sum rate 

    %outage probability
    po_conv(mi) = sum1/ct;
    po_conv_block(mi) = sum2/ct;
    po_pin(mi) = sum3/ct;
    po_pin_block(mi) = sum4/ct;
 
    %analytical results
    term1 = 0;
    stepsize = D/2/1000;
    yvec = [-D/2:stepsize:D/2];
    for iy = 1 : length(yvec)
        y = yvec(iy);
        term1 = term1 + (1-exp(-blockage*sqrt(y^2+height^2)))*stepsize/D;
    end
    term2 = 0; 
    tau2 = max(-D/2,-sqrt(tau1^2-height^2));
    stepsize = (tau2+D/2)/1000;
    yvec = [-D/2:stepsize:tau2];
    for iy = 1 : length(yvec)
        y = yvec(iy);
        term2 = term2 + exp(-blockage*sqrt(y^2+height^2))*stepsize/D;
    end
    term3 = 0; 
    tau3 = min(D/2,sqrt(tau1^2-height^2));
    stepsize = (D/2-tau3)/1000;
    yvec = [tau3:stepsize:D/2];
    for iy = 1 : length(yvec)
        y = yvec(iy);
        term3 = term3 + exp(-blockage*sqrt(y^2+height^2))*stepsize/D;
    end
    Pana(mi) = term1+term2+term3;

    %approximation
    Papp(mi) = term1;

    % %conventional
    % term4 = 0;
    % stepsizey = D/2/1000;
    % stepsizex = D_leng/2/1000;
    % yvec = [0:stepsizey:D/2];
    % xvec = [0:stepsizex:D_leng/2];
    % for iy = 1 : length(yvec)
    %     for ix = 1 : length(xvec)
    %         y = yvec(iy);
    %         x = xvec(ix);
    %         term4 = term4 + (1-exp(-blockage*sqrt(x^2+y^2+height^2)))*stepsizey/D*stepsizex/D_leng;
    %     end
    % end
    % Pappconv(mi) = term4;
end
 
%plot(snrvec, r_ave_sum) 
%plot(snrvec, r_ave_sum,  snrvec, r_ave_pin_sum,snrvec, r_ave_sum_blockage,  snrvec, r_ave_pin_sum_blockage ) 
%[(r_ave_pin_sum-r_ave_sum)./r_ave_sum; (r_ave_pin_sum_blockage-r_ave_sum_blockage)./r_ave_sum_blockage]
semilogy(Dlvec, po_conv_block, Dlvec, po_pin_block, Dlvec, Pana)
 