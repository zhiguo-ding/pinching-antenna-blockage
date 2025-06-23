clear
%close all

ct = 1000000;
snr=10000;
Mvec=[1:5];
D = 10;
height = 3; %d 
snrvec = [100:5 : 120 ];
M=2;
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % free space wavelength 
eta = (c/4/pi/f)^2;
D_leng =D*4;
leaky = 1;
blockage = 0.1;
neff = 1.4;%low-index materials like Teflon
lambda_g = lambda/neff;
feedpoint = -D_leng/2; %the feeder is at one side
for m = 1 : M %the y-axis locations of the waveguides (pinching antennas)
    betam(m,1) = -D/2+(m-1)*D/M+D/2/M;
end

parfor mi = 1 : length(snrvec) 
    snr = 10^(snrvec(mi)/10); 
    %conventional antenna
    r_all_low = zeros(M,ct);
    r_all_pin_low = zeros(M,ct);
    for i = 1 : ct    
        %M users, each in a strip, length D_length, width D/M, centered at betam(m)
        loc = zeros(M,2);
        for m = 1 : M
            loc(m,1) = D_leng*rand(1,1)-D_leng/2; %length   same as before     
            loc(m,2) = betam(m); %right underneath of wavegudie 
        end     
        
        %pinching antenna
        antenna_location = [loc(:,1) betam]; %x-corrdinates decided by users, y-corridantes decided by waveguides
        %find the distances - M columns, each column for one user, from user m to antenna k
        dist_all=zeros(M,M);
        for m = 1 : M
            dist_temp = loc(m,:)-antenna_location;%user m's location - loc(m,:)
            dist_all(:,m) = sqrt(sum(abs(dist_temp).^2,2) +height^2);
        end
        %use the distances to get the blockage coefficients, alpha
        blockage_probability_all = exp(-blockage*dist_all.^2); % the smaller the distance is, the larger this probability is
        alpha_all = rand(M,M) <= blockage_probability_all; %the larger the probability, closer to 1
        
        % %test
        % alpha11(i) = alpha_all(1,1);
        % alpha12(i) = alpha_all(2,1);

        %alpha_all = ones(M,M);
        %case III
        wd_loss = 10.^(0.08*(abs(loc(:,1)+D_leng/2))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %case II: waveguide loss
        %wd_loss = 10.^(0.08*(abs(antenna_location(:,1)))/10); %0.1dB/m, so 0.1x dB, P_ini/Pnew = 10^(0.1x/10)
        %wd_loss = ones(M,1);
        
        %build the channel matrix: each column for each user
        channel_temp = zeros(M,M);
        for um = 1 : M  % for each user m, it has the following M channels
            for ak = 1 : M  %from user m to antenna k         
                %in-waveguide phase shift
                theta_temp = 2*pi* (antenna_location(ak,1)-feedpoint)/lambda_g; %just assume 100 as the location of the feed
                %free space phase shift
                theta_channel = 2*pi*dist_all(ak,um)/lambda;
                %channel without blockage
                channel_temp(ak,um) =  sqrt(eta)*exp(-complex(0,1)*(theta_channel+theta_temp))/dist_all(ak,um)/sqrt(wd_loss(ak));
            end
        end
        H_channl = channel_temp.*alpha_all;% the composite channel with all effects
 

        %low complexity algorithms
        for m = 1 : M
            interference_blockage = 0;
            for im = 1 : M
                if im == m
                    continue
                else
                    interference_blockage = interference_blockage + abs(H_channl(im,m))^2;
                end
            end
            r_all_pin_low(m,i) = log2(1+abs(H_channl(m,m))^2*snr/(M+snr*interference_blockage));
        end

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conventional case: exactly same as before, just different antenna locations 
        antenna_location_c = [[1:M]'*lambda/2 zeros(M,1) ];
        %find the distances - M columns, each column for one user, from user m to antenna k
        dist_all_c = zeros(M,M);
        for m = 1 : M
            dist_temp = loc(m,:)-antenna_location_c;
            dist_all_c(:,m) = sqrt(sum(abs(dist_temp).^2,2) +height^2);
        end
        blockage_probability_all_c = exp(-blockage*dist_all_c.^2); % the smaller the distance is, the larger this probability is
        alpha_all_c = rand(M,M) <= blockage_probability_all_c;
        for m = 1 : M
            alpha_all_c(:,m) = alpha_all_c(m,m); %due to small antenna spacing, one blockage parameter per user
        end
     

        %alpha_all_c = ones(M,M);% <= blockage_probability_all_c;

        %build the channel matrix
        channel_temp_c=zeros(M,M);
        for um = 1 : M  % for each user m, it has the following M channels
            for ak = 1 : M  %from user m to antenna k        
                %only free-space phase shift. 
                theta_channel = 2*pi*dist_all_c(ak,um)/lambda;
                channel_temp_c(ak,um) =  sqrt(eta)*exp(-complex(0,1)*(theta_channel))/dist_all_c(ak,um);
            end
        end
        H_channl_c = channel_temp_c.*alpha_all_c;% the composite channel with all effects
 
        %low-complexity                
        for m = 1 : M
            interference_blockage = 0;
            for im = 1 : M
                if im == m
                    continue
                else
                    interference_blockage = interference_blockage + abs(H_channl_c(im,m))^2;
                end
            end
            %r_all_pin(m,i) = log2(1+eta*snr/dist_all(m,m)^2/wd_loss(m)/(1+interference));
            r_all_low(m,i) = log2(1+abs(H_channl_c(m,m))^2*snr/(M+snr*interference_blockage));
        end
    end

    %conventional antenna
    r_ave_sum_low(mi) = mean(r_all_low(1,:));%  sum rate
    %pinching antenna  
    r_ave_pin_sum_low(mi) = mean(r_all_pin_low(1,:)); 

    %analytical 
    tau4 = (betam(1)-betam(2))^2 + height^2;
    r_ana(mi) = 2*log2(1+eta*snr/M/(height^2))*exp(-blockage*(height^2))...
        *(1/2-exp(-blockage*tau4)/D_leng*sqrt(pi)/2/sqrt(blockage)*erf(sqrt(blockage)*D_leng)...
        +exp(-blockage*tau4)/D_leng^2/2/(blockage)*(1-exp(-blockage*D_leng^2)));

    % %acurate estimations
    % stepdl = D_leng/500;
    % dlvec = [-D_leng/2 : stepdl : D_leng/2];
    % sum1=0;
    % for i1 = 1 : length(dlvec)
    %     for i2 = 1 : length(dlvec)
    %         x1 = dlvec(i1);x2 = dlvec(i2);
    %         sum1 = sum1 + (1 - exp(-blockage* ((x1-x2)^2+tau4)))*stepdl^2/D_leng^2;
    %     end
    % end
    % r_test(mi) = log2(1+eta*snr/M/( height^2))...
    %     *exp(-blockage* (height^2))*sum1;

     
end
 
%plot(snrvec, r_ave_sum) 
%plot(snrvec, r_ave_sum_low, snrvec, r_ave_pin_sum_low, snrvec, r_ana  ) 
plot( snrvec, r_ave_pin_sum_low   ) 
 