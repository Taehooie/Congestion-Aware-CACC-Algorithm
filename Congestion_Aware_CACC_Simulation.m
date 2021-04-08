clc;clear all;
%Given parameters for CACC, ACC, and human driven vehicles.
L = 5000;    %length of road
l = 5;       %effective length
vc = 1;      %cluster velocity
hc = 1;      %cluster headway
hf = 100;    %free headway for ACC and Human
chf = 300;   %free headway for CACC
vf = 25;     %free flow velocity
human = 0.4; %Driver sensitivity
Tr = 1.35;    %Truncation ratio
kh = (vf)/((hf)^(human)); %state dependent
t_leave = 4.5; %Time taken to leave a cluster

rho_steps = 100;
rho = linspace(0.1,0.8,rho_steps);
% Loop for CACC Penetration Rate (1->0 percent cacc, 11->100 percent cacc)
for pro=1:1:11
    % Loop for Density 
    for y=1:1:rho_steps        
        %Size of clusters
        N = round(rho(y)*L/l); %Total number of vehicles
        n = linspace(1,N,N);
        
        %Sigmoid parameters and function for driver sensitivity
        A = 0.4;     %The lower asymptote
        K = 0.8;     %The upper asymptote
        B = 0.06;     %growth rate
        v = 0.03;    %Starting point of growth
        Al(n) = A + (K-A)./(1+exp(-B.*n)).^(1/v);
        
        %Loop for h_cri values and K values, which is state dependent value
        for i=1:1:length(n)
            hcri(i) = chf.*(power((vc./vf),(1./Al(i))));
            KC(i)  = vf./(chf.^Al(i));
            %headway numerical value dependent on size of cluster
            hvalue(i,1:2000) = linspace(chf,1.001*hcri(i),2000);
        end
        
        %Loop for t_join
        m=1;%initial powerterm
        powerTerm=39;
        % Truncation ratio dependent on function of clustersize.
        for sizen=1:1:length(n)
            while m<=powerTerm
                t(m,1:length(hvalue)) = (1/vc).*(1./(1-m.*Al(sizen)).*(vc./KC(sizen))...
                    .^m.*(chf.^(1-m.*Al(sizen))-(hvalue(sizen,1:2000).^(1-m.*Al(sizen)))));
                m=m+1;
            end
            T40(sizen,1:length(hvalue)) = sum(t); %Last iteration of powerseries
            T1(sizen,1:length(hvalue)) = t(1,1:length(hvalue)); %First iteration of powerseries
            TR(sizen)=T40(sizen,2000)/T1(sizen,2000);
            m=1;
        end
        
        %Loop for proportion rate of CACC
        P=linspace(0,1,11);
        for p=1:1:11
            %Loop for the number of vehicles we count
            i=1;
            while i<=length(n)
                h(p,i)=(L - l*N - (n(i)-1)*hc)/(N - n(i) + 1);
                wplus_human(p,i) = kh*(1-human)/Tr*(1./(h(p,i).^(1-human)-hc^(1-human)));
                wplus_CACC(p,i) = KC(i).*(1-Al(i))./TR(i).*(1./(h(p,i).^(1-Al(i))-hc^(1-Al(i))));
                wplus(p,i) = (1-P(p)).*wplus_human(p,i) + P(p).*wplus_CACC(p,i);
                i=i+1;
            end
        end
        
        ze=find(wplus<0);
        wplus(ze)=0;
        one = find(wplus>1);
        wplus(one) = 1;
        wminus = 1/t_leave;
        %Monte-Carlo simulation and Random cluster size
        monte_steps = 5000;
        clustersize(1) = round(N*rand(1));
        RandomIN(1)=rand(1);
        RandomOUT(1)=rand(1);
        maximum_size = N;
        % Loop for Monte-Carlo mean
        for k=1:1:10
            % Loop for Monte-Carlo Simulation Steps
            for time=1:1:monte_steps
                if clustersize(time)==0 %lower bound
                    clustersize(time)=1;
                end
                %wplus is defined by clustersize at time
                if clustersize(time)>= maximum_size % upper bound
                    clustersize(time+1) = clustersize(time);
                elseif RandomIN(time)<wplus(pro,clustersize(time)) && RandomOUT(time)<wminus
                    clustersize(time+1) = clustersize(time);
                elseif RandomIN(time)>wplus(pro,clustersize(time)) && RandomOUT(time)>wminus
                    clustersize(time + 1) = clustersize(time);
                elseif RandomIN(time)>wplus(pro,clustersize(time)) && RandomOUT(time)<wminus
                    clustersize(time + 1) = clustersize(time) - 1;
                elseif RandomIN(time)<wplus(pro,clustersize(time)) && RandomOUT(time)>wminus
                    clustersize(time + 1) = clustersize(time) + 1;
                end
                RandomIN(time+1) = rand(1);
                RandomOUT(time+1)= rand(1);
            end
            mean_value(y,k) = mean(clustersize*l/L);
        end       
    end
    convert = transpose(mean_value);
    Result(:,pro) = convert(:);
end  
xlswrite('Monte_plotFile', Result)

[data] = xlsread('Monte_plotFile.xls');
x1 = data(:,1); % 0 percent CACC
x3 = data(:,3); % 20 percent CACC
x5 = data(:,5); % 40 percent CACC
x7 = data(:,7); % 60 percent CACC
x9 = data(:,9); % 80 percent CACC
x11 = data(:,11); % 100 percent CACC
rho = linspace(0,0.8,length(x1));
pe1 = linspace(0,0,length(x1));
pe3 = linspace(0.2,0.2,length(x3));
pe5 = linspace(0.4,0.4,length(x5));
pe7 = linspace(0.6,0.6,length(x7));
pe9 = linspace(0.8,0.8,length(x9));
pe11 = linspace(1,1,length(x11));
plot(rho,x1,rho,x3,rho,x5,rho,x7,rho,x9,rho,x11);
xlabel('Normalized density')
ylabel('Normalized cluster size')
grid on;