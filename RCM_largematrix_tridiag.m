% *************************************************************************
% This program creates the ensemble of Normalized z from the Random Matrix
% Hamiltonian for apertures radiating into 3D chaotic cavities
% *************************************************************************

clc;
clear all;
%close all;

%%% IOTA: the imaginary unit
j = sqrt(-1);

%%% Important for the critical parameter \mu = \rho / \alpha
%alpha =10;
rho = 1.0;

%%% ***********************************************************************
%%% Definition of the parameters used for generating the universal
%%% fluctuation term
%%% ***********************************************************************
 %%% Size of cavity Hamiltonian (number of chaotic modes)
 Ns = 1000; 
 for nn=1:1
     nn
    %%% IOTA: the imaginary unit
    j = sqrt(-1);

    %%% Important for the critical parameter \mu = \rho / \alpha
%     alpha = 0.1^nn;
    alpha = 10000;
    rho = 1.0;

    %%% ***********************************************************************
    %%% Definition of the parameters used for generating the universal
    %%% fluctuation term
    %%% ***********************************************************************
     %%% Size of cavity Hamiltonian (number of chaotic modes)
     Ns = 50000; 
     %%% Number of radiating elements (antennas/ports)
     Nt = 1;
     Nr = 1;

     samples=1; % Number of runs per iteration. 
     runs=10000; % Number of iterations. Final size of output will contain "runs*samples" number of znorm entries. 

    %%% Radiation impedance matrix antennas/ports, loads, generator **********************
    Z_rad_t = 1*eye(Nt);
    Z_rad_r = 1*eye(Nr);
    Z_L = 1*eye(Nr);
    Z_G = 1*eye(Nt);
    %%% ***********************************************************************

    %%% ***********************************************************************
    %%% Calculation of fluctuating cavity impedance 
    %%% ***********************************************************************

    eigen = zeros(Ns,samples); % Initializing the eigen matrix for cavity 1 - 8 aperture modes

    mask1=zeros(Ns);% Creating hamiltonian mask for diagonal terms
    mask2=zeros(Ns);% Creating hamiltonian mask for off-diagonal terms
    for i=1:Ns
        for ind=1:Ns
            if(i==ind)
              mask1(i,ind)=1/2;
            elseif(i<ind)
              mask2(i,ind)=1;
            elseif (i>ind)
              0;
            end
        end
    end

    %%% ***********************************************************************
    %%% ***********************************************************************
    inds = 1;
    for p=1:runs
            p
            for cnt=1:samples
                beta=1;      	%GOE
                d0=normrnd(0,sqrt(2),1,Ns)';                                     %Diagonal Term
                d1=sqrt(chi2rnd(beta*((Ns-1):-1:1)))';                           %Off-diagonal Terms
                hamilt=spdiags([[d1;0],d0,[0;d1]],[-1,0,1],Ns,Ns)/sqrt(2);        %Sparse Matrix
                eigen(:,cnt) = eig(hamilt);% computing eigenvalues of hamiltonian for cavity 
            end

            %%% Mapping the eigenvalues from the semi-circle to have uniform spacing
            Elevel=(Ns/(2*pi))*(pi+2*asin(eigen./sqrt(2*Ns))+2.*(eigen./sqrt(2*Ns)).*sqrt(2*Ns-eigen.^2)/sqrt(2*Ns))-Ns/2;
            ElevAppo(p,:)=eigen;

            CoupleMatrixtN=normrnd(0,1,Ns,Nt); % setting up cavity mode coupling matrix
            CoupleMatrixrN=normrnd(0,1,Ns,Nr); % setting up cavity mode coupling matrix

            % setting up cavity Green's function
            green=(j/pi)./(Elevel+(j*alpha));
            greend = diag(green);

            % setting up cavity normalized impedances
            XI_tt = transpose(CoupleMatrixtN) * greend * CoupleMatrixtN;
            XI_tr = transpose(CoupleMatrixtN) * greend * CoupleMatrixrN;
            XI_rt = transpose(CoupleMatrixrN) * greend * CoupleMatrixtN;
            XI_rr = transpose(CoupleMatrixrN) * greend * CoupleMatrixrN;

            % setting up antenna array impedance matrices
            Z_tt(p,:,:) = j .* imag(Z_rad_t) + sqrtm(real(Z_rad_t)) * XI_tt * sqrtm(real(Z_rad_t));
            Z_tr(p,:,:) = sqrtm(real(Z_rad_t)) * XI_tr * sqrtm(real(Z_rad_r));
            Z_rt(p,:,:) = sqrtm(real(Z_rad_r)) * XI_rt * sqrtm(real(Z_rad_t));
            Z_rr(p,:,:) = j .* imag(Z_rad_r) + sqrtm(real(Z_rad_r)) * XI_rr * sqrtm(real(Z_rad_r));


    end

    fname=['znorm_mn_' num2str(Nt) 'a' num2str(Nr) 'p_' num2str(floor(alpha)) 'd' num2str(10*mod(alpha,1)) 'hs' num2str(Ns) '.mat'];
    save(fname, 'Z_tt', 'Z_tr','Z_rt', 'Z_rr', 'alpha','-v7.3');
    clc;
    clear all;
 end