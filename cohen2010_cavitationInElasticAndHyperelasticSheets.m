%-------------------------------------------------------------
% Developed by: 
% Shivam Agarwal (ME14B143)
% Dual Degree Student
% Department of Mechanical Engineering
% Indian Institute of Technology Madras
% Chennai, India, 600036
% 12-7-2019
%-------------------------------------------------------------

% For nomenclature, please read the article:
% Cavitation in elastic and hyperelastic sheets
% Tal Cohen and David Durban
% International Journal of Engineering Science
% Volume 48, Issue 1, January 2010, Pages 52-66

%inputs:
n = 30000;        %number of data points on S axis
m=60;             %number of Einf vals
Smin = 0;         %minimum value of S on S-axis
Smax = 6;         %maximum value of S on S-axis
Einfmax = 1.2;    %maximum value of remote stress
nu=[0,0.1,0.2,0.33,0.5];    %values of poisson's ratio
%--------------------------------------------------------------

Einf=linspace(0,Einfmax,m)';
o = length(nu);
Sinc = (Smax-Smin)/(n-1);
S = linspace(Smin,Smax,n)';

Er = zeros(n,m,o);
Erd = Er;
for i=1:o
    Er(1,:,i)= Einf(:);      %boundary condition
end
Erd(1,:,:)= -1/2;            %limit S->0 applied on differential equation
for k=1:o
    for j=1:m
        for i=2:n
            Er(i,j,k)=Er(i-1,j,k)+ Erd(i-1,j,k)*Sinc;
            Erd(i,j,k)=-((S(i)/(exp((1+nu(k))*S(i))-1))-nu(k)*Er(i,j,k))/...
                ((1-nu(k))*(S(i)/(exp((1+nu(k))*S(i))-1))+1-2*nu(k)*Er(i,j,k));
        end
    end
end

%to find where S becomes zero:
Erabs = Er;
minval = zeros(m,o);
minpos = zeros(m,o);
SatErzero = zeros(m,o)+7;

for k=1:o
    for j=1:m
        Erabs(:,j,k) = abs(Er(:,j,k));
        [minval(j,k),minpos(j,k)] = min(Erabs(:,j,k));
        if minval(j,k)<1e-3
            SatErzero(j,k) = S(minpos(j,k));
        end
    end
end

% % plotting
% hold on;
% for j=1:m
%     plot(S,Er(:,j,1),'linewidth',1.5);
% end
% plot([0,6],[0,0],'linewidth',2);
% hold off;
% axis([0 6 -2 2]);


% plotting
figure
hold on
for k=1:o
    plot(SatErzero(:,k),Einf(:),'linewidth',1.5);
end
legend('\nu=0','\nu=0.1','\nu=0.2','\nu=0.33','\nu=0.5')
xlabel('\epsilon_{\theta}|_{r=a_{o}}','fontweight','bold')
ylabel('\Sigma_{\infty}','fontweight','bold')
axis([0 6 0 1.2])
hold off
%-----------------------------------------------------------