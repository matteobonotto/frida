function [psibar,FdF,dP] = fun_ArmonicSmoothing(psibar,FdF,dP,mm,m_tresh)

weigths = [ones(1,m_tresh) linspace(1,0,mm(end)-m_tresh)];
weigths = [fliplr(weigths) 1 weigths]'; 

theta = pi*[-fliplr(psibar(2:end)) psibar]';

%% smoothing of FFpsime
fun = [fliplr(FdF(2:end)) FdF]';


% % for ii = 25:2:80
% %     mm = -ii:ii;
% %     
% %     m_tresh = 10;
% %     weigths = [ones(1,m_tresh) linspace(1,0,mm(end)-m_tresh)];
% %     weigths = [fliplr(weigths) 1 weigths]';
% %     
% %     [c_n]=fun_FFT(fun,mm,theta');
% %     TT = exp(1i*theta*mm);
% %     figure
% %     plot(theta,fun,'k'); hold on; xlim([-pi,pi]); plot(theta,real(TT*(c_n)),'r'); plot(theta,real(TT*(c_n.*weigths)),'b--'); title(ii)
% % 
% %     % % figure; semilogy(mm,abs(real(c_n)))
% %     pause
% % end


[c_n]=fun_FFT(fun,mm,theta');
TT = exp(1i*theta*mm);
% % figure
% % plot(theta,fun); hold on; xlim([-pi,pi]); plot(theta,real(TT*(c_n))); plot(theta,real(TT*(c_n.*weigths))); title(ii)
% % figure; semilogy(mm,abs(real(c_n)),'k'); hold on; semilogy(mm,abs(real(c_n.*weigths)));
fun_smooth = real(TT*(c_n.*weigths));
FdF = fun_smooth(numel(psibar):end)';
FdF(end) = 0;

%% smoothing of Ppsime
% % mm = -35:35;
fun = [fliplr(dP(2:end)) dP]';
[c_n]=fun_FFT(fun,mm,theta');
TT = exp(1i*theta*mm);
% % figure
% % plot(theta,fun); hold on; xlim([-pi,pi]); plot(theta,real(TT*c_n))
% % figure; semilogy(mm,abs(real(c_n)))
fun_smooth = real(TT*(c_n.*weigths));
dP = fun_smooth(numel(psibar):end)';
dP(end) = 0;

%% Double sampling points (with pchip interpolation)
psibar_int = unique([psibar .5*(psibar(2:end)+psibar(1:end-1))]);

FdF_int = interp1(psibar,FdF,psibar_int,'PCHIP');
dP_int = interp1(psibar,dP,psibar_int,'PCHIP');


psibar = psibar_int;
FdF = FdF_int;
dP = dP_int;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c_n]=fun_FFT(fun,mm,x)

%% BEST VERSION
[z,w] = MacGaussQuad1D(10);
x0 = (x(1:end-1)+x(2:end))*0.5;
h2 = diff(x)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

yy    = interp1(x,fun,xx').*wh';
expmt = exp(-1i*xx*mm)/(2*pi);
c_n  = (yy*expmt).';


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = MacGaussQuad1D(n)

%% Compute rule for 1D for interval [-1,1], with the order = n

eps = 2.0e-15;
xstart = -1.0;  xstop = 1.0;

x = zeros(1,n);
w = zeros(1,n);

m  = round((n+1)/2);
xm = 0.5*(xstop+xstart);
xl = 0.5*(xstop-xstart);

if n>1
    for j1=1:m
        z = cos(pi*(j1-0.25)/(n+0.5));
        z1 = 2*z;
        while (abs(z-z1)>eps)
            p1=1.0;
            p2=0.0;
            for j2=1:n
                p3=p2;
                p2=p1;
                p1=((2.0*j2-1.0)*z*p2-(j2-1.0)*p3)/j2;
            end
            pp=n*(z*p1-p2)/(z*z-1);
            z1=z;
            z=z1-p1/pp;
        end
        k=n+1-j1;
        x(j1)=xm-xl*z;
        x(k)=xm+xl*z;
        w(j1)=2*xl/((1-z*z)*pp*pp);
        w(k)=w(j1);
    end
end

if n==1
    x(1) = 0.0;
    w(1) = 2.0;
end

%disp(['   GaussQuad1D: x=[' num2str(x) '], w=[',num2str(w) ']'])
end





