printc='bw';

%if usejava('jvm') || MATLABVERSION <= 7.12
%% Garbage collection and initialization
format compact %remove blank lines from output
format long e %lots of digits
%clearvars %clear all variables
close all %close all figures
ColorOrder=get(gca,'ColorOrder'); close all
set(0,'defaultaxesfontsize',21,'defaulttextfontsize',21) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%% Initialize parameters
mmax=22;
mplot=12;
mlag=4; % the coefficients to shrink wrt mplot

testfun=@(x) exp(-3*x).*sin(3*pi*x.^2);

r=1;
lambda_i=@(i,r) lambda_k(floor(i/2),r);
pos_v=@(k,x,r) sqrt(2)*sin(2*pi*bsxfun(@times,k,x));
neg_v=@(k,x,r) sqrt(2)*cos(2*pi*bsxfun(@times,k,x));
v_i=@(i,x,r) bsxfun(@times,mod(i,2)==0,neg_v(-i/2,x,r))+bsxfun(@times,and(mod(i,2)==1,i>1),pos_v((i-1)/2,x,r))+(i==1);
u_i=@(i,x,r) bsxfun(@times,lambda_i(i,r),v_i(i,x,r));

i_max=(0:2^mmax-1)';
y = testfun(i_max/2^mmax);
yfft_ = fft(y);
yfft = yfft_;
yfft(2:2:end) = sqrt(2)*real(yfft_(2:2^(mmax-1)+1))./lambda_i(2:2:size(yfft,1),r)';
yfft(3:2:end) = -sqrt(2)*imag(yfft_(2:2^(mmax-1)))./lambda_i(3:2:size(yfft,1),r)';

yticks=10.^(-10:5:10);
ymin = 10^-11;
ymax = 10^8;


%% Plot function
figure
xplot=(0:0.002:1);
yplot=testfun(xplot);
if strcmp(printc,'color')
    plot(xplot,yplot,'-');
else
    plot(xplot,yplot,'-k');
end
yminf=1.1*min(yplot);
ymaxf=1.1*max(yplot);
axis([0 1 yminf ymaxf])
xlabel('\(x\)','interpreter','latex')
ylabel('\(f(x)\)','interpreter','latex')
set(gca,'Xtick',[0 0.5 1.0])
print('Function','-depsc');

%% Plot scaled function
figure
scale = 100.0;
if strcmp(printc,'color')
    plot(xplot,scale * yplot,'-');
else
    plot(xplot,scale * yplot,'-k');
end
axis([0 1 1.1*min(scale * yplot) 1.1*max(scale * yplot)])
xlabel('\(x\)','interpreter','latex')
ylabel('\(f_{\textnormal{big}}(x)\)','interpreter','latex')
set(gca,'Xtick',[0 0.5 1.0])
print('Scaled','-depsc');


%% Plot FW coefficients
nplot=2^mplot;
i_plot=1:nplot;
ymap=yfft;
yfwtabs=abs(ymap(i_plot));
figure
switch printc
    case 'color'
        h=loglog(i_plot, yfwtabs,'.','MarkerSize',20);
    case 'bw'
        h=loglog(i_plot, yfwtabs,'.k','MarkerSize',20);
end
maxexp=floor(log10(nplot-1));
set(gca,'Xtick',10.^(0:maxexp))
set(gca,'Ytick',yticks)
axis([1 nplot-1 ymin ymax])
xlabel('\(i\)','interpreter','latex')
ylabel('\(|\widehat{f}_{i}|\)','interpreter','latex')
set(gca,'Position',[0.2 0.155 0.75 0.77])
print('FunctionWalshFourierCoeffDecay.eps', '-depsc');


%% Plot scaled function Fourier coefficients
figure
switch printc
    case 'color'
        h=loglog(i_plot, scale*yfwtabs,'.','MarkerSize',20);
    case 'bw'
        h=loglog(i_plot, scale*yfwtabs,'.k','MarkerSize',20);
end
maxexp=floor(log10(nplot-1));
set(gca,'Xtick',10.^(0:maxexp))
set(gca,'Ytick',yticks)
axis([1 nplot-1 ymin ymax])
xlabel('\(i\)','interpreter','latex')
ylabel('\(|\widehat{f_{\textnormal{big}}}_{i}|\)','interpreter','latex')
set(gca,'Position',[0.2 0.155 0.75 0.77])
print('ScaledWalshFourierCoeffDecay.eps', '-depsc');


%% Plot filtered functions Fourier coefficients
whkill=(2^(mplot-mlag)+1):2^(mplot-mlag+1); %choose ones to shrink
ymap(whkill)=1e-10*ymap(whkill); %shrink them
yfwtabs=abs(ymap);
figure
loglog(i_plot,yfwtabs(i_plot),'k.','MarkerSize',20);
hold on
loglog(whkill,yfwtabs(whkill),'.','MarkerSize',20);
hold off
maxexp=floor(log10(nplot-1));
set(gca,'Xtick',10.^(0:maxexp))
set(gca,'Ytick',yticks)
axis([1 nplot-1 ymin ymax])
xlabel('\(i\)','interpreter','latex')
ylabel('\(|\widehat{f_{\textnormal{fuzzy}}}_{i}|\)','interpreter','latex')
set(gca,'Position',[0.2 0.155 0.75 0.77])
print('FilteredWalshFourierCoeffDecay.eps', '-depsc');

%% Plot filtered function
% ------ Old usual basis method, two lines below
% hat_f=ymap;
% yplot = sum(bsxfun(@times,hat_f,u_i(i_max+1,xplot,r)),1) / 2^mmax;
% ------ Inverse FFT method, 6 lines below
yfft_ = ymap;
yfft_real = ymap(2:2:end).*lambda_i(2:2:size(ymap,1),r)';
yfft_imag = ymap(3:2:end).*lambda_i(3:2:size(ymap,1),r)';
yfft_(2:end/2+1) = (yfft_real - [sqrt(-1)*yfft_imag;0])/sqrt(2);
yfft_(end:-1:end/2+2) = (yfft_real(1:end-1) + sqrt(-1)*yfft_imag)/sqrt(2);
yplot = ifft(yfft_);

figure
if strcmp(printc,'color')
    plot(i_max/2^mmax,yplot,'-');
else
    plot(i_max/2^mmax,yplot,'-k');
end
axis([0 1 yminf ymaxf])
xlabel('\(x\)','interpreter','latex')
ylabel('\(f_{\textnormal{fuzzy}}\)','interpreter','latex')
set(gca,'Xtick',[0 0.5 1.0])
print('Filtered','-depsc');



