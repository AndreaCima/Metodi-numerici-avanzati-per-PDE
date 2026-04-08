% MNAPDEdrawwaves.m     Andrea Moiola, Pavia, MNAPDE 2019-2026
%
% Draw several Helmholtz solutions in 2D:
% draw real and imaginary part, magnitude,
% gif animation of time-harmonic evolution.
%
% If Matlab opens empty figures and freezes, try typing 
% opengl('save','software') at the command prompt and restart Matlab.


function MNAPDEdrawwaves
%% Set all parameters here
PlotLim = 1;               % Plot functions on square max{|x|,|y|}<PlotLim
k = 30;                    % Wavenumber (>0)
np = 300;                  % Number of points in x & y directions for field plots, affects speed and size of png/gif

anglePW = pi/6;            % Angle defining direction of plane wave
EvanescentParam = 0.3;     % Parameter for definition of evanescent wave
CWindex = 3;               % Index of Fourier-Bessel smooth circular wave
OutCWindex = 5;            % Index of Fourier-Hankel outgoing circular wave
OutCWtruncation = 2;       % Truncate field at this magnitude to improve plot (keep ~2)
nquad = 10;                % Number of Gauss quadrature points in angular coordinate for Herglotz wave
th0 = 0;                   % Herglotz density is = 1 on interval [th0,th1] and 0 outside
th1 = pi/6;
MieRadius = 1/4;           % Radius of disc for Mie scattering

WithSaveFig = 0;           % Save all figures as .png, choose 1 or 0
WithGif = 0;               % Generates time-harmonic animated gifs
% Warning: gif-making is slow, closes all pics, don't use PC while running it!

%% Set up  ----------------------------------------------------------
tic;
close all;
x = linspace(-PlotLim,PlotLim,np);
[xx,yy] = meshgrid(x,x);                % Matrices of x and y coordinates where to evaluate and plot fields
zz = xx + 1i * yy;                      % Treat as complex points
TitleString = sprintf('on $(%g,%g)^2$, $k=%g$', -PlotLim,PlotLim, k);

%% Plane wave
PW_Val = exp(1i*k*(xx*cos(anglePW)+yy*sin(anglePW)));
FH_PW = MyFieldPlot(xx,yy,PW_Val, ['Travelling plane wave $\mathrm e^{\mathrm ik\mathbf x\cdot\mathbf d}$ ', TitleString, ', angle = ',num2str(anglePW)]); 

%% Stationary wave, sum of two plane waves
SW_Val = cos(k*(xx*cos(anglePW)+yy*sin(anglePW)));
FH_SW = MyFieldPlot(xx,yy,SW_Val, ['Stationary plane wave $\cos(k\mathbf x\cdot\mathbf d)=(\mathrm e^{\mathrm ik\mathbf x\cdot\mathbf d}+\mathrm e^{-\mathrm ik\mathbf x\cdot\mathbf d})/2$ ',...
                                   TitleString, ', angle = ',num2str(anglePW)]);

%% Evanescent wave
EW_Val = exp(1i*sqrt(EvanescentParam^2+1)*k*xx-EvanescentParam*k*yy);
FH_EW = MyFieldPlot(xx,yy,EW_Val, ['Evanescent plane wave $\mathrm e^{\mathrm i \mathbf x\cdot\mathbf k}$ ',...
                                    TitleString,', $\mathbf k=(',num2str(sqrt(EvanescentParam^2+1)*k),',',num2str(EvanescentParam*k),')$']);

%% Fundamental solution
FS_Val = (1i/4)*besselh(0,1,k*abs(zz));
FH_FS = MyFieldPlot(xx,yy,FS_Val,['Fundamental solution ',TitleString]);

%% Smooth circular wave
CW_Val = besselj(CWindex, k*abs(zz)) .* exp(1i*CWindex*angle(zz)) ;
FH_CW = MyFieldPlot(xx,yy,CW_Val, ['Smooth circular wave $J_\ell(kr)\mathrm e^{\mathrm i\ell\theta}$ (Fourier--Bessel) ',TitleString,', index $\ell=',num2str(CWindex),'$']);

%% Outgoing circular wave
OutCW_Val = besselh(OutCWindex,1,k*abs(zz)) .* exp(1i*OutCWindex*angle(zz));
OutCW_Val(abs(OutCW_Val)>OutCWtruncation) = nan + 1i*nan;   % Set to nan where magnitude is >OutCWtruncation
FH_OutCW = MyFieldPlot(xx,yy,OutCW_Val, ['Outgoing circular wave $H^{(1)}_\ell(kr)\mathrm e^{\mathrm i\ell\theta}$ (Fourier--Hankel) ',...
                                         TitleString,', index $\ell=',num2str(OutCWindex),'$, truncated at $|u|<',num2str(OutCWtruncation),'$']);

%% Sum of two outgoing circular waves with opposite index
OutCW2_Val = besselh(OutCWindex,1,k*abs(zz)) .* exp(1i*OutCWindex*angle(zz)) ...
             + besselh(-OutCWindex,1,k*abs(zz)) .* exp(-1i*OutCWindex*angle(zz));
OutCW2_Val(abs(OutCW2_Val)>OutCWtruncation) = nan+1i*nan;    % Set to nan where magnitude is >OutCWtruncation
FH_OutCW2 = MyFieldPlot(xx,yy,OutCW2_Val, ['Sum of 2 opposite-$\ell$ outgoing circular waves $H^{(1)}_{\pm\ell}(kr)\mathrm e^{\pm\mathrm i\ell\theta}$ (Fourier--Hankel) ',...
                                           TitleString,', index $\ell=',num2str(OutCWindex),'$, truncated at $|u|<',num2str(OutCWtruncation),'$']);

%% Sum of two outgoing circular waves with different centres
x0 = .5-.5i;
OutCWindexB = 0;
Out2CW_Val = besselh(OutCWindex,1,k*abs(zz)) .* exp(1i*OutCWindex*angle(zz)) ...
             + besselh(OutCWindexB,1,k*abs(zz-x0)) .* exp(1i*OutCWindexB*angle(zz-x0));
Out2CW_Val(abs(Out2CW_Val)>OutCWtruncation) = nan+1i*nan;    % Set to nan where magnitude is >OutCWtruncation
FH_Out2C = MyFieldPlot(xx,yy,Out2CW_Val, ['Sum of 2 circular waves $H^{(1)}_\ell(kr)\mathrm e^{\mathrm i\ell\theta}+H^{(1)}_0(k|\mathbf x-\mathbf x_0|)$ ',TitleString, ...
                                           ', centre $\mathbf x_0=(',num2str(real(x0)),',',num2str(imag(x0)),')$, index $\ell=',num2str(OutCWindex),'$, truncated at $|u|<',num2str(OutCWtruncation),'$']);

%% Herglotz function with piecewise-constant kernel
[xth,wth] = gaussquad(nquad);  % Gauss quadrature in angular coordinate
xth = th0 + xth*(th1-th0);
wth = wth * ( th1-th0 );
Hergl_Val = zeros(size(xx));
for QuadratureIndex = 1:nquad
    Hergl_Val = Hergl_Val + wth(QuadratureIndex)* exp(1i*k*(xx*cos(xth(QuadratureIndex))+yy*sin(xth(QuadratureIndex))));
end
FH_Hergl = MyFieldPlot(xx,yy,Hergl_Val, sprintf('Herglotz function %s, constant kernel supported on $(%.3g,%.3g)$',TitleString,th0,th1));

%% Phased array: equispaced point sources on segment (-0.5,0.5)x{-0.5} emitting with phase shift/delay
PAnumSources = 50;
PAsources = linspace(-0.5,0.5,PAnumSources+1) ; 
PAdelay = 0.3;
PA_Val = zeros(np);
for j = 0:PAnumSources
    PA_Val = PA_Val + exp(1i*j*PAdelay)*besselh(0,1,k*abs(zz-PAsources(j+1)+0.5i)) ;
end
FH_PA = MyFieldPlot(xx,yy,PA_Val, ['Phased array $u(\mathbf x)=\sum_{j=0}^{',num2str(PAnumSources),...
        '}\mathrm e^{\mathrm i',num2str(PAdelay),'j} H^{(1)}_0(k|\mathbf x-\mathbf x_j|)$, $\mathbf x_j=(-0.5+\frac{j}{',num2str(PAnumSources),'},-0.5)$, ', TitleString]);

%% PW scattering by sound-soft disc, Mie series / multipole expansion
Mie_SS_Val =  zeros(size(xx));
MaxL = 2*k*MieRadius + 10;   % Cautious choice for truncation of Mie series
for lMie = -MaxL:MaxL    
    Mie_SS_Val = Mie_SS_Val - (1i^lMie*besselj(lMie,k*MieRadius)*exp(-1i*lMie*anglePW)/besselh(lMie,1,k*MieRadius)) ...
        * besselh(lMie,1,k*abs(zz)) .* exp(1i*lMie*angle(zz));
end
Mie_SS_Val(abs(zz) < MieRadius) = nan + 1i*nan; % chop the field in the scatterer to NaN
Mie_SStot_Val = Mie_SS_Val + PW_Val ;
ScatteringString = [' disc with Mie series ',TitleString,', incoming PW angle $',num2str(anglePW/pi),'\pi$, radius ',num2str(MieRadius)];
FH_SSscat = MyFieldPlot(xx,yy,Mie_SS_Val,  strcat('Field scattered by sound-soft (Dirichlet)', ScatteringString));
FH_SStot = MyFieldPlot(xx,yy,Mie_SStot_Val,  strcat('Total field for scattering by sound-soft (Dirichlet)', ScatteringString));

%% PW scattering by sound-hard disc, Mie series / multipole expansion
Mie_SH_Val =  zeros(size(xx));
for lMie = -MaxL:MaxL    
    Mie_SH_Val = Mie_SH_Val - ...
        (  1i^lMie*(besselj(lMie-1,k*MieRadius)-besselj(lMie+1,k*MieRadius))*exp(-1i*lMie*anglePW)...
        /(besselh(lMie-1,1,k*MieRadius)-besselh(lMie+1,1,k*MieRadius))  ) ...
        * besselh(lMie,1,k*abs(zz)) .* exp(1i*lMie*angle(zz));
end
Mie_SH_Val(abs(zz) < MieRadius) = nan + 1i*nan;
Mie_SHtot_Val = Mie_SH_Val + PW_Val;
FH_SHscat = MyFieldPlot(xx,yy,Mie_SH_Val,  strcat('Field scattered by sound-hard (Neumann)', ScatteringString));
FH_SHtot = MyFieldPlot(xx,yy,Mie_SHtot_Val,  strcat('Total field for scattering by sound-hard (Neumann)', ScatteringString));

%% Scattering of fundamental solution by sound-soft disc, Mie series
FScentre = -1.2-.7i;  % Center of the "incoming" fundamental solution as complex number, should be out of plot box for decent figures
Mie_FS_Val =  zeros(size(xx));
for lMie = -MaxL:MaxL    
    Mie_FS_Val = Mie_FS_Val - ...
        (1i/4)* besselj(lMie,k*MieRadius) * besselh(lMie,1,k*abs(FScentre)) * exp(-1i*lMie*angle(FScentre))  ...
        / besselh(lMie,1,k*MieRadius) * besselh(lMie,1,k*abs(zz)) .* exp(1i*lMie*angle(zz));
end
Mie_FS_Val(abs(zz) < MieRadius) = nan + 1i*nan;    % chop field in the scatterer to 0
Mie_FStot_Val = Mie_FS_Val + (1i/4)*besselh(0,1,k*abs(zz-FScentre));
Mie_FStot_Val = Mie_FStot_Val .* (abs(zz) > MieRadius) ;
Mie_FSstring = ['by sound-soft (Dirichlet) disc with Mie series ', TitleString, ', point source centred at $(',...
                num2str(real(FScentre)),',',num2str(imag(FScentre)),')$, radius ',num2str(MieRadius)];
FH_FSscat = MyFieldPlot(xx,yy,Mie_FS_Val,  ['Field scattered ',Mie_FSstring]);
FH_FStot = MyFieldPlot(xx,yy,Mie_FStot_Val,  ['Total field for scattering ',Mie_FSstring]);

%% Dirichlet eigenfunction in unit disc
IndexDisc = 9; kDisc = 13.354300477435331;
EigenDisc_Val = besselj(IndexDisc,kDisc*abs(zz)) .* exp(1i*IndexDisc*angle(zz)) .* (abs(zz)<1);
EigenDisc_Val(abs(zz)>1) = nan+1i*nan;
FH_EigenDisc = MyFieldPlot(xx,yy,EigenDisc_Val, ['Unit disc eigenfunction, $k=',num2str(kDisc),'$: $u(\mathbf x)=J_{',num2str(IndexDisc),'}(kr)\mathrm e^{\mathrm i',num2str(IndexDisc),'\theta}$']);

%% Random eigenfunction of the square
EigenIndex = ceil(rand(1,2)*10);
k_Eigenv = .5*pi*norm(EigenIndex)/PlotLim;
Eigen_Val = sin((.5*pi*EigenIndex(1)/PlotLim)*(xx+PlotLim)).* sin((.5*pi*EigenIndex(2)/PlotLim)*(yy+PlotLim));
FH_Eigen = MyFieldPlot(xx,yy,Eigen_Val,  sprintf('Dirichlet eigenfunction on $(-%g,%g)^2$, index (%g,%g), $k=%g$',...
    PlotLim,PlotLim,EigenIndex(1),EigenIndex(2),k_Eigenv));

%% Downward plane wave reflected by horizontal line with Dirichlet/Neumann/Impedance BCs
d1 = abs(cos(anglePW));
d2 = abs(sin(anglePW));
IncomingVal = exp(1i*k*(xx*d1-yy*d2));
SnellReflVal = exp(1i*k*(xx*d1+yy*d2));
ReflDir_Val = IncomingVal - SnellReflVal ;
ReflNeu_Val = IncomingVal + SnellReflVal ;
% Impedance BCon {x_2=const}:  du/dn-ikZu = -du/dx_2-ikZu = 0  
ImpParamZ = 1;
%ImpParamZ = d2;   % For this impedance the incoming wave is completely absorbed
ReflImp_Val = IncomingVal + (d2-ImpParamZ)/(d2+ImpParamZ)* SnellReflVal;
%FH_ReflDir = MyFieldPlot(xx,yy,ReflDir_Val,sprintf('Plane wave reflected by $(x_2=-%g)$, %s, Dirichlet BCs',PlotLim,TitleString));
FH_ReflDir = MyFieldPlot(xx,yy,ReflDir_Val,['Plane wave reflected by $\{x_2=',num2str(-PlotLim),'\}$, ',TitleString,', Dirichlet BCs']); 
FH_ReflNeu = MyFieldPlot(xx,yy,ReflNeu_Val,['Plane wave reflected by $\{x_2=',num2str(-PlotLim),'\}$, ',TitleString,', Neumann BCs']);
FH_ReflImp = MyFieldPlot(xx,yy,ReflImp_Val,['Plane wave reflected by $\{x_2=',num2str(-PlotLim),'\}$, ',TitleString,', impedance BCs, $\vartheta=',num2str(ImpParamZ),'$']);

%% Reflection of Herglotz function on horizontal line
th0R = -pi/3;                   % Choose downward-propagating Herglotz function (-pi<th0R<th1R<0)
th1R = -pi/6;   
yyR = yy+1;                     % Translate coordinates upwards
[xthR,wthR] = gaussquad(nquad);  
xthR = th0R + xthR*(th1R-th0R);
wthR = wthR * (th1R-th0R);
HerglR_Val = zeros(size(xx));
for QuadratureIndex = 1:nquad
    HerglR_Val = HerglR_Val + wthR(QuadratureIndex)* exp(1i*k*(xx*cos(xthR(QuadratureIndex))+yyR*sin(xthR(QuadratureIndex))));
end
% Use uI(-x,-y)=conj{uI(x,y)} so uR(x,y):=uI(x,-y)=conj{uI(-x,y)} and
% symmetry around x=0 of domain:
HerglR_Val = HerglR_Val - conj(fliplr(HerglR_Val));
FH_HerglR = MyFieldPlot(xx,yyR,HerglR_Val, ['Reflection of Herglotz function with constant kernel supported on $(',...
                        num2str(th0R),',',num2str(th1R),')$, on sound-soft $\{x_2=0\}$, on $(-1,1)\times(0,2)$, $k=',num2str(k),'$']);

%% Transmission problem: downward plane wave hitting penetrable half space
d1 = cos(anglePW);
d2 = -abs(sin(anglePW));        % make sure PW points downward
kbelow = k*2;                   % wavenumber in lower half plane
% kbelow = k*0.8;               %    kbelow < k can give total internal reflection
d1below = k/kbelow *d1;         % transmitted PW direction vector
d2below = -sqrt(1-d1below^2);
A_transmission = 1;             % normal derivative jump parameter in the transmission conditions
R_transmission = (k*d2-A_transmission*kbelow*d2below)/(k*d2+A_transmission*kbelow*d2below);     % Reflection coefficient R
T_transmission = 1 + R_transmission;                                                            % Transmission coefficient T
uFlatTransmissionVal = ( exp(1i*k*(d1*xx-d2*yy)) +  R_transmission* exp(1i*k*(d1*xx+d2*yy)) ) .*(yy>0)...
                       + T_transmission .* (exp(1i*kbelow*(d1below*xx+d2below*yy)))  .*(yy<=0); % Evaluate solution
FH_FlatTransm = MyFieldPlot(xx,yy,uFlatTransmissionVal,...
    ['Transmission problem: $k^+=',num2str(k),', k^-=',num2str(kbelow),', A=',num2str(A_transmission),...
     ', \mathbf d^+=(',num2str(d1),',',num2str(d2),'), \mathbf d^-=(',num2str(d1below),',',num2str(d2below),')$, on $(-',num2str(PlotLim),',',num2str(PlotLim),')^2$']);  
subplot(1,3,1), hold on  % draw arrows
quiver(0,0,d1/2,-d2/2,0,'r','linewidth',2,'maxheadsize',.7)
quiver(-d1/2,-d2/2,d1/2,d2/2,0,'r','linewidth',2,'maxheadsize',.7)
if isreal(d2below), quiver(0,0,d1below/2,d2below/2,0,'r','linewidth',2,'maxheadsize',.7), end

%% Save figure and make GIF animations, if requested:
if WithSaveFig || WithGif
    dirpath = 'Figs';            % Name of folder where to save figures
    % if absent, make folder where to save figures and gifs
    if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end
end
if WithSaveFig
    MyPrintFig = @(FigHandle,FigName)  exportgraphics(FigHandle,[dirpath,'/FieldPlot',FigName,'.png']);
    %MyPrintFig = @(FigHandle,FigName)  exportgraphics(FigHandle,[dirpath,'/FieldPlot',FigName,'.pdf']);
    %MyPrintFig = @(FigHandle,FigName)  print(FigHandle,'-dpng',[dirpath,'/FieldPlot',FigName,'.png']);
    MyPrintFig(FH_PW,'PW');
    MyPrintFig(FH_SW,'SW');
    MyPrintFig(FH_EW,'EW');    
    MyPrintFig(FH_FS,'FS');
    MyPrintFig(FH_CW,'CW');
    MyPrintFig(FH_OutCW,'OutCW');
    MyPrintFig(FH_OutCW2,'OutCW2');
    MyPrintFig(FH_Out2C,'Out2C');
    MyPrintFig(FH_Hergl,'Hergl');
    MyPrintFig(FH_PA,'PhasedArray');
    MyPrintFig(FH_SSscat,'SoundSoftDiscScat');
    MyPrintFig(FH_SStot,'SoundSoftDiscTot');    
    MyPrintFig(FH_SHscat,'SoundHardDiscScat');
    MyPrintFig(FH_SHtot,'SoundHardDiscTot');
    MyPrintFig(FH_FSscat,'FSDiscScat');    
    MyPrintFig(FH_FStot,'FSDiscTot');
    MyPrintFig(FH_EigenDisc,'EigenfunctionDisc');
    MyPrintFig(FH_Eigen,'Eigenfunction');
    MyPrintFig(FH_HerglR,'HerglRefl');
    MyPrintFig(FH_FlatTransm,['FlatTransmissionK',num2str(kbelow)]);
    disp('Figures saved.')
    close all
end
if WithGif 
    close all  % Close all figure otherwise sometimes Matlab messes up
    MyGif = @(Val,Str) MyGiffer(xx,yy,Val,strcat(dirpath,'/AnimFieldPlot',Str,'.gif'));
    MyGif(PW_Val,'PW');
    MyGif(SW_Val,'SW');
    MyGif(EW_Val,'EW');
    MyGif(FS_Val,'FS');    
    MyGif(CW_Val,'CW');
    MyGif(OutCW_Val,'OutCW');    
    MyGif(OutCW2_Val,'OutCW2');  
    MyGif(Out2CW_Val,'Out2C');
    MyGif(Hergl_Val,'Hergl');
    MyGif(PA_Val,'PhasedArray');
    MyGif(Mie_SS_Val,'SoundSoftDiscScat');
    MyGif(Mie_SStot_Val,'SoundSoftDiscTot');
    MyGif(Mie_SH_Val,'SoundHardDiscScat');
    MyGif(Mie_SHtot_Val,'SoundHardDiscTot');
    MyGif(Mie_FS_Val,'FSDiscScat');
    MyGif(Mie_FStot_Val,'FSDiscTot');
    MyGif(EigenDisc_Val,'EigenfunctionDisc'); 
    MyGif(Eigen_Val,'Eigenfunction');
    MyGif(ReflDir_Val,'ReflDir');
    MyGif(ReflNeu_Val,'ReflNeu');
    MyGif(ReflImp_Val,'ReflImp'); 
    MyGif(HerglR_Val,'HerglRefl');  
    MyGif(uFlatTransmissionVal,['FlatTransmissionK',num2str(kbelow)]);  
    disp('Gif animations saved.')
end
toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions
function FH = MyFieldPlot(xx,yy,Field,PlotTitle)
% Plot real part, imaginary part and magnitude of Field in points xx/yy
% in three subplots of a figure, superimposing PlotTitle string as title.
% Returns figure handle FH.

FH = figure;
set(gcf,'position',[50, 100, 1400,400]);

subplot(1,3,1)
pcolor(xx,yy,real(Field));
axis off;   axis equal; shading flat;   colorbar;
title('Real part')

subplot(1,3,2)
pcolor(xx,yy,imag(Field));
axis off;   axis equal; shading flat;   colorbar;
title('Imaginary part')

sp3 = subplot(1,3,3);
pcolor(xx,yy,abs(Field)); title('Absolute value');
%pcolor(xx,yy,log10(abs(Field))); title('log10(Absolute value)');
%pcolor(xx,yy,abs(Field).^2); title('(Absolute value)^2');
%pcolor(xx,yy,angle(Field)); title('Complex argument');
axis off;   axis equal; shading flat;
colormap(sp3,hot);
colorbar;
sgtitle(PlotTitle,'interpreter','latex'); % add main title
addToolbarExplorationButtons(gcf);  % Matlab old style zoom buttons   
end

%% ------------------------------------------------------------------
function MyGiffer(xx,yy,M,GifName)
% Generate time-harmonic animation of complex field M on points xx/yy, 
% export as gif in file GifName --- not very reliable
% Do not touch mouse & keyboard while running!
figure
pcolor(xx,yy,(real(M)));
colormax = max(max(abs(M)));
clim([-colormax,colormax]);
axis tight;  axis equal;  shading interp;
set(gca,'nextplot','replacechildren','visible','off')
f = getframe; %(gcf); % gcf sometimes avoids errors but gives ugly gray box
[im,map] = rgb2ind(f.cdata,256,'nodither');
for t_count = 1:25  % loop over 25 frames
    t = t_count*2*pi/25;
    pcolor((real(exp(-t*1i)*M)));
    axis tight;    axis square;   shading interp;
    f = getframe; %(gcf);
    im(:,:,1,t_count) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,GifName,'gif','DelayTime',0,'LoopCount',inf);
hold off;
close;
end     

%% ------------------------------------------------------------------
function [x,w] = gaussquad(q)
% Computes weights and nodes for Gaussian quadrature on [0,1],
% exact for polynomials up to degree 2p-1
%
% Input:  p : number of quadrature nodes
% Output: x : nodes in [0,1] (column vector)
%         w : weights        (column vector)

b = (1:(q-1)) ./ sqrt(4*(1:(q-1)).^2-1) ;
[ev,ew] = eig( diag(b,-1) + diag(b,1) );
x = (diag(ew)+1)/2;
w = (ev(1,:).*ev(1,:))';
end
