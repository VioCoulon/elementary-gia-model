% geoid model + isostatic adjustment based on spatially-varying ELRA
% Coulon et al. 2021
% violaine.coulon@ulb.be
% september 2023

function Elementary_GIA_Model

clear;
% close all;

%%% INPUT DATA %%%

load Antarctica25km.mat; % Load input data
ctr.imax=225; % Maximum number of grid points in the x-direction
ctr.jmax=225; % Maximum number of grid points in the y-direction
ctr.delta=25e3; % Resolution in meters

%%% EXPERIMENT SET-UP %%%

forcing=3; % Experiment type: 1 - Uniform mass loss, 2 - Uniform mass gain, 3 - WAIS Collapse

ctr.nsteps=50; % Number of time steps
ctr.dt=10; % Time step duration (years)
time=linspace(ctr.dt,ctr.dt*ctr.nsteps,ctr.nsteps)'; %  Time vector

ctr.GeoidCalc=1; % Geoid calculation: 0 - No geoid changes, 1 - Geoid changes
ctr.BedAdj=1; % Bedrock adjustment: 0 - No isostatic adjustment, 1 - ELRA model (uniform or spatially-varying), 2 - Local isostatic adjustment
var_ELRA=1; % ELRA model type: 0 - Uniform ELRA model (Green's function formalism), 1 - Spatially-varying ELRA model
par.FlexRigid=1e25; % Uniform flexural rigidity (if var_ELRA = 0)
par.bedrelax=3000.; % Uniform asthenosphere relaxation time (if var_ELRA = 0)

%%% CONSTANTS AND PARAMETERS %%%

par.g=9.81; % gravitational acceleration
par.rho=917.; % ice density
par.rhow=1027.; % sea water density
par.rhom=3370.; % lithosphere density
par.rhof=1000; % fresh water density
par.nuB=0.25; % Poisson ratio in flexural rigidity
par.Re=6.3781e6; % Radius Earth
par.Me=5.972e24; % Mass Earth
par.Aoc=3.618e14; % ocean surface
par.geoidist=5000e3; % size of convolution filter for geoid changes
par.SLref=0; % reference sea level (Goelzer et al. 2020)

%%% INITIALISATION %%%

H0 = H; % Initial ice thickness
B0 = B; % Initial bedrock elevation
Hn = H; % Ice thickness at the current time step
Bn = B; % Bedrock elevation at the current time step

fc.DeltaSL=zeros(ctr.nsteps,1); % Background sea-level change (here fixed to zero)

% Definition of matrix sizes
VAF = zeros(ctr.nsteps, 1); % Volume above floatation
VA0 = zeros(ctr.nsteps, 1); % Volume change relative to reference sea level
POV = zeros(ctr.nsteps, 1); % Potential ocean volume
DENScorr = zeros(ctr.nsteps, 1); % Density correction for transformation from ice to freshwater
SLC = zeros(ctr.nsteps, 1); % Sea level change

% Node count for sparse matrix for sparse_solver_bedrock

nodes=ctr.imax*ctr.jmax;
node=linspace(1,nodes,nodes)';
node=reshape(node,[ctr.imax ctr.jmax]);
VM=ones(nodes,1);

% Initial Volume Above Floatation

dSLR=zeros(ctr.imax,ctr.jmax); % spatially varying component of sea level (geoid effect)
SLR=dSLR+fc.DeltaSL(1); % total sea level (background + fingerprint)
SLR0=SLR;
VAF0=max(0,H0+min(B0-SLR,0)*(par.rhow/par.rho)); % Initial volume above floatation
POV0=max(0,par.SLref-B0); % Initial potential ocean volume (Goelzer et al. 2020)

% Initialisation of geoid mass

if ctr.GeoidCalc==1
    frg=round(par.geoidist/ctr.delta-0.5);
    Pg0=zeros(ctr.imax,ctr.jmax);
    Pg0(MASK==1)=par.rho*H0(MASK==1)*ctr.delta^2.;
    Pg0(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
    Pg0=Pg0+par.rhom*B0*ctr.delta^2.; % addition of bed change (if applicable)
    [gkx,gky]=meshgrid(-par.geoidist:ctr.delta:par.geoidist,-par.geoidist:ctr.delta:par.geoidist);
    dist=sqrt(gkx.^2+gky.^2);
    geoide=par.Re/par.Me*(1./(2*sin(dist/(2*par.Re))));
    geoidmax=par.Re/par.Me*(1./(2*sin(10/(2*par.Re))));
    geoide(geoide>geoidmax)=0;
    geoide(geoide<1.5e-18)=0;
end

% Initialization of bedrock loads (considering initial bedrock in isostatic
% equilibrium) for constant and variable flexural rigidity

if ctr.BedAdj==1 % ELRA
    if var_ELRA==1
        % Varying flexural rigidity 
        load0=zeros(ctr.imax,ctr.jmax);
        load0(MASK==1)=par.rho*par.g*H(MASK==1);
        load0(MASK==0)=(-par.rhow*par.g)*(B(MASK==0)-SLR(MASK==0));
        B0=B;
    else
        % Constant flexural rigidity
        Ll=(par.FlexRigid/par.rhom/par.g)^0.25;
        frb=round(6*Ll/ctr.delta);
        P=zeros(ctr.imax,ctr.jmax);
        P(MASK==1)=par.rho*par.g*H(MASK==1)*ctr.delta^2.;
        P(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
        P0=NaN(ctr.imax+2*frb,ctr.jmax+2*frb);
        P0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=P;
        P0(isnan(P0))=P(ctr.imax,ctr.jmax);
        [kerx,kery]=meshgrid(-frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll, ...
            -frb*ctr.delta/Ll:ctr.delta/Ll:frb*ctr.delta/Ll);
        ker=sqrt((kerx.^2+kery.^2));
        ker(ker<1e-8)=1e-8;
        kei=imag(besselk(0,ker*(1+1i)/sqrt(2)));
        load0=-xcorr2(P0,kei)*Ll^2./(2*pi*par.FlexRigid);
        load0=load0(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb);
        B0=B;
    end
end
if ctr.BedAdj==2 % Local bedrock adjustment
    load0=zeros(ctr.imax,ctr.jmax);
    load0(MASK==1)=par.rho*H(MASK==1)/par.rhom;
    load0(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))/par.rhom;
    B0=B;
end

%%% Forcing %%%

if forcing==1
    H=H-200; % uniform mass loss
    H(H<0)=0;
elseif forcing==2
    H=H+100; % uniform mass gain
elseif forcing==3 
    H(ZB<2 | ZB>17)=0; % WAIS collapse
    H(H<0)=0;
end
Hn=H;

%%%%%%% ITERATION IN TIME %%%%%%%

for cnt=1:ctr.nsteps
    
    % Copying H and B on old values
    if cnt>1
        H=Hn; % If H evolves through time (i.e. whithin an ice sheet model)
        if ctr.BedAdj>0
            B=Bn; % assigning new bed elevation to previous ones
        end
    end
    
    %------------------------------------------------------
    % Adjusting MASK based on floating condition
    %------------------------------------------------------

    SLR=dSLR+fc.DeltaSL(cnt);
    HAF=B-SLR+H*par.rho/par.rhow;
    MASK(HAF<0)=0;
    MASK(HAF>=0)=1;
    HB=max(SLR-par.rho/par.rhow*H,B);
    sn=H+HB; % ice elevation

    if ctr.BedAdj==1
        % Bedrock loads for isostatic adjustment for constant and 
        % variable flexural rigidity
        if var_ELRA==1
            % variable flexural rigidity 
            loadB=zeros(ctr.imax,ctr.jmax);
            loadB(MASK==1)=par.rho*par.g*Hn(MASK==1)-load0(MASK==1);
            loadB(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))-load0(MASK==0);
            bload=sparse_solver_bedrock(Db,par,ctr,VM,node,nodes,loadB);
        else
            % constant flexural rigidity
            P=zeros(ctr.imax,ctr.jmax);
            P(MASK==1)=par.rho*par.g*Hn(MASK==1)*ctr.delta^2.;
            P(MASK==0)=(-par.rhow)*par.g*(B(MASK==0)-SLR(MASK==0))*ctr.delta^2.;
            P0=NaN(ctr.imax+2*frb,ctr.jmax+2*frb);
            P0(frb+1:frb+ctr.imax,frb+1:frb+ctr.jmax)=P;
            P0(isnan(P0))=P(ctr.imax,ctr.jmax);
            bload=-xcorr2(P0,kei)*Ll^2./(2*pi*par.FlexRigid);
            bload=bload(2*frb+1:ctr.imax+2*frb,2*frb+1:ctr.jmax+2*frb)-load0;
        end
    end
    if ctr.BedAdj==2
        % Local isostasy
        bload=zeros(ctr.imax,ctr.jmax);
        bload(MASK==1)=par.rho*Hn(MASK==1)/par.rhom-load0(MASK==1);
        bload(MASK==0)=(-par.rhow)*(B(MASK==0)-SLR(MASK==0))/par.rhom-load0(MASK==0);
    end

    % Bedrock adjustment
    if ctr.BedAdj>0
        if var_ELRA==1
            Bn=(B-B0+bload)*ctr.dt./(-Btau)+B;
        else
            Bn=(B-B0+bload)*ctr.dt/(-par.bedrelax)+B;
        end
    end

%-----------------------------------------------
% Volume above floatation and fingerprinting
%-----------------------------------------------
    VAFi=max(0,H+min(B-SLR,0)*(par.rhow/par.rho));
    % VAF variation (SL equivalent in ocean water)
    VAF(cnt)=sum(sum((VAF0-VAFi)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    VA0i=max(0,H+min(B-par.SLref,0)*(par.rhow/par.rho));
    % VAreferenceSL variation (SL equivalent in ocean water)
    VA0(cnt)=sum(sum((VAF0-VA0i)*ctr.delta^2.))*par.rho/(par.Aoc*par.rhow);
    POVi=max(0,par.SLref-B);
    % Potential Ocean volume variation (SL equivalent in ocean water)
    POV(cnt)=sum(sum((POV0-POVi)*ctr.delta^2.))/par.Aoc;
    % Density correction for transformation from ice to freshwater and 
    % not ocean water (SL equivalent)
    DENScorr(cnt)=sum(sum((H0-H)*ctr.delta^2.))*((par.rho/par.rhow)-(par.rho/par.rhof))/par.Aoc;
    SLC(cnt)=VA0(cnt)+POV(cnt)-DENScorr(cnt); % SL calculation following Goelzer et al. 2020 (TC)
    
    if ctr.GeoidCalc==1

        % Kori-ULB
        % Calculation of geoid change due to mass changes of the ice sheet. The
        % effect is stored in SLR (spatial variability of global sea level)

        % use new ice thickess and bed elevation
        HAF=Bn-SLR+Hn*par.rho/par.rhow;
        MASK(HAF<0)=0;
        MASK(HAF>=0)=1;
        Pg1=zeros(ctr.imax,ctr.jmax);
        Pg1(MASK==1)=par.rho*Hn(MASK==1)*ctr.delta^2.;
        Pg1(MASK==0)=(-par.rhow)*(Bn(MASK==0)-SLR(MASK==0)+fc.DeltaSL(cnt))*ctr.delta^2.; % Takes into account previous geoid change but not DeltaSL (external SL forcing) as ocean mass
        Pg1=Pg1+par.rhom*Bn*ctr.delta^2.; % Takes into account land deformation
        DPg=Pg1-Pg0;
        % Fingerprint
        loadH=conv2fft(DPg,geoide); % FFT convolution
        MaP=zeros(ctr.imax+2*frg,ctr.jmax+2*frg);
        MaP(frg+1:frg+ctr.imax,frg+1:frg+ctr.jmax)=MASK;
        % Determine local sea-level change
        Vd=sum(loadH(MaP==0))*ctr.delta^2;  % Volume change in domain without mass addition
        loadG=loadH(frg+1:frg+ctr.imax,frg+1:frg+ctr.jmax);
        dSLR=loadG+SLC(cnt)-Vd/par.Aoc; 
    end
    
    fprintf('%d %f %f\n',cnt, mean(B(:)-B0(:)), mean(SLR(:)));

    % Plot the results (customize as needed)
    figure(1);
    subplot(1, 3, 1);
    imagesc(SLR);
    axis xy;
    axis equal;
    axis tight;
    title('Geoid Change (m)');
    colorbar;
    axis equal tight;
    
    subplot(1, 3, 2);
    imagesc(Bn-B0);
    axis xy;
    axis equal;
    axis tight;
    title('Bedrock Elevation Change (m)');
    colorbar;
    axis equal tight;

    subplot(1, 3, 3);
    imagesc(B0-B+SLR);
    axis xy;
    axis equal;
    axis tight;
    title('Relative Sea-Level Change (m)');
    colorbar;
    axis equal tight;
    
    pause(0.01); % Pause to allow time for figure to update (adjust as needed)
end

SLR=SLR0+dSLR+fc.DeltaSL(cnt); % SLR = Initial SLR + Geoid perturbation + External SL forcing

save toto;

fprintf('Simulation completed.\n');

end

%----------------------------------------------------------------------
% Solving thin-plate equation with spatially-varying flexural rigidity
%----------------------------------------------------------------------

function bload=sparse_solver_bedrock(Db,par,ctr,VM,node,nodes,loadB)

% Kori-ULB
% Solving thin-plate equation with spatially-varying flexural rigidity for
% isostatic adjustment

Db1=circshift(Db,[0 -1]); % d(i,j+1)
Db2=circshift(Db,[0  1]); % d(i,j-1)
Db3=circshift(Db,[-1 0]); % d(i+1,j)
Db4=circshift(Db,[1 0]); % d(i-1,j)
Db5=circshift(Db,[-1 -1]); % d(i+1,j+1)
Db6=circshift(Db,[1  -1]); % d(i-1,j+1)
Db7=circshift(Db,[-1 1]); % d(i+1,j-1)
Db8=circshift(Db,[1 1]); % d(i-1,j-1)

MASKb=zeros(ctr.imax,ctr.jmax);
MASKb(1:2,:)=1;
MASKb(ctr.imax-1:ctr.imax,:)=1;
MASKb(:,1:2)=1;
MASKb(:,ctr.jmax-1:ctr.jmax)=1;

nabla2Db = (-4*Db+Db1+Db2+Db3+Db4);
dDbx = (Db1-Db2);
dDby = (Db3-Db4);
dDbx2 = (-2*Db+Db1+Db2);
dDbxy = (Db5-Db6-Db7+Db8)/4;
dDby2 = (-2*Db+Db3+Db4);

V0 = (20*Db-4*nabla2Db-(1-par.nuB)*(-2*dDbx2-2*dDby2))/ctr.delta^4+par.rhom*par.g;
V1 = (-8*Db-dDbx-dDbx+nabla2Db-(1-par.nuB)*dDby2)/ctr.delta^4;
V2 = (-8*Db+dDbx+dDbx+nabla2Db-(1-par.nuB)*dDby2)/ctr.delta^4;
V3 = (-8*Db-dDby-dDby+nabla2Db-(1-par.nuB)*dDbx2)/ctr.delta^4; 
V4 = (-8*Db+dDby+dDby+nabla2Db-(1-par.nuB)*dDbx2)/ctr.delta^4;         
V5 = (2*Db+0.5*dDbx+0.5*dDby-0.5*(1-par.nuB)*(-dDbxy))/ctr.delta^4;
V6 = (2*Db+0.5*dDbx-0.5*dDby-0.5*(1-par.nuB)*(dDbxy))/ctr.delta^4;
V7 = (2*Db-0.5*dDbx+0.5*dDby-0.5*(1-par.nuB)*(dDbxy))/ctr.delta^4;
V8 = (2*Db-0.5*dDbx-0.5*dDby-0.5*(1-par.nuB)*(-dDbxy))/ctr.delta^4;
V9 = (Db+0.5*dDbx)/ctr.delta^4;
V10 = (Db-0.5*dDbx)/ctr.delta^4;
V11 = (Db+0.5*dDby)/ctr.delta^4;
V12 = (Db-0.5*dDby)/ctr.delta^4;

R0 = loadB;

V0(MASKb==1)=par.rhom*par.g;
V1(MASKb==1)=0;
V2(MASKb==1)=0;
V3(MASKb==1)=0;
V4(MASKb==1)=0;
V5(MASKb==1)=0;
V6(MASKb==1)=0;
V7(MASKb==1)=0;
V8(MASKb==1)=0;
V9(MASKb==1)=0;
V9(MASKb==1)=0;
V10(MASKb==1)=0;
V11(MASKb==1)=0;
V12(MASKb==1)=0;
R0(MASKb==1)=loadB(MASKb==1);

V=[reshape(V0(VM==1),nodes,1)
    V1(V1~=0)
    V2(V2~=0)
    V3(V3~=0)
    V4(V4~=0)
    V5(V5~=0)
    V6(V6~=0)
    V7(V7~=0)
    V8(V8~=0)
    V9(V9~=0)
    V10(V10~=0)
    V11(V11~=0)
    V12(V12~=0)
    ];

row=[reshape(node(VM==1),nodes,1)
    node(V1~=0)
    node(V2~=0)
    node(V3~=0)
    node(V4~=0)
    node(V5~=0)
    node(V6~=0)
    node(V7~=0)
    node(V8~=0)
    node(V9~=0)
    node(V10~=0)
    node(V11~=0)
    node(V12~=0)
    ];

nodeV1=circshift(node,[0 -1]); % i,j+1
nodeV2=circshift(node,[0 1]); % i,j-1
nodeV3=circshift(node,[-1 0]); % i+1,j
nodeV4=circshift(node,[1 0]); % i-1,j
nodeV5=circshift(node,[-1 -1]); % i+1,j+1
nodeV6=circshift(node,[1 -1]); % i-1,j+1
nodeV7=circshift(node,[-1 1]); % i+1,j-1
nodeV8=circshift(node,[1 1]); % i-1,j-1
nodeV9=circshift(node,[0 -2]); % i,j+2
nodeV10=circshift(node,[0 2]); % i,j-2
nodeV11=circshift(node,[-2 0]); % i+2,j
nodeV12=circshift(node,[2 0]); % i-2,j

col=[reshape(node(VM==1),nodes,1)
    nodeV1(V1~=0)
    nodeV2(V2~=0)
    nodeV3(V3~=0)
    nodeV4(V4~=0)
    nodeV5(V5~=0)
    nodeV6(V6~=0)
    nodeV7(V7~=0)
    nodeV8(V8~=0)
    nodeV9(V9~=0)
    nodeV10(V10~=0)
    nodeV11(V11~=0)
    nodeV12(V12~=0)
    ];

R=reshape(R0(VM==1),nodes,1);

% construct sparse matrix

A=sparse(row,col,V);

% solve

s=A\R;
bload = zeros(ctr.imax,ctr.jmax);
bload(node>0)=s(node(node>0));

end

function C = conv2fft(varargin)
% C = conv2fft(A, B)
% C = conv2fft(H1, H2, A)
%
%   C = CONV2FFT(A, B) performs the 2-D convolution of matrices A and B.
%   If [ma,na] = size(A), [mb,nb] = size(B), and [mc,nc] = size(C), then
%   mc = max([ma+mb-1,ma,mb]) and nc = max([na+nb-1,na,nb]).
%
%   C = CONV2FFT(H1, H2, A) first convolves each column of A with the vector
%   H1 and then convolves each row of the result with the vector H2.  If
%   n1 = length(H1), n2 = length(H2), and [mc,nc] = size(C) then
%   mc = max([ma+n1-1,ma,n1]) and nc = max([na+n2-1,na,n2]).
%   CONV2(H1, H2, A) is equivalent to CONV2FFT(H1(:)*H2(:).', A) up to
%   round-off.
%
%   C = CONV2FFT(..., SHAPE) returns a subsection of the 2-D
%   convolution with size specified by SHAPE:
%     'full'  - (default) returns the full 2-D convolution,
%     'same'  - returns the central part of the convolution
%               that is the same size as A.
%     'valid' - returns only those parts of the convolution
%               that are computed without the zero-padded edges.
%               size(C) = max([ma-max(0,mb-1),na-max(0,nb-1)],0).
%
%   See also CONV2, CONVN, CONVNFFT
%
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-April-2014
if length(varargin)>=3 && isnumeric(varargin{3})
    [H1, H2, A] = deal(varargin{1:3});
    varargin = varargin(4:end);
    C = convnfft(H1(:), A, varargin{:});
    C = convnfft(H2(:).', C, varargin{:});
else
    C = convnfft(varargin{:});
end
end % conv2fft

function A = convnfft(A, B, shape, dims, options)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
% 
%   C = CONVNFFT(A, B, SHAPE) controls the size of the answer C:
%       'full'   - (default) returns the full N-D convolution
%       'same'   - returns the central part of the convolution that
%                  is the same size as A.
%       'valid'  - returns only the part of the result that can be
%                  computed without assuming zero-padded arrays.
%                  size(C,k) = max([nak-max(0,nbk-1)],0).
%
%   C = CONVNFFT(..., SHAPE, DIMS) with DIMS is vector of dimensions where
%       the convolution will be carried out. By default DIMS is
%       [1:max(ndims(A),ndims(B))] (all dimensions). A and B must have the
%       same lengths on other dimensions.
%   C = CONVNFFT(..., SHAPE, DIMS, GPU)
%       GPU is boolean flag, see next
%
%   C = CONVNFFT(..., SHAPE, DIMS, OPTIONS)
%
%   OPTIONS is structure with following optional fields
%       - 'GPU', boolean. If GPU is TRUE Jacket/GPU FFT engine will be used
%       By default GPU is FALSE.
%       - 'Power2Flag', boolean. If it is TRUE, use FFT with length rounded
%       to the next power-two. It is faster but requires more memory.
%       Default value is TRUE.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-Jun-2009
%       23-Jun-2009: correct bug when ndims(A)<ndims(B)
%       02-Sep-2009: GPU/JACKET option
%       04-Sep-2009: options structure
%       16-Sep-2009: inplace product
if nargin<3 || isempty(shape)
    shape = 'full';
end
if nargin<5 || isempty(options)
    options = struct();
elseif ~isstruct(options) % GPU options
    options = struct('GPU', options);
end
nd = max(ndims(A),ndims(B));
% work on all dimensions by default
if nargin<4 || isempty(dims)
    dims = 1:nd;
end
dims = reshape(dims, 1, []); % row (needed for for-loop index)
% GPU enable flag
GPU = getoption(options, 'GPU', false);
% Check if Jacket is installed
GPU = GPU && ~isempty(which('ginfo'));
% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
switch lower(shape)
    case 'full',
        ifun = @(m,n) 1:m+n-1;
    case 'same',
        ifun = @(m,n) ceil((n-1)/2)+(1:m);
    case 'valid',
        ifun = @(m,n) n:m;
    otherwise
        error('convnfft: unknown shape %s', shape);
end
classA = class(A);
classB = class(B);
ABreal = isreal(A) && isreal(B);
% Special case, empty convolution, try to follow MATLAB CONVN convention
if any(size(A)==0) || any(size(B)==0)
    szA = zeros(1,nd); szA(1:ndims(A))=size(A);
    szB = zeros(1,nd); szB(1:ndims(B))=size(B);
    % Matlab wants these:
    szA = max(szA,1); szB = max(szB,1);
    szC = szA;
    for dim=dims
        szC(dim) = length(ifun(szA(dim),szB(dim)));
    end
    A = zeros(szC,classA); % empty -> return zeros
    return
end
power2flag = getoption(options, 'Power2Flag', true);
if power2flag
    % faster FFT if the dimension is power of 2
    lfftfun = @(l) 2^nextpow2(l);
else
    % slower, but smaller temporary arrays
    lfftfun = @(l) l;
end
if GPU % GPU/Jacket FFT
    if strcmp(classA,'single')
        A = gsingle(A);
    else
        A = gdouble(A);
    end
    if strcmp(classB,'single')
        B = gsingle(B);
    else
        B = gdouble(B);
    end
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        % compute the FFT length
        l = lfftfun(m+n-1);
        % We need to swap dimensions because GPU FFT works along the
        % first dimension
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
            B = permute(B, swap);
        end
        A = fft(A,l);
        B = fft(B,l);
        subs{dim} = ifun(m,n);
    end
else % Matlab FFT
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        % compute the FFT length
        l = lfftfun(m+n-1);
        A = fft(A,l,dim);
        B = fft(B,l,dim);
        subs{dim} = ifun(m,n);
    end
end
 
if GPU
    A = A.*B;
    clear B
else
%     % inplace product to save 1/3 of the memory
%     inplaceprod(A,B);
    A = A.*B;
end
% Back to the non-Fourier space
if GPU % GPU/Jacket FFT
    for dim=dims(end:-1:1) % reverse loop
        A = ifft(A,[]);
        % Swap back the dimensions
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
        end        
    end   
else % Matlab IFFT  
    for dim=dims
        A = ifft(A,[],dim);
    end
end
% Truncate the results
if ABreal
    % Make sure the result is real
    A = real(A(subs{:}));
else
    A = A(subs{:});
end
% GPU/Jacket
if GPU
    % Cast the type back
    if strcmp(class(A),'gsingle')
        A = single(A);
    else
        A = double(A);
    end
end
end % convnfft
% Get defaut option
function value = getoption(options, name, defaultvalue)
% function value = getoption(options, name, defaultvalue)
    value = defaultvalue;
    fields = fieldnames(options);
    found = strcmpi(name,fields);
    if any(found)
        i = find(found,1,'first');
        if ~isempty(options.(fields{i}))
            value = options.(fields{i});
        end
    end
end

