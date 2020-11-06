%% IQE vs Thickness plot generation

%Use this for comparing EQE for different thicknesses of the active layer
clear;
[J,J1,J2,R]=Stack();



function [Jsc,Jsc1,Jsc_Si,R] = Stack()
%% new instance of TMat Script

spectrum=370:1:800;            %examined spectrum for R,T, in nm

%% Solar Panel Model
% standard model = 5 layers: glass, ITO, PEDOT/PSS, the active layer (PCPDTBT/PC60BM) and the aluminum 
% (PEDOT/PSS is poly(3,4-ethylenedioxythiophene)/poly(styrene sulfonate) and ITO is indium tin oxide)

% initiate layer parameters
layers = {'Air' 'Glass' 'ITO2' 'lukas_triple' 'Si'}; % Names of layers of materials starting from side light is incident from
                                                                            % Names must match the strings in the excel n,k source file
thicknesses = [100000 120 150];                                            % thickness of each corresponding layer in nm
                                                                                                                                         
active_layer=3; %in reference to thicknesses

% Load in index of refraction for each material
ntotal = zeros(size(layers,2),size(spectrum,2));
for index = 1:size(layers,2)
    ntotal(index,:) = LoadRefrIndex(layers{index},spectrum);
end

% Constants
h = 6.62606957e-34; 	% Js Planck's constant
c = 2.99792458e8;	% m/s speed of light
q = 1.60217657e-19;	% C electric charge 

%% Calculate the R,T for the spectrum

t = thicknesses;
ls=length(spectrum);
i=0;

%calculate total thickness for x-axis, as well as layer positions on it
dTOT=0;
pos(1)=0;
for i=2:length(t)
    if(t(i)>99999)
        break
    else
        dTOT=dTOT+t(i);
        pos(i)=pos(i-1)+t(i);
    end

end


%% Absorbance and Generation
% absorption coefficient a gives the fraction of incident
% radiant energy absorbed per unit thickness (usually cm-1 but here taken m-1)
% The absorbed intensity in one layer is therefore int(a*I(x)dx)
%I=some constant*n*E2 but we use n*E2(x) instead of I(x) since absorbance is a
%relative value (I/I0)

% Load in 1sun AM 1.5 solar global tilt in mW/cm2
AM=xlsread('AM15.xls');                                 

i=0;
Abs_lay=zeros(length(t),ls);
G=zeros(length(t),ls);
for lam=spectrum
    i=i+1;
    n = ntotal(:,i:i).';
    k_Si=imag(n(length(n)-1));
    t_Si=t(length(t));
    [R(i),T(i),E2_lam]=Tmat(lam,n,t,0);
    %assumption that Si absorbs all
    A_Si(i)=T(i);
    A_metal(i)=0;
    
    AM15=AM(find(AM==lam),2);  %load in AM1.5, in mW/cm2/nm  
    G2_L(lam)=AM15*A_Si(i) *lam/h/c; %generation in Si incoh layer
    
    for lay=2:length(t) %Si and metal absorption added later

        alpha=zeros(1,dTOT+1); %make sure to reset alpha for each layer
        a=4*pi*imag(n(lay+1))/lam/1e-9; %in metres (usually in cm)
        for z=pos(lay-1):pos(lay)-1
            alpha(z+1)=a*real(n(lay+1))/real(n(1))*E2_lam(z+1);  %calculate Qj for a certain position/layer
                                                                %note that alpha=Iabs/I0 means we can leave out constants, 
                                                                %but we still need to include n of I0. |E+0|=1
            if(lay==active_layer)
                G(z+1,lam)=a*real(n(lay+1))*AM15*E2_lam(z+1) *lam*1e-9/h/c; %# of generated pairs in the active layer per m2, at zi
                                                                                  %G=a*AM15*E2*nl /(hf) = a*AM15*E2*nl /(hc)*lam
                                                                                  %to see why I=AM15*n*E2, see notes p.13
            end

        end
        
        Abs_lay(lay,i)=sum(alpha)*1e-9;  %stepsize for the integral must match metres unit       
    end

    %add other absorptances
    Abs_lay(length(t)+1,i)=A_Si(i);
    %Abs_lay(length(t)+1,i)=A_metal(i); always 0
    
    %total absorptance
    Absorption(i)=sum(Abs_lay(:,i));
end
parasitic=1-R-sum(Abs_lay(active_layer:active_layer+1,:));
%next plot Absorption
figure(1)
plot(spectrum,R,spectrum,parasitic,spectrum,Abs_lay(:,:),'LineWidth',2)
stringl={'Reflection','Parasitic Abs'};
for i=3:length(layers)+1
    stringl(i)=layers(i-1);
end
legend(stringl);

%% Current density
%Calculate G(x) from I(x), J(x) from G(x)

Jsc1=sum(G(:))*1e-9*q  %in mA/cm^2 (q = 1.60217657e-19;	% C electric charge) Perovskite current
Jsc_Si=sum(G2_L(:))*1e-9*q   %incoherent layer (silicon) current

Jsc=Jsc1+Jsc_Si;

end



%% Function LoadRefrIndex 
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library2.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Optical_Constants.xls');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch(strcat(name,'_lambda'),IndRefr_names));

n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));   %has to be minus for some reason
Nan=find(isnan(file_wavelengths));
if(~isempty(Nan))
    file_wavelengths=file_wavelengths(1:Nan(1)-1);
    n=n(1:Nan(1)-1);
    k=k(1:Nan(1)-1);
end
% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
end

function [EQEe] = LoadEQE(name,wavelengths)

[IndRefr,IndRefr_names]=xlsread('EQE.xls');
file_wavelengths=IndRefr(:,strmatch(strcat(name,'_lambda'),IndRefr_names));

p=IndRefr(:,strmatch(strcat(name,'_p'),IndRefr_names));

% Interpolate/Extrapolate data linearly to desired wavelengths
p_interp=interp1(file_wavelengths, p, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
EQEe = p_interp.';
end
