clear;
clc;


c=3e8;% light speed
wavelength=488e-9; %wavelength of illumination light in nanometers
k0=2*pi/wavelength;%wave vector of light wave in the vacuum
resize_time=4; % scaling factor to speed up the simulation
diameter_lens=2e-3/resize_time; %diameter of the illumination objective lens
radius=diameter_lens/2; %radius of the illumination objective lens
NA=0.65; %NA of the objective lens
focal_length=radius/tan(asin(NA));%focal length of the objective lens

%%%%%%%%%%  defination of the electic field (xz plane) after the annular ring mask
%%%%%%%%%%
x=linspace(-radius,radius,radius*2/500e-9);  % pixel: 500 nm
z=linspace(-radius,radius,radius*2/500e-9);

NA_outer=0.35;%outer NA of the annular ring
radius_outer=focal_length*tan(asin(NA_outer));%outer annular ring radius
NA_inner=0.07;%inner NA of the annular ring
radius_inner=focal_length*tan(asin(NA_inner));%inner annular ring radius
NA_x_lines=0.0;%NA location of the Bessel line beams at the rear pupil plane
x_lines=focal_length*tan(asin(NA_x_lines)); % x postion of the Bessel line beams at the rear pupil plane

length_1=sqrt(radius_outer^2-x_lines^2)-sqrt(radius_inner^2-x_lines^2);%length of the a1 line beams: L1
z_c1=(sqrt(radius_outer^2-x_lines^2)+sqrt(radius_inner^2-x_lines^2))/2;% distances from the center points of the a1 line beamlets to the x-axis
ratio=1/2; % ratio of z_c2/z_c1 obtained by the theoretical derivation
z_c2=z_c1*ratio;
length_2=(z_c2-sqrt(radius_inner^2-x_lines^2))*2;%length of the a2 line beams: L2
width_1=1.25e-5/resize_time; %width of the a1 line beams: W1
width_2=width_1*length_1/length_2; %width of the a2 line beams: W2

amplitude_incident=zeros(length(z),length(x)); 
phase_incident=zeros(length(z),length(x));
condition=zeros(length(z),length(x));
for i=1:length(x)
    for j=1:length(z)
%        condition(j,i)=((x(i)^2+z(j)^2)<=radius_outer^2)&&((x(i)^2/1^2+z(j)^2/1)>=radius_inner^2)&&(abs(abs(x(i))-x_lines)<=width_1/2);% a1-line beams
%        condition(j,i)=((x(i)^2+z(j)^2)<=radius_outer^2)&&((x(i)^2/1^2+z(j)^2/1)>=radius_inner^2)&&(abs(abs(x(i))-x_lines)<=width_2/2)&&abs(abs(z(j))-z_c2)<=length_2/2;    % a2-line beams
        condition(j,i)=((x(i)^2+z(j)^2)<=radius_outer^2)&&((x(i)^2/1^2+z(j)^2/1)>=radius_inner^2)&&((abs(abs(x(i))-x_lines)<=width_1/2)||(abs(abs(x(i))-x_lines)<=(width_1+width_2)/2)&&abs(abs(z(j))-z_c2)<=length_2/2);    % a3-line beams.
    
       if  condition(j,i)
           amplitude_incident(j,i)=1;
       end
    end
end
 U_incident=amplitude_incident.*exp(1i*phase_incident); % electric field after the annular ring mask
 figure(1); 
 surf(x,z,amplitude_incident.^2); 
 xlabel('x (m)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 ylabel('z (m)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 view(2);
 shading interp;
 axis  equal
 colormap(hot)
  

%%%%%%%%%%  calculation of the electic field (xz plane) before the objective lens
%%%%%%%%%%  
tic
x_lens=x;
z_lens=z;
y_lens=focal_length;
E_xz_before_lens=zeros(length(z_lens),length(x_lens));
[X,Z]=meshgrid(x,z);

%%% simulation method:  straightforward diffraction integral approach, according to formula 3-51 in book "Introduction to fourier optics" (J. W. Goodman,3rd) 
inner_n=length(x_lens);
  parfor m=1:length(z_lens)
        E_xz_inner=zeros(1,inner_n);
    for n=1:inner_n
              distance=sqrt((X-x_lens(n)).^2+(Z-z_lens(m)).^2+(0-y_lens)^2);                 
              E=(-1j/(wavelength)).*(y_lens./distance)./distance.*exp(1i.*(distance./wavelength.*2.*pi)).*U_incident;
              E_lens=sum(sum(E));
              E_xz_inner(n)=E_lens;
      end
        E_xz_before_lens(m,:)=E_xz_inner; %E_xz_before_lens: electric field at the plane just before the objective lens
 end
 toc 
 figure(2);
 surf(x_lens,z_lens,abs(E_xz_before_lens).^2/max(max(abs(E_xz_before_lens).^2))); 
 xlabel('x (m)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8); 
 ylabel('z (m)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 view(2);
 shading interp;
 axis  equal
 colormap(hot)
 
 
%%%%%%%%%%% calculation of the phase shift due to the objective lens
%%%%%%%%%%%
phase_lens=zeros(length(x_lens),length(z_lens));
for i=1:length(x_lens)
    for j=1:length(z_lens)
            light_path_difference=sqrt((radius)^2+focal_length^2)-sqrt(x_lens(i)^2+z_lens(j)^2+focal_length^2);
            phase_lens(i,j)=light_path_difference/wavelength*2*pi;
    end
end
U_after_lens=E_xz_before_lens.*exp(1i*phase_lens);  % U_after_lens: the electric field at the plane right behind the objective lens

 
%%%%%%%%%%% focal electric field: yz plane (x = 0)
%%%%%%%%%%%
tic
y_imaging=linspace(focal_length*2-15*wavelength, focal_length*2+15*wavelength,2*15*10);
z_imaging=linspace(-8*wavelength, 8*wavelength,2*8*10);
x_imaging=0;
     
E_yz_at_focal=zeros(length(z_imaging),length(y_imaging));
[X_lens,Z_lens]=meshgrid(x_lens,z_lens);
          
inner_n=length(y_imaging);
  parfor m=1:length(z_imaging)
      E_yz_inner=zeros(1,inner_n);
    for n=1:inner_n         
                distance=sqrt((X_lens-x_imaging).^2+(Z_lens-z_imaging(m)).^2+(y_imaging(n)-y_lens)^2);
                E=(-1i./(wavelength)).*((y_imaging(n)-y_lens)./distance)./distance.*exp(1i.*(distance./wavelength.*2.*pi)).*U_after_lens;
                E_focal=sum(sum(E));
                E_yz_inner(n)=E_focal;
    end
    E_yz_at_focal(m,:)=E_yz_inner; % E_yz_at_focal: focal electric field in yz plane 
 end
 toc
 
 figure(3);
 surf((y_imaging-2*focal_length)/1e-6,z_imaging/1e-6,abs(E_yz_at_focal).^2/max(max(abs(E_yz_at_focal).^2)));
 view(2);
 xlim([min((y_imaging-2*focal_length)/1e-6),max((y_imaging-2*focal_length)/1e-6)]);
 ylim([min(z_imaging/1e-6),max(z_imaging/1e-6)]);
 xlabel('distance from the focal point \mum','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 ylabel('z \mum','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 shading interp;
 axis  equal
 colormap(hot)
 
 %%%%%%% normalized field intensity along y axis (x=0,z=0): y-profile
 E_y_intensity=abs(E_yz_at_focal(length(z_imaging)/2,:)).^2/max(max(abs(E_yz_at_focal(length(z_imaging)/2,:)).^2)); 
 figure(4);
 plot((y_imaging-2*focal_length)/1e-6,E_y_intensity);
 xlim([min((y_imaging-2*focal_length)/1e-6),max((y_imaging-2*focal_length)/1e-6)]);
 ylim([0 1.2]);
 xlabel('distance from the focal point \mum','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 ylabel('normalized intensity','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 
%%%%%%% calculation of illumination uniformity within the effective propagation distance: +-2.5 um from the focal point 
for m=1:length(y_imaging)
    if (y_imaging(m)-2*focal_length)>-2.5e-6
           bound_left_point_num=m;
         break
    end
end

for m=1:length(y_imaging)
    if (y_imaging(m)-2*focal_length)>2.5e-6
           bound_right_point_num=m;
         break
    end
end

E_effective_distance=E_y_intensity(bound_left_point_num:bound_right_point_num);
illumination_variation=(max(E_effective_distance)-min(E_effective_distance))/max(E_effective_distance);

%%%%%%% calculation of propagation length: FWHM of the electric field y-profile
y_imaging_interp=linspace(min(y_imaging), max(y_imaging),5000);
E_y_intensity_interp=interp1(y_imaging,E_y_intensity,y_imaging_interp,'linear');
half_y=y_imaging_interp(length(y_imaging_interp)/2:length(y_imaging_interp));
half_E_z_intensity =E_y_intensity_interp(length(y_imaging_interp)/2:length(y_imaging_interp));

max_diff=max(abs(diff(half_E_z_intensity)));
FWHM_point_num=find(abs(half_E_z_intensity-0.5)<=max_diff);
FWHM_y=(half_y(FWHM_point_num(1))-2*focal_length)*2;
 
%%%%%%%%%%% focal electric field: xz plane 
%%%%%%%%%%%
tic
x_imaging=linspace(-10*wavelength, 10*wavelength,2*10*10);
z_imaging=linspace(-10*wavelength, 10*wavelength,2*10*10);
y_imaging=focal_length*2;     
E_xz_at_focal=zeros(length(z_imaging),length(x_imaging));
[X_lens,Z_lens]=meshgrid(x_lens,z_lens);
          
 inner_n=length(x_imaging); 
 parfor m=1:length(z_imaging)
        E_xz_inner=zeros(1,inner_n);        
    for n=1:inner_n                                                                   
                 distance=sqrt((X_lens-x_imaging(n)).^2+(Z_lens-z_imaging(m)).^2+(y_imaging-y_lens)^2);                 
                 E_focal=(-1i./(wavelength)).*((y_imaging-y_lens)./distance)./distance.*exp(1i.*(distance./wavelength.*2.*pi)).*U_after_lens;                 
                 E=sum(sum(E_focal));                                        
                 E_xz_inner(n)=E;
    end        
    E_xz_at_focal(m,:)=E_xz_inner; % E_xz_at_focal: focal electric field in xz plane
 end 
 toc
 
 figure(5);
 surf(x_imaging/1e-6,z_imaging/1e-6,abs(E_xz_at_focal).^2/max(max(abs(E_xz_at_focal).^2)));
 view(2);
 xlim([min(x_imaging/1e-6),max(x_imaging/1e-6)]);
 ylim([min(z_imaging/1e-6),max(z_imaging/1e-6)]);
 xlabel('x \mum','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8); 
 ylabel('z \mum','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
 shading interp;
 axis  equal 
 colormap(hot)


%%%%%%%%%% z-profile of the dithered light sheet
%%%%%%%%%%
dithered_intensity = mean(abs(E_xz_at_focal).^2/max(max(abs(E_xz_at_focal).^2)),2);%%% dithered intensity calculation by taking the mean value of electric field intensity in x direction
dithered_intensity =dithered_intensity ./max(dithered_intensity );
dithered_intensity =dithered_intensity';

figure(6);
plot(z_imaging/1e-6,dithered_intensity,'w','LineWidth',2)
set(gcf,'color','k');
set(gca,'xcolor','w','ycolor','w','color','k');
xlim([min(z_imaging/1e-6) max(z_imaging/1e-6)]);
ylim([0 1.2]);
set(gca,'xtick',-floor(max(z_imaging/1e-6)):1:floor(max(z_imaging/1e-6)));
set(gca,'FontName','Times New Roman','FontSize',22);
xlabel('z (\mum)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8); 
grid on
set(gca,'gridcolor','w','LineWidth',0.5,'GridAlpha',1,'GridLineStyle',':');


%%%%%%%%%%% calculation of P_side/P_main and FWHM using interpolated values
%%%%%%%%%%% of the dithered intensity profile 
z_imaging_interp=linspace(min(z_imaging), max(z_imaging),5000);
intensity_interp=interp1(z_imaging,dithered_intensity,z_imaging_interp,'linear');
half_z=z_imaging_interp(length(z_imaging_interp)/2:length(z_imaging_interp));
half_z_intensity =intensity_interp(length(z_imaging_interp)/2:length(z_imaging_interp));

%find the first dip point of the intensity profile 
for n=1:length(half_z)-2
          if sign(half_z_intensity(n+2)-half_z_intensity(n+1))*sign(half_z_intensity(n+1)-half_z_intensity(n))<0
              dip_piont_num=n+1;
              break
          end
end

%calculation of P_side/P_main
half_P_side=trapz(half_z(dip_piont_num:end),half_z_intensity(dip_piont_num:end));
half_P_main=trapz(half_z(1:dip_piont_num),half_z_intensity(1:dip_piont_num));
T=half_P_side/half_P_main;

%calculation of FWHM
max_diff=max(abs(diff(half_z_intensity)));
FWHM_point_num=find(abs(half_z_intensity-0.5)<=max_diff);
FWHM=half_z(FWHM_point_num(1))*2;

% load  relative_P_1st_side_Bessel_line.mat   % P_1st_side for the a1 line beams
% half_P_1st_side=trapz(half_z(point_num_start:point_num_end),half_z_intensity(point_num_start:point_num_end));
% %here 'point_num_start' is corresponding to the left dip edge of the first side lobe,  and 'point_num_end' the right dip edge  
% P_1st_side_area=half_P_1st_side/half_P_1st_side_a1;



%%%%%%%%%% Fourier transform calculation for z-profile  (to verify the results of z-profile obtained by 
%%%%%%%%%% diffraction integration method, no corresponding figures in the manuscript)
%%%%%%%%%%

%%%%  kz spatial frequency space array
k01=1/wavelength;%spatial frequency of light wave in the vacuum
kz=linspace(-20*k01,20*k01,50000);

%%%%% kz distribution and focal electric field z-profile calculation: a1 line beam
kz_a1_max=k01*sqrt(radius_outer^2-x_lines^2)/(sqrt(focal_length^2+radius_outer^2)); 
kz_a1_min=k01*sqrt(radius_inner^2-x_lines^2)/(sqrt(focal_length^2+radius_inner^2)); 
for i=1:length(kz)
       if (abs(kz(i))>=kz_a1_min)&&(abs(kz(i))<=kz_a1_max)
           F_kz_a1(i)=1; % F_kz_a1: kz distribution of a1 line beams
       else
           F_kz_a1(i)=0;
       end
end
F_z_a1=fft(F_kz_a1);
F_z_a1=fftshift(F_z_a1);%F_z_a1: focal electric field z-profile in real space
z=1/(40*k01)*(1:length(F_z_a1));
z=z-max(z)/2;

%%%%% kz distribution and focal electric field z-profile calculation: a2 line beam
kz_c1=(kz_a1_max+kz_a1_min)/2; 
kz_c2=kz_c1*ratio;
kz_range_1=(kz_a1_max-kz_a1_min); 
kz_range_2=(kz_c2-kz_a1_min)*2;
kz_a2_max=kz_c2+kz_range_2/2; 
kz_a2_min=kz_a1_min;

for i=1:length(kz)
       if (abs(kz(i))<=kz_a2_max)&&(abs(kz(i))>=kz_a2_min)
           F_kz_a2(i)=1;% F_kz_a1: kz distribution of a2 line beams
       else
           F_kz_a2(i)=0;
       end
end
F_z_a2=fft(F_kz_a2);
F_z_a2=fftshift(F_z_a2);%F_z_a2: focal electric field z-profile in real space

%%%%% kz distribution and focal electric field z-profile calculation: a3 line beam
F_z_a3=(F_z_a1+F_z_a2*kz_range_1/kz_range_2); %F_z_a3: focal electric field z-profile in real space


%%%%%% compare the z-profiles obtained of a1 (or a2, a3) lines by diffraction integration method and Fourier transform method  
figure(7)
plot(z/1e-6,abs(F_z_a3.^2)/max(abs(F_z_a3.^2)),'LineWidth',2,'Color','k');
set(gca,'FontName','Times New Roman','FontSize',22);
z_imaging_range=max(z_imaging);
xlim([-z_imaging_range/1e-6  z_imaging_range/1e-6]);
ylim([0 1.4]);
xlabel('z (\mum)');
ylabel('normalized field intensity');
set(gca,'xtick',-floor(z_imaging_range/1e-6):1:floor(z_imaging_range/1e-6));
grid on

hold on
plot(z_imaging./1e-6,dithered_intensity,'^m')
xlabel('z (\mum)','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
ylabel('Normalized Intensity ','FontSize',22,'FontName','Times New Roman','FontWeight','bold','LineWidth',8);
whitebg(gcf,[1,1,1]);
h=legend('fft calculation','integration simulation');
set(h, 'Box', 'off')




 

