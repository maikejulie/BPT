
clear

% First load data. The data is downloaded directly from the LAS server, then I
% just calculated the time-mean.   
load time_mean_press
load time_mean_pressbot 
load ecco_params 

%------------------------------------------------------------

% establish some standard values
lon_diff=1.111774765625000e+05;
Omega=7.2921e-5;%2*pi/T; 
f=2*Omega*sin(([min(lat_c):(lat_c(2)-lat_c(1)):max(lat_c)])*pi./180); % coriolis parameter at cell centre
beta_fgrid=nan(1,length(lat_c)); 
beta_fgrid(:,2:end)=((f(2:end)-f(1:end-1))./(lon_diff)); % Gives beta at the cell corner (the vorticity point)
beta_fgrid=(beta_fgrid'*ones(1,length(lon_u)));  % Just makes it 2D

% grid cell widths and heights etc. 
grid_width_u=(lon_diff.*cos(lat_c.*pi./180))*ones(1,length(lon_c));
grid_width_v=(lon_diff.*cos(lat_v*pi./180))*ones(1,length(lon_c));
depth_l_diff=[depth_l(1);diff(depth_l)];
depth_l_diff_glob=repmat(depth_l_diff,[1 length(lat_c) length(lon_c)]);
grid_width_u_glob=permute(repmat(grid_width_u,[1 1 length(depth_l)]),[3 1 2]);
grid_width_v_glob=permute(repmat(grid_width_v,[1 1 length(depth_l)]),[3 1 2]);

%------------------------------------------------------

%%% determine the pressure gradient term in the depth-integrated momentum equation. 

% calculate gradients. 
gradp_i=nan(23,160,360);
gradp_j=nan(23,160,360);
gradp_i(:,:,2:end)=-((hydpress_mean(:,:,2:end)-hydpress_mean(:,:,1:end-1)));
gradp_j(:,2:end,:)=-((hydpress_mean(:,2:end,:)-hydpress_mean(:,1:end-1,:)));
gradp_i(:,:,1)=-((hydpress_mean(:,:,1)-hydpress_mean(:,:,end)));
gradp_i=gradp_i./grid_width_u_glob;
gradp_j=gradp_j./lon_diff;

% vertically integrate grad(p) to bottom. 
gradp_int_i=squeeze(nansum(gradp_i(1:23,:,:).*depth_l_diff_glob(1:23,:,:),1));
gradp_int_j=squeeze(nansum(gradp_j(1:23,:,:).*depth_l_diff_glob(1:23,:,:),1));

%----------------------------------------------------------
%%% Calculate the form stress term in the depth-integrated momentum equation.

% vertically integrate p to bottom. 
p_int=squeeze(nansum(hydpress_mean(1:23,:,:).*depth_l_diff_glob(1:23,:,:),1));

% take gradient
grad_intp_i=nan(160,360);
grad_intp_j=nan(160,360);
grad_intp_i(:,2:end)=-((p_int(:,2:end)-p_int(:,1:end-1)));
grad_intp_j(2:end,:)=-((p_int(2:end,:)-p_int(1:end-1,:)));
grad_intp_i(:,1)=-((p_int(:,1)-p_int(:,end)));
grad_intp_i=grad_intp_i./grid_width_u;
grad_intp_j=grad_intp_j./lon_diff;

% The form stress term:
pb_gradH_i=gradp_int_i-grad_intp_i;
pb_gradH_j=gradp_int_j-grad_intp_j;

%----------------------------------------------------------

%%% Compute curl (in spherical polar coordinates) to get the depth-integrated vorticity equation. 

coslatc_glob=permute(repmat(cos(lat_c.*pi./180)*ones(1,360),[1 1 23]),[3 1 2]); % cosine of latitude (on the c-point) everywhere on the globe
coslatv_glob=permute(repmat(cos(lat_v.*pi./180)*ones(1,360),[1 1 23]),[3 1 2]); % cosine of latitude (on the v-point) everywhere on the globe

% first scale the zonal terms with cos(lat). Needed for the spherical coordinate calculation. 
% Spherical coordinate calculations are based on this page: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
gradp_i_scale=gradp_i.*coslatc_glob;
gradp_int_i_scale=gradp_int_i.*squeeze(coslatc_glob(1,:,:));
grad_intp_i_scale=grad_intp_i.*squeeze(coslatc_glob(1,:,:));
pb_gradH_i_scale=pb_gradH_i.*squeeze(coslatc_glob(1,:,:));


% Find the curl of the depth integrated pressure gradient. This is the term used for BPT. 
curl_gradp_int_kcomp=(padarray((gradp_int_i_scale(2:end,:)-gradp_int_i_scale(1:end-1,:)),[1 0],nan,'pre')./(squeeze(coslatv_glob(1,:,:)).*lon_diff))...
    -(padarray((gradp_int_j(:,2:end)-gradp_int_j(:,1:end-1)),[0 1],nan,'pre')./grid_width_v);
%and calculate the first longitude index
curl_gradp_int_kcomp(:,1)=(padarray((gradp_int_i_scale(2:end,1)-gradp_int_i_scale(1:end-1,1)),[1 0],nan,'pre')./(coslatv_glob(1,:,1)'.*lon_diff))...
    -((gradp_int_j(:,1)-gradp_int_j(:,end))./grid_width_v(:,1));

% Curl of the gradient of pressure. This is just a verification step to
% verify I am taking the curl correctly. The curl of a gradient is zero everywhere. 
curl_gradp_kcomp=(padarray((gradp_i_scale(:,2:end,:)-gradp_i_scale(:,1:end-1,:)),[0 1 0],nan,'pre')./(coslatv_glob.*lon_diff))...
    -(padarray((gradp_j(:,:,2:end)-gradp_j(:,:,1:end-1)),[0 0 1],nan,'pre')./grid_width_v_glob);
%and calculate the first longitude index
curl_gradp_kcomp(:,:,1)=(padarray((gradp_i_scale(:,2:end,1)-gradp_i_scale(:,1:end-1,1)),[0 1 0],nan,'pre')./(coslatv_glob(:,:,1).*lon_diff))...
    -((gradp_j(:,:,1)-gradp_j(:,:,end))./grid_width_v_glob(:,:,1));

% Another verification step to show that the curl of the gradient of
% depth-integrated pressure is also zero, as it should be. 
curl_grad_intp_kcomp=(padarray((grad_intp_i_scale(2:end,:)-grad_intp_i_scale(1:end-1,:)),[1 0],nan,'pre')./(squeeze(coslatv_glob(1,:,:)).*lon_diff))...
    -(padarray((grad_intp_j(:,2:end)-grad_intp_j(:,1:end-1)),[0 1],nan,'pre')./grid_width_v);

% A final verification step to show that the curl of the form stress term
% equals the curl of the depth integrated pressure gradient. 
curl_pb_gradH_kcomp=(padarray((pb_gradH_i_scale(2:end,:)-pb_gradH_i_scale(1:end-1,:)),[1 0],nan,'pre')./(squeeze(coslatv_glob(1,:,:)).*lon_diff))...
    -(padarray((pb_gradH_j(:,2:end)-pb_gradH_j(:,1:end-1)),[0 1],nan,'pre')./grid_width_v);



%-------------------------------------------------

%%% Calculate BPT (in units m^2/s)

BPT=-curl_gradp_int_kcomp./beta_fgrid; 
figure, pcolor(BPT), shading flat, colorbar, caxis([-10 10])



