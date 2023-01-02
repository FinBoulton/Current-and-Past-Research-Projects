% RAY-TRACING RK4fn method for C(x) = (r,z)
%% data for ssp over range

clear

%location: is path connecting two points
%current: south pacific 
lat=[-45 -45];
lon=[-160 -158];

%range discretisation number
rdn = 500 ;

%func to generate profile, quite high sample of 50 
[r, glat, glon]=dist(lat,lon,rdn);
[P1,z]=get_lev(glat,glon,'c');

%plot of ssp
%{

figure
plot(P1,-z/1000); xlabel('Sound Speed (m/s)'); ylabel('Depth (km)')
meanSSP = mean(P1,2);
hold on
plot(meanSSP,-z/1000,'kx','LineWidth',2)
hold off
%}

%plot of variation in SSP along path
%{
%figure
%surf(r,-z/1000,P1)
%xlabel('Range')
%ylabel('Depth (km)'); title('Sound speed over range')
%}

%array with interpolated ssp
%figure
C = zeros(5501,500) ;
for i = 1: rdn
    zq = 0:1:5500 ;
    c = interp1(z, P1(:,i), zq, 'linear') ;
    C(:,i) = c ;
    %plot(c,-zq/1000)
    %hold on
end
%xlabel('Sound Speed (m/s)'); ylabel('Depth (km)')

%data variables 
depth_lim = 5500.0 ;

%bottom fn symbolic 
%{  

syms x
y = 0.02*x -5.5 ;
figure
fplot(y,'k')
%}

%boundary points positive gradient 
%{
nb=140000;
m = 0.02; %gradient of bottom
bz    = zeros(1, nb) ; %pre-allocate array
br    = zeros(1, nb) ;
bz(1) = 5500 ; %array first value
br(1) = 0 ;
for i= 1:nb-1 %loop to create boundary points 
    bz(i+1) = bz(i) - m ;
    br(i+1) = br(i) + 1;
end
depth_min = min(abs(bz)) ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}

%alt Boundary negative gradient
%{
%boundary points
nb=140000;
m = -0.02; %gradient of bottom
bz    = zeros(1, nb) ; %pre-allocate array
br    = zeros(1, nb) ;
bz(1) = 2700 ; %array first value
br(1) = 0 ;
for i= 1:nb-1 %loop to create boundary points 
    bz(i+1) = bz(i) - m ;
    br(i+1) = br(i) + 1;
end
depth_min = min(abs(bz)) ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}
 
%alt Boundary zero gradient
%{
%boundary points
nb=140000;
m = 0.0 ; %gradient of bottom
bz    = zeros(1, nb) ; %pre-allocate array
br    = zeros(1, nb) ;
bz(1) = 5500 ; %array first value
br(1) = 0 ;
for i= 1:nb-1 %loop to create boundary points 
    bz(i+1) = bz(i) - m ;
    br(i+1) = br(i) + 1 ;
end
depth_min = min(abs(bz)) ;
%depth_min = 2000;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}

%alt Boundary steep gradient
%{
%boundary points
nb=140000;
m = 0.04; %gradient of bottom
bz    = zeros(1, nb) ; %pre-allocate array
br    = zeros(1, nb) ;
bz(1) = 5500 ; %array first value
br(1) = 0 ;
for i= 1:nb-1 %loop to create boundary points 
    bz(i+1) = bz(i) - m ;
    br(i+1) = br(i) + 1;
end
depth_min = min(abs(bz)) ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}

%seamount
%{
n=140000;
m = 0.09; %gradient of bottom
bz    = zeros(1, n) ; %pre-allocate array
br    = zeros(1, n) ;
bz(1) = 5500 ; %array first value
br(1) = 0 ;
for i = 1:20000 %loop to create boundary points
    bz(i+1) = bz(i) ;
    br(i+1) = br(i) + 1;
end
for i = 20000:50000
    bz(i+1) = bz(i) - m ;
    br(i+1) = br(i) + 1 ;
end
for i = 50000:55000
    bz(i+1) = bz(i) ;
    br(i+1) = br(i) + 1;
end
for i = 55000:85000
    bz(i+1) = bz(i) + m ;
    br(i+1) = br(i) + 1 ;
end
for i = 85000:n
    bz(i+1) = bz(i) ;
    br(i+1) = br(i) + 1;
end
depth_min = min(abs(bz)) ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}

%normal distrobution
%{%}
nb = 140000 ;
mu = 65000 ;
sigma = 20000 ;
br=0:nb ;
bz = -100000000*0.9*pdf('normal',br,mu,sigma) + 4000;
depth_min = min(abs(bz)) - 10 ;
%test plot
%figure
%plot(br/1000,-bz/1000)


%normal dist but trough not peak 
%{
nb = 140000 ;
mu = 65000 ;
sigma = 20000 ;
br=0:nb ;
bz = 100000000*1.2*pdf('normal',br,mu,sigma) + 2000;
depth_min = min(abs(bz)) - 10 ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}

%steep normal distrobution
%{
nb = 140000 ;
mu = 65000 ;
sigma = 20000 ;
br=0:nb ;
bz = -100000000*1.7*pdf('normal',br,mu,sigma) + 4000;
depth_min = min(abs(bz)) - 10 ;
%test plot
%figure
%plot(br/1000,-bz/1000)
%}


%% Perform the angle iteration

%iteration conditions
ds = 1.0 ;
n = round(140000.0/ds) ;

%for many angle
a1 = 0.01 ;
a2 = 8.10 ;
a3 = 15.0 ;
a4 = 11.6 ;
a5 = 13.0 ;

angle_arr = [a1, a2, a3] ; %angle array


figure
plot(br/1000,-bz/1000,'k.','DisplayName', 'Bottom') ; hold on
for j = 1:length(angle_arr)
    
    k=1;
    
    % Set initial conditions
    r0 = 0.0 ; 
    z0 = 0.0 ;
    xi0   = cosd(angle_arr(j)) / C(1,1) ;
    zeta0 = sind(angle_arr(j)) / C(1,1) ;
    
    %function
    f = [r0, z0, xi0, zeta0] ;

    % Arrays to save solution for plotting
    r_arr    = [f(1)] ; z_arr    = [f(2)] ;
    xi_arr   = [f(3)] ; zeta_arr = [f(4)] ;
    
    %sound speed
    c_arr    = zeros(1, n) ;
    c_arr(1) = C(1,1) ;
    
    %distance arrarys
    dr = zeros(1,n);
    dz = zeros(1,n);
    
    % Calculate for i = 2 via Euler Method
    r_arr(2)    = r_arr(1) + ds * c_arr(1) * xi_arr(1) ;
    z_arr(2)    = z_arr(1) + ds * c_arr(1) * zeta_arr(1) ;

    c_arr(2)    = interp1(zq, C(:,1), z_arr(2), 'linear') ;

    xi_arr(2)   = xi_arr(1)   - ds * (1/c_arr(1)^2 * (c_arr(2)-c_arr(1))/(r_arr(2)-r_arr(1))) ;
    zeta_arr(2) = zeta_arr(1) - ds * (1/c_arr(1)^2 * (c_arr(2)-c_arr(1))/(z_arr(2)-z_arr(1))) ;

    %f = [r_arr(2), z_arr(2), xi_arr(2), zeta_arr(2)] ;
   
    yes = 0 ; %flag start condition
    
    % RK4 iteration from seecond step and backwards
    for i = 2:n-1

        f  = [r_arr(i)  , z_arr(i)  , xi_arr(i)  , zeta_arr(i)] ;
        fi = [r_arr(i-1), z_arr(i-1), xi_arr(i-1), zeta_arr(i-1)] ;

        eta1 = F(f            , fi            , c_arr(i), c_arr(i-1)) ;
        eta2 = F(f + ds/2*eta1, fi + ds/2*eta1, c_arr(i), c_arr(i-1)) ;
        eta3 = F(f + ds/2*eta2, fi + ds/2*eta2, c_arr(i), c_arr(i-1)) ;
        eta4 = F(f + ds  *eta3, fi + ds  *eta3, c_arr(i), c_arr(i-1)) ;
 
        % f = [r(i+1), z(i+1), xi(i+1), z(i+1)]
        f = f + ds/6*(eta1 + 2*eta2 + 2*eta3 + eta4) ;

        %{%}
        %angle of reflection check, needds = ang incidence
        if yes == 1
            beta = acosd(( xi_arr(i)*(-tangent_z) + zeta_arr(i)*tangent_r )/(sqrt(xi_arr(i)^2 + zeta_arr(i)^2)*sqrt((-tangent_z)^2 + tangent_r^2))) ;
            disp(beta-90)
        end
        
        
        r_arr(i+1)    = f(1) ;
        z_arr(i+1)    = f(2) ;
        xi_arr(i+1)   = f(3) ;
        zeta_arr(i+1) = f(4) ;
        
        %{%}
        %backward ray prevention
        if r_arr(i+1)<r_arr(i)
            disp('Backward Ray')
            for y=i:n-1
               r_arr(y+1) = r_arr(i) ; %fill r,z arrays with current point
               z_arr(y+1) = z_arr(i) ;
            end
            %r_tst(i+1) = r_tst(i) ;
            %z_tst(i+1) = z_tst(i) ;
            break
        end
        
        
        %BR surf z = 0.
        if z_arr(i+1) < 0
            r_arr(i+1) = (r_arr(i)*z_arr(i+1) - r_arr(i+1)*z_arr(i)) / (z_arr(i+1) - z_arr(i)) ;
            z_arr(i+1) = 0 ;
            zeta_arr(i+1) = -zeta_arr(i+1) ;
            %disp('Func')
            %disp(r_arr(i+1))
            %c_arr(i+1) = interp1(depth, sspd, z_arr(i+1), 'makima') ;
        end
        
        %BR bottom
        yes = 0 ;
        ind = 1 ;
        flag = 0 ; %reset flag parameter 
        %ind  = 0 ; %for saving value of q for closest boundary point
        if z_arr(i+1) > depth_min % only run below known min depth. boundary detection has large perfomance penalty
            for q=1:nb-1 %loop overall boundary points 
                dr(q) = br(q) - r_arr(i+1) ; % vectors pointing from x to boundary 
                dz(q) = bz(q) - z_arr(i+1) ;
                if q >= 2 %need to test q-1 < q as q+1 will always be zero(from array pre-allocation)
                    if ((dr(q-1))^2 + (dz(q-1))^2) < ((dr(q))^2 + (dz(q))^2) % test if distance of boundary to x(i+1) is minimum
                        ind = q ; %save index of boundary points
                        flag=1 ; %not sure if flag or break method is best 
                        %break
                    end
                end
                %{ %}
                if flag ==1 %to stop for loop
                    break
                end   
            end
            %distances 
            D0r =  r_arr(i) - br(ind); %vector pointing from boundary to ray previous point
            D0z =  z_arr(i) - bz(ind);
            Dr =  r_arr(i+1) - br(ind); %vector point from boundary to ray next point
            Dz =  z_arr(i+1) - bz(ind);
            %boundary tangent at nearest point
            tangent_r = br(ind+1) - br(ind) ; %tangent to boundary defined as vector from nearest boundary point to next boundary point
            tangent_z = bz(ind+1) - bz(ind) ; %used for calculating the boundary normal
            %deltas
            delta0 = (D0r)*(-tangent_z) + (D0z)*(tangent_r); %D0x scalar product with normal to boundary rotate tangent vector by pi/2
            delta = (Dr)*(-tangent_z) + (Dz)*(tangent_r); % same as previous line but for Dx
            %{%}
            if delta0<=0 %this test method ensure only reflect rays comming into boundary
                if delta>=0
                   %new step size
                   nds = -(delta0*ds)/(-delta0 + delta) ; %given equation for new step size nds
                   %nds2 = delta0/(-xi(i)*Bz + zeta(i)*Br) ; %alt method for nds; doesnt work (maybe will with updated tangent thing)

                   %new ray location on boundary
                   r_arr(i+1) = r_arr(i) + nds * c_arr(i) * xi_arr(i) ; % new ray postion so it is on boundary
                   z_arr(i+1) = z_arr(i) + nds * c_arr(i) * zeta_arr(i) ;

                   %ray_r = r_arr(i+1) - r_arr(i) ; %vector for ray
                   %ray_z = z_arr(i+1) - z_arr(i) ; % can use tangent components instead. 

                   %angle to normal calculation
                   alpha = acosd(( xi_arr(i)*(-tangent_z) + zeta_arr(i)*tangent_r )/(sqrt(xi_arr(i)^2 + zeta_arr(i)^2)*sqrt((-tangent_z)^2 + tangent_r^2))) ; %angle from ray to normal

                   %test stuff
                   disp(90-alpha) %insidence angle 
                   %disp(i+1)
                   %zeta(i) = -zeta(i); %for flat bottom comparison test

                   %gradient of tangent
                   m1 = (bz(ind+1) - bz(ind))/(br(ind+1) - br(ind)) ;
                   %disp(m1) 

                   %new tangent components
                   xi_arr(i+1) = sind(alpha + atand(m1))/c_arr(i) ; %for xi,zeta with angle reflection
                   zeta_arr(i+1) = -cosd(alpha + atand(m1))/c_arr(i); %used 180-alpha and account for bottom gradient
                    
                   yes = 1 ;
                end
            end
        end        
        
        %{%}
        %BR bottom backup
        if z_arr(i+1) > depth_lim
            r_arr(i+1) = (depth_lim*(r_arr(i+1)-r_arr(i)) + r_arr(i)*z_arr(i+1) - r_arr(i+1)*z_arr(i)) / (z_arr(i+1) - z_arr(i)) ;
            z_arr(i+1) = depth_lim ;
            zeta_arr(i+1) = -zeta_arr(i+1) ;
            %c_arr(i+1) = interp1(depth, sspd, z_arr(i+1), 'makima') ;
        end
        
        
        %for range update in ssp
        if abs(r_arr(i+1) - r(k+1)) < 1
            C(:,k) = C(:,k+1); 
            %test stuff
            %{
            %disp('k')
            %disp(k+1)
            %disp('r')
            %disp(r_arr(i+1))
            %disp(r(k+1))
            %disp('subtraction')
            %disp(abs(r_arr(i+1) - r(k+1)))
            %}
            k=k+1 ;
        end
        
        c_arr(i+1) = interp1(zq, C(:,k), z_arr(i+1), 'linear') ;      
   
    end
    disp(k)
    %plot
    plot(r_arr/1000, -z_arr/1000,".",'DisplayName',strcat(num2str(angle_arr(j)),'º')) ; hold on
    xlabel('Range (km)') ; 
    ylabel('Depth (km)') ;
    legend('show');
end
%plot(br/1000, zeros(1,n), 'b.','DisplayName', 'Surface') ; hold on 
%%
function res = F(fi, fi1, ci, ci1) 
    res = [ci*fi(3), ci*fi(4), -1/ci^2 * (ci-ci1)/(fi(1)-fi1(1)), -1/ci^2 * (ci-ci1)/(fi(2)-fi1(2))] ;
end