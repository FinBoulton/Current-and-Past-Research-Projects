% Euler method for non-flat bottom
%% data for ssp over range

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
figure
plot(P1,-z/1000); xlabel('Sound Speed (m/s)'); ylabel('Depth (km)')
meanSSP = mean(P1,2);
hold on
plot(meanSSP,-z/1000,'kx','LineWidth',2)
hold off


%plot of variation in SSP along path
%figure
%surf(r,-z/1000,P1)
%xlabel('Range')
%ylabel('Depth (km)'); title('Sound speed over range')



figure

C = zeros(5501,500) ;
for i = 1: rdn
    zq = 0:1:5500 ;
    c = interp1(z, P1(:,i), zq, 'linear') ;
    C(:,i) = c ;
    plot(c,-zq/1000)
    hold on
end

%data variables 
depth_lim = 5500.0 ;




%% Euler method

%bottom fn
syms x
y = piecewise((0<=x)&(x<25), -5.5,(25<=x)&(x<50), 1/10*x - 8,(50<=x)&(x<75), -1/10*x + 2, (75<=x),-5.5  ) ;
figure
fplot(y,'k')
hold on

%IC
r0 = 0.0 ; z0 = 0.0 ; %Source point

k=1;

%Solution interval and discretization step
ds = 1.0 ;
n=140000;

%for many anglep
angle_arr = [1.0, 13.0, 16.0] ;


for j = 1:length(angle_arr)
    
    k=1;
    
    %aux varibles IC
    xi0   = cosd(angle_arr(j)) ./ C(1,1) ;
    zeta0 = sind(angle_arr(j)) ./ C(1,1) ;

    %Create array that will contain the numerical solution
    r_tst = zeros(1, n) ; %range
    z_tst = zeros(1, n) ; %depth
    s     = zeros(1, n) ; %arclength parameter
    xi    = zeros(1, n) ; %aux
    zeta  = zeros(1, n) ;
    c_arr = zeros(1, n) ; %sound speed

    %Set the initial values
    r_tst(1) = r0 ;
    z_tst(1) = z0 ;
    s(1)     = ds ;
    xi(1)    = xi0 ;
    zeta(1)  = zeta0 ;
    c_arr(1) = C(1,1) ;
        
        % Perform the iteration
        for i = 1:n-1
            r_tst(i+1) = r_tst(i) + ds * c_arr(i) * xi(i) ;
            z_tst(i+1) = z_tst(i) + ds * c_arr(i) * zeta(i) ;
            
            %BR surf
            if z_tst(i+1) < 0.0
                z_tst(i+1) = 0.0;
                r_tst(i+1) = r_tst(i) - (z_tst(i) ./ zeta(i)) * xi(i); 
                zeta(i) = -zeta(i);
                %disp('Euler')
                %disp(r_tst(i+1))
            end
            
            %BR seamount rise
            if (25000 < r_tst(i+1)) && (r_tst(i+1) < 50000)
                ray = [r_tst(i), z_tst(i)];
                rayy = [r_tst(i+1), z_tst(i+1)];
                d0 = abs(ray - [25000,5500]) ;
                d = abs(rayy - [25000,5500]) ;
                norm = [37375,3000] ;
                delta0 = dot(d0,norm) ;
                delta = dot(d,norm) ;
                if delta0 < 0
                    if (delta > 0)
                        ds = ((-delta0)/(-delta0 - delta))*ds ;
                        r_tst(i+1) = r_tst(i) + ds * c_arr(i) * xi(i) ;
                        z_tst(i+1) = z_tst(i) + ds * c_arr(i) * zeta(i) ;
                        zeta(i) = -zeta(i) ;
                        disp('yes')
                    end
                end
            end
            
            %BR seamount fall
            if (50000 < r_tst(i+1)) && (r_tst(i+1) < 75000)
                ray = [r_tst(i), z_tst(i)];
                rayy = [r_tst(i+1), z_tst(i+1)];
                d0 = abs(ray - [75000,5500]) ;
                d = abs(rayy - [75000,5500]) ;
                norm = [61775,3000] ;
                delta0 = dot(d0,norm) ;
                delta = dot(d,norm) ;
                if delta0 < 0
                    if (delta > 0)
                        ds = ((-delta0)/(-delta0 - delta))*ds ;
                        r_tst(i+1) = r_tst(i) + ds * c_arr(i) * xi(i) ;
                        z_tst(i+1) = z_tst(i) + ds * c_arr(i) * zeta(i) ;
                        zeta(i) = -zeta(i) ;
                        disp('yes')
                    end
                end
            end
            
            %BR bottom outside of seamount
            if z_tst(i+1) > depth_lim
                if r_tst(i+1) < 25000
                    z_tst(i+1) = depth_lim;
                    r_tst(i+1) = r_tst(i) + ((depth_lim - z_tst(i)) ./ zeta(i)) * xi(i);
                    zeta(i) = -zeta(i);
                    %disp(r_tst(i+1))
                    %disp(z_tst(i+1))
                end
                if r_tst(i+1) > 75000
                    z_tst(i+1) = depth_lim;
                    r_tst(i+1) = r_tst(i) + ((depth_lim - z_tst(i)) ./ zeta(i)) * xi(i);
                    zeta(i) = -zeta(i);
                    %disp(r_tst(i+1))
                    %disp(z_tst(i+1))
                end
            end
            
            %{
            %BR bottom norm implimentation
            if z_tst(i+1) > depth_lim
                z_tst(i+1) = depth_lim;
                r_tst(i+1) = r_tst(i) + ((depth_lim - z_tst(i)) ./ zeta(i)) * xi(i);
                zeta(i) = -zeta(i);
                %disp(r_tst(i+1))
                %disp(z_tst(i+1))
            end
            %}
            
            %for range update in ssp
            if mod(round(r_tst(i+1)),round(r(2))) == 0
                C(:,k) = C(:,k+1);
                k=k+1 ;
                %if k == rdn
                %    k=1 ;
                %end
                %disp(round(r_tst(i+1)))
                %disp(k)
            end
            
            c_arr(i+1) = interp1(zq, C(:,k), z_tst(i+1), 'linear') ;

            xi(i+1)    = xi(i)    - ds * (1./c_arr(i)^2 * (c_arr(i+1)-c_arr(i))./(r_tst(i+1)-r_tst(i))) ;
            zeta(i+1)  = zeta(i)  - ds * (1./c_arr(i)^2 * (c_arr(i+1)-c_arr(i))./(z_tst(i+1)-z_tst(i))) ;
            %s(i+1)     = (i+1) * ds ;
        end
 
    %plot
    plot(r_tst/1000, -z_tst/1000) ; hold on
    xlabel('Range (km)') ; 
    ylabel('Depth (km)') ;
    legend('Bottom','1º','13º','16º');
    
end
%plot([0, n/1000], [0, 0], 'k-') ;
