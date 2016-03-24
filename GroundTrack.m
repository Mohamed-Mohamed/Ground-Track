function [ RA, DEC ] = GroundTrack( h, muo, J2, R, e, T, i0, w0, Omega0, theta0, we, dt )
% This function is used to get right ascension (longitude east of x_dash) and declination (latitude) relative to the rotating plant
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% h              : specific angular momentum vector in km^2/s
% muo         :  Gravitational Parameter
% e              : eccentricity
% a              : semimajor axis in km
% T              : orbital period in seconds
% i0              : inclination angle in degreesa
% w0            : argument of perigee in degree
% Omega0   : right ascension of the ascending node in degree
% theta0       : true anomaly in 
% J2              :Second Zonal Harmonics
% R               : raduis of plane in km
% we             : rate of rotation on the plante
%% OUTPUTS:
% RA                  : Right ascension angle in degree
% DEC                : Declination angle  in degree
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
if e > 1
    a=h^2/muo/(e^2-1)
elseif e == 0
    a=h^2/muo
elseif e > 0 && e < 1
    a=h^2/muo/(1-e^2)
end
Omega_dot=-(3/2*sqrt(muo)*J2*R^2/(1-e^2)^2/a^(7/2))*cosd(i0)
w_dot=-(3/2*sqrt(muo)*J2*R^2/(1-e^2)^2/a^(7/2))*(5/2*(sind(i0))^2-2)
% eccentric anomaly
E0=2*atand(sqrt((1-e)/(1+e))*tand(theta0/2))
% mean anomaly
M0=E0-e*sind(E0)
% initail time
t0=M0/2/pi*T
% update time
t=t0+dt 
% mean anomaly
n=floor(t/T);
Me = 2*pi/T*(t-n*T)
%


f0=Me-e/2
f=f0-e*sind(f0)-Me
f_dash=1-e*cos(f0)
sol=f0
ratio=f/f_dash
while abs(ratio) >= 1e-8
    sol=(sol-ratio);
    ratio=(sol-e*sind(sol)-Me)/(1-e*cos(sol));
end
E=sol

% 

theta=2*atand(tan(E/2)*sqrt((1+e)/(1-e)))
if theta < 0 
    theta=theta+360
end
Omega=Omega0+dt*Omega_dot
w=w0+dt*w_dot

    function [ r, v ] = OrbitalElements2rvGeo( h, muo, e, i, Omega, w, theta )
        % r @ perifocal coordinates
        r_xyz_bar=h^2/muo/(1+e*cosd(theta))*[cosd(theta);sind(theta);0];
        % v @ perifocal coordinates
        v_xyz_bar=muo/h*[-sind(theta);e+cosd(theta);0];
        % transformation matrix from perifocal to geocentric equatorial coordinates
        QxX=[-sind(Omega)*cosd(i)*sind(w)+cosd(Omega)*cosd(w),-sind(Omega)*cosd(i)*cosd(w)-cosd(Omega)*sind(w),sind(Omega)*sind(i);...
            cosd(Omega)*cosd(i)*sind(w)+sind(Omega)*cosd(w),cosd(Omega)*cosd(i)*cosd(w)-sind(Omega)*sind(w),-cosd(Omega)*sind(i);...
            sind(i)*sind(w),sind(i)*cosd(w),cosd(i)];
        % geocentric  r
        r=QxX*r_xyz_bar;
        % geocentric  v
        v=QxX*v_xyz_bar;
    end
% geocentric  r v
[ r_XYZ, v_XYZ ] = OrbitalElements2rvGeo( h, muo, e, i0, Omega, w, theta )
theta1=we*(t-t0)
R3=[cos(theta1),sin(theta1),0;-sin(theta1),cos(theta1),0;0,0,1]
% perifocal  r
r_xyz=R3*r_XYZ


    function [ ra, dec ] = lmn2RaDec( r )
        mag_r=norm(r);
        l=r(1)/mag_r;
        m=r(2)/mag_r;
        n=r(3)/mag_r;
        dec=asind(n);
        if m > 0
            ra=acosd(l/cosd(dec));
        elseif m <= 0
            ra=360-acosd(l/cosd(dec));
        end
    end
[ RA, DEC ] = lmn2RaDec( r_xyz )        
end

