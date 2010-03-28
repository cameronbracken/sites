

function report_copy(nxgrd,ntgrd,az,el,a,b,map)

%damped_vibrating_membrane(nxgrd,ntgrd,az,el,a,b,map)
%damped_vibrating_membrane(20,500,60,20,.01,10,'jet')
%solve the wave equation on a rectangle(explicit finite difference scheme)
%Author Cameron Bracken
%
%nxgrd   =number of x grid points
%ntgrd   = how long to run for
%deltat  = time step (set to max allowable by stability conditions)
%a       = magnitude of impact (not a real value just guess ~0.01)
%b       = damping coefficient 0:no damping 1:very floppy
%                             10:moderately damped 100:VERY damped
%az,el   = azumuth and elevation for viewing movie
%u(x,y,t)= 3D array displacement of membrane size(u)=nxgrd x nygrd x ntgrd
%i,j,n   = x,y,time grid increments respectively
%s,r     = random numbers used to hit the drum in random places
%p,alpha = coefficeints from derivation of algorithm
%hit     = array of timesteps at which membrane is hit
%map     = colormap options include: jet,summer,hot,[0 0 0] is black
%
%note:unstabe if c*deltat/deltax<=1/sqrt(2)
%deltax=1/(nxgrd-1)


deltax=1/(nxgrd-1);                  %calculate x grid spacing
nygrd=nxgrd;                         %number of gridpoints sam in x and y
deltay=deltax;                       %y grid spacing same as x spacing (square)
x=(0:deltax:1);                      %dimensionless grid
y=(0:deltay:1);

c=10;                                %this is a total guess of wave number, 
                                     %program blows up if c>50
deltat=(1/sqrt(2))*(deltax/c)-.001;  %%set timestep to max allowed by stabiliy
p=c*deltat/deltax;                   %   
alpha=2-4*p^2-deltat*b/2;            %coefficient from
%1/sqrt(2)



u=zeros(nxgrd,nygrd,ntgrd);          %initial conditions zero everywhere

[X,Y]=meshgrid(x,y);                 %create mesh points

%aviobj = avifile('dampedwave.avi'); %create avi file to make movie 

hit=[0,2,10:10:340];                   %timesteps at which membrane is hit
k=2;

 for n=2:(ntgrd-1)                   %time loop
     if (n==hit(k))
                  if or(k+1>size(hit),k==1)
                      k=1;
                  else
                      k=k+1;
                  end 
                  s=round((nxgrd-5)*rand(1,1)+3);  %generate random spot to hit away from boundary
                  r=round((nxgrd-5)*rand(1,1)+3);
                  u(s-2:s+2,r-2:r+2,n-1)=u(s-2:s+2,r-2:r+2,n-1)-a;  %Hit that shit
     end
     for i=2:(nxgrd-1)              %x direction loop
         for j=2:(nygrd-1)          %y direction loop
            u(i,j,n+1)=(2/(2+deltat*b))*(alpha*u(i,j,n)+p^2*(u(i-1,j,n)...
                       +u(i+1,j,n)+u(i,j-1,n)+u(i,j+1,n))-u(i,j,n-1));
         end
     end
        
     surf(X,Y,u(:,:,n))             %plot membrane at evey timestep
     colormap(map)                  %set plot color
     if n==2
          climits = [-max(max(u(:,:,3))), max(max(u(:,:,3))) ];
          caxis(climits);
     end 
     hold off;                      %dont plot on top of old mesh
     axis([0 1 0 1 -.5 .5])         %set axis
     view(az,el)                    %set view angle and elevation
     %frame=getframe(gcf);              %grab the current frame as a picture
     %aviobj = addframe(aviobj,frame);  %put the current frame in the avi file 
     if n==2
         pause(1.5)                 %stop at the begining 
     end
     pause(.05)                     %pause in between timesteps to see plot
 end
 
 %aviobj = close(aviobj);               %finish making the movie file
        