module odec
  double precision::umax,kp1,kp2,kp3,ka1
end module odec


program maxeff
  use dislin
  use odec
  implicit none
  integer::neq,npts,exitflag,i,j
  double precision::deltat,h,hmin,hmax,tstart,tend,eps1,eps2,ans
  double precision,allocatable,dimension(:)::y,P,A,T
  character(len=8)::fn,xlabel1,xlabel2,ylabel1,ylabel2,plottitle1,plottitle2,pictype
  integer, parameter::maxn=100
  double precision::xa,xe,xor,xstep,ya,ye,yor,ystep
  character (len=10)::fmt="(a)"
  character (len=200)::legendstring

  interface
    subroutine rkf(tstart,tend,n,y,h,hmin,hmax,eps1,eps2,f,exitflag)
      double precision,dimension(:),intent(inout)::y
      double precision,intent(in)::hmin,hmax,eps1,eps2
      double precision,intent(inout)::tstart,tend,h
      integer,intent(inout)::exitflag
      integer,intent(in)::n
      interface
        subroutine f(t,y,Dy,exitflag)
          integer,intent(inout)::exitflag
          double precision,intent(in)::t
          double precision,dimension(:),intent(in)::y
          double precision,dimension(:),intent(out)::Dy
        end subroutine f
      end interface
    end subroutine rkf
    subroutine ode(t,y,Dy,exitflag)
      integer,intent(inout)::exitflag
      double precision,intent(in)::t
      double precision,dimension(:),intent(in)::y
      double precision,dimension(:),intent(out)::Dy
    end subroutine ode
  end interface

  !This program will solve a sysyem of  first order ODE's and 
  !plot them using the dislin graphics package
  !also see lab5_dislin.f90
  !variable list
  !neq=       number of equations
  !npts=      number of output points
  !deltat=    distace between output points
  !tstart=    starting point for rkf routine
  !tend=      ending point for rkf routine
  !h=         step size for rkf routine
  !hmin=      minimum allowable step size
  !hmax=      maximum allowable step size
  !y=         solution vector at tend
  !P=         vector of phosphorus concentrations
  !A=         vector of chlorophyll-A concentrations
  !fn=        name of input file
  !exitflag=  error checking
  !eps1=      minimum allowable error
  !eps2=      maximum allowable error
 
  fn="vars.dat"
  open(11,file=fn)
  write(*,*)"How many points?"
  read(*,*)npts
  read(11,*)neq,h,hmin,hmax,eps1,eps2,umax,kp1,kp2,kp3,ka1
  allocate(y(neq),P(npts+2),A(npts+2),T(npts+2))
  deltat=200d0/dble(npts)                !fix intervals
  !deltat=60d0/dble(npts)                !for test ode
  y(1)=0.85d0
  y(2)=2d-3
  !y(1)=300                              !for test ode
  tend=0
  P(1)=y(1)
  A(1)=y(2)
  T(1)=tstart
  
  xa=0d0 ! xa is the lower limit of the x-axis.
  xe=69.11d0 ! xe is the upper limit of the x-axis.
  xor=0 ! xor is the first x-axis label.
  xstep=8.64 ! xstep is the step between x-axis labels.
  ya=0d0 ! ya is the lower limit of the y-axis.
  ye=.85d0 ! ye is the upper limit of the y-axis.
  yor=0 ! yor is the first y-axis label.
  ystep=.1d0 ! ystep is the step between y-axis labels.
  !Plot data using DISLIN
  call metafl("PS") ! or "PS", "EPS", "PDF", "WMF" "BMP"
  call setpag("USAL") !"USAL" is US size A landscape, "USAP" is portrait
  call scrmod("REVERS") !sets black on white background
  call disini() !Initialize dislin
  call complx ! Sets the font
  call name("Distance downstream (km)","X") ! Set label for x-axis
  call name("Phosphorus (mg/l)","Y") ! Set label for y-axis
  !call titlin("Phosphorus Concentrations",1) ! Set 1st line of plot title
  call psfont("Helvetica")
  call graf (xa, xe, xor, xstep, ya, ye, yor, ystep) ! sets up axis
  call title ! Actually draw the title in over the axis
  call grid(1,2)
  call legini(legendstring,7,8) ! Store 2 lines of legend text, max20 characters/line
  CALL LEGPOS(2205,0) !defines a global position for the legend where NX and NY are the 
                  !plot coordinates of the upper left corner. After a call to LEGPOS, 
                  !the second parameter in LEGEND will be ignored.
  call FRAME(5)
  call legtit("") ! set legend title (default="legend")
  call leglin(legendstring,"Po=0.850",1) ! Specify the legend text for curve 1
  call leglin(legendstring,"Po=0.735",2)
  call leglin(legendstring,"Po=0.620",3)
  call leglin(legendstring,"Po=0.505",4)
  call leglin(legendstring,"Po=0.390",5)
  call leglin(legendstring,"Po=0.275",7)
  call leglin(legendstring,"Po=0.160",6)
  
  do j=1,7
    tstart=0
    tend=0
    A(1)=y(2)
    P(1)=y(1)
    T(1)=tstart
    do i=1,npts+1
      tstart=tend                         !advance interval
      tend=tend+deltat
      call rkf(tstart,tend,neq,y,h,hmin,hmax,eps1,eps2,ode,exitflag)
      if(exitflag/=0)then
        write(*,*)"error"
        stop
      end if
      P(i+1)=y(1)
      A(i+1)=y(2)
      T(i+1)=4.32*(tend/10d0)
      if (T(i+1)>=43.2d0 .and. j==1)then
        ans=P(npts/2+1)
      end if
    end do
    call curve(T,P,npts+1) ! draw the x-y curve
    call lintyp(j) ! Change the line style (values are from 1 to 7)
    if(j<=5)then
      call lintyp(j)
    else if(j==6)then 
      call lintyp(7)
    else if(j==7)then
      call lintyp(6)
    end if
    write(*,*)"Po=",.85-dble(j-1)*.115d0
    y(1)=.85-dble(j)*.115d0
    y(2)=2d-3
    write(*,*)"T=",T(npts/2+1),j
    write(*,*)"P=",P(npts/2+1),j
   end do
  call legend(legendstring,3) ! draw legend in 7 (upper right inside axis)
  !draw legend in location 1-8. 1-4=page corner, 5-8=axis corner,1 and 5=lowerleft
  call disfin ! finish off the plot
  write(*,*)"***",ans,"***"
end program maxeff


subroutine rkf(tstart,tend,n,y,h,hmin,hmax,eps1,eps2,f,exitflag)
  implicit none
  double precision,dimension(:),intent(inout)::y
  double precision,intent(in)::hmin,hmax,eps1,eps2
  double precision,intent(inout)::tstart,tend,h
  integer,intent(inout)::exitflag
  integer,intent(in)::n
  double precision,dimension(n)::K1,K2,K3,K4,K5,K6,Dy,y4,ysave
  double precision::t,hsave,emax
  double precision,parameter::c1=1d0/5d0,c2=3d0/10d0,c3=3d0/40d0,c4=9d0/40d0,&
            c5=3d0/5d0,c6=3d0/10d0,c7=-9d0/10d0,c8=6d0/5d0,c9=11d0/54d0,&
            c10=5d0/2d0,c11=-70d0/27d0,c12=35d0/27d0,c13=7d0/8d0,&
            c14=1631d0/55296d0,c15=175d0/512d0,c16=575d0/13824d0,c17=44275d0/110592d0,&
            c18=253d0/4096d0,c19=37d0/378d0,c20=250d0/621d0,c21=125d0/594d0,&
            c22=512d0/1771d0,c23=2825d0/27648d0,c24=18575d0/48384d0,&
            c25=13525d0/55296d0,c26=277d0/14336d0,c27=1d0/4d0

  interface
    subroutine f(t,y,Dy,exitflag)
      integer,intent(out)::exitflag
      double precision,intent(in)::t
      double precision,dimension(:),intent(in)::y
      double precision,dimension(:),intent(out)::Dy
    end subroutine f
    subroutine test(t,y,Dy,exitflag)
      integer,intent(inout)::exitflag
      double precision,intent(in)::t
      double precision,dimension(:),intent(in)::y
      double precision,dimension(:),intent(out)::Dy
    end subroutine test
  end interface

  !variable list
  !y=          solution vector
  !h=          step size
  !hmin=       minimum step size
  !hmax=       maximum step size
  !hsave=      save step size for end of interval
  !t=          current time
  !tstart=     start of time interval
  !tend=       end of time interval
  !n=          number of equations
  !Kn=         Runge Kutta constants
  !Dy=         Derivative estimate at current time
  !y4=         fourth order runge-kutta 
  !emax=       error between fourth and fifth order estimates
  !exitflag=   error checking
  !eps1=       minimum allowable error
  !eps2=       maximum allowable error

  t=tstart
  exitflag=0
  if((t+h)>tend)then
    hsave=h
    h=tend-t
  end if
  do
    ysave=y
    call f(t,y,Dy,exitflag)
      if(exitflag/=0)return
      K1=h*Dy
    call f(t+c1*h,y+c1*K1*h,Dy,exitflag)
      if(exitflag/=0)return
      K2=h*Dy
    call f(t+c2*h,y+c3*K1*h+c4*K2*h,Dy,exitflag)
      if(exitflag/=0)return
      K3=h*Dy
    call f(t+c5*h,y+c6*K1*h+c7*K2*h+c8*K3*h,Dy,exitflag)
      if(exitflag/=0)return
      K4=h*Dy
    call f(t+h,y+c9*K1*h+c10*K2*h+c11*K3*h+c12*K4*h,Dy,exitflag)
      if(exitflag/=0)return
      K5=h*Dy
    call f(t+c13*h,y+c14*K1*h+c15*K2*h+c16*K3*h+c17*K4*h+c18*K5*h,Dy,exitflag)
      if(exitflag/=0)return
      K6=h*Dy
    y4=y+(c19*K1+c20*K3+c21*K4+c22*K6)*h
    y=y+(c23*K1+c24*K3+c25*K4+c26*K5+c27*K6)*h
    emax=maxval(abs((y-y4)/y))          !max relative truncation error
    if(emax>eps2 .and. abs(h-hmin)>10d-6)then
      h=h/2                             !large error,reduce step size and try again
      y=ysave
    else
      t=t+h                             !advance time, accept solution
      if(emax>eps2)exitflag=1           ! big error but h=hmin
      if(emax<eps1 .and. h<hmax)then
        h=h*2d0                         !small error increase step size
        if(h>hmax)h=hmax
      end if
      hsave=h                           !save the step size we are on
      if(t>=tend)exit                   !are we done?
      if((t+h)>tend)h=tend-t            !will the next step be beyond the end
    end if
  end do
  h=hsave
  return
end subroutine rkf

subroutine ode(t,y,Dy,exitflag)
  use odec
  implicit none
  integer,intent(inout)::exitflag
  double precision,intent(in)::t
  double precision,dimension(:),intent(in)::y
  double precision,dimension(:),intent(out)::Dy
  double precision::u
  
  u=umax*(y(1)/(kp3+y(1)))
  Dy(1)=-kp1*y(1)-kp2*u*y(2)
  Dy(2)=-ka1*y(2)+u*y(2)
  
end subroutine ode
