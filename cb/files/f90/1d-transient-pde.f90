program pde2
  use dislin
  implicit none
  double precision,dimension(:),allocatable::a,b,c,d,x,t
  double precision::tstep,xstep,tini,phi,th,k,thk,ccon,p,h,Ts,time
  double precision::xa,xe,xor,xastep,ya,ye,yor,ystep
  integer::nxgrd,ntgrd,i,j
  character(len=6)::legendstring

  !Variable list

  interface
    subroutine tomas(a,b,c,d,x,neq)
      double precision,dimension(:),intent(out)::x
      double precision,dimension(:),intent(in)::a,b,c,d
      integer,intent(in)::neq
    end subroutine
  end interface


  open(11,file="prm.dat")
  open(12,file="heat.out")
  read(11,*)ntgrd,nxgrd,tini,thk,k,ccon,p,h,Ts,time
  write(*,*)nxgrd

  tstep=time/dble(ntgrd)
  write(*,*)tstep
  xstep=thk/nxgrd
  phi=(k*tstep)/(ccon*p*xstep**2)
  th=(h*tstep)/(ccon*p*xstep)

  allocate(a(nxgrd),b(nxgrd),c(nxgrd),d(nxgrd),x(nxgrd),t(nxgrd))
  write(*,*)"here"
  d(1:nxgrd-1)=30d0
  write(12,"(100f10.5)")(d(j),j=1,nxgrd)
  d(nxgrd)=30d0+2d0*th*Ts
  a(1)=0
  a(2:nxgrd-1)=-phi
  a(nxgrd)=-2d0*phi
  b(1:nxgrd-1)=1+2d0*phi
  b(nxgrd)=1d0+2d0*phi+2d0*th
  c(1)=-2d0*phi
  c(2:nxgrd-1)=-phi
  c(nxgrd)=0

  t(1)=0d0
  do i=2,nxgrd
    t(i)=t(i-1)+xstep
  end do
  
  xa=0d0 ! xa is the lower limit of the x-axis.
  xe=30d0 ! xe is the upper limit of the x-axis.
  xor=0 ! xor is the first x-axis label.
  xastep=2d0 ! xstep is the step between x-axis labels.
  ya=23d0 ! ya is the lower limit of the y-axis.
  ye=30.01d0 ! ye is the upper limit of the y-axis.
  yor=23d0 ! yor is the first y-axis label.
  ystep=1d0 ! ystep is the step between y-axis labels.
  !Plot data using DISLIN
  call metafl("xwin") ! or "PS", "EPS", "PDF", "WMF" "BMP"
  call setpag("USAL") !"USAL" is US size A landscape, "USAP" is portrait
  call scrmod("REVERS") !sets black on white background
  call disini() !Initialize dislin
  call complx ! Sets the font
  call name("Distance From Back of Slab (cm)","X") ! Set label for x-axis
  call name("Temperature (degrees C)","Y") ! Set label for y-axis
  !call titlin("Phosphorus Concentrations",1) ! Set 1st line of plot title
  call psfont("Helvetica")
  call graf (xa, xe, xor, xastep, ya, ye, yor, ystep) ! sets up axis
  call title ! Actually draw the title in over the axis
  !call grid(1,2)
  call legini(legendstring,12,8) ! Store 2 lines of legend text, max20 characters/line
  !CALL LEGPOS(2405,500) !defines a global position for the legend where NX and NY are the
                  !plot coordinates of the upper left corner. After a call to LEGPOS,
                  !the second parameter in LEGEND will be ignored.
  call FRAME(3)
  call legtit("") ! set legend title (default="legend")
  call leglin(legendstring,"1 hr.",1) ! Specify the legend text for curve 1
  call leglin(legendstring,"2 hr.",2)
  call leglin(legendstring,"3 hr.",3)
  call leglin(legendstring,"4 hr.",4)
  call leglin(legendstring,"5 hr.",5)
  call leglin(legendstring,"6 hr.",6)
  call leglin(legendstring,"7 hr.",7)
  call leglin(legendstring,"8 hr.",8)
  call leglin(legendstring,"9 hr.",9)
  call leglin(legendstring,"10 hr.",10)
  call leglin(legendstring,"11 hr.",11)
  call leglin(legendstring,"12 hr.",12)
  
  do i=1,ntgrd
    call tomas(a,b,c,d,x,nxgrd)
    write(12,"(100f10.5)")(x(j),j=1,nxgrd)
    write(12,*)""
    call curve(t,x,nxgrd) ! draw the x-y curve
    if(i>=8)then
      call lintyp(i-7) ! Change the line style (values are from 1 to 7)
    else
      call lintyp(i)
    end if
    d=x
    d(nxgrd)=x(nxgrd)+2d0*th*Ts
  end do

  call legend(legendstring,5) ! draw legend in 7 (upper right inside axis)
  !draw legend in location 1-8. 1-4=page corner, 5-8=axis corner,1 and 5=lowerleft
  call disfin ! finish off the plot
  
  rewind(12)
  rewind(11)
  stop
end program pde2


subroutine tomas(a,b,c,d,x,neq)
  implicit none
  double precision,dimension(:),intent(out)::x
  double precision,dimension(:),intent(in)::a,b,c,d
  integer,intent(in)::neq
  double precision,dimension(neq)::e,f
  integer::i,j

  e(1)=d(1)/b(1)
  f(1)=c(1)/b(1)

  do i=2,neq
    e(i)=(d(i)-a(i)*e(i-1))/(b(i)-a(i)*f(i-1))
    f(i)=c(i)/(b(i)-a(i)*f(i-1))
  end do

  x(neq)=e(neq)

  do i=neq-1,1,-1
    x(i)=e(i)-f(i)*x(i+1)
  end do

end subroutine
