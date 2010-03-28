program pde1
  use dislin
  implicit none
  double precision,allocatable,dimension(:,:)::P,Q
  double precision,allocatable,dimension(:)::x
  double precision::ub,lb,eps,T,w,psave,prate
  integer::nxgrd,prow,pcol,i,j,numit,maxit,exit,npump

  !Variable list
  !P           =Hydraulic head matrix
  !Q           =pumping rate matrix
  !x           =x and y array for shaded surface
  !ub,lb       =upper and lower boundry conditions
  !eps         =convergence tolerance
  !T           =transmisivity
  !psave       =saves value of current node for convergence check
  !prate       =pumping rate
  !nxgrd       =number of nodes in x and y direction nxgrd^2 total nodes
  !prow,pcol   =pump row and column number
  !numit       =number of iterations
  !maxit       =maximum number of iterations
  !exit        =convercence criteria exit==nxgrd**2-2*nxgrd
  !npump       =number of pumps

  open(11,file="vars.dat")
  open(12,file="head.out")
  open(13,file="pump.dat")
  read(11,*)nxgrd,ub,lb,eps,T,w,maxit
  allocate(P(nxgrd,nxgrd),Q(nxgrd,nxgrd),x(nxgrd))
  do i=1,nxgrd
    x(i)=x(i-1)+100d0
  end do
  Q=0
  read(13,*)npump
  do i=1,npump
    read(13,*)prow,pcol,prate
    !write(*,*)prow,pcol,prate
    Q(prow,pcol)=prate
  end do 
  P=0
  P(1,1:nxgrd)=ub
  P(nxgrd,1:nxgrd)=lb
  numit=0
  do
    numit=numit+1
    exit=0
    do i=2,nxgrd-1                !dont include
      do j=1,nxgrd
        psave=P(i,j)
        if(j==1)then             !left edge no flow
          P(i,1)=P(i,1)+(w/4d0)*(P(i+1,1)+P(i-1,1)+2d0*P(i,2)-4d0*P(i,1))
        else if(j==nxgrd)then   !right edge no flow
          P(i,nxgrd)=P(i,nxgrd)+(w/4d0)*(P(i+1,nxgrd)+P(i-1,nxgrd)+2d0*P(i,nxgrd-1)-4d0*P(i,j))
        else                       !pump site
          P(i,j)=P(i,j)+(w/4d0)*(P(i+1,j)+P(i-1,j)+P(i,j-1)+P(i,j+1)-4d0*P(i,j)+Q(i,j)/T)
        end if
        if (abs(psave-P(i,j))<=eps)then
          exit=exit+1            !node isnt changing
        end if
      end do
    end do
    if(numit>=maxit)then
      write(*,*)(nxgrd**2-2*nxgrd),exit,numit,maxit,"maxit"
      do i=1,nxgrd
        write(12,"(10000f10.5)")(P(i,j),j=1,nxgrd)
      end do
      stop
    else if(exit>=(nxgrd**2-2*nxgrd))then
      do i=1,nxgrd
        write(12,"(10000f10.5)")(P(i,j),j=1,nxgrd)
      end do
      write(*,*)numit
      write(*,*)"This is the number",P(13,15)
      call setpag("USAP") !"USAL" is US size A landscape, "USAP" is portrait
      call scrmod("REVERS") !sets black on white background
      call metafl("xwin") ! or "PS", "EPS", "PDF", "WMF" "BMP"
      call disini() !Initialize dislin 
      call disalf
      call psfont("helvetica")
      call name("Distance (m)","XY") ! Set label for x-axis
      !call name("Phosphorus (mg/l)","Y") ! Set label for y-axis
      call name("Head (m)","Z") ! Set label for z-axis
      call axis3d(2d0,2d0,2d0)           !set up axis in absolute plot coordianates
      call view3d(5d0,-4.5d0,2d0,"abs")  !set view point in absolute coordinates
      call vfoc3d(0d0,0d0,0d0,"abs")     !set location to look at in absolute coordiates
      
      !these two for surface mesh
      call labl3d("HORIZONTAL")  ! 'STANDARD', 'HORIZONTAL', 'PARALLEL' and 'OTHER'.
      call graf3d(0d0,3000d0,0d0,600d0,0d0,3000d0,0d0,600d0,45d0,50d0,45d0,1d0)
      call surmat(P,nxgrd,nxgrd,2,2)
  
      !these for contour plot
      !!!!!call nobar
      !call axspos(320,2400)
      !call labtyp("vert","X")
      !CALL NAMDIS (5,"Y")
      !CALL NAMDIS (0,"z")
      !call graf3(3000d0,0d0,3000d0,-200d0,0d0,3000d0,0d0,200d0,50d0,43.8d0,50d0,-1d0)   
                  !set up axis: label range, first label,step betwen labels, X-Y-Z
      !!!!call zaxis(45d0,50d0,45d0,1d0,1800,"Head",0,1,480,200)   !custop color bAR
      !call colran(1,254)                      !color range
      !call crvmat(transpose(P),nxgrd,nxgrd,20,30)        !draw 2d color contour plot
      
      !these two for color surface plot
      !call labl3d("HORIZONTAL")  ! 'STANDARD', 'HORIZONTAL', 'PARALLEL' and 'OTHER'.
      !call shdmod("smooth","surface") !plot a "smooth" or "flat" interpolated graph
      !call surshd(x,nxgrd,x,nxgrd,P)  !plot the surface

      call title     !draw on title
      call disfin
      stop
    end if
  end do
  rewind(12)
  rewind(11)
  stop
end program pde1
