module con
  double precision,parameter::e=2.718281828,pi=3.1415926,m=45000,D=.086,k=.153,Ts=250
end module con

program heatex
  use con
  implicit none
  integer::numtrap,maxit,numit,s,j
  double precision::T1,T2,ans,eps
  character(len=30)::func
  double precision,allocatable,dimension(:,:)::I
  interface
    subroutine romtrap(nt,a,b,ans,maxit,eps,numit,IG,f)
      implicit none
      integer,intent(in)::nt,maxit
      integer,intent(out)::numit                
      double precision,intent(in)::b,a,eps
      double precision,intent(out)::ans
      double precision,dimension(:,:),intent(out)::IG
      interface
        function f(x)
          double precision::x
          double precision::f
        end function f
      end interface
    end subroutine romtrap
    function temp(x)
      double precision::x
      double precision::temp
    end function temp
  end interface

    !this program will integrate a pre set flow equation using an integration 
    !subroutine. The subroutine uses the trapezoid rule given upper and lower 
    !limits, a number of trapezoids and a pipe radius.
    !variable list:
    !T1=lower limit of integration
    !T2=upper limit of integration
    !numtrap=number of trapezoids

  write(*,*)"This program will integrate an equation ",&
            "using a trapezoid rule integration subroutine(answer in feet)."
  write(*,*)"Enter the lower limit of integration."
  read(*,*)T1
  write(*,*)"Enter the upper limit of integration." 
  read(*,*)T2
  numtrap=1
  !write(*,*)"Enter the number of trapezoids to use."
  !read(*,*)numtrap
  maxit=20
  allocate(I(maxit,maxit))
  eps=.0000000001
  call romtrap(numtrap,T1,T2,ans,maxit,eps,numit,I,temp)
  write(*,"(/,a,x,f14.10)")"The best estimate is:",ans
  write(*,"(5x,a,i4,/)")"# of iterations:",numit+1
  write(*,"(a)")"The Romberg extrapolation table is (first col. is stepsize):"
  write(*,"(7x,a,13x)",advance="no")"h"
  do j=1,numit+1
    write(*,"(a,i2,11x)",advance="no")"I",j
  end do 
  write(*,"(a)",advance="yes")" "
  do s=1,numit+1
      write(*,"(100f14.10)")(I(s,j),j=1,numit+2)
  end do
  stop
end program heatex

subroutine romtrap(nt,a,b,ans,maxit,eps,numit,IG,f)
  implicit none
  integer,intent(inout)::nt,maxit
  integer,intent(out)::numit                
  double precision,intent(in)::b,a,eps
  double precision,intent(out)::ans
  double precision,dimension(:,:),intent(out)::IG
  double precision::h,intp,endpts,newp,r,tol        !local variables
  integer::i,j,k               
  
    !variable list:
    !h=        step size
    !a,b=      upper/lower limit of integration
    !interior= area under the curve without the endpoints
    !f=        function name
    
  interface
    function f(x)
      double precision::x
      double precision::f
    end function f
  end interface 
  
  h=(b-a)/nt
    !write(*,*)h
  IG(1,1)=h
  intp=0
  do i=2,nt
    intp=intp+f(a+(i-1)*h)
  end do
  endpts=f(b)/(2d0)+f(a)/(2d0)
  IG(1,2)=(h)*(endpts+intp)
  
  r=2d0
  j=0
  do 
    j=j+1
    h=h/r
    nt=2d0*nt
    IG(j+1,1)=h
    do i=1,nt/2
      newp=newp+f(a+(i-1)*(2d0*h)+h)
    end do
    IG(j+1,2)=(h)*(endpts+intp+newp)
    k=1
    do i=1,j
      k=2*k
      IG(j+1,i+2)=((r**(k)*IG(j+1,i+1))-IG(j,i+1))/(r**(k)-1)
    end do
    tol=abs(IG(j+1,j+2)-IG(j+1,j+1))
    if(j>maxit-2)then
      write(*,*)"Reached max iteration"
      numit=j
      return
    else if(tol<eps)then
      ans=IG(j+1,j+2)
      numit=j
      exit
      return
    end if
  end do
  
  return
end subroutine romtrap

function temp(T)
  use con
  double precision::temp,T,cp,lnu,u,ht
  cp=.53d0+.00065*T
  lnu=.00002d0*(T*T)-((.0225d0)*T)+5.4874d0
  u=e**lnu
  ht=((.023d0*k)/d)*(((4d0*m)/(pi*D*u))**(.8d0))*(((u*cp)/k)**(.4d0))
  temp=(m/(pi*D))*(cp/(ht*((Ts-T)**(2d0))))
  return
end function temp

