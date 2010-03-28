module constants
  double precision,parameter::pi=3.1415926,Q=.3,L=95,h=10,specw=998.2,&
                              visc=.000001007,g=9.81,k=.0002591
  double precision::hf
end module constants

Program pipe_dim
  use constants
  implicit none
  double precision::di,x1,x2,hp,e
  integer::tries
  logical::success
  interface
    subroutine secant(xold,xolder,maxit,epsi1,epsi2,root,numit,exitflag,f)
      double precision,intent(inout)::xold,xolder
      double precision,intent(in)::epsi1,epsi2
      double precision,intent(out)::root
      double precision::xnew,fxnew,fxold,fxolder
      integer,intent(in)::maxit
      integer,intent(out)::numit
      logical,intent(out)::exitflag
      interface
        function f(x)
          double precision::x
          double precision::f
        end function f
      end interface 
    end subroutine secant
  end interface
  interface
    function f(x)
      double precision::x
      double precision::f
    end function f
  end interface

  !This program will find the optimal pipe diameter between two water 
  !storage tanks at different heights with a pump in between as described
  !in engr 325 lab assignment 7.  To find the pipe diameter this program will call
  !a root finding subroutine which will call a function subprogram to evaluate
  !the given pipe diameter function.
  ! 
  !variable list:
  ! local variables:
  !hf      =friction loss(pump energy-change of potential) (m)
  !specw   =specific weight of water (kg/m^3)
  !Q       =flow rate(m^3/s)     
  ! inputs:
  !hp      =horsepower of pump (horsepower)
  !e       =pump efficiency
  !x1,x2   =initial pipe diameter guesses (m)
  ! outputs:
  !di      =pipe diameter (m)
  !tries   =number of iterations the rootfinding subroutine went through to determine the root
  !success =value is .true. if root successfully found,else .false.

  write(*,*)"This Program will calculate..."
  write(*,*)"Enter the horsepower of the pump"
  read(*,*)hp
  write(*,*)"Enter the efficiency of the pump"
  read(*,*)e
  hf=(76.04d0*e*hp)/(specw*Q)-h
  !write(*,*)"hf=",hf
  write(*,*)"Enter the first root estimate"
  read(*,*)x1
  write(*,*)"Enter the second root estimate"
  read(*,*)x2
  call secant(x1,x2,20,.00001d0,.00000001d0,di,tries,success,f)
  if(success)then
    write(*,"(a18,x,f10.5)")"The pipe diameter=",di
    write(*,"(a7,x,i2,x,a27)")"It took",tries,"iterations to find the root"
  end if
  stop
end program pipe_dim


subroutine secant(xold,xolder,maxit,epsi1,epsi2,root,numit,exitflag,f)
  implicit none
  double precision,intent(inout)::xold,xolder  
  double precision,intent(in)::epsi1,epsi2
  double precision,intent(out)::root
  double precision::xnew,fxnew,fxold,fxolder     !local variables
  integer,intent(in)::maxit
  integer,intent(out)::numit
  logical,intent(out)::exitflag
  interface
    function f(x)
      double precision::x
      double precision::f
    end function f
  end interface 
  
  !This is a general rootfinding subroutine that employs the secant method
  !it must be used in conjunction with an external function subprogram that
  !will evaluate the function in question
  !
  !variable list:
  ! local variables:
  !xnew      =updated root estimate
  !fxnew     =function evaluated at updated root estimate
  !fxold     =function evaluated at old root estimate
  !fxolder   =function evaluated at older root estimate
  ! inputs:
  !xold      =old root estimate
  !xolder    =older root estimate
  !epsi1     =stopping criteria for root found, if f(xnew)<epsi1 then root found 
  !epsi2     =stopping criteria for slow progress, if abs(xold-xolder)<epsi2) then too slow
  !maxit     =maxumum allowable iterations subroutine will perform until it exits
  ! outputs:
  !root      =value of particular found root
  !numit     =records number of iterations done by subroutine
  !exitflag  =value is .true. if root successfully found,else .false.

  numit=0
  fxold=f(xold)
  fxolder=f(xolder)
  do
    xnew=xold-(fxold*((xold-xolder)/(fxold-fxolder)))  
              !fxold*((xold-xolder)/(fxold-fxolder)) is root update
    fxnew=f(xnew)
    !write(*,*)"fxnew=",fxnew
    numit=numit+1
    if(abs(fxnew)<epsi1)then
      write(*,*)"root found"
      root=xnew
      exitflag=.true.
      return
      exit
    else if(numit>maxit)then
      write(*,*)"No root found (max iterations exceeded), try a better guess."
      exitflag=.false.
      return
      exit
    else if(abs(xold-xolder)<epsi2)then
      write(*,*)"No root found (progress too slow),try better guess."
      exitflag=.false.
      return
    else if ((abs(xnew-xold))>=(abs(xold-xolder)) .and. numit/=1)then
      write(*,"(a38,x,i2,x,a28)")"No root found (solution diverged after",numit,"iterations),&
                                  try better guess."
      exitflag=.false.
      return
    end if
    fxolder=fxold    !swap values to reduce number of functional evaluations
    fxold=fxnew
    xolder=xold
    xold=xnew
    !write(*,*)"fxolder=",fxolder
    !write(*,*)"fxold=",fxold
    !write(*,*)"xolder=",xolder
    !write(*,*)"xold=",xold
  end do
end subroutine secant


function f(d)
  use constants
  double precision,intent(in)::d
  double precision::A,v,fdw,Re,f    !local variables except f
  
  !This function subprogram evaluates the diameter function below
  !variable list:
  ! local variables
  !A     =cross sectional area of the pipe (m^2)
  !v     =fluid velocity (m/s)
  !fdw   =Darcy-Weisbach friction coefficient (m^2/s)
  !Re    =Reynolds number (dimensionless)
  !pi    =pi (3.1415926)
  !Q     =flow rate(m^3/s)
  !hf    =friction loss(pump energy-change of potential) (m)
  !g     =gravitational constant (9.81 m/s^2)
  !L     =pipe length (m)
  !visc  =fluid viscocity (m^2/s)
  ! inputs:
  !d     =pipe diameter (m)
  ! outputs:
  !f     =function value at given pipe diameter
  
  A=(pi*d**2d0)/4d0
  v=(4d0*Q)/(pi*d**2d0)
  fdw=(2d0*hf*g*d)/(L*v**2d0)
  Re=(v*d)/(visc)
   !write(*,*)"A=",A
   !write(*,*)"v=",v
   !write(*,*)"fdw=",fdw
   !write(*,*)"re=",Re
  f=-1d0/sqrt(fdw)-2d0*log10(k/(3.7d0*d)+2.51d0/(Re*sqrt(fdw)))
   !f=-d**(2d0)+1  !check function root=1
end function f
