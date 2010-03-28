program cmts
  implicit none
  integer::ntanks,itr,maxit,i,j
  double precision::Q,Qr,V,k,cin,w,offdiagsum,eps,check,done
  double precision,allocatable,dimension(:,:)::A,AB
  double precision,allocatable,dimension(:)::b,x
  logical::diagdom,singular
  character(len=1)::ans

  interface
    subroutine gsswap(a,neqn,diagdom,singular)
      implicit none
      integer, intent(in)::neqn
      double precision, dimension(:,:), intent(inout)::a
      logical, intent(out)::diagdom, singular
    end subroutine gsswap
  end interface

  !This program will compute the effluent concentration values of 
  !a pond using the complete-mix tank-in series approach using
  !the Gauss-Seidell method to solve for the concentration values.
  !~~variable list~~
  !--input variables--
  !ntanks     =number of tanks to use as an estimate
  !Q          =influent flowrate (m^3/day)
  !Qr         =recycle flow rate (m^3/day)
  !V          =tank or compartment volume
  !k          =first order reaction rate constant (day^-1)
  !cin        =influent polution concentration
  !w          =relaxation coefficent
  !eps        =solution found tolerance
  !maxit      =maximum number of iterations
  !x          =solution vector needs initial guess
  !--output variables--
  !x          =solution vector of final concentration values
  !itr        =number of iteration carried out
  !--internal variables--
  !i,j        =loop vaiables
  !A          =holds initial coefficient matrix
  !b          =initial right hand side vector
  !ans        =answer to "do you want to continue?"

  write(*,*)"This program will compute the effluent concentration values"
  write(*,*)"of a pond using the complete-mix tank-in series approach using"
  write(*,*)"the Gauss-Seidell method to solve for the concentration values."
  write(*,*)" "
  write(*,*)"How many tanks?"
  read(*,*)ntanks
  if(ntanks<=0)then
    write(*,"(a21,i2,a6)")"Can't do operation on",ntanks," tanks"
    stop
  end if
  allocate(A(ntanks,ntanks),x(ntanks),b(ntanks),AB(ntanks,ntanks+1))
  A=0
  b=0
  Q=10000
  Qr=.2*Q
  V=50000/ntanks
  k=.093
  cin=30
  write(*,*)"Enter the relaxation coefficient 0<w<2"
  read(*,*)w
  eps=.00001
  maxit=4999
  do i=1,ntanks
    write(*,"(a28,i3)")"Enter concentration for tank",i
    read(*,*)x(i)
  end do
  A(1,1)=-((Q+Qr)/V+k)
  A(1,2)=Qr/V
  b(1)=(-Q*cin)/V
  if(ntanks>1)then
    do i=2,ntanks-1
      A(i,i-1)=(Q+Qr)/V
      A(i,i)=-((Q+2*Qr)/V+k)
      A(i,i+1)=Qr/V
    end do
    A(ntanks,ntanks-1)=(Q+Qr)/V
    A(ntanks,ntanks)=-((Q+Qr)/V+k)
  end if
  AB=A                   !create augmented matrix instead of modifying originals
  AB(1:ntanks,ntanks+1)=b(1:ntanks)   !augment A with b
  write(*,*)"[A|b]="
  do i=1,ntanks
    write(*,"(a)",advance="no")"|"
    write(*,"(1000f8.3)",advance="no")(AB(i,j),j=1,ntanks+1)
    write(*,"(a)")"|"
  end do
  
  do i=1,ntanks
    offdiagsum=0
    offdiagsum=sum(abs(A(i,1:ntanks)))-abs(A(i,i))
    if(offdiagsum>abs(A(i,i)))then
      write(*,*)"Answer may not converge."
      write(*,*)"Continue? (y,n)"
      read(*,*)ans
      do
        if(ans=="n")then
          stop
        else if(ans/="y")then
          write(*,*)"not an option"
        else if(ans=="y")then
          exit
        end if
      end do
    end if
  end do
  
  call gsswap(AB,ntanks,diagdom,singular)
  !do i=1,ntanks
  !  write(*,"(a)",advance="no")"|"
  !  write(*,"(1000f8.3)",advance="no")(AB(i,j),j=1,ntanks+1)
  !  write(*,"(a)")"|"
  !end do
  if(singular)then
    write(*,*)"The coefficient matrix is singular."
    stop
  end if
  if(.not.diagdom)then
    write(*,*)"The augmented matrix is not diagonally dominant, may not converge."
    write(*,*)"Continue? (y,n)"
    read(*,*)ans
    do
      if(ans=="n")then
        stop
      else if(ans/="y")then
        write(*,*)"not an option"
      else if(ans=="y")then
        exit
      end if
    end do
  end if
  
  !write(*,"(a11,f8.5)")"offdiagsum=",offdiagsum
  itr=0
  do
    itr=itr+1
    done=0
    do i=1,ntanks
      check=0
      check=x(i)
      x(i)=x(i)+w*((AB(i,ntanks+1)-sum(AB(i,1:i-1)*x(1:i-1))-sum(AB(i,i:ntanks)*x(i:ntanks)))/AB(i,i))
      if(abs(check-x(i))<eps)then
        done=done+1
      end if
    end do
       !write(*,*)"x="
       !do i=1,ntanks
       !  write(*,"(f10.5)")x(i)
       !end do
    if(done==ntanks)then
      write(*,*)" "
      write(*,*)"Solution found"
      write(*,*)"x="
      do i=1,ntanks
        write(*,"(f10.5)")x(i)
      end do
      exit
    else if(itr>maxit)then
      write(*,*)"No solution,max iterations exceeded"
      exit
    end if
  end do
  write(*,*)" "
  write(*,"(a12,i2,a9,f9.5)")"average for ",ntanks," tanks = ",sum(x(1:ntanks))/ntanks
  write(*,"(a5,i4,a12)")"took ",itr," iterations."
  write(*,*)" "
  stop
end program
