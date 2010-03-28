program pipe_network
  implicit none
  double precision,allocatable,dimension(:,:)::A,b,x
  double precision,allocatable,dimension(:,:)::rowef
  double precision::det,b01,b12,b13,b23,b24,b34,b45,p0,tol
  integer::neqn,nrhs,i,j,check
  character(len=1)::ans
  character(len=30)::file
  logical::scale,solution,pivot

  interface
    subroutine gauss(A,b,x,ref,neqn,nrhs,scale,solution,det,pivot,tol)
      double precision,dimension(:,:),intent(in)::A,b
      double precision,dimension(:,:),intent(out)::x
      double precision,dimension(neqn,neqn+nrhs),intent(out)::ref
      double precision,intent(in)::tol
      double precision,intent(out)::det
      integer,intent(in)::neqn,nrhs
      logical,intent(in)::scale,pivot
      logical,intent(out)::solution
    end subroutine gauss
  end interface

  !This Program will read in a square coefficient matrix A and a right 
  !hand side vector b and call a gauss elimination subroutine to find 
  !the solution vector x.
  !~~It can handle b and x vector with rank>1~~
  !variable list:
  !--inputs:
  !A       =coefficient matrix(pressure values in this case)
  !         must be square
  !b       =right hand side vector
  !neqn    =number of equations/unknowns in matrix A (number of rows)
  !nrhs    =number of columns in right hand side vector b
  !ans     =answer to "do you want to scale"
  !file    =name of file containing data
  !pivot   =if ans is yes then pivot=.true. else pivot=.false.
  !tol     =singularity tolerance for diagonal elements
  !--outputs:
  !x       =solution vector
  !det     =determinant of matrix A
  !ref     =row eschelon form of marrix A|b
  !--internal variables:
  !i,j     =loop variables
  !check   =iostat check for read ans
  !scale   =if ans is yes then scale=.true. else scale=.false.
  !solution=if the system has a solution then solution=.true.

  write(*,*)"This program will read in a square coefficient matrix A"
  write(*,*)"of pressure values and a right hand side vector b and call a"
  write(*,*)"gauss elimination subroutine to find the solution vector x."
  write(*,*)" "
  write(*,*)"Enter the name of the file that contains the initial pressure"
  write(*,*)"and friction coefficients of matrix A and vector b."
  read(*,"(a)")file
  open(11,file=file)
  neqn=4
  nrhs=1
  allocate(A(neqn,neqn),b(neqn,nrhs),x(neqn,nrhs),rowef(neqn,neqn+nrhs))
  do
    read(11,*,iostat=check)p0,b01,b12,b13,b23,b24,b34,b45
    if(check/=0)then
      write(*,*)"something is wrong, check your file and try again."
      exit
      stop
    else if(check==0)then
      exit
    end if
  end do
  
  A(1,1)=(b01+b12+b13)
  A(1,2)=-b12
  A(1,3)=-b13
  A(1,4)=0
  A(2,1)=-b12
  A(2,2)=(b12+b23+b24)
  A(2,3)=-b23
  A(2,4)=-b24
  A(3,1)=-b13
  A(3,2)=-b23
  A(3,3)=(b13+b23+b34)
  A(3,4)=-b34
  A(4,1)=0
  A(4,2)=-b24
  A(4,3)=-b34
  A(4,4)=(b24+b34+b45)
  b(1,1)=b01*p0
  b(2,1)=0
  b(3,1)=0
  b(4,1)=0

  write(*,*)"Enter the matrix singularity tolerance"
  read(*,*)tol
  do
    write(*,*)"Do you want to scale matrix A|B? (y or n)"
    read(*,"(a)",iostat=check)ans
    if(ans=="y" .and. check==0)then
      scale=.true.
      exit
    else if(ans=="n")then
      scale=.false.
      exit
    else if(check/=0)then
      write(*,*)"Not an option."
    end if
  end do
  do
    write(*,*)"Do you want to use pivoting? (y or n)"
    read(*,"(a)",iostat=check)ans
    if(ans=="y" .and. check==0)then
      pivot=.true.
      exit
    else if(ans=="n")then
      pivot=.false.
      exit
    else if(check/=0)then
      write(*,*)"Not an option."
    end if
  end do
  
  call gauss(A,b,x,rowef,neqn,nrhs,scale,solution,det,pivot,tol)

  if(solution)then
    write(*,*)"ref(A|b)="
    do i=1,neqn
      write(*,"(a)",advance="no")"|"
      write(*,"(100f8.3)",advance="no")(rowef(i,j),j=1,nrhs+neqn)
      write(*,"(a)")"|"
    end do
    write(*,*)" "
    write(*,*)"The solution vector:"
    write(*,*)"x="
    do i=1,neqn
      write(*,"(100f10.3)")(x(i,j),j=1,nrhs)
    end do
    write(*,"(/,a8,f8.5,/)")" det(A)=",det
  end if
  stop
end program pipe_network


subroutine gauss(A,b,x,ref,neqn,nrhs,scale,solution,det,pivot,tol)
  implicit none
  double precision,dimension(:,:),intent(in)::A,b
  double precision,dimension(:,:),intent(out)::x
  double precision,dimension(neqn,neqn+nrhs),intent(out)::ref
  double precision,intent(out)::det
  double precision,intent(in)::tol
  integer,intent(in)::neqn,nrhs
  logical,intent(in)::scale,pivot
  logical,intent(out)::solution
  double precision,dimension(neqn,neqn+nrhs)::AB
  double precision,dimension(neqn+nrhs)::dummy
  double precision::max,sum,hold
  integer::maxpos,row,i,j,k,l,imax

  !This subroutine will input an coefficient matrix A and a right hand 
  !side matrix b from the form Ax=b, augment one with the other to get 
  !a new matrix AB, scale if desired, find the ref of the AB, and using 
  !back substitution, find the solution matrix x associated with A and b.
  !~~b may be a vector or a matrix~~
  !variable list:
  !--inputs:
  !A        =coefficient matrix(pressure values in this case)
  !          must be square
  !b        =right hand side vector
  !neqn     =number of equations/unknowns in matrix A (number of rows)
  !nrhs     =number of columns in right hand side vector b
  !scale    =if ans is yes then scale=.true. else scale=.false.
  !pivot    =if ans is yes then pivot=.true. else pivot=.false.
  !tol      =matrix singularity tolerance
  !--outputs:
  !x        =solution vector
  !det      =determinant of matrix A
  !ref      =row eschelon form of matrix A|b
  !solution =if the system has a solution then solution=.true.
  !--internal variables:
  !AB       =matrix A augmented with matrix b
  !max      =max of row, used for scaling
  !row      =working row for augmentation
  !sum      =sum of terms used in back substitution
  !hold     =temporary hold of values in row reduction and row interchange
  !imax     =row with max value for each working column (used to indicate if 
  !          row interchange is necessary
  !i,j,k,l  =loop variables
 

  AB=A
  row=0
  do j=neqn+1,neqn+nrhs    !Augment A with b
    row=row+1
    do i=1,neqn
      AB(i,j)=b(i,row)
    end do
  end do

  det=1d0
  if(scale)then       !This if,then statement is to scale
    do i=1,neqn
      max=0d0
      maxpos=1
      do j=1,neqn+nrhs
        if(abs(A(i,j))>abs(max))then
          max=A(i,j)
          maxpos=j
        end if
      end do
      AB(i,1:neqn+nrhs)=AB(i,1:neqn+nrhs)/max
      det=det*max
    end do
  end if

  
  do j=1,neqn                        !this whole loop for gauss ref
    imax=1 
    if(pivot)then
      do i=j,neqn                    !check for row max
        if(abs(AB(i,j))>abs(AB(imax,j)))then
          imax=i
        end if
      end do
      if(imax/=j)then
        do k=j,neqn+nrhs
          hold=AB(j,k)
          AB(j,k)=AB(imax,k)        !row interchange
          AB(imax,k)=hold
               !write(*,*)"  ",k
          end do
        end if
    end if
    solution=.true.
    if(abs(AB(j,j))<tol)then
      solution=.false.
      write(*,*)"NO SOLUTION"
      return
      stop
    end if
    do i=j+1,neqn                    !row reduction
      hold=AB(i,j)/AB(j,j)
      do k=j,neqn+nrhs            
        AB(i,k)=AB(i,k)-(hold)*AB(j,k)      
      end do
    end do
  end do

  do k=1,nrhs                    !back substitution, last row first
    x(neqn,k)=AB(neqn,neqn+k)/AB(neqn,neqn)
    do i=neqn-1,1,-1             !rest of rows working up
      sum=0
          !sum(AB(i,i+1:neqn)*x(i+1:neqn,k))
      do j=i+1,neqn
        sum=sum+x(j,k)*AB(i,j)
      end do
      x(i,k)=(AB(i,neqn+k)-sum)/AB(i,i)
    end do
  end do

  do i=1,neqn
    det=det*AB(i,i)        !calculate the determinate by multiplying
  end do                   !diagonal elements because AB is triangular
  
  ref=0
  do i=1,neqn              !divide each row by leading entry to get ref
    ref(i,1:neqn+nrhs)=AB(i,1:neqn+nrhs)/AB(i,i)
  end do
  
  return
end subroutine gauss
