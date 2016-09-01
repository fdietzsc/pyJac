program test_pyjac

  use pyjac_mod
  use, intrinsic :: iso_c_binding

  implicit none

  integer, parameter :: nspec=10
  real(8), allocatable :: hi(:)
  real(8), allocatable :: cpi(:)
  real(8), allocatable :: X(:)
  real(8), allocatable :: Y(:)
  real(8), allocatable :: omega(:)
  real(8) :: T = 1000.0D0
  real(8) :: P = 101325.0
  real(8) :: rho = 0.0
  integer :: i

  allocate(hi(nspec))
  allocate(cpi(nspec))
  allocate(X(nspec))
  allocate(Y(nspec))
  allocate(omega(nspec))

  ! get hi and cpi
  call compute_hi(T,hi)
  call compute_cpi(T,cpi)
  do i=1,nspec
    print'(2es12.5)',hi(i),cpi(i)
  end do

  print*,''

  ! forward and backward conversion of X and Y
  ! density is also computed
  X = 0.1
  call convert_x2y(X,Y)
  do i=1,nspec
    print'(2es12.5)',X(i),Y(i)
  end do
  print*,''
  rho =  computeDensity(T,P,X) 
  print*,'density: ',rho
  print*,''
  call convert_y2x(Y,x)
  do i=1,nspec
    print'(2es12.5)',Y(i),X(i)
  end do
  print*,''
  rho =  computeDensity(T,P,X) 
  print*,'density: ',rho

  print*,''
  ! forward and backward conversion of X and Y
  Y = 0.1
  call convert_y2x(Y,x)
  do i=1,nspec
    print'(2es12.5)',Y(i),X(i)
  end do
  print*,''
  call convert_x2y(X,Y)
  do i=1,nspec
    print'(2es12.5)',X(i),Y(i)
  end do
  ! source terms
  print*,''
  Y(1) = T ! first entry should be temperature
  call reaction_rates(T,P,Y,omega)
  do i=1,nspec
    print'(2es12.5)',Y(i),omega(i)
  end do


end program test_pyjac
