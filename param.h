!!!!!!!!!!!!!!!!!!!!!!
!
! Parameter file
!
!!!!!!!!!!!!!!!!!!!!!!

! parameters !
integer, parameter :: n_x = 64, n_y = 64, n_z = 64
real*8, parameter :: small_num = 1.0d-50

! parameters !
real*8, parameter :: pi = 3.1415926535897932384626433832795D0
integer, parameter :: max_iter = 10000000

! density !
real*8, parameter :: rho_a = 0.0d0
real*8, parameter :: rho_0 = 1.0d0
real*8, parameter :: r_surf = 0.5
real*8, parameter :: mass = (4.0d0/3.0d0)*pi*rho_0*r_surf**3

! parameters !
real*8, parameter :: x2_start = -1.0d0 
real*8, parameter :: x2_end = 1.0d0 
real*8, parameter :: y2_start = -1.0d0 
real*8, parameter :: y2_end = 1.0d0 
real*8, parameter :: z2_start = -1.0d0 
real*8, parameter :: z2_end = 1.0d0 

! parameters !
REAL*8, PARAMETER :: dx2 = (x2_end - x2_start)/DBLE(n_x)	
REAL*8, PARAMETER :: dy2 = (y2_end - y2_start)/DBLE(n_y)	
REAL*8, PARAMETER :: dz2 = (z2_end - z2_start)/DBLE(n_z)	

! strecth grid? !
LOGICAL, PARAMETER :: strecth = .true.
