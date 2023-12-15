!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! testing uneven-mesh poisson solver in cartesian coordinates 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program poisson_solver
USE HDF5
implicit none
include "param.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Arrays !
real*8, dimension(-2:n_x+3, -2:n_y+3, -2:n_z+3) :: V, rho, V_old
real*8, dimension(-2:n_x+3) :: x2
real*8, dimension(-3:n_x+3) :: xf2
real*8, dimension(-2:n_y+3) :: y2
real*8, dimension(-3:n_y+3) :: yf2
real*8, dimension(-2:n_z+3) :: z2
real*8, dimension(-3:n_z+3) :: zf2

! For poisson solver !
real*8, dimension(1:n_x) :: ajp1, ajm1
real*8, dimension(1:n_y) :: bkp1, bkm1
real*8, dimension(1:n_z) :: clp1, clm1
real*8, dimension(1:n_x,1:n_y,1:n_z) :: epsc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Integer !
INTEGER :: i, j, k, n

! Real !
REAL*8 :: dx_big
REAL*8 :: dx_small
REAL*8 :: ratio
REAL*8 :: tol
REAL*8 :: rad
REAL*8 :: err, tmp
REAL*8 :: rho_in, rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for HDF5 !
character(len=99) :: filename
integer :: error, space_rank
integer(HSIZE_T) :: data_dims(3), scalar_dims(1)
integer(HID_T) :: file_id, dspace_id, dset_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: time_start, time_end
INTEGER :: cr
REAL*8 :: rate
CALL system_clock(count_rate=cr)
rate = REAL(cr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set number of threads !
CALL OMP_SET_NUM_THREADS(56)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! machine epsilon !
tol = 1.0d-10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize mesh at cell interface !

! strecth grid !
IF(strecth) THEN

	! Solve for the dx !
	ratio = 1.01d0
	dx_small = (ratio-1.0d0)/(ratio**(FLOAT(n_x)/2) - 1.0d0)	
	
	! get cell-interface !
	DO i = FLOAT(n_x)/2 + 1, n_x + 3
		xf2(i) = xf2(i-1) + dx_small*ratio**(i - (FLOAT(n_x)/2 + 1))
	END DO
	DO i = FLOAT(n_x)/2 - 1, -3, -1
		xf2(i) = -xf2(n_x - i)
	END DO

	! Get y and z !
	yf2(:) = xf2(:)
	zf2(:) = xf2(:)

	! Get dx_big !
	dx_big = dx_small*(ratio**(FLOAT(n_x)/2-1))

ELSE

	! Radius !
	Do i = -3, n_x + 3
		xf2(i) = x2_start + DBLE(i)*dx2
	END DO
	Do j = -3, n_y + 3
		yf2(j) = y2_start + DBLE(j)*dy2
	END DO
	Do k = -3, n_z + 3
		zf2(k) = z2_start + DBLE(k)*dz2
	END DO

	! dx big and dx small !
	dx_big = dx2
	dx_small = dx2

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mesh at cell center !

! Radius !
Do i = -2, n_x + 3
  x2(i) = 0.5d0*(xf2(i) + xf2(i-1))
END DO
Do j = -2, n_y + 3
	y2(j) = 0.5d0*(yf2(j) + yf2(j-1))
END DO
Do k = -2, n_z + 3
	z2(k) = 0.5d0*(zf2(k) + zf2(k-1))
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get poisson coefficient !
Do i = 1, n_x
	Do j = 1, n_y
		Do k = 1, n_z
			CALL poisson_coef(x2(i-1), x2(i), x2(i+1), &
												y2(j-1), y2(j), y2(j+1), &
												z2(k-1), z2(k), z2(k+1), &
												ajp1(i), ajm1(i), bkp1(j), bkm1(j), &
												clp1(k), clm1(k), epsc(i,j,k))
		end do
	end do
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set initial density !
DO i = 1, n_x
	DO j = 1, n_y
		DO k = 1, n_z
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			if(rad < r_surf) then
			  rho(i,j,k) = 1.0d0
			else
				rho(i,j,k) = 0.0d0
			end if
		end do
	end do
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set initial guess !
V(:,:,:) = -1.0d0

! Set boundary conditions !
do k = 0, n_z+1
	do j = 0, n_y+1 
		do i = n_x+1, n_x+1
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do
do k = 0, n_z+1
	do j = 0, n_y+1 
		do i = 0, 0
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do

! Set boundary conditions !
do k = 0, n_z+1
	do j = n_y+1, n_y+1 
		do i = 0, n_x+1
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do
do k = 0, n_z+1
	do j = 0, 0
		do i = 0, n_x+1
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do

! Set boundary conditions !
do k = n_z+1, n_z+1
	do j = 0, n_y+1
		do i = 0, n_x+1
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do
do k = 0, 0
	do j = 0, n_y+1
		do i = 0, n_x+1
			rad = DSQRT(x2(i)**2+y2(j)**2+z2(k)**2)
			V(i,j,k) = -mass/rad
		end do
	end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$ACC DATA COPY(v_old, V, rho, ajp1, ajm1, bkp1, bkm1, clp1, clm1, epsc)

CALL system_clock(time_start)

! Solve the poisson equation !
! Loop over !
do n = 1, max_iter

	!$OMP PARALLEL PRIVATE(rho_in)
  ! Back up !
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
  !$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT)
	do k = 0, n_z + 1
		do j = 0, n_y + 1
			do i = 0, n_x + 1
          v_old(i,j,k) = V(i,j,k)
        END DO
      END DO
    END DO
  !$ACC END PARALLEL
	!$OMP END DO

	! Minimum error !
	!$OMP SINGLE
	err = small_num
  !$OMP END SINGLE

	! Loop over !
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) PRIVATE (rho_in)
	do k = 1, n_z
		do j = 1, n_y 
			do i = 1, n_x 
				IF(rho(i,j,k) > rho_a) THEN
					rho_in = rho(i,j,k)
				ELSE
					rho_in = 0.0d0
				END IF
				V(i,j,k) = (ajp1(i)*v_old(i+1,j,k) + ajm1(i)*v_old(i-1,j,k) + & 
							 	  	bkp1(j)*v_old(i,j+1,k) + bkm1(j)*v_old(i,j-1,k) + & 
							 			clp1(k)*v_old(i,j,k+1) + clm1(k)*v_old(i,j,k-1) - 4.0d0*pi*rho_in)/epsc(i,j,k)
			end do
		end do
	end do
	!$ACC END PARALLEL
	!$OMP END DO

	! Update error 
	!$OMP DO COLLAPSE(3) SCHEDULE(STATIC) Reduction(MAX:err)
	!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(3) DEFAULT(PRESENT) Reduction(MAX:err)
	do k = 1, n_z
		do j = 1, n_y 
			do i = 1, n_x 
				err = max(err, abs(V(i, j, k) - v_old(i, j, k)))
			end do
		end do
	end do
  !$ACC END PARALLEL
  !$OMP END DO
  
	!$OMP END PARALLEL
	! Check maximum error !
	!WRITE (*,*) n, err

	! Check for convergence
	if (err < tol) exit
  
end do

CALL system_clock(time_end)
WRITE(*,*) 'total time = ', REAL(time_end - time_start) / rate

!$ACC END DATA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE (*,*) V(30, 30, 30)

! assign !
filename = 'results.hdf5'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! create interface !
call h5open_f(error)

! open the file !
call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
scalar_dims(1) = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,scalar_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"n_iter",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,n,scalar_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,scalar_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"dx_small",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,dx_small,scalar_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,scalar_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"dx_big",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,dx_big,scalar_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,scalar_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"r_surf",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,r_surf,scalar_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,scalar_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"mass",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,mass,scalar_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
data_dims(1) = n_x

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"x2",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,x2(1:n_x),data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
data_dims(1) = n_y

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"y2",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,y2(1:n_y),data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 1
data_dims(1) = n_z

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"z2",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,z2(1:n_z),data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define DIMENSION !
space_rank = 3
data_dims(1) = n_x
data_dims(2) = n_y
data_dims(3) = n_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! open dataspace !
call h5screate_simple_f(space_rank,data_dims,dspace_id,error)

! create dataset !
call h5dcreate_f(file_id,"phi",H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)

! write dataset !
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,V(1:n_x,1:n_y,1:n_z),data_dims,error)

! close dataset !
call h5dclose_f(dset_id,error)

! close data space !
call h5sclose_f(dspace_id,error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! close the file !
call h5fclose_f(file_id,error)

! close interface !
call h5close_f(error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program poisson_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get left hand side of the discrete poisson equation 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE poisson_coef(xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1, &
	alphajp1, alphajm1, betakp1, betakm1, gammalp1, gammalm1, epsc)


implicit none
include "param.h"

! input !
real*8, INTENT(IN) :: xm1, xc, xp1, ym1, yc, yp1, zm1, zc, zp1

! output !
real*8, INTENT(OUT) :: epsc
real*8, INTENT(OUT) :: alphajp1, alphajm1
real*8, INTENT(OUT) :: betakp1, betakm1
real*8, INTENT(OUT) :: gammalp1, gammalm1

! assign !
epsc = 2.0d0*(1.0d0/(xp1-xc)/(xc-xm1) + 1.0d0/(yp1-yc)/(yc-ym1) + 1.0d0/(zp1-zc)/(zc-zm1))

! assign !
alphajp1 = 2.0d0/(xp1-xm1)/(xp1-xc)
alphajm1 = 2.0d0/(xp1-xm1)/(xc-xm1)

! assign !
betakp1 = 2.0d0/(yp1-yc)/(yp1-ym1)
betakm1 = 2.0d0/(yp1-ym1)/(yc-ym1)

! assign !
gammalp1 = 2.0d0/(zp1-zm1)/(zp1-zc)
gammalm1 = 2.0d0/(zp1-zm1)/(zc-zm1) 

END SUBROUTINE
