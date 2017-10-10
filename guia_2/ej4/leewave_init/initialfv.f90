! External mechanical forcing.
! This file contains the expression used for the external
! mechanical forcing. You can use temporary real arrays R1-R3
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate
! computations. The variable f0 should control the global
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz.

! Small vertical perturbation in a narrow vertical band

      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               fx(k,j,i) = 0.
               fy(k,j,i) = 0.
            END DO
         END DO
      END DO

      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
               IF (i.le.(nx/20)) THEN
                  R1(i,j,k) = f0
               ELSE
                  R1(i,j,k) = 0.
               ENDIF
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,fz,MPI_COMM_WORLD)
