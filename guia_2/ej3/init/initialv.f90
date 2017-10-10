! Initial condition for the velocity.
! This file contains the expression used for the initial 
! velocity field. You can use temporary real arrays R1-R3 
! of size (1:n,1:n,ksta:kend) and temporary complex arrays 
! C1-C8 of size (n,n,ista:iend) to do intermediate 
! computations. The variable u0 should control the global 
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Superposition of harmonic modes with random phases in a
! Fourier shell [kdn,kup]. We create three random components
! for the initial velocity field in Fourier space. Modes with
! the same wave number has the same amplitude but different
! random phases. The amplitudes of the modes with different
! wave numbers decay as k^(-2/5).
!     kdn    : minimum wave number
!     kup    : maximum wave number
  
      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))**5
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kk2(k,1,1).le.kup**2).and.(kk2(k,1,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,1,1))**5
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.kup**2).and.(kk2(k,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(k,j,1))**5
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**5
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.kup**2).and.(kk2(k,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(k,j,i))**5
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

! C1,C2,C3 are the components of the random field. We multiply
! it by zeros outside a thin slab to confine the perturbation
! to a horizontal band. Note we must first go to real space.

      CALL fftp3d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C2,R2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,C3,R3,MPI_COMM_WORLD)

      DO k = ksta,kend
         rmp = 1.0_GP
         IF ((k.lt.real(nz)/3).or.(k.gt.2*real(nz)/3)) rmp = 0.0_GP
         DO j = 1,ny
            DO i = 1,nx
               R1(i,j,k) = R1(i,j,k)*rmp
               R2(i,j,k) = R3(i,j,k)*rmp
               R3(i,j,k) = R3(i,j,k)*rmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,C1,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,C2,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,C3,MPI_COMM_WORLD)

! We go back to Fourier space. Then we project the field into
! its solenoidal component (div.v = 0) and normalize its amplitude.
      
      CALL nonlhd3(C1,C2,C3,vx,1)
      CALL nonlhd3(C1,C2,C3,vy,2)
      CALL nonlhd3(C1,C2,C3,vz,3)      
      CALL normalize(vx,vy,vz,u0,1,MPI_COMM_WORLD)
