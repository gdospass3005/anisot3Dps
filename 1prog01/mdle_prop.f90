MODULE mdle_prop
CONTAINS
  SUBROUTINE prop_anisot3D(Vx,Vy,Vz,Sxx,Sxy,Sxz,Syy,Syz,Szz,&
       C11,C12,C13,C33,C44,C66,rhox,dx,dt,Nx,Ny,Nz,&
       trig_x,ifax_x,nfac_x,rkx,trig_y,ifax_y,nfac_y,rky,&
       trig_z,ifax_z,nfac_z,rkz)
    USE mdle_fft
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Vx,Vz,Sxx,Sxz,Szz
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: Vy,Syy,Sxy,Syz
    REAL, DIMENSION(:,:,:), INTENT(IN)    :: C11,C12,C13,C33,C44,C66,rhox
    REAL,    INTENT(IN) :: dx,dt
    INTEGER, INTENT(IN) :: Nx,Ny,Nz
    INTEGER :: i,j,k
    REAL, DIMENSION(2*Nx), INTENT(IN)  :: trig_x
    REAL, DIMENSION(2*Ny), INTENT(IN)  :: trig_y
    REAL, DIMENSION(2*Nz), INTENT(IN)  :: trig_z
    INTEGER, DIMENSION(Nx), INTENT(IN)  :: ifax_x
    INTEGER, DIMENSION(Ny), INTENT(IN)  :: ifax_y
    INTEGER, DIMENSION(Nz), INTENT(IN)  :: ifax_z
    INTEGER, INTENT(IN)  :: nfac_x,nfac_y,nfac_z
    REAL, DIMENSION(Nx), INTENT(IN)  :: rkx
    REAL, DIMENSION(Ny), INTENT(IN)  :: rky
    REAL, DIMENSION(Nz), INTENT(IN)  :: rkz
    REAL, DIMENSION(Nx) :: aux1
    REAL, DIMENSION(Ny) :: aux2
    REAL, DIMENSION(Nz) :: aux3

    ! Space derivatives of velocities => stresses
    DO k=1,Nz
       DO j=1,Ny
          aux1=Vx(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Sxx(:,j,k) = Sxx(:,j,k) + (dt * C11(:,j,k) * aux1)
          Syy(:,j,k) = Syy(:,j,k) + (dt * C12(:,j,k) * aux1)
          Szz(:,j,k) = Szz(:,j,k) + (dt * C13(:,j,k) * aux1)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Vx(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Sxx(i,:,k) = Sxx(i,:,k) + (dt * C12(i,:,k) * aux2)
          Syy(i,:,k) = Syy(i,:,k) + (dt * C11(i,:,k) * aux2)
          Szz(i,:,k) = Szz(i,:,k) + (dt * C13(i,:,k) * aux2)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Vz(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Sxx(i,j,:) = Sxx(i,j,:) + (dt * C13(i,j,:) * aux3)
          Syy(i,j,:) = Syy(i,j,:) + (dt * C13(i,j,:) * aux3)
          Szz(i,j,:) = Szz(i,j,:) + (dt * C33(i,j,:) * aux3)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Vz(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Syz(i,:,k) = Syz(i,:,k) + (dt * C44(i,:,k) * aux2)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Vy(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Syz(i,j,:) = Syz(i,j,:) + (dt * C44(i,j,:) * aux3)
       END DO
    END DO
    DO k=1,Nz
       DO j=1,Ny
          aux1=Vz(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Sxz(:,j,k) = Sxz(:,j,k) + (dt * C44(:,j,k) * aux1)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Vx(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Sxz(i,j,:) = Sxz(i,j,:) + (dt * C44(i,j,:) * aux3)
       END DO
    END DO
    DO k=1,Nz
       DO j=1,Ny
          aux1=Vy(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Sxy(:,j,k) = Sxy(:,j,k) + (dt * C66(:,j,k) * aux1)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Vx(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Sxy(i,:,k) = Sxy(i,:,k) + (dt * C66(i,:,k) * aux2)
       END DO
    END DO

    ! Space derivatives of stresses => velocities
    DO k=1,Nz
       DO j=1,Ny
          aux1=Sxx(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Vx(:,j,k) = Vx(:,j,k) + (dt * rhox(:,j,k) * aux1)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Sxy(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Vx(i,:,k) = Vx(i,:,k) + (dt * rhox(i,:,k) * aux2)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Sxz(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Vx(i,j,:) = Vx(i,j,:) + (dt * rhox(i,j,:) * aux3)
       END DO
    END DO
    DO k=1,Nz
       DO j=1,Ny
          aux1=Sxy(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Vy(:,j,k) = Vy(:,j,k) + (dt * rhox(:,j,k) * aux1)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Syy(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Vy(i,:,k) = Vy(i,:,k) + (dt * rhox(i,:,k) * aux2)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Syz(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Vy(i,j,:) = Vy(i,j,:) + (dt * rhox(i,j,:) * aux3)
       END DO
    END DO
    DO k=1,Nz
       DO j=1,Ny
          aux1=Sxz(:,j,k)
          call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
          Vz(:,j,k) = Vz(:,j,k) + (dt * rhox(:,j,k) * aux1)
       END DO
    END DO
    DO k=1,Nz
       DO i=1,Nx
          aux2=Syz(i,:,k)
          call ddx(aux2,Ny,trig_y,ifax_y,nfac_y,rky)
          Vz(i,:,k) = Vz(i,:,k) + (dt * rhox(i,:,k) * aux2)
       END DO
    END DO
    DO j=1,Ny
       DO i=1,Nx
          aux3=Szz(i,j,:)
          call ddx(aux3,Nz,trig_z,ifax_z,nfac_z,rkz)
          Vz(i,j,:) = Vz(i,j,:) + (dt * rhox(i,j,:) * aux3)
       END DO
    END DO

  END SUBROUTINE prop_anisot3D
END MODULE mdle_prop
