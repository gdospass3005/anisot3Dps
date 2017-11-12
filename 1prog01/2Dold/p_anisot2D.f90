PROGRAM anisot

  USE mdle_source
  USE mdle_prop
  USE mdle_utils
  USE mdle_taper
  USE mdle_io_utils                                                  
  USE mdle_fft

  IMPLICIT NONE
  REAL, DIMENSION(:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz

  INTEGER :: Nx,Nz
  REAL    :: dx,dz,dt,fpeak
  INTEGER :: itmax
  INTEGER :: lx, lz
  REAL    :: r
  REAL, DIMENSION(:,:), ALLOCATABLE :: C11,C13,C33,C44
  REAL, DIMENSION(:,:), ALLOCATABLE :: rhox

  REAL, DIMENSION(:,:), ALLOCATABLE :: source
  REAL, DIMENSION(:),   ALLOCATABLE :: trs
  REAL    :: tdelay

  REAL, DIMENSION(:), ALLOCATABLE :: taper  
  REAL :: F
  INTEGER :: nb,fsf

  REAL    :: t
  INTEGER :: it,i,j
  INTEGER :: aux,ir

  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_Vx,csg_Vz,csg_P,csg_S
  INTEGER :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER :: isnap,nsnaps,snapmin,dsnap

  REAL, DIMENSION(:), ALLOCATABLE  :: trig_x
  REAL, DIMENSION(:), ALLOCATABLE  :: trig_z
  INTEGER, DIMENSION(:), ALLOCATABLE  :: ifax_x
  INTEGER, DIMENSION(:), ALLOCATABLE  :: ifax_z
  INTEGER  :: nfac_x,nfac_z
  REAL, DIMENSION(:), ALLOCATABLE  :: rkx
  REAL, DIMENSION(:), ALLOCATABLE  :: rkz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,r,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
       nb,F,fsf)

  dz=dx

  ALLOCATE(C11(Nx,Nz),C13(Nx,Nz),C33(Nx,Nz),C44(Nx,Nz))
  ALLOCATE(rhox(Nx,Nz))
  !C11=1.181; C13=0.711; C33=0.968; C44=0.170; rhox=2250.;

  CALL inputmodel(Nx,Nz,C11,C13,C33,C44,rhox)
  C11 = C11*1.e10
  C13 = C13*1.e10
  C33 = C33*1.e10
  C44 = C44*1.e10

  PRINT*,'anisot - Modeling seismic waves in anisotropic media'
  PRINT*,''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL quads(itmax,aux)
  ALLOCATE(trs(aux))
  CALL source_dervgauss_init(trs,dt,aux,fpeak,tdelay)
  ir = CEILING(r)
  ALLOCATE(source(2*ir+1,2*ir+1))

  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

  ! Initialize fft vectors
  ALLOCATE(trig_x(2*Nx),ifax_x(Nx),rkx(Nx))
  ALLOCATE(trig_z(2*Nz),ifax_z(Nz),rkz(Nz))
  call fft_init(NX,trig_x,ifax_x,nfac_x,dx,rkx)
  call fft_init(NZ,trig_z,ifax_z,nfac_z,dz,rkz)

  ALLOCATE(csg_Vx(itmax,nphones),csg_Vz(itmax,nphones))
  ALLOCATE(csg_P(itmax,nphones),csg_S(itmax,nphones))

  n1=31; n2=32; n3=33; n4=34; n5=35; n6=36; n7=37; n8=38
  OPEN(n1,FILE="csg_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n2,FILE="csg_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)

  OPEN(n3,FILE="snapshots_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(n4,FILE="snapshots_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

  OPEN(n5,FILE="csg_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n6,FILE="snapshots_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

  ! OPEN(n7,FILE="csg_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
  !      ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  ! OPEN(n8,FILE="snapshots_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
  !      ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(Vx(Nx,Nz),Vz(Nx,Nz),Sxx(Nx,Nz),Sxz(Nx,Nz),Szz(Nx,Nz))
  Sxx=0.; Sxz=0.; Szz=0.; Vx=0.; Vz=0.
  rhox=1./rhox;

  isnap=0
  DO it=1,itmax
     t=it*dt
     OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
     WRITE(77,*) &
          'anisot - Modeling seismic waves in anisotropic media'
     WRITE(77,*) &
          'Iteracao',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)
     CLOSE(77)

     ! Insert source function 
     CALL source_UU1_2Dc(t,it,trs,source,r,tdelay)
     Vx(lx-ir:lx+ir,lz-ir:lz+ir) = & 
          Vx(lx-ir:lx+ir,lz-ir:lz+ir) + source    
     CALL source_UU3_2Dc(t,it,trs,source,r,tdelay)
     Vz(lx-ir:lx+ir,lz-ir:lz+ir) = &
          Vz(lx-ir:lx+ir,lz-ir:lz+ir) + source

    ! Sxx(lx,lz)=trs(it)
    ! Szz(lx,lz)=trs(it)

     CALL prop_anisot(Vx,Vz,Sxx,Sxz,Szz,&
          C11,C13,C33,C44,rhox,dx,dt,Nx,Nz,&
             trig_x,ifax_x,nfac_x,rkx,trig_z,ifax_z,nfac_z,rkz)

     PRINT*,'it',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)

     CALL save_shotgather_n_snapshots(csg_Vx,csg_Vz,&
          csg_P,csg_S,Vx,Vz,&
          Sxx,Sxz,Szz,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
          isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
          Nx,Nz,it,itmax)

     CALL bt_apply_multiple(Vx,Vz,Sxx,Sxz,Szz,Nx,Nz,nb,taper,fsf)

  END DO

  PRINT*,'Successful run.'

END PROGRAM anisot
