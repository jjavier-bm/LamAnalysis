! **********************************************************************
! *** LAMELLAR DI-BLOCK COPOLYMER + NANODIMERS *************************
! *** Structural Data Analysis                 *************************
! *** Author: J. Javier Burgos Mármol -- 2020  *************************
! **********************************************************************
PROGRAM LamellarStructure
IMPLICIT NONE
!---DATA------------------------------------
CHARACTER(64) :: inputfilename='stdyn.pure.lammpstrj' 
!CONSTANTS
REAL*8, PARAMETER :: PI=3.141592654d0, TWOPI=6.283185307d0, HALFPI=1.570796327d0
!PARAMETERS
INTEGER, PARAMETER :: Nchains = 12696, LchainsA = 15, LchainsB = 15	!Number of Chains and length of each of the blocks
REAL*8, PARAMETER :: Mbeadchains=1.d0 ! Mass of each bead in the chains
INTEGER, PARAMETER :: Npart = 0 !Number of nano-dimers
REAL*8, PARAMETER :: NPdiam=3.d0, RbondEqNP=2.5d0 ! NP bead diameter, NP bond eq. distance (mass derived from here)
INTEGER, PARAMETER :: Ntimesteps=100000000, period=500000, Nconfigsperfile=3 !Number of configurations and timesteps between each consecutive pair of configurations
INTEGER, PARAMETER :: Nconfigs=Ntimesteps/period+1
REAL*8, PARAMETER :: delta_t=0.001d0
INTEGER, PARAMETER :: Nslabs=180		   				!Number of slabs in which the box is divided
INTEGER, PARAMETER :: NAngleSlabsPer180=36 ! Number of slabs per 180 degrees
!OTHER PARAMETERS
INTEGER, PARAMETER :: Lchains=LchainsA+LchainsB
REAL*8, PARAMETER :: Vlens=(PI*(NPdiam-RbondEqNP)**2*(RbondEqNP**2+2.d0*RbondEqNP*NPdiam)/(12.d0*RbondEqNP)) ! Overlapping volume (lens shaped)
REAL*8, PARAMETER :: MbeadNpart=0.74d0*(NPdiam**3/6.d0-Vlens/2.d0) ! Mass of Nano-dimer's beads (scalates as FCC packed)
!OTHER VARIABLES
REAL*8, DIMENSION(:) :: RChain(Nchains,0:Lchains,3) !Positions of chains ! 0 is for the CoM
REAL*8, DIMENSION(:) :: Rpart(Npart,0:2,3) !Positions of Nano-dimers ! 0 is for the CoM
REAL*8, DIMENSION(:) :: Lbox(Nconfigs,0:3)
REAL*8, DIMENSION(:) :: L0box(0:3) ! SIMULATION BOX, 0=volume
REAL*8, DIMENSION(:) :: pos(3), vel(3), jvol(6) ! jvol(stresstensor*volume): xx,yy,zz,xy,xz,yz
INTEGER, DIMENSION(:) :: ir(3)
!REAL*8 :: AvDensMonomersA,AvDensMonomersB,AvDensAllMonomers,AvDensChains,AvDensJNDs
INTEGER ::  i, j, k, conf, kconf, d, d2, d3, nsl, nsl2, id, idmol, idbead, tipo
REAL*8, DIMENSION(:) :: acdensProf(6,3,Nslabs)!, acdensProf2(6,3,Nslabs) !Density profile accumulation (ChA=1,ChB=2,ChCoM=3,NPA=4,NPB=5,NPCoM=6)
REAL*8, DIMENSION(:) :: acdirProf(3,3,0:Nslabs), acdirProf2(3,3,0:Nslabs) !Angle profile accumulation (Polar=1,Azim=2,abs(polar)=3)
REAL*8, DIMENSION(:) :: slabthickness(3), AvDens(4) ! DensProf Slab Thickness, AvDensity (ChA=1,ChB=2,ChCoM=3,NPA=NPB=NPCoM=4)
REAL*8, DIMENSION(:) :: NematicV(3), NematicV2(3), vnem(0:3)	!Nematic directors.
REAL*8, DIMENSION(:) :: NematicT(3,3) !Nematic tensor
REAL*8, PARAMETER :: slabthicknessAngles=180.d0/dfloat(NAngleSlabsPer180)
REAL*8 :: P2, P4,zen,azim,sgn
REAL*8, DIMENSION(:) :: acVdist(2,3,2*NAngleSlabsPer180)	!Angle distribution accumulation,  dims: 1=Theta(polar), 2=Phi(Azim)
INTEGER, DIMENSION(:) :: counterAngles(2,3,NSlabs)

!LAMMPS NOTE: As explained on the dump doc page, atom coordinates in a dump file may be slightly outside the simulation box. This is because periodic boundary conditions are enforced only on timesteps when neighbor lists are rebuilt, which will not typically coincide with the timesteps dump snapshots are written.

!***********************************************************************************************************************************************
!-- DENSITY PROFILES

open(1111,file=trim(inputfilename))
IF (Npart>0) open(777,file='NematicDirEvolution.dat')
acdensProf=0.d0
L0box(0)=1.d0
AvDens=0.d0
conf=0
counterAngles=0
acVdist=0.d0
acdirProf=0.d0
acdirProf2=0.d0
P2=0.d0
P4=0.d0
NematicV=0.d0
NematicV2=0.d0
DO kconf=0,NconfigsperFile*Ntimesteps/((NconfigsperFile-1)*period)-1
    IF (kconf>0 .and. MOD(kconf,NconfigsperFile)==0) THEN
        DO i=1,9+Nchains*Lchains+Npart*2
            READ(1111,*)
        END DO !i
    ELSE
        NematicT=0.d0
        conf=conf+1
        !Read coordinates and accumulate values for beads
        READ(1111,*)
        READ(1111,*)
        READ(1111,*)
        READ(1111,*)
        READ(1111,*)
        Lbox(conf,0)=1.d0
        DO d=1,3
            READ(1111,*) L0box(d), Lbox(conf,d)
            Lbox(conf,d)=Lbox(conf,d)-L0box(d)
            Lbox(conf,0)=Lbox(conf,0)*Lbox(conf,d)
            IF (conf==1) L0box(0)=L0box(0)*Lbox(conf,d) !initial volume
            slabthickness(d)=Lbox(conf,d)/dfloat(Nslabs)
        END DO !dims
        AvDens(1)=AvDens(1)+dfloat(LchainsA*Nchains)/Lbox(conf,0)
        AvDens(2)=AvDens(2)+dfloat(LchainsB*Nchains)/Lbox(conf,0)
        AvDens(3)=AvDens(3)+dfloat(Nchains)/Lbox(conf,0)
        AvDens(4)=AvDens(4)+dfloat(Npart)/Lbox(conf,0)
        Lbox(conf,0)=Lbox(conf,0)/L0box(0) !Volume normalisation V/V_0
        READ(1111,*)
        DO i=1,Nchains*Lchains+Npart*2
            READ(1111,*) id, pos(1), pos(2), pos(3), ir(1), ir(2), ir(3),&
                vel(1), vel(2), vel(3), jvol(1), jvol(2), jvol(3), jvol(4), jvol(5), jvol(6)

            IF (id<=Nchains*Lchains) THEN
                idmol=int((id-1)/Lchains)+1
                idbead=id-(idmol-1)*Lchains
                tipo=1
                IF (idbead>LchainsA) tipo=2
                DO d=1,3
                    Rchain(idmol,idbead,d)=pos(d)-L0box(d)
                    !Put bead back in the box. See LAMMPS note above.
                    IF (Rchain(idmol,idbead,d)<0.d0) Rchain(idmol,idbead,d)=Rchain(idmol,idbead,d)+Lbox(conf,d)
                    IF (Rchain(idmol,idbead,d)>=Lbox(conf,d)) Rchain(idmol,idbead,d)=Rchain(idmol,idbead,d)-Lbox(conf,d)
                    !Accumulator
                    nsl=int(Rchain(idmol,idbead,d)/slabthickness(d))+1
                    acdensProf(tipo,d,nsl)=acdensProf(tipo,d,nsl)+1.d0/Lbox(conf,0)
                END DO !dims
            ELSE
                idmol=int((id-Nchains*Lchains-1)/2)+1
                idbead=id-Nchains*Lchains-(idmol-1)*2
                tipo=4
                IF (idbead==2) tipo=5
                DO d=1,3
                    Rpart(idmol,idbead,d)=pos(d)-L0box(d)
                    !Put bead back in the box. See LAMMPS note above.
                    IF (Rpart(idmol,idbead,d)<0.d0) Rpart(idmol,idbead,d)=Rpart(idmol,idbead,d)+Lbox(conf,d)
                    IF (Rpart(idmol,idbead,d)>=Lbox(conf,d)) Rpart(idmol,idbead,d)=Rpart(idmol,idbead,d)-Lbox(conf,d)
                    !Accumulator
                    nsl=int(Rpart(idmol,idbead,d)/slabthickness(d))+1
                    acdensProf(tipo,d,nsl)=acdensProf(tipo,d,nsl)+1.d0/Lbox(conf,0)
                END DO !dims       
            END IF
        END DO !Total Nbeads
        !Compute Properties for CoMs
        tipo=3
        DO i=1,Nchains
            DO d=1,3
                ir(d)=0
                RChain(i,0,d)=0.d0
                DO j=1,Lchains
                    IF (j>1) THEN
                        pos(d)=Rchain(i,j,d)-Rchain(i,j-1,d)
                        IF (pos(d)>Lbox(conf,d)/2.d0) ir(d)=ir(d)-1
                        IF (pos(d)<-Lbox(conf,d)/2.d0) ir(d)=ir(d)+1
                    END IF
                    RChain(i,0,d)=RChain(i,0,d)+Mbeadchains*(Rchain(i,j,d)+Lbox(conf,d)*ir(d))
                END DO !Lchains
                RChain(i,0,d)=RChain(i,0,d)/(dfloat(Lchains)*Mbeadchains)
                IF (RChain(i,0,d)<0.d0) RChain(i,0,d)=RChain(i,0,d)-(int(RChain(i,0,d)/Lbox(conf,d))-1)*Lbox(conf,d)
                IF (RChain(i,0,d)>=Lbox(conf,d)) RChain(i,0,d)=RChain(i,0,d)-int(RChain(i,0,d)/Lbox(conf,d))*Lbox(conf,d)
                nsl=int(Rchain(i,0,d)/slabthickness(d))+1
                acdensProf(tipo,d,nsl)=acdensProf(tipo,d,nsl)+1.d0/Lbox(conf,0)
            END DO !dims
        END DO !Nchains
        tipo=6
        DO i=1,Npart
            Call DistAB(Rpart(i,2,1),Rpart(i,2,2),Rpart(i,2,3),&
                Rpart(i,1,1),Rpart(i,1,2),Rpart(i,1,3),&
                Lbox(conf,1),Lbox(conf,2),Lbox(conf,3),&
                vnem(1),vnem(2),vnem(3),vnem(0))
            DO d=1,3
                !density profiles
                ir(d)=0
                Rpart(i,0,d)=MbeadNpart*Rpart(i,1,d)
                pos(d)=Rpart(i,2,d)-Rpart(i,1,d)
                IF (pos(d)>Lbox(conf,d)/2.d0) ir(d)=ir(d)-1
                IF (pos(d)<-Lbox(conf,d)/2.d0) ir(d)=ir(d)+1
                Rpart(i,0,d)=Rpart(i,0,d)+MbeadNpart*(Rpart(i,2,d)+Lbox(conf,d)*ir(d))
                Rpart(i,0,d)=Rpart(i,0,d)/(2.d0*MbeadNpart)
                IF (Rpart(i,0,d)<0.d0) Rpart(i,0,d)=Rpart(i,0,d)-(int(Rpart(i,0,d)/Lbox(conf,d))-1)*Lbox(conf,d)
                IF (Rpart(i,0,d)>=Lbox(conf,d)) Rpart(i,0,d)=Rpart(i,0,d)-int(Rpart(i,0,d)/Lbox(conf,d))*Lbox(conf,d)

                nsl=int(Rpart(i,0,d)/slabthickness(d))+1
                acdensProf(tipo,d,nsl)=acdensProf(tipo,d,nsl)+1.d0/Lbox(conf,0)

               !orientational behaviour
                DO d2=1,3
                    NematicT(d,d2)=NematicT(d,d2)+vnem(d)*vnem(d2)/vnem(0)
                END DO !dims2
                
                zen=180.d0*dacos(vnem(d)/dsqrt(vnem(0)))/PI-90.d0
                nsl2=int((zen+90.d0)/slabthicknessAngles)+1
                IF (nsl2 == int(180.d0/slabthicknessAngles)+1) nsl2=nsl2-1
                acVdist(1,d,nsl2)=acVdist(1,d,nsl2)+1.d0

                d2=d-1
                d3=d+1
                IF (d2==0) d2=3    
                IF (d3==4) d3=1
                sgn=dfloat(nint(vnem(d2)/abs(vnem(d2))))
                azim=180.d0*sgn*dacos(vnem(d3)/dsqrt((vnem(d2)**2+vnem(d3)**2)))/PI
                nsl2=int((azim+180.d0)/slabthicknessAngles)+1
                IF (nsl2 == int(360.d0/slabthicknessAngles)+1) nsl2=nsl2-1
                acVdist(2,d,nsl2)=acVdist(2,d,nsl2)+1.d0

                counterAngles(1,d,nsl)=counterAngles(1,d,nsl)+1
                counterAngles(2,d,nsl)=counterAngles(2,d,nsl)+1
                acdirProf(1,d,nsl)=acdirProf(1,d,nsl)+zen
                acdirProf2(1,d,nsl)=acdirProf2(1,d,nsl)+zen**2
                acdirProf(2,d,nsl)=acdirProf(2,d,nsl)+azim
                acdirProf2(2,d,nsl)=acdirProf2(2,d,nsl)+azim**2
                acdirProf(3,d,nsl)=acdirProf(3,d,nsl)+dabs(zen)
                acdirProf2(3,d,nsl)=acdirProf2(3,d,nsl)+zen**2

            END DO !dims
        END DO !Npart
        IF (Npart>0) THEN
            NematicT=NematicT/dfloat(Npart)
            Call Order(NematicT,vnem(1),vnem(2),vnem(3),vnem(0))
            P2=P2+vnem(0)
            P4=P4+vnem(0)**2
            DO d=1,3
                NematicV(d)=NematicV(d)+vnem(d)
                NematicV2(d)=NematicV2(d)+vnem(d)**2
            END DO !dims

            WRITE(777,777) dfloat(period*(conf-1))*delta_t, vnem(0), vnem(1), vnem(2), vnem(3)
        END IF !Npart>0
    END IF !ignore every
END DO !kconf Nconfigs
close(1111)
IF (Npart>0) close(777)
!-- OUTPUT DENSITY PROFILES
open(500,file='densprof_monomers.dat')
open(501,file='densprof_monomersH.dat')
open(502,file='densprof_monomersT.dat')
open(503,file='densprof_chainsCM.dat')
open(600,file='densprof_monomers_notNorm.dat')
open(601,file='densprof_monomersH_notNorm.dat')
open(602,file='densprof_monomersT_notNorm.dat')
IF (Npart>0) THEN
    open(504,file='densprof_NPsA.dat')!0, for 00, H for the rest
    open(505,file='densprof_NPsB.dat')!0,H or T
    open(506,file='densprof_NPsCM.dat')
END IF

acdensProf=dfloat(Nslabs)*acdensProf/(dfloat(Nconfigs)*L0box(0))
AvDens=AvDens/dfloat(Nconfigs)

DO i=1,3
    DO k=1,Nslabs
        IF (i==1) THEN
            WRITE (500, 500) real((dfloat(k)-0.5)/Nslabs),&
                (real(acdensProf(1,1,k))+real(acdensProf(2,1,k)))/(AvDens(1)+AvDens(2)),&
                (real(acdensProf(1,2,k))+real(acdensProf(2,2,k)))/(AvDens(1)+AvDens(2)),&
                (real(acdensProf(1,3,k))+real(acdensProf(2,3,k)))/(AvDens(1)+AvDens(2))
            WRITE (600, 500) real((dfloat(k)-0.5)/Nslabs),&
                real(acdensProf(1,1,k))+real(acdensProf(2,1,k)),&
                real(acdensProf(1,2,k))+real(acdensProf(2,2,k)),&
                real(acdensProf(1,3,k))+real(acdensProf(2,3,k))
        END IF
        IF (i==1 .or. i==2) THEN
            WRITE (500+i, 500) real((dfloat(k)-0.5)/Nslabs), real(acdensProf(i,1,k))/(AvDens(1)+AvDens(2)),&
                real(acdensProf(i,2,k))/(AvDens(1)+AvDens(2)), real(acdensProf(i,3,k))/(AvDens(1)+AvDens(2))
            WRITE (600+i, 500) real((dfloat(k)-0.5)/Nslabs), real(acdensProf(i,1,k)),&
                real(acdensProf(i,2,k)), real(acdensProf(i,3,k))
        ELSE IF (i==3) THEN
            WRITE (500+i, 500) real((dfloat(k)-0.5)/Nslabs), real(acdensProf(i,1,k))/AvDens(i),&
                 real(acdensProf(i,2,k))/AvDens(i), real(acdensProf(i,3,k))/AvDens(i)
        END IF
        IF (Npart>0) THEN
            WRITE (500+i+3, 500) real((dfloat(k)-0.5)/Nslabs), real(acdensProf(i+3,1,k))/AvDens(4),&
                real(acdensProf(i+3,2,k))/AvDens(4), real(acdensProf(i+3,3,k))/AvDens(4)
        END IF
    END DO !k Nslabs
END DO !i type
close(500)
DO i=1,3
    close(500+i)
    IF (Npart>0) close(500+i+3)
END DO !output files
!***********************************************************************************************************************************************
!-- ORIENTATIONAL BEHAVIOUR

IF (Npart>0) THEN
    !-- OUTPUT ORIENTATIONAL DATA
    open(400,file='Nematic_average.out')

    WRITE (400,*) 
    WRITE (400,*) 'Ens. Average Nematic Director of Janus Nanodimers'
    WRITE (400,974) '<n_x>= ', NematicV(1)/dfloat(Nconfigs),' +- ',&
        dsqrt(NematicV2(1)/dfloat(Nconfigs)-(NematicV(1)/dfloat(Nconfigs))**2)
    WRITE (400,974) '<n_y>= ', NematicV(2)/dfloat(Nconfigs),' +- ',&
        dsqrt(NematicV2(2)/dfloat(Nconfigs)-(NematicV(2)/dfloat(Nconfigs))**2)
    WRITE (400,974) '<n_z>= ', NematicV(3)/dfloat(Nconfigs),' +- ',&
        dsqrt(NematicV2(3)/dfloat(Nconfigs)-(NematicV(3)/dfloat(Nconfigs))**2)
    WRITE (400,*) 'Nematic Order Parameter of Janus Nanodimers'
    WRITE (400,974) '<P_2>= ', P2/dfloat(Nconfigs),' +- ',&
        dsqrt(P4/dfloat(Nconfigs)-(P2/dfloat(Nconfigs))**2)


    close(400)

    open(100,file='Angle_average.dat')
    open(110,file='polar_X_prof.dat')
    open(120,file='polar_Y_prof.dat')
    open(130,file='polar_Z_prof.dat')
    open(210,file='azimutal_YZ_prof.dat')
    open(220,file='azimutal_ZX_prof.dat')
    open(230,file='azimutal_XY_prof.dat')
    open(310,file='abs_polar_X_prof.dat')
    open(320,file='abs_polar_Y_prof.dat')
    open(330,file='abs_polar_Z_prof.dat')

    DO i=1,3 !type
        IF (i==1) WRITE(100,*) '# dir(x=1,y=2,z=3) Polar angles(º) stdev(º)'
        IF (i==2) WRITE(100,*) '# plane(yz=1,zx=2,xy=3) Azimuthal angles(º) stdev(º)'
        IF (i==3) WRITE(100,*) '# dir(x=1,y=2,z=3) ABS(Polar angles)(º) stdev(º)'
        DO d=1,3 ! dims
            DO j=1,Nslabs
                acdirProf(i,d,0)=acdirProf(i,d,0)+acdirProf(i,d,j)
                acdirProf2(i,d,0)=acdirProf2(i,d,0)+acdirProf2(i,d,j)
                nsl=counterAngles(1,d,j)
                IF (i==2) nsl=counterAngles(2,d,j)
                IF (nsl>0) WRITE (100*i+10*d,600) &
                    real((dfloat(j)-0.5d0)/dfloat(Nslabs)),&
                    acdirProf(i,d,j)/dfloat(nsl),&
                    dsqrt(acdirProf2(i,d,j)/dfloat(nsl)-&
                    (acdirProf(i,d,j)/dfloat(nsl))**2)
            END DO !j Nslabs
            close(100*i+10*d)
            WRITE (100,601) &
                 d, acdirProf(i,d,0)/dfloat(Nconfigs*Npart),&
                 dsqrt(acdirProf2(i,d,0)/dfloat(Nconfigs*Npart)-&
                 (acdirProf(i,d,0)/dfloat(Nconfigs*Npart))**2)
         END DO !d dims
    END DO !i type

!*************
    acVdist=acVdist/dfloat(Npart*Nconfigs)
    open(100,file='polar_dist.dat')
    open(200,file='azimutal_dist.dat') 
    DO i=1,2
        DO j=1, int(i*180.d0/slabthicknessAngles)
            WRITE(100*i,602) real((dfloat(j)-0.5d0)*slabthicknessAngles)-i*90.d0,&
                real(acVdist(i,1,j)),real(acVdist(i,2,j)),real(acVdist(i,3,j))
        END DO !j Nslabs
        close(100*d)
    END DO !i Angle type
END IF !NPart>0
!***********************************************************************************************************************************************
!-- FORMATS
500 FORMAT (F8.5,2X,F10.5,2X,F10.5,2X,F10.5)
600 FORMAT (F10.5,2X,F12.5,2X,F12.5)
601 FORMAT (I2,2X,F12.5,2X,F12.5)
602 FORMAT (F7.1,2X,F10.5,2X,F10.5,2X,F10.5)
777 FORMAT(F14.6,2X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5)
974 FORMAT(A7,F8.5,A4,F8.5)
!***********************************************************************************************************************************************
!-- SUBROUTINES
CONTAINS
!-- DISTANCE BETWEEN TWO POINTS
SUBROUTINE DistAB(Ax,Ay,Az,Bx,By,Bz,Lx,Ly,Lz,ABx,ABy,ABz,AB2)
IMPLICIT NONE
REAL*8, INTENT(IN) :: Ax,Ay,Az,Bx,By,Bz,Lx,Ly,Lz
REAL*8, INTENT(OUT) :: ABx,ABy,ABz,AB2

ABx= Ax - Bx
 IF (ABx>Lx/2.d0)  ABx=ABx-Lx
 IF (ABx<-Lx/2.d0) ABx=ABx+Lx
ABy= Ay - By
 IF (ABy>Ly/2.d0)  ABy=ABy-Ly
 IF (ABy<-Ly/2.d0) ABy=ABy+Ly
ABz= Az - Bz
 IF (ABz>Lz/2.d0)  ABz=ABz-Lz
 IF (ABz<-Lz/2.d0) ABz=ABz+Lz
AB2=ABx*ABx+ABy*ABy+ABz*ABz
RETURN
END SUBROUTINE DistAB
!-- NEMATIC DIRECTOR COMPUTATION
SUBROUTINE Order(Tensorin,ex,ey,ez,OPAR)
IMPLICIT NONE
REAL*8, DIMENSION(3,3), INTENT(IN) :: Tensorin
REAL*8, INTENT(OUT) :: ex,ey,ez, OPAR
REAL*8 :: EIGVA1,EIGVA2,EIGVA3
REAL*8, DIMENSION(3,3) :: Eigenvectors, Tensor
INTEGER :: IP, NP
INTEGER :: IQ, NROT
REAL*8 ::  C,GOL, H, P, SM,SO, T,TAU, TRESH, THETA
REAL*8, DIMENSION(3) :: Dee
REAL*8, DIMENSION(100) :: B, Z

Eigenvectors=0.d0
Eigenvectors(1,1)=1.d0
Eigenvectors(2,2)=1.d0
Eigenvectors(3,3)=1.d0
Tensor=Tensorin
! Diagonalization

!     INITIALIZE D AND B TO THE DIAGONAL OF A
NP=3
DO IP=1,NP
    B(IP)=Tensor(IP,IP)
    Dee(IP)=B(IP)
    Z(IP)=0.0D00     
END DO !IP
NROT=0
DO I=1,50
    SM=0.0D00
!    SUM OFF-DIAGONAL ELEMENTS
    DO IP=1,NP-1
        DO IQ=IP+1,NP
            SM=SM+DABS(Tensor(IP,IQ))
        END DO !IQ
    END DO !IP 
!    CONDITION FOR MACHINE CONVERGENCE AND EXIT
    IF (SM.EQ.0.0D00) THEN
        EXIT
    END IF
    IF (I.LT.4) THEN
        TRESH=0.2D00*SM/DFLOAT(NP)**2
    ELSE
        TRESH=0.0D00
    END IF
    DO IP=1,NP-1
        DO IQ=IP+1,NP
            GOL=100.D00*DABS(Tensor(IP,IQ))
            IF ((I.GT.4).AND.(DABS(Dee(IP))+GOL.EQ.DABS(Dee(IP)))&
                .AND.(DABS(Dee(IQ))+GOL.EQ.DABS(Dee(IQ)))) THEN
                Tensor(IP,IQ)=0.0D00
            ELSE IF (DABS(Tensor(IP,IQ)).GT.TRESH) THEN
                H=Dee(IQ)-Dee(IP)
                IF (DABS(H)+GOL.EQ.DABS(H)) THEN
                    T=Tensor(IP,IQ)/H
                ELSE
                    THETA=0.5D00*H/Tensor(IP,IQ)
                    T=1.0D00/(DABS(THETA)+DSQRT(1.0D00+THETA**2))
                    IF (THETA.LT.0.0D00) T=-T
                END IF
                C=1.0D00/DSQRT(1.0D00+T*T)
                SO=T*C
                TAU=SO/(1.0D00+C)
                H=T*Tensor(IP,IQ)
                Z(IP)=Z(IP)-H
                Z(IQ)=Z(IQ)+H
                Dee(IP)=Dee(IP)-H
                Dee(IQ)=Dee(IQ)+H
                Tensor(IP,IQ)=0.0D00
                DO J=1,IP-1
                    GOL=Tensor(J,IP)
                    H=Tensor(J,IQ)
                    Tensor(J,IP)=GOL-SO*(H+GOL*TAU)
                    Tensor(J,IQ)=H+SO*(GOL-H*TAU)
                END DO !J
                DO J=IP+1,IQ-1
                    GOL=Tensor(IP,J)
                    H=Tensor(J,IQ)
                    Tensor(IP,J)=GOL-SO*(H+GOL*TAU)
                    Tensor(J,IQ)=H+SO*(GOL-H*TAU)
                END DO !J
                DO J=IQ+1,NP
                    GOL=Tensor(IP,J)
                    H=Tensor(IQ,J)
                    Tensor(IP,J)=GOL-SO*(H+GOL*TAU)
                    Tensor(IQ,J)=H+SO*(GOL-H*TAU)
                END DO !J
                DO J=1,NP
                    GOL=Eigenvectors(J,IP)
                    H=Eigenvectors(J,IQ)
                    Eigenvectors(J,IP)=GOL-SO*(H+GOL*TAU)
                    Eigenvectors(J,IQ)=H+SO*(GOL-H*TAU)
                END DO !J
                NROT=NROT+1
            END IF
        END DO !IQ
    END DO !IP
    DO IP=1,NP
        B(IP)=B(IP)+Z(IP)
        Dee(IP)=B(IP)
        Z(IP)=0.0D00
    END DO !IP
END DO !I

!      PAUSE '50 ITERATIONS SHOULD NEVER HAPPEN'
!    SORTS TO PUT THE EIGENVALUES INTO ASCENDING ORDER, AND TO 
!     REARRANGE THE COLUMNS OF V
DO I=1,NP-1
    K=I
    P=Dee(I)
    DO J=I+1,NP
        IF (Dee(J).GE.P) THEN
            K=J
            P=Dee(J)
        END IF
    END DO !J
    IF (K.NE.I) THEN
        Dee(K)=Dee(I)
        Dee(I)=P
        DO J=1,NP
            P=Eigenvectors(J,I)
            Eigenvectors(J,I)=Eigenvectors(J,K)
            Eigenvectors(J,K)=P
        END DO !J
    END IF
END DO !I
!     EIGENVALUES OF MATRIX Q THE ORDERING TENSOR (IN DECREASING ORDER)
!     Q=A-1/2
EIGVA1=3.0D00/2.0D00*Dee(1)-0.5D00
EIGVA2=3.0D00/2.0D00*Dee(2)-0.5D00
EIGVA3=3.0D00/2.0D00*Dee(3)-0.5D00
! OUTPUT
OPAR=EIGVA1
      !DIRECTOR (EIGENVECTOR Eigenvectors(I,1)), ASSOCIATED WITH EIGENVALUE 1
ex=Eigenvectors(1,1)
ey=Eigenvectors(2,1)
ez=Eigenvectors(3,1)

RETURN
END SUBROUTINE Order
!******************************************************************************

END PROGRAM LamellarStructure
