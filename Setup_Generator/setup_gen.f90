! **********************************************************************
! *** LAMELLAR DI-BLOCK COPOLYMER + NANODIMERS *************************
! *** Setup LAMMPS Script generator            *************************
! *** Author: J. Javier Burgos Mármol -- 2021  *************************
! **********************************************************************
!BSD 3-Clause License
!
!Copyright (c) 2021, J. Javier Burgos-Mármol (Github username: jjavier-bm)
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
!3. Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
PROGRAM Setup_script_gen
IMPLICIT NONE
!---DATA------------------------------------
!CONSTANTS
REAL*8, PARAMETER :: PI=3.141592654d0, TWOPI=6.283185307d0, HALFPI=1.570796327d0, PI6TH=0.523598776d0
!PARAMETERS
INTEGER, PARAMETER :: Nchains = 12696, LchainsA = 15, LchainsB = 15	!Number of Chains and length of each of the blocks
!REAL*8, PARAMETER :: DensChains=1.02d0  !Number density of chains desired in the final box (discounting the volume occupied by NPs)
REAL*8, PARAMETER :: Mbeadchains=1.d0 ! Mass of each bead in the chains
REAL*8, PARAMETER :: wpcpart=0.d0, NPdiam=3.d0 ! % Weight of Nano-dimers, Nano-dimer bead diameter 
REAL*8, PARAMETER :: epsRepulsive=52.d0, epsNeutral=1.d0 !Lennard-Jones energy parameters used
REAL*8, PARAMETER :: kFENEchains=30.d0, R0FENEchains=1.5d0, RbondEqChains=0.960897d0  !Chains stretching potential (FENE) parameters
REAL*8, PARAMETER :: Theta0=PI, ktheta=0.d0 !Chains bending potential. The form of the potential is 0.5*k*(x-x0)^2. !Note the 0.5 term
REAL*8, PARAMETER :: Rcutoff=2.d0**(1.d0/6.d0) ! Lennard-Jones cutoff distance for chain and nano-dimer beads
REAL*8, PARAMETER :: kFENENP=15.d0, R0FENENP=3.d0, RbondEqNP=2.5d0  !Nano-dimers stretching potential (FENE) parameters
!OTHER PARAMETERS
REAL*8, PARAMETER :: Vlens=(PI*(NPdiam-RbondEqNP)**2*(RbondEqNP**2+2.d0*RbondEqNP*NPdiam)/(12.d0*RbondEqNP)) ! Overlapping volume (lens shaped)
REAL*8, PARAMETER :: MbeadNpart=0.74d0*(PI6TH*NPdiam**3-Vlens/2.d0)/PI6TH ! Mass of Nano-dimer's beads (scalates as FCC packed), in Mbeadchains units.
INTEGER, PARAMETER :: Lchains=LchainsA+LchainsB
INTEGER, PARAMETER :: Npartopt=NINT(Nchains*Lchains*Mbeadchains*wpcpart/(2.d0*MbeadNpart*(100.d0-wpcpart))) !Number of nanodimers
!OTHER VARIABLES
REAL*8, DIMENSION(:) :: RChain(Nchains,Lchains,3) !Positions of chains
REAL*8, allocatable, DIMENSION(:) :: Rpart(:,:,:) !Positions of Nano-dimers
INTEGER, DIMENSION(:) :: iChain(Nchains,Lchains,3) !Box image of chains
INTEGER, allocatable, DIMENSION(:) :: ipart(:,:,:) !Box image of Nano-dimers
REAL*8, DIMENSION(:) :: Lbox(3) ! SIMULATION BOX
INTEGER, DIMENSION(:) :: Lfactors(3), Lpfactors(3) ! Multiplicative factors(nx,ny,nz) with nx*ny*nz=Lchains, nz even
REAL*8, DIMENSION(:) :: eps(4,4), REV(4,4) !Lennard-Jones energy and hard-core parameters
REAL*8, DIMENSION(:) :: sigma(4), Mass(4)
INTEGER :: ntypes,i,j,k,Npart,n,seed,z,Lprest
REAL*8 :: randnum
CHARACTER(len=64) :: fdata,fscript,fxyz
LOGICAL :: ok
LOGICAL, allocatable, DIMENSION(:) :: PartSlot(:,:)
INTEGER, allocatable, DIMENSION(:) :: blanks(:)
!---INIT------------------------------------
!print*, 'init'
sigma(1:2)=1.d0
sigma(3:4)=NPdiam
! These values must fulfil Lfactors(1)*Lfactors(2)*Lfactors(3)=Nchains
Lfactors(1)=46
Lfactors(2)=92
Lfactors(3)=3

! These values will be included in data and lammps script files
eps=epsNeutral
eps(1,2)=epsRepulsive
eps(3,3)=epsRepulsive
eps(4,4)=epsRepulsive
eps(3,4)=epsRepulsive
! The total number of particles will be adjusted to be distributed uniformally
!   in all interphases. If Lfactors(3) is too high, the final %wpc might be
!   significantly different from input wpcpart
IF (Npartopt>0) THEN
    IF (MOD(Npartopt,2*Lfactors(3))<2*Lfactors(3)) THEN
        Npart=Npartopt-MOD(Npartopt,2*Lfactors(3))
    ELSE
        Npart=Npartopt+2*Lfactors(3)-MOD(Npartopt,2*Lfactors(3))
    END IF
    !print*, 0.5d0*Vlens/(PI6TH*NPDiam**3), MbeadNpart
    allocate(Rpart(Npart,2,3))
    allocate(ipart(Npart,2,3))
    Lpfactors(3)=2*Lfactors(3)
    Lpfactors(1)=nint(dsqrt(dfloat(Npart*Lfactors(1))/dfloat(Lpfactors(3)*Lfactors(2))))
    Lpfactors(2)=nint(dsqrt(dfloat(Npart*Lfactors(1))/dfloat(Lpfactors(3)*Lfactors(2)))*&
                     dfloat(Lfactors(2))/dfloat(Lfactors(1)))
    DO WHILE (Lpfactors(1)*Lpfactors(2)*Lpfactors(3)<Npart)
        Lpfactors(2)=Lpfactors(2)+1
    END DO
! Since it is likely that not all slots in the lattice are to be filled with particles,
!   we need to randomly choose which ones to leave blank
    Lprest=Lpfactors(1)*Lpfactors(2)-Npart/Lpfactors(3)
    allocate(PartSlot(Lpfactors(3),Lpfactors(1)*Lpfactors(2)))
    PartSlot=.True.
    IF (Lprest>0) THEN
        allocate(blanks(Lprest))
        !print*, 'Lprest=',Lprest,', Lpf1*Lpf2=',Lpfactors(1)*Lpfactors(2)
        DO z=1,Lpfactors(3)
            DO i=1,Lprest
                seed=((z-1)*Lprest+i)*17+53
                Call Ran01(seed,randnum)
                blanks(i)=int(randnum*Lpfactors(1)*Lpfactors(2))+1
                ok=.False.
                DO WHILE (.not.ok)           
                    ok=.True.
                    DO j=1, i-1
                        IF (blanks(i)==blanks(j)) THEN
                            !print*, 'repeated',z,blanks(i)
                            ok=.False.
                            seed=seed+19
                            Call Ran01(seed,randnum)
                            blanks(i)=int(randnum*Lpfactors(1)*Lpfactors(2))+1
                        END IF
                    END DO !j
                END DO !ok
                !print*, 'Partslot ', z,blanks(i), randnum
                PartSlot(z,blanks(i))=.False.
            END DO !i
        END DO !z 
    END IF
ELSE
    Npart=0
END IF
!print*, "Number of particles = ",Npart, Npartopt
IF (Npart==0) THEN
    ntypes=2
ELSE
    ntypes=4
END IF

DO i=1,ntypes
    IF (i<=2) THEN
        Mass(i)=Mbeadchains
    ELSE
        Mass(i)=MbeadNpart
    END IF
    DO j=1,ntypes
        IF (j>i) eps(j,i)=eps(i,j)
        REV(i,j)=(sigma(i)+sigma(j))/2.d0-1.d0
    END DO
END DO
Call files(wpcpart,fdata,fscript,fxyz)
!---PLACING---------------------------------
!print*, 'placing'
Call ChainPlacing(Nchains,Lchains,RbondEqChains,REV(1,1)+Rcutoff,Lfactors,RChain,ichain)
!print*, 'chains done'
Lbox(1)=RChain(Lfactors(3)*Lfactors(1),Lchains,1)+RbondEqChains/2.d0
Lbox(2)=RChain(Nchains,Lchains,2)-RChain(1,1,2)+REV(1,1)+Rcutoff
Lbox(3)=RChain(Lfactors(3),Lchains,3)+RbondEqChains/2.d0
IF (Npart>0) THEN
    Lbox(2)=Lbox(2)-(REV(1,1)+Rcutoff)+2.d0*(REV(1,3)+Rcutoff)+&
                    (Lpfactors(2)-1)*(REV(3,3)+Rcutoff)
    Call PartPlacing(Npart,Lpfactors, Lfactors,PartSlot, Lbox, &
                     RbondEqNP,REV(3,3)+Rcutoff,REV(1,3)+Rcutoff, Rpart,ipart)
    Call ChainRealloc(int(Lfactors(2)/Lpfactors(2)),REV(3,3)+Rcutoff,&
                      Lfactors(2),Lfactors(1)*Lfactors(3),Lpfactors(2))
    !print*, 'LboxY: ', Lbox(2), Rpart(Npart,2,2)
    !print*, 'particles done'
END IF
!---LAMMPS DATA FILE------------------------
open(9000, file=trim(fdata))

WRITE(9000,*) 'LAMMPS Description'
WRITE(9000,*)
WRITE(9000,*) Nchains*Lchains+2*Npart, 'atoms'
WRITE(9000,*) Nchains*(Lchains-1)+Npart, 'bonds'
IF (ktheta>0.d0) WRITE(9000,*) Nchains*(Lchains-2), 'angles'
WRITE(9000,*)
WRITE(9000,*) ntypes, 'atom types'
WRITE(9000,*) ntypes/2, 'bond types'
IF (ktheta>0.d0) WRITE(9000,*) 1, 'angle types'
WRITE(9000,*)
WRITE(9000,1001) -Lbox(1)/2.d0,Lbox(1)/2.d0,'xlo xhi'
WRITE(9000,1001) -Lbox(2)/2.d0,Lbox(2)/2.d0,'ylo yhi'
WRITE(9000,1001) -Lbox(3)/2.d0,Lbox(3)/2.d0,'zlo zhi'
WRITE(9000,*)
WRITE(9000,*) ' Masses'
WRITE(9000,*)
DO i=1,ntypes
    WRITE(9000,3001) i, Mass(i)
END DO
WRITE(9000,*)
WRITE(9000,*) ' Pair Coeffs'
WRITE(9000,*)
DO i=1,ntypes
    WRITE(9000,1002) i, eps(i,i), 1.d0, REV(i,i), Rcutoff
END DO
WRITE(9000,*)
WRITE(9000,*) ' PairIJ Coeffs'
WRITE(9000,*)
DO i=1,ntypes
    DO j=i,ntypes
        WRITE(9000,1003) i, j, eps(i,j), 1.d0, REV(i,j), Rcutoff
    END DO
END DO
WRITE(9000,*)
WRITE(9000,*) ' Bond Coeffs'
WRITE(9000,*)
WRITE(9000,1004) 1, kFENEchains, R0FENEchains, 1.d0, 1.d0
IF (Npart>0) WRITE(9000,1004) 2, kFENENP, R0FENENP, 1.d0, NPdiam
WRITE(9000,*)
IF (ktheta>0.d0) THEN
    WRITE(9000,*)
    WRITE(9000,*) ' Angle Coeffs'
    WRITE(9000,*)
    WRITE(9000,1005) 1, 0.5d0*ktheta, Theta0
END IF
WRITE(9000,*)
WRITE(9000,*) ' Atoms'
WRITE(9000,*)
DO i=1, Nchains
    DO j=1,Lchains
        IF (j<=LchainsA) THEN            
            WRITE(9000,1006) (i-1)*Lchains +j, i,1,&
                 Rchain(i,j,1)-Lbox(1)/2.d0,Rchain(i,j,2)-Lbox(2)/2.d0,&
                 Rchain(i,j,3)-Lbox(3)/2.d0,ichain(i,j,1),&
                 ichain(i,j,2),ichain(i,j,3)
        ELSE
            WRITE(9000,1006) (i-1)*Lchains +j, i,2,&
                 Rchain(i,j,1)-Lbox(1)/2.d0,Rchain(i,j,2)-Lbox(2)/2.d0,&
                 Rchain(i,j,3)-Lbox(3)/2.d0,ichain(i,j,1),&
                 ichain(i,j,2),ichain(i,j,3)
        END IF
    END DO
END DO
IF (Npart>0) THEN
    DO i=1,Npart
        DO j=1,2     
            WRITE(9000,1006) Nchains*Lchains+(i-1)*2 +j, Nchains+i,j+2,&
                 Rpart(i,j,1)-Lbox(1)/2.d0,Rpart(i,j,2)-Lbox(2)/2.d0,&
                 Rpart(i,j,3)-Lbox(3)/2.d0,ipart(i,j,1),&
                 ipart(i,j,2),ipart(i,j,3)
        END DO
    END DO
END IF
WRITE(9000,*)
WRITE(9000,*) ' Bonds'
WRITE(9000,*)
DO i=1,Nchains
     DO j=1,Lchains-1
         WRITE(9000,1007) (i-1)*(Lchains-1) +j, 1,(i-1)*Lchains+j,(i-1)*Lchains+j+1
     END DO
END DO
IF (Npart>0) THEN
    DO i=1,Npart
        WRITE(9000,1007) Nchains*(Lchains-1)+i, 2,Nchains*Lchains+(i-1)*2 +1,Nchains*Lchains+(i-1)*2 +2
    END DO
END IF
IF (ktheta>0.d0) THEN
    WRITE(9000,*)
    WRITE(9000,*) ' Angles'
    WRITE(9000,*)
    DO i=1,Nchains
        DO j=1,Lchains-2
            WRITE(9000,1008) (i-1)*(Lchains-2) +j, 1,&
                (i-1)*Lchains+j,(i-1)*Lchains+j+1,(i-1)*Lchains+j+2
        END DO
    END DO
END IF
WRITE(9000,*)
WRITE(9000,*) ' Velocities'
WRITE(9000,*)
DO i=1, NChains*Lchains+Npart*2
    WRITE(9000,1009) i, 0.d0,0.d0,0.d0
END DO
 close(9000)

!---LAMMPS SCRIPT TEMPLATE FILE-------------
open(9001, file=trim(fscript))

WRITE(9001,*) '# Polymer nanocomposite'
WRITE(9001,*) '# Number of block-copolimer chains', Nchains
WRITE(9001,*) '# Number of beads in block A(1)', LchainsA
WRITE(9001,*) '# Number of beads in block B(2)', LchainsB
WRITE(9001,*) '# Number of Janus Nano-dimers', Npart
WRITE(9001,*)  
WRITE(9001,*) '# TEMPLATE SCRIPT'
WRITE(9001,*)
WRITE(9001,*) 'units           lj'
WRITE(9001,*) 'atom_style      molecular'
WRITE(9001,*) 'pair_style      lj/expand ', Rcutoff
WRITE(9001,*) 'pair_modify     shift yes'
WRITE(9001,*) 'bond_style      fene'
IF (KTheta==0.d0) WRITE(9001,*) 'angle_style     none'
IF (KTheta>0.d0) WRITE(9001,*) 'angle_style     harmonic'
WRITE(9001,*) 'dihedral_style  none'
WRITE(9001,*) 'improper_style  none'
WRITE(9001,*) 'boundary        p p p'
WRITE(9001,*) '# ----------'
WRITE(9001,*) 'read_data       ',trim(fdata)
WRITE(9001,*) 'special_bonds   fene'
WRITE(9001,*)
DO i=1,ntypes
    DO j=i,ntypes
        WRITE(9001,2001) 'pair_coeff      ',i,j,eps(i,j), 1.d0, REV(i,j), Rcutoff
    END DO
END DO
WRITE(9001,2002) 'bond_coeff      ',1, kFENEchains, R0FENEchains, 1.d0, 1.d0
IF (Npart>0) WRITE(9001,2002) 'bond_coeff      ',2, kFENENP, R0FENENP, 1.d0, NPdiam
IF (KTheta>0.d0) WRITE(9001,2003) 'angle_coeff     ',1, 0.5d0*ktheta, Theta0

WRITE(9001,*)
WRITE(9001,*) 'timestep        0.001'
WRITE(9001,*) 'thermo          1000'
WRITE(9001,*) 'thermo_style    custom step temp etotal density press lx ly lz vol'
WRITE(9001,*) 'group           chains type 1 2'
IF (Npart>0) WRITE(9001,*) 'group           nps type 3 4'
WRITE(9001,*)
WRITE(9001,*) 'dump            snapshot all xyz 10000 template.xyz'
IF (Npart==0) WRITE(9001,*) 'dump_modify     snapshot element C N'
IF (Npart>0) WRITE(9001,*) 'dump_modify     snapshot element C N Rb Ba'
WRITE(9001,*) '###################################'
WRITE(9001,*) 'fix             1 all nvt temp 0.01 0.01 0.1 tchain 1'
WRITE(9001,*) 'run             0'
WRITE(9001,*) '#unfix           1'
WRITE(9001,*) '###################################'
WRITE(9001,*) '#fix             1 all npt temp 0.01 1.0 0.1 tchain 1 aniso 14.25 14.25 1.0 couple xy'
WRITE(9001,*) '#run             1000000'
 close(9001)
!---XYZ INITIAL COORDINATES-------------
open(9002, file=trim(fxyz))
WRITE(9002,*) Nchains*Lchains+Npart*2
WRITE(9002,*) 'Atoms. Timestep: 0'

DO i=1,Nchains
    DO j=1,Lchains
        IF (j<=LchainsA) THEN
            WRITE(9002,4001) 'C', RChain(i,j,1)-Lbox(1)/2.d0,&
                                  RChain(i,j,2)-Lbox(2)/2.d0,&
                                  RChain(i,j,3)-Lbox(3)/2.d0 
        ELSE
            WRITE(9002,4001) 'N', RChain(i,j,1)-Lbox(1)/2.d0,&
                                  RChain(i,j,2)-Lbox(2)/2.d0,&
                                  RChain(i,j,3)-Lbox(3)/2.d0
        END IF
    END DO !j
END DO !i
DO i=1,Npart
    DO j=1,2
        IF (j==1) THEN
            WRITE(9002,4002) 'Rb', Rpart(i,j,1)-Lbox(1)/2.d0,&
                                   Rpart(i,j,2)-Lbox(2)/2.d0,&
                                   Rpart(i,j,3)-Lbox(3)/2.d0
        ELSE
            WRITE(9002,4002) 'Ba', Rpart(i,j,1)-Lbox(1)/2.d0,&
                                   Rpart(i,j,2)-Lbox(2)/2.d0,&
                                   Rpart(i,j,3)-Lbox(3)/2.d0
        END IF
    END DO !j
END DO !i
 close(9002)
!---FORMATS-------------
1001 FORMAT(2X,F12.8,2X,F12.8,2X,A7)
1002 FORMAT(2X,I2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F11.8)
1003 FORMAT(2X,I2,2X,I2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F11.8)
1004 FORMAT(2X,I2,2X,F6.2,2X,F5.2,2X,F5.2,2X,F5.2)
1005 FORMAT(2X,I2,2X,F6.2,2X,F6.2)
1006 FORMAT(2X,I9,2X,I9,2X,I2,2X,F12.8,2X,F12.8,2X,F12.8,2X,I2,2X,I2,2X,I2)
1007 FORMAT(2X,I9,2X,I2,2X,I9,2X,I9)
1008 FORMAT(2X,I9,2X,I2,2X,I9,2X,I9,2X,I9)
1009 FORMAT(2X,I9,2X,F11.8,2X,F11.8,2X,F11.8)
2001 FORMAT(A16,2X,I2,2X,I2,2X,F5.2,2X,F5.2,2X,F5.2,2X,F11.8)
2002 FORMAT(A16,2X,I2,2X,F6.2,2X,F5.2,2X,F5.2,2X,F5.2)
2003 FORMAT(A16,2X,I2,2X,F6.2,2X,F6.2)
3001 FORMAT(2X,I2,2X,F5.2)
4001 FORMAT(A1,1X,F9.4,1X,F9.4,1X,F9.4)
4002 FORMAT(A2,1X,F9.4,1X,F9.4,1X,F9.4)
!---SUBROUTINES-------------
Contains
!*****************************************
SUBROUTINE ChainPlacing(N,L,bond,rdist,f,R,iR)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, L
REAL*8, INTENT(IN) :: bond, rdist
INTEGER, DIMENSION(:), INTENT(IN) :: f(3)
REAL*8, DIMENSION(:), INTENT(OUT) :: R(N,L,3)
INTEGER, DIMENSION(:), INTENT(OUT) :: iR(N,L,3)
INTEGER :: dir, nx,ny, nz, nc,i,j,k
REAL*8, DIMENSION(:) :: r0(3)

dir=1
r0=rdist/2.d0
ir=0
nz=0
DO i=1,f(3)
    DO j=1,L
        DO k=1,3
            R(i,j,k)=r0(k)
        END DO
        IF (j<L) THEN
            r0(3)=r0(3)+dir*bond
        ELSE
            r0(3)=r0(3)+dir*rdist
        END IF
    END DO
END DO
!print*, 'First ', f(3), ' chains placed'
dir=-dir
r0(1)=r0(1)+rdist
r0(3)=R(f(3),L/2,3)
DO i=f(3)+1, 2*f(3)
    DO j=1,L
        DO k=1,3
            R(i,j,k)=r0(k)
            IF (k==3) iR(i,j,k)=nz
        END DO
        IF (j<L) THEN
            r0(3)=r0(3)+dir*bond
        ELSE
            r0(3)=r0(3)+dir*rdist 
        END IF
        IF (r0(3)<0.d0) THEN
            r0(3)=R(f(3),L,3)
            nz=nz-1
        END IF
    END DO
END DO
!print*, 'Second 3 chains placed'

nc=2*f(3)
DO nx=1,f(1)/2-1
    DO i=1,2*f(3)
        nc=nc+1 
        DO j=1,L
            R(nc,j,1)=R(i,j,1)+2*nx*rdist
            R(nc,j,2)=R(i,j,2)
            R(nc,j,3)=R(i,j,3)
            iR(nc,j,3)=iR(i,j,3)
        END DO
    END DO
END DO
!print*, 'expansion on X dir Finished', nc
DO ny=1,f(2)-1
    DO i=1,f(3)*f(1)
        nc=nc+1
        DO j=1,L
            R(nc,j,1)=R(i,j,1)
            R(nc,j,2)=R(i,j,2)+ny*rdist
            R(nc,j,3)=R(i,j,3)
            iR(nc,j,3)=iR(i,j,3)
        END DO
    END DO
END DO
!print*, 'expansion on Y dir Finished', nc

RETURN
END SUBROUTINE ChainPlacing
!*****************************************
SUBROUTINE PartPlacing(N,f,fch,place,Lb,bond,distNP,distChNP,R,iR)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER, DIMENSION(:), INTENT(IN) :: f(3),fch(3)
LOGICAL, DIMENSION(:), INTENT(IN) :: place(f(3),f(1)*f(2))
REAL*8, DIMENSION(:), INTENT(IN) :: Lb(3)
REAL*8, INTENT(IN) :: bond,distNP,distChNP
REAL*8, DIMENSION(:), INTENT(OUT) :: R(N,2,3)
INTEGER, DIMENSION(:), INTENT(OUT) :: iR(N,2,3)
INTEGER :: nf,i,j,k,ninter,np,nx,ny,nz, nn,nchperlayer
REAL*8, DIMENSION(:) :: r0(3)
REAL*8 :: distCh

!print*, 'particle placing'

distCh=2.d0*distChNp-distNP
iR=0
r0(3)=0.001d0
r0(1)=distNP/2.d0
nchperlayer=int(fch(2)/f(2))
r0(2)=(nchperlayer-0.5d0)*distCh+distChNP
IF (nint(dfloat(f(2))/dfloat(MOD(fch(2),f(2))))==1) r0(2)=r0(2)+distCh

!print*, 'Ly=',Lb(2),', f(2)=',f(2)
!print*, f

np=0
DO ny=1,f(2)
    DO nx=1,f(1)
        nn=(ny-1)*f(1)+nx
        DO nz=1,f(3)
            IF (place(nz,nn)) THEN
                np=np+1
                DO j=1,2
                    R(np,j,1)=r0(1)
                    R(np,j,2)=r0(2)
                    IF (j==1) THEN 
                        R(np,j,3)=r0(3)-bond/2.d0
                        IF (R(np,j,3)<0.d0) THEN
                            R(np,j,3)=R(np,j,3)+Lb(3)
                            iR(np,j,3)=-1
                        END IF
                    ELSE
                        R(np,j,3)=r0(3)+bond/2.d0
                    END IF
                END DO !j
            END IF
            IF (nz<f(3)) THEN
                r0(3)=r0(3)+Lb(3)/f(3)
            ELSE
                r0(3)=0.001d0
                IF (nx<f(1)) THEN
                    r0(1)=r0(1)+Lb(1)/f(1)
                ELSE
                    r0(1)=distNP/2.d0
                    r0(2)=r0(2)+2.d0*distChNP+(nchperlayer-1)*distCh
                    DO k=1,f(2)
                        IF (ny==nint(k*dfloat(f(2))/dfloat(MOD(fch(2),f(2))))-1) r0(2)=r0(2)+distCh
                    END DO !k
                END IF
            END IF
        END DO !nz
    END DO !nx
END DO !ny

RETURN
END SUBROUTINE PartPlacing
!*****************************************
SUBROUTINE ChainRealloc(ncpl,dnp,fcy,fcxz,fpy)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ncpl,fcy, fcxz, fpy
REAL*8, INTENT(IN) :: dnp
INTEGER :: a,b,nc,lc,np,nl
np=0
!print*, 'chains per y-layer = ',ncpl

nl=0
DO a=1,fpy
    nl=nl+ncpl*fcxz
    DO b=1,fpy
        IF (a==nint(dfloat(b*fpy)/dfloat(MOD(fcy,fpy)))) nl=nl+fcxz
    END DO
    DO nc=nl+1,Nchains
        DO lc=1,Lchains
            RChain(nc,lc,2)=RChain(nc,lc,2)+dnp
        END DO !lc
    END DO !nc
END DO !a
    
RETURN
END SUBROUTINE ChainRealloc
!*****************************************
SUBROUTINE files(w,f1,f2,f3)
IMPLICIT NONE
REAL*8, INTENT(IN) :: w
CHARACTER(len=64), INTENT(OUT) :: f1,f2,f3

IF (w==0.d0) f1='data.pure.initial.lammpsdata'
IF (w<1.d0 .and. w>0.d0) write (f1, "(A9,I1,A23)") "data.wpc0", INT(w*10),".N00.initial.lammpsdata"
IF (w<10.d0 .and. w>=1.d0) write (f1, "(A8,I1,A23)") "data.wpc", INT(w),".N00.initial.lammpsdata"
IF (w>=10.d0) write (f1, "(A8,I2,A23)") "data.wpc", INT(w),".N00.initial.lammpsdata"

IF (w==0.d0) f2='in.pure.initial.lammpsscript'
IF (w<1.d0 .and. w>0.d0) write (f2, "(A7,I1,A24)") "in.wpc0", INT(w*10),".N00.initial.lammpsscript"
IF (w<10.d0 .and. w>=1.d0) write (f2, "(A6,I1,A24)") "in.wpc", INT(w),".N00.initial.lammpsscript"
IF (w>=10.d0) write (f2, "(A6,I2,A24)") "in.wpc", INT(w),".N00.initial.lammpsscript"

IF (w==0.d0) f3='traj.pure.initial.xyz'
IF (w<1.d0 .and. w>0.d0) write (f3, "(A9,I1,A16)") "traj.wpc0", INT(w*10),".N00.initial.xyz"
IF (w<10.d0 .and. w>=1.d0) write (f3, "(A8,I1,A16)") "traj.wpc", INT(w),".N00.initial.xyz"
IF (w>=10.d0) write (f3, "(A8,I2,A16)") "traj.wpc", INT(w),".N00.initial.xyz"

RETURN
END SUBROUTINE
!*****************************************
  SUBROUTINE Ran01(inseed,variable)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: inseed
  REAL*8, INTENT(OUT) :: variable
  integer(kind=4), dimension(8) :: val
  INTEGER :: idummy

   call date_and_time(values=val)  ! random numbers are linked to time
   idummy=val(3)+val(6)+val(7)+val(8)+int(inseed,4)
   variable=ran2(idummy)
  RETURN
  END SUBROUTINE Ran01
!-----------------------------------------
  real function ran2(idum) ! NUMERICAL RECIPES

    integer  idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real  AM,EPS,RNMX

    parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, & 
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,  &
         IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

    integer  idum2,j,k,iv(NTAB),iy
    save iv,iy,idum2
    data idum2/123456789/,iv/NTAB*0/,iy/0/

    if(idum.le.0) then
       idum=max(-idum,1)
       idum2=idum

       do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if(idum.lt.0) idum=idum+IM1
          if(j.le.NTAB) iv(j)=idum
       end do

       iy=iv(1)
    end if

    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if(idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if(idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1) iy=iy+IMM1
    ran2=min(AM*iy,RNMX)

    return
  end function ran2
END PROGRAM Setup_script_gen
