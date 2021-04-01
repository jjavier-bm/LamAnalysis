! **********************************************************************
! *** LAMELLAR DI-BLOCK COPOLYMER + NANODIMERS *************************
! *** Dynamical Data Analysis                  *************************
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
PROGRAM LamellarDynamics
IMPLICIT NONE
!---DATA------------------------------------
!CONSTANTS
REAL*8, PARAMETER :: PI=3.141592654d0, TWOPI=6.283185307d0, HALFPI=1.570796327d0
!PARAMETERS
INTEGER, PARAMETER :: Nchains = 12696, LchainsA = 15, LchainsB = 15	!Number of Chains and length of each of the blocks
REAL*8, PARAMETER :: Mbeadchains=1.d0 ! Mass of each bead in the chains
INTEGER, PARAMETER :: Npart = 1080 !Number of nano-dimers
REAL*8, PARAMETER :: NPdiam=3.d0, RbondEqNP=2.5d0 ! NP bead diameter, NP bond eq. distance (mass derived from here)
REAL*8, PARAMETER :: delta_t=0.001d0
INTEGER, PARAMETER :: normaldir=3
INTEGER, PARAMETER :: NtimestepsperFile=1000000                 !Number of timesteps per file (power of 10)
INTEGER, PARAMETER :: Nfiles=100                                !Number of statistical repetitions (files, power of 10)
!OTHER PARAMETERS
INTEGER, PARAMETER :: Lchains=LchainsA+LchainsB
REAL*8, PARAMETER :: Vlens=(PI*(NPdiam-RbondEqNP)**2*(RbondEqNP**2+2.d0*RbondEqNP*NPdiam)/(12.d0*RbondEqNP)) ! Overlapping volume (lens shaped)
REAL*8, PARAMETER :: MbeadNpart=0.74d0*(NPdiam**3/6.d0-Vlens/2.d0) ! Mass of Nano-dimer's beads (scalates as FCC packed)
REAL*8, PARAMETER :: Mass=dfloat(NChains*Lchains)*Mbeadchains+dfloat(2*Npart)*MbeadNpart
INTEGER, PARAMETER :: finaltimestep=NtimestepsperFile*Nfiles !Final timestep printed (in a log scale)
INTEGER, PARAMETER :: Nconfigs=int(dlog10(dfloat(finaltimestep)))*9&
                               +finaltimestep/(10**(int(dlog10(dfloat(finaltimestep))))) !Total number of configurations
INTEGER, PARAMETER :: NconfigsperFile=int(dlog10(dfloat(NtimestepsperFile)))*9&
                                      +NtimestepsperFile/(10**(int(dlog10(dfloat(NtimestepsperFile))))) !Number of configurations per file
!OTHER VARIABLES
REAL*8, DIMENSION(:) :: RChain(Nchains,0:Lchains,3) !Positions of chains ! 0 is for the CoM
REAL*8, DIMENSION(:) :: Rpart(Npart,0:2,3) !Positions of Nano-dimers ! 0 is for the CoM
REAL*8, DIMENSION(:) :: R0Chain(Nchains,0:Lchains,3,0:Nfiles) !Positions of chains at t=0! 0 is for the CoM
REAL*8, DIMENSION(:) :: R0part(Npart,0:2,3,0:Nfiles) !Positions of Nano-dimers at t=0! 0 is for the CoM
REAL*8, DIMENSION(:) :: RCoM(3) !Position of Centre of Mass
REAL*8, DIMENSION(:) :: R0CoM(3,0:Nfiles) !Position of Centre of Mass at t=0
REAL*8, DIMENSION(:) :: pos(3), vel(3), jvol(6) ! jvol(stresstensor*volume): xx,yy,zz,xy,xz,yz
INTEGER, DIMENSION(:) :: t(0:Nconfigs), avconf(Nconfigs)
INTEGER ::  i, j, k, l, conf, conf0, d, d2, f, f2, id, idmol, idbead, tipo
REAL*8, DIMENSION(:) :: msd(Nconfigs,6,0:4) ! 1=ChainMonomers, 2=Chain Centre beads, 3=Chain End beads, 4=Chain CoMs, 5=Particle Beads, 6=Particle CoMs; 0=3D, 1=x,2=y,3=z,4=2D plane
REAL*8, DIMENSION(:) :: avdir(0:Nconfigs), integral(Nconfigs)
!***********************************************************************************************************************************************
R0chain=0.d0
R0part=0.d0
R0CoM=0.d0
msd=0.d0
avdir=0.d0
avdir(0)=1.d0
!print*, Nconfigs, NConfigsPerFile
DO conf=0,Nconfigs
    IF (conf<10) THEN
        t(conf)=conf
    ELSE
        t(conf)=t(conf-9)*10
    END IF
    IF (conf>=0 .and. conf<=NconfigsPerFile) avconf(conf)=Nfiles
    IF (conf>NconfigsperFile) avconf(conf)=0
END DO !i N configs
DO f=1,Nfiles
    open(100+f)
    !print*, 'Reading file fort.',100+f
    IF (f==1) conf0=0
    IF (f>1) conf0=1
    DO conf=conf0, NconfigsperFile
        !print*, conf
        Rchain=0.d0
        Rpart=0.d0
        RCoM=0.d0
        DO i=1,9
            READ(100+f,*)
        END DO !i
        ! READ THE COORDINATES
        DO i=1,Nchains*Lchains+Npart*2
            READ(100+f,*) id, pos(1), pos(2), pos(3), vel(1), vel(2), vel(3),&
                         jvol(1), jvol(2), jvol(3), jvol(4), jvol(5), jvol(6)
            IF (id<=Nchains*Lchains) THEN
                idmol=int((id-1)/Lchains)+1
                idbead=id-(idmol-1)*Lchains
                DO d=1,3
                    Rchain(idmol,idbead,d)=pos(d)
                    Rchain(idmol,0,d)=Rchain(idmol,0,d)+pos(d)*MbeadChains
                    RCoM(d)=RCoM(d)+pos(d)*MbeadChains
                END DO !dims
            ELSE
                idmol=int((id-Nchains*Lchains-1)/2)+1
                idbead=id-Nchains*Lchains-(idmol-1)*2
                DO d=1,3
                    Rpart(idmol,idbead,d)=pos(d)
                    Rpart(idmol,0,d)=Rpart(idmol,0,d)+pos(d)*MbeadNpart
                    RCoM(d)=RCoM(d)+pos(d)*MbeadNpart
                END DO !dims
            END IF !id         
        END DO !i
        ! COMPUTE CoM
        DO d=1,3
            DO i=1,Nchains
                Rchain(i,0,d)=Rchain(i,0,d)/(dfloat(Lchains)*MbeadChains)
                IF (conf==0 .or. conf==NconfigsperFile) THEN
                    DO j=0,Lchains
                        IF (conf==0) R0chain(i,j,d,0)=Rchain(i,j,d)
                        IF (conf==NconfigsperFile) R0chain(i,j,d,f)=Rchain(i,j,d)
                    END DO !j Lchains
                END IF !terminal confs
            END DO !i Nchains
            DO i=1,Npart
                Rpart(i,0,d)=Rpart(i,0,d)/(2.d0*MbeadNpart)
                IF (conf==0 .or. conf==NconfigsperFile) THEN
                    DO j=0,2
                        IF (conf==0) R0part(i,j,d,0)=Rpart(i,j,d)
                        IF (conf==NconfigsperFile) R0part(i,j,d,f)=Rpart(i,j,d)
                    END DO !j dimer beads and CoM
                END IF !terminal confs
            END DO !i Npart
            RCoM(d)=RCoM(d)/Mass
            IF (conf==0) R0CoM(d,0)=RCoM(d)
            IF (conf==NconfigsperFile) R0CoM(d,f)=RCoM(d)
        END DO !dims    
        IF (conf>0) THEN
            ! Compute MSD and director
            DO i=1,Nchains
                DO d=1,3
                    DO j=1,Lchains
                        msd(conf,1,d)=msd(conf,1,d)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !All chain beads, 1D
                        msd(conf,1,0)=msd(conf,1,0)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 ! 3D
                        IF (d/=normaldir) msd(conf,1,4)=msd(conf,1,4)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                        IF (j==LchainsA .or. j==LchainsA+1) THEN !Central beads
                            msd(conf,2,d)=msd(conf,2,d)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !1D
                            msd(conf,2,0)=msd(conf,2,0)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !3D
                            IF (d/=normaldir) msd(conf,2,4)=msd(conf,2,4)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                        ELSE IF (j==1 .or. j==Lchains) THEN !End beads
                            msd(conf,3,d)=msd(conf,3,d)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !1D
                            msd(conf,3,0)=msd(conf,3,0)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !3D
                            IF (d/=normaldir) msd(conf,3,4)=msd(conf,3,4)+(Rchain(i,j,d)-R0chain(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                        END IF!bead type
                    END DO !j Lchains
                    msd(conf,4,d)=msd(conf,4,d)+(Rchain(i,0,d)-R0chain(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !Chain CoMs, 1D
                    msd(conf,4,0)=msd(conf,4,0)+(Rchain(i,0,d)-R0chain(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !3D
                    IF (d/=normaldir) msd(conf,4,4)=msd(conf,4,4)+(Rchain(i,0,d)-R0chain(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                END DO !dims
            END DO !i Nchains
            DO i=1,Npart
                DO d=1,3
                    DO j=1,2
                        msd(conf,5,d)=msd(conf,5,d)+(Rpart(i,j,d)-R0part(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !All ND beads, 1D
                        msd(conf,5,0)=msd(conf,5,0)+(Rpart(i,j,d)-R0part(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 ! 3D
                        IF (d/=normaldir) msd(conf,5,4)=msd(conf,5,4)+(Rpart(i,j,d)-R0part(i,j,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                    END DO !dimer bead
                    msd(conf,6,d)=msd(conf,6,d)+(Rpart(i,0,d)-R0part(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !NP CoMs, 1D
                    msd(conf,6,0)=msd(conf,6,0)+(Rpart(i,0,d)-R0part(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !3D
                    IF (d/=normaldir) msd(conf,6,4)=msd(conf,6,4)+(Rpart(i,0,d)-R0part(i,0,d,f-1)-RCoM(d)+R0CoM(d,f-1))**2 !2D
                    avdir(conf)=avdir(conf)+(Rpart(i,2,d)-Rpart(i,1,d))*(R0part(i,2,d,f-1)-R0part(i,1,d,f-1))&
                                           /(dsqrt((Rpart(i,2,1)-Rpart(i,1,1))**2&
                                                  +(Rpart(i,2,2)-Rpart(i,1,2))**2&
                                                  +(Rpart(i,2,3)-Rpart(i,1,3))**2)&
                                            *dsqrt((R0part(i,2,1,f-1)-R0part(i,1,1,f-1))**2&
                                                  +(R0part(i,2,2,f-1)-R0part(i,1,2,f-1))**2&
                                                  +(R0part(i,2,3,f-1)-R0part(i,1,3,f-1))**2))
                END DO !dims
            END DO !Npart
        END IF !conf
    END DO !conf
    close(100+f)
    !print*, '    --> file read'
    !print*, 'Accumulators'
    ! Multi-file statistics (only for conf>NconfigsperFile)
    DO f2=1,f-1
        conf=0
        DO l=NconfigsperFile+1,Nconfigs
            IF ((f-(f2-1))*NtimestepsperFile==t(l)) conf=l
        END DO !l loop
        IF (conf>0) THEN
            avconf(conf)=avconf(conf)+1
            DO d=1,3
                DO i=1,Nchains 
                    DO j=1,Lchains
                        msd(conf,1,d)=msd(conf,1,d)+&
                                    +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                    -RCoM(d)+R0CoM(d,f2-1))**2 !All chain beads, 1D
                        msd(conf,1,0)=msd(conf,1,0)+&
                                    +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                    -RCoM(d)+R0CoM(d,f2-1))**2 !All chain beads, 3D
                        IF (d/=normaldir) msd(conf,1,4)=msd(conf,1,4)+&
                                    +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                    -RCoM(d)+R0CoM(d,f2-1))**2 !All chain beads, 2D
                        IF (j==LchainsA .or. j==LchainsA+1) THEN !Central beads
                            msd(conf,2,d)=msd(conf,2,d)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !1D
                            msd(conf,2,0)=msd(conf,2,0)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !3D
                            IF (d/=normaldir) msd(conf,2,4)=msd(conf,2,4)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !2D
                        ELSE IF (j==1 .or. j==Lchains) THEN !End beads
                            msd(conf,3,d)=msd(conf,3,d)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !1D
                            msd(conf,3,0)=msd(conf,3,0)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !3D
                            IF (d/=normaldir) msd(conf,3,4)=msd(conf,3,4)+&
                                        +(Rchain(i,j,d)-R0chain(i,j,d,f2-1)&
                                        -RCoM(d)+R0CoM(d,f2-1))**2 !2D
                        END IF !bead type
                    END DO !j Lchains
                    msd(conf,4,d)=msd(conf,4,d)+&
                                +(Rchain(i,0,d)-R0chain(i,0,d,f2-1)&
                                -RCoM(d)+R0CoM(d,f2-1))**2 !Chain CoMs, 1D
                    msd(conf,4,0)=msd(conf,4,0)+&
                                +(Rchain(i,0,d)-R0chain(i,0,d,f2-1)&
                                -RCoM(d)+R0CoM(d,f2-1))**2 !Chain CoMs, 3D
                    IF (d/=normaldir) msd(conf,4,4)=msd(conf,4,4)+&
                                +(Rchain(i,0,d)-R0chain(i,0,d,f2-1)&
                                -RCoM(d)+R0CoM(d,f2-1))**2 !Chain CoMs, 2D
                END DO !i Nchains
                DO i=1, Npart
                    DO j=1, 2
                        msd(conf,5,d)=msd(conf,5,d)&
                                +(Rpart(i,j,d)-R0part(i,j,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 !All ND beads, 1D
                        msd(conf,5,0)=msd(conf,5,0)&
                                +(Rpart(i,j,d)-R0part(i,j,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 ! 3D
                        IF (d/=normaldir) msd(conf,5,4)=msd(conf,5,4)&
                                +(Rpart(i,j,d)-R0part(i,j,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 !2D
                    END DO !j dimer bead
                    msd(conf,6,d)=msd(conf,6,d)&
                            +(Rpart(i,0,d)-R0part(i,0,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 !NP CoMs, 1D
                    msd(conf,6,0)=msd(conf,6,0)&
                            +(Rpart(i,0,d)-R0part(i,0,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 !3D
                    IF (d/=normaldir) msd(conf,6,4)=msd(conf,6,4)&
                            +(Rpart(i,0,d)-R0part(i,0,d,f2-1)-RCoM(d)+R0CoM(d,f2-1))**2 !2D
                    avdir(conf)=avdir(conf)&
                            +(Rpart(i,2,d)-Rpart(i,1,d))*(R0part(i,2,d,f2-1)-R0part(i,1,d,f2-1))&
                            /(dsqrt((Rpart(i,2,1)-Rpart(i,1,1))**2&
                                   +(Rpart(i,2,2)-Rpart(i,1,2))**2&
                                   +(Rpart(i,2,3)-Rpart(i,1,3))**2)&    
                             *dsqrt((R0part(i,2,1,f2-1)-R0part(i,1,1,f2-1))**2&
                                   +(R0part(i,2,2,f2-1)-R0part(i,1,2,f2-1))**2&
                                   +(R0part(i,2,3,f2-1)-R0part(i,1,3,f2-1))**2))    
                END DO !i Npart
            END DO !dims
        END IF !conf
    END DO !f2 files
    !print*, '    --> accumulation finished'
END DO !f

!print*, 'Averaging'
! DIVIDE BY NUMBER OF TIMES COMPUTED
DO i=1,6
    IF (i==1) d2=Nchains*Lchains
    IF (i==2 .or. i==3) d2=Nchains*2
    IF (i==4) d2=Nchains
    IF (i==5) d2=Npart*2
    IF (i==6) d2=Npart
    DO conf=1,Nconfigs
        IF (Npart>0 .and. i==6) avdir(conf)=avdir(conf)/dfloat(d2*avconf(conf))
        DO d=0,4
            msd(conf,i,d)=msd(conf,i,d)/dfloat(d2*avconf(conf))
        END DO !dims
    END DO !conf Nconfigs
END DO !i types
!print*, '    --> average done'
!print*, 'Integrating autocorrelation'
! Trapezoid rule
IF (Npart>0) THEN
    integral=0.d0
    DO i=1,Nconfigs
        DO j=1,i
            integral(i)=integral(i)+0.5d0*(avdir(j)+avdir(j-1))*delta_t*dfloat(t(j)-t(j-1))
        END DO !j 1,i
    END DO !i Nconfigs
!print*, '    --> integration done'
!print*, 'Printing outputs'
!OUTPUT MSD, ROT
open(501,file='msd_allmonomers.dat')
open(502,file='msd_centralmonomers.dat')
open(503,file='msd_endingmonomers.dat')
open(504,file='msd_chainsCM.dat')
IF (Npart>0) THEN
    open(505,file='msd_NPbeads.dat')
    open(506,file='msd_NPCM.dat')
END IF
tipo=4
IF (Npart>0) tipo=6
DO i=1,tipo
    WRITE(500+i,*) '# time(tau) msd_x msd_y msd_z msd_plane msd_3d'    
    DO conf=1, Nconfigs
        WRITE(500+i,500) dfloat(t(conf))*delta_t, msd(conf,i,1), msd(conf,i,2),&
            msd(conf,i,3), msd(conf,i,4), msd(conf,i,0)        
    END DO !Nconfigs
    close(500+i)
END DO !Ntypes
500 FORMAT (F14.6,2X,F20.12,2X,F20.12,2X,F20.12,2X,F20.12,2X,F20.12)
open(601,file='avdir_autocorr_Drot.dat')
WRITE(601,*) '# time(tau) avdir ln_avdir integral_avdir(tau) D_rot(tau^-1)'   
DO i=1,Nconfigs
    WRITE(601,600) dfloat(t(i))*delta_t, avdir(i), dlog(avdir(i)), integral(i), 0.5d0/integral(i)
END DO
600 FORMAT (F14.3,2X,F20.12,2X,F20.12,2X,F20.12,2X,F20.12)
close(601)

END PROGRAM LamellarDynamics
