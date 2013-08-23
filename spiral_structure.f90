!> @file
!! Describing the data structure, Physics related parameter and other machine
!! related parameters here.
!! @author Chien-Chang Feng

!> Define constans
MODULE NUM
USE,INTRINSIC                           ::iso_fortran_env
INTEGER,PARAMETER                       ::qp            = REAL128       
!< Presicion of Real Numbers
REAL(REAL128),PARAMETER                 ::zerolimit     = 1.d-6         
!< Numbers smaller than this number  will be considered zero.
REAL(REAL128),PARAMETER                 ::GravConst     = 4.3d-6 
!< Gravitational Constant with the unit in kpc
REAL(REAL128),PARAMETER                 ::g             = 4.3d0
!< Gravitational Constant with the unit in pc
REAL(REAL128),PARAMETER                 ::pi            = 4.d0*atan(1.d0)
!< Pi
ENDMODULE NUM

!> Modul: STELLARDISK_MODEL
!> Define data structure only. Related procedures other than io attached to them.
MODULE STELLARDISK_MODEL
USE NUM
!> Define type typgagaxy_para for store parameters describing rotation curve.
TYPE   typgalaxy_para
        !> Parameters relate to rotation curve.
        REAL(REAL128)                    ::para(14)
        !> Flag to show if para been initialized.
        LOGICAL                          ::inited=.false.
        CONTAINS 
        !> Print para on the screen.
        PROCEDURE                        ::print
        !> Read para from the file defined by #parafnm.
        PROCEDURE                        ::readstd
ENDTYPE
!> Define type typspiralmodel which contains all 1-D calculation results,
!!  pattern speed and special locations like CO and Four to One
!! and other program controlling flags.
TYPE,EXTENDS(typgalaxy_para)            ::typspiralmodel
        !> Solved eigen value.
        COMPLEX(REAL128)                 ::w
        !> Mode number
        INTEGER                          ::mode
        !> Resolution of RK4 and the length of u, h1, k3 and r.
        INTEGER                          ::N
        !> @brief Solved u by RK4
        !! - u(1,:) is data point of r.
        !! - u(2,:) is u.
        !! - u(3,:) is du.
        COMPLEX(REAL128),ALLOCATABLE     ::u(:,:) 
        !> First perturbation of enthalpy.
        !! Coordinate is in u(1,:) or r(:)
        COMPLEX(REAL128),ALLOCATABLE     ::h1(:)
        !> Data points of k3.
        !! I don't even remebmer I have set this.
        !! Coordinate is in u(1,:) or r(:)
        COMPLEX(REAL128),ALLOCATABLE     ::k3(:)
        !> The same as u(1,:) but in REAL.
        REAL(REAL128),ALLOCATABLE        ::r(:)
        !> End of RK4 calculation. This should larger than four to one.
        REAL(REAL128)                    ::rmax
        !> Begining of RK4 calculation. This may deviated from zero to avoid
        !! singularity
        REAL(REAL128)                    ::rmin
        !> Position of corotation
        REAL(REAL128)                    ::co              
        !> position of outer four to one
        REAL(REAL128)                    ::fortoone        
        !> Exact error value at the end of calculation.
        COMPLEX(REAL128)                      ::error
        !> Starting angle for 2D plot
        !! @todo check if phase here is using.
        REAL(REAL128)                    ::phase     = 0.d0
        !> Separation between two points.
        REAL(REAL128)                    ::dr              
        !> Pattern speed.
        REAL(REAL128)                    ::pspd
        !> Controll flag. Show if u is calculated.
        LOGICAL                          ::ucaled    = .false.
        !> Controll flag. Show if co is calculated.
        LOGICAL                          ::cocaled   = .false.
        !> Controll flag. Show if h1 is calculated.
        LOGICAL                          ::h1caled   = .false.
        !> Controll flag. Show which boundary condition at r=0.
        LOGICAL                          ::bndu0     = .true.
        !> Controll flag. Show if eigen value is calculated.
        LOGICAL                          ::winit     = .false.
        CONTAINS
        !> Initalizing typspiral.
        PROCEDURE,PASS                   ::init
        !> Free typspiral.
        PROCEDURE,PASS                   ::free
        !> Read eigen value from file.
        !! @todo try to merge this with setw.
        PROCEDURE,PASS                   ::readw
        !> Set the value of eigen value by hand.
        PROCEDURE,PASS                   ::setw
ENDTYPE

!> Standard stellar variables.
TYPE(typgalaxy_para),TARGET,SAVE        ::stdpara
!> File name where para be load. Default is 'para.list'
CHARACTER(*),PARAMETER                  ::parafnm = 'para.list'
CONTAINS

!> Output parameters on the stdio
SUBROUTINE print(this)
IMPLICIT NONE
class(typgalaxy_para),intent(in)           ::this
!Halo
REAL(REAL128)                ::Lh,rhoh,gHalo,VHalo
!bulge
REAL(REAL128)                ::rb,Mb,gBulge,VBulge
!disk
REAL(REAL128)                ::dM,da,db,VDisk
!Toomre Q
REAL(REAL128)                ::Q,Qod,rq
!Lau Disk
REAL(REAL128)                ::a1,a2,M1,M2
!pspd from readin
REAL(REAL128)                ::w(4)
!NAME LIST
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2
if(.not.this.inited)then
        write(0,*)'para not init,print failed'
else
      Lh        = this.para(1)
      rhoh      = this.para(2)
      Mb        = this.para(3)
      rb        = this.para(4)
      dM        = this.para(5)
      da        = this.para(6)
      db        = this.para(7)
      Qod       = this.para(8)
      q         = this.para(9)
      rq        = this.para(10)
      a1        = this.para(11)
      a2        = this.para(12)
      M1        = this.para(13)
      M2        = this.para(14)
      write(6,*)'para print:'
      write(6,nml=paralist)
endif
ENDSUBROUTINE

!> Read parameters from #parafnm
SUBROUTINE readstd(this)
IMPLICIT NONE
CLASS(typgalaxy_para)           ::this
!Halo
REAL(REAL128)                ::Lh,rhoh,gHalo,VHalo
!bulge
REAL(REAL128)                ::rb,Mb,gBulge,VBulge
!disk
REAL(REAL128)                ::dM,da,db,VDisk
!Toomre Q
REAL(REAL128)                ::Q,Qod,rq
!Lau Disk
REAL(REAL128)                ::a1,a2,M1,M2
!NAME LIST
namelist /paralist/ Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2

!$OMP CRITICAL
open(10,file='para.list')
read(10,nml=paralist)
close(10)
!$OMP END CRITICAL

stdpara.para = (/Lh,rhoh,Mb,rb,dM,da,db,Qod,q,rq,a1,a2,M1,M2/)
stdpara.inited = .true.
ENDSUBROUTINE

!> Read eigen value from file.
!! @param[in] mode Mode Number
!! @todo remove?
SUBROUTINE readw(this,mode)
IMPLICIT NONE
class(typspiralmodel)                   ::this
!pspd from readin
REAL(REAL128)                        ::w(4)
INTEGER                                 ::mode
LOGICAL                                 ::bndu0(2)
namelist /spiralnml/ w,bndu0

this.mode = mode
!$OMP CRITICAL
open(10,file=parafnm)
read(10,nml=spiralnml)
close(10)
!$OMP END CRITICAL
!this.bndu0 = bndu0(mode)

this.w = dcmplx(w(mode*2-1),w(mode*2))
this.winit = .true.
this.pspd = real(this.w)/2.d0
ENDSUBROUTINE
 
!> Set the value of eigen value by hand.
!! @param[in] w DOUBLE COMPLEX, value of eigenvalue.
!! @param[in] mode INTEGER, mode number.
SUBROUTINE setw(this,w,mode)
IMPLICIT NONE
class(typspiralmodel)                   ::this
COMPLEX(REAL128)                          ::w
INTEGER                                 ::mode

this.mode = mode
this.w = w 
this.winit = .true.
this.pspd = real(this.w)/2.d0
ENDSUBROUTINE

!> Initalizing typspiral.
!! @param[in] n Resolution of solving u.
!! @param[in] domain Right boundary of gas simulation. 
!!  The outter bondary of RK4 will be domain*1.5.
!! @param[in] para Parameters template for typspiral.
!! @param[in] mode Mode Number.
SUBROUTINE init(this,n,domain,para,mode)
IMPLICIT NONE
class(typspiralmodel)                   ::this
type(typgalaxy_para)                    ::para
INTEGER                                 ::n,mode
REAL(REAL128)                        ::domain
IF(.NOT.ALLOCATED(this.u))ALLOCATE(this.u(3,2*n))
IF(.NOT.ALLOCATED(this.h1))ALLOCATE(this.h1(2*n))
IF(.NOT.ALLOCATED(this.r))ALLOCATE(this.r(2*n))
IF(.NOT.ALLOCATED(this.k3))ALLOCATE(this.k3(2*n))
this.rmin = 0.d0 
this.rmax = 1.5d0*domain
this.n    = 2*n
this.para = para.para
this.inited = .true.
this.dr   = (this.rmax - this.rmin)/dble(this.n)
ENDSUBROUTINE

!> Free typspiral.
SUBROUTINE free(this)
IMPLICIT NONE
class(typspiralmodel)                   ::this
IF(ALLOCATED(this.u))DEALLOCATE(this.u)
IF(ALLOCATED(this.h1))DEALLOCATE(this.h1)
IF(ALLOCATED(this.r))DEALLOCATE(this.r)
IF(ALLOCATED(this.k3))DEALLOCATE(this.k3)
ENDSUBROUTINE

ENDMODULE

