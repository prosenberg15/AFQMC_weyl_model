
!This parameter contains the useful parameters
!---------------------------------------------
module param
implicit none
 complex(kind=8),parameter:: Xi=dcmplx(0d0,1d0)
 complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
 complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)
 real(kind=8)::Pi=2.d0*asin(1.d0)
 !real(kind=8)::Pi=3.1415926535898d0
end module param



!-------------------------------------------
!This parameter contains all the model_param
!-------------------------------------------
module model_param
implicit none
 real(kind=8)::onsitU      !Hubbard U interaction on the same site.
 integer:: Ntot            !the tot number of electrons
 character(len=1):: dtype  !for the determinate type: d decouple, c couple.
 integer:: Nspin(2)        !Nup and Ndn
end module model_param



!-------------------------------------------------
!This parameter contains all the QMC project param
!-------------------------------------------------
module project_param
implicit none
 real(kind=8):: dt        ! each slice of imagine time
 integer::kcrn      !Decide the different kind of free projection
 complex(kind=8),allocatable::exp_halfK(:,:)   ! DNsite,DNsite.
 complex(kind=8),allocatable::exp_mhalfK(:,:)   ! DNsite,DNsite.
 complex(kind=8),allocatable::exp_K(:,:)   ! DNsite,DNsite.
 complex(kind=8),allocatable::exp_mK(:,:)   ! DNsite,DNsite.
 real(kind=8),allocatable::exp_he(:),exp_mhe(:),exp_e(:),exp_me(:) !DNsite
 complex(kind=8),allocatable::uk(:),vk(:) !Nsite,for the unitary transformation
 complex(kind=8)::gamaf     !the gama for the free-projection
 complex(kind=8),allocatable::ng(:) !Nsite, the back ground for free projection
 integer::bgset !0. mean field 1. dynamic background walker 
 integer::pfft  !0. do not use fftw  1. use fftw for the code
 integer::diagm !0. get all the matrix 
                !1. only measure diaganal matrix. pfft must be 1.
end module project_param


!----------------------------------------------
!This parameter contains all the QMC loop param
!----------------------------------------------
module mc_loop_param
implicit none
 integer::Ntherm          ! number of therm
 integer::Nmeas           ! number of measurement
 integer::StepforGram     ! number of steps when we do modified GS
 integer::meastep         ! measure after each meastep
 integer::thermstep       ! thermstep therm for non-commute observables.
end module mc_loop_param



!----------------------------------------------------------
!In FPMC and RCPMC: |psi>= sum_i rx(i) weight(i) phi(:,:,i)
!----------------------------------------------------------
module phi_x_param
implicit none
 integer::PP   !0 phi and aux from code, 1 aux from code read phi, 2 read aux and phi
 integer::Nlen !beta length of the projection
 integer::Nblk,blk !Nblk: number of block, blk: size of each block
 complex(kind=8),allocatable::phi(:,:,:)       !DNsite,Ntot,2 determinate.
 real(kind=8),allocatable::aux(:,:) !Nsite, Nlen
 complex(kind=8),allocatable::phi_Nblk(:,:,:) !DNsite,Ntot,Nblk+3
 complex(kind=8),allocatable::phi_blk(:,:,:)  !DNsite,Ntot,blk
 real(kind=8),allocatable::rphi(:) !blk: the rescaling after MGS
 complex(kind=8)::imp_save !save the imp function
end module phi_x_param


!----------------------------------------
!The meas parameter in measure subroutine
!----------------------------------------
module meas_param
implicit none
integer(kind=8)::i_energy
integer(kind=8)::i_observ
complex(kind=8)::K_one,K_c
complex(kind=8)::V_one,V_c
complex(kind=8)::E_one,E_c
complex(kind=8)::Nup_one,Nup_c,Ndn_one,Ndn_c
complex(kind=8)::nconds_one,ncondt_one
complex(kind=8)::nconds_c,ncondt_c
complex(kind=8),allocatable::didj_one(:),didj_c(:),dk_one(:)
complex(kind=8),allocatable::didj_true_site_one(:),didj_true_site_c(:)
complex(kind=8),allocatable::sisj_one(:),sisj_c(:),sk_one(:)
complex(kind=8),allocatable::ninj_one(:),ninj_c(:),nnk_one(:)
complex(kind=8),allocatable::ninj_true_site_one(:),ninj_true_site_c(:)
complex(kind=8),allocatable::d_one(:),d_c(:)
complex(kind=8),allocatable::edgecup_one(:,:),edgecup_c(:,:)
complex(kind=8),allocatable::edgecdn_one(:,:),edgecdn_c(:,:)
complex(kind=8),allocatable::niup_one(:),niup_c(:)
complex(kind=8),allocatable::nidn_one(:),nidn_c(:)
complex(kind=8),allocatable::skx_one(:),skx_c(:)
complex(kind=8),allocatable::sky_one(:),sky_c(:)
complex(kind=8),allocatable::ck_one(:),ck_c(:)
complex(kind=8),allocatable::pair_full(:,:),onebody(:,:)
character(len=300)::basename,EnergyName,numName,scorrName,skName,cksName,didjName,dkName
character(len=300)::ninjName,nnkName,pairmName,edgecName,niName,bondName
integer(kind=8)::Nupdate,Naccept
end module meas_param

