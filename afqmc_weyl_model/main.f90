!main program of code
program main
use all_param
implicit none

!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
 call start_code()

!-----------------------------------------------
!BeginTiming and EndTiming get the running time.
!-----------------------------------------------
 if(rank.eq.0) write(*,*) "Start to run the routine."
 call BeginTiming()
 call PrintTiming()


!----------------------------------------------
!initial parameter:random number and read param
!----------------------------------------------
 call init_genrand()
 call readparam()


!-----------------------------------------------------------
!Get the lattice information,Get Hzero matrix of the lattice
!-----------------------------------------------------------
 call set_lattice()

!--------------------------------------------------
!Set the mpi htype and phtype after set the lattice
!--------------------------------------------------
 call init_mpi_type(dtype,Nsite,Ntot,Nspin)


!------------------------
!initial the plan in fftw
!------------------------
 call init_fftw()


!-------------------------------------
!Inital the projection part of K and V
!-------------------------------------
 call inital_k_v()


!-------------------------------------
!Initial the QMC,phil and phir and aux
!------------------------------------- 
 call set_qmc()

!-----------------------
!run the QMC and meausre
!-----------------------
 call run_qmc()


!-------------------
!manipulate the data
!-------------------
call data_mani()


!-------------------------------
!write the wave function and aux
!-------------------------------
!call write_phi_aux_serial() 

!--------------------
!Prepare for the stop
!--------------------
 call clean
end program main


!use serial code to generate the phi and aux for mpi code
subroutine write_phi_aux_serial()
use all_param
implicit none
integer::i

if(Nsize.NE.1) then
  write(*,*) "Something is wrong, this routine only used for serial version:",Nsize
  call mystop
end if
call write_phi(DNsite,Ntot,phi_Nblk(1,1,0),phi_Nblk(1,1,Nblk+2))

do i=0,95,1
   write(*,*) "writing aux:",i
   call one_step()
   call write_aux_rank(Nsite,Nlen,aux,i)
end do
end subroutine write_phi_aux_serial


!write the phi and aux for re-run the mpi code
subroutine write_phi_aux_mpi()
use all_param
implicit none
if(rank.eq.0) then
  call write_phi(DNsite,Ntot,phi_Nblk(1,1,0),phi_Nblk(1,1,Nblk+2))
end if
call write_aux_rank(Nsite,Nlen,aux,rank)
end subroutine write_phi_aux_mpi
