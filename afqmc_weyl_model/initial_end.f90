module all_param
use param
use model_param
use project_param
use mc_loop_param
use phi_x_param
use meas_param
use lattice_param
use rand_num
use timing_module
use io_module
use mpi_serial_param
use fftw_param
implicit none
end module all_param


!----------------------------------------------------
!This subroutine is used to read parameter from param
!at the beginning of running the code++++++++++++++++
!----------------------------------------------------
subroutine readparam
use all_param
implicit none
open(unit=10,file='param',status='old')

!read the lattice parameter
read(10,*) set
read(10,*) Nsite
read(10,*) Nhop
read(10,*) Nl(1)
read(10,*) Nl(2)
read(10,*) Nl(3)
read(10,*) kbound(1)
read(10,*) kbound(2)
read(10,*) kbound(3)
read(10,*) openbcx
read(10,*) openbcy
read(10,*) t1
read(10,*) vhop
read(10,*) whop
read(10,*) tyhop
read(10,*) tdhop
read(10,*) lamda

!read the model parameter
read(10,*) onsitU
read(10,*) Ntot
read(10,*) dtype
read(10,*) Nspin(1)
read(10,*) Nspin(2)

if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
  Ntot=Nspin(1)+Nspin(2)
end if


!Project parameter
read(10,*) dt
read(10,*) kcrn
read(10,*) bgset
read(10,*) pfft
read(10,*) diagm
!if we use diagm for measure, we must set pfft.eq.1
if(openbcx.ne.1.and.openbcy.ne.1)then
   if(pfft.NE.1.or.diagm.NE.1) then
      write(*,*) "pfft and diagm must be 1"
      stop
   end if
endif

!MC loop parameter
read(10,*) Ntherm
read(10,*) Nmeas
read(10,*) StepforGram


!Phi parameter
read(10,*) PP
read(10,*) Nlen
if(mod(Nlen,2).NE.0) then
  write(*,*) "Nlen must be even:",Nlen
  call mystop  
end if

read(10,*) blk
if(mod(Nlen,blk).NE.0) then
  write(*,*) "Nlen must be devided by blk:",Nlen,blk
  call mystop
end if
Nblk=Nlen/blk
if(Nblk.NE.blk) then
  if(rank.eq.0) write(*,*) Nblk,blk,Nlen
  if(rank.eq.0) write(*,*) "suggestion: blk*blk=Nlen will save the memory!"
end if
if(mod(blk,StepforGram).NE.0) then
  if(rank.eq.0) write(*,*) "warning: blk need to be devided by StepforGram!:",blk,StepforGram
  call mystop
end if

read(10,*) meastep
if(mod(Nlen/2,meastep).NE.0) then
  if(rank.eq.0) write(*,*) "We want Nlen/2 to be devided by meastep to reach middle point:"
  if(rank.eq.0) write(*,*) Nlen/2,meastep
  call mystop
end if

read(10,*) thermstep
if(thermstep.GT.Nlen/2) then
  if(rank.eq.0) write(*,*) "Thermstep is too big, only energy will be produced!"
end if
close(10)
end subroutine readparam



!-----------------------------------------------------------------
!This subroutine deallocate all the arrays at the end of the code.
!-----------------------------------------------------------------
subroutine deallocatearray()
implicit none
call deallocate_lattice_array()
call deallocate_project_array()
call deallocate_qmc() 
call deallocate_one_meas()
end subroutine deallocatearray



!-------------------------------------------------
!This subroutine clean call the things before stop
!-------------------------------------------------
subroutine clean
use all_param
implicit none
if(rank.eq.0) write(*,*) "Clean for the stop."
call end_fftw()
call end_genrand()
call deallocatearray()
call EndTiming()
call end_mpi()
end subroutine clean



!----------------------------------------------------------------
!This subroutine end the program before deallocate all the arrays
!----------------------------------------------------------------
subroutine mystop
implicit none
call clean()
stop
end subroutine mystop
