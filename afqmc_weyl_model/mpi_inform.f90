!--------------------------
!mpi or serial  parammeters
!--------------------------
module mpi_serial_param
implicit none
 integer::ierr
 integer::rank
 integer::Nsize
 integer::htype
 integer::phtype
end module mpi_serial_param


!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
subroutine start_code()
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

#ifdef MPI
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,Nsize,ierr)
#else
 rank=0
 Nsize=1
#endif
end subroutine start_code



!--------------------------------------------
!This subroutine set the mpi htype and phtype
!depend on the dtype: c or d, if dtype=c then
!htype: 2*Nsite*2Nsite, phtype: 2*Nsite*Ntot
!--------------------------------------------
subroutine init_mpi_type(dtype,Nsite,Ntot,Nspin)
use mpi_serial_param
implicit none
character(len=1),intent(IN)::dtype
integer,intent(IN)::Nsite
integer,intent(IN)::Ntot
integer,intent(IN)::Nspin(2)

integer::bl(2*Nsite)
integer::disp(2*Nsite)
integer::i
#ifdef MPI
include "mpif.h"
#endif


#ifdef MPI
 if(dtype.EQ.'c') then
   call MPI_TYPE_CONTIGUOUS(4*Nsite*Nsite,MPI_DOUBLE_COMPLEX,htype,ierr)
   call MPI_TYPE_CONTIGUOUS(2*Nsite*Ntot,MPI_DOUBLE_COMPLEX,phtype,ierr)
 else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then

   bl=Nsite
   do i=1,Nsite,1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nsite,1
      disp(i+Nsite)=(Nsite+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(2*Nsite,bl,disp,MPI_DOUBLE_COMPLEX,htype,ierr)


   bl=Nsite
   do i=1,Nspin(1),1
      disp(i)=(i-1)*2*Nsite
   end do
   do i=1,Nspin(2),1
      disp(i+Nspin(1))=(Nspin(1)+i-1)*2*Nsite+Nsite
   end do
   call MPI_TYPE_INDEXED(Ntot,bl,disp,MPI_DOUBLE_COMPLEX,phtype,ierr)

 end if
 call MPI_TYPE_COMMIT(htype,ierr)
 call MPI_TYPE_COMMIT(phtype,ierr)
#else
 return
#endif

end subroutine init_mpi_type



!End the all the mpi part
subroutine end_mpi()
use mpi_serial_param 
implicit none
#ifdef MPI
include "mpif.h"
#endif

#ifdef MPI
 call MPI_TYPE_FREE(htype,ierr)
 call MPI_TYPE_FREE(phtype,ierr)
 call MPI_Finalize(ierr)
#else
 return
#endif
end subroutine end_mpi
