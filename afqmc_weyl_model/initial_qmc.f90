subroutine set_qmc()
use all_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer i,j
complex(kind=8)::phtmp(DNsite,Ntot)

call allocate_qmc()

if(PP.eq.0) then

   if((openbcx.eq.1).or.(openbcy.eq.1)) then
      if(rank.eq.0) then
         write(*,*) "Generate phi and aux from code."
         if(dtype.EQ.'c') then
            call input_phiT(DNsite,Ntot,Hzero,phtmp)
         else if(dtype.EQ.'d') then
            call input_phiT_d(phtmp)
         else if(dtype.eq.'w') then
            call input_phiT_d(phtmp)
         end if
      end if
   else
      if(rank.eq.0) then
         call sort_phiT(phtmp)
      end if
   end if
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phtmp(1,1),1,phtype,0,MPI_COMM_WORLD,IERR)
#endif


!CHANGE TO MAKE SURE UP=CONJ(DN)
!uncomment if lambda=0 (i.e. no SOC)
do i=1, Nsite, 1
   do j=1, Nspin(1),1
      phtmp(i,j) = conjg(phtmp(i+Nsite,j+Nspin(1)))
   end do
end do

  phi(:,:,1)=phtmp(:,:)
  phi(:,:,2)=phtmp(:,:)
  !do i = 1, DNsite, 1
  !   do j= 1, DNsite, 1
  !      phi(i,j,2)=conjg(phtmp(j,i))
  !   enddo
  !enddo

!if (rank.eq.0) then
!   do i=1, DNsite, 1
!      do j=1, Ntot, 1
!         write(*,'(2E26.16)',advance="no") dble(phtmp(i,j)),dimag(phtmp(i,j))
!      end do
!      write(*,*)
!   end do
!end if

else if(PP.eq.1) then
  if(rank.eq.0) then
    write(*,*) "Generate aux from code and read phi."
    call read_phi(DNsite,Ntot,phi)
  end if
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phi(1,1,1),2,phtype,0,MPI_COMM_WORLD,IERR)
#endif

else if(PP.eq.2) then

  if(rank.eq.0) then
    write(*,*) "Read aux and phi."
    call read_phi(DNsite,Ntot,phi)
  end if
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call MPI_BCAST(phi(1,1,1),2,phtype,0,MPI_COMM_WORLD,IERR)
#endif

  call read_aux_rank(Nsite,Nlen,aux,rank)

else
   write(*,*) "PP must be 0,1,2:",PP
   call mystop
end if

!initial aux and the path
call init_path()

end subroutine set_qmc


!Get the phiT by sort the eigenvalue and fftw
subroutine sort_phiT(phtmp)
use all_param
implicit none
complex(kind=8),intent(OUT)::phtmp(DNsite,Ntot)
complex(kind=8)::phtmpk(DNsite,Ntot)
real(kind=8)::env(DNsite)
integer::posit(DNsite)
integer::mi,mj,singlezerocount,doublezerocount
integer::singlecount,doublecount,count,countminus,countplus

if(dtype.ne.'w') then
   do mi=1,DNsite,1
      env(mi)=-dlog(exp_e(mi))/dt
   end do

   ! uncomment below for half-filling

   !do mi=1,Nsite,1
   !  if(abs(env(mi))<1e-12) exit
   !end do

   !env(mi)=env(mi)-0.000001
   !env(mi+Nsite)=env(mi+Nsite)-0.000001

   ! break degeneracy
   singlezerocount=0
   doublezerocount=0
   do mi=1,Nsite,1
      if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))>1e-12))then
         singlezerocount=singlezerocount+1
      end if
      if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))<1e-12))then
         doublezerocount=doublezerocount+1
      end if
   end do

   if (rank.eq.0) then
      write(*,*) singlezerocount
      write(*,*) doublezerocount
   endif

   singlecount=1
   doublecount=1
   !do while (count.le.(zerocount/2))
   do mi=1,Nsite,1
      if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))>1e-12).and.(singlecount.le.singlezerocount))then
         env(mi)=env(mi)-0.000001
         singlecount=singlecount+1
      end if
      if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))<1e-12).and.(doublecount.le.(doublezerocount/2)))then
         env(mi)=env(mi)-0.000001
         env(mi+Nsite)=env(mi+Nsite)-0.000001
         doublecount=doublecount+1
      end if
   end do
   !end do

!countminus=1
!countplus=1
!do mi=1,Nsite,1
!  if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))<1e-12).and.(countplus.eq.1)) then
!     if(rank.eq.0) then
!        write(*,*) mi, env(mi), env(mi+Nsite)
!     end if
!     exit
!     env(mi)=env(mi)-0.000001
!     env(mi+Nsite)=env(mi+Nsite)-0.000001
!     countplus=countplus+1
!  end if
!  if((abs(env(mi))<1e-12).and.(abs(env(mi+Nsite))>1e-12))then!.and.(countminus<3)) then
!     env(mi)=env(mi)-0.000001
!     countminus=countminus+1
!  end if
  !if((abs(env(mi+Nsite))<1e-12).and.(abs(env(mi))>1e-12).and.(countplus<3)) then
  !   env(mi+Nsite)=env(mi+Nsite)-0.000001
  !   countplus=countplus+1
  !end if
!     count=count+1
!  end if
!end do

!env(mi)=env(mi)-0.000001
!env(mi+Nsite)=env(mi+Nsite)-0.000001

!if (rank.eq.0) then
!   do mi=1,DNsite,1
!      write(*,*) mi, env(mi), env(mi+Nsite)
!   end do
!end if
!stop
elseif(dtype.eq.'w') then
   do mi=1,DNsite,1
      env(mi)=-dlog(exp_e(mi))/dt
   end do

   ! break degeneracy
   !singlezerocount=0
   !doublezerocount=0
   !do mi=1,Nbravais,1
   !   if((abs(env(mi))<1e-12).and.(abs(env(mi+Nbravais))>1e-12))then
   !      singlezerocount=singlezerocount+1
   !   end if
   !   if((abs(env(mi))<1e-12).and.(abs(env(mi+Nbravais))<1e-12))then
   !      doublezerocount=doublezerocount+1
   !   end if
   !end do

   !if (rank.eq.0) then
   !   write(*,*) singlezerocount
   !   write(*,*) doublezerocount
   !endif

   !singlecount=1
   !doublecount=1
   !do while (count.le.(zerocount/2))
   !do mi=1,Nsite,1
   !   if((abs(env(mi))<1e-12).and.(abs(env(mi+Nbravais))>1e-12).and.(singlecount.le.singlezerocount))then
   !      env(mi)=env(mi)-0.000001
   !      singlecount=singlecount+1
   !   end if
   !   if((abs(env(mi))<1e-12).and.(abs(env(mi+Nbravais))<1e-12).and.(doublecount.le.(doublezerocount/2)))then
   !      env(mi)=env(mi)-0.000001
   !      env(mi+Nbravais)=env(mi+Nbravais)-0.000001
   !      doublecount=doublecount+1
   !   end if
   !end do

endif
   
if(dtype.EQ.'c') then
  call sort2(DNsite,env,posit)
  phtmpk=zero
!  if (rank.eq.0) then
!     do mi=1,DNsite,1
!        write(*,*) mi, env(mi)
!     end do
!  end if
  do mi=1,Ntot,1
     mj=posit(mi)
     phtmpk(mj,mi)=one
  end do
else if(dtype.EQ.'d') then
  call sort2(Nsite,env,posit)
  phtmpk=zero
  do mi=1,Nspin(1),1
     mj=posit(mi)
     phtmpk(mj,mi)=one
  end do
  do mi=1,Nspin(2),1
     mj=posit(mi)
     phtmpk(mj+Nsite,mi+Nspin(1))=one
  end do
else if(dtype.eq.'w') then
   call sort2(Nsite,env,posit)
   !call sort2(Nsite,env(1+Nsite:DNsite),posit(1+Nsite:DNsite))
   if (rank.eq.0) then
      do mi=1,Nsite,1
         write(*,*) mi, env(mi)!, env(mi+Nsite)
      end do
   end if
   phtmpk=zero
   do mi=1,Nspin(1),1
      mj=posit(mi)
      phtmpk(mj,mi)=one
   end do
   do mi=1,Nspin(2),1
      mj=posit(mi)
      phtmpk(mj+Nsite,mi+Nspin(1))=one
   end do
end if
call ph_dc_fftw_b(phtmpk(1:DNsite,1:Ntot),phtmp(1:DNsite,1:Ntot),uk,vk)
end subroutine sort_phiT



!Diagonialize h0, give the lowest eigenstate to ph
subroutine input_phiT(nl,nt,h0,ph)
implicit none
integer,intent(IN)::nl
integer,intent(IN)::nt
complex(kind=8),intent(IN)::h0(nl,nl)
complex(kind=8),intent(OUT)::ph(nl,nt)
complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)
integer:: i

call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)
call zcopy(nl*nt,hu,1,ph,1)
write(*,*) 'phi input correctly'
do i=1,nl
   write(*,*) i, ev(i)
enddo
end subroutine input_phiT



!Used to call input_phiT in decouple condition
subroutine input_phiT_d(phtmp)
use all_param
implicit none
complex(kind=8),intent(OUT)::phtmp(DNsite,Ntot)
complex(kind=8)::h0(Nsite,Nsite)
complex(kind=8)::ph(Nsite,Ntot)

h0(1:Nsite,1:Nsite)=Hzero(1:Nsite,1:Nsite)
call input_phiT(Nsite,Nspin(1),h0,ph(1,1))
phtmp(1:Nsite,1:Nspin(1))=ph(1:Nsite,1:Nspin(1))

h0(1:Nsite,1:Nsite)=Hzero((Nsite+1):(DNsite),(Nsite+1):(DNsite))
call input_phiT(Nsite,Nspin(2),h0,ph(1,1))
phtmp((Nsite+1):(DNsite),(Nspin(1)+1):Ntot)=ph(1:Nsite,1:Nspin(2))
end subroutine input_phiT_d


!read the phi
subroutine read_phi(DNsite,Ntot,phi)
use io_module
implicit none
integer,intent(IN)::DNsite
integer,intent(IN)::Ntot
complex(kind=8),intent(OUT)::phi(DNsite,Ntot,2)
integer::j,k

open(unit=10,file='phi.dat',status='old')
  do k=1,Ntot,1
     do j=1,DNsite,1
        read(10,*) phi(j,k,1)
     end do
  end do
  do k=1,Ntot,1
     do j=1,DNsite,1
        read(10,*) phi(j,k,2)
     end do
  end do
close(10)
end subroutine read_phi


!write the phi
subroutine write_phi(DNsite,Ntot,phi1,phi2)
use io_module
implicit none
integer,intent(IN)::DNsite
integer,intent(IN)::Ntot
complex(kind=8),intent(IN)::phi1(DNsite,Ntot),phi2(DNsite,Ntot)
integer::j,k

open(unit=10,file='phi.dat',status='replace')
  do k=1,Ntot,1
     do j=1,DNsite,1
        write(10,*) phi1(j,k)
     end do
  end do
  do k=1,Ntot,1
     do j=1,DNsite,1
        write(10,*) phi2(j,k)
     end do
  end do
close(10)
end subroutine write_phi


!read the aux for PP.eq.2
subroutine read_aux_rank(Nsite,Nlen,aux,rank)
use io_module
implicit none
integer,intent(IN)::Nsite
integer,intent(IN)::Nlen
real(kind=8),intent(OUT)::aux(Nsite,Nlen)
integer,intent(IN)::rank
character(len=300)::phi_name
integer::j,k

call createFileName(phi_name,'aux')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call openUnit(phi_name,10,'B')
read(10) aux(1:Nsite,1:Nlen)
close(10)
end subroutine read_aux_rank


!write the aux
subroutine write_aux_rank(Nsite,Nlen,aux,rank)
use io_module
implicit none
integer,intent(IN)::Nsite
integer,intent(IN)::Nlen
real(kind=8),intent(IN)::aux(Nsite,Nlen)
integer,intent(IN)::rank
character(len=300)::phi_name
integer::j,k

call createFileName(phi_name,'aux')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call openUnit(phi_name,10,'C')
write(10) aux(1:Nsite,1:Nlen)
close(10)
end subroutine write_aux_rank


subroutine init_path
use all_param
implicit none
complex(kind=8)::temp(DNsite,Ntot)
complex(kind=8)::explr(DNsite)
real(kind=8)::anm,anm1,anm2
integer::i,i_GS,i_nblk,i_blk


call copy_wf_dc(phi(1,1,1),phi_Nblk(1,1,0))

!apply exp(-dt*K/2) to phi_r 
if(pfft.eq.0) then
   call copy_wf_dc(phi(1,1,1),temp(1,1))
   call k_to_ph_dc(exp_halfK,temp(1,1),phi(1,1,1))
else if(pfft.eq.1) then
   call k_to_ph_dc_fftw(exp_he,phi(1,1,1),uk,vk)
end if



!apply exp(-dt*K)exp(-dt*V) to phi_r(initial aux)
i_GS=0
do i_nblk=1,Nblk,1
   call copy_wf_dc(phi(1,1,1),phi_Nblk(1,1,i_nblk))
   if(i_GS.NE.0) write(*,*) "Warning!! phi_Nblk is not saved before MGS."
   do i_blk=1,blk,1
      i=(i_nblk-1)*blk+i_blk

      !We do not need to get aux if PP.EQ.2 
      if(PP.NE.2) then
        call cal_n_back(phi(1,1,2),phi(1,1,1))
        call sample_field(aux(1,i))
      end if
      call diag_op(aux(1,i),explr)
      call d_to_phi_dc(explr,phi(1,1,1))
 
     !apply exp(-dt*K) to phi_r
      if(pfft.eq.0) then
         call copy_wf_dc(phi(1,1,1),temp(1,1))
         call k_to_ph_dc(exp_K,temp(1,1),phi(1,1,1))
      else if(pfft.eq.1) then
         call k_to_ph_dc_fftw(exp_e,phi(1,1,1),uk,vk)
      end if 


      i_GS=i_GS+1
      if(i_GS.EQ.StepforGram) then
        i_GS=0
        if(dtype.EQ.'c') then
          call modGS(phi(1,1,1),DNsite,Ntot,anm)
        else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
          call modGS(phi(1:Nsite,1:Nspin(1),1),Nsite,Nspin(1),anm1)
          call modGS(phi((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,1),Nsite,Nspin(2),anm2)
        end if
      end if
   end do
end do

call copy_wf_dc(phi(1,1,2),phi_Nblk(1,1,Nblk+2))

!apply Dagger[exp(dt*K/2)] to phi_l
if(pfft.eq.0) then
   call copy_wf_dc(phi(1,1,2),temp(1,1))
   call k_to_ph_dc(exp_mhalfK,temp(1,1),phi(1,1,2))
else if(pfft.eq.1) then
   call k_to_ph_dc_fftw(exp_mhe,phi(1,1,2),uk,vk)
end if

call copy_wf_dc(phi(1,1,2),phi_Nblk(1,1,Nblk+1))

end subroutine init_path


!--------------------------------------
!This subroutine allocate arrays in phi
!--------------------------------------
subroutine allocate_qmc()
use lattice_param
use model_param
use phi_x_param
implicit none
allocate(phi(DNsite,Ntot,2))
allocate(aux(Nsite,Nlen))
allocate(phi_Nblk(DNsite,Ntot,0:Nblk+2))
allocate(phi_blk(DNsite,Ntot,blk))
allocate(rphi(blk))
end subroutine allocate_qmc


!----------------------------------------
!This subroutine deallocate arrays in phi
!----------------------------------------
subroutine deallocate_qmc()
use phi_x_param
implicit none
if(allocated(phi)) deallocate(phi)
if(allocated(aux)) deallocate(aux)
if(allocated(phi_Nblk)) deallocate(phi_Nblk)
if(allocated(phi_blk)) deallocate(phi_blk)
if(allocated(rphi)) deallocate(rphi)
end subroutine deallocate_qmc
