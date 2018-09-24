module fftw_param
implicit none
complex(kind=8),allocatable::aforwin(:),aforwout(:)
complex(kind=8),allocatable::abackin(:),abackout(:)
integer(kind=8)::planf,planb
!complex(kind=8),allocatable::aforwin_w(:),aforwout_w(:)
!complex(kind=8),allocatable::abackin_w(:),abackout_w(:)
!integer(kind=8)::planf_w,planb_w
end module fftw_param


subroutine init_fftw
use fftw_param
use lattice_param
use model_param
implicit none
allocate(aforwin(Nsite),aforwout(Nsite))
allocate(abackin(Nsite),abackout(Nsite))
!choose different plan here
call dfftw_plan_dft_3d(planf,Nl(1),Nl(2),Nl(3),aforwin,aforwout,-1,0)
call dfftw_plan_dft_3d(planb,Nl(1),Nl(2),Nl(3),abackin,abackout,1,0)
!if(dtype.eq.'w') then
!   allocate(aforwin_w(Nbravais),aforwout_w(Nbravais))
!   allocate(abackin_w(Nbravais),abackout_w(Nbravais))
!   !choose different plan here
!   call dfftw_plan_dft_3d(planf_w,Nl(1),Nl(2),Nl(3),aforwin_w,aforwout_w,-1,0)
!   call dfftw_plan_dft_3d(planb_w,Nl(1),Nl(2),Nl(3),abackin_w,abackout_w,1,0)
!endif
end subroutine init_fftw


subroutine end_fftw
use fftw_param
use model_param
implicit none
deallocate(aforwin,aforwout,abackin,abackout)
call dfftw_destroy_plan(planf)
call dfftw_destroy_plan(planb)
!if(dtype.eq.'w') then
!   call dfftw_destroy_plan(planf_w)
!   call dfftw_destroy_plan(planb_w)
!endif
end subroutine end_fftw


subroutine k_to_ph_dc_fftw(eig,phn,uk,vk)
use param
use lattice_param
use model_param
implicit none
real(kind=8),intent(IN)::eig(DNsite)
complex(kind=8),intent(INOUT)::phn(DNsite,Ntot)
complex(kind=8),intent(IN)::uk(Nsite),vk(Nsite)

if(dtype.EQ.'c') then
   call k_to_ph_fftw_c(DNsite,Ntot,eig(1),phn(1:DNsite,1:Ntot),uk(1),vk(1))
else if(dtype.EQ.'d') then
   call k_to_ph_fftw(Nsite,Nspin(1),eig(1),phn(1:Nsite,1:Nspin(1)))
   call k_to_ph_fftw(Nsite,Nspin(2),eig(Nsite+1),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot))
else if(dtype.EQ.'w') then
   call k_to_ph_fftw_c(Nsite,Nspin(1),eig(1),phn(1:Nsite,1:Nspin(1)),uk(1),vk(1))
   call k_to_ph_fftw_c(Nsite,Nspin(2),eig(Nsite+1),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot),uk(1),vk(1))
end if
end subroutine k_to_ph_dc_fftw


subroutine k_to_ph_fftw_c(L,N,exp_e,ph,uk,vk)
use fftw_param
implicit none
integer,intent(IN)::L,N
real(kind=8),intent(IN)::exp_e(L)
complex(kind=8),intent(INOUT)::ph(L,N)
complex(kind=8),intent(IN)::uk(L/2),vk(L/2)
complex(kind=8)::da(L/2),db(L/2),dc(L/2),dd(L/2)
integer::i,j,k

!L must equal the DNsite in the lattice
do i=1,N,1
   !to momentum space
   aforwin(1:L/2)=ph(1:L/2,i)
   call dfftw_execute(planf,aforwin,aforwout)
   da(1:L/2)=aforwout(1:L/2)

   aforwin(1:L/2)=ph((1+L/2):L,i)
   call dfftw_execute(planf,aforwin,aforwout)
   db(1:L/2)=aforwout(1:L/2)

   !apply K matrix 
   ! (must transform to natural orbital (helicity) basis using 
   ! unitary transformation defined by uk, vk)
   do j=1,L/2,1
      dc(j)=(conjg(uk(j))*da(j)+conjg(vk(j))*db(j))*dcmplx(exp_e(j))/dble(L/2)
      dd(j)=(-conjg(uk(j))*da(j)+conjg(vk(j))*db(j))*dcmplx(exp_e(j+L/2))/dble(L/2)
   end do

   !tranform back to spin basis 
   do j=1,L/2,1
      da(j)=uk(j)*(dc(j)-dd(j))
      db(j)=vk(j)*(dc(j)+dd(j))
   end do

   !to real space
   abackin(1:L/2)=da(1:L/2)
   call dfftw_execute(planb,abackin,abackout)
   ph(1:L/2,i)=abackout(1:L/2)

   abackin(1:L/2)=db(1:L/2)
   call dfftw_execute(planb,abackin,abackout)
   ph((1+L/2):L,i)=abackout(1:L/2)

end do
end subroutine k_to_ph_fftw_c


subroutine k_to_ph_fftw(L,N,exp_e,ph)
use fftw_param
implicit none
integer,intent(IN)::L,N
real(kind=8),intent(IN)::exp_e(L)
complex(kind=8),intent(INOUT)::ph(L,N)
integer::i,j,k

!L must equal the Nsite in the lattice
do i=1,N,1
   aforwin(1:L)=ph(1:L,i)
   call dfftw_execute(planf,aforwin,aforwout)
   do j=1,L,1
      abackin(j)=aforwout(j)*dcmplx(exp_e(j))/dble(L)
   end do
   call dfftw_execute(planb,abackin,abackout)
   ph(1:L,i)=abackout(1:L)
end do
end subroutine k_to_ph_fftw



subroutine ph_dc_fftw_f(ph,phn,uk,vk)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::ph(DNsite,Ntot)
complex(kind=8),intent(OUT)::phn(DNsite,Ntot)
complex(kind=8),intent(IN)::uk(Nsite),vk(Nsite)

if(dtype.EQ.'c') then
   call ph_fftw_f_c(DNsite,Ntot,ph(1:DNsite,1:Ntot),phn(1:DNsite,1:Ntot),uk,vk)
else if(dtype.EQ.'d') then
   call ph_fftw_f(Nsite,Nspin(1),ph(1:Nsite,1:Nspin(1)),phn(1:Nsite,1:Nspin(1)))
   call ph_fftw_f(Nsite,Nspin(2),ph((Nsite+1):DNsite,(Nspin(1)+1):Ntot),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot))
else if(dtype.EQ.'w') then
   call ph_fftw_f_c(Nsite,Nspin(1),ph(1:Nsite,1:Nspin(1)),phn(1:Nsite,1:Nspin(1)),uk,vk)
   call ph_fftw_f_c(Nsite,Nspin(2),ph((Nsite+1):DNsite,(Nspin(1)+1):Ntot),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot),uk,vk)
end if
end subroutine ph_dc_fftw_f



subroutine ph_fftw_f_c(L,N,ph,phn,uk,vk)
use fftw_param
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(IN)::ph(L,N)
complex(kind=8),intent(OUT)::phn(L,N)
complex(kind=8),intent(IN)::uk(L/2),vk(L/2)
complex(kind=8)::da(L/2),db(L/2)
integer::i,j,k

!L must equal the Nsite in the lattice
do i=1,N,1
   
   aforwin(1:L/2)=ph(1:L/2,i)
   call dfftw_execute(planf,aforwin,aforwout)
   da(1:L/2)=aforwout(1:L/2)

   aforwin(1:L/2)=ph((1+L/2):L,i)
   call dfftw_execute(planf,aforwin,aforwout)
   db(1:L/2)=aforwout(1:L/2)


   do j=1,L/2,1
      phn(j,i)=(conjg(uk(j))*da(j)+conjg(vk(j))*db(j))/sqrt(L/2.d0)
      phn(j+L/2,i)=(-conjg(uk(j))*da(j)+conjg(vk(j))*db(j))/sqrt(L/2.d0)
   end do 
end do
end subroutine ph_fftw_f_c


subroutine ph_fftw_f(L,N,ph,phn)
use fftw_param
use model_param
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(IN)::ph(L,N)
complex(kind=8),intent(OUT)::phn(L,N)
complex(kind=8)::da(L/2),db(L/2)
integer::i,j,k

if(dtype.ne.'w') then
   !L must equal the Nsite in the lattice
   do i=1,N,1
      aforwin(1:L)=ph(1:L,i)
      call dfftw_execute(planf,aforwin,aforwout)
      do j=1,L,1
         phn(j,i)=aforwout(j)/sqrt(dble(L))
      end do
   end do
else if(dtype.eq.'w') then
   !L must equal the Nsite in the lattice
   do i=1,N,1
      aforwin(1:L/2)=ph(1:L/2,i)
      call dfftw_execute(planf,aforwin,aforwout)
      da(1:L/2)=aforwout(1:L/2)

      aforwin(1:L/2)=ph((1+L/2):L,i)
      call dfftw_execute(planf,aforwin,aforwout)
      db(1:L/2)=aforwout(1:L/2)
      
      do j=1,L/2,1
         phn(j,i)=da(j)/sqrt(dble(L/2.d0))
         phn(j+L/2,i)=db(j)/sqrt(dble(L/2.d0))
      end do
   end do
end if
end subroutine ph_fftw_f



subroutine ph_dc_fftw_b(ph,phn,uk,vk)
use param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(IN)::ph(DNsite,Ntot)
complex(kind=8),intent(OUT)::phn(DNsite,Ntot)
complex(kind=8),intent(IN)::uk(Nsite),vk(Nsite)


if(dtype.EQ.'c') then
   call ph_fftw_b_c(DNsite,Ntot,ph(1:DNsite,1:Ntot),phn(1:DNsite,1:Ntot),uk,vk)
else if(dtype.EQ.'d') then
   call ph_fftw_b(Nsite,Nspin(1),ph(1:Nsite,1:Nspin(1)),phn(1:Nsite,1:Nspin(1)))
   call ph_fftw_b(Nsite,Nspin(2),ph((Nsite+1):DNsite,(Nspin(1)+1):Ntot),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot))
else if(dtype.eq.'w') then
   call ph_fftw_b_c(Nsite,Nspin(1),ph(1:Nsite,1:Nspin(1)),phn(1:Nsite,1:Nspin(1)),uk,vk)
   call ph_fftw_b_c(Nsite,Nspin(2),ph((Nsite+1):DNsite,(Nspin(1)+1):Ntot),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot),uk,vk)
end if
end subroutine ph_dc_fftw_b


subroutine ph_fftw_b_c(L,N,ph,phn,uk,vk)
use fftw_param
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(IN)::ph(L,N)
complex(kind=8),intent(OUT)::phn(L,N)
complex(kind=8),intent(IN)::uk(L/2),vk(L/2)
complex(kind=8)::da(L/2),db(L/2)
integer::i,j,k

!L must equal the Nsite in the lattice
do i=1,N,1

   do j=1,L/2,1
      da(j)=uk(j)*(ph(j,i)-ph(j+L/2,i))/sqrt(L/2.d0)
      db(j)=vk(j)*(ph(j,i)+ph(j+L/2,i))/sqrt(L/2.d0)
   end do

   abackin(1:L/2)=da(1:L/2)
   call dfftw_execute(planb,abackin,abackout)
   phn(1:L/2,i)=abackout(1:L/2) 

   abackin(1:L/2)=db(1:L/2)
   call dfftw_execute(planb,abackin,abackout)
   phn((1+L/2):L,i)=abackout(1:L/2)

end do
end subroutine ph_fftw_b_c


subroutine ph_fftw_b(L,N,ph,phn)
use fftw_param
use model_param
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(IN)::ph(L,N)
complex(kind=8),intent(OUT)::phn(L,N)
complex(kind=8)::da(L/2),db(L/2)
integer::i,j,k

!L must equal the Nsite in the lattice
if(dtype.ne.'w') then
   do i=1,N,1
      abackin(1:L)=ph(1:L,i)
      call dfftw_execute(planb,abackin,abackout)
      do j=1,L,1
         phn(j,i)=abackout(j)/sqrt(dble(L))
      end do
   end do
elseif(dtype.eq.'w') then
   !L must equal the Nsite in the lattice
   do i=1,N,1
      abackin(1:L/2)=ph(1:L/2,i)
      call dfftw_execute(planb,abackin,abackout)
      da(1:L/2)=abackout(1:L/2)

      abackin(1:L/2)=ph((1+L/2):L,i)
      call dfftw_execute(planb,abackin,abackout)
      db(1:L/2)=abackout(1:L/2)
      
      do j=1,L/2,1
         phn(j,i)=da(j)/sqrt(dble(L/2.d0))
         phn(j+L/2,i)=db(j)/sqrt(dble(L/2.d0))
      end do
   end do
end if
end subroutine ph_fftw_b



!!Here we do the fourier transformation to the k space instead of nature orbital space
!subroutine ph_dc_fftw_f_true_k(ph,phn,uk,vk)
!use param
!use lattice_param
!use model_param
!implicit none
!complex(kind=8),intent(IN)::ph(DNsite,Ntot)
!complex(kind=8),intent(OUT)::phn(DNsite,Ntot)
!complex(kind=8),intent(IN)::uk(Nsite),vk(Nsite)
!
!if(dtype.EQ.'c') then
!   call ph_fftw_f_c_true_k(DNsite,Ntot,ph(1:DNsite,1:Ntot),phn(1:DNsite,1:Ntot),uk,vk)
!else if(dtype.EQ.'d') then
!   call ph_fftw_f(Nsite,Nspin(1),ph(1:Nsite,1:Nspin(1)),phn(1:Nsite,1:Nspin(1)))
!   call ph_fftw_f(Nsite,Nspin(2),ph((Nsite+1):DNsite,(Nspin(1)+1):Ntot),phn((Nsite+1):DNsite,(Nspin(1)+1):Ntot))
!end if
!end subroutine ph_dc_fftw_f_true_k
!
!
!subroutine ph_fftw_f_c_true_k(L,N,ph,phn,uk,vk)
!use fftw_param
!implicit none
!integer,intent(IN)::L,N
!complex(kind=8),intent(IN)::ph(L,N)
!complex(kind=8),intent(OUT)::phn(L,N)
!complex(kind=8),intent(IN)::uk(L/2),vk(L/2)
!complex(kind=8)::da(L/2),db(L/2)
!integer::i,j,k
!
!!L must equal the Nsite in the lattice
!do i=1,N,1
!
!   aforwin(1:L/2)=ph(1:L/2,i)
!   call dfftw_execute(planf,aforwin,aforwout)
!   phn(1:L/2,i)=aforwout(1:L/2)/sqrt(L/2.d0)
!
!   aforwin(1:L/2)=ph((1+L/2):L,i)
!   call dfftw_execute(planf,aforwin,aforwout)
!   phn((j+L/2):L,i)=aforwout(1:L/2)/sqrt(L/2.d0)
!
!end do
!end subroutine ph_fftw_f_c_true_k
