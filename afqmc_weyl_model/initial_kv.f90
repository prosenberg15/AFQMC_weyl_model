subroutine inital_k_v()
use all_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

!allocated exp_K,exp_halfK,exp_mhalfK and ng
call allocate_project_array()

if(pfft.eq.0) then
  !since different thread might get different result due to the degeneracy in
  !Hzero, we only do rank.eq.0 and Bcast the result.
  if(dtype.EQ.'c') then
    if(rank.eq.0) call initial_K(DNsite,Ntot,Hzero,dt,exp_K,exp_mK,exp_halfK,exp_mhalfK,exp_he,exp_mhe,exp_e,exp_me,rank)
  else if(dtype.EQ.'d') then
     if(rank.eq.0) call initial_K_d()
  else if(dtype.EQ.'w') then
     if(rank.eq.0) call initial_K_d()
  else
    write(*,*) "Something is wrong with dtype input:",dtype
    call mystop
  end if


#ifdef MPI
   call MPI_BARRIER(MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_K(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_mK(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_halfK(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_mhalfK(1,1),1,htype,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_he(1),DNsite,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_mhe(1),DNsite,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_e(1),DNsite,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(exp_me(1),DNsite,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
#endif

else if(pfft.eq.1) then
   call initial_K_fourier()
else
   write(*,*) "Something is wrong with fft input:",pfft
   call mystop
end if

!get the gamaf and ng(Nsite)
call initial_V() 

end subroutine inital_k_v



!-----------------------------------------------------------
!Input the h0,dt
!Output the exp(-dt*h0/2.d0),exp(-dt*h0) and exp(dt*h0/2.d0)
!-----------------------------------------------------------
subroutine initial_K(nl,np,h0,dt,exp_h0,exp_mh0,exp_half_h0,exp_mhalf_h0,exp_he,exp_mhe,exp_e,exp_me,rank)
implicit none
integer,intent(IN)::nl,np,rank
complex(kind=8),intent(IN)::h0(nl,nl)
real(kind=8),intent(IN)::dt
complex(kind=8),intent(OUT)::exp_h0(nl,nl)
complex(kind=8),intent(OUT)::exp_mh0(nl,nl)
complex(kind=8),intent(OUT)::exp_half_h0(nl,nl)
complex(kind=8),intent(OUT)::exp_mhalf_h0(nl,nl)
real(kind=8),intent(OUT)::exp_he(nl),exp_mhe(nl),exp_e(nl),exp_me(nl)

complex(kind=8)::hu(nl,nl)
real(kind=8)::ev(nl)
integer::i,j,k

complex(kind=8)::htmp(nl,nl)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)
complex(kind=8)::Xi=dcmplx(0.d0,1.d0)

!complex(kind=8)::Amat_loc(nl,nl)
!complex(kind=8)::Sx_loc(nl/2),Sy_loc(nl/2), Sz_loc(nl/2)
!complex(kind=8)::sxsx_loc(nl/2,nl/2), sxsy_loc(nl/2,nl/2), sysy_loc(nl/2,nl/2)
!complex(kind=8)::edgecup(nl/2,nl/2)
!complex(kind=8)::edgecdn(nl/2,nl/2)

call check_Hermite_c(h0,nl)
call zcopy(nl*nl,h0,1,hu,1)
call eigen(Hu,nl,ev)

! calculate Green's functions
! G_ij=[(Psi_0)*(Psi_0)+]_ji

!call cal_Amat(nl,np,Hu,Hu,Amat_loc)
!call cal_Amat_step_dc(Hu,Hu,Amat_loc)

! measure edge current : j_nm = sum_{sigma} -i*(c+_{m,sigma}c_{n,sigma})
!do j=1,nl/2,1
!   do i=1,j,1
!      if (abs(h0(i,j)).gt.1d-12) then         
!         edgecup(i,j)=-Xi*(Amat_loc(j,i)-Amat_loc(i,j))
!         edgecup(j,i)=conjg(edgecup(i,j))
!         edgecdn(i,j)=-Xi*(Amat_loc(j+nl/2,i+nl/2)-Amat_loc(i+nl/2,j+nl/2))
!         edgecdn(j,i)=conjg(edgecdn(i,j))
!      end if
!   end do
!end do

!measure <Sx_i Sx_j>,<Sx_i Sy_j>,<Sy_i Sy_j> 
!do i=1,nl/2,1
!   do j=1,nl/2,1
!      ! Sx Sx
!      sxsx_loc(i,j)=Amat_loc(i,i+nl/2)*Amat_loc(j,j+nl/2) -&
!           &Amat_loc(i,j+nl/2)*Amat_loc(j,i+nl/2) +&
!       &Amat_loc(i,i+nl/2)*Amat_loc(j+nl/2,j) -&
!           &Amat_loc(i,j)*Amat_loc(j+nl/2,i+nl/2) +&
!       &Amat_loc(i+nl/2,i)*Amat_loc(j,j+nl/2) -&
!           &Amat_loc(i+nl/2,j+nl/2)*Amat_loc(j,i) +&
!       &Amat_loc(i+nl/2,i)*Amat_loc(j+nl/2,j) -&
!           &Amat_loc(i+nl/2,j)*Amat_loc(j+nl/2,i)
!      sxsx_loc(i,j)=sxsx_loc(i,j)/4.d0
!      ! Sx Sy
!      sxsy_loc(i,j)=Amat_loc(i,i+nl/2)*Amat_loc(j,j+nl/2) -&
!           &Amat_loc(i,j+nl/2)*Amat_loc(j,i+nl/2) -&
!           Amat_loc(i+nl/2,i)*Amat_loc(j+nl/2,j) +&
!           &Amat_loc(i+nl/2,j)*Amat_loc(j+nl/2,i)
!      sxsy_loc(i,j)=-Xi*sxsy_loc(i,j)/4.d0
!      ! Sy Sy
!      sysy_loc(i,j)=Amat_loc(i,i+nl/2)*Amat_loc(j,j+nl/2) -&
!           &Amat_loc(i,j+nl/2)*Amat_loc(j,i+nl/2) -&
!       &Amat_loc(i,i+nl/2)*Amat_loc(j+nl/2,j) +&
!           &Amat_loc(i,j)*Amat_loc(j+nl/2,i+nl/2) -&
!       &Amat_loc(i+nl/2,i)*Amat_loc(j,j+nl/2) +&
!           &Amat_loc(i+nl/2,j+nl/2)*Amat_loc(j,i) +&
!       &Amat_loc(i+nl/2,i)*Amat_loc(j+nl/2,j) -&
!           &Amat_loc(i+nl/2,j)*Amat_loc(j+nl/2,i)
!      sysy_loc(i,j)=-1.0*sysy_loc(i,j)/4.d0
!   end do
!end do
! measure S_x, S_y
!do i=1,nl/2,1
!   Sx_loc(i)=0.5*(Amat_loc(i,i+nl/2)+Amat_loc(i+nl/2,i))
!   Sy_loc(i)=-0.5*Xi*(Amat_loc(i,i+nl/2)-Amat_loc(i+nl/2,i))
!   Sz_loc(i)=0.5*(Amat_loc(i,i)-Amat_loc(i+nl/2,i+nl/2))
!end do

!if (rank.eq.0) then
!   do i=1,nl/2,1
!      write(*,'(100(A1,F8.5,A1,F8.5,A1))') ('(',real(edgecup(i,j)),',',imag(edgecup(i,j)),')',j=1,nl/2)
!   end do
!end if

!if (rank.eq.0) then
!   do i=1,nl/2,1
!      do j=1,nl/2,1
!         write(*,'(2I4,12F8.5)') i, j, real(Amat_loc(j,i)), imag(Amat_loc(j,i)), real(Amat_loc(i,j)), imag(Amat_loc(i,j)), real(Amat_loc(j+nl/2,i+nl/2)), imag(Amat_loc(j+nl/2,i+nl/2)), real(Amat_loc(i+nl/2,j+nl/2)), imag(Amat_loc(i+nl/2,j+nl/2)), real(edgecup(i,j)), imag(edgecup(i,j)), real(edgecdn(i,j)), imag(edgecdn(i,j))
!      end do
!   end do
!end if
!
!if (rank.eq.0) then
!   do i=1,nl/2,1
!!      do j=1, nl/2,1
!!         write(*,'(2I4,12F8.5)') i, j, real(Amat_loc(i,j)), imag(Amat_loc(i,j))
!      write(*,*) i, Sx_loc(i), Sy_loc(i), Sz_loc(i)
!!      end do
!   end do
!end if
!if (rank.eq.0) then
!   do i=1,nl/2,1
!      do j=1, nl/2,1
!!         write(*,'(2I4,12F8.5)') i, j, real(Amat_loc(i,j)), imag(Amat_loc(i,j))
!      write(*,*) i, j, nl*nl*sxsx_loc(i,j)/(np*(np-1)), nl*nl*sxsy_loc(i,j)/(np*(np-1)), nl*nl*sysy_loc(i,j)/(np*(np-1))
!      end do
!   end do
!end if


!if (rank.eq.0) then
!   do i=1,nl/2,1
!      write(*,'(100F8.5)') (real(edgec(i,j)*conjg(edgec(i,j))),j=1,nl/2)
!   end do
!end if

!-------------------------------------------------------
!Try to find the blacs code or use openmp to parallelize
!-------------------------------------------------------
!exp_half_h0=dcmplx(0.d0,0.d0)
!exp_mhalf_h0=dcmplx(0.d0,0.d0)
!exp_h0=dcmplx(0.d0,0.d0)
!
!do i=1,nl,1
!   do j=1,nl,1
!      do k=1,nl,1
!         exp_half_h0(i,j)=exp_half_h0(i,j)+hu(i,k)*dcmplx(exp(-dt/2.d0*ev(k)))*conjg(hu(j,k))
!         exp_mhalf_h0(i,j)=exp_mhalf_h0(i,j)+hu(i,k)*dcmplx(exp(dt/2.d0*ev(k)))*conjg(hu(j,k))
!         exp_h0(i,j)=exp_h0(i,j)+hu(i,k)*dcmplx(exp(-dt*ev(k)))*conjg(hu(j,k))
!      end do
!   end do
!end do
!
!write(*,*) exp_h0(7,8)

do k=1,nl,1
   exp_he(k)=exp(-dt/2.d0*ev(k))
end do
do i=1,nl,1
   do k=1,nl,1
      htmp(i,k)=hu(i,k)*dcmplx(exp_he(k))
   end do
end do
call zgemm('N','C',nl,nl,nl,one,htmp,nl,hu,nl,zero,exp_half_h0,nl)

do k=1,nl,1
   exp_mhe(k)=exp(dt/2.d0*ev(k))
end do
do i=1,nl,1
   do k=1,nl,1
      htmp(i,k)=hu(i,k)*dcmplx(exp_mhe(k))
   end do
end do
call zgemm('N','C',nl,nl,nl,one,htmp,nl,hu,nl,zero,exp_mhalf_h0,nl)


do k=1,nl,1
   exp_e(k)=exp(-dt*ev(k))
end do
do i=1,nl,1
   do k=1,nl,1
      htmp(i,k)=hu(i,k)*dcmplx(exp_e(k))
   end do
end do
call zgemm('N','C',nl,nl,nl,one,htmp,nl,hu,nl,zero,exp_h0,nl)


do k=1,nl,1
   exp_me(k)=exp(dt*ev(k))
end do
do i=1,nl,1
   do k=1,nl,1
      htmp(i,k)=hu(i,k)*dcmplx(exp_me(k))
   end do
end do
call zgemm('N','C',nl,nl,nl,one,htmp,nl,hu,nl,zero,exp_mh0,nl)

!write(*,*) exp_h0(7,8)

end subroutine initial_K


!--------------------------------------------------
!used to call initial_K for the decoupled condition
!--------------------------------------------------
subroutine initial_K_d
use all_param
implicit none
complex(kind=8)::h0(Nsite,Nsite)
complex(kind=8)::exp_h0(Nsite,Nsite)
complex(kind=8)::exp_mh0(Nsite,Nsite)
complex(kind=8)::exp_half_h0(Nsite,Nsite)
complex(kind=8)::exp_mhalf_h0(Nsite,Nsite)
real(kind=8)::exp_he0(Nsite),exp_mhe0(Nsite),exp_e0(Nsite),exp_me0(Nsite)
integer::i,j,k


!allocate(h0(Nsite,Nsite),exp_h0(Nsite,Nsite),exp_half_h0(Nsite,Nsite),exp_mhalf_h0(Nsite,Nsite))

h0(1:Nsite,1:Nsite)=Hzero(1:Nsite,1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(1,1),2*Nsite,h0(1,1),Nsite)
call initial_K(Nsite,Ntot,h0,dt,exp_h0,exp_mh0,exp_half_h0,exp_mhalf_h0,exp_he0,exp_mhe0,exp_e0,exp_me0,rank)
exp_K(1:Nsite,1:Nsite)=exp_h0(1:Nsite,1:Nsite)
exp_mK(1:Nsite,1:Nsite)=exp_mh0(1:Nsite,1:Nsite)
exp_halfK(1:Nsite,1:Nsite)=exp_half_h0(1:Nsite,1:Nsite)
exp_mhalfK(1:Nsite,1:Nsite)=exp_mhalf_h0(1:Nsite,1:Nsite)
exp_he(1:Nsite)=exp_he0(1:Nsite)
exp_mhe(1:Nsite)=exp_mhe0(1:Nsite)
exp_e(1:Nsite)=exp_e0(1:Nsite)
exp_me(1:Nsite)=exp_me0(1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_h0(1,1),Nsite,exp_K(1,1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_half_h0(1,1),Nsite,exp_halfK(1,1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_mhalf_h0(1,1),Nsite,exp_mhalfK(1,1),2*Nsite)


h0(1:Nsite,1:Nsite)=Hzero((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),Hzero(Nsite+1,Nsite+1),2*Nsite,h0(1,1),Nsite)
call initial_K(Nsite,Ntot,h0,dt,exp_h0,exp_mh0,exp_half_h0,exp_mhalf_h0,exp_he0,exp_mhe0,exp_e0,exp_me0,rank)
exp_K((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_h0(1:Nsite,1:Nsite)
exp_mK((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_mh0(1:Nsite,1:Nsite)
exp_halfK((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_half_h0(1:Nsite,1:Nsite)
exp_mhalfK((Nsite+1):(2*Nsite),(Nsite+1):(2*Nsite))=exp_mhalf_h0(1:Nsite,1:Nsite)
exp_he((Nsite+1):(2*Nsite))=exp_he0(1:Nsite)
exp_mhe((Nsite+1):(2*Nsite))=exp_mhe0(1:Nsite)
exp_e((Nsite+1):(2*Nsite))=exp_e0(1:Nsite)
exp_me((Nsite+1):(2*Nsite))=exp_me0(1:Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_h0(1,1),Nsite,exp_K(Nsite+1,Nsite+1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_half_h0(1,1),Nsite,exp_halfK(Nsite+1,Nsite+1),2*Nsite)
!call mkl_zomatcopy('C','N',Nsite,Nsite,(1.d0,0.d0),exp_mhalf_h0(1,1),Nsite,exp_mhalfK(Nsite+1,Nsite+1),2*Nsite)


!deallocate(h0,exp_h0,exp_half_h0,exp_mhalf_h0)
end subroutine initial_K_d

!----------------------
!test with mathematica
!subroutine test()
!implicit none
!integer,parameter::nl=2
!complex(kind=8)::h0(2*nl,2*nl)
!real(kind=8)::dt=0.01d0
!complex(kind=8)::exp_h0(2*nl,2*nl),exp_half_h0(2*nl,2*nl),exp_mhalf_h0(2*nl,2*nl)
!h0(1,1)=dcmplx(1.d0,0.d0);h0(1,2)=dcmplx(0.d0,0.d0);h0(1,3)=dcmplx(1.d0,-1.d0);h0(1,4)=dcmplx(2.d0,-1.d0)
!h0(2,1)=dcmplx(0.d0);h0(2,2)=dcmplx(2.d0);h0(2,3)=dcmplx(0);h0(2,4)=dcmplx(3.d0,1.d0)
!h0(3,1)=dcmplx(1.d0,1.d0);h0(3,2)=dcmplx(0);h0(3,3)=dcmplx(3.d0);h0(3,4)=dcmplx(2.d0)
!h0(4,1)=dcmplx(2.d0,1.d0);h0(4,2)=dcmplx(3.d0,-1.d0);h0(4,3)=dcmplx(2.d0);h0(4,4)=dcmplx(4.d0)
!call initial_K(2*nl,h0,dt,exp_h0,exp_half_h0,exp_mhalf_h0)
!write(*,*) exp_mhalf_h0;call mystop
!end subroutine test
!end test
!-----------------------


!use the fourier tranformation to get the K eigenvalue matrix
subroutine initial_K_fourier()
use all_param
implicit none
integer::i,j,k,Ns
real(kind=8)::tmp
real(kind=8)::eig,esoc
real(kind=8)::kx,ky
real(kind=8)::hkx,hky
real(kind=8)::lit=1d-10

if(dtype.ne.'w') then
   Ns=Nsite
elseif(dtype.eq.'w') then
   Ns=Nbravais
endif
  
do i=1,Ns,1

   !Get the momentum k
   kx=(dble(coor(i,1)-1)+kbound(1))*2*Pi/Nl(1)
   ky=(dble(coor(i,2)-1)+kbound(2))*2*Pi/Nl(2)
   if(kx.GT.Pi) kx=kx-2.d0*Pi
   if(ky.GT.Pi) ky=ky-2.d0*Pi

   if(dtype.ne.'w') then
      eig=0.d0
      do j=1,Dimen,1
         eig=eig+cos((dble(coor(i,j)-1)+kbound(j))*2*Pi/Nl(j))
      end do
      eig=2.d0*dble(t1)*eig
      !eig=kx*kx+ky*ky

      tmp=0.d0
      do j=1,Dimen,1
         tmp=tmp+sin((dble(coor(i,j)-1)+kbound(j))*2*Pi/Nl(j))**2
      end do
      esoc=2.d0*lamda*sqrt(tmp) 
      !esoc=lamda*sqrt(kx*kx+ky*ky)

      if(abs(kx).LT.lit.AND.abs(ky).LT.lit) then
         uk(i)=-1.d0/sqrt(2.d0)
      else
         uk(i)=-1.d0*(sin(ky)+Xi*sin(kx))/(sqrt(2.d0)*abs((sin(ky)-Xi*sin(kx))))
      end if
      vk(i)=1.d0/sqrt(2.d0)

      !if(abs(kx).LT.lit.AND.abs(ky).LT.lit) then
      !  uk(i)=-1.d0/sqrt(2.d0)
      !else
      !  uk(i)=-1.d0*(ky+Xi*kx)/(sqrt(2.d0)*abs(ky-Xi*kx))
      !end if
      !vk(i)=1.d0/sqrt(2.d0)
   elseif(dtype.eq.'w') then
      eig=2.d0*dble(tyhop)*cos(ky)

      hkx=-1.d0*(dble(vhop)+2.d0*dble(tdhop)*cos(ky)) &
           & -1.d0*(dble(whop)+2.d0*dble(tdhop)*cos(ky))*cos(kx)
      hky=-1.d0*(dble(whop)+2.d0*dble(tdhop)*cos(ky))*sin(kx)
      esoc=sqrt(hkx**2+hky**2)

      if(abs(kx).LT.lit.AND.abs(ky).LT.lit) then
         uk(i)=-1.d0/sqrt(2.d0)
      else
         uk(i)=-1.d0*(hkx-Xi*hky)/(sqrt(2.d0*(hkx**2+hky**2)))
      end if
      vk(i)=1.d0/sqrt(2.d0)

      if((abs(hkx).lt.1d-8).and.(abs(hky).lt.1d-8)) then
         uk(i)=1/sqrt(2.d0)
         vk(i)=1/sqrt(2.d0)
      endif
   endif
      

   !set the diagnoal matrix
   if(dtype.ne.'w') then
      exp_he(i)=exp(-dt/2.d0*(eig-esoc))
      exp_he(i+Nsite)=exp(-dt/2.d0*(eig+esoc))
      exp_mhe(i)=exp(dt/2.d0*(eig-esoc))
      exp_mhe(i+Nsite)=exp(dt/2.d0*(eig+esoc))
      exp_e(i)=exp(-dt*(eig-esoc))
      exp_e(i+Nsite)=exp(-dt*(eig+esoc))
      exp_me(i)=exp(dt*(eig-esoc))
      exp_me(i+Nsite)=exp(dt*(eig+esoc))
   elseif(dtype.eq.'w') then
      exp_he(i)=exp(-dt/2.d0*(eig-esoc))
      exp_he(i+Nbravais)=exp(-dt/2.d0*(eig+esoc))
      exp_he(i+Nsite)=exp(-dt/2.d0*(eig-esoc))
      exp_he(i+Nsite+Nbravais)=exp(-dt/2.d0*(eig+esoc))
      
      exp_mhe(i)=exp(dt/2.d0*(eig-esoc))
      exp_mhe(i+Nbravais)=exp(dt/2.d0*(eig+esoc))
      exp_mhe(i+Nsite)=exp(dt/2.d0*(eig-esoc))
      exp_mhe(i+Nsite+Nbravais)=exp(dt/2.d0*(eig+esoc))
      
      exp_e(i)=exp(-dt*(eig-esoc))
      exp_e(i+Nbravais)=exp(-dt*(eig+esoc))
      exp_e(i+Nsite)=exp(-dt*(eig-esoc))
      exp_e(i+Nsite+Nbravais)=exp(-dt*(eig+esoc))
      
      exp_me(i)=exp(dt*(eig-esoc))
      exp_me(i+Nbravais)=exp(dt*(eig+esoc))
      exp_me(i+Nsite)=exp(dt*(eig-esoc))
      exp_me(i+Nsite+Nbravais)=exp(dt*(eig+esoc))
   endif
      
end do
end subroutine initial_K_fourier



!--------------------------------------------------------
!We initial the gamaf which is used in H-S transformation
!Also initial the ng(Nsite) in the free projection
!--------------------------------------------------------
subroutine initial_V()
use all_param
implicit none
real(kind=8)::y
integer::i,j


!For the gamaf in Free projection
if(kcrn.eq.1) then
  y=dt*onsitU/2.d0
else if(kcrn.eq.2) then
  y=-1.d0*dt*onsitU/2.d0
else if(kcrn.eq.3.or.kcrn.eq.4) then
  y=0.d0
else
  write(*,*) "Do not know what kind of decouple method in gama free projection."
  call mystop
end if
call get_gama(gamaf,y)


!For the ng background in Free projection
if(kcrn.eq.3.or.kcrn.eq.1) then
  do i=1,Nsite,1
     ng(i)=dble(Nspin(1)-Nspin(2))/dble(Nsite)
  end do
else if(kcrn.eq.4.or.kcrn.eq.2) then
  do i=1,Nsite,1
     ng(i)=dble(Ntot)/dble(Nsite)
  end do
end if
end subroutine initial_V



!---------------------------------
!solve the equation cosh(x)=exp(y)
!---------------------------------
subroutine get_gama(x,y)
implicit none
real(kind=8),intent(IN)::y
complex(kind=8),intent(OUT)::x
complex(kind=8)::ey,rdummy
if(abs(y)<1d-10) then
   x=dcmplx(0.d0,0.d0)
else
   ey=dcmplx(exp(y))
   x=log(ey-sqrt(ey**2-1.d0))
endif

!Check x
rdummy=log((exp(x)+exp(-x))/2d0)
if(abs(rdummy-y)<1d-10) then
  !print *, 'gamma checked!'
else
  print *, 'wrong x!'
  call mystop
endif
end subroutine get_gama


!-------------------------------------------------------------
!This subroutine allocate the arrays in the projection K and V
!------------------------------------------------------------- 
subroutine allocate_project_array()
use lattice_param
use project_param
implicit none
if(pfft.eq.0) then
  allocate(exp_halfK(DNsite,DNsite))
  allocate(exp_mhalfK(DNsite,DNsite))
  allocate(exp_K(DNsite,DNsite))
  allocate(exp_mK(DNsite,DNsite))
end if
allocate(exp_he(DNsite))
allocate(exp_mhe(DNsite))
allocate(exp_e(DNsite))
allocate(exp_me(DNsite))
allocate(ng(Nsite))
allocate(uk(Nsite))
allocate(vk(Nsite))
end subroutine allocate_project_array


!---------------------------------------------------------------
!This subroutine deallocate the arrays in the projection K and V
!--------------------------------------------------------------- 
subroutine deallocate_project_array()
use project_param
implicit none
if(allocated(exp_halfK)) deallocate(exp_halfK)
if(allocated(exp_mhalfK)) deallocate(exp_mhalfK)
if(allocated(exp_K)) deallocate(exp_K)
if(allocated(exp_mK)) deallocate(exp_mK)
if(allocated(exp_he)) deallocate(exp_he)
if(allocated(exp_mhe)) deallocate(exp_mhe)
if(allocated(exp_e)) deallocate(exp_e)
if(allocated(exp_me)) deallocate(exp_me)
if(allocated(ng)) deallocate(ng)
if(allocated(uk)) deallocate(uk)
if(allocated(vk)) deallocate(vk)
end subroutine deallocate_project_array
