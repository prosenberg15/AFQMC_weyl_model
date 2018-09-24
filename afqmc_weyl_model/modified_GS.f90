!This subroutine do the Modified GS periodically.
!subroutine Modified_GS()
!use all_param
!implicit none
!integer::i
!do i=1,Nwalkers,1
!   call modGS_i(i)
!end do
!end subroutine Modified_GS


!subroutine modGS_i(i)
!use all_param
!implicit none
!integer,intent(IN)::i
!real(kind=8)::anm,anm1,anm2
!complex(kind=8)::Rmat(Ntot,Ntot)
!
!!For the right phi
!if(dtype.EQ.'c') then
!  call modGS(phi(1,1,i,1),2*Nsite,Ntot,anm)
!else if(dtype.EQ.'d') then
!  call modGS(phi(1:Nsite,1:Nspin(1),i,1),Nsite,Nspin(1),anm1)
!  call modGS(phi((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,i,1),Nsite,Nspin(2),anm2)
!  anm=anm1*anm2
!end if
!if(anm.LT.0.d0) then
!  write(*,*) "Something is wrong in Modified_GS."
!  write(*,*) "anm:",anm
!  call mystop
!end if
!if(weight(i).le.0.d0) write(*,*) "negative weight!"
!weight(i)=weight(i)/anm
!dlogw(i)=dlogw(i)-dlog(anm)
!
!
!!For the left phi
!if(dtype.EQ.'c') then
!  call modGS(phi(1,1,i,2),2*Nsite,Ntot,anm)
!else if(dtype.EQ.'d') then
!  call modGS(phi(1:Nsite,1:Nspin(1),i,2),Nsite,Nspin(1),anm1)
!  call modGS(phi((Nsite+1):(2*Nsite),(Nspin(1)+1):Ntot,i,2),Nsite,Nspin(2),anm2)
!  anm=anm1*anm2
!end if
!if(anm.LT.0.d0) then
!  write(*,*) "Something is wrong in Modified_GS."
!  write(*,*) "anm:",anm
!  call mystop
!end if
!if(weight(i).le.0.d0) write(*,*) "negative weight!"
!weight(i)=weight(i)/anm
!dlogw(i)=dlogw(i)-dlog(anm)
!end subroutine modGS_i


!use QR to do the modified GS
subroutine modGS(ph,L,N,det,Rmax)
implicit none
integer,intent(IN)::L,N
complex(kind=8),intent(INOUT)::ph(L,N)
real(kind=8),intent(OUT)::det
complex(kind=8),optional::Rmax(N,N)

complex(kind=8)::tau(N)
complex(kind=8),allocatable::work(:)
integer::lwork
integer::info

integer::i


allocate(work(1))
call zgeqrf(L,N,ph,L,tau,work,-1,info)
lwork=work(1)
deallocate(work)
if(info.NE.0) then
  write(*,*) "Something is wrong in QR:",info
end if


allocate(work(lwork))
call zgeqrf(L,N,ph,L,tau,work,lwork,info)
if(info.NE.0) then
  write(*,*) "Something is wrong in QR:",info
end if


det=1.d0
do i=1,N,1
   det=det/ph(i,i)
end do
call zungqr(L,N,N,ph,L,tau,work,lwork,info)
deallocate(work)


if(det.LT.0.d0) then
  det=-1.d0*det
  do i=1,L,1
     ph(i,1)=ph(i,1)*(-1.d0)
  end do
end if
end subroutine modGS
