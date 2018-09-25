subroutine meas_energy(phL_r,phR_r)
use all_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

complex(kind=8),intent(IN)::phL_r(DNsite,Ntot),phR_r(DNsite,Ntot)
complex(kind=8)::phL_k(DNsite,Ntot),phR_k(DNsite,Ntot)
complex(kind=8)::Amat_local_dr(DNsite),Amat_local_dk(DNsite),Amat_local_ofdr(DNsite)
complex(kind=8),allocatable::Amat_local(:,:)
complex(kind=8)::ovp_local(Ntot,Ntot)
complex(kind=8)::imp_local
complex(kind=8)::Kinm_local,Vinm_local,E_local
complex(kind=8)::Nup_local,Ndn_local
complex(kind=8)::tmp
integer::mi,mj

if(diagm.eq.0) then
  allocate(Amat_local(DNsite,DNsite))
end if

i_energy=i_energy+1


!Get the inverse,imp_local,tot_local and Amat
 call over_lap_dc(phL_r(1,1),phR_r(1,1),ovp_local(1,1))
 call linear_equation_d_dc(ovp_local(1:Ntot,1:Ntot),imp_local,phR_r(1:2*Nsite,1:Ntot))
 if(diagm.eq.0) then
    call cal_Amat_step_dc(phL_r(1,1),phR_r(1,1),Amat_local(1,1))
 else
    call cal_Amat_diag_dc(phL_r(1,1),phR_r(1,1),Amat_local_dr(1))
    !natural orbital space
    call ph_dc_fftw_f(phR_r(1:DNsite,1:Ntot),phR_k(1:DNsite,1:Ntot),uk,vk)
    call ph_dc_fftw_f(phl_r(1:DNsite,1:Ntot),phl_k(1:DNsite,1:Ntot),uk,vk)
    call cal_Amat_diag_dc(phL_k(1,1),phR_k(1,1),Amat_local_dk(1))
 end if


 !Kinetic energy
 if(diagm.eq.0) then
    if(dtype.ne.'w') then
       Kinm_local=zero
       do mi=1,Nhop,1
          Kinm_local=Kinm_local+hopt(mi)*Amat_local(sit(mi,1),sit(mi,2))
       end do
    elseif(dtype.eq.'w') then
       Kinm_local=zero
       do mi=1,Nhop,1
          Kinm_local=Kinm_local+hopt(mi)*Amat_local(sit(mi,1),sit(mi,2))
          Kinm_local=Kinm_local+hopt(mi)*Amat_local(sit(mi,1)+Nsite,sit(mi,2)+Nsite)
       end do
    endif
 else
   Kinm_local=zero
   do mi=1,DNsite,1
      Kinm_local=Kinm_local-dlog(exp_e(mi))*Amat_local_dk(mi)/dt
      !if(rank.eq.0) then
      !   write(*,*) mi, -dlog(exp_e(mi))/dt, dble(Amat_local_dk(mi)), Kinm_local
      !endif
   end do
 end if


 !Potential energy
 if(diagm.eq.0) then
   Vinm_local=zero
   if(dtype.EQ.'c') then
     do mi=1,Nsite,1
        Vinm_local=Vinm_local+Amat_local(mi,mi)*Amat_local(mi+Nsite,mi+Nsite) &
                   & -Amat_local(mi,mi+Nsite)*Amat_local(mi+Nsite,mi)
     end do
   else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
     do mi=1,Nsite,1
        Vinm_local=Vinm_local+Amat_local(mi,mi)*Amat_local(mi+Nsite,mi+Nsite)
     end do
   end if
 else
   Vinm_local=zero
   if(dtype.EQ.'c') then
     !Set off diagonal term first
     Amat_local_ofdr=zero
     do mi=1,Nsite,1
       do mj=1,Ntot,1
          Amat_local_ofdr(mi)=Amat_local_ofdr(mi)+phR_r(mi,mj)*conjg(phL_r(mi+Nsite,mj))
          Amat_local_ofdr(mi+Nsite)=Amat_local_ofdr(mi+Nsite)+phR_r(mi+Nsite,mj)*conjg(phL_r(mi,mj))
       end do
     end do
     do mi=1,Nsite,1
        Vinm_local=Vinm_local+Amat_local_dr(mi)*Amat_local_dr(mi+Nsite) &
                  & -Amat_local_ofdr(mi)*Amat_local_ofdr(mi+Nsite)
     end do
   else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
     do mi=1,Nsite,1
        Vinm_local=Vinm_local+Amat_local_dr(mi)*Amat_local_dr(mi+Nsite)
     end do
   end if
 end if
 Vinm_local=dcmplx(onsitU)*Vinm_local

 E_local=Kinm_local+Vinm_local

 call kahan_sum_c(E_local,E_c,E_one)

 !if(rank.eq.0) write(*,*) Kinm_local,Vinm_local, i_energy, E_one

 !Nup_local and Ndn_local
 !if(diagm.eq.0) then
 !  Nup_local=zero
 !  Ndn_local=zero
 !  do mi=1,Nsite,1
 !     Nup_local=Nup_local+Amat_local(mi,mi)
 !     Ndn_local=Ndn_local+Amat_local(mi+Nsite,mi+Nsite)
 !  end do
 !else
 !  Nup_local=zero
 !  Ndn_local=zero
 !  do mi=1,Nsite,1
 !     Nup_local=Nup_local+Amat_local_dk(mi)
 !     Ndn_local=Ndn_local+Amat_local_dk(mi+Nsite)
 !  end do
 !end if
 !call kahan_sum_c(Nup_local,Nup_c,Nup_one)
 !call kahan_sum_c(Ndn_local,Ndn_c,Ndn_one)


if(diagm.eq.0) then
  deallocate(Amat_local)
end if
end subroutine meas_energy



subroutine meas_observ(phL_r,phR_r)
use all_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

complex(kind=8),intent(IN)::phL_r(DNsite,Ntot),phR_r(DNsite,Ntot)
complex(kind=8)::phL_k(DNsite,Ntot),phR_k(DNsite,Ntot)
complex(kind=8),allocatable::Amat_local(:,:),Amat_local_k(:,:)
complex(kind=8)::ovp_local(Ntot,Ntot)
complex(kind=8)::imp_local
complex(kind=8)::Kinm_local,Vinm_local,E_local
complex(kind=8)::Nup_local,Ndn_local
complex(kind=8)::hkx,hky
real(kind=8)::kx,ky
integer::mi,mj,Ns

allocate(Amat_local(DNsite,DNsite),Amat_local_k(DNsite,DNsite))

i_energy=i_energy+1
i_observ=i_observ+1

!Get the imp_local,tot_local and Amat
!For real space
 call over_lap_dc(phL_r(1,1),phR_r(1,1),ovp_local(1,1))
 call linear_equation_d_dc(ovp_local(1:Ntot,1:Ntot),imp_local,phR_r(1:DNsite,1:Ntot))
 call cal_Amat_step_dc(phL_r(1,1),phR_r(1,1),Amat_local(1,1))
!For K space true k space, not natural orbital space
 if((openbcx.eq.0).and.(openbcy.eq.0))then
       call ph_fftw_f(Nsite,Ntot,phR_r(1:Nsite,1:Ntot),phR_k(1:Nsite,1:Ntot))
       call ph_fftw_f(Nsite,Ntot,phR_r((Nsite+1):DNsite,1:Ntot),phR_k((Nsite+1):DNsite,1:Ntot))
       call ph_fftw_f(Nsite,Ntot,phl_r(1:Nsite,1:Ntot),phl_k(1:Nsite,1:Ntot))
       call ph_fftw_f(Nsite,Ntot,phl_r((Nsite+1):DNsite,1:Ntot),phl_k((Nsite+1):DNsite,1:Ntot))
       call cal_Amat_step_dc(phL_k(1,1),phR_k(1,1),Amat_local_k(1,1))
 end if

 !Kinectic energy
 !Kinm_local=zero
 !do mi=1,DNsite,1
 !   Kinm_local=Kinm_local-dlog(exp_e(mi))*Amat_local_k(mi,mi)/dt
 !end do
 !call kahan_sum_c(Kinm_local,K_c,K_one)
 if((openbcx.eq.0).and.(openbcy.eq.0))then
    Kinm_local=zero
    if(dtype.ne.'w') then
       Ns=Nsite
    elseif(dtype.eq.'w') then
       Ns=Nbravais
    endif
    do mi=1,Ns,1
       !Get the momentum k
       kx=(dble(coor(mi,1)-1)+kbound(1))*2*Pi/Nl(1)
       ky=(dble(coor(mi,2)-1)+kbound(2))*2*Pi/Nl(2)
       if(kx.GT.Pi) kx=kx-2.d0*Pi
       if(ky.GT.Pi) ky=ky-2.d0*Pi
       !Kinm_local=Kinm_local+(kx*kx+ky*ky)*(Amat_local_k(mi,mi)+Amat_local_k(mi+Nsite,mi+Nsite))
       !Kinm_local=Kinm_local+(ky+Xi*kx)*lamda*Amat_local_k(mi,mi+Nsite)
       !Kinm_local=Kinm_local+(ky-Xi*kx)*lamda*Amat_local_k(mi+Nsite,mi)
       if(dtype.ne.'w') then
          Kinm_local=Kinm_local+2.d0*dble(t1)*(cos(kx)+cos(ky))*(Amat_local_k(mi,mi)+Amat_local_k(mi+Nsite,mi+Nsite))
          Kinm_local=Kinm_local+2*(sin(ky)+Xi*sin(kx))*lamda*Amat_local_k(mi,mi+Nsite)
          Kinm_local=Kinm_local+2*(sin(ky)-Xi*sin(kx))*lamda*Amat_local_k(mi+Nsite,mi)
       else if(dtype.eq.'w') then
          hkx=-1.d0*(dble(vhop)+2.d0*dble(tdhop)*cos(ky)) &
           & -1.d0*(dble(whop)+2.d0*dble(tdhop)*cos(ky))*cos(kx)
          hky=-1.d0*(dble(whop)+2.d0*dble(tdhop)*cos(ky))*sin(kx)
          Kinm_local=Kinm_local+2.d0*dble(tyhop)*cos(ky)*(Amat_local_k(mi,mi)+Amat_local_k(mi+Nsite,mi+Nsite) &
               & + Amat_local_k(mi+Nbravais,mi+Nbravais)+Amat_local_k(mi+Nbravais+Nsite,mi+Nbravais+Nsite))
          Kinm_local=Kinm_local+(hkx-Xi*hky)*(Amat_local_k(mi,mi+Nbravais)+Amat_local_k(mi+Nsite,mi+Nbravais+Nsite))
          Kinm_local=Kinm_local+(hkx+Xi*hky)*(Amat_local_k(mi+Nbravais,mi)+Amat_local_k(mi+Nbravais+Nsite,mi+Nsite))
          !if(rank.eq.0) then
          !   write(*,*) kx,ky,dble(Amat_local_k(mi,mi)),dble(Amat_local_k(mi+Nsite,mi+Nsite)),&
          !        & dble(Amat_local_k(mi+Nbravais,mi+Nbravais)),dble(Amat_local_k(mi+Nbravais+Nsite,mi+Nbravais+Nsite))
          !   write(*,*) kx,ky,hkx-Xi*hky,Amat_local_k(mi,mi+Nbravais),Amat_local_k(mi+Nsite,mi+Nbravais+Nsite)
          !   write(*,*) kx,ky,hkx+Xi*hky,Amat_local_k(mi+Nbravais,mi),Amat_local_k(mi+Nbravais+Nsite,mi+Nsite)
          !endif
       endif
    end do
    !if(rank.eq.0) write(*,*) 'total K: ', Kinm_local
    call kahan_sum_c(Kinm_local,K_c,K_one)
 else if((openbcx.eq.1).or.(openbcy.eq.1)) then
   Kinm_local=zero
   do mi=1,Nhop,1
      Kinm_local=Kinm_local+hopt(mi)*Amat_local(sit(mi,1),sit(mi,2))
      Kinm_local=Kinm_local+hopt(mi)*Amat_local(sit(mi,1)+Nsite,sit(mi,2)+Nsite)
   end do
   if(rank.eq.0) write(*,*) 'total K: ', Kinm_local
   call kahan_sum_c(Kinm_local,K_c,K_one)
end if

 !Potential energy
 Vinm_local=zero
 if(dtype.EQ.'c') then
   do mi=1,Nsite,1
      Vinm_local=Vinm_local+Amat_local(mi,mi)*Amat_local(mi+Nsite,mi+Nsite) &
                 & -Amat_local(mi,mi+Nsite)*Amat_local(mi+Nsite,mi)
   end do
 else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
   do mi=1,Nsite,1
      Vinm_local=Vinm_local+Amat_local(mi,mi)*Amat_local(mi+Nsite,mi+Nsite)
   end do
 end if
 Vinm_local=dcmplx(onsitU)*Vinm_local
 call kahan_sum_c(Vinm_local,V_c,V_one)


 E_local=Kinm_local+Vinm_local
 call kahan_sum_c(E_local,E_c,E_one)


 !Nup_local and Ndn_local
 Nup_local=zero
 Ndn_local=zero
 do mi=1,Nsite,1
    Nup_local=Nup_local+Amat_local(mi,mi)
    Ndn_local=Ndn_local+Amat_local(mi+Nsite,mi+Nsite)
 end do
 call kahan_sum_c(Nup_local,Nup_c,Nup_one)
 call kahan_sum_c(Ndn_local,Ndn_c,Ndn_one)

 if((openbcx.eq.1).or.(openbcy.eq.1)) then
    call add_edgec(Amat_local)

    !density
    call add_ni(Amat_local)
 end if

 if(((openbcx.eq.1).and.(openbcy.eq.1)) & 
      .or. ((openbcx.eq.0).and.(openbcy.eq.0))) then
    !singlet pairing
    call add_didj(Amat_local)

    !spin spin correlation
    call add_sisj(Amat_local)

    !charge charge correlation
    call add_ninj(Amat_local)
 endif

 if((openbcx.eq.0).and.(openbcy.eq.0)) then
    !momentum distribution
    call add_ck(Amat_local_k)

    !spin direction
    call add_sk(Amat_local_k)

    !for the pairing wf 
    call add_pair_one_body(Amat_local_k)
    deallocate(Amat_local,Amat_local_k)
 end if
end subroutine meas_observ


subroutine add_didj(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::didj_local
integer,external::latt_label
integer,external::bound
integer::cc(1:Dimen),ctmp
integer::i,j,n,m

if((openbcx.eq.0).and.(openbcy.eq.0)) then
   if(dtype.EQ.'c') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite) &
                 & -Amat_local(i,j+Nsite)*Amat_local(i+Nsite,j)
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))
            call kahan_sum_c(didj_local,didj_c(n),didj_one(n))  
         end do
      end do
   else if(dtype.EQ.'d') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite)
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))
            call kahan_sum_c(didj_local,didj_c(n),didj_one(n))
         end do
      end do
   else if(dtype.EQ.'w') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite)
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))
            !d^{dagger,A}_i d^A_j
            if((i.le.Nbravais).and.(j.le.Nbravais)) then
               call kahan_sum_c(didj_local,didj_c(n),didj_one(n))
            !d^{dagger,A}_i d^B_j
            else if((i.le.Nbravais).and.(j.gt.Nbravais)) then
               call kahan_sum_c(didj_local,didj_c(n+Nbravais),didj_one(n+Nbravais))
            !d^{dagger,B}_i d^A_j
            else if((i.gt.Nbravais).and.(j.le.Nbravais)) then
               call kahan_sum_c(didj_local,didj_c(n+2*Nbravais),didj_one(n+2*Nbravais))
            !d^{dagger,B}_i d^B_j
            else if((i.gt.Nbravais).and.(j.gt.Nbravais)) then
               call kahan_sum_c(didj_local,didj_c(n+3*Nbravais),didj_one(n+3*Nbravais))
            endif
         end do
      end do
   end if
! OPEN BC ALONG X, PERIODIC ALONG Y
else if((openbcx.eq.1).and.(openbcy.eq.0)) then
   if(dtype.EQ.'c') then
      i=1
      do j=1,Nsite,1
         didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite) &
              & -Amat_local(i,j+Nsite)*Amat_local(i+Nsite,j)
         call kahan_sum_c(didj_local,didj_c(j),didj_one(j))  
      end do

      !do i=1,Nsite,1
      !   do j=1,Nsite,1
      !      didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite) &
      !           & -Amat_local(i,j+Nsite)*Amat_local(i+Nsite,j)
      !      do m=1,Dimen,1
      !         !We calculate didj==>d1d(j-i+1)
      !         !so we need to focus on j-i+1
      !         if(m.eq.1) then
      !            ctmp=coor(j,m)-coor(i,m)+1
      !            cc(m)=bound(ctmp,Nl(m))
      !         else
      !            cc(m)=coor(j,m)-coor(i,m)+1
      !         end if
      !      end do
      !      n=latt_label(cc(1:Dimen))
      !      call kahan_sum_c(didj_local,didj_c(n),didj_one(n))  
      !   end do
      !end do
   else if(dtype.EQ.'d') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite)
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               if(m.eq.1) then
                  ctmp=coor(j,m)-coor(i,m)+1
                  cc(m)=bound(ctmp,Nl(m))
               else
                  cc(m)=coor(j,m)-coor(i,m)+1
               end if
            end do
            n=latt_label(cc(1:Dimen))
            call kahan_sum_c(didj_local,didj_c(n),didj_one(n))
         end do
      end do
   end if
else if ((openbcx.eq.1).and.(openbcy.eq.1)) then
      if(dtype.EQ.'c') then
         i=1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite) &
                 & -Amat_local(i,j+Nsite)*Amat_local(i+Nsite,j)
            call kahan_sum_c(didj_local,didj_c(j),didj_one(j))  
         end do
      else if(dtype.EQ.'d') then
         i=1
         do j=1,Nsite,1
            didj_local=Amat_local(i,j)*Amat_local(i+Nsite,j+Nsite)
            call kahan_sum_c(didj_local,didj_c(j),didj_one(j))
         end do
      end if
   end if
end subroutine add_didj


subroutine add_sisj(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::sisj_local
integer,external::latt_label
integer,external::bound
integer::cc(1:Dimen),ctmp
integer::i,j,n,m

if ((openbcx.eq.0).and.(openbcy.eq.0)) then
   if(dtype.EQ.'c') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            sisj_local=zero
            !term 1
            if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
            sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            !term 2
            if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
            sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)- &
                 & Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
            !term 3
            sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)+ &
                 &  Amat_local(i+Nsite,j)*Amat_local(j,i+Nsite)
            !term 4
            sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)+ &
                 &  Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
            !term 5
            if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
            sisj_local=sisj_local+2.d0*Amat_local(i,i+Nsite)*Amat_local(j+Nsite,j)- &
                 &  2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
            !term 6
            if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
            sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i)*Amat_local(j,j+Nsite)- &
                 &  2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
            sisj_local=sisj_local/4.d0
            do m=1,Dimen,1
               !We calculate sisj==>s1s(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))
            call kahan_sum_c(sisj_local,sisj_c(n),sisj_one(n))
         end do
      end do
      
   else if(dtype.EQ.'d') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            sisj_local=zero
            !term 1
            if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
            sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            !term 2
            if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
            sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)-&
                 &Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
            !term 3
            sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)
            
            !term 4
            sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
            
            !term 5
            if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
            sisj_local=sisj_local-2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
            
            !term 6
            if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
            sisj_local=sisj_local-2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
            sisj_local=sisj_local/4.d0
            do m=1,Dimen,1
               !We calculate sisj==>s1s(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))
            call kahan_sum_c(sisj_local,sisj_c(n),sisj_one(n))
         end do
      end do
   end if
! OPEN BC ALONG Y, PBC AlONG X
else if ((openbcx.eq.0).and.(openbcy.eq.1)) then
   if(dtype.EQ.'c') then
      i=1
      do j=1,Nsite,1
         sisj_local=zero
         !term 1
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
         sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         !term 2
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)- &
              & Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
         !term 3
         sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)+ &
              &  Amat_local(i+Nsite,j)*Amat_local(j,i+Nsite)
         !term 4
         sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)+ &
              &  Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
         !term 5
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
         sisj_local=sisj_local+2.d0*Amat_local(i,i+Nsite)*Amat_local(j+Nsite,j)- &
              &  2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
         !term 6
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i)*Amat_local(j,j+Nsite)- &
              &  2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
         sisj_local=sisj_local/4.d0
         do m=1,Dimen,1
            !We calculate sisj==>s1s(j-i+1)
            !so we need to focus on j-i+1
            if(m.eq.1) then
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            else
               cc(m)=coor(j,m)-coor(i,m)+1
            end if
         end do
         n=latt_label(cc(1:Dimen))
         call kahan_sum_c(sisj_local,sisj_c(n),sisj_one(n))
      end do
   else if(dtype.EQ.'d') then
      i=1
      do j=1,Nsite,1
         sisj_local=zero
         !term 1
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
         sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         !term 2
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)-&
              &Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
         !term 3
         sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)
         
         !term 4
         sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
         
         !term 5
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
         sisj_local=sisj_local-2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
         
         !term 6
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local-2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
         sisj_local=sisj_local/4.d0
         do m=1,Dimen,1
            !We calculate sisj==>s1s(j-i+1)
            !so we need to focus on j-i+1
            if(m.eq.1) then
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            else
               cc(m)=coor(j,m)-coor(i,m)+1
            end if
         end do
         n=latt_label(cc(1:Dimen))
         call kahan_sum_c(sisj_local,sisj_c(n),sisj_one(n))
      end do
   end if
else if ((openbcx.eq.1).and.(openbcy.eq.1)) then
   if(dtype.EQ.'c') then
      i=1
      do j=1,Nsite,1
         sisj_local=zero
         !term 1
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
         sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         !term 2
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)- &
              & Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
         !term 3
         sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)+ &
              &  Amat_local(i+Nsite,j)*Amat_local(j,i+Nsite)
         !term 4
         sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)+ &
              &  Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
         !term 5
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
         sisj_local=sisj_local+2.d0*Amat_local(i,i+Nsite)*Amat_local(j+Nsite,j)- &
              &  2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
         !term 6
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i)*Amat_local(j,j+Nsite)- &
              &  2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
         sisj_local=sisj_local/4.d0
         call kahan_sum_c(sisj_local,sisj_c(j),sisj_one(j))
      end do
   else if(dtype.EQ.'d') then
      i=1
      do j=1,Nsite,1
         sisj_local=zero
         !term 1
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i,i)
         sisj_local=sisj_local+Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         !term 2
         if(i.eq.j) sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local+Amat_local(i+Nsite,i+Nsite)*Amat_local(j+Nsite,j+Nsite)-&
              &Amat_local(i+Nsite,j+Nsite)*Amat_local(j+Nsite,i+Nsite)
         !term 3
         sisj_local=sisj_local-Amat_local(i+Nsite,i+Nsite)*Amat_local(j,j)
         
         !term 4
         sisj_local=sisj_local-Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
         
         !term 5
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i,i)
         sisj_local=sisj_local-2.d0*Amat_local(i,j)*Amat_local(j+Nsite,i+Nsite)
         
         !term 6
         if(i.eq.j) sisj_local=sisj_local+2.d0*Amat_local(i+Nsite,i+Nsite)
         sisj_local=sisj_local-2.d0*Amat_local(i+Nsite,j+Nsite)*Amat_local(j,i)
         sisj_local=sisj_local/4.d0
         call kahan_sum_c(sisj_local,sisj_c(j),sisj_one(j))
      end do
   end if
end if

end subroutine add_sisj


subroutine add_ninj(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::ninj_local
integer,external::latt_label
integer,external::bound
integer::cc(1:Dimen),ctmp
integer::i,j,n,m

if ((openbcx.eq.0).and.(openbcy.eq.0)) then
   if(dtype.EQ.'c') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))

            !Niu Nju
            if(i.eq.j) then
               ninj_local=Amat_local(i,i)
            else
               ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            end if
            call kahan_sum_c(ninj_local,ninj_c(n),ninj_one(n))

            !Niu Njd
            ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)-Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
            call kahan_sum_c(ninj_local,ninj_c(n+Nsite),ninj_one(n+Nsite))
         end do
      end do
   else if(dtype.EQ.'d') then
      do i=1,Nsite,1
         do j=1,Nsite,1

            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))

            !Niu Nju
            if(i.eq.j) then
               ninj_local=Amat_local(i,i)
            else
               ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            end if
            call kahan_sum_c(ninj_local,ninj_c(n),ninj_one(n))

            !Niu Njd
            ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
            call kahan_sum_c(ninj_local,ninj_c(n+Nsite),ninj_one(n+Nsite))
         end do
      end do
   else if(dtype.EQ.'w') then
      do i=1,Nsite,1
         do j=1,Nsite,1

            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               ctmp=coor(j,m)-coor(i,m)+1
               cc(m)=bound(ctmp,Nl(m))
            end do
            n=latt_label(cc(1:Dimen))

            !AA
            if((i.le.Nbravais).and.(j.le.Nbravais)) then
               !Niu Nju
               if(i.eq.j) then
                  ninj_local=Amat_local(i,i)
               else
                  ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
               end if
               call kahan_sum_c(ninj_local,ninj_c(n),ninj_one(n))

               !Niu Njd
               ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
               call kahan_sum_c(ninj_local,ninj_c(n+2*Nsite),ninj_one(n+2*Nsite))
            !AB   
            else if((i.le.Nbravais).and.(j.gt.Nbravais)) then
               !Niu Nju
               if(i.eq.j) then
                  ninj_local=Amat_local(i,i)
               else
                  ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
               end if
               call kahan_sum_c(ninj_local,ninj_c(n+Nbravais),ninj_one(n+Nbravais))

               !Niu Njd
               ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
               call kahan_sum_c(ninj_local,ninj_c(n+Nbravais+2*Nsite),ninj_one(n+Nbravais+2*Nsite))
            !BA
            else if((i.gt.Nbravais).and.(j.le.Nbravais)) then
               !Niu Nju
               if(i.eq.j) then
                  ninj_local=Amat_local(i,i)
               else
                  ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
               end if
               call kahan_sum_c(ninj_local,ninj_c(n+2*Nbravais),ninj_one(n+2*Nbravais))

               !Niu Njd
               ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
               call kahan_sum_c(ninj_local,ninj_c(n+2*Nbravais+2*Nsite),ninj_one(n+2*Nbravais+2*Nsite))
            !BB
            else if((i.gt.Nbravais).and.(j.gt.Nbravais)) then
               !Niu Nju
               if(i.eq.j) then
                  ninj_local=Amat_local(i,i)
               else
                  ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
               end if
               call kahan_sum_c(ninj_local,ninj_c(n+3*Nbravais),ninj_one(n+3*Nbravais))

               !Niu Njd
               ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
               call kahan_sum_c(ninj_local,ninj_c(n+3*Nbravais+2*Nsite),ninj_one(n+3*Nbravais+2*Nsite))
            end if
         end do
      end do
   end if
! OPEN BC ALONG Y, PBC ALONG X
else if ((openbcx.eq.0).and.(openbcy.eq.1)) then
   if(dtype.EQ.'c') then
      do i=1,Nsite,1
         do j=1,Nsite,1
            
            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               if(m.eq.1) then
                  ctmp=coor(j,m)-coor(i,m)+1
                  cc(m)=bound(ctmp,Nl(m))
               else
                 cc(m)=coor(j,m)-coor(i,m)+1
              end if
            end do
            n=latt_label(cc(1:Dimen))

            !Niu Nju
            if(i.eq.j) then
               ninj_local=Amat_local(i,i)
            else
               ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            end if
            call kahan_sum_c(ninj_local,ninj_c(n),ninj_one(n))

            !Niu Njd
            ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)-Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
            call kahan_sum_c(ninj_local,ninj_c(n+Nsite),ninj_one(n+Nsite))
         end do
      end do
   else if(dtype.EQ.'d') then
      do i=1,Nsite,1
         do j=1,Nsite,1

            do m=1,Dimen,1
               !We calculate didj==>d1d(j-i+1)
               !so we need to focus on j-i+1
               if(m.eq.1) then
                  ctmp=coor(j,m)-coor(i,m)+1
                  cc(m)=bound(ctmp,Nl(m))
               else
                  cc(m)=coor(j,m)-coor(i,m)+1
               end if
            end do
            n=latt_label(cc(1:Dimen))

            !Niu Nju
            if(i.eq.j) then
               ninj_local=Amat_local(i,i)
            else
               ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
            end if
            call kahan_sum_c(ninj_local,ninj_c(n),ninj_one(n))

            !Niu Njd
            ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
            call kahan_sum_c(ninj_local,ninj_c(n+Nsite),ninj_one(n+Nsite))
         end do
      end do
   end if
else if ((openbcx.eq.1).and.(openbcy.eq.1)) then
   if(dtype.EQ.'c') then
      i=1
      do j=1,Nsite,1
         !Niu Nju
         if(i.eq.j) then
            ninj_local=Amat_local(i,i)
         else
            ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         end if
         call kahan_sum_c(ninj_local,ninj_c(j),ninj_one(j))
         
         !Niu Njd
         ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)-Amat_local(i,j+Nsite)*Amat_local(j+Nsite,i)
         call kahan_sum_c(ninj_local,ninj_c(j+Nsite),ninj_one(j+Nsite))
      end do
   else if(dtype.EQ.'d') then
      i=1
      do j=1,Nsite,1
         !Niu Nju
         if(i.eq.j) then
            ninj_local=Amat_local(i,i)
         else
            ninj_local=Amat_local(i,i)*Amat_local(j,j)-Amat_local(i,j)*Amat_local(j,i)
         end if
         call kahan_sum_c(ninj_local,ninj_c(j),ninj_one(j))
         
         !Niu Njd
         ninj_local=Amat_local(i,i)*Amat_local(j+Nsite,j+Nsite)
         call kahan_sum_c(ninj_local,ninj_c(j+Nsite),ninj_one(j+Nsite))
      end do
   end if
end if

end subroutine add_ninj


subroutine add_edgec(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::edgecup_local(Nsite,Nsite)
complex(kind=8)::edgecdn_local(Nsite,Nsite)
integer::i,j,mi

edgecup_local(:,:)=0.d0
edgecdn_local(:,:)=0.d0

do j=1,Nsite,1
   do i=1,j,1
      if (abs(Hzero(i,j)).gt.1d-12) then         
         edgecup_local(i,j)=-Xi*(Amat_local(j,i)-Amat_local(i,j))
         !if(rank.eq.0) then
         !   write(*,*), Amat_local(i,j), Amat_local(j,i)
         !endif
         !edgecup_local(j,i)=conjg(edgecup(i,j))
         edgecdn_local(i,j)=-Xi*(Amat_local(j+Nsite,i+Nsite)-Amat_local(i+Nsite,j+Nsite))
         !edgecdn_local(j,i)=conjg(edgecdn(i,j))
         call kahan_sum_c(edgecup_local(i,j),edgecup_c(i,j),edgecup_one(i,j))
         call kahan_sum_c(edgecdn_local(i,j),edgecdn_c(i,j),edgecdn_one(i,j))
      end if
   end do
end do
end subroutine add_edgec

subroutine add_ni(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::niup_local(Nsite)
complex(kind=8)::nidn_local(Nsite)
integer::i

   do i=1,Nsite,1
      niup_local(i)=Amat_local(i,i)
      nidn_local(i)=Amat_local(i+Nsite,i+Nsite)
      call kahan_sum_c(niup_local(i),niup_c(i),niup_one(i))
      call kahan_sum_c(nidn_local(i),nidn_c(i),nidn_one(i))
   end do

 end subroutine add_ni


subroutine add_ck(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::ck_local
integer::n

do n=1,Nsite,1
   !total n(k)
   ck_local=Amat_local(n,n)+Amat_local(n+Nsite,n+Nsite)
   call kahan_sum_c(ck_local,ck_c(n),ck_one(n))
end do
end subroutine add_ck



subroutine add_sk(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::sk_local
integer::n

do n=1,Nsite,1
   !total S^x(k)
   sk_local=(Amat_local(n,n+Nsite)+Amat_local(n+Nsite,n))/2.d0
   call kahan_sum_c(sk_local,skx_c(n),skx_one(n))
   !total S^y(k)
   sk_local=(Amat_local(n,n+Nsite)-Amat_local(n+Nsite,n))/(2.d0*Xi)
   call kahan_sum_c(sk_local,sky_c(n),sky_one(n))
end do
end subroutine add_sk


subroutine add_pair_one_body(Amat_local)
use all_param
implicit none
complex(kind=8),intent(IN)::Amat_local(DNsite,DNsite)
complex(kind=8)::nconds_local,ncondt_local
integer::k,q,mk,mq
!add one-body density matrix
do k=1,DNsite,1
   do q=1,DNsite,1
      onebody(k,q)=onebody(k,q)+Amat_local(k,q)
   end do
end do

!build up the full matrix, only build half of the matrix
do k=1,Nsite,1
   call inverse_momentum(k,mk)
   do q=1,Nsite,1
      call inverse_momentum(q,mq)

      !triplet condensate fraction
      ncondt_local=Amat_local(k,q)*Amat_local(mk,mq)-Amat_local(k,mq)*Amat_local(mk,q)
      if (k.eq.q) then
         call kahan_sum_c(ncondt_local,ncondt_c,ncondt_one)
      endif

      !11
      pair_full(k,q)=pair_full(k,q)+ncondt_local!+Amat_local(k,q)*Amat_local(mk,mq)-Amat_local(k,mq)*Amat_local(mk,q)
      
      !12
      pair_full(k,q+Nsite)=pair_full(k,q+Nsite)+Amat_local(k,q+Nsite)*Amat_local(mk,mq+Nsite)- &
                        &  Amat_local(k,mq+Nsite)*Amat_local(mk,q+Nsite)
      !13
      pair_full(k,q+DNsite)=pair_full(k,q+DNsite)+0.5d0*(Amat_local(k,q)*Amat_local(mk,mq+Nsite)- &
       & Amat_local(k,mq+Nsite)*Amat_local(mk,q))-0.5d0*(Amat_local(k,q+Nsite)*Amat_local(mk,mq)- &
       & Amat_local(k,mq)*Amat_local(mk,q+Nsite))

      !22
      pair_full(k+Nsite,q+Nsite)=pair_full(k+Nsite,q+Nsite)+ &
     &  Amat_local(k+Nsite,q+Nsite)*Amat_local(mk+Nsite,mq+Nsite)-Amat_local(k+Nsite,mq+Nsite)*Amat_local(mk+Nsite,q+Nsite)
      !23
      pair_full(k+Nsite,q+DNsite)=pair_full(k+Nsite,q+DNsite)+0.5d0*(Amat_local(k+Nsite,q)*Amat_local(mk+Nsite,mq+Nsite)- &
       & Amat_local(k+Nsite,mq+Nsite)*Amat_local(mk+Nsite,q))-0.5d0*(Amat_local(k+Nsite,q+Nsite)*Amat_local(mk+Nsite,mq)- &
       & Amat_local(k+Nsite,mq)*Amat_local(mk+Nsite,q+Nsite))

      !singlet condensate fraction
      nconds_local=Amat_local(k,q)*Amat_local(mk+Nsite,mq+Nsite)-Amat_local(k,mq+Nsite)*Amat_local(mk+Nsite,q)
      if (k.eq.q) then
         call kahan_sum_c(nconds_local,nconds_c,nconds_one)
      endif

      !33
      pair_full(k+DNsite,q+DNsite)=pair_full(k+DNsite,q+DNsite)+ & 
      0.25d0*(nconds_local &
             +(Amat_local(k+Nsite,q+Nsite)*Amat_local(mk,mq)-Amat_local(k+Nsite,mq)*Amat_local(mk,q+Nsite)) &
             -(Amat_local(k,q+Nsite)*Amat_local(mk+Nsite,mq)-Amat_local(k,mq)*Amat_local(mk+Nsite,q+Nsite)) &
             -(Amat_local(k+Nsite,q)*Amat_local(mk,mq+Nsite)-Amat_local(k+Nsite,mq+Nsite)*Amat_local(mk,q)))
   end do
end do
end subroutine add_pair_one_body


!--------------------------------------------------
!This subroutine allocate the arrays in one measure
!--------------------------------------------------
subroutine allocate_one_meas()
use all_param
implicit none
if(dtype.ne.'w') then
   allocate(didj_one(Nsite))
   allocate(didj_c(Nsite))
else if(dtype.eq.'w') then
   allocate(didj_one(DNsite))
   allocate(didj_c(DNsite))
endif
allocate(dk_one(Nsite))
allocate(sisj_one(Nsite))
allocate(sisj_c(Nsite))
allocate(sk_one(Nsite))
if(dtype.ne.'w') then
   allocate(ninj_one(DNsite))
   allocate(ninj_c(DNsite))
else if(dtype.eq.'w') then
   allocate(ninj_one(2*DNsite))
   allocate(ninj_c(2*DNsite))
end if
allocate(edgecup_one(Nsite,Nsite))
allocate(edgecup_c(Nsite,Nsite))
allocate(edgecdn_one(Nsite,Nsite))
allocate(edgecdn_c(Nsite,Nsite))
allocate(niup_one(Nsite))
allocate(niup_c(Nsite))
allocate(nidn_one(Nsite))
allocate(nidn_c(Nsite))
allocate(nnk_one(DNsite))
allocate(skx_one(Nsite))
allocate(skx_c(Nsite))
allocate(sky_one(Nsite))
allocate(sky_c(Nsite))
allocate(ck_one(Nsite))
allocate(ck_c(Nsite))
allocate(pair_full(3*Nsite,3*Nsite))
allocate(onebody(DNsite,DNsite))
end subroutine allocate_one_meas


!----------------------------------------------------
!This subroutine deallocate the arrays in one measure
!----------------------------------------------------
subroutine deallocate_one_meas()
use all_param
implicit none
if(allocated(didj_one)) deallocate(didj_one)
if(allocated(didj_c)) deallocate(didj_c)
if(allocated(dk_one)) deallocate(dk_one)
if(allocated(sisj_one)) deallocate(sisj_one)
if(allocated(sisj_c)) deallocate(sisj_c)
if(allocated(sk_one)) deallocate(sk_one)
if(allocated(ninj_one)) deallocate(ninj_one)
if(allocated(ninj_c)) deallocate(ninj_c)
if(allocated(edgecup_c)) deallocate(edgecup_c)
if(allocated(edgecup_one)) deallocate(edgecup_one)
if(allocated(edgecdn_c)) deallocate(edgecdn_c)
if(allocated(edgecdn_one)) deallocate(edgecdn_one)
if(allocated(niup_one)) deallocate(niup_one)
if(allocated(niup_c)) deallocate(niup_c)
if(allocated(nidn_one)) deallocate(nidn_one)
if(allocated(nidn_c)) deallocate(nidn_c)
if(allocated(nnk_one)) deallocate(nnk_one)
if(allocated(skx_one)) deallocate(skx_one)
if(allocated(skx_c)) deallocate(skx_c)
if(allocated(sky_one)) deallocate(sky_one)
if(allocated(sky_c)) deallocate(sky_c)
if(allocated(ck_one)) deallocate(ck_one)
if(allocated(ck_c)) deallocate(ck_c)
if(allocated(pair_full)) deallocate(pair_full)
if(allocated(onebody)) deallocate(onebody)
end subroutine deallocate_one_meas