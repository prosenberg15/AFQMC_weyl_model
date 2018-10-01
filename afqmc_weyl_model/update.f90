subroutine run_qmc()
use all_param
implicit none
integer::i_qmc


Nupdate=0;Naccept=0

!energy and commute quantity
i_energy=0
E_one=zero;E_c=zero

!other obersverables
if(thermstep.LE.Nlen/2) then
  i_observ=0
  K_one=zero;K_c=zero
  V_one=zero;V_c=zero
  Nup_one=zero;Nup_c=zero
  Ndn_one=zero;Ndn_c=zero
  nconds_one=zero;nconds_c=zero
  ncondt_one=zero;ncondt_c=zero
  call allocate_one_meas()
  niup_one=zero;niup_c=zero
  nidn_one=zero;nidn_c=zero
  if(dtype.eq.'w') then
     d_one=zero;d_c=zero
     ninj_true_site_one=zero
     ninj_true_site_c=zero
     didj_true_site_one=zero
     didj_true_site_c=zero
  end if
  didj_one=zero;didj_c=zero
  sisj_one=zero;sisj_c=zero
  ninj_one=zero;ninj_c=zero
  edgecup_one=zero;edgecup_c=zero
  edgecdn_one=zero;edgecdn_c=zero
  skx_one=zero;skx_c=zero
  sky_one=zero;sky_c=zero
  ck_one=zero;ck_c=zero
  pair_full=zero;onebody=zero
end if

if(rank.eq.0) write(*,*) "Start therm:"
do i_qmc=1,Ntherm,1
   if(rank.eq.0) write(*,*) i_qmc
   call one_step()
end do
if(rank.eq.0) write(*,*) "Therm done!"


if(rank.eq.0) write(*,*) "Start measure:"
do i_qmc=1,Nmeas,1
   if(rank.eq.0) write(*,*) i_qmc
   call one_step_meas()
end do

end subroutine run_qmc


!only update one step
subroutine one_step()
use all_param
implicit none
integer::i_GS
integer::i,i_nblk,i_blk

i_GS=0
do i_nblk=Nblk,1,-1
   call set_blk_left(i_nblk)
   do i_blk=blk,1,-1
       i=(i_nblk-1)*blk+i_blk
       call v_k_update_right(i,i_nblk,i_blk,i_GS)
   end do
   if(i_nblk.ne.1) call copy_wf_dc(phi(1,1,2),phi_Nblk(1,1,i_nblk))
   if(i_GS.NE.0) write(*,*) "Warning!! phi_Nblk is not saved before MGS."
end do

i_GS=0
do i_nblk=1,Nblk,1
   call set_blk_right(i_nblk)
   do i_blk=1,blk,1
       i=(i_nblk-1)*blk+i_blk
       call k_v_update_left(i,i_nblk,i_blk,i_GS)
   end do
   if(i_nblk.ne.Nblk) call copy_wf_dc(phi(1,1,1),phi_Nblk(1,1,i_nblk+1))
   if(i_GS.NE.0) write(*,*) "Warning!! phi_Nblk is not saved before MGS."
end do
!set the right phi(1,1,2)
 call copy_wf_dc(phi_Nblk(1,1,Nblk+1),phi(1,1,2))
end subroutine one_step


!update one step and measure
subroutine one_step_meas()
use all_param
implicit none
integer::i_GS
complex(kind=8)::templ(DNsite,Ntot),tempr(DNsite,Ntot)
integer::i,i_nblk,i_blk

i_GS=0
do i_nblk=Nblk,1,-1
   call set_blk_left(i_nblk)
   do i_blk=blk,1,-1
       i=(i_nblk-1)*blk+i_blk
       call v_k_update_right(i,i_nblk,i_blk,i_GS)

       !-------
       !measure
       !-------
       if(mod(i-1,meastep).EQ.0.or.i.eq.1) then
          call prepare_meas_right(templ,tempr)
          if((Nlen-i+1).GE.thermstep.and.(i-1).GE.thermstep) then
            call meas_observ(templ,tempr)
          else
            call meas_energy(templ,tempr)
          end if
       end if


   end do
   if(i_nblk.ne.1) call copy_wf_dc(phi(1,1,2),phi_Nblk(1,1,i_nblk))
   if(i_GS.NE.0) write(*,*) "Warning!! phi_Nblk is not saved before MGS."
end do


i_GS=0
do i_nblk=1,Nblk,1
   call set_blk_right(i_nblk)
   do i_blk=1,blk,1
       i=(i_nblk-1)*blk+i_blk
       call k_v_update_left(i,i_nblk,i_blk,i_GS)

       !-------
       !measure
       !-------
       if(mod(i,meastep).EQ.0.or.i.eq.Nlen) then
         call prepare_meas_left(templ,tempr)
         if(i.GE.thermstep.and.(Nlen-i).GE.thermstep) then
           call meas_observ(templ,tempr)
         else
           call meas_energy(templ,tempr)
         end if
       end if

   end do
   if(i_nblk.ne.Nblk) call copy_wf_dc(phi(1,1,1),phi_Nblk(1,1,i_nblk+1))
   if(i_GS.NE.0) write(*,*) "Warning!! phi_Nblk is not saved before MGS."
end do

!set the right phi(1,1,2)
call copy_wf_dc(phi_Nblk(1,1,Nblk+1),phi(1,1,2))

end subroutine one_step_meas


!prepare before measure
subroutine prepare_meas_right(templ,tempr)
use all_param
implicit none
complex(kind=8),intent(OUT)::templ(DNsite,Ntot),tempr(DNsite,Ntot)
!apply Dagger[exp(-dt*K/2)] to phi_l
!apply exp(dt*K/2) to phi_r
if(pfft.eq.0) then
   call k_to_ph_dc(exp_halfK,phi(1,1,2),templ(1,1))
   call k_to_ph_dc(exp_mhalfK,phi(1,1,1),tempr(1,1))
else if(pfft.eq.1) then
   call copy_wf_dc(phi(1,1,2),templ(1,1))
   call k_to_ph_dc_fftw(exp_he,templ(1,1),uk,vk)
   call copy_wf_dc(phi(1,1,1),tempr(1,1))
   call k_to_ph_dc_fftw(exp_mhe,tempr(1,1),uk,vk)
end if
end subroutine prepare_meas_right


!prepare before measure
subroutine prepare_meas_left(templ,tempr)
use all_param
implicit none
complex(kind=8),intent(OUT)::templ(DNsite,Ntot),tempr(DNsite,Ntot)
!apply Dagger[exp(-dt*K/2)] to phi_l
!apply exp(dt*K/2) to phi_r
if(pfft.eq.0) then
   call k_to_ph_dc(exp_mhalfK,phi(1,1,2),templ(1,1))
   call k_to_ph_dc(exp_mhalfK,phi(1,1,1),tempr(1,1))
else if(pfft.eq.1) then
   call copy_wf_dc(phi(1,1,2),templ(1,1))
   call k_to_ph_dc_fftw(exp_mhe,templ(1,1),uk,vk)
   call copy_wf_dc(phi(1,1,1),tempr(1,1))
   call k_to_ph_dc_fftw(exp_mhe,tempr(1,1),uk,vk)
end if
end subroutine prepare_meas_left



subroutine k_v_update_left(i,i_nblk,i_blk,i_GS)
use all_param
implicit none
integer,intent(IN)::i,i_nblk,i_blk
integer,intent(INOUT)::i_GS
complex(kind=8)::temp(DNsite,Ntot)
complex(kind=8)::explr(DNsite),explr_tmp(DNsite)
complex(kind=8)::impn,impo
real(kind=8)::alpha
complex(kind=8)::an,ao
real(kind=8)::anm,anm1,anm2
real(kind=8)::naux(Nsite)


!apply exp(-dt*V) to phi_r
 call diag_op(aux(1,i),explr)
 call d_to_phi_dc_n(explr,phi(1,1,1),temp)

!we need two aux fields for weyl model
 
!get the new phi_l
 call copy_wf_dc(phi_blk(1,1,i_blk),phi(1,1,2))


!get the new field 
 call cal_n_back(phi(1,1,2),phi(1,1,1))
 call sample_field(naux)
 call diag_op(naux,explr)
 call d_to_phi_dc(explr,phi(1,1,1))


!get the scale by a number
 call get_scale(ao,aux(1,i))
 call get_scale(an,naux)

!Get the importance function
 impo=imp_save
 call imp_fun_dc(phi(1,1,2),phi(1,1,1),impn)

!correct by metropolis
 alpha=dble((impn/impo)/(ao/an))
 if(alpha.LT.0.d0) then
   write(*,*) "sign problem!"
   write(*,*) "alpha:",alpha
   write(*,*) "impn:",impn
   write(*,*) "impo:",impo
   call mystop
 end if

 Nupdate=Nupdate+1
 if(rndm().LT.alpha) then
   Naccept=Naccept+1
   aux(1:Nsite,i)=naux(1:Nsite)
   imp_save=impn
 else
   call copy_wf_dc(temp,phi(1,1,1)) 
 end if


!apply exp(-dt*K) to phi_r
 if(pfft.eq.0) then
    call copy_wf_dc(phi(1,1,1),temp(1,1))
    call k_to_ph_dc(exp_K,temp(1,1),phi(1,1,1))
 else if(pfft.eq.1) then
    call k_to_ph_dc_fftw(exp_e,phi(1,1,1),uk,vk)
 end if



!add MGS
 i_GS=i_GS+1
 if(i_GS.EQ.StepforGram) then
    i_GS=0
    if(dtype.EQ.'c') then
       call modGS(phi(1,1,1),DNsite,Ntot,anm)
       imp_save=imp_save*anm
    else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
       call modGS(phi(1:Nsite,1:Nspin(1),1),Nsite,Nspin(1),anm1)
       call modGS(phi((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,1),Nsite,Nspin(2),anm2)
       imp_save=imp_save*anm1*anm2
    end if
 end if
 imp_save=imp_save/rphi(i_blk)
end subroutine k_v_update_left



subroutine v_k_update_right(i,i_nblk,i_blk,i_GS)
use all_param
implicit none
integer,intent(IN)::i,i_nblk,i_blk
integer,intent(INOUT)::i_GS
complex(kind=8)::temp(DNsite,Ntot)
complex(kind=8)::explr(DNsite),explr_tmp(DNsite)
complex(kind=8)::impn,impo
real(kind=8)::alpha
complex(kind=8)::an,ao
real(kind=8)::anm,anm1,anm2
real(kind=8)::naux(Nsite)


!apply Dagger[exp(-dt*K)] to phi_l
 if(pfft.eq.0) then
    call copy_wf_dc(phi(1,1,2),temp(1,1))
    call k_to_ph_dc(exp_K,temp(1,1),phi(1,1,2))
 else if(pfft.eq.1) then
    call k_to_ph_dc_fftw(exp_e,phi(1,1,2),uk,vk)
 end if

!backup the phi_l
 call copy_wf_dc(phi(1,1,2),temp(1,1))


!get the new phi_r
 call copy_wf_dc(phi_blk(1,1,i_blk),phi(1,1,1))


!get the new field 
 call cal_n_back(phi(1,1,2),phi(1,1,1))
 call sample_field(naux)
 call diag_op(naux,explr)
 explr_tmp(1:DNsite)=conjg(explr(1:DNsite))
 call d_to_phi_dc(explr_tmp,phi(1,1,2))


!get the scale by a number
 call get_scale(ao,aux(1,i))
 call get_scale(an,naux)


!Get the importance function
 impo=imp_save
 call imp_fun_dc(phi(1,1,2),phi(1,1,1),impn)


!correct by metropolis
 alpha=dble((impn/impo)/(ao/an))
 if(alpha.LT.0.d0) then
   write(*,*) "sign problem!"
   write(*,*) "alpha:",alpha
   write(*,*) "impn:",impn
   write(*,*) "impo:",impo
   call mystop
 end if

 Nupdate=Nupdate+1
 if(rndm().LT.alpha) then
   Naccept=Naccept+1
   aux(1:Nsite,i)=naux(1:Nsite)
   imp_save=impn
 else
   call diag_op(aux(1,i),explr)
   explr_tmp(1:DNsite)=conjg(explr(1:DNsite))
   call d_to_phi_dc_n(explr_tmp,temp,phi(1,1,2))
 end if

!add MGS
 i_GS=i_GS+1
 if(i_GS.EQ.StepforGram) then
   i_GS=0
   if(dtype.EQ.'c') then
     call modGS(phi(1,1,2),DNsite,Ntot,anm)
     imp_save=imp_save*anm
   else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
     call modGS(phi(1:Nsite,1:Nspin(1),2),Nsite,Nspin(1),anm1)
     call modGS(phi((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,2),Nsite,Nspin(2),anm2)
     imp_save=imp_save*anm1*anm2
   end if
 end if
 imp_save=imp_save/rphi(i_blk)
end subroutine v_k_update_right


!For propagation to the left
subroutine set_blk_left(i_nblk)
use all_param
implicit none
integer,intent(IN)::i_nblk
complex(kind=8)::temp(DNsite,Ntot)
complex(kind=8)::explr(DNsite),explr_tmp(DNsite)
real(kind=8)::anm,anm1,anm2,tmp_scale
integer::mi,i,i_GS

rphi=one

call copy_wf_dc(phi_Nblk(1,1,i_nblk),phi_blk(1,1,1))

i_GS=0
do mi=2,blk,1
   !apply exp(-dt*V) to phi_r
   i=mi-1+(i_nblk-1)*blk
   call diag_op(aux(1,i),explr)
   call d_to_phi_dc_n(explr,phi_blk(1,1,mi-1),phi_blk(1,1,mi))

   !apply exp(-dt*K) to phi_r
   if(pfft.eq.0) then
      call copy_wf_dc(phi_blk(1,1,mi),temp(1,1))
      call k_to_ph_dc(exp_K,temp(1,1),phi_blk(1,1,mi))
   else if(pfft.eq.1) then
      call k_to_ph_dc_fftw(exp_e,phi_blk(1,1,mi),uk,vk)
   end if

   i_GS=i_GS+1
   if(i_GS.EQ.StepforGram) then
     i_GS=0
     if(dtype.EQ.'c') then
       call modGS(phi_blk(1:DNsite,1:Ntot,mi),DNsite,Ntot,anm)
       rphi(mi)=anm
     else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
       call modGS(phi_blk(1:Nsite,1:Nspin(1),mi),Nsite,Nspin(1),anm1)
       call modGS(phi_blk((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,mi),Nsite,Nspin(2),anm2)
       rphi(mi)=anm1*anm2
     end if
   end if

end do


!set the phir and imp_save
i=i+1
call diag_op(aux(1,i),explr)
call d_to_phi_dc_n(explr,phi_blk(1,1,blk),phi(1,1,1))
!apply exp(-dt*K) to phi_r
if(pfft.eq.0) then
   call copy_wf_dc(phi(1,1,1),temp(1,1))
   call k_to_ph_dc(exp_K,temp(1,1),phi(1,1,1))
else if(pfft.eq.1) then
   call k_to_ph_dc_fftw(exp_e,phi(1,1,1),uk,vk)
end if
!add MGS
i_GS=i_GS+1
tmp_scale=1.d0
if(i_GS.EQ.StepforGram) then
  i_GS=0
  if(dtype.EQ.'c') then
    call modGS(phi(1,1,1),DNsite,Ntot,anm)
    tmp_scale=tmp_scale*anm
  else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
    call modGS(phi(1:Nsite,1:Nspin(1),1),Nsite,Nspin(1),anm1)
    call modGS(phi((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,1),Nsite,Nspin(2),anm2)
    tmp_scale=tmp_scale*anm1*anm2
  end if
end if
call imp_fun_dc(phi(1,1,2),phi(1,1,1),imp_save)
imp_save=imp_save/tmp_scale
end subroutine set_blk_left


!For propagation to the right
subroutine set_blk_right(i_nblk)
use all_param
implicit none
integer,intent(IN)::i_nblk
complex(kind=8)::temp(DNsite,Ntot)
complex(kind=8)::explr(DNsite),explr_tmp(DNsite)
real(kind=8)::anm,anm1,anm2
integer::mi,i,i_GS

rphi=one

!apply exp(-dt*K) to phi_l
if(pfft.eq.0) then
   call k_to_ph_dc(exp_K,phi_Nblk(1,1,i_nblk+1),phi_blk(1,1,blk))
else if(pfft.eq.1) then
   call copy_wf_dc(phi_Nblk(1,1,i_nblk+1),phi_blk(1,1,blk))
   call k_to_ph_dc_fftw(exp_e,phi_blk(1,1,blk),uk,vk)
end if

i_GS=0
do mi=blk-1,1,-1

   !apply exp(-dt*V) to phi_l
   i=(i_nblk-1)*blk+mi+1
   call diag_op(aux(1,i),explr)
   explr_tmp(1:DNsite)=conjg(explr(1:DNsite))
   call d_to_phi_dc_n(explr_tmp,phi_blk(1,1,mi+1),phi_blk(1,1,mi))

   !apply exp(-dt*K) to phi_l
   if(pfft.eq.0) then
      call copy_wf_dc(phi_blk(1,1,mi),temp(1,1))
      call k_to_ph_dc(exp_K,temp(1,1),phi_blk(1,1,mi))
   else if(pfft.eq.1) then
      call k_to_ph_dc_fftw(exp_e,phi_blk(1,1,mi),uk,vk)
   end if

   i_GS=i_GS+1
   if(i_GS.EQ.StepforGram) then
     i_GS=0
     if(dtype.EQ.'c') then
       call modGS(phi_blk(1:DNsite,1:Ntot,mi),DNsite,Ntot,anm)
       rphi(mi)=anm
     else if((dtype.EQ.'d').or.(dtype.EQ.'w')) then
       call modGS(phi_blk(1:Nsite,1:Nspin(1),mi),Nsite,Nspin(1),anm1)
       call modGS(phi_blk((Nsite+1):(DNsite),(Nspin(1)+1):Ntot,mi),Nsite,Nspin(2),anm2)
       rphi(mi)=anm1*anm2
     end if
   end if
end do

!set the imp_save
i=i-1
call diag_op(aux(1,i),explr)
explr_tmp(1:DNsite)=conjg(explr(1:DNsite))
call d_to_phi_dc_n(explr_tmp,phi_blk(1,1,1),temp(1,1))
call imp_fun_dc(temp(1,1),phi(1,1,1),imp_save)
end subroutine set_blk_right


subroutine sample_field(naux)
use all_param
implicit none
real(kind=8),intent(OUT)::naux(Nsite)
real(kind=8)::x,p(2),Normp
integer::j

if(kcrn.eq.1) then !Use discrete spin decouple
  do j=1,Nsite,1

     x=rndm()
     p(1)=abs(exp(-gamaf*ng(j)))
     p(2)=abs(exp(gamaf*ng(j)))
     Normp=p(1)+p(2)

     if(x.LT.p(1)/Normp) then
       naux(j)=-1.d0
     else
       naux(j)=1.d0
     end if

  end do
else if(kcrn.eq.2) then !Use discrete charge decouple
  do j=1,Nsite,1

     x=rndm()
     p(1)=abs(exp((1-ng(j))*gamaf))
     p(2)=abs(exp((ng(j)-1)*gamaf))
     Normp=p(1)+p(2)

     if(x.LT.p(1)/Normp) then
       naux(j)=-1.d0
     else
       naux(j)=1.d0
     end if

  end do
else
  write(*,*) "Something is wrong with kcrn input:",kcrn
  call mystop
end if
end subroutine sample_field


subroutine diag_op(naux,explr)
use all_param
implicit none
real(kind=8),intent(IN)::naux(Nsite)
complex(kind=8),intent(OUT)::explr(DNsite)
integer::j
if(kcrn.eq.1) then !Use discrete spin decouple
   do j=1,Nsite,1
      explr(j)=exp(-1.d0*(dt*onsitU*0.5d0-naux(j)*gamaf))
      explr(j+Nsite)=exp(-1.d0*(dt*onsitU*0.5d0+naux(j)*gamaf))
   end do
else if(kcrn.eq.2) then !Use discrete charge decouple
   do j=1,Nsite,1
      explr(j)=exp(-1.d0*(dt*onsitU*0.5d0-naux(j)*gamaf))
      explr(j+Nsite)=exp(-1.d0*(dt*onsitU*0.5d0-naux(j)*gamaf))
   end do
else
   write(*,*) "Something is wrong with kcrn input:",kcrn
   call mystop
end if
end subroutine diag_op 


subroutine get_scale(alpha,naux)
use all_param
implicit none
complex(kind=8),intent(OUT)::alpha
real(kind=8),intent(IN)::naux(Nsite)
complex(kind=8)::en
integer::j

if(kcrn.eq.1) then !Use discrete spin decouple
  en=one
  do j=1,Nsite,1
     en=en/abs(exp(gamaf*naux(j)*ng(j)))
  end do
  alpha=en
else if(kcrn.eq.2) then !Use discrete charge decouple
  en=one
  do j=1,Nsite,1
     en=en*exp(-gamaf*naux(j))/abs(exp(gamaf*naux(j)*(ng(j)-1.d0)))
  end do
  alpha=en
else
  write(*,*) "Something is wrong with kcrn input:",kcrn
  call mystop
end if
end subroutine get_scale


subroutine cal_n_back(phl,phr)
use all_param
implicit none
complex(kind=8),intent(IN)::phl(DNsite,Ntot),phr(DNsite,Ntot)
complex(kind=8)::ovp_local(Ntot,Ntot)
complex(kind=8)::phr_tmp(DNsite,Ntot)
complex(kind=8)::imp_local
complex(kind=8)::avgn(DNsite)
integer::j
!set the background
if(bgset.eq.0) then !use mean field
  return
else if(bgset.eq.1) then !use dynamic field from walker
  call over_lap_dc(phl(1,1),phr(1,1),ovp_local(1,1))
  call copy_wf_dc(phr(1:DNsite,1:Ntot),phr_tmp(1:DNsite,1:Ntot))
  call linear_equation_d_dc(ovp_local(1:Ntot,1:Ntot),imp_local,phr_tmp(1:DNsite,1:Ntot))
  call cal_Amat_diag_dc(phl(1,1),phr_tmp(1,1),avgn(1))
  do j=1,Nsite,1
     if(kcrn.eq.3.or.kcrn.eq.1) then
        ng(j)=avgn(j)-avgn(j+Nsite)
     else if(kcrn.eq.4.or.kcrn.eq.2) then
        ng(j)=avgn(j)+avgn(j+Nsite)
     end if
  end do
  !write(*,*) ng;pause
else
   write(*,*) "something is wrong with bgset:",bgset
   call mystop
end if
end subroutine cal_n_back
