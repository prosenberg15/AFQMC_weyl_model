subroutine data_mani()
use all_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer(kind=8)::num_tmp,Ns
complex(kind=8)::temp_array(Nsize)
real(kind=8)::rtemp_array(Nsize)
complex(kind=8)::mean_c,mean_cup,mean_cdn
real(kind=8)::mean,error,errorup,errordn
real(kind=8)::kx,ky
real(kind=8)::m_n,m_n_e,m_np,m_np_e,m_nm,m_nm_e,m_sx,m_sx_e,m_sy,m_sy_e
integer::i,j

call get_filename()

!print the accept ratio:
#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(Nupdate,num_tmp,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
 Nupdate=num_tmp
 call MPI_ALLREDUCE(Naccept,num_tmp,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
 Naccept=num_tmp
#endif
 if(rank.eq.0) then
   write(*,*) "Total update number is:",Nupdate
   write(*,*) "Total accept number is:",Naccept
   write(*,*) "The accept ratio is:",dble(Naccept)/dble(Nupdate)
   write(*,*) "Measure for energy in each thread:",i_energy
   write(*,*) "Measure for observable in each thread:",i_observ
 end if

 E_one=E_one/dble(i_energy)

#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(E_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 rtemp_array(:)=dble(temp_array(:)) 
 call err_anal(rtemp_array,Nsize,mean,error)
#else
 mean=E_one
 error=0.d0
#endif
 if(rank.eq.0) then
   write(*,*) "The energy is:",mean,error
 end if


!print the observables
if(thermstep.LE.Nlen/2) then

 !write edge current
 if(rank.eq.0) call openUnit(edgecName,16,'R')
 do i=1,Nsite,1
    do j=1,Nsite,1
       edgecup_one(i,j)=edgecup_one(i,j)/dble(i_observ)
       edgecdn_one(i,j)=edgecdn_one(i,j)/dble(i_observ)
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call MPI_GATHER(edgecup_one(i,j),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
    call err_anal_c(temp_array(1),Nsize,mean_cup,errorup)
    call MPI_GATHER(edgecdn_one(i,j),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
    call err_anal_c(temp_array(1),Nsize,mean_cdn,errordn)
#else
    mean_cup=edgecup_one(i,j)
    errorup=0.d0
    mean_cdn=edgecdn_one(i,j)
    errordn=0.d0
#endif
    if(rank.eq.0) write(16,'(2I4,6E26.16)') i,j,dble(mean_cup),dimag(mean_cup),errorup,dble(mean_cdn),dimag(mean_cdn),errordn
 end do
end do
 if(rank.eq.0) close(16)

 K_one=K_one/dble(i_observ)
#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(K_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 rtemp_array(:)=dble(temp_array(:))
 call err_anal(rtemp_array,Nsize,mean,error)
#else
 mean=K_one
 error=0.d0
#endif
 if(rank.eq.0) then
   write(*,*) "The kinetic energy is:",mean,error
 end if



 V_one=V_one/dble(i_observ)
#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(V_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 rtemp_array(:)=dble(temp_array(:))
 call err_anal(rtemp_array,Nsize,mean,error)
#else
 mean=V_one
 error=0.d0
#endif
 if(rank.eq.0) then
   write(*,*) "The potential energy is:",mean,error
 end if


 Nup_one=Nup_one/dble(i_observ)
#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(Nup_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 rtemp_array(:)=dble(temp_array(:))
 call err_anal(rtemp_array,Nsize,mean,error)
#else
 mean=Nup_one
 error=0.d0
#endif
 if(rank.eq.0) then
   write(*,*) "Nup is:",mean,error
 end if


 Ndn_one=Ndn_one/dble(i_observ)
#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_GATHER(Ndn_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 rtemp_array(:)=dble(temp_array(:))
 call err_anal(rtemp_array,Nsize,mean,error)
#else
 mean=Ndn_one
 error=0.d0
#endif
 if(rank.eq.0) then
   write(*,*) "Ndn is:",mean,error
 end if

    !write pairing correlation didj
 if(rank.eq.0) call openUnit(didjName,16,'R')
 if(dtype.ne.'w') then
    Ns=Nsite
 else if(dtype.eq.'w') then
    Ns=DNsite
 end if
    do i=1,Ns,1
       if((openbcx.eq.0).and.(openbcy.eq.0)) then
          if(dtype.ne.'w') then
             didj_one(i)=didj_one(i)/dble(i_observ*Nsite)
          else if (dtype.eq.'w') then
             didj_one(i)=didj_one(i)/dble(i_observ*Nbravais)
          end if
       else if((openbcx.eq.1).or.(openbcy.eq.1)) then
          didj_one(i)=didj_one(i)/dble(i_observ)
       end if
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(didj_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=didj_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)

    if((openbcx.eq.0).and.(openbcy.eq.0)) then
    !write pairing stucture factor
    if(rank.eq.0) call openUnit(dkName,16,'R')
    aforwin(1:Nsite)=didj_one(1:Nsite)
    call dfftw_execute(planf,aforwin,aforwout)
    dk_one(1:Nsite)=aforwout(1:Nsite)/dble(Nsite)
    do i=1,Nsite,1
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(dk_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=dk_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)

 end if
    
    !write spin spin correlation
    if(rank.eq.0) call openUnit(ScorrName,16,'R')
    do i=1,Nsite,1
       if((openbcx.eq.0).and.(openbcy.eq.0)) then
          sisj_one(i)=sisj_one(i)/dble(i_observ*Nsite)
       else if((openbcx.eq.1).or.(openbcy.eq.1)) then
          sisj_one(i)=sisj_one(i)/dble(i_observ)
       end if
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(sisj_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=sisj_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)

 if((openbcx.eq.0).and.(openbcy.eq.0)) then
    !write spin stucture factor
    if(rank.eq.0) call openUnit(SkName,16,'R')
    aforwin(1:Nsite)=sisj_one(1:Nsite)
    call dfftw_execute(planf,aforwin,aforwout)
    sk_one(1:Nsite)=aforwout(1:Nsite)/dble(Nsite)
    do i=1,Nsite,1
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(sk_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=sk_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)
    
 end if

 !write density
    if(rank.eq.0) call openUnit(niName,16,'R')
    do i=1,Nsite,1
       if((openbcx.eq.0).and.(openbcy.eq.0)) then
          niup_one(i)=niup_one(i)/dble(i_observ*Nsite)
          nidn_one(i)=nidn_one(i)/dble(i_observ*Nsite)
       else if((openbcx.eq.1).or.(openbcy.eq.1)) then
          niup_one(i)=niup_one(i)/dble(i_observ)
          nidn_one(i)=nidn_one(i)/dble(i_observ)
       end if
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(niup_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_cup,errorup)
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(nidn_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_cdn,errordn)
#else
       mean_cup=niup_one(i)
       mean_cdn=nidn_one(i)
       errorup=0.d0
       errordn=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,6E26.16)') i,dble(mean_cup),dimag(mean_cup),errorup,dble(mean_cdn),dimag(mean_cdn),errordn
    end do
    if(rank.eq.0) close(16)

 !write charge charge correlation
    if(rank.eq.0) call openUnit(ninjName,16,'R')
    if(dtype.ne.'w') then
       Ns=DNsite
    else if(dtype.eq.'w') then
       Ns=2*DNsite
    end if
    do i=1,Ns,1
       if((openbcx.eq.0).and.(openbcy.eq.0)) then
          if(dtype.ne.'w') then
             ninj_one(i)=ninj_one(i)/dble(i_observ*Nsite)
          else if(dtype.eq.'w') then
             ninj_one(i)=ninj_one(i)/dble(i_observ*Nbravais)
          end if
       else if((openbcx.eq.1).or.(openbcy.eq.1)) then
          ninj_one(i)=ninj_one(i)/dble(i_observ)
       end if
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(ninj_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=ninj_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)

 if((openbcx.eq.0).and.(openbcy.eq.0)) then
    !write charge stucture factor
    if(rank.eq.0) call openUnit(nnkName,16,'R')
    aforwin(1:Nsite)=ninj_one(1:Nsite)
    call dfftw_execute(planf,aforwin,aforwout)
    nnk_one(1:Nsite)=aforwout(1:Nsite)/dble(Nsite)
    
    aforwin(1:Nsite)=ninj_one(1+Nsite:DNsite)
    call dfftw_execute(planf,aforwin,aforwout)
    nnk_one(1+Nsite:DNsite)=aforwout(1:Nsite)/dble(Nsite)
    do i=1,DNsite,1
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(nnk_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=nnk_one(i)
       error=0.d0
#endif
       if(rank.eq.0) write(16,'(1I4,3E26.16)') i,dble(mean_c),dimag(mean_c),error
    end do
    if(rank.eq.0) close(16)

    
    !write momentum distribution
    if(rank.eq.0) call openUnit(cksName,16,'R')
    do i=1,Nsite,1
       ck_one(i)=ck_one(i)/dble(i_observ)
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(ck_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=ck_one(i)
       error=0.d0
#endif
       m_n=dble(mean_c);m_n_e=error

       skx_one(i)=skx_one(i)/dble(i_observ)
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(skx_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=skx_one(i)
       error=0.d0
#endif
       m_sx=dble(mean_c);m_sx_e=error

       sky_one(i)=sky_one(i)/dble(i_observ)
#ifdef MPI
       call MPI_BARRIER(MPI_COMM_WORLD,IERR)
       call MPI_GATHER(sky_one(i),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
       call err_anal_c(temp_array(1),Nsize,mean_c,error)
#else
       mean_c=sky_one(i)
       error=0.d0
#endif
       m_sy=dble(mean_c);m_sy_e=error

       m_np=m_n/2.d0+sqrt(m_sx**2+m_sy**2)
       m_np_e=sqrt(m_n_e**2/2.d0+m_sx_e**2+m_sy_e**2)
       
       m_nm=m_n/2.d0-sqrt(m_sx**2+m_sy**2)
       m_nm_e=sqrt(m_n_e**2/2.d0+m_sx_e**2+m_sy_e**2)
       
       if(rank.eq.0) write(16,'(1I4,10E26.16)') i,m_n,m_n_e,m_np,m_np_e,m_nm,m_nm_e,m_sx,m_sx_e,m_sy,m_sy_e
    end do
    if(rank.eq.0) close(16)


    ! condensate fraction
    nconds_one=nconds_one/dble(i_observ)
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call MPI_GATHER(nconds_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
    call err_anal_c(temp_array,Nsize,mean_c,error)
#else
    mean_c=nconds_one
    error=0.d0
#endif
    if(rank.eq.0) then
       write(*,*) "singlet condensate:",mean_c,error
    end if

    ncondt_one=ncondt_one/dble(i_observ)
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call MPI_GATHER(ncondt_one,1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
    call err_anal_c(temp_array,Nsize,mean_c,error)
#else
    mean_c=ncondt_one
    error=0.d0
#endif
    if(rank.eq.0) then
       write(*,*) "singlet condensate:",mean_c,error
    end if

    !deal with pairing matrix
    do i=1,3*Nsite,1
       do j=i,3*Nsite,1 ! j count from i, only use upper triangular matrix
#ifdef MPI
          pair_full(i,j)=pair_full(i,j)/dble(i_observ)
          call MPI_BARRIER(MPI_COMM_WORLD,IERR)
          call MPI_GATHER(pair_full(i,j),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
          call err_anal_c(temp_array(1),Nsize,pair_full(i,j),error)
#else
          pair_full(i,j)=pair_full(i,j)/dble(i_observ)
#endif
       end do
    end do

    do i=1,DNsite,1
       do j=1,DNsite,1
#ifdef MPI
          onebody(i,j)=onebody(i,j)/dble(i_observ)
          call MPI_BARRIER(MPI_COMM_WORLD,IERR)
          call MPI_GATHER(onebody(i,j),1,MPI_DOUBLE_COMPLEX,temp_array(1),1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
          call err_anal_c(temp_array(1),Nsize,onebody(i,j),error)
#else
          onebody(i,j)=onebody(i,j)/dble(i_observ)
#endif
       end do
    end do

    if(rank.eq.0) then
       call minus_pairing_bg()
       call openUnit(pairmName,16,'R')
       do i=1,3*Nsite,1
          do j=1,3*Nsite,1
             write(16,'(2I4,2E26.16)') i,j,dble(pair_full(i,j)),dimag(pair_full(i,j))
          end do
       end do
       close(16)
       call diag_pair(pair_full,3*Nsite)
    end if


 end if
end if

!Write the useful label
if(rank.eq.0) then
  open(unit=16,file='lattice-label.dat',status='replace')
      do i=1,Nsite,1
         write(16,*) i,coor(i,1:Dimen)
      end do
  close(16)

  open(unit=16,file='momentum-label.dat',status='replace')
      do i=1,Nsite,1
         kx=(dble(coor(i,1)-1)+kbound(1))*2*Pi/Nl(1)
         ky=(dble(coor(i,2)-1)+kbound(2))*2*Pi/Nl(2)
         if(kx.GT.Pi) kx=kx-2.d0*Pi
         if(ky.GT.Pi) ky=ky-2.d0*Pi
         write(16,'(1I4,2E26.16)') i,kx,ky
      end do
  close(16)
end if
end subroutine data_mani



!sumk=sumk+input, avoid the low digits error
subroutine kahan_sum_c(input,c,sumk)
implicit none
complex(kind=8),intent(IN)::input
complex(kind=8),intent(INOUT)::c,sumk
complex(kind=8)::y,t

y=input-c
t=sumk+y
c=(t-sumk)-y
sumk=t
end subroutine kahan_sum_c


!---------------------------
!Get the error bar of dat(N)
!---------------------------
subroutine err_anal(dat,N,m,er)
implicit none
integer,intent(IN)::N
real(kind=8),intent(IN)::dat(N)
real(kind=8),intent(OUT)::m,er
real(kind=8)::erray(N)
integer::i

if(N.LE.0) then
  write(*,*) "N should not be smaller than or EQ 0", N
  !call mystop
  stop
else if(N.EQ.1) then
  !write(*,*) "N eq 1 warning", N
  m=dat(1)
  er=0.d0
  return
end if

!Get the mean
call kahan_sum_r_array(N,dat,m)
m=m/dble(N)

do i=1,N,1
   erray(i)=(dat(i)-m)**2
end do
call kahan_sum_r_array(N,erray,er)
er=er/dble(N)
er=sqrt(er/dble(N-1))
end subroutine err_anal


!sum a real array
subroutine kahan_sum_r_array(N,input,sumk)
implicit none
integer,intent(IN)::N
real(kind=8),intent(IN)::input(N)
real(kind=8),intent(OUT)::sumk
real(kind=8)::y,c,t
integer(kind=8)::i

sumk=0.d0;c=0.d0
do i=1,N,1
  y=input(i)-c
  t=sumk+y
  c=(t-sumk)-y
  sumk=t
end do
end subroutine kahan_sum_r_array


!-----------------------------------
!Get the error bar of dat(N) complex
!-----------------------------------
subroutine err_anal_c(dat,N,m,er)
implicit none
integer,intent(IN)::N
complex(kind=8),intent(IN)::dat(N)
complex(kind=8),intent(OUT)::m
real(kind=8),intent(OUT)::er
real(kind=8)::erray(N)
integer::i

if(N.LE.0) then
  write(*,*) "N should not be smaller than or EQ 0", N
  !call mystop
  stop
else if(N.EQ.1) then
  !write(*,*) "N eq 1 warning", N
  m=dat(1)
  er=0.d0
  return
end if

!Get the mean
call kahan_sum_c_array(N,dat,m)
m=m/dble(N)

do i=1,N,1
   erray(i)=abs(dat(i)-m)**2
end do
call kahan_sum_r_array(N,erray,er)
er=er/dble(N)
er=sqrt(er/dble(N-1))
end subroutine err_anal_c


!sum a complex array
subroutine kahan_sum_c_array(N,input,sumk)
implicit none
integer,intent(IN)::N
complex(kind=8),intent(IN)::input(N)
complex(kind=8),intent(OUT)::sumk
complex(kind=8)::y,c,t
integer(kind=8)::i

sumk=dcmplx(0.d0,0.d0);c=dcmplx(0.d0,0.d0)
do i=1,N,1
  y=input(i)-c
  t=sumk+y
  c=(t-sumk)-y
  sumk=t
end do
end subroutine kahan_sum_c_array

subroutine minus_pairing_bg()
use all_param
implicit none
integer::k,q,mk,mq
do k=1,Nsite,1
   call inverse_momentum(k,mk)
   do q=1,Nsite,1
      call inverse_momentum(q,mq)

      if(q.GE.k) then
        !11
         pair_full(k,q)=pair_full(k,q)-(onebody(k,q)*onebody(mk,mq)-onebody(k,mq)*onebody(mk,q))
        !22
         pair_full(k+Nsite,q+Nsite)=pair_full(k+Nsite,q+Nsite)- &
     & (onebody(k+Nsite,q+Nsite)*onebody(mk+Nsite,mq+Nsite)-onebody(k+Nsite,mq+Nsite)*onebody(mk+Nsite,q+Nsite))
        !33
         pair_full(k+DNsite,q+DNsite)=pair_full(k+DNsite,q+DNsite)- &
     & 0.25d0*((onebody(k,q)*onebody(mk+Nsite,mq+Nsite)-onebody(k,mq+Nsite)*onebody(mk+Nsite,q)) &  
     &        +(onebody(k+Nsite,q+Nsite)*onebody(mk,mq)-onebody(k+Nsite,mq)*onebody(mk,q+Nsite)) &
     &        -(onebody(k,q+Nsite)*onebody(mk+Nsite,mq)-onebody(k,mq)*onebody(mk+Nsite,q+Nsite)) &
     &        -(onebody(k+Nsite,q)*onebody(mk,mq+Nsite)-onebody(k+Nsite,mq+Nsite)*onebody(mk,q)))
      end if
         !12
         pair_full(k,q+Nsite)=pair_full(k,q+Nsite)-(onebody(k,q+Nsite)*onebody(mk,mq+Nsite)- &
                        &  onebody(k,mq+Nsite)*onebody(mk,q+Nsite))
         !13
         pair_full(k,q+DNsite)=pair_full(k,q+DNsite)-(0.5d0*(onebody(k,q)*onebody(mk,mq+Nsite)- &
       & onebody(k,mq+Nsite)*onebody(mk,q))-0.5d0*(onebody(k,q+Nsite)*onebody(mk,mq)- &
       & onebody(k,mq)*onebody(mk,q+Nsite)))
         !23
         pair_full(k+Nsite,q+DNsite)=pair_full(k+Nsite,q+DNsite)-(0.5d0*(onebody(k+Nsite,q)*onebody(mk+Nsite,mq+Nsite)- &
       & onebody(k+Nsite,mq+Nsite)*onebody(mk+Nsite,q))-0.5d0*(onebody(k+Nsite,q+Nsite)*onebody(mk+Nsite,mq)- &
       & onebody(k+Nsite,mq)*onebody(mk+Nsite,q+Nsite))) 

   end do
end do
end subroutine minus_pairing_bg


subroutine diag_pair(M,N)
implicit none
integer,intent(IN)::N
complex(kind=8),intent(IN)::M(N,N)
complex(kind=8)::U(N,N)
real(kind=8)::V(N)
integer::i,j

U=M
call eigen(U,N,V)
write(*,*) "pair value:",V(N),V(1)
open(unit=33,file='pair_wf.dat',status='replace')
do i=1,N,1
   write(33,'(2E26.16)') dble(U(i,N)),dimag(U(i,N))
end do
close(33)
end subroutine diag_pair
