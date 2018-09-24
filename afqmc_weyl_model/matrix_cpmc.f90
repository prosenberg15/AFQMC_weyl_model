!------------------------------------------------------
!SUBROUTINE: QMC matrix manipulation
!AUTHOR:  Hao Shi
!VERSION: 10-Jan-2014
!NOTICE: This subroutine should compile with the module 
!caldet_module and matrixcal routines
!TYPE: Serial code
!EMAIL: boruoshihao@gmail.com
!COMMENT:
!REFERENCE:
!------------------------------------------------------


!This subroutine calculate <phi_l|phi_r>,Tanspose[Nsite*N(sigma)].Nsite*N(sigma)
subroutine deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: ovlpinv_temp(m,m)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)

 call zgemm('C','N',m,m,n,one,phi_l,n,phi_r,n,zero,ovlpinv_temp,m)
 !ovlpinv_temp= matmul(transpose(conjg(phi_l)), phi_r)
 !write(*,*) ovlpinv_temp
end subroutine



!This subroutine calculate Det[<phi_l|phi_r>]
subroutine deter_overlap_imp(n,m,phi_l,phi_r,imp)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT)::imp
complex(kind=8):: ovlpinv_temp(m,m)
complex(kind=8)::one=dcmplx(1.d0,0.d0)
complex(kind=8)::zero=dcmplx(0.d0,0.d0)

 call zgemm('C','N',m,m,n,one,phi_l,n,phi_r,n,zero,ovlpinv_temp,m)
 call caldet(m,ovlpinv_temp(1:m,1:m),imp)
end subroutine deter_overlap_imp




!This subroutine calcuate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r>,all the i,j matrix
!is wrote into Amax(n,n),<phi_l|phi_r>/=0
subroutine cal_Amat(n,m,phi_l,phi_r,Amat)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat(n,n)

complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)


complex(kind=8):: tempmat(n,m)
complex(kind=8):: ovlpinv_temp(m,m)
real(kind=8)::rdummy

  !ovlpinv_temp= matmul(transpose(conjg(phi_l)), phi_r)
  call deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)

  call inverse(ovlpinv_temp,m)

  call ZGEMM('N','N',n,m,m,one,phi_r,n,ovlpinv_temp,m,zero,tempmat,n)
  call ZGEMM('N','c',n,n,m,one,tempmat,n,phi_l,n,zero,Amat,n)

  Amat=transpose(Amat)

end subroutine cal_Amat


!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r>,here we use
!the information of ovlpinv, which is the inverse of
!<phi_l|phi_r>, <phi_l|phi_r>/=0
subroutine cal_Amat_withovlpinv(n,m,phi_l,phi_r,ovlpinv,Amat)
implicit none

integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m), ovlpinv(m,m)
complex(kind=8),intent(OUT):: Amat(n,n)

complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)


integer:: i,j
integer:: k,l

complex(kind=8):: tempmat(n,m)


   call ZGEMM('N','N',n,m,m,one,phi_r,n,ovlpinv,m,zero,tempmat,n)
   call ZGEMM('N','c',n,n,m,one,tempmat,n,phi_l,n,zero,Amat,n)

   Amat=transpose(Amat)

  ! G=RO^(-1)L^{dagger}
end subroutine cal_Amat_withovlpinv


!phi_r= phi.(phiT^{+}. phi)^{-1}
!phi_l= phiT
subroutine cal_Amat_step(n,m,phi_l,phi_r,Amat)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat(n,n)
complex(kind=8),parameter:: one=dcmplx(1d0,0d0)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)

 call ZGEMM('N','C',n,n,m,one,phi_r,n,phi_l,n,zero,Amat,n)

 Amat=transpose(Amat)

end subroutine cal_Amat_step


!phi_r= phi.(phiT^{+}. phi)^{-1}
!phi_l= phiT, only get the diagnal matrix
subroutine cal_Amat_diag(n,m,phi_l,phi_r,Amat)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat(n)
complex(kind=8),parameter:: zero=dcmplx(0d0,0d0)
integer::i,j
do i=1,n,1
   Amat(i)=zero
   do j=1,m,1
      Amat(i)=Amat(i)+phi_r(i,j)*conjg(phi_l(i,j))
   end do
end do
end subroutine cal_Amat_diag



!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> element,
!<phi_l|phi_r>/=0
subroutine cal_cidcj(n,m,phi_l,phi_r,Amat,i,j)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
complex(kind=8),intent(OUT):: Amat
integer,intent(IN)::i,j
complex(kind=8):: ovlpinv_temp(m,m)
integer::k,l
 !Get the inverse of phi_l over lap phi_r
  call deter_overlap(n,m,phi_l,phi_r,ovlpinv_temp)
  call inverse(ovlpinv_temp,m)

 !Get the element.
  Amat=dcmplx(0.d0,0.d0)
  do k=1,m,1
     do l=1,m,1
        Amat=Amat+phi_r(j,k)*ovlpinv_temp(k,l)*conjg(phi_l(i,l))
     end do
  end do
end subroutine cal_cidcj



!This subroutine also calculate <phi_l|ci^+ cj|phi_r>/<phi_l|phi_r> element.
!We use the ovlpinv to faster the simulation, <phi_l|phi_r>/=0
subroutine cal_cidcj_withovlpinv(n,m,phi_l,phi_r,ovlpinv,Amat,i,j)
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m),ovlpinv(m,m)
complex(kind=8),intent(OUT):: Amat
integer,intent(IN)::i,j
integer::k,l

 !Get the element.
  Amat=dcmplx(0.d0,0.d0)
  do k=1,m,1
     do l=1,m,1
        Amat=Amat+phi_r(j,k)*ovlpinv(k,l)*conjg(phi_l(i,l))
     end do
  end do
end subroutine cal_cidcj_withovlpinv




!Green function element, calculate <phi_l|ci^+ cj|phi_r>, it is a exact
!expand method, which calculate <phi_l|exp(ci^+cj)|phi_r>, can get the green
!function element by transfer, it can be used when <phi_l|phi_r>=0
subroutine cal_ij_exp(n,m,phi_l,phi_r,det,i,j,Amat)
use caldet_module
implicit none
integer,intent(IN)::n,m
complex(kind=8),intent(IN):: phi_l(n,m),phi_r(n,m)
character(len=1),intent(IN)::det
integer,intent(IN)::i,j
complex(kind=8),intent(OUT):: Amat
complex(kind=8)::phi_rn(n,m),ovp,ovpexp,ovlp(m,m) !ovp: <phi_l|phi_r> 
                                                  !ovpexp:<phi_l|exp(ni)|phi_r>
integer::k,p,q

!Get the overlap of phi_l.phi_r
if(det.eq.'Y') then !zero overlap
  ovp=dcmplx(0.d0,0.d0)
else if(det.eq.'N') then !None zero overlap
  call deter_overlap(n,m,phi_l,phi_r,ovlp)
  call caldet(m,ovlp,ovp)
else
write(*,*) "Something is wrong with det input in cal_ij_derivation"
stop
end if


if(i.eq.j) then
  !Get exp(ni).phi_r
  phi_rn(1:n,1:m)=phi_r(1:n,1:m)
  do p=1,m,1
     phi_rn(i,p)=phi_r(i,p)*exp(1.d0)
  end do
  !Get det(phi_l.exp(ni).phi_r)
  call deter_overlap(n,m,phi_l,phi_rn,ovlp)
  call caldet(m,ovlp,ovpexp)
  
  !Get the <phi_l|ci^+ ci|phi_r>
  Amat=dcmplx(1.d0/(exp(1.d0)-1.d0))*(ovpexp-ovp)
else
  !Get exp(ci^+cj).phi_r
  phi_rn(1:n,1:m)=phi_r(1:n,1:m)
  do p=1,m,1
     phi_rn(i,p)=phi_rn(i,p)+phi_r(j,p)
  end do
  !Get det(phi_l.exp(ci^+cj).phi_r)
  call deter_overlap(n,m,phi_l,phi_rn,ovlp)
  call caldet(m,ovlp,ovpexp)

  !Get the <phi_l|ci^+ cj|phi_r>
  Amat=ovpexp-ovp
end if
end subroutine cal_ij_exp
