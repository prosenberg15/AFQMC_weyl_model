!----------------------------------------------------------------
!Subroutine: lattice subroutine
!AUTHOR:  Hao Shi
!VERSION: 11-Jan-2014
!NOTICE: This routine set the hopt,sit,Hzero,coor,Tmatrix by
!the lattice size. It will need:set,Nsite,Nhop,Dimen,Nl(3),
!kbound(3),t1. If set=1, we will also need the input hop file, 
!We will need the matrixcal.f90 file to check the Hermition.
!TYPE: Serial
!EMAIL: boruoshihao@gmail.com
!COMMENT:
!REFERENCE:
!----------------------------------------------------------------

!This parameter contains all the one with lattice and Hzero
module lattice_param
implicit none
integer::set !set the hopt by yourself or by computer
integer::Nsite !the number of the whole sites
integer::DNsite !double of the whole sites
integer::Nbravais !number of sites in bravais lattice (i.e. number of unit cells)
integer::Nhop ! the number of the hoping terms need to consider
integer::openbcx, openbcy ! 1 for open BCs, 0 for other BCs
integer,parameter::Dimen=2 !the dimension
integer::Nl(3) !the number of different axis
real(kind=8)::kbound(3) !the twist boundary condition number
complex(kind=8)::t1      !Hubbard hopping t1 in nearest direction(-1).
complex(kind=8)::vhop,whop !v, w parameter from SSH model
complex(kind=8)::tyhop,tdhop !vertical and diagonal hopping b/w 1D SSH chains
real(kind=8)::lamda   !spin orbit coupling parameters-some in the note

complex(kind=8),allocatable::hopt(:) !hopt and sit record the information
integer,allocatable::sit(:,:)        !of hoping term in different sites
integer,allocatable::coor(:,:) ! we label the site, it record the coordinate
integer,allocatable::Tmatrix(:,:) !Use to store the nearest hopping in different direction.
complex(kind=8),allocatable::Hzero(:,:)    !2*Nsite,2*Nsite,the Hezo Hamiltonian of the lattice
end module lattice_param



!This subroutine get the Hamilton of the no interaction lattice
subroutine set_lattice()
use lattice_param
use model_param
use mpi_serial_param
implicit none
integer::mi,mj
integer::i,j

 !Set the number of lattice and Nhop
if(set.NE.1) then
   Nsite=1
   do mi=1,Dimen,1
      Nsite=Nsite*Nl(mi)
   end do
   Nhop=Nsite*2*Dimen*4
   if(dtype.eq.'w') then
      Nbravais=1
      do mi=1,Dimen,1
         Nbravais=Nbravais*Nl(mi)
      end do
      !Nbravais unit cells, with 2 sites (i.e. 2 bands) per unit cell
      Nsite=2*Nbravais
      Nhop=16*Nbravais
   endif

   !Set Nl to 1 for non-relevant dimension
   if(Dimen.eq.1) then
     Nl(2)=1;Nl(3)=1
   else if(Dimen.eq.2) then
     Nl(3)=1
   end if
 end if
 DNsite=2*Nsite

if(rank.eq.0) write(*,*) 'Nbravias: ', Nbravais
if(rank.eq.0) write(*,*) 'Nsites: ', Nsite

 call allocate_lattice_array()
 call set_lattice_hop()

 Hzero=(0.d0,0.d0)
 do mi=1,Nhop,1
    Hzero(sit(mi,1),sit(mi,2))=Hzero(sit(mi,1),sit(mi,2))+hopt(mi)
    Hzero(sit(mi,1)+Nsite,sit(mi,2)+Nsite)=Hzero(sit(mi,1)+Nsite,sit(mi,2)+Nsite)+hopt(mi)
 end do

!if (rank.eq.0) then
!    do i=1,Nsite,1
!       do j=1,Nsite,1
!          !write(*,'(100(A1,F5.2,A1,F5.2,A1))') ('(',real(Hzero(i,j)),',',imag(Hzero(i,j)),')',j=1,DNsite)
!          write(*,*) i,j,real(Hzero(i,j)),imag(Hzero(i,j)), real(Hzero(j,i)), imag(Hzero(j,i))
!!       !write(*,'(100(F5.2,A3,F5.2,A3))') (real(Hzero(i,j)),'+i*',imag(Hzero(i,j)),' ',j=1,DNsite)
!       end do
!    end do
! end if

 call check_Hermite_c(Hzero,DNsite)

 if(set.NE.1) then
   call set_Tmatrix()
 end if

 !------------------------------------------------
 !If (set.EQ.1): We need to pay attention that:
 !coor(:,:) and Tmatrix(:,:) will not be set here.
 !------------------------------------------------
 
end subroutine set_lattice



!-----------------------------------------------------------
!This subroutine get the sit(Nhop,2) and hopt(Nhop) and coor
!-----------------------------------------------------------
subroutine set_lattice_hop()
use model_param
use lattice_param
implicit none
integer::mi,mj

if(set.EQ.1) then
  open(unit=10,file='hop',status='old')
    do mi=1,Nhop,1
       read(10,*) sit(mi,1),sit(mi,2),hopt(mi)
    end do
  close(10)
else if(set.EQ.2) then
   call setnumber()
   if(dtype.ne.'w') then
      call sethopTBC()
   elseif(dtype.eq.'w') then
      call sethopTBC_multi_open()
   endif
else
  write(*,*) "something is wrong with the set"
  call mystop
end if

end subroutine set_lattice_hop


!-----------------------------------------
!we label the lattice by 1~Nsite,set
!by (x,y,z) the coor(Nsite,Dimen) is
!dependent on the dimension of the lattice
!It can be think as the coordinate of the
!lattice point is from 1~Nl(:)
!-----------------------------------------
subroutine setnumber
use model_param
use lattice_param
implicit none
integer::i,j,k
integer::ntemp,den

if(dtype.ne.'w') then
   do i=1,Nsite,1
      ntemp=i-1
      do j=Dimen,1,-1
         den=1
         do k=1,j-1,1
            den=den*Nl(k)
         end do !den is z coor Nl(1)*Nl(2)
                !       y coor Nl(1)
                !       x coor 1
         coor(i,j)=ntemp/den
         ntemp=ntemp-coor(i,j)*den
         coor(i,j)=coor(i,j)+1
      end do
   end do
elseif(dtype.eq.'w') then
   do i=1,Nbravais,1
      ntemp=i-1
      do j=Dimen,1,-1
         den=1
         do k=1,j-1,1
            den=den*Nl(k)
         end do !den is z coor Nl(1)*Nl(2)
                !       y coor Nl(1)
                !       x coor 1
         coor(i,j)=ntemp/den
         ntemp=ntemp-coor(i,j)*den
         coor(i,j)=coor(i,j)+1
      end do
   end do
   do i=1,Nbravais,1
      do j=Dimen,1,-1
         coor(i+Nbravais,j)=coor(i,j)
      enddo
   enddo
endif

 
end subroutine setnumber


subroutine sethopTBC_multi
use param
use lattice_param
implicit none
integer,external::latt_label
integer,external::bound
integer::coord(Dimen)
integer::i,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

   !inside unit cell hopping
   Nh=Nh+1
   sit(Nh,1)=i
   sit(Nh,2)=i+Nbravais   ! A_nm --> B_nm
   hopt(Nh)=vhop

   Nh=Nh+1
   sit(Nh,1)=i+Nbravais
   sit(Nh,2)=i            ! B_nm --> A_nm
   hopt(Nh)=vhop

   !outside unit cell hopping
   do j=1,Dimen,1

      den=1
      do k=1,j-1,1
         den=den*Nl(k)
      end do

      if(coor(i,j).EQ.Nl(j)) then
         ntemp=(1-Nl(j))*den+i
      else
         ntemp=i+den
      endif

      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp            ! B_nm --> A_{n+1}m
         hopt(Nh)=whop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp            ! A_nm --> A_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais   ! B_nm --> B_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais   ! A_nm --> B_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m+1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m+1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp   ! B_nm --> A_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      endif
      
      if(coor(i,j).EQ.1) then
         ntemp=(Nl(j)-1)*den+i
      else
         ntemp=i-den
      endif

      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_{n+1}m --> B_nm
         hopt(Nh)=whop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp                ! A_n{m+1} --> A_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais       ! B_n{m+1} --> B_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_nm --> B_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))

         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m-1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m-1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp              ! B_nm --> A_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
      endif
      
   enddo
enddo
if(Nh.ne.Nhop) then
   write(*,*) "Something is wrong with Nhop PBC"
   call mystop
end if

end subroutine sethopTBC_multi

subroutine sethopTBC_multi_open
use param
use lattice_param
implicit none
integer,external::latt_label
integer,external::bound
integer::coord(Dimen)
integer::i,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

   !inside unit cell hopping
   Nh=Nh+1
   sit(Nh,1)=i
   sit(Nh,2)=i+Nbravais   ! A_nm --> B_nm
   hopt(Nh)=vhop

   Nh=Nh+1
   sit(Nh,1)=i+Nbravais
   sit(Nh,2)=i            ! B_nm --> A_nm
   hopt(Nh)=vhop

   !outside unit cell hopping
   do j=1,Dimen,1

      den=1
      do k=1,j-1,1
         den=den*Nl(k)
      end do

      if(coor(i,j).EQ.Nl(j)) then
         ntemp=(1-Nl(j))*den+i
      else
         ntemp=i+den
      endif

      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp            ! B_nm --> A_{n+1}m
         hopt(Nh)=whop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcx.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp            ! A_nm --> A_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais   ! B_nm --> B_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais   ! A_nm --> B_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m+1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if((openbcx.eq.1).and.(coor(i,1).eq.1)) hopt(Nh)=zero
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m+1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if((openbcx.eq.1).and.(coor(i,1).eq.Nl(1))) hopt(Nh)=zero
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp   ! B_nm --> A_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.Nl(j))) hopt(Nh)=zero
         
      endif
      
      if(coor(i,j).EQ.1) then
         ntemp=(Nl(j)-1)*den+i
      else
         ntemp=i-den
      endif

      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_{n+1}m --> B_nm
         hopt(Nh)=whop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcx.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp                ! A_n{m+1} --> A_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais       ! B_n{m+1} --> B_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_nm --> B_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m-1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if((openbcx.eq.1).and.(coor(i,1).eq.1)) hopt(Nh)=zero
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m-1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if((openbcx.eq.1).and.(coor(i,1).eq.Nl(1))) hopt(Nh)=zero
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp              ! B_nm --> A_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if((openbcy.eq.1).and.(coor(i,j).eq.1)) hopt(Nh)=zero
      endif
      
   enddo
enddo
if(Nh.ne.Nhop) then
   write(*,*) "Something is wrong with Nhop PBC"
   call mystop
end if

end subroutine sethopTBC_multi_open

subroutine sethopTBC_multi_open_v2
use param
use lattice_param
implicit none
integer,external::latt_label
integer,external::bound
integer::coord(Dimen)
integer::i,j,k,ntemp,den
integer::Nh

Nh=0
do i=1,Nbravais,1

   !inside unit cell hopping
   Nh=Nh+1
   sit(Nh,1)=i
   sit(Nh,2)=i+Nbravais   ! A_nm --> B_nm
   hopt(Nh)=vhop

   Nh=Nh+1
   sit(Nh,1)=i+Nbravais
   sit(Nh,2)=i            ! B_nm --> A_nm
   hopt(Nh)=vhop

   !outside unit cell hopping
   do j=1,Dimen,1

      den=1
      do k=1,j-1,1
         den=den*Nl(k)
      end do

      if(coor(i,j).EQ.Nl(j)) then
         ntemp=(1-Nl(j))*den+i
      else
         ntemp=i+den
      endif

      if((j.eq.1).and.(coor(i,j).eq.Nl(j))) then
         hopt(Nh)=zero
      endif
      if((j.eq.1).and.(coor(i,j).eq.Nl(j))) then
         hopt(Nh)=zero
      endif
      
      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp            ! B_nm --> A_{n+1}m
         hopt(Nh)=whop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp            ! A_nm --> A_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais   ! B_nm --> B_n{m+1}
         hopt(Nh)=tyhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais   ! A_nm --> B_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m+1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if(coor(i,1).eq.1) hopt(Nh)=zero
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m+1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if(coor(i,1).eq.Nl(1)) hopt(Nh)=zero
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp   ! B_nm --> A_n{m+1}
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.Nl(j)) hopt(Nh)=zero
         
      endif
      
      if(coor(i,j).EQ.1) then
         ntemp=(Nl(j)-1)*den+i
      else
         ntemp=i-den
      endif

      if(j.eq.1) then
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_{n+1}m --> B_nm
         hopt(Nh)=whop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.1) hopt(Nh)=zero
      else
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp                ! A_n{m+1} --> A_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.1) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp+Nbravais       ! B_n{m+1} --> B_nm
         hopt(Nh)=tyhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.1) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=ntemp+Nbravais       ! A_nm --> B_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.1) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i
         if(mod(ntemp-1,Nl(1)).eq.0)then
            sit(Nh,2)=(ntemp-1)+Nl(1)+Nbravais
         else
            sit(Nh,2)=ntemp+Nbravais-1 ! A_nm --> B_{n-1}{m-1}
         endif
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))+kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if(coor(i,1).eq.1) hopt(Nh)=zero
         if(coor(i,j).eq.1) hopt(Nh)=zero

         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         if(mod(ntemp,Nl(1)).eq.0)then
            sit(Nh,2)=ntemp+1-Nl(1)
         else
            sit(Nh,2)=ntemp+1         ! B_nm --> A_{n+1}{m-1}
         endif
         !sit(Nh,2)=ntemp+Nbravais   
         hopt(Nh)=tdhop*exp((0.d0,1.d0)*(kbound(1)*(2.d0*Pi/dble(Nl(1)))-kbound(2)*(2.d0*Pi/dble(Nl(2)))))
         if(coor(i,1).eq.Nl(1)) hopt(Nh)=zero
         if(coor(i,j).eq.1) hopt(Nh)=zero
         
         Nh=Nh+1
         sit(Nh,1)=i+Nbravais
         sit(Nh,2)=ntemp              ! B_nm --> A_n{m-1}
         hopt(Nh)=tdhop*exp((0.d0,-1.d0)*kbound(j)*2.d0*Pi/dble(Nl(j)))
         if(coor(i,j).eq.1) hopt(Nh)=zero

      endif
      
   enddo
enddo
if(Nh.ne.Nhop) then
   write(*,*) "Something is wrong with Nhop PBC"
   call mystop
end if

end subroutine sethopTBC_multi_open_v2


!----------------------------------------------------------------------
!we get the hoping matrix by Twist boundary condition,site(:,2),hopt(:)
!----------------------------------------------------------------------
subroutine sethopTBC
use param
use lattice_param
implicit none
integer,external::latt_label
integer,external::bound
integer::coord(Dimen)
integer::i,j,k,ntemp,den
integer::Nh

!write(*,*) 'openbcx: ', openbcx
!write(*,*) 'openbcy: ', openbcy

Nh=0
do i=1,Nsite,1

   if ((openbcx.eq.1).and.(openbcy.eq.1)) then
      !Ci{+}C{i+x}
      coord(1)=bound(coor(i,1)+1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.1) then
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1))) 

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      
         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=-lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.1) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      
         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i-x}
      coord(1)=bound(coor(i,1)-1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.1) then
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=-lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.1.) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i+y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)+1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.Nl(1)) then
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.Nl(1)) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i-y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)-1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.Nl(1)) then
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.Nl(1)) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

   elseif ((openbcx.eq.1).and.(openbcy.ne.1)) then
      !Ci{+}C{i+x}
      coord(1)=bound(coor(i,1)+1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.1) then
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1))) 

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      
         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=-lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.1) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      
         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i-x}
      coord(1)=bound(coor(i,1)-1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.1) then
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=-lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.1.) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i+y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)+1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      !Ci{+}C{i-y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)-1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

   elseif ((openbcx.ne.1).and.(openbcy.eq.1)) then
      !Ci{+}C{i+x}
      coord(1)=bound(coor(i,1)+1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1))) 

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=-lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      !Ci{+}C{i-x}
      coord(1)=bound(coor(i,1)-1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=-lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      !Ci{+}C{i+y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)+1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.Nl(1)) then
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.Nl(1)) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

      !Ci{+}C{i-y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)-1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      if (abs(sit(Nh,1)-sit(Nh,2)).eq.Nl(1)) then
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))
      else if (abs(sit(Nh,1)-sit(Nh,2)).ne.Nl(1)) then
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i+Nsite
         sit(Nh,2)=latt_label(coord(1:2))
         hopt(Nh)=0.d0

         Nh=Nh+1
         sit(Nh,1)=i
         sit(Nh,2)=latt_label(coord(1:2))+Nsite
         hopt(Nh)=0.d0
      end if

   else if((openbcx.ne.1).and.(openbcy.ne.1)) then

      !Ci{+}C{i+x}
      coord(1)=bound(coor(i,1)+1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1))) 

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=-lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=lamda*exp((0.d0,1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      !Ci{+}C{i-x}
      coord(1)=bound(coor(i,1)-1,Nl(1))
      coord(2)=coor(i,2)

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=-lamda*exp((0.d0,-1.d0)*kbound(1)*2.d0*Pi/dble(Nl(1)))

      !Ci{+}C{i+y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)+1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=-lamda*Xi*exp((0.d0,1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      !Ci{+}C{i-y}
      coord(1)=coor(i,1)
      coord(2)=bound(coor(i,2)-1,Nl(2))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=t1*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i+Nsite
      sit(Nh,2)=latt_label(coord(1:2))
      hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))

      Nh=Nh+1
      sit(Nh,1)=i
      sit(Nh,2)=latt_label(coord(1:2))+Nsite
      hopt(Nh)=lamda*Xi*exp((0.d0,-1.d0)*kbound(2)*2.d0*Pi/dble(Nl(2)))
   end if
end do


if(Nh.ne.Nhop) then
  write(*,*) "Something is wrong with Nhop PBC"
  call mystop
end if
end subroutine sethopTBC


!-------------------------------------------
!get the shift matrix in different direction
!-------------------------------------------
subroutine set_Tmatrix
use lattice_param
implicit none
integer,external::bound
integer::i,j,k
integer::NT(Dimen)
integer::tmp,dem
do i=1,Nsite,1
   do j=1,Dimen,1

      tmp=1;dem=1
      do k=1,Dimen,1
         NT(k)=coor(i,k)
         if(k.eq.j) then
           NT(k)=bound(NT(k)+1,Nl(k))
         end if

         tmp=tmp+(NT(k)-1)*dem
         dem=dem*Nl(k)
      end do

      Tmatrix(i,j)=tmp

      if(tmp.GT.Nsite.OR.tmp.LT.1) then
        write(*,*) "Something is wrong with the Tmatrix",tmp
        call mystop
      end if
   end do
end do
end subroutine set_Tmatrix


!--------------------------------------------------------
!It is a bound function use to fit the boundary condtion.
!--------------------------------------------------------
integer function bound(i,ni)
integer::i,ni
bound=i
if(i.GT.ni) bound=i-ni
if(i.LT.1) bound=i+ni
if(bound.GT.ni.OR.bound.LT.1) write(*,*) "Something is wrong is the boundary condition"
end function


!-------------------------------------------------------
!give in coord(Dimen), give out the number label 1~Nsite
!-------------------------------------------------------
integer function latt_label(coord)
use lattice_param
implicit none
integer,intent(IN)::coord(1:Dimen)
integer::den,i

latt_label=1
den=1
do i=1,Dimen,1
   latt_label=latt_label+(coord(i)-1)*den
   den=den*Nl(i)
end do
if(latt_label.LT.1.OR.latt_label.GT.Nsite) then
  write(*,*) "Something is wrong with latt_label output:",latt_label
  call mystop
end if
end function


!input k, get q=-k
subroutine inverse_momentum(k,q)
use lattice_param
implicit none
integer,intent(IN)::k
integer,intent(OUT)::q
integer,external::latt_label
integer,external::bound
integer::coord(1:Dimen)
integer::i
do i=1,Dimen,1
   coord(i)=bound(2-coor(k,i),Nl(i))
end do
q=latt_label(coord)
end subroutine inverse_momentum


!-----------------------------------------------------------------
!This subroutine allocate the arrays we need to use in set lattice
!-----------------------------------------------------------------
subroutine allocate_lattice_array()
use lattice_param
implicit none
allocate(hopt(Nhop))
allocate(sit(Nhop,2))
allocate(coor(Nsite,Dimen))
allocate(Hzero(DNsite,DNsite))
allocate(Tmatrix(Nsite,Dimen))
end subroutine allocate_lattice_array


!-------------------------------------------------------------------
!This subroutine deallocate the arrays we need to use in set lattice
!-------------------------------------------------------------------
subroutine deallocate_lattice_array()
use lattice_param
implicit none
if(allocated(hopt)) deallocate(hopt)
if(allocated(sit)) deallocate(sit)
if(allocated(coor)) deallocate(coor)
if(allocated(Hzero)) deallocate(Hzero)
if(allocated(Tmatrix)) deallocate(Tmatrix)
end subroutine deallocate_lattice_array




!test
!subroutine test()
!use lattice_param
!implicit none
!integer,external::latt_label
!integer,external::bound
!integer::i,j,m,n
!integer::cc(1:Dimen),ctmp
!do i=1,Nsite,1
!   do j=1,Nsite,1
!
!      do m=1,Dimen,1
!         ctmp=coor(j,m)-coor(i,m)+1
!         cc(m)=bound(ctmp,Nl(m))
!      end do
!      n=latt_label(cc(1:Dimen))
!      write(*,*) i,j,n;pause
!   end do
!end do
!write(*,*) coor(1,1:Dimen)
!write(*,*) latt_label(coor(1,1:Dimen))
!call mystop
!end subroutine test
!end test
