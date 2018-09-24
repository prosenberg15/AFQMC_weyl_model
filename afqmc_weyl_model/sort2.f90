!This subroutine use Heapsort method to sort the array of ra.
!n: input,integer, the length of the array ra and rb
!ra: input,output, real(kind=8)::ra(n) 
!rb: input,output, integer::rb(n)
!Output ra: the new array which is arranged from smaller value to large value.
!Output rb: rb(i) is to the number of new ra(i) in old ra.

subroutine sort2(n,ra,rb)
implicit none
integer,intent(IN)::n
real(kind=8),intent(INOUT)::ra(n)
integer,intent(OUT)::rb(n)
real(kind=8)::rra
integer::rrb
integer::i,j,l,ir

do i=1,n,1
   rb(i)=i
end do

l=n/2+1

ir=n

do

  if(l>1) then

    l=l-1

    rra=ra(l)

    rrb=rb(l)

  else

    rra=ra(ir)

    rrb=rb(ir)

    ra(ir)=ra(1)

    rb(ir)=rb(1)

    ir=ir-1

    if(ir==1) then

      ra(1)=rra

      rb(1)=rrb

      return

    endif

  endif

  i=l

  j=l+l

  do while(j<=ir) 

    if(j.lt.ir) then

      if(ra(j)<ra(j+1)) j=j+1

    endif

    if(rra<ra(j)) then

      ra(i)=ra(j)

      rb(i)=rb(j)

      i=j

      j=j+j

    else

      j=ir+1

    endif

  end do

  ra(i)=rra

  rb(i)=rrb

end do

END SUBROUTINE sort2
