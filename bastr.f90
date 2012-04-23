subroutine bastr(DIR,REC,VOL,iflag)
  implicit none
!...........................................................................
!  Using vectors DIR(#,xyz), it calculates:
! (1) vectors REC(#,xyz) which are reciprocal to DIR in the following sense:
!
!         scalar product DIR(#1)*REC(#2) = Kroneker's delta(#1,#2) * factor
!
!     where factor=2*pi if iflag=1 and factor=1 otherwise.
! (2) volume VOL of the cell built using DIR
!     (i.e. without 2*pi).
!...........................................................................
  real*8 DIR(3,3),REC(3,3),factor,VOL
  integer iflag
  real*8, parameter :: TwoPi=6.2831853072,tiny=0.00001

  factor=1.0
  if(iflag.eq.1) factor=TwoPi
!
  REC(1,1)=DIR(2,2)*DIR(3,3)-DIR(3,2)*DIR(2,3)
  REC(1,2)=DIR(2,3)*DIR(3,1)-DIR(3,3)*DIR(2,1)
  REC(1,3)=DIR(2,1)*DIR(3,2)-DIR(3,1)*DIR(2,2)
  REC(2,1)=DIR(3,2)*DIR(1,3)-DIR(1,2)*DIR(3,3)
  REC(2,2)=DIR(3,3)*DIR(1,1)-DIR(1,3)*DIR(3,1)
  REC(2,3)=DIR(3,1)*DIR(1,2)-DIR(1,1)*DIR(3,2)
  REC(3,1)=DIR(1,2)*DIR(2,3)-DIR(2,2)*DIR(1,3)
  REC(3,2)=DIR(1,3)*DIR(2,1)-DIR(2,3)*DIR(1,1)
  REC(3,3)=DIR(1,1)*DIR(2,2)-DIR(2,1)*DIR(1,2)
  VOL=DIR(1,1)*REC(1,1)+DIR(1,2)*REC(1,2)+DIR(1,3)*REC(1,3)
  if(abs(vol).lt.tiny) then
     write(*,*)'ERROR! Collinear vectors!'
     stop
  end if
  REC=factor*REC/VOL
  VOL=abs(VOL)
  return
end subroutine bastr
