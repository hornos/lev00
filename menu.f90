module menu
!.....................................................................
! 'icase' gives the method used to specify the plane:
!  icase = 1 - by a normal vector;
!          2 - by 3 points lying in the plane.
!  angstr - gives a method of specifying coordinates
!.....................................................................
integer :: nresol=160,nresol3=160,nresols=30,nresold=30,Nrad=0
integer :: nresol_prv=30,nclasses=30,icase=1,icase1=1
real*8  :: NumbEl=8.0
real*8  :: vers0x=0.0,vers0y=0.0,vers0z=0.0
real*8  :: vers1x=0.0,vers1y=0.0,vers1z=0.0,vers2x=0.0,vers2y=0.0,vers2z=0.0
real*8  :: vers3x=0.0,vers3y=0.0,vers3z=0.0,aCENTX=0.0,aCENTY=0.0,aCENTZ=0.0
real*8  :: width1=0.0,width2=0.0
real*8  :: multcon=1.0,lochop=0.0,hichop=0.0
real*8,dimension(3) ::    Ra=(/0.0,0.0,0.0/),Rb,Rc
character(len=7)  :: type_prv='3d-col-'
character(len=1)  :: way_res
character(len=10) :: method='conserving'
character(len=12) :: angstr='<Angstroms> '
logical :: central_p=.false.
!
! Point() is used also in [module dos]; should not be a problem
! as these two options (dos and density) are never called together
!
integer,parameter :: Num_Sph=5
real*8 :: Point(3,Num_Sph),Shrink(Num_Sph),RadiusS=0.0,RadiusL=0.0

CONTAINS

subroutine choose1()
!.......... for a 2-dimensional plot ...............................
implicit none
character enable*8
integer item
!.............. show present setting
 1    write(*,*)'..............CUSTOMISATION ......................'
      write(*,*)'..... Change these parameters if necessary:...'
      write(*,*)
      write(*,'(a31,i5)')'   1. Resolution for the plot = ',nresol
      if(lochop.eq.0.0) then
           enable='disabled'
      else
           enable=' enabled'
      end if
      write(*,'(a18,f6.2,a11)') &
                       '   2. Low chop  = ',lochop,' {'//enable//'}'
      if(hichop.eq.0.0) then
           enable='disabled'
      else
           enable=' enabled'
      end if
      write(*,'(a18,f6.2,a11)') &
                       '   3. High chop = ',hichop,' {'//enable//'}'
      write(*,'(a30,f6.2)')'   4. Multiplication factor = ',multcon
      write(*,'(a)')'   5. Return to the previous menu.'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ give the resolution in either direction
      if(item.eq.1) then
 15      write(*,*)'Specify this number:'
         read(*,*,err=16) nresol
         if(nresol.lt.2) go to 16
         go to 1
 16      write(*,*)'Incorrect or too small resolution! Try again!'
         go to 15
!
!__________ give low chop value
      else if(item.eq.2) then
 20      write(*,*)'Specify this number:'
         read(*,*,err=23) lochop
         go to 1
 23      write(*,*)'Error! Try again!'
         go to 20
!
!__________ give high chop value
      else if(item.eq.3) then
 25      write(*,*)'Specify this number:'
         read(*,*,err=26) hichop
         go to 1
 26      write(*,*)'Error! Try again!'
         go to 25
!
!__________ give a multiplication factor for the density
      else if(item.eq.4) then
 35      write(*,*)'Specify this number:'
         read(*,*,err=37) multcon
         if(multcon.lt.0.0) go to 37
         go to 1
 37      write(*,*)'Error! Try again!'
         go to 35
      else if(item.eq.5) then
         return
      else
         go to 100
      end if
!............. error
 100  write(*,*)'Incorrect item number! Try again!'
      go to 1
end subroutine choose1

subroutine choose3()
!.......... for a 3-dimensional plot ...............................
!...................................................................
implicit none
character enable*8,cha2*2

!.............. show present setting
 1    write(*,*)'..............MENU ONE .......................'
      write(*,*)'..... Change these parameters if necessary:...'
      write(*,*)
      write(*,'(a31,i5)')'   1. Resolution for the plot = ',nresol3
      write(*,'(a38,i5)') &
            '   2. Resolution for the previewing = ',nresol_prv
      if(lochop.eq.0.0) then
           enable='disabled'
      else
           enable=' enabled'
      end if
      write(*,'(a18,f6.2,a11)') &
                       '   4. Low chop  = ',lochop,' {'//enable//'}'
      if(hichop.eq.0.0) then
           enable='disabled'
      else
           enable=' enabled'
      end if
      write(*,'(a18,f6.2,a11)') &
                       '   5. High chop = ',hichop,' {'//enable//'}'
      write(*,'(a30,f6.2)')'   6. Multiplication factor    = ',multcon

      if(type_prv.eq.'2d-old-') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 2D + contours'
         write(*,*) ' cl. PREVIEW/Ps: number of contour levels = ',nclasses
      else if(type_prv.eq.'2d-col-') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 2D + colours'
      else if(type_prv.eq.'2d-gray') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 2D + gray palette'
      else if(type_prv.eq.'3d-old-') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 3D + grid lines'
      else if(type_prv.eq.'3d-c-2d') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 3D + 2D underneath; in colour'
      else if(type_prv.eq.'3d-g-2d') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 3D + 2D underneath; in gray'
      else if(type_prv.eq.'3d-col-') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 3D; in colour'
      else if(type_prv.eq.'3d-col-') then
         write(*,'(a,a7)') '  3D. PREVIEW/Ps type: 3D; in colour'
      end if

      write(*,'(a)')'   Q. Return to the previous menu.'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
      IF(trim(cha2).eq.'cl' .and. type_prv.eq.'2d-old-')  THEN
 75      write(*,*)'Specify this number:'
         read(*,*,err=75) nclasses
         if(nclasses.lt.2) go to 75
         go to 1
!
!__________ give the resolution in either direction
       ELSE IF(trim(cha2).eq.'1') THEN
 15      write(*,*)'Specify this number:'
         read(*,*,err=16) nresol3
         if(nresol3.lt.2) go to 16
         go to 1
 16      write(*,*)'Incorrect or too small resolution! Try again!'
         go to 15
!
!__________ give the resolution in either direction
      ELSE IF(trim(cha2).eq.'2') THEN
 45      write(*,*)'Specify this number:'
         read(*,*,err=46) nresol_prv
         if(nresol_prv.lt.2) go to 46
         go to 1
 46      write(*,*)'Incorrect or too small resolution! Try again!'
         go to 45
!
!__________ give the type of the previewing
      ELSE IF(trim(cha2).eq.'3D') THEN
         if(type_prv.eq.'2d-old-') then
            type_prv='2d-col-'
         else if(type_prv.eq.'2d-col-') then
            type_prv='2d-gray'
         else if(type_prv.eq.'2d-gray') then
            type_prv='3d-old-'
         else if(type_prv.eq.'3d-old-') then
            type_prv='3d-c-2d'
         else if(type_prv.eq.'3d-c-2d') then
            type_prv='3d-g-2d-'
         else if(type_prv.eq.'3d-g-2d') then
            type_prv='3d-col-'
         else if(type_prv.eq.'3d-col-') then
            type_prv='2d-old-'
         end if
         go to 1
!
!__________ give low chop value
      ELSE IF(trim(cha2).eq.'4') THEN
 20      write(*,*)'Specify this number:'
         read(*,*,err=23) lochop
         go to 1
 23      write(*,*)'Error! Try again!'
         go to 20
!
!__________ give high chop value
      ELSE IF(trim(cha2).eq.'5') THEN
 25      write(*,*)'Specify this number:'
         read(*,*,err=26) hichop
         go to 1
 26      write(*,*)'Error! Try again!'
         go to 25
!
!__________ give a multiplication factor for the density
      ELSE IF(trim(cha2).eq.'6') THEN
 35      write(*,*)'Specify this number:'
         read(*,*,err=37) multcon
         if(multcon.lt.0.0) go to 37
         go to 1
 37      write(*,*)'Error! Try again!'
         go to 35
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         return
      ELSE 
         go to 100      
      END IF
!............. error
 100  write(*,*)'Incorrect item number! Try again!'
      go to 1
end subroutine choose3

end module menu

