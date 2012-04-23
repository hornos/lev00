subroutine Tag(NumLin,NIONS,listat,list,len)
!.......................................................................
! analyses the Line with NumLin words (by LinPos and LinEnd), each word
! may contain either a number or two numbers separated by a dash, and
! compiles a list for tagging, namely, it sets:
!
!    listat(atoms in the list)     = .true.
!    listat(atoms not in the list) = .false.
!    list(1:len)  - contains a nicely written list of atoms
!.......................................................................
implicit none
integer LinEnd(100),LinPos(100)
character Line*200
logical listat(NIONS)
character(*) list
integer NumLin,len,i,i1,i2,i0,NIONS,l1,l2,iErr

!........... initialise the listat first

   listat=.false.

!............ now nicely ask

44 write(*,*) 'Specify atoms using space and/or comma to separate lists'
   write(*,*)'and a dash (without space) to indicate from/to'
   write(*,*)'within the same list, for example: 1-5,9, 7-8 15 25'
   read(*,'(a)') line
   call CutStr(Line,NumLin,LinPos,LinEnd,0,1,iErr)
   if(iErr.ne.0) go to 44
   if(NumLin.eq.0) then
      write(*,*)'WARNING! The tagging is disabled!'
      return
   end if

!............ produce the listat
   
   do i=1,NumLin
      i0=index(line(LinPos(i):LinEnd(i)),'-')
      if(i0.eq.0) then
         read(line(LinPos(i):LinEnd(i)),*) i1
         i2=i1
      else
         read(line(LinPos(i):LinPos(i)+i0-2),*) i1
         read(line(LinPos(i)+i0:LinEnd(i)),*) i2
      end if
      if(i1.gt.i2 .or. i1.lt.1 .or. i2.gt.NIONS) go to 44
      listat(i1:i2)=.true.
   end do

!__________ compose the nice list

   l1=1
   do i=1,NumLin
      l2=l1+LinEnd(i)-LinPos(i)
      if(i.ne.NumLin) then
         list(l1:l2+1)=line(LinPos(i):LinEnd(i))//','
      else
         list(l1:l2)=line(LinPos(i):LinEnd(i))
      end if
      l1=l2+2
   end do
   len=l2
   
   return
 end subroutine Tag

 subroutine TTag(TagYes,NIONS,TI,Labels,listat,list,len)
!.......................................................................
! compiles a list for tagging, namely, it sets:
!
!    listat(atoms in the list)     = .true.
!    listat(atoms not in the list) = .false.
!    list(1:len)  - contains a nicely written list of atoms
!.......................................................................
! A number of options are offered in a little menu, including:
!   - by coordinates
!   - by numbers
!.......................................................................
! If asked by numbers, then it analyses the Line with NumLin words
! (by LinPos and LinEnd), each word may contain either a number or two
! numbers separated by a dash, and
!.......................................................................
integer LinEnd(100),LinPos(100),i,i1,i2,j,iErr,NumLin,nat,NIONS,len,i0
logical listat(NIONS),TagYes
character(*) list
character Line*200,cha*10,chb*10,cha2*2,Labels(*)*2,ch,spec*2
real*8 TI(3,NIONS)
real*8,dimension(3) :: cp=(/0.0,0.0,0.0/)
real*8 :: Rad=1.0,x1,x2,y1,y2,z1,z2,d

!........... initialise the listat first
   listat=.false.
   TagYes=.false.
   x1=-1.e10 ; x2= 1.e10 ; y1=-1.e10 ; y2= 1.e10 ; z1=-1.e10 ; z2= 1.e10
!
!.......................................................................
!.............. menu starts here .......................................
!.......................................................................
!
21 write(*,'(/a/)')'========= CHOOSE how to tag atoms'
   if(TagYes) then
      write(*,*) '     >> Current set of tagged atoms: '
      call do_nice_list(listat,NIONS,list,len)
      do i=1,len,100
         i2=min(i+99,len)
         write(*,'(a)') list(i:i2)
      end do
   end if
   write(*,*)
   write(*,'(2x,3a)') char(47), '========================================',char(92)
   write(*,'(a)')' | You can execute each Produce option many |'
   write(*,'(a)')' |  times until the final list is obtained  |'
   write(*,'(2x,3a)') char(92), '========================================',char(47)
   write(*,*)
   write(*,'(a)')'  N. Produce: by atomic numbers'
   write(*,'(a)')' Sp. Produce: by species'
   write(*,'(a)')' Si. Produce: inside the sphere'
   write(*,'(a)')' So. Produce: outside the sphere'
   write(*,'(a)')' Cr. Produce: by coordinates X,Y,Z intervals'
   write(*,'(a)')'  E. Empty the list; you can start again!'

!_____ show settings
   write(*,*) ' Co. Show atoms with current tagging'
   write(*,'(a)') &
        '-------------- g e n e r a l  s e t t i n g s ----------------'
   write(*,'(3(a,f10.5),a)') &
        ' Cp. The central point: (',cp(1),',',cp(2),',',cp(3),')'
   write(*,'(a,f10.5)')'  R. The sphere radius: ',Rad
   call interval(x1,x2,cha,chb)
   write(*,'(a)')'  X. X-interval:         ('//cha//','//chb//')'
   call interval(y1,y2,cha,chb)
   write(*,'(a)')'  Y. Y-interval:         ('//cha//','//chb//')'
   call interval(z1,z2,cha,chb)
   write(*,'(a)')'  Z. Z-interval:         ('//cha//','//chb//')'
   write(*,'(a)')'  Q. Quit with the empty list'
   write(*,'(a)')'  P. Exit: quit with the current list'
   write(*,*)
   write(*,'(a)')'----------> Choose an appropriate option:'
   read (*,'(a)',err=21) cha2

![N]...... choose atoms by numbers

   IF(trim(cha2).eq.'N') THEN

44    write(*,*) &
           'Specify atoms using space and/or comma to separate lists'
      write(*,*)'and a dash (without space) to indicate from/to'
      write(*,*)'within the same list, for example: 1-5,9, 7-8 15 25'
      read(*,'(a)') line
      call CutStr(Line,NumLin,LinPos,LinEnd,0,1,iErr)
      if(iErr.ne.0) go to 44
      if(NumLin.eq.0) go to 21

!________________ produce the listat
      do i=1,NumLin
         i0=index(line(LinPos(i):LinEnd(i)),'-')
         if(i0.eq.0) then
            read(line(LinPos(i):LinEnd(i)),*) i1
            i2=i1
         else
            read(line(LinPos(i):LinPos(i)+i0-2),*) i1
            read(line(LinPos(i)+i0:LinEnd(i)),*) i2
         end if
         if(i1.gt.i2 .or. i1.lt.1 .or. i2.gt.NIONS) go to 44
         listat(i1:i2)=.true.
      end do
      TagYes=.true.

![Sp]...... choose atoms of the same species

   ELSE IF(trim(cha2).eq.'Sp') THEN

      write(*,*)'Choose the species [all atoms of this species will be selected]:'
      read(*,'(a)') spec
      i=0
      do nat=1,NIONS
         if(trim(adjustl(spec)).eq.trim(adjustl(Labels(nat)))) then
            i=i+1
            listat(nat)=.true.
         end if
      end do
      if(i.eq.0) then
         write(*,*)'ERROR: wrong species!'
      else
         write(*,*)'... Number of selected atoms = ',i
      end if
      if(i.gt.0) TagYes=.true.

![Si]...... choose atoms inside the sphere

   ELSE IF(trim(cha2).eq.'Si') THEN
      
      i=0
      do nat=1,NIONS
         d=sqrt( (TI(1,nat)-Cp(1))**2 + (TI(2,nat)-Cp(2))**2 + &
              (TI(3,nat)-Cp(3))**2 )
         if(d.le.Rad) then
            i=i+1
            listat(nat)=.true.
         end if
      end do
      if(i.gt.0) TagYes=.true.

![So]...... choose atoms outside the sphere

   ELSE IF(trim(cha2).eq.'So') THEN
      
      i=0
      do nat=1,NIONS
         d=sqrt( (TI(1,nat)-Cp(1))**2 + (TI(2,nat)-Cp(2))**2 + &
              (TI(3,nat)-Cp(3))**2 )
         if(d.gt.Rad) then
            i=i+1
            listat(nat)=.true.
         end if
      end do
      if(i.gt.0) TagYes=.true.

![Cr]...... choose atoms by coordinates

   ELSE IF(trim(cha2).eq.'Cr') THEN
      
      i=0
      do nat=1,NIONS
         if((TI(1,nat).gt.x1 .and. TI(1,nat).lt.x2) .and. &
              (TI(2,nat).gt.y1 .and. TI(2,nat).lt.y2) .and. &
              (TI(3,nat).gt.z1 .and. TI(3,nat).lt.z2) ) then
            i=i+1
            listat(nat)=.true.
         end if
      end do
      if(i.gt.0) TagYes=.true.

![E]...... empty the list

   ELSE IF(trim(cha2).eq.'E') THEN
      
      listat=.false.
      TagYes=.false.
      
![Co]...... show atoms with tagging

   ELSE IF(trim(cha2).eq.'Co') THEN

      write(*,*)' #  Species               Position            Tag'
      do i=1,NIONS
         ch=' '
         if(listat(i)) ch='Y'
         write(*,'(i5,x,a,4x,3(f10.5,x),2x,a)') i,Labels(i),(TI(j,i),j=1,3),ch
      end do
      write(*,*)'Hit ENTER when done ...'
      read(*,*)

![X]...... choose X interval

   ELSE IF(trim(cha2).eq.'X') THEN

      write(*,*) 'Specify the interval of atomic coordinates: A < X < B.'
1     write(*,*) &
           'Give the smallest (left) boundary of X (hit ENTER for -Infty):'
      call read_number_or_infty(x1,-1,iErr)
      if(iErr.ne.0) go to 1
2     write(*,*) &
           'Give the largest (right) boundary of X (hit ENTER for Infty):'
      call read_number_or_infty(x2,1,iErr)
      if(iErr.ne.0) go to 2

![Y]...... choose Y interval
      
   ELSE IF(trim(cha2).eq.'Y') THEN

      write(*,*) 'Specify the interval of atomic coordinates: A < Y < B.'
3     write(*,*) &
           'Give the smallest (left) boundary of Y (hit ENTER for -Infty):'
      call read_number_or_infty(y1,-1,iErr)
      if(iErr.ne.0) go to 3
4     write(*,*) &
           'Give the largest (right) boundary of Y (hit ENTER for Infty):'
      call read_number_or_infty(y2,1,iErr)
      if(iErr.ne.0) go to 4

![Z]...... choose Z interval

   ELSE IF(trim(cha2).eq.'Z') THEN

      write(*,*) 'Specify the interval of atomic coordinates: A < Z < B.'
5     write(*,*) &
           'Give the smallest (left) boundary of Z (hit ENTER for -Infty):'
      call read_number_or_infty(z1,-1,iErr)
      if(iErr.ne.0) go to 5
6     write(*,*) &
           'Give the largest (right) boundary of Z (hit ENTER for Infty):'
      call read_number_or_infty(z2,1,iErr)
      if(iErr.ne.0) go to 6
      
![Cp]...... choose Central point

   ELSE IF(trim(cha2).eq.'Cp') THEN

7     write(*,*)'Specify the Central Point:'
      read(*,*,err=7) Cp

![R]...... choose radius

   ELSE IF(trim(cha2).eq.'R') THEN

8     write(*,*)'Specify the radius of the sphere:'
      read(*,*,err=8) Rad

![Q]...... quit and do nothing

   ELSE IF(trim(cha2).eq.'Q') THEN

      listat=.false. ; TagYes=.false.
      return

![P]...... exit with the current list

   ELSE IF(trim(cha2).eq.'P') THEN
      return
   ELSE
      write(*,*)'Error! Try again!'
   END IF
   go to 21
end subroutine TTag

subroutine interval(x1,x2,cha,chb)
!......................................................................
! creates interval (x1,x2) in the character form as: (cha,chb)
!......................................................................
implicit none
real*8 x1,x2
character cha*10,chb*10
 
   if(x1.le.-1.0e10) then
      cha='-Infinity '
   else
      write(cha,'(f10.5)') x1
   end if
   if(x2.ge. 1.0e10) then
      chb='+Infinity '
   else
      write(chb,'(f10.5)') x2
   end if
end subroutine interval

subroutine read_number_or_infty(a,icase,iErr)
!......................................................................
!     reads a number or replaces it with:
!  icase=-1:  -infinity
!  icase= 1:  +infinity
!......................................................................
implicit none
integer iErr,icase,i0
real*8 a
character Line*20

  iErr=0
  if(icase.eq.-1) then
     write(*,*) 'Enter number or I for -Infinity:'
  else if(icase.eq.1) then
     write(*,*) 'Enter number or I for Infinity:'
  end if
  read(*,'(a)') line
  i0=index(line,'I')
  if(i0.ne.0) then
     if(icase.eq.-1) a=-1.0e10
     if(icase.eq. 1) a= 1.0e10
  else
     read(line,*,err=1) a
  end if
  return
1 iErr=1
end subroutine read_number_or_infty

subroutine do_sets_from_list(list_at,NIONS,n,set_beg,set_end)
!........................................................................
! create a nice list (of length len) out of a logical list of atoms 
! contained in list_at
!........................................................................
implicit none
logical list_at(*)
integer set_beg(100),set_end(100),n,i_prev,i,NIONS

!_____________ find continuous sets of numbers 
   n=0
   i_prev=-1
   do i=1,NIONS
      if(list_at(i)) then
         
         if(i.eq.i_prev+1) then
            set_end(n)=i
            i_prev=i
         else
            n=n+1
            if(n.gt.100) then
               write(*,*) 'FATAL! Too many sets > 100 !'
               stop 'in do_sets_from_list'
            end if
            set_beg(n)=i
            set_end(n)=i
            i_prev=i
         end if
         
      end if
   end do
end subroutine do_sets_from_list

subroutine do_nice_list(list_at,NIONS,list,len)
!........................................................................
! create a nice list (of length len) out of a logical list of atoms 
! contained in list_at
!........................................................................
implicit none
character(*)  list
character cha1*10,cha2*10
logical list_at(*)
integer set_beg(100),set_end(100),len,i1,i,n,i2,len1,NIONS,len2

!_____________ find continuous sets of numbers 

   call do_sets_from_list(list_at,NIONS,n,set_beg,set_end)

!_____________ compile the nice list
   if(n.eq.0) then
      len=1
      list(1:1)=' '
   else
      i1=0
      do i=1,n
         call number_to_char(set_beg(i),cha1,len1)
         if(set_end(i).eq.set_beg(i)) then
            i2=i1+len1
            list(i1+1:i2)=cha1(1:len1)
         else
            call number_to_char(set_end(i),cha2,len2)
            i2=i1+len1+len2+1
            list(i1+1:i2)=cha1(1:len1)//'-'//cha2(1:len2)
         end if
         i1=i2+1
         list(i1:i1)=','
      end do
      len=i2
   end if
 end subroutine do_nice_list

subroutine number_to_char(num,cha,len)
implicit none
character cha*10
integer i,num,len,l

   write(cha,'(i10)') num
   do i=1,10
      if(cha(i:i).ne.' ') then
         l=i
         go to 10
      end if
   end do
   len=10
   return
10 len=11-l
   do i=l,10
      cha(i-l+1:i-l+1)=cha(i:i)
      cha(i:i)=' '
   end do
 end subroutine number_to_char
