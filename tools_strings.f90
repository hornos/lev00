     subroutine CutStr(Line,NumLin,LinPos,LinEnd,NFIL,iCom,iErr)
!.................................................................
!....It reads the line from the file UNIT=NFIL as a string   .....
!.... Line*80 (if NFIL.ne.0) and replies with a set of       .....
!.... substrings. Their number is NumLin, while starting and .....
!.... ending positions are LinPos and LinEnd.                .....
!.................................................................
! Up to 100 substrings is implied to be in the initial string Line.
!.................................................................
!  If NFIL=0 it is assumed that the string Line*80 exists and
!  must not be read from the external source, but must be taken
!  as an input from the parent program
!.................................................................
!  iErr=1 in the case of error input; oterwise iErr=0
!.................................................................
!  Commas are also interpreted as string division symbols 
!  alongside spaces if iCom=1; not, if iCom=/=1.
!.................................................................
!
      implicit none
      character Line*200,cha
      integer LinEnd(100),LinPos(100),NFIL,iCom,iErr,i,i0,icase,NumLin
      iErr=0
      if(NFIL.ne.0)  then
         call SkipL(NFIL,iErr)
         if(iErr.ne.0) go to 10
         read (NFIL,'(a)',err=10,end=10) Line
      end if
      NumLin=0
      icase=0
!............ replace commas, if present, by spaces
      if(iCom.eq.1) then
         do i=1,200
            if(Line(i:i).eq.',') Line(i:i)=' '
         end do
      end if
!............ Loop over all characters in the Line(1:200)............
      do 5 i=1,200
         i0=i
         cha=Line(i:i)
!________ if the 1st not blank character is found
         if(cha.ne.' '.and.cha.ne.'	'.and.icase.eq.0) then
            NumLin=NumLin+1
            if(NumLin.gt.100) go to 10
            LinPos(NumLin)=i
            icase=1
         end if
!_________ if the 1st blank character is found after the substring
         if(cha.eq.' '.and.icase.eq.1) then
            LinEnd(NumLin)=i-1
            icase=0
         end if
 5    end do
      if(icase.eq.1) LinEnd(NumLin)=200
      return
 10   iErr=1
      return
    end subroutine CutStr

    subroutine SkipL(NFIL,iErr)
!.......................................................................x
!...... Program skips all "empty" string containing comments, etc.......
!...... All lines, containing **,--,==,*-, ## annd -* symbols are.......
!...... recognized as the comments lines and are skipped. The    .......
!...... driver is set up at the beginning of the first "nonempty".......
!...... record                                                   .......
!.......................................................................
!  iErr=1 in the case of error input; oterwise iErr=0
!.................................................................
!
      implicit none
      integer, parameter :: nComm=6
      character Line*200,Comm(nComm)*2
      integer NFIL,iErr,i,i0
      data Comm/'*-','-*','--','==','**','##'/

!.......... read the current line from the driver NFIL
1     read (NFIL,'(a)',err=10,end=10) Line
!_________ analyze Line if it is empty
      if(Line(1:1).eq.CHAR(0)) go to 1
!_________ analyze Line if it is filled by simple blanks
      do i=1,200
         if(Line(i:i).ne.CHAR(0).and.Line(i:i).ne.' ') go to 2
      end do
      go to 1
!_________ analyze Line for the comment symbols; in the case if any
!          one from the list Comm() was found it reads the next line
 2    do i=1,nComm
         i0=INDEX(Line,Comm(i))
         if(i0.ne.0) go to 1
      end do
      iErr=0
!_________ return to the beginning of the "nonempty" line (record)
      backspace NFIL
      return
!......... in the case of an error
10    iErr=1
    end subroutine SkipL

    subroutine find_string(string,l,line,NFIL,speak,iErr)
!................................................................
! finds a string(1:l) in a file NFIL and returns with iErr=0
! and the whole Line*100 from the input. Otherwise, iErr=1
!.................................................................
      implicit none
      character*(*) string,line*200
      logical speak
      integer l,NFIL,iErr,i0
!      rewind(NFIL)
      iErr=0
10    read(NFIL,'(a)',err=20,end=20) line
      i0=index(line,string(1:l))
      if(i0.ne.0) then
         if(speak) then
            write(*,'(a)')'... string:: ['//string(1:l)//'] found'
!            write(*,'(a)') line
         end if
         return
      end if
      go to 10
20    iErr=1
    end subroutine find_string

    subroutine find_2strings(string1,l1,string2,l2,line,NFIL,speak,iErr)
!.................................................................
! finds 2 strings  string1(1:l1) and string2(1:l2) in the file NFIL
! and returns iErr=0 and the Line*100 itself from the input.
! Otherwise, iErr=1.
!.................................................................
      implicit none
      character*(*) string1,string2,line*100
      logical speak
      integer l1,l2,NFIL,iErr,i1,i2
      iErr=0
!      rewind(NFIL)
 10   read(NFIL,'(a)',err=20,end=20) line
      i1=index(line,string1(1:l1))
      i2=index(line,string2(1:l2))
      if(i1.ne.0 .and. i2.ne.0) then
         if(speak) then
            write(*,'(a)')'...2 string:: '//string1(1:l1)//' found'
            write(*,'(a)')'...2 string:: '//string2(1:l2)//' found'
         end if
         return
      end if
      go to 10
 20   iErr=1
      return
    end subroutine find_2strings

    subroutine find_3strings(string1,l1,string2,l2,string3,l3,line,NFIL,speak,iErr)
!.................................................................
! finds 3 strings  string1(1:l1) - string3(1:l3) in the file NFIL 
! and returns iErr=0 and the Line*100 itself from the input. 
! Otherwise, iErr=1.
!.................................................................
      implicit none
      character*(*) string1,string2,string3,line*100
      logical speak
      integer l1,l2,l3,NFIL,iErr,i1,i2,i3
!      rewind(NFIL)
      iErr=0
 10   read(NFIL,'(a)',err=20,end=20) line
      i1=index(line,string1(1:l1))
      i2=index(line,string2(1:l2))
      i3=index(line,string3(1:l3))
      if(i1.ne.0 .and. i2.ne.0 .and.i3.ne.0) then
         if(speak) then
            write(*,'(a)')'...3 string:: '//string1(1:l1)//' found'
            write(*,'(a)')'...3 string:: '//string2(1:l2)//' found'
            write(*,'(a)')'...3 string:: '//string3(1:l3)//' found'
         end if
         return
      end if
      go to 10
 20   iErr=1
      return
    end subroutine find_3strings

subroutine position_quote(line,pos,NFIL,iErr)
!............................................................
! reads the next line and finds the two quotes "". 
!  pos(1) - starting position just after the 1st quote "
!  pos(2) - ending position before the 2nd "
!  iErr=1 if quote was not found
!............................................................
  implicit none
  character*(*) line*200
  integer LinEnd(100),LinPos(100),NFIL,NumLin,iErr,i,j,pos(2)
   call CutStr(Line,NumLin,LinPos,LinEnd,NFIL,0,iErr)
   j=0
   do i=1,LinEnd(NumLin)
      if(line(i:i).eq.'"') then
         j=j+1 ; pos(j)=i ; if(j.eq.2) exit
      end if
   end do
   pos(1)=pos(1)+1 ; pos(2)=pos(2)-1 
   return
20 iErr=1
end subroutine position_quote

subroutine print_long_char(list,len,long_list)
!....................................................................
! Prints nicely a list (of length len) indicating if it has been cut.
!....................................................................
implicit none
character(*) list
logical long_list
integer i,len,i2
  do i=1,len,100
     i2=min(i+99,len)
     if(i2.eq.len .and. long_list) then
        write(*,'(a)') list(i:i2)//' ... (cut)'
     else
        write(*,'(a)') list(i:i2)
     end if
  end do
end subroutine print_long_char
