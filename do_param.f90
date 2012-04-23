subroutine do_param(Err)
use param
use atoms
use code
use kpoints
implicit none
!.................................................................
!  reads OUTCAR file and creates all parameters from module PARAM
!.................................................................
character Line*200,cha
integer LinEnd(100),LinPos(100),i0
logical                             :: found=.false.,Err1=.true.,Err
integer                             :: isp,NumLin,iErr,i,NKP,jNB,jNB1,NB,k
real*8, dimension(:,:), allocatable :: E
!
   Err=.false.
!
!...... read in all the parameters from the VASP output OUTCAR
!
   write(*,*)'======================================================'
   write(*,*)'===== VASP output files found (my guess) =============' 
   write(*,*)'======================================================'
   call system('ls *OUTCAR*')
   write(*,*)'======================================================'
   do
      write(*,*)' Give the name of the VASP output file:'
      read(*,'(a)') nameout
!
      open(1,file=trim(nameout),form='formatted',status='old',err=1)
      write(*,'(a)') '# Analyzing VASP file <'//trim(nameout)//'>...'
      exit
1     write(*,'(a)')'ERROR openining/reading file '//trim(nameout)//' !'
   end do
!
!___________ determine the number of species first
!
   call find_string('INCAR:',6,line,1,.true.,iErr)
   if(iErr.eq.1) then
      write(*,*)'FATAL! Cannot find pseudopotential information'
      go to 100
   end if
   isp=0 
   do 
      read(1,'(a)') line
      if(index(line,'POTCAR:').ne.0) then
         isp=isp+1
      else
         exit
      end if
   end do
   NSPEC=isp-1
   write(*,*)'found NSPEC = ',NSPEC
!
!___________ count the number of k-points and bands
!
   call find_string(' NKPTS = ',9,line,1,.true.,iErr)
!20 read(1,'(a)') line
!   if(index(line,'k-Points           NKPTS =').eq.0) go to 20
   call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(line(LinPos(4):LinEnd(4)),*) NKPTS
   write(*,*)'found NKPTS = ',NKPTS
!___________________ different format in vasp version 5.2: the NBANDS appears
!                    as word 15; in earlier versions it is word 9; the lines
!                    below account for that (28.05.2010, LK)
   write(*,'(a)') line
   read(line(LinPos(9):LinEnd(9)),*,err=18) NBANDS
   go to 19
18 read(line(LinPos(15):LinEnd(15)),*) NBANDS
19 write(*,*)'found NBANDS = ',NBANDS
!
!___________ count the number of ions
!
   call find_string(' NIONS = ',9,line,1,.true.,iErr)
!22 read(1,'(a)') line
!   if(index(line,'number of ions     NIONS =').eq.0) go to 22
   call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(line(LinPos(12):LinEnd(12)),*) NIONS
   write(*,*)'found NIONS = ',NIONS
!
!__________ check for the charge density grid used
!
25 read(1,'(a)') line
   if(index(line,'dimension x,y,z NGXF=').eq.0) go to 25
   call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(line(LinPos(4):LinEnd(4)),*) NGX
   write(*,*)'found NGX = ',NGX
   read(line(LinPos(6):LinEnd(6)),*) NGY
   write(*,*)'found NGY = ',NGY
   read(line(LinPos(8):LinEnd(8)),*) NGZ
   write(*,*)'found NGZ = ',NGZ
!
!___________ check if this is a spin-polarised run
!
10 read(1,'(a)') line
   if(index(line,'ISPIN').eq.0) go to 10
   call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(line(LinPos(3):LinEnd(3)),*) ISPIN
   if(ispin.eq.1) then
      write(*,*)'This is NOT spin-polarised calculation'
   else if(ISPIN.eq.2) then
      write(*,*)'This IS spin-polarised calculation'
   else
      go to 10
   end if
!
!.................. read eigenvalues from OUTCAR
!
   allocate(E(NKPTS,NBANDS))
   do isp=1,ISPIN
      if(isp.eq.1) then
         if(ISPIN.eq.1) then
            open(3,file='band.out',form='formatted')
         else
            open(3,file='band.out.1',form='formatted')
         end if
         write(*,'(/a)')'# Producing file <band.out.1> ...'
      else
         open(3,file='band.out.2',form='formatted')
         write(*,'(/a)')'# Producing file <band.out.2> ...'
      end if
30    read(1,'(a)',end=37) line
      if(index(line,'k-point   1 :').eq.0) go to 30
      read(1,'(a)',end=37) line
      do NKP=1,NKPTS
         do NB=1,NBANDS
            read(1,*,err=37,end=37) i,E(NKP,NB)
            if(i.ne.NB) go to 100
         end do
         read(1,'(a/a/a)',end=37) line,line,line
      end do
!____________________write energies
      do NKP=1,NKPTS
         do jNB=1,NBANDS,5
            jNB1=jNB+4
            if(jNB1.gt.NBANDS) jNB1=NBANDS
            write(3,'(i10,3x,5(1x,f10.5))') NKP,(E(NKP,NB),NB=jNB,jNB1)
         end do
      end do
   end do
   close (3)
   write(*,'(a)')'# Done band.out[1,2] !'
   yesen=.true.
   go to 38
!
37 write(*,'(/4x,a)') 'WARNING! Eigenvalues are not available in '//trim(nameout)//' file!'
   write(*,'(a)') '=========================================================='
   write(*,'(a)')  ' ===>         DOS option will not be available!       <==='
   write(*,'(a/)') '=========================================================='
38 deallocate(E)
!
!...................... now determine species themselves ...............
!
   allocate(Species(NSPEC))
   allocate(Mendel(NSPEC))
   rewind(1)
   do
      read(1,'(a)') line
      if(index(line,'INCAR:').ne.0) exit
   end do
   do isp=1,NSPEC

      call CutStr(line,NumLin,LinPos,LinEnd,1,0,iErr)
      i0=index(Line(LinPos(3):LinEnd(3)),'_')
      if(i0.ne.0) then
         LinEnd(3)=LinPos(3)+i0-2
         write(*,*)'WARNING: symbols before underscore taken!'
      end if
      i=LinEnd(3)-LinPos(3)+1
      if(i.eq.1) then
         Species(isp)=' '//Line(LinPos(3):LinEnd(3))
      else if(i.eq.2) then 
         Species(isp)=Line(LinPos(3):LinEnd(3))
      else
         Species(isp)=Line(LinPos(3):LinPos(3)+1)
         write(*,*)'WARNING: found species '//Line(LinPos(3):LinEnd(3))//&
                   ' containing more than 2 characters; two 1st characters taken'
      end if

!___________ checking from NAZV
      call check_Mendel(Species(isp),k,Err1,.true.)
      if(Err1) then
         write(*,*)'WARNING! Species '//Species(isp)//' not recognised!'
         write(*,*)'         Attempting to proceed ...'
      end if
      mendel(isp)=k

   end do
   write(*,*) 'Species found are: ',(Species(isp),' ',isp=1,NSPEC)
   close (1)
   return
!
!____________ errors
!
100 write(*,'(a)')'ERROR reading file '//trim(nameout)//' !'
200 Err=.true. ; deallocate(E) ; close(1)
    return
end subroutine do_param

