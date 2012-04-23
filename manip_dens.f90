subroutine manipulate()
!...............................................................
!   Read charge two densities and manipulate them using the 
!   following recipe:
!
!       final(grid) = c1*A(grid) + c2*B(grid)
!
!   Write result to another density file.
!...............................................................
!use param
use code
implicit none
character filenA*100,filenB*100,filen*101
character file*24,item*2
real*8 :: c1=1.0,c2=-1.0,totdens,totdensA,totdensB
integer  iErr,NPLWV,NGXa,NGYa,NGZa,NGXb,NGYb,NGZb,NPLWVa,NPLWVb
real*8  TotA,TotB
logical Yes_Do,goodA,goodB
!
!......... defaults for file names
!
      if(Which_Code.eq.'  VASP') then
         file='CHGCAR'
      else if(Which_Code.eq.'SIESTA') then
         file='siesta.RHO'
      end if
      filenA=trim(file)//'_a'
      filenB=trim(file)//'_b'
      filen =trim(file)//'_ab'
!
!......... ask what to do:
!
      Yes_Do=.false.
      goodA=.false.
      goodB=.false.
 10   write(*,*)' ....... DENSITIES MANIPULATION MENU .................'
      write(*,*)' ... CALCULATE the density from A,B-densities as: ....'
      if(c2.lt.0.0) then
         write(*,'(/a,f8.3,a,f9.3,a)') &
           '    final(grid) = ',c1,' * A(grid) ',c2,' * B(grid)'
      else
         write(*,'(/a,f8.3,a,f9.3,a)') &
           '    final(grid) = ',c1,' * A(grid) + ',c2,' * B(grid)'
      end if
      write(*,*)
      write(*,*)'   <<< active code = '//Which_Code//' >>>'
      write(*,*)' A. file name of the A-density = ',trim(filenA)
      write(*,*)'rA. read the file [',trim(filenA),'] <= only CHECK'
      if(goodA) then
         write(*,'(6x,3(a,i5))')'Grid => ',NGXa,' x ',NGYa,' x ',NGZa
         NPLWVa=NGXa*NGYa*NGZa
         write(*,'(6x,a,f10.5,a)') '...> Total charge in ['//trim(filenA)// &
              '] = ',TotA/NPLWVa,' <...'
      end if
      write(*,*)' B. file name of the B-density = ',trim(filenB)
      write(*,*)'rB. read the file [',trim(filenB),'] <= only CHECK'
      if(goodB) then
         write(*,'(6x,3(a,i5))')'Grid => ',NGXb,' x ',NGYb,' x ',NGZb
         NPLWVb=NGXb*NGYb*NGZb
         write(*,'(6x,a,f10.5,a)') '...> Total charge in ['//trim(filenB)// &
              '] = ',TotB/NPLWVb,' <...'
      end if
      write(*,*)'AB. file name of the final density = ',trim(filen)
      write(*,*)'C1. coefficient c1= ',c1
      write(*,*)'C2. coefficient c2= ',c2
      if(Yes_Do) then
         write(*,*)' F. calculate and write the final density <= DONE'
         write(*,'(a,f10.5,a)') '...> Total charge in ['//trim(filenA)// &
              '] = ',totdensA/NPLWV,' <...'
         write(*,'(a,f10.5,a)') '...> Total charge in ['//trim(filenB)// &
              '] = ',totdensB/NPLWV,' <...'
         write(*,'(a,f10.5,a)') '...> Final charge in ['//trim(filen)// &
              '] = ',totdens/NPLWV,' <...'
      else
         write(*,*)' F. calculate and write the final density'
      end if
      write(*,*)' Q. quit'
      write(*,*)
      write(*,*)'------------>'
      read(*,'(a)',ERR=11) item

![A]..... A-density filename
      if(trim(item).eq.'A') then
 13      write(*,*)'Give the file name for the A-density:'
         read(*,'(a)',err=13) filenA

![rA]..... read A-density to check the file
      else if(trim(item).eq.'rA') then
         if(Which_Code.eq.'  VASP') then
            call check_dens_vasp(19,filenA,NGXa,NGYa,NGZa,TotA,iErr)
         else if(Which_Code.eq.'SIESTA') then
            call check_dens_siesta(19,filenA,NGXa,NGYa,NGZa,TotA,iErr)
         end if
         if(iErr.eq.0) goodA=.true.

![B]..... B-density filename
      else if(trim(item).eq.'B') then
 14      write(*,*)'Give the file name for the B-density:'
         read(*,'(a)',err=14) filenB

![rB]..... read B-density to check the file
      else if(trim(item).eq.'rB') then
         if(Which_Code.eq.'  VASP') then
            call check_dens_vasp(19,filenB,NGXb,NGYb,NGZb,TotB,iErr)
         else if(Which_Code.eq.'SIESTA') then
            call check_dens_siesta(19,filenB,NGXb,NGYb,NGZb,TotB,iErr)
         end if
         if(iErr.eq.0) goodB=.true.

![AB]..... final density filename
      else if(trim(item).eq.'AB') then
 15      write(*,*)'Give the file name for the final density:'
         read(*,'(a)',err=15) filen

![C1]..... c1
      else if(trim(item).eq.'C1') then
 16      write(*,*)'Give c1:'
         read(*,*,err=16) c1

![C2]..... c2
      else if(trim(item).eq.'C2') then
 17      write(*,*)'Give c2:'
         read(*,*,err=17) c2

![F]..... calculate + write the final density
!  The total density 'totdens' is calculated in each case to
!  take care of the calculation
!
      else if(trim(item).eq.'F') then

!===================== open files for VASP

         if(Which_Code.eq.'  VASP') then

            write(*,*)'Reading in density from '//trim(filenA)//' ...'
            open(1,file=trim(filenA),status='old',form='formatted', &
                 err=150)
            write(*,*)'Reading in density from '//trim(filenB)//' ...'
            open(2,file=trim(filenB),status='old',form='formatted', &
                 err=151)
            write(*,*)'Writing density to '//trim(filen)//' ...'
            open(3,file=trim(filen),status='unknown',form='formatted')

            call vasp_dens_manip(1,2,3,c1,c2,totdensA,totdensB,totdens,NPLWV,iErr)
            close (1)
            close (2)
            if(iErr.eq.2) go to 153
            if(iErr.eq.1) go to 172
            close (3)

!===================== open files for SIESTA

         else if(Which_Code.eq.'SIESTA') then
            call siesta_dens_manip(1,2,3,c1,c2,filenA, &
                 filenB,filen,totdensA,totdensB,totdens,NPLWV,iErr)
            close (1)
            close (2)
            if(iErr.eq.2) go to 153
            if(iErr.eq.1) go to 172
            close (3)
         end if
         Yes_Do=.true.

![7]..... quit
      else if(trim(item).eq.'Q') then
         return
      else
         go to 11
      end if
      go to 10
 11   write(*,*) "Incorrect menu option! Try again!"
      go to 10

!..... errors

 150  write(*,*)'FATAL! Error while opening '//trim(filenA)//' file!'
      go to 10
 151  write(*,*)'FATAL! Error while opening '//trim(filenB)//' file!'
      close (1)
      go to 10
 153  write(*,*) 'FATAL! Files '//trim(filenA)//' and '// &
           trim(filenB)//' are inconsistent'
      go to 172
 172  close (3,status='delete')
      go to 10
end subroutine manipulate

subroutine check_dens_vasp(unit1,filen,NGX1,NGY1,NGZ1,Tot,iErr)
implicit none
integer, parameter ::  NCol0=10
integer iErr,NGX1,NGY1,NGZ1,unit1,i,j,ijk5,ijk,j0,ix,iy,iz,last,Ncol,jj
real*8 Tot,DIR1(3,3),factor,x,y,z,dens(NCol0)
character filen*100,line*200,title*10
integer LinEnd(100),LinPos(100),NSPEC1,NPLWV
integer ix5(NCol0),iy5(NCol0),iz5(NCol0)
integer, dimension(:), allocatable :: NatSpec1

  iErr=0
  ncol=5
  open(unit1,file=trim(filen),status='old',form='formatted',err=10)
!
!.......... going through the heading
!
  write(*,*) '_______> (A) density file: <________'
  
  read(unit1,'(a)',err=10,end=10) title
  read(unit1,*,err=10,end=10) factor
  do i=1,3
     read(unit1,*,err=10,end=10) (DIR1(i,j), j=1,3)
  end do
  write(*,'(a)')'(A) ... found lattice vectors:'
  write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR1(1,j),j=1,3)
  write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR1(2,j),j=1,3)
  write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR1(3,j),j=1,3)
  call CutStr(Line,NSPEC1,LinPos,LinEnd,unit1,0,iErr)
  allocate(NatSpec1(NSPEC1))
  j=0
  do i=1,NSPEC1
     read(Line(LinPos(i):LinEnd(i)),*,err=10) NatSpec1(i)
     j=j+NatSpec1(i)
  end do
  write(*,*) '...  found number of species = ',NSPEC1, ' with species numbers: '
  write(*,'(20(i4,x))') (NatSpec1(i),i=1,NSPEC1)
  write(*,*) '...  found total number of atoms = ',j
!
  read(unit1,'(a)',err=10,end=10) title
!
!_____________Skipping coordinates

  do i=1,NSPEC1
     do j=1,NatSpec1(i)
        read(unit1,*,err=10,end=10) x,y,z
      end do
  end do
 
!.......... grid
  read(unit1,*,err=10,end=10) NGX1,NGY1,NGZ1
  NPLWV=NGX1*NGY1*NGZ1

!.......... reading density and calculating the total charge Tot

  ijk5=0
  ijk=0
  j0=NPLWV/10
  Tot=0.0
  do iz=1,NGZ1
     do iy=1,NGY1
        do ix=1,NGX1
           ijk=ijk+1
           ijk5=ijk5+1
           last=NCol
           if(ijk.eq.NPLWV) last=ijk5
           if(ijk5.le.last) then
              ix5(ijk5)=ix
              iy5(ijk5)=iy
              iz5(ijk5)=iz
           end if
           if(ijk5.eq.last) then
              read(unit1,*,err=10,end=10) (dens(jj),jj=1,last)
              do jj=1,last
                 Tot=Tot+dens(jj)
              end do
              ijk5=0
           end if
           if(ijk/j0*j0.eq.ijk) &
                write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
        end do
     end do
  end do
  write(*,*)'Done!'
  go to 20


!........... error messages
10 write(*,*)'FATAL! The density file is bad!'
  iErr=1
  go to 20

!............ finish

20 deallocate(NatSpec1)
   close (unit1)
end subroutine check_dens_vasp

subroutine check_dens_siesta(unit1,filen,NGX1,NGY1,NGZ1,Tot,iErr)
implicit none
real*8,parameter :: A_au=1.889725989,au_A=1.0/A_au
integer unit1,iErr,i,j,NGX1,NGY1,NGZ1,nspins
integer j0,ijk,iz,iy,ix,NPLWV
real*8 DIR1(3,3),RECC(3,3),VOLC,Tot,factor
real  ,dimension(:),allocatable :: dens1
character filen*100
logical formA, alloc

   iErr=0
   alloc=.false.

!................. open density file

   write(*,*) 'Trying reading density as FORMATTED ...'
   open(unit1,file=trim(filen),status='old',form='formatted',err=1)
   do i=1,3
      read(unit1,*,err=1,end=1) (DIR1(i,j),j=1,3)
   end do
   formA=.true.
   write(*,'(a)')'Opened '//trim(filen)//' as FORMATTED'
   go to 5
   
1  close (unit1)
   write(*,*) 'Trying reading density as UNFORMATTED ...'
   open(unit1,file=trim(filen),status='old',form='unformatted',err=10)
   read(unit1,err=10,end=10) ((DIR1(i,j),j=1,3),i=1,3)
   formA=.false.
   write(*,'(a)') 'Opened '//trim(filen)//' as UNFORMATTED'
   
5  write(*,'(a)')'... found lattice vectors:'
   write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR1(1,j)*au_A,j=1,3)
   write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR1(2,j)*au_A,j=1,3)
   write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR1(3,j)*au_A,j=1,3)
!
   if(formA) then
      read(unit1,*,err=10,end=10) NGX1,NGY1,NGZ1,nspins
   else
      read(unit1,err=10,end=10) NGX1,NGY1,NGZ1,nspins
   end if
   write(*,'(a,3(i5,x))')  '... found grid: ',NGX1,NGY1,NGZ1
   write(*,'(a,i2)')       '... found spin: ',nspins
!
   allocate(dens1(NGX1))
   alloc=.true.
   
   DIR1=DIR1*au_A
   call bastr(DIR1,RECC,VOLC,1)
   factor=VOLC*A_au**3
   Tot = 0.0

   ijk=0
   j0=NGX1*NGY1*NGZ1/10
   do iz=1,NGZ1
      do iy=1,NGY1
         ijk=ijk+NGX1
         
         if(formA) then
            read(unit1,*,err=10,end=10) (dens1(ix),ix=1,NGX1)
         else
            read(unit1,err=10,end=10) (dens1(ix),ix=1,NGX1)
         end if
         
         do ix=1,NGX1
            Tot=Tot+dens1(ix)*factor
         end do
         
         if(ijk/j0*j0.eq.ijk) &
              write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
      end do
   end do
   write(*,*)'Done!'
   go to 20
!
!.......errors
10 write(*,*)'FATAL! File '//trim(filen)//' not found or bad!'
   iErr=1
   go to 20
   
!............ finish
   
20 if(alloc) deallocate(dens1)
   close (unit1)
end subroutine check_dens_siesta
