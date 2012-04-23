subroutine vasp_dens(grid,charge,name,spin,iErr)
!..............................................................
!    read density from name=CHGCAR/LOCPOT/... file of VASP
!
!    ncol=5  in the case of CHGCAR (density)
!    ncol=5  in the case of LOCPOT (potential)
!    ncol=10 in the case of PARCHG (partial density)
!..............................................................
!    spin =.false. - total charge \____ for all files
!         = .true. - spin density /
!..............................................................
use param
use atoms
implicit none
integer, parameter :: NCol0=10
real*8, parameter  :: tiny=0.0001
integer, dimension(:), allocatable :: NatSpec 
real*8 GRID(NGX,NGY,NGZ),DIR(3,3),charge,factor,x,y,z,dens(NCol0)
integer ix5(NCol0),iy5(NCol0),iz5(NCol0),iErr,ncol,i,j
!integer icase
integer NGX1,NGY1,NGZ1,ijk,ijk5,iz,iy,ix,last,jj,i0,j0
character name*6,title*10,line*100
logical spin
!
!......... which case
!
      allocate(NatSpec(NSPEC))
      iErr=0
      if(name.eq.'CHGCAR' .or. name.eq.'LOCPOT') then
!         icase=1
         ncol=5
      else if(name.eq.'PARCHG') then
!         icase=2
         ncol=10
      else
         write(*,*)'Not yet implemented!'
         iErr=1
         go to 300
      end if
      open(19,file=name,status='old',form='formatted',err=100)
!
!.......... skipping stuff and checking
!
      read(19,'(a)',err=100,end=100) title
      read(19,*,err=100,end=100) factor
      read(19,*,err=100,end=100) DIR(1,1), DIR(1,2), DIR(1,3)
      read(19,*,err=100,end=100) DIR(2,1), DIR(2,2), DIR(2,3)
      read(19,*,err=100,end=100) DIR(3,1), DIR(3,2), DIR(3,3)
      write(*,*) DIR(1,1), DIR(1,2), DIR(1,3)
      write(*,*) DIR(2,1), DIR(2,2), DIR(2,3)
      write(*,*) DIR(3,1), DIR(3,2), DIR(3,3)
      do i=1,3
         do j=1,3
            if(abs(DIRC(i,j)-factor*DIR(i,j)).gt.tiny) then
               write(*,*)'FATAL! Found different lattice vectors!'
               go to 100
            end if
         end do
      end do
!______________ attempt to read other lines; in version 5.2 we have a new line
!               with atomic species; read it if possible, otherwise skip for
!               older versions (29.05.2010, LK)
      read(19,*,err=10) (NatSpec(i),i=1,NSPEC)
      go to 15
10    read(19,*,err=100,end=100) (NatSpec(i),i=1,NSPEC)
15    write(*,*) (NatSpec(i),i=1,NSPEC)
      do i=1,NSPEC
         if(NatSpec(i).ne.NspN(i)) then
            write(*,*)'FATAL! Found different species numbers!'
            go to 100
         end if
      end do
!      write(*,*) (NatSpec(i),i=1,NSPEC)
      read(19,'(a)',err=100,end=100) title
      write(*,'(a)') title
!_____________Skipping coordinates: not checking. It could be worth doing
!             though...
      do i=1,NSPEC
         do j=1,NspN(i)
            read(19,*,err=100,end=100) x,y,z
         end do
      end do
!_____________checking grid
      read(19,*,err=100,end=100) NGX1,NGY1,NGZ1
      write(*,*)  NGX1,NGY1,NGZ1
      if(NGX.ne.NGX1 .or.NGY.ne.NGY1 .or.NGZ.ne.NGZ1) then
         write(*,*)'FATAL! Found different grid dimensions!'
         go to 100
      end if
!
!..........reading the charge density
!     It was written in NCol=5 or 10 columns using something like:
!
!     WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      if(spin) then
         write(*,*)'Skipping the charge density, please wait ...'
      else
         if(name.eq.'CHGCAR') then
            write(*,*)'Reading in charge density, please wait ...'
         else if(name.eq.'PARCHG') then
            write(*,*)'Reading in partial density, please wait ...'
         else if(name.eq.'LOCPOT') then
            write(*,*)'Reading in potential, please wait ...'
         end if
      end if
!
      ijk5=0
      ijk=0
      j0=NPLWV/10
      charge=0.0
      do iz=1,NGZ
         do iy=1,NGY
            do ix=1,NGX
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
                  read(19,*,err=100,end=100) (dens(jj),jj=1,last)
                  do jj=1,last
!                    if(dens(jj).lt.0.0) dens(jj)=0.0
                    if(.not.spin) grid(ix5(jj),iy5(jj),iz5(jj))=dens(jj)
                    charge=charge+dens(jj)
                  end do
                  ijk5=0
               end if
               if(ijk/j0*j0.eq.ijk) &
                    write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
            end do
         end do
      end do
      write(*,*)'Done!'
      write(*,'(a,e12.6)')'The charge is = ',charge/NPLWV
      if(.not.spin) go to 70
!
!.......... reading the spin density if spin=.true.
!     It was written in 5 columns using something like:
!
!  WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      write(*,*)'... jumping to spin density'
 51   read(19,'(a)',end=200) line
!      write(*,'(a)') line
      i0=index(line,'.')
      if(i0.ne.0) go to 51
      read(line,*,err=51,end=51) NGX1,NGY1,NGZ1
      write(*,*) NGX1,NGY1,NGZ1
      if(NGX.ne.NGX1 .or.NGY.ne.NGY1 .or.NGZ.ne.NGZ1) go to 51
      write(*,*)'Reading in spin density, please wait ...'
!
      ijk5=0
      ijk=0
      j0=NPLWV/10
      charge=0.0
      do iz=1,NGZ
         do iy=1,NGY
            do ix=1,NGX
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
                  read(19,*,err=200,end=200) (dens(jj),jj=1,last)
                  do jj=1,last
                    grid(ix5(jj),iy5(jj),iz5(jj))=dens(jj)
                    charge=charge+dens(jj)
                  end do
                  ijk5=0
               end if
               if(ijk/j0*j0.eq.ijk) &
                    write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
            end do
         end do
      end do
      write(*,*)'Done!'
      write(*,'(a,e12.6)')'The spin charge is = ',charge/NPLWV
 70   close (19)
      go to 300
!
!.......errors
 100  write(*,*)'FATAL! File '//name//' does not exist, bad or wrong!'
      iErr=1
      go to 300
 200  write(*,*)'WARNING! There is no spin density found!'
      iErr=1

!......... finish
300   deallocate(NatSpec)
      return
end subroutine vasp_dens

subroutine siesta_dens(grid,charge,spin,iErr)
!..............................................................
!    read density from name=seed.RHO... file of siesta
!..............................................................
!    spin = .false. - total charge
!         = .true.  - spin density
!..............................................................
use param
use atoms
use code
implicit none
real*8,parameter :: A_au=1.889725989,au_A=1.0/A_au
real*8,parameter :: tiny=0.001
real*8 GRID(NGX,NGY,NGZ),DIR(3,3)
real*8 charge,factor
logical spin,form
integer iErr,i,j,NGX1,NGY1,NGZ1,nspins,ijk,j0,iz,iy,ix
!.......... should be real, not real*8, otherwise it gives an error
!           when reading an unformatted density file
real dens(NGX)

      iErr=0
      write(*,'(a)') &
              'Trying '//trim(seed)//'.RHO as FORMATTED ...'
      open(19,file=trim(seed)//'.RHO',status='old',form='formatted', &
                                    err=10)
      write(*,'(a)') 'Opened '//trim(seed)//'.RHO as FORMATTED'
      do i=1,3
         read(19,*,err=10,end=10) (DIR(i,j),j=1,3)
      end do
      write(*,'(a)')'Reading lattice vectors done'
      form=.true.
      go to 20

 10   close (19)
      write(*,'(a)') 'Trying '//trim(seed)//'.RHO as UNFORMATTED ...'
      open(19,file=trim(seed)//'.RHO',status='old', &
                                           form='unformatted',err=100)
      write(*,'(a)')'Opened '//trim(seed)//'.RHO as UNFORMATTED'
      read(19,err=100,end=100) ((DIR(i,j),j=1,3),i=1,3)
      write(*,'(a)')'Reading lattice vectors done'
      form=.false.
!
!.......... skipping stuff and checking
!
 20   factor=VOLC*A_au**3
      write(*,'(a)')'..........in the density .........|........from input...............'
      do j=1,3
         write(*,'(3(f10.5,x),a,3(f10.5,x))') (DIR(j,i)*au_A, i=1,3),' | ',(DIRC(j,i), i=1,3)
      end do
      do i=1,3
         do j=1,3
            if(abs(DIRC(i,j)-au_A*DIR(i,j)).gt.tiny) then
               write(*,*) 'ERROR! Mismatch of lattice vectors!'
               write(*,'(2i2,3(x,e12.6))') i,j,DIRC(i,j),au_A*DIR(i,j),abs(DIRC(i,j)-au_A*DIR(i,j))
               go to 100
            end if
         end do
      end do

!_____________checking grid and spin

      if(form) then
         read(19,*,err=100,end=100) NGX1,NGY1,NGZ1,nspins
      else
         read(19,err=100,end=100) NGX1,NGY1,NGZ1,nspins
      end if
      write(*,*)  NGX1,NGY1,NGZ1,nspins
      if(NGX.ne.NGX1 .or.NGY.ne.NGY1 .or.NGZ.ne.NGZ1) go to 100
      if(spin .and. nspins.eq.1) go to 200

!
!..........reading the charge density
!
!     WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      if(spin) then
         write(*,*)'Skipping the charge density, please wait ...'
      else
         write(*,*)'Reading in charge density, please wait ...'
      end if
!
      charge=0.0
      ijk=0
      j0=NPLWV/10
      do iz=1,NGZ
         do iy=1,NGY
            ijk=ijk+NGX
            if(form) then
               read(19,*,err=100,end=100) (dens(ix),ix=1,NGX)
            else
               read(19,err=100,end=100) (dens(ix),ix=1,NGX)
               write(111,'(a,x,2(i5,x),5(f10.5,x))') '19 ',iz,iy,dens(1:5)
            end if
            do ix=1,NGX
               if(.not.spin) then
                  grid(ix,iy,iz)=dens(ix)*factor
                  charge=charge+grid(ix,iy,iz)
               else
                  charge=charge+dens(ix)*factor
               end if
            end do
            if(ijk/j0*j0.eq.ijk) &
                 write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
         end do
      end do
      write(*,*)'Done!'
      write(*,'(a,e12.6)')'The charge is = ',charge/NPLWV
      if(.not.spin) go to 70
!
!.......... reading the spin density if spin=.true.
!
!  WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      write(*,*)'... jumping to spin density'
      write(*,*)'Reading in spin density, please wait ...'
!
      charge=0.0
      ijk=0
      j0=NPLWV/10
      do iz=1,NGZ
         do iy=1,NGY
            ijk=ijk+NGX
            if(form) then
               read(19,*,err=100,end=100) (dens(ix),ix=1,NGX)
            else
               read(19,err=100,end=100) (dens(ix),ix=1,NGX)
            end if
            do ix=1,NGX
               grid(ix,iy,iz)=dens(ix)*factor
               charge=charge+grid(ix,iy,iz)
            end do
            if(ijk/j0*j0.eq.ijk) &
                 write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
         end do
      end do
      write(*,*)'Done!'
      write(*,'(a,e12.6)')'The spin charge is = ',charge/NPLWV
 70   close (19)
      return
!
!.......errors
 100  write(*,*)'FATAL! File '//trim(seed)//'.RHO not found or bad!'
      iErr=1
      return
 200  write(*,*)'WARNING! There is no spin density found!'
      iErr=1
      return
end subroutine siesta_dens

subroutine siesta_dens_manip(unit1,unit2,unit3,c1,c2, &
                      filenA,filenB,filen,totdensA,totdensB,totdens,NPLWV,iErr)
!..............................................................
!    reading density from unit1, unit2 of SIESTA and
!    manipulaiting into:
!                         density=c1*den1+c2*den2
!..............................................................
!use param
use atoms
use code
implicit none
real*8,parameter :: tiny=0.001
real*8,parameter :: A_au=1.889725989,au_A=1.0/A_au
integer unit1,unit2,unit3,iErr,i,j,NGX1,NGY1,NGZ1,nspins,NGX2,NGY2,NGZ2
integer j0,ijk,iz,iy,ix,NPLWV
real*8 DIR1(3,3),DIR2(3,3),c1,c2
real*8,dimension(:),allocatable :: dens
real  ,dimension(:),allocatable :: dens1,dens2
real*8 totdensA,totdensB,totdens,factor
character filenA*100,filenB*100,filen*101
logical formA, formB, alloc

      iErr=0
      alloc=.false.

!................. open A-density file

      write(*,*) 'Trying reading A-density as FORMATTED ...'
      open(unit1,file=trim(filenA),status='old',form='formatted',err=1)
      do i=1,3
         read(unit1,*,err=1,end=1) (DIR1(i,j),j=1,3)
      end do
      formA=.true.
      write(*,'(a)')'Opened '//trim(filenA)//' as FORMATTED'
      go to 5

 1    close (unit1)
      write(*,*) 'Trying reading A-density as UNFORMATTED ...'
      open(unit1,file=trim(filenA),status='old',form='unformatted',err=100)
      read(unit1,err=100,end=100) ((DIR1(i,j),j=1,3),i=1,3)
      formA=.false.
      write(*,'(a)') 'Opened '//trim(filenA)//' as UNFORMATTED'

5     write(*,'(a)')'(A) ... found lattice vectors:'
      write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR1(1,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR1(2,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR1(3,j),j=1,3)

!................. open B-density file

      write(*,*) 'Trying reading B-density as FORMATTED ...'
      open(unit2,file=trim(filenB),status='old',form='formatted',err=2)
      do i=1,3
         read(unit2,*,err=2,end=2) (DIR2(i,j),j=1,3)
      end do
      formB=.true.
      write(*,'(a)')'Opened '//trim(filenB)//' as FORMATTED'
      go to 6

 2    close (unit2)
      write(*,*) 'Trying reading B-density as UNFORMATTED ...'
      open(unit2,file=trim(filenB),status='old',form='unformatted',err=101)
      read(unit2,err=101,end=101) ((DIR2(i,j),j=1,3),i=1,3)
      formB=.false.
      write(*,'(a)') 'Opened '//trim(filenB)//' as UNFORMATTED'

6     write(*,'(a)')'(B) ... found lattice vectors:'
      write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR2(1,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR2(2,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR2(3,j),j=1,3)

!................. open the final file as formatted

      write(*,*)'Writing density to '//trim(filen)//' as FORMATTED ...'
      open(unit3,file=trim(filen),status='unknown',form='formatted',err=102)

!................. checking lattice vectors and the grid

!__________ checking lattice vectors

      do i=1,3
         do j=1,3
            if(abs(DIR2(i,j)-DIR1(i,j)).gt.tiny) then
               iErr=2
               return
            end if
         end do
      end do
      write(*,*)'... Lattice vectors: check passed!'

!_____________checking grid

      if(formA) then
         read(unit1,*,err=100,end=100) NGX1,NGY1,NGZ1,nspins
      else
         read(unit1,err=100,end=100) NGX1,NGY1,NGZ1,nspins
      end if
      write(*,'(a,3(i5,x))')  '(A) ... found grid: ',NGX1,NGY1,NGZ1
      write(*,'(a,i2)')  '(A) ... found spin: ',nspins
      if(formB) then
         read(unit2,*,err=101,end=101) NGX2,NGY2,NGZ2,nspins
      else
         read(unit2,err=101,end=101) NGX2,NGY2,NGZ2,nspins
      end if
      write(*,'(a,3(i5,x))')  '(B) ... found grid: ',NGX2,NGY2,NGZ2
      write(*,'(a,i2)')  '(B) ... found spin: ',nspins
      if(NGX2.ne.NGX1 .or.NGY2.ne.NGY1 .or.NGZ2.ne.NGZ1) then
         iErr=2
         return
      else
         write(*,*)'... Grid: check passed!'
         NPLWV=NGX1*NGY1*NGZ1
      end if

!................. start READING/WRITING

      do i=1,3
         write(unit3,*) (DIR1(i,j),j=1,3)
         do j=1,3
            DIR1(i,j)=DIR1(i,j)*au_A
         end do
      end do
      write(unit3,*) NGX1,NGY1,NGZ1,' 1'

      call bastr(DIR1,RECC,VOLC,1)
      factor=VOLC*A_au**3
!
!..........reading the charge density
!
!     WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      allocate(dens(NGX1))
      allocate(dens1(NGX1))
      allocate(dens2(NGX1))
      alloc=.true.

      totdensA = 0.0
      totdensB = 0.0
      totdens  = 0.0
      ijk=0
      j0=NGX1*NGY1*NGZ1/10
      do iz=1,NGZ1
         do iy=1,NGY1
            ijk=ijk+NGX1
            
            if(formA) then
               read(unit1,*,err=100,end=100) (dens1(ix),ix=1,NGX1)
            else
               read(unit1,err=100,end=100) (dens1(ix),ix=1,NGX1)
            end if
            if(formB) then
               read(unit2,*,err=101,end=101) (dens2(ix),ix=1,NGX1)
            else
               read(unit2,err=101,end=101) (dens2(ix),ix=1,NGX1)
            end if

            do ix=1,NGX1
               dens(ix)=c1*dens1(ix)+c2*dens2(ix)
               totdensA=totdensA+dens1(ix)*factor
               totdensB=totdensB+dens2(ix)*factor
               totdens =totdens + dens(ix)*factor
            end do
            write(unit3,*) (dens(ix),ix=1,NGX1)

            if(ijk/j0*j0.eq.ijk) &
                 write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
         end do
      end do
      write(*,*)'Done!'
      go to 200
!
!.......errors
 100  write(*,*)'FATAL! File '//trim(filenA)//' not found or bad!'
      iErr=1
      go to 200
 101  write(*,*)'FATAL! File '//trim(filenB)//' not found or bad!'
      iErr=1
      go to 200
 102  write(*,*)'FATAL! File '//trim(filen)//' cannlot be opened!'
      iErr=1

!............ finish

200   if(alloc) then
         deallocate(dens)
         deallocate(dens1)
         deallocate(dens2)
      end if
end subroutine siesta_dens_manip

subroutine vasp_dens_manip(unit1,unit2,unit3,c1,c2, &
                           totdensA,totdensB,totdens,NPLWV,iErr)
!.................................................................
!    reading density from unit1, unit2 of VASP and manipulaiting
!    into:
!                density=c1*den1+c2*den2
!.................................................................
!use param
use atoms
implicit none
real*8, parameter  :: pi=3.1415927d0,tiny=0.0001
integer, parameter ::  NCol0=10
integer unit1,unit2,unit3,ncol,iErr,i,j,j0,NPLWV
real*8 DIR1(3,3),DIR2(3,3),dens(NCol0),dens1(NCol0),dens2(NCol0),factor
real*8 x,y,z,totdensA,totdensB,totdens,c1,c2
integer, dimension(:), allocatable :: NatSpec1, NatSpec2
integer ix5(NCol0),iy5(NCol0),iz5(NCol0),NGX1,NGY1,NGZ1,NGX2,NGY2,NGZ2
integer ijk5,ijk,iz,iy,ix,last,jj
character title*10,line*200
logical check
integer LinEnd(100),LinPos(100),NSPEC1,NSPEC2

      iErr=0
      ncol=5
!
!.......... going through the heading
!
!..... A-file

      write(*,*) '_______> (A) density file: <________'

      read(unit1,'(a)',err=100,end=100) title
      read(unit1,*,err=100,end=100) factor
      do i=1,3
         read(unit1,*,err=100,end=100) (DIR1(i,j), j=1,3)
      end do
      write(*,'(a)')'(A) ... found lattice vectors:'
      write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR1(1,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR1(2,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR1(3,j),j=1,3)
      call CutStr(Line,NSPEC1,LinPos,LinEnd,unit1,0,iErr)
      allocate(NatSpec1(NSPEC1))
      j=0
      do i=1,NSPEC1
         read(Line(LinPos(i):LinEnd(i)),*,err=100) NatSpec1(i)
         j=j+NatSpec1(i)
      end do
      write(*,*) '(A) ...  found number of species = ',NSPEC1, ' with species numbers: '
      write(*,'(20(i4,x))') (NatSpec1(i),i=1,NSPEC1)
      write(*,*) '(A) ...  found total number of atoms = ',j

!..... B-file

      write(*,*) '_______> (B) density file: <________'
      read(unit2,'(a)',err=101,end=101) title
      read(unit2,*,err=101,end=101) factor
      do i=1,3
         read(unit2,*,err=101,end=101) (DIR2(i,j), j=1,3)
      end do
      write(*,'(a)')'(B) ... found lattice vectors:'
      write(*,'(a,3(f10.5,x))')'   DIR(1) ',(DIR2(1,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(2) ',(DIR2(2,j),j=1,3)
      write(*,'(a,3(f10.5,x))')'   DIR(3) ',(DIR2(3,j),j=1,3)
      call CutStr(Line,NSPEC2,LinPos,LinEnd,unit2,0,iErr)
      allocate(NatSpec2(NSPEC2))
      j=0
      do i=1,NSPEC2
         read(Line(LinPos(i):LinEnd(i)),*,err=100) NatSpec2(i)
         j=j+NatSpec2(i)
      end do
      write(*,*) '(B) ...  found number of species = ',NSPEC2, ' with species numbers: '
      write(*,'(20(i4,x))') (NatSpec2(i),i=1,NSPEC2)
      write(*,*) '(B) ...  found total number of atoms = ',j

!.... checks

      check=.true.
      if(NSPEC1.ne.NSPEC2) then
         check=.false.
      else
         do i=1,NSPEC1
            if(NatSpec1(i).ne.NatSpec2(i)) check=.false.
         end do
      end if

      if(.not.check) then
         write(*,*)'WARNING! Species are different in the two densities!'
         write(*,*)'WARNING! The file to be created will have information from file A'
         write(*,*)'         Please, EDIT BEFORE USE! '
         write(*,*)'         Press ENTER when ready ...'
         read(*,*)
      end if

!.......... check consistency of the two files: lattice vectors
      do i=1,3
         do j=1,3
            if(abs(DIR1(i,j)-DIR2(i,j)).gt.tiny) then
               iErr=2
               go to 200
            end if
         end do
      end do
      write(*,*)'... Lattice vectors: check passed!'
!
!.......... writing this to unit3
      write(unit3,'(a)') title
      write(unit3,*) factor
      write(unit3,'(x,3(2x,f10.6))') DIR1(1,1), DIR1(1,2), DIR1(1,3)
      write(unit3,'(x,3(2x,f10.6))') DIR1(2,1), DIR1(2,2), DIR1(2,3)
      write(unit3,'(x,3(2x,f10.6))') DIR1(3,1), DIR1(3,2), DIR1(3,3)
      write(unit3,*) (NatSpec1(i),i=1,NSPEC1)
!
      read(unit1,'(a)',err=100,end=100) title
      read(unit2,'(a)',err=101,end=101) title
      write(unit3,'(a)') title
!
!_____________Skipping coordinates: not checking; writing coordinates to unit3

      do i=1,NSPEC1
         do j=1,NatSpec1(i)
            read(unit1,*,err=100,end=100) x,y,z
            write(unit3,'(3(x,f9.6))') x,y,z
         end do
      end do
      do i=1,NSPEC2
         do j=1,NatSpec2(i)
            read(unit2,*,err=101,end=101) x,y,z
         end do
      end do

!.......... check consistency of the two files: grid
      read(unit1,*,err=100,end=100) NGX1,NGY1,NGZ1
      read(unit2,*,err=101,end=101) NGX2,NGY2,NGZ2
      if(NGX2.ne.NGX1 .or.NGY2.ne.NGY1 .or.NGZ2.ne.NGZ1) then
         iErr=2
         go to 200
      else
         write(*,*)'... Grid: check passed!'
         NPLWV=NGX1*NGY1*NGZ1
      end if
      write(unit3,*) NGX1,NGY1,NGZ1
!
!..........reading the charge density
!     It was written in NCol=5 columns using something like:
!
!     WRITE(..,..) (((DENS(I,J,K),I=1,NGX), J=1,NGY), K=1,NGZ)
!
      write(*,*)'Reading in charge densities, please wait ...'
!
      ijk5=0
      ijk=0
      j0=NPLWV/10

      totdensA = 0.0
      totdensB = 0.0
      totdens  = 0.0
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
                  read(unit1,*,err=100,end=100) (dens(jj),jj=1,last)
                  read(unit2,*,err=101,end=101) (dens1(jj),jj=1,last)
                  do jj=1,last
                    dens2(jj)=c1*dens(jj)+c2*dens1(jj)
                    totdensA=totdensA+dens(jj)
                    totdensB=totdensB+dens1(jj)
                    totdens=totdens+dens2(jj)
                  end do
                  write(unit3,'(5(e18.11,x))') (dens2(jj),jj=1,last)
                  ijk5=0
               end if
               if(ijk/j0*j0.eq.ijk) &
                    write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
            end do
         end do
      end do
      write(*,*)'Done!'
      go to 200
!
!........... error messages
 100  write(*,*)'FATAL! The A-density file is bad!'
      iErr=1
      go to 200
 101  write(*,*)'FATAL! The B-density file is bad!'
      iErr=1

!............ finish

200   deallocate(NatSpec1)
      deallocate(NatSpec2)
      return
end subroutine vasp_dens_manip

