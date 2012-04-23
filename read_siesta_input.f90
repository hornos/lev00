subroutine read_siesta_geom(Err)
!.................................................................
!  reading SIESTA output geometry file seed.XV
!.................................................................
!     It is assumed that coordinates are given in SIESTA in a.u.,
! so that we convert them back to A.
!.................................................................
use param
use code
use atoms
use kpoints
use mendeleev
implicit none
real*8,parameter :: A_au=1.889725989,au_A=1.0/A_au
character Line*200
integer k,i,ni,nat,j,nkp,i1,i2,i0
integer LinEnd(100),LinPos(100),NumLin,iErr
logical good,Err
real*8,dimension(:,:),allocatable ::  TI1
integer,dimension(:),allocatable  ::  isp,mendl
real*8 r(3),s,vscale

   good=.true.
   write(*,*) 'Reading from '//trim(seedfdf)//'  ...'
   open(14,file=trim(seedfdf),form='formatted',status='old',err=29)

!_________ reading lattice vectors
   write(*,*) 'Reading from '//trim(seedfdf)//' lattice vectors ...'
   call find_string('LatticeConstant',15,line,14,.false.,iErr)
   write(*,*) line
   if(iErr.eq.1) then
      vscale=1.0
   else
      call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
      read(Line(LinPos(2):LinEnd(2)),*,err=29,end=29) vscale
   end if
   write(*,*)'SIESTA: found LatticeConstant = ',vscale
   
   rewind(14)
   call find_string('%block LatticeVectors',21,Line,14,.false.,iErr)
   if(iErr.eq.1) go to 29
   read (14,*,err=29,end=29) DIRC(1,1), DIRC(1,2), DIRC(1,3)
   read (14,*,err=29,end=29) DIRC(2,1), DIRC(2,2), DIRC(2,3)
   read (14,*,err=29,end=29) DIRC(3,1), DIRC(3,2), DIRC(3,3)
   DIRC=DIRC*vscale
   write(*,*)'SIESTA: found lattice vectors:'
   do i=1,3
      write(*,'(a,i1,a,3(f10.5,x))') '[',i,']', (DIRC(i,j),j=1,3)
   end do
   call bastr(DIRC,RECC,VOLC,1)
   call bastr(DIRC,BCELL,VOLC,0)

!_________ reading Mendeleev's numbers and setting species
   allocate(Species(NSPEC)) ; allocate(Mendel(NSPEC))

   rewind(14)
   call find_string('%block ChemicalSpeciesLabel',27,Line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,*)'ERROR: no ChemicalSpeciesLabel found'
      go to 29
   end if
   do i=1,NSPEC
      read(14,*,err=29) i1,Mendel(i)
      Species(i)=NAZV(Mendel(i))
   end do
   write(*,'(a,20(i2,3a))')'SIESTA: found species: ',(Mendel(i),&
        '[',Species(i),'] ',i=1,NSPEC)

!_________ reading geometry as it is in the input

   allocate(isp(NIONS)) ; allocate(mendl(NIONS)) ; allocate(TI1(3,NIONS))
   allocate(old_to_new(NIONS))

!_______________________ reading format
   rewind(14)
   call find_string('AtomicCoordinatesFormat',22,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,*)'SIESTA: have not found format, implying Ang'
      i1=1
   else
      call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
      if(iErr.eq.1.or.NumLin.lt.2) go to 29
      i0=index(Line(LinPos(2):LinEnd(2)),'Cartesian')
      i1=index(Line(LinPos(2):LinEnd(2)),'Ang')
      i2=index(Line(LinPos(2):LinEnd(2)),'Bohr')
   end if
   if(i0.ne.0 .or. i1.ne.0 .or. i2.ne.0) then
      write(*,'(a)') 'SIESTA: implying Cartesian format, reading in ...'
   else
      write(*,'(a)') 'SIESTA: implying fractional format, reading in ...'
   end if

!_________________________ reading geometry 
   rewind(14)
   call find_string('%block AtomicCoordinatesAndAtomicSpecies',40,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,*)'ERROR: no atomic coordinates'
      go to 29
   end if
   do nat=1,NIONS
      if(i0.ne.0 .or. i1.ne.0 .or. i2.ne.0) then
         read(14,*,err=29,end=29) (TI(j,nat),j=1,3),isp(nat)
         if(i1.ne.0) then
            TI(1:3,nat)=TI(1:3,nat)*vscale
         else if(i2.ne.0) then
            TI(1:3,nat)=TI(1:3,nat)*vscale*au_A
         end if
      else
         read(14,*,err=29,end=29) (r(j),j=1,3),isp(nat)
         TI(1,nat)=r(1)*DIRC(1,1)+r(2)*DIRC(2,1)+r(3)*DIRC(3,1)
         TI(2,nat)=r(1)*DIRC(1,2)+r(2)*DIRC(2,2)+r(3)*DIRC(3,2)
         TI(3,nat)=r(1)*DIRC(1,3)+r(2)*DIRC(2,3)+r(3)*DIRC(3,3)
      end if
      mendl(nat)=Mendel(isp(nat))
   end do

!__________ sorting atoms by species
!           ni - counts atoms in the current species i
!            k - counts all atoms
   k=0
   do i=1,NSPEC
      ni=0
      do nat=1,NIONS
         if(isp(nat).eq.i) then
            ni=ni+1 ; k=k+1 ; TI1(1:3,k)=TI(1:3,nat)
            old_to_new(nat)=k
        end if
      end do
      NspN(i)=ni
   end do
   if(k.ne.NIONS) good=.false.
   TI=TI1
   write(*,'(a,10(i5,x))')'SIESTA: found NspN = ',(NspN(i),i=1,NSPEC)
   write(*,*)'Done geometry!'
  
!............. finishing
!
    if(.not. good) go to 29
    write(*,*) '........> Total number of atoms in the cell is ',NIONS
!
!___________ read relaxation flags
!
    relax_flag='T'
    rewind (14)
    call get_relax_flags_SIESTA(14)
    yesrlx=.true.
!
!___________ show atoms 
!
    nat=0
    do i=1,NSPEC
       write(*,'(a,i3,a)')'______> Atoms in species ',i,' ('//trim(Species(i))//') <______'
       do k=1,NspN(i)
          nat=nat+1
          write(*,'(i5,5x,3(f10.5,1x),x,3a)') nat,(TI(j,nat),j=1,3),(relax_flag(j,nat),j=1,3)
       end do
    end do
!
!............. read the k-points in fract coordinates (if available)
!
!___________ default option: DOS not available
!
    if(NKPTS.eq.1) then
       VKPT=0.0 ; WTKPT=1.0
    else
!
!___________ check if k-points are available (if NKPTS > 1)
!
       rewind(14)
       call find_string('BandPoints',10,line,14,.false.,iErr)
       if(iErr.eq.1) go to 31
       s=0.0
       do nkp=1,NKPTS
          read(14,*,err=31,end=31) (VKPT(i,nkp),i=1,3),WTKPT(nkp)  
          s=s+WTKPT(nkp)
       end do
       WTKPT=WTKPT/s
       read(14,'(a)',end=31,err=31) line      
       if(index(line,'endblock').eq.0) go to 31
       yeskp=.true.
    end if
    write(*,*) 'O.K.! This '//trim(seedfdf)//' is fine! Let us proceed!' 

    close (14)
!
!.......... showing k points
!
    do nkp = 1 , NKPTS
       write (*,'(a,3(1x,f10.5),5x,f10.5)') 'K-vectors (fract): ', &
                         VKPT(1,nkp),VKPT(2,nkp),VKPT(3,nkp), WTKPT(nkp)  
    end do
    Err=.false.
    return
     
!......... errors
29  write(*,*)'WARNING! Bad or absent SIESTA '//trim(seedfdf)//' file!'
    go to 31
30  write(*,*)'WARNING! k-points NOT available - no DOS possible!'
31  close (14)
    return
end subroutine read_siesta_geom

subroutine get_relax_flags_SIESTA(un)
!...........................................................................
! read geometry constraints from SIESTA fdf input file;
! the input file at unit = un  should be opened  already.
!...........................................................................
use param
use atoms
implicit none
character Line*200,seed0*40,cha(3)
integer LinEnd(100),LinPos(100),un,iErr,NUMLIN,m,n1,n2,i,i1,m1
real*8 x(3)

   call find_string('%block GeometryConstraints',26,Line,un,.true.,iErr)
   if(iErr.eq.1) return
   write(*,*)'FOUND block Geometry Constraints ...'
15 call CutStr(Line,NumLin,LinPos,LinEnd,un,0,iErr)
   write(*,'(a)') Line
   if(index(Line(1:9),'%endblock').ne.0) go to 17
   if(index(Line,'position').eq.0) go to 18
   if(index(Line(LinPos(2):LinEnd(2)),'from').ne.0) then
      cha(1:3)='F'
      if(Numlin.eq.8) then
         read(Line(LinPos(6):LinEnd(8)),*,err=18) x
         if(x(1).eq.0.0) cha(1)='T'
         if(x(2).eq.0.0) cha(2)='T'
         if(x(3).eq.0.0) cha(3)='T'
      else if(Numlin.ne.5) then
         go to 18
      end if
      read(Line(LinPos(3):LinEnd(3)),*,err=18) n1
      read(Line(LinPos(5):LinEnd(5)),*,err=18) n2
      if(n1.lt.0 .and. n2.lt.0) then
         m=NIONS+n1+1 ; n1=NIONS+n2+1 ; n2=m
      end if
      if(n1.gt.n2 .or. n1.lt.1.or.n2.gt.NIONS) go to 18
      do i=n1,n2
         i1=old_to_new(i)
         relax_flag(1:3,i1)=cha(1:3)
      end do
   else
      if(Numlin.lt.2) go to 18
      read(Line(LinPos(2):LinEnd(2)),*,err=18) m
      if(m.lt.1 .or. m.gt.NIONS) go to 18
      m1=old_to_new(m)
      relax_flag(1:3,m1)='F'
      if(Numlin.eq.5) then       
         read(Line(LinPos(3):LinEnd(5)),*,err=18) x
         if(x(1).eq.0.0) relax_flag(1,m1)='T'
         if(x(2).eq.0.0) relax_flag(2,m1)='T'
         if(x(3).eq.0.0) relax_flag(3,m1)='T'
      end if
   end if
   go to 15
17 write(*,*)'FINE: block Geometry Constraints read in successfully!'
   return
18 write(*,*)'ERROR: constraints not recognised! Ignored!'
   do i=1,NIONS
      relax_flag(1:3,i)='T' 
   end do
 end subroutine get_relax_flags_SIESTA
