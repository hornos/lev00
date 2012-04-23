subroutine read_vasp_geom(Err)
use param
use atoms
use kpoints
use code
implicit none
character Line*200
integer LinEnd(100),LinPos(100)
!..................................................................
!  AJ  - direct lattice vectors;
!  TI  - sublattices: they are stored in one array for all species,
!        one after another;
!..................................................................
character               :: how_geometry*5
logical good,Err,sel_dyn
integer i0,i1,i2,j0,j1,j2,i,nat,n,k,NumLin,iErr,j,NKP1,nkp,ITI
real*8 vscale,r(3),s
!
!.................. read VASP geometry
!
    Err=.false.
    write(*,*)'======================================================'
    write(*,*)'===== VASP geometry files found (my guess) ===========' 
    write(*,*)'======================================================'
    call system('ls *POSCAR* *CONTCAR*')
    write(*,*)'======================================================'
    write(*,*)' Give the name of the VASP geometry file:'
    read(*,'(a)') nameg
!
!__________________________ reading geometry
!
    good=.true.
    write(*,*)'Reading geometry from '//trim(nameg)//' ...'
    open (15,file=trim(nameg),status='old',FORM='FORMATTED',err=200)
    read (15,*)
    read (15,*) vscale
    read (15,*,err=200,end=200) DIRC(1,1), DIRC(1,2), DIRC(1,3)
    read (15,*,err=200,end=200) DIRC(2,1), DIRC(2,2), DIRC(2,3)
    read (15,*,err=200,end=200) DIRC(3,1), DIRC(3,2), DIRC(3,3)
    DIRC=DIRC*vscale
!
!_________ reciprocal vectors BCELL without 2*pi
    call bastr(DIRC,BCELL,VOLC,0)
!_________ reciprocal vectors RECC with 2*pi
    call bastr(DIRC,RECC,VOLC,1)

!_________ numbers of atoms NspN in every species; in v. 5.2 atomic species
!          symbols may go before the NspN numbers; we shuld identify and skip them
!
    read(15,*,err=21) (NspN(i),i=1,NSPEC)
    go to 25
21  read(15,*,err=200) (NspN(i),i=1,NSPEC)

!_________ checking against NIONS
!
25  n=0
    do i=1,NSPEC
       n=n+NspN(i)
    end do
    if(n.ne.NIONS) then
       Err=.true.
       write(*,*) 'Wrong number of atoms in species!'
       return
    else
       write(*,*)'O.K.! We have got NspN to be:'
       write(*,'(15(i3,1x))') (NspN(i),i=1,NSPEC)
       write(*,*)'Proceeding ...'
    end if

!__________ setting up atomic species numbers
!
    n=0
    do i=1,NSPEC
       do k=1,NspN(i)
          n=n+1 ; atom_species(n)=i
       end do
    end do
    
!_________ read actual coordinates
!
    nat=0
10  read(15,'(a)',err=200,end=200) line 
    i0=index(line,'Cart')
    i1=index(line,'cart')
    i2=index(line,'CART')
    j0=index(line,'direct')
    j1=index(line,'Direct')
    j2=index(line,'DIRECT')
    if(i0.ne.0 .or. i1.ne.0 .or. i2.ne.0 ) then
       how_geometry='carte'
       write(*,'(a)') '... implying Cartesian format, reading in ...'
       do i=1,NSPEC
          do k=1,NspN(i)
             nat=nat+1
             call CutStr(Line,NumLin,LinPos,LinEnd,15,0,iErr)
             if(iErr.eq.1.or.NumLin.lt.3) go to 200 
             read(Line(LinPos(1):LinEnd(3)),*,err=200) (TI(j,nat),j=1,3)
             TI(1:3,nat)=vscale*TI(1:3,nat)
             call check_atoms(nat,good)
             if(NumLin.eq.6) then
                sel_dyn=.true.
                read(Line(LinPos(4):LinEnd(4)),'(a)',err=200) relax_flag(1,nat)
                read(Line(LinPos(5):LinEnd(5)),'(a)',err=200) relax_flag(2,nat)
                read(Line(LinPos(6):LinEnd(6)),'(a)',err=200) relax_flag(3,nat)
             else
                sel_dyn=.false. ; relax_flag(1:3,nat)='F'
             end if
          end do
       end do
    else if(j0.ne.0 .or. j1.ne.0 .or. j2.ne.0 ) then
       how_geometry='fract'
       write(*,'(a)') '... implying fractional format, reading in ...'
       do i=1,NSPEC
          do k=1,NspN(i)
             nat=nat+1
             call CutStr(Line,NumLin,LinPos,LinEnd,15,0,iErr)
             if(iErr.eq.1.or.NumLin.lt.3) go to 200 
             read(Line(LinPos(1):LinEnd(3)),*,err=200) (r(j),j=1,3)
             TI(1,nat)=r(1)*DIRC(1,1)+r(2)*DIRC(2,1)+r(3)*DIRC(3,1)
             TI(2,nat)=r(1)*DIRC(1,2)+r(2)*DIRC(2,2)+r(3)*DIRC(3,2)
             TI(3,nat)=r(1)*DIRC(1,3)+r(2)*DIRC(2,3)+r(3)*DIRC(3,3)
             call check_atoms(nat,good)
             if(NumLin.eq.6) then
                sel_dyn=.true.
                read(Line(LinPos(4):LinEnd(4)),'(a)',err=200) relax_flag(1,nat)
                read(Line(LinPos(5):LinEnd(5)),'(a)',err=200) relax_flag(2,nat)
                read(Line(LinPos(6):LinEnd(6)),'(a)',err=200) relax_flag(3,nat)
             else
                sel_dyn=.false. ; relax_flag(1:3,nat)='F'
             end if
          end do
       end do
    else
       go to 10
    end if
!
    write(*,*)'Done geometry!'
    close (15)
    ITI=nat
    if(good) then
       write(*,*)'O.K.! This '//trim(nameg)//' is fine! Let us proceed!' 
       write(*,*)'........> Total number of atoms in the cell is ',ITI
    else
       write(*,*)'FATAL! Equivalent atoms in geometry file!' 
       Err=.true.
       return
    end if
!
!___________ show atoms 
!
    nat=0 ; yesrlx=.true.
    do i=1,NSPEC
       write(*,'(a25,i3,a)')'______> Atoms in species ',i,' ('//&
                              trim(Species(i))//') <______'
       do k=1,NspN(i)
          nat=nat+1
          write(*,'(i5,5x,3(f10.5,1x),3x,3(a,x))')  &
               nat,(TI(j,nat),j=1,3),(relax_flag(j,nat),j=1,3)
       end do
    end do
    do i=1,3
       write(*,'(a,i1,a,3(f10.5,x))') '[',i,']  ',(DIRC(i,j),j=1,3)
    end do
    close (15)
    write(*,*)'Done geometry!'
!
!................. read k-points from KPOINTS file
!
    write(*,*)'Trying to read k-points from KPOINTS ...'
    open (15,file='KPOINTS',status='old',FORM='FORMATTED',err=20)
    read(15,'(a)',err=20,end=20) line
    read(15,*,err=20,end=20) NKP1
    if(NKP1.ne.NKPTS) go to 20
    read(15,'(a)',err=20,end=20) line
    i0=index(line,'recip') ; i1=index(line,'Recip') ; i2=index(line,'RECIP')
    if(i0.ne.0 .or. i1.ne.0 .or. i2.ne.0 ) then
       do nkp=1,NKPTS
          read(15,*,err=20,end=20) (VKPT(i,nkp),i=1,3),WTKPT(nkp)  
       end do
       write(*,*) '... K-points have been read in successfully from KPOINTS!'
    else
       write(*,*)'Unknown format in KPOINTS.'
       go to 20
    end if
    close (15)
    go to 50  
!
!................. read k-points from OUTCAR file
!
20  write(*,*)'Trying to read k-points from '//trim(nameout)//' ...'
    s=0.0
    open (15,file=trim(nameout),status='old',FORM='FORMATTED',err=300)
11  read(15,'(a)',err=300,end=300) line
    i0=index(Line,'Subroutine IBZKPT') 
    i1=index(Line,'k-points in reciprocal lattice and weights')
    if(i0.ne.0) then
       read(15,'(//a)',err=300,end=300) line
       read(line,'(6x,i7)',err=300) NKP1
       if(NKP1.ne.NKPTS) go to 300
       read(15,'(//)',err=300,end=300)
    else if(i1.ne.0) then
       go to 12
    else
       go to 11
    end if
12  do nkp=1,NKPTS
       read(15,'(a)',err=300,end=300) line
       read(line,*,err=300,end=300) (VKPT(i,nkp),i=1,3),WTKPT(nkp)  
       s=s+WTKPT(nkp)
    end do
    do nkp=1,NKPTS
       WTKPT(nkp)=WTKPT(nkp)/s
       write(*,'(3(f10.5,x),5x,f10.5)') (VKPT(i,nkp),i=1,3),WTKPT(nkp)  
    end do
    close(15)
    write(*,*) '... K-points have been read in successfully from '//trim(nameout)
!
50  yeskp=.true.
    return
!
!......... errors
200 write(*,*)'FATAL! File '//trim(nameg)//' is bad or absent!'
    Err=.true.
    return
300 write(*,*)'FATAL! File '//trim(nameout)//' is bad or absent!'
    yeskp=.false.
    return
end subroutine read_vasp_geom

subroutine get_rwigs(RWIGS)
!.................................................................
!  VASP: read Wigner-Seits atomic radii from OUTCAR 
!.................................................................
use param
use code
implicit none
character Line*200
integer LinEnd(100),LinPos(100)
integer iErr,isp1,ii,NumLin,isp
real*8 RWIGS(NSPEC)
!
!........ default values for RWIGS from pseudopotentials reading from OUTCAR
!   RWIGS=1.0
!
   write(*,'(a)')'Trying to read RWIGS from '//trim(nameout)//' ...'
   open (15,file=trim(nameout),status='old',FORM='FORMATTED',err=30)
   do isp=1,NSPEC
1     read(15,'(a)') line
      if(index(line,'RWIGS').eq.0) go to 1
      call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
      if(NumLin.ge.6) then
         read(line(LinPos(6):LinEnd(6)),*,err=2,end=2) RWIGS(isp)
      else
         go to 2
      end if
   end do
   write(*,'(a/(10(f10.5,x)))') '... default values for RWIGS: ',(RWIGS(isp),isp=1,NSPEC)
   go to 3
2  write(*,'(a/(10(f10.5,x)))') 'WARNING: no default values for RWIGS found'
!
!........ reading entered values (if any)
!
   rewind(15)
3  call find_3strings('Atomic',6,'Wigner-Seitz',12,'radii',5,line,15,.false.,iErr)
   isp1=1
65 call CutStr(line,NumLin,LinPos,LinEnd,15,0,iErr)
   if(iErr.eq.1.or.NumLin.le.2) go to 30
   ii=NumLin-2
   read(line(LinPos(3):),*) (RWIGS(isp),isp=isp1,isp1+ii-1)
   if(isp1+ii-1.lt.NSPEC) then
      isp1=isp1+ii
      go to 65
   end if
   go to 31
30 write(*,*)'WARNING! The file OUTCAR is bad or absent!'
   write(*,*)'         Using defaults for RWIGS' 
31 RWIGS=abs(RWIGS)
   write(*,'(a,6(x,f10.5)/(18x,6(x,f10.5)))')  &
                       '... using RWIGS = ',(RWIGS(isp),isp=1,NSPEC)
   close (15)
   return
end subroutine get_rwigs

subroutine check_Mendel(Spec,k,Err,verbose)
  use mendeleev
  implicit none
  logical Err,verbose
  character Spec*2
  integer k
  Err=.false.

!___________ checking from NAZV
  do k=1,112
     if(NAZV(k).eq.Spec) then
        if(verbose) write(*,'(a,i2,a)')'Species '//Spec// &  
         ' => recognised as #',k,' in Mendeleev''s Table'
        return
     end if
  end do
  write(*,*)'ERROR! Species '//Spec//' has not been recognised!'
  Err=.true.
end subroutine check_Mendel
