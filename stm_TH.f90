subroutine stm_TH(grid)
!....................................................................
!   A box is chosen inside the cell and the density in the box is
! written into a file in the cube gOpenMol format.
!....................................................................
use param
use atoms
use menu
implicit none
real*8,parameter :: tiny=0.0001
real*8 GRID(NGX,NGY,NGZ)
real*8,dimension(3) :: Center=(/0.0,0.0,0.0/),Face=(/1.0,1.0,1.0/)
integer,dimension(3) :: ngrid=(/20,20,20/)
real*8 Direct(3,3),da(3),vect(3,3),Bot,Top,corner(3),LDOS,Zmin,Zmax
integer nChrg,j,i,ix,iy,iz,i1,i2,i3,nat,j1,Nltrp
logical Yes_Do,Yes_Z,Yes_Curr,Yes_XY,Yes_Spctr,Yes_Empty_lines
character cha2*2,Title*50,title_pl(9)*7,cha1
data Title/'                                                  '/
real*8,allocatable,dimension(:,:) :: Jmin,Jmax
real*8,allocatable,dimension(:,:) :: LatPos
real*8,allocatable,dimension(:)   :: Zldos
real*8 Imin,Imax,z,x(3),denval,z1,z2,zz,R(3),R1(3),R2(3),spctr(9)
!
!................. setting for the box orientation (fixed here)
     Direct=0.0 ; Direct(1,1)=1.0 ; Direct(2,2)=1.0 ; Direct(3,3)=1.0
!................. initial setting for the Z-window for atoms to be shown in [Sh]
     Zmax=0.0 ; Zmin=-2.0 
!................. initialisation for spectrosopy
     Nltrp=1 ; allocate(LatPos(3,Nltrp))

!................. start the main menu
!
      Yes_Do=.false. ; Yes_Z=.false. ; Yes_Curr=.false. ; Yes_XY=.false. 
      Yes_Spctr=.false. ; Yes_Empty_lines=.true.
1     write(*,*) 
      write(*,*)'............MENU for STM/TH ...........................'
      write(*,*)'***[ Surface MUST be perpendicular to the Z axis ! ]***'
      write(*,*)'......... Change these parameters if necessary:........'
      write(*,*)

      write(*,'(a)') '---------- T H - S T M   i m a g e  -------------'
      write(*,'(a,2(f10.5,a))') '   1. Surface area center is at: (', &
                Center(1),',',Center(2),')'

      write(*,'(a)')'      Directions of the surface area are fixed along:'
      write(*,'(10x,a,3(f10.5,a))') &
          '1   <',Direct(1,1),',',Direct(1,2),',',Direct(1,3),'  >'
      write(*,'(10x,a,3(f10.5,a))') &
          '2   <',Direct(2,1),',',Direct(2,2),',',Direct(2,3),'  >'
      write(*,'(a,2(1x,f10.5))') &
          '   2. Lengths of the surface area (in A) are: ',Face(1),Face(2)
      corner(1:2)=Center(1:2)-0.5*Face(1:2)

      if(Yes_Z) then
         Face(3)=Top-Bot ; corner(3)=Bot 
         write(*,'(5x,a,3(1x,f10.5))') &
              '>>>> the starting point (corner) of the box is at ', &
              (corner(i),i=1,3)
         if(Yes_Curr) then
            write(*,'(a,f10.5)') '   I. Value of LDOS within the window = ',LDOS
         else
            write(*,'(a)') '   I. Value of LDOS within the window = UNDERFINED'
         end if
      end if

      write(*,'(a,3(1x,i3))') &
      '   4. The grid  in the box in 3 directions: ',(ngrid(i),i=1,3)
      nChrg=ngrid(1)*ngrid(2)*ngrid(3)
      write(*,'(5x,a,i10)')'>>>> total number of grid points = ',nChrg
      if(Yes_Z) then
         do i=1,3
            da(i)=Face(i)/ngrid(i) ; vect(i,1:3)=Direct(i,1:3)*da(i)
         end do
         write(*,'(5x,a,3(x,f10.5))')'>>>> steps along 3 directions are ', &
              (da(i),i=1,3)
         if(Yes_Curr) then
            if(.not.Yes_Do) then
               write(*,'(a)') '  Ci. Calculate IMAGE of constant LDOS topography'
            else
               write(*,'(a)') '  Ci. Calculate constant LDOS topography - DONE!'
               write(*,'(a)') '  Vi. PREVIEW the LDOS' 
               write(*,'(a)') '  Ps. Create the Postscript of the LDOS' 
            end if
         end if
      end if
      write(*,'(a)') '  Sh. SHOW (X,Y)-coordinates of surface atoms within the PLOT'
      write(*,'(a)') '------------ S p e c t r o s c o p y ------------'
      write(*,'(a,i3)') '  Np. Number of points along Z axis = ',ngrid(3)
      write(*,'(a,i3)') '  Nl. Number of lateral points to consider = ',Nltrp
      if(Yes_XY) then
         write(*,'(a)') '  XY. Lateral points - KNOWN'
         write(*,'(a)') '  Sl. Show Lateral points'
         if(Yes_Z) then
            if(Yes_Spctr) then
               write(*,'(a)') '  Cs. Calculate spectroscopy curves- DONE!'
               write(*,'(a)') '  Vs. PREVIEW the spectroscopy' 
            else
               write(*,'(a)') '  Cs. Calculate spectroscopy curves'
            end if
         end if
      else
         write(*,'(a,i3)') '  XY. Specify lateral points - yet UNDEFINED'
      end if

      write(*,'(a)') '-------  G e n e r a l  s e t t i n g s ---------'
      if(Yes_Z) then
         write(*,'(2(a,f10.5))') &
              '  BT. Bottom and top of the box (along Z, in A) are: ',Bot,',',Top
         write(*,'(5x,2(a,f10.5),a)') '>>>> LDOS window found = [ ',Imin,',',Imax,' ]'
      else
         write(*,'(a)') &
              '  BT. Bottom and top of the box (along Z, in A) are UNDERFINED'
      end if
      
      write(*,'(a)') '  An. Coordinates are specified in: '//angstr
      write(*,'(a)') '  Co. Show current atomic positions in fractional/Cartesian'
      write(*,'(2(a,f10.5),a)') &
           '  zW. Z-window to SHOW (X,Y)-positions of atoms: [ ',Zmin,',',Zmax,' ]'
      if(Yes_Empty_lines) then
         write(*,*)'  El. Add empty lines in the constant LDOS file <--- YES'
      else
         write(*,*)'  El. Add empty lines in the constant LDOS file <--- NO'
      end if
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

      write(*,'(a)') '------------ C a l c u l a t o r ----------------'
      write(*,'(a)') '  cT. Centre of a triangle of points/atoms'
      write(*,'(a)') '  mD. Middle distance between two points/atoms'
      write(*,'(a)') '------- L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![1]__________ give the central point
!
      IF(trim(cha2).eq.'1') THEN
         write(*,*)'Give the center of your surface area :'
         write(*,*)'*** Z coordinate will be ignored ***'
         call givepoint(Center(1),Center(2),Center(3),angstr)
         Yes_Do=.false.
!
![2]__________ choose lengths of the area sides
!
      ELSE IF(trim(cha2).eq.'2') THEN
 12      WRITE(*,'(a)')'Give lengths of the surface area (in A)'
         read(*,*,err=12) (Face(j),j=1,2)
         Face(1)=abs(Face(1)) ; if(Face(1).lt.tiny) go to 12
         Face(2)=abs(Face(2)) ; if(Face(2).lt.tiny) go to 12
         Yes_Do=.false.
!
![BT]__________ choose bottom and top of the box (along Z)
!
      ELSE IF(trim(cha2).eq.'BT') THEN
 14      WRITE(*,'(a)')'Give BOTTOM and TOP of the surface area (in A)'
         read(*,*,err=14) Bot,Top
         if(Top-Bot .lt. tiny) go to 14
!
!___________________ analyse the LDOS: find the window of currents for the
!                    heights between Bot and Top
         write(*,*)'... Analysing the LDOS ...'
         allocate (Jmin(NGX,NGY)) ; allocate (Jmax(NGX,NGY))
         do ix=1,NGX
            do iy=1,NGY
               Imin=1.e10 ; Imax=-1.e10
               do iz=1,NGZ
                  z = (ix-1)*DIRC(1,3)/NGX + (iy-1)*DIRC(2,3)/NGY + (iz-1)*DIRC(3,3)/NGZ
                  if(z.ge.Bot .and. z.le.Top) then
                     if(grid(ix,iy,iz) .gt. Imax) Imax=grid(ix,iy,iz)
                     if(grid(ix,iy,iz) .lt. Imin) Imin=grid(ix,iy,iz)
                  endif
               end do
               Jmin(ix,iy)=Imin
               Jmax(ix,iy)=Imax
            end do
         end do
!_______________________ find the lagest minimum value Imin and the smallest maximum
!                        value Imax across all lateral grid points
         Imax= 1.e10 ; Imin=-1.e10
         do ix=1,NGX
            do iy=1,NGY
               if(Jmin(ix,iy) .ge. Imin) Imin=Jmin(ix,iy)
               if(Jmax(ix,iy) .lt. Imax) Imax=Jmax(ix,iy)
            end do
         end do
         deallocate (Jmin) ; deallocate (Jmax)
         write(*,'(2(a,f10.5),a)') '... LDOS window found = [ ',Imin,',',Imax,' ]'
         Yes_Z=.true. ; Yes_Do=.false. ; Yes_Spctr=.false. ; Yes_Curr=.false.
!
![4]__________ choose the grid in the box
!
      ELSE IF(trim(cha2).eq.'4') THEN
 13      write(*,'(a)') 'Specify the number of grid points in all 3 directions:'
         read(*,*,err=13) (ngrid(i),i=1,3)
         do i=1,3
           if(ngrid(i).lt.1) go to 13
         end do
         Yes_Do=.false.
!
![I]__________ choose LDOS value within the window
!
      ELSE IF(trim(cha2).eq.'I' .and. Yes_Z) THEN
21       write(*,'(2(a,f10.5),a)')'Specify the LDOS value within [ ',Imin,',',Imax,' ]'
         read(*,*,err=21) LDOS
         if(LDOS.lt.Imin .or. LDOS.gt.Imax) go to 21
         Yes_Curr=.true. ; Yes_Do=.false.
!
![Ci]__________ calculate the current topography
!
      ELSE IF(trim(cha2).eq.'Ci' .and. Yes_Z .and. Yes_Curr) THEN
         allocate(Zldos(ngrid(3)))
         write(*,'(a)') 'Opening the file for the constant LDOS ....'
         open(31,file='ldos.dat',form='formatted',status='unknown')
         do i1=0,ngrid(1)-1
            do i2=0,ngrid(2)-1

               do i3=0,ngrid(3)-1
!_________________________ Cartesian position of the grid point in the box
                  x(1)= i1*vect(1,1)+i2*vect(2,1)+i3*vect(3,1)+corner(1)
                  x(2)= i1*vect(1,2)+i2*vect(2,2)+i3*vect(3,2)+corner(2)
                  x(3)= i1*vect(1,3)+i2*vect(2,3)+i3*vect(3,3)+corner(3)
!_________________________ interpolate the density at x()
                  call reducn(x,DIRC,BCELL)
                  call interpolate(x,BCELL,denval,grid)
                  Zldos(i3+1)=denval
               end do
!_________________________ determine the height corresponding to the given LDOS
               do i3=1,ngrid(3)-1
                  if( (LDOS.ge.Zldos(i3) .and. LDOS.lt.Zldos(i3+1)) .or. &
                          (LDOS.ge.Zldos(i3+1) .and. LDOS.lt.Zldos(i3)) ) then
                     x(1)= i1*vect(1,1)+i2*vect(2,1)+i3*vect(3,1)+corner(1)-Center(1)
                     x(2)= i1*vect(1,2)+i2*vect(2,2)+i3*vect(3,2)+corner(2)-Center(2)
                     z1=i1*vect(1,3)+i2*vect(2,3)+i3*vect(3,3)+corner(3)
                     z2=i1*vect(1,3)+i2*vect(2,3)+(i3+1)*vect(3,3)+corner(3)
                     zz=z1+(LDOS-Zldos(i3))*(z2-z1)/(Zldos(i3+1)-Zldos(i3))
                     write(31,'(3(e13.5,x))') x(1),x(2),zz 
                     exit
                  end if
               end do

            end do
            if(Yes_Empty_lines) write(31,*)
         end do
         close (31)
         write(*,*)'Done! The file <ldos.dat> created!'
         deallocate(Zldos)
         Yes_Do=.true.
!
![Vi]__________ preview LDOS image
!
      ELSE IF(trim(cha2).eq.'Vi' .and. Yes_Do) THEN
         call Plot3d('ldos.dat',8,Title, &
              'X-coordinate (A)    ','Y-coordinate (A)    ', &
              'Z-coordinate        ','Screen', 33,nclasses,type_prv)
!
![Ps]__________ create a Postscript file for the LDOS
!
      ELSE IF(trim(cha2).eq.'Ps' .and. Yes_Do) THEN
         call Plot3d('ldos.dat',8,Title, &
              'X-coordinate (A)    ','Y-coordinate (A)    ', &
              'Z-coordinate        ','Postsc', 33,nclasses,type_prv)
         write(*,*)'... Postscript file <ldos.dat.ps> has been created' 
!
![Sh]__________ show atomic positions via their X,Y coordinates on the PLOT [V];
!               only atoms within the specified Z-window [Zmin,Zmax] are given 
!
      ELSE IF(trim(cha2).eq.'Sh') THEN
         write(*,'(a)')'*****************************************************'
         write(*,'(a)')'*** Only atoms within the Z-window are shown.     ***'
         write(*,'(a)')'*** Only atoms within the surface area are shown. ***'
         write(*,'(2(a,f10.5),a)') &
                       '*** Z-window = [',Zmin,',',Zmax,']            ***'
         write(*,*)    '*****************************************************'
         nat=0
         do i=1,NSPEC
            write(*,'(3a,i5)')'.... Species <',Species(i),'> ......> ',i
            write(*,'(a)')' Tot   Sp#      |-X-in-PLOT-| |-Y-in-PLOT-| [--actual-Z-]'
            do j=1,NspN(i)
               nat=nat+1
               if(TI(3,nat).ge.Zmin .and. TI(3,nat).le.Zmax) then
                  x(1:2)=TI(1:2,nat)-Center(1:2)
                  if(abs(x(1)).le.0.5*Face(1) .and. abs(x(2)).le.0.5*Face(2)) then 
                     write(*,'(2(x,i5),2(4x,f10.5),2x,a,f10.5,a)') &
                          nat,j,(x(j1),j1=1,2),'  [',TI(3,nat),' ]'
                  end if
               end if
            end do 
         end do
         write(*,*)'Hit ENTER when done ...' ; read(*,*)
!
![Np]__________ number of points along Z for spectroscopy
!
      ELSE IF(trim(cha2).eq.'Np') THEN
68       write(*,*)'Specify number of points for 1D plots:'
         read(*,*,err=68) ngrid(3)
         if(ngrid(3).lt.2) go to 68
         Yes_Spctr=.false.
!
![Nl]__________ number of lateral points for spectroscopy
!
      ELSE IF(trim(cha2).eq.'Nl') THEN
         deallocate(LatPos)
69       write(*,*)'Specify number of lateral points for 1D plots (<10):'
         read(*,*,err=69) Nltrp
         if(Nltrp.lt.1 .or. Nltrp.gt.9) go to 69
         allocate(LatPos(3,Nltrp)) 
         Yes_Spctr=.false.
!
![XY]__________ specify lateral points for spectroscopy
!
      ELSE IF(trim(cha2).eq.'XY') THEN
         do i=1,Nltrp
            write(*,'(a,i2,a)') &
                 'Specify X,Y,Z of the lateral point [',i,'] (only X,Y are kept)'
            call givepoint(LatPos(1,i),LatPos(2,i),LatPos(3,i),angstr)
            write(cha1,'(i1)') i
            title_pl(i)='point '//cha1
         end do
         Yes_XY=.true. ;  Yes_Spctr=.false.

!
![Sl]__________ show the lateral points for spectroscopy
!
      ELSE IF(trim(cha2).eq.'Sl' .and. Yes_XY) THEN
         write(*,'(a)')' #  --actual-X- --actual-Y- |  X-in-2D-PLOT Y-in-2D-PLOT'
         do i=1,Nltrp
            write(*,'(i2,2(x,f10.5),3x,a,2(x,f10.5))') &
     i,LatPos(1,i),LatPos(2,i),' | ',LatPos(1,i)-Center(1),LatPos(2,i)-Center(2)
         end do
         write(*,*)'Hit ENTER when done ...' ; read(*,*)
!
![Cs]__________ calculate spectroscopy curves
!
      ELSE IF(trim(cha2).eq.'Cs' .and. Yes_XY .and. Yes_Z) THEN
         open(31,file='spectr.dat',form='formatted',status='unknown')
         do i3=0,ngrid(3)-1
            do i=1,Nltrp
               x(1)=LatPos(1,i) ; x(2)=LatPos(2,i)  ; x(3)=corner(3)+(i3-1)*da(3)
!_________________________ interpolate the density at R() (will be changed!)
               R=x
               call reducn(R,DIRC,BCELL)
               call interpolate(R,BCELL,spctr(i),grid)
            end do
            write(31,'(f10.5,9(x,f10.5))') x(3),(spctr(i),i=1,Nltrp)
         end do
         close(31)
         Yes_Spctr=.true.
!
![Vs]__________ preview spectroscopy curves
!
      ELSE IF(trim(cha2).eq.'Vs' .and. Yes_Spctr) THEN
         call Plot_bunch(Nltrp,'spectr.dat',10,Title,title_pl, &
                'Z (A)               ',&
                'LDOS (arb. units)   ', 'Screen', 33)
!
![An]__________ choose the way how the coordinates are given
!
      ELSE IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
![cT]__________ centre of a triangle
!
      ELSE IF(trim(cha2).eq.'cT') THEN
         write(*,*)'Give the 1st point/atom:'
         call givepoint(R(1),R(2),R(3),angstr)
         write(*,*)'Give the 2nd point/atom:'
         call givepoint(R1(1),R1(2),R1(3),angstr)
         write(*,*)'Give the 3rd point/atom:'
         call givepoint(R2(1),R2(2),R2(3),angstr)
         x=(R+R1+R2)/3.0
         write(*,*)'The centre of the triangle is at:',(x(i),i=1,3)

!
![mD]__________ middle distance between two points
!
      ELSE IF(trim(cha2).eq.'mD') THEN
         write(*,*)'Give the 1st point/atom:'
         call givepoint(R(1),R(2),R(3),angstr)
         write(*,*)'Give the 2nd point/atom:'
         call givepoint(R1(1),R1(2),R1(3),angstr)
         x=(R+R1)/2.0
         write(*,*)'The middle distance between the points is at:', &
              (x(i),i=1,3)
!
![3D].... contour or 3-dimensional for preview
!
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
!
![cl].... number of contour levels for preview
!
      ELSE IF(trim(cha2).eq.'cl' .and. type_prv.eq.'2d-old-') THEN
75       write(*,*)'Specify the number of contour levels (2D only; must be > 2):'
         read(*,*,err=75) nclasses
         if(nclasses.lt.2) go to 75
!
![zW].... Z-window: z-values of atoms to be shown in option 'Sh'
!
      ELSE IF(trim(cha2).eq.'zW') THEN
76       write(*,*)'Specify the Z window of coordinates for atoms to be shown in [Sh].'
         write(*,*)'Give the low and high Z coordinates (A):'
         read(*,*,err=76) Zmin,Zmax
         if(Zmax-Zmin.lt.tiny) go to 76
!
![El].... whether add or not empty lines when saving the constant LDOS plot 
!
      ELSE IF(trim(cha2).eq.'El') THEN
         if(Yes_Empty_lines) then
            Yes_Empty_lines=.false.
         else
            Yes_Empty_lines=.true.
         end if
!
![Q]__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         deallocate(LatPos)
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine stm_TH
