subroutine line(grid,DIRC,BCELL,VOLC,filen,lenght)
!...................................................................
!  Line Calculation of the MAP.
!  31 - unit number for the file filen(1:lenght) with output data.
!  iQuit = 0 - all parameters are properly defined; can plot
!          1 - there are undefined parameters; cannot plot
!....................................................................
use param
use menu
implicit none 
real*8 GRID(NGX,NGY,NGZ),DIRC(3,3),BCELL(3,3),R(3),x(3),dzero
real*8 bCENTX,bCENTY,bCENTZ,denval
character filen*12,Title*50,title_pl*7,cha2*2
data Title/'                                                  '/
data title_pl/'       '/,dzero/0.0/
logical Yes_Do
real*8 fCENTX,fCENTZ,fCENTY,a,xcoord,absden,VOLC
integer iQuit,k2,lenght
!......................................................................
!....................... LINE MENU ....................................
!......................................................................
!_____ choose the vector along the line and normalize it;
!      give starting point; give length.
!......................................................................
      Yes_Do=.false.
1     iQuit=0
      write(*,*)'.............. LINE MENU .......................'
      write(*,*)'...... Change these parameters if necessary:....'
      write(*,*)
      write(*,'(a35,f10.5,2(a1,f10.5),a1)') &
      '   1. Starting point (Angstroms): (', &
                             aCENTX,',',aCENTY,',',aCENTZ,')'
      fCENTX=BCELL(1,1)*aCENTX+BCELL(1,2)*aCENTY+BCELL(1,3)*aCENTZ
      fCENTY=BCELL(2,1)*aCENTX+BCELL(2,2)*aCENTY+BCELL(2,3)*aCENTZ
      fCENTZ=BCELL(3,1)*aCENTX+BCELL(3,2)*aCENTY+BCELL(3,3)*aCENTZ
      write(*,'(a36,f10.5,2(a1,f10.5),a1)') &
      '      Starting point (fractional): (', &
                             fCENTX,',',fCENTY,',',fCENTZ,')'
      a=vers0x*vers0x + vers0y*vers0y + vers0z*vers0z
      if(a.eq.dzero) then
         iQuit=1
         write(*,'(a)') '   2. Vector along the line: ....... undefined .......'
      else
         write(*,'(a30,f10.5,2(a1,f10.5),a1)') &
      '   2. Vector along the line: (', &
                             vers0x,',',vers0y,',',vers0z,')'
      end if

      write(*,'(a44,f10.5)') &
        '   3. Lenght along the line (in Angstroms): ',width1

      write(*,'(a)')'   4. Parameters for the plotting'

      if(Yes_Do) then
        write(*,'(a)') &
         '   5. Perform calculation for the plotting: file '//filen &
                                             //' <= DONE!'
      else
        write(*,'(a)') '   5. Perform calculation for the plotting: file '//filen
      end if

      write(*,'(a)')'   6. Preview the density'
      write(*,'(a)')'   7. Create a postscript file '// &
                                trim(filen)//'.ps for the plot'

      write(*,'(a)') &
                    '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
      write(*,'(a)') &
          '  Co. Show current atomic positions in fractional/Cartesian'
      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![An]__________ choose a way how the coordinates are given
!
      IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![1]__________ give starting point for the line
!
      ELSE IF(trim(cha2).eq.'1') THEN
         if(angstr.eq.'<AtomNumber>') &
               write(*,*)'Specify the 1st atom to be started from.'
         call givepoint(aCENTX,aCENTY,aCENTZ,angstr)
         Yes_Do=.false.
!
![2]__________ give a vector along the line
!
      ELSE IF(trim(cha2).eq.'2') THEN
         if(angstr.eq.'<AtomNumber>') then
           write(*,*)'Specify the 2nd atom to be connected with.'
           call givepoint(bCENTX,bCENTY,bCENTZ,angstr)
           vers0x=bCENTX-aCENTX
           vers0y=bCENTY-aCENTY
           vers0z=bCENTZ-aCENTZ
           WIDTH1=sqrt( vers0x**2+vers0y**2+vers0z**2 )
         else
 7         write(*,*)'Give a vector (x,y,z) along your line:'
           read (*,*,err=7)  vers0x, vers0y, vers0z
         end if
         call normalize(vers0x,vers0y,vers0z)
         Yes_Do=.false.
!
![3]__________ give length along the line
!
      ELSE IF(trim(cha2).eq.'3') THEN
   10    WRITE(*,*) 'Enter length (in Angstroms):'
         READ(*,*,err=10) WIDTH1
         if(width1.lt.dzero) go to 10
         Yes_Do=.false.
!
![4]__________ give the resolution in either direction, chop values
!  and a multiplication factor for the density, etc.
!
      ELSE IF(trim(cha2).eq.'4') THEN
         call choose1()
         Yes_Do=.false.
!
![5]__________ perform calculation of the density along the line
!
      ELSE IF(trim(cha2).eq.'5') THEN
         if(iQuit.ne.0) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
         open(31,file=filen(:lenght),status='unknown',form='formatted')
         write(*,*)'The file '//filen//' has been opened for the PLOT.'
         write (*,*)'Writing to the file '//filen(1:lenght)//' ...'
         do k2=0,NRESOL
            R(1)=acentx+k2*vers0x*width1/nresol
            R(2)=acenty+k2*vers0y*width1/nresol
            R(3)=acentz+k2*vers0z*width1/nresol
            call reducn(R,DIRC,BCELL)
            call interpolate(R,BCELL,denval,grid)
            xcoord=k2*width1/nresol
            absden=(denval*multcon)/VOLC
            if(lochop.ne.hichop) then
               if(absden.gt.hichop) then
                  absden=hichop
               else if(absden.lt.lochop) then
                  absden=lochop
               end if
            end if
            write(31,*) xcoord,absden
         end do
         close (31)
         write(*,*)'.... File '//trim(filen)//' has been created! ....'
         Yes_Do=.true.
!
![6]__________ preview the file just created
!
      ELSE IF(trim(cha2).eq.'6') THEN
         if(Yes_Do) then
           call Plot1(filen,lenght,Title,title_pl, &
                'Coordinate (A)      ', &
                'Charge density      ',  'Screen', 33,0, &
                 'N',.false.,dzero,dzero)
         else
           write(*,*) &
             'IGNORED! You have to accomplish the item 5 first!'
         end if
!
![7]__________ create a PostScript file of the plot
!
      ELSE IF(trim(cha2).eq.'7') THEN
         if(Yes_Do) then
           write(*,*)'Give the title:'
           read(*,'(a)') Title
           call Plot1(filen,lenght,Title,title_pl, &
                'Coordinate (A)      ', &
                'Charge density      ', 'Postsc', 33,0, &
                 'N',.false.,dzero,dzero)
         else
           write(*,*) &
             'IGNORED! You have to accomplish the item 5 first!'
         end if
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
!__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine line

subroutine plane(grid,DIRC,BCELL,VOLC,filen,lenght)
!....................................................................
!  Plane Calculation of the MAP.
!  nfile - unit number for the file filen(1:lenght) with output data.
!  iQuit = 0 - all parameters are properly defined; can plot
!          1 - there are undefined parameters; cannot plot
!....................................................................
use param
use menu
implicit none
real*8 DIRC(3,3),BCELL(3,3),pA(2),pB(2),pC(2),x(3),GRID(NGX,NGY,NGZ),R(3)
real*8,parameter :: tiny=0.00001,dzero=0.0
character filen*12, Title*50,cha2*2
data Title/'                                                  '/
logical Yes_Do
integer iQuit,k3,k2,lenght3,lenght
real*8 a,fCENTX,fCENTY,fCENTZ,displ,xcoord,ycoord,VOLC,denval,absden
!......................................................................
!....................... PLANE MENU ...................................
!......................................................................
!_____ choose the vector along the plane normal and normalize it;
!      other two vectors lying in the plane are then generated
!      (rather arbitrarily though);
!      give center point of the plane; give widths.
!......................................................................
      Yes_Do=.false.
1     iQuit=0
      write(*,*)'..............MENU for PLANE .......................'
      write(*,*)'........ Change these parameters if necessary:......'
      write(*,*)
      a=vers1x*vers1x + vers1y*vers1y + vers1z*vers1z
      if(a.lt.tiny) then
         iQuit=1
         write(*,'(a)') &
      ' / 1. Normal vector to the plane: ....... undefined .......'
         write(*,'(a)') &
      ' |    X1 vector in the plane:     ....... undefined .......'
         write(*,'(a)') &
      ' |    Y1 vector in the plane:     ....... undefined .......'
      else
         iQuit=0
         write(*,'(a35,f10.5,2(a1,f10.5),a1)') &
      ' / 1. Normal vector to the plane: (', &
                             vers1x,',',vers1y,',',vers1z,')'
         write(*,'(a31,f10.5,2(a1,f10.5),a1)') &
      ' |    X1 vector in the plane: (', &
                             vers2x,',',vers2y,',',vers2z,')'
         write(*,'(a31,f10.5,2(a1,f10.5),a1)') &
      ' |    Y1 vector in the plane: (', &
                             vers3x,',',vers3y,',',vers3z,')'
      end if
      if(icase.eq.1) then
         write(*,'(a)') &
        ' \\ 2. The plane has been specified by 3 points: NO'
      else if(icase.eq.2) then
         write(*,'(a)') &
        ' \\ 2. The plane has been specified by 3 points: YES'
      end if
      if(icase.eq.2.and.icase1.eq.1) then
         write(*,'(a)') &
          '   3. Central point => the center of the triangle: NO'
      else if(icase.eq.2.and.icase1.eq.2) then
         write(*,'(a)') &
          '   3. Central point => the center of the triangle: YES'
         aCENTX=(Ra(1)+Rb(1)+Rc(1))/3.
         aCENTY=(Ra(2)+Rb(2)+Rc(2))/3.
         aCENTZ=(Ra(3)+Rb(3)+Rc(3))/3.
         central_p=.true.
      end if
      if(central_p) then
         if(icase.eq.2.and.icase1.eq.2) then
            write(*,'(a)') '      Central point on the plane: '
         else if(icase1.eq.1) then
            write(*,'(a)') '   4. Central point on the plane: '
         end if
         write(*,'(a29,f10.5,2(a1,f10.5),a1)') &
                            '           in Angstroms  => (', &
                             aCENTX,',',aCENTY,',',aCENTZ,')'
         fCENTX=BCELL(1,1)*aCENTX+BCELL(1,2)*aCENTY+BCELL(1,3)*aCENTZ
         fCENTY=BCELL(2,1)*aCENTX+BCELL(2,2)*aCENTY+BCELL(2,3)*aCENTZ
         fCENTZ=BCELL(3,1)*aCENTX+BCELL(3,2)*aCENTY+BCELL(3,3)*aCENTZ
         write(*,'(a29,f10.5,2(a1,f10.5),a1)') &
                              '           in fractional => (', &
                             fCENTX,',',fCENTY,',',fCENTZ,')'
         if(icase.eq.2) then
            write(*,'(11x,(a))') &
                 'The reference points A,B,C in (X1,Y1) are given as:'
            pA(1)=vers2x*(Ra(1)-aCENTX)+vers2y*(Ra(2)-aCENTY)+ &
                                                vers2z*(Ra(3)-aCENTZ)
            pB(1)=vers2x*(Rb(1)-aCENTX)+vers2y*(Rb(2)-aCENTY)+ &
                                                vers2z*(Rb(3)-aCENTZ)
            pC(1)=vers2x*(Rc(1)-aCENTX)+vers2y*(Rc(2)-aCENTY)+ &
                                                vers2z*(Rc(3)-aCENTZ)
            pA(2)=vers3x*(Ra(1)-aCENTX)+vers3y*(Ra(2)-aCENTY)+ &
                                                vers3z*(Ra(3)-aCENTZ)
            pB(2)=vers3x*(Rb(1)-aCENTX)+vers3y*(Rb(2)-aCENTY)+ &
                                                vers3z*(Rb(3)-aCENTZ)
            pC(2)=vers3x*(Rc(1)-aCENTX)+vers3y*(Rc(2)-aCENTY)+ &
                                                vers3z*(Rc(3)-aCENTZ)
            write(*,13) 'A = (',pA(1),',',pA(2),')'
            write(*,13) 'B = (',pB(1),',',pB(2),')'
            write(*,13) 'C = (',pC(1),',',pC(2),')'
 13         format(15x,a5,f10.5,a1,f10.5,a1)
         end if
         write(*,'(a)') &
               '  44. Move the plane parallel to the existing one '
      else
         write(*,'(a)') &
          '   4. Central point on the plane: ....... undefined .......'
      end if
      write(*,'(a39,f10.5)') '   5. Width along X1 axis (Angstroms): ',width1
      write(*,'(a39,f10.5)') '   6. Width along Y1 axis (Angstroms): ',width2
      write(*,'(a)')'   7. Parameters for the plotting'

      write(*,'(a)')'   8. Preview the density'
      if(Yes_Do) then
        write(*,'(a)') &
         '   9. Perform calculation for the plotting: file '//filen &
                                             //' <= DONE!'
      else
        write(*,'(a)') &
         '   9. Perform calculation for the plotting: file '//filen
      end if

      write(*,'(a)') &
                    '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
      write(*,'(a)') &
          '  Co. Show current atomic positions in fractional/Cartesian'
      if(central_p .and. a.gt.tiny) write(*,'(a)') &
          '  Ap. Show atomic positions with respect to the plane'
      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![An]__________ choose the way how the coordinates are given
!
      IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![1]__________ give a normal vector to the plane and generate two
!           others in the plane
!
      ELSE IF(trim(cha2).eq.'1') THEN
         icase=1
         icase1=1
         call vector3(DIRC)
         Yes_Do=.false.
!
![2]__________ specify the plane by 3 points
!
      ELSE IF(trim(cha2).eq.'2') THEN
         icase=2
         if(NIONS.le.2 .and. angstr.eq.'<AtomNumber>') then
            write(*,*)'WARNING! Not enough atoms for this option!'
            write(*,*)'Change to <Fractional> or <Angstroms> using An!'
            write(*,*)'Hit ENTER when ready ...'
            read(*,*) 
         else
            call vector3(DIRC)
         end if
         Yes_Do=.false.
!
![3]__________ give a method to choose the central point on the plane
!           in the case of 3 points (icase=2)
!
      ELSE IF(trim(cha2).eq.'3' .and. icase.eq.2) THEN
         if(icase1.eq.1) then
           icase1=2
           aCENTX=(Ra(1)+Rb(1)+Rc(1))/3.
           aCENTY=(Ra(2)+Rb(2)+Rc(2))/3.
           aCENTZ=(Ra(3)+Rb(3)+Rc(3))/3.
           central_p=.true.
           Yes_Do=.false.
         else if(icase1.eq.2) then
           icase1=1
         end if
!
![4]__________ give central point on the plane in a general way
!
      ELSE IF(trim(cha2).eq.'4' .and. icase1.eq.1) THEN
         if(iQuit.eq.1) then
            write(*,*)'ERROR! You must acomplish the item 1 first!'
         else
            call centralP(DIRC)
            central_p=.true.
            Yes_Do=.false.
         end if
!
!
![44]__________ give central point on the plane parallel to the existing one
!
      ELSE IF(trim(cha2).eq.'44'.and.central_p) THEN
         if(iQuit.eq.1) then
            write(*,*)'ERROR! You must acomplish the item 1 first!'
         else
 20         write(*,*)'Specify the displacement of the existing plane:'
            read(*,*,err=20) displ
            aCENTX=aCENTX+displ*vers1x
            aCENTY=aCENTY+displ*vers1y
            aCENTZ=aCENTZ+displ*vers1z
            icase1=1
            Yes_Do=.false.
         end if
!
![5,6]__________ give length along the X1,Y1 axes
!
      ELSE IF(trim(cha2).eq.'5') THEN
 10      WRITE(*,*) 'Enter length along X1 axis (in Angstroms):'
         READ(*,*,err=10) WIDTH1
         if(width1.lt.dzero) go to 10
         Yes_Do=.false.

      ELSE IF(trim(cha2).eq.'6') THEN
 11      WRITE(*,*) 'Enter length along Y1 axis (in Angstroms):'
         READ(*,*,err=11) WIDTH2
         if(width2.lt.dzero) go to 11
         Yes_Do=.false.
!
![7]__________ give the resolution in either direction, chop values
!  and a multiplication factor for the density, etc.
!
      ELSE IF(trim(cha2).eq.'7') THEN
         call choose3()
         Yes_Do=.false.
!
![8]__________ preview
!
      ELSE IF(trim(cha2).eq.'8') THEN
         if(iQuit.eq.0 .and. central_p .and. &
                          width1.gt.dzero .and. width2.gt.dzero) then
            open(32,file='test.dat',status='unknown',form='formatted')
            write(*,*)'The file test.dat has been opened to preview.'
            write(*,*)'Working on previewing ...'
            DO K3=-NRESOL_PRV/2,NRESOL_PRV/2
               DO K2=-NRESOL_PRV/2,NRESOL_PRV/2
                  R(1)=acentx+(k2*vers2x*width1 + &
                                   k3*vers3x*width2)/nresol_prv
                  R(2)=acenty+(k2*vers2y*width1 + &
                                   k3*vers3y*width2)/nresol_prv
                  R(3)=acentz+(k2*vers2z*width1 + &
                                   k3*vers3z*width2)/nresol_prv
                  call reducn(R,DIRC,BCELL)
                  call interpolate(R,BCELL,denval,grid)
                  xcoord=k2*width1/nresol_prv
                  ycoord=k3*width2/nresol_prv
                  absden=(denval*multcon)/VOLC
                  if(lochop.ne.hichop) then
                     if(absden.gt.hichop) then
                        absden=hichop
                     else if(absden.lt.lochop) then
                        absden=lochop
                     end if
                  end if
                  write(32,'(3(1x,e11.5))') xcoord,ycoord,absden
               END DO
               write(32,*)
            END DO
            close(32)
!______ plot the density: previewing
            lenght3=8
            call Plot3d('test.dat',lenght3,Title, &
                 'X-coordinate (A)    ','Y-coordinate (A)    ', &
                 'Charge density      ',  'Screen', 33, &
                 nclasses,nresol_prv,type_prv)
         else
            write(*,*)'ERROR! You still have undefined parameters!'
         end if
!
![9]__________ real calculation
!
      ELSE IF(trim(cha2).eq.'9') THEN
         if(iQuit.eq.0 .and. central_p .and. &
                             width1.gt.dzero .and. width2.gt.dzero) then
            open(31,file=filen(:lenght),status='unknown',form='formatted')
            write(*,*) &
                 'The file '//filen//' has been opened for the PLOT.'
            write(*,*)'Working on the real plot: writing to '// &
                 filen(1:lenght)//' ...'
            DO K2=-NRESOL3/2,NRESOL3/2
               DO K3=-NRESOL3/2,NRESOL3/2
                 R(1)=acentx+(k2*vers2x*width1+k3*vers3x*width2)/nresol3
                 R(2)=acenty+(k2*vers2y*width1+k3*vers3y*width2)/nresol3
                 R(3)=acentz+(k2*vers2z*width1+k3*vers3z*width2)/nresol3
                 call reducn(R,DIRC,BCELL)
                 call interpolate(R,BCELL,denval,grid)
                 xcoord=k2*width1/nresol3
                 ycoord=k3*width2/nresol3
                 absden=(denval*multcon)/VOLC
                 if(lochop.ne.hichop) then
                    if(absden.gt.hichop) then
                       absden=hichop
                    else if(absden.lt.lochop) then
                       absden=lochop
                    end if
                 end if
                 write(31,'(3(1x,e11.5))') xcoord,ycoord,absden
               END DO
            END DO
            close(31)
            write(*,*) &
               '.... File '//filen(1:lenght)//' has been created! ....'
            Yes_Do=.true.
         else
            write(*,*)'ERROR! You still have undefined parameters!'
         end if
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
![Ap].... display atomic positions with respect to the local coordinate system
!         fixed with the chosen plane
!
      ELSE IF(trim(cha2).eq.'Ap' .and. central_p .and. iQuit.eq.0) THEN
         call show_atoms_in_plane()
!
![Q]__________ quit option
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         return
      else
         write(*,*)'ERROR! Try again!'
      end if
      go to 1
end subroutine plane

subroutine charge_sph(grid,DIRC,BCELL,VOLC,totdens,filen,lenght)
!....................................................................
!  Charge inside a sphere.
!  nfile - unit number for the file filen(1:lenght) with output data.
!....................................................................
use param
use menu
implicit none
real*8 GRID(NGX,NGY,NGZ),DIRC(3,3),BCELL(3,3),R(3),x(3),factor,rCharge
real*8 Radius(Num_Sph),NumbEl_max,fCENTX,fCENTY,fCENTZ,drad,Rad_Fit
character iask,cha,filen*12,Title*50,title_pl*7,cha2*2
logical Yes_Obtain
data Title/'                                                  '/
data title_pl/'       '/
real*8,parameter :: dzero=0.0
logical Yes_Do
integer i,N_Sp,iQuit,ii,k1,k2,k3,ij,j1,j2,j3,iPnt,iX,iY,iZ,lenght,j
real*8 dV,dX,aR,denval,Rad,c1,c2,rad1,factor_way,totdens,VOLC,charg,rad2

      do i=1,Num_Sph
         Point(1,i)=dzero
         Point(2,i)=dzero
         Point(3,i)=dzero
      end do
      N_Sp=1
      Shrink(1)=1.0
!......................................................................
!_____ choose the starting point, the smallest and the largest radii,
!      the number of points in between and the grid inside the sphere.
! iQuit = 0 - not quit, proceed with plotting in the parent program;
!         1 - quit, do not proceed with plotting.
! method='conserving' - for a "conserving" charge method (each grid point
!                       is counted strictly once);
! method='nonconserv' - for a "non-conserving" algorithm when we scan the
!                       sphere rather than the UC so that each point
!                       may enter several times.
!......................................................................
      Yes_Do=.false.
      Yes_Obtain=.false.
1     iQuit=0
      write(*,*)'..............MENU for CHARGE ........................'
      write(*,*)'......... Change these parameters if necessary:.......'
      write(*,*)
      write(*,'(a51,i1)') &
        '   2. The number of (possibly) overlaping spheres: ',N_Sp
      write(*,'(a)')'   3. The center(s) of your sphere(s):'
      do i=1,N_Sp
        fCENTX=BCELL(1,1)*Point(1,i)+BCELL(1,2)*Point(2,i)+ &
             BCELL(1,3)*Point(3,i)
        fCENTY=BCELL(2,1)*Point(1,i)+BCELL(2,2)*Point(2,i)+ &
             BCELL(2,3)*Point(3,i)
        fCENTZ=BCELL(3,1)*Point(1,i)+BCELL(3,2)*Point(2,i)+ &
             BCELL(3,3)*Point(3,i)
        write(*,'(a2,i1,a7,f8.3,2(a1,f8.3),a9,f8.3,2(a1,f8.3),a1)') &
        ' [',i,'] A=> (',Point(1,i),',',Point(2,i),',',Point(3,i), &
        '), fr=> (', fCENTX,',',fCENTY,',',fCENTZ,')'
      end do

      write(*,'(a39,f10.5)') &
           '   4. The smallest radius (Angstroms): ', RadiusS
      write(*,'(a38,f10.5)') &
           '   5. The largest radius (Angstroms): ',RadiusL
      if(N_Sp.gt.1) write(*,'(a25,5(1x,f8.3))') &
               '  55. Shrinking factors: ',(Shrink(i),i=1,N_Sp)

      if(Nrad.eq.0) then
         iQuit=1
         write(*,'(a)') '   6. The number of points between '// &
                    'these radii: ... undefined ...'
      else
         write(*,'(a48,i5)') &
              '   6. The number of points between these radii: ',Nrad
         if(Nrad.eq.1) then
            dRad=RadiusS
         else
            dRad=(RadiusL-RadiusS)/(Nrad-1)
         end if
      end if

      if(Yes_Do) then
         write(*,'(a)') '   9. Perform calculation for the plotting: file '//filen &
                                             //' <= DONE!'
         write(*,'(a,f10.5)') '    98. Number of electrons to fit the sphere in: ',NumbEl
         if(Yes_Obtain) then
            write(*,'(a,f10.5)')'      >>>>> Radius = ',Rad_Fit
            write(*,'(a)') '    99. Obtain the radius for the given # of electrons:' &
                 //' <= DONE!'
         else
            write(*,'(a)') '    99. Obtain the radius for the given # of electrons:'
         end if
      else
         write(*,'(a)') '   9. Perform calculation for the plotting: file '//filen
      end if
      
      write(*,'(a)') '  10. Preview the dependence of charge versus Radius'
      write(*,'(a)') '  11. Create a PostScript file '//filen(:lenght)//'.ps'

      write(*,'(a)') &
                    '-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
      write(*,'(a)')'   7. Algorithm for the charge integration: <'//method//'>'
      if(method.eq.'conserving') then
         write(*,'(a)') '   8. X,Y,Z integration grid inside the sphere: IRRELEVANT'
      else
         write(*,'(a48,i5)') '   8. X,Y,Z integration grid inside the sphere: ',NRESOLs
         if(NRESOLs.le.1) iQuit=1
      end if
      if(way_res.eq.'%') then
          write(*,'(a)') '  Ne. Representation of results: through % to the density'
          cha='%'
      else
          write(*,'(a)') '  Ne. Representation of results: through number of electrons'
          cha='e'
      end if
      write(*,'(a)')'  Co. Show current atomic positions in fractional/Cartesian'
      write(*,'(a)')'------ L e a v e   t h e   m e n u -------------'
      write(*,'(a)')'   Q. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read (*,'(a)',err=1) cha2
!
![An]__________ choose the way how the coordinates are given
!
      IF(trim(cha2).eq.'An') THEN
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
![Ne]__________ choose a way to present the results
!
      ELSE IF(trim(cha2).eq.'Ne') THEN
         if(way_res.eq.'%') then
            way_res='e'
         else
            way_res='%'
         end if
!
![2]__________ give number of spheres
!
      ELSE IF(trim(cha2).eq.'2') THEN
 12      WRITE(*,*)'Give number of (ppossibly) overlapping spheres:'
         read(*,*,err=12) N_Sp
         if(N_Sp.gt.Num_Sph) then
           write(*,*)'IGNORED! It is too big (see <menu.inc> file)!'
           go to 12
         end if
         do i=1,N_Sp
           Shrink(i)=1.0
         end do
         Yes_Do=.false.
!
![3]__________ give central point on the plane
!
      ELSE IF(trim(cha2).eq.'3') THEN
         WRITE(*,*)'Give positions of your sphere(s)'
         do i=1,N_Sp
           write(*,*)'The sphere number ',i
           call givepoint(Point(1,i),Point(2,i),Point(3,i),angstr)
         end do
         Yes_Do=.false.
!
![4,5,55]__________ give radii
!
      ELSE IF(trim(cha2).eq.'4') THEN
 10      WRITE(*,*) 'Enter the smallest radius (in Angstroms):'
         READ(*,*,err=10) RadiusS
         if(RadiusS.lt.dzero) go to 10
         Yes_Do=.false.

      ELSE IF(trim(cha2).eq.'5') THEN
 11      WRITE(*,*) 'Enter the largest radius (in Angstroms):'
         READ(*,*,err=11) RadiusL
         if(RadiusL.lt.dzero) go to 11
         Yes_Do=.false.

      ELSE IF(trim(cha2).eq.'55' .and. N_Sp.gt.1) THEN
 13      WRITE(*,*) 'Give shrinking factors for every sphere:'
         READ(*,*,err=13) (Shrink(i),i=1,N_Sp)
         do i=1,N_Sp
           if(Shrink(i).lt.dzero) go to 13
         end do
         Yes_Do=.false.
!
![6]__________ number of points between RadiusS and RadiusL
!
      ELSE IF(trim(cha2).eq.'6') THEN
997      write(*,'(a31,f10.5,a4,f10.5)') 'Give the number of '// &
                  'different radii: '
         read(*,*,err=997) Nrad
         if(Nrad.lt.1) go to 997
         Yes_Do=.false.
!
![7]__________ algorithm: charge conserving or not
!
      ELSE IF(trim(cha2).eq.'7') THEN
         if(method.eq.'conserving') then
           method='nonconserv'
         else
           method='conserving'
         end if
         Yes_Do=.false.
!
![8]__________ number of grid points inside the sphere
!
      ELSE IF(trim(cha2).eq.'8') THEN
         if(method.eq.'conserving') then

           write(*,*)'This option is IRRELEVANT for this method!'
         else
 996       write(*,*)'Give this number:'
           read(*,*,err=996) NRESOLs
           if(NRESOLs.lt.2) go to 996
           Yes_Do=.false.
         end if
!
![9]__________ calculation
!   (loop over values of Radius from RadiusS till RadiusL)
!    NumbEl_max - maximal number of electrons found
!
      ELSE IF(trim(cha2).eq.'9') THEN
         if(iQuit.ne.0) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
!
!_________ if in number of electrons: unchanged
         factor_way=1.0
!_________ if in % to totdens - total density
         if(way_res.eq.'%') factor_way=NPLWV/totdens*100
!
         open(31,file=filen(:lenght),status='unknown',form='formatted')
         write(*,*)'The file '//filen//' has been opened for the PLOT.'
         write(*,*)'Working on the charge: writing to '// &
                                       filen(1:lenght)//' ...'
         NumbEl_max=dzero
         do i=1,Nrad
            Rad = dRad * (i-1) + RadiusS
            do j=1,N_Sp
               Radius(j)=Rad*Shrink(j)
            end do
            write(*,*)'..... Radius = ',Rad,' ......'
            rCharge=dzero
            if(method.eq.'nonconserv') then
!
!_________  charge "non-conserving" algorithm: since each point of the
!            UC may be met more than once, the charge is not properly
!            normalised. Scan a net of points inside the sphere of
!            Radius using NRESOLs and calculate the amount of charge
!            inside Radius
!
              DO ii=1,N_Sp
                 dX=Radius(ii)/(NRESOLs/2)
                 dV=dX*dX*dX
                 factor=dV/VOLC
                 if(Radius(ii).ne.dzero) then
                    do k1=-NRESOLs/2,NRESOLs/2
                       do k2=-NRESOLs/2,NRESOLs/2
                          do k3=-NRESOLs/2,NRESOLs/2
                             R(1)= dX*k1
                             R(2)= dX*k2
                             R(3)= dX*k3
                             aR=sqrt(R(1)*R(1)+R(2)*R(2)+R(3)*R(3))
                             if(aR.le.Radius(ii)) then
                                R(1)=R(1)+Point(1,ii)
                                R(2)=R(2)+Point(2,ii)
                                R(3)=R(3)+Point(3,ii)
                                
                                if(ii.ne.1) then
                                   do ij=1,ii-1
                                      do j1=-1,1
                                         do j2=-1,1
                                            do j3=-1,1
                                               x(1)=j1*DIRC(1,1)+j2*DIRC(2,1)+ &
                                                    j3*DIRC(3,1)+R(1)-Point(1,ij)
                                               x(2)=j1*DIRC(1,2)+j2*DIRC(2,2)+ &
                                                    j3*DIRC(3,2)+R(2)-Point(2,ij)
                                               x(3)=j1*DIRC(1,3)+j2*DIRC(2,3)+ &
                                                    j3*DIRC(3,3)+R(3)-Point(3,ij)
                                               aR=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
                                               if(aR.le.Radius(ij)) go to 21
                                            end do
                                         end do
                                      end do
                                   end do
                                end if
                                
                                call reducn(R,DIRC,BCELL)
                                call interpolate(R,BCELL,denval,grid)
                                rCharge = rCharge + denval*factor
                             end if
21                           continue
                          end do
                       end do
                    end do
                 end if
              END DO
              write(*,'(2(f10.5,5x))') Rad, rCharge*factor_way
              write(31,'(2(f10.5,5x))') Rad, rCharge*factor_way
              
           else
!
!____________ charge "conserving" algorithm: each point can be met
!             only once since we scan UC rather than the sphere
!
              iPnt=0
              do iZ=0,NGZ-1
                 do iY=0,NGY-1
                    do iX=0,NGX-1
!_____________ ask whether the point (iX,iY,iZ) is inside the sphere
                       do ii=1,N_Sp
                          call ask(iX,iY,iZ,Point(1,ii),Radius(ii),DIRC,iask)
                          if(iask.eq.'y') then
                             rCharge = rCharge + grid(iX+1,iY+1,iZ+1)
                             iPnt=iPnt+1
                             go to 22
                          end if
                       end do
22                     continue
                    end do
                 end do
              end do
              if(way_res.eq.'%') then
                 write(*,'(f10.5,5x,i10,5x,f10.5)') &
                      Rad, iPnt, rCharge/totdens*100
                 write(31,'(2(f10.5,5x))') Rad, rCharge/totdens*100
              else
                 write(*,'(f10.5,5x,i10,5x,f10.5)') &
                      Rad, iPnt, rCharge/NPLWV
                 write(31,'(2(f10.5,5x))') Rad, rCharge/NPLWV
              end if
              rCharge=rCharge/NPLWV
           end if
           if(rCharge.gt.NumbEl_max) NumbEl_max=rCharge
        end do
        close (31)
        write(*,*)'.... File '//filen(1:lenght)//' has been created! ....'
        write(*,*)'... The maximum # of electrons found = ',NumbEl_max
        Yes_Do=.true.
!
![98]__________ number of electrons to fit the sphere in
!
     ELSE IF(Yes_Do .and. trim(cha2).eq.'98') THEN
33      WRITE(*,*) 'Give the number of electrons:'
        READ(*,*,err=33) NumbEl
        if(NumbEl.le.dzero .or. NumbEl.gt.NumbEl_max) go to 33
        Yes_Obtain=.false.
!
![99]__________ number of electrons to fit the sphere in (using
!           linear interpolation)
!
      ELSE IF(Yes_Do .and. trim(cha2).eq.'99') THEN
        rad1=dzero
        c1=dzero
        open(31,file=filen(:lenght),status='old', &
                                         form='formatted',err=150)
 35     read(31,*,end=45) rad,charg
        if(NumbEl.ge.c1.and.NumbEl.le.charg) then
          rad2=rad
          c2=charg
          go to 40
        else
          rad1=rad
          c1=charg
        end if
        go to 35
 40     Rad_Fit=rad1+(NumbEl-c1)/(c2-c1)*(rad2-rad1)
        Yes_Obtain=.true.
        close (31)
        go to 1
 45     write(*,*)'ERROR! Cannot find the fit!'
        close (31)
        go to 1
150     write(*,*)'FATAL! Cannot open '//filen(:lenght)//' file!'
        go to 1
!
![10]__________ preview the file just created
!
      ELSE IF(trim(cha2).eq.'10') THEN
         if(Yes_Do) then
           call Plot1(filen,lenght,Title,title_pl, &
              'Sphere radius (A)   ', &
              'Charge involved ('//cha//') ',  'Screen', 33,0, &
                 'N',.false.,dzero,dzero)
         else
           write(*,*) &
             'IGNORED! You have to accomplish the item 8 first!'
         end if
!
![11]__________ create a PostScript file of the plot
!
      ELSE IF(trim(cha2).eq.'11') THEN
         if(Yes_Do) then
           write(*,*)'Give the title:'
           read(*,'(a)') Title
           call Plot1(filen,lenght,Title,title_pl, &
              'Sphere radius (A)   ', &
              'Charge involved ('//cha//') ', 'Postsc', 33,0, &
                 'N',.false.,dzero,dzero)
         else
           write(*,*) &
             'IGNORED! You have to accomplish the item 8 first!'
         end if
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         call show_atoms()
!
![Q]__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
end subroutine charge_sph


subroutine max_dens(grid,DIRC,BCELL,VOLC,totdens)
!....................................................................
!     Points in the cell with the maximal values of the density are
!  generated by scaning the density in grid(iX,iY,iZ)
!....................................................................
!  Npeak - number of peaks found
!  Then, for every i=1,...,Npeak:
!  x(3,i)  - positions of the peaks in Cartesian coord.
!  High(i) - values on peaks tip;
!  Rad(i)  - radius of a sphere for the peak i which separates it
!            from the rest of the cell
!  rCharg(i) - charge (in %) involved in the sphere associated with
!             the peak i
!....................................................................
use param
use explore
implicit none
real*8 :: multcon=1.0
real*8 absden(3),xcoord(3),width(3),VOLC,totdens
real*8 GRID(NGX,NGY,NGZ),DIRC(3,3),BCELL(3,3),denval,rCharge,R(3)
character iask,Title*50,answer,show
integer iPeak,iX,iY,iZ,i,iXr,iYr,iZr,k2,j,iPnt
data Title/'                                                  '/
!
!______ lenghts along cell axes
!
   width(1)=sqrt(DIRC(1,1)**2+DIRC(1,2)**2+DIRC(1,3)**2)
   width(2)=sqrt(DIRC(2,1)**2+DIRC(2,2)**2+DIRC(2,3)**2)
   width(3)=sqrt(DIRC(3,1)**2+DIRC(3,2)**2+DIRC(3,3)**2)
   write(*,*) (width(j),j=1,3)
!
!_____ dialog: results in % or in # of electrons
!
   show='e'
   write(*,*) &
           'Results in % of charge or in # of electrons (%/any cha)?'
   read(*,'(a)') answer
   if(answer.eq.'%') show='%'
   if(show.eq.'%') then
      write(*,'(a)')'MIND: Results will be given via % of charge'
   else
      write(*,'(a)')'MIND: Results will be given via # of electrons'
   end if
!
!.....................................................................
!.................  iPeak - is used to accumulate peaks ..............
!.....................................................................
!
   iPeak = 0
1  iPeak=iPeak+1
   if(iPeak.gt.Mpeak) then
      write(*,*)'WARNING! No more peaks can be generated!'
      write(*,*)'Increase Mpeak in the file ~/TOOLS/explore.inc'
      return
   end if
   write(*,*)'........Searching for the ',iPeak,' peak center ......'
!
!............. loop over grid points: look for the point with the
!       highest density, High(iPeak), for the current peak, iPeak
!
   High(iPeak)=0.0
   do iZ=0,NGZ-1
      do iY=0,NGY-1
         do 10 iX=0,NGX-1
!
!_________ ask whether the point (iX,iY,iZ) is inside the spheres
!          associated with other peaks; if so, reject the point
!
            if(iPeak.gt.1) then
               do i=1,iPeak-1
                  call ask(iX,iY,iZ,x(1,i),Rad(i),DIRC,iask)
                  if(iask.eq.'y') go to 10
               end do
            end if
!
!_________ check the density at the point
!
            if(grid(iX+1,iY+1,iZ+1).gt.High(iPeak)) then
               High(iPeak)=grid(iX+1,iY+1,iZ+1)
               iXr=iX
               iYr=iY
               iZr=iZ
            end if
10       end do
      end do
   end do
   High(iPeak)=High(iPeak)*multcon/VOLC
!
!....... give the point (iXr,iYr,iZr) in Cartesian coordinates
!
   x(1,iPeak)=iXr*DIRC(1,1)/NGX+iYr*DIRC(2,1)/NGY + &
        iZr*DIRC(3,1)/NGZ
   x(2,iPeak)=iXr*DIRC(1,2)/NGX+iYr*DIRC(2,2)/NGY + &
        iZr*DIRC(3,2)/NGZ
   x(3,iPeak)=iXr*DIRC(1,3)/NGX+iYr*DIRC(2,3)/NGY + &
        iZr*DIRC(3,3)/NGZ
   write(*,*)'The point with the highest density ',High(iPeak),' is:'
   write(*,'(a16,3(i5,a1))') &
                ' fractional => (',iXr,',',iYr,',',iZr,')'
   write(*,'(a16,3(f10.5,a1))') &
        '  Cartesian => (',x(1,iPeak),',',x(2,iPeak),',',x(3,iPeak),')'
   write(*,*)'Press ENTER to proceed ...'
   read(*,*)
!
!....... explore the density around and plot it ....................
!
!__________ scanning along axes to choose the radius, Rad(iPeak)
!
   open(32,file='test.dat',status='unknown',form='formatted')
   write(*,*)'The file test.dat has been opened to preview.'
   do k2=-NRESOL/2,NRESOL/2
      do j=1,3
!______ along j-th axis
         R(1)=x(1,iPeak)+k2*DIRC(j,1)/nresol
         R(2)=x(2,iPeak)+k2*DIRC(j,2)/nresol
         R(3)=x(3,iPeak)+k2*DIRC(j,3)/nresol
         call reducn(R,DIRC,BCELL)
         call interpolate(R,BCELL,denval,grid)
         xcoord(j)=k2*width(j)/nresol
         absden(j)=(denval*multcon)/VOLC
      end do
      write(32,'(6(e12.6,1x))') (xcoord(j),absden(j),j=1,3)
   end do
   close (32)
!
!........... show density along axes on one plot simultaneously
!
   call Plot_Explr('test.dat',8,Title, &
        'Coordinate (A)      ', &
        'Charge density      ',  'Screen', 33)
!
!........... choose the radius around this point
!
20 write(*,*)'Choose the radius (in Angstrem):'
   read(*,*,err=20) Rad(iPeak)
!
!.......... loop over grid points: match all points that fit this
!          radius and calculate the amount of charge inside sphere
!          (in % to totdens - total density);
!
   write(*,*)'Calculating the charge involved ...'
   rCharge=0.0
   iPnt=0
   do iZ=0,NGZ-1
      do iY=0,NGY-1
         do 30 iX=0,NGX-1
!
!_________ ask whether the point (iX,iY,iZ) is inside the current
!          sphere and outside all previous ones
!
            do i=1,iPeak
               call ask(iX,iY,iZ,x(1,i),Rad(i),DIRC,iask)
               if(i.lt.iPeak .and. iask.eq.'y') go to 30
            end do
            if(iask.eq.'y') then
               rCharge = rCharge + grid(iX+1,iY+1,iZ+1)
               iPnt=iPnt+1
            end if
30       end do
      end do
   end do

   if(show.eq.'%') then
!________________ if in units of % of the total electronic charge
      rCharg(iPeak)=rCharge/totdens*100
      write(*,'(a9,f6.3,a26,f6.3)') 'There is ',rCharg(iPeak), &
           '% of charge inside radius=',Rad(iPeak)
   else
!_________________ if in units of the # of electrons
      rCharg(iPeak)=rCharge/NPLWV
      write(*,'(a9,f6.3,a26,f6.3)') 'There is ',rCharg(iPeak), &
           'electrons inside radius=',Rad(iPeak)
   end if
!
   write(*,*)'Try another radius (y,Y/other char)?'
   read(*,'(a)') answer
   if(answer.eq.'y'.or.answer.eq.'Y') go to 20
!
!....... finishing by printing
   write(*,*)'........> Final distribution of the density:'
   if(show.eq.'%') then
      write(*,*) &
           'Peak  Hight        _________(x,y,z)_________'// &
           '      Charge(%)    Radius'
   else
      write(*,*) &
           'Peak  Hight        _________(x,y,z)_________'// &
           '      Charge(e)    Radius'
   end if
   do i=1,iPeak
      write(*,'(i2,1x,f10.5,a2,3(f10.5,a1),2(1x,f10.5))') &
           i,High(i),' (',x(1,i),',',x(2,i),',',x(3,i),')',rCharg(i),Rad(i)
   end do
!
!........ whether to repeat
   write(*,*) &
        'Would you like to repeat exploration (y,Y/other char)?'
   read(*,'(a)') answer
   if(answer.eq.'y'.or.answer.eq.'Y') go to 1
   open(32,file='explor.dat',status='unknown',form='formatted')
   write(32,*)'........> Final distribution of the density:'
   write(32,*) &
        'Peak  Hight        _________(x,y,z)_________'// &
        '      Charge(%)    Radius'
   do i=1,iPeak
      write(32,'(i2,1x,f10.5,a2,3(f10.5,a1),2(1x,f10.5))') &
           i,High(i),' (',x(1,i),',',x(2,i),',',x(3,i),')',rCharg(i),Rad(i)
   end do
   close (32)
   write(*,*)'WARNING! The table with peaks - see explor.dat file!'
end subroutine max_dens
