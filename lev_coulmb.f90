subroutine lev_coulmb()
!.........................................................................
!  Calculates the Coulomb potential (in Ev unless factor=/=1.0) from a
! point-ion lattice in an arbitrary point inside the UC, as well as along
! a line or in a plane.
!.........................................................................
! Charge(species) - ionic charge on each species;
! Q(ion) - ionic charge on each ion in the cell.
! QradI(ion) - atomic radius for the ion used to eliminate
!              a singularity in the Madelung potential within the
!              sphere of this radius near the ion (in menu.inc);
! QradS(species) - atomic radius for the species (is used to build up QradI)
!.........................................................................
! DIRC - direct lattice vectors
! RECC - reciprocal lattice vectors (with 2*pi)
! VOLC - unit cell volume
!..................................................................
! BCELL - reciprocal lattice vectors (without 2*pi)
!..................................................................
use param
use atoms
use code
use menu
implicit none
real*8 a(3),x(3),tot,dist,distm,dd
real*8,dimension(:),allocatable :: qCharge,Q,QradS
character cha,cha2*2,filen*12
real*8 :: pi=3.141592654,EPSew=1.0E-10,factor=1.0,tiny=0.00001,gEwald,EPSx
integer :: iCharge=0,iChargeS=0,mCharge=0,iOver=0
integer i,j,ijk,lenght,iQuit, iCheck,i1,item,ion
integer iOverlp,ion1,k

!................... memory

      allocate(qCharge(NIONS))
      allocate(Q(NIONS))
      allocate(QradS(NSPEC))

!........... do BCELL = RECC/(2*pi) - it is used in transforming
!            grid point coordinates
!
       do i=1,3
          do j=1,3
             BCELL(i,j)=RECC(i,j)/(2*pi)
          end do
       end do
!
!..........  The "best" Ewald constant gEwald is estimated as:
!
       call best_Ewald(DIRC,BCELL,gEwald)
!
!...................... RUN a MENU here ...........................
!....... set up charges, precision, the point of interest (Pnt)
!        and a multiplication factor
!..................................................................
!
!....................................................................
!............ General part: let us plot just once ...................
!....................................................................
!.... ijk - counts different cycles of calculations (not more than 9),
!           i.e. different plots
!
      ijk=1
!
!............ name of the file for the output
!
2     if(ijk.le.9) then
         write(cha,'(i1)') ijk
         filen='mad.dat.'//cha
         lenght=9
      else if(ijk.le.99) then
         write(cha2,'(i2)') ijk
         filen='mad.dat.'//cha2
         lenght=10
      else
         write(*,*)'LEV_COULMB: You cannot trial my patience so much!'
         go to 200
      end if
      iQuit=0
      iCheck=0
!
!____________ choose between a line, plane or charge
!
      write(*,*)'..............MENU for Madelung potential ...........'
      write(*,*)'......... Change these parameters if necessary: .....'
      write(*,*)
      write(*,'(a33,i2)')'     NUMBER OF THE CURRENT PLOT: ',ijk
      write(*,*)
      write(*,'(a)')'   0. Coordinates are specified in: '//angstr
      if(mCharge.eq.0) then
         write(*,'(a)') '   1. Charges on all species are the same: YES'
      else
         write(*,'(a)') '   1. Charges on all species are the same: NO'
      end if
      if(iChargeS.eq.0) then
         if(mCharge.eq.0) then
            iQuit=1
            write(*,'(a)') '   2. Charges on species: ....... undefined .......'
         else
            write(*,'(a)') '   2. Charges on species (not used): UNKNOWN'
         end if
      else
         write(*,'(a)') '   2. Charges on species:   KNOWN'
      end if
      if(iCharge.eq.0) then
         iQuit=1
         write(*,'(a)') '   3. Charges on ions: ....... undefined .......'
      else
         write(*,'(a)') '   3. Charges on ions:   KNOWN'
      end if
      write(*,'(a35,e12.6)')'   4. Precision of the summations: ',EPSew
      EPSx=-log(EPSew)
      write(*,'(a29,e12.6)')'   5. Multiplication factor: ',factor

      write(*,'(a)') &
      '   6. Calculate the potential at a single point (Volts*factor)'
      if(iOver.eq.0) then
         write(*,'(a)') &
              '   7. Atomic radii (only line/plane): ..... undefined .....'
      else
         write(*,'(a)') &
              '   7. Atomic radii (only line/plane) are given for species:'
         do i=1,NSPEC,6
            i1=i+5
            if(i1.gt.NSPEC) i1=NSPEC
            write(*,'(7x,6(f8.5,1x))') (QradS(j),j=i,i1)
         end do
      end if

      write(*,'(a)') '   8. Calculate the potential along a line (in eV*factor)'
      write(*,'(a)') '   9. Calculate the potential in a plane (in eV*factor)'
      write(*,'(a)') '  10. Quit'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ choose the way how the coordinates are given
!
      if(item.eq.0) then
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ choose the method: individually for all ions or through
!           the species
!
      else if(item.eq.1) then
         if(mCharge.eq.0) then
            mCharge=1
         else
            mCharge=0
         end if
!
!__________ specify charges on species and assign these charges to
!           all ions in Q(ion):
!
      else if(item.eq.2) then
 20      if(mCharge.eq.0) then
             write(*,'(a22,i2,a9)') 'Give charges for each ',NSPEC,' species:'
             write(*,'(10(a,x))') (Species(i),i=1,NSPEC)
             read(*,*,err=20) (qCharge(i),i=1,NSPEC)
             iChargeS=1
             write(*,*)'______ charges on ions: ________'
             ion=0
             tot=0.0
             do i=1,NSPEC
                do j=1,NspN(i)
                   ion = ion +1
                   Q(ion)=qCharge(i)
                   tot=tot+Q(ion)
                   write(*,'(a,a,i2,a,i3,a,f10.5)') &
                                Species(i),' #=',j,' Q(',ion,')=',Q(ion)
                end do
             end do
             write(*,*)'________________________________'
             if(abs(tot).gt.tiny) then
                write(*,*)'ERROR! Your unit cell is not neutral!'
                write(*,*)'Total charge = ', tot
                go to 20
             end if
             iCharge=1
          else
             write(*,*)'WANRNING! You do NOT need them unless the option 2 is YES!'
          end if
!
!__________ specify charges for all ions in the cell individually
!
       else if(item.eq.3) then
29        if(mCharge.eq.1) then
              write(*,*) 'Give ionic charges individually (in the order):'
              ion=1
              do i=1,NSPEC
30               write(*,'(a17,i2,a2,i2,a6)') &
                      '____ The species ',i,': ',NspN(i),' ions:'
                 read(*,*,err=30) (Q(k),k=ion,ion+NspN(i)-1)
                 ion=ion+NspN(i)
              end do
              iCharge=1
              write(*,*)'______ charges on ions: ________'
              ion=0
              tot=0.0
              do i=1,NSPEC
                 do j=1,NspN(i)
                    ion = ion +1
                    write(*,'(a4,i2,a3,i2,a3,i3,a2,f10.5)') &
                         ' sp=',i,' #=',j,' Q(',ion,')=',Q(ion)
                    tot=tot+Q(ion)
                 end do
              end do
              write(*,*)'________________________________'
              if(abs(tot).gt.tiny) then
                 write(*,*)'ERROR! Your unit cell is not neutral!'
                 go to 29
              end if
           else
              write(*,*) 'WANRNING! You do NOT need them unless the option 2 is NO!'
           end if
!
!__________ specify the precision of the summations in the Ewald method
!
      else if(item.eq.4) then
40       write(*,*)'Give the precision:'
         read(*,*,err=40) EPSew
         if(EPSew .lt. 1.0e-13) go to 40
!
!__________ specify the multiplication factor
!
      else if(item.eq.5) then
45       write(*,*)'Give the multiplication factor:'
         read(*,*,err=45) factor
!
!__________ Madelung potential at a point
!
      else if(item.eq.6) then
         if(iQuit.eq.1) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 2
         end if
         call pointM(Q,factor,gEwald,EPSx)
         ijk=ijk+1
!
!__________ specify atomic radii for every species to eliminate
!           singularities near atoms while making plots. Then we check,
!           whether the spheres overlap: they should not.
!
      else if(item.eq.7) then
46       write(*,*)'Give atomic radii (in Angstroms) for every ', &
              NSPEC,' species:'
         read(*,*,err=46) (QradS(j),j=1,NSPEC)
         ion=0
         do i=1,NSPEC
            do j=1,NspN(i)
               ion = ion +1
               QradI(ion)=QradS(i)
               write(*,'(a4,i2,a3,i2,a7,i3,a2,f10.5)') &
                    ' sp=',i,' #=',j,' QradI(',ion,')=',QradI(ion)
            end do
         end do
         write(*,*)'Checking if the spheres overlap ...'
         iOverlp=0
         do ion=1,NIONS
            do 50 ion1=ion,NIONS
               if(ion.eq.ion1) go to 50
               a(1)=TI(1,ion)-TI(1,ion1)
               a(2)=TI(2,ion)-TI(2,ion1)
               a(3)=TI(3,ion)-TI(3,ion1)
               distm=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
               do i=-1,1
                  do j=-1,1
                     do k=-1,1
                        x(1) = i*DIRC(1,1)+j*DIRC(2,1)+k*DIRC(3,1) + a(1)
                        x(2) = i*DIRC(1,2)+j*DIRC(2,2)+k*DIRC(3,2) + a(2)
                        x(3) = i*DIRC(1,3)+j*DIRC(2,3)+k*DIRC(3,3) + a(3)
                        dist=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
                        if(dist.lt.distm) distm=dist
                     end do
                  end do
               end do
               dd=distm-(QradI(ion)+QradI(ion1))
               if(dd.lt.0.0) then
                  iOverlp=1
                  write(*,'(a24,i3,a5,i3,a12,f5.2,a2)') &
                       'ERROR! Spheres of atoms ',ion,' and ',ion1, &
                       ' overlap by ',abs(dd),' !'
               else
                  iCheck=1
                  iOver=1
               end if
50          end do
         end do
         if(iOverlp.eq.1) go to 46
!
!__________ Madelung potential along a line
!
      else if(item.eq.8) then
         if(iOver.eq.0 .or. iQuit.eq.1) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 2
         end if
         call lineM(Q,filen,lenght,factor,gEwald,EPSx)
         ijk=ijk+1
!
!__________ Madelung potential in a plane
!
      else if(item.eq.9) then
         if(iOver.eq.0 .or. iQuit.eq.1) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 2
         end if
         call planeM(Q,filen,lenght,factor,gEwald,EPSx)
         ijk=ijk+1
!
!__________ Quit or skip
!
      else if(item.eq.10) then
         go to 200
      else
         go to 100
      end if
      go to 2
100   write(*,*)'Incorrect item number! Try again!'
      go to 2
!
!............. finish
200   deallocate(qCharge)
      deallocate(Q)
      deallocate(QradS)
end subroutine lev_coulmb

subroutine pointM(Q,factor,gEwald,EPSx)
!...................................................................
!  Madelung potential at a single point
!...................................................................
use param
use menu
use atoms 
implicit none
real*8 :: Q(NIONS),vMad,factor,gEwald,EPSx,bCENTX,bCENTY,bCENTZ,pot
real*8, dimension(3) :: Pnt=(/0.0,0.0,0.0/),fPnt(3)
integer :: iCoord=0,iQuit,item,iCheck
!
1     iQuit=0
      write(*,*)'..............MENU for Madelung (Point) .............'
      write(*,*)'......... Change these parameters if necessary: .....'
      write(*,*)
      write(*,'(a)')'   0. Coordinates are specified in: '//angstr
      if(iCoord.eq.0) then
         iQuit=1
         write(*,'(a)') '   1. Point of interest: ....... undefined .......'
      else
         write(*,'(a38,f10.5,2(a1,f10.5),a1)') &
              '   1. Point of interest (Angstroms): (', &
              Pnt(1),',',Pnt(2),',',Pnt(3),')'
         fPnt(1)=BCELL(1,1)*Pnt(1)+BCELL(1,2)*Pnt(2)+BCELL(1,3)*Pnt(3)
         fPnt(2)=BCELL(2,1)*Pnt(1)+BCELL(2,2)*Pnt(2)+BCELL(2,3)*Pnt(3)
         fPnt(3)=BCELL(3,1)*Pnt(1)+BCELL(3,2)*Pnt(2)+BCELL(3,3)*Pnt(3)
         write(*,'(a38,f10.5,2(a1,f10.5),a1)') &
         '     Point of interest (fractional): (', &
                           fPnt(1),',',fPnt(2),',',fPnt(3),')'
      end if
      write(*,'(a)')'   2. Calculate the potential'
      write(*,'(a)')'   3. Quit'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ choose the way how the coordinates are given
!
      if(item.eq.0) then
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ give the point
!
      else if(item.eq.1) then
         call givepoint(Pnt(1),Pnt(2),Pnt(3),angstr)
         iCoord=1
!
!__________ calculate the potential
!
      else if(item.eq.2) then
         if(iQuit.eq.1) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
         write(*,*)'Please, wait ...'
         pot=vMad(Pnt,Q,gEwald,EPSx,iCheck)
         write(*,'(a,e12.6)')'Madelung potential: ',pot*factor
!
!__________ Quit
!
      else if(item.eq.3) then
         return
      else
         go to 100
      end if
      go to 1
100   write(*,*)'ERROR! Try again!'
      go to 1
end subroutine pointM

subroutine lineM(Q,filen,lenght,factor,gEwald,EPSx)
!....................................................................
!  Line Calculation of the Madelung potential
!  31 - unit number for the file filen(1:lenght) with output data.
!  iQuit = 0 - all parameters are properly defined; can plot
!          1 - there are undefined parameters; cannot plot
!....................................................................
! Note: the potential of the atom within its radius is calculated
! properly so that there is no discontinuity neither in the potential
! nor in its derivative at the sphere surface.
!....................................................................
use param
use menu
use atoms
implicit none
real*8 :: R(3),Q(NIONS),tiny=0.00001,dzero=0.0,vMad,factor,gEwald,EPSx
real*8 :: fCENTX,fCENTY,fCENTZ,xcoord,a,absden,bCENTX,bCENTY,bCENTZ
character filen*12,Title*50,title_pl*7
integer lenght,iQuit,item,k2,lenght3,iCheck
data Title/'                                                  '/
data title_pl/'       '/
logical Yes_Do
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
      write(*,'(a)')'   0. Coordinates are specified in: '//angstr
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
      '   2. Vector along the line: (',vers0x,',',vers0y,',',vers0z,')'
      end if

      write(*,'(a44,f10.5)') '   3. Lenght along the line (in Angstroms): ',width1

      write(*,'(a)')'   4. Parameters for the plotting'

      if(Yes_Do) then
        write(*,'(a)') '   5. Perform calculation of the potential: file '//filen &
                                             //' <= DONE!'
      else
        write(*,'(a)') '   5. Perform calculation of the potential: file '//filen
      end if

      write(*,'(a)')'   6. Preview the potential'
      write(*,'(a)')'   7. Create a postscript file '// &
                                filen(:lenght)//'.ps for the plot'
      write(*,'(a)')'   8. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'

      read(*,*,err=100) item
!
!__________ choose a way how the coordinates are given
!
      if(item.eq.0) then
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ give starting point for the line
!
      else if(item.eq.1) then
         if(angstr.eq.'<AtomNumber>') &
              write(*,*)'Specify the 1st atom to be started from.'
         call givepoint(aCENTX,aCENTY,aCENTZ,angstr)
         Yes_Do=.false.
!
!__________ give a vector along the line
!
      else if(item.eq.2) then
         if(angstr.eq.'<AtomNumber>') then
            write(*,*)'Specify the 2nd atom to be connected with.'
            call givepoint(bCENTX,bCENTY,bCENTZ,angstr)
            vers0x=bCENTX-aCENTX
            vers0y=bCENTY-aCENTY
            vers0z=bCENTZ-aCENTZ
            WIDTH1=sqrt( vers0x**2+vers0y**2+vers0z**2 )
         else
7           write(*,*)'Give a vector (x,y,z) along your line:'
            read (*,*,err=7)  vers0x, vers0y, vers0z
         end if
         call normalize(vers0x,vers0y,vers0z)
         Yes_Do=.false.
         !
!__________ give length along the line
!
      else if(item.eq.3) then
10       write(*,*) 'Enter length (in Angstroms):'
         read(*,*,err=10) width1
         if(width1.lt.dzero) go to 10
         Yes_Do=.false.
!
!__________ give the resolution in either direction, chop values
!  and a multiplication factor for the potential, etc.
!
      else if(item.eq.4) then
         multcon=factor
         call choose1()
         Yes_Do=.false.
!
!__________ perform calculation of the potential along the line;
!
      else if(item.eq.5) then
         if(iQuit.ne.0) then
            write(*,*)'ERROR! You still have undefined parameters!'
            go to 1
         end if
         write(*,*)'Please, wait ...'
         open(31,file=filen(:lenght),status='unknown',form='formatted')
         write(*,*)'The file '//filen//' has been opened ...'
         write (*,*)'Writing to the file '//filen(1:lenght)//' ...'
         do K2=0,NRESOL
            R(1)=acentx+k2*vers0x*width1/nresol
            R(2)=acenty+k2*vers0y*width1/nresol
            R(3)=acentz+k2*vers0z*width1/nresol
            call reducn(R,DIRC,BCELL)
            absden=vMad(R,Q,gEwald,EPSx,iCheck)
            xcoord=k2*width1/nresol
            if(lochop.ne.hichop) then
               if(absden.gt.hichop) then
                  absden=hichop
               else if(absden.lt.lochop) then
                  absden=lochop
               end if
            end if
            write(31,*) xcoord,absden*multcon
         end do
         close (31)
         write(*,*)'.... File '//filen(1:lenght)//' has been created! ....'
         Yes_Do=.true.
!
!__________ preview the file just created
!
      else if(item.eq.6) then
         if(Yes_Do) then
            call Plot1(filen,lenght,Title,title_pl, &
                 'Coordinate (A)      ', &
                 'Coulomb potential   ',  'Screen', 33,0, &
                 'N',.false.,dzero,dzero)
         else
            write(*,*) 'IGNORED! You have to accomplish the item 5 first!'
         end if
         !
!__________ create a PostScript file of the plot
!
      else if(item.eq.7) then
         if(Yes_Do) then
            write(*,*)'Give the title:'
            read(*,'(a)') Title
            call Plot1(filen,lenght,Title,title_pl, &
                'Coordinate (A)      ', &
                'Coulomb potential   ', 'Postsc', 33,0, &
                'N',.false.,dzero,dzero)
         else
            write(*,*) 'IGNORED! You have to accomplish the item 5 first!'
         end if
!
      else if(item.eq.8) then
         return
      else
         go to 100
      end if
      go to 1
100   write(*,*)'ERROR! Try again!'
      go to 1
end subroutine lineM

subroutine planeM(Q,filen,lenght,factor,gEwald,EPSx)
!....................................................................
!  Plane Calculation of the Madelung potential.
!  nfile - unit number for the file filen(1:lenght) with output data.
!  iQuit = 0 - all parameters are properly defined; can plot
!          1 - there are undefined parameters; cannot plot
!....................................................................
! Note: the potential of the atom within its radius is calculated
! properly so that there is no discontinuity neither in the potential
! nor in its derivative at the sphere surface.
!....................................................................
use param
use menu
use atoms
implicit none
real*8 pA(2),pB(2),pC(2),R(3),Q(NIONS),vMad,factor,gEwald,EPSx,absden
integer lenght,iQuit,item,k3,k2,k1,lenght3,iCheck
character filen*12, Title*50
data Title/'                                                  '/
real*8 :: tiny=0.00001,dzero=0.0,a,fCENTX,fCENTY,fCENTZ,xcoord,ycoord
logical Yes_Do
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
      write(*,'(a)')'   0. Coordinates are specified in: '//angstr
      a=vers1x*vers1x + vers1y*vers1y + vers1z*vers1z
      if(a.lt.tiny) then
         iQuit=1
         write(*,'(a)') ' / 1. Normal vector to the plane: ....... undefined .......'
         write(*,'(a)') ' |    X1 vector in the plane:     ....... undefined .......'
         write(*,'(a)') ' |    Y1 vector in the plane:     ....... undefined .......'
      else
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
         write(*,'(a)') ' \\ 2. The plane has been specified by 3 points: NO'
      else if(icase.eq.2) then
         write(*,'(a)') ' \\ 2. The plane has been specified by 3 points: YES'
      end if
      if(icase.eq.2.and.icase1.eq.1) then
         write(*,'(a)') '   3. Central point => the center of the triangle: NO'
      else if(icase.eq.2.and.icase1.eq.2) then
         write(*,'(a)') '   3. Central point => the center of the triangle: YES'
         aCENTX=(Ra(1)+Rb(1)+Rc(1))/3.
         aCENTY=(Ra(2)+Rb(2)+Rc(2))/3.
         aCENTZ=(Ra(3)+Rb(3)+Rc(3))/3.
         central_p=.true.
      end if
      if(central_p) then
         if(icase.eq.2.and.icase1.eq.2) then
            write(*,'(a)') '      Central point on the plane: '
         else
            write(*,'(a)') '   4. Central point on the plane: '
         end if
         write(*,'(a29,f10.5,2(a1,f10.5),a1)') &
              '           in Angstroms  => (',aCENTX,',',aCENTY,',',aCENTZ,')'
         fCENTX=BCELL(1,1)*aCENTX+BCELL(1,2)*aCENTY+BCELL(1,3)*aCENTZ
         fCENTY=BCELL(2,1)*aCENTX+BCELL(2,2)*aCENTY+BCELL(2,3)*aCENTZ
         fCENTZ=BCELL(3,1)*aCENTX+BCELL(3,2)*aCENTY+BCELL(3,3)*aCENTZ
         write(*,'(a29,f10.5,2(a1,f10.5),a1)')'           in fractional => (', &
                             fCENTX,',',fCENTY,',',fCENTZ,')'
         if(icase.eq.2) then
            write(*,'(11x,(a))') 'The reference points A,B,C in (X1,Y1) are given as:'
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
13          format(15x,a5,f10.5,a1,f10.5,a1)
         end if
      else
         write(*,'(a)') '   4. Central point on the plane: ....... undefined .......'
      end if
      write(*,'(a39,f10.5)') '   5. Width along X1 axis (Angstroms): ',width1
      write(*,'(a39,f10.5)') '   6. Width along Y1 axis (Angstroms): ',width2
      write(*,'(a)')'   7. Parameters for the plotting'

      write(*,'(a)')'   8. Preview the potential'
      if(Yes_Do) then
         write(*,'(a)')'   9. Perform calculation for the potential: file '//filen &
                                             //' <= DONE!'
      else
        write(*,'(a)') '   9. Perform calculation for the potential: file '//filen
      end if

      write(*,'(a)')'  10. Return to the previous menu'
      write(*,*)
      write(*,*)'------> Choose the item and press ENTER:'
      read(*,*,err=100) item
!
!__________ choose the way how the coordinates are given
!
      if(item.eq.0) then
         if(angstr.eq.'<Fractional>') then
            angstr='<Angstroms> '
         else if(angstr.eq.'<Angstroms> ') then
            angstr='<AtomNumber>'
         else if(angstr.eq.'<AtomNumber>') then
            angstr='<Fractional>'
         end if
!
!__________ give a normal vector to the plane and generate two
!           others in the plane
!
      else if(item.eq.1) then
         icase=1
         call vector3(DIRC)
         Yes_Do=.false.
!
!__________ specify the plane by 3 points
!
      else if(item.eq.2) then
         icase=2
         if(NIONS.le.2 .and. angstr.eq.'<AtomNumber>' ) then
            write(*,*)'ERROR! Not enough atoms for this option!'
            write(*,*)'Change to <Fractional> or <Angstroms> using 0'
            write(*,*)'Hit ENTER when ready ...'
            read(*,*) 
         else
            call vector3(DIRC)
         end if
         Yes_Do=.false.
!
!__________ give a method to choose the central point on the plane
!           in the case of 3 points (icase=2)
!
      else if(item.eq.3) then
         if(icase.eq.2.and.icase1.eq.1) then
            icase1=2
            aCENTX=(Ra(1)+Rb(1)+Rc(1))/3.
            aCENTY=(Ra(2)+Rb(2)+Rc(2))/3.
            aCENTZ=(Ra(3)+Rb(3)+Rc(3))/3.
            central_p=.true.
            Yes_Do=.false.
         else if(icase.eq.2.and.icase1.eq.2) then
            icase1=1
         end if
!
!__________ give central point on the plane in a general way
!
      else if(item.eq.4) then
         if(iQuit.eq.1) then
            write(*,*)'ERROR! You must acomplish the item 1 first!'
         else
            call centralP(DIRC)
            central_p=.true.
            Yes_Do=.false.
         end if
!
!__________ give length along the X1,Y1 axes
!
      else if(item.eq.5) then
10       write(*,*) 'Enter length along X1 axis (in Angstroms):'
         read(*,*,err=10) width1
         if(width1.lt.dzero) go to 10
         Yes_Do=.false.
      else if(item.eq.6) then
 11      write(*,*) 'Enter length along Y1 axis (in Angstroms):'
         read(*,*,err=11) width2
         if(width2.lt.dzero) go to 11
         Yes_Do=.false.
!
!__________ give the resolution in either direction, chop values
!  and a multiplication factor for the density, etc.
!
      else if(item.eq.7) then
         multcon=factor
         call choose3()
         Yes_Do=.false.
!
!__________ preview
!
      else if(item.eq.8) then
         if(iQuit.eq.0.and.central_p) then
            open(32,file='test.dat',status='unknown',form='formatted')
            write(*,*)'The file test.dat has been opened to preview.'
            write(*,*)'Working on previewing. Please, wait ...'
            DO K3=-NRESOL_PRV/2,NRESOL_PRV/2
               DO K2=-NRESOL_PRV/2,NRESOL_PRV/2
                  R(1)=acentx+(k2*vers2x*width1 + &
                       k3*vers3x*width2)/nresol_prv
                  R(2)=acenty+(k2*vers2y*width1 + &
                       k3*vers3y*width2)/nresol_prv
                  R(3)=acentz+(k2*vers2z*width1 + &
                       k3*vers3z*width2)/nresol_prv
                  call reducn(R,DIRC,BCELL)
                  absden=vMad(R,Q,gEwald,EPSx,iCheck)
                  xcoord=k2*width1/nresol_prv
                  ycoord=k3*width2/nresol_prv
                  if(lochop.ne.hichop) then
                     if(absden.gt.hichop) then
                        absden=hichop
                     else if(absden.lt.lochop) then
                        absden=lochop
                     end if
                  end if
                  write(32,*) xcoord,ycoord,absden*multcon
               END DO
               write(32,*)
            END DO
            close(32)
!______ plot the density: previewing
            lenght3=8
            call Plot3d('test.dat',lenght3,Title, &
                 'X-coordinate (A)    ','Y-coordinate (A)    ', &
                'Coulomb potential   ',  'Screen', 33, &
                nclasses,nresol_prv,type_prv)
         else
            write(*,*)'ERROR! You still have undefined parameters!'
         end if
!
!__________ real calculation
!
     else if(item.eq.9) then
        if(iQuit.eq.0.and.central_p) then
           open(31,file=filen(:lenght),status='unknown',form='formatted')
           write(*,*)'The file '//filen//' has been opened for the PLOT.'
           write(*,*)'Working on the real plot: writing to '// &
                                       filen(1:lenght)//' ...'
           DO K2=-NRESOL/2,NRESOL/2
              DO K3=-NRESOL/2,NRESOL/2
                 R(1)=acentx+(k2*vers2x*width1+k3*vers3x*width2)/nresol
                 R(2)=acenty+(k2*vers2y*width1+k3*vers3y*width2)/nresol
                 R(3)=acentz+(k2*vers2z*width1+k3*vers3z*width2)/nresol
                 call reducn(R,DIRC,BCELL)
                 absden=vMad(R,Q,gEwald,EPSx,iCheck)
                 xcoord=k2*width1/nresol
                 ycoord=k3*width2/nresol
                 if(lochop.ne.hichop) then
                    if(absden.gt.hichop) then
                       absden=hichop
                    else if(absden.lt.lochop) then
                       absden=lochop
                    end if
                 end if
                 write(31,*) xcoord,ycoord,absden*multcon
              END DO
           END DO
           close(31)
           write(*,*) '.... File '//filen(1:lenght)//' has been created! ....'
           Yes_Do=.true.
        else
           write(*,*)'ERROR! You still have undefined parameters!'
        end if
!
!__________ quit option
!
     else if(item.eq.10) then
        return
     else
        go to 100
     end if
     go to 1
100  write(*,*)'ERROR! Try again!'
     go to 1
end subroutine planeM

real*8 function vMad(Pnt,Q,gEwald,EPSx,iCheck)
!.....................................................................
! vMad - Madelung potential at the point Pnt from the whole lattice
!        of charges Q(ion).
! iCheck=0 - atomic radii are ignored, i.e. atoms are meant to be
!            point charges; the potential near any atom is therefore
!            very large; however, exactly at the atomic site it is
!            defined properly;
! iCheck=1 - atomic radii play their role: if Pnt happens to be inside
!            any atomic sphere (spheres are not allowed to overlap!),
!            then the atomic charge is assumed to be uniformly spread
!            over the sphere, so that the potential inside the sphere
!            produced by this atom is well defined and is calculated
!            properly (i.e. from the inside and the outside parts).
!            With this definition, the potential exactly at the lattice
!            site is NOT the same as the convential Ewald method gives!
!.....................................................................
use param
use atoms
implicit none
real*8 :: Q(NIONS),Pnt(3),a(3),dzero=0.0,EPSx,gEwald,Ew
integer ion,iCheck
      do ion=1,NIONS
         a(1)=TI(1,ion)-Pnt(1)
         a(2)=TI(2,ion)-Pnt(2)
         a(3)=TI(3,ion)-Pnt(3)
         call Madelung(a,gEwald,EPSx,Ew,iCheck,ion)
         vMad=vMad+Q(ion)*Ew
      end do
!________ convert to eV (this is the case, however, if factor=1)
!      vMad=vMad*51.42322361
!________ convert to Volts (this is the case, however, if factor=1)
!         (as in CETEP, see ewaltr.f)
      vMad=vMad*14.39976868
end function vMad

subroutine Madelung(X,gEwald,EPSx,Ew,iCheck,ion)
!.......................................................................
! The Coulomb potential (Ew) at the point X from point-ion lattice.
! The Ewald's method is used here.
! Besides, the Evien's idea of organizing the summations over the
! direct and the reciprocal lattices ("by shells") is implemented.
!.......................................................................
! EPSx - the precision of the lattice Ewald's summation for x2:
!       EPSx=-ln(EPSew)
!.......................................................................
use param
use atoms
implicit none
real*8 X(3),Xc(3),Y(3),gEwald,EPSx,Ew,gE2,Em0D,Em0I,Em
integer iCheck,ion,N,iDir,iInv,N1,N2,N3
logical FlagD,FlagI,Singul
real*8 :: pi=3.141592654,tiny=0.00001,dzero=0.0,urfc9
real*8 y2,y1,x2,dist,x1,g2,qu,Xg,cosXg

      gE2=gEwald*gEwald
      Singul=.false.
!.......... put summands and Em to zero:
      Em0D=dzero
      Em0I=dzero
      Em=dzero
!
!....... The both lattices are built by shells numbered using the.......
!   index N=0,1,2,... where N=0 belongs to the 0 site; inside every
!   shell the lattice vectors are computed as N1*a1+N2*a2+N3*a3,
!   where a1,a2,a3 are basic translations (AI for the direct and BI for
!   the invers lattices, respectively), and N1,N2,N3 - indices for the
!   shell, at least one of them is +N or -N.
!.......................................................................
!
!....... The contributions from all shells N. Construction of the shells
!        N by means of the Evien's method.
      FlagD=.True.
      FlagI=.True.
!....... FlagD and FlagI are logical variables for interrupting of the
!   summations over the direct and inverse lattices, respectively. At the
!   beginning they are .true. and the summations are allowed. But, if for
!   every term in the shell the corresponding contribution is small
!   enough, then the variable becomes .false. and suppresses the corresp.
!   summation for all sequential shells. If both are .false., the both
!   summations are stopped. The property is obtained by checking the
!   variables iDir and iInv for 0 values at the end of every shell.
!   iDir and iInv are the numbers of sites in the shell which give a
!   nonzero contribution).
!
      N=-1
30    N=N+1
      iDir=0
      iInv=0
      do N3=-N,N
         do N2=-N,N
            do 250 N1=-N,N
               if(N3.ne.N.and.N3.ne.-N) then
                  if(N2.ne.N.and.N2.ne.-N) then
                     if(N1.ne.N.and.N1.ne.-N) go to 250
                  end if
               end if
!
!____________ summation over the direct lattice; if iCheck=1 and the
! distance to the atom y2 < QradI, then the atom 'ion' is first of all
! removed here; then, at the end of the routine, it is added back but
! with a proper contribution.
               if( .not.FlagD ) go to 100
               Y(1)=N1*DIRC(1,1)+N2*DIRC(2,1)+N3*DIRC(3,1)+X(1)
               Y(2)=N1*DIRC(1,2)+N2*DIRC(2,2)+N3*DIRC(3,2)+X(2)
               Y(3)=N1*DIRC(1,3)+N2*DIRC(2,3)+N3*DIRC(3,3)+X(3)
               Xc(1)=gEwald*Y(1)
               Xc(2)=gEwald*Y(2)
               Xc(3)=gEwald*Y(3)
               y2= Y(1)*Y(1) + Y(2)*Y(2) + Y(3)*Y(3)
               y1=sqrt(y2)
               x2= y2*gEwald*gEwald
               if( x2.le.EPSx ) then
                  iDir=iDir + 1
                  if(iCheck.eq.1 .and. y1.lt.QradI(ion)) then
                     dist=y1
                     Singul=.true.
                  end if
                  if(x2.lt.tiny .or. Singul) then
                     Em=-gEwald*1.128379167*(1.0-x2*(0.33333333-0.1*x2))
                  else
                     x1= sqrt(x2)
                     Em0D=Em0D+urfc9(x1)*exp(-x2)
                  end if
               end if
!
!____________summation over the reciprocal lattice
 100           if(N.eq.0) go to 30
               if( .not.FlagI ) go to 250
               Xc(1)=(N1*BCELL(1,1)+ N2*BCELL(2,1)+ N3*BCELL(3,1))*2*PI
               Xc(2)=(N1*BCELL(1,2)+ N2*BCELL(2,2)+ N3*BCELL(3,2))*2*PI
               Xc(3)=(N1*BCELL(1,3)+ N2*BCELL(2,3)+ N3*BCELL(3,3))*2*PI
               g2= Xc(1)*Xc(1) + Xc(2)*Xc(2) + Xc(3)*Xc(3)
               x2= g2/(4.0*gE2)
               if( x2.le.EPSx ) then
                  iInv=iInv + 1
                  qu = exp(-x2)/g2
                  Xg = Xc(1)*X(1) + Xc(2)*X(2) + Xc(3)*X(3)
                  cosXg = cos(Xg)
                  Em0I=Em0I+qu*cosXg
               end if
250         end do
         end do
      end do
      if(iDir.eq.0) FlagD = .False.
      if(iInv.eq.0) FlagI = .False.
      if(FlagD.or.FlagI) go to 30
!
!.......... compile the whole thing from Dir. and INv. contributions
!
      Ew=4*PI/VOLC*Em0I+gEwald*Em0D+Em
!
!.......... this is because the real-space part does not give (as it is)
!           zero while averaging over the unit cell (is needed only for
!           charged cells wrt point charges, i.e. when the non-point part
!           is also added)
!
!      Ew=Ew-PI/(gE2*VOLC)
!
!_________ add back the proper contribution from the 'ion'
!      if(Singul) Ew=Ew-dist*dist/(2*QradI(ion)**3)+1.5/QradI(ion)
      if(Singul) Ew=Ew+dist*dist/QradI(ion)**3
end subroutine Madelung

real*8 function URFC9(Y)
!......................................................................
!                 urfc9(y)=erfc(y)*exp(y**2)/y ,
!  where erfc(y) is the usual error function, erf(y)=1-erfc(y).
!......................................................................
!  We use here the interpolation given in the Abramovitz and Steagan
!  book for the erfc(y).
!......................................................................
implicit none
real*8 T,Y
      T=1./(1.+0.3275911*Y)
      URFC9 = T*(0.2548296+T*(-0.2844967+T*(1.4214137+T*(-1.45315203+ &
      T*1.0614054))))/Y
end function URFC9

subroutine best_Ewald(DIRC,BCELL,gEwald)
!.......................................................................
!       The "best" Ewald constant gEwald is estimated as:
!
!              gEwald**2 = BCELL_min/(2*DIRC_min)
!
!  where BCELL_min and DIRC_min are the minimal inverse and direct lattice
!  vectors, respectively.
!.......................................................................
! DIRC - direct lattice vectors
! BCELL - reciprocal lattice vectors (without 2*pi)
!..................................................................
implicit none
real*8 :: a(3),DIRC(3,3),BCELL(3,3),pi=3.141592654
real*8 DIRC_min,BCELL_min,aa,gEwald
integer n1,n2,n3
      DIRC_min=100000.0
      BCELL_min=100000.0
      do n1=-1,1
         do n2=-1,1
            do 10 n3=-1,1
               if(n1.eq.0 .and. n2.eq.0 .and. n3.eq.0) go to 10
!_______ find DIRC_min
               a(1)=n1*DIRC(1,1)+n2*DIRC(2,1)+n3*DIRC(3,1)
               a(2)=n1*DIRC(1,2)+n2*DIRC(2,2)+n3*DIRC(3,2)
               a(3)=n1*DIRC(1,3)+n2*DIRC(2,3)+n3*DIRC(3,3)
               aa=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
               if(aa.le.DIRC_min) DIRC_min=aa
!_______ find BCELL_min
               a(1)=n1*BCELL(1,1)+n2*BCELL(2,1)+n3*BCELL(3,1)
               a(2)=n1*BCELL(1,2)+n2*BCELL(2,2)+n3*BCELL(3,2)
               a(3)=n1*BCELL(1,3)+n2*BCELL(2,3)+n3*BCELL(3,3)
               aa=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))*2*PI
               if(aa.le.BCELL_min) BCELL_min=aa
10          end do
         end do
      end do
      gEwald=sqrt( 0.5*BCELL_min/DIRC_min )
      write(*,*)'... The "best" gEwald = ',gEwald,'...'
end subroutine best_Ewald
