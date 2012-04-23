subroutine simulate(grid,Dip,Quadr,filen,lenght)
!.......................................................................
!   We simulate the charge density in the cell by a collection of point
! charges which are positioned in a predefined parallelepiped box in the
! position of a grid inside the box. The charges are chosen to match
! the exact potential grid(i,j,k) in the area in the cell which is
! OUTSIDE the molecule (i.e.outside the box). Then, the positions and
! the values of charges are gievn as the output to the file filen.
!
!_____Input:
!  grid(i,j,k) - place in the memory for the exact potential (including
!                the nuclei) within the cell
!  Dip(xyz) - dipole moment of the molecule (e*A)
!  Quadr(xyz,xyz) - quadrupole moment of the molecule (e*A^2)
!.......................................................................
use param
use menu
use atoms
implicit none
integer, parameter :: N_Simul0=200
real*8 GRID(NGX,NGY,NGZ),cntr(3),Box(3)
integer,dimension(3) :: ngrid=(/5,5,5/)
real*8 :: tiny=0.0001,EPSew=1.0E-10
real*8 :: Dip(3),Quadr(3,3),Bdir(3,3),Rdir(3,3),x(3),EPSx,D(N_Simul0)
character filen*12,answer,Poten_units*3,cha2*2
logical Yes_Do,Box_Def,Poten,Box_Dir,Err,Yes_Check
integer iQuit,ifreq,i,j,N_Simul,iGrid,ix,iy,iz,k,i1,j1,ixyz,lenght
real*8 fCENTX,fCENTY,fCENTZ,w,VOL,error

      cntr=0.0
      ifreq=5
!......................................................................
! iQuit = 0 - not quit, proceed with calculation
!         1 - quit, do not proceed with calculation.
!......................................................................
      Yes_Do=.false.
      Yes_Check=.false.
      Box_Def=.false.
      Box_Dir=.false.
      Poten=.false.
      Poten_units='e/A'
1     iQuit=0
      write(*,*)'..............MENU for SIMULATE ......................'
      write(*,*)'......... Change these parameters if necessary:.......'

      write(*,'(a,3(f10.5,1x))') &
                          '>>>> Dipole moment (e*A): ',(Dip(j),j=1,3)
      write(*,'(a)') '>>>> Quadrupole moment tensor (e*A*A):'
      do i=1,3
        write(*,'(10x,3(1x,f10.5))') (Quadr(i,j),j=1,3)
      end do

      write(*,*)

      write(*,'(a35,e12.6)')'   1. Precision of the summations: ',EPSew
      EPSx=-log(EPSew)

      write(*,'(a)')'   2. The center of your box:'
      fCENTX=BCELL(1,1)*Cntr(1)+BCELL(1,2)*Cntr(2)+BCELL(1,3)*Cntr(3)
      fCENTY=BCELL(2,1)*Cntr(1)+BCELL(2,2)*Cntr(2)+BCELL(2,3)*Cntr(3)
      fCENTZ=BCELL(3,1)*Cntr(1)+BCELL(3,2)*Cntr(2)+BCELL(3,3)*Cntr(3)
      write(*,'(a,f8.3,2(a1,f8.3),a9,f8.3,2(a1,f8.3),a1)') &
        ' A=> (',Cntr(1),',',Cntr(2),',',Cntr(3), &
        '), fr=> (', fCENTX,',',fCENTY,',',fCENTZ,')'

      if(Box_Dir) then
        write(*,'(a)')'   3. Directions of the box sides are:'
        write(*,'(7x,a,i1,a,f10.5,1x,f10.5,1x,f10.5)') &
                       ('vector ',i,' is ',(Bdir(i,j),j=1,3),i=1,3)
      else
        write(*,'(a)')'   3. Directions of the box sides are: undefined'
        iQuit=1
      end if

      if(Box_Def) then
        write(*,'(a,3(f10.5,1x))') &
         '   4. Lengthes of the sides of the box are:',(Box(i),i=1,3)
      else
        write(*,'(a,3(f10.5,1x))') &
         '   4. Lengthes of the sides of the box are: undefined'
        iQuit=1
      end if

      write(*,'(a,3i5)')'   5. Grid for the simulating charge: ', &
                                                 (ngrid(i),i=1,3)
      N_Simul=(ngrid(1)-2)*(ngrid(2)-2)*(ngrid(3)-2)
      write(*,'(a,i5)')'      Total number of simulating charges = ',N_Simul

      if(Poten) then
        write(*,'(a)') '   6. Read the exact potential on the grid <=== DONE!'
        write(*,'(a)') '   7. Units of the potential: '//Poten_units
      else
        write(*,'(a)') '   6. Read the exact potential on the grid <=== undefined!'
        iQuit=1
      end if

      write(*,'(a,i5,a)') &
       '   8. We fit to every ',ifreq,'-th grid point in the potential'
      if(Yes_Check) write(*,'(6x,a,i7,a,i7)') 'There are ',iGrid, &
         ' outside grid points in the cell to use out of ',NPLWV

      if(Yes_Do) then
        write(*,'(a)') &
         '   9. Calculate simulating charges and put them in the file ' &
                                            //filen//' <= DONE!'
        write(*,'(a,e12.6)') '      After fitting the average error is ',error
      else
        write(*,'(a)') &
         '   9. Calculate simulating charges and put them in the file '//filen
      end if

      write(*,'(a)')'-------  G e n e r a l  s e t t i n g s ---------'
      write(*,'(a)')'  An. Coordinates are specified in: '//angstr
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
![1]__________ specify the precision of the summations in the Ewald method
!
      ELSE IF(trim(cha2).eq.'1') THEN
 40      write(*,*)'Give the precision for the Ewald summation:'
         read(*,*,err=40) EPSew
         if(EPSew .lt. 1.0e-13) go to 40
!
![2]__________ give the box position in the cell
!
      ELSE IF(trim(cha2).eq.'2') THEN
         WRITE(*,*)'Give position of your box'
         call givepoint(Cntr(1),Cntr(2),Cntr(3),angstr)
         Yes_Do=.false.
!
![3]__________ specify directions of the box sides (3 vectors)
!
      ELSE IF(trim(cha2).eq.'3') THEN
 12      write(*,*)'Specify 3 directions of the box sides:'
         write(*,*) 'Give manually (m,M) or along the basic vectors (other)?'
         read(*,'(a)',err=12) answer
         if(answer.eq.'m' .or. answer.eq.'M') then
13          write(*,*)' Specify the 1st direction:'
            read(*,*,err=13) (Bdir(1,j),j=1,3)
            w=sqrt(Bdir(1,1)**2+Bdir(1,2)**2+Bdir(1,3)**2)
            if(w.lt.tiny) go to 13
14          write(*,*)' Specify the 2nd direction:'
            read(*,*,err=14) (Bdir(2,j),j=1,3)
            w=sqrt(Bdir(2,1)**2+Bdir(2,2)**2+Bdir(2,3)**2)
            if(w.lt.tiny) go to 14
15          write(*,*)' Specify the 3rd direction:'
            read(*,*,err=15) (Bdir(3,j),j=1,3)
            w=sqrt(Bdir(3,1)**2+Bdir(3,2)**2+Bdir(3,3)**2)
            if(w.lt.tiny) go to 15
         else
            Bdir=DIRC
         end if
         do i=1,3
            w=sqrt(Bdir(i,1)**2+Bdir(i,2)**2+Bdir(i,3)**2)
            do j=1,3
               Bdir(i,j)=Bdir(i,j)/w
            end do
         end do
         call BASTR(Bdir,Rdir,VOL,0)
         if(VOL.le.tiny) then
            write(*,*)'ERROR! Vectors are collinear! Try again!'
            go to 12
         end if
         Yes_Do=.false.
         Box_Dir=.true.
!
![4]__________ lengthes of the box sides
!
      ELSE IF(trim(cha2).eq.'4') THEN
 10      WRITE(*,*)'Enter the length of every side of the box (in Ang):'
         READ(*,*,err=10) (Box(i),i=1,3)
         do i=1,3
            if(Box(i).le.0.0) go to 10
         end do
         Box_Def=.true.
         Yes_Do=.false.
!
![5]__________ number of grid points along every direction
!
      ELSE IF(trim(cha2).eq.'5') THEN
997      write(*,*) 'Give numbers of grids along every box side (>=3):'
         read(*,*,err=997) (ngrid(i),i=1,3)
         do i=1,3
            if(ngrid(i).le.2) go to 997
         end do
         N_Simul=(ngrid(1)-2)*(ngrid(2)-2)*(ngrid(3)-2)
         if(N_Simul.gt.N_Simul0) then
           write(*,'(a,i5)')'ERROR! Too fine grid! N_Simul > ',N_Simul0
           go to 997
         else if(N_Simul.lt.9) then
           write(*,'(a,i5)')'ERROR! The grid is too scarse! Not enough!'
           go to 997
         end if
         Yes_Do=.false.
!
![6]__________ read the exact potential from poten.200 file
!
      ELSE IF(trim(cha2).eq.'6') THEN
         write(*,*)'Reading in the exact potential from poten.200 ...'
         open(33,file='poten.200',status='old',form='formatted',err=150)
         do iXYZ=1,NPLWV
            read(33,*,err=170,end=170) iX,iY,iZ,grid(iX,iY,iZ)
         end do
         close (33)
         write(*,*)'Done!'
         Yes_Do=.false.
         Poten=.true.
!
![7]__________ choose the units in which potential (poten.200) is given
!
      ELSE IF(trim(cha2).eq.'7') THEN
         if(Poten_units.eq.' eV') then
            Poten_units=' au'
         else if(Poten_units.eq.' au') then
            Poten_units='e/A'
         else if(Poten_units.eq.'e/A') then
            Poten_units=' eV'
         end if
!
![8]__________ specify the frequency at which we scan the potential on the grid
!           and also calculate the total number of grid points outside the cell
!           which will be used with this frequency and the box specified
!
      ELSE IF(trim(cha2).eq.'8') THEN
        if(.not.(Box_Dir.and.Box_Def)) then
           write(*,*)'ERROR! Specify the box completely!'
           go to 1
        end if
 17     write(*,*) 'Specify the frequency at which to scan the '// &
                                           'potential grid (1,2,...)'
        read(*,*,err=17) ifreq
        if(ifreq.lt.1 .or. ifreq.ge.min(NGX,NGY,NGZ)) go to 17
        call check_grid(iGrid,ifreq,DIRC,Cntr,Box,Bdir)
        Yes_Check=.true.
!
!
![9]__________ calculation: fitting the potential outside the box by a defined
!                        above grid of point charges
!
      ELSE IF(trim(cha2).eq.'9') THEN
        if(iQuit.ne.0) then
           write(*,*)'ERROR! You still have undefined parameters!'
           go to 1
        end if
        open(31,file=filen(:lenght),status='unknown', &
                                             form='formatted',err=160)
        write(*,*)'The file '//filen//' has been opened for the simulation.'
        write(*,*)'Working on the simulation ...'
        call do_simulate(31,grid,Box,ngrid,Cntr, &
                           N_Simul,Dip,Quadr,EPSx,Bdir,Poten_units,D, &
                           ifreq,error,Err)
        if(Err) then
           write(*,*)'FATAL! The algorithm failed!'
           close (31,status='delete')
           write(*,*)'.... File '//filen(1:lenght)//' has been deleted! ....'
        else
           close (31)
           write(*,*)'.... File '//filen(1:lenght)//' has been created! ....'
           Yes_Do=.true.
        end if
!
![Co].... display atomic positions
!
      ELSE IF(trim(cha2).eq.'Co') THEN
         k=0
         do i=1,NSPEC
            write(*,'(3a,i5)')'.... Species <',Species(i),'> ......> ',i
            write(*,'(a)')' Tot   # |----------  Fractional ---------|' &
                 //'-------------  Cartesian ---------|'
            do j=1,NspN(i)
               k=k+1
               do i1=1,3
                  x(i1)=0.0
                  do j1=1,3
                     x(i1)=x(i1)+TI(j1,k)*BCELL(i1,j1)
                  end do
               end do
               write(*,'(2(x,i3),x,3(x,f10.5),2x,3(x,f10.5),3x,a,4x,a)') &
                   k,j,(x(i1),i1=1,3),(TI(j1,k),j1=1,3)
            end do
         end do
         write(*,*)'Hit ENTER when done ...'
         read(*,*)
!
![Q]__________ return to the previous menu
!
      ELSE IF(trim(cha2).eq.'Q') THEN
         return
      ELSE
         write(*,*)'ERROR! Try again!'
      END IF
      go to 1
!
!....... error
 170  write(*,*)'FATAL! Error while reading poten.200 file!'
      go to 1
 150  write(*,*)'FATAL! Error while opening poten.200 file!'
      go to 1
 160  write(*,*)'FATAL! Error while opening '//filen(1:lenght)//' file'
      go to 1
end subroutine simulate

subroutine do_simulate(nunit,grid,Box,ngrid,Cntr, &
                           N_Simul,Dip,Quadr,EPSx,Bdir,Poten_units,D, &
                           ifreq,error,Err)
!.......................................................................
!   We simulate the charge density in the cell by a collection of point
! charges within a box on a grid by minimising (analytically) the error
! between the exact potential in the OUTSIDE of the box and the potential
! there due to the simulating charges.
!
!_____Input:
!  grid(i,j,k) - the exact potential (including the nuclei) within the cell
!  Dip(xyz) - dipole moment of the molecule (e*A)
!  Quadr(xyz,xyz) - quadrupole moment of the molecule (e*A^2)
!_____Box information:
!  Bdir(#,xyz) - directions of the box sides
!  Box(#) - lengthes of every side of the box
!  ngrid(#) - number of grid points along every side of the box (in fact,
!             boundary points will not be used)
!  Cntr(xyz) - position of the box center
!.......................................................................

use param
use atoms
implicit none
integer, parameter :: N_Simul0=200,N_Simul9=N_simul0-9
real*8 GRID(NGX,NGY,NGZ),R(3),x(3)
real*8 Cntr(3),Box(3),Charges(N_Simul0),R_Ch(3,N_Simul0),VV(N_Simul0)
integer ngrid(3),ng(3),INDX(N_Simul0),i,j,ic,n1,n2,n3,ic1,id,icc,icc1
integer ijk,j0,ix,iy,iz,num,N_simul,ifreq,idc,nunit
real*8 Dip(3),Quadr(3,3),e(3,3),corner(3),Bdir(3,3),BoxD(3,3)
real*8 g(9),C(9,9),Cm1(9,9),Gtlda(9),Ftlda(9,N_Simul9),Fi_Tlda
real*8 A(N_Simul9,N_Simul9),B(N_Simul9),D(N_Simul0),EPSx,Ew
real*8 Dtlda(N_Simul9),BoxR(3,3),unfact,dx,VOL,s,gEwald,error
character Poten_units*3
logical out,Err
      Err=.false.
!
!______ working out e/A units for the potential in grid(iX,iY,iZ) from the
!       units Poten_units it is in
!
      if(Poten_units.eq.'e/A') unfact=1.0
      if(Poten_units.eq.' au') unfact=1.0/0.529177249
      if(Poten_units.eq.' eV') unfact=1.0/(0.529177249*27.212)
!
!
!______ box sides are given by vectors BoxD(#,xyz) which are constructed out
!       of Bdir(#,xyz).
!       The vectors e(#,xyz) give basic vectors of the box grid.
      do i=1,3
         ng(i)=ngrid(i)-1
         dx=Box(i)/ng(i)
         do j=1,3
            BoxD(i,j)=Bdir(i,j)*Box(i)
            e(i,j)=Bdir(i,j)*dx
         end do
      end do
!______ position of the 0-th corner (n1=n2=n3=0) of the box
      do i=1,3
         corner(i)=Cntr(i)-0.5*( BoxD(1,i)+BoxD(2,i)+BoxD(3,i) )
      end do
!______ reciprocal vectors BoxR(#,xyz) (without 2*pi) associated with
!       the box grid
      call BASTR(BoxD,BoxR,VOL,0)
!
!......... positions of N_Simul simulating charges; the box surface
!          is ignored so that the charges are only inside the box
!
      write(*,*)'Calculating the positions of charges ...'
      ic=0
      do n1=1,ng(1)-1
         do n2=1,ng(2)-1
            do n3=1,ng(3)-1
               ic=ic+1
               R_Ch(1,ic)=n1*e(1,1)+n2*e(2,1)+n3*e(3,1)+corner(1)
               R_Ch(2,ic)=n1*e(1,2)+n2*e(2,2)+n3*e(3,2)+corner(2)
               R_Ch(3,ic)=n1*e(1,3)+n2*e(2,3)+n3*e(3,3)+corner(3)
            end do
         end do
      end do
!
!......... preparing 9-dim vectors and matrices of conditions
!
      write(*,*)'Preparing 9D vectors and arrays ...'
      g(1)=0.0
      g(2)=Dip(1)
      g(3)=Dip(2)
      g(4)=Dip(3)
      g(5)=Quadr(1,1)
      g(6)=Quadr(2,2)
      g(7)=Quadr(1,2)
      g(8)=Quadr(1,3)
      g(9)=Quadr(2,3)
      do ic=1,9
         call do_F_cond(R_Ch(1,ic),C(1,ic))
      end do
      call invers(C,Cm1,INDX,VV,9,9,Err)
      if(Err) return
!______________ G-tilda
      do id=1,9
         s=0.0
         do i=1,9
            s=s+Cm1(id,i)*g(i)
         end do
         Gtlda(id)=s
      end do
!______________ F-tilda
      do ic1=10,N_Simul
         ic=ic1-9
         call do_F_cond(R_Ch(1,ic1),g)
         do id=1,9
            s=0.0
            do i=1,9
               s=s+Cm1(id,i)*g(i)
            end do
            Ftlda(id,ic)=s
        end do
     end do
!
!.......... calculating the Ewald summations: running the cell grid points
!           which are outside the box. We eventually calculate an array
!           A and a vector B for the equation A*Q=B, where Q are charges
!           10,11,... (designated by ic). After that, we will be able to
!           calculate other charges (designated by id=1,...,9)
!
     call best_Ewald(DIRC,BCELL,gEwald)
!
!________ set A and B to zero before the summation over the grid points
!
     do ic1=10,N_Simul
        ic=ic1-9
        B(ic)=0.0
        do icc1=ic1,N_Simul
           icc=icc1-9
           A(ic,icc)=0.0
           A(icc,ic)=A(ic,icc)
        end do
     end do
!
!________ in the loop over original grid points. We keep the summations
!         - the matrix D - in the file d.tmp for future reference, i.e.
!         calculation of the approximation error
!
      write(*,'(a,i3,a)')'Scanning the grid points with the frequency ',ifreq,' ...'
      QradI(1)=1.0
      open(32,file='d.tmp',status='unknown',form='unformatted')
      ijk=0
      j0=NPLWV/(ifreq**3*10)
      do iZ=0,NGZ-1,ifreq
         do iY=0,NGY-1,ifreq
            do 50 iX=0,NGX-1,ifreq
               !
!_______ get the point (iX,iY,iZ) in Cartesian coordinates
               x(1) = iX*DIRC(1,1)/NGX + iY*DIRC(2,1)/NGY + iZ*DIRC(3,1)/NGZ
               x(2) = iX*DIRC(1,2)/NGX + iY*DIRC(2,2)/NGY + iZ*DIRC(3,2)/NGZ
               x(3) = iX*DIRC(1,3)/NGX + iY*DIRC(2,3)/NGY + iZ*DIRC(3,3)/NGZ
!_____________ ask whether the point (iX,iY,iZ) is outside the box
               ijk=ijk+1
               if(ijk/j0*j0.eq.ijk) write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
               call ask_outside(DIRC,x,corner,BoxR,out)
               if(.not.out) go to 50
!
!_____________ calculate Ewald at the point x from every charge; note that
!              the point X is outside the box so that there is no singularity
               do idc=1,N_Simul
                  R(1)=x(1)-R_Ch(1,idc)
                  R(2)=x(2)-R_Ch(2,idc)
                  R(3)=x(3)-R_Ch(3,idc)
                  call Madelung(R,gEwald,EPSx,Ew,0,1)
                  D(idc)=Ew
               end do
               write(32)  D
!           write(32,'(8(e12.6,1x))') (D(ic),ic=1,N_Simul)
!
!_________________ get Potential-tilda
               s=0.0
               do id=1,9
                  s=s+D(id)*Gtlda(id)
               end do
               Fi_Tlda=grid(iX+1,iY+1,iZ+1)*unfact - s
!
!_________________ get all D(ic)-tilda in Dtlda
               do ic1=10,N_Simul
                  ic=ic1-9
                  s=0.0
                  do id=1,9
                     s=s+D(id)*Ftlda(id,ic)
                  end do
                  Dtlda(ic)=D(ic1)-s
               end do
!
!____________ sum contributions for this grid point to get A and B
!
               do ic1=10,N_Simul
                  ic=ic1-9
                  B(ic)=B(ic)+Fi_Tlda*Dtlda(ic)
                  do icc1=ic1,N_Simul
                     icc=icc1-9
                     A(ic,icc)=A(ic,icc)+Dtlda(ic)*Dtlda(icc)
                     A(icc,ic)=A(ic,icc)
                  end do
               end do
!
50          end do
         end do
      end do
      close (32)
      write(*,*)'Done A and B !'
!
!.......... solve for charges 10, 11, ... by solving A*Q=B. The charges
!           will be in B
!
      write(*,*)'Solving for charges 10,11,...'
      call linear(A,B,INDX,VV,N_Simul-9,N_Simul9,Err)
      if(Err) return
!
!.......... calculate the first 9 charges
!
      write(*,*)'Calculating first 9 charges  ..'
      do id=1,9
         s=0.0
         do ic1=10,N_Simul
            ic=ic1-9
            s=s+Ftlda(id,ic)*B(ic)
         end do
         Charges(id)=Gtlda(id)-s
         write(nunit,'(i5,3(1x,f10.5),5x,f10.5)') &
              id, (R_Ch(j,id),j=1,3),Charges(id)
      end do
      do ic1=10,N_Simul
         ic=ic1-9
         Charges(ic1)=B(ic)
         write(nunit,'(i5,3(1x,f10.5),5x,f10.5)') &
                          ic1, (R_Ch(j,ic1),j=1,3),Charges(ic1)
      end do
      write(*,*)'Done!'
!
!.................... check: calculation of the error in simulating the
!                            potential in the outside area
!
      write(*,*)'Doing the error ...'
      error=0.0
      num=0
      open(32,file='d.tmp',status='old',form='unformatted')
      do iZ=0,NGZ-1,ifreq
         do iY=0,NGY-1,ifreq
            do 60 iX=0,NGX-1,ifreq
!
!_______ get the point (iX,iY,iZ) in Cartesian coordinates
               x(1) = iX*DIRC(1,1)/NGX + iY*DIRC(2,1)/NGY + iZ*DIRC(3,1)/NGZ
               x(2) = iX*DIRC(1,2)/NGX + iY*DIRC(2,2)/NGY + iZ*DIRC(3,2)/NGZ
               x(3) = iX*DIRC(1,3)/NGX + iY*DIRC(2,3)/NGY + iZ*DIRC(3,3)/NGZ
!_____________ ask whether the point (iX,iY,iZ) is outside the box
               call ask_outside(DIRC,x,corner,BoxR,out)
               if(.not.out) go to 60
               num=num+1
!
!_____________ read the D-matrix for the given grid point and all charges
!           read(32,'(8(e12.6,1x))') (D(ic),ic=1,N_Simul)
               read(32) D
!
!_____________ calculate the potential from all charges and the error
               s=0.0
               do ic1=1,N_Simul
                  s=s+Charges(ic1)*D(ic1)
               end do
               error=error+( s-grid(iX+1,iY+1,iZ+1)*unfact )**2
!
60          end do
         end do
      end do
      error=sqrt(error)/num
      close (32)
      write(*,*)'Done!'
end subroutine do_simulate

subroutine ask_outside(DIRC,r,cor,BoxR,out)
!...................................................................
!   Checks if the vector r(xyz) is inside (out=.false.) or outside
! (out=.true.) the box which is defined by:
!  (1) cor(xyz) - position of the box 1st corner
!  (2) vectors BoxR(#,xyz) - "reciprocal" vectors associated with
!                            box sides
!...................................................................
!  Note that all 27 images of the box are taken into account.
!...................................................................
implicit none
real*8 cor(3),r(3),BoxR(3,3),x(3),DIRC(3,3),a
logical out
integer n1,n2,n3,ijk,j,i
      out=.false.
      do n1=-1,1
         do n2=-1,1
            do n3=-1,1
               x(1)=r(1)-cor(1)-n1*DIRC(1,1)-n2*DIRC(2,1)-n3*DIRC(3,1)
               x(2)=r(2)-cor(2)-n1*DIRC(1,2)-n2*DIRC(2,2)-n3*DIRC(3,2)
               x(3)=r(3)-cor(3)-n1*DIRC(1,3)-n2*DIRC(2,3)-n3*DIRC(3,3)
               ijk=0
               do i=1,3
                  a=x(1)*BoxR(i,1)+x(2)*BoxR(i,2)+x(3)*BoxR(i,3)
                  if(a.ge.0.0 .and. a.le.1) ijk=ijk+1
               end do
               if(ijk.eq.3) return
            end do
         end do
      end do
      out=.true.
end subroutine ask_outside

subroutine check_grid(iGrid,ifreq,DIRC,Cntr,Box,Bdir)
!........................................................................
!   We run the scan of the grid outside the box to counter the # of grid
! points outside the box.
!........................................................................
use param
implicit none
real*8 DIRC(3,3),x(3),Cntr(3),Box(3)
integer ngrid(3),ng(3),ifreq,i,j,ijk,ix,iy,iz,j0,iGrid
real*8 corner(3),Bdir(3,3),BoxD(3,3),BoxR(3,3),VOL
logical out
!
!______ box sides are given by vectors BoxD(#,xyz) which are constructed out
!       of Bdir(#,xyz).
!
      do i=1,3
         do j=1,3
            BoxD(i,j)=Bdir(i,j)*Box(i)
         end do
      end do
!______ position of the 0-th corner (n1=n2=n3=0) of the box
      do i=1,3
         corner(i)=Cntr(i)-0.5*( BoxD(1,i)+BoxD(2,i)+BoxD(3,i) )
      end do
!______ reciprocal vectors BoxR(#,xyz) (without 2*pi) associated with
!       the box grid
      call BASTR(BoxD,BoxR,VOL,0)
!
!........ In the loop over original grid points, we scan only the points
!         outside the box.
!
      ijk=0
      j0=NPLWV/(ifreq**3*10)
      do iZ=0,NGZ-1,ifreq
         do iY=0,NGY-1,ifreq
            do 50 iX=0,NGX-1,ifreq
!
!_______ get the point (iX,iY,iZ) in Cartesian coordinates
               x(1) = iX*DIRC(1,1)/NGX + iY*DIRC(2,1)/NGY + iZ*DIRC(3,1)/NGZ
               x(2) = iX*DIRC(1,2)/NGX + iY*DIRC(2,2)/NGY + iZ*DIRC(3,2)/NGZ
               x(3) = iX*DIRC(1,3)/NGX + iY*DIRC(2,3)/NGY + iZ*DIRC(3,3)/NGZ
!_____________ ask whether the point (iX,iY,iZ) is outside the box
               ijk=ijk+1
               if(ijk/j0*j0.eq.ijk) &
                    write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
               call ask_outside(DIRC,x,corner,BoxR,out)
               if(.not.out) go to 50
50          end do
         end do
      end do
      iGrid=ijk
end subroutine check_grid

subroutine do_F_cond(R,F)
!..............................................................
!    The 9D vector of conditions is calculated for the given
! charge position R.
!..............................................................
implicit none
real*8 R(3),F(9),ar2
      ar2=R(1)*R(1)+R(2)*R(2)+R(3)*R(3)
      F(1)=1.0
      F(2)=R(1)
      F(3)=R(2)
      F(4)=R(3)
      F(5)=3*R(1)*R(1)-ar2
      F(6)=3*R(2)*R(2)-ar2
      F(7)=3*R(1)*R(2)
      F(8)=3*R(1)*R(3)
      F(9)=3*R(2)*R(3)
end subroutine do_F_cond
    
