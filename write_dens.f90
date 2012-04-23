subroutine write_dens(grid)
!..................................................................
!    Write density from GRID into a file with the extension .new
!..................................................................
use param
use atoms
use code
implicit none
real*8 GRID(NGX,NGY,NGZ),factor
real*8,parameter :: A_au=1.889725989,au_A=1.0/A_au
character filen*12
integer i,ix,iy,iz,j,nspins

!_______________ VASP input

      if(Which_Code.eq.'  VASP') then
         write(*,*)'Writing the density to CHGCAR.new ...'
         filen='  CHGCAR.new'
         open(11,file='CHGCAR.new',status='unknown',form='formatted')
         call write_vasp_dens(grid)

!_______________ SIESTA input
      else if(Which_Code.eq.'SIESTA') then
         write(*,*) 'Writing density to job.RHO.new ...'
         filen=' job.RHO.new'
         open(11,file='job.RHO.new',status='unknown',form='formatted')
         do i=1,3
            write(11,*) (DIRC(i,j)/au_A,j=1,3)
         end do
         nspins=1
         write(11,*) NGX,NGY,NGZ,nspins
         factor=VOLC*A_au**3
         do iZ=1,NGZ
            do iY=1,NGY
               do iX=1,NGX
                  write(11,*) grid(iX,iY,iZ)/factor
               end do
            end do
         end do
      end if
      close (11)
      write(*,*) '.....> Charge density reaas been created!'
      write(*,*)'Done!'
end subroutine write_dens

subroutine write_vasp_dens(grid)
!..................................................................................
! the charge density GRID is written down into the file CHGCAR.new ( in VASP form.)
!..................................................................................
use param
use atoms
implicit none
integer, parameter :: NCol0=10
integer ix5(NCol0),iy5(NCol0),iz5(NCol0),i,j,Ncol,ijk,j0,ix,iy,iz,ijk5,last
integer nat,jj
real*8 GRID(NGX,NGY,NGZ),r(3)
!
!.......... info about the system
!
      write(11,'(a)') 'created density by lev00 (see write_dens.f)'
      write(11,*) ' 1.0 '
      write(11,'(x,3(2x,f10.6))') DIRc(1,1), DIRc(1,2), DIRc(1,3)
      write(11,'(x,3(2x,f10.6))') DIRc(2,1), DIRc(2,2), DIRc(2,3)
      write(11,'(x,3(2x,f10.6))') DIRc(3,1), DIRc(3,2), DIRc(3,3)
      write(11,*) (NspN(i),i=1,NSPEC)

!_____________ writing coordinates (in fractional)

      write(11,'(a)') 'Direct'
      nat=0
      do i=1,NSPEC
         do j=1,NspN(i)
            nat=nat+1
            r(1)=BCELL(1,1)*TI(1,nat)+BCELL(1,2)*TI(2,nat)+BCELL(1,3)*TI(3,nat)
            r(2)=BCELL(2,1)*TI(1,nat)+BCELL(2,2)*TI(2,nat)+BCELL(2,3)*TI(3,nat)
            r(3)=BCELL(3,1)*TI(1,nat)+BCELL(3,2)*TI(2,nat)+BCELL(3,3)*TI(3,nat)
            write(11,'(3(x,f9.6))') r(1),r(2),r(3)
         end do
      end do

!_____________ writing grid
      write(11,'(/(3i5,x))') NGX,NGY,NGZ
!
!..........writing the charge density in 5 columns
!
      NCol=5
      ijk5=0
      ijk=0
      j0=NPLWV/10

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
                  write(11,'(5(e18.11,x))') (grid(ix5(jj),iy5(jj),iz5(jj)),jj=1,last)
                  ijk5=0
               end if
               if(ijk/j0*j0.eq.ijk) &
                    write(*,'(a,i3,a)') '... done ',ijk/j0*10,' %'
            end do
         end do
      end do
end subroutine write_vasp_dens



