subroutine Plot(ng,file,lenght,title,title_pl, &
                           xaxis,yaxis,where,l3,iFlag, &
                               Smeared,SmearedPrev,Up_lim,Up_lim_E)
!......................................................................
!   Plot curves eighter on the screen or on to a ps-file
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
!......................................................................
! ng - number of data files to be used simultaneously
! file(i) - name of the data file for the GNUPLOT, i=1,...,ng;
! lenght(i) - lenght of the data files names;
! title_pl(i) - subtitle for each curve
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!    When Smeared='S' and SmearedPrev=.true., we show column 3 for
! the DOS and column 5 for the LDOS (i.e. smeared densities); otherwise,
! i.e. if SmearedPrev=.false., we show columns 2 and 4, respectively.
!    If Smeared='N', we show columns 2 and 3 respectively as there
! are no smeared results at all.
!    If Smeared='F' and SmearedPrev=.true., we show column 2 for the
! smeared DOS and column 3 for the smeared LDOS in file//'_' ; otherwise,
! i.e. if SmearedPrev=.false., we show columns 2 and 3 in file.
!.....................................................................
!
implicit none
character title*50,xaxis*20,yaxis*20,where*6,Smeared
character file(ng)*12,f_gr*12,title_pl(ng)*7,cha1,cha2,cha*2
integer lenght(ng),ng,l3,iFlag,i
logical SmearedPrev
real*8 Up_lim,Up_lim_E

!.......... name of the file for the group (is used if ng.ne.1)
      f_gr='dos.dat0'//file(1)(9:lenght(1))
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis
300   format ('set data style lines'/ &
              'set format y "%.4f"'/ &
              'set key'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')

!.......... axes limits
      if(Up_lim.ne.0.0) then
         write(l3,306) Up_lim
306      format('set yrange [0:',f10.5,']')
      end if
      if(Up_lim_E.ne.0.0) then
         write(l3,307) Up_lim_E
307      format('set xrange [:',f10.5,']')
      end if
!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
         write (l3,401)
401      format ('#set terminal post color "Times-Roman" 14')
         if(ng.eq.1) then
            write(l3,'(a)')'#set output "'//file(1)(1:lenght(1))//'.ps"'
         else
            write(l3,'(a)')'#set output "'//f_gr(1:lenght(1))//'.ps"'
         end if
      else
         write (l3,301)
301      format ('set terminal post color "Times-Roman" 14')
         if(ng.eq.1) then
            write (l3,'(a)')'set output "'//file(1)(1:lenght(1))//'.ps"'
         else
            write (l3,'(a)')'set output "'//f_gr(1:lenght(1))//'.ps"'
         end if
      end if
!.......... common part: .........................................
!____ we show columns cha1 and cha2 for DOS and LDOS, respectively,
!     depending on Smeared='N','F','S':
      cha='" '
      if(Smeared.eq.'F'.and.SmearedPrev) cha='_"'
      if(Smeared.eq.'N' .or. Smeared.eq.'F') then
         cha1='2'
         cha2='3'
      else if(Smeared.eq.'S') then
         if(SmearedPrev) then
            cha1='3'
            cha2='5'
         else
            cha1='2'
            cha2='4'
         end if
      end if
!_____ plot commands
      if(iFlag.eq.1) then
         write (l3,*)  'plot "'//file(1)(1:lenght(1))//cha// &
              ' u 1:'//cha1//' t "'//title_pl(1)//'"'//char(92)
         write (l3,*)', "'//file(1)(1:lenght(1))//cha//' u 1:'// &
              cha2//' t "'//title_pl(1)//'(loc)"'//char(92)
      else
         write (l3,*)'plot "'//file(1)(1:lenght(1))//cha//' u 1:'// &
              cha1//' t "'//title_pl(1)//'"'//char(92)
      end if
      if(ng.eq.1) go to 10
      do i=2,ng
         if(iFlag.eq.1) then
            write (l3,*)', "'//file(i)(1:lenght(i))//cha//' u 1:'// &
                 cha1//' t "'//title_pl(i)//'"'//char(92)
            write (l3,*)', "'//file(i)(1:lenght(i))//cha//' u 1:'// &
                 cha2//' t "'//title_pl(i)//'(loc)"'//char(92)
         else
            write (l3,*)', "'//file(i)(1:lenght(i))//cha//' u 1:'// &
                 cha1//' t "'//title_pl(i)//'"'//char(92)
         end if
      end do
10    write(l3,*)' '
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
!.......... Check postscript file if it was specified ................
      if(where.ne.'Screen') then
         if(ng.ne.1) then
            write(6,*)'The postscript file '// &
                 f_gr(1:lenght(1))//'.ps" has been created!'
         else
            write(6,*)'The postscript file '// &
                 file(1)(1:lenght(1))//'.ps" has been created!'
         end if
      end if
end subroutine Plot

subroutine Plot1(file,lenght,title,title_pl, &
                                   xaxis,yaxis,where,l3,iFlag, &
                               Smeared,SmearedPrev,Up_lim,Up_lim_E)
!......................................................................
!   Plot curves eighter on the screen or on to a ps-file
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
!......................................................................
! file - name of the data file for the GNUPLOT
! lenght - lenght of the data files names;
! title_pl - subtitle for each curve
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!    When Smeared='S' and SmearedPrev=.true., we show column 3 for
! the DOS and column 5 for the LDOS (i.e. smeared densities); otherwise,
! i.e. if SmearedPrev=.false., we show columns 2 and 4, respectively.
!    If Smeared='N', we show columns 2 and 3 respectively as there
! are no smeared results at all.
!    If Smeared='F' and SmearedPrev=.true., we show column 2 for the
! smeared DOS and column 3 for the smeared LDOS in file//'_' ; otherwise,
! i.e. if SmearedPrev=.false., we show columns 2 and 3 in file.
!.....................................................................
!
implicit none
character title*50,xaxis*20,yaxis*20,where*6,Smeared
character file*12,f_gr*12,title_pl*7,cha1,cha2,cha*2
integer lenght,iFlag,l3
logical SmearedPrev
real*8 Up_lim,Up_lim_E

!.......... name of the file for the group (is used if ng.ne.1)
      f_gr='dos.dat0'//file(9:lenght)
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis
300   format ('set data style lines'/ &
              'set format y "%.4f"'/ &
              'set key'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')

!.......... exes limits
      if(Up_lim.ne.0.0) then
        write(l3,306) Up_lim
 306    format('set yrange [0:',f10.5,']')
      end if
      if(Up_lim_E.ne.0.0) then
        write(l3,307) Up_lim_E
 307    format('set xrange [:',f10.5,']')
      end if

!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
         write (l3,401)
 401     format ('#set terminal post color "Times-Roman" 14')
         write(l3,'(a)')'#set output "'//file(1:lenght)//'.ps"'
      else
         write (l3,301)
 301     format ('set terminal post color "Times-Roman" 14')
         write (l3,'(a)')'set output "'//file(1:lenght)//'.ps"'
      end if
!.......... common part: .........................................
!____ we show columns cha1 and cha2 for DOS and LDOS, respectively,
!     depending on Smeared='N','F','S':
      cha='" '
      if(Smeared.eq.'F'.and.SmearedPrev) cha='_"'
      if(Smeared.eq.'N' .or. Smeared.eq.'F') then
        cha1='2'
        cha2='3'
      else if(Smeared.eq.'S') then
        if(SmearedPrev) then
          cha1='3'
          cha2='5'
        else
          cha1='2'
          cha2='4'
        end if
      end if
!_____ plot commands
      if(iFlag.eq.1) then
         write (l3,*)  'plot "'//file(1:lenght)//cha// &
                ' u 1:'//cha1//' t "'//title_pl//'"'//char(92)
         write (l3,*)', "'//file(1:lenght)//cha//' u 1:'// &
                cha2//' t "'//title_pl//'(loc)"'//char(92)
      else
         write (l3,*)'plot "'//file(1:lenght)//cha//' u 1:'// &
                cha1//' t "'//title_pl//'"'//char(92)
      end if

      write(l3,*)' '
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
!.......... Check postscript file if it was specified ................
      if(where.ne.'Screen') then
         write(6,*)'The postscript file '// &
                      file(1:lenght)//'.ps" has been created!'
      end if
end subroutine Plot1

subroutine Plot2(ng,file,lenght,title,title_pl, &
                     xaxis,yaxis,where,l3,Up_lim,Up_lim_E,ShowTotDos)
!......................................................................
!   Plot curves eighter on the screen or on to a ps-file
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
!......................................................................
! ng - number of data files to be used simultaneously
! file(i) - name of the data file for the GNUPLOT, i=1,...,ng;
! lenght(i) - lenght of the data files names;
! title_pl(i) - subtitle for each curve
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!  ShowTotDos = 'O' - only total DOS
!  ShowTotDos = 'Y' - total DOS shown alongside the PDOS
!  ShowTotDos = 'N' - total not DOS shown alongside the PDOS
!......................................................................
!
implicit none
character title*50,xaxis*20,yaxis*20,where*6
character file(ng)*12,f_gr*12,title_pl(ng)*8,cha1,cha2,cha*2
integer lenght(ng),ng,l3,i
real*8 Up_lim,Up_lim_E
character ShowTotDos

!.......... name of the file for the group (is used if ng.ne.1)
      f_gr='dos.dat0'//file(1)(9:lenght(1))
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis
300   format ('set data style lines'/ &
              'set format y "%.4f"'/ &
              'set key'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')

!.......... axes limits
      if(Up_lim.ne.0.0) then
         write(l3,306) Up_lim
306      format('set yrange [0:',f10.5,']')
      end if
      if(Up_lim_E.ne.0.0) then
         write(l3,307) Up_lim_E
307      format('set xrange [:',f10.5,']')
      end if
!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
         write (l3,401)
401      format ('#set terminal post color "Times-Roman" 14')
         if(ng.eq.1) then
            write(l3,'(a)')'#set output "'//file(1)(1:lenght(1))//'.ps"'
         else
            write(l3,'(a)')'#set output "'//f_gr(1:lenght(1))//'.ps"'
         end if
      else
         write (l3,301)
301      format ('set terminal post color "Times-Roman" 14')
         if(ng.eq.1) then
            write (l3,'(a)')'set output "'//file(1)(1:lenght(1))//'.ps"'
         else
            write (l3,'(a)')'set output "'//f_gr(1:lenght(1))//'.ps"'
         end if
      end if
!.......... common part: .........................................
!____ we show columns cha1 and cha2 for DOS and LDOS, respectively,

      cha='" '
      cha1='2'
      cha2='3'
!_____ plot commands: show DOS only for the first file
      if(ShowTotDos.eq.'O') then
         write (l3,*)  'plot "'//file(1)(1:lenght(1))//cha// &
              ' u 1:'//cha1//' t "total"'//char(92)
      else if(ShowTotDos.eq.'Y') then
         write (l3,*)  'plot "'//file(1)(1:lenght(1))//cha// &
              ' u 1:'//cha1//' t "total"'//char(92)
         write (l3,*)', "'//file(1)(1:lenght(1))//cha//' u 1:'// &
              cha2//' t "'//title_pl(1)//'(loc)"'//char(92)
      else if(ShowTotDos.eq.'N') then
         write (l3,*)'plot "'//file(1)(1:lenght(1))//cha//' u 1:'// &
              cha2//' t "'//title_pl(1)//'(loc)"'//char(92)
      end if
      if(ng.eq.1) go to 10
      do i=2,ng
         write (l3,*)', "'//file(i)(1:lenght(i))//cha//' u 1:'// &
              cha2//' t "'//title_pl(i)//'(loc)"'//char(92)
      end do
10    write(l3,*)' '
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
!.......... Check postscript file if it was specified ................
      if(where.ne.'Screen') then
         if(ng.ne.1) then
            write(6,*)'The postscript file '// &
                 f_gr(1:lenght(1))//'.ps" has been created!'
         else
            write(6,*)'The postscript file '// &
                 file(1)(1:lenght(1))//'.ps" has been created!'
         end if
      end if
end subroutine Plot2

subroutine Plot_sm(ng,file1,lenght1,file,lenght,title,title_pl, &
                        xaxis,yaxis,where,l3,Smeared,Up_lim,Up_lim_E)
!......................................................................
!   Plot curves eighter on the screen or on to a ps-file
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
!......................................................................
! ng - number of data files to be used simultaneously
! file(i) - name of the data file for the GNUPLOT, i=1,...,ng;
! lenght(i) - lenght of the data files names;
! title_pl(i) - subtitle for each curve
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!   If Smeared='F', we show columns 2 for the raw DOS from file and,
! at the same time, all smeared DOS from file//'_'.
!   If Smeared='S', then we show columns 2 and 3 for the raw and smeared
! DOS, respectively.
!......................................................................
!
implicit none
character title*50,xaxis*20,yaxis*20,where*6,Smeared
character file(ng)*12,f_gr*12,title_pl(ng)*7,cha1,cha2,cha*2
character file1*12
integer lenght(ng),ng,lenght1,l3,i
real*8 Up_lim,Up_lim_E
!.......... name of the file for the group (is used if ng.ne.1)
      f_gr='dos.dat0'//file(1)(9:lenght(1))
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis
300   format ('set data style lines'/ &
              'set format y "%.4f"'/ &
              'set key'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')

!.......... axes limits
      if(Up_lim.ne.0.0) then
        write(l3,306) Up_lim
 306    format('set yrange [0:',f10.5,']')
      end if
      if(Up_lim_E.ne.0.0) then
        write(l3,307) Up_lim_E
 307    format('set xrange [:',f10.5,']')
      end if

!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
      else
         write (l3,301)
 301     format ('set terminal post  "Times-Roman" 14')
! 301     format ('set terminal post portrait "Times-Roman" 14'/
!     &        'set size 0.7,1.4')
         if(ng.eq.1) then
           write (l3,'(a)')'set output "'//file(1)(1:lenght(1))//'.ps"'
         else
           write (l3,'(a)')'set output "'//f_gr(1:lenght(1))//'.ps"'
         end if
      end if
!.......... common part: .........................................
!____ we show columns cha1 and cha2 depending on Smeared='F','S'
      if(Smeared.eq.'F') then
        cha='_"'
        cha1='2'
        cha2='2'
      else
        cha='" '
        cha1='2'
        cha2='3'
      end if
!_____ plot commands
      write (l3,*)  'plot "'//file1(1:lenght1)//'" u 1:'// &
                cha1//' t "original"'//char(92)
      write (l3,*)', "'//file(1)(1:lenght(1))//cha//' u 1:'// &
                cha2//' t "'//title_pl(1)//'(smr)"'//char(92)
      if(ng.eq.1) go to 10
      do i=2,ng
        write (l3,*)', "'//file(i)(1:lenght(i))//cha//' u 1:'// &
                cha2//' t "'//title_pl(i)//'(smr)"'//char(92)
      end do
 10   write(l3,*)' '
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
!.......... Check postscript file if it was specified ................
      if(where.ne.'Screen') then
        if(ng.ne.1) then
          write(6,*)'The postscript file '// &
                      f_gr(1:lenght(1))//'.ps" has been created!'
        else
          write(6,*)'The postscript file '// &
                      file(1)(1:lenght(1))//'.ps" has been created!'
        end if
      end if
!     call system ('gs '//file(1:lenght)//'.ps')
      return
end subroutine Plot_sm

subroutine Plot3d(file,lenght,title,xaxis,yaxis,zaxis,where,l3, &
                                   nclasses,nresol_prv,type_prv)
!......................................................................
!   Plot the 3D density
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! zaxis - the title of the z-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
! type_prv = 'contour' - it makes a contour plot;
!          = '3dimens' - it makes a 3D plot;
! nresol_prv - number of grid lines in X and Y directions in 'test.dat'
!......................................................................
! file - name of the data file for the GNUPLOT;
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!
implicit none
integer lenght,l3,nclasses,nresol_prv
character title*50,xaxis*20,yaxis*20,where*6,zaxis*20,type_prv*7
character file*8
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis,zaxis
300   format ('set data style lines'/ &
              'set parametric'/ &
              '#set nokey'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0'/ &
              'set zlabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')
!
!.......... write differently for the contour plot
      if(type_prv.eq.'contour') then
         write (l3,200) nclasses
200      format ('set view 0,0,1,1'/ &
                 'set nosurface'/ &
                 'set contour'/&
                 'set cntrparam level ',i5 )
      end if
!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
         write (l3,302)
 302     format ('#set terminal post color "Times-Roman" 14')
      else
         write (l3,301)
 301     format ('set terminal post color "Times-Roman" 14')
         write (l3,'(a)')'set output "'//file(1:lenght)//'.ps"'
      end if
!.......... common part: plot commands ..............................
      write (l3,*)'splot "'//file(1:lenght)//'" u 1:2:3 '
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
!.......... Check postscript file if it was specified ................
      if(where.ne.'Screen') call system ('gs '//file(1:lenght)//'.ps')
end subroutine Plot3d

subroutine Plot_Explr(file,lenght,title,xaxis,yaxis,where,l3)
!......................................................................
!   Plot the DOS curves eighter on the screen or on to a ps-file
!......................................................................
! title - the title of the plot (char*50)
! xaxis - the title of the x-axis (char*20)
! yaxis - the title of the y-axis (char*20)
! where - where to print: on the 'Screen' or on the 'Postsc' (char*6)
!......................................................................
! file - name of the data file for the GNUPLOT;
! test.gnu - standart name of the file  used by the GNUPLOT;
! file.ps  - name of the postscript file prepared by GNUPLOT.
! l3 - reserved UNIT for the files used here.
!......................................................................
!
implicit none
integer lenght,l3
character title*50,xaxis*20,yaxis*20,where*6
character file*8
!.......... common part for both regimes ............................
      call system('rm test.gnu')
      open (l3, file='test.gnu', status='new')
      write (l3,300) xaxis,yaxis
300   format ('set data style lines'/ &
              'set format y "%.4f"'/ &
              'set key'/ &
              'set xlabel "',a,'" 0,0'/ &
              'set ylabel "',a,'" 0,0')
      write (l3,305) title
305   format ('set title "', a, '" 0,0')
!.......... distinct part ...........................................
      if(where.eq.'Screen') then
         write (l3,*) 'set terminal x11'
         write (l3,'(a)')'#set output "'//file(1:lenght)//'.eps"'
         write(l3,*)'#set terminal post color "Times-Roman" 14'
      else
         write (l3,301)
 301     format ('set terminal post color "Times-Roman" 14')
! 301     format ('set terminal post portrait "Times-Roman" 14'/
!     &        'set size 0.7,1.4')
         write (l3,'(a)')'set output "'//file(1:lenght)//'.eps"'
      end if
!.......... common part: plot commands ..............................
      write (l3,*)'plot "'//file(1:lenght)//'" u 1:2 t "A1"'//char(92)
      write (l3,*)', "'//file(1:lenght)//'" u 3:4 t "A2"'//char(92)
      write (l3,*)', "'//file(1:lenght)//'" u 5:6 t "A3"'
      if(where.eq.'Screen') write (l3,*)'pause -1'
      close (l3)
      call system ('gnuplot test.gnu')
end subroutine Plot_Explr
