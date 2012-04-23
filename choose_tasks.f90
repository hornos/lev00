subroutine choose_tasks(iTask,tlist,tlen)
use param
use dos
use atoms
implicit none
integer iTask,tlen,len,it,nat,i,ii
character(*) tlist
character type*1,item*2,list*200,chal*2
logical at_done,it_done,TagYes,Yes_Done
character(len=2),dimension(:),allocatable :: Labels
logical,dimension(:),allocatable :: listat
!
!........... species by atom number
   allocate(Labels(NIONS)) ; allocate(listat(NIONS))
   nat=0
   do i=1,NSPEC
      do ii=1,NspN(i)
         nat=nat+1 ; Labels(nat)=Species(i)
      end do
   end do
!
!.....................................................................
!............. START MAIN MENU HERE ..................................
!.....................................................................
!
    type='t' ; listat=.false. 
    Yes_Done=.false. ; at_done=.false. ; it_done=.false.
21  write(*,*)'...................Choose PDOS Tasks .......................'
    write(*,*)'... IMPORTNAT: PDOS for multiple atoms will be summed up ...'
    write(*,*)'............................................................'
    if(.not.at_done) then
       write(*,'(a)')'   A. Choose atom(s): ... undefined ..'
    else
       write(*,'(a)')'   A. Chosen atom(s):'
       write(*,'(a)')'         '//list(1:len)
    end if
    write(*,'(a,a)') '  Tp. Current type of PDOS: ',type
    if(at_done)  then
       if(Yes_Done) then
          write(*,'(a)')   '   L. Produce the list -> Done!'
          write(*,'(7x,a)') tlist(1:tlen)
       else
          write(*,'(a)')   '   L. Produce the list'
       end if
    end if
    if(Yes_Done) write(*,'(a)')   '   E. Empty the list of task: start again!'
    write(*,'(a)')   '  ST. Show all tasks indicating the chosen ones (if any)'
    write(*,'(a)')   '   Q. Quit with the empty list'
    write(*,'(a)')   '   P. Exit: quit with the current list'
    write(*,*)
    write(*,*)'------------>'
    read(*,'(a)') item

![A].............. choose atoms

    IF(trim(item).eq.'A') THEN
       
       call TTag(TagYes,NIONS,TI,Labels,listat,list,len)
       if(TagYes) at_done=.true.

![Tp].............. choose between: t,s,p,d

    ELSE IF(trim(item).eq.'Tp') THEN
       
       write(*,*)'Choose the type of PDOS: t,s,p,d'
1      read(*,'(a1)') type
       if(type.ne.'t' .and. type.ne.'s' .and. type.ne.'p' .and. type.ne.'d') go to 1
       
![L].............. produce the list of tasks

    ELSE IF(trim(item).eq.'L' .and. at_done) THEN
       
       do nat=1,NIONS
          if(listat(nat)) then
             if(type.eq.'t') then
                it=(nat-1)*4+1
             else if(type.eq.'s') then
                it=(nat-1)*4+2
             else if(type.eq.'p') then
                it=(nat-1)*4+3
             else if(type.eq.'d') then
                it=(nat-1)*4+4
             end if
             which_task(it)=.true.
          end if
       end do
       call do_nice_list(which_task,Ntask,tlist,tlen)
       Yes_Done=.true. 
 
![E].............. empty the current list

    ELSE IF(trim(item).eq.'E' .and. Yes_Done) THEN
       
       which_task=.false. ; Yes_Done=.false.

![ST].............. show all tasks

    ELSE IF(trim(item).eq.'ST') THEN

        do it=1,Ntask
           chal=' N' ; if(which_task(it)) chal=' Y'
           if(mod(it+3,4).eq.0) then
              write(*,'(2(a,i5),a,a,a,f5.2,a,3(f7.3,1x),a)') &
                   'at= ',v_atm(it),' '//Species(v_spec(it))//' ',&
                   it,' task {'//Flag(it)//'}', &
                   chal,' Radius= ',Radius(it), ' center= (', &
                   (Point(i,it),i=1,3),') '//method(it)
           else
              write(*,'(13x,i5,2a)') it,' task {'//Flag(it)//'}',chal
           end if
        end do
        write(*,*)'Press ENTER when done ...' ; read(*,*)

![Q].............. quit with the empty list

    ELSE IF(trim(item).eq.'Q') THEN

       which_task=.false. ; iTask=0 ; return

![P].............. exit with the current list

   ELSE IF(trim(item).eq.'P') THEN
      if(Yes_Done) iTask=1 ; return
   ELSE
      write(*,*)'Error! Try again!'
   END IF
   go to 21
end subroutine choose_tasks







