program ew_eventgen
   use event_module
   use mc_module
   use dirac_matrices
   use mathtool
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=938.0d0,xmpi=139.d0,xmmu=105.658357
   real*8 :: progress
   integer*4 :: nw,nev,i,nZ,xA,i_fg,np0,np,ne,j,ilept,gen_events,num_events
   real*8 :: wmax,enu,thetalept,xpf,Q2min,hw,sig,sig_err,xmlept,start,finish,total_sig,total_sig_err
   integer*8 :: irn
   character*50 :: fname,intf_char,en_char  
   character*40 :: nk_fname
   type(event_container_t) :: saved_events

   read(5,*) gen_events
   read(5,*) nev
   read(5,*) enu 
   read(5,*) irn
   read(5,*) ilept
   read(5,*) xpf
   read(5,*) Q2min
   read(5,*) nZ,xA
   read(5,*) i_fg
   read(5,*) np0,ne 
   read(5,*) np

   write(en_char,'(i5)') int(enu)
   en_char=adjustl(en_char)
   en_char=trim(en_char)  

   !fname='test.out'
   fname='C12_CCintf_'//trim(en_char)//'.out'
   fname=trim(fname)
   open(unit=7, file=fname)

   if(ilept.eq.0) then
      xmlept = 0.0d0
   else
      xmlept = xmmu
   endif

   call dirac_matrices_in(xmd,xmn,xmpi,0.0d0,xmlept)

   !Initialize currents and spinors

   call mc_init(i_fg,nev,xpf,xmlept,xA,nZ,np,np0,ne,Q2min)
   num_events = 0

   write(6,*) 'Computing total cross section for Ev = ', enu, ' MeV'
   i=0
   total_sig = 0

   do while(num_events.le.gen_events)

      call update_progress_bar(num_events, gen_events)

      call event_container_init(saved_events,nev)
      irn = irn+1
      call mc_eval(enu,sig,sig_err,saved_events,irn)
      call event_unweight(saved_events)
      call print_unweighted_events(saved_events,7)


      total_sig = total_sig + sig
      i = i + 1
      deallocate(saved_events%events)
      num_events = num_events + saved_events%num_gen_events
   enddo

   total_sig = total_sig/dble(i)
   write(6,*)'total cross section = ', total_sig
   write(6,*)'total number of events generated = ', num_events
   
   call cpu_time(finish)


   print '("Time = ",f8.2," seconds.")',finish-start

   contains
      subroutine update_progress_bar(current_step, total_steps)
          integer*4, intent(in) :: current_step, total_steps
          real*8 :: percent_done
          integer*4 :: bar_width, num_hashes

          ! Calculate the progress percentage
          percent_done = real(current_step) / real(total_steps) * 100.0

          ! Calculate the number of hashes to display in the progress bar
          bar_width = 100
          num_hashes = int(percent_done * real(bar_width) / 100.0)

          ! Clear the line and print the progress bar
          write(*, "('Progress: [', A, A, '] ', F3.0, '%')", advance="no") &
            repeat("=", num_hashes), repeat(" ", bar_width - num_hashes), percent_done
          ! Move the cursor to the beginning of the line
          write(*, '(A1)', advance="no") char(13)
          

      end subroutine update_progress_bar


end program
