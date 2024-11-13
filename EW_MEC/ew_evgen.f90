program ew_eventgen
   use event_module
   use mc_module
   use dirac_matrices
   use mathtool
   use mympi
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=938.0d0,xmpi=139.d0,xmmu=105.658357
   integer :: clocks(2), count_rate, seeds(2)
   real*8 :: progress,ti,tf, xsec_acc
   integer*4 :: nw,nZ,xA,i_fg,j,ilept,gen_events,num_events,nwlk,isospin
   integer*4 :: gen_events_perproc, ierr
   real*8 :: wmax,enu,thetalept,xpf,hw,sig,sig_err
   real*8 :: xmlept,start,finish,total_sig,total_sig_err
   integer*8, allocatable :: irn_int(:),irn_event(:),irn_int0(:),irn_event0(:)
   integer*8 :: i
   character*50 :: fname,isospin_char,en_char,temp_fname
   character*40 :: nk_fname
   character*200 :: command

   type(event_container_t) :: saved_events

   call init0()

   if (myrank().eq.0) then

      !Getting a random seed from my computer
      call system_clock(count=clocks(1),count_rate=count_rate)
      call sleep(1)
      call system_clock(count=clocks(2),count_rate=count_rate)

      seeds(1)=mod(clocks(1),100000)
      seeds(2)=mod(clocks(2),100000)
      print*,'Seed 1: ', seeds(1)
      print*,'Seed 2: ', seeds(2)

      read(5,*) gen_events
      read(5,*) nwlk
      read(5,*) enu 
      read(5,*) thetalept
      read(5,*) isospin
      read(5,*) xsec_acc
      read(5,*) ilept
      read(5,*) xpf
      read(5,*) nZ,xA
      read(5,*) i_fg

      write(en_char,'(i5)') int(enu)
      en_char=adjustl(en_char)
      en_char=trim(en_char)  

      write(isospin_char, '(I0)') isospin

      !fname='test_'//trim(isospin_char)//'_FG_700_31p79_w100_RMF.out'
      fname='test_PS_FG_700_31p79_w100_RMF.out'
      !fname='C12_CC_'//trim(en_char)//'_parallel.out'
      fname=trim(fname)
      print*,'Outputfile: ', fname
      !open(unit=7, file=fname, status='replace')
   endif

   write(temp_fname,'(A,I0,A)') 'process_', myrank(), '.out'
   open(unit=11+myrank(), file=temp_fname, status='replace')

   call bcast(gen_events)
   call bcast(nwlk)
   call bcast(enu)
   call bcast(thetalept)
   call bcast(seeds(1))
   call bcast(seeds(2))
   call bcast(isospin)
   call bcast(xsec_acc)
   call bcast(ilept)
   call bcast(xpf)
   call bcast(nZ)
   call bcast(xA)
   call bcast(i_fg)

   ti=MPI_Wtime()

   allocate(irn_int0(nwlk),irn_event0(nwlk))
   do i=1,nwlk
       irn_int0(i)=seeds(1) + i
       irn_event0(i)=seeds(2) +i
    enddo
    if (myrank().eq.0) then
       print*,'theta lept read in = ', thetalept
       write (6,'(''number of cpus ='',t50,i10)') nproc()
       if (mod(nwlk,nproc()).ne.0) then
          write(6,*)'Error: nwalk must me a multiple of nproc'
          stop
       endif
    endif

    !Convert to radians
    thetalept=thetalept/180.0d0*pi

    nwlk=nwlk/nproc()
    allocate(irn_int(nwlk),irn_event(nwlk))
    irn_int(:)=irn_int0(myrank()*nwlk+1:myrank()*nwlk+nwlk)
    irn_event(:)=irn_event0(myrank()*nwlk+1:myrank()*nwlk+nwlk)

   if(ilept.eq.0) then
      xmlept = 0.0d0
   else
      xmlept = xmmu
   endif

   gen_events_perproc = gen_events/nproc()

   !Including scalar self energy in nucleon mass
   call dirac_matrices_in(xmd,0.8d0*xmn,xmpi,0.0d0,xmlept)

   !Initialize currents and spinors
   call mc_init(gen_events_perproc,xsec_acc,i_fg,irn_int,irn_event, &
         &  nwlk,xpf,xmlept,xA,nZ,isospin)
   num_events = 0

   if(myrank().eq.0) then
      write(6,*) 'Computing total cross section for Ev = ', enu, ' MeV'
   endif

   !Compute the cross section and generate events
   call mc_eval(enu,thetalept,sig,sig_err,saved_events)

   !Print events to temp files
   call print_unweighted_events(saved_events,11+myrank())
   call flush(11+myrank())
   call MPI_Barrier(mpi_comm_world,ierr)
   
   !Concatenate each of the temp files together into output
   if(myrank().eq.0) then
      command = 'cat'
      command = trim(command) // ' process* > ' // trim(fname) 
      print*,command
      call execute_command_line(command)
   endif

   call MPI_Barrier(mpi_comm_world,ierr)

   ! !Delete each of the temp files
   ! if(myrank().eq.0) then
   !    command = 'rm '
   !    do i = 0,nproc()-1
   !       write(temp_fname,'(A,I0,A)') 'process_', i, '.txt'
   !       command = trim(command) // ' ' // trim(temp_fname)
   !    enddo 
   !    call execute_command_line(command)
   ! endif
   
   tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()
   
   !call cpu_time(finish)

   !print '("Time = ",f8.2," seconds.")',finish-start

   contains
      

end program
