program ew_eventgen
   use event_module
   use mc_module
   use dirac_matrices
   use mathtool
   use mympi
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=938.0d0,xmpi=139.d0,xmmu=105.658357
   real*8 :: progress,ti,tf, xsec_acc  
   integer :: clocks(2), count_rate, seeds(2)
   integer*4 :: nw,nZ,xA,i_fg,Deltapropfull,DeltaPot,j,ilept,gen_events,num_events,nwlk
   integer*4 :: gen_events_perproc, ierr
   integer*4 :: local_trials, global_trials, local_events, global_events
   real*8 :: wmax,enu,thetalept,Eshift,xpf,hw,sig,sig_err,omega,qval,thetaproton,phiproton
   real*8 :: xmlept,start,finish,total_sig,total_sig_err
   integer*8, allocatable :: irn_int(:),irn_event(:),irn_int0(:),irn_event0(:)
   integer*8 :: ran1,ran2,i,idx
   character*50 :: intf_char,temp_fname
   character*50 :: nk_fname,int_string,FG_string,Delta_string
   character*200 :: command, sig_char, fname
   logical :: CC

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
      read(5,*) qval
      read(5,*) omega
      read(5,*) thetaproton
      read(5,*) phiproton
      read(5,*) xpf
      read(5,*) Eshift
      read(5,*) Deltapropfull
      read(5,*) DeltaPot
      read(5,*) nZ,xA
      read(5,*) i_fg
      read(5,*) CC

      if(CC.eqv..true.) then
         int_string = 'EW'
      else
         int_string = 'EM'
      endif

      if(i_fg.eq.1) then
         FG_string = 'FG'
      else
         FG_string = 'SF'
      endif

      if(Deltapropfull.eq.1) then 
         Delta_string = 'FullDeltaProp_'
      else if(Deltapropfull.eq.1 .and.DeltaPot.eq.1) then 
         Delta_string = 'FullDeltaProp_DeltaPot_'
      else if(Deltapropfull.eq.0 .and.DeltaPot.eq.1) then 
         Delta_string = 'DeltaPot_'
      else
         Delta_string = ''
      endif



      write(fname,'(A,A,A,A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'test_resp_',trim(Delta_string),trim(int_string), &
      &  '_',trim(FG_string),'_Ebeam_', int(enu),'_qval_',int(qval), &
      &  '_w_',int(omega),'_ptheta_',int(thetaproton),'_phi_',int(phiproton), '.out'
      if (myrank().eq.0) then
         print*, 'Output file: ', fname
      endif

      fname=trim(fname)
   endif

   !Temporary files for each processor
   write(temp_fname,'(A,I0,A)') 'process_', myrank(), '.out'
   open(unit=11+myrank(), file=temp_fname, status='replace')

   call bcast(gen_events)
   call bcast(nwlk)
   call bcast(enu)
   call bcast(qval)
   call bcast(omega)
   call bcast(thetaproton)
   call bcast(phiproton)
   call bcast(seeds(1))
   call bcast(seeds(2))
   call bcast(xpf)
   call bcast(Eshift)
   call bcast(Deltapropfull)
   call bcast(DeltaPot)
   call bcast(nZ)
   call bcast(xA)
   call bcast(i_fg)
   call bcast(CC)

   ti=MPI_Wtime()
   
   !Do random seed allocation
   allocate(irn_int0(nwlk),irn_event0(nwlk))
   do i=1,nwlk
       irn_int0(i)=seeds(1) + i
       irn_event0(i)=seeds(2) +i
    enddo
    if (myrank().eq.0) then
       write (6,'(''number of cpus ='',t50,i10)') nproc()
       if (mod(nwlk,nproc()).ne.0) then
          write(6,*)'Error: nwalk must me a multiple of nproc'
          stop
       endif
    endif

    nwlk=nwlk/nproc()
    allocate(irn_int(nwlk),irn_event(nwlk))
    irn_int(:)=irn_int0(myrank()*nwlk+1:myrank()*nwlk+nwlk)
    irn_event(:)=irn_event0(myrank()*nwlk+1:myrank()*nwlk+nwlk)

   !Number of events each processor should generate
   gen_events_perproc = gen_events/nproc()

   if(CC.eqv..true.) then
      xmlept = xmmu
   else
      xmlept = 0.0d0
   endif

   if (myrank().eq.0) then
      write(6,*) 'Using mlept = ', xmlept
   endif

   !Initialize currents module
   call dirac_matrices_in(xmd,xmn,xmpi,0.0d0,xmlept,CC,Deltapropfull,DeltaPot)

   !Initialize spectral function and other necessary inputs
   call mc_init(gen_events_perproc,i_fg,irn_int,irn_event, &
         &  nwlk,xpf,Eshift,xmlept,xA,nZ,CC)
   num_events = 0

   if(myrank().eq.0) then
      write(6,*) 'Computing total cross section for Ev = ', enu, ' MeV', &
      & ', omega = ', omega, ', q = ', qval, &
      & ', thetaproton = ', thetaproton, ', phiproton = ', phiproton
   endif

   !Compute the cross section and generate events
   call mc_eval(enu,qval,thetaproton,phiproton,omega,sig,sig_err,saved_events)

   !Print events to temp files
   call print_unweighted_events(saved_events,11+myrank())

   !Get the total number of trials
   local_trials = saved_events%trials
   call addall(local_trials,global_trials)
   call bcast(global_trials)
   local_events = saved_events%num_gen_events
   call addall(local_events,global_events)
   call bcast(global_events)
   call sleep(1)

   close(11+myrank())
   call MPI_Barrier(mpi_comm_world,ierr)

   !Concatenate each of the temp files together into output
   if(myrank().eq.0) then

      open(unit=1, file=fname, status='replace', action='write')
      write(1,'(ES24.16,1X,I0,1X,I0)') sig, global_trials, global_events
      close(1)
      command = 'cat'
      command = trim(command) // ' process* >> ' // trim(fname) 
      call execute_command_line(command)
   endif

   call MPI_Barrier(mpi_comm_world,ierr)
   
   tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()  
      
end program
