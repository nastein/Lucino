program ew_eventgen
   use event_module
   use constants
   use fast_flux_module
   use mc_module
   use dirac_matrices
   use mathtool
   use mympi
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0)
   real*8 :: progress,ti,tf, xsec_acc  
   integer :: clocks(2), count_rate, seeds(2)
   integer*4 :: nw,nZ,xA,i_mode,j,ilept,gen_events,num_events,nwlk,isospin
   integer*4 :: unweight_mode
   integer*4 :: DeltaPropFull,DeltaProp3half,DeltaPot,intfsign,np_del
   integer*4 :: gen_events_perproc, ierr,nenu,bin
   integer*4 :: local_trials, global_trials, local_events, global_events
   real*8 :: wmax,thetalept,xpf,Eshift,hw,sig,sig_err,enu_1,enu_2,flux_norm
   real*8 :: xmlept,start,finish,total_sig,total_sig_err, enu_max,henu
   real*8, allocatable :: enu_v(:),flux_v(:),enu_width(:),pdel(:),pot_del(:)
   integer*8, allocatable :: irn_int(:),irn_event(:),irn_int0(:),irn_event0(:)
   integer*8 :: ran1,ran2,i,idx
   character*50 :: intf_string,temp_fname
   character*40 :: nk_fname,int_string,FG_string
   character*200 :: command,sig_char,theta_str,fname,flux_file
   logical :: CC, rotate_beam_along_z

   type(event_container_t) :: saved_events

   call init0()

   call read_physics_constants('MEC_constants.dat')

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
      read(5,*) unweight_mode
      read(5,*) xpf
      read(5,*) Eshift
      read(5,*) DeltaPropFull
      read(5,*) DeltaProp3half
      read(5,*) DeltaPot
      read(5,*) intfsign
      read(5,*) nZ,xA
      read(5,*) i_mode
      read(5,*) rotate_beam_along_z
      read(5,*) CC
      read(5,*) flux_file
      close(5)

      
      if(CC.eqv..true.) then
         int_string = 'EW'
      else
         int_string = 'EM'
      endif

      if(i_mode.eq.0) then
         FG_string = 'FG'
      else if(i_mode.eq.1) then
         FG_string = 'SFp1p2'
      else
         FG_string = 'SFqQ'
      endif

      if(intfsign.eq.1) then
         intf_string = 'Noemi'
      else
         intf_string = 'Amaro'
      endif

      write(fname,'(A,A,A,A,A,A,A)') 'test_',trim(int_string),'_',trim(FG_string),'_',trim(intf_string),'_T2K.out'
      
      print*, 'Output file: ', fname

      !write(fname, '(A)') 'test.out'
      

      fname=trim(fname)

         !.......flux folded cross section
      open(unit=5,file=flux_file,status='unknown',form='formatted') 
      read(5,*) nenu 

      allocate(enu_v(nenu),flux_v(nenu))
      do i=1,nenu
         read(5,*) bin, enu_1, enu_2,flux_v(i)
         enu_v(i)= 0.5d0*(enu_1+enu_2)
      enddo   
      close(5)
      enu_v=enu_v*1.e3
      flux_v=flux_v*1.e-3
      henu=enu_v(2)-enu_v(1)
      flux_norm=sum(flux_v(:))*henu
      write(6,*) 'The normalization of the flux [10^-5/m^2] is', flux_norm
      enu_max=enu_v(nenu)

      write(6,*)'Reading in Delta Potential'
      open(10, file='rho_1.dat')
      read(10,*) np_del
      allocate(pdel(np_del),pot_del(np_del))
      do i=1,np_del
         read(10,*) pdel(i),pot_del(i)
      enddo
   endif

   !Temporary files for each processor
   write(temp_fname,'(A,I0,A)') 'process_', myrank(), '.out'
   open(unit=11+myrank(), file=temp_fname, status='replace')

   call bcast(gen_events)
   call bcast(nwlk)
   call bcast(unweight_mode)
   call bcast(seeds(1))
   call bcast(seeds(2))
   call bcast(ilept)
   call bcast(xpf)
   call bcast(Eshift)
   call bcast(DeltaPropFull)
   call bcast(DeltaProp3half)
   call bcast(DeltaPot)
   call bcast(intfsign)
   call bcast(nZ)
   call bcast(xA)
   call bcast(i_mode)
   call bcast(rotate_beam_along_z)
   call bcast(CC)
   call bcast(flux_norm)
   call bcast(nenu)
   if(myrank().ne.0) then
      allocate(enu_v(nenu),flux_v(nenu))
   endif
   call bcast(enu_v)
   call bcast(flux_v)
   call bcast(enu_max)
   call bcast(np_del)
   if(myrank().ne.0) then
      allocate(pot_del(np_del),pdel(np_del))
   endif
   call bcast(pot_del)
   call bcast(pdel)

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

   call init_flux_lookup(enu_v,flux_v,nenu)

   if(CC.eqv..true.) then
      xmlept = xmmu
   else
      xmlept = 0.0d0
   endif

   if (myrank().eq.0) then
      write(6,*) 'Using mlept = ', xmlept
   endif

   !Initialize currents module
   call dirac_matrices_in(xmd,xmn,xmpi,0.0d0,xmlept,CC,DeltaPropFull,DeltaProp3half,DeltaPot,&
      &  intfsign,np_del,pdel,pot_del)

   !Initialize spectral function and other necessary inputs
   call mc_init(gen_events_perproc,unweight_mode,xsec_acc,i_mode,irn_int, &
         &  irn_event,nwlk,xpf,Eshift,xmlept,xA,nZ,CC,rotate_beam_along_z)

   call set_up_flux(flux_v,enu_v,enu_max,flux_norm,nenu)
   num_events = 0

   !Compute the cross section and generate events
   call mc_eval(sig,sig_err,saved_events)

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

      ! Delete temporary per-process files
      command = 'rm -f process_*.out'
      call execute_command_line(command)
   endif

   call MPI_Barrier(mpi_comm_world,ierr)
   
   tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()  
      
end program
