program ew_eventgen
   use event_module
   use mc_module
   use dirac_matrices
   use mathtool
   !use mympi
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=938.0d0,xmpi=139.d0,xmmu=105.658357
   real*8 :: progress
   integer*4 :: nw,nev,nZ,xA,i_fg,np0,np,ne,j,ilept,gen_events,num_events,nwlk
   real*8 :: wmax,enu,thetalept,xpf,Q2min,hw,sig,sig_err,xmlept,start,finish,total_sig,total_sig_err
   integer*8, allocatable :: irn_int(:),irn_event(:)
   integer*8 :: irn0,irn1,i
   character*50 :: fname,intf_char,en_char  
   character*40 :: nk_fname
   type(event_container_t) :: saved_events

   !call init0()

   !if (myrank().eq.0) then
      read(5,*) gen_events
      read(5,*) nev
      read(5,*) nwlk
      read(5,*) enu 
      read(5,*) irn0
      read(5,*) irn1
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
      fname='C12_CC_'//trim(en_char)//'new.out'
      fname=trim(fname)
      open(unit=7, file=fname)
   !endif

   !call bcast(gen_events)
   !call bcast(nev)
   !call bcast(enu)
   !call bcast(irn)
   !call bcast(ilept)
   !call bcast(xpf)
   !call bcast(Q2min)
   !call bcast(nZ,xA)
   !call bcast(i_fg)
   !call bcast(np0)
   !call bcast(np)
   !call bcast(ne)

   !allocate(irn0(nwlk))
   !do i=1,nwlk
   !    irn0(i)=19+i
   ! enddo
   ! if (myrank().eq.0) then
   !    write (6,'(''number of cpus ='',t50,i10)') nproc()
   !    if (mod(nwlk,nproc()).ne.0) then
   !       write(6,*)'Error: nwalk must me a multiple of nproc'
   !       stop
   !    endif
   ! endif
   ! nwlk=nwlk/nproc()
    allocate(irn_int(nwlk),irn_event(nwlk))
    do i=1,nwlk
      irn_int(i)=irn0 + i
      irn_event(i)=irn1 + i
   enddo

   if(ilept.eq.0) then
      xmlept = 0.0d0
   else
      xmlept = xmmu
   endif


   call dirac_matrices_in(xmn, 0.0d0, xmlept)

   !Initialize currents and spinors
   call mc_init(gen_events,i_fg,irn_int,irn_event, &
         &  nev,nwlk,xpf,xmlept,xA,nZ,np,np0,ne,Q2min)
   num_events = 0

   write(6,*) 'Computing total cross section for Ev = ', enu, ' MeV'
   call mc_eval(enu,sig,sig_err,saved_events)
   call print_unweighted_events(saved_events,7)
   num_events = saved_events%num_gen_events

   write(6,*)'Total cross section = ', sig,' +/- ', sig_err
   write(6,*)'Total number of events generated = ', num_events
   
   call cpu_time(finish)

   print '("Time = ",f8.2," seconds.")',finish-start

   contains
      

end program
