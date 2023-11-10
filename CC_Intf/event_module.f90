module event_module
    use mathtool
    implicit none

    type :: particle_t
        real*8 :: p4(4)
        integer*4 :: pdgcode
    endtype 

    type, extends(particle_t) :: event_t 
        type(particle_t), allocatable :: particles(:)
        real*8 :: weight=0
        logical :: unweighted=.FALSE.
    endtype

    type, extends(event_t) :: event_container_t
        type(event_t), allocatable :: events(:)
        real*8 :: max_weight=0
        integer*4 :: num_accepted_events=0
        integer*4 :: num_gen_events=0
    endtype

    contains

        subroutine event_container_init(event_container, numEvents)
            type(event_container_t), intent(inout) :: event_container
            integer*4 :: numEvents

            allocate(event_container%events(numEvents))
        end subroutine event_container_init

        subroutine print_event_container(event_container)
            type(event_container_t), intent(in) :: event_container
            integer*4 :: n,i 
            
            do i=1,event_container%num_accepted_events
                if(event_container%events(i)%weight.ne.0.0d0) then
                    call print_event(event_container%events(i))
                endif 
            enddo 
            write(6,*) 'Maximum weight = ', event_container%max_weight
            write(6,*) 'Accepted events = ', event_container%num_accepted_events

        end subroutine print_event_container

        subroutine print_unweighted_events(event_container,fileunit)
            type(event_container_t), intent(inout) :: event_container
            integer*4 :: n,i,fileunit,eventnum

            eventnum=0
            n = event_container%num_accepted_events
            
            do i=1,n
                if(event_container%events(i)%unweighted) then
                    !call rotate_event(event_container%events(i))
                    call print_unweighted_event(event_container%events(i),&
                        &   eventnum, fileunit)
                    eventnum = eventnum+1
                endif 
            enddo 

            event_container%num_gen_events = eventnum

        end subroutine print_unweighted_events

        subroutine rotate_event(event)
            type(event_t), intent(inout) :: event 
            integer*4 :: n,i 
            real*8 :: angle, rot_matrix(3,3), probe_mag, probe_z

            n = size(event%particles)

            probe_z = event%particles(1)%p4(4)
            probe_mag = sqrt(sum(event%particles(1)%p4(2:4)**2))

            angle = -acos(probe_z/probe_mag)
            !write(6,*)' rotation angle = ', angle*180.0d0/acos(-1.0d0)

            rot_matrix = reshape((/cos(angle), 0.0d0, -sin(angle), 0.0d0, &
                &   1.0d0, 0.0d0, sin(angle), 0.0d0, cos(angle)/),shape(rot_matrix))

            do i=1,n
                event%particles(i)%p4(2:4) = matmul(rot_matrix, event%particles(i)%p4(2:4))
            enddo



        end subroutine rotate_event

        subroutine event_unweight(event_container)
            type(event_container_t), intent(inout) :: event_container
            integer*4 :: i
            real*8 :: ratio,r
            real*8 :: max

            max = event_container%max_weight

            do i=1,event_container%num_accepted_events
                ratio = event_container%events(i)%weight/max
                r = ran()
                if(r.le.ratio) event_container%events(i)%unweighted=.TRUE.
            enddo

        end subroutine event_unweight

        subroutine event_init(event,numPart)
            type(event_t), intent(inout) :: event
            integer*4 :: numPart 

            allocate(event%particles(numPart))
        end subroutine event_init

        subroutine print_event(event)
            type(event_t), intent(in) :: event 
            integer*4 :: n,i

            n = size(event%particles)

            write(6,*) 'weight = ', event%weight
            write(6,*) 'unweighted = ', event%unweighted
            do i=1,n
                !write(6,*) 'i = ', i
                write(6,*) '    P4 = ', event%particles(i)%p4 
            enddo
        end subroutine print_event

        subroutine print_unweighted_event(event,eventnum,fileunit)
            type(event_t), intent(in) :: event 
            integer*4 :: n,i,fileunit,eventnum

            n = size(event%particles)

            !write(fileunit,*) eventnum
            do i=1,n
                write(fileunit,*) event%particles(i)%pdgcode, event%particles(i)%p4 
            enddo
        end subroutine print_unweighted_event
end module

! program test
!   use event_module
!   implicit none

!   type(event_container_t) my_events 
!   type(event_t) :: my_event
!   type(particle_t) :: particles(4)

!   call event_container_init(my_events, 2)

!   call event_init(my_event, 4)
!   particles(1)%p4 = (/1.0, 0.0, 2.0, 3.0/)
!   particles(2)%p4 = (/3.0, 0.0, 2.0, 0.0/)
!   particles(3)%p4 = (/1.0, -2.0, 2.0, 3.0/)
!   particles(4)%p4 = (/1.0, 7.0, 2.0, 5.0/)
    
!   my_event%particles = particles
!   my_event%weight = 1.0e-26

!   my_events%events(1) = my_event


!   particles(1)%p4 = (/-1.0, 0.0, -2.0, -3.0/)
!   particles(2)%p4 = (/-3.0, 0.0, -2.0, -0.0/)
!   particles(3)%p4 = (/-1.0, -2.0, -2.0, -3.0/)
!   particles(4)%p4 = (/-1.0, -7.0, -2.0, -5.0/)
    
!   my_event%particles = particles
!   my_event%weight = 1.0e-26

!   my_events%events(2) = my_event

!   call print_event_container(my_events)

! end program test


