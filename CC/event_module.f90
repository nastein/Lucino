module event_module
    use mathtool
    implicit none

    type :: particle_t
        real*8 :: p4(4)
    endtype 

    type :: event_t 
        type(particle_t), allocatable :: particles(:)
        real*8 :: weight=0
        logical :: unweighted=.FALSE.
    endtype

    type :: event_container_t
        type(event_t), allocatable :: events(:)
        real*8 :: max_weight=0
        integer*4 :: size = 0
        integer*4 :: num_gen_events=0
        integer*4 :: capacity = 0
    contains
        procedure :: add_event => vector_add_event
    endtype

    contains

        subroutine vector_add_event(this, event)
            class(event_container_t), intent(inout) :: this 
            type(event_t), intent(in) :: event 

            if (this%size.eq.this%capacity) then
                call increase_event_container_capacity(this)
            endif
            !print*,'Added event'
            this%size = this%size + 1
            this%events(this%size) = event
        end subroutine vector_add_event

        subroutine increase_event_container_capacity(this)
            type(event_container_t), intent(inout) :: this 
            type(event_t), allocatable :: temp(:)
            integer*4 :: new_capacity

            if(this%capacity.eq.0) then
                new_capacity = 1
            else
                new_capacity = this%capacity * 2
            endif

            allocate(temp(new_capacity))
            if (this%size.gt.0) then
                temp(1:this%size) = this%events(1:this%size)
            endif

            call move_alloc(temp, this%events)
            this%capacity = new_capacity
        end subroutine increase_event_container_capacity

        subroutine move_alloc(source, dest)
            type(event_t), allocatable, intent(inout) :: source(:)
            type(event_t), allocatable, intent(out) :: dest(:)
            dest = source
            deallocate(source)
        end subroutine move_alloc

        subroutine print_unweighted_events(this,fileunit)
            type(event_container_t), intent(inout) :: this
            integer*4 :: i,fileunit,eventnum

            eventnum=0
            
            do i=1,this%size
                if(this%events(i)%unweighted) then
                    call print_event(this%events(i),fileunit)
                    eventnum = eventnum+1
                endif 
            enddo 

            this%num_gen_events = eventnum

        end subroutine print_unweighted_events

        subroutine unweight_events(this)
            type(event_container_t), intent(inout) :: this
            integer*4 :: i
            real*8 :: ratio,r

            do i=1,this%size
                ratio = this%events(i)%weight/this%max_weight
                r = ran()

                if(r.le.ratio) then
                    this%events(i)%unweighted=.TRUE.
                    print*,'r = ', r
                    print*,'ratio = ', ratio
                endif
            enddo

        end subroutine unweight_events

        subroutine event_init(event,numPart)
            type(event_t), intent(inout) :: event
            integer*4 :: numPart 

            allocate(event%particles(numPart))
        end subroutine event_init

        subroutine print_event(event,fileunit)
            type(event_t), intent(in) :: event 
            integer*4 :: n,i,fileunit

            n = size(event%particles)

            do i=1,n
                write(fileunit,*) event%particles(i)%p4 
            enddo
        end subroutine print_event

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


