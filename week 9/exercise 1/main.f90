program main

    !! The "main" routine of the program
    !!
    !! This is where execution starts when running the program

    use processor
    use fea

    implicit none

    ! Read model data (in processor)
    call input

    ! Initialize problem (in fea)
    call initial

    ! Calculate displacements
    !call displ

    ! eigen stuff
    call eigen

end program main
