subroutine Close_IO
    implicit none

    close( 20, status = 'KEEP' )
    CLOSE( 25, STATUS = 'KEEP' )
    CLOSE( 30, STATUS = 'KEEP' )
    CLOSE( 35, STATUS = 'KEEP' )

end subroutine Close_IO