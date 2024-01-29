subroutine Open_IO
    implicit none
    
    character(LEN=10),PARAMETER:: u = 'unknown', f = 'formatted'

    integer:: i

    open( 20, file = 'OUTPUT.DAT',       status=u, form=f )
    OPEN( 25, FILE = 'TENSOR_RATES.DAT', STATUS=U, FORM=F )
    OPEN( 30, FILE = 'JUNK.DAT'  ,       STATUS=U, FORM=F )
    OPEN( 35, FILE = 'STM_HISTORY.DAT'  ,       STATUS=U, FORM=F )

end SUBROUTINE OPEN_IO

!=======================================================

