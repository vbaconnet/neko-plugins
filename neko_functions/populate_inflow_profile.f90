module populate_inflow_profile

use neko
use matrix, only: matrix_t

contains

  !> Reads in a csv file and counts lines+columns to automatically
  !! allocate the corresponding matrix_t object.
  subroutine populate_inflow(file_name, INFLOW_PROFILE)
    character(len=*), intent(in) :: file_name
    type(matrix_t), intent(inout) :: INFLOW_PROFILE

    type(file_t) :: f
    integer :: unit, i, ierr, num_columns, num_lines, ios
    character(len=1000) :: line
    character(len=1) :: delimiter
    character(len=LOG_SIZE) :: log_buf
    delimiter = ','


    ! 
    ! Count lines and columns for initialization
    !
    open(unit=unit, file=trim(file_name), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      call neko_error("Error opening file " // trim(file_name))
    end if
    
    num_columns = 1
    num_lines = 0

    ! Read the file line by line
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        ! If it's the first line, count the columns
        if (num_columns .eq. 1) then

          ! Count the number of delimiters in the line
          do i = 1, len_trim(line)
              if (line(i:i) == delimiter) then
                  num_columns = num_columns + 1
              end if
          end do

        end if ! if num_columns .eq. 1

        num_lines = num_lines + 1
    end do
    
    ! Close the file
    close(unit)

    write (log_buf, '(A,A,A,I5,I5)') "Reading file ", trim(file_name), " of size ", num_lines, num_columns
    call neko_log%message(log_buf)

    !
    ! Reading of the inflow profile
    !
    call INFLOW_PROFILE%init(num_lines, num_columns)

    ! Read csv file and broadcast to all ranks since csv is only read in serial  
    f = file_t(trim(file_name))
    call f%read(INFLOW_PROFILE)
    call MPI_Bcast(INFLOW_PROFILE%x , INFLOW_PROFILE%n , MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)
    call file_free(f)

  end subroutine populate_inflow
 end module populate_inflow_profile
