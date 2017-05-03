module co2
    !--------------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------------

    use constants
    use namelist_define
    use mapinfo_define
    use obs_define
    use mpi_module
    use map_utils
    use netcdf
    use wrf_tools

    implicit none

contains
    !--------------------------------------------------------------------------
    ! get_co2_tower_obs
    !
    ! Read CO2 tower files for the current time into raw%co2_tower variables.
    ! 
    ! Input files:
    !   tower_%Y%m%d%H%M%S.dat (e.g. tower_20170503103800.dat)
    !
    !--------------------------------------------------------------------------
    subroutine get_co2_tower_obs(times, ix, jx, kx, p, ph, proj)
        implicit none

        character (len=80), intent(in)              :: times
        integer, intent(in)                         :: ix, jx, kx
        real, dimension(ix, jx, kx), intent(in)     :: p
        real, dimension(ix, jx, kx+1), intent(in)   :: ph
        type(proj_info), intent(in)                 :: proj

        integer                                     :: i
        character (len=80)                          :: input_file
        integer                                     :: ierr
        integer                                     :: nlines

        character (len=80)                          :: tower_name
        real                                        :: lat, lon, elev, co2

        input_file = 'tower_' // times(1:4) // times(6:7) // times(9:10) // &
                     times(12:13) // times(15:16) // times(18:19) // ".dat"
        open(10, file=trim(input_file), status='old', form='formatted', &
             iostat=ierr)
        if (ierr /= 0) then
            if (my_proc_id == 0) then
                write(*, *) trim(input_file), ' does not exsit.'
            end if
        end if

        if (my_proc_id == 0) then
            write(*, *) ''
            write(*, *) '---------------------------------------------------'
            write(*, *) '.... Getting tower CO2 data ....'
        end if

        nlines = 0
        do
            read(10, *, iostat=ierr)
            if (ierr /= 0) then
                exit
            end if
            nlines = nlines + 1
        end do
        raw%co2_tower%num = nlines

        allocate(raw%co2_tower%latitude(raw%co2_tower%num))
        allocate(raw%co2_tower%longitude(raw%co2_tower%num))
        allocate(raw%co2_tower%elevation(raw%co2_tower%num))
        allocate(raw%co2_tower%co2(raw%co2_tower%num))

        rewind(10)
        do i = 1, raw%co2_tower%num
            read(10, '(a12, 4f12.3)', iostat=ierr) tower_name, lat, lon, &
                                                   elev, co2

            raw%co2_tower%latitude(i) = lat
            raw%co2_tower%longitude(i) = lon
            raw%co2_tower%elevation(i) = elev
            raw%co2_tower%co2(i) = co2

            if (ierr /= 0) then
                exit
            end if
        end do
        close(10)

    end subroutine get_co2_tower_obs

end module co2
