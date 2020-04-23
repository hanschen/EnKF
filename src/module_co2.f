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
    ! height_to_k
    !
    ! Index of closest model level corresponding to height.
    !
    ! Parameters
    ! ----------
    ! ph:
    !   Array with model geopotential heights (staggered) with the same units
    !   as ``height``.
    ! height:
    !   Desired height above sea level.
    !
    ! Return
    ! ------
    ! height_to_k:
    !   k at height.
    !
    !   Special return values:
    !   - Negative number: k is below the lowest model level (number indicates
    !     the height difference in the units of ``height``).
    !   - 0: k is above the highest model level.
    !
    !--------------------------------------------------------------------------
    function height_to_k(ph, height_msl)
        implicit none

        real                                :: height_to_k
        real, dimension(:), intent(in)      :: ph
        real, intent(in)                    :: height_msl

        integer                             :: k, kx, k1
        real, dimension(:), allocatable     :: ph_unstaggered
        real                                :: height_lower, height_upper
        real                                :: alpha

        kx = size(ph) - 1
        allocate(ph_unstaggered(kx))

        ph_unstaggered = (ph(1:kx) + ph(2:kx+1)) / 2.

        ! Find closest vertical level below 'height_msl'
        do k = 1, kx
            if (ph_unstaggered(k) < height_msl) then
                k1 = k
            end if
        end do

        ! Linearly interpolate to find k at 'height_msl'
        height_lower = ph_unstaggered(k1)
        height_upper = ph_unstaggered(k1+1)

        if (height_msl < height_lower) then
            height_to_k = height_msl - height_lower
        else if (height_msl > height_upper) then
            height_to_k = 0
        else
            alpha = (height_msl - height_lower)/(height_upper - height_lower)
            height_to_k = k1 + alpha
        end if

        ! write(*, *) height_msl, height_lower, height_upper, k1, height_to_k

        return

    end function height_to_k

    !--------------------------------------------------------------------------
    ! get_co2_tower_obs
    !
    ! Read CO2 tower files for the current time into raw%co2_tower variables.
    !
    ! Parameters
    ! ----------
    ! ix:
    !   x (longitude) dimension size.
    ! jx:
    !   y (latitude) dimension size.
    ! kx:
    !   z (vertical) dimension size.
    ! proj:
    !   Model projection (of type proj_info).
    !
    ! Input files
    ! -----------
    ! co2_tower_{cycle}.dat (e.g. co2_tower_0000.dat)
    !
    !--------------------------------------------------------------------------
    subroutine get_co2_tower_obs(ix, jx, kx, proj)
        implicit none

        integer, intent(in)                         :: ix, jx, kx
        type(proj_info), intent(in)                 :: proj

        character (len=80)                          :: wrf_file
        real, dimension(ix,jx,kx+1)                 :: ph, phb
        integer                                     :: icycle
        character (len=4)                           :: cycle_num
        integer                                     :: n, nn
        character (len=80)                          :: input_file
        integer                                     :: ierr
        integer                                     :: nlines, nlines_total

        character (len=80)                          :: tower_name
        real                                        :: lat, lon, height_msl, co2

        real                                        :: aio, ajo, ako
        integer                                     :: io, jo, ko
        integer                                     :: time_window

        time_window = time_window_length
        if (time_window == 0) then
            ! include current time if time_window_length is 0
            time_window = 1
        end if

        allocate(raw%co2_tower%nlines(0:time_window-1))

        nlines_total = 0
        do icycle = 0, time_window-1
            write(cycle_num, '(i4.4)') icycle
            input_file = 'co2_tower_' // cycle_num // ".dat"

            open(10, file=trim(input_file), status='old', form='formatted', &
                 iostat=ierr)

            if (ierr /= 0) then
                if (my_proc_id == 0) then
                    write(*, *) trim(input_file), ' does not exist.'
                    ! TODO: Fix error message
                    write(*, *) 'Tower CO2 concentration data will not be &
                                &assimilated!'
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

            raw%co2_tower%nlines(icycle) = nlines
            nlines_total = nlines_total + nlines

            close(10)
        end do

        raw%co2_tower%num = nlines_total

        allocate(raw%co2_tower%icycle(raw%co2_tower%num))
        allocate(raw%co2_tower%tower_name(raw%co2_tower%num))
        allocate(raw%co2_tower%latitude(raw%co2_tower%num))
        allocate(raw%co2_tower%longitude(raw%co2_tower%num))
        allocate(raw%co2_tower%height_msl(raw%co2_tower%num))
        allocate(raw%co2_tower%co2(raw%co2_tower%num))
        allocate(raw%co2_tower%ii(raw%co2_tower%num))
        allocate(raw%co2_tower%jj(raw%co2_tower%num))
        allocate(raw%co2_tower%kk(raw%co2_tower%num))

        nn = 1

        do icycle = 0, time_window-1
            write(cycle_num, '(i4.4)') icycle
            input_file = 'co2_tower_' // cycle_num // ".dat"

            ! Get staggered heights for each cycle time.
            ! The wrfinput filename is currently hard-coded, we may want to
            ! make more general later.
            wrf_file = 'wrfinput_d01_input_' // cycle_num
            call get_variable3d(trim(wrf_file), 'PH        ', &
                                ix, jx, kx+1, 1, ph )
            call get_variable3d(trim(wrf_file), 'PHB       ', &
                                ix, jx, kx+1, 1, phb)
            ph = (ph + phb)/9.8

            open(10, file=trim(input_file), status='old', form='formatted', &
                 iostat=ierr)

            do n = 1, raw%co2_tower%nlines(icycle)
                read(10, '(a20, f10.4, f10.4, f7.1, f11.6)', &
                     iostat=ierr) tower_name, lat, lon, height_msl, co2

                raw%co2_tower%icycle(nn) = icycle
                raw%co2_tower%tower_name(nn) = tower_name
                raw%co2_tower%latitude(nn) = lat
                raw%co2_tower%longitude(nn) = lon
                raw%co2_tower%height_msl(nn) = height_msl
                raw%co2_tower%co2(nn) = co2

                call latlon_to_ij(proj, lat, lon, aio, ajo)
                io = nint(aio)
                jo = nint(ajo)
                ako = height_to_k(ph(io,jo,:), height_msl)

                raw%co2_tower%ii(nn) = aio
                raw%co2_tower%jj(nn) = ajo
                raw%co2_tower%kk(nn) = ako

                ! Diagnostics
                ! write(*, *) trim(tower_name)
                ! write(*, *) lat
                ! write(*, *) lon
                ! write(*, *) height_msl
                ! write(*, *) co2

                nn = nn + 1
            end do
            close(10)
        end do

        return

    end subroutine get_co2_tower_obs

    !--------------------------------------------------------------------------
    ! sort_co2_tower_data
    !
    ! Prepare CO2 tower data to be assimilated.
    !
    ! Parameters
    ! ----------
    ! ix:
    !   x (longitude) dimension size.
    ! jx:
    !   y (latitude) dimension size.
    ! datathin:
    !   Ignore every ``datathin``th tower (NOT IMPLEMENTED).
    !   For example, if ``datathin`` is 2, then every second tower is ignored.
    ! hroi:
    !   Horizontal radius of influence in gridpoints.
    ! vroi:
    !   Vertical radius of influence in gridpoints.
    !
    !--------------------------------------------------------------------------
    subroutine sort_co2_tower_data(ix, jx, datathin, hroi, vroi)
        implicit none

        real, parameter                      :: HEIGHT_TOLERENCE = 100
        integer, intent(in)                  :: ix, jx
        integer, intent(in)                  :: datathin, hroi, vroi

        character (len=4)                    :: cycle_num
        real                                 :: ii, jj, kk
        integer                              :: n

        do n = 1, raw%co2_tower%num
            write(cycle_num, '(i4.4)') raw%co2_tower%icycle(n)

            ii = raw%co2_tower%ii(n)
            jj = raw%co2_tower%jj(n)
            kk = raw%co2_tower%kk(n)

            ! Check if the observation is within the domain
            if ((ii < 1. .or. ii > real(ix)) .or. &
                    (jj < 1. .or. jj > real(jx))) then
                if (my_proc_id == 0) then
                    write(*, *) 'Tower ', trim(raw%co2_tower%tower_name(n)), &
                                ' is outside the domain, skipping observation'
                end if

                cycle
            end if

            ! Check if the vertical position of the observation is reasonable
            if (kk == 0) then
                if (my_proc_id == 0) then
                    write(*, *) 'Tower ', trim(raw%co2_tower%tower_name(n)), &
                                ' is above the highest model level', &
                                ' (how is that even possible?), ', &
                                ' skipping observation'
                end if

                cycle
            else if (kk < 0) then
                if (abs(kk) < HEIGHT_TOLERENCE) then
                    kk = 1
                else
                    if (my_proc_id == 0) then
                        write(*, *) 'Tower ', &
                                    trim(raw%co2_tower%tower_name(n)), &
                                    ' is too far below the lowest model ', &
                                    ' level, skipping observation'
                    end if

                    cycle
                end if
            end if

            if (raw%co2_tower%co2(n) > 0) then
                obs%num                 = obs%num + 1
                obs%dat     (obs%num  ) = raw%co2_tower%co2(n)
                obs%type    (obs%num  ) = 'co2t_' // cycle_num
                obs%err     (obs%num  ) = co2_error_tower
                obs%position(obs%num,1) = ii
                obs%position(obs%num,2) = jj
                obs%position(obs%num,3) = kk
                obs%roi     (obs%num,1) = hroi
                obs%roi     (obs%num,2) = vroi
            endif
        end do

        return

    end subroutine sort_co2_tower_data

    !--------------------------------------------------------------------------
    ! xb_to_co2_tower
    !
    ! Translate background state vector to CO2 concentration value at tower
    ! locations.
    !
    ! The tower value is obtained through linear interpolation horizontally and
    ! vertically.
    !
    ! Parameters
    ! ----------
    ! inputfile:
    !   Path to wrfout file with the background. Used only if 'CO2_{cycle}' is
    !   not a state variable.
    ! xb:
    !   State vector (background or prior).
    ! cycle_num:
    !   Cycle number as a string, e.g. '0000'.
    ! ix:
    !   x (longitude) dimension size.
    ! jx:
    !   y (latitude) dimension size.
    ! kx:
    !   z (vertical) dimension size.
    ! nv:
    !   Number of state variables.
    ! iob:
    !   Number of the current observation.
    !
    ! Returns
    ! -------
    ! hxb:
    !   h(xb) for the current observation.
    !
    !--------------------------------------------------------------------------
    subroutine xb_to_co2_tower(inputfile, xb, cycle_num, ix, jx, kx, nv, iob, hxb)
        implicit none
        character(len=10), intent(in)           :: inputfile
        character(len=4)                        :: cycle_num
        integer, intent(in)                     :: ix, jx, kx, nv, iob
        real, dimension(3,3,kx+1,nv), intent(in) :: xb
        real, intent(out)                       :: hxb

        character(len=8)                        :: varname
        real, dimension(ix,jx,kx)               :: co2
        real, dimension(kx)                     :: co2_vert
        integer                                 :: m
        integer                                 :: i_co2 = 0
        real                                    :: obs_ii, obs_jj, obs_kk
        integer                                 :: i1, j1, k1
        real                                    :: dx, dxm, dy, dym, dz, dzm

        obs_ii = obs%position(iob,1)
        obs_jj = obs%position(iob,2)
        obs_kk = obs%position(iob,3)

        i1 = int(obs_ii)
        j1 = int(obs_jj)
        k1 = int(obs_kk)

        dx = obs_ii - real(i1)
        dy = obs_jj - real(j1)
        dz = obs_kk - real(k1)
        dxm = real(i1+1) - obs_ii
        dym = real(j1+1) - obs_jj
        dzm = real(k1+1) - obs_kk

        varname = 'CO2_' // cycle_num

        do m = 1, nv
            if (trim(enkfvar(m)) == varname) then
                i_co2 = m
            end if
        end do

        if (i_co2 > 0) then
            co2(i1:i1+1,j1:j1+1,1:kx) = xb(1:2,1:2,1:kx,i_co2)
        else
            call get_variable3d(inputfile, varname, ix, jx, kx, 1, co2)
        end if

       co2_vert(:) = dym*(dx*co2(i1+1,j1,:) + dxm*co2(i1,j1,:)) + &
                     dy*(dx*co2(i1+1,j1+1,:) + dxm*co2(i1,j1+1,:))

        hxb = dzm*co2_vert(k1) + dz*co2_vert(k1+1)

        return

    end subroutine xb_to_co2_tower

    !--------------------------------------------------------------------------
    ! get_co2_airborne_obs
    !
    ! Read CO2 airborne files for the current time into raw%co2_airborne
    ! variables.
    !
    ! Parameters
    ! ----------
    ! times:
    !   Time of the observation formatted as %Y-%m-%d_%H:%M:%S, e.g.
    !   2017-05_03_10:38:00.
    ! ph:
    !   Unstaggered model geopotential heights.
    ! proj:
    !   Model projection (of type proj_info).
    !
    ! Input files
    ! -----------
    ! airborne_%Y%m%d%H%M%S.dat (e.g. airborne_20170503103800.dat)
    !
    !--------------------------------------------------------------------------
    subroutine get_co2_airborne_obs(times, ph, proj)
        implicit none

        character(len=80), intent(in)               :: times
        real, dimension(:,:,:), intent(in)          :: ph
        type(proj_info), intent(in)                 :: proj

        integer                                     :: n
        character (len=80)                          :: input_file
        integer                                     :: ierr
        integer                                     :: nlines

        real                                        :: lat, lon, height_msl, co2

        real                                        :: aio, ajo, ako
        integer                                     :: io, jo, ko

        input_file = 'airborne_' // times(1:4) // times(6:7) // &
                     times(9:10) // & times(12:13) // times(15:16) // &
                     times(18:19) // ".dat"
        open(10, file=trim(input_file), status='old', form='formatted', &
             iostat=ierr)
        if (ierr /= 0) then
            if (my_proc_id == 0) then
                write(*, *) trim(input_file), ' does not exist.'
                write(*, *) 'Airborne CO2 concentration data will not be &
                            &assimilated!'
            end if
        end if

        if (my_proc_id == 0) then
            write(*, *) ''
            write(*, *) '---------------------------------------------------'
            write(*, *) '.... Getting airborne CO2 data ....'
        end if

        nlines = 0
        do
            read(10, *, iostat=ierr)
            if (ierr /= 0) then
                exit
            end if
            nlines = nlines + 1
        end do
        raw%co2_airborne%num = nlines

        allocate(raw%co2_airborne%latitude(raw%co2_airborne%num))
        allocate(raw%co2_airborne%longitude(raw%co2_airborne%num))
        allocate(raw%co2_airborne%height_msl(raw%co2_airborne%num))
        allocate(raw%co2_airborne%co2(raw%co2_airborne%num))
        allocate(raw%co2_airborne%ii(raw%co2_airborne%num))
        allocate(raw%co2_airborne%jj(raw%co2_airborne%num))
        allocate(raw%co2_airborne%kk(raw%co2_airborne%num))

        rewind(10)
        do n = 1, raw%co2_airborne%num
            read(10, '(f9.4, f10.4, f8.1, f9.4)', iostat=ierr) &
                lat, lon, height_msl, co2

            raw%co2_airborne%latitude(n) = lat
            raw%co2_airborne%longitude(n) = lon
            raw%co2_airborne%height_msl(n) = height_msl
            raw%co2_airborne%co2(n) = co2

            call latlon_to_ij(proj, lat, lon, aio, ajo)
            io = nint(aio)
            jo = nint(ajo)
            ako = height_to_k(ph(io,jo,:), height_msl)

            raw%co2_airborne%ii(n) = aio
            raw%co2_airborne%jj(n) = ajo
            raw%co2_airborne%kk(n) = ako
        end do
        close(10)

        return

    end subroutine get_co2_airborne_obs

    !--------------------------------------------------------------------------
    ! sort_co2_airborne_data
    !
    ! Prepare CO2 airborne data to be assimilated.
    !
    ! Parameters
    ! ----------
    ! ix:
    !   x (longitude) dimension size.
    ! jx:
    !   y (latitude) dimension size.
    ! datathin:
    !   Ignore every ``datathin``th airborne (NOT IMPLEMENTED).
    !   For example, if ``datathin`` is 2, then every second airborne is
    !   ignored.
    ! hroi:
    !   Horizontal radius of influence in gridpoints.
    ! vroi:
    !   Vertical radius of influence in gridpoints.
    !
    !--------------------------------------------------------------------------
    subroutine sort_co2_airborne_data(ix, jx, datathin, hroi, vroi)
        implicit none

        real, parameter                      :: HEIGHT_TOLERENCE = 100
        integer, intent(in)                  :: ix, jx
        integer, intent(in)                  :: datathin, hroi, vroi

        real                                 :: ii, jj, kk
        integer                              :: n

        do n = 1, raw%co2_airborne%num
            ii = raw%co2_airborne%ii(n)
            jj = raw%co2_airborne%jj(n)
            kk = raw%co2_airborne%kk(n)

            ! Check if the observation is within the domain
            if ((ii < 1. .or. ii > real(ix)) .or. &
                    (jj < 1. .or. jj > real(jx))) then
                if (my_proc_id == 0) then
                    write(*, *) 'Airborne obs is outside the domain, ', &
                                'skipping observation'
                end if

                cycle
            end if

            ! Check if the vertical position of the observation is reasonable
            if (kk == 0) then
                if (my_proc_id == 0) then
                    write(*, *) 'Airborne is above the highest model level', &
                                ' (how is that even possible?), ', &
                                ' skipping observation'
                end if

                cycle
            else if (kk < 0) then
                if (abs(kk) < HEIGHT_TOLERENCE) then
                    kk = 1
                else
                    if (my_proc_id == 0) then
                        write(*, *) 'Airborne is too far below the lowest ', &
                                    'model level, skipping observation'
                    end if

                    cycle
                end if
            end if

            if (raw%co2_airborne%co2(n) > 0) then
                obs%num                 = obs%num + 1
                obs%dat     (obs%num  ) = raw%co2_airborne%co2(n)
                obs%type    (obs%num  ) = 'co2a_   '
                obs%err     (obs%num  ) = co2_error_airborne
                obs%position(obs%num,1) = ii
                obs%position(obs%num,2) = jj
                obs%position(obs%num,3) = kk
                obs%roi     (obs%num,1) = hroi
                obs%roi     (obs%num,2) = vroi
            endif
        end do

        return

    end subroutine sort_co2_airborne_data

    !--------------------------------------------------------------------------
    ! xb_to_co2_airborne
    !
    ! Translate background state vector to CO2 concentration value at airborne
    ! locations.
    !
    ! The airborne value is obtained through linear interpolation horizontally
    ! and vertically.
    !
    ! Parameters
    ! ----------
    ! inputfile:
    !   Path to wrfout file with the background. Used only if 'CO2' is not
    !   a state variable.
    ! xb:
    !   State vector (background or prior).
    ! ix:
    !   x (longitude) dimension size.
    ! jx:
    !   y (latitude) dimension size.
    ! kx:
    !   z (vertical) dimension size.
    ! nv:
    !   Number of state variables.
    ! iob:
    !   Number of the current observation.
    !
    ! Returns
    ! -------
    ! hxb:
    !   h(xb) for the current observation.
    !
    !--------------------------------------------------------------------------
    subroutine xb_to_co2_airborne(inputfile, xb, ix, jx, kx, nv, iob, hxb)
        implicit none
        character(len=10), intent(in)           :: inputfile
        integer, intent(in)                     :: ix, jx, kx, nv, iob
        real, dimension(3,3,kx+1,nv), intent(in) :: xb
        real, intent(out)                       :: hxb

        real, dimension(ix,jx,kx)               :: co2
        real, dimension(kx)                     :: co2_vert
        integer                                 :: m
        integer                                 :: i_co2 = 0
        real                                    :: obs_ii, obs_jj, obs_kk
        integer                                 :: i1, j1, k1
        real                                    :: dx, dxm, dy, dym, dz, dzm

        obs_ii = obs%position(iob,1)
        obs_jj = obs%position(iob,2)
        obs_kk = obs%position(iob,3)

        i1 = int(obs_ii)
        j1 = int(obs_jj)
        k1 = int(obs_kk)

        dx = obs_ii - real(i1)
        dy = obs_jj - real(j1)
        dz = obs_kk - real(k1)
        dxm = real(i1+1) - obs_ii
        dym = real(j1+1) - obs_jj
        dzm = real(k1+1) - obs_kk

        do m = 1, nv
            if (trim(enkfvar(m)) == 'CO2') then
                i_co2 = m
            end if
        end do

        if (i_co2 > 0) then
            co2(i1:i1+1,j1:j1+1,1:kx) = xb(1:2,1:2,1:kx,i_co2)
        else
            call get_variable3d(inputfile, 'CO2', ix, jx, kx, 1, co2)
        end if

        co2_vert(:) = dym*(dx*co2(i1+1,j1,:) + dxm*co2(i1,j1,:)) + &
                      dy*(dx*co2(i1+1,j1+1,:) + dxm*co2(i1,j1+1,:))

        hxb = dzm*co2_vert(k1) + dz*co2_vert(k1+1)

        return

    end subroutine xb_to_co2_airborne

end module co2
