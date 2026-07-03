
    program test_genData
!---^^^^^^^^^^^^^^^^^^^^
!*      generates a test file with known mu,sig
!*      Daniel Mason
!*      (c) UKAEA July 2026

        use Lib_Lognormal
        use Lib_LogisticFunction
        use Lib_CommandLineArguments
        use iso_fortran_env
        implicit none


    !---    command line args
        type(CommandLineArguments)          ::      cla
        real(kind=real64)                   ::      mu = 4.0            !   mean
        real(kind=real64)                   ::      sigma = 2.0         !   stdev
        real(kind=real64)                   ::      d0 = LIB_CLA_NODEFAULT_R            !   visibility characteristic size
        real(kind=real64)                   ::      w = -1                              !   visibility characteristic width: -1 = all visible
        integer                             ::      n = 1000            !   count
        character(len=256)                  ::      filename = ""       !   output file name
        logical                             ::      voids = .true.
        real(kind=real64)                   ::      omega0 = LIB_CLA_NODEFAULT_R                      !   volume per atom
    !---

        type(Lognormal)                     ::      ln
        type(LogisticFunction)              ::      lf
        integer                             ::      ii,nn
        real(kind=real64)                   ::      d,dbar,d2bar,zeta,npd,npd_obs
        real(kind=real64)                   ::      mm,ss



    !---    
        cla = CommandLineArguments_ctor(30)
        call setProgramDescription( cla, "test_genData"//new_line("A")//" generate data with known lognormal distribution." )
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"         output filename")          
        call get( cla,"mu",mu ,LIB_CLA_OPTIONAL,"      mean")          
        call get( cla,"sigma",sigma ,LIB_CLA_OPTIONAL,"   stdev")          
        call get( cla,"n",n ,LIB_CLA_OPTIONAL,"       count")          
        call get( cla,"d0",d0 ,LIB_CLA_OPTIONAL,"      visibility characteristic size")          
        call get( cla,"w",w ,LIB_CLA_OPTIONAL,"       visibility characteristic width - set -1 for all visible")          
        call get( cla,"voids",voids ,LIB_CLA_OPTIONAL,"     distribution is for voids not loops")          
        call get( cla,"omega0",omega0 ,LIB_CLA_OPTIONAL,"  volume per atom")  

        call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop "fitLogNormal error - required arguments unset"
        call delete(cla)
    !---        

    

    !---    write file
        d = sigma*sigma + mu*mu             !   mean square
        d = log( d/(mu*mu) )                !    = sig^2, lognormal shape func
        mm = mu / exp( d/2 )                !   lognormal scale
        ss = sqrt( d )                      !   lognormal shape func
        ln = Lognormal_ctor(mm,ss)
        print *,"generation function"
        call report(ln)        
        print *,"68% confidence interval   ",confidenceLevel(ln,0.16d0),":",moment(ln,1),":",confidenceLevel(ln,0.84d0)
        print *,"<mean>   : ",moment(ln,1)
        print *,"<stdev>  : ",sqrt( moment(ln,2) - moment(ln,1)**2 )

        print *,"visibility function"
        lf = LogisticFunction_ctor(d0,w)
        call report(lf)        
        print *,""


        open(unit=600,file=trim(filename),action="write")
            npd = 0
            if (isUnset(lf)) then
                do ii = 1,n
                    d = variate(ln)
                    write(unit=600,fmt=*) d
                    dbar = dbar + d
                    d2bar = d2bar + d*d
                    npd = npd + pointDefectCount(d,omega0)
                end do
                npd_obs = npd
            else
                nn = 0
                do ii = 1,n
                    do
                        nn = nn + 1 
                        d = variate(ln)
                        npd = npd + pointDefectCount(d,omega0)
                        call random_number(zeta)
                        if (zeta < func(lf,d)) then
                            write(unit=600,fmt=*) d
                            dbar = dbar + d
                            d2bar = d2bar + d*d
                            npd_obs = npd_obs + pointDefectCount(d,omega0)
                            exit
                        end if
                    end do
                end do
            end if
        close(unit=600)
        dbar = dbar/n
        d2bar = d2bar/n
        print *,"n real   : ",nn
        print *,"n seen   : ",n
        if (omega0 /= LIB_CLA_NODEFAULT_R) &
        print *,"npd real : ",npd
        print *,"npd obs  : ",npd_obs
        print *,"<d>      : ",dbar
        print *,"<d²>     : ",d2bar
        print *,"stdev    : ",sqrt(d2bar  - dbar*dbar) 
    !---



        print *,""
        print *,"pass"
        print *,""


    contains
!---^^^^^^^^

    
        pure real(kind=real64) function pointDefectCount(d,omega0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return count of point defects given vol per atom 
            real(kind=real64),intent(in)            ::      d   !   diam
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            real(kind=real64),parameter     ::  PI = 3.141592654d0
            
            if (voids) then
                pointDefectCount = (PI/6) * d**3 / omega0
            else
                pointDefectCount = (PI/4) * d**2 / (4 * omega0)
            end if

            return
        end function pointDefectCount

    end program test_genData
