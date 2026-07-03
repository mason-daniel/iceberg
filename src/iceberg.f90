!   A program to read in a list of defect diameters, and offer possible lognormal distributions
!
!   Daniel Mason    
!   (c) UKAEA 2026
!


!   version history
!   vol.0.0.1     July 2026       First working version

    program iceberg
!---^^^^^^^^^^^^^^^
        use Lib_ColouredTerminal
        use Lib_CommandLineArguments
        use Lib_LogNormal
        use Lib_ChiSquare
        use Lib_LogisticFunction
        use Lib_DownhillSimplex
        use Lib_RandomSeed
        use Lib_Quicksort
        use iso_fortran_env
        implicit none

        character(len=8),parameter          ::      VERSION = "0.0.1"

    !---    command line args
        type(CommandLineArguments)          ::      cla
        real(kind=real64)                   ::      d0 = LIB_CLA_NODEFAULT_R            !   visibility characteristic size
        real(kind=real64)                   ::      w = -1                              !   visibility characteristic width: -1 = all visible
        character(len=256)                  ::      filename = ""                       !   input file name
        real(kind=real64)                   ::      rho = LIB_CLA_NODEFAULT_R          !   expected point defect density
        real(kind=real64)                   ::      omega0 = 1.0d0                      !   volume per atom
        real(kind=real64)                   ::      vol = LIB_CLA_NODEFAULT_R           !   observed volume
        logical                             ::      nDefectsObservedFixed = .true.             !   use observed defect count
        logical                             ::      nPointDefectsFixed = .false.               !   use rho to determine expected defect count
        logical                             ::      logisticFuncFixed = .false.         !   use a logistic func for invisibility
        logical                             ::      voids = .true.                      !  distribution is for voids not loops
        real(kind=real64)                   ::      b = 1                               !    burgers vector magnitude
    !---


    !---    input file
        integer                             ::      n                           !   number of data entries
        real(kind=real64)                   ::      x,xbar,x2bar
        real(kind=real64),dimension(10000)  ::      dat
        real(kind=real64)                   ::      observed_pointdefect_count,expected_pointdefect_count       
        real(kind=real64)                   ::      observed_defect_count,defect_count      
    !---

    !---    output
        real(kind=real64)                   ::      mu,sig
        type(LogNormal)                     ::      ln
        type(LogisticFunction)              ::      lf
        integer                             ::      NBINS = 20
        integer,dimension(:),allocatable        ::      observed    
        real(kind=real64),dimension(:),allocatable      ::      histogram_knot
        real(kind=real64),dimension(:),allocatable      ::      histogram_knot_fine
        real(kind=real64),dimension(:),allocatable      ::      expected,expected_nolf
    !---

    !---    fit
        type(DownhillSimplex)               ::      simp
        real(kind=real64),dimension(:),allocatable      ::      ss      !  a simplex point
        real(kind=real64),dimension(:),allocatable      ::      bestss
        integer                             ::      nTrials = 100
    !---    

    !---    dummy variable
        integer                             ::      ii,jj,trial
        logical                             ::      ok
    !---



    


    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)
        call setProgramDescription( cla, "iceberg"//new_line("A")//" A program to read in a list of defect diameters, and offer possible lognormal distributions." )
        call setProgramVersion( cla, VERSION )
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"         input filename")          
        call get( cla,"d0",d0 ,LIB_CLA_OPTIONAL,"      visibility characteristic size")          
        call get( cla,"w",w ,LIB_CLA_OPTIONAL,"       visibility characteristic width - set -1 for all visible")          
        !if (hasArgument(cla,"d0") .neqv. hasArgument(cla,"w")) stop "iceberg error - must set both -d0 and -w or neither"
        logisticFuncFixed = (hasArgument(cla,"d0") .and. hasArgument(cla,"w")) .or. (hasArgument(cla,"w") .and. (w==-1))


        call get( cla,"rho",rho ,LIB_CLA_OPTIONAL,"    expected point defect density (at fr)")          
        call get( cla,"omega0",omega0 ,LIB_CLA_OPTIONAL,"  volume per atom")          
        call get( cla,"vol",vol ,LIB_CLA_OPTIONAL,"       observed volume")          
        if (hasArgument(cla,"vol")) then
            if (.not. hasArgument(cla,"omega0")) stop "iceberg error - can't use observed -vol without volume per atom -omega0"
            
        else 
            if (hasArgument(cla,"rho")) stop "iceberg error - can't use desired -rho without volume -vol"
        end if 
        call get( cla,"voids",voids ,LIB_CLA_OPTIONAL,"     distribution is for voids not loops")          
        call get( cla,"b",b ,LIB_CLA_OPTIONAL,"         Burgers vector magnitude ( for loop point defect count )")          
        if ((.not. voids).and.(.not. hasArgument(cla,"voids"))) stop "iceberg error - can't do a loops estimate without -b"


        nPointDefectsFixed = (hasArgument(cla,"rho")) !   fixed by expected point defect count

        nDefectsObservedFixed = (.not. nPointDefectsFixed) .and. (logisticFuncFixed.and.(w==-1))
        

        call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop "iceberg error - required arguments unset"
        call delete(cla)
    !---
 

    !---
        print *,""
        print *,colour(LIGHT_AQUA,"iceberg")
        print *,""
        print *,"   use observed defect count? ",nDefectsObservedFixed
        print *,"   number of pd fixed?        ",nPointDefectsFixed
        print *,"   logistic func fixed?       ",logisticFuncFixed
        print *,""
    !---




    !---    read input file
        print *,"iceberg info - reading from """//trim(filename)//""""
        inquire(file=trim(filename),exist=ok)
        if (.not. ok) stop "iceberg error - could not find input file"
        open(unit=500,file=trim(filename),action="read")
            do n = 1,size(dat)
                read(unit=500,fmt=*,iostat=ii) dat(n)
                if (ii/=0) exit     !   failed to parse a number
            end do        
            n = n - 1   !   because last read failed            
        close(unit=500)
    !---    

        NBINS = ceiling( sqrt(real(n)) )
        allocate(observed(0:NBINS-1))
        allocate(histogram_knot(0:NBINS))
        allocate(histogram_knot_fine(0:NBINS*NBINS))
        allocate(expected(0:NBINS-1))
        allocate(expected_nolf(0:NBINS-1))

    !---

    
    !---    report success of read
        print *,"iceberg info - read ",n," lines"
        call quickSort(dat(1:n))
        observed_defect_count = n
        xbar = 0 ; x2bar = 0
        do ii = 1,n
            x = dat(ii)
            xbar = xbar + x
            x2bar = x2bar + x*x
        end do
        if (n<1) stop "iceberg error - zero lines read"
        xbar = xbar/n
        x2bar = x2bar/n
        if (n<10) stop "iceberg warning - very few lines read - expect to fail"
        sig = sqrt( max(0.0d0,x2bar - xbar*xbar) )      !   true stdev, not lognormal shape function
        print *,"iceberg info - <d>,stdev(d) = ",xbar,sig
        print *,""
        print *,""
    !---

    !---    fit data - first pass
        print *,"Naive fitting - lognormal with mean and stdev matching input data"
        x = log( x2bar/(xbar*xbar) )        !    = sig^2, lognormal shape func
        mu = xbar / exp( x/2 )              !   lognormal scale
        sig = sqrt( x )                     !   lognormal shape func
        ln = LogNormal_ctor(mu,sig)        
        call report(ln,o=1)
        print *,"observed defect count        ",observed_defect_count
        print *,"68% confidence interval  <d> ",confidenceLevel(ln,0.16d0),":",moment(ln,1),":",confidenceLevel(ln,0.84d0)
        print *,""
    !---    



    !---    construct the logistic function
        print *,"construct initial logistic function"
        if (.not. logisticFuncFixed) then            
            d0 = confidenceLevel(ln,0.1d0)
            w = confidenceLevel(ln,0.2d0) - d0
        end if
        lf = LogisticFunction_ctor(d0,w)
        call report(lf,o=1)        
        print *,""
    !---



    !---    construct a histogram to bin the results. Put about 10% of the results into each bin.
        do jj = 0,NBINS
            x = jj*(1.0d0-1.0d-8)/NBINS      !   from 0 to <1
            histogram_knot(jj) = confidenceLevel(ln,x)
        end do

        do jj = 0,NBINS*NBINS
            x = jj*(1.0d0-1.0d-8)/(NBINS*NBINS)      !   from 0 to 1
            histogram_knot_fine(jj) = confidenceLevel(ln,x)
        end do

        observed = 0
        do ii = 1,n
            x = dat(ii)
            do jj = 0,NBINS-1
                if (x < histogram_knot(jj+1)) then
                    observed(jj) = observed(jj) + 1
                    exit
                end if
            end do
        end do
        call expectedHistogram( ln,lf,real(n,kind=real64),histogram_knot,histogram_knot_fine,expected )
        print *,"observed results histogram"
        write(*,fmt='(a16,a6,a16)') "min diam","count","expected"
        do jj = 0,NBINS-1
            write (*,fmt='(g16.8,i6,g16.8)') histogram_knot(jj),observed(jj),expected(jj)
        end do
        print *,"chi-square value ",chi_square(observed,expected)
        print *,""
    !---
     

    !---    find point defect density if requested
        observed_pointdefect_count = 0
        !if (vol /= LIB_CLA_NODEFAULT_R) then
            xbar = 0 ; x2bar = 0
            do ii = 1,n
                x = pointDefectCount( dat(ii),omega0 )
                xbar = xbar + x
                x2bar = x2bar + x*x
            end do
            xbar = xbar/n
            x2bar = x2bar/n
            print *,"point defects observed per defect <n>,stdev ",xbar,sqrt(max(0.0d0,x2bar-xbar*xbar))
            if (vol /= LIB_CLA_NODEFAULT_R) &
            print *,"point defect frac ",xbar/vol,sqrt(max(0.0d0,x2bar-xbar*xbar))/vol
            observed_pointdefect_count = xbar * n
            if (rho == LIB_CLA_NODEFAULT_R) then
                expected_pointdefect_count = observed_pointdefect_count 
            else 
                expected_pointdefect_count = rho * vol
                print *,"expected_pointdefect_count ",rho * vol
            end if
            print *,""

        !end if
    !---


    !---    fit
        print *,"Fitting to data"
        call constructInitialSimplex( ln,lf,simp )
        select case(getD(simp))
            case (2)        !   mu,sig
                print *,"fitting lognormal mu,sig"
            case (3)        !   mu,sig,n
                print *,"fitting lognormal mu,sig + defect count"                
            case (4)        !   mu,sig,d0,w
                print *,"fitting lognormal mu,sig + logistic function d0,w"                
            case (5)        !   mu,sig,n,d0,w
                print *,"fitting lognormal mu,sig + defect count + logistic function d0,w"                
        end select
        call report(simp,o=4)
        allocate(ss(getD(simp)))
        allocate(bestss(0:getD(simp)))       
        bestss(0) = huge(1.0) 
        do trial = 1,nTrials
            do ii = 1,10000
                call suggestNextPoint( simp,ss )
                x = getChiSquared( ss )
                call addSimplexPoint( simp,ss,x )
            !    if (mod(ii,1000)==0) then
            !        write(*,fmt='(a,i6,a,100g16.8)',advance="no") "step ",ii," point ",ss
            !        write(*,fmt='(a,g16.8,a)',advance="no") " chisq ",x,"   "
            !        call report(simp)
            !    end if
                if (isConverged( simp, 1d-8,1d-8 )) exit
            end do
        !    call report(simp)
            call getBestSimplexPoint(simp,ss,x)
            if (x < bestss(0)) then 
                bestss(1:) = ss
                bestss(0) = x
                call report(simp,o=4)
            end if
            
            if (trial < nTrials) then
                call delete(simp)
                call constructInitialSimplex( ln,lf,simp )
            end if
        end do
        ss = bestss(1:)
        !print *,bestss
        print *,""
    !---

    !---    report
        select case(getD(simp))
            case (2)        !   mu,sig
                mu = ss(1) 
                sig = ss(2)                 
                defect_count = expectedDefectCount( LogNormal_ctor(mu,sig),omega0,expected_pointdefect_count )
            case (3)        !   mu,sig,n
                mu = ss(1) 
                sig = ss(2)                 
                defect_count = ss(3)                        
            case (4)        !   mu,sig,d0,w
                mu = ss(1) 
                sig = ss(2)               
                d0 = ss(3)
                w = ss(4)  
                defect_count = expectedDefectCount( LogNormal_ctor(mu,sig),omega0,expected_pointdefect_count )
            case (5)        !   mu,sig,n,d0,w
                mu = ss(1) 
                sig = ss(2)                 
                defect_count = ss(3)        
                d0 = ss(4)
                w = ss(5)  
        end select
        !print *,"mu,sig,defect_count,d0,w = ",mu,sig,defect_count,d0,w
        ln = LogNormal_ctor(mu,sig)        
        print *,"Fitting - lognormal"
        call report(ln,o=1)
        if ( .not. isUnset(lf)) then
            if (logisticFuncFixed) then
                print *,"using fixed visibility function"
                call report(lf,o=1)
            else
                print *,"fitted visibility function"
                lf = LogisticFunction_ctor(d0,w)
                call report(lf,o=1)
            end if
        else
            print *,"no visibility function"
        end if
        call expectedHistogram( ln,lf,defect_count,histogram_knot,histogram_knot_fine,expected )
        print *,"expected results histogram"
        write(*,fmt='(a16,a6,a16)') "min diam","count","expected"
        do jj = 0,NBINS-1
            write (*,fmt='(g16.8,i6,g16.8)') histogram_knot(jj),observed(jj),expected(jj)
        end do
        print *,"chi-square value            ",chi_square(observed,expected)
        print *,"alpha value                 ",alpha_score( chi_square(observed,expected), df=nBins-getD(simp) )
        print *,"observed defect count       ",observed_defect_count
        print *,"expected defect count       ",defect_count
        print *,"observed pointdefect count  ",observed_pointdefect_count
        print *,"expected point defect count ",expectedPointDefectCount( ln,omega0,defect_count )
        print *,"68% confidence interval <d> ",confidenceLevel(ln,0.16d0),":",moment(ln,1),":",confidenceLevel(ln,0.84d0)
        
        print *,""



        print *,"even bin histogram for output"
       ! ii = int( observed_defect_count*0.99 )      !   99% percentile
        d0 = roundedBin( dat(n) )        
        do jj = 0,NBINS
            x = jj*d0/NBINS
            histogram_knot(jj) = x
        end do
        do jj = 0,NBINS*NBINS
            x = jj*d0/(NBINS*NBINS)
            histogram_knot_fine(jj) = x
        end do
        observed = 0
        do ii = 1,n
            x = dat(ii)
            do jj = 0,NBINS-1
                if (x < histogram_knot(jj+1)) then
                    observed(jj) = observed(jj) + 1
                    exit
                end if
            end do
        end do
        call expectedHistogram( ln,lf,defect_count,histogram_knot,histogram_knot_fine,expected )
        if (isUnset(lf)) then
            write(*,fmt='(a16,a6,2a16)') "mean diam","count","expected"
            do jj = 0,NBINS-1
                write (*,fmt='(g16.8,i6,2g16.8)') (histogram_knot(jj)+histogram_knot(jj+1))/2,observed(jj),expected(jj)
            end do
        else
            call expectedHistogram( ln,LogisticFunction_ctor(),defect_count,histogram_knot,histogram_knot_fine,expected_nolf )    
            write(*,fmt='(a16,a6,2a16)') "mean diam","count","expected","if all vis"
            do jj = 0,NBINS-1
                write (*,fmt='(g16.8,i6,2g16.8)') (histogram_knot(jj)+histogram_knot(jj+1))/2,observed(jj),expected(jj),expected_nolf(jj)
            end do
        end if

    !---    bye bye
        print *,""
        print *,"done"
        print *,""
        


    contains
!---^^^^^^^^

        real(kind=real64) function roundedBin(dmax)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      choose a nice looking end point bin
            real(kind=real64),intent(in)        ::  dmax
            real(kind=real64),dimension(31)  ::  binends = (/ 0.1,0.2,0.4,0.5,0.6,0.8,   &
                                                             1.0,2.0,4.0,5.0,6.0,8.0,   &
                                                             10.,20.,40.,50.,60.,80.,   &
                                                             1e2,2e2,4e2,5e2,6e2,8e2,   &
                                                             1e3,2e3,4e3,5e3,6e3,8e3,   &
                                                             1e4 /)
            integer ::  ii
            if (dmax < binends(1)) then
                roundedBin = dmax
            else
                do ii = 2,size(binends)
                    if (dmax < binends(ii)) then
                        roundedBin = binends(ii)
                        return
                    end if
                end do
                roundedBin = dmax
            end if
            return
        end function roundedBin


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


        subroutine expectedHistogram( logn,logf,defect_count,histogram_knot,histogram_knot_fine,expected )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given an input lognormal probability distribution, 
    !*      and an expected logistic function
    !*      and a total void ount,
    !*      and a set of knot points,
    !*      return the expected count in each bin
            type(LogNormal),intent(in)          ::      logn
            type(LogisticFunction),intent(in)   ::      logf
            real(kind=real64),intent(in)        ::      defect_count
            real(kind=real64),dimension(0:),intent(in)  ::      histogram_knot,histogram_knot_fine
            real(kind=real64),dimension(0:),intent(out) ::      expected
            integer             ::      nBins,jj,kk
            real(kind=real64)   ::      c_old,cc,xx,gg,x_old,intgrl

            nBins = size(histogram_knot) - 1
            c_old = 0.0d0


            if (isUnset(logf)) then

                do jj = 0,nBins-1
    
                    xx = histogram_knot(jj+1)
                    cc = cdf( logn,xx )                !   cdf is integral
                    expected(jj) = (cc - c_old) 
                    c_old = cc

                end do

            else
                x_old = 0.0d0
                do jj = 0,nBins-1

                    intgrl = 0.0d0
                    do kk = jj*nBins,jj*nBins+nBins-1
                        xx = histogram_knot_fine(kk+1)
                        gg = func( logf,(x_old + xx)/2 )        !   compute logistic func at midpoint
                        cc = cdf( logn,xx )                     !   cdf is integral
                        intgrl = intgrl + gg * (cc - c_old) 
                        c_old = cc
                        x_old = xx
                    end do
                    expected(jj) = intgrl

                end do
            end if


            expected = expected * defect_count
            return
        end subroutine expectedHistogram
 


        real(kind=real64) function expectedDefectCount( logn,omega0,npd )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given we expect to have npd point defects, return the expected defect count
    !*          npd = Nvoid (pi/6) (1/omega0) <d³>
            type(LogNormal),intent(in)              ::      logn
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            real(kind=real64),intent(in)            ::      npd     !  expected pd count
            real(kind=real64),parameter             ::  PI = 3.141592654d0
            real(kind=real64)           ::      xx

            expectedDefectCount = observed_defect_count
            if (nPointDefectsFixed .and. .not. nDefectsObservedFixed) then
                if (voids) then
                    xx = moment(logn,3)
                    expectedDefectCount = npd * 6 * omega0 / ( PI * xx )
                else
                    xx = moment(logn,2)
                    expectedDefectCount = npd * 4 * omega0 / ( PI * b * xx )
                end if
            end if

            return
        end function expectedDefectCount


        real(kind=real64) function expectedPointDefectCount( logn,omega0,defect_count )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogNormal),intent(in)              ::      logn
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            real(kind=real64),intent(in)            ::      defect_count     !  expected pd count
            real(kind=real64),parameter             ::  PI = 3.141592654d0
            real(kind=real64)           ::      xx

            if (voids) then
                xx = moment(logn,3)
                expectedPointDefectCount = ( PI * xx )/( 6 * omega0 )
            else
                xx = moment(logn,2)
                expectedPointDefectCount = ( PI * b * xx )/( 4 * omega0 )
            end if
            expectedPointDefectCount = expectedPointDefectCount * defect_count

            return
        end function expectedPointDefectCount        


        real(kind=real64) function getChiSquared( x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      get the chi-square value at search point x
            real(kind=real64),dimension(:),intent(in)       ::      x
            integer     ::      dd
            type(Lognormal)         ::  logn
            type(LogisticFunction)  ::  logf
            real(kind=real64),dimension(0:NBINS-1)  ::      expected
            real(kind=real64)       ::  defect_count
            real(kind=real64),dimension(size(x))    ::      sanitisedx

            dd = size(x)
            sanitisedx = sanitise(x)
            select case(dd)
                case (2)        !   mu,sig
                    logn = Lognormal_ctor(sanitisedx(1),sanitisedx(2))
                    logf = LogisticFunction_ctor()
                    defect_count = expectedDefectCount( logn,omega0,expected_pointdefect_count )
                case (3)        !   mu,sig,n
                    logn = Lognormal_ctor(sanitisedx(1),sanitisedx(2))
                    logf = LogisticFunction_ctor()
                    defect_count = sanitisedx(3)
                case (4)        !   mu,sig,d0,w
                    logn = Lognormal_ctor(sanitisedx(1),sanitisedx(2))
                    logf = LogisticFunction_ctor(sanitisedx(3),sanitisedx(4))
                    defect_count = expectedDefectCount( logn,omega0,expected_pointdefect_count )
                case (5)        !   mu,sig,n,d0,w
                    logn = Lognormal_ctor(sanitisedx(1),sanitisedx(2))
                    defect_count = sanitisedx(3)
                    logf = LogisticFunction_ctor(sanitisedx(4),sanitisedx(5))
            end select

            call expectedHistogram( logn,logf,defect_count,histogram_knot,histogram_knot_fine,expected )
            getChiSquared = chi_square(observed,expected)

            return


        end function getChiSquared

        elemental real(kind=real64) function sanitise(x)
            real(kind=real64),intent(in)        ::      x
            sanitise = max( 1.0d-8,min(1d8,x) )
        end function sanitise



        subroutine constructInitialSimplex( logn,logf,simp )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(Lognormal),intent(in)          ::      logn
            type(LogisticFunction),intent(in)   ::      logf
            type(DownhillSimplex),intent(out)   ::      simp

            integer         ::      dd,nn       !   search space dimensionality, number of simplex points
            real(kind=real64),dimension(:),allocatable      ::      xx      !   simplex point
            real(kind=real64)                               ::      cc      !   chi-square value at this point
            integer         ::      ii


            if (logisticFuncFixed) then
                if (nPointDefectsFixed .or. nDefectsObservedFixed) then
                    !   only adjusting lognormal mu,sig
                    dd = 2
                else 
                    !   adjusting mu,sig, and defect count 
                    dd = 3
                end if
            else 
                if (nPointDefectsFixed .or. nDefectsObservedFixed) then
                    !   adjusting lognormal mu,sig + logistic d0,w
                    dd = 4
                else 
                    !   adjusting mu,sig, and defect count + logistic d0,w
                    dd = 5
                end if
            end if

            nn = dd * 2
            simp = DownhillSimplex_ctor(nn,dd)
            allocate(xx(dd))

            do ii = 1,nn
                xx(1:dd) = 1 + gaussianVariate(dd)/4
                xx = max(1.0d-3,xx)
                select case(dd)
                    case (2)        !   mu,sig
                        xx(1) = getMu(logn)*xx(1)
                        xx(2) = getSigma(logn)*xx(2)
                    case (3)        !   mu,sig,n
                        xx(1) = getMu(logn)*xx(1)
                        xx(2) = getSigma(logn)*xx(2)
                        xx(3) = expectedDefectCount( logn,omega0,expected_pointdefect_count )*xx(3)
                    case (4)        !   mu,sig,d0,w
                        xx(1) = getMu(logn)*xx(1) 
                        xx(2) = getSigma(logn)*xx(2) 
                        xx(3) = getx0(logf)*xx(3) 
                        xx(4) = getw(logf)*xx(4) 
                    case (5)        !   mu,sig,n,d0,w
                        xx(1) = getMu(logn)*xx(1)
                        xx(2) = getSigma(logn)*xx(2)
                        xx(3) = expectedDefectCount( logn,omega0,expected_pointdefect_count )*xx(3)
                        xx(4) = getx0(logf)*xx(4)
                        xx(5) = getw(logf)*xx(5)
                end select

                cc = getChiSquared( xx )
                !print *,"xx = ",xx," c = ",cc
                call addSimplexPoint(simp,ii,xx,cc)

            end do
       

            return
        end subroutine constructInitialSimplex
    
        


    end program iceberg