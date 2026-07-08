!   A program to read in a list of defect diameters, and offer possible lognormal distributions
!
!   Daniel Mason    
!   (c) UKAEA 2026
!


!   version history
!   v.0.0.1     July 2026       First working version
!   v.0.0.2     July 2026       Option ( on by default ) to eliminate tiny (zero point defect) diameters from analysis

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
        real(kind=real64),parameter         ::      PI = 3.141592653590d0

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
        integer                             ::      NBINS = LIB_CLA_NODEFAULT_I
        real(kind=real64)                   ::      maxBin = LIB_CLA_NODEFAULT_R
        integer                             ::      nTrials = 500
        logical                             ::      zero = .true.                       !   zero size defect correction - expected count for pd [0:0.5] = 0
    !---


    !---    input file
        integer                             ::      n                           !   number of data entries
        real(kind=real64)                   ::      x,xbar,x2bar
        real(kind=real64),dimension(10000)  ::      dat
        real(kind=real64)                   ::      observed_pointdefect_count,expected_pointdefect_count               
        real(kind=real64)                   ::      observed_defect_count,defect_count      
        real(kind=real64)                   ::      observed_defect_mean,observed_defect_stdev    
    !---

    !---    output
        real(kind=real64)                   ::      mu,sig
        type(LogNormal)                     ::      ln
        type(LogisticFunction)              ::      lf
        
        integer,dimension(:),allocatable        ::      observed    
        real(kind=real64),dimension(:),allocatable      ::      histogram_knot
        real(kind=real64),dimension(:),allocatable      ::      histogram_knot_fine
        real(kind=real64),dimension(:),allocatable      ::      expected,expected_nolf
        integer                             ::      oversized
        real(kind=real64)                   ::      alpha16,alpha84
    !---

    !---    fit
        real(kind=real64)                   ::      dzero       !   diameter which gives npd = 0.5
        type(DownhillSimplex)               ::      simp
        real(kind=real64),dimension(:),allocatable      ::      ss      !  a simplex point
        real(kind=real64),dimension(:),allocatable      ::      bestss
    !---    

    !---    dummy variable
        integer                             ::      ii,jj,trial
        logical                             ::      ok
    !---



    


    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)
        call setProgramDescription( cla, "iceberg"//new_line("A")//" A program to read in a list of defect diameters, and offer possible lognormal distributions.",LIGHT_AQUA )
        call setProgramVersion( cla, VERSION )
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"         input filename")          
        call get( cla,"d0",d0 ,LIB_CLA_OPTIONAL,"      visibility characteristic size")          
        call get( cla,"w",w ,LIB_CLA_OPTIONAL,"       visibility characteristic width - set -1 for all visible")          
        !if (hasArgument(cla,"d0") .neqv. hasArgument(cla,"w")) stop "iceberg error - must set both -d0 and -w or neither"
        logisticFuncFixed = (hasArgument(cla,"d0") .and. hasArgument(cla,"w")) .or. (hasArgument(cla,"w") .and. (w==-1))


        call get( cla,"rho",rho ,LIB_CLA_OPTIONAL,"     expected point defect density (at fr)")          
        call get( cla,"omega0",omega0 ,LIB_CLA_OPTIONAL,"  volume per atom")          
        call get( cla,"vol",vol ,LIB_CLA_OPTIONAL,"     observed volume")          
        if (hasArgument(cla,"vol")) then
            if (.not. hasArgument(cla,"omega0")) stop "iceberg error - can't use observed -vol without volume per atom -omega0"
            
        else 
            if (hasArgument(cla,"rho")) stop "iceberg error - can't use desired -rho without volume -vol"
        end if 
        call get( cla,"voids",voids ,LIB_CLA_OPTIONAL,"   distribution is for voids not loops")          
        call get( cla,"b",b ,LIB_CLA_OPTIONAL,"       Burgers vector magnitude ( for loop point defect count )")          
        if ((.not. voids).and.(.not. hasArgument(cla,"voids"))) stop "iceberg error - can't do a loops estimate without -b"

        call get( cla,"n",nBins ,LIB_CLA_OPTIONAL,"       number of histogram bins")          
        call get( cla,"dmax",maxBin ,LIB_CLA_OPTIONAL,"    max histogram diameter")          
        call get( cla,"zero",zero ,LIB_CLA_OPTIONAL,"    zero size defect correction")          
        call get( cla,"ntrials",ntrials ,LIB_CLA_OPTIONAL,"  number of random start points to use")



        nPointDefectsFixed = (hasArgument(cla,"rho")) !   fixed by expected point defect count

        nDefectsObservedFixed = (.not. nPointDefectsFixed) .and. (logisticFuncFixed.and.(w==-1))
        

        call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop "iceberg error - required arguments unset"
        call delete(cla)
    !---
 

    !---
        print *,""
!        print *,colour(LIGHT_AQUA,"iceberg")
        print *,""
        print *,"   void diameter-volume?      ",voids
        print *,"   use observed defect count? ",nDefectsObservedFixed
        print *,"   number of pd fixed?        ",nPointDefectsFixed
        print *,"   logistic func fixed?       ",logisticFuncFixed
        print *,""
    !---


!------------------------------------------------------------------------------
!
!       read input file
!
!------------------------------------------------------------------------------
      
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
        if ( NBINS == LIB_CLA_NODEFAULT_I )  NBINS = ceiling( sqrt(real(n)) )
        allocate(observed(0:NBINS-1))
        allocate(histogram_knot(0:NBINS))
        allocate(histogram_knot_fine(0:NBINS*NBINS))
        allocate(expected(0:NBINS-1))
        allocate(expected_nolf(0:NBINS-1))
        print *,""
        print *,""
    !---



!------------------------------------------------------------------------------
!
!       report basic stats about read to prove correct input
!
!------------------------------------------------------------------------------
        print *,"iceberg info - basic data: read ",n," lines"
        dzero = 0.0d0
        if (zero) then
            dzero = defectDiameter( 1.0d0,omega0 )
            print *,"iceberg info - zero size defect correction- expected count zero for d < ",dzero," npd = 1"
        end if

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
        observed_defect_mean = xbar
        observed_defect_stdev = sqrt( max(0.0d0,x2bar - xbar*xbar) )      !   true stdev, not lognormal shape function
        print *,"iceberg info - <d>,stdev(d) = ",observed_defect_mean,observed_defect_stdev
        print *,"iceberg info - naive fitting lognormal with mean and stdev matching input data"
        x = log( x2bar/(xbar*xbar) )        !    = sig^2, lognormal shape func
        mu = xbar / exp( x/2 )              !   lognormal scale
        sig = sqrt( x )                     !   lognormal shape func
        ln = LogNormal_ctor(mu,sig)        
        call report(ln,o=1)
    !---    construct the logistic function
        print *,"iceberg info - initial logistic function"
        if (.not. logisticFuncFixed) then            
            d0 = confidenceLevel(ln,0.1d0)
            w = confidenceLevel(ln,0.2d0) - d0
        end if
        lf = LogisticFunction_ctor(d0,w)
        call report(lf,o=1)        
        print *,""
        print *,""
    !---


!------------------------------------------------------------------------------
!
!       make an estimate of the observed point defect density
!
!------------------------------------------------------------------------------                
        observed_pointdefect_count = 0
        xbar = 0 ; x2bar = 0
        do ii = 1,n
            x = pointDefectCount( dat(ii),omega0 )
            xbar = xbar + x
            x2bar = x2bar + x*x
        end do
        xbar = xbar/n
        x2bar = x2bar/n
        print *,"iceberg info - point defects observed per defect <n>,stdev ",xbar,sqrt(max(0.0d0,x2bar-xbar*xbar))
        observed_pointdefect_count = xbar * n
        if (vol /= LIB_CLA_NODEFAULT_R) &
        print *,"iceberg info - observed point defect volume fraction ",observed_pointdefect_count*omega0/vol 
        
        if (rho == LIB_CLA_NODEFAULT_R) then
            expected_pointdefect_count = observed_pointdefect_count 
!            print *,"iceberg info - expected point defect volume fraction ",expected_pointdefect_count            
        else 
            expected_pointdefect_count = rho * vol / omega0
            print *,"iceberg info - expected point defect volume fraction ",rho
            print *,"iceberg info - expected point defect count ",expected_pointdefect_count
        end if
        print *,""
        print *,""




!------------------------------------------------------------------------------
!
!       construct a histogram to bin the results, such that 1/NBINS of the observed results are in each bin
!
!------------------------------------------------------------------------------       
        print *,"iceberg info - constructing histogram before fitting with approximately even filling of the observed results" 
        do jj = 0,NBINS
            x = jj*(1.0d0-1.0d-8)/NBINS             !   from 0 to <1
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

        call expectedHistogram( ln,lf,histogram_knot,histogram_knot_fine,expected )
        if (nPointDefectsFixed) then
            defect_count = expectedDefectCount( ln,omega0,expected_pointdefect_count )  
        else
            defect_count = observed_defect_count / sum(expected)
        end if
        expected = expected * defect_count

        print *,"iceberg info - histogram before fitting"
        write(*,fmt='(a16,a6,a16)') "min diam","count","expected"
        do jj = 0,NBINS-1
            write (*,fmt='(g16.8,i6,g16.8)') histogram_knot(jj),observed(jj),expected(jj)
        end do
        print *,"chi-square value            ",chi_square(observed,expected)
        print *,"alpha value                 ",alpha_score( chi_square(observed,expected), df=nBins-2 )
        print *,"observed defect count       ",observed_defect_count
        print *,"unfitted defect count       ",defect_count
        print *,"observed point defect count ",observed_pointdefect_count
        print *,"unfitted point defect count ",expected_pointdefect_count
        print *,"observed <d>,+/-,stdev(d)   ",observed_defect_mean,observed_defect_stdev/sqrt(observed_defect_count-1),observed_defect_stdev
        x = moment(ln,2,dzero) - moment(ln,1,dzero)**2        !   variance
        x = sqrt( defect_count*x/(defect_count-1) )         !   unbiassed stdev  
        print *,"unfitted <d>,+/-            ",moment(ln,1,dzero),x/sqrt(defect_count)
        print *,"68% confidence unfitted d   ",confidenceLevel(ln,0.16d0),":",confidenceLevel(ln,0.84d0)
        print *,""
        print *,""
    !---
     



    !---


!------------------------------------------------------------------------------
!
!       fit the histogram using downhill simplex method
!
!------------------------------------------------------------------------------       
        print *,"iceberg info - Fitting to data"
        call constructInitialSimplex( ln,lf,simp )
        select case(getD(simp))
            case (2)        !   mu,sig
                print *,"fitting lognormal mu,sig"
            case (4)        !   mu,sig,d0,w
                print *,"fitting lognormal mu,sig + logistic function d0,w"                
        end select
        write (*,fmt='(a,i8)',advance="no") "   trial: ",0
        call report(simp,o=2)
        allocate(ss(getD(simp)))
        allocate(bestss(0:getD(simp)))       
        bestss(0) = huge(1.0) 
        do trial = 1,nTrials
            do ii = 1,10000
                call suggestNextPoint( simp,ss )
                x = getChiSquared( ss )
                call addSimplexPoint( simp,ss,x )
                if (isConverged( simp, 1d-8,1d-8 )) exit
            end do
        !    call report(simp)
            call getBestSimplexPoint(simp,ss,x)
            if (x < bestss(0)) then 
                bestss(1:) = ss
                bestss(0) = x
                write (*,fmt='(a,i8)',advance="no") "   trial: ",trial
                call report(simp,o=2)
            end if
            
            if (trial < nTrials) then
                call delete(simp)
                call constructInitialSimplex( ln,lf,simp )
            end if
        end do
        ss = bestss(1:)
        print *,""
        print *,""
    !---


!------------------------------------------------------------------------------
!
!       report best solution
!
!------------------------------------------------------------------------------            
        select case(getD(simp))
            case (2)        !   mu,sig
                mu = ss(1) 
                sig = ss(2)                 
            case (4)        !   mu,sig,d0,w
                mu = ss(1) 
                sig = ss(2)               
                d0 = ss(3)
                w = ss(4)  
        end select
        ln = LogNormal_ctor(mu,sig)        
        print *,"iceberg info - best fit lognormal"
        call report(ln,o=1)
        if ( .not. isUnset(lf)) then
            if (logisticFuncFixed) then
                print *,"iceberg info - using fixed visibility function"
                call report(lf,o=1)
            else
                print *,"iceberg info - best fit visibility function"
                lf = LogisticFunction_ctor(d0,w)
                call report(lf,o=1)
            end if
        else
            print *,"iceberg info - no visibility function"
        end if
        call expectedHistogram( ln,lf,histogram_knot,histogram_knot_fine,expected )
        if (nPointDefectsFixed) then
            defect_count = expectedDefectCount( ln,omega0,expected_pointdefect_count )  
        else
            defect_count = observed_defect_count / sum(expected)
        end if
        expected = expected * defect_count
        ! print *,"iceberg info - best fit results histogram"
        ! write(*,fmt='(a16,a6,a16)') "min diam","count","expected"
        ! do jj = 0,NBINS-1
        !     write (*,fmt='(g16.8,i6,g16.8)') histogram_knot(jj),observed(jj),expected(jj)
        ! end do
        print *,"chi-square value            ",chi_square(observed,expected)
        print *,"alpha value                 ",alpha_score( chi_square(observed,expected), df=nBins-getD(simp) )
        print *,"observed defect count       ",observed_defect_count
        print *,"fitted defect count         ",defect_count
        print *,"observed point defect count ",observed_pointdefect_count
        print *,"fitted point defect count   ",expectedPointDefectCount( ln,omega0,defect_count )
        print *,"observed <d>,+/-,stdev(d)   ",observed_defect_mean,observed_defect_stdev/sqrt(observed_defect_count-1),observed_defect_stdev
        x = moment(ln,2,dzero) - moment(ln,1,dzero)**2       !   variance
        x = sqrt( defect_count*x/(defect_count-1) )          !   unbiassed stdev  
        print *,"fitted <d>,+/-              ",moment(ln,1,dzero),x/sqrt(defect_count)
        print *,"68% confidence fitted d     ",confidenceLevel(ln,0.16d0),":",confidenceLevel(ln,0.84d0)

        print *,""
        print *,""




!------------------------------------------------------------------------------
!
!       reconstruct histogram with even bins for output
!
!------------------------------------------------------------------------------                   
       print *,"iceberg info - best fit even bin histogram for output"
       ! ii = int( observed_defect_count*0.99 )      !   99% percentile
        if (maxBin == LIB_CLA_NODEFAULT_R) then
            maxBin = roundedBin( dat(n) )        
        end if
        do jj = 0,NBINS
            x = jj*maxBin/NBINS
            histogram_knot(jj) = x
        end do

        !print *,""
        do jj = 0,NBINS*NBINS
            x = jj*maxBin/(NBINS*NBINS)
            histogram_knot_fine(jj) = x
        !    write(*,fmt='(100f16.8)') x,func(ln,x),cdf(ln,x),cdf(ln,dzero,x)
        end do
        !print *,""
 
        observed = 0
        oversized = 0
        do ii = 1,n
            x = dat(ii)
            ok = .false.
            do jj = 0,NBINS-1
                if (x < histogram_knot(jj+1)) then
                    observed(jj) = observed(jj) + 1
                    ok = .true.
                    exit
                end if
            end do
            if (.not. ok) oversized = oversized + 1
        end do
        call expectedHistogram( ln,lf,histogram_knot,histogram_knot_fine,expected )
        if (nPointDefectsFixed) then
            defect_count = expectedDefectCount( ln,omega0,expected_pointdefect_count )  
        else
            defect_count = observed_defect_count / sum(expected)
        end if
        expected = expected * defect_count

        if (isUnset(lf)) then
            write(*,fmt='(6a16)') "#mean diam       ","npd  ","count  ","expected  "
            do jj = 0,NBINS-1
                x = ( histogram_knot(jj) + histogram_knot(jj+1) )/2
                write (*,fmt='(g16.8,g16.8,i16,3g16.8)') x,pointdefectcount(x,omega0),observed(jj),expected(jj)
            end do
            write (*,fmt='(a16,i16,2g16.8)') "#over"//repeat(" ",11),oversized!
            
        else
            call expectedHistogram( ln,LogisticFunction_ctor(),histogram_knot,histogram_knot_fine,expected_nolf )   
            !defect_count = expectedDefectCount( ln,omega0,expected_pointdefect_count ) / sum(expected)             
            expected_nolf = expected_nolf * defect_count        !   expected_nolf should give integral 1
            write(*,fmt='(6a16)') "#mean diam       ","npd  ","count  ","expected  ","vis func  ","exp if all vis"
            do jj = 0,NBINS-1
                x = ( histogram_knot(jj) + histogram_knot(jj+1) )/2
                write (*,fmt='(g16.8,g16.8,i16,4g16.8)') x,pointdefectcount(x,omega0),observed(jj),expected(jj),func(lf,x),expected_nolf(jj)
            end do            

            write (*,fmt='(a16,i6,2g16.8)') "#over"//repeat(" ",11),oversized!,sum(expected)*(1 - x),defect_count - sum(expected_nolf)
        end if
        print *,""
        print *,""
        


!------------------------------------------------------------------------------
!
!       bye bye
!
!------------------------------------------------------------------------------                           
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

 

        subroutine expectedHistogram( logn,logf,histogram_knot,histogram_knot_fine,expected )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given an input lognormal probability distribution, 
    !*      and an expected logistic function
    !*      and a set of knot points,
    !*      return the expected probability in each bin- note this needs to be scaled to give counts in each bin
            type(LogNormal),intent(in)          ::      logn
            type(LogisticFunction),intent(in)   ::      logf
            !real(kind=real64),intent(in)        ::      defect_count
            real(kind=real64),dimension(0:),intent(in)  ::      histogram_knot,histogram_knot_fine
            real(kind=real64),dimension(0:),intent(out) ::      expected
            integer             ::      nBins,jj,kk
            real(kind=real64)   ::      c_old,cc,xx,gg,x_old,intgrl,cdfdzero

            nBins = size(histogram_knot) - 1
            

            cdfdzero = cdf(logn,dzero)
           
            if (isUnset(logf)) then

                c_old = - cdfdzero
                do jj = 0,nBins-1
        
                    xx = histogram_knot(jj+1)
                    cc = cdf( logn,xx ) - cdfdzero                !   cdf is integral
                    expected(jj) = (cc - c_old) 
                    c_old = cc

                end do

            else
                x_old = 0.0d0
                c_old = - cdfdzero
                do jj = 0,nBins-1

                    intgrl = 0.0d0
                    do kk = jj*nBins,jj*nBins+nBins-1
                        xx = histogram_knot_fine(kk+1)
                        gg = func( logf,(x_old + xx)/2 )            !   compute logistic func at midpoint
                        cc = cdf( logn,xx ) - cdfdzero              !   cdf is integral
                        intgrl = intgrl + gg * (cc - c_old) 
                        c_old = cc
                        x_old = xx
                    end do
                    expected(jj) = intgrl

                end do
            end if
             
            return
        end subroutine expectedHistogram
 


        pure real(kind=real64) function pointDefectCount(d,omega0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return count of point defects given vol per atom 
            real(kind=real64),intent(in)            ::      d   !   diam
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            !real(kind=real64),parameter     ::  PI = 3.141592654d0
          
            if (voids) then
                pointDefectCount = ( PI * d**3 ) / ( 6 * omega0 )
            else
                pointDefectCount = ( PI * b * d**2 ) / (4 * omega0)
            end if

            return
        end function pointDefectCount



        pure real(kind=real64) function defectDiameter(n,omega0)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return diameter given pd count
    !*      inverse func of pointDefectCount()
            real(kind=real64),intent(in)            ::      n   !   diam
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            !real(kind=real64),parameter     ::  PI = 3.141592654d0
          
            if (voids) then
                defectDiameter = ( 6 * n * omega0 / PI )**(1.0d0/3)
            else
                defectDiameter = sqrt( 4 * n * omega0 / ( b * PI ) )
            end if

            return
        end function defectDiameter        

        real(kind=real64) function expectedDefectCount( logn,omega0,npd )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given we expect to have npd point defects, return the expected defect count
    !*          npd = Nvoid (pi/6) (1/omega0) <d³>
            type(LogNormal),intent(in)              ::      logn
            real(kind=real64),intent(in)            ::      omega0  !  volume per atom
            real(kind=real64),intent(in)            ::      npd     !  expected pd count
            !real(kind=real64),parameter             ::  PI = 3.141592654d0
            real(kind=real64)           ::      xx

            if (voids) then
                xx = moment(logn,3,dzero)
                expectedDefectCount = npd * 6 * omega0 / ( PI * xx )
            else
                xx = moment(logn,2,dzero)
                expectedDefectCount = npd * 4 * omega0 / ( PI * b * xx )
            end if

            return
        end function expectedDefectCount


        real(kind=real64) function expectedPointDefectCount( logn,omega0,defect_count )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      reverse function of expectedDefectCount() - gets point defecvts from total defect count.
            type(LogNormal),intent(in)              ::      logn
            real(kind=real64),intent(in)            ::      omega0          !  volume per atom
            real(kind=real64),intent(in)            ::      defect_count     !  expected defect count
            !real(kind=real64),parameter             ::  PI = 3.141592654d0
            real(kind=real64)           ::      xx

            if (voids) then
                xx = moment(logn,3,dzero)
                expectedPointDefectCount = ( defect_count * PI * xx )/( 6 * omega0 )
            else
                xx = moment(logn,2,dzero)
                expectedPointDefectCount = ( defect_count * PI * b * xx )/( 4 * omega0 )
            end if
           

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
                    if ( isUnset(lf)) then
                        logf = LogisticFunction_ctor()
                    else
                        logf = LogisticFunction_ctor(d0,w)
                    end if
                case (4)        !   mu,sig,d0,w
                    logn = Lognormal_ctor(sanitisedx(1),sanitisedx(2))
                    logf = LogisticFunction_ctor(sanitisedx(3),sanitisedx(4))
            end select

            call expectedHistogram( logn,logf,histogram_knot,histogram_knot_fine,expected )
            if (nPointDefectsFixed) then
                defect_count = expectedDefectCount( logn,omega0,expected_pointdefect_count )  
            else
                defect_count = observed_defect_count / sum(expected)
            end if
            expected = expected * defect_count

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
                !if (nPointDefectsFixed .or. nDefectsObservedFixed) then
                    !   only adjusting lognormal mu,sig
                    dd = 2
                ! else 
                !     !   adjusting mu,sig, and defect count 
                !     dd = 3
                !end if
            else 
                !if (nPointDefectsFixed .or. nDefectsObservedFixed) then
                    !   adjusting lognormal mu,sig + logistic d0,w
                    dd = 4
                ! else 
                !     !   adjusting mu,sig, and defect count + logistic d0,w
                !     dd = 5
                !end if
            end if

            nn = dd * 2
            simp = DownhillSimplex_ctor(nn,dd)
            allocate(xx(dd))

            do ii = 1,nn

                call random_number(xx)
                xx = max( 1.0d-3, xx * 4 )
                select case(dd)
                    case (2)        !   mu,sig
                        xx(1) = getMu(logn)*xx(1)
                        xx(2) = getSigma(logn)*xx(2)
                    case (4)        !   mu,sig,d0,w
                        xx(1) = getMu(logn)*xx(1) 
                        xx(2) = getSigma(logn)*xx(2) 
                        xx(3) = getx0(logf)*xx(3) 
                        xx(4) = getw(logf)*xx(4) 
                end select

                cc = getChiSquared( xx )
                call addSimplexPoint(simp,ii,xx,cc)

            end do
       

            return
        end subroutine constructInitialSimplex
    
        


    end program iceberg