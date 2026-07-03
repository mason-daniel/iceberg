!*  A module containing a single simple function, an exponential that won't overflow
!*      y = safeExp(x) 
!*                          = exp(x)        if x is within bounds
!*                          = 0             if x too negative
!*                          = huge(1.0)     if x too positive       ( not huge(1.0d0) so that we can find safe_exp(x1) + safe_exp(x2)
!*  
!*  Set PASS_INVALID_ARG = .true. if you want to see what happens when an invalid argument is passed to exponential - could help debugging?
!*  
!*  Daniel Mason
!*  (c) UKAEA March 2025
!*
!*  version history
!*      0.0.1       March 2025      First working version
!*

    module Lib_SafeExp
!---^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private

        real(kind=real64),private,parameter             ::      MIN_ARG64 = log( tiny(1.0_real64) )
        real(kind=real64),private,parameter             ::      MAX_ARG64 = log( real( huge(1.0_real32),kind=real64 ) )
        
        real(kind=real128),private,parameter             ::      MIN_ARG128 = log( tiny(1.0_real128) )
        real(kind=real128),private,parameter             ::      MAX_ARG128 = log( real( huge(1.0_real64),kind=real64 ) )
        
        public          ::      safeExp


#ifdef DEBUG        
        logical,private,parameter   ::      PASS_INVALID_ARG = .false.
        integer,private             ::      static_underflow_count = 0
        integer,private             ::      static_overflow_count = 0
#endif

        interface   safeExp
            module procedure    safeExp64
            module procedure    safeExp128
#ifdef DEBUG
            module procedure    safeExp64a
            module procedure    safeExp128a
#endif
        end interface


    contains
!---^^^^^^^^

#ifdef DEBUG

    real(kind=real64) function safeExp64(x)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real64),intent(in)            ::      x
        real(kind=real64)                       ::      safeArg

        safeArg = x
        if (x<MIN_ARG64) then
            static_underflow_count = static_underflow_count + 1
            if (static_underflow_count == 1) then
                print *,"Lib_SafeExp::safeExp64 warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_underflow_count )) &
                        write(0,fmt='(a,g20.6,a,i20)') "Lib_SafeExp::safeExp64 WARNING - arg = ",x,", underflow count = ",static_underflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MIN_ARG64
        else if (x>MAX_ARG64) then
            static_overflow_count = static_overflow_count + 1
            if (static_overflow_count == 1) then
                print *,"Lib_SafeExp::safeExp64 warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_overflow_count )) &
                    write(0,fmt='(a,g20.6,a,i20)') "Lib_SafeExp::safeExp64 warning - arg = ",x,", overflow count = ",static_overflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MAX_ARG64
        else if (x/=x) then
            if (.not. PASS_INVALID_ARG) stop "Lib_SafeExp::safeExp64 ERROR - arg = NaN"
        end if
            
        safeExp64 = exp( safeArg )

        return
    end function safeExp64

    real(kind=real128) function safeExp128(x)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real128),intent(in)            ::      x
        real(kind=real128)                       ::      safeArg

        safeArg = x
        if (x<MIN_ARG128) then
            static_underflow_count = static_underflow_count + 1
            if (static_underflow_count == 1) then
                print *,"Lib_SafeExp::safeExp128 warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_underflow_count )) &
                    write(0,fmt='(a,g20.6,a,i20)') "Lib_SafeExp::safeExp128 WARNING - arg = ",x,", underflow count = ",static_underflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MIN_ARG128
        else if (x>MAX_ARG128) then
            static_overflow_count = static_overflow_count + 1
            if (static_overflow_count == 1) then
                print *,"Lib_SafeExp::safeExp128 warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_overflow_count )) &
                    write(0,fmt='(a,g20.6,a,i20)') "Lib_SafeExp::safeExp128 warning - arg = ",x,", overflow count = ",static_overflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MAX_ARG128
        else if (x/=x) then
            if (.not. PASS_INVALID_ARG) stop "Lib_SafeExp::safeExp128 ERROR - arg = NaN"
        end if
            
        safeExp128 = exp( safeArg )

        return
    end function safeExp128

    function safeExp64a(x) result(se)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real64),dimension(:),intent(in)            ::      x
        real(kind=real64),dimension(size(x))                 ::      se
        integer     ::      ii
        do ii = 1,size(x)
            se(ii) = safeExp64(x(ii))
        end do        
        return
    end function safeExp64a

    function safeExp128a(x) result(se)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real128),dimension(:),intent(in)            ::      x
        real(kind=real128),dimension(size(x))                 ::      se
        integer     ::      ii
        do ii = 1,size(x)
            se(ii) = safeExp128(x(ii))
        end do        
        return
    end function safeExp128a

#else

        elemental real(kind=real64) function safeExp64(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)            ::      x
            real(kind=real64)                       ::      safeArg

            safeArg = max( MIN_ARG64, min( MAX_ARG64,x ) )

            safeExp64 = exp( safeArg )

            return
        end function safeExp64

        elemental real(kind=real128) function safeExp128(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real128),intent(in)            ::      x
            real(kind=real128)                       ::      safeArg

            safeArg = max( MIN_ARG128, min( MAX_ARG128,x ) )

            safeExp128 = exp( safeArg )

            return
        end function safeExp128

#endif


    end module Lib_SafeExp
        

