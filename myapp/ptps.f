C FILE: ptps.f
      subroutine ptps(spikes, n_electrodes, ec_traces, vscale_HP,
     / max_ptp, max_ch)
                    
        integer spikes(:,:)
        real*4 vscale_HP
        real*4 ptp(1000)
        real*4 max_ptp(1000)
        integer n_electrodes, u, i, k
        integer max_ch(1000,1)
        integer*2 ec_traces(:,:)

Cf2py intent(in) spikes, n_electrodes, n_vd_samples, f_name
Cf2py intent(in) ec_traces, vscale_HP
Cf2py intent(out)  max_ptp, max_ch


        !write(*,*) spikes(28,168)
        !write(*,*) spikes(28,169)
        !write(*,*) spikes(28,170)
        !write(*,*) spikes(28,171)
        !return
        do u = 1, size(spikes,1)
            ptp = 0
            do i = 1,n_electrodes
            !NB: Must set n_spikes to max to catch largest unit
                n_spikes = size(spikes(u,:)) 
                do k = 1, size(spikes(u,:))
                    if (spikes(u,k).eq.0) then
                        n_spikes = k-1
                        exit
                    endif
                    if((spikes(u,k).ge.10).and.((spikes(u,k)-10).le.
     &               size(ec_traces(i,:)))) then
                        ptp(i) = ptp(i) + 
     &   vscale_HP*(maxval(ec_traces(i,spikes(u,k)-10:spikes(u,k)+10)) -
     &              minval(ec_traces(i,spikes(u,k)-10:spikes(u,k)+10)))
                    endif
                enddo
            enddo
            max_ptp(u) = maxval(ptp)/n_spikes
            max_ch(u,:) = maxloc(ptp)
            write(*,*) "cell: ", u, " n_spikes: ", n_spikes, 
     &       max_ptp(u), max_ch(u,:)
        enddo
        
      end subroutine ptps
C END FILE ptps.f
