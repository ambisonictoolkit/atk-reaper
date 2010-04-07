def fir_lp(N, Wn, width=pi):
    """fir_lp(N, Wn, width=pi)

    Lowpass FIR Filter Design using windowed ideal filter method.
    (A wrapper for firwin)
    
    Inputs:
    
      N  -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    return firwin(N, Wn, window = ('kaiser', width))


def fir_hp(N, Wn, width=pi, method = 'hilbert'):
    """fir_hp(N, Wn, width=pi, method = 'hilbert')

    Highpass FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
      method  -- highpass methods: 'hilbert', 'sinc', qmf'. For 'hilbert' and 'sinc', 
                  fir_hp and fir_lp filters are matched for allpass reconstruction
                  with matching parameters. 'sinc' gives better performance at high
                  frequencies, 'hilbert' gives better performance at low frequencies.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if method is 'qmf' or method is 'hilbert':

        lp_Wn = 1. - Wn
        hp_b = qmf(fir_lp(N, lp_Wn, width))

        if method is 'hilbert':
            
            if (N % 2) is 0:
                hp_h = hilbert(hp_b)
                hp_b = cos(pi * Wn) * real(hp_h) + sin(pi * Wn) * imag(hp_h)

        if ((N - 1)) / 2 % 2 is 1:
            hp_b *= -1.

    else:
        lp_b = fir_lp(N, Wn, width)

        if (N % 2) is 0:
            ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)
        else:
            ap_b = zeros(N)
            ap_b[N / 2] = 1.

        hp_b = ap_b - lp_b

    return hp_b


def fir_ls(N, Wn, k, width=pi, method = 'hilbert'):
    """fir_ls(N, Wn, k, width=pi, method = 'hilbert')

    Low shelf FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)

      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at low frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
      method  -- highpass methods: 'hilbert' or 'sinc'. 'sinc' gives better
                  performance at high frequencies, 'hilbert' gives better
                  performance at low frequencies.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    lp_b = fir_lp(N, Wn, width)

    if method is 'hilbert':
        hp_b = fir_hp(N, Wn, width, method)

        ls_b = k * lp_b + hp_b

    else:
        if (N % 2) is 0:
            ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

        else:
            ap_b = zeros(N)
            ap_b[N / 2] = 1.

        ls_b = ap_b + (k - 1.) * lp_b

    return ls_b


def fir_hs(N, Wn, k, width=pi, method = 'hilbert'):
    """fir_hs(N, Wn, k, width=pi, method = 'hilbert')

    High shelf FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N  -- order of filter (number of taps)

      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at high frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
    
      method  -- highpass methods: 'hilbert' or 'sinc'. 'sinc' gives better
                  performance at high frequencies, 'hilbert' gives better
                  performance at low frequencies.
    
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    lp_b = fir_lp(N, Wn, width)

    if method is 'hilbert':
        hp_b = fir_hp(N, Wn, width, method)

        hs_b = lp_b + k * hp_b

    else:
        if (N % 2) is 0:
            ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

        else:
            ap_b = zeros(N)
            ap_b[N / 2] = 1.

        hs_b = lp_b + k * (ap_b - lp_b)

    return hs_b


def fir_bp(N, Wn, bw, width=pi, method = 'cos'):
    """fir_bp(N, Wn, bw, width=pi, method = 'cos')

    Bandpass FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- center frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

      method -- Bandpass methods: 'cos' or 'dual'.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if method is 'cos':

        lp_Wn = bw / 2
        lp_b = fir_lp(N, lp_Wn, width)

        phi = ((1 - N)  / 2.) * pi * Wn

        Wp_b = 2 * fcososc(Wn, phi, N)

        bp_b = Wp_b * lp_b

    else:
        lp1_Wn = Wn + bw / 2
        lp0_Wn = Wn - bw / 2

        bp_b = fir_lp(N, lp1_Wn, width) - fir_lp(N, lp0_Wn, width)

    return bp_b


def fir_bs(N, Wn, bw, width=pi, method = 'cos'):
    """fir_bs(N, Wn, bw, width=pi, method = 'cos')

    Bandstop FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.
        
      method -- Bandpass methods: 'cos', 'hilbert' or 'sinc'. The last two
                  methods determine allpass technique for even N.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if method is 'cos':
        bp_b = fir_bp(N, Wn, bw, width)

        if N % 2 is 1:
            ap_b = zeros(N)
            ap_b[N / 2] = 1.
        else:
            ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

        bs_b = ap_b - (bp_b)

    else:
        hp_Wn = Wn + bw / 2
        lp_Wn = Wn - bw / 2

        bs_b = fir_lp(N, lp_Wn, width) + fir_hp(N, hp_Wn, width, method)

    return bs_b


def fir_pk(N, Wn, bw, k, width=pi, method = 'cos'):
    """fir_pk(N, Wn, bw, k, width=pi, method = 'cos')

    Peaking FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at peaking frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

      method -- Peaking methods: 'cos', 'hilbert' or 'sinc'. The last two
                  methods determine allpass technique for even N.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bp_b = fir_bp(N, Wn, bw, width, method)

    if (N % 2) is 0:
        ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

    else:
        ap_b = zeros(N)
        ap_b[N / 2] = 1.

    pk_b = ap_b + (k - 1.) * bp_b

    return pk_b


def fir_sk(N, Wn, bw, k, width=pi, method = 'cos'):
    """fir_sk(N, Wn, bw, k, width=pi, method = 'cos')

    Skirting FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- cutoff frequency of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      bw -- bandwidth of filter (normalized so that 1 corresponds to
                  Nyquist or pi radians / sample)
    
      k -- scale at peaking frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

      method -- Skirting methods: 'cos', 'dual' or 'hilbert'.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bs_b = fir_bs(N, Wn, bw, width, method)

    if (N % 2) is 0:
        ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

    else:
        ap_b = zeros(N)
        ap_b[N / 2] = 1.

    sk_b = ap_b + (k - 1.) * bs_b

    return sk_b


def fir_bps(N, Wn, bw, k, width=pi, method = 'cos'):
    """fir_bps(N, Wn, bw, k, width=pi, method = 'cos')

    Bandpass bank FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- an array of center frequencies of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      bw -- an array of cutoff bandwidths of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      k -- an array of scales at center frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

      method -- Bandpass methods: 'cos' or 'dual'.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    bps_b = zeros(N)

    for n in range(len(Wn)):
        bps_b += k[n] * fir_bp(N, Wn[n], bw[n], width, method)

    return bps_b


def fir_pks(N, Wn, bw, k, width=pi, method = 'cos'):
    """fir_pks(N, Wn, bw, k, width=pi, method = 'cos')

    Peaking bank FIR Filter Design using windowed ideal filter method.
    
    Inputs:
    
      N      -- order of filter (number of taps)
      Wn -- an array of center frequencies of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      bw -- an array of cutoff bandwidths of filter (normalized so that
                 1 corresponds to Nyquist or pi radians / sample)
    
      k -- an array of scales at center frequencies
    
      width -- beta for Kaiser window FIR design.
                  pi = minimum ripple for steepest cutoff.

      method -- Bandpass methods: 'cos' or 'dual'.
        
    Outputs:
    
      b      -- coefficients of length N FIR filter.
    """
    if (N % 2) is 0:
        ap_b = kaiser(N, pi) * sinc(arange(N) - N / 2 + .5)

    else:
        ap_b = zeros(N)
        ap_b[N / 2] = 1.

    pks_b = ap_b

    for n in range(len(Wn)):
        pks_b += (k[n] - 1.) * fir_bp(N, Wn[n], bw[n], width, method)

    return pks_b
