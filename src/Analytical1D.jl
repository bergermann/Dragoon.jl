###		1-dimensional boostfactor calculation, matrix formalism

function transform_surface(freq::AbstractArray, gamma::AbstractArray;
    axion_out::AbstractArray = [],
    eps_i::Union{AbstractFloat, Complex} = 1.,
    eps_m::Union{AbstractFloat, Complex} = 1.,
    l::AbstractFloat = 0.,
    surfaceloss::AbstractFloat = 0.,
    USE_RAY = true)

    # gamma same length as freq

    c = 299792458.
    #println("STARTING TRANSFORM_SURFACE")
    if USE_RAY
        # Ray Tracing

        r_surface = (sqrt(eps_m) - sqrt(eps_i)) / (sqrt(eps_m) + sqrt(eps_i))
        r_internal = gamma

		#println("r_surface is: ")
        #println(r_surface)
        #println("r_internal is: ")
        #println(r_internal)

        # Any transmission or reflection factor simply gets multiplied with this loss factor
        g = (1. - surfaceloss)

		#println("g is: ")
        #println(g)
        # The full relflectivity in front of the surface
        gamma = g^2 * (1. + r_surface) * (1. - r_surface) * r_internal ./ (1. .+ r_internal * r_surface * g) .+ g * r_surface
		#println("gamma is: ")
        #println(gamma)
	else
        # Impedance Transformation

        # Calculate impedance from new reflection coefficient
        z = (1. .+ gamma) ./ (1. .- gamma)
		#println("z is: ")
        #println(z)

        Z_i = 1. / sqrt(eps_i)
        Z_m = 1. / sqrt(eps_m)

		#println("Z_i is: ")
        #println(Z_i)
        #println("Z_m is: ")
        #println(Z_m)

        # Change basis to new medium
        z *= Z_i / Z_m
		#println("z is: ")
        #println(z)
        # Calculate reflection coefficient
        # According to ray calculation that is not fully correct, only if gamma = 1
        gamma = ((z .- 1.) ./ (z .+ 1.)) .* (1. - surfaceloss)
		#println("gamma is: ")
        #println(gamma)
    end

    # Propagation (bidirectional)
    wavel = c ./ (freq .* sqrt(eps_m))
	#println("wavel is: ")
    #println(wavel)

    gamma = complex(gamma)

    for i in eachindex(wavel)
        gamma[i] *= exp(-1.0im * 4. * pi * l / wavel[i])
    end

	#println("gamma is: ")
    #println(gamma)

    # Calculate axion contribution
    if !isempty(axion_out)
        if !USE_RAY
            throw(ArgumentError("Need to use ray approx."))
        end

        denominator = eps_i * sqrt(eps_m) + eps_m * sqrt(eps_i)
		#println("denominator is: ")
        #println(denominator)
        ax_i_t = sqrt(eps_i) * (1 - eps_m / eps_i) / denominator
		#println("ax_i_t is: ")
        #println(ax_i_t)
        ax_m = sqrt(eps_m) * (1 - eps_i / eps_m) / denominator
		#println("ax_m is: ")
        #println(ax_m)

		ax_i = fill!(similar(r_internal), ax_i_t)


        #Propagation of the internal new axion field inside and through the surface
        ax_i .= ax_i_t * g * r_internal * (1. - r_surface) ./ (r_internal * r_surface * g .+ 1.)
		#println("ax_i is: ")
        #println(ax_i)
		#ax_i *= -g*(1-r_surface)*(-g*r_surface*r_internal)/((-g*r_surface*r_internal)-1)

        # Propagation of the present axion_out through the surface

		#println("---------------------------------------------------------------------------------")
		#println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		#println("g is:")
        #println(g)
        #println("r_surface is:")
        #println(r_surface)
        #println("r_internal is:")
        #println(r_internal)
		#println("PREVIOUS axion_out is:")
        #println(axion_out)
		@. axion_out *= g * (1 - r_surface) / (r_internal * r_surface * g + 1.)
		#println("POST axion_out is: ")
        #println(axion_out)

        # Total outgoing axion amplitude
        axion_out .+= ax_i .+ ax_m
		#println("axion_out is: ")
        #println(axion_out)

        # display(axion_out)
        # Propagate the axion wave
        # Propagates in the direction of the reflected wave therefore the leading sign in the exponent is minus

		# axion_out = complex(axion_out)

        for i in eachindex(wavel)
            axion_out[i] *= exp(-1.0im * 2. * pi * l / wavel[i])
        end
        
        # display(axion_out)
		#println("axion_out is: ")
        #println(axion_out)
    end

	#println("The final result of this function is: ")
    #println(gamma)
    #println("ENDING TRANSFORM_SURFACE")

    gamma
end


function lossy_eps(freq::AbstractArray,
    eps::Union{AbstractFloat, Complex},
    tand::AbstractFloat;
    SvenssonDjordjevic = false)

	#println("STARTING LOSSY_EPS")

	if !SvenssonDjordjevic
        result =  (eps  -1.0im * eps * tand)
    else
        FreqForEpsrTanD = 1e9 #Hz
        HighFreqForTanD = 1e12 #Hz
        LowFreqForTanD = 1e3 #Hz
        # Svensson/Djordjevic Model
        # http://edadocs.software.keysight.com/display/ads2009/About+Dielectric+Loss+Models
        L = log((HighFreqForTanD + 1.0im * FreqForEpsrTanD) / (LowFreqForTanD + 1.0im * FreqForEpsrTanD))
		#println("L is: ")
        #println(L)
		a = - (eps * tand) / imag(L)
		#println("a is: ")
        #println(a)
		Einf = eps - a * real(L)

		#println("Einf is: ")
        #println(Einf)


		result = complex(fill!(similar(freq), 0.))
		for i in eachindex(freq)
			result[i] = (Einf  + a * log((HighFreqForTanD + 1.0im * freq[i]) / (LowFreqForTanD + 1.0im * freq[i])))
		end
    end

	#println("The final result of this function is: ")
	#println(result)
	#println("ENDING LOSSY_EPS")
	result
end


function disk_system(freq::AbstractArray;
    tand::AbstractFloat = 0.,
    num_disk::Integer = 5,
    non_uniform_surfaceloss::AbstractArray = [],
    spacings::AbstractArray = [],
    mirror = true,
    disk_thickness::AbstractFloat = 0.001,
    disk_epsilon::AbstractFloat = 9.,
    kwargs...)

	#println("STARTING DISK_SYSTEM")

    if isempty(spacings)
        spacings = 8e-3 * ones(num_disk + 1)
		#println("spacings is: ")
        #println(spacings)
    end

    Z_0 = 0.001
    eps_1 = lossy_eps(freq, disk_epsilon, 0.; SvenssonDjordjevic = false) #9.
	#println("eps_1 is: ")
    #println(eps_1)
	l_1 = 0.001
    # Traditional Loss Model
    eps_2 = 1.
    eps_2_tr = lossy_eps(freq, 1.0, tand; SvenssonDjordjevic = false)
	#println("eps_2_tr is: ")
    #println(eps_2_tr)
	#eps_2_sd = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic = true)
    l_2 = 0.

    # Air 0
    #gamma = transform_surface(freq, 0., eps_i=1., eps_m=1., l=0, **kwargs)

    #if mirror:
    #    axion_out=np.ones(freq.shape) + 0*1j
    #else:
    # axion_out = fill!(similar(freq), 0.) .+ 0.0 * 1.0im
    axion_out = zeros(ComplexF64,length(freq))

    # Metal Disk
    if mirror
        gamma = transform_surface(freq, zeros(ComplexF64,length(freq)); axion_out = axion_out, eps_i = 1., eps_m = 1e20, l = 100e-3, kwargs...)
		#println("gamma is: ")
        #println(gamma)
	else
        gamma = fill!(similar(freq), 0.)
		#println("gamma is: ")
        #println(gamma)
    end

    ##printlnln(axion_out)

    if isempty(non_uniform_surfaceloss)
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = (mirror ? 1e20 : 1.), eps_m = eps_2_tr, l = spacings[1], kwargs...)
		#println("gamma is: ")
        #println(gamma)
        for i in collect(1:num_disk)
            # Disk i+1
            # display(axion_out)
            gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = eps_2_tr, eps_m = eps_1, l = disk_thickness, kwargs...)
            # Air i+2
            # display(axion_out)
            gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = eps_1, eps_m = eps_2_tr, l = spacings[i+1], kwargs...)
        end
		#println("gamma is: ")
        #println(gamma)
    else
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = 1e20, eps_m = eps_2_tr, l = spacings[1], surfaceloss = non_uniform_surfaceloss[1], kwargs...)

        for i in collect(1:num_disk)
            # Disk i+1
            gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = eps_2_tr, eps_m = eps_1, l = disk_thickness,  surfaceloss = non_uniform_surfaceloss[2 * i ], kwargs...)
			#println("gamma is: ")
	        #println(gamma)
			# Air i+2
            gamma = transform_surface(freq, gamma; axion_out = axion_out, eps_i = eps_1, eps_m = eps_2_tr, l = spacings[i+1],  surfaceloss = non_uniform_surfaceloss[2 * i + 1], kwargs...)
			#println("gamma is: ")
	        #println(gamma)
		end
    end

	#println("ENDING DISK_SYSTEM")

    return gamma, axion_out


end


function disk_system_phase_depths(d_air::AbstractFloat;
    d_disk::AbstractFloat = pi / 2.,
    tand::AbstractFloat = 0.,
    num_disk::Integer = 5,
    non_uniform_surfaceloss::AbstractArray = [],
    kwargs...)


    c = 299792458.

    # The result should be independent of frequency now...
    freq = [20e9]
    wavel = c ./ (freq)
    l_air  = d_air .* wavel ./ (2. * pi)
    l_disk = d_disk .* wavel ./ (2. * pi * sqrt(9.4))


    Z_0 = 0.001
    eps_1 = lossy_eps(freq, 9.4, 0., SvenssonDjordjevic = false) #9.
    l_1 = 0.001
    # Traditional Loss Model
    eps_2 = 1.
    eps_2_tr = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic = false)
    #eps_2_sd = lossy_eps(freq, 1.0, tand, SvenssonDjordjevic = true)
    l_2 = 0.


    # Air 0
    #gamma = transform_surface(freq, 0., eps_i=1., eps_m=1., l=0, **kwargs)

    axion_out = fill!(similar(freq), 0.) .+ 0.0 * 1.0im
    # Metal Disk
    gamma = transform_surface(freq, fill!(similar(freq), 0.), axion_out = axion_out, eps_i = 1., eps_m = 1e20, l = 100e-3, kwargs)

    if isemtpy(non_uniform_surfaceloss)
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out = axion_out, eps_i = 1e20, eps_m = eps_2_tr, l = l_air, kwargs)

        for i in collect(0:num_disk - 1)
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out = axion_out, eps_i = eps_2_tr, eps_m = eps_1, l = l_disk, kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out = axion_out, eps_i = eps_1, eps_m = eps_2_tr, l = l_air, kwargs)
        end
    else
        # Air 1
        # 1. = Z_0 / Z_0 (at the beginning the normalized impedance is always 1)
        gamma = transform_surface(freq, gamma, axion_out=axion_out, eps_i = 1e20, eps_m = eps_2_tr, l = 8e-3, surfaceloss = non_uniform_surfaceloss[0], kwargs)

        for i in collect(0:num_disk - 1)
            # Disk i+1
            gamma = transform_surface(freq, gamma, axion_out = axion_out, eps_i = eps_2_tr, eps_m = eps_1, l = l_disk,  surfaceloss = non_uniform_surfaceloss[2*i+1], kwargs)
            # Air i+2
            gamma = transform_surface(freq, gamma, axion_out = axion_out, eps_i = eps_1, eps_m = eps_2_tr, l = l_air,  surfaceloss = non_uniform_surfaceloss[2*i+2], kwargs)
        end
    end

    return gamma, axion_out
end
