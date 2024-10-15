module eps
   !
   ! ==============================================================
   ! This module contains wavelength-dependnet dielectric functions
   ! epsXX(lambda) for various systems XX. These functions use the
   ! analytical expression given in Eric and Pablo's book.
   ! ==============================================================
   !
   implicit none
   !
   private
   !
   public :: epsAu, epsAg, epsPd, epsPt, epsSi, epsAl, epsCr, epsWater, interp1, epsDrude
   !
   real(8), parameter :: pi = 3.141592653589793d0 ! pi
   real(8), parameter :: c = 2.99792458d+17 ! nm/s
   real(8), parameter :: h = 4.135667662d-15 ! eV*s
   real(8), parameter :: hc = 1.2398419738620932d+3 ! nm*eV
   complex(8), parameter :: imu = (0, 1) ! the imaginary unit sqrt(-1)
   !
contains

!
   function epsAu(wavelength) result(eps)
      !
      ! Returns the wavelength-dependent relative dielectric function of gold
      ! This function uses the analytical expression given in Eq. (E.2).
      ! The exp(-i omega t) convention is assumed.
      !
      ! Input:
      ! - wavelength: real scalar,
      !               wavelength in NANOMETERS (nm)
      !
      ! Returns:
      ! - eps: complex scalar,
      !        epsilon(wavelength) as a complex number.
      !
      real(8), intent(in) :: wavelength
      real(8), parameter :: eps_infty = 1.54d0
      real(8), parameter :: lambda_p = 177.5d0
      real(8), parameter :: mu_p = 14500.0d0
      real(8), parameter :: A1 = 1.27d0
      real(8), parameter :: lambda1 = 470.0d0
      real(8), parameter :: mu_p1 = 1900.0d0
      real(8), parameter :: A2 = 1.1d0
      real(8), parameter :: lambda2 = 325.0d0
      real(8), parameter :: mu_p2 = 1060.0d0
      real(8), parameter :: phi = -pi/4.0d0
      complex(8) :: eps, etp, etm
      !
      etp = exp(imu*phi)
      etm = 1.0d0/etp ! =exp(-imu*phi)
      !
      eps = eps_infty*(1 - 1/(lambda_p**2*((1/wavelength)**2 + imu/(mu_p*wavelength)))) &
            + A1/lambda1*(etp/(1/lambda1 - 1/wavelength - imu/mu_p1) + &
                          etm/(1/lambda1 + 1/wavelength + imu/mu_p1)) &
            + A2/lambda2*(etp/(1/lambda2 - 1/wavelength - imu/mu_p2) + &
                          etm/(1/lambda2 + 1/wavelength + imu/mu_p2))
      !
   end function epsAu
   !
   function epsAg(wavelength) result(eps)
      !
      ! Returns the wavelength-dependent relative dielectric function of silver
      ! This function uses the analytical expression given in Eq. (E.1).
      !
      ! Input:
      ! - wavelength: real scalar,
      !               wavelength in NANOMETERS (nm)
      !
      ! Returns:
      ! - eps:   complex scalar,
      !          epsilon(wavelength) as a complex number.
      !
      real(8), intent(in) :: wavelength
      real(8), parameter :: eps_infty = 4.0d0
      real(8), parameter :: lambda_p = 282.0d0
      real(8), parameter :: mu_p = 17000.0d0
      complex(8) :: eps
      !
      eps = eps_infty*(1 - 1/(lambda_p**2*((1/wavelength)**2 + imu/(mu_p*wavelength))))
      !
   end function epsAg
   !
   function epsPt(wavelength) result(eps)
      !
      ! Returns the wavelength-dependent relative dielectric function
      ! of a Lorentz-Drude metal, with the parameters for Pt from:
      ! Rakic et al. 1998, https://doi.org/10.1364/AO.37.005271
      !
      real(8), intent(in) :: wavelength ! in nm
      real(8), parameter :: &
         wp = 9.59d0, & ! eV
         f0 = 0.333d0, &
         G0 = 0.080d0 ! eV
      real(8), parameter :: &
         f1 = 0.191d0, &
         G1 = 0.517d0, & ! eV
         w1 = 0.780d0 ! eV
      real(8), parameter :: &
         f2 = 0.659d0, &
         G2 = 1.838d0, & ! eV
         w2 = 1.314d0 ! eV
      real(8), parameter :: &
         f3 = 0.547d0, &
         G3 = 3.668d0, & ! eV
         w3 = 3.141d0 ! eV
      real(8), parameter :: &
         f4 = 3.576d0, &
         G4 = 8.517d0, & ! eV
         w4 = 9.249d0 ! eV
      real(8) :: Op, w ! omega
      complex(8) :: eps
      !
      Op = sqrt(f0)*wp
      w = hc/wavelength
      !
      eps = 1 - Op**2/(w*(w + cmplx(0, G0, kind(eps))))
      eps = eps + f1*wp**2/((w1**2 - w**2) - cmplx(0, w*G1, kind(eps)))
      eps = eps + f2*wp**2/((w2**2 - w**2) - cmplx(0, w*G2, kind(eps)))
      eps = eps + f3*wp**2/((w3**2 - w**2) - cmplx(0, w*G3, kind(eps)))
      eps = eps + f4*wp**2/((w4**2 - w**2) - cmplx(0, w*G4, kind(eps)))
      !
   end function epsPt
   !
   function epsPd(wavelength) result(eps)
      !
      ! Returns the wavelength-dependent relative dielectric function
      ! of a Lorentz-Drude metal, with the parameters for Pd from:
      ! Rakic et al. 1998, https://doi.org/10.1364/AO.37.005271
      !
      real(8), intent(in) :: wavelength ! in nm
      real(8), parameter :: &
         wp = 9.72d0, & ! eV
         f0 = 0.330d0, &
         G0 = 0.008d0   ! eV
      real(8), parameter :: &
         f1 = 0.649d0, &
         G1 = 2.950d0, & ! eV
         w1 = 0.336d0    ! eV
      real(8), parameter :: &
         f2 = 0.121d0, &
         G2 = 0.555d0, & ! eV
         w2 = 0.501d0    ! eV
      real(8), parameter :: &
         f3 = 0.638d0, &
         G3 = 4.621d0, & ! eV
         w3 = 1.659d0    ! eV
      real(8), parameter :: &
         f4 = 0.453d0, &
         G4 = 3.236d0, & ! eV
         w4 = 5.715d0    ! eV
      real(8) :: Op, w     ! omega
      complex(8) :: eps
      !
      Op = sqrt(f0)*wp ! eV
      w = hc/wavelength ! convert lambda to omega
      !
      eps = 1 - Op**2/(w*(w + cmplx(0, G0, kind(eps))))
      eps = eps + f1*wp**2/((w1**2 - w**2) - cmplx(0, w*G1, kind(eps)))
      eps = eps + f2*wp**2/((w2**2 - w**2) - cmplx(0, w*G2, kind(eps)))
      eps = eps + f3*wp**2/((w3**2 - w**2) - cmplx(0, w*G3, kind(eps)))
      eps = eps + f4*wp**2/((w4**2 - w**2) - cmplx(0, w*G4, kind(eps)))
      !
   end function epsPd
   !
   function epsDrude(wavelength, eps_infty, lambda_p, mu_p) result(eps)
      !
      ! Returns the wavelength-dependent relative dielectric function of a Drude metal
      ! The analytical expression is given in Eq. (3.2).
      ! The exp(-i omega t) convention is assumed.
      !
      ! Input parameters:
      ! - wavelength: real scalar,
      !               wavelength in NANOMETERS (nm).
      ! - eps_infty: scalar
      !              eps_infty (see Eq. 3.2).
      ! - lambda_p: scalar
      !             wavelength in NANOMETERS (nm) corresponding to
      !             plasma frequency omega_p [rad/s]: lambda_p=2*pi*c/(omega_p)
      ! - mu_p: scalar
      !         wavelength in NANOMETERS (nm) corresponding to
      !         damping term gamma_0 [rad/s]: mu_p=2*pi*c/(gamma_0)
      !
      ! Returns:
      ! - eps: complex scalar,
      !        epsilon(wavelength) as a complex number.
      !
      real(8), intent(in) :: wavelength, eps_infty, lambda_p, mu_p
      complex(8) :: eps
      !
      eps = eps_infty*(1 - 1/(lambda_p**2*((1/wavelength)**2 + imu/(mu_p*wavelength))))
      !
      return
      !
   end function epsDrude

   function epsSi(wavelength) result(eps)
      ! epsSi
      ! Dielectric function of silicon in the UV-vis-nearIR region from tabulated
      ! values in palik's handbook )
      real(8), intent(in) :: wavelength(:)
      complex(8) :: eps(size(wavelength))
      complex(8) :: eps_(236) !
      real(8)    :: re_eps(236), im_eps(236)
      real(8), dimension(236) :: lambda
      integer :: i

      lambda = (/ &
           & 206.6, 207.3, 208.0, 208.7, 209.4, 210.1, 210.9, &
           & 211.6, 212.3, 213.0, 213.8, 214.5, 215.3, 216.0, &
           & 216.8, 217.5, 218.3, 219.1, 219.8, 220.6, 221.4, &
           & 222.2, 223.0, 223.8, 224.6, 225.4, 226.3, 227.1, &
           & 227.9, 228.8, 229.6, 230.5, 231.3, 232.2, 233.1, &
           & 233.9, 234.8, 235.7, 236.6, 237.5, 238.4, 239.4, &
           & 240.3, 241.2, 242.2, 243.1, 244.1, 245.0, 246.0, &
           & 247.0, 248.0, 249.0, 250.0, 251.0, 252.0, 253.0, &
           & 254.1, 255.1, 256.2, 257.2, 258.3, 259.4, 260.5, &
           & 261.6, 262.7, 263.8, 264.9, 266.1, 267.2, 268.4, &
           & 269.5, 270.7, 271.9, 273.1, 274.3, 275.5, 276.8, &
           & 278.0, 279.2, 280.5, 281.8, 283.1, 284.4, 285.7, &
           & 287.0, 288.3, 289.7, 291.0, 292.4, 293.8, 295.2, &
           & 296.6, 298.0, 299.5, 300.9, 302.4, 303.9, 305.4, &
           & 306.9, 308.4, 310.0, 311.5, 313.1, 314.7, 316.3, &
           & 317.9, 319.5, 321.2, 322.9, 324.6, 326.3, 328.0, &
           & 329.7, 331.5, 333.3, 335.1, 336.9, 338.8, 340.6, &
           & 342.5, 344.4, 346.3, 348.3, 350.2, 352.2, 354.2, &
           & 356.3, 358.3, 360.4, 362.5, 364.7, 366.8, 369.0, &
           & 371.2, 373.4, 375.7, 378.0, 380.3, 382.7, 385.0, &
           & 387.5, 389.9, 392.4, 394.9, 397.4, 400.0, 402.5, &
           & 405.2, 407.8, 410.5, 413.3, 416.1, 418.9, 421.7, &
           & 424.6, 427.5, 430.5, 433.5, 436.6, 439.7, 442.8, &
           & 446.0, 449.2, 452.5, 455.8, 459.2, 462.6, 466.1, &
           & 469.6, 473.2, 476.9, 480.6, 484.3, 488.1, 492.0, &
           & 495.9, 499.9, 504.0, 508.1, 512.3, 516.6, 520.9, &
           & 525.4, 529.9, 534.4, 539.1, 543.8, 548.6, 553.5, &
           & 558.5, 563.6, 568.7, 574.0, 579.4, 584.8, 590.4, &
           & 596.1, 601.9, 607.8, 613.8, 619.9, 626.2, 632.6, &
           & 639.1, 645.8, 652.6, 659.5, 666.6, 673.8, 681.2, &
           & 688.8, 696.5, 704.5, 712.6, 720.8, 729.3, 738.0, &
           & 746.9, 756.0, 765.3, 774.9, 784.7, 794.8, 805.1, &
           & 815.7, 826.6, 1120.0, 1144.0, 1200.0, 1372.0, 1400.0, &
           & 1532.0, 1600.0, 1696.0, 1800.0, 2000.0/)

      re_eps = (/ &
          & -7.44218, -7.49989, -7.57302, -7.48961, -7.63447, -7.71944, -7.73843, &
          & -7.81562, -7.86034, -7.89846, -7.98834, -8.0724, -8.1093, -8.1685, &
          & -8.24178, -8.29214, -8.4002, -8.45598, -8.54928, -8.65087, -8.72343, &
          & -8.81976, -8.89162, -8.98589, -9.05146, -9.1076, -9.16072, -9.18623, &
          & -9.21744, -9.19486, -9.16612, -9.08742, -9.01475, -8.92019, -8.82381, &
          & -8.74937, -8.68349, -8.65424, -8.65105, -8.66794, -8.7244, -8.79437, &
          & -8.89465, -8.99208, -9.14014, -9.29, -9.4445, -9.62777, -9.81939, &
          & -10.0203, -10.2443, -10.465, -10.695, -10.9598, -11.2254, -11.5046, &
          & -11.7709, -12.0893, -12.4008, -12.7297, -13.0835, -13.45, -13.8489, &
          & -14.2798, -14.7154, -15.1896, -15.7041, -16.2369, -16.7922, -17.3548, &
          & -17.9315, -18.4562, -18.9332, -19.3437, -19.6392, -19.8193, -19.8867, &
          & -19.8124, -19.6102, -19.2877, -18.8239, -18.2164, -17.4503, -16.3333, &
          & -14.7679, -12.4106, -9.45548, -6.12934, -2.93152, -0.067319, 2.37222, &
          & 4.34308, 5.9783, 7.31638, 8.42475, 9.36796, 10.1172, 10.7756, &
          & 11.3419, 11.8277, 12.2335, 12.6362, 13.0106, 13.332, 13.6551, &
          & 13.9645, 14.2543, 14.5652, 14.891, 15.2122, 15.5291, 15.8485, &
          & 16.1746, 16.5346, 16.8777, 17.233, 17.5865, 17.9552, 18.3213, &
          & 18.7081, 19.1254, 19.5746, 20.0724, 20.6812, 21.4212, 22.3879, &
          & 23.7106, 25.6007, 28.1836, 31.4907, 35.2195, 38.7911, 41.4811, &
          & 43.1383, 43.7395, 43.2656, 42.1301, 40.7353, 39.2276, 37.7444, &
          & 36.3509, 35.064, 33.8739, 32.7856, 31.7947, 30.8751, 30.047, &
          & 29.2682, 28.5138, 27.836, 27.1969, 26.6019, 26.0423, 25.5314, &
          & 25.0456, 24.5703, 24.1294, 23.7022, 23.3043, 22.9248, 22.5644, &
          & 22.2373, 21.8989, 21.5862, 21.2811, 20.987, 20.7126, 20.4305, &
          & 20.1906, 19.9308, 19.7233, 19.4922, 19.274, 19.0644, 18.8557, &
          & 18.6571, 18.4675, 18.2884, 18.0998, 17.9316, 17.7626, 17.6033, &
          & 17.4445, 17.2954, 17.1376, 16.9968, 16.8573, 16.718, 16.5883, &
          & 16.4578, 16.3367, 16.2075, 16.0952, 15.9753, 15.8634, 15.7521, &
          & 15.6492, 15.5466, 15.4521, 15.3501, 15.2564, 15.155, 15.0696, &
          & 14.9766, 14.8839, 14.7992, 14.7223, 14.6381, 14.554, 14.4779, &
          & 14.4094, 14.3412, 14.2731, 14.1977, 14.145, 14.0774, 14.0249, &
          & 13.9576, 13.8979, 13.8458, 13.7937, 13.727, 13.6678, 13.6013, &
          & 13.5497, 13.4909, 12.504, 12.4574, 12.3855, 12.2549, 12.1634, &
          & 12.0993, 12.0478, 12.0021, 11.9564, 11.8956/)

      im_eps = (/ &
        & 5.87618, 6.06682, 6.15885, 6.26168, 6.34082, 6.45901, 6.49971, &
          & 6.62302, 6.68727, 6.76995, 6.89997, 6.97296, 7.09863, 7.18421, &
          & 7.2897, 7.34432, 7.49265, 7.6293, 7.74504, 7.8793, 7.99576, &
          & 8.16684, 8.3072, 8.48767, 8.66583, 8.84936, 9.04096, 9.26185, &
          & 9.4872, 9.70751, 9.90277, 10.1175, 10.2791, 10.4149, 10.5173, &
          & 10.5888, 10.6068, 10.6339, 10.6406, 10.6568, 10.659, 10.6674, &
          & 10.6817, 10.7025, 10.7216, 10.7739, 10.8361, 10.9108, 10.9885, &
          & 11.0866, 11.1941, 11.3337, 11.4771, 11.6266, 11.802, 11.9743, &
          & 12.1854, 12.4101, 12.641, 12.9074, 13.1944, 13.4877, 13.8338, &
          & 14.2145, 14.629, 15.0928, 15.6078, 16.208, 16.8885, 17.6899, &
          & 18.5997, 19.6182, 20.7537, 22.0406, 23.4415, 24.912, 26.4813, &
          & 28.1124, 29.7862, 31.5482, 33.3466, 35.2671, 37.2916, 39.5016, &
          & 41.8694, 44.0879, 45.7794, 46.6799, 46.7569, 46.2433, 45.3509, &
          & 44.2714, 43.1506, 42.0316, 40.9581, 39.9492, 39.0132, 38.1319, &
          & 37.3339, 36.6095, 35.9417, 35.3535, 34.8395, 34.3517, 33.9413, &
          & 33.5671, 33.239, 32.9399, 32.6794, 32.4439, 32.2337, 32.0383, &
          & 31.8743, 31.7437, 31.6152, 31.5341, 31.478, 31.4426, 31.4592, &
          & 31.5134, 31.6383, 31.8346, 32.1257, 32.5323, 33.079, 33.8171, &
          & 34.6961, 35.6351, 36.3148, 36.3467, 35.284, 32.8858, 29.481, &
          & 25.5405, 21.444, 17.7252, 14.6187, 12.1943, 10.2951, 8.83218, &
          & 7.63812, 6.67366, 5.90042, 5.23853, 4.70413, 4.31118, 3.90003, &
          & 3.56636, 3.34847, 3.07529, 2.80944, 2.63364, 2.49319, 2.30645, &
          & 2.1138, 2.01417, 1.90741, 1.80264, 1.78747, 1.62894, 1.54948, &
          & 1.40596, 1.39524, 1.23637, 1.20913, 1.19158, 1.19289, 1.2119, &
          & 1.0788, 1.07184, 0.79956, 0.830208, 0.728906, 0.689986, 0.668822, &
          & 0.63072, 0.627508, 0.564564, 0.61272, 0.5082, 0.5058, 0.469952, &
          & 0.442762, 0.357674, 0.3726, 0.395808, 0.361328, 0.359832, 0.260672, &
          & 0.308332, 0.258688, 0.273768, 0.24072, 0.215838, 0.23898, 0.23814, &
          & 0.213624, 0.19715, 0.19655, 0.188064, 0.171864, 0.171292, 0.147516, &
          & 0.13932, 0.131172, 0.123104, 0.122784, 0.11478, 0.10682, 0.09893, &
          & 0.098696, 0.098462, 0.090672, 0.082896, 0.082742, 0.07504, 0.0749, &
          & 0.067248, 0.067104, 0.059536, 0.059424, 0.05187, 0.051758, 0.044256, &
          & 0.044172, 0.03673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

      do i = 1, 236
         eps_(i) = complex(re_eps(i), im_eps(i))
      end do

      call interp1(lambda, eps_, wavelength, eps)
   end function epsSi

   function epsAl(wavelength) result(eps)
!--------------------------------------------------------
      !eps_Aluminum. ref:
      !Algorithm for the determination of
      !intrinsic optical constants of
      !metal films: application to aluminum, Aleksandar D. Rakić,
      !Applied Optics Vol. 34, Issue 22, pp. 4755-4767 (1995)
!----------------------------------------------------------
      real(8), intent(in) :: wavelength(:)
      complex(8) :: eps(size(wavelength))
      complex(8) :: eps_(38)
      real(8)    :: re_eps(38), im_eps(38)
      real(8), dimension(38) :: lambda
      integer :: i

      lambda = (/103.32, 112.71, 123.99, 137.76, 154.98, 177.12, 206.64, &
      247.97, 309.96, 326.28, 364.66, 413.28, 442.8, 476.87, 516.6, 563.57, &
      & 619.93, 652.25, 688.81, 729.32, 774.91, 794.78, 815.69, 837.74, &
      & 885.61, 911.66, 939.28, 968.63, 999.88, 1033.2, 1127.1, 1239.9, &
      & 1377.6, 1549.8, 1771.2, 2066.4, 2479.7, 2755.2/)

      re_eps = (/-0.594135, -0.913929, -1.33304, -1.89424, -2.6732, -3.80103, &
      & -5.53608, -8.39393, -13.6716, -15.2399, -19.1632, -24.7362, -28.4416, &
      & -32.8567, -38.1972, -44.8704, -52.971, -57.364, -61.6159, -65.1629, &
      &-65.2636, -62.676, -60.5278, -59.7797, -60.6279, -65.09, -71.0914, &
      & -79.973, -88.0723, -95.8804, -118.555, -148.209, -188.067, -242.62, &
      & -319.988, -434.125, -614.063, -745.149/)

      im_eps = (/0.0551762, 0.0736101, 0.107009, 0.157495, 0.237323, 0.367878, &
       & 0.597416, 1.053, 2.07676, 2.46536, 3.50575, 5.21433, 6.52593, 8.4106, &
       & 10.9032, 14.5555, 20.231, 24.3263, 29.5016, 36.1101, 44.4168, 46.4198, &
       & 45.6926, 44.1241, 37.0003, 32.7896, 28.8584, 26.9554, 27.2646, 27.692, &
       & 29.1359, 32.2215, 38.3168, 49.4166, 69.1034, 103.811, 166.887, 217.22/)

      do i = 1, 38
         eps_(i) = complex(re_eps(i), im_eps(i))
      end do

      call interp1(lambda, eps_, wavelength, eps)
   end function epsAl
!---------------------------------------------------------------------
   function epsCr(wavelength) result(eps)
!--------------------------------------------------------
      !eps_Chromium. ref:
      !Handbook of Palik, page :382-385
!----------------------------------------------------------
      real(8), intent(in) :: wavelength(:)
      complex(8) :: eps(size(wavelength))
      complex(8) :: eps_(125)
      real(8)    :: re_eps(125), im_eps(125)
      real(8), dimension(125) :: lambda
      integer :: i

      lambda = (/100.8, 103.3, 106.0, 107.8, 109.7, 111.7, 114.8, 118.1, 120.4, 124.0, &
       & 127.0, 130.0, 134.0, 137.0, 139.0, 142.0, 144.0, 146.0, 149.0, 151.0, 153.1, &
       & 155.0, 156.9, 159.0, 161.0, 162.9, 164.9, 168.0, 170.1, 172.0, 173.9, 176.9, 179.9, &
       & 182.1, 185.1, 187.9, 191.0, 194.0, 197.1, 200.0, 202.9, 207.0, 210.1, 214.1, &
       & 217.1, 221.0, 225.0, 230.0, 233.9, 238.0, 243.1, 248.0, 253.0, 258.3, 263.8, 269.5, &
       & 275.5, 281.8, 287.0, 293.8, 300.2, 307.7, 316.3, 323.7, 333.3, 341.6, 351.2, &
       & 362.5, 372.3, 385.0, 396.1, 409.2, 424.6, 438.1, 455.8, 471.4, 490.1, 512.3, &
       & 532.1, 558.5, 582.1, 610.8, 700.5, 815.7, 826.6, 849.2, 861.0, 885.6, 911.6, &
       & 939.3, 968.6, 999.9, 1033.0, 1069.0, 1107.0, 1148.0, 1192.0, 1240.0, 1292.0, &
       & 1348.0, 1378.0, 1442.0, 1512.0, 1590.0, 1675.0, 1771.0, 1879.0, 2000.0, &
       & 2138.0, 2296.0, 2480.0, 2695.0, 2952.0, 3263.0, 3647.0, 4133.0, 4769.0, 5636.0, &
       & 6888.0, 8856.0, 12400.0, 13780.0, 15500.0, 20660.0, 31000.0/)

      re_eps = (/0.7955, 0.7869, 0.7339, 0.712, 0.6688, 0.6264, 0.5709, 0.5017, &
        & 0.4816, 0.4275, 0.3696, 0.2988, 0.1467, 0.0328, -0.0328, -0.1485, -0.2145, &
        & -0.2988, -0.4225, -0.4959, -0.6055, -0.6825, -0.7697, -0.8736, -0.9672, &
        & -1.0773, -1.164, -1.2928, -1.3596, -1.428, -1.498, -1.5549, -1.6425, &
        & -1.7176, -1.8565, -1.8881, -1.9845, -2.0331, -2.048, -2.064, -2.112, &
        & -2.1093, -2.0713, -2.0867, -2.1216, -2.1251, -2.214, -2.3751, -2.5389, &
        & -2.74, -3.024, -3.3176, -3.5453, -3.7973, -4.1445, -4.424, -4.7304, &
        &-5.0464, -5.372, -5.7728, -6.1685, -6.5772, -6.9989, -7.4481, -7.8492, &
        & -8.1468, -8.3435, -8.5655, -8.9112, -9.4864, -10.3412, -11.3925, &
        & -12.4096, -13.2436, -13.8483, -14.0812, -13.5135, -12.3291, -10.9221, &
        & -9.3357, -8.0288, -6.8992, -4.3513, -0.9427, -0.516, -0.0863, 0.0865, &
        & 0.6083, 1.0464, 1.5768, 1.8417, 1.9316, 2.0263, 2.0309, 1.9448, 1.6853, &
        & 1.3305, 0.356, -0.2691, -0.9911, -1.5351, -2.7931, -5.1585, -8.244, &
        & -11.6064, -15.168, -19.7209, -23.736, -28.1424, -34.4772, -42.7875, &
        & -53.4016, -68.42, -94.4919, -135.648, -194.567, -277.122, -396.27, &
        & -568.947, -462.87, -748.8, -1009.97, -1595.88, -1314.56, -4029.03/)

      im_eps = (/1.6188, 1.582, 1.518, 1.5042, 1.4766, 1.449, 1.442, 1.4544, &
       & 1.44, 1.4308, 1.387, 1.3616, 1.3244, 1.3446, 1.3446, 1.3572, 1.3528, &
       & 1.3616, 1.3968, 1.42, 1.4352, 1.4552, 1.5096, 1.541, 1.5946, 1.6236, &
       & 1.7018, 1.8354, 1.904, 1.9738, 2.0448, 2.146, 2.2648, 2.265, 2.4492, &
       & 2.544, 2.6732, 2.822, 2.9568, 3.0082, 3.1648, 3.2524, 3.3216, 3.3756, &
       & 3.395, 3.306, 3.3088, 3.348, 3.294, 3.2538, 3.3368, 3.417, 3.5604, &
       & 3.6636, 3.8012, 4.0128, 4.183, 4.356, 4.5318, 4.8504, 5.2332, 5.6304, &
       & 6.042, 6.608, 7.1744, 7.8624, 8.4588, 9.0072, 9.4666, 9.792, 10.4784, &
       & 11.4268, 12.837, 14.616, 16.7956, 19.3584, 22.1112, 24.53, 26.522, &
       & 28.0476, 29.2584, 30.3456, 33.5616, 36.7164, 36.9782, 37.2384, 37.4112, &
       & 37.7556, 38.012, 38.3526, 38.4344, 38.52, 38.7816, 38.958, 39.0486, &
       & 39.3204, 39.3272, 39.6042, 40.23, 40.584, 40.756, 40.542, 40.7888, &
       & 41.5478, 43.036, 44.8318, 47.124, 50.6062, 53.333, 55.3504, 58.0628, &
       & 61.008, 62.2518, 63.036, 69.36, 90.09, 133.722, 218.36, 443.484, 807.84, &
       & 703.28, 800.04, 1147.84, 1780.8, 1942.96/)

      do i = 1, 125
         eps_(i) = complex(re_eps(i), im_eps(i))
      end do

      call interp1(lambda, eps_, wavelength, eps)

   end function epsCr
!--------------------------------------------------------------------------
   function epsWater(wavelength) result(eps)
!--------------------------------------------------------
!after Measurement of the refractive index of distilled water
! from the near-infrared region to the ultraviolet region
!Appl Opt. 2007 Jun 20;46(18):3811-20 (for temp=20)

!----------------------------------------------------------
      real(8), intent(in) :: wavelength(:)
      complex(8) :: eps(size(wavelength))
      complex(8) :: eps_(281)
      real(8)    :: re_eps(281)
      real(8), dimension(281) :: lambda
      integer :: i

      do i = 1, 281
         lambda(i) = 200 + 10*(i - 1)
      end do
      re_eps = (/2.02847, 1.98776, 1.95711, 1.9331, 1.91374, 1.89779, 1.88442, &
        & 1.87307, 1.86332, 1.85485, 1.84744, 1.84091, 1.83511, 1.82993, &
        & 1.82529, 1.8211, 1.8173, 1.81384, 1.81069, 1.8078, 1.80514, &
        & 1.80269, 1.80043, 1.79833, 1.79637, 1.79455, 1.79285, 1.79125, &
        & 1.78976, 1.78835, 1.78702, 1.78577, 1.78459, 1.78346, 1.78239, &
        & 1.78138, 1.78041, 1.77949, 1.77861, 1.77776, 1.77695, 1.77617, &
        & 1.77542, 1.7747, 1.77401, 1.77334, 1.77269, 1.77206, 1.77145, &
        & 1.77085, 1.77028, 1.76972, 1.76917, 1.76864, 1.76811, 1.7676, &
        & 1.76711, 1.76662, 1.76614, 1.76566, 1.7652, 1.76475, 1.7643, &
        & 1.76385, 1.76342, 1.76298, 1.76256, 1.76213, 1.76172, 1.7613, &
        & 1.76089, 1.76048, 1.76008, 1.75967, 1.75927, 1.75887, 1.75848, &
        & 1.75808, 1.75769, 1.75729, 1.7569, 1.75651, 1.75612, 1.75572, &
        & 1.75533, 1.75494, 1.75455, 1.75415, 1.75376, 1.75336, 1.75297, &
        & 1.75257, 1.75217, 1.75177, 1.75137, 1.75097, 1.75056, 1.75015, &
        & 1.74974, 1.74933, 1.74891, 1.7485, 1.74807, 1.74765, 1.74723, &
        & 1.7468, 1.74636, 1.74593, 1.74549, 1.74504, 1.7446, 1.74415, &
        & 1.74369, 1.74324, 1.74277, 1.74231, 1.74183, 1.74136, 1.74088, &
        & 1.74039, 1.7399, 1.73941, 1.73891, 1.7384, 1.73789, 1.73738, &
        & 1.73686, 1.73633, 1.7358, 1.73526, 1.73471, 1.73416, 1.7336, &
        & 1.73304, 1.73247, 1.73189, 1.73131, 1.73072, 1.73012, 1.72951, &
        & 1.7289, 1.72828, 1.72765, 1.72701, 1.72637, 1.72572, 1.72505, &
        & 1.72438, 1.7237, 1.72302, 1.72232, 1.72161, 1.7209, 1.72017, &
        & 1.71944, 1.71869, 1.71793, 1.71717, 1.71639, 1.7156, 1.7148, &
        & 1.71399, 1.71317, 1.71233, 1.71148, 1.71062, 1.70975, 1.70886, &
        & 1.70796, 1.70705, 1.70612, 1.70518, 1.70423, 1.70325, 1.70227, &
        & 1.70127, 1.70025, 1.69921, 1.69816, 1.69709, 1.69601, 1.6949, &
        & 1.69378, 1.69264, 1.69148, 1.6903, 1.6891, 1.68787, 1.68663, &
        & 1.68536, 1.68408, 1.68277, 1.68143, 1.68007, 1.67869, 1.67728, &
        & 1.67584, 1.67438, 1.67289, 1.67137, 1.66982, 1.66824, 1.66663, &
        & 1.66499, 1.66332, 1.66161, 1.65987, 1.65809, 1.65627, 1.65442, &
        & 1.65252, 1.65059, 1.64861, 1.64659, 1.64453, 1.64242, 1.64026, &
        & 1.63806, 1.6358, 1.63349, 1.63113, 1.62871, 1.62623, 1.6237, &
        & 1.6211, 1.61843, 1.6157, 1.6129, 1.61003, 1.60708, 1.60405, &
        & 1.60095, 1.59775, 1.59448, 1.59111, 1.58764, 1.58408, 1.58041, &
        & 1.57664, 1.57275, 1.56875, 1.56462, 1.56037, 1.55598, 1.55145, &
        & 1.54677, 1.54194, 1.53695, 1.53178, 1.52644, 1.52091, 1.51517, &
        & 1.50923, 1.50307, 1.49667, 1.49002, 1.48312, 1.47593, 1.46845, &
        & 1.46066, 1.45253, 1.44405, 1.43519, 1.42593, 1.41623, 1.40607, &
        & 1.39542, 1.38423, 1.37246, 1.36008, 1.34702, 1.33324, 1.31867, &
        & 1.30325, 1.28688, 1.2695, 1.25099, 1.23125, 1.21014, 1.18754, 1.16325/)
      do i = 1, 281
         eps_(i) = complex(re_eps(i), 0.0)
      end do

      call interp1(lambda, eps_, wavelength, eps)

   end function epsWater
!-----------------------------------------------------------------------

   subroutine interp1(x1, y1, x2, y2)
      ! Inputs:
      ! x1 = vector of x-values for the data to be interpolated
      ! y1 = vector of y-values for the data to be interpolated
      ! x2  = vector of x-values for the interpolated data
      ! Output:
      ! y2  = vector of y-values for the interpolated data
      !
      implicit none
      !
      real(8), intent(in) :: x1(:), x2(:)
      complex(8), intent(in) :: y1(size(x1))
      complex(8), intent(out) :: y2(size(x2))
      character(*), parameter :: myname = 'interp1'
      integer :: i1, i2, j
      real(8) ::  w
      !
      if (minval(x1) > minval(x2) .or. maxval(x1) < maxval(x2)) then
         write (*, '(A,A)') myname, &
            '> ERROR: x2-range not inside x1-range'
         STOP
      end if
      !
      do i1 = 2, size(x1)
         if (x1(i1) <= x1(i1 - 1)) then
            write (*, '(A,A)') myname, &
               '> ERROR: x1 values not monotonically increasing'
            STOP
         end if
      end do
      !
      do i2 = 2, size(x2)
         if (x2(i2) <= x2(i2 - 1)) then
            write (*, '(A,A)') myname, &
               '> ERROR: x2 values not monotonically increasing'
            STOP
         end if
      end do
      !
      j = 1
      do i2 = 1, size(x2)
         do i1 = j, size(x1)
            if (x1(i1) > x2(i2)) exit
         end do
         w = (x2(i2) - x1(i1 - 1))/(x1(i1) - x1(i1 - 1))
         y2(i2) = (real(1, kind(w)) - w)*y1(i1 - 1) + w*y1(i1)
         j = i1 ! to avoid unnecessary re-scanning
      end do
      !
   end subroutine interp1
   !
end module eps