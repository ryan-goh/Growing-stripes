## Overview of core-far field codes for  stripe-forming fronts in quenched Swift-Hohenberg equation,
As explained in Sec. 6 of _Growing stripes, with and without wrinkles_
Directory overviews and brief explanations. Finds heteroclinic fronts in the elliptic modulated traveling wave problem:
c k_x u_y = -( 1+\partial_x^2 + k_y\partial y^2)^2 u + \mu(x) u + c u_x - u^3.

General outline and process to create data in Figures. See headers below for descriptions of each directory.

To understand skeleton of codes, go to /GeneralContinuationCodes/c-cont for continuation in speed c, and /GeneralContinuation/ky-cont for continuation in ky.

  To obtain data, run
  * /GeneralContinuationCodes/c-cont/moduli_cont_perp.m
  * /GeneralContinuationCodes/ky-cont/moduli_cont_perp.m
  * /GeneralContinuationCodes/ky-cont/moduli_cont.m
  * In /GeneralContinuationCodes/ky-cont call function in matlab:  cont(1,0,0.001,1.02,-1)
    - copy data file 'ContData0.mat' to /FullBubble, /Catastrophe, and /Zig-ZagBubble
  * /Zig-ZagBubble/moduli_cont_zz.m
  * /Catastrophe/moduli_cont_zz_full.m  -> can edit this to only run a few values near the catastrophe
  * /FullBubble/moduli_cont_zz_full.m
  * /FullBubble/Process.m  - Creates file 'moduli_dat_ZZBubble.mat', move this into  /GeneralContinuationCodes/c-cont/moduli_cont_perp.m
  * For data of isolas, see discussion below on /GeneralContinuationCodes/c-cont/BarbaPlot.m

After all data is created, see each individual directory below for detail of how to plot figures.  General key for figure files:
* Run /GeneralContinuationCodes/ky-cont/plot_sol.m  to obtain the figures
  - cprof.eps
  - cslice.eps
  - csmallprof.eps
  - csmallslice.eps
  - coreff.eps

* Run /GeneralContinuationCodes/ky-cont/plot_fig.m to obtain the figures
  - moduli_wire.eps
  - moduli_topsurf.eps
  - detachment.eps

* Run /GeneralContinuationCodes/c-cont/plot_fig.m  to obtain the figures
  - moduli_topdown.eps
  - moduli_topdownzz.eps
  - moduli_v1.eps
  - moduli_v2.eps
  - moduli_v3.eps
  - moduli.avi  -> movie giving a 360degree view of moduli space, cutoff at the perp. stripes domain (set makevideo = 1)
  - moduli.eps
  - perp_comp.eps



* Run /GeneralContinuationCodes/c-cont/plot_sol.m to obtain the figures
  - zzprof.eps
  - zzslice.eps
  - obprofiles.eps
  - obslice.eps
  - movie_ky*.avi (set makevideo = 1 and set kychoose to desired value)



* Run /GeneralContinuationCodes/c-cont/BarbaPlot.m to obtain the figures
  - isola8.eps
  - isolap8plain.eps
  - isola8profiles.eps
  - perp8.eps
  - barba_ky*.avi (set makevideo=1);

* Run /Zig-ZagBubble/plot_fig.m
  - shbubble1.eps
  - shbubble2.eps
  - shbubbleprofiles.eps
  - shbubbleslice.eps
  - Also various projections of the saddle node curve sn_cxkx.eps, sn_kycx.eps, sn_kykx.eps
* Run /Catastrophe/CatastrophePlot.m
  - afterprofiles.eps
  - afterslice.eps
  - beforeprofiles.eps
  - beforeslice.eps
  - Catastrophe_ky8.420000e-01.avi
  - Catastrophe_ky8.471200e-01.avi




### /GeneralContinuationCodes
The skeleton continuation codes, which continue fronts in either c or ky solving for the core perturbation w and wavenumber k_x see /GeneralCode

  #### /c-cont
  moduli_cont_perp.m - main continuation script, continues in c, for a range of fixed ky values, starting with perpendicular stripes, continuing past the fold and back down in c towards c = 0, Then continuation restarts near fold, to find branch of oblique stripes bifurcating and continue to all-stripe detachment curves. Saves data in moduli_dat_perp.mat, but also in .mat solution files, which also contain solution profiles.




  cont.m - Main contuation function, makes one run for fixed ky, takes argument (dtcont,cc,kx_init,ky_init,dir)

      dtcont - continuation stepsize, 0.5 to 1 seems to work well
      cc - initial speed c to start continuation at
      kx_init - initial kx value for initial Newton guess
      ky_init - initial ky value for initial Newton guess, fixed for whole run
      dir - direction of continuation in ky: 1 = up, -1=down, [1,-1] tells continuation to first go up in ky and then restart from (kx_init, ky_init) going down in ky

##### Function files
  shwf.m — builds nonlinear function for fsolve, also outputs certain derivatives of input

  shw.m — wrapper for shwf.m,  also calculates Jacobean: in w explicitly, in kx and ky with finite difference quotient

  shu.m — Builds nonlinear function and Jacobian for far-field solver getroll.m

  Lw.m — for input kx,ky, calculates the constant coefficient linear operator, mu>0

  inhom.m — includes the quenching inhomogeneity

  getroll.m — given a bulk wavenumber k, uses fsolve to find a 1-d roll solution

##### Spectral script
    spec_analysis.m - Loads in "core "solutions (w,kx,c) of a continuation run, contructs the full solution u = w+ \chi u_p, and then studies spectrum near the origin, finds saddle-node and pitchfork locations

##### Plotting scripts
  plot_fig.m -  Plots numerous figures and movies for paper.
  Loads data from
    - moduli_dat_perp.mat,
    - moduli_dat_ZZBubble.mat (from /FullBubble),
    -  moduli_dat_perp1.mat (was another set of continuation runs, made to fill in a few gaps for better surface interpolation, this shouldn't been needed now and associated  code can be commented out)
    -  pf_data.mat, sn_cont_4th_0.25.mat, sn_cont_small_ky.mat (From NWS continuation, codes not provided here, only data)

    In general, this script, takes individual continuation runs in c (i.e. sets of points in (ky,c,kx) in R^3, along a slice of the moduli space), interpolates them onto a uniform parameterized grid, then surf can be used to plot an interpolated surface over the data.

    moduli_topdown.eps
    moduli_topdownzz.eps
    moduli_v1.eps
    Fig. 15 moduli_v2.eps
    Fig. 15 moduli_v3.eps
    Footnote4 moduli.avi  -> movie giving a 360degree view of moduli space, cutoff
    at the perp. stripes domain (set makevideo = 1 )
    moduli.eps
    Figure 20 perp_comp.eps







  plot_sol.m - Plots a slice of moduli space with ky fixed, along with a few solution profiles along the curve.
    Fig. 19 	obprofiles.eps
	            obslice.eps
              zzprof.eps
              zzslice.eps
              movie_ky*.avi (set makevideo = 1 and set kychoose to desired value)



  sprof.m - input (x,y,kky,z), where x,y are the spatial mesh, kky is the fixed ky value, and z is a core solutions with bifurcatoin parameters (w,kx,c)

  BarbaPlot.m -  Script to plot isolas near perpendicular stripes.
        isola8.eps
      isolap8plain.eps
      isola8profiles.eps
      perp8.eps
      barba_ky*.avi (set makevideo=1)

      - To obtain the data for this plotting script run the additional four runs of the function cont.m -


        dtcont = 0.1;
        cc = 0.03;
        kx_init = 0;
        ky_init = [.7779,0.778,0.779,0.78];
        dir = 1;,
        cont(dtcont,cc,kx_init,ky_init,dir))
        also reduced the STEPS variable to 1000

        and one additional run with (c,kx_init) = (0.03068,0.02857) and ky_init = 0.78 to get the small hump of oblique stripes.










  #### /ky-cont
  moduli_cont_perp.m  - Continuation script to continue perpendicular stripes, saves data in moduli_dat_full.mat

  moduli_cont.m  - Continuation script to continue oblique stripes, saves data in moduli_dat_full_perp.mat



  cont.m  - Main continuation function takes argument (dtcont,cc,kx_init,ky_init,dir)

          dtcont - continuation step size, 1 or 2 seems to work well
		      cc  - quenching speed c, fixed throughout continuation run
		      kx_init - initial kx value for initial Newton guess
		      ky_init - initial ky value for initial Newton guess
		      dir - direction of continuation in ky: 1 = up, -1=down, [1,-1] tells continuation to first go up in ky and then restart from (kx_init, ky_init) going down in ky

  ##### Function files
    shwf.m — builds nonlinear function for fsolve, also outputs certain derivatives of input

    shw.m — wrapper for shwf.m,  also calculates Jacobean: in w explicitly, in kx and ky with finite difference quotient

    shu.m — Builds nonlinear function and Jacobian for far-field solver getroll.m

    Lw.m — for input kx,ky, calculates the constant coefficient linear operator, mu>0

    inhom.m — includes the quenching inhomogeneity

    getroll.m — given a bulk wavenumber k, uses fsolve to find a 1-d roll solution

    sprof.m - function which takes in (x,y,cc,PROF, pl), outputs a full core+far-field solution, in matrix form, for plotting using imagesc.  x, y - are spatial meshes, PROF gives a solution snapshot (w,kx,ky), cc is the speed.


 ##### Spectral script
  spec_analysis.m - Loads in "core "solutions (w,kx,ky) of a continuation run, contructs the full solution u = w+ \chi u_p, and then studies spectrum near the origin, finds saddle-node and pitchfork locations

  ##### Plotting scripts
    plot_fig.m  - Combines data from perp and oblique stripe continuation scripts, and plots them
      plots figures 	Fig. 15, top left moduli_wire.eps
	                    Fig. 15, top right moduli_topsurf.eps
	                    Fig. 20 right detachment.eps

        In general, sthis script takes individual continuation runs in ky (i.e. slices of the moduli space (ky,c,kx) in R^3), interpolates them onto a uniform parameterized grid, then surf can be used to plot an interpolated surface over the data.



    plot_sol.m - Script to plot a slice of moduli space, with solution profiles for 4 points along bifurcation curves
      plots from Figure 17: 	cprof.eps    with "cchoose" set to 2.948980e-01
	                   cslice.eps
	                   csmallprof.eps with "cchose" set to 9.081633e-02
	                    csmallslice.eps
	                Figure 16: coreff.eps
    isola.m  - plots a few specific continuation runs with c ~0.53, where the solution lies on an isola, bounded by two folds in c




### /Zig-ZagBubble
Takes a c = 0 continuation run from /GeneralContinuationCodes/ky-cont, then continues solutions in c, for various ky fixed

  cont.m - similar continuation function as in /c-cont
  moduli_cont_zz.m - main continuation script.
    To run:
        * first go to /GeneralContinuationCodes/ky-cont directory and run cont(1,0,0.001,1.02,-1) with SAVE = true (this catches the pitchfork and continues onto the zig-zag critical branch with kx\neq0 for c = 0.
        * Move the saved data file ContData0.mat into  /Zig-Zag Bubble Directory
        * then run the script moduli_cont_zz.m
    This starts a continuation runs in c for initial guesses on the c = 0 curve near the bifurcation point, at (ky,c,kx) = (kzz,0,0), tracing out the kink-dragging bubble.
  #### plotting scripts
      plot_fig.m -
          Yields the figures, some of which are depicted in Figure 18. Requires data set Bubbleplots.mat (included), from 1-D Cahn-Hilliard continuation (codes not yet included)

          shbubble1.eps - Figure 18 upper left
	        shbubble2.eps -
	        shbubbleprofiles.eps -Fig. 18 lower right
        	shbubbleslice.eps - Fig. 18 lower left
        	sn_cxkx.eps
        	sn_kycx.eps   - Fig. 18 upper right
        	sn_kykx.eps

      plot_sol.m    - not needed for figures in paper, just extra plotting scripts

### /Catastrophe

A few runs in c (ky fixed) to resolve the hyperbolic catastrophe near ky = 0.84
  moduli_cont_zz_full.m makes a series of c-continuations off of zig-zag critical curve, with   increased resolution near ky = 0.842

  CatastrophePlot.m  - plots moduli slices, and solution profiles before and after Catastrophe,
  with ky = 8.471268e-01; and ky = 8.419514e-01; Before running, must run /GeneralContinuationCodes/c-cont/moduli_cont_perp.m continuation script, and copy over corresponding oblique stripe reattachment curves, 'ContData%ky.mat' with %ky equal to values above.





### /FullBubble
Doesn’t  make any figures, creates data used in plotting code /GeneralContinuationCodes/c-cont/plot_fig.m (in particular the figures which include the "delicate arch" from the zig-zag bubble)
Similar code to /Zig-ZagBubble, must copy over ContData0.mat from /GeneralContinuationCodes/ky-cont before running

moduli_cont_zz_full.m - Continue in c off of the zig-zag critical curve at c = 0, for all values of ky\leq kzz

Process.m complies data from continuation script, and makes the data file moduli_dat_ZZBubble.mat  (contains continuation data and surface interpolation data) which is required for plotting 	in	/GeneralContinuationCodes/plot_fig.m.
  This script also plots some parts of the moduli space.

  It looks for data created from moduli_cont_zz_full.m, had to re-run a few continuations (note the variables 'rerun' and 'touchup'.)
