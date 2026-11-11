pro moments_to_props $
   , indata = indata $
   , infile = infile $
   , inbootfile = inbootfile $
   , hdrfile = hdrfile $
   , inhdr = inhdr $
   , dist = dist $
   , alpha = alpha $
   , rmstorad = rmstorad $
   , chantosig = chantosig $
   , vircoeff = vircoeff $
   , verbose = verbose $
   , text_file = text_file $
   , idl_file = idl_file $
   , fits_file = fits_file

  compile_opt idl2

;   20170213  LJ  add bootstrap uncertainties based on cprops_dist
;   20170223  LJ  add gaussian-fitted velocity dispersion which gets rid of
;                 large scale motions (more details see Utomo+ 2015)


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS & DEFINITIONS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; INPUT MOMENT DATA
  if n_elements(indata) eq 0 then begin

     if n_elements(infile) gt 0 then begin
        restore, infile
     endif else begin
        message, "Need input moments data or file", /info
     endelse

  endif  


; INPUT BOOT-MOMENT DATA
  if n_elements(inbootfile) gt 0 then begin
     restore, inbootfile   ; bootmomra[bootstrap, nclouds]
	 bootstrap=(size(bootmomra))[1]
     message, "Calculate uncertainties by bootstrap resampling", /info
  endif else begin
     message, "Do not calculate uncertainties", /info
  endelse

; INPUT HEADER FILE
  if n_elements(inhdr) eq 0 then begin

     if n_elements(hdrfile) gt 0 then begin
        hdr = headfits(hdrfile)
     endif else begin
        message, "Need input header data or file to extract header.", /info
     endelse

  endif  
  
; DISTANCE
  if n_elements(dist) eq 0 then begin
     dist = sxpar(hdr, 'DIST', count=dist_found)
     if dist_found eq 0 then begin
        message, "Need a distance to derive properties.", /info
        return
     endif
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; LOOP OVER MOMENTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
 
  sz = size(moments)

  n_mom = n_elements(moments)
  
  for i = 0, n_mom-1 do begin
       
     if keyword_set(verbose) then begin
        counter, i+1, n_mom, "Properties for moment "
     endif

;    Calculate 
     this_prop = calc_props_from_moments( $
                 moments[i] $
                 , empty_props = empty_props $
                 , hdr = hdr $
                 , astr = astr $
                 , vaxis = vaxis $
                 , dist = dist $
                 , alpha = alpha $
                 , bmaj = bmaj $
                 , bmin = bmin $
                 , bpa = bpa $
                 , rmstorad = rmstorad $
                 , chantosig = chantosig $
                 , vircoeff = vircoeff)

     if n_elements(props) eq 0 then begin

        if sz[0] eq 1 then begin
           props = replicate(empty_props, sz[1])
        endif

        if sz[0] eq 2 then begin
           props = replicate(empty_props, sz[1], sz[2])
        endif

     endif
	 
	 if n_elements(inbootfile) gt 0 then begin  
	   ; THE FACTOR TO SCALE THE ERRORS BY TO ACCOUNT FOR CORRELATION
	   ; IN DATA DUE TO THE OVERSAMPLED BEAM
	   ppbeam = !pi*(this_prop.beamfwhm_pix/2.)^2./alog(2.)
	   indfac = sqrt(ppbeam)

	   bootpropra=[]
	   for jj = 0, bootstrap -1 do begin
         bootprop = calc_props_from_moments( $
                     bootmomra[jj,i] $
                     , empty_props = empty_props $
                     , hdr = hdr $
                     , astr = astr $
                     , vaxis = vaxis $
                     , dist = dist $
                     , alpha = alpha $
                     , bmaj = bmaj $
                     , bmin = bmin $
                     , bpa = bpa $
                     , rmstorad = rmstorad $
                     , chantosig = chantosig $
                     , vircoeff = vircoeff)
		bootpropra=[bootpropra,bootprop]
		endfor		 

		this_prop.radell_deconv_uc= indfac*mad(bootpropra.radell_deconv) $
			/median(bootpropra.radell_deconv)
		this_prop.radrms_extrap_deconv_uc= indfac*mad(bootpropra.radrms_extrap_deconv) $
			/median(bootpropra.radrms_extrap_deconv)
		this_prop.vrms_extrap_deconv_uc= indfac*mad(bootpropra.vrms_extrap_deconv) $
			/median(bootpropra.vrms_extrap_deconv)
		res_ind = where((bootpropra.vrms_gauss_deconv) ne 0., res_ct)
		if res_ct gt 0 then $
    	  this_prop.vrms_gauss_deconv_uc= indfac* $
			  mad((bootpropra.vrms_gauss_deconv)[res_ind]) $
			  /median((bootpropra.vrms_gauss_deconv)[res_ind])
		res_ind = where((bootpropra.vrms_gauss_extrap_deconv) ne 0., res_ct)
		if res_ct gt 0 then $
    	  this_prop.vrms_gauss_extrap_deconv_uc= indfac* $
			  mad((bootpropra.vrms_gauss_extrap_deconv)[res_ind]) $
			  /median((bootpropra.vrms_gauss_extrap_deconv)[res_ind])
		this_prop.lum_extrap_uc= indfac*mad(bootpropra.lum_extrap) $
			/median(bootpropra.lum_extrap)
		this_prop.mass_extrap_uc= indfac*mad(bootpropra.mass_extrap) $
			/median(bootpropra.mass_extrap)
		this_prop.virmass_radell_deconv_uc= indfac*$
			    mad(bootpropra.virmass_radell_deconv) $
			    /median(bootpropra.virmass_radell_deconv)
		res_ind = where((bootpropra.virmass_extrap_deconv) ne 0., res_ct)
		if res_ct gt 0 then $
    	    this_prop.virmass_extrap_deconv_uc= indfac*$
    		    mad((bootpropra.virmass_extrap_deconv)[res_ind]) $
    		    /median((bootpropra.virmass_extrap_deconv)[res_ind])

	 endif	 ; end bootstrap


     
     props[i] = this_prop

  endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; OUTPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if n_elements(idl_file) gt 0 then begin
     if n_elements(levels) gt 0 then begin
        save, file=idl_file, props, levels, kernel_ind, hdr, merger_matrix
     endif else begin
        save, file=idl_file, props
     endelse
  endif

  if n_elements(text_file) gt 0 then begin
     ; write to csv is super easy 
     ; need to tweak this
     write_csv, text_file, props, header=tag_names(props)
  endif
  
  if n_elements(fits_file) gt 0 then begin 
     mwrfits, props, fits_file, /create
  endif 

end                             ; OF MOMENTS_TO_PROPS
