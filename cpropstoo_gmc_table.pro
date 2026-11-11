PRO cpropstoo_gmc_table, infile=infile, $
						inprop=inprop, $
						inpmom=inpmom, $
						ingalpar=ingalpar, $
	                    outtable=outtable

;+
; NAME: 
;   cpropstoo_GMC_table
;
; PURPOSE:
;   Generate a table for cloud properties based on output of cpropstoo code
;
; INPUTS:
;   infile       --  the '.fits' file for data cube
;   inprop      --  the '.idl' file for output of cpropstoo-'moments_to_props'
;   inpmom       --  the '.idl' file for output of cpropstoo-'cube_to_moments'
;   ingalpar     --  the dat file of galaxy parameters
;    
; OUTPUTS:
;   outtable   --  the '.csv' file of output table
;        column 1: ID -- ID of cloud
;        column 2: RA(2000) -- RA in [h:m:s]
;        column 3: Dec(2000) -- Dec in [d:m:s]
;        column 4: Vlsr -- line velocity center [km/s]
;        column 5: R  - the cloud physcial radius (after extrap and deconv) in [pc]
;        column 6: vrms -- the cloud sigma linewidth (after extrap 
;                  and deconv) [km/s]
;        column 7: L -- total luminosity of cloud in [10^5 K km/s pc^2]
;                  (after extrap)
;        column 8: Mlum -- total luminous molecular mass in [10^6 Msun]
;                  (after extrap)
;        column 9: S/N -- SNR at peak location of cloud
;        column 10: Tmax -- peak flux in [K]
;        column 11: Oshear -- galactic shear in [km/s/pc]
;        column 12: d -- distance of cloud to galaxy center in [pc]
; 
; HISTORY:
;   20170208  LJ   introduced


; HEAD MESSAGE
  print, '================================'
  print, 'CPROPSTOO: GENERATE GMC TABLE'
  print, '================================'

; READ IN THE IDLFILE TO VARIABLE PROPS/MOMENTS
; props - a list of "this_prop" for all clouds
  restore, inprop
  restore, inpmom
  ; numbers - number of clouds
  numbers=(size(props))[1]
  print, '==========================='
  print, 'Cloud Numbers:', numbers
  print, '==========================='

; EXTRACT GALAXY NAME AND ASSOCIATED PROPERTIES
  galaxy=(strsplit(infile,'/',/extract))[-1]
  if (strmatch(galaxy,'*_*',/fold_case) eq 1) then begin
    galaxy=(strsplit(galaxy,'_',/extract))[0]
  endif
  if (strmatch(galaxy,'*.*',/fold_case) eq 1) then begin
    galaxy=(strsplit(galaxy,'.',/extract))[0]
  endif
  print, 'Galaxy Name: ', galaxy

; AXIS_RATIO AND PA. PA is the position angle in [degree] of major axis
; of galaxy disk from the North to the East: PA = 0 corresponds to a galaxy
; whose major axis is oriented North-South
  pa=0.0
  inc = 0.0
  ra_c = ''
  dec_c = ''
  dist_pc = 0.0
  blank = ''
  vsys = 0.0

  openr,lun,ingalpar,/get_lun
  readf,lun,blank
  readf,lun,blank
  readf,lun,pa
  readf,lun,blank
  readf,lun,inc
  readf,lun,blank
  readf,lun,ra_c
  readf,lun,blank
  readf,lun,dec_c
  readf,lun,blank
  readf,lun,dist_pc
  readf,lun,blank
  readf,lun,vsys

  close,lun
  free_lun,lun


  if ((pa eq 0.) and (ra_c eq 0.) and (dec_c eq 0.)) then begin
      print, 'Error: Wrong galaxy name, or No records for the galaxy !'
	  return 
  endif
 
	  
; GET THE POSITION IN [PIXEL] OF GALAXY CENTER  
  cube=readfits(infile,hdr)
  noise=mad(cube,/finite)
  extast, hdr, astr

  ihr = float((strsplit(ra_c,'h',/extract))[0])
  imin = float((strsplit((strsplit(ra_c,'h',/extract))[1], 'm', /extract))[0])
  xsec = float((strsplit((strsplit(ra_c,'m',/extract))[1], 's', /extract))[0])
  ideg = float((strsplit(dec_c,'d',/extract))[0])
  imn = float((strsplit((strsplit(dec_c,'d',/extract))[1], 'm', /extract))[0])
  xsc = float((strsplit((strsplit(dec_c,'m',/extract))[1], 's', /extract))[0])
  ra_c = ihr * (360./24.) + imin * (360./(24.*60.)) + $
	  xsec * (360./(24. * 60. * 60))     ; in [degree]
  dec_c = abs(ideg) + imn * (1./60.) + xsc * (1./3600.)   ; in [degree]
  dec_c = dec_c * ideg/abs(ideg)

  ad2xy, ra_c, dec_c, astr, xc, yc
  ;print,'===================================='
  print,'Center of Galaxy in Pixel:', xc, yc
  ;print,'===================================='


  rdhd, hdr, s = hdstruct
  vaxis = hdstruct.v
  if strmatch(sxpar(hdr,'CUNIT3'),'*KM/S*') or $
      strmatch(sxpar(hdr,'CUNIT3'),'*km/s*') then begin
      vaxis = vaxis *1e3    ; in [km/s]
  endif
  chanwidth = vaxis[1]-vaxis[0]

; GET DEGPERPIX
  xy2ad, [0,1], [0,0], astr, ra, dec
  degperpix = sphdist(ra[0], dec[0], ra[1], dec[1], /deg)
  pcperpix= degperpix*!dtor*(props[0].dist_pc)


; VECTORIFY ORIGINAL DATA CUBE 
  vectorify, cube, x = x, y = y, v = v, $
	  t = t



; INITIATE ALL VARIABLE ARRAY
  ID=make_array(numbers,/integer)
  RA=make_array(numbers,/str)
  Dec=make_array(numbers,/str)
  Vlsr=make_array(numbers,/float)         ; [km/s]
  R=make_array(numbers,/float)
  d_R=make_array(numbers,/float)
  vrms=make_array(numbers,/float)
  d_vrms=make_array(numbers,/float)
  vturb=make_array(numbers,/float)
  d_vturb=make_array(numbers,/float)
  L=make_array(numbers,/float)
  d_L=make_array(numbers,/float)
  Mlum=make_array(numbers,/float)
  d_Mlum=make_array(numbers,/float)
  SNR=make_array(numbers,/float)
  Tmax=make_array(numbers,/float)
  Omega_shear=make_array(numbers,/float)
  Theta_shear=make_array(numbers,/float)
  shear_vel=make_array(numbers,/float)
  d=make_array(numbers,/float)

  ; resolves = props.resolved_extrap
  ; We consider a cloud resolved when its deconvolved diameter is at 
  ; least one beam wide (its radius is greater than or equal to half 
  ; of the beam size) and its deconvolved velocity width is at least 
  ; as wide as half of one (Hanning smoothed) velocity channel
  degperpix = sqrt(abs(sxpar(hdr, "cdelt1"))*abs(sxpar(hdr, "cdelt1")))
  bmaj_deg = sxpar(hdr, "BMAJ")
  bmin_deg = sxpar(hdr, "BMIN")
  beamfwhm_deg = sqrt(bmaj_deg*bmin_deg)
  beamfwhm_pc = beamfwhm_deg * !dtor * dist_pc    ; in pc
  beamfwhm_pix = beamfwhm_deg/degperpix
  pixperbeam = (beamfwhm_pix/2.0)^2*!pi/alog(2.)
  ;print,'beamfwhm_pc',beamfwhm_pc
  dv =  abs(sxpar(hdr, 'CDELT3'))  ; in 'km/s'
  dv2 = dv
  if galaxy eq 'NGC4429' then begin
	  beamfwhm_pc = 14. * sqrt(alog(2.))
	  dv = 2.
	  dv2 = 0.0
  endif
  resolves = (props.radrms_extrap_deconv ge beamfwhm_pc/2./sqrt(alog(2.)) $
	  and (props.vrms_extrap_deconv ge dv/2.) $
	  and (props.vrms_gauss_extrap_deconv ge dv2/2.) and $
	  (props.area ge pixperbeam))
  res_ind = where(resolves eq 1)

  planefit_chisq = props.velchisq
  print, '========================================================='
  print, 'mean(planefit_chisq) for resolved clouds:', $
	  mean(planefit_chisq[res_ind], /nan)
  theind = where(planefit_chisq[res_ind] le 1, thect)
  print, 'number or resolved clouds with chisq <=1:', thect
  print, '========================================================='

; ATTRIBUE VALUES TO EACH VARIABLE
  for i = 0, numbers-1 do begin
    this_prop=props[i]
	this_mom=moments[i]
    ; 1. Cloud ID
    ID[i]=this_mom.peaknum
    ; 2. RA and Dec
    this_ra=this_prop.xpos       ; [degree]
    this_dec=this_prop.ypos      ; [degree]
    radec, this_ra, this_dec, ihr, imin, xsec, ideg, imn, xsc 
    RA[i]=string(ihr,FORMAT='(I3)')+':'+ $
	   	string(imin,FORMAT='(I02)')+':'+ $
		string(xsec,FORMAT='(F06.4)')
    Dec[i]=string(ideg,FORMAT='(I3)')+':'+ $
	   	string(imn,FORMAT='(I02)')+':'+ $
		string(xsc,FORMAT='(F06.3)')
    ; 3. Vlsr - line velocity center in [km/s]
    Vlsr[i]=this_prop.vpos
    ; 4. R/d_R - the cloud physcial radius (after extrap and deconv) in [pc]
    resolved=resolves[i]
	if resolved eq 1  then begin
	  R[i]=this_prop.radrms_extrap_deconv
	  d_R[i]=this_prop.radrms_extrap_deconv_uc*R[i]
    endif else begin
	  R[i]=0.0
	  d_R[i]=0.0
    endelse  
	; 5. vrms - the cloud linewidth in sigma (after extrap and deconv) in
	;    [km/s]
	vrms[i]=this_prop.vrms_extrap_deconv
	d_vrms[i]=this_prop.vrms_extrap_deconv_uc*vrms[i]
	; 6. vturb - the cloud linewidth in sigma (after extrap and deconv) in
	;    [km/s]
	vturb[i]=this_prop.vrms_gauss_extrap_deconv
	; LJ 20190123, this is to consider the velocity gradient caused by
    ; turbulence (Burkert & Bodenheimer 2000)
    ; vturb[i] = sqrt(vturb[i]^2 + 0.08^2.*R[i]^2/5.)
	d_vturb[i]=this_prop.vrms_gauss_extrap_deconv_uc*vturb[i]
	; 7. L - total luminosity of cloud in [10^5 K km/s pc^2] after extrap
	L[i]=this_prop.lum_extrap/1e4
	d_L[i]=this_prop.lum_extrap_uc*L[i]
	; 8. Mlum - total molecular mass in [10^6 Msun] after extrap
	Mlum[i]=this_prop.mass_extrap/1e5
	d_Mlum[i]=this_prop.mass_extrap_uc*Mlum[i]
	;; 9. S/N - SNR at peak location of cloud
	;SNR[i]=(this_prop.maxval)/(this_prop.noise)
	; 9. Tmax - peak flux in [k]
	Tmax[i]=this_prop.maxval
	; 10. CaseI: Omega_shear - galactic shear in [km/s/pc]
	;  	Note omgea_shear is calculated by Equation (15) rather than 
	;	Equation (16) in Utomo+ 2015
	Omega_shear[i]=this_prop.velgrad
	; 11. angle of rotation axis
    coef = this_prop.planefit_coef
	Theta_shear[i]= atan(coef[1],-coef[2]) * !radeg
	; angle of rotation axis relative to North direction, turning positive
    ; into the direction of the RA
	Theta_shear[i] = (Theta_shear[i] + 270) mod 360
	; 12. d - distance from the galaxy center in [pc]
	x2=this_mom.mom1x        ; in [pixel]
	y2=this_mom.mom1y        ; in [pixel]
	d[i]=calc_distance(x1=xc, y1=yc, x2=x2, y2=y2, $
		    in_pa=pa, in_inc=inc, pcperpix=pcperpix)

  endfor
  
  print, "============================================="
  print, 'min(Omega_shear[res_ind]):',min(Omega_shear[res_ind], /nan)
  print, 'max(Omega_shear[res_ind]):',max(Omega_shear[res_ind], /nan)
  print, 'mean(Omega_shear[res_ind]):',mean(Omega_shear[res_ind], /nan)
  print, 'min(R[res_ind]), max(R[res_ind]):',min(R[res_ind]), max(R[res_ind])
  print, 'min(vrms[res_ind]), max(vrms[res_ind]):',min(vrms[res_ind]), max(vrms[res_ind])
  print, 'min(vturb[res_ind]), max(vturb[res_ind]):',min(vturb[res_ind]), max(vturb[res_ind])
  print, 'min(Mlum[res_ind]), max(Mlum[res_ind]):',min(Mlum[res_ind]), max(Mlum[res_ind])
  print, "============================================="

; Write out csv file
;  write_csv, outtable, ID, RA, Dec, Vlsr, R, vrms, L, Mlum, SNR, Tmax, Omega_shear, d
  if (n_elements(incircVfield) le 0) then begin
    openw, lun, outtable, /get_lun
    printf, lun, ';ID',',','RA(2000)',',','Dec(2000)',',','Vlsr',',','R',',',$
		'd_R', ',', 'vrms',',','d_vrms',',','vturb',',','d_vturb',',','L',',',$
		'd_L',',', 'Mlum',',','d_Mlum',',','Tmax',',','Omega',',','Theta',',','d', $
        format='(A3,a1,A11,a1,A12,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A9,a1,A9,a1,A9,a1,A9,a1,A7,a1,A9,a1,A9,a1,A7)'
    printf,lun,'; ',',','[h:m:s]',',','[d:m:s]',',','[km/s]',',','[pc]',',','[pc]',$
        ',','[km/s]',',','[km/s]',',','[km/s]',',','[km/s]',',','[10^4*]',',',$
		'[10^4*]',',','[10^5*]',',','[10^5*]',',','[K]',',','[km/s/pc]',',',$
		'[degree]',',','[pc]', $
        format='(A3,a1,A11,a1,A12,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A9,a1,A9,a1,A9,a1,A9,a1,A9,a1,A9,a1,A7,a1,A9,a1,A9,a1,A7)'
    for i = 0, numbers-1 do begin
      printf,lun,ID[i],',',RA[i],',',Dec[i],',',Vlsr[i],',',R[i],',',d_R[i],',',$
          vrms[i],',',d_vrms[i],',',vturb[i],',',d_vturb[i],',',L[i],',',d_L[i],$
		  ',', Mlum[i],',',d_Mlum[i],',',Tmax[i],',',Omega_shear[i],$
		  ',',Theta_shear[i],',',d[i], $
          format='(I4,a1,A11,a1,A12,a1,F7.1,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F9.3,a1,F9.3,a1,F9.3,a1,F9.3,a1,F7.1,a1,F9.2,a1,F9.1,a1,I7)'
    endfor
    close, lun
    free_lun, lun
  endif else begin
    openw, lun, outtable, /get_lun
    printf, lun, ';ID',',','RA(2000)',',','Dec(2000)',',','Vlsr',',','R',',',$
		'd_R', ',', 'vrms',',','d_vrms',',','vturb',',','d_vturb',',','L',',',$
		'd_L',',', 'Mlum',',','d_Mlum',',','Tmax',',','Vshear',',','d', $
        format='(A3,a1,A11,a1,A12,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A9,a1,A9,a1,A9,a1,A9,a1,A7,a1,A9,a1,A7)'
    printf,lun,'; ',',','[h:m:s]',',','[d:m:s]',',','[km/s]',',','[pc]',',','[pc]',$
        ',','[km/s]',',','[km/s]',',','[km/s]',',','[km/s]',',','[10^4*]',',',$
		'[10^4*]',',','[10^5*]',',','[10^5*]',',','[K]',',','[km/s]',',','[pc]', $
        format='(A3,a1,A11,a1,A12,a1,A7,a1,A7,a1,A7,a1,A7,a1,A7,a1,A9,a1,A9,a1,A9,a1,A9,a1,A9,a1,A9,a1,A7,a1,A9,a1,A7)'
    for i = 0, numbers-1 do begin
      printf,lun,ID[i],',',RA[i],',',Dec[i],',',Vlsr[i],',',R[i],',',d_R[i],',',$
          vrms[i],',',d_vrms[i],',',vturb[i],',',d_vturb[i],',',L[i],',',d_L[i],$
		  ',', Mlum[i],',',d_Mlum[i],',',Tmax[i],',',shear_vel[i],$
		  ',',d[i], $
          format='(I4,a1,A11,a1,A12,a1,F7.1,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F7.2,a1,F9.2,a1,F9.2,a1,F9.2,a1,F9.2,a1,F7.1,a1,F9.2,a1,I7)'
    endfor
    close, lun
    free_lun, lun
  endelse



END
