PRO NGC5806_GMC_output, text1 = text1, text2 = text2, fn = fn

text1=strtrim(string(text1),2) ; hithreshold
text2=strtrim(string(text2),2) ; minconvexity
fn = strtrim(string(fn),2)
print, text1
print, text2
print, fn
infile='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_correct.fits'
cpropstoo_gmc_figure, infile='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_correct.fits', $
             inmask='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_mask.fits', $
             inpmom='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_moments_clfriendtoo.idl', $
             inprop='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_props_clfriendtoo.idl', $
             inassign='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_assign_clfriendtoo.fits', $
             outfile='./output' + fn + '/NGC5806_gmc_figure_clfriendtoo_new_hi'+text1+'_conv'+text2+'_all.eps', $
             ;/plcloudNr, $
			 ;plellradius = [300,600], $
			 ingalpar='./data/NGC5806_galpar.dat', $
             impad = [-50,-50, -50, -50], $
             ;impad = [-200,-200,-200,-200], $
		     scale_label = '200 pc'

print, infile

cpropstoo_gmc_figure, infile='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_correct.fits', $
             inmask='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_mask.fits', $
             inpmom='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_moments_clfriend.idl', $
             inprop='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_props_clfriend.idl', $
             inassign='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_assign_clfriend.fits', $
             outfile='./output' + fn + '/NGC5806_gmc_figure_clfriend_new_hi'+text1+'_conv'+text2+'_all.eps', $
             ;/plcloudNr, $
			 ;plellradius = [300,600], $
			 ingalpar='./data/NGC5806_galpar.dat', $
             impad = [-50,-50, -50, -50], $
             ;impad = [-200,-200,-200,-200], $
		     scale_label = '200 pc'

cpropstoo_gmc_figure, infile='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_correct.fits', $
             inmask='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_mask.fits', $
             inpmom='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_moments_cpropstoo.idl', $
             inprop='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_props_cpropstoo.idl', $
             inassign='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_assign_cpropstoo.fits', $
             outfile='./output' + fn + '/NGC5806_gmc_figure_cpropstoo_new_hi'+text1+'_conv'+text2+'_all.eps', $
             ;/plcloudNr, $
			 ;plellradius = [300,600], $
			 ingalpar='./data/NGC5806_galpar.dat', $
             impad = [-50,-50, -50, -50], $
             ;impad = [-200,-200,-200,-200], $
		     scale_label = '200 pc'

cpropstoo_gmc_figure, infile='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_correct.fits', $
             inmask='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_mask.fits', $
             inpmom='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_moments_cprops.idl', $
             inprop='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_props_cprops.idl', $
             inassign='./measurements' + fn + '/5806.CO21_rob0.5_b13_sub_assign_cprops.fits', $
             outfile='./output' + fn + '/NGC5806_gmc_figure_cprops_new_hi'+text1+'_conv'+text2+'_all.eps', $
             ;/plcloudNr, $
			 ;plellradius = [300,600], $
			 ingalpar='./data/NGC5806_galpar.dat', $
             impad = [-50,-50, -50, -50], $
             ;impad = [-200,-200,-200,-200], $
		     scale_label = '200 pc'

END
