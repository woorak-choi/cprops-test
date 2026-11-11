pro cap_diagonal, val1, VAL=val, _EXTRA=extra
;
; Plots a line parallel to the one-to-one line in existing plot.
; In logarithmic coordinates the input values represent
; fractional differences from the diagonal line.
;
; V1.0.0: Harald Kuntchner, secret safe location, 2001 (?)
; V1.1.0: Use Coyote Graphics. Michele Cappellari, Oxford, 9 October 2011
; V1.2.0: Allow input values to be give either as parameter or
;   as keyword to keep the procedure backward compatible.
;   Removed unnecessary /LOG keyword. MC, Oxford 9 November 2011
; V2.0.0: Adapted for diagonal line. MC, Oxford 9 November 2011
; V2.1.0: Works for logarithmic plots. MC, 27 October 2012
; V2.1.1: Make it compatible with IDL < 8. MC, Oxford, 2 February 2017

if n_elements(linestyle) eq 0 then linestyle = 2
if n_elements(val) eq 0 then val = N_PARAMS() ? val1 : !x.type

xx = !x.type ? 10^!x.crange : !x.crange
if !x.type then $
    for j=0,n_elements(val)-1 do $
        cgplot, xx, xx * val[[j,j]], _EXTRA=extra, LINESTYLE=linestyle, /OVERPLOT $
else $
    for j=0,n_elements(val)-1 do $
        cgplot, xx, xx + val[[j,j]], _EXTRA=extra, LINESTYLE=linestyle, /OVERPLOT

end
