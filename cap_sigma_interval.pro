;-----------------------------------------------------------------
function cap_sigma_interval, vec, int, CONF=p
COMPILE_OPT idl2
;
; Half of the interval between the 16th and 84th percentiles
; percentiles are given by 50 + [-1,1]*68.27/2
; 16% percentile (0.158655)
; 84% percentile (0.841345)
; http://en.wikipedia.org/wiki/Percentile
; Linear interpolation between closest ranks
; 
; V1.0: Michele Cappellari, Oxford, 31 March 2011
; V1.1: Use interpolation. MC, Oxford, 16 April 2011

if n_elements(vec) lt 3 then begin
    int = [0,0]
    return, 0
endif    

if n_elements(p) eq 0 then p = 0.6827d ; Default 1sigma
n = n_elements(vec)
v = vec[sort(vec)]
p = (1d + [-p,p])/2d
pn = (dindgen(n) + 0.5d)/n
int = interpol(v,pn,p)
sig = (int[1] - int[0])/2d

return, sig
end
;-----------------------------------------------------------------
